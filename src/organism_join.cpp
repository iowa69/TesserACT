#include "organism_join.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <unordered_map>

#include "graph.h"

namespace ts {

namespace {

// How far back from a contig end markers are collected. Twenty kilobases holds
// roughly forty markers at the sampling density, so a handful of markers the
// assembly happens to miss does not silence the query.
constexpr uint32_t kWindow = 20000;

// A chromosomal join needs this many panel genomes behind its strongest marker
// adjacency. Twenty out of hundreds is a low bar in absolute terms and a high
// one in kind: it cannot be met by a coincidence between two genomes.
constexpr uint32_t kMinPanelChr = 20;
// Plasmids are mosaic and each one is carried by far fewer isolates, so the
// same absolute bar would silence the table entirely. The share test below
// does the real work here.
constexpr uint32_t kMinPanelPls = 5;

// Share of the panel members carrying both markers that must support the join.
// A pair present in only a few genomes cannot claim unanimity on thin evidence.
constexpr double kMinFractionChr = 0.5;
// Held higher for plasmids precisely because the absolute bar is lower: a
// plasmid adjacency is worth acting on only when it is near-universal among
// the plasmids that carry both flanks.
constexpr double kMinFractionPls = 0.7;

// Independent marker pairs that must agree on the same gap. One pair agreeing
// with itself is not agreement.
constexpr size_t kMinAgreeingPairs = 3;

constexpr int32_t kGapTolerance = 500;    // gaps within this are the same cluster
constexpr int32_t kMaxGap = 20000;

struct MarkerAt {
    uint64_t oriented;   // (marker id << 1) | orientation as walked
    int32_t offset;      // bases from the port to the marker's start
};

struct Vote {
    int32_t gap;
    uint32_t support;
};

struct Hit { size_t contig; uint32_t pos; int orient; uint32_t id; };

// One pass of model joining, restricted to a single replicon class. Returns
// the number of joins made and rewrites `contigs` in place.
size_t joinPass(const OrganismModel& model, Replicon cls, std::vector<std::string>& contigs,
                std::vector<double>& covs, int k, OrganismJoinStats& st) {
    const size_t n = contigs.size();
    if (n < 2) return 0;

    const uint32_t minPanel = cls == Replicon::Chromosome ? kMinPanelChr : kMinPanelPls;
    const double minFraction = cls == Replicon::Chromosome ? kMinFractionChr : kMinFractionPls;

    // ---- locate model markers in the assembly ----------------------------
    // A marker occurring more than once across the assembly names a repeat
    // here even if it is single-copy in the panel, and cannot identify a
    // locus, so it is dropped.
    std::vector<Hit> hits;
    std::unordered_map<uint32_t, uint32_t> seen;
    for (size_t c = 0; c < n; ++c) {
        forEachMarkerKmer(contigs[c], [&](uint64_t km, uint32_t pos, int orient) {
            const uint32_t id = model.markerOf(km);
            if (id == UINT32_MAX) return;
            hits.push_back({c, pos, orient, id});
            ++seen[id];
        });
    }

    // ---- classify each contig by the evidence it carries -------------------
    // A contig whose markers are seen on panel plasmids and not on panel
    // chromosomes is a plasmid contig, and vice versa. Contigs that say
    // nothing either way take part in both passes.
    std::vector<uint32_t> chrVotes(n, 0), plsVotes(n, 0);
    for (const Hit& h : hits) {
        const uint32_t gc = model.markerGenomes(h.id, Replicon::Chromosome);
        const uint32_t gp = model.markerGenomes(h.id, Replicon::Plasmid);
        if (gc > 0 && gp == 0) ++chrVotes[h.contig];
        else if (gp > 0 && gc == 0) ++plsVotes[h.contig];
    }
    std::vector<char> participates(n, 1);
    for (size_t c = 0; c < n; ++c) {
        const uint32_t cv = chrVotes[c], pv = plsVotes[c];
        if (cv == 0 && pv == 0) continue;                 // undecided: both passes
        const bool chromosomal = cv >= pv * 3;
        const bool plasmidic = pv >= cv * 3;
        if (cls == Replicon::Chromosome) participates[c] = chromosomal || !plasmidic;
        else participates[c] = plasmidic || !chromosomal;
    }

    std::vector<Hit> unique;
    unique.reserve(hits.size());
    for (const Hit& h : hits) {
        if (seen[h.id] != 1) continue;
        if (!participates[h.contig]) continue;
        if (model.markerGenomes(h.id, cls) == 0) continue;   // nothing to say in this class
        unique.push_back(h);
    }
    st.markersFound += unique.size();
    if (unique.empty()) return 0;

    // ---- gather markers at each port --------------------------------------
    // Port 2c+1 leaves contig c at its right end walking forward; port 2c+0
    // leaves at its left end, which is the same as walking forward along the
    // reverse complement. Entry ports mirror this: entering at the left end
    // means the contig is used as written.
    const size_t nports = n * 2;
    std::vector<std::vector<MarkerAt>> exitM(nports), entryM(nports);
    for (const Hit& h : unique) {
        const int32_t L = static_cast<int32_t>(contigs[h.contig].size());
        const int32_t pos = static_cast<int32_t>(h.pos);
        const uint64_t fwd = (static_cast<uint64_t>(h.id) << 1) | static_cast<uint64_t>(h.orient);
        const uint64_t rev = (static_cast<uint64_t>(h.id) << 1) | static_cast<uint64_t>(1 - h.orient);

        const int32_t dRight = L - pos;
        if (dRight <= static_cast<int32_t>(kWindow)) exitM[h.contig * 2 + 1].push_back({fwd, dRight});
        // Leaving to the left: the walk runs on the reverse complement, where
        // this marker sits `pos + k` from that end and reads flipped.
        const int32_t dLeft = pos + kMarkerK;
        if (dLeft <= static_cast<int32_t>(kWindow)) exitM[h.contig * 2 + 0].push_back({rev, dLeft});

        if (pos <= static_cast<int32_t>(kWindow)) entryM[h.contig * 2 + 0].push_back({fwd, pos});
        // Entering at the right end means the contig arrives reverse
        // complemented, putting this marker `L - pos - k` from its new start.
        const int32_t dIn = L - pos - kMarkerK;
        if (dIn >= 0 && dIn <= static_cast<int32_t>(kWindow))
            entryM[h.contig * 2 + 1].push_back({rev, dIn});
    }

    struct EntryRef { uint64_t oriented; int32_t offset; uint32_t port; };
    std::vector<EntryRef> allEntries;
    for (uint32_t p = 0; p < nports; ++p) {
        for (const MarkerAt& m : entryM[p]) allEntries.push_back({m.oriented, m.offset, p});
    }

    // ---- score every exit port against every entry port -------------------
    struct Best { uint32_t port = UINT32_MAX; size_t pairs = 0; uint32_t support = 0; int32_t gap = 0; };
    std::vector<Best> best(nports);

    std::unordered_map<uint32_t, std::vector<Vote>> votesByPort;
    for (uint32_t p = 0; p < nports; ++p) {
        if (exitM[p].empty()) continue;
        votesByPort.clear();

        for (const MarkerAt& a : exitM[p]) {
            const uint32_t bothA = model.markerGenomes(static_cast<uint32_t>(a.oriented >> 1), cls);
            if (bothA == 0) continue;
            for (const EntryRef& b : allEntries) {
                if (b.port / 2 == p / 2) continue;            // same contig
                const MarkerEdge* e = model.edge(a.oriented, b.oriented, cls);
                if (!e) continue;
                // The panel measures marker-start to marker-start. Here that
                // distance is spent on the tail of this contig, the gap, and
                // the head of the next.
                const int32_t gap = e->medianDist - a.offset - b.offset;
                if (gap < -(k - 1) || gap > kMaxGap) continue;
                const uint32_t both = std::min(
                    bothA, model.markerGenomes(static_cast<uint32_t>(b.oriented >> 1), cls));
                if (both == 0) continue;
                if (static_cast<double>(e->support) < minFraction * static_cast<double>(both)) continue;
                votesByPort[b.port].push_back({gap, e->support});
            }
        }

        for (auto& kv : votesByPort) {
            std::vector<Vote>& v = kv.second;
            if (v.empty()) continue;
            ++st.candidates;
            std::sort(v.begin(), v.end(), [](const Vote& x, const Vote& y) { return x.gap < y.gap; });
            // Largest run of votes agreeing on a gap. Marker pairs that
            // disagree about the distance are disagreeing about the join.
            size_t bestBegin = 0, bestLen = 0;
            for (size_t i = 0, j = 0; i < v.size(); ++i) {
                if (j < i) j = i;
                while (j < v.size() && v[j].gap - v[i].gap <= kGapTolerance) ++j;
                if (j - i > bestLen) { bestLen = j - i; bestBegin = i; }
            }
            uint32_t strongest = 0;
            for (size_t i = bestBegin; i < bestBegin + bestLen; ++i)
                strongest = std::max(strongest, v[i].support);
            if (bestLen < kMinAgreeingPairs) { ++st.rejectedInconsistent; continue; }
            if (strongest < minPanel) { ++st.rejectedWeak; continue; }
            if (bestLen > best[p].pairs ||
                (bestLen == best[p].pairs && strongest > best[p].support)) {
                best[p].port = kv.first;
                best[p].pairs = bestLen;
                best[p].support = strongest;
                best[p].gap = v[bestBegin + bestLen / 2].gap;
            }
        }
    }

    // ---- accept mutual bests ---------------------------------------------
    // Exit port 2c+e and entry port 2c+e are the same physical end of the same
    // contig -- leaving contig c at its right end and arriving into it at that
    // end are one junction described from either side. So the partner's exit
    // port is its entry port, and a mutual join is simply best[best[p]] == p.
    std::vector<uint32_t> joinTo(nports, UINT32_MAX);
    std::vector<int32_t> joinGap(nports, 0);
    size_t made = 0;
    for (uint32_t p = 0; p < nports; ++p) {
        const Best& b = best[p];
        if (b.port == UINT32_MAX) continue;
        const uint32_t mirror = b.port;
        if (best[mirror].port == UINT32_MAX) { ++st.rejectedNotMutual; continue; }
        if (best[mirror].port != p) { ++st.rejectedNotMutual; continue; }
        if (joinTo[p] != UINT32_MAX || joinTo[mirror] != UINT32_MAX) continue;
        joinTo[p] = mirror;
        joinTo[mirror] = p;
        joinGap[p] = joinGap[mirror] = std::max(1, b.gap);
        ++made;
    }
    if (made == 0) return 0;

    // ---- emit --------------------------------------------------------------
    std::vector<char> done(n, 0);
    std::vector<std::string> out;
    std::vector<double> outCov;

    for (size_t c = 0; c < n; ++c) {
        if (done[c]) continue;
        // Start where a chain of joins can only run one way, so each path is
        // walked from an end rather than from its middle.
        if (joinTo[c * 2 + 0] != UINT32_MAX && joinTo[c * 2 + 1] != UINT32_MAX) continue;

        std::string seq;
        double covWeighted = 0;
        size_t covLen = 0;
        size_t cur = c;
        int inEnd = (joinTo[c * 2 + 0] == UINT32_MAX) ? 0 : 1;
        size_t steps = 0;
        while (true) {
            done[cur] = 1;
            const std::string piece = inEnd == 0 ? contigs[cur] : reverseComplement(contigs[cur]);
            seq += piece;
            covWeighted += covs[cur] * static_cast<double>(piece.size());
            covLen += piece.size();

            const uint32_t outPort = static_cast<uint32_t>(cur * 2 + (inEnd == 0 ? 1 : 0));
            const uint32_t partner = joinTo[outPort];
            if (partner == UINT32_MAX || ++steps > n) break;
            const size_t nextC = partner / 2;
            if (done[nextC]) break;
            seq.append(static_cast<size_t>(joinGap[outPort]), 'N');
            st.gapBases += static_cast<size_t>(joinGap[outPort]);
            ++st.joins;
            cur = nextC;
            // Arriving through that port means entering at the opposite end.
            inEnd = (partner & 1u) ? 1 : 0;
        }
        out.push_back(std::move(seq));
        outCov.push_back(covLen ? covWeighted / static_cast<double>(covLen) : 0.0);
    }
    // Anything caught in a cycle of joins is emitted untouched.
    for (size_t c = 0; c < n; ++c) {
        if (done[c]) continue;
        out.push_back(contigs[c]);
        outCov.push_back(covs[c]);
    }

    contigs.swap(out);
    covs.swap(outCov);
    return made;
}

}  // namespace

OrganismJoinStats joinByModel(const OrganismModel& model, std::vector<std::string>& contigs,
                              std::vector<double>& covs, int k, bool verbose) {
    OrganismJoinStats st;
    st.contigsIn = contigs.size();
    st.contigsOut = contigs.size();
    if (!model.loaded() || contigs.size() < 2) return st;

    // The chromosome first: its evidence is the stronger of the two, and
    // resolving it removes chromosomal ends from contention before the weaker
    // plasmid table is consulted at all.
    st.chromosomeJoins = joinPass(model, Replicon::Chromosome, contigs, covs, k, st);
    st.plasmidJoins = joinPass(model, Replicon::Plasmid, contigs, covs, k, st);
    st.contigsOut = contigs.size();

    if (verbose) {
        std::fprintf(stderr,
                     "      %zu markers placed, %zu joins (%zu chromosomal, %zu plasmid) "
                     "spanning %zu gap bases\n",
                     st.markersFound, st.joins, st.chromosomeJoins, st.plasmidJoins, st.gapBases);
        if (std::getenv("TESSERA_DEBUG_ORGANISM")) {
            std::fprintf(stderr,
                         "      [debug] organism join: candidates=%zu weak=%zu "
                         "inconsistent=%zu not-mutual=%zu\n",
                         st.candidates, st.rejectedWeak, st.rejectedInconsistent,
                         st.rejectedNotMutual);
        }
    }
    return st;
}

}  // namespace ts
