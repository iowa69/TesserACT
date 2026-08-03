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

// How much of a contig end is inspected for insertion-sequence content, and
// how much of that window must match before the end is considered to be inside
// a mobile element rather than merely near one.
constexpr size_t kIsWindow = 600;
constexpr double kIsDensity = 0.30;

constexpr int32_t kGapTolerance = 500;    // gaps within this are the same cluster
constexpr int32_t kMaxGap = 20000;

// Canonical 31-mer of a sequence's first k bases, or UINT64_MAX if it contains
// anything but ACGT.
uint64_t canonicalHead(const std::string& seq) {
    if (seq.size() < static_cast<size_t>(kMarkerK)) return UINT64_MAX;
    const uint64_t mask = (1ULL << (2 * kMarkerK)) - 1;
    uint64_t fwd = 0, rev = 0;
    for (int i = 0; i < kMarkerK; ++i) {
        const int b = markerBaseCode(seq[static_cast<size_t>(i)]);
        if (b < 0) return UINT64_MAX;
        fwd = ((fwd << 2) | static_cast<uint64_t>(b)) & mask;
        rev = (rev >> 2) | (static_cast<uint64_t>(3 - b) << (2 * (kMarkerK - 1)));
    }
    return fwd <= rev ? fwd : rev;
}

uint64_t canonicalTail(const std::string& seq) {
    if (seq.size() < static_cast<size_t>(kMarkerK)) return UINT64_MAX;
    return canonicalHead(seq.substr(seq.size() - static_cast<size_t>(kMarkerK)));
}



// Every canonical 31-mer of `seq`, with its start position.
template <typename F>
void forEachKmer31(const std::string& seq, F&& fn) {
    if (seq.size() < static_cast<size_t>(kMarkerK)) return;
    const uint64_t mask = (1ULL << (2 * kMarkerK)) - 1;
    uint64_t fwd = 0, rev = 0;
    int valid = 0;
    for (size_t i = 0; i < seq.size(); ++i) {
        const int b = markerBaseCode(seq[i]);
        if (b < 0) { valid = 0; fwd = rev = 0; continue; }
        fwd = ((fwd << 2) | static_cast<uint64_t>(b)) & mask;
        rev = (rev >> 2) | (static_cast<uint64_t>(3 - b) << (2 * (kMarkerK - 1)));
        if (++valid < kMarkerK) continue;
        fn(fwd <= rev ? fwd : rev, static_cast<int32_t>(i + 1 - static_cast<size_t>(kMarkerK)));
    }
}

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
                std::vector<double>& covs, int k, OrganismJoinStats& st,
                const IsPanel* isPanel, std::vector<uint32_t>* source) {
    const size_t n = contigs.size();
    // `source` maps each emitted contig back to the single input it came from,
    // or UINT32_MAX when it was built by joining several. Callers holding data
    // parallel to `contigs` need it to stay in step. Filled on every return
    // path, identity included, so the caller never has to guess.
    auto identity = [&]() {
        if (!source) return;
        source->resize(n);
        for (size_t i = 0; i < n; ++i) (*source)[i] = static_cast<uint32_t>(i);
    };
    if (n < 2) { identity(); return 0; }

    // Tunable for the threshold sweep; the shipped values are the measured
    // ones, not the first ones that worked on a single isolate.
    auto envNum = [](const char* name, double dflt) {
        const char* e = std::getenv(name);
        return e ? std::atof(e) : dflt;
    };
    const uint32_t minPanel = static_cast<uint32_t>(envNum(
        "TESSERA_MODEL_MIN_PANEL", cls == Replicon::Chromosome ? kMinPanelChr : kMinPanelPls));
    const double minFraction = envNum(
        "TESSERA_MODEL_MIN_FRACTION", cls == Replicon::Chromosome ? kMinFractionChr : kMinFractionPls);
    const size_t minPairs = static_cast<size_t>(envNum("TESSERA_MODEL_MIN_PAIRS", kMinAgreeingPairs));

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

    // An end lying inside a known insertion sequence is where it is *because*
    // this isolate carries an element the panel need not have. The panel's
    // adjacency across that point describes a genome without the insertion, so
    // acting on it deletes the element. Silence those ports entirely.
    // The two replicons are not the same problem. Chromosomal insertion sites
    // recur across the panel, so an end lying on one can be *placed*: the
    // element's length is known and the join can be sized. Plasmids are mosaic,
    // their insertion sites are not stable, and there the only safe move is to
    // leave the end alone. So placement runs on the chromosome and the veto on
    // the plasmid.
    std::vector<char> isFlanked(nports, 0);
    if (isPanel && isPanel->loaded() && cls == Replicon::Plasmid) {
        for (size_t c = 0; c < n; ++c) {
            const std::string& seq = contigs[c];
            const size_t w = std::min<size_t>(kIsWindow, seq.size());
            if (isPanel->density(seq, 0, w) >= kIsDensity) isFlanked[c * 2 + 0] = 1;
            if (isPanel->density(seq, seq.size() - w, w) >= kIsDensity) isFlanked[c * 2 + 1] = 1;
        }
        for (size_t p = 0; p < nports; ++p) {
            if (isFlanked[p]) { exitM[p].clear(); entryM[p].clear(); ++st.rejectedInsertionSeq; }
        }
    }

    // ---- place elements at known insertion sites (chromosome only) --------
    // A contig ending on the left flank of a recurrent site, and another
    // starting on its right flank, are the two halves of a locus the panel has
    // seen intact. The element between them has a known length, so the join can
    // be made and sized rather than refused.
    std::vector<uint32_t> siteJoin(nports, UINT32_MAX);
    std::vector<int32_t> siteGap(nports, 0);
    if (isPanel && isPanel->siteCount() && cls == Replicon::Chromosome) {
        // A contig stops wherever the graph broke, which is rarely the flank
        // boundary itself, so the signature is looked for across a window and
        // its offset from the end is carried into the gap arithmetic.
        constexpr size_t kSiteWindow = 3000;
        struct Occ { uint64_t kmer; int32_t offset; uint32_t port; };
        std::vector<Occ> tails, heads;
        for (size_t c = 0; c < n; ++c) {
            const std::string& seq = contigs[c];
            if (seq.size() < static_cast<size_t>(kMarkerK)) continue;
            const std::string rc = reverseComplement(seq);
            const int32_t L = static_cast<int32_t>(seq.size());
            forEachKmer31(seq, [&](uint64_t km, int32_t pos) {
                const int32_t fromRight = L - pos - kMarkerK;
                if (fromRight >= 0 && fromRight <= static_cast<int32_t>(kSiteWindow))
                    tails.push_back({km, fromRight, static_cast<uint32_t>(c * 2 + 1)});
                if (pos <= static_cast<int32_t>(kSiteWindow))
                    heads.push_back({km, pos, static_cast<uint32_t>(c * 2 + 0)});
            });
            forEachKmer31(rc, [&](uint64_t km, int32_t pos) {
                const int32_t fromRight = L - pos - kMarkerK;
                if (fromRight >= 0 && fromRight <= static_cast<int32_t>(kSiteWindow))
                    tails.push_back({km, fromRight, static_cast<uint32_t>(c * 2 + 0)});
                if (pos <= static_cast<int32_t>(kSiteWindow))
                    heads.push_back({km, pos, static_cast<uint32_t>(c * 2 + 1)});
            });
        }
        std::unordered_map<uint64_t, std::vector<const Occ*>> headIdx;
        for (const Occ& h : heads) headIdx[h.kmer].push_back(&h);

        struct Cand { uint32_t port = UINT32_MAX; int32_t gap = 0; uint32_t support = 0; uint32_t hits = 0; };
        std::vector<Cand> cand(nports);
        for (const Occ& t : tails) {
            const std::vector<IsSite>* sites = isPanel->sitesFor(t.kmer);
            if (!sites) continue;
            for (const IsSite& site : *sites) {
                auto it = headIdx.find(site.rightKmer);
                if (it == headIdx.end()) continue;
                for (const Occ* h : it->second) {
                    if (h->port / 2 == t.port / 2) continue;
                    // Panel distance runs flank-signature to flank-signature
                    // through the element; here it is spent on this contig's
                    // tail, the gap, and the next contig's head.
                    const int32_t gap = site.elementLen - t.offset - h->offset;
                    if (gap < 1 || gap > 20000) continue;
                    Cand& c = cand[t.port];
                    ++c.hits;
                    if (site.genomes > c.support) {
                        c.support = site.genomes;
                        c.port = h->port;
                        c.gap = gap;
                    }
                }
            }
        }
        for (uint32_t p = 0; p < nports; ++p) {
            // One partner only: several means the site does not name a locus here.
            if (cand[p].port != UINT32_MAX && cand[p].hits >= 1) {
                bool unique = true;
                for (uint32_t q = 0; q < nports && unique; ++q) {
                    if (q != p && cand[q].port == cand[p].port && cand[q].support > cand[p].support)
                        unique = false;
                }
                if (unique) { siteJoin[p] = cand[p].port; siteGap[p] = cand[p].gap; }
            }
        }
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
            if (bestLen < minPairs) { ++st.rejectedInconsistent; continue; }
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
    // Insertion-site placements first: they are the more specific evidence,
    // naming a locus rather than an adjacency, and both ends must agree.
    for (uint32_t p = 0; p < nports; ++p) {
        const uint32_t q = siteJoin[p];
        if (q == UINT32_MAX || siteJoin[q] != p) continue;
        if (joinTo[p] != UINT32_MAX || joinTo[q] != UINT32_MAX) continue;
        joinTo[p] = q;
        joinTo[q] = p;
        joinGap[p] = joinGap[q] = std::max(1, siteGap[p]);
        ++made;
        ++st.insertionSiteJoins;
    }
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
    if (made == 0) { identity(); return 0; }

    // ---- emit --------------------------------------------------------------
    std::vector<char> done(n, 0);
    std::vector<std::string> out;
    std::vector<double> outCov;
    std::vector<uint32_t> outSrc;

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
        size_t members = 0;
        while (true) {
            done[cur] = 1;
            ++members;
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
        outSrc.push_back(members == 1 ? static_cast<uint32_t>(c) : UINT32_MAX);
    }
    // Anything caught in a cycle of joins is emitted untouched.
    for (size_t c = 0; c < n; ++c) {
        if (done[c]) continue;
        out.push_back(contigs[c]);
        outCov.push_back(covs[c]);
        outSrc.push_back(static_cast<uint32_t>(c));
    }

    contigs.swap(out);
    covs.swap(outCov);
    if (source) source->swap(outSrc);
    return made;
}


}  // namespace

bool IsPanel::loadSites(const std::string& tsvPath, std::string& error) {
    std::FILE* f = std::fopen(tsvPath.c_str(), "r");
    if (!f) { error = "cannot open insertion-site table: " + tsvPath; return false; }
    char buf[1 << 16];
    bool header = true;
    size_t kept = 0;
    while (std::fgets(buf, sizeof(buf), f)) {
        if (header) { header = false; continue; }
        // family, genomes_sharing_site, median_element_len, left_sig, right_sig
        std::string line(buf);
        std::vector<std::string> col;
        size_t pos = 0;
        while (col.size() < 5) {
            const size_t tab = line.find('\t', pos);
            if (tab == std::string::npos) { col.push_back(line.substr(pos)); break; }
            col.push_back(line.substr(pos, tab - pos));
            pos = tab + 1;
        }
        if (col.size() < 5) continue;
        while (!col[4].empty() && (col[4].back() == '\n' || col[4].back() == '\r')) col[4].pop_back();

        const uint32_t genomes = static_cast<uint32_t>(std::atoi(col[1].c_str()));
        // A site seen in one genome is that genome's accident, not a site.
        if (genomes < 2) continue;
        const int32_t elen = static_cast<int32_t>(std::atoi(col[2].c_str()));
        if (elen <= 0 || elen > 20000) continue;

        const uint64_t left = canonicalTail(col[3]);
        const uint64_t right = canonicalHead(col[4]);
        if (left == UINT64_MAX || right == UINT64_MAX) continue;
        sites_[left].push_back(IsSite{right, elen, genomes});
        // The mirrored walk: arriving from the other side, the roles swap.
        sites_[right].push_back(IsSite{left, elen, genomes});
        ++kept;
    }
    std::fclose(f);
    if (kept == 0) { error = "no recurrent insertion sites in " + tsvPath; return false; }
    return true;
}

bool IsPanel::load(const std::string& fastaPath, std::string& error) {
    std::FILE* f = std::fopen(fastaPath.c_str(), "r");
    if (!f) { error = "cannot open IS panel: " + fastaPath; return false; }
    kmers_.clear();
    copies_ = 0;
    std::string cur;
    char buf[1 << 16];
    auto absorb = [&]() {
        if (cur.empty()) return;
        ++copies_;
        // Every 31-mer of the element, not a sample: the query window is short
        // and a missed k-mer is a missed veto.
        const uint64_t mask = (1ULL << (2 * kMarkerK)) - 1;
        uint64_t fwd = 0, rev = 0;
        int valid = 0;
        for (char c : cur) {
            const int b = markerBaseCode(c);
            if (b < 0) { valid = 0; fwd = rev = 0; continue; }
            fwd = ((fwd << 2) | static_cast<uint64_t>(b)) & mask;
            rev = (rev >> 2) | (static_cast<uint64_t>(3 - b) << (2 * (kMarkerK - 1)));
            if (++valid < kMarkerK) continue;
            kmers_.insert(fwd <= rev ? fwd : rev);
        }
        cur.clear();
    };
    while (std::fgets(buf, sizeof(buf), f)) {
        if (buf[0] == '>') { absorb(); continue; }
        for (char* p = buf; *p && *p != '\n' && *p != '\r'; ++p) cur.push_back(*p);
    }
    absorb();
    std::fclose(f);
    if (kmers_.empty()) { error = "IS panel contained no usable sequence: " + fastaPath; return false; }
    return true;
}

double IsPanel::density(const std::string& seq, size_t from, size_t len) const {
    if (kmers_.empty() || seq.size() < static_cast<size_t>(kMarkerK)) return 0.0;
    const size_t end = std::min(seq.size(), from + len);
    if (end < from + static_cast<size_t>(kMarkerK)) return 0.0;
    const uint64_t mask = (1ULL << (2 * kMarkerK)) - 1;
    uint64_t fwd = 0, rev = 0;
    int valid = 0;
    size_t total = 0, hit = 0;
    for (size_t i = from; i < end; ++i) {
        const int b = markerBaseCode(seq[i]);
        if (b < 0) { valid = 0; fwd = rev = 0; continue; }
        fwd = ((fwd << 2) | static_cast<uint64_t>(b)) & mask;
        rev = (rev >> 2) | (static_cast<uint64_t>(3 - b) << (2 * (kMarkerK - 1)));
        if (++valid < kMarkerK) continue;
        ++total;
        if (kmers_.count(fwd <= rev ? fwd : rev)) ++hit;
    }
    return total ? static_cast<double>(hit) / static_cast<double>(total) : 0.0;
}

OrganismJoinStats joinByModel(const OrganismModel& model, std::vector<std::string>& contigs,
                              std::vector<double>& covs, int k, bool verbose,
                              const IsPanel* isPanel, std::vector<uint32_t>* source) {
    OrganismJoinStats st;
    st.contigsIn = contigs.size();
    st.contigsOut = contigs.size();
    if (source) {
        source->resize(contigs.size());
        for (size_t i = 0; i < contigs.size(); ++i) (*source)[i] = static_cast<uint32_t>(i);
    }
    if (!model.loaded() || contigs.size() < 2) return st;

    // The chromosome first: its evidence is the stronger of the two, and
    // resolving it removes chromosomal ends from contention before the weaker
    // plasmid table is consulted at all.
    std::vector<uint32_t> m1, m2;
    st.chromosomeJoins =
        joinPass(model, Replicon::Chromosome, contigs, covs, k, st, isPanel, source ? &m1 : nullptr);
    st.plasmidJoins =
        joinPass(model, Replicon::Plasmid, contigs, covs, k, st, isPanel, source ? &m2 : nullptr);
    st.contigsOut = contigs.size();

    // Compose the two passes so the caller sees one input->output mapping.
    if (source) {
        source->assign(m2.size(), UINT32_MAX);
        for (size_t i = 0; i < m2.size(); ++i) {
            if (m2[i] != UINT32_MAX && m2[i] < m1.size()) (*source)[i] = m1[m2[i]];
        }
    }

    if (verbose) {
        std::fprintf(stderr,
                     "      %zu markers placed, %zu joins (%zu chromosomal, %zu plasmid, "
                     "%zu placed at known insertion sites) spanning %zu gap bases\n",
                     st.markersFound, st.joins, st.chromosomeJoins, st.plasmidJoins,
                     st.insertionSiteJoins, st.gapBases);
        if (std::getenv("TESSERA_DEBUG_ORGANISM")) {
            std::fprintf(stderr,
                         "      [debug] organism join: candidates=%zu weak=%zu "
                         "inconsistent=%zu not-mutual=%zu is-flanked=%zu\n",
                         st.candidates, st.rejectedWeak, st.rejectedInconsistent,
                         st.rejectedNotMutual, st.rejectedInsertionSeq);
        }
    }
    return st;
}

}  // namespace ts
