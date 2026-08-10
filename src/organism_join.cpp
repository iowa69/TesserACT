#include "organism_join.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
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

// How far two contig ends may OVERLAP and still be joined.
//
// The graph breaks at a repeat, but it does not always break cleanly on either side of
// it: often both contigs run the whole way through the shared element, so the two ends
// carry the same sequence and the distance the panel implies between them is negative
// by exactly the element's length. Measured on this cohort, such junctions overlap by
// a median of 638 bp and up to about 3 kb, and the recurring lengths -- 1202, 1445,
// 1914 -- are the IS families the panel already knows.
//
// Admitting them is OFF by default, and the history is worth keeping because the obvious
// reading of it is wrong.
//
// Raising the bound to 5000 doubled the share of the chromosome in the largest contig,
// 27% to 57%, which looked decisive. Requiring the overlapping sequence to match exactly
// looked like a strong safeguard. It is not: two different *copies* of the same repeat
// match exactly along their whole length, so the test confirms a correct merge and a
// collapsed repeat equally, and no minimum length separates them -- measured, 50 bp and
// 300 bp give identical results.
//
// Against the panel layout the whole change is dominated. Laying out the pre-overlap
// contigs gives the same closure (95.6% of isolates at >=90% against 95.2%), FEWER
// misassemblies (median 4 against 6; better on 72 isolates of 106, worse on 16) and a
// longer NGA50 (3.03 Mb against 2.07 Mb). The extra aligned bases the merge contributes
// arrive with misassemblies, and NGA50 breaks the scaffold at every one of them.
//
// So the bound returns to -(k-1), the largest overlap the de Bruijn graph itself can
// produce. TESSERA_MODEL_MAX_OVERLAP still raises it, because the measurement that
// retired it should stay reproducible.
constexpr int32_t kMaxOverlap = 0;   // 0 = fall back to -(k-1)

// Smallest exact suffix/prefix match accepted as an overlap. Shorter than this and a
// coincidental match becomes plausible in a 5 Mb genome.
constexpr size_t kMinOverlapMerge = 50;

// The panel's median distance is only accurate to about kGapTolerance, so the exact
// overlap is looked for in a window around the predicted length rather than at it.
constexpr int32_t kOverlapSearchSlack = 600;

// Largest L within the search window such that `left` ends with the same L bases that
// `right` begins with. Returns 0 when no overlap of at least kMinOverlapMerge is exact.
//
// Exactness is the whole safeguard. The model proposes the join; the sequence has to
// confirm it. A predicted overlap the sequence does not corroborate is not trimmed,
// because trimming it would delete real bases.
size_t mergeOnOverlap(const std::string& left, const std::string& right, int32_t predicted,
                      size_t minMerge) {
    if (predicted <= 0) return 0;
    const size_t hi = std::min<size_t>(
        {static_cast<size_t>(predicted) + static_cast<size_t>(kOverlapSearchSlack),
         left.size(), right.size()});
    if (hi < minMerge) return 0;
    const size_t lo = std::max<size_t>(
        minMerge,
        predicted > kOverlapSearchSlack ? static_cast<size_t>(predicted - kOverlapSearchSlack)
                                        : minMerge);

    // Seed and verify, rather than testing every candidate length.
    //
    // An overlap of length L means `left` ends with right[0..L). So the first
    // kMinOverlapMerge bases of `right` must appear in `left` exactly L from its end.
    // Finding the few places that seed occurs at all costs one pass over the window;
    // testing every L instead costs a full comparison per candidate length, which on a
    // 5 kb window with no overlap present is about 31 million character comparisons for
    // a single junction -- and the no-overlap case is precisely the one that has to be
    // ruled out before the join can be butted with an N.
    const std::string seed = right.substr(0, minMerge);
    const size_t windowStart = left.size() - hi;
    size_t best = 0;
    size_t from = windowStart;
    while (true) {
        const size_t at = left.find(seed, from);
        if (at == std::string::npos || at > left.size() - minMerge) break;
        const size_t L = left.size() - at;
        if (L >= lo && L <= hi &&
            std::equal(right.begin(), right.begin() + static_cast<long>(L),
                       left.begin() + static_cast<long>(at))) {
            // Longest wins, and positions are scanned left to right, so the first
            // qualifying hit is already the longest.
            best = L;
            break;
        }
        from = at + 1;
    }
    return best;
}

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
    // 0 means "the graph's own bound", which depends on k and so cannot be a constant.
    const int32_t maxOverlapCfg =
        static_cast<int32_t>(envNum("TESSERA_MODEL_MAX_OVERLAP", kMaxOverlap));
    const int32_t maxOverlap = maxOverlapCfg > 0 ? maxOverlapCfg : (k - 1);
    const size_t minOverlapMerge =
        static_cast<size_t>(envNum("TESSERA_MODEL_MIN_OVERLAP", kMinOverlapMerge));

    // Optional per-port diagnostic. The aggregate counters say how many candidates
    // were rejected but not whether the *correct* partner was among them, which is
    // the only question that tells you which threshold is worth moving. Ports are
    // identified by the terminal 31-mer of the end they leave from, so the caller
    // can locate them on a reference and line them up against the true junctions.
    std::FILE* dump = nullptr;
    if (const char* dp = std::getenv("TESSERA_JOIN_DUMP")) {
        dump = std::fopen(dp, "a");
        // The header marks every pass, not just the first. Both passes append to one
        // file and the plasmid pass runs on contigs the chromosome pass already
        // rewrote, so without a separator its port signatures would be read as if
        // they described the same contigs.
        if (dump) {
            std::fprintf(dump, "# pass=%s ncontigs=%zu\n",
                         cls == Replicon::Chromosome ? "chr" : "pls", n);
        }
    }
    // Terminal 31-mer of the sequence a port leaves from, oriented so the port's
    // own end comes last -- the same convention the site table uses.
    auto portSignature = [&](uint32_t port) -> std::string {
        const std::string& s = contigs[port / 2];
        if (s.size() < static_cast<size_t>(kMarkerK)) return std::string("NA");
        return (port & 1u) ? s.substr(s.size() - kMarkerK)
                           : reverseComplement(s.substr(0, kMarkerK));
    };

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
        }, model.markerDenom());
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
    if (unique.empty()) {
        if (dump) std::fclose(dump);
        identity();
        return 0;
    }

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
    // The veto now applies to the CHROMOSOME as well, which is a change of position that a
    // census of 57,119 IS copies across 2,910 closed chromosomes forced.
    //
    // The old reasoning was that chromosomal insertion sites recur across the panel, so an
    // end lying on one can be placed rather than refused. Measured, they recur as *sites*
    // and not as *insertions*: the median locus is occupied in 0.008 of the panel genomes
    // that carry both its flanks, and 88-99.9% of chromosomal copies sit where the panel
    // lacks the element entirely. So the panel's adjacency across such a point describes a
    // genome without the insertion, and acting on it deletes this isolate's element -- the
    // exact failure the veto exists to prevent, on the replicon where it was switched off.
    //
    // Copy-number stability does not rescue any family. IS1182 has a copy-number CV among
    // carriers of 0.05, the most stable in the census, and the panel still lacks it at its
    // locus 85% of the time, because it is a lineage marker rather than species furniture.
    //
    // Ports the site table can place are exempted below, which is what the old
    // chromosome-only site machinery was for: a locus the panel has actually seen intact
    // can be sized, and refusing it would discard the better evidence.
    const bool vetoChromosome = [] {
        const char* e = std::getenv("TESSERA_IS_VETO_CHR");
        return !(e && *e == '0');
    }();
    std::vector<char> isFlanked(nports, 0);
    if (isPanel && isPanel->loaded() &&
        (cls == Replicon::Plasmid || vetoChromosome)) {
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

    // Every candidate that passed both thresholds, kept so acceptance can consider
    // the whole field at once rather than each port's own favourite in isolation.
    struct Cand2 {
        uint32_t exitPort, entryPort;
        size_t pairs;
        uint32_t support;
        int32_t gap;
    };
    std::vector<Cand2> eligible;

    std::unordered_map<uint32_t, std::vector<Vote>> votesByPort;
    for (uint32_t p = 0; p < nports; ++p) {
        // Logged before the empty check: a port with no markers in its window is a
        // different failure from a port whose candidates were all rejected, and
        // conflating the two would point at the thresholds when the fix is coverage.
        if (dump) {
            std::fprintf(dump, "port\t%s\t%zu\t%zu\t%zu\n", portSignature(p).c_str(),
                         exitM[p].size(), entryM[p].size(), contigs[p / 2].size());
        }
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
                if (gap < -maxOverlap || gap > kMaxGap) continue;
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
            if (dump) {
                std::fprintf(dump, "cand\t%s\t%s\t%zu\t%u\t%d\t%zu\t%s\n",
                             portSignature(p).c_str(), portSignature(kv.first).c_str(),
                             bestLen, strongest, v[bestBegin + bestLen / 2].gap, v.size(),
                             bestLen < minPairs ? "inconsistent"
                                 : (strongest < minPanel ? "weak" : "eligible"));
            }
            if (bestLen < minPairs) { ++st.rejectedInconsistent; continue; }
            if (strongest < minPanel) { ++st.rejectedWeak; continue; }
            const int32_t gapv = v[bestBegin + bestLen / 2].gap;
            eligible.push_back({p, kv.first, bestLen, strongest, gapv});
            if (bestLen > best[p].pairs ||
                (bestLen == best[p].pairs && strongest > best[p].support)) {
                best[p].port = kv.first;
                best[p].pairs = bestLen;
                best[p].support = strongest;
                best[p].gap = gapv;
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
    // Mutual best has a specific, visible failure: port p's true partner q is discarded
    // whenever q happens to prefer some third contig r, even when p-q rests on far more
    // evidence than q-r. On this cohort the junctions lost that way carried a median of
    // 19 agreeing marker pairs and 239 genomes of panel support -- not marginal calls.
    //
    // Greedy matching was the obvious fix: take the whole field in strength order and
    // accept a pair whenever both ends are free, so a correct pair can only be beaten by
    // a stronger claim on one of its own ends. Measured over all 666 isolates it is
    // WORSE -- median chromosome-in-one-contig 55.0% against 57.3%, and 3 wins to 21
    // losses. Requiring both ends to prefer each other is not merely a tie-break; it is
    // a consistency check, and strength alone lets a well-supported but wrong pair
    // consume an end that a better join needed. So mutual best stays the default and
    // greedy is kept, off, because the negative result is worth being able to reproduce.
    const bool greedyMatch = [] {
        const char* e = std::getenv("TESSERA_MODEL_MATCH");
        return e && std::strcmp(e, "greedy") == 0;
    }();

    if (greedyMatch) {
        std::sort(eligible.begin(), eligible.end(), [](const Cand2& x, const Cand2& y) {
            if (x.pairs != y.pairs) return x.pairs > y.pairs;
            if (x.support != y.support) return x.support > y.support;
            // Deterministic tie-break: without it the order depends on the hash map's
            // iteration order and the same input could join differently between runs.
            if (x.exitPort != y.exitPort) return x.exitPort < y.exitPort;
            return x.entryPort < y.entryPort;
        });
        for (const Cand2& c : eligible) {
            if (c.exitPort / 2 == c.entryPort / 2) continue;
            if (joinTo[c.exitPort] != UINT32_MAX || joinTo[c.entryPort] != UINT32_MAX) {
                ++st.rejectedNotMutual;   // outcompeted for one of its ends
                continue;
            }
            joinTo[c.exitPort] = c.entryPort;
            joinTo[c.entryPort] = c.exitPort;
            // The raw gap, sign included. A negative one is an overlap to be trimmed at
            // emission, not a gap to be filled; clamping it to 1 here is what duplicated
            // the overlapping bases into the output.
            joinGap[c.exitPort] = joinGap[c.entryPort] = c.gap;
            ++made;
        }
    } else {
        for (uint32_t p = 0; p < nports; ++p) {
            const Best& b = best[p];
            if (b.port == UINT32_MAX) continue;
            const uint32_t mirror = b.port;
            if (best[mirror].port == UINT32_MAX) { ++st.rejectedNotMutual; continue; }
            if (best[mirror].port != p) { ++st.rejectedNotMutual; continue; }
            if (joinTo[p] != UINT32_MAX || joinTo[mirror] != UINT32_MAX) continue;
            joinTo[p] = mirror;
            joinTo[mirror] = p;
            joinGap[p] = joinGap[mirror] = b.gap;
            ++made;
        }
    }
    if (dump) {
        for (uint32_t p = 0; p < nports; ++p) {
            if (joinTo[p] == UINT32_MAX) continue;
            std::fprintf(dump, "join\t%s\t%s\t%d\n", portSignature(p).c_str(),
                         portSignature(joinTo[p]).c_str(), joinGap[p]);
        }
        std::fclose(dump);
        dump = nullptr;
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
        // The gap belonging to the join that brought us to `cur`, carried forward
        // because a negative one can only be resolved against the incoming sequence.
        int32_t pendingGap = 0;
        bool havePending = false;
        while (true) {
            done[cur] = 1;
            ++members;
            std::string piece = inEnd == 0 ? contigs[cur] : reverseComplement(contigs[cur]);
            if (havePending) {
                if (pendingGap < 0) {
                    const size_t trim = mergeOnOverlap(seq, piece, -pendingGap, minOverlapMerge);
                    if (trim) {
                        piece.erase(0, trim);
                        ++st.overlapMerges;
                    } else {
                        // The panel says these ends overlap and the sequence does not
                        // agree. Rather than delete bases on the panel's word alone,
                        // the join is still made but butted with a single N, which
                        // marks the seam without inventing or destroying anything.
                        ++st.overlapUnconfirmed;
                        seq.append(1, 'N');
                        st.gapBases += 1;
                    }
                } else {
                    const size_t g = static_cast<size_t>(std::max(1, pendingGap));
                    seq.append(g, 'N');
                    st.gapBases += g;
                }
            }
            seq += piece;
            covWeighted += covs[cur] * static_cast<double>(piece.size());
            covLen += piece.size();

            const uint32_t outPort = static_cast<uint32_t>(cur * 2 + (inEnd == 0 ? 1 : 0));
            const uint32_t partner = joinTo[outPort];
            if (partner == UINT32_MAX || ++steps > n) break;
            const size_t nextC = partner / 2;
            if (done[nextC]) break;
            pendingGap = joinGap[outPort];
            havePending = true;
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

        // Read signed, then checked: a negative count cast straight to uint32_t
        // became about four billion and sailed past the "< 2" filter below.
        const long genomesSigned = std::atol(col[1].c_str());
        if (genomesSigned < 2 || genomesSigned > 100000000L) continue;
        const uint32_t genomes = static_cast<uint32_t>(genomesSigned);
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

    // A pass accepts only mutual bests, so a port whose true partner happens to
    // prefer a third contig gets nothing -- and once that third contig is consumed,
    // the preference that blocked the join no longer exists. Joining also lengthens
    // contigs, which brings different markers inside the window at the ends that
    // remain. Both effects mean a second pass has strictly more to work with than
    // the first, so passes repeat until one makes no progress.
    static const int kMaxRounds = [] {
        const char* e = std::getenv("TESSERA_MODEL_ROUNDS");
        const int v = e ? std::atoi(e) : 4;
        return v < 1 ? 1 : (v > 32 ? 32 : v);
    }();

    // Composed input->output mapping, rebuilt after every pass. `acc[i]` is the
    // index in the ORIGINAL contig vector that output contig i came from, or
    // UINT32_MAX once it is a join of several.
    std::vector<uint32_t> acc;
    if (source) {
        acc.resize(contigs.size());
        for (size_t i = 0; i < acc.size(); ++i) acc[i] = static_cast<uint32_t>(i);
    }
    auto compose = [&](const std::vector<uint32_t>& step) {
        if (!source || step.empty()) return;
        std::vector<uint32_t> next(step.size(), UINT32_MAX);
        for (size_t i = 0; i < step.size(); ++i) {
            if (step[i] != UINT32_MAX && step[i] < acc.size()) next[i] = acc[step[i]];
        }
        acc.swap(next);
    };

    // The chromosome first in every round: its evidence is the stronger of the two,
    // and resolving it removes chromosomal ends from contention before the weaker
    // plasmid table is consulted at all.
    for (int round = 0; round < kMaxRounds; ++round) {
        std::vector<uint32_t> m1, m2;
        const size_t chr = joinPass(model, Replicon::Chromosome, contigs, covs, k, st,
                                    isPanel, source ? &m1 : nullptr);
        compose(m1);
        const size_t pls = joinPass(model, Replicon::Plasmid, contigs, covs, k, st,
                                    isPanel, source ? &m2 : nullptr);
        compose(m2);
        st.chromosomeJoins += chr;
        st.plasmidJoins += pls;
        st.rounds = static_cast<size_t>(round) + 1;
        if (chr + pls == 0) break;
    }
    st.contigsOut = contigs.size();
    if (source) source->swap(acc);

    if (verbose) {
        std::fprintf(stderr,
                     "      %zu markers placed, %zu joins (%zu chromosomal, %zu plasmid, "
                     "%zu placed at known insertion sites) spanning %zu gap bases "
                     "in %zu rounds\n",
                     st.markersFound, st.joins, st.chromosomeJoins, st.plasmidJoins,
                     st.insertionSiteJoins, st.gapBases, st.rounds);
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
