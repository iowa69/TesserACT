#include "resolve.h"
#include "emit_fixflags.h"
#include "envflags.h"
#include "graph_coverage.h"
#include "resolve_evidence.h"
#include "read_thread_evidence.h"

#include "organism.h"   // forEachMarkerKmer, for the plasmid vouch below

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <mutex>
#include <cstdlib>
#include <thread>

namespace ts {

namespace {

// ---- build_v3 G-resolve fixes ----------------------------------------------------------
// Every output-changing fix below is default OFF and is read per resolve() call, never
// cached in a function-local static, so a driver that runs several configurations in one
// process sees each one. Resolution: TESSERACT_FIX_<NAME>=1 on, =0 off; unset follows the
// umbrella TESSERACT_FIXES (=1 on, =0 off); both unset = off = release 1.3.0 behaviour.
// Any other value is an error (exit status 2), never a silent default: a flag that is set
// but mis-typed must not look like an arm that ran.
// build_v3: the values are read through the envflags table (T16), whose Binary kind is
// exactly this contract. The campaign package's TESSERACT_GAP_KEEP_FLANK (build_v2_merge
// PKG-GAP) fixed the same k-1 deletion as T03 with a render-whole design; build_v3 keeps ONE
// implementation, T03's record-and-restore, and the package flag is its alias. Precedence:
// the fix's own variable, then the package alias, then the umbrella.
bool fixEnabled(const char* name) {
    if (env::isSet(name)) return env::on(name, false);
    if (std::strcmp(name, "TESSERACT_FIX_GAP_FLANK") == 0 && env::isSet("TESSERACT_GAP_KEEP_FLANK"))
        return env::on("TESSERACT_GAP_KEEP_FLANK", false);
    return env::on("TESSERACT_FIXES", false);
}

// The resolver's repeat test: a unitig deeper than this multiple of the single-copy depth
// (theta) is treated as a collapsed repeat. One constant, so the value reported in the
// [resolver] line cannot drift from the value used.
constexpr double kRepeatMultiplier = 1.6;

// T03: the shortest exact suffix/prefix overlap between the two flanks of a scaffold gap
// that is taken as a real overlap and dropped once. A chance match of l bases has
// probability ~(3/4)4^-l; measured on 59 reference-placed K2 gaps (combo3/verify/T03,
// gapflank_rows.tsv), true overlaps had L = 1..108 (18/20 at L >= 4) and every chance
// match at D >= 0 had L <= 3 (one L=45 row is a misplaced join).
constexpr size_t kMinFlankOverlap = 4;
// T20: a join whose flanks the estimate says overlap, but which share no verified
// overlap, has no measurable length. It is written as this many Ns (AGP type U length).
constexpr int kUnknownGapN = 100;
// T03: bases of flank context on each side of an N-run that identify a gap record when
// restoreGapFlanks() runs after gap filling (both flanks are at least k+1 >= 32 bases).
constexpr size_t kGapFlankContext = 32;

// The always-printed resolver counter lines (OBJECTIVE amendment A2: every flag prints
// its counter on every run, zeros included). One struct and one printer, so a run in
// which the resolver never runs prints exactly the same lines with run=0. Decision
// counts from bestContinuation() are EVALUATIONS: every chain end is re-evaluated each
// round (up to 24), so a side refused in three rounds counts three times.
struct ResolverCounters {
    int run = 0;
    // [resolveflags] -- release flags in force (the values actually used) and outcomes (T14).
    int requireSupportSingle = 0, excludeSharedRepeat = 0, sharedAudit = 0, noUnspannedFallback = 0;
    size_t minFallbackDest = 0;
    double linkBar = 0, tieRatio = 0;
    long long evals = 0, ok = 0, noCand = 0, noPick = 0, tie = 0, midChain = 0, byCoverage = 0;
    long long matched = 0, shortDest = 0, unspanned = 0, loneRepeat = 0;
    long long sharedExcludedAll = 0, esrsChanged = 0;
    // [enumtrunc] (T21)
    int truncGuard = 0;
    long long trSides = 0, trBudget = 0, trMaxNodes = 0, trCap = 0, trLone = 0, trLoneBelowBar = 0, trRefused = 0;
    // [revisit] (T07)
    int revisitGuard = 0;
    long long rvPrunedSides = 0, rvPicked = 0, rvRefused = 0;
    // [scafcycle] (T13)
    int cycleBreak = 0;
    size_t scCycles = 0, scChains = 0, scJoins = 0, scDropped = 0, scBases = 0, scBroken = 0;
    // [gapest] (T20)
    int gapEstimate = 0;
    size_t geJoins = 0, geEstimated = 0, geChanged = 0, geFloor1 = 0, geOverlapEst = 0, geSkippedNoMean = 0;
    int geMedianShift = 0;
    // [gapflank] (T03)
    int gapFlank = 0;
    size_t gfGaps = 0, gfKept = 0, gfOverlapGaps = 0, gfOverlapBases = 0, gfOverlapRejected = 0;
    size_t gfFloor1 = 0, gfUnknown = 0;
    // [cov_contrib] (T24)
    int covContrib = 0;
    size_t cvPaths = 0, cvOff10 = 0;
    // [routeorder] (T31)
    int routeOrder = 0;
    long long roDecisions = 0, roSorted = 0;
    // [mirrorroute] (R6, housekeeping)
    int mirrorRoute = 0;
    long long mrDisagree = 0, mrChoseBack = 0;
};

void printResolverCounters(const ResolverCounters& c) {
    std::fprintf(stderr,
        "[resolveflags] run=%d requireSupportSingle=%d minFallbackDest=%zu excludeSharedRepeat=%d "
        "sharedAudit=%d noUnspannedFallback=%d linkBar=%.3f tieRatio=%.3f evals=%lld ok=%lld noCand=%lld "
        "noPick=%lld tie=%lld midChain=%lld byCoverage=%lld matched=%lld shortDest=%lld unspanned=%lld "
        "loneRepeat=%lld sharedExcludedAll=%lld esrsChanged=%lld\n",
        c.run, c.requireSupportSingle, c.minFallbackDest, c.excludeSharedRepeat, c.sharedAudit,
        c.noUnspannedFallback, c.linkBar, c.tieRatio, c.evals, c.ok, c.noCand, c.noPick, c.tie, c.midChain,
        c.byCoverage, c.matched, c.shortDest, c.unspanned, c.loneRepeat, c.sharedExcludedAll, c.esrsChanged);
    std::fprintf(stderr,
        "[enumtrunc] run=%d guard=%d sides=%lld budget=%lld maxNodes=%lld cap64=%lld lone=%lld "
        "loneBelowBar=%lld refused=%lld\n",
        c.run, c.truncGuard, c.trSides, c.trBudget, c.trMaxNodes, c.trCap, c.trLone, c.trLoneBelowBar,
        c.trRefused);
    std::fprintf(stderr, "[revisit] run=%d guard=%d prunedSides=%lld pickedSkipped=%lld refused=%lld\n",
        c.run, c.revisitGuard, c.rvPrunedSides, c.rvPicked, c.rvRefused);
    std::fprintf(stderr,
        "[scafcycle] run=%d enabled=%d cycles=%zu chains=%zu joinsInCycles=%zu joinsDropped=%zu "
        "cycleBases=%zu broken=%zu\n",
        c.run, c.cycleBreak, c.scCycles, c.scChains, c.scJoins, c.scDropped, c.scBases, c.scBroken);
    std::fprintf(stderr,
        "[gapest] run=%d enabled=%d joins=%zu estimated=%zu changed=%zu floor1=%zu overlapEst=%zu "
        "skippedNoMean=%zu medianShift=%d\n",
        c.run, c.gapEstimate, c.geJoins, c.geEstimated, c.geChanged, c.geFloor1, c.geOverlapEst,
        c.geSkippedNoMean, c.geMedianShift);
    std::fprintf(stderr,
        "[gapflank] run=%d enabled=%d gaps=%zu restorableBases=%zu overlapGaps=%zu overlapBases=%zu "
        "pLineOverspell=%zu overlapRejected=%zu floor1=%zu unknownLen=%zu\n",
        c.run, c.gapFlank, c.gfGaps, c.gfKept, c.gfOverlapGaps, c.gfOverlapBases, c.gfOverlapBases,
        c.gfOverlapRejected, c.gfFloor1, c.gfUnknown);
    std::fprintf(stderr, "[cov_contrib] run=%d enabled=%d paths=%zu full_vs_contributed_off10=%zu\n",
        c.run, c.covContrib, c.cvPaths, c.cvOff10);
    std::fprintf(stderr, "[routeorder] run=%d enabled=%d decisions=%lld sorted=%lld\n",
        c.run, c.routeOrder, c.roDecisions, c.roSorted);
    std::fprintf(stderr, "[mirrorroute] run=%d enabled=%d disagreements=%lld choseBack=%lld\n",
        c.run, c.mirrorRoute, c.mrDisagree, c.mrChoseBack);
}

// The distance enumeration may search when it is only establishing which
// continuations exist. Repeats in an enterobacterial genome run to a few
// kilobases (IS elements ~0.8-2.5 kb, rRNA operons ~5 kb); this covers them
// without letting the search wander across the chromosome.
constexpr double kTopologyReach = 6000.0;

// Coverage-continuity fallback, used only where the paired reads gave nothing.
// A candidate counts as on-depth within kDepthMatch of the chain's own depth,
// and the decision is only taken when every rival is at least kDepthSeparation
// away -- a 2x repeat sits at 1.0, so the bar clears real alternatives while
// refusing anything that merely looks different.
constexpr double kDepthMatch = 0.25;
constexpr double kDepthSeparation = 0.60;
constexpr double kDepthWindow = 20000.0;   // how much of the chain end to average

// Repeat resolution by paired matching: only for repeats short enough that no
// pair can span them (which is why they are unresolved), and the intended
// assignment must beat the crossed one by this factor.
// A vouched group must total at least this much sequence. At a sampling of one k-mer in
// 512 this is about three sampled k-mers, which is the least that can carry two marker
// hits; it is also the length below which replicon.cpp declines to classify a contig.
constexpr size_t kMinVouchLen = 1500;
// How many plasmid-exclusive markers a group must carry before its depth is excused.
constexpr uint32_t kMinPlasmidVouch = 2;

constexpr size_t kMaxMatchedRepeat = 3000;
constexpr double kMatchDominance = 3.0;

constexpr int kMaxProbes = 12;        // k-mer lookups per read before deciding
constexpr int kMinVotes = 2;
constexpr int kAnchorK = 31;
constexpr uint64_t kAmbiguous = UINT64_MAX;

inline uint64_t packIndex(uint32_t unitig, uint32_t pos, int strand) {
    return (static_cast<uint64_t>(unitig) << 33) | (static_cast<uint64_t>(pos) << 1) |
           static_cast<uint64_t>(strand & 1);
}
inline uint32_t idxUnitig(uint64_t v) { return static_cast<uint32_t>(v >> 33); }
inline uint32_t idxPos(uint64_t v) { return static_cast<uint32_t>((v >> 1) & 0xFFFFFFFFULL); }
inline int idxStrand(uint64_t v) { return static_cast<int>(v & 1); }

// combo2/ends (default OFF). TESSERACT_PAIR_ANCHORED_PREFIX=<m>|bar: a terminal walk is kept only as far
// as >= m read pairs vouch for it (-1 = the resolver's own linkBar). build_v3: read per call
// through the envflags table (an integer in [0, INT_MAX] or "bar"; anything else exits 2).
long pairAnchoredPrefixSetting() {
    const char* e = env::text("TESSERACT_PAIR_ANCHORED_PREFIX");
    if (!e) return 0L;
    if (std::string(e) == "bar") return -1L;
    return std::strtol(e, nullptr, 10);
}

}  // namespace

bool gapFlankFixEnabled() { return fixEnabled("TESSERACT_FIX_GAP_FLANK"); }

int checkFixCoupling(std::FILE* log) {
    if (gapFlankFixEnabled() && !emitfix::enabled(emitfix::kPolishSkipN)) {
        std::fprintf(log, "error: TESSERACT_FIX_GAP_FLANK needs TESSERACT_FIX_POLISH_SKIP_N "
                          "(see DEFECTS.md T02/T03)\n");
        return 2;
    }
    return 0;
}

PairedResolver::PairedResolver(const UnitigGraph& graph, const SequenceStore& reads, int threads,
                               int minLinkSupport, double tieRatio, double linkSupportPerX,
                               int minScaffoldSupport)
    : g_(graph), reads_(reads), threads_(threads > 0 ? threads : 1), k_(graph.k()),
      kMap_(std::min(kAnchorK, graph.k())),
      minLinkSupport_(minLinkSupport), tieRatio_(tieRatio),
      linkSupportPerX_(linkSupportPerX >= 0 ? linkSupportPerX : 0.10),
      minScaffoldSupport_(minScaffoldSupport) {
    medianCoverage_ = graph.medianCoverage();
    const double legacyMedian = medianCoverage_;
    const char* estimator = "unweighted";
    const bool useEligibleCoverage = env::on("TESSERACT_WEIGHTED_ELIGIBLE_COVERAGE", false);
    // Explicit precedence: eligible-only calibration wins if both experiments
    // are enabled. They are alternative estimators, never cumulative weights.
    // Length-weighted by default since 1.3.0. The unweighted node median lets a cloud of
    // short carry-floor fragments set the single-copy baseline: on senterica GCA006365335v1
    // 327/980 eligible unitigs sit at ~4x and carry 2% of the length, the median reads 5.73x
    // for a genome sequenced at 35x, every unitig above 9.2x becomes a "repeat" and the
    // resolver joins almost nothing. 146 isolates: NGA50 tie -> WIN, genome-fraction losses
    // 21 -> 15. TESSERACT_WEIGHTED_RESOLVER_COVERAGE=0 restores the unweighted median.
    const bool useWeightedCoverage = env::on("TESSERACT_WEIGHTED_RESOLVER_COVERAGE", true);
    if (useEligibleCoverage || useWeightedCoverage) {
        const double calibrated = useEligibleCoverage
            ? graphEligibleLengthWeightedMedianCoverage(graph)
            : graphLengthWeightedMedianCoverage(graph);
        if (calibrated > 0.0) medianCoverage_ = calibrated;
        estimator = calibrated > 0.0 ? (useEligibleCoverage ? "eligible_weighted" : "weighted")
                                     : "unweighted_fallback";
    }

    // T35: theta (the single-copy depth every repeat decision is taken against) and the
    // repeat threshold were never reported; report.json's median_coverage is the legacy
    // unweighted node median, which on 17/248 K2 isolates is >25% away from theta.
    // Output-neutral: stats and one stderr line, printed on every construction.
    // Population = the nodes the chosen estimator actually weighed.
    const std::string est(estimator);
    const bool massWeighted = est == "weighted" || est == "eligible_weighted";
    const size_t minLen = static_cast<size_t>(std::max(graph.k(), 0)) * (est == "weighted" ? 1 : 2);
    size_t population = 0;
    for (const Unitig& u : graph.nodes) {
        if (u.deleted || u.seq.size() < minLen) continue;
        if (massWeighted && (!std::isfinite(u.coverage) || u.coverage <= 0.0)) continue;
        ++population;
    }
    stats_.theta = medianCoverage_;
    stats_.repeatThreshold = medianCoverage_ * kRepeatMultiplier;
    stats_.legacyMedian = legacyMedian;
    stats_.thetaEstimator = estimator;
    stats_.thetaPopulation = population;
    std::fprintf(stderr,
        "[resolver] run=1 theta=%.4f repeatThreshold=%.4f estimator=%s legacyMedian=%.4f population=%zu\n",
        stats_.theta, stats_.repeatThreshold, estimator, legacyMedian, population);
}

void PairedResolver::buildIndex() {
    size_t totalKmers = 0;
    for (const Unitig& u : g_.nodes) {
        if (!u.deleted && u.seq.size() >= static_cast<size_t>(kMap_)) {
            totalKmers += u.seq.size() - static_cast<size_t>(kMap_) + 1;
        }
    }
    index_.reserve(totalKmers * 2);

    for (uint32_t i = 0; i < g_.nodes.size(); ++i) {
        const Unitig& u = g_.nodes[i];
        if (u.deleted || u.seq.size() < static_cast<size_t>(kMap_)) continue;
        Kmer fwd = 0, rc = 0;
        int valid = 0;
        for (uint32_t p = 0; p < u.seq.size(); ++p) {
            int c = baseCode(u.seq[p]);
            if (c < 0) { valid = 0; continue; }
            fwd = pushBack(fwd, c, kMap_);
            rc = pushFrontRc(rc, c, kMap_);
            if (++valid < kMap_) continue;
            const uint32_t start = p + 1 - static_cast<uint32_t>(kMap_);
            const Kmer canon = fwd < rc ? fwd : rc;
            // strand flag records whether the unitig's forward k-mer is the
            // canonical one, so a read hit can be resolved to an orientation.
            const int strand = (fwd == canon) ? 0 : 1;
            auto it = index_.find(canon);
            if (it == index_.end()) index_.emplace(canon, packIndex(i, start, strand));
            else it->second = kAmbiguous;   // occurs elsewhere; cannot locate a read
        }
    }
}

Anchor PairedResolver::anchorRead(size_t read) const {
    Anchor best;
    const int len = static_cast<int>(reads_.length(read));
    if (len < kMap_) return best;

    const int span = len - kMap_;
    const int probes = std::min(kMaxProbes, span + 1);

    struct Vote { uint32_t unitig; int32_t pos; uint8_t orient; int count; };
    Vote votes[kMaxProbes];
    int distinct = 0;

    // Probes are spread across the read because a sequencing error invalidates
    // every k-mer overlapping it; one error must not veto the whole anchor.
    for (int t = 0; t < probes; ++t) {
        const int rp = probes == 1 ? 0 : span * t / (probes - 1);
        Kmer fwd = 0, rcv = 0;
        bool ok = true;
        for (int j = 0; j < kMap_; ++j) {
            int c = reads_.baseAt(read, static_cast<uint32_t>(rp + j));
            if (c < 0) { ok = false; break; }
            fwd = pushBack(fwd, c, kMap_);
            rcv = pushFrontRc(rcv, c, kMap_);
        }
        if (!ok) continue;

        const Kmer canon = fwd < rcv ? fwd : rcv;
        auto it = index_.find(canon);
        if (it == index_.end() || it->second == kAmbiguous) continue;

        const uint32_t u = idxUnitig(it->second);
        const int up = static_cast<int>(idxPos(it->second));
        const int sflag = idxStrand(it->second);
        const int rflag = (fwd == canon) ? 0 : 1;
        const int orient = rflag ^ sflag;

        // A read spanning a junction legitimately hangs off the unitig, so the
        // position is deliberately allowed to fall outside [0, unitigLen).
        const int32_t startPos = (orient == 0)
                                     ? static_cast<int32_t>(up - rp)
                                     : static_cast<int32_t>(up - (len - rp - kMap_));

        int found = -1;
        for (int q = 0; q < distinct; ++q) {
            if (votes[q].unitig == u && votes[q].pos == startPos &&
                votes[q].orient == static_cast<uint8_t>(orient)) { found = q; break; }
        }
        if (found >= 0) ++votes[found].count;
        else if (distinct < kMaxProbes) {
            votes[distinct++] = {u, startPos, static_cast<uint8_t>(orient), 1};
        }
    }

    int bestIdx = -1, bestVotes = 0;
    for (int q = 0; q < distinct; ++q) {
        if (votes[q].count > bestVotes) { bestVotes = votes[q].count; bestIdx = q; }
    }
    if (bestIdx < 0 || bestVotes < kMinVotes) return best;

    best.unitig = votes[bestIdx].unitig;
    best.pos = votes[bestIdx].pos;
    best.orient = votes[bestIdx].orient;
    return best;
}

void PairedResolver::buildSupport() {
    buildIndex();
    if (!reads_.paired()) return;

    routeDistance_ = env::on("TESSERACT_ROUTE_DISTANCE", false);
    routeDensity_ = {};
    const size_t pairs = reads_.pairCount();
    std::vector<std::vector<detail::RouteTrainingInterval>> routeIntervals(routeDistance_ ? g_.nodes.size() : 0);
    std::vector<detail::RouteTrainingSequence> routePopulation;
    size_t routeTrainingNodes = 0;
    if (routeDistance_) {
        std::vector<std::vector<uint8_t>> uniqueStarts(g_.nodes.size());
        for (size_t u = 0; u < g_.nodes.size(); ++u) {
            const auto& node = g_.nodes[u];
            if (!node.deleted && node.seq.size() >= size_t(2 * k_) &&
                std::isfinite(node.coverage) && node.coverage > 0 &&
                node.coverage <= kRepeatMultiplier * medianCoverage_)
                uniqueStarts[u].resize(node.seq.size() - size_t(kMap_) + 1, 0);
        }
        // Reuse the actual ambiguity-aware anchoring index. Nominal low depth
        // alone does not imply a mappable training opportunity.
        for (const auto& entry : index_) {
            if (entry.second == kAmbiguous) continue;
            const uint32_t u = idxUnitig(entry.second), pos = idxPos(entry.second);
            if (pos < uniqueStarts[u].size()) uniqueStarts[u][pos] = 1;
        }
        for (size_t u = 0; u < g_.nodes.size(); ++u) {
            routeIntervals[u] = detail::uniqueSeedIntervals(uniqueStarts[u], size_t(kMap_));
            if (!routeIntervals[u].empty()) ++routeTrainingNodes;
            for (const auto& interval : routeIntervals[u])
                routePopulation.push_back({interval.end - interval.begin, g_.nodes[u].coverage});
        }
    }
    std::vector<std::vector<int>> routeSamples(static_cast<size_t>(threads_));
    std::vector<std::vector<int>> insertSamples(static_cast<size_t>(threads_));
    std::vector<std::unordered_map<uint64_t, std::unordered_map<uint64_t, std::vector<int32_t>>>>
        localSupport(static_cast<size_t>(threads_));
    std::atomic<size_t> mappedCount{0}, linkCount{0};
    const bool keepAnchors = pairAnchoredPrefixSetting() != 0;
    if (keepAnchors) readAnchors_.assign(reads_.size(), Anchor{});

    auto worker = [&](int tid) {
        auto& samples = insertSamples[static_cast<size_t>(tid)];
        auto& sup = localSupport[static_cast<size_t>(tid)];
        size_t localMapped = 0, localLinks = 0;

        for (size_t p = static_cast<size_t>(tid); p < pairs; p += static_cast<size_t>(threads_)) {
            const size_t i1 = p * 2, i2 = p * 2 + 1;
            const Anchor a1 = anchorRead(i1);
            if (keepAnchors) readAnchors_[i1] = a1;
            if (!a1.mapped()) {
                if (keepAnchors) readAnchors_[i2] = anchorRead(i2);
                continue;
            }
            const Anchor a2 = anchorRead(i2);
            if (keepAnchors) readAnchors_[i2] = a2;
            if (!a2.mapped()) continue;
            localMapped += 2;

            const int len1 = static_cast<int>(reads_.length(i1));
            const int len2 = static_cast<int>(reads_.length(i2));

            if (a1.unitig == a2.unitig) {
                // Same unitig: a properly oriented pair measures the fragment
                // length directly, which is how the insert model is learned.
                if (a1.orient == 0 && a2.orient == 1 && a2.pos + len2 > a1.pos) {
                    samples.push_back(a2.pos + len2 - a1.pos);
                } else if (a2.orient == 0 && a1.orient == 1 && a1.pos + len1 > a2.pos) {
                    samples.push_back(a1.pos + len1 - a2.pos);
                }
                if (routeDistance_ && !routeIntervals[a1.unitig].empty()) {
                    const int nodeLength = int(g_.nodes[a1.unitig].seq.size());
                    // Both entire read footprints must belong to this unitig;
                    // overhanging anchors cannot measure contained opportunity.
                    if (a1.pos >= 0 && a2.pos >= 0 && a1.pos + len1 <= nodeLength &&
                        a2.pos + len2 <= nodeLength) {
                        int begin = -1, end = -1;
                        if (a1.orient == 0 && a2.orient == 1 && a2.pos >= a1.pos && a1.pos + len1 <= a2.pos + len2) {
                            begin = a1.pos; end = a2.pos + len2;
                        } else if (a2.orient == 0 && a1.orient == 1 && a1.pos >= a2.pos && a2.pos + len2 <= a1.pos + len1) {
                            begin = a2.pos; end = a1.pos + len1;
                        }
                        if (begin >= 0 && end - begin >= kMap_ + 1 &&
                            detail::containedTrainingFragment(routeIntervals[a1.unitig], size_t(begin), size_t(end)))
                            routeSamples[size_t(tid)].push_back(end - begin);
                    }
                }
                continue;
            }

            const int eu = (a1.orient == 0) ? 1 : 0;
            const int ev = (a2.orient == 1) ? 0 : 1;
            // Orientation of each unitig as traversed by the fragment.
            const int dU = (eu == 1) ? 0 : 1;
            const int dV = (ev == 0) ? 0 : 1;

            const int lenU = static_cast<int>(g_.nodes[a1.unitig].seq.size());
            const int lenV = static_cast<int>(g_.nodes[a2.unitig].seq.size());
            // Bases of the fragment that fall inside each unitig, measured from
            // the read's 5' end to the boundary the fragment crosses.
            const int consumedU = (a1.orient == 0) ? (lenU - a1.pos) : (a1.pos + len1);
            const int consumedV = (dV == 0) ? (a2.pos + len2) : (lenV - a2.pos);
            // Concatenating two unitigs merges their shared (k-1) overlap.
            const int32_t span = static_cast<int32_t>(consumedU + consumedV - (k_ - 1));

            sup[orientedId(a1.unitig, dU)][orientedId(a2.unitig, dV)].push_back(span);
            // The mirrored traversal must be recorded too, so a path arriving
            // from the other side sees the same evidence.
            sup[orientedId(a2.unitig, 1 - dV)][orientedId(a1.unitig, 1 - dU)].push_back(span);
            ++localLinks;
        }
        mappedCount += localMapped;
        linkCount += localLinks;
    };

    std::vector<std::thread> pool;
    for (int t = 0; t < threads_; ++t) pool.emplace_back(worker, t);
    for (auto& th : pool) th.join();

    for (auto& sup : localSupport) {
        for (auto& kv : sup) {
            auto& dst = support_[kv.first];
            for (auto& kv2 : kv.second) {
                auto& v = dst[kv2.first];
                v.insert(v.end(), kv2.second.begin(), kv2.second.end());
            }
        }
    }

    std::vector<int> all;
    for (auto& s : insertSamples) all.insert(all.end(), s.begin(), s.end());
    for (int v : all) {
        if (v < 0) continue;
        const size_t bin = static_cast<size_t>(v);
        if (bin >= insertHistogram_.size()) insertHistogram_.resize(bin + 1, 0);
        ++insertHistogram_[bin];
    }
    if (all.size() >= 1000) {
        std::sort(all.begin(), all.end());
        // Trim the tails before fitting: chimeric pairs and mismapped reads sit
        // far out and would inflate the standard deviation badly.
        const size_t lo = all.size() / 100;
        const size_t hi = all.size() - all.size() / 100;
        double sum = 0;
        for (size_t i = lo; i < hi; ++i) sum += all[i];
        const double n = static_cast<double>(hi - lo);
        insert_.mean = sum / n;
        double var = 0;
        for (size_t i = lo; i < hi; ++i) {
            const double d = all[i] - insert_.mean;
            var += d * d;
        }
        insert_.stddev = std::sqrt(var / n);
        insert_.observations = all.size();
        insert_.minPlausible = std::max(0, static_cast<int>(insert_.mean - 4 * insert_.stddev));
        insert_.maxPlausible = static_cast<int>(insert_.mean + 4 * insert_.stddev);
        insert_.usable = true;
    }
    // An externally measured window wins over the fitted one, and -- unlike the fit --
    // does not need 1000 same-unitig pairs to exist before it can be used. On a badly
    // fragmented graph those pairs are scarce precisely when the resolver needs the
    // constraint most, so this also rescues the case where no model could be fitted here
    // at all.
    if (haveForcedBounds_) {
        insert_.minPlausible = forcedMin_;
        insert_.maxPlausible = forcedMax_;
        insert_.usable = true;
    }

    if (routeDistance_ && insert_.usable) {
        std::vector<uint64_t> histogram;
        for (const auto& samples : routeSamples) for (int length : samples) {
            if (length < insert_.minPlausible || length > insert_.maxPlausible) continue;
            if (size_t(length) >= histogram.size()) histogram.resize(size_t(length) + 1, 0);
            ++histogram[size_t(length)];
        }
        const size_t trainingIntervals = routePopulation.size();
        routeDensity_ = detail::fitRouteInsertDensity(histogram, std::move(routePopulation),
                                                      std::max(insert_.minPlausible, kMap_ + 1), insert_.maxPlausible);
        std::fprintf(stderr,
            "[routedensity] nodes=%zu intervals=%zu pairs=%llu effective=%.1f bandwidth=%d usable=%d bounds=%d:%d\n",
            routeTrainingNodes, trainingIntervals, static_cast<unsigned long long>(routeDensity_.observations),
            routeDensity_.effectiveObservations, routeDensity_.bandwidth, routeDensity_.usable ? 1 : 0,
            std::max(insert_.minPlausible, kMap_ + 1), insert_.maxPlausible);
    }

    size_t distinct = 0;
    for (auto& kv : support_) distinct += kv.second.size();
    stats_.readsMapped = mappedCount.load();
    stats_.pairsLinking = linkCount.load();
    stats_.distinctLinks = distinct;
    stats_.insert = insert_;
}


double PairedResolver::scoreCandidate(uint64_t from, uint64_t terminal, int interLen) const {
    auto it = support_.find(from);
    if (it == support_.end()) return 0;
    auto jt = it->second.find(terminal);
    if (jt == it->second.end()) return 0;
    if (!insert_.usable) return static_cast<double>(jt->second.size());

    // A pair supports this particular path only if the fragment length it
    // implies once the path's intermediate sequence is inserted is one the
    // library could actually have produced.
    double score = 0;
    for (int32_t span : jt->second) {
        const int implied = span + interLen;
        if (implied >= insert_.minPlausible && implied <= insert_.maxPlausible) score += 1;
    }
    return score;
}

void PairedResolver::resolve(std::vector<std::string>& contigs, std::vector<double>& covs) {
    const size_t n = g_.nodes.size();
    const size_t ov = static_cast<size_t>(k_ - 1);

    // build_v3 G-resolve fixes: all default off (see fixEnabled). Each prints its counter
    // line at the end of this function on every run, zero included.
    const bool fixGapFlank = fixEnabled("TESSERACT_FIX_GAP_FLANK");            // T03
    const bool fixRevisitGuard = fixEnabled("TESSERACT_FIX_REVISIT_GUARD");    // T07 refuse
    const bool fixCycleBreak = fixEnabled("TESSERACT_FIX_SCAFFOLD_CYCLE");     // T13
    const bool fixGapEstimate = fixEnabled("TESSERACT_FIX_GAP_ESTIMATE");      // T20
    const bool fixTruncGuard = fixEnabled("TESSERACT_FIX_TRUNC_GUARD");        // T21
    const bool fixCovContrib = fixEnabled("TESSERACT_FIX_COV_CONTRIB");        // T24
    const bool fixRouteOrder = fixEnabled("TESSERACT_FIX_ROUTE_ORDER");        // T31
    const bool fixMirrorRoute = fixEnabled("TESSERACT_FIX_MIRROR_ROUTE");      // R6 (housekeeping)
    gapFlank_.clear();

    // A collapsed repeat may be traversed more than once and carries no unique
    // paired evidence, so it cannot seed a chain, terminate a path, or be
    // traversed -- it is emitted verbatim and walls off both neighbours.
    //
    // Depth alone does not identify one. A plasmid at two copies per cell sits
    // at the same ~2x median as a two-copy repeat while its k-mers are
    // perfectly unique, and a bacterial isolate usually carries several
    // plasmids; on one panel isolate 4.8% of the assembly was excluded this
    // way, including a 68 kb unitig at 1.78x median that was the whole of that
    // assembly's NGA50. Raising the depth cutoff only trades that error for the
    // opposite one, and stops two-copy repeats resolving.
    //
    // Two refinements were tried and measured on the closed-reference panel,
    // and both were rejected. Raising the multiplier to 2.4 frees the plasmids
    // but stops genuine two-copy repeats resolving. Additionally requiring a
    // branching end -- on the theory that distinct copies of a repeat are
    // flanked by different sequence -- left contig NGA50 unchanged to the base
    // pair on both datasets tested while adding misassemblies (3 -> 5 and
    // 2 -> 3), because a collapsed tandem array has a single link at each end
    // and was then traversed once instead of many times.
    //
    // So this stays depth-only. The plasmid cost is real but is the cheaper
    // error, and the contiguity it was meant to buy came from the scaffolding
    // and tie-break changes instead.
    const double repeatThreshold = medianCoverage_ * kRepeatMultiplier;
    // Vouching is done over a CONNECTED GROUP of high-depth unitigs, not one unitig at a
    // time, because the marker sampling will not support the per-unitig question. The
    // model samples one canonical k-mer in `markerDenom_` -- 512 for the shipped S. aureus
    // model -- so a 1,000 bp unitig contributes about two sampled k-mers in total, and
    // asking it for two plasmid-exclusive markers asks for near enough all of them. The
    // first version of this rule did exactly that and vouched nothing on either isolate
    // tested.
    //
    // A multi-copy plasmid is anyway not one deep unitig: it is a group of them, joined to
    // each other and cut off from the chromosome by the same depth test. Summing markers
    // over the group both matches that structure and gives the sampling enough sequence to
    // work with -- a 20 kb plasmid contributes about 40 sampled k-mers however many pieces
    // it is in.
    std::vector<char> plasmidVouched(n, 0);
    size_t vouchedCount = 0, vouchedGroups = 0;
    if (exclusiveMarkers_ && !exclusiveMarkers_->empty() && medianCoverage_ > 0) {
        std::vector<char> deep(n, 0);
        for (uint32_t u = 0; u < n; ++u) {
            if (!g_.nodes[u].deleted && g_.nodes[u].coverage > repeatThreshold) deep[u] = 1;
        }
        // Connected groups among the deep unitigs only. Union-find would do; the graph is
        // small enough here that a flood fill is clearer.
        std::vector<uint32_t> group(n, UINT32_MAX);
        std::vector<std::vector<uint32_t>> groups;
        std::vector<uint32_t> stack;
        for (uint32_t seed = 0; seed < n; ++seed) {
            if (!deep[seed] || group[seed] != UINT32_MAX) continue;
            const uint32_t gid = static_cast<uint32_t>(groups.size());
            groups.emplace_back();
            stack.assign(1, seed);
            group[seed] = gid;
            while (!stack.empty()) {
                const uint32_t u = stack.back(); stack.pop_back();
                groups[gid].push_back(u);
                for (int o = 0; o < 2; ++o) {
                    for (const Link& l : g_.exits(u, o)) {
                        if (g_.nodes[l.to].deleted || !deep[l.to]) continue;
                        if (group[l.to] != UINT32_MAX) continue;
                        group[l.to] = gid;
                        stack.push_back(l.to);
                    }
                }
            }
        }
        for (const std::vector<uint32_t>& grp : groups) {
            size_t len = 0;
            uint32_t pls = 0, chr = 0;
            for (uint32_t u : grp) {
                len += g_.nodes[u].seq.size();
                forEachMarkerKmer(g_.nodes[u].seq, [&](uint64_t km, uint32_t, int) {
                    auto it = exclusiveMarkers_->find(km);
                    if (it == exclusiveMarkers_->end()) return;
                    if (it->second == 1) ++pls; else ++chr;
                }, markerDenom_);
            }
            // One-sided on purpose. A wrong exemption lets a real repeat be traversed, so a
            // single chromosome-exclusive marker anywhere in the group vetoes the whole
            // group, and a repeat family shared between chromosome and plasmid -- which
            // contributes no exclusive markers either way -- is never exempted.
            if (chr != 0 || pls < kMinPlasmidVouch || len < kMinVouchLen) continue;
            ++vouchedGroups;
            for (uint32_t u : grp) { plasmidVouched[u] = 1; ++vouchedCount; }
        }
    }
    if (env::present("TESSERACT_DEBUG_RESOLVE") && exclusiveMarkers_) {
        size_t hi = 0, hiLen = 0, hiMax = 0;
        for (uint32_t u = 0; u < n; ++u) {
            if (g_.nodes[u].deleted || medianCoverage_ <= 0) continue;
            if (g_.nodes[u].coverage <= repeatThreshold) continue;
            ++hi; hiLen += g_.nodes[u].seq.size();
            hiMax = std::max(hiMax, g_.nodes[u].seq.size());
        }
        std::fprintf(stderr,
            "      [debug] high-depth unitigs=%zu (mean len %zu, max %zu); vouched %zu in "
            "%zu group(s) (denom %u, exclusive markers %zu)\n",
            hi, hi ? hiLen / hi : 0, hiMax, vouchedCount, vouchedGroups, markerDenom_,
            exclusiveMarkers_->size());
    }

    const bool noPlasmidVouch = env::present("TESSERACT_NO_PLASMID_VOUCH");
    auto isRepeat = [&](uint32_t u) {
        if (medianCoverage_ <= 0) return false;
        if (g_.nodes[u].coverage <= repeatThreshold) return false;
        if (plasmidVouched[u] && !noPlasmidVouch) return false;
        return true;
    };

    // Reads only anchor uniquely inside a unitig long enough to own k-mers its
    // neighbours do not share across the (k-1) junction overlaps. Shorter
    // unitigs carry no paired evidence, so extension has to look *through* them
    // to the next unitig that does.
    const size_t anchorableLen = static_cast<size_t>(2 * k_);
    auto anchorable = [&](uint32_t u) {
        return !g_.nodes[u].deleted && g_.nodes[u].seq.size() >= anchorableLen && !isRepeat(u);
    };

    const double reach = insert_.usable ? insert_.maxPlausible : 1000.0;
    constexpr size_t kMaxNodes = 10;
    constexpr double kForcedReachFactor = 1.5;
    constexpr size_t kMaxCandidates = 64;

    // How far enumeration may look when it is only asking *what continuations
    // exist*, as opposed to asking the read pairs to choose between them.
    //
    // Tying the search itself to fragment reach conflates the two. On libraries
    // whose mates overlap -- most of the closed-reference panel, where the
    // fragment is ~230 bp against ~190 bp reads -- reach is a few hundred bases,
    // so any repeat longer than that ends the search before a single candidate
    // is generated: on one isolate 2,620 of 3,288 chain ends returned no
    // candidate at all, against 269 successful joins. The paired evidence never
    // got the chance to be insufficient; it was never consulted.
    //
    // Searching further does not weaken any decision. A candidate beyond
    // fragment reach collects no paired support and so cannot win a contested
    // choice -- it can only be taken when it is the sole way through, which is
    // a statement about the graph's topology that needs no reads to back it.
    const double searchBudget = std::max(reach, kTopologyReach);

    // How much paired support a contested join needs. A flat count cannot mean
    // the same thing at every depth: the number of pairs crossing a junction
    // scales with coverage, so a fixed 2 is a real bar on a 40x library and
    // almost none on a 100x one. That showed up as soon as the trimming stopped
    // discarding a fifth of the bases -- contig NGA50 rose sharply and so did
    // long-range chimeric joins, which is what a threshold too weak for the
    // depth looks like. Scaling keeps the bar constant in the units that
    // matter, with the configured value as the floor so a shallow library is
    // never asked for less than it was before.
    const double linkBar = std::max(static_cast<double>(minLinkSupport_),
                                    std::min(6.0, medianCoverage_ * linkSupportPerX_));

    // ON by default since 1.3.0. Paired against the 1.3.0 resolver on 146 isolates: NGA50
    // 27 better / 1 worse, contigs 65 / 4, genome fraction 39 / 10, size accuracy 37 / 15;
    // misassemblies 1 better / 4 worse (+5 events net). TESSERACT_EXACT_READ_THREADS=0 disables.
    const bool exactReadThreads = reads_.paired() && env::on("TESSERACT_EXACT_READ_THREADS", true);
    ReadThreadEvidence threadEvidence;
    // Preserve observed routes even if all their molecules already supplied
    // paired evidence: contrary exact observations must remain visible vetoes.
    std::vector<size_t> freshThreadSupport;
    std::unordered_map<uint64_t, std::vector<std::pair<size_t, bool>>> threadPorts;
    size_t threadPairedExcluded = 0, threadFreshFragments = 0;
    size_t threadNominations = 0, threadConflictingPorts = 0, threadJoins = 0;
    if (exactReadThreads) {
        std::vector<uint8_t> eligible(n, 0);
        for (uint32_t u = 0; u < n; ++u) eligible[u] = anchorable(u);
        threadEvidence = collectReadThreadEvidence(g_, reads_, eligible);
        freshThreadSupport.resize(threadEvidence.routes.size(), 0);
        for (size_t i = 0; i < threadEvidence.routes.size(); ++i) {
            const auto& route = threadEvidence.routes[i];
            for (size_t fragment : route.fragments) {
                bool alreadyPaired = false;
                if (fragment < reads_.pairCount()) {
                    const Anchor first = anchorRead(fragment * 2);
                    const Anchor second = anchorRead(fragment * 2 + 1);
                    if (first.mapped() && second.mapped() && first.unitig != second.unitig) {
                        const uint64_t from = orientedId(first.unitig, first.orient);
                        const uint64_t to = orientedId(second.unitig, 1 - second.orient);
                        alreadyPaired = (from == route.oriented.front() && to == route.oriented.back()) ||
                            ((from ^ 1) == route.oriented.back() && (to ^ 1) == route.oriented.front());
                    }
                }
                if (alreadyPaired) ++threadPairedExcluded;
                else { ++freshThreadSupport[i]; ++threadFreshFragments; }
            }
            threadPorts[route.oriented.front()].push_back({i, false});
            threadPorts[route.oriented.back() ^ 1].push_back({i, true});
        }
        const auto& stats = threadEvidence.stats;
        std::fprintf(stderr,
            "[readthread] reads=%zu acceptedReads=%zu molecules=%zu fresh=%zu pairedExcluded=%zu "
            "routes=%zu unknown=%zu noFlanks=%zu noExact=%zu ambiguous=%zu searchLimited=%zu "
            "mateDuplicates=%zu conflictingMates=%zu\n",
            stats.readsExamined, stats.readsAccepted, stats.fragmentsAccepted, threadFreshFragments,
            threadPairedExcluded, threadEvidence.routes.size(), stats.readsUnknown, stats.readsNoFlanks,
            stats.readsNoExactPath, stats.readsAmbiguous, stats.readsSearchLimited,
            stats.matesDeduplicated, stats.fragmentsConflicting);
    }

    auto flip = [](uint64_t oid) { return orientedId(unitigOf(oid), 1 - orientOf(oid)); };
    auto addedLen = [&](uint64_t oid) { return g_.nodes[unitigOf(oid)].seq.size() - ov; };

    // Depth of a chain, weighted by how much sequence each unitig contributes,
    // and only over the part near the end being extended -- a 400 kb chain's
    // far end says nothing about the coverage where it is about to continue.
    auto chainDepth = [&](const std::vector<uint64_t>& tailFirst) {
        double num = 0, den = 0, acc = 0;
        for (size_t i = tailFirst.size(); i-- > 0;) {
            const uint32_t u = unitigOf(tailFirst[i]);
            const double len = static_cast<double>(addedLen(tailFirst[i]));
            num += g_.nodes[u].coverage * len;
            den += len;
            acc += len;
            if (acc > kDepthWindow) break;
        }
        return den > 0 ? num / den : 0.0;
    };

    // True when every candidate ends on the same oriented unitig, so the only
    // thing in doubt is which way through the repeat, not where it comes out.
    auto allSameDestination = [](const std::vector<std::vector<uint64_t>>& cands) {
        if (cands.empty()) return false;
        const uint64_t t = cands.front().back();
        for (const auto& c : cands) {
            if (c.back() != t) return false;
        }
        return true;
    };

    // Among routes to the same destination, the one whose weakest interior
    // unitig is best covered -- a route through sequence the reads support all
    // the way is likelier than one that dips through a thinly covered branch.
    auto pickByInterior = [&](const std::vector<std::vector<uint64_t>>& cands) {
        int pickIdx = 0;
        double bestInterior = -1.0;
        for (size_t i = 0; i < cands.size(); ++i) {
            double minCov = std::numeric_limits<double>::max();
            for (size_t j = 0; j + 1 < cands[i].size(); ++j) {
                minCov = std::min(minCov, g_.nodes[unitigOf(cands[i][j])].coverage);
            }
            // A path with no interior is the most direct route.
            if (minCov > bestInterior) { bestInterior = minCov; pickIdx = static_cast<int>(i); }
        }
        return pickIdx;
    };

    // The candidate whose terminal unitig runs at the chain's own depth, when
    // exactly one does and the rest are clearly off it. Returns -1 when the
    // answer is not clean enough to act on -- which is most of the time, and
    // deliberately so: this fires only where the paired reads had nothing to
    // say, so it has no second opinion to check itself against.
    auto pickByCoverage = [&](const std::vector<uint64_t>& tailFirst,
                              const std::vector<std::vector<uint64_t>>& cands) -> int {
        const double depth = chainDepth(tailFirst);
        if (depth <= 0) return -1;
        int match = -1;
        double worstOther = 1e9;
        for (size_t i = 0; i < cands.size(); ++i) {
            const double cov = g_.nodes[unitigOf(cands[i].back())].coverage;
            const double rel = std::fabs(cov - depth) / depth;
            if (rel <= kDepthMatch) {
                if (match >= 0) return -1;      // two candidates fit; no answer
                match = static_cast<int>(i);
            } else {
                worstOther = std::min(worstOther, rel);
            }
        }
        if (match < 0) return -1;
        // The runner-up has to be properly off-depth, not just outside the
        // matching band, or this is a tie dressed up as a decision.
        if (cands.size() > 1 && worstOther < kDepthSeparation) return -1;
        return match;
    };

    // What enumerate() discarded without completing (T21). A discarded frame is a
    // route whose destination was never learned, not proof that none exists; the
    // "sole way through" argument above holds only when nothing was discarded.
    // Dead ends and revisits are not truncation.
    struct EnumTrunc {
        size_t budget = 0, maxNodes = 0, cap = 0;
        bool any() const { return budget || maxNodes || cap; }
    };
    // Revisits enumerate() pruned (T07): the frame reached `prefix.back()` again,
    // i.e. the graph offers a cycle through it that no enumerated route traverses
    // more than once.
    struct EnumLoops {
        std::vector<std::vector<uint64_t>> prefix;
    };

    // Every way out of `from` that terminates on an anchorable unitig within
    // fragment reach; intermediate nodes are unanchorable repeats.
    auto enumerate = [&](uint64_t from, std::vector<std::vector<uint64_t>>& out,
                         EnumTrunc* trunc = nullptr, EnumLoops* loops = nullptr) {
        // `forced` marks a path that has never had an alternative. Paired reads
        // cannot vouch for anything beyond one fragment length, but where the
        // graph offers no other way through, length is irrelevant -- so a forced
        // path may run past `reach` and cross repeats no pair could span.
        struct Frame { std::vector<uint64_t> nodes; size_t len; bool forced; };
        std::vector<Frame> work;
        const auto& firstExits = g_.exits(unitigOf(from), orientOf(from));
        for (const Link& l : firstExits) {
            if (g_.nodes[l.to].deleted) continue;
            work.push_back({{orientedId(l.to, UnitigGraph::enterOrient(l))},
                            g_.nodes[l.to].seq.size() - ov, firstExits.size() == 1});
        }
        while (!work.empty() && out.size() < kMaxCandidates) {
            Frame f = std::move(work.back());
            work.pop_back();
            const uint64_t tail = f.nodes.back();
            const uint32_t tu = unitigOf(tail);
            if (anchorable(tu)) { out.push_back(std::move(f.nodes)); continue; }
            if (f.nodes.size() >= kMaxNodes) { if (trunc) ++trunc->maxNodes; continue; }
            // A forced path (one the graph offered no alternative to) may run a
            // little past what a fragment can span, but not far: measured on the
            // benchmark panel, letting it run unbounded bought only a lower
            // contig count and cost a misassembly on K. pneumoniae.
            const double budget = f.forced ? searchBudget * kForcedReachFactor : searchBudget;
            if (static_cast<double>(f.len) > budget) { if (trunc) ++trunc->budget; continue; }
            const auto& exits = g_.exits(tu, orientOf(tail));
            for (const Link& l : exits) {
                if (g_.nodes[l.to].deleted) continue;
                const uint64_t nxt = orientedId(l.to, UnitigGraph::enterOrient(l));
                size_t seenAt = SIZE_MAX;
                for (size_t q = 0; q < f.nodes.size(); ++q) if (f.nodes[q] == nxt) { seenAt = q; break; }
                if (seenAt != SIZE_MAX) {
                    if (loops)
                        loops->prefix.emplace_back(f.nodes.begin(), f.nodes.begin() + static_cast<std::ptrdiff_t>(seenAt) + 1);
                    continue;
                }
                Frame g2 = f;
                g2.nodes.push_back(nxt);
                g2.len += g_.nodes[l.to].seq.size() - ov;
                g2.forced = f.forced && exits.size() == 1;
                work.push_back(std::move(g2));
            }
        }
        if (trunc && out.size() >= kMaxCandidates) trunc->cap += work.size();
    };

    // Where SPAdes was ahead, the junction is a short branch node -- 219 to
    // 398 bp, in-degree 2 and out-degree 2, confirmed from the emitted GFA --
    // sitting between two long chains. No single read pair can span it: with a
    // 355 +/- 119 fragment and ~200 bp reads, crossing a 308 bp node needs
    // read + node + read, about 708 bp, three standard deviations out. Each
    // junction therefore collects one to three supporting pairs against a bar
    // of four and is refused, and the sequence leaves in the output as a
    // contig of its own.
    //
    // The bar is not wrong; the evidence being weighed is incomplete. Two ways
    // in and two ways out means both traversals of the node are used, which
    // turns four independent guesses into a choice between two perfect
    // matchings: if our tail pairs with terminal S, the node's *other*
    // predecessor must pair with the other successor. Scoring both assignments
    // and requiring the intended one to dominate is much stronger than asking
    // a single join to clear a threshold alone, because a wrong pairing has to
    // beat the right one twice over.
    //
    // Measured against the depth-scaled bar on ten isolates: mean NGA50
    // +6,924 and genome fraction +0.02 pp for one extra misassembly, which is
    // roughly 69,000 bases of contiguity per misassembly spent. Simply
    // lowering the bar buys about 2,900. The rule is kept for that ratio, not
    // for the size of the gain, which is modest.
    //
    // Returns true when the matching is decisive. `interLen` is the sequence
    // the connector adds, as scoreCandidate counts it.
    // Sequence a connector adds ahead of its terminal, counted as
    // scoreCandidate expects.
    auto interLenOf = [&](const std::vector<uint64_t>& cand) {
        int len = 0;
        for (size_t j = 0; j + 1 < cand.size(); ++j) {
            len += static_cast<int>(g_.nodes[unitigOf(cand[j])].seq.size()) - (k_ - 1);
        }
        return len;
    };

    const double dom = env::real("TESSERACT_MATCH_DOMINANCE", kMatchDominance);
    auto matchingAgrees = [&](uint64_t from, const std::vector<uint64_t>& connector,
                              uint64_t terminal, int interLen) {
        if (connector.size() != 2) return false;      // one repeat between, then the terminal
        const uint64_t rep = connector.front();
        const uint32_t ru = unitigOf(rep);
        if (!isRepeat(ru)) return false;
        if (g_.nodes[ru].seq.size() > kMaxMatchedRepeat) return false;

        // The repeat's two ways in and two ways out, in the frame it is
        // traversed. Anything other than exactly two of each is a different
        // problem and is left alone.
        const auto& outs = g_.exits(ru, orientOf(rep));
        const auto& ins = g_.exits(ru, 1 - orientOf(rep));
        if (outs.size() != 2 || ins.size() != 2) return false;

        uint64_t succ[2], pred[2];
        for (int i = 0; i < 2; ++i) {
            succ[i] = orientedId(outs[i].to, UnitigGraph::enterOrient(outs[i]));
            pred[i] = flip(orientedId(ins[i].to, UnitigGraph::enterOrient(ins[i])));
        }
        // Identify which way in is ours and which way out we are proposing.
        int mine = -1, want = -1;
        for (int i = 0; i < 2; ++i) {
            if (pred[i] == from) mine = i;
            if (succ[i] == terminal) want = i;
        }
        if (mine < 0 || want < 0) return false;
        const uint64_t other = pred[1 - mine];
        const uint64_t alt = succ[1 - want];
        if (other == from || alt == terminal) return false;

        const double matched = scoreCandidate(from, terminal, interLen) +
                               scoreCandidate(other, alt, interLen);
        const double crossed = scoreCandidate(from, alt, interLen) +
                               scoreCandidate(other, terminal, interLen);
        return matched >= 1.0 && matched >= dom * crossed;
    };

    // Complete terminal repeats across a single-substitution bubble only when
    // explicitly requested. This extends sequence; it does not relax the anchor
    // join chooser or make a choice between distinct flanking destinations.
    const bool completePrefixBubbles = env::on("TESSERACT_PREFIX_SNP_BUBBLES", false);

    // Start with every anchorable unitig as a chain of one, then repeatedly
    // join chains whose paired evidence mutually prefers each other. Growing
    // chains lets support accumulate over a whole contig tail rather than just
    // the last unitig, which is what carries a walk past a repeat.
    std::vector<std::vector<uint64_t>> chains;
    for (uint32_t u = 0; u < n; ++u) {
        if (anchorable(u)) chains.push_back({orientedId(u, 0)});
    }

    std::vector<uint32_t> ownerChain(n, UINT32_MAX);
    std::vector<uint32_t> ownerPos(n, 0);
    auto reindex = [&]() {
        std::fill(ownerChain.begin(), ownerChain.end(), UINT32_MAX);
        for (uint32_t c = 0; c < chains.size(); ++c) {
            for (uint32_t j = 0; j < chains[c].size(); ++j) {
                const uint32_t u = unitigOf(chains[c][j]);
                ownerChain[u] = c;
                ownerPos[u] = j;
            }
        }
    };

    struct Cont {
        uint32_t chainB = UINT32_MAX;
        int endB = 0;
        std::vector<uint64_t> connector;
        std::vector<uint64_t> legacyConnector;
        bool distanceChanged = false;
        bool exactThread = false;
        double score = 0;
        bool ok = false;
        // TESSERACT_ONE_SIDED_TIE (default off): `weak` marks a nomination made from
        // inside the tie zone, which may only join a `decisive` partner (a lone
        // candidate with best >= linkBar, or several candidates and no rival support).
        bool weak = false;
        bool decisive = false;
    };

    // Best continuation off one end of a chain, scored with every member of
    // that chain that still lies within fragment reach of the boundary.
    long long dbgNoCand = 0, dbgLowSupport = 0, dbgTie = 0, dbgMidChain = 0, dbgOk = 0;
    long long dbgCoverage = 0, dbgMatched = 0, dbgShortDest = 0, dbgUnspanned = 0;
    long long dbgLoneRepeat = 0;
    // T14 (always printed): what the shared-repeat exclusion did to the decisions.
    long long dbgEsrsAllExcluded = 0, dbgEsrsChanged = 0, dbgEvals = 0;
    // T21 (always printed): evaluations whose enumeration discarded frames.
    long long trSides = 0, trBudget = 0, trMaxNodes = 0, trCap = 0, trLone = 0;
    long long trLoneBelowBar = 0, trRefused = 0;
    // T07 (always printed): evaluations with a pruned revisit, picks that pass a pruned
    // loop, and refusals.
    long long rvPrunedSides = 0, rvPicked = 0, rvRefused = 0;
    // T31 (always printed): route-distance decisions whose spans were put in canonical order.
    long long roDecisions = 0, roSorted = 0;
    // R6 housekeeping (always printed): mutual legacy joins whose two ends chose different
    // interior routes.
    long long mrDisagree = 0, mrChoseBack = 0;
    const bool requireSupportSingle_ = env::on("TESSERACT_REQUIRE_SUPPORT_SINGLE", false);
    const bool joinTrace_ = env::on("TESSERACT_JOIN_TRACE", false);
    const bool noUnspannedFallback_ = env::on("TESSERACT_NO_UNSPANNED_FALLBACK", false);
    long long dbgTieBest = 0, dbgTieSecond = 0;

    // Withdraw the fallbacks when the continuation they would take lands on a unitig
    // shorter than this. 0 disables, which is the shipped behaviour until measured.
    const size_t minFallbackDest_ = static_cast<size_t>(env::integer("TESSERACT_MIN_FALLBACK_DEST", 0));

    const bool excludeSharedRepeatSupport = env::on("TESSERACT_EXCLUDE_SHARED_REPEAT_SUPPORT", false);
    const bool auditSharedRepeatSupport = env::on("TESSERACT_SHARED_SUPPORT_AUDIT", false);

    // Default-off experiment (combo2/joinrule, L4): a near-tied pick is nominated
    // "weak" instead of refused, and joins only when the reverse nomination had no
    // rival; long-long joins at a branching chain end stay refused.  Measured on
    // reference-labelled traces, the tie-zone ratio carries no error signal while
    // the partner's lack of a rival does (one-sided ~1.3% wrong, two-sided ~25%).
    const bool oneSidedTie = env::on("TESSERACT_ONE_SIDED_TIE", false);
    const size_t oneSidedLongBp = static_cast<size_t>(env::integer("TESSERACT_ONE_SIDED_TIE_LONG", 2000));
    // Default-off experiment (FB-M): of the three fallbacks keep only the perfect
    // matching (matchingAgrees); withdraw same-destination and coverage picks.
    // Default-off experiment (XO1, combo2/joinrule L6): a near-tied pick at a chain
    // end whose unitig has ONE live exit nominates normally.  The candidates then
    // diverge downstream of a connector, so every chain member's pairs are specific
    // to one terminal; at a branching end the end's own pairs reach every terminal.
    const bool tieSingleExit = env::on("TESSERACT_TIE_SINGLE_EXIT", false);
    // Length gate for XO1 (default: none).  L7 = XO1 only when the shorter of the
    // two chains is below this, weak nominations only below ONE_SIDED_TIE_LONG.
    const long long tieSingleExitLongSet = env::integer("TESSERACT_TIE_SINGLE_EXIT_LONG", -1);   // -1: unset
    const size_t tieSingleExitLong = tieSingleExitLongSet < 0 ? std::numeric_limits<size_t>::max()
                                                              : static_cast<size_t>(tieSingleExitLongSet);
    long long dbgSingleExit = 0, dbgSingleExitLong = 0;
    const bool fallbackMatchingOnly = env::on("TESSERACT_FALLBACK_MATCHING_ONLY", false);
    long long dbgWeakNom = 0, dbgWeakLongBranch = 0, dbgWeakTwoSided = 0, dbgWeakPartner = 0;
    long long dbgWeakJoin = 0, dbgWeakThread = 0, dbgFbmKept = 0, dbgFbmRefused = 0;
    long long dbgJtExt = 0, dbgFbTraced = 0, dbgFbTracedRefused = 0;   // JOIN_TRACE extensions

    auto bestContinuation = [&](uint32_t c, int end) {
        Cont result;
        const std::vector<uint64_t>& ch = chains[c];
        if (ch.empty()) return result;

        // Orient the chain so the end under consideration is the tail.
        std::vector<uint64_t> tailFirst;
        if (end == 1) tailFirst.assign(ch.begin(), ch.end());
        else {
            tailFirst.reserve(ch.size());
            for (size_t i = ch.size(); i-- > 0;) tailFirst.push_back(flip(ch[i]));
        }

        // Distance from the far edge of each member to the chain boundary.
        std::vector<size_t> distToEnd(tailFirst.size(), 0);
        size_t acc = 0;
        for (size_t i = tailFirst.size(); i-- > 0;) {
            distToEnd[i] = acc;
            acc += addedLen(tailFirst[i]);
        }

        std::vector<std::vector<uint64_t>> cands;
        EnumTrunc trunc;
        EnumLoops loops;
        enumerate(tailFirst.back(), cands, &trunc, &loops);
        ++dbgEvals;
        if (trunc.any()) {
            ++trSides;
            if (trunc.budget) ++trBudget;
            if (trunc.maxNodes) ++trMaxNodes;
            if (trunc.cap) ++trCap;
            if (cands.size() == 1) ++trLone;
        }
        if (!loops.prefix.empty()) ++rvPrunedSides;
        if (cands.empty()) { ++dbgNoCand; return result; }

        double best = -1, second = -1;
        size_t pick = 0;
        std::vector<double> scores(cands.size(), 0.0);
        for (size_t i = 0; i < cands.size(); ++i) {
            int interLen = 0;
            for (size_t j = 0; j + 1 < cands[i].size(); ++j) {
                interLen += static_cast<int>(g_.nodes[unitigOf(cands[i][j])].seq.size()) - (k_ - 1);
            }
            const uint64_t terminal = cands[i].back();
            double sc = 0;
            for (size_t m = tailFirst.size(); m-- > 0;) {
                if (static_cast<double>(distToEnd[m]) > reach) break;
                sc += scoreCandidate(tailFirst[m], terminal,
                                     interLen + static_cast<int>(distToEnd[m]));
            }
            scores[i] = sc;
            if (sc > best) { second = best; best = sc; pick = i; }
            else if (sc > second) second = sc;
        }

        if (excludeSharedRepeatSupport || auditSharedRepeatSupport) {
            // Compare a source's support across distinct destinations before
            // summing it into a branch score. This follows SPAdes' exclusion
            // of common history evidence, restricted here to collapsed repeats
            // so distinctive low-depth anchor support remains unchanged.
            std::vector<uint64_t> destinations;
            for (const auto& candidate : cands) destinations.push_back(candidate.back());
            std::vector<detail::BranchEvidence> evidence;
            for (size_t m = tailFirst.size(); m-- > 0;) {
                if (static_cast<double>(distToEnd[m]) > reach) break;
                detail::BranchEvidence row;
                row.source = tailFirst[m];
                row.repeat = isRepeat(unitigOf(tailFirst[m]));
                for (size_t i = 0; i < cands.size(); ++i) {
                    row.scores.push_back(scoreCandidate(tailFirst[m], destinations[i],
                        interLenOf(cands[i]) + static_cast<int>(distToEnd[m])));
                }
                evidence.push_back(std::move(row));
            }
            const auto adjusted = detail::distinctiveBranchScores(destinations, evidence);
            if (auditSharedRepeatSupport || joinTrace_) {
                const auto shared = detail::sharedRepeatSources(destinations, evidence);
                for (size_t m = 0; m < evidence.size(); ++m) {
                    if (!shared[m]) continue;
                    std::fprintf(stderr,
                        "[sharedsupport] chain=%u end=%d source=%llu active=%d destinations=",
                        c, end, static_cast<unsigned long long>(evidence[m].source),
                        excludeSharedRepeatSupport ? 1 : 0);
                    for (size_t i = 0; i < cands.size(); ++i)
                        std::fprintf(stderr, "%s%llu:%.1f:%.1f", i ? "," : "",
                            static_cast<unsigned long long>(destinations[i]),
                            evidence[m].scores[i], adjusted[i]);
                    std::fprintf(stderr, "\n");
                }
            }
            if (excludeSharedRepeatSupport) {
                const double rawBest = best;
                const size_t rawPick = pick;
                scores = adjusted;
                best = second = -1.0;
                pick = 0;
                for (size_t i = 0; i < scores.size(); ++i) {
                    const double sc = scores[i];
                    if (sc > best) { second = best; best = sc; pick = i; }
                    else if (sc > second) second = sc;
                }
                // T14 counters only; no decision reads them.
                if (rawBest > 0 && best <= 0) ++dbgEsrsAllExcluded;
                if (pick != rawPick || ((rawBest >= linkBar) != (best >= linkBar))) ++dbgEsrsChanged;
            }
        }

        // Gating this bar on whether the route passes through a repeat was
        // tried and does not work. The reasoning was sound -- every misassembly
        // the bar was introduced to stop was a relocation through a repeat, and
        // a continuation that reaches its destination through unique sequence
        // has no other copy to be confused with, so charging it the same bar
        // should buy nothing. Measured across eight isolates it bought a mean
        // +15,153 of NGA50 and four extra misassemblies, which is no better
        // than simply lowering the factor and is worse than leaving it alone.
        //
        // The reason is the proxy, not the idea: isRepeat is depth-only
        // (coverage > 1.6x median), and the unitigs that actually carry a
        // chimeric join through are frequently not flagged by it -- a long
        // repeat's unitig coverage is diluted, and a two-copy repeat with one
        // copy in a thin region never reaches the threshold. Identifying the
        // risky joins needs a structural signal rather than a depth one.

        // Retained, default-off negative experiment: require paired support
        // for a lone candidate routed through a repeat. The historical trace
        // counted 381 PROPOSALS, not accepted joins: its print preceded tie,
        // terminal, and mutual-choice checks. The claim that 25 unsupported
        // proposals caused accepted misjoins was therefore unwarranted; this
        // guard's measured output was unchanged. Actual arbitration now emits
        // [acceptedjoin], and terminal-repeat additions emit [prefixtrace].
        //
        // T21 (TESSERACT_FIX_TRUNC_GUARD, default off): a lone candidate is "the sole
        // way through" only when enumerate() discarded nothing. With a discarded rival
        // the side is contested by a route whose destination is unknown; counted, that
        // rival scores 0 (a budget-dropped route is longer than any fragment), so the
        // side needs the paired bar a contested side needs, and the fallbacks cannot
        // weigh a rival they never saw.
        if (cands.size() == 1 && trunc.any() && best < linkBar) {
            ++trLoneBelowBar;
            if (fixTruncGuard) { ++trRefused; return result; }
        }
        if (requireSupportSingle_ && cands.size() == 1 && best <= 0.0) {
            bool throughRepeat = false;
            for (size_t j = 0; j + 1 < cands[0].size(); ++j) {
                if (isRepeat(unitigOf(cands[0][j]))) { throughRepeat = true; break; }
            }
            if (throughRepeat) {
                ++dbgLoneRepeat;
                // (H33 housekeeping: the inline lambda that summed interLen here was
                // indented as if the for loop guarded its return; interLenOf computes the
                // same value.)
                if (joinTrace_)
                    std::fprintf(stderr, "[refused] from=%u dest=%u interLen=%d\n",
                                 unitigOf(tailFirst.back()), unitigOf(cands[0].back()),
                                 interLenOf(cands[0]));
                return result;
            }
        }

        // Output-neutral JOIN_TRACE extension (combo2 PKG-MERGE): the chain and end a candidate
        // terminal would attach to, mirroring the terminal check below (-1 = no chain / mid-chain).
        auto termPort = [&](uint64_t term, long long& chainB, int& endB) {
            chainB = -1;
            endB = -1;
            const uint32_t tu2 = unitigOf(term);
            if (ownerChain[tu2] == UINT32_MAX) return;
            chainB = ownerChain[tu2];
            if (ownerChain[tu2] == c) return;
            const std::vector<uint64_t>& oth = chains[ownerChain[tu2]];
            const uint32_t jx = ownerPos[tu2];
            if (orientOf(oth[jx]) == orientOf(term) && jx == 0) endB = 0;
            else if (orientOf(oth[jx]) != orientOf(term) && jx + 1 == oth.size()) endB = 1;
        };
        // One [fallbackproposal] line per fallback decision, refused ones included: a refused
        // fallback returns before [jointrace] and would otherwise leave no record.
        auto fbTrace = [&](const char* type, int chosen, const char* outcome) {
            if (!joinTrace_) return;
            ++dbgFbTraced;
            if (outcome[0] == 'r') ++dbgFbTracedRefused;   // "refused_*"
            long long chainB = -1;
            int endB = -1;
            if (chosen < 0) {
                std::fprintf(stderr,
                    "[fallbackproposal] chain=%u end=%d from=%u cands=%zu best=%.1f second=%.1f "
                    "linkBar=%.2f type=%s dest=-1 destOrient=-1 destLen=0 destCov=0.0 interLen=-1 "
                    "destChain=-1 destEnd=-1 outcome=%s\n",
                    c, end, unitigOf(tailFirst.back()), cands.size(), best, second, linkBar,
                    type, outcome);
                return;
            }
            const uint64_t term = cands[static_cast<size_t>(chosen)].back();
            termPort(term, chainB, endB);
            std::fprintf(stderr,
                "[fallbackproposal] chain=%u end=%d from=%u cands=%zu best=%.1f second=%.1f "
                "linkBar=%.2f type=%s dest=%u destOrient=%d destLen=%zu destCov=%.1f interLen=%d "
                "destChain=%lld destEnd=%d outcome=%s\n",
                c, end, unitigOf(tailFirst.back()), cands.size(), best, second, linkBar, type,
                unitigOf(term), orientOf(term), g_.nodes[unitigOf(term)].seq.size(),
                g_.nodes[unitigOf(term)].coverage, interLenOf(cands[static_cast<size_t>(chosen)]),
                chainB, endB, outcome);
        };

        // Whether the paired reads decided this, or coverage had to.
        bool byCoverage = false;
        if (cands.size() > 1 && best < linkBar) {
            // No paired evidence to choose with -- which on a library whose
            // mates overlap is the normal case, not the exception. With a
            // ~230 bp fragment against ~190 bp reads only about 0.3% of pairs
            // land on two different unitigs, and once the permissive cutoff
            // connected the graph this became the single largest reason a chain
            // stops: on one panel isolate 616 chain ends had candidates and no
            // support, against 239 successful joins.
            //
            // Two things can still settle it. First, and much the commoner:
            // every candidate may end on the same unitig in the same
            // orientation, differing only in the route taken through the repeat
            // between. Then the destination is not in question at all and the
            // absent paired evidence was never needed -- only the arm is
            // uncertain, and interior coverage picks that. This is the same
            // reasoning the tie branch below already applies; it simply has to
            // apply before support is consulted, not only after.
            //
            // Failing that, coverage continuity: a chain running at 45x
            // continues into sequence at 45x, while a repeat it merely passes
            // through sits at a multiple of that.
            int chosen = -1;
            const char* fbType = "coverage";
            if (matchingAgrees(tailFirst.back(), cands[pick], cands[pick].back(),
                               interLenOf(cands[pick]))) {
                chosen = static_cast<int>(pick);
                fbType = "matching";
                ++dbgMatched;
                if (fallbackMatchingOnly) ++dbgFbmKept;
            } else if (fallbackMatchingOnly) {
                ++dbgFbmRefused;
                ++dbgLowSupport;
                if (joinTrace_) {   // what the withdrawn fallback would have picked (pure helpers)
                    const bool sd = allSameDestination(cands);
                    fbTrace(sd ? "samedest" : "coverage",
                            sd ? pickByInterior(cands) : pickByCoverage(tailFirst, cands),
                            "refused_matchingonly");
                }
                return result;
            } else if (allSameDestination(cands)) {
                chosen = pickByInterior(cands);
                fbType = "samedest";
            } else {
                chosen = pickByCoverage(tailFirst, cands);
            }
            if (chosen < 0) { fbTrace(fbType, chosen, "refused_nopick"); ++dbgLowSupport; return result; }
            // A short destination reached on fallback evidence is where our
            // misassemblies actually come from. Measured over 177 *S. aureus*
            // isolates, splitting every extensive misassembly by the length of
            // the SHORTER alignment block flanking its breakpoint:
            //
            //   shorter flank    model   base   SPAdes   excess     p
            //   < 2 kb             192    185      100      +92   1.9e-04
            //   2 - 20 kb          137    120      140       -3   0.74
            //   >= 20 kb            71     59       40      +31   3.3e-04
            //
            // The short class is identical in `base` and `model` (192 v 185,
            // 5/7, p=0.50), so it is not the organism model, and it is
            // significant even in the 29 isolates where SPAdes is the MORE
            // contiguous assembler (+11, 2/10, p=0.031), so it is not bought
            // contiguity either. 81% of those short blocks sit at a contig
            // terminus, their lengths cluster at IS-element scale, and the
            // sequence they belong to is a median 209 kb away -- a short
            // dispersed repeat appended from the wrong copy.
            //
            // Paired support that clears `linkBar` is not affected: this only
            // withdraws the three fallbacks when what they would append is both
            // short and unvouched. The cost in contiguity should be near zero,
            // because a sub-2 kb tail is not contiguity.
            if (minFallbackDest_ > 0) {
                const uint32_t destU = unitigOf(cands[static_cast<size_t>(chosen)].back());
                if (g_.nodes[destU].seq.size() < minFallbackDest_) {
                    fbTrace(fbType, chosen, "refused_shortdest");
                    ++dbgShortDest;
                    return result;
                }
            }
            // Refuse a fallback whose route no read pair could have spanned.
            //
            // Measured: the <2 kb blocks flanking our excess misassemblies are dispersed
            // repeats with 5-7 copies in the reference, sitting at 6-10x median depth --
            // far ABOVE the 1.6x repeat guard, so they are not escaping it by dilution.
            // They sit in the INTERIOR of the route, between two unique flanks. On
            // GCF003184985v1 the interior is 1,063 bp against an insert window of 204+/-140,
            // so maxPlausible is 764: no pair can reach across it. Every candidate therefore
            // scores zero, `best < linkBar` holds, and the decision falls to the fallbacks --
            // which choose by coverage, with nothing to check themselves against.
            //
            // This is why the three threshold sweeps all came back empty. Raising `linkBar`
            // pushes MORE joins into the fallback path rather than fewer; `tieRatio` never
            // applies there at all; and the first version of this guard tested the
            // DESTINATION unitig's length, which is long unique sequence -- the repeat is
            // interior -- so it withdrew 2 joins out of 166.
            //
            // The test is structural, not depth-based, so the dilution that defeats
            // `isRepeat` does not defeat it: if the sequence between here and the
            // destination is longer than the library's longest plausible fragment, then no
            // paired evidence about this junction exists, and picking anyway is a guess.
            if (noUnspannedFallback_ && insert_.usable) {
                int interLenChosen = 0;
                const std::vector<uint64_t>& route = cands[static_cast<size_t>(chosen)];
                for (size_t j = 0; j + 1 < route.size(); ++j) {
                    interLenChosen += static_cast<int>(g_.nodes[unitigOf(route[j])].seq.size())
                                    - (k_ - 1);
                }
                if (interLenChosen > insert_.maxPlausible) {
                    fbTrace(fbType, chosen, "refused_unspanned");
                    ++dbgUnspanned;
                    return result;
                }
            }
            fbTrace(fbType, chosen, "proposed");
            pick = static_cast<size_t>(chosen);
            byCoverage = true;
            ++dbgCoverage;
        }

        // Candidate diagnostic only: the tie, terminal-end, and mutual-choice
        // checks below can still reject this proposal. Actual accepted joins are
        // recorded separately as [acceptedjoin] after arbitration.
        if (joinTrace_) {
            const uint32_t fromU = unitigOf(tailFirst.back());
            const uint32_t destU = unitigOf(cands[pick].back());
            int interLen = 0;
            for (size_t j = 0; j + 1 < cands[pick].size(); ++j)
                interLen += static_cast<int>(g_.nodes[unitigOf(cands[pick][j])].seq.size()) - (k_ - 1);
            double interDepth = 0; size_t interN = 0;
            for (size_t j = 0; j + 1 < cands[pick].size(); ++j) {
                interDepth += g_.nodes[unitigOf(cands[pick][j])].coverage; ++interN;
            }
            // Trailing fields (combo2 PKG-MERGE): chain, end, pick index and, per candidate,
            // terminal unitig:orient:score:interLen:destChain:destEnd.
            std::string candScores;
            for (size_t i = 0; i < cands.size(); ++i) {
                long long chainB = -1;
                int endB = -1;
                termPort(cands[i].back(), chainB, endB);
                char buf[160];
                std::snprintf(buf, sizeof(buf), "%s%u:%d:%.1f:%d:%lld:%d", i ? "," : "",
                              unitigOf(cands[i].back()), orientOf(cands[i].back()), scores[i],
                              interLenOf(cands[i]), chainB, endB);
                candScores += buf;
            }
            ++dbgJtExt;
            std::fprintf(stderr,
                "[jointrace] from=%u fromLen=%zu fromCov=%.1f dest=%u destLen=%zu destCov=%.1f "
                "cands=%zu best=%.1f second=%.1f interLen=%d interN=%zu interCov=%.1f "
                "byCov=%d med=%.1f maxPl=%d chain=%u end=%d pick=%zu candScores=%s\n",
                fromU, g_.nodes[fromU].seq.size(), g_.nodes[fromU].coverage,
                destU, g_.nodes[destU].seq.size(), g_.nodes[destU].coverage,
                cands.size(), best, second, interLen, interN,
                interN ? interDepth / interN : 0.0, byCoverage ? 1 : 0,
                medianCoverage_, insert_.usable ? insert_.maxPlausible : -1,
                c, end, pick, candScores.c_str());
        }

        if (cands.size() > 1 && !byCoverage) {
            bool tieExempt = false;   // L4 weak or XO1 nomination: skip the same-destination re-pick
            // A near tie means the repeat is genuinely unresolved; guessing
            // would manufacture a misassembly.
            if (second > 0 && best < tieRatio_ * second) {
                // Unless every near-tied path ends on the same unitig in the
                // same orientation. Then the destination is not in doubt at
                // all. The broad-window pair COUNTS can tie even when insert
                // lengths distinguish the interiors. Preserve the established
                // endpoint join; legacy behavior settles the route by coverage.
                // The default-off distance model below can use that otherwise
                // discarded information without changing endpoint support.
                const uint64_t bestTerm = cands[pick].back();
                bool sameDestination = true;
                for (size_t i = 0; i < cands.size(); ++i) {
                    if (scores[i] * tieRatio_ < best) continue;   // not near-tied
                    if (cands[i].back() != bestTerm) { sameDestination = false; break; }
                }
                if (!sameDestination) {
                    ++dbgTie;
                    dbgTieBest += static_cast<long long>(best);
                    dbgTieSecond += static_cast<long long>(second);
                    if (!oneSidedTie && !tieSingleExit) return result;
                    const uint64_t tailOid = tailFirst.back();
                    size_t liveExits = 0;
                    for (const Link& l : g_.exits(unitigOf(tailOid), orientOf(tailOid)))
                        if (!g_.nodes[l.to].deleted) ++liveExits;
                    const size_t lenA = acc + ov;
                    size_t lenB = 0;
                    const uint32_t cbw = ownerChain[unitigOf(cands[pick].back())];
                    if (cbw != UINT32_MAX) {
                        for (uint64_t o : chains[cbw]) lenB += addedLen(o);
                        lenB += ov;
                    }
                    const size_t shorter = std::min(lenA, lenB);
                    if (tieSingleExit && liveExits == 1 && shorter < tieSingleExitLong) {
                        tieExempt = true;
                        ++dbgSingleExit;
                    } else {
                    if (tieSingleExit && liveExits == 1) ++dbgSingleExitLong;
                    if (!oneSidedTie) return result;
                    // Weak nomination.  Refused when the shorter chain is >= the gate, except
                    // (L4 alone, XO1 off) at a single-exit end, which L4 never gated.
                    const bool gateExempt = liveExits < 2 && !tieSingleExit;
                    if (shorter >= oneSidedLongBp && !gateExempt) {
                        ++dbgWeakLongBranch;
                        return result;
                    }
                    result.weak = true;
                    tieExempt = true;
                    ++dbgWeakNom;
                    }
                }
                if (!tieExempt) {
                double bestInterior = -1.0;
                for (size_t i = 0; i < cands.size(); ++i) {
                    if (scores[i] * tieRatio_ < best) continue;
                    double minCov = std::numeric_limits<double>::max();
                    for (size_t j = 0; j + 1 < cands[i].size(); ++j) {
                        minCov = std::min(minCov, g_.nodes[unitigOf(cands[i][j])].coverage);
                    }
                    // A path with no interior is the most direct route.
                    if (cands[i].size() < 2) minCov = std::numeric_limits<double>::max();
                    if (minCov > bestInterior) { bestInterior = minCov; pick = i; }
                }
                if (routeDistance_ && routeDensity_.usable) {
                    // Endpoint nomination and its count score are settled above.
                    // Only distinguish near-tied routes to that SAME endpoint.
                    // Same-length routes remain one hypothesis, so duplicating
                    // an enumerated path cannot manufacture paired evidence.
                    std::vector<int> lengths;
                    for (size_t i = 0; i < cands.size(); ++i) {
                        if (scores[i] * tieRatio_ >= best) lengths.push_back(interLenOf(cands[i]));
                    }
                    std::sort(lengths.begin(), lengths.end());
                    lengths.erase(std::unique(lengths.begin(), lengths.end()), lengths.end());
                    if (lengths.size() > 1) {
                        std::vector<int> spans;
                        std::vector<uint64_t> seenSources;
                        for (size_t m = tailFirst.size(); m-- > 0;) {
                            if (double(distToEnd[m]) > reach) break;
                            if (isRepeat(unitigOf(tailFirst[m])) ||
                                std::find(seenSources.begin(), seenSources.end(), tailFirst[m]) != seenSources.end())
                                continue;
                            seenSources.push_back(tailFirst[m]);
                            const auto source = support_.find(tailFirst[m]);
                            if (source == support_.end()) continue;
                            const auto target = source->second.find(bestTerm);
                            if (target == source->second.end()) continue;
                            for (int span : target->second) {
                                const int64_t adjusted = int64_t(span) + int64_t(distToEnd[m]);
                                if (adjusted >= std::numeric_limits<int>::min() &&
                                    adjusted <= std::numeric_limits<int>::max()) spans.push_back(int(adjusted));
                            }
                        }
                        std::vector<double> contrast;
                        // T31 (TESSERACT_FIX_ROUTE_ORDER, default off): support_ lists are
                        // merged in thread order, so the SAME multiset of spans arrives in a
                        // -t dependent order and the floating-point totals below differ in the
                        // last bits with the thread count. A canonical order makes them exact.
                        ++roDecisions;
                        if (fixRouteOrder) { std::sort(spans.begin(), spans.end()); ++roSorted; }
                        const auto allocated = detail::routeDistanceAllocation(routeDensity_, spans, lengths, &contrast, linkBar);
                        if (joinTrace_) {
                            // Exact bits (%a): the [routedistance] line's %.6f hides -t differences.
                            std::fprintf(stderr, "[routeorder-trace] chain=%u end=%d sorted=%d pairs=%zu totals=",
                                         c, end, fixRouteOrder ? 1 : 0, spans.size());
                            for (size_t i = 0; i < allocated.size(); ++i)
                                std::fprintf(stderr, "%s%a", i ? "," : "", allocated[i]);
                            std::fprintf(stderr, " contrasts=");
                            for (size_t i = 0; i < contrast.size(); ++i)
                                std::fprintf(stderr, "%s%a", i ? "," : "", contrast[i]);
                            std::fprintf(stderr, "\n");
                        }
                        // The existing bar applies to distinguishing route mass,
                        // not common/flat pairs. Fractions are never rounded up:
                        // two near-certain pairs may need a third to clear bar2.
                        const int winner = detail::supportedRouteLength(allocated, contrast, tieRatio_, linkBar);
                        if (winner >= 0) {
                            const size_t legacy = pick;
                            double interior = -1;
                            for (size_t i = 0; i < cands.size(); ++i) {
                                if (scores[i] * tieRatio_ < best || interLenOf(cands[i]) != lengths[size_t(winner)]) continue;
                                double minCov = std::numeric_limits<double>::max();
                                for (size_t j = 0; j + 1 < cands[i].size(); ++j)
                                    minCov = std::min(minCov, g_.nodes[unitigOf(cands[i][j])].coverage);
                                if (minCov > interior) { interior = minCov; pick = i; }
                            }
                            if (pick != legacy) {
                                result.legacyConnector.assign(cands[legacy].begin(), cands[legacy].end() - 1);
                                result.distanceChanged = true;
                            }
                            if (joinTrace_) {
                                std::fprintf(stderr,
                                    "[routedistance] chain=%u end=%d target=%llu pairs=%zu oldLen=%d newLen=%d changed=%d allocations=",
                                    c, end, static_cast<unsigned long long>(bestTerm), spans.size(),
                                    interLenOf(cands[legacy]), interLenOf(cands[pick]), legacy != pick ? 1 : 0);
                                for (size_t i = 0; i < lengths.size(); ++i)
                                    std::fprintf(stderr, "%s%d:%.6f", i ? "," : "", lengths[i], allocated[i]);
                                std::fprintf(stderr, " contrasts=");
                                for (size_t i = 0; i < lengths.size(); ++i)
                                    std::fprintf(stderr, "%s%d:%.6f", i ? "," : "", lengths[i], contrast[i]);
                                std::fprintf(stderr, "\n");
                            }
                        }
                    }
                }
                }   // !tieExempt
            }
        }

        // T07 (TESSERACT_FIX_REVISIT_GUARD, default off). enumerate() never extends a frame to
        // a node it already holds, so every route through a tandem R -> body -> R carries ONE
        // copy of the unit, whatever the genome has. When the pick passes a node that its own
        // history could have looped back to, its copy number is 1 by construction, not by
        // evidence: QUAST finds whole tandem units missing (21 K2 local misassemblies at +1/2/3x
        // the loop length) and shorter ones as indels. The guard refuses such a pick; the end
        // stays open. The exact read-thread fallback can still nominate the join, because a
        // read thread may repeat a node and so carries the real copy number. (Deciding the
        // copy number from paired distances was prototyped and dropped: the distance model
        // was 24-129 bp off on real libraries -- see combo3/fix_G-resolve/README.md.)
        if (!loops.prefix.empty()) {
            const std::vector<uint64_t>& route = cands[pick];
            bool onRoute = false;   // a loop entered with the pick's own history
            for (const auto& pre : loops.prefix)
                if (pre.size() < route.size() && std::equal(pre.begin(), pre.end(), route.begin())) {
                    onRoute = true;
                    break;
                }
            if (onRoute) {
                ++rvPicked;
                if (fixRevisitGuard) { ++rvRefused; return result; }
            }
        }

        // The terminal has to be an *end* of its chain; landing in the middle
        // would mean cutting an already-supported contig in half.
        const uint64_t terminal = cands[pick].back();
        const uint32_t tu = unitigOf(terminal);
        const uint32_t cb = ownerChain[tu];
        if (cb == UINT32_MAX || cb == c) return result;
        const std::vector<uint64_t>& other = chains[cb];
        const uint32_t j = ownerPos[tu];
        const int storedOrient = orientOf(other[j]);

        if (storedOrient == orientOf(terminal) && j == 0) result.endB = 0;
        else if (storedOrient != orientOf(terminal) && j + 1 == other.size()) result.endB = 1;
        else { ++dbgMidChain; return result; }
        ++dbgOk;

        result.chainB = cb;
        result.connector.assign(cands[pick].begin(), cands[pick].end() - 1);
        result.score = best;
        result.decisive = (cands.size() == 1 && best >= linkBar) ||
                          (cands.size() > 1 && !byCoverage && second <= 0);
        result.ok = true;
        return result;
    };

    // A separate fallback preserves every accepted legacy nomination. This
    // first experiment uses routes beginning at the actual chain port, without
    // inferring upstream history or skipping other anchorable chain members.
    auto exactThreadContinuation = [&](uint32_t c, int end) {
        Cont result;
        if (chains[c].empty()) return result;
        const uint64_t source = end == 1 ? chains[c].back() : flip(chains[c].front());
        const auto found = threadPorts.find(source);
        if (found == threadPorts.end()) return result;
        // No count dominance rule: even a single contrary observed full route
        // means the short reads do not give a mutually consistent continuation.
        if (found->second.size() != 1) { ++threadConflictingPorts; return result; }
        const auto entry = found->second.front();
        if (static_cast<double>(freshThreadSupport[entry.first]) < linkBar) return result;
        auto route = threadEvidence.routes[entry.first].oriented;
        if (entry.second) route = reverseReadThreadPath(route);
        for (size_t j = 1; j + 1 < route.size(); ++j)
            if (anchorable(unitigOf(route[j]))) return result;
        const uint64_t target = route.back();
        const uint32_t targetNode = unitigOf(target);
        const uint32_t otherChain = ownerChain[targetNode];
        if (otherChain == UINT32_MAX || otherChain == c) return result;
        const auto& other = chains[otherChain];
        const size_t pos = ownerPos[targetNode];
        if (pos >= other.size()) return result;
        if (other[pos] == target && pos == 0) result.endB = 0;
        else if (other[pos] == flip(target) && pos + 1 == other.size()) result.endB = 1;
        else return result;
        result.chainB = otherChain;
        result.connector.assign(route.begin() + 1, route.end() - 1);
        result.score = static_cast<double>(freshThreadSupport[entry.first]);
        result.exactThread = true;
        result.ok = true;
        ++threadNominations;
        return result;
    };

    for (int round = 0; round < 24; ++round) {
        reindex();
        const size_t nc = chains.size();
        std::vector<Cont> cont(nc * 2);
        for (uint32_t c = 0; c < nc; ++c) {
            if (chains[c].empty()) continue;
            cont[c * 2 + 0] = bestContinuation(c, 0);
            cont[c * 2 + 1] = bestContinuation(c, 1);
        }

        if (routeDistance_) {
            // Both nominations retain their original endpoint. Unequal chain
            // histories can prefer different interiors; in that case restore
            // their legacy routes before arbitration rather than making the
            // chosen route depend on which chain happens to be visited first.
            // This does not reject a join or alter any endpoint/support score.
            for (uint32_t c = 0; c < nc; ++c) for (int e = 0; e < 2; ++e) {
                Cont& f = cont[c * 2 + e];
                if (!f.ok) continue;
                Cont& back = cont[f.chainB * 2 + f.endB];
                if (!back.ok || back.chainB != c || back.endB != e ||
                    (!f.distanceChanged && !back.distanceChanged)) continue;
                bool same = f.connector.size() == back.connector.size();
                for (size_t j = 0; same && j < f.connector.size(); ++j)
                    same = f.connector[j] == flip(back.connector[back.connector.size() - 1 - j]);
                if (same) continue;
                if (f.distanceChanged) f.connector = f.legacyConnector;
                if (back.distanceChanged) back.connector = back.legacyConnector;
                f.distanceChanged = back.distanceChanged = false;
                if (joinTrace_) std::fprintf(stderr,
                    "[routedistance-fallback] round=%d chainA=%u endA=%d chainB=%u endB=%d reason=mirror-route-disagreement\n",
                    round, c, e, f.chainB, f.endB);
            }
        }

        if (exactReadThreads) {
            for (uint32_t c = 0; c < nc; ++c) for (int e = 0; e < 2; ++e) {
                if (!cont[c * 2 + e].ok) cont[c * 2 + e] = exactThreadContinuation(c, e);
                else if (cont[c * 2 + e].weak) {
                    // Exact read threads outrank a weak (near-tied) nomination.
                    Cont t = exactThreadContinuation(c, e);
                    if (t.ok) { cont[c * 2 + e] = t; ++dbgWeakThread; }
                }
            }
        }

        std::vector<char> merged(nc, 0);
        size_t joins = 0;
        // Existing reciprocal choices have first claim on chains each round.
        // New fallback joins cannot consume a chain before such a legacy join.
        for (int pass = 0; pass < (exactReadThreads ? 2 : 1); ++pass) {
            for (uint32_t c = 0; c < nc; ++c) {
                if (merged[c] || chains[c].empty()) continue;
                for (int e = 0; e < 2; ++e) {
                    const Cont& f = cont[c * 2 + e];
                    if (!f.ok || merged[f.chainB] || merged[c]) continue;
                    // Only join when both ends independently chose each other.
                    const Cont& back = cont[f.chainB * 2 + f.endB];
                    if (!back.ok || back.chainB != c || back.endB != e) continue;
                    if (f.weak || back.weak) {
                        // Only reachable with TESSERACT_ONE_SIDED_TIE.
                        if (f.weak && back.weak) { ++dbgWeakTwoSided; continue; }
                        const Cont& strong = f.weak ? back : f;
                        if (!strong.decisive || strong.exactThread) { ++dbgWeakPartner; continue; }
                    }
                    const bool threadJoin = f.exactThread || back.exactThread;
                    if (exactReadThreads && threadJoin != (pass == 1)) continue;
                    if (threadJoin) {
                        bool same = f.connector.size() == back.connector.size();
                        for (size_t j = 0; same && j < f.connector.size(); ++j)
                            same = f.connector[j] == flip(back.connector[back.connector.size() - 1 - j]);
                        if (!same) continue;
                    }

                    // R6 (TESSERACT_FIX_MIRROR_ROUTE, default off). A legacy join takes the
                    // interior route of whichever end is visited first (the lower chain index),
                    // so when the two ends chose different interiors to the same endpoints the
                    // emitted arm depends on chain numbering. Choose symmetrically instead: the
                    // route whose weakest interior unitig is better covered (the rule the tie
                    // branch uses), then the smaller orientation-independent SEQUENCE (node ids
                    // are numbering, so they cannot break the tie).
                    std::vector<uint64_t> connector = f.connector;
                    if (!threadJoin) {
                        std::vector<uint64_t> other;
                        for (size_t j = back.connector.size(); j-- > 0;) other.push_back(flip(back.connector[j]));
                        if (other != f.connector) {
                            ++mrDisagree;
                            if (fixMirrorRoute) {
                                auto minCov = [&](const std::vector<uint64_t>& r) {
                                    double mc = std::numeric_limits<double>::max();
                                    for (uint64_t o : r) mc = std::min(mc, g_.nodes[unitigOf(o)].coverage);
                                    return mc;
                                };
                                auto canon = [&](const std::vector<uint64_t>& r) {
                                    std::string sp;
                                    for (uint64_t o : r) {
                                        const std::string piece = g_.oriented(unitigOf(o), orientOf(o));
                                        sp += sp.empty() ? piece : piece.substr(ov);
                                    }
                                    return std::min(sp, reverseComplement(sp));
                                };
                                const double cf = minCov(f.connector), cb = minCov(other);
                                const bool takeBack = cb > cf || (cb == cf && canon(other) < canon(f.connector));
                                if (takeBack) { connector = other; ++mrChoseBack; }
                            }
                        }
                    }

                    std::vector<uint64_t> a;
                    if (e == 1) a = chains[c];
                    else {
                        a.reserve(chains[c].size());
                        for (size_t i = chains[c].size(); i-- > 0;) a.push_back(flip(chains[c][i]));
                    }
                    std::vector<uint64_t> b;
                    if (f.endB == 0) b = chains[f.chainB];
                    else {
                        b.reserve(chains[f.chainB].size());
                        for (size_t i = chains[f.chainB].size(); i-- > 0;) {
                            b.push_back(flip(chains[f.chainB][i]));
                        }
                    }

                    if (joinTrace_) {
                        bool sameRoute = f.connector.size() == back.connector.size();
                        if (sameRoute) {
                            for (size_t j = 0; j < f.connector.size(); ++j) {
                                if (f.connector[j] != flip(back.connector[back.connector.size() - 1 - j])) {
                                    sameRoute = false;
                                    break;
                                }
                            }
                        }
                        std::fprintf(stderr,
                            "[acceptedjoin] round=%d chainA=%u endA=%d chainB=%u endB=%d "
                            "score=%.1f backScore=%.1f sameRoute=%d route=%llu",
                            round, c, e, f.chainB, f.endB, f.score, back.score,
                            sameRoute ? 1 : 0, static_cast<unsigned long long>(a.back()));
                        for (uint64_t oid : f.connector)
                            std::fprintf(stderr, ",%llu", static_cast<unsigned long long>(oid));
                        std::fprintf(stderr, ",%llu\n", static_cast<unsigned long long>(b.front()));
                    }

                    a.insert(a.end(), connector.begin(), connector.end());
                    a.insert(a.end(), b.begin(), b.end());
                    chains[c].swap(a);
                    chains[f.chainB].clear();
                    merged[c] = 1;
                    merged[f.chainB] = 1;
                    ++joins;
                    if (f.weak || back.weak) {
                        ++dbgWeakJoin;
                        if (joinTrace_) std::fprintf(stderr,
                            "[weakjoin] round=%d chainA=%u endA=%d chainB=%u endB=%d\n",
                            round, c, e, f.chainB, f.endB);
                    }
                    if (threadJoin) {
                        ++threadJoins;
                        if (joinTrace_) std::fprintf(stderr,
                            "[readthread-join] round=%d chainA=%u endA=%d chainB=%u endB=%d "
                            "molecules=%.0f backMolecules=%.0f\n",
                            round, c, e, f.chainB, f.endB, f.score, back.score);
                    }
                    break;
                }
            }
        }
        if (joins == 0) break;
    }
    if (exactReadThreads) std::fprintf(stderr,
        "[readthread-result] nominations=%zu conflictingPorts=%zu joins=%zu\n",
        threadNominations, threadConflictingPorts, threadJoins);
    // build_v3: the package counters are printed on every run, zeros included, with a trailing
    // enabled= field (OBJECTIVE A2); the leading fields are the package's own.
    std::fprintf(stderr,
        "[onesidedtie] weakNominations=%lld weakJoins=%lld refusedLongBranching=%lld "
        "refusedTwoSided=%lld refusedPartnerNotDecisive=%lld threadPreferred=%lld long=%zu enabled=%d\n",
        dbgWeakNom, dbgWeakJoin, dbgWeakLongBranch, dbgWeakTwoSided, dbgWeakPartner, dbgWeakThread,
        oneSidedLongBp, oneSidedTie ? 1 : 0);
    std::fprintf(stderr, "[tiesingleexit] exempted=%lld refusedLong=%lld long=%zu enabled=%d\n",
                 dbgSingleExit, dbgSingleExitLong,
                 tieSingleExitLong == std::numeric_limits<size_t>::max() ? size_t(0) : tieSingleExitLong,
                 tieSingleExit ? 1 : 0);
    std::fprintf(stderr, "[fallbackmatching] kept=%lld refused=%lld enabled=%d\n", dbgFbmKept, dbgFbmRefused,
                 fallbackMatchingOnly ? 1 : 0);
    if (joinTrace_) std::fprintf(stderr,
        "[jointrace-ext] jointraceWithCands=%lld fallbackProposals=%lld fallbackRefused=%lld\n",
        dbgJtExt, dbgFbTraced, dbgFbTracedRefused);

    if (env::present("TESSERACT_DEBUG_RESOLVE")) {
        std::fprintf(stderr,
                     "      [debug] continuation outcomes: ok=%lld no-candidate=%lld "
                     "low-support=%lld tie=%lld mid-chain=%lld by-coverage=%lld matched=%lld "
                     "short-dest=%lld unspanned=%lld lone-repeat=%lld  (tie mean best=%.1f second=%.1f)\n",
                     dbgOk, dbgNoCand, dbgLowSupport, dbgTie, dbgMidChain, dbgCoverage, dbgMatched,
                     dbgShortDest, dbgUnspanned, dbgLoneRepeat,
                     dbgTie ? static_cast<double>(dbgTieBest) / static_cast<double>(dbgTie) : 0.0,
                     dbgTie ? static_cast<double>(dbgTieSecond) / static_cast<double>(dbgTie) : 0.0);
    }

    // ---- scaffolding -----------------------------------------------------
    // Chain ends with strong paired support but no path through the graph are
    // joined across a gap of Ns whose length comes from the fragment model.
    // Each chain has two ports; a join consumes one port at each end, so the
    // accepted joins form simple paths that are then walked in order.
    const size_t nc = chains.size();
    std::vector<uint32_t> joinTo(nc * 2, UINT32_MAX);   // port -> port
    std::vector<int> joinGap(nc * 2, 0);
    // Per accepted join, both ports: the unfloored estimate behind joinGap (release frame:
    // true gap + k-1), its standard error, whether the T20 estimate says the flanks
    // overlap or abut (true gap <= 0), and its support (T13's cycle break).
    std::vector<int> joinEst(nc * 2, 0);
    std::vector<double> joinSigma(nc * 2, 0.0);
    std::vector<char> joinOverlapEst(nc * 2, 0);
    std::vector<double> joinScore(nc * 2, 0.0);
    // T20 counters (always printed).
    size_t geJoins = 0, geEstimated = 0, geChanged = 0, geFloor1 = 0, geOverlapEst = 0, geSkippedNoMean = 0;
    std::vector<int> geShift;

    // T20: forced insert bounds (TESSERACT_QC_INSERT) without a fitted model leave
    // insert_.mean at 0, so every gap would be -span-dist: skip scaffolding instead.
    const bool noFittedMean = haveForcedBounds_ && insert_.observations == 0;
    if (scaffolding_ && insert_.usable && fixGapEstimate && noFittedMean) geSkippedNoMean = 1;

    if (scaffolding_ && insert_.usable && !(fixGapEstimate && noFittedMean)) {
        reindex();
        struct End { uint64_t oriented; std::vector<uint64_t> members; std::vector<size_t> dist; };
        std::vector<End> ends(nc * 2);
        for (uint32_t c = 0; c < nc; ++c) {
            if (chains[c].empty()) continue;
            for (int e = 0; e < 2; ++e) {
                std::vector<uint64_t> tailFirst;
                if (e == 1) tailFirst.assign(chains[c].begin(), chains[c].end());
                else {
                    for (size_t i = chains[c].size(); i-- > 0;) tailFirst.push_back(flip(chains[c][i]));
                }
                End& E = ends[c * 2 + e];
                E.oriented = tailFirst.back();
                size_t acc = 0;
                for (size_t i = tailFirst.size(); i-- > 0;) {
                    if (static_cast<double>(acc) > reach) break;
                    E.members.push_back(tailFirst[i]);
                    E.dist.push_back(acc);
                    acc += addedLen(tailFirst[i]);
                }
            }
        }

        // The port an oriented unitig represents when a fragment enters it.
        std::unordered_map<uint64_t, uint32_t> entryOf;
        for (uint32_t c = 0; c < nc; ++c) {
            if (chains[c].empty()) continue;
            entryOf[chains[c].front()] = c * 2 + 0;
            entryOf[flip(chains[c].back())] = c * 2 + 1;
        }

        struct Best { uint32_t partner = UINT32_MAX; double score = 0; int gap = 0; };
        std::vector<Best> best(ends.size());
        for (size_t idx = 0; idx < ends.size(); ++idx) {
            const End& E = ends[idx];
            if (E.members.empty()) continue;
            // Anything reachable through the graph was already handled by chain
            // extension; scaffolding only spans true gaps.
            std::vector<std::vector<uint64_t>> reachable;
            enumerate(E.oriented, reachable);
            std::unordered_map<uint64_t, char> viaGraph;
            for (const auto& r : reachable) viaGraph[r.back()] = 1;

            std::unordered_map<uint32_t, std::vector<int>> gaps;
            for (size_t m = 0; m < E.members.size(); ++m) {
                auto it = support_.find(E.members[m]);
                if (it == support_.end()) continue;
                for (const auto& kv : it->second) {
                    if (viaGraph.count(kv.first)) continue;
                    auto pt = entryOf.find(kv.first);
                    if (pt == entryOf.end()) continue;
                    if (pt->second / 2 == idx / 2) continue;   // same chain
                    for (int32_t span : kv.second) {
                        const int gap = static_cast<int>(insert_.mean) - span -
                                        static_cast<int>(E.dist[m]);
                        if (gap < -static_cast<int>(ov) || gap > insert_.maxPlausible) continue;
                        gaps[pt->second].push_back(gap);
                    }
                }
            }
            // `gaps` is an unordered_map fed from support_, whose insertion
            // order depends on which thread saw a read pair first. A strict >
            // therefore let two partners tied on supporting-pair count be
            // settled by iteration order, changing both the partner and the
            // gap with the thread count. The port breaks the tie instead.
            for (auto& kv : gaps) {
                const double score = static_cast<double>(kv.second.size());
                if (score > best[idx].score ||
                    (score == best[idx].score && kv.first < best[idx].partner)) {
                    std::sort(kv.second.begin(), kv.second.end());
                    best[idx].score = score;
                    best[idx].partner = kv.first;
                    best[idx].gap = kv.second[kv.second.size() / 2];
                }
            }
        }

        // How many spanning pairs a join needs. A fixed count cannot be right
        // at both ends of the depth range: the number of pairs crossing a
        // junction scales with coverage, so a flat 5 is nearly unreachable on a
        // 40x library and trivially met on a 200x one. On one panel isolate 153
        // of the 173 dead ends that found a partner were rejected on this
        // threshold alone -- the mutual-best rule turned away none of them.
        // Scaling with the observed depth keeps the evidence bar constant in
        // the units that matter, with a floor of 3 so a shallow library still
        // needs corroboration.
        const double scafOverride = env::real("TESSERACT_SCAF_SUPPORT", 0.0);
        const double minScaffoldSupport =
            scafOverride > 0 ? scafOverride
            : minScaffoldSupport_ > 0
                ? static_cast<double>(minScaffoldSupport_)
                : std::min(10.0, std::max(3.0, medianCoverage_ * 0.06));
        for (uint32_t idx = 0; idx < best.size(); ++idx) {
            const Best& b = best[idx];
            if (b.partner == UINT32_MAX || b.score < minScaffoldSupport) continue;
            if (best[b.partner].partner != idx) continue;      // must be mutual
            if (joinTo[idx] != UINT32_MAX || joinTo[b.partner] != UINT32_MAX) continue;
            joinTo[idx] = b.partner;
            joinTo[b.partner] = idx;
            joinGap[idx] = std::max(1, b.gap);
            joinGap[b.partner] = joinGap[idx];
            joinEst[idx] = joinEst[b.partner] = b.gap;
            joinSigma[idx] = joinSigma[b.partner] =
                1.2533 * insert_.stddev / std::sqrt(std::max(1.0, b.score));
            joinScore[idx] = joinScore[b.partner] = std::min(b.score, best[b.partner].score);
            ++geJoins;
            if (b.gap < 1) ++geFloor1;
        }

        // T20 (TESSERACT_FIX_GAP_ESTIMATE, default off). The release gap of a join is the
        // estimate of whichever port has the lower index -- the median over ITS members'
        // pairs to the partner's entry unitig -- so renumbering the nodes changes the N-run
        // (47/48 fixture replicates). Here both ports' estimates are averaged, which no
        // numbering can change; each is still the release estimate, so the release's
        // accuracy is kept. (The release estimate is also biased short where the insert
        // distribution is wide -- linking pairs are length-biased -- but that bias did not
        // follow any model on real libraries: MEASURED against the references, an analytic
        // spanning-conditioned correction moved efaecalis GCA053117665v1 from -220 to +73 bp
        // but abaumannii GCF046097175v1 from +3 to +336 bp, and a correction calibrated on
        // graph junctions does not transfer to scaffold gaps. Not corrected; see README.)
        // With T03, an estimate <= 0 without a verified flank overlap is written as an
        // unknown-length gap instead of a 1-N adjacency claim.
        if (fixGapEstimate) {
            for (uint32_t idx = 0; idx < joinTo.size(); ++idx) {
                const uint32_t q = joinTo[idx];
                if (q == UINT32_MAX || q < idx) continue;
                ++geEstimated;
                const double sum = static_cast<double>(best[idx].gap) + static_cast<double>(best[q].gap);
                const int est = static_cast<int>(std::floor(sum / 2.0 + 0.5));   // release frame
                if (est != joinEst[idx]) ++geChanged;
                geShift.push_back(est - joinEst[idx]);
                joinEst[idx] = joinEst[q] = est;
                joinGap[idx] = joinGap[q] = std::max(1, est);
                joinSigma[idx] = joinSigma[q] =
                    1.2533 * insert_.stddev / std::sqrt(std::max(1.0, std::max(best[idx].score, best[q].score)));
                const bool overlapEst = est - static_cast<int>(ov) <= 0;
                joinOverlapEst[idx] = joinOverlapEst[q] = overlapEst ? 1 : 0;
                if (overlapEst) ++geOverlapEst;
            }
        }
    }

    // ---- scaffold cycles (T13) ---------------------------------------------
    // Accepted joins give each port at most one partner and never pair a chain with
    // itself, so they form simple paths AND simple cycles (a circular element split at two
    // or more scaffold gaps). The emit walk below starts only at a free port, so a cycle
    // has no start, and the fallback after it renders the cycle's chains one by one: every
    // join in the cycle was accepted and is then silently dropped (no N-run, no
    // scaffoldJoins, no gap for the closer). Counted on every run (stderr only);
    // TESSERACT_FIX_SCAFFOLD_CYCLE=1 breaks each cycle once, at its weakest join (lowest
    // support; tie -> lowest port), so the rest is walked as one linear scaffold.
    size_t scCycles = 0, scChains = 0, scJoins = 0, scBases = 0, scBroken = 0;
    {
        std::vector<char> seenChain(nc, 0);
        for (uint32_t c = 0; c < nc; ++c) {
            if (seenChain[c] || chains[c].empty()) continue;
            if (joinTo[c * 2 + 0] == UINT32_MAX || joinTo[c * 2 + 1] == UINT32_MAX) continue;
            // Walk out of port 1; a cycle comes back to c through port 0. A walk that meets
            // a chain already seen, or a free port, is on a path (a cycle's chains are
            // reachable only from inside it, and the first visit walks all of it).
            std::vector<uint32_t> outPorts, members;
            uint32_t cur = c;
            int inPort = 0;
            bool isCycle = false;
            for (size_t steps = 0; steps <= nc; ++steps) {
                members.push_back(cur);
                seenChain[cur] = 1;
                const uint32_t out = cur * 2 + (inPort == 0 ? 1 : 0);
                const uint32_t partner = joinTo[out];
                if (partner == UINT32_MAX) break;
                outPorts.push_back(out);
                cur = partner / 2;
                inPort = static_cast<int>(partner % 2);
                if (cur == c) { isCycle = inPort == 0; break; }
                if (seenChain[cur]) break;
            }
            if (!isCycle) continue;
            ++scCycles;
            scChains += members.size();
            scJoins += outPorts.size();
            for (uint32_t m : members)
                for (size_t i = 0; i < chains[m].size(); ++i)
                    scBases += g_.nodes[unitigOf(chains[m][i])].seq.size() - (i ? ov : 0);
            if (joinTrace_) {
                std::fprintf(stderr, "[scafcycle-detail] cycle=%zu chains=", scCycles);
                for (size_t i = 0; i < members.size(); ++i)
                    std::fprintf(stderr, "%s%u:u%u", i ? "," : "", members[i], unitigOf(chains[members[i]].front()));
                std::fprintf(stderr, " joins=");
                for (size_t i = 0; i < outPorts.size(); ++i)
                    std::fprintf(stderr, "%s%u-%u:s%.0f:g%d", i ? "," : "", outPorts[i], joinTo[outPorts[i]],
                                 joinScore[outPorts[i]], joinGap[outPorts[i]]);
                std::fprintf(stderr, "\n");
            }
            if (!fixCycleBreak) continue;
            uint32_t weakest = outPorts.front();
            for (uint32_t p : outPorts) {
                const uint32_t lo = std::min(p, joinTo[p]);
                const uint32_t wlo = std::min(weakest, joinTo[weakest]);
                if (joinScore[p] < joinScore[weakest] || (joinScore[p] == joinScore[weakest] && lo < wlo))
                    weakest = p;
            }
            const uint32_t other = joinTo[weakest];
            joinTo[weakest] = joinTo[other] = UINT32_MAX;
            joinGap[weakest] = joinGap[other] = 0;
            ++scBroken;
        }
    }

    // ---- emit ------------------------------------------------------------
    std::vector<char> placed(n, 0);
    std::vector<char> done(nc, 0);

    // T03 (TESSERACT_FIX_GAP_FLANK, default off). The first unitig after a scaffold N-gap
    // is not graph-adjacent to what precedes it -- the Ns are -- yet the release renders
    // it like a graph successor and drops its first k-1 bases, while joinGap estimates the
    // true gap PLUS k-1 (span = consumedU + consumedV - (k-1)), so the Ns stand in for
    // those real bases and the GFA P-line (gN + the WHOLE segment) spells k-1 more than
    // the FASTA. A gap the gap filler closes gets those bases back (it spells everything
    // between its seed and a target k-mer that lies wholly in full-depth sequence); a gap
    // it leaves open loses them. Rendering the unitig whole here instead moved the
    // filler's target onto the dropout tip and cost closures (MEASURED, ecloacae
    // GCF053691565v1: 15 -> 10 of 32 closed with G-emit's back-off on). So the resolver
    // keeps the release layout, which is what the filler sees, and records per gap what an
    // OPEN gap must look like; restoreGapFlanks() applies it after gap filling: the
    // unitig's first k-1 bases back, and the N-run = joinGap - (k-1), the estimated true
    // gap. Two flanks may genuinely overlap by fewer than k-1 bases (a dropout of a few
    // k-mers); an exact suffix/prefix overlap of L >= kMinFlankOverlap bases that the
    // estimate does not contradict is restored once and the gap is one N. An estimate
    // <= 0 without such an overlap is one N, or, with T20, an unknown-length gap of
    // kUnknownGapN Ns. curPath.gaps carries the open-gap count, so for an open gap FASTA,
    // AGP and the GFA P-line agree (an overlap-case trim of L bases cannot be expressed in
    // a GFA1 P-line; counted as pLineOverspell).
    int flankJoin = 0;            // > 0: the next piece follows the gap of port flankPort
    uint32_t flankPort = 0;
    int flankOpen = 0;            // Ns an open gap gets
    size_t flankNStart = 0;       // where the Ns of that gap begin in seq
    size_t gfGaps = 0, gfKept = 0, gfOverlapGaps = 0, gfOverlapBases = 0, gfFloor1 = 0;
    size_t gfUnknown = 0, gfOverlapRejected = 0;
    // T24 (TESSERACT_FIX_COV_CONTRIB, default off): shadow sums of both weightings.
    double cvRelW = 0, cvConW = 0;
    size_t cvRelL = 0, cvConL = 0, cvPaths = 0, cvOff10 = 0;
    auto notePathCov = [&]() {
        if (cvRelL && cvConL) {
            const double a = cvRelW / static_cast<double>(cvRelL), b = cvConW / static_cast<double>(cvConL);
            ++cvPaths;
            if (a > 0 && b > 0 && std::max(a / b, b / a) > 1.10) ++cvOff10;
        }
        cvRelW = cvConW = 0;
        cvRelL = cvConL = 0;
    };

    ResolvedPath curPath;
    // PKG-MERGE (GAP x ENDS), build_v3: the overlap bases T03 (TESSERACT_FIX_GAP_FLANK, alias
    // TESSERACT_GAP_KEEP_FLANK) trims from the first unitig after a gap when restoreGapFlanks()
    // runs, by curPath index, so endBody() below measures the N-free bases finally emitted.
    std::vector<std::pair<size_t, size_t>> flankDrops;
    auto renderChain = [&](uint32_t c, bool flipped, std::string& seq,
                           double& covWeighted, size_t& covLen) {
        const size_t sz = chains[c].size();
        for (size_t i = 0; i < sz; ++i) {
            const uint64_t oid = flipped ? flip(chains[c][sz - 1 - i]) : chains[c][i];
            curPath.oriented.push_back(oid);
            curPath.gaps.push_back(0);
            const uint32_t u = unitigOf(oid);
            const std::string piece = g_.oriented(u, orientOf(oid));
            size_t added = 0;   // bases of this piece appended to seq (Ns excluded)
            if (seq.empty()) { seq = piece; added = piece.size(); }
            else if (flankJoin > 0) {
                // T03: the first piece after an N-gap (reached only with the fix on). The
                // A side is seq up to the Ns.
                size_t L = 0;
                const size_t lmax = std::min(ov, std::min(flankNStart, piece.size()));
                for (size_t l = lmax; l > 0; --l)
                    if (seq.compare(flankNStart - l, l, piece, 0, l) == 0) { L = l; break; }
                const long dEst = static_cast<long>(joinEst[flankPort]) - static_cast<long>(ov);
                bool trim = false;
                if (L >= kMinFlankOverlap) {
                    // An overlap the estimate does not contradict (3 standard errors, plus
                    // slack for the release estimator's short bias).
                    if (static_cast<double>(dEst + static_cast<long>(L)) <= 3.0 * joinSigma[flankPort] + 20.0)
                        trim = true;
                    else ++gfOverlapRejected;
                }
                int nN = 1;
                if (trim) {
                    ++gfOverlapGaps;
                    gfOverlapBases += L;
                    // PKG-MERGE GAP x ENDS: after restoreGapFlanks() this piece is emitted less
                    // its L overlap bases; endBody() below subtracts them (flankDrops).
                    flankDrops.emplace_back(curPath.oriented.size() - 1, L);
                }
                else if (dEst >= 1) nN = static_cast<int>(dEst);
                else if (fixGapEstimate) { nN = kUnknownGapN; ++gfUnknown; }
                else ++gfFloor1;
                const size_t from = trim ? L : 0;
                GapFlankRecord rec;
                const size_t lc = std::min(kGapFlankContext, flankNStart);
                rec.left = seq.substr(flankNStart - lc, lc);
                rec.nWritten = flankJoin;
                rec.nOpen = nN;
                rec.restore = piece.substr(from, ov - from);
                gapFlank_.push_back(std::move(rec));   // `right` is filled once the chain is rendered
                // The release layout: the Ns already written stand for this piece's first k-1.
                seq += piece.substr(ov);
                added = piece.size() - ov;
                ++gfGaps;
                gfKept += ov - from;
                flankOpen = nN;
                flankJoin = 0;
            } else {
                seq += piece.substr(ov);
                added = piece.size() - ov;
            }
            const double cov = g_.nodes[u].coverage;
            cvRelW += cov * static_cast<double>(piece.size());
            cvRelL += piece.size();
            cvConW += cov * static_cast<double>(added);
            cvConL += added;
            // T24 (TESSERACT_FIX_COV_CONTRIB): weight by the bases the piece contributes
            // (whole first piece, piece-(k-1) after, whole after an N-gap under T03), as
            // extendByCommonPrefix already does. The release weights by the full length,
            // so a k-long connector that adds one base carries k bases of weight
            // (nmeningitidis GCA002073235v2 NODE_98: 668x reported vs 62x contributed).
            const size_t w = fixCovContrib ? added : piece.size();
            covWeighted += cov * static_cast<double>(w);
            covLen += w;
            placed[u] = 1;
        }
        if (sz > 1) stats_.unitigsJoined += sz - 1;
    };

    // A contig that stops because its only exit enters a multi-entrance repeat
    // can absorb that repeat regardless of which copy it belongs to: every
    // continuation the graph offers begins with the same sequence, so appending
    // it decides nothing. This recovers the repeat copy at this locus -- which
    // would otherwise be represented once and counted as missing everywhere
    // else -- and lengthens the contig, without choosing between candidates.
    // On by default: measured over 37 closed-reference isolates this raises
    // genome fraction from 98.83% to 99.02% with the duplication ratio
    // unchanged at 1.000, because the sequence it adds is sequence every
    // candidate continuation agreed on. Set to 0 to disable.
    const double kPrefixBudget = env::real("TESSERACT_COMMON_PREFIX", 3000.0);
    // A chain owns each eligible anchor exactly once. Terminal context must
    // not make another copy of that anchor after chain arbitration declined
    // to join it. The owner is computed before any sequence is rendered, so
    // this decision does not depend on output order. Default-off experiment.
    const bool ownedAnchorPrefix = env::on("TESSERACT_OWNED_ANCHOR_PREFIX", false);
    // ---- evidence-backed terminal extension (combo2/ends, 2026-09-24; all default OFF) ----
    // TESSERACT_PREFIX_MIN_BODY=<bp>|reach  certified starting context: walk an end only when the
    //   N-free path segment it continues carries at least that many bases ('reach' = the insert
    //   model's maxPlausible, i.e. the context is long enough to anchor the fragments that would
    //   vouch for it). Measured on the 250-bp walk: 8 of the 11 Q-clean misassemblies it added
    //   sit in contigs <= 1.4 kb that the walk lengthened, not in walked sequence.
    // TESSERACT_CERTIFIED_PREFIX=1  a step into an anchorable unitig is never taken (entering an
    //   anchor is a join claim arbitration owns and declined), and a step out of a walked unitig
    //   needs depth(next) >= TESSERACT_CERTIFIED_PREFIX_FLOW (default 0.8) * depth(current): fewer
    //   copies leaving a repeat than entering it means the other copies' exits are not in the graph.
    // build_v3: read per call through the envflags table (T16), never cached in a static.
    const long kPrefixMinBody = [] {           // 0 = off, -1 = fragment reach
        const char* e = env::text("TESSERACT_PREFIX_MIN_BODY");
        if (!e) return 0L;
        if (std::string(e) == "reach") return -1L;
        return std::strtol(e, nullptr, 10);   // validated: an integer in [0, INT_MAX]
    }();
    const bool kCertifiedPrefix = env::on("TESSERACT_CERTIFIED_PREFIX", false);
    const double kCertifiedFlow = env::real("TESSERACT_CERTIFIED_PREFIX_FLOW", 0.8);
    // TESSERACT_PREFIX_NO_REVISIT=1: a walk never enters a unitig (either orientation) already on this
    //   contig's path -- re-entering emitted sequence is a hairpin/cycle claim, not a common prefix.
    const bool kPrefixNoRevisit = env::on("TESSERACT_PREFIX_NO_REVISIT", false);
    size_t noRevisitStops = 0;
    const size_t minBodyBp = kPrefixMinBody < 0 ? static_cast<size_t>(reach)
                                                : static_cast<size_t>(kPrefixMinBody);
    size_t ctpAnchorStops = 0, ctpFlowStops = 0, minBodySkips = 0, minBodyChecked = 0;
    // PKG-MERGE self-check, JOIN_TRACE + MIN_BODY only (stderr): endBody(tail) vs the N-free tail of seq.
    size_t ebTails = 0, ebExact = 0, ebReleaseGapCut = 0, ebOther = 0;
    size_t pathBeforeTailWalk = 0;
    // N-free bases of curPath adjacent to one end, before any walk was appended.
    // PKG-MERGE: the sum counts one element in full and the rest at addedLen. Behind a scaffold gap
    // that is the N-free length finally emitted when T03 (TESSERACT_FIX_GAP_FLANK, alias
    // TESSERACT_GAP_KEEP_FLANK) gives the first unitig after the gap its k-1 bases back
    // (restoreGapFlanks, after gap filling), less the verified overlap T03 trims there
    // (flankDrops). With the fix off the formula is unchanged, so it still counts the k-1 bases
    // release drops after a gap, and ENDS alone is byte-identical.
    auto endBody = [&](bool head) -> size_t {
        const std::vector<uint64_t>& P = curPath.oriented;
        const std::vector<int>& Gp = curPath.gaps;
        const size_t n = std::min(pathBeforeTailWalk, P.size());
        size_t bp = 0;
        if (n == 0) return 0;
        if (!head) {
            for (size_t i = n; i-- > 0;) {
                bp += bp == 0 ? g_.nodes[unitigOf(P[i])].seq.size() : addedLen(P[i]);
                if (Gp[i] > 0) {                          // an N-gap precedes element i
                    for (const auto& d : flankDrops) if (d.first == i) bp -= d.second;
                    break;
                }
            }
        } else {
            for (size_t i = 0; i < n; ++i) {
                if (i > 0 && Gp[i] > 0) break;
                bp += bp == 0 ? g_.nodes[unitigOf(P[i])].seq.size() : addedLen(P[i]);
            }
        }
        return bp;
    };
    auto bodyAllowsWalk = [&](bool head) {
        if (kPrefixMinBody == 0 || kPrefixBudget <= 0) return true;
        ++minBodyChecked;
        if (endBody(head) >= minBodyBp) return true;
        ++minBodySkips;
        return false;
    };
    // TESSERACT_PAIR_ANCHORED_PREFIX=<m>|bar : a walk is kept only as far as read pairs vouch for it. A pair
    //   vouches when one mate is anchored (buildSupport) on a non-repeat unitig of the N-free contig body within
    //   fragment reach of the end, runs outward, and its mate -- reverse-complemented onto the outward strand --
    //   matches the body+walk window exactly past the body end, with a plausible fragment length.
    //   TESSERACT_PAIR_ANCHORED_MODE=node (default): a walked unitig is kept whole iff >= m mates reach into its
    //   novel bases and every unitig before it was kept (a unitig does not branch inside). =pos: the walk is cut
    //   at the m-th largest mate end, possibly inside a unitig (that unitig is then left out of the GFA path).
    const long paSetting = pairAnchoredPrefixSetting();
    const int paM = paSetting < 0 ? std::max(1, static_cast<int>(std::lround(linkBar)))
                                  : static_cast<int>(paSetting);
    const bool paPosMode = [] {
        const char* e = env::text("TESSERACT_PAIR_ANCHORED_MODE");   // node (default) | pos
        return e && std::string(e) == "pos";
    }();
    std::vector<uint32_t> paStart, paReads;       // CSR: reads anchored on each unitig
    if (paSetting != 0 && !readAnchors_.empty()) {
        paStart.assign(n + 1, 0);
        for (const Anchor& a : readAnchors_) if (a.mapped() && a.unitig < n) ++paStart[a.unitig + 1];
        for (size_t u = 0; u < n; ++u) paStart[u + 1] += paStart[u];
        paReads.resize(paStart[n]);
        std::vector<uint32_t> fill(paStart.begin(), paStart.end() - 1);
        for (size_t r = 0; r < readAnchors_.size(); ++r) {
            const Anchor& a = readAnchors_[r];
            if (a.mapped() && a.unitig < n) paReads[fill[a.unitig]++] = static_cast<uint32_t>(r);
        }
    }
    size_t paEnds = 0, paCapped = 0, paZero = 0, paBpWalk = 0, paBpKept = 0, paMates = 0;
    // body: outward-oriented N-free body (the last bases are adjacent to the walk); elems: outward oriented
    // unitigs of that body with their start offset in `body`; novel: the walked bases; cum: cumulative novel
    // bases per walked unitig. Returns the number of novel bases to keep.
    auto paKeep = [&](const std::string& body, const std::vector<std::pair<uint64_t, long>>& elems,
                      const std::string& novel, const std::vector<size_t>& cum) -> size_t {
        ++paEnds;
        paBpWalk += novel.size();
        const long reachL = static_cast<long>(reach);
        const size_t bw = std::min(body.size(), static_cast<size_t>(reachL) + 200);
        const long shift = static_cast<long>(body.size() - bw);          // body offset -> window offset
        const std::string W = body.substr(body.size() - bw) + novel;
        std::unordered_map<Kmer, int32_t, KmerHasher> wk;
        {
            Kmer f = 0; int valid = 0;
            for (size_t p = 0; p < W.size(); ++p) {
                const int c = baseCode(W[p]);
                if (c < 0) { valid = 0; continue; }
                f = pushBack(f, c, kMap_);
                if (++valid < kMap_) continue;
                const int32_t start = static_cast<int32_t>(p + 1 - static_cast<size_t>(kMap_));
                auto it = wk.find(f);
                if (it == wk.end()) wk.emplace(f, start); else it->second = -1;
            }
        }
        std::vector<long> ends;
        std::string mate;
        for (const auto& el : elems) {
            const uint32_t u = unitigOf(el.first);
            if (isRepeat(u) || u + 1 >= paStart.size()) continue;
            const long uLen = static_cast<long>(g_.nodes[u].seq.size());
            const long elStartW = el.second - shift;
            if (elStartW + uLen < 0) continue;
            for (uint32_t x = paStart[u]; x < paStart[u + 1]; ++x) {
                const uint32_t r = paReads[x];
                if (r >= 2 * reads_.pairCount()) continue;
                const Anchor& a = readAnchors_[r];
                const int d = orientOf(el.first);
                if (static_cast<int>(a.orient) != d) continue;            // must run outward
                const long rl = static_cast<long>(reads_.length(r));
                const long posOut = d == 0 ? a.pos : uLen - (a.pos + rl);
                const long aStart = elStartW + posOut;
                if (aStart < 0 || aStart >= static_cast<long>(bw)) continue;
                reads_.decode(r ^ 1, mate);
                const std::string m = reverseComplement(mate);         // outward strand
                const int ml = static_cast<int>(m.size());
                if (ml < kMap_) continue;
                // diagonal vote over spread probes, then exact extension from the last agreeing probe
                const int span = ml - kMap_;
                const int probes = std::min(kMaxProbes, span + 1);
                long bestDiag = 0; int bestVotes = 0, lastProbe = -1;
                long diags[kMaxProbes]; int votes[kMaxProbes]; int lastp[kMaxProbes]; int nd = 0;
                for (int t = 0; t < probes; ++t) {
                    const int rp = probes == 1 ? 0 : span * t / (probes - 1);
                    Kmer f = 0; bool ok = true;
                    for (int j = 0; j < kMap_; ++j) {
                        const int c = baseCode(m[static_cast<size_t>(rp + j)]);
                        if (c < 0) { ok = false; break; }
                        f = pushBack(f, c, kMap_);
                    }
                    if (!ok) continue;
                    auto it = wk.find(f);
                    if (it == wk.end() || it->second < 0) continue;
                    const long dg = static_cast<long>(it->second) - rp;
                    int q = 0;
                    while (q < nd && diags[q] != dg) ++q;
                    if (q == nd) { diags[nd] = dg; votes[nd] = 0; lastp[nd] = rp; ++nd; }
                    ++votes[q]; lastp[q] = std::max(lastp[q], rp);
                    if (votes[q] > bestVotes) { bestVotes = votes[q]; bestDiag = dg; lastProbe = lastp[q]; }
                }
                if (bestVotes < kMinVotes) continue;
                const long frag = bestDiag + ml - aStart;
                if (insert_.usable && (frag < insert_.minPlausible || frag > insert_.maxPlausible)) continue;
                long i = lastProbe + kMap_;
                while (i < ml && bestDiag + i < static_cast<long>(W.size()) &&
                       m[static_cast<size_t>(i)] == W[static_cast<size_t>(bestDiag + i)]) ++i;
                const long bEnd = bestDiag + i;
                if (bEnd <= static_cast<long>(bw)) continue;
                ends.push_back(bEnd - static_cast<long>(bw));
                ++paMates;
            }
        }
        std::sort(ends.begin(), ends.end(), std::greater<long>());
        size_t keep = 0;
        if (static_cast<int>(ends.size()) >= paM && paM > 0) {
            if (paPosMode) {
                keep = std::min(novel.size(), static_cast<size_t>(ends[static_cast<size_t>(paM - 1)]));
            } else {
                size_t acc = 0;
                for (size_t c : cum) {
                    size_t reachIn = 0;
                    for (long e : ends) { if (e > static_cast<long>(acc)) ++reachIn; else break; }
                    if (static_cast<int>(reachIn) < paM) break;
                    acc = c;
                }
                keep = acc;
            }
        }
        if (keep < novel.size()) ++paCapped;
        if (keep == 0) ++paZero;
        paBpKept += keep;
        return keep;
    };
    // Truncates a walk that extendByCommonPrefix appended: `ext` holds `novelBefore` bases before the walk,
    // curPath holds `pathBefore` elements before it.
    auto paApply = [&](bool head, std::string& ext, size_t novelBefore, size_t pathBefore,
                       double& covW, size_t& covL, const std::string& scaffold, size_t scaffoldEnd) {
        if (paSetting == 0 || readAnchors_.empty()) return;
        if (curPath.oriented.size() <= pathBefore) return;
        std::vector<uint64_t> walked(curPath.oriented.begin() + static_cast<std::ptrdiff_t>(pathBefore),
                                     curPath.oriented.end());
        std::vector<size_t> cum;
        size_t c = 0;
        for (uint64_t o : walked) { c += addedLen(o); cum.push_back(c); }
        const std::string novel = ext.substr(novelBefore);
        // outward N-free body and its outward-oriented elements with start offsets in it
        std::string body;
        std::vector<std::pair<uint64_t, long>> elems;
        const std::vector<uint64_t>& P = curPath.oriented;
        const std::vector<int>& Gp = curPath.gaps;
        if (!head) {
            size_t lastN = scaffold.find_last_of("Nn", scaffoldEnd ? scaffoldEnd - 1 : 0);
            const size_t from = (lastN == std::string::npos || scaffoldEnd == 0) ? 0 : lastN + 1;
            body = scaffold.substr(from, scaffoldEnd - from);
            long end = static_cast<long>(body.size());
            for (size_t i = pathBefore; i-- > 0;) {
                const long len = static_cast<long>(g_.nodes[unitigOf(P[i])].seq.size());
                elems.push_back({P[i], end - len});
                end -= len - static_cast<long>(ov);
                if (Gp[i] > 0 || end - static_cast<long>(ov) < 0) break;
            }
        } else {
            size_t firstN = scaffold.find_first_of("Nn");
            const size_t to = std::min(firstN == std::string::npos ? scaffold.size() : firstN, scaffoldEnd);
            body = reverseComplement(scaffold.substr(0, to));
            long startFwd = 0;
            for (size_t i = 0; i < pathBefore; ++i) {
                if (i > 0 && Gp[i] > 0) break;
                const long len = static_cast<long>(g_.nodes[unitigOf(P[i])].seq.size());
                if (startFwd + len > static_cast<long>(to)) break;
                elems.push_back({flip(P[i]), static_cast<long>(to) - (startFwd + len)});
                startFwd += len - static_cast<long>(ov);
            }
        }
        const size_t keep = paKeep(body, elems, novel, cum);
        if (keep >= novel.size()) return;
        // drop the walked unitigs (and their coverage) beyond `keep`; a partial unitig leaves the path
        size_t acc = 0, keepNodes = 0;
        for (size_t i = 0; i < walked.size(); ++i) {
            const size_t add = addedLen(walked[i]);
            if (acc + add <= keep) { acc += add; keepNodes = i + 1; continue; }
            const size_t removed = acc + add - std::max(acc, keep);
            covW -= g_.nodes[unitigOf(walked[i])].coverage * static_cast<double>(removed);
            covL -= removed;
            acc += add;
        }
        ext.resize(novelBefore + keep);
        curPath.oriented.resize(pathBefore + keepNodes);
        curPath.gaps.resize(pathBefore + keepNodes);
    };
    auto extendByCommonPrefix = [&](uint64_t tail, std::string& seq,
                                    double& covWeighted, size_t& covLen, bool head) {
        if (kPrefixBudget <= 0) return;
        size_t added = 0;
        uint64_t cur = tail;
        std::vector<uint64_t> visited{tail};
        for (int step = 0; step < 8; ++step) {
            const auto& exits = g_.exits(unitigOf(cur), orientOf(cur));
            std::vector<uint64_t> next;
            bool bubble = false;
            if (exits.size() == 1) {
                const Link& l = exits[0];
                if (g_.nodes[l.to].deleted) break;
                next.push_back(orientedId(l.to, UnitigGraph::enterOrient(l)));
            } else if (completePrefixBubbles && exits.size() == 2 && isRepeat(unitigOf(cur))) {
                // A substitution inside a collapsed repeat creates two arms of
                // length 2k-1 differing at one base, with the same destination.
                // Stopping here emits an arbitrary half-repeat. Completing the
                // allele bubble retains the repeat context without picking its
                // unique exit. The higher k-mer-depth allele is a consensus
                // choice, not a phased claim; evaluate mismatches separately.
                uint64_t arm[2], dest[2];
                std::string allele[2];
                bool valid = true;
                for (size_t i = 0; i < 2; ++i) {
                    const Link& l = exits[i];
                    if (g_.nodes[l.to].deleted || !isRepeat(l.to)) { valid = false; break; }
                    arm[i] = orientedId(l.to, UnitigGraph::enterOrient(l));
                    const auto& out = g_.exits(l.to, orientOf(arm[i]));
                    if (out.size() != 1 || g_.nodes[out[0].to].deleted ||
                        !isRepeat(out[0].to)) { valid = false; break; }
                    dest[i] = orientedId(out[0].to, UnitigGraph::enterOrient(out[0]));
                    allele[i] = g_.oriented(l.to, orientOf(arm[i]));
                    if (allele[i].size() != static_cast<size_t>(2 * k_ - 1)) {
                        valid = false;
                        break;
                    }
                }
                if (!valid || dest[0] != dest[1] || arm[0] == arm[1]) break;
                size_t differences = 0;
                for (size_t i = 0; i < allele[0].size(); ++i)
                    differences += allele[0][i] != allele[1][i];
                if (differences != 1) break;
                const double c0 = g_.nodes[unitigOf(arm[0])].coverage;
                const double c1 = g_.nodes[unitigOf(arm[1])].coverage;
                const size_t pick = c1 > c0 || (c1 == c0 && arm[1] < arm[0]) ? 1 : 0;
                next = {arm[pick], dest[pick]};
                bubble = true;
            } else break;

            if (kPrefixNoRevisit) {
                bool again = false;
                for (uint64_t x : next) {
                    const uint32_t xu = unitigOf(x);
                    for (uint64_t y : curPath.oriented) {
                        if (unitigOf(y) == xu) { again = true; break; }
                    }
                    if (again) break;
                }
                if (again) { ++noRevisitStops; break; }
            }
            if (kCertifiedPrefix) {
                bool certified = true;
                uint64_t pred = cur;
                for (size_t i = 0; i < next.size(); ++i) {
                    const uint32_t nu = unitigOf(next[i]);
                    if (anchorable(nu)) { ++ctpAnchorStops; certified = false; break; }
                    const bool predWalked = step > 0 || i > 0;
                    if (predWalked &&
                        g_.nodes[nu].coverage < kCertifiedFlow * g_.nodes[unitigOf(pred)].coverage) {
                        ++ctpFlowStops;
                        certified = false;
                        break;
                    }
                    pred = next[i];
                }
                if (!certified) break;
            }

            size_t extra = 0;
            bool valid = true;
            for (uint64_t nxt : next) {
                const size_t length = g_.nodes[unitigOf(nxt)].seq.size();
                if (length <= ov || unitigOf(nxt) == unitigOf(cur)) { valid = false; break; }
                if (ownedAnchorPrefix && anchorable(unitigOf(nxt)) &&
                    ownerChain[unitigOf(nxt)] != UINT32_MAX) {
                    valid = false;
                    break;
                }
                // Keep the baseline path unchanged when the experiment is off.
                if (completePrefixBubbles &&
                    std::find(visited.begin(), visited.end(), nxt) != visited.end()) {
                    valid = false;
                    break;
                }
                extra += length - ov;
            }
            if (!valid || added + extra > static_cast<size_t>(kPrefixBudget)) break;
            if (joinTrace_) {
                std::fprintf(stderr, "[prefixtrace] head=%d from=%llu addedBefore=%zu bubble=%d route=",
                    head ? 1 : 0, static_cast<unsigned long long>(cur), added, bubble ? 1 : 0);
                for (size_t i = 0; i < next.size(); ++i)
                    std::fprintf(stderr, "%s%llu", i ? "," : "",
                                 static_cast<unsigned long long>(next[i]));
                std::fprintf(stderr, " added=%zu\n", extra);
            }
            for (uint64_t nxt : next) {
                const std::string piece = g_.oriented(unitigOf(nxt), orientOf(nxt));
                seq += piece.substr(ov);
                covWeighted += g_.nodes[unitigOf(nxt)].coverage *
                               static_cast<double>(piece.size() - ov);
                covLen += piece.size() - ov;
                added += piece.size() - ov;
                curPath.oriented.push_back(nxt);
                curPath.gaps.push_back(0);
                visited.push_back(nxt);
                cur = nxt;
            }
            if (!completePrefixBubbles && g_.exits(unitigOf(cur), orientOf(cur)).size() != 1) break;
        }
    };

    for (uint32_t c = 0; c < nc; ++c) {
        if (chains[c].empty() || done[c]) continue;
        // Start scaffolds at a free port so each path is walked from one end.
        if (joinTo[c * 2 + 0] != UINT32_MAX && joinTo[c * 2 + 1] != UINT32_MAX) continue;

        std::string seq;
        double covWeighted = 0;
        size_t covLen = 0;
        curPath = ResolvedPath();
        flankJoin = 0;
        cvRelW = cvConW = 0;
        cvRelL = cvConL = 0;
        flankDrops.clear();
        uint32_t cur = c;
        // Entering through the port that has no join leaves the other free to
        // continue; a chain joined only at port 1 is therefore used forward.
        int inPort = (joinTo[c * 2 + 0] == UINT32_MAX) ? 0 : 1;
        size_t steps = 0;

        int pendingGap = 0;
        while (true) {
            done[cur] = 1;
            const size_t before = curPath.oriented.size();
            renderChain(cur, inPort == 1, seq, covWeighted, covLen);
            if (fixGapFlank && pendingGap > 0) {
                // T03: the path describes the open-gap layout restoreGapFlanks() gives.
                pendingGap = flankOpen;
                const size_t rs = flankNStart + static_cast<size_t>(gapFlank_.back().nWritten);
                gapFlank_.back().right = seq.substr(rs, std::min(kGapFlankContext, seq.size() - rs));
            }
            if (pendingGap > 0 && before < curPath.gaps.size()) curPath.gaps[before] = pendingGap;
            pendingGap = 0;
            const uint32_t outPort = cur * 2 + (inPort == 0 ? 1 : 0);
            const uint32_t partner = joinTo[outPort];
            if (partner == UINT32_MAX || ++steps > nc) {
                pathBeforeTailWalk = curPath.oriented.size();
                if (joinTrace_ && kPrefixMinBody != 0 && !curPath.oriented.empty()) {
                    const size_t lastN = seq.find_last_of("Nn");
                    size_t actual = lastN == std::string::npos ? seq.size() : seq.size() - lastN - 1;
                    // T03 keeps the release layout here and restores the flank after gap filling:
                    // compare against the tail as restoreGapFlanks() will emit it.
                    if (fixGapFlank && lastN != std::string::npos && !gapFlank_.empty())
                        actual += gapFlank_.back().restore.size();
                    const size_t eb = endBody(false);
                    ++ebTails;
                    if (eb == actual) ++ebExact;
                    else if (eb == actual + ov) ++ebReleaseGapCut;   // T03 off, behind a gap
                    else ++ebOther;
                }
                if (!curPath.oriented.empty() && bodyAllowsWalk(false)) {
                    const size_t seqBefore = seq.size();
                    const std::string scaffoldCopy = paSetting != 0 ? seq : std::string();
                    extendByCommonPrefix(curPath.oriented.back(), seq, covWeighted, covLen, false);
                    paApply(false, seq, seqBefore, pathBeforeTailWalk, covWeighted, covLen,
                            scaffoldCopy, seqBefore);
                }
                break;
            }
            const uint32_t nextC = partner / 2;
            if (done[nextC] || chains[nextC].empty()) break;
            if (fixGapFlank) {
                // T03: the next chain's first piece records how this gap is restored.
                flankJoin = joinGap[outPort];
                flankPort = outPort;
                flankNStart = seq.size();
            }
            seq.append(static_cast<size_t>(joinGap[outPort]), 'N');
            stats_.gapBases += static_cast<size_t>(joinGap[outPort]);
            ++stats_.scaffoldJoins;
            pendingGap = joinGap[outPort];
            cur = nextC;
            inPort = static_cast<int>(partner % 2);
        }

        if (seq.empty()) continue;
        // Both ends of a contig stop for the same reason, so both deserve the
        // same treatment. The walk extends the tail as it finishes; the head is
        // reached by walking outward from its flipped first unitig and
        // prepending what comes back. On one plasmid, seventeen gaps -- every
        // base it was short -- were repeat copies sitting just off a contig end.
        if (!curPath.oriented.empty()) {
            std::string headExt;
            double hCov = 0;
            size_t hLen = 0;
            const size_t beforeHead = curPath.oriented.size();
            if (bodyAllowsWalk(true)) {
                extendByCommonPrefix(flip(curPath.oriented.front()), headExt, hCov, hLen, true);
                paApply(true, headExt, 0, beforeHead, hCov, hLen, seq, seq.size());
            }
            if (!headExt.empty()) {
                // extendByCommonPrefix already drops each piece's overlap as it
                // appends, so headExt is novel sequence only and is prepended
                // whole -- exactly as the tail case appends whole. Subtracting
                // the overlap a second time here would delete k-1 bases that are
                // really in the genome.
                seq.insert(0, reverseComplement(headExt));
                covWeighted += hCov;
                covLen += hLen;
            }
            // That walk ran outward from the front in the flipped frame, so the
            // unitigs it recorded landed at the end of the path facing the wrong
            // way. Move them to the front, reversed and flipped, or the GFA walk
            // describes a contig the FASTA does not contain.
            if (curPath.oriented.size() > beforeHead) {
                std::vector<uint64_t> head(curPath.oriented.begin() +
                                               static_cast<std::ptrdiff_t>(beforeHead),
                                           curPath.oriented.end());
                curPath.oriented.resize(beforeHead);
                curPath.gaps.resize(beforeHead);
                for (uint64_t& u : head) u = flip(u);
                std::reverse(head.begin(), head.end());
                curPath.oriented.insert(curPath.oriented.begin(), head.begin(), head.end());
                curPath.gaps.insert(curPath.gaps.begin(), head.size(), 0);
            }
        }
        contigs.push_back(std::move(seq));
        covs.push_back(covLen ? covWeighted / static_cast<double>(covLen) : 0);
        notePathCov();
        paths_.push_back(curPath);
        ++stats_.pathsBuilt;
    }

    // Chains inside a scaffold cycle have no free port. Release 1.3.0 renders them here
    // one by one, dropping every join of the cycle (T13; counted in [scafcycle], and
    // TESSERACT_FIX_SCAFFOLD_CYCLE breaks each cycle before the walk above instead).
    for (uint32_t c = 0; c < nc; ++c) {
        if (chains[c].empty() || done[c]) continue;
        std::string seq;
        double covWeighted = 0;
        size_t covLen = 0;
        done[c] = 1;
        curPath = ResolvedPath();
        flankJoin = 0;
        cvRelW = cvConW = 0;
        cvRelL = cvConL = 0;
        flankDrops.clear();
        renderChain(c, false, seq, covWeighted, covLen);
        if (seq.empty()) continue;
        contigs.push_back(std::move(seq));
        covs.push_back(covLen ? covWeighted / static_cast<double>(covLen) : 0);
        notePathCov();
        paths_.push_back(curPath);
        ++stats_.pathsBuilt;
    }

    // Anything no chain walked through would otherwise be lost.
    for (uint32_t u = 0; u < n; ++u) {
        if (g_.nodes[u].deleted || placed[u]) continue;
        contigs.push_back(g_.nodes[u].seq);
        covs.push_back(g_.nodes[u].coverage);
        ResolvedPath solo;
        solo.oriented.push_back(orientedId(u, 0));
        solo.gaps.push_back(0);
        paths_.push_back(solo);
    }

    ResolverCounters rc;
    rc.run = 1;
    rc.requireSupportSingle = requireSupportSingle_ ? 1 : 0;
    rc.minFallbackDest = minFallbackDest_;
    rc.excludeSharedRepeat = excludeSharedRepeatSupport ? 1 : 0;
    rc.sharedAudit = auditSharedRepeatSupport ? 1 : 0;
    rc.noUnspannedFallback = noUnspannedFallback_ ? 1 : 0;
    rc.linkBar = linkBar;
    rc.tieRatio = tieRatio_;
    rc.evals = dbgEvals; rc.ok = dbgOk; rc.noCand = dbgNoCand; rc.noPick = dbgLowSupport; rc.tie = dbgTie;
    rc.midChain = dbgMidChain; rc.byCoverage = dbgCoverage; rc.matched = dbgMatched;
    rc.shortDest = dbgShortDest; rc.unspanned = dbgUnspanned; rc.loneRepeat = dbgLoneRepeat;
    rc.sharedExcludedAll = dbgEsrsAllExcluded; rc.esrsChanged = dbgEsrsChanged;
    rc.truncGuard = fixTruncGuard ? 1 : 0;
    rc.trSides = trSides; rc.trBudget = trBudget; rc.trMaxNodes = trMaxNodes; rc.trCap = trCap;
    rc.trLone = trLone; rc.trLoneBelowBar = trLoneBelowBar; rc.trRefused = trRefused;
    rc.revisitGuard = fixRevisitGuard ? 1 : 0;
    rc.rvPrunedSides = rvPrunedSides; rc.rvPicked = rvPicked; rc.rvRefused = rvRefused;
    rc.cycleBreak = fixCycleBreak ? 1 : 0;
    rc.scCycles = scCycles; rc.scChains = scChains; rc.scJoins = scJoins;
    rc.scDropped = fixCycleBreak ? scBroken : scJoins; rc.scBases = scBases; rc.scBroken = scBroken;
    rc.gapEstimate = fixGapEstimate ? 1 : 0;
    rc.geJoins = geJoins; rc.geEstimated = geEstimated; rc.geChanged = geChanged; rc.geFloor1 = geFloor1;
    rc.geOverlapEst = geOverlapEst; rc.geSkippedNoMean = geSkippedNoMean;
    if (!geShift.empty()) {
        std::nth_element(geShift.begin(), geShift.begin() + static_cast<std::ptrdiff_t>(geShift.size() / 2), geShift.end());
        rc.geMedianShift = geShift[geShift.size() / 2];
    }
    rc.gapFlank = fixGapFlank ? 1 : 0;
    rc.gfGaps = gfGaps; rc.gfKept = gfKept; rc.gfOverlapGaps = gfOverlapGaps; rc.gfOverlapBases = gfOverlapBases;
    rc.gfOverlapRejected = gfOverlapRejected; rc.gfFloor1 = gfFloor1; rc.gfUnknown = gfUnknown;
    rc.covContrib = fixCovContrib ? 1 : 0;
    rc.cvPaths = cvPaths; rc.cvOff10 = cvOff10;
    rc.routeOrder = fixRouteOrder ? 1 : 0;
    rc.roDecisions = roDecisions; rc.roSorted = roSorted;
    rc.mirrorRoute = fixMirrorRoute ? 1 : 0;
    rc.mrDisagree = mrDisagree; rc.mrChoseBack = mrChoseBack;
    printResolverCounters(rc);
    // build_v3: the ENDS counters are printed on every run, zeros included, with a trailing
    // enabled= field (OBJECTIVE A2); the leading fields are the package's own.
    std::fprintf(stderr,
        "[pairprefix] m=%d mode=%s anchors=%zu ends=%zu capped=%zu zero=%zu bp_walk=%zu bp_kept=%zu mates=%zu "
        "enabled=%d\n",
        paSetting != 0 ? paM : 0, paPosMode ? "pos" : "node", readAnchors_.size(), paEnds, paCapped, paZero,
        paBpWalk, paBpKept, paMates, paSetting != 0 ? 1 : 0);
    std::fprintf(stderr,
        "[endext] budget=%.0f min_body=%zu ends_checked=%zu ends_skipped_short_body=%zu "
        "certified=%d flow=%.2f anchor_stops=%zu flow_stops=%zu no_revisit=%d revisit_stops=%zu enabled=%d\n",
        kPrefixBudget, kPrefixMinBody != 0 ? minBodyBp : 0, minBodyChecked, minBodySkips,
        kCertifiedPrefix ? 1 : 0, kCertifiedFlow, ctpAnchorStops, ctpFlowStops,
        kPrefixNoRevisit ? 1 : 0, noRevisitStops, (kPrefixMinBody != 0 || kCertifiedPrefix || kPrefixNoRevisit) ? 1 : 0);
    if (joinTrace_ && kPrefixMinBody != 0) {
        std::fprintf(stderr, "[endbody-check] tails=%zu exact=%zu releaseGapCut=%zu other=%zu\n",
                     ebTails, ebExact, ebReleaseGapCut, ebOther);
    }
}

void restoreGapFlanks(std::vector<std::string>& seqs, const std::vector<GapFlankRecord>& records) {
    // Key: left|right|nWritten in the resolver's orientation (forward) and the reverse
    // complement of the whole gap (reverse). A key two records share is ambiguous and
    // restores nothing.
    struct Hit { size_t rec; bool reverse; bool ambiguous; };
    std::unordered_map<std::string, Hit> index;
    auto add = [&](const std::string& key, size_t r, bool rev) {
        auto it = index.find(key);
        if (it == index.end()) index.emplace(key, Hit{r, rev, false});
        else if (it->second.rec != r || it->second.reverse != rev) it->second.ambiguous = true;
    };
    for (size_t r = 0; r < records.size(); ++r) {
        const GapFlankRecord& g = records[r];
        if (g.left.size() != kGapFlankContext || g.right.size() != kGapFlankContext) continue;
        const std::string n = std::to_string(g.nWritten);
        add(g.left + "|" + g.right + "|" + n, r, false);
        add(reverseComplement(g.right) + "|" + reverseComplement(g.left) + "|" + n, r, true);
    }
    size_t runs = 0, restored = 0, ambiguous = 0, bases = 0;
    long long nDelta = 0;
    std::vector<char> used(records.size(), 0);
    for (std::string& s : seqs) {
        if (s.find('N') == std::string::npos) continue;
        std::string out;
        out.reserve(s.size() + 256);
        size_t pos = 0;
        while (pos < s.size()) {
            if (s[pos] != 'N') { out.push_back(s[pos++]); continue; }
            size_t e = pos;
            while (e < s.size() && s[e] == 'N') ++e;
            ++runs;
            bool done = false;
            if (pos >= kGapFlankContext && e + kGapFlankContext <= s.size()) {
                const std::string key = s.substr(pos - kGapFlankContext, kGapFlankContext) + "|" +
                                        s.substr(e, kGapFlankContext) + "|" + std::to_string(e - pos);
                auto it = index.find(key);
                if (it != index.end() && it->second.ambiguous) ++ambiguous;
                else if (it != index.end()) {
                    const GapFlankRecord& g = records[it->second.rec];
                    if (it->second.reverse) {
                        out += reverseComplement(g.restore);
                        out.append(static_cast<size_t>(g.nOpen), 'N');
                    } else {
                        out.append(static_cast<size_t>(g.nOpen), 'N');
                        out += g.restore;
                    }
                    ++restored;
                    used[it->second.rec] = 1;
                    bases += g.restore.size();
                    nDelta += static_cast<long long>(g.nOpen) - static_cast<long long>(e - pos);
                    done = true;
                }
            }
            if (!done) out.append(s, pos, e - pos);
            pos = e;
        }
        s.swap(out);
    }
    size_t unmatched = 0;
    for (char u : used) unmatched += u ? 0 : 1;
    std::fprintf(stderr,
        "[gapflank-restore] records=%zu nRuns=%zu restored=%zu restoredBases=%zu nDelta=%lld "
        "unmatchedRecords=%zu ambiguous=%zu\n",
        records.size(), runs, restored, bases, nDelta, unmatched, ambiguous);
}

void printResolverCountersNotRun() {
    // Same lines as a resolver run, run=0, zero counts; the flag values are the ones the
    // environment resolves to (release flags parsed as release parses them).
    ResolverCounters rc;
    rc.run = 0;
    rc.requireSupportSingle = env::on("TESSERACT_REQUIRE_SUPPORT_SINGLE", false) ? 1 : 0;
    rc.minFallbackDest = static_cast<size_t>(env::integer("TESSERACT_MIN_FALLBACK_DEST", 0));
    rc.excludeSharedRepeat = env::on("TESSERACT_EXCLUDE_SHARED_REPEAT_SUPPORT", false) ? 1 : 0;
    rc.sharedAudit = env::on("TESSERACT_SHARED_SUPPORT_AUDIT", false) ? 1 : 0;
    rc.noUnspannedFallback = env::on("TESSERACT_NO_UNSPANNED_FALLBACK", false) ? 1 : 0;
    rc.truncGuard = fixEnabled("TESSERACT_FIX_TRUNC_GUARD") ? 1 : 0;
    rc.revisitGuard = fixEnabled("TESSERACT_FIX_REVISIT_GUARD") ? 1 : 0;
    rc.cycleBreak = fixEnabled("TESSERACT_FIX_SCAFFOLD_CYCLE") ? 1 : 0;
    rc.gapEstimate = fixEnabled("TESSERACT_FIX_GAP_ESTIMATE") ? 1 : 0;
    rc.gapFlank = fixEnabled("TESSERACT_FIX_GAP_FLANK") ? 1 : 0;
    rc.covContrib = fixEnabled("TESSERACT_FIX_COV_CONTRIB") ? 1 : 0;
    rc.routeOrder = fixEnabled("TESSERACT_FIX_ROUTE_ORDER") ? 1 : 0;
    rc.mirrorRoute = fixEnabled("TESSERACT_FIX_MIRROR_ROUTE") ? 1 : 0;
    std::fprintf(stderr,
        "[resolver] run=0 theta=0.0000 repeatThreshold=0.0000 estimator=none legacyMedian=0.0000 population=0\n");
    printResolverCounters(rc);
    // build_v3: the package resolver counters, run=0 form (same lines, zeros, flag state).
    const long long oneSidedLong = env::integer("TESSERACT_ONE_SIDED_TIE_LONG", 2000);
    const long long singleExitLong = env::integer("TESSERACT_TIE_SINGLE_EXIT_LONG", 0);
    std::fprintf(stderr,
        "[onesidedtie] weakNominations=0 weakJoins=0 refusedLongBranching=0 refusedTwoSided=0 "
        "refusedPartnerNotDecisive=0 threadPreferred=0 long=%lld enabled=%d\n",
        oneSidedLong, env::on("TESSERACT_ONE_SIDED_TIE", false) ? 1 : 0);
    std::fprintf(stderr, "[tiesingleexit] exempted=0 refusedLong=0 long=%lld enabled=%d\n", singleExitLong,
                 env::on("TESSERACT_TIE_SINGLE_EXIT", false) ? 1 : 0);
    std::fprintf(stderr, "[fallbackmatching] kept=0 refused=0 enabled=%d\n",
                 env::on("TESSERACT_FALLBACK_MATCHING_ONLY", false) ? 1 : 0);
    const bool paOn = pairAnchoredPrefixSetting() != 0;
    const char* paMode = env::text("TESSERACT_PAIR_ANCHORED_MODE");
    std::fprintf(stderr,
        "[pairprefix] m=0 mode=%s anchors=0 ends=0 capped=0 zero=0 bp_walk=0 bp_kept=0 mates=0 enabled=%d\n",
        paMode && std::string(paMode) == "pos" ? "pos" : "node", paOn ? 1 : 0);
    const bool minBodyOn = env::text("TESSERACT_PREFIX_MIN_BODY") != nullptr &&
                           std::string(env::text("TESSERACT_PREFIX_MIN_BODY")) != "0";
    const bool certified = env::on("TESSERACT_CERTIFIED_PREFIX", false);
    const bool noRevisit = env::on("TESSERACT_PREFIX_NO_REVISIT", false);
    std::fprintf(stderr,
        "[endext] budget=%.0f min_body=0 ends_checked=0 ends_skipped_short_body=0 certified=%d flow=%.2f "
        "anchor_stops=0 flow_stops=0 no_revisit=%d revisit_stops=0 enabled=%d\n",
        env::real("TESSERACT_COMMON_PREFIX", 3000.0), certified ? 1 : 0,
        env::real("TESSERACT_CERTIFIED_PREFIX_FLOW", 0.8), noRevisit ? 1 : 0,
        (minBodyOn || certified || noRevisit) ? 1 : 0);
}

}  // namespace ts
