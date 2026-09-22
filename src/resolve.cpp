#include "resolve.h"
#include "graph_coverage.h"
#include "resolve_evidence.h"
#include "read_thread_evidence.h"

#include "organism.h"   // forEachMarkerKmer, for the plasmid vouch below

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <mutex>
#include <cstdlib>
#include <thread>

namespace ts {

namespace {

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

}  // namespace

PairedResolver::PairedResolver(const UnitigGraph& graph, const SequenceStore& reads, int threads,
                               int minLinkSupport, double tieRatio, double linkSupportPerX,
                               int minScaffoldSupport)
    : g_(graph), reads_(reads), threads_(threads > 0 ? threads : 1), k_(graph.k()),
      kMap_(std::min(kAnchorK, graph.k())),
      minLinkSupport_(minLinkSupport), tieRatio_(tieRatio),
      linkSupportPerX_(linkSupportPerX >= 0 ? linkSupportPerX : 0.10),
      minScaffoldSupport_(minScaffoldSupport) {
    medianCoverage_ = graph.medianCoverage();
    const char* weightedCoverage = std::getenv("TESSERACT_WEIGHTED_RESOLVER_COVERAGE");
    const char* eligibleCoverage = std::getenv("TESSERACT_WEIGHTED_ELIGIBLE_COVERAGE");
    const bool useEligibleCoverage = eligibleCoverage && std::string(eligibleCoverage) == "1";
    // Explicit precedence: eligible-only calibration wins if both experiments
    // are enabled. They are alternative estimators, never cumulative weights.
    // Length-weighted by default since 1.3.0. The unweighted node median lets a cloud of
    // short carry-floor fragments set the single-copy baseline: on senterica GCA006365335v1
    // 327/980 eligible unitigs sit at ~4x and carry 2% of the length, the median reads 5.73x
    // for a genome sequenced at 35x, every unitig above 9.2x becomes a "repeat" and the
    // resolver joins almost nothing. 146 isolates: NGA50 tie -> WIN, genome-fraction losses
    // 21 -> 15. TESSERACT_WEIGHTED_RESOLVER_COVERAGE=0 restores the unweighted median.
    const bool useWeightedCoverage = !weightedCoverage || std::string(weightedCoverage) != "0";
    if (useEligibleCoverage || useWeightedCoverage) {
        const double calibrated = useEligibleCoverage
            ? graphEligibleLengthWeightedMedianCoverage(graph)
            : graphLengthWeightedMedianCoverage(graph);
        if (calibrated > 0.0) medianCoverage_ = calibrated;
    }
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

    const char* routeFlag = std::getenv("TESSERACT_ROUTE_DISTANCE");
    routeDistance_ = routeFlag && std::string(routeFlag) == "1";
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
                node.coverage <= 1.6 * medianCoverage_)
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

    auto worker = [&](int tid) {
        auto& samples = insertSamples[static_cast<size_t>(tid)];
        auto& sup = localSupport[static_cast<size_t>(tid)];
        size_t localMapped = 0, localLinks = 0;

        for (size_t p = static_cast<size_t>(tid); p < pairs; p += static_cast<size_t>(threads_)) {
            const size_t i1 = p * 2, i2 = p * 2 + 1;
            const Anchor a1 = anchorRead(i1);
            if (!a1.mapped()) continue;
            const Anchor a2 = anchorRead(i2);
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
    const double repeatThreshold = medianCoverage_ * 1.6;
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
    if (getenv("TESSERACT_DEBUG_RESOLVE") && exclusiveMarkers_) {
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

    auto isRepeat = [&](uint32_t u) {
        if (medianCoverage_ <= 0) return false;
        if (g_.nodes[u].coverage <= repeatThreshold) return false;
        if (plasmidVouched[u] && !getenv("TESSERACT_NO_PLASMID_VOUCH")) return false;
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

    const bool exactReadThreads = reads_.paired() && [] {
        const char* flag = std::getenv("TESSERACT_EXACT_READ_THREADS");
        return flag && std::string(flag) == "1";
    }();
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

    // Every way out of `from` that terminates on an anchorable unitig within
    // fragment reach; intermediate nodes are unanchorable repeats.
    auto enumerate = [&](uint64_t from, std::vector<std::vector<uint64_t>>& out) {
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
            if (f.nodes.size() >= kMaxNodes) continue;
            // A forced path (one the graph offered no alternative to) may run a
            // little past what a fragment can span, but not far: measured on the
            // benchmark panel, letting it run unbounded bought only a lower
            // contig count and cost a misassembly on K. pneumoniae.
            const double budget = f.forced ? searchBudget * kForcedReachFactor : searchBudget;
            if (static_cast<double>(f.len) > budget) continue;
            const auto& exits = g_.exits(tu, orientOf(tail));
            for (const Link& l : exits) {
                if (g_.nodes[l.to].deleted) continue;
                const uint64_t nxt = orientedId(l.to, UnitigGraph::enterOrient(l));
                bool seen = false;
                for (uint64_t o : f.nodes) if (o == nxt) { seen = true; break; }
                if (seen) continue;
                Frame g2 = f;
                g2.nodes.push_back(nxt);
                g2.len += g_.nodes[l.to].seq.size() - ov;
                g2.forced = f.forced && exits.size() == 1;
                work.push_back(std::move(g2));
            }
        }
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
        static const double dom = [] {
            const char* e = std::getenv("TESSERACT_MATCH_DOMINANCE");
            return e ? std::atof(e) : kMatchDominance;
        }();
        return matched >= 1.0 && matched >= dom * crossed;
    };

    // Complete terminal repeats across a single-substitution bubble only when
    // explicitly requested. This extends sequence; it does not relax the anchor
    // join chooser or make a choice between distinct flanking destinations.
    const bool completePrefixBubbles = [] {
        const char* e = std::getenv("TESSERACT_PREFIX_SNP_BUBBLES");
        return e && std::atoi(e) != 0;
    }();

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
    };

    // Best continuation off one end of a chain, scored with every member of
    // that chain that still lies within fragment reach of the boundary.
    long long dbgNoCand = 0, dbgLowSupport = 0, dbgTie = 0, dbgMidChain = 0, dbgOk = 0;
    long long dbgCoverage = 0, dbgMatched = 0, dbgShortDest = 0, dbgUnspanned = 0;
    long long dbgLoneRepeat = 0;
    static const bool requireSupportSingle_ = [] {
        const char* e = std::getenv("TESSERACT_REQUIRE_SUPPORT_SINGLE");
        return e && std::atoi(e) != 0;
    }();
    static const bool joinTrace_ = [] {
        const char* e = std::getenv("TESSERACT_JOIN_TRACE");
        return e && std::atoi(e) != 0;
    }();
    static const bool noUnspannedFallback_ = [] {
        const char* e = std::getenv("TESSERACT_NO_UNSPANNED_FALLBACK");
        return e && std::atoi(e) != 0;
    }();
    long long dbgTieBest = 0, dbgTieSecond = 0;

    // Withdraw the fallbacks when the continuation they would take lands on a unitig
    // shorter than this. 0 disables, which is the shipped behaviour until measured.
    static const size_t minFallbackDest_ = [] {
        const char* e = std::getenv("TESSERACT_MIN_FALLBACK_DEST");
        return e ? static_cast<size_t>(std::atoi(e)) : 0u;
    }();

    const bool excludeSharedRepeatSupport = [] {
        const char* e = std::getenv("TESSERACT_EXCLUDE_SHARED_REPEAT_SUPPORT");
        return e && std::atoi(e) != 0;
    }();
    const bool auditSharedRepeatSupport = [] {
        const char* e = std::getenv("TESSERACT_SHARED_SUPPORT_AUDIT");
        return e && std::atoi(e) != 0;
    }();

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
        enumerate(tailFirst.back(), cands);
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
                scores = adjusted;
                best = second = -1.0;
                pick = 0;
                for (size_t i = 0; i < scores.size(); ++i) {
                    const double sc = scores[i];
                    if (sc > best) { second = best; best = sc; pick = i; }
                    else if (sc > second) second = sc;
                }
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
        if (requireSupportSingle_ && cands.size() == 1 && best <= 0.0) {
            bool throughRepeat = false;
            for (size_t j = 0; j + 1 < cands[0].size(); ++j) {
                if (isRepeat(unitigOf(cands[0][j]))) { throughRepeat = true; break; }
            }
            if (throughRepeat) {
                ++dbgLoneRepeat;
                if (joinTrace_)
                    std::fprintf(stderr, "[refused] from=%u dest=%u interLen=%d\n",
                                 unitigOf(tailFirst.back()), unitigOf(cands[0].back()),
                                 [&]{ int L=0; for (size_t j=0;j+1<cands[0].size();++j)
                                      L+=static_cast<int>(g_.nodes[unitigOf(cands[0][j])].seq.size())-(k_-1);
                                      return L; }());
                return result;
            }
        }

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
            if (matchingAgrees(tailFirst.back(), cands[pick], cands[pick].back(),
                               interLenOf(cands[pick]))) {
                chosen = static_cast<int>(pick);
                ++dbgMatched;
            } else if (allSameDestination(cands)) {
                chosen = pickByInterior(cands);
            } else {
                chosen = pickByCoverage(tailFirst, cands);
            }
            if (chosen < 0) { ++dbgLowSupport; return result; }
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
                    ++dbgUnspanned;
                    return result;
                }
            }
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
            std::fprintf(stderr,
                "[jointrace] from=%u fromLen=%zu fromCov=%.1f dest=%u destLen=%zu destCov=%.1f "
                "cands=%zu best=%.1f second=%.1f interLen=%d interN=%zu interCov=%.1f "
                "byCov=%d med=%.1f maxPl=%d\n",
                fromU, g_.nodes[fromU].seq.size(), g_.nodes[fromU].coverage,
                destU, g_.nodes[destU].seq.size(), g_.nodes[destU].coverage,
                cands.size(), best, second, interLen, interN,
                interN ? interDepth / interN : 0.0, byCoverage ? 1 : 0,
                medianCoverage_, insert_.usable ? insert_.maxPlausible : -1);
        }

        if (cands.size() > 1 && !byCoverage) {
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
                    return result;
                }
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
                        const auto allocated = detail::routeDistanceAllocation(routeDensity_, spans, lengths, &contrast, linkBar);
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
                    const bool threadJoin = f.exactThread || back.exactThread;
                    if (exactReadThreads && threadJoin != (pass == 1)) continue;
                    if (threadJoin) {
                        bool same = f.connector.size() == back.connector.size();
                        for (size_t j = 0; same && j < f.connector.size(); ++j)
                            same = f.connector[j] == flip(back.connector[back.connector.size() - 1 - j]);
                        if (!same) continue;
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

                    a.insert(a.end(), f.connector.begin(), f.connector.end());
                    a.insert(a.end(), b.begin(), b.end());
                    chains[c].swap(a);
                    chains[f.chainB].clear();
                    merged[c] = 1;
                    merged[f.chainB] = 1;
                    ++joins;
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

    if (getenv("TESSERACT_DEBUG_RESOLVE")) {
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

    if (scaffolding_ && insert_.usable) {
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
        static const double scafOverride = [] {
            const char* e = std::getenv("TESSERACT_SCAF_SUPPORT");
            return e ? std::atof(e) : 0.0;
        }();
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
        }
    }

    // ---- emit ------------------------------------------------------------
    std::vector<char> placed(n, 0);
    std::vector<char> done(nc, 0);

    ResolvedPath curPath;
    auto renderChain = [&](uint32_t c, bool flipped, std::string& seq,
                           double& covWeighted, size_t& covLen) {
        const size_t sz = chains[c].size();
        for (size_t i = 0; i < sz; ++i) {
            const uint64_t oid = flipped ? flip(chains[c][sz - 1 - i]) : chains[c][i];
            curPath.oriented.push_back(oid);
            curPath.gaps.push_back(0);
            const uint32_t u = unitigOf(oid);
            const std::string piece = g_.oriented(u, orientOf(oid));
            if (seq.empty()) seq = piece;
            else seq += piece.substr(ov);
            covWeighted += g_.nodes[u].coverage * static_cast<double>(piece.size());
            covLen += piece.size();
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
    static const double kPrefixBudget = [] {
        const char* e = std::getenv("TESSERACT_COMMON_PREFIX");
        return e ? std::atof(e) : 3000.0;
    }();
    // A chain owns each eligible anchor exactly once. Terminal context must
    // not make another copy of that anchor after chain arbitration declined
    // to join it. The owner is computed before any sequence is rendered, so
    // this decision does not depend on output order. Default-off experiment.
    const bool ownedAnchorPrefix = [] {
        const char* flag = std::getenv("TESSERACT_OWNED_ANCHOR_PREFIX");
        return flag && std::string(flag) == "1";
    }();
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
            if (pendingGap > 0 && before < curPath.gaps.size()) curPath.gaps[before] = pendingGap;
            pendingGap = 0;
            const uint32_t outPort = cur * 2 + (inPort == 0 ? 1 : 0);
            const uint32_t partner = joinTo[outPort];
            if (partner == UINT32_MAX || ++steps > nc) {
                if (!curPath.oriented.empty())
                    extendByCommonPrefix(curPath.oriented.back(), seq, covWeighted, covLen, false);
                break;
            }
            const uint32_t nextC = partner / 2;
            if (done[nextC] || chains[nextC].empty()) break;
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
            extendByCommonPrefix(flip(curPath.oriented.front()), headExt, hCov, hLen, true);
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
        paths_.push_back(curPath);
        ++stats_.pathsBuilt;
    }

    // Chains inside a scaffold cycle have no free port; break them arbitrarily.
    for (uint32_t c = 0; c < nc; ++c) {
        if (chains[c].empty() || done[c]) continue;
        std::string seq;
        double covWeighted = 0;
        size_t covLen = 0;
        done[c] = 1;
        curPath = ResolvedPath();
        renderChain(c, false, seq, covWeighted, covLen);
        if (seq.empty()) continue;
        contigs.push_back(std::move(seq));
        covs.push_back(covLen ? covWeighted / static_cast<double>(covLen) : 0);
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
}

}  // namespace ts
