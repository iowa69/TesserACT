#include "polish.h"

#include <algorithm>
#include <atomic>
#include <cstdio>
#include <thread>
#include <unordered_map>
#include <vector>

#include "kmer.h"
#include "polish_quality.h"

namespace ts {

namespace {

constexpr int kMaxProbes = 12;
constexpr int kMinVotes = 2;
constexpr uint64_t kAmbiguous = UINT64_MAX;
// A read that disagrees with the contig this often is mismapped, and letting it
// vote would move the consensus toward a repeat copy that does not belong here.
constexpr double kMinReadAgreement = 0.95;

inline uint64_t packIndex(uint32_t contig, uint32_t pos, int strand) {
    return (static_cast<uint64_t>(contig) << 33) | (static_cast<uint64_t>(pos) << 1) |
           static_cast<uint64_t>(strand & 1);
}

struct AcceptedPlacement {
    uint32_t contig = UINT32_MAX;
    int32_t start = 0;
    uint8_t orientation = 0;
};

struct QualitySite {
    size_t flat = 0;
    int originalContigBase = -1;
    quality_consensus::Pileup evidence;
};

void applyOriginalQuality(std::vector<std::string>& contigs, const SequenceStore& reads,
                          const std::vector<size_t>& offset,
                          const std::vector<AcceptedPlacement>& placements,
                          std::vector<QualitySite>& sites, int minDepth, double minFraction,
                          QualityPolishStats& stats) {
    stats.candidatePositions = sites.size();
    stats.candidateBytes = sites.capacity() * sizeof(QualitySite);
    if (sites.empty()) return;
    auto lower = [&](size_t flat) {
        return static_cast<size_t>(std::lower_bound(sites.begin(), sites.end(), flat,
            [](const QualitySite& s, size_t p) { return s.flat < p; }) - sites.begin());
    };
    auto range = [&](size_t read) {
        const auto& p = placements[read];
        if (p.contig == UINT32_MAX) return std::pair<size_t, size_t>{0, 0};
        const long long begin = std::max(0LL, static_cast<long long>(p.start));
        const long long end = std::min(static_cast<long long>(contigs[p.contig].size()),
                                      static_cast<long long>(p.start) + reads.length(read));
        if (end <= begin) return std::pair<size_t, size_t>{0, 0};
        return std::pair<size_t, size_t>{lower(offset[p.contig] + static_cast<size_t>(begin)),
                                        lower(offset[p.contig] + static_cast<size_t>(end))};
    };
    auto observe = [&](size_t read, const QualitySite& site, quality_consensus::Observation& out) {
        const auto& p = placements[read];
        const long long decoded = static_cast<long long>(site.flat - offset[p.contig]) - p.start;
        const uint32_t raw = static_cast<uint32_t>(p.orientation ? reads.length(read) - 1 - decoded : decoded);
        const int original = reads.originalBaseAt(read, raw);
        if (original < 0) { ++stats.excludedOriginalAmbiguous; return false; }
        const int current = reads.baseAt(read, raw);
        if (current < 0) { ++stats.excludedMasked; return false; }
        if (current != original) { ++stats.excludedModified; return false; }
        const int phred = reads.originalQualityAt(read, raw);
        if (phred < 0) { ++stats.excludedMissingQuality; return false; }
        ++stats.eligibleObservations;
        if (!quality_consensus::logOdds(static_cast<unsigned>(phred))) ++stats.uniformObservations;
        if (phred > static_cast<int>(quality_consensus::kMaxPhred)) ++stats.cappedObservations;
        out = {static_cast<uint8_t>(p.orientation ? 3 - original : original),
               static_cast<uint8_t>(phred), p.orientation};
        return true;
    };

    // One stable traversal of physical fragments. Each mate retains its one
    // already accepted native placement; no extra alignment or sequence vote.
    for (size_t read = 0; read < reads.size();) {
        const unsigned members = reads.hasMate(read) ? 2 : 1;
        auto a = range(read);
        auto b = members == 2 ? range(read + 1) : std::pair<size_t, size_t>{0, 0};
        while (a.first < a.second || b.first < b.second) {
            const size_t ia = a.first < a.second ? a.first : sites.size();
            const size_t ib = b.first < b.second ? b.first : sites.size();
            const size_t index = std::min(ia, ib);
            quality_consensus::Observation observations[2];
            unsigned count = 0;
            if (ia == index) { if (observe(read, sites[index], observations[count])) ++count; ++a.first; }
            if (ib == index) { if (observe(read + 1, sites[index], observations[count])) ++count; ++b.first; }
            if (!count) continue;
            const auto result = quality_consensus::addFragment(sites[index].evidence, observations, count);
            if (result == quality_consensus::AddResult::Informative) ++stats.informativeFragmentSites;
            else if (result == quality_consensus::AddResult::Uniform) ++stats.uniformFragmentSites;
            else if (result == quality_consensus::AddResult::Invalid) sites[index].evidence.overflow = true;
        }
        read += members;
    }
    for (const auto& site : sites) {
        const auto& pile = site.evidence;
        stats.maxFragmentDepth = std::max(stats.maxFragmentDepth, static_cast<size_t>(pile.fragments));
        const auto d = quality_consensus::decide(pile, site.originalContigBase, minDepth, minFraction);
        if (d.best >= 0 && d.best != site.originalContigBase) ++stats.alternateProposals;
        using Reason = quality_consensus::DecisionReason;
        if (d.reason == Reason::BelowDepth) ++stats.belowDepth;
        else if (d.reason == Reason::OneOrientation) ++stats.oneOrientation;
        else if (d.reason == Reason::BelowPosterior) ++stats.belowPosterior;
        else if (d.reason == Reason::Tie) ++stats.tied;
        else if (d.reason == Reason::Overflow) ++stats.overflow;
        if (d.reason != Reason::Accept) continue;
        const size_t c = static_cast<size_t>(std::upper_bound(offset.begin(), offset.end(), site.flat) - offset.begin() - 1);
        const size_t p = site.flat - offset[c];
        const char previous = contigs[c][p];
        contigs[c][p] = codeBase(d.best);
        ++stats.extraChanges;
        std::fprintf(stderr, "[qualitypolish] proposal contig=%zu position=%zu from=%c to=%c fragments=%llu observations=%llu orientations=%u posterior_score=%.12g scores=%llu,%llu,%llu,%llu\n",
                     c, p, previous, contigs[c][p], static_cast<unsigned long long>(pile.fragments),
                     static_cast<unsigned long long>(pile.observations), unsigned(pile.orientations[d.best]), d.posterior,
                     static_cast<unsigned long long>(pile.scores[0]), static_cast<unsigned long long>(pile.scores[1]),
                     static_cast<unsigned long long>(pile.scores[2]), static_cast<unsigned long long>(pile.scores[3]));
    }
}

void logQualitySummary(const QualityPolishStats& s) {
    std::fprintf(stderr, "[qualitypolish] enabled=1 provenance=%d phred_cap=%u score_scale=%.0f original_bytes=%zu placement_bytes=%zu candidate_bytes=%zu candidates=%zu alternate_proposals=%zu extra_changes=%zu informative_fragment_sites=%zu uniform_fragment_sites=%zu max_fragment_depth=%zu eligible_observations=%zu excluded_original_ambiguous=%zu excluded_masked=%zu excluded_net_modified=%zu excluded_missing_quality=%zu uniform_observations=%zu capped_observations=%zu below_depth=%zu one_orientation=%zu below_posterior=%zu ties=%zu overflow=%zu\n",
                 int(s.provenanceAvailable), quality_consensus::kMaxPhred, quality_consensus::kScoreScale,
                 s.originalBytes, s.placementBytes, s.candidateBytes, s.candidatePositions, s.alternateProposals,
                 s.extraChanges, s.informativeFragmentSites, s.uniformFragmentSites, s.maxFragmentDepth,
                 s.eligibleObservations, s.excludedOriginalAmbiguous, s.excludedMasked, s.excludedModified,
                 s.excludedMissingQuality, s.uniformObservations, s.cappedObservations, s.belowDepth,
                 s.oneOrientation, s.belowPosterior, s.tied, s.overflow);
}

}  // namespace

PolishStats polishContigs(std::vector<std::string>& contigs, const SequenceStore& reads,
                          int threads, int anchorK, int minDepth, double minFraction) {
    PolishStats stats;
    stats.quality.enabled = quality_consensus::enabled();
    stats.quality.provenanceAvailable = reads.hasOriginalQualities();
    stats.quality.originalBytes = reads.originalQualityBytes();
    const bool qualityEnabled = stats.quality.enabled && stats.quality.provenanceAvailable;
    if (contigs.empty() || reads.size() == 0) return stats;
    if (threads <= 0) threads = 1;

    const int k = std::max(15, std::min(anchorK, 31));

    // Flatten the contigs so one counter array can cover all of them.
    std::vector<size_t> offset(contigs.size() + 1, 0);
    for (size_t i = 0; i < contigs.size(); ++i) offset[i + 1] = offset[i] + contigs[i].size();
    const size_t totalBases = offset.back();
    if (totalBases == 0) return stats;

    std::unordered_map<Kmer, uint64_t, KmerHasher> index;
    index.reserve(totalBases * 2);
    for (uint32_t c = 0; c < contigs.size(); ++c) {
        const std::string& s = contigs[c];
        if (s.size() < static_cast<size_t>(k)) continue;
        Kmer fwd = 0, rc = 0;
        int valid = 0;
        for (uint32_t p = 0; p < s.size(); ++p) {
            const int b = baseCode(s[p]);
            if (b < 0) { valid = 0; continue; }
            fwd = pushBack(fwd, b, k);
            rc = pushFrontRc(rc, b, k);
            if (++valid < k) continue;
            const Kmer canon = fwd < rc ? fwd : rc;
            auto it = index.find(canon);
            if (it == index.end()) {
                index.emplace(canon, packIndex(c, p + 1 - static_cast<uint32_t>(k),
                                               (fwd == canon) ? 0 : 1));
            } else {
                it->second = kAmbiguous;   // repeated k-mer cannot place a read
            }
        }
    }

    // counts[4 * position + base]; incremented atomically from every worker.
    std::vector<uint32_t> counts(totalBases * 4, 0);
    std::atomic<size_t> used{0};
    std::vector<AcceptedPlacement> placements;
    std::vector<QualitySite> qualitySites;
    if (qualityEnabled) {
        placements.resize(reads.size());
        stats.quality.placementBytes = placements.capacity() * sizeof(AcceptedPlacement);
    }

    auto worker = [&](int tid) {
        size_t localUsed = 0;
        std::string decoded;
        for (size_t r = static_cast<size_t>(tid); r < reads.size(); r += static_cast<size_t>(threads)) {
            const int len = static_cast<int>(reads.length(r));
            if (len < k) continue;

            const int span = len - k;
            const int probes = std::min(kMaxProbes, span + 1);
            struct Vote { uint32_t contig; int32_t pos; uint8_t orient; int count; };
            Vote votes[kMaxProbes];
            int distinct = 0;

            for (int t = 0; t < probes; ++t) {
                const int rp = probes == 1 ? 0 : span * t / (probes - 1);
                Kmer fwd = 0, rcv = 0;
                bool ok = true;
                for (int j = 0; j < k; ++j) {
                    const int b = reads.baseAt(r, static_cast<uint32_t>(rp + j));
                    if (b < 0) { ok = false; break; }
                    fwd = pushBack(fwd, b, k);
                    rcv = pushFrontRc(rcv, b, k);
                }
                if (!ok) continue;
                const Kmer canon = fwd < rcv ? fwd : rcv;
                auto it = index.find(canon);
                if (it == index.end() || it->second == kAmbiguous) continue;

                const uint32_t c = static_cast<uint32_t>(it->second >> 33);
                const int up = static_cast<int>((it->second >> 1) & 0xFFFFFFFFULL);
                const int sflag = static_cast<int>(it->second & 1);
                const int orient = ((fwd == canon) ? 0 : 1) ^ sflag;
                const int32_t startPos = (orient == 0)
                                             ? static_cast<int32_t>(up - rp)
                                             : static_cast<int32_t>(up - (len - rp - k));

                int found = -1;
                for (int q = 0; q < distinct; ++q) {
                    if (votes[q].contig == c && votes[q].pos == startPos &&
                        votes[q].orient == static_cast<uint8_t>(orient)) { found = q; break; }
                }
                if (found >= 0) ++votes[found].count;
                else if (distinct < kMaxProbes) {
                    votes[distinct++] = {c, startPos, static_cast<uint8_t>(orient), 1};
                }
            }

            int bestIdx = -1, bestVotes = 0;
            for (int q = 0; q < distinct; ++q) {
                if (votes[q].count > bestVotes) { bestVotes = votes[q].count; bestIdx = q; }
            }
            if (bestIdx < 0 || bestVotes < kMinVotes) continue;

            const uint32_t c = votes[bestIdx].contig;
            const int32_t pos = votes[bestIdx].pos;
            const int orient = votes[bestIdx].orient;
            const std::string& ctg = contigs[c];

            // Reconstruct the read as it lies on the contig's forward strand.
            decoded.assign(static_cast<size_t>(len), 'N');
            for (int i = 0; i < len; ++i) {
                const int b = reads.baseAt(r, static_cast<uint32_t>(orient == 0 ? i : len - 1 - i));
                if (b < 0) { decoded[static_cast<size_t>(i)] = 'N'; continue; }
                decoded[static_cast<size_t>(i)] = codeBase(orient == 0 ? b : 3 - b);
            }

            int matches = 0, compared = 0;
            for (int i = 0; i < len; ++i) {
                const long long cp = static_cast<long long>(pos) + i;
                if (cp < 0 || cp >= static_cast<long long>(ctg.size())) continue;
                const char rb = decoded[static_cast<size_t>(i)];
                if (rb == 'N') continue;
                ++compared;
                if (rb == ctg[static_cast<size_t>(cp)]) ++matches;
            }
            if (compared == 0) continue;
            if (static_cast<double>(matches) / static_cast<double>(compared) < kMinReadAgreement) continue;

            for (int i = 0; i < len; ++i) {
                const long long cp = static_cast<long long>(pos) + i;
                if (cp < 0 || cp >= static_cast<long long>(ctg.size())) continue;
                const int b = baseCode(decoded[static_cast<size_t>(i)]);
                if (b < 0) continue;
                const size_t slot = (offset[c] + static_cast<size_t>(cp)) * 4 + static_cast<size_t>(b);
                __atomic_fetch_add(&counts[slot], 1u, __ATOMIC_RELAXED);
            }
            if (qualityEnabled) placements[r] = {c, pos, static_cast<uint8_t>(orient)};
            ++localUsed;
        }
        used += localUsed;
    };

    std::vector<std::thread> pool;
    pool.reserve(static_cast<size_t>(threads));
    for (int t = 0; t < threads; ++t) pool.emplace_back(worker, t);
    for (auto& th : pool) th.join();

    size_t changed = 0, covered = 0, lowCov = 0;
    double depthSum = 0;
    for (size_t c = 0; c < contigs.size(); ++c) {
        std::string& s = contigs[c];
        for (size_t p = 0; p < s.size(); ++p) {
            const uint32_t* q = &counts[(offset[c] + p) * 4];
            const uint32_t depth = q[0] + q[1] + q[2] + q[3];
            if (depth == 0) continue;
            ++covered;
            depthSum += depth;
            // A fragment contributes at most one depth unit, so positions
            // below the existing read-depth floor cannot qualify. Cache only
            // mixed positions where the finalized native pile makes no edit.
            if (qualityEnabled && static_cast<int>(depth) >= minDepth) {
                int best = 0;
                for (int b = 1; b < 4; ++b) if (q[b] > q[best]) best = b;
                const double nativeFraction = static_cast<double>(q[best]) / depth;
                const int current = baseCode(s[p]);
                const bool nativeChanges = nativeFraction >= minFraction && s[p] != codeBase(best);
                const bool hasAlternate = current < 0 || q[current] != depth;
                if (!nativeChanges && hasAlternate) qualitySites.push_back({offset[c] + p, current, {}});
            }
            if (static_cast<int>(depth) < minDepth) { ++lowCov; continue; }

            int bi = 0;
            for (int b = 1; b < 4; ++b) if (q[b] > q[bi]) bi = b;
            const double frac = static_cast<double>(q[bi]) / static_cast<double>(depth);
            if (frac < minFraction) continue;
            // A read that carries the true base where the contig is wrong still
            // agrees with the contig everywhere else, so it does anchor and it
            // does vote -- the reason nothing was ever corrected is that the
            // winning fraction has to clear `minFraction`, and at a real error
            // the split is nothing like unanimous.
            const char want = codeBase(bi);
            if (s[p] != want) { s[p] = want; ++changed; }
        }
    }

    if (qualityEnabled) {
        applyOriginalQuality(contigs, reads, offset, placements, qualitySites, minDepth, minFraction, stats.quality);
        changed += stats.quality.extraChanges;
    }
    if (stats.quality.enabled) logQualitySummary(stats.quality);

    stats.readsUsed = used.load();
    stats.basesChanged = changed;
    stats.positionsCovered = covered;
    stats.lowCoveragePositions = lowCov;
    stats.meanDepth = covered ? depthSum / static_cast<double>(covered) : 0;
    return stats;
}

}  // namespace ts
