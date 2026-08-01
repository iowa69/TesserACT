#include "polish.h"

#include <algorithm>
#include <atomic>
#include <thread>
#include <unordered_map>
#include <vector>

#include "kmer.h"

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

}  // namespace

PolishStats polishContigs(std::vector<std::string>& contigs, const SequenceStore& reads,
                          int threads, int anchorK, int minDepth, double minFraction) {
    PolishStats stats;
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

    stats.readsUsed = used.load();
    stats.basesChanged = changed;
    stats.positionsCovered = covered;
    stats.lowCoveragePositions = lowCov;
    stats.meanDepth = covered ? depthSum / static_cast<double>(covered) : 0;
    return stats;
}

}  // namespace ts
