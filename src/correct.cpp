#include "correct.h"

#include <algorithm>
#include <atomic>
#include <thread>
#include <vector>

namespace ts {

namespace {

// How far past a candidate substitution the k-mers must stay trusted before the
// fix is believed. A real correction restores a long run; a coincidental one
// usually restores only the single k-mer it was chosen to fix.
constexpr int kLookahead = 8;
constexpr int kMinAnchorRun = 2;

struct Fix {
    uint32_t read;
    uint32_t pos;
    uint8_t code;
};

}  // namespace

CorrectionStats correctReads(SequenceStore& reads, const KmerTable& solid, int k, int threads) {
    CorrectionStats stats;
    if (threads <= 0) threads = 1;
    if (solid.size() == 0) return stats;

    std::vector<std::vector<Fix>> perThread(static_cast<size_t>(threads));
    std::atomic<size_t> examined{0}, corrected{0}, uncorrectable{0};

    auto worker = [&](int tid) {
        std::vector<Fix>& fixes = perThread[static_cast<size_t>(tid)];
        std::vector<int> codes;
        std::vector<Kmer> fwd;
        std::vector<char> ok;
        size_t localExamined = 0, localCorrected = 0, localUncorrectable = 0;

        for (size_t r = static_cast<size_t>(tid); r < reads.size(); r += static_cast<size_t>(threads)) {
            const int len = static_cast<int>(reads.length(r));
            if (len < k + 1) continue;
            ++localExamined;

            codes.assign(static_cast<size_t>(len), 0);
            bool ambiguous = false;
            for (int i = 0; i < len; ++i) {
                const int c = reads.baseAt(r, static_cast<uint32_t>(i));
                if (c < 0) { ambiguous = true; break; }
                codes[static_cast<size_t>(i)] = c;
            }
            if (ambiguous) continue;

            const int m = len - k + 1;
            fwd.assign(static_cast<size_t>(m), Kmer());
            ok.assign(static_cast<size_t>(m), 0);
            {
                Kmer f = 0;
                for (int i = 0; i < len; ++i) {
                    f = pushBack(f, codes[static_cast<size_t>(i)], k);
                    if (i >= k - 1) {
                        const int p = i - k + 1;
                        fwd[static_cast<size_t>(p)] = f;
                        ok[static_cast<size_t>(p)] = solid.contains(canonical(f, k)) ? 1 : 0;
                    }
                }
            }

            // Longest trusted run anchors the correction; without one there is
            // nothing to extend from and the read is left untouched.
            int bestLo = -1, bestLen = 0, curLo = -1, curLen = 0;
            for (int i = 0; i < m; ++i) {
                if (ok[static_cast<size_t>(i)]) {
                    if (curLo < 0) { curLo = i; curLen = 0; }
                    ++curLen;
                    if (curLen > bestLen) { bestLen = curLen; bestLo = curLo; }
                } else {
                    curLo = -1;
                    curLen = 0;
                }
            }
            if (bestLo < 0 || bestLen < kMinAnchorRun) { ++localUncorrectable; continue; }

            const int hi = bestLo + bestLen - 1;
            const int maxFixes = std::max(2, len / 10);
            int applied = 0;

            // Extend right: the only base the next k-mer adds is at pos+k.
            Kmer cur = fwd[static_cast<size_t>(hi)];
            for (int pos = hi; pos + 1 < m && applied < maxFixes; ++pos) {
                const int idx = pos + k;
                const int obs = codes[static_cast<size_t>(idx)];
                Kmer cand = pushBack(cur, obs, k);
                if (solid.contains(canonical(cand, k))) { cur = cand; continue; }

                int bestBase = -1, bestScore = 0;
                for (int b = 0; b < 4; ++b) {
                    if (b == obs) continue;
                    Kmer t = pushBack(cur, b, k);
                    if (!solid.contains(canonical(t, k))) continue;
                    int score = 1;
                    Kmer tmp = t;
                    for (int j = 1; j <= kLookahead && pos + 1 + j < m; ++j) {
                        tmp = pushBack(tmp, codes[static_cast<size_t>(idx + j)], k);
                        if (!solid.contains(canonical(tmp, k))) break;
                        ++score;
                    }
                    if (score > bestScore) { bestScore = score; bestBase = b; }
                }
                if (bestBase < 0) break;
                codes[static_cast<size_t>(idx)] = bestBase;
                fixes.push_back({static_cast<uint32_t>(r), static_cast<uint32_t>(idx),
                                 static_cast<uint8_t>(bestBase)});
                ++applied;
                cur = pushBack(cur, bestBase, k);
            }

            // Extend left: the k-mer one position earlier adds the base at lo-1.
            cur = fwd[static_cast<size_t>(bestLo)];
            for (int pos = bestLo; pos > 0 && applied < maxFixes; --pos) {
                const int idx = pos - 1;
                const int obs = codes[static_cast<size_t>(idx)];
                Kmer cand = pushFront(cur, obs, k);
                if (solid.contains(canonical(cand, k))) { cur = cand; continue; }

                int bestBase = -1, bestScore = 0;
                for (int b = 0; b < 4; ++b) {
                    if (b == obs) continue;
                    Kmer t = pushFront(cur, b, k);
                    if (!solid.contains(canonical(t, k))) continue;
                    int score = 1;
                    Kmer tmp = t;
                    for (int j = 1; j <= kLookahead && idx - j >= 0; ++j) {
                        tmp = pushFront(tmp, codes[static_cast<size_t>(idx - j)], k);
                        if (!solid.contains(canonical(tmp, k))) break;
                        ++score;
                    }
                    if (score > bestScore) { bestScore = score; bestBase = b; }
                }
                if (bestBase < 0) break;
                codes[static_cast<size_t>(idx)] = bestBase;
                fixes.push_back({static_cast<uint32_t>(r), static_cast<uint32_t>(idx),
                                 static_cast<uint8_t>(bestBase)});
                ++applied;
                cur = pushFront(cur, bestBase, k);
            }

            if (applied) ++localCorrected;
        }
        examined += localExamined;
        corrected += localCorrected;
        uncorrectable += localUncorrectable;
    };

    std::vector<std::thread> pool;
    pool.reserve(static_cast<size_t>(threads));
    for (int t = 0; t < threads; ++t) pool.emplace_back(worker, t);
    for (auto& th : pool) th.join();

    // Applied serially: neighbouring reads share 64-bit words in the packed
    // store, so concurrent writes would race on the boundary words.
    for (const auto& v : perThread) {
        for (const Fix& f : v) {
            reads.setBase(f.read, f.pos, f.code);
            ++stats.basesCorrected;
        }
    }

    stats.readsExamined = examined.load();
    stats.readsCorrected = corrected.load();
    stats.readsUncorrectable = uncorrectable.load();
    return stats;
}

}  // namespace ts
