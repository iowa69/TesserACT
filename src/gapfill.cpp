#include "gapfill.h"

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "kmer.h"
#include "util.h"

namespace ts {

namespace {

// A gap needs at least this much sequence on both sides to be worked on: the
// flank is what recruits the reads and anchors the walk.
constexpr uint32_t kMinFlank = 40;

// Search limits. A gap is a small, local problem -- if it does not fall out
// within this budget the answer would not be trustworthy anyway.
constexpr int kMaxExpansions = 200000;
constexpr size_t kMaxSolutions = 24;  // enough to tell "one clear way" from "many"
// How far past the estimated gap to look. The estimate comes from a fragment
// model whose spread is tens of bases and whose mean is pulled around by the
// chimeric pairs in any real library, so it is a hint about the order of
// magnitude, not a length.
constexpr int kSlackBases = 900;
constexpr uint16_t kCountCeil = 60000;

// A gap this long was estimated from the fragment model with little to go on;
// closing it from 150 bp reads is not something to trust.
constexpr uint32_t kMaxGapLen = 3000;

// A flank that pulls in more reads than a few hundred x of local depth is not
// a flank, it is a repeat.
constexpr size_t kMaxRecruits = 40000;

// How far the best path must outrun the runner-up on its weakest step
// before it is written as sequence rather than left as Ns.
constexpr double kDominanceRatio = 3.0;

constexpr uint32_t kAmbiguousGap = UINT32_MAX;

struct Gap {
    uint32_t contig = 0;
    uint32_t start = 0;    // first N
    uint32_t len = 0;      // number of Ns
};

// Encoded 2-bit sequence of one recruited read, so the per-gap reassembly does
// not have to go back through the store's ambiguity bitmap for every k-mer.
struct Local {
    std::unordered_map<Kmer, uint16_t, KmerHasher> counts;

    uint16_t at(const Kmer& km, int k) const {
        auto it = counts.find(canonical(km, k));
        return it == counts.end() ? 0 : it->second;
    }
};

// One way through the gap, with the evidence behind its weakest step.
struct Solution {
    std::string seq;
    uint16_t weakest = 0;   // smallest k-mer count anywhere along the path
};

// Depth-first walk from the left anchor towards the right one, appending one
// base at a time. Every distinct way through is collected rather than the
// first, because a gap almost always sits where the graph was hard, and which
// of two candidate paths is real is a question about their coverage -- not
// about which the search happened to reach first.
struct Walker {
    const Local& local;
    int k;
    Kmer target;
    uint32_t minLen;
    uint32_t maxLen;
    uint16_t minCount;
    int expansions = 0;
    std::vector<Solution> solutions;

    void walk(Kmer cur, std::string& out, uint16_t weakest) {
        if (solutions.size() >= kMaxSolutions || expansions > kMaxExpansions) return;
        ++expansions;

        // Candidate next bases, strongest first: a real continuation is
        // normally the deepest one, so the best answer is found early and the
        // rest of the search only has to weigh what else could fit.
        struct Cand { int base; uint16_t count; };
        Cand cands[4];
        int nc = 0;
        for (int b = 0; b < 4; ++b) {
            const Kmer nxt = pushBack(cur, b, k);
            const uint16_t c = local.at(nxt, k);
            if (c >= minCount) cands[nc++] = {b, c};
        }
        // Insertion sort by hand: with at most four candidates std::sort's
        // threshold machinery costs more than the sort, and it makes GCC's
        // bounds analysis complain about a fixed array it cannot overrun.
        for (int i = 1; i < nc; ++i) {
            const Cand key = cands[i];
            int j = i - 1;
            while (j >= 0 && cands[j].count < key.count) { cands[j + 1] = cands[j]; --j; }
            cands[j + 1] = key;
        }

        for (int i = 0; i < nc; ++i) {
            const Kmer nxt = pushBack(cur, cands[i].base, k);
            const uint16_t w = std::min(weakest, cands[i].count);
            out.push_back(codeBase(cands[i].base));
            // `minLen` is k: the walk has to have written the whole target
            // k-mer itself before the prefix in front of it is the gap. A
            // shorter hit would mean the target overlaps the left flank, which
            // is a claim about the join, not about the gap, and is not this
            // stage's to make.
            if (nxt == target && out.size() >= minLen) {
                // The target k-mer is the first k-mer of the right flank, so
                // the bases that fill the gap are everything before it.
                solutions.push_back({out.substr(0, out.size() - static_cast<size_t>(k)), w});
                out.pop_back();
                if (solutions.size() >= kMaxSolutions) return;
                continue;
            }
            if (out.size() < maxLen) walk(nxt, out, w);
            out.pop_back();
            if (solutions.size() >= kMaxSolutions || expansions > kMaxExpansions) return;
        }
    }
};

}  // namespace

GapFillStats closeGaps(std::vector<std::string>& contigs, const SequenceStore& reads,
                       int threads, int k, int flank) {
    GapFillStats stats;
    util::Timer timer;
    if (contigs.empty() || reads.size() == 0) return stats;
    if (threads <= 0) threads = 1;
    k = std::max(19, std::min(k, 31));
    if (flank < static_cast<int>(kMinFlank)) flank = static_cast<int>(kMinFlank);

    // ---- 1. locate the gaps ---------------------------------------------
    std::vector<Gap> gaps;
    for (uint32_t c = 0; c < contigs.size(); ++c) {
        const std::string& s = contigs[c];
        uint32_t i = 0;
        while (i < s.size()) {
            if (s[i] != 'N') { ++i; continue; }
            uint32_t j = i;
            while (j < s.size() && s[j] == 'N') ++j;
            ++stats.gapsSeen;
            const uint32_t len = j - i;
            const bool roomLeft = i >= std::max<uint32_t>(kMinFlank, static_cast<uint32_t>(k));
            const bool roomRight = s.size() - j >= std::max<uint32_t>(kMinFlank, static_cast<uint32_t>(k));
            if (len <= kMaxGapLen && roomLeft && roomRight) gaps.push_back({c, i, len});
            i = j;
        }
    }
    if (gaps.empty()) {
        stats.seconds = timer.elapsed();
        return stats;
    }

    // ---- 2. index the flanks --------------------------------------------
    // A k-mer shared by two gaps identifies neither, so it is dropped rather
    // than dragging one gap's reads into the other's reassembly.
    std::unordered_map<Kmer, uint32_t, KmerHasher> flankIndex;
    flankIndex.reserve(gaps.size() * static_cast<size_t>(flank) * 4);
    auto indexRange = [&](const std::string& s, size_t from, size_t to, uint32_t gid) {
        if (to <= from || to - from < static_cast<size_t>(k)) return;
        Kmer fwd = 0, rc = 0;
        int valid = 0;
        for (size_t p = from; p < to; ++p) {
            const int b = baseCode(s[p]);
            if (b < 0) { valid = 0; continue; }
            fwd = pushBack(fwd, b, k);
            rc = pushFrontRc(rc, b, k);
            if (++valid < k) continue;
            const Kmer canon = fwd < rc ? fwd : rc;
            auto it = flankIndex.find(canon);
            if (it == flankIndex.end()) flankIndex.emplace(canon, gid);
            else if (it->second != gid) it->second = kAmbiguousGap;
        }
    };
    for (uint32_t g = 0; g < gaps.size(); ++g) {
        const std::string& s = contigs[gaps[g].contig];
        const size_t ls = gaps[g].start > static_cast<uint32_t>(flank)
                              ? gaps[g].start - static_cast<uint32_t>(flank) : 0;
        indexRange(s, ls, gaps[g].start, g);
        const size_t re = std::min(s.size(), static_cast<size_t>(gaps[g].start + gaps[g].len + flank));
        indexRange(s, gaps[g].start + gaps[g].len, re, g);
    }

    // ---- 3. recruit reads ------------------------------------------------
    // A read is recruited by its own k-mers or by its mate's: the reads that
    // actually sit *inside* the gap share no k-mer with either flank, and the
    // mate is the only thing that can place them.
    std::vector<std::vector<uint32_t>> bucket(gaps.size());
    {
        std::vector<std::vector<std::vector<uint32_t>>> local(
            static_cast<size_t>(threads), std::vector<std::vector<uint32_t>>(gaps.size()));
        const bool paired = reads.paired();
        auto worker = [&](int tid) {
            auto& mine = local[static_cast<size_t>(tid)];
            for (size_t r = static_cast<size_t>(tid); r < reads.size();
                 r += static_cast<size_t>(threads)) {
                uint32_t hit = kAmbiguousGap;
                bool conflict = false;
                forEachKmerRaw(reads, r, k, [&](const Kmer& km, uint32_t) {
                    if (conflict) return;
                    auto it = flankIndex.find(km);
                    if (it == flankIndex.end() || it->second == kAmbiguousGap) return;
                    if (hit == kAmbiguousGap) hit = it->second;
                    else if (hit != it->second) conflict = true;
                });
                if (conflict || hit == kAmbiguousGap) continue;
                mine[hit].push_back(static_cast<uint32_t>(r));
                if (paired) mine[hit].push_back(static_cast<uint32_t>(reads.mateOf(r)));
            }
        };
        std::vector<std::thread> pool;
        pool.reserve(static_cast<size_t>(threads));
        for (int t = 0; t < threads; ++t) pool.emplace_back(worker, t);
        for (auto& th : pool) th.join();

        for (uint32_t g = 0; g < gaps.size(); ++g) {
            size_t total = 0;
            for (int t = 0; t < threads; ++t) total += local[static_cast<size_t>(t)][g].size();
            bucket[g].reserve(total);
            for (int t = 0; t < threads; ++t) {
                auto& v = local[static_cast<size_t>(t)][g];
                bucket[g].insert(bucket[g].end(), v.begin(), v.end());
                std::vector<uint32_t>().swap(v);
            }
            std::sort(bucket[g].begin(), bucket[g].end());
            bucket[g].erase(std::unique(bucket[g].begin(), bucket[g].end()), bucket[g].end());
            stats.readsRecruited += bucket[g].size();
        }
    }

    // ---- 4. one local reassembly per gap ---------------------------------
    std::vector<std::string> fill(gaps.size());
    std::vector<uint8_t> closed(gaps.size(), 0);
    std::atomic<size_t> cursor{0};
    std::atomic<size_t> nAmbiguous{0}, nNoPath{0}, nThinPool{0}, nBudget{0};
    std::atomic<size_t> nSeedGone{0}, nTargetGone{0};
    std::atomic<size_t> dbgDepth{0}, dbgFloor{0};
    const bool debug = std::getenv("TESSERA_GF_DEBUG") != nullptr;

    auto solve = [&](int) {
        Local local;
        std::string path;
        for (;;) {
            const size_t g = cursor.fetch_add(1);
            if (g >= gaps.size()) break;
            const Gap& gp = gaps[g];
            const std::string& s = contigs[gp.contig];
            // Too few reads is no evidence; absurdly many means the flank is
            // itself repetitive and the pool is a mixture of loci, where a
            // "unique" path would be an artefact of the abundance floor.
            if (bucket[g].size() < 4 || bucket[g].size() > kMaxRecruits) { ++nThinPool; continue; }

            local.counts.clear();
            uint64_t recruitedBases = 0;
            for (uint32_t r : bucket[g]) {
                recruitedBases += reads.length(r);
                forEachKmerRaw(reads, r, k, [&](const Kmer& km, uint32_t) {
                    uint16_t& c = local.counts[km];
                    if (c < kCountCeil) ++c;
                });
            }

            // Abundance floor for the local graph. The reads here cover roughly
            // 2*flank + gap bases, so their depth is knowable, and an error
            // k-mer sits far below it. Reads from a repeat copy elsewhere in
            // the genome also land in this pool, which is exactly why the floor
            // must scale rather than sit at a constant 2.
            const double span = static_cast<double>(2 * flank + gp.len);
            const double depth = span > 0 ? static_cast<double>(recruitedBases) / span : 0;
            const uint16_t minCount = static_cast<uint16_t>(
                std::max(2.0, std::min(12.0, depth * 0.12)));

            bool ok = false;
            const Kmer seed = stringToKmer(
                s.substr(gp.start - static_cast<size_t>(k), static_cast<size_t>(k)), k, ok);
            if (!ok) { ++nNoPath; continue; }
            const Kmer target = stringToKmer(
                s.substr(gp.start + gp.len, static_cast<size_t>(k)), k, ok);
            if (!ok) { ++nNoPath; continue; }

            // An anchor missing from the pool its own flank recruited means the
            // recruitment or the floor is wrong, not the data; worth counting
            // separately from a gap nothing spans.
            if (local.at(seed, k) < minCount) ++nSeedGone;
            if (local.at(target, k) < minCount) ++nTargetGone;
            dbgDepth += static_cast<size_t>(depth);
            dbgFloor += minCount;

            if (debug && g < 12) {
                char buf[256];
                int off = std::snprintf(buf, sizeof(buf),
                    "      [gf] gap%u len=%u reads=%zu depth=%.0f floor=%u seedNext=",
                    static_cast<unsigned>(g), gp.len, bucket[g].size(), depth,
                    static_cast<unsigned>(minCount));
                for (int b = 0; b < 4; ++b)
                    off += std::snprintf(buf + off, sizeof(buf) - off, "%u,",
                                         static_cast<unsigned>(local.at(pushBack(seed, b, k), k)));
                std::snprintf(buf + off, sizeof(buf) - off, " seedCnt=%u tgtCnt=%u",
                              static_cast<unsigned>(local.at(seed, k)),
                              static_cast<unsigned>(local.at(target, k)));
                std::fprintf(stderr, "%s\n", buf);
            }

            // A gap exists because coverage there is poor -- that is usually
            // why the graph stopped in the first place. Starting at the floor
            // the flanks justify and stepping down only when nothing spans the
            // gap keeps the strict answer when there is one, and still reaches
            // the sequence sitting under a dropout. Stepping down after an
            // *ambiguous* result would be pointless: a lower floor can only add
            // paths, never remove them.
            const uint16_t floors[3] = {minCount, static_cast<uint16_t>(3), static_cast<uint16_t>(2)};
            size_t solutions = 0;
            int lastExpansions = 0;
            for (int f = 0; f < 3; ++f) {
                if (f > 0 && floors[f] >= floors[f - 1]) continue;
                // The estimated gap is only an estimate -- the true distance
                // can be shorter (the model over-shot) or longer.
                Walker w{local, k, target, static_cast<uint32_t>(k),
                         gp.len + static_cast<uint32_t>(kSlackBases) + static_cast<uint32_t>(k),
                         floors[f], 0, {}};
                path.clear();
                w.walk(seed, path, kCountCeil);
                solutions = w.solutions.size();
                lastExpansions = w.expansions;
                if (solutions == 0) continue;   // nothing spanned it; try a lower floor

                // Rank by the weakest link: a path that never drops below 40x
                // is read-backed along its whole length, while one that dips to
                // 3x is a chain of coincidences that happens to end in the
                // right place. Ties on that go to the length the fragment model
                // predicted.
                size_t bestI = 0;
                for (size_t i = 1; i < w.solutions.size(); ++i) {
                    const Solution& a = w.solutions[i];
                    const Solution& b = w.solutions[bestI];
                    const long da = std::labs(static_cast<long>(a.seq.size()) - static_cast<long>(gp.len));
                    const long db = std::labs(static_cast<long>(b.seq.size()) - static_cast<long>(gp.len));
                    if (a.weakest > b.weakest || (a.weakest == b.weakest && da < db)) bestI = i;
                }
                uint16_t runnerUp = 0;
                for (size_t i = 0; i < w.solutions.size(); ++i) {
                    if (i != bestI && w.solutions[i].weakest > runnerUp)
                        runnerUp = w.solutions[i].weakest;
                }
                // One path, or one path that dominates everything else by a
                // clear margin. Anything closer than that is a repeat the reads
                // cannot separate, and guessing there is how a gap becomes a
                // misassembly.
                if (w.solutions.size() == 1 ||
                    static_cast<double>(w.solutions[bestI].weakest) >=
                        kDominanceRatio * static_cast<double>(runnerUp)) {
                    fill[g] = w.solutions[bestI].seq;
                    closed[g] = 1;
                }
                break;
            }

            if (!closed[g]) {
                if (solutions > 0) ++nAmbiguous;
                else {
                    ++nNoPath;
                    if (lastExpansions > kMaxExpansions) ++nBudget;
                }
            }
        }
    };
    {
        std::vector<std::thread> pool;
        pool.reserve(static_cast<size_t>(threads));
        for (int t = 0; t < threads; ++t) pool.emplace_back(solve, t);
        for (auto& th : pool) th.join();
    }
    stats.gapsAmbiguous = nAmbiguous.load();
    stats.gapsNoPath = nNoPath.load() + nThinPool.load();
    stats.gapsThinPool = nThinPool.load();
    stats.gapsOutOfBudget = nBudget.load();
    stats.seedBelowFloor = nSeedGone.load();
    stats.targetBelowFloor = nTargetGone.load();
    const size_t nWorked = gaps.size() - nThinPool.load();
    stats.meanLocalDepth = nWorked ? static_cast<double>(dbgDepth.load()) / static_cast<double>(nWorked) : 0;
    stats.meanFloor = nWorked ? static_cast<double>(dbgFloor.load()) / static_cast<double>(nWorked) : 0;

    // ---- 5. splice -------------------------------------------------------
    // Right to left, so an edit never moves a gap that has not been done yet.
    for (size_t gi = gaps.size(); gi-- > 0;) {
        if (!closed[gi]) continue;
        Gap& gp = gaps[gi];
        contigs[gp.contig].replace(gp.start, gp.len, fill[gi]);
        ++stats.gapsClosed;
        stats.nBasesRemoved += gp.len;
        stats.basesInserted += fill[gi].size();
    }

    stats.seconds = timer.elapsed();
    return stats;
}

}  // namespace ts
