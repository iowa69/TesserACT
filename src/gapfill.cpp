#include "gapfill.h"
#include "envflags.h"

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

#include "emit_fixflags.h"
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

// T11 back-off limits: moved anchors tried per side, and the expansion budget all back-off
// ladders of one gap may spend together (the release ladder alone may spend 3 budgets).
constexpr size_t kBackoffCandidatesPerSide = 8;
constexpr long long kBackoffExpansionBudget = 8LL * kMaxExpansions;

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
    // T09: set when the expansion budget stopped the search while some branch was still
    // unexplored. The walk returns at exactly the release points; this only records it.
    bool truncated = false;

    void walk(Kmer cur, std::string& out, uint16_t weakest) {
        if (solutions.size() >= kMaxSolutions) return;
        if (expansions > kMaxExpansions) { truncated = true; return; }
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
            if (solutions.size() >= kMaxSolutions) return;
            if (expansions > kMaxExpansions) {
                if (i + 1 < nc) truncated = true;
                return;
            }
        }
    }
};

// Outcome of one floor ladder between one pair of anchors.
struct LadderResult {
    bool closed = false;
    size_t solutions = 0;       // of the deciding search (the last one run)
    int expansions = 0;         // of the deciding search
    long long totalExpansions = 0;
    bool truncated = false;     // deciding search stopped at the budget with work left
    bool capped = false;        // deciding search stopped at kMaxSolutions
    bool refusedTruncated = false;  // would have closed; withheld by the strict flag
};

// A closure after back-off replaces flank bases as well as Ns; two such closures a few
// dozen bases apart could claim the same bases.
struct Replacement {
    uint32_t contig = 0;
    uint32_t start = 0;
    uint32_t len = 0;
    bool backedOff = false;
    bool closed = false;
};

// Walking each contig left to right, a closure's replaced stretch must end at or before
// the next closure's seed k-mer starts -- which also keeps the next stretch clear of this
// closure's target k-mer, the two conditions being the same inequality. Release closures
// always pass (their anchors hold no N, so a seed never starts inside the Ns before it);
// where a clash involves a back-off closure, the back-off closure is dropped.
size_t cancelBackoffClashes(std::vector<Replacement>& reps, int k) {
    size_t cancelled = 0;
    size_t prev = SIZE_MAX;
    for (size_t gi = 0; gi < reps.size(); ++gi) {
        Replacement& c = reps[gi];
        if (!c.closed) continue;
        if (prev != SIZE_MAX && reps[prev].contig == c.contig &&
            static_cast<uint64_t>(c.start) <
                static_cast<uint64_t>(reps[prev].start) + reps[prev].len + static_cast<uint64_t>(k)) {
            if (c.backedOff) { c.closed = false; ++cancelled; continue; }
            if (reps[prev].backedOff) { reps[prev].closed = false; ++cancelled; }
        }
        prev = gi;
    }
    return cancelled;
}

}  // namespace

size_t cancelBackoffClashes(std::vector<GapClosure>& closures, int k) {
    std::vector<Replacement> reps(closures.size());
    for (size_t i = 0; i < closures.size(); ++i) {
        reps[i].contig = closures[i].contig;
        reps[i].start = closures[i].repStart;
        reps[i].len = closures[i].repLen;
        reps[i].backedOff = closures[i].backedOff;
        reps[i].closed = closures[i].closed;
    }
    const size_t cancelled = cancelBackoffClashes(reps, k);
    for (size_t i = 0; i < closures.size(); ++i) closures[i].closed = reps[i].closed;
    return cancelled;
}

std::vector<uint8_t> readsWithInputAmbiguity(const SequenceStore& reads, int threads) {
    std::vector<uint8_t> out(reads.size(), 0);
    if (threads <= 0) threads = 1;
    auto worker = [&](int tid) {
        for (size_t r = static_cast<size_t>(tid); r < reads.size(); r += static_cast<size_t>(threads)) {
            const uint32_t len = reads.length(r);
            for (uint32_t p = 0; p < len; ++p) {
                if (reads.baseAt(r, p) < 0) { out[r] = 1; break; }
            }
        }
    };
    std::vector<std::thread> pool;
    pool.reserve(static_cast<size_t>(threads));
    for (int t = 0; t < threads; ++t) pool.emplace_back(worker, t);
    for (auto& th : pool) th.join();
    return out;
}

GapFillStats gapFillIdleStats() {
    GapFillStats s;
    s.strictBudget = emitfix::enabled(emitfix::kGapfillStrict);
    s.backoff = emitfix::backoffBases();
    s.skipInputN = emitfix::enabled(emitfix::kGapfillSkipInputN);
    return s;
}

void logGapFillCounters(const GapFillStats& s) {
    std::fprintf(stderr,
                 "[gapfill-budget] strict=%d truncated=%zu truncAccepted=%zu truncRefused=%zu "
                 "capped=%zu cappedAccepted=%zu searches=%zu\n",
                 s.strictBudget ? 1 : 0, s.gapsTruncated, s.gapsTruncatedAccepted,
                 s.gapsTruncatedRefused, s.gapsCapped, s.gapsCappedAccepted, s.searches);
    std::fprintf(stderr,
                 "[gapfill-backoff] B=%d seedBelowFloor=%zu targetBelowFloor=%zu seedBackedOff=%zu "
                 "targetBackedOff=%zu attempts=%zu closedAfterBackoff=%zu trimmedBp=%zu cancelled=%zu\n",
                 s.backoff, s.seedBelowFloor, s.targetBelowFloor, s.seedBackedOff, s.targetBackedOff,
                 s.backoffAttempts, s.closedAfterBackoff, s.trimmedBp, s.backoffCancelled);
    std::fprintf(stderr,
                 "[gapfill-inputn] enabled=%d provenance=%d reads_input_n=%zu recruited_input_n=%zu "
                 "kmers_skipped_input_n=%zu\n",
                 s.skipInputN ? 1 : 0, s.inputNProvenance ? 1 : 0, s.readsInputN, s.recruitedInputN,
                 s.kmersSkippedInputN);
    if (s.skipInputN && !s.inputNProvenance) {
        std::fprintf(stderr, "[gapfill-inputn] WARNING: %s is on but no input-N provenance was "
                             "supplied; reads were read as in release\n", emitfix::kGapfillSkipInputN);
    }
}

GapFillStats closeGaps(std::vector<std::string>& contigs, const SequenceStore& reads,
                       int threads, int k, int flank,
                       const std::vector<uint8_t>* inputAmbiguousReads) {
    return closeGaps(contigs, reads, threads, k, flank, inputAmbiguousReads, nullptr, nullptr);
}

GapFillStats closeGaps(std::vector<std::string>& contigs, const SequenceStore& reads,
                       int threads, int k, int flank,
                       const std::vector<uint8_t>* inputAmbiguousReads,
                       const std::vector<uint8_t>* allow, size_t* blocked) {
    if (blocked) *blocked = 0;
    GapFillStats stats;
    util::Timer timer;
    // Fix switches, read per call (never cached): see emit_fixflags.h.
    const bool strictBudget = emitfix::enabled(emitfix::kGapfillStrict);
    const int backoff = emitfix::backoffBases();
    const bool skipInputN = emitfix::enabled(emitfix::kGapfillSkipInputN);
    stats.strictBudget = strictBudget;
    stats.backoff = backoff;
    stats.skipInputN = skipInputN;
    const bool provenance = inputAmbiguousReads != nullptr && inputAmbiguousReads->size() == reads.size();
    stats.inputNProvenance = provenance;
    // T27: only a read that carried a non-ACGT base in the input is read with its ambiguous
    // positions skipped. Every other read keeps the release raw read-through of correction
    // masks, which is what this stage relies on at coverage dropouts.
    const bool useInputN = skipInputN && provenance;
    if (useInputN) {
        for (uint8_t v : *inputAmbiguousReads) stats.readsInputN += v ? 1 : 0;
    }
    auto hadInputN = [&](size_t r) { return useInputN && (*inputAmbiguousReads)[r] != 0; };

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
            // Organism Model 2.0 C1e: a gap the junction ledger does not allow is left as N.
            if (allow && stats.gapsSeen - 1 < allow->size() && !(*allow)[stats.gapsSeen - 1]) {
                if (blocked) ++*blocked;
                i = j;
                continue;
            }
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
                auto probe = [&](const Kmer& km, uint32_t) {
                    if (conflict) return;
                    auto it = flankIndex.find(km);
                    if (it == flankIndex.end() || it->second == kAmbiguousGap) return;
                    if (hit == kAmbiguousGap) hit = it->second;
                    else if (hit != it->second) conflict = true;
                };
                if (hadInputN(r)) forEachKmer(reads, r, k, probe);
                else forEachKmerRaw(reads, r, k, probe);
                if (conflict || hit == kAmbiguousGap) continue;
                mine[hit].push_back(static_cast<uint32_t>(r));
                if (paired && reads.hasMate(r)) {
                    mine[hit].push_back(static_cast<uint32_t>(reads.mateOf(r)));
                }
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
            if (useInputN) {
                for (uint32_t r : bucket[g]) stats.recruitedInputN += (*inputAmbiguousReads)[r] ? 1 : 0;
            }
        }
    }

    // ---- 4. one local reassembly per gap ---------------------------------
    std::vector<std::string> fill(gaps.size());
    std::vector<uint8_t> closed(gaps.size(), 0);
    // What a closure replaces: the Ns, or -- after a back-off -- the Ns plus the flank tip
    // bases between the moved anchors. Read only where closed[g].
    std::vector<uint32_t> repStart(gaps.size(), 0), repLen(gaps.size(), 0);
    // 0 = release anchors; else the release attempt's outcome, kept so a cancelled
    // back-off closure is counted as the release would count it: 1 ambiguous,
    // 2 unspanned, 3 unspanned and out of budget.
    std::vector<uint8_t> backedOff(gaps.size(), 0);
    std::atomic<size_t> cursor{0};
    std::atomic<size_t> nAmbiguous{0}, nNoPath{0}, nThinPool{0}, nBudget{0};
    std::atomic<size_t> nSeedGone{0}, nTargetGone{0};
    std::atomic<size_t> dbgDepth{0}, dbgFloor{0};
    std::atomic<size_t> nTruncated{0}, nTruncAccepted{0}, nTruncRefused{0}, nCapped{0}, nCappedAccepted{0};
    std::atomic<size_t> nSeedBackedOff{0}, nTargetBackedOff{0}, nAttempts{0};
    std::atomic<size_t> nKmersSkippedInputN{0}, nSearches{0};
    const bool debug = env::present("TESSERACT_GF_DEBUG");

    auto solve = [&](int) {
        Local local;
        std::string path;
        // The floor ladder from `seed` to `target`. `expect` is the length the replaced
        // stretch is predicted to have (the tie-break) and `maxLen` the walk ceiling.
        // Writes `fillOut` only when it closes.
        auto ladder = [&](const Kmer& seed, const Kmer& target, uint32_t expect, uint32_t maxLen,
                          uint16_t minCount, std::string& fillOut) -> LadderResult {
            LadderResult res;
            // A gap exists because coverage there is poor -- that is usually
            // why the graph stopped in the first place. Starting at the floor
            // the flanks justify and stepping down only when nothing spans the
            // gap keeps the strict answer when there is one, and still reaches
            // the sequence sitting under a dropout. Stepping down after an
            // *ambiguous* result would be pointless: a lower floor can only add
            // paths, never remove them.
            const uint16_t floors[3] = {minCount, static_cast<uint16_t>(3), static_cast<uint16_t>(2)};
            // T38: each floor runs at most once, and only below every floor already tried.
            // The release test compared a floor with the array entry before it only, so with
            // minCount == 2 the floor-2 walk ran twice; a walk depends only on its inputs, so
            // the repeat returned the same nothing at up to a full budget of expansions.
            uint16_t lowestTried = UINT16_MAX;
            for (int f = 0; f < 3; ++f) {
                if (floors[f] >= lowestTried) continue;
                lowestTried = floors[f];
                // The estimated gap is only an estimate -- the true distance
                // can be shorter (the model over-shot) or longer.
                Walker w{local, k, target, static_cast<uint32_t>(k), maxLen, floors[f], 0, {}};
                path.clear();
                w.walk(seed, path, kCountCeil);
                ++nSearches;
                res.solutions = w.solutions.size();
                res.expansions = w.expansions;
                res.totalExpansions += w.expansions;
                res.truncated = w.truncated;
                res.capped = w.solutions.size() >= kMaxSolutions;
                if (res.solutions == 0) continue;   // nothing spanned it; try a lower floor

                // Rank by the weakest link: a path that never drops below 40x
                // is read-backed along its whole length, while one that dips to
                // 3x is a chain of coincidences that happens to end in the
                // right place. Ties on that go to the length the fragment model
                // predicted.
                size_t bestI = 0;
                for (size_t i = 1; i < w.solutions.size(); ++i) {
                    const Solution& a = w.solutions[i];
                    const Solution& b = w.solutions[bestI];
                    const long da = std::labs(static_cast<long>(a.seq.size()) - static_cast<long>(expect));
                    const long db = std::labs(static_cast<long>(b.seq.size()) - static_cast<long>(expect));
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
                    // T09: a truncated search saw only a prefix of the paths, so "one path"
                    // and "dominates the rest" are unproven. Under the strict flag that is
                    // ambiguity, not a closure (and the ladder does not step down after it).
                    if (w.truncated && strictBudget) {
                        res.refusedTruncated = true;
                        ++nTruncRefused;
                    } else {
                        fillOut = w.solutions[bestI].seq;
                        res.closed = true;
                        if (w.truncated) ++nTruncAccepted;
                        if (res.capped) ++nCappedAccepted;
                    }
                }
                break;
            }
            return res;
        };

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
            size_t skippedHere = 0;
            for (uint32_t r : bucket[g]) {
                recruitedBases += reads.length(r);
                auto count = [&](const Kmer& km, uint32_t) {
                    uint16_t& c = local.counts[km];
                    if (c < kCountCeil) ++c;
                };
                if (hadInputN(r)) {
                    // T27: an input N read as 'A' is a base no read observed; restart the
                    // window there instead. Correction never masks such a read, so its
                    // ambiguous positions are exactly its input Ns.
                    const uint32_t len = reads.length(r);
                    size_t emitted = 0;
                    forEachKmer(reads, r, k, [&](const Kmer& km, uint32_t p) { ++emitted; count(km, p); });
                    const size_t windows = len >= static_cast<uint32_t>(k)
                                               ? static_cast<size_t>(len) - static_cast<size_t>(k) + 1 : 0;
                    skippedHere += windows - emitted;
                } else {
                    forEachKmerRaw(reads, r, k, count);
                }
            }
            nKmersSkippedInputN += skippedHere;

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
            const bool seedLow = local.at(seed, k) < minCount;
            const bool targetLow = local.at(target, k) < minCount;
            if (seedLow) ++nSeedGone;
            if (targetLow) ++nTargetGone;
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

            const LadderResult rel = ladder(
                seed, target, gp.len,
                gp.len + static_cast<uint32_t>(kSlackBases) + static_cast<uint32_t>(k), minCount, fill[g]);
            if (rel.truncated) ++nTruncated;
            if (rel.capped) ++nCapped;
            if (rel.closed) {
                closed[g] = 1;
                repStart[g] = gp.start;
                repLen[g] = gp.len;
            }

            // T11 (TESSERACT_FIX_GAPFILL_BACKOFF=<B>): tried only where the release anchors
            // left the gap open and an anchor is below the floor, so every release closure
            // stays exactly as it was. A flank tip the reads were already thinning out over
            // is often spelled by a couple of reads only; when it is not what the deeper
            // reads spell, no walk from it reaches the other side. The moved anchor is a
            // supported flank k-mer at most B bases from the gap, inside the recruited flank
            // window and never across an N. Branch points -- a supported k-mer with a
            // supported successor that leaves the contig -- are tried first, nearest first,
            // then the other supported k-mers: the first supported k-mer is often where the
            // tip joins a paralog, and a walk from there follows the paralog. A closure
            // re-spells the stepped-over tip bases from the same local reads, with the ladder
            // and the dominance rule unchanged; the first closure wins.
            if (!closed[g] && backoff > 0 && !rel.refusedTruncated && (seedLow || targetLow)) {
                const size_t ls = gp.start > static_cast<uint32_t>(flank)
                                      ? gp.start - static_cast<uint32_t>(flank) : 0;
                const size_t re = std::min(s.size(), static_cast<size_t>(gp.start + gp.len + flank));
                const size_t seedPos0 = gp.start - static_cast<size_t>(k);
                const size_t targetPos0 = gp.start + gp.len;
                struct Anchor { size_t pos; Kmer km; bool branch; };
                auto order = [](std::vector<Anchor>& v) {
                    std::stable_sort(v.begin(), v.end(), [](const Anchor& a, const Anchor& b) {
                        return a.branch && !b.branch;
                    });
                    if (v.size() > kBackoffCandidatesPerSide) v.resize(kBackoffCandidatesPerSide);
                };
                std::vector<Anchor> seedC, targetC;
                if (seedLow) {
                    for (int j = 1; j <= backoff; ++j) {
                        if (seedPos0 < static_cast<size_t>(j)) break;
                        const size_t pos = seedPos0 - static_cast<size_t>(j);
                        if (pos < ls || baseCode(s[pos]) < 0) break;
                        bool okk = false;
                        const Kmer km = stringToKmer(s.substr(pos, static_cast<size_t>(k)), k, okk);
                        if (!okk) break;
                        if (local.at(km, k) < minCount) continue;
                        const int next = baseCode(s[pos + static_cast<size_t>(k)]);
                        bool branch = false;
                        for (int b = 0; b < 4 && !branch; ++b)
                            if (b != next && local.at(pushBack(km, b, k), k) >= minCount) branch = true;
                        seedC.push_back({pos, km, branch});
                    }
                    order(seedC);
                }
                if (targetLow) {
                    for (int j = 1; j <= backoff; ++j) {
                        const size_t pos = targetPos0 + static_cast<size_t>(j);
                        if (pos + static_cast<size_t>(k) > re) break;
                        if (baseCode(s[pos + static_cast<size_t>(k) - 1]) < 0) break;
                        bool okk = false;
                        const Kmer km = stringToKmer(s.substr(pos, static_cast<size_t>(k)), k, okk);
                        if (!okk) break;
                        if (local.at(km, k) < minCount) continue;
                        const int prev = baseCode(s[pos - 1]);
                        bool branch = false;
                        for (int b = 0; b < 4 && !branch; ++b)
                            if (b != prev && local.at(pushFront(km, b, k), k) >= minCount) branch = true;
                        targetC.push_back({pos, km, branch});
                    }
                    order(targetC);
                }
                if (!seedC.empty()) ++nSeedBackedOff;
                if (!targetC.empty()) ++nTargetBackedOff;

                struct Attempt { size_t seedPos; Kmer seed; size_t targetPos; Kmer target; };
                std::vector<Attempt> attempts;
                for (const Anchor& a : seedC) attempts.push_back({a.pos, a.km, targetPos0, target});
                for (const Anchor& a : targetC) attempts.push_back({seedPos0, seed, a.pos, a.km});
                if (!seedC.empty() && !targetC.empty())
                    attempts.push_back({seedC.front().pos, seedC.front().km, targetC.front().pos,
                                        targetC.front().km});
                long long spent = 0;
                for (const Attempt& at : attempts) {
                    if (spent > kBackoffExpansionBudget) break;
                    const uint32_t shiftL = static_cast<uint32_t>(seedPos0 - at.seedPos);
                    const uint32_t shiftR = static_cast<uint32_t>(at.targetPos - targetPos0);
                    const uint32_t expect = gp.len + shiftL + shiftR;
                    ++nAttempts;
                    const LadderResult r = ladder(
                        at.seed, at.target, expect,
                        expect + static_cast<uint32_t>(kSlackBases) + static_cast<uint32_t>(k),
                        minCount, fill[g]);
                    spent += r.totalExpansions;
                    if (r.closed) {
                        closed[g] = 1;
                        repStart[g] = static_cast<uint32_t>(at.seedPos) + static_cast<uint32_t>(k);
                        repLen[g] = expect;
                        backedOff[g] = rel.solutions > 0 ? 1 : (rel.expansions > kMaxExpansions ? 3 : 2);
                        break;
                    }
                    // A strict-budget refusal is ambiguity; no other anchor overrides it.
                    if (r.refusedTruncated) break;
                }
            }

            if (!closed[g]) {
                if (rel.solutions > 0) ++nAmbiguous;
                else {
                    ++nNoPath;
                    if (rel.expansions > kMaxExpansions) ++nBudget;
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
    // A back-off closure replaces flank bases as well as Ns, so two gaps a few dozen bases
    // apart could claim the same bases. A dropped closure is counted as the release attempt
    // left that gap.
    if (backoff > 0) {
        std::vector<Replacement> reps(gaps.size());
        for (size_t gi = 0; gi < gaps.size(); ++gi) {
            reps[gi] = {gaps[gi].contig, repStart[gi], repLen[gi], backedOff[gi] != 0, closed[gi] != 0};
        }
        stats.backoffCancelled = cancelBackoffClashes(reps, k);
        for (size_t gi = 0; gi < gaps.size(); ++gi) {
            if (!closed[gi] || reps[gi].closed) continue;
            closed[gi] = 0;
            if (backedOff[gi] == 1) ++nAmbiguous;
            else {
                ++nNoPath;
                if (backedOff[gi] == 3) ++nBudget;
            }
            backedOff[gi] = 0;
        }
    }
    stats.gapsAmbiguous = nAmbiguous.load();
    stats.gapsNoPath = nNoPath.load() + nThinPool.load();
    stats.gapsThinPool = nThinPool.load();
    stats.gapsOutOfBudget = nBudget.load();
    stats.seedBelowFloor = nSeedGone.load();
    stats.targetBelowFloor = nTargetGone.load();
    stats.gapsTruncated = nTruncated.load();
    stats.gapsTruncatedAccepted = nTruncAccepted.load();
    stats.gapsTruncatedRefused = nTruncRefused.load();
    stats.gapsCapped = nCapped.load();
    stats.gapsCappedAccepted = nCappedAccepted.load();
    stats.seedBackedOff = nSeedBackedOff.load();
    stats.targetBackedOff = nTargetBackedOff.load();
    stats.backoffAttempts = nAttempts.load();
    stats.kmersSkippedInputN = nKmersSkippedInputN.load();
    stats.searches = nSearches.load();
    const size_t nWorked = gaps.size() - nThinPool.load();
    stats.meanLocalDepth = nWorked ? static_cast<double>(dbgDepth.load()) / static_cast<double>(nWorked) : 0;
    stats.meanFloor = nWorked ? static_cast<double>(dbgFloor.load()) / static_cast<double>(nWorked) : 0;

    // ---- 5. splice -------------------------------------------------------
    // Right to left, so an edit never moves a gap that has not been done yet.
    for (size_t gi = gaps.size(); gi-- > 0;) {
        if (!closed[gi]) continue;
        Gap& gp = gaps[gi];
        contigs[gp.contig].replace(repStart[gi], repLen[gi], fill[gi]);
        ++stats.gapsClosed;
        stats.nBasesRemoved += gp.len;
        stats.basesInserted += fill[gi].size();
        if (backedOff[gi]) {
            stats.trimmedBp += repLen[gi] - gp.len;
            ++stats.closedAfterBackoff;
        }
    }

    stats.seconds = timer.elapsed();
    return stats;
}

}  // namespace ts
