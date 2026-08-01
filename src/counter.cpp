#include "counter.h"

#include "util.h"

#include <algorithm>
#include <atomic>
#include <mutex>
#include <thread>

namespace ts {
namespace {

constexpr double kMaxLoad = 0.6;
constexpr size_t kMinCapacity = 1024;
constexpr size_t kHistMax = 100000;

// Per-thread, per-shard staging buffer. Sized so a whole thread's staging area
// stays around 2 MiB with the default 256 shards: big enough that the shard
// mutex is taken once per 512 k-mers, small enough to stay out of the way.
constexpr size_t kShardBuffer = 512;

void rehashTable(std::vector<Kmer>& keys, std::vector<uint32_t>& counts, size_t& mask,
                 size_t newCapacity) {
    std::vector<Kmer> nk(newCapacity);
    std::vector<uint32_t> nc(newCapacity, 0);
    const size_t nmask = newCapacity - 1;
    for (size_t i = 0; i < counts.size(); ++i) {
        if (counts[i] == 0) continue;
        size_t j = kmerHash(keys[i]) & nmask;
        while (nc[j] != 0) j = (j + 1) & nmask;
        nk[j] = keys[i];
        nc[j] = counts[i];
    }
    keys.swap(nk);
    counts.swap(nc);
    mask = nmask;
}

size_t capacityFor(size_t expectedEntries) {
    size_t cap = kMinCapacity;
    if (expectedEntries) {
        const double target = static_cast<double>(expectedEntries) / kMaxLoad;
        while (static_cast<double>(cap) < target && cap <= (~static_cast<size_t>(0) >> 1)) cap <<= 1;
    }
    return cap;
}

}  // namespace

void KmerTable::reserve(size_t expectedEntries) {
    const size_t cap = capacityFor(expectedEntries);
    if (cap <= capacity()) return;
    rehashTable(keys_, counts_, mask_, cap);
}

void KmerTable::growIfNeeded() {
    const size_t cap = capacity();
    if (cap == 0) return;
    if (size_ * 5 > cap * 3) rehashTable(keys_, counts_, mask_, cap * 2);
}

void KmerTable::put(Kmer key, uint32_t count) {
    if (count == 0) {
        erase(key);
        return;
    }
    if (capacity() == 0) rehashTable(keys_, counts_, mask_, kMinCapacity);
    size_t i = kmerHash(key) & mask_;
    while (counts_[i] != 0) {
        if (keys_[i] == key) {
            counts_[i] = count;
            return;
        }
        i = (i + 1) & mask_;
    }
    keys_[i] = key;
    counts_[i] = count;
    ++size_;
    growIfNeeded();
}

void KmerTable::erase(Kmer key) {
    if (capacity() == 0) return;
    size_t hole = kmerHash(key) & mask_;
    while (true) {
        if (counts_[hole] == 0) return;
        if (keys_[hole] == key) break;
        hole = (hole + 1) & mask_;
    }
    // Backward shift: pull every following entry of the same probe run back into
    // the hole if the hole sits inside its [ideal, current] window. The load
    // factor guarantees an empty slot terminates the walk.
    size_t j = hole;
    for (;;) {
        j = (j + 1) & mask_;
        if (counts_[j] == 0) break;
        const size_t ideal = kmerHash(keys_[j]) & mask_;
        if (((j - hole) & mask_) <= ((j - ideal) & mask_)) {
            keys_[hole] = keys_[j];
            counts_[hole] = counts_[j];
            hole = j;
        }
    }
    counts_[hole] = 0;
    keys_[hole] = 0;
    --size_;
}

void KmerTable::clear() {
    std::fill(counts_.begin(), counts_.end(), 0u);
    size_ = 0;
}

struct KmerCounter::Shard {
    KmerTable table;
    std::mutex mu;
};

KmerCounter::~KmerCounter() = default;

KmerCounter::KmerCounter(int k, int threads) : k_(k), threads_(threads > 0 ? threads : 1) {
    size_t target = std::max<size_t>(256, static_cast<size_t>(8) * static_cast<size_t>(threads_));
    size_t nShards = 1;
    while (nShards < target) nShards <<= 1;
    shards_.resize(nShards);
    for (size_t i = 0; i < nShards; ++i) {
        shards_[i].reset(new Shard());
        shards_[i]->table.reserve(4096);
    }
}

void KmerCounter::count(const SequenceStore& store, const std::vector<std::string>& extraSeqs,
                        uint32_t extraWeight) {
    if (shards_.empty()) return;
    const size_t nShards = shards_.size();
    const size_t shardMask = nShards - 1;
    const size_t nReads = store.size();
    const int k = k_;

    // Blocks amortise the atomic hand-out; reads inside a block are contiguous
    // so the packed sequence buffer is walked linearly.
    constexpr size_t kBlock = 1024;
    const size_t nBlocks = (nReads + kBlock - 1) / kBlock;
    std::atomic<size_t> nextBlock{0};
    std::atomic<bool> stop{false};
    std::atomic<size_t> nextExtra{0};
    std::atomic<uint64_t> totalKmers{0};

    size_t nThreads = static_cast<size_t>(threads_);
    if (nThreads > nBlocks + extraSeqs.size()) nThreads = std::max<size_t>(1, nBlocks + extraSeqs.size());

    auto worker = [&]() {
        std::vector<Kmer> buf(nShards * kShardBuffer);
        std::vector<uint32_t> fill(nShards, 0);
        uint64_t local = 0;
        uint32_t weight = 1;

        auto flushShard = [&](size_t s) {
            const uint32_t n = fill[s];
            if (n == 0) return;
            const Kmer* p = &buf[s * kShardBuffer];
            Shard& shard = *shards_[s];
            std::lock_guard<std::mutex> lock(shard.mu);
            for (uint32_t i = 0; i < n; ++i) {
                const uint64_t c = static_cast<uint64_t>(shard.table.get(p[i])) + weight;
                shard.table.put(p[i], c > 0xFFFFFFFFull ? 0xFFFFFFFFu : static_cast<uint32_t>(c));
            }
            fill[s] = 0;
        };

        // Shard selection uses high hash bits; the tables index with low bits.
        auto push = [&](Kmer key) {
            const size_t s = static_cast<size_t>(kmerHash(key) >> 40) & shardMask;
            buf[s * kShardBuffer + fill[s]] = key;
            if (++fill[s] == kShardBuffer) flushShard(s);
            local += weight;
        };

        for (;;) {
            const size_t b = nextBlock.fetch_add(1, std::memory_order_relaxed);
            if (b >= nBlocks) break;
            // Checked once per block rather than per k-mer: reading /proc is a
            // syscall, and a block is large enough that overshoot is bounded.
            if (memoryLimit_ > 0 && !stop.load(std::memory_order_relaxed)) {
                const long long rss = util::currentMemoryBytes();
                if (rss > 0 && rss > memoryLimit_) stop.store(true, std::memory_order_relaxed);
            }
            if (stop.load(std::memory_order_relaxed)) break;
            const size_t begin = b * kBlock;
            const size_t end = std::min(nReads, begin + kBlock);
            for (size_t r = begin; r < end; ++r) {
                forEachKmer(store, r, k, [&](Kmer km, uint32_t) { push(km); });
            }
        }
        for (size_t s = 0; s < nShards; ++s) flushShard(s);
        if (stop.load(std::memory_order_relaxed)) return;

        // Contigs carried over from a smaller k are trusted evidence, so each of
        // their k-mers is worth `extraWeight` observations.
        weight = extraWeight;
        if (extraWeight > 0) {
            for (;;) {
                const size_t i = nextExtra.fetch_add(1, std::memory_order_relaxed);
                if (i >= extraSeqs.size()) break;
                const std::string& s = extraSeqs[i];
                Kmer fwd = 0, rc = 0;
                int valid = 0;
                for (size_t p = 0; p < s.size(); ++p) {
                    const int c = baseCode(s[p]);
                    if (c < 0) {
                        valid = 0;
                        fwd = 0;
                        rc = 0;
                        continue;
                    }
                    fwd = pushBack(fwd, c, k);
                    rc = pushFrontRc(rc, c, k);
                    if (++valid >= k) push(fwd < rc ? fwd : rc);
                }
            }
            for (size_t s = 0; s < nShards; ++s) flushShard(s);
        }

        totalKmers.fetch_add(local, std::memory_order_relaxed);
    };

    if (nThreads <= 1) {
        worker();
        if (stop.load(std::memory_order_relaxed)) exceeded_ = true;
    } else {
        std::vector<std::thread> pool;
        pool.reserve(nThreads);
        for (size_t t = 0; t < nThreads; ++t) pool.emplace_back(worker);
        for (auto& th : pool) th.join();
        if (stop.load(std::memory_order_relaxed)) exceeded_ = true;
    }

    uint64_t distinct = 0;
    for (const auto& s : shards_) distinct += s->table.size();
    stats_.totalKmers = totalKmers.load(std::memory_order_relaxed);
    stats_.distinctKmers = distinct;
}

uint32_t KmerCounter::chooseCutoff(const std::vector<uint64_t>& histogram, double& peakOut) {
    peakOut = 0;
    const size_t n = histogram.size();
    if (n <= 3) return 2;

    // Weighting by c picks the coverage mode rather than the error shoulder,
    // which always dominates by raw distinct-k-mer count.
    size_t peak = 0;
    unsigned __int128 best = 0;
    for (size_t c = 3; c < n; ++c) {
        const unsigned __int128 v = static_cast<unsigned __int128>(histogram[c]) * c;
        if (v > best) {
            best = v;
            peak = c;
        }
    }
    peakOut = static_cast<double>(peak);
    if (peak < 5) return 2;

    size_t valley = 1;
    uint64_t lowest = ~static_cast<uint64_t>(0);
    for (size_t c = 1; c < peak; ++c) {
        if (histogram[c] < lowest) {
            lowest = histogram[c];
            valley = c;
        }
    }

    uint32_t cutoff = static_cast<uint32_t>(valley);
    const uint32_t ceiling = std::max<uint32_t>(3, static_cast<uint32_t>(peak / 4));
    if (cutoff < 2) cutoff = 2;
    if (cutoff > ceiling) cutoff = ceiling;
    return cutoff;
}

void KmerCounter::extractSolid(uint32_t forcedCutoff, KmerTable& out) {
    std::vector<uint64_t> hist(kHistMax + 1, 0);
    const size_t nShards = shards_.size();

    if (nShards) {
        std::atomic<size_t> nextShard{0};
        std::mutex mergeMu;
        size_t nThreads = std::min<size_t>(static_cast<size_t>(threads_), nShards);
        if (nThreads < 1) nThreads = 1;
        auto worker = [&]() {
            std::vector<uint64_t> localHist(kHistMax + 1, 0);
            for (;;) {
                const size_t s = nextShard.fetch_add(1, std::memory_order_relaxed);
                if (s >= nShards) break;
                shards_[s]->table.forEach([&](Kmer, uint32_t c) {
                    localHist[c < kHistMax ? c : kHistMax] += 1;
                });
            }
            std::lock_guard<std::mutex> lock(mergeMu);
            for (size_t i = 0; i <= kHistMax; ++i) hist[i] += localHist[i];
        };
        if (nThreads <= 1) {
            worker();
        } else {
            std::vector<std::thread> pool;
            pool.reserve(nThreads);
            for (size_t t = 0; t < nThreads; ++t) pool.emplace_back(worker);
            for (auto& th : pool) th.join();
        }
    }

    double peak = 0;
    const uint32_t autoCutoff = chooseCutoff(hist, peak);
    const uint32_t cutoff = forcedCutoff ? forcedCutoff : autoCutoff;

    uint64_t solid = 0;
    for (size_t c = cutoff; c <= kHistMax; ++c) solid += hist[c];

    out.clear();
    out.reserve(solid);
    for (size_t s = 0; s < nShards; ++s) {
        shards_[s]->table.forEach([&](Kmer key, uint32_t c) {
            if (c >= cutoff) out.put(key, c);
        });
        // Released one shard at a time: the compact table and the sharded ones
        // are never both fully resident.
        shards_[s].reset();
    }
    shards_.clear();
    shards_.shrink_to_fit();

    size_t used = kHistMax + 1;
    while (used > 0 && hist[used - 1] == 0) --used;
    hist.resize(used);
    hist.shrink_to_fit();

    stats_.histogram = std::move(hist);
    stats_.solidKmers = solid;
    stats_.cutoff = cutoff;
    stats_.peakCoverage = peak;
    stats_.estimatedGenomeSize = static_cast<double>(solid);
}

}  // namespace ts
