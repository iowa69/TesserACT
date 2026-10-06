// Phase 2, emit-B: canonical 31-mer helpers and the two count tables the end-of-run read pass fills.
//
// Encoding (identical to the kqa prototype, backward/novel/tools/kqa.cpp, and to forEachMarkerKmer in
// organism.h): A=0 C=1 G=2 T=3, the newest base in the low bits, canonical = min(forward, reverse
// complement); any other character (N, IUPAC) breaks the run. Lower-case bases count as their upper-case
// base, so a contigs.fasta record with lower-case tips gives the same k-mers as the upper-case record.
#pragma once

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

namespace ts {
namespace p2 {

constexpr int kK = 31;
constexpr uint64_t kKMask = (1ULL << (2 * kK)) - 1;

// splitmix64 finalizer: the same function as kqa's mix() and organism.h's markerHash().
inline uint64_t mix64(uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

// kqa's 1-in-S hash sampling of canonical k-mers (kqa.cpp `sampled`).
inline bool sqaSampled(uint64_t canon, uint64_t s) { return (mix64(canon ^ 0x5bd1e995ULL) % s) == 0; }

inline int baseCode2(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return -1;
    }
}

// f(pos, canon) for every valid canonical 31-mer of s[0, len); pos is the k-mer's start.
template <class F>
inline void forEachCanon31(const char* s, size_t len, F&& f) {
    uint64_t fw = 0, rv = 0;
    int valid = 0;
    constexpr int shift = 2 * (kK - 1);
    for (size_t i = 0; i < len; ++i) {
        const int c = baseCode2(s[i]);
        if (c < 0) { valid = 0; fw = rv = 0; continue; }
        fw = ((fw << 2) | static_cast<uint64_t>(c)) & kKMask;
        rv = (rv >> 2) | (static_cast<uint64_t>(3 - c) << shift);
        if (++valid >= kK) f(i + 1 - kK, fw < rv ? fw : rv);
    }
}

// Canonical 31-mer of s[pos, pos+31), or UINT64_MAX when it holds a non-ACGT character.
uint64_t canon31At(const std::string& s, size_t pos);
std::string reverseComplement(const std::string& s);

// Exact counts over a key set fixed before counting: keys are added single-threaded, then frozen;
// hit() is then thread-safe (relaxed atomic increments) and count() reads the total.
class FixedCountTable {
public:
    void addKey(uint64_t canon);
    // `filter`: also build the 2 MiB negative-test bitset (worth it for a table probed with every read
    // k-mer; small per-end tables skip it).
    void freeze(bool filter = true);
    void resetCounts();
    size_t size() const { return n_; }
    bool contains(uint64_t canon) const;
    void hit(uint64_t canon) { hitHashed(canon, mix64(canon)); }   // thread-safe after freeze()
    // The same with h = mix64(canon) already computed by the caller (the read pass shares it).
    void hitHashed(uint64_t canon, uint64_t h);
    void hitSerial(uint64_t canon);           // after freeze(); single-threaded use
    uint32_t count(uint64_t canon) const;     // 0 for an unknown key

private:
    size_t slotOf(uint64_t canon) const { return slotOfH(canon, mix64(canon)); }
    size_t slotOfH(uint64_t canon, uint64_t h) const;
    std::vector<uint64_t> pending_;
    std::vector<uint64_t> keys_;              // key + 1; 0 = empty
    std::unique_ptr<std::atomic<uint32_t>[]> counts_;
    std::vector<uint64_t> filter_;            // one bit per top-24-bit hash value: a cheap negative test
    size_t mask_ = 0, n_ = 0;
    bool frozen_ = false;
};

// The kqa table, sharded so several threads can insert: per canonical k-mer, the read count `rc`
// (saturating at 2^32-1) and the assembly multiplicity `ac` (saturating at 65535). Totals are sums,
// so the result does not depend on the thread count or the insertion order.
class ShardedCountTable {
public:
    static constexpr int kShardBits = 6;
    static constexpr size_t kShards = size_t(1) << kShardBits;
    struct Slot { uint64_t key; uint32_t rc; uint16_t ac; uint16_t pad; };

    ShardedCountTable();
    // Single-threaded, before the read pass.
    void addAssembly(uint64_t canon);
    // Thread-safe batch insert of read k-mers (each +1 to rc). `shard` = shardOf(k) of every key.
    void addReadBatch(size_t shard, const uint64_t* keys, size_t n);
    static size_t shardOf(uint64_t canon) { return static_cast<size_t>(mix64(canon) >> (64 - kShardBits)); }
    size_t distinct() const;
    size_t memoryBytes() const;
    // f(canon, rc, ac) over every stored k-mer; serial, after the pass (order is unspecified: every
    // consumer either sums or sorts what it collects).
    template <class F>
    void forEach(F&& f) const {
        for (const Shard& s : shards_)
            for (const Slot& x : s.slots)
                if (x.key) f(x.key - 1, x.rc, x.ac);
    }
    // rc of one k-mer (0 when absent).
    uint32_t readCount(uint64_t canon) const;

private:
    struct Shard {
        std::mutex m;
        std::vector<Slot> slots;
        size_t n = 0, mask = 0;
    };
    static size_t slotHash(uint64_t canon) { return static_cast<size_t>(mix64(canon)); }
    Slot& insert(Shard& s, uint64_t canon);
    void grow(Shard& s);
    std::vector<Shard> shards_;
};

}  // namespace p2
}  // namespace ts
