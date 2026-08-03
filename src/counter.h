// Parallel k-mer counting.
//
// Counting runs into sharded hash tables so threads almost never contend, then
// the solid (above-threshold) k-mers are compacted into a single dense table
// that the graph builder queries. The intermediate tables are freed as soon as
// the compact one is built, which is what keeps peak memory near the size of
// the genome rather than the size of the error cloud.
#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "kmer.h"
#include "seqio.h"

namespace ts {

// Dense open-addressed table mapping k-mer -> count. A count of 0 marks an
// empty slot, so a k-mer must have a positive count to be stored.
class KmerTable {
public:
    KmerTable() = default;
    explicit KmerTable(size_t expectedEntries) { reserve(expectedEntries); }

    void reserve(size_t expectedEntries);

    // Inserts or overwrites. Not thread-safe.
    void put(Kmer key, uint32_t count);

    uint32_t get(Kmer key) const {
        if (mask_ == 0) return 0;
        size_t i = kmerHash(key) & mask_;
        while (true) {
            if (counts_[i] == 0) return 0;
            if (keys_[i] == key) return counts_[i];
            i = (i + 1) & mask_;
        }
    }

    bool contains(Kmer key) const { return get(key) != 0; }

    // Marks a k-mer as removed. Uses backward-shift deletion so lookups stay
    // correct without tombstones.
    void erase(Kmer key);

    size_t size() const { return size_; }
    size_t capacity() const { return mask_ ? mask_ + 1 : 0; }

    template <typename Fn>
    void forEach(Fn&& fn) const {
        if (counts_.empty()) return;   // mask_ == 0 is also the empty state
        for (size_t i = 0; i <= mask_; ++i) {
            if (counts_[i] != 0) fn(keys_[i], counts_[i]);
        }
    }

    void clear();

private:
    void growIfNeeded();

    std::vector<Kmer> keys_;
    std::vector<uint32_t> counts_;
    size_t mask_ = 0;
    size_t size_ = 0;
};

struct CountingStats {
    uint64_t totalKmers = 0;        // every k-mer instance observed
    uint64_t distinctKmers = 0;     // distinct before filtering
    uint64_t solidKmers = 0;        // distinct after the abundance cutoff
    uint32_t cutoff = 0;            // the threshold actually applied
    double peakCoverage = 0;        // coverage mode of solid k-mers
    double estimatedGenomeSize = 0;
    std::vector<uint64_t> histogram;   // count -> distinct k-mers with it
};

class KmerCounter {
public:
    KmerCounter(int k, int threads);

    // Abort counting if resident memory exceeds this many bytes (0 = no limit).
    // Counting is where memory peaks, so this is the only place worth guarding.
    void setMemoryLimit(long long bytes) { memoryLimit_ = bytes; }
    bool exceededMemory() const { return exceeded_; }
    // Declared here and defaulted in the .cpp so `Shard` may stay incomplete in
    // every other translation unit.
    ~KmerCounter();

    // Counts all k-mers in `store`. `extraSeqs` lets earlier assembly stages
    // feed contigs back in as additional evidence for the next k.
    void count(const SequenceStore& store, const std::vector<std::string>& extraSeqs,
               uint32_t extraWeight);

    // Chooses an abundance cutoff from the count histogram and builds the
    // compact table of solid k-mers. Pass 0 for `forcedCutoff` to auto-detect.
    void extractSolid(uint32_t forcedCutoff, KmerTable& out);

    const CountingStats& stats() const { return stats_; }

    // Cutoff chosen by locating the error/signal valley in the histogram.
    static uint32_t chooseCutoff(const std::vector<uint64_t>& histogram, double& peakOut);

private:
    struct Shard;

    int k_;
    int threads_;
    std::vector<std::unique_ptr<Shard>> shards_;
    CountingStats stats_;
    long long memoryLimit_ = 0;
    bool exceeded_ = false;
};

}  // namespace ts
