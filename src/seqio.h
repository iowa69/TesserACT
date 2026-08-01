// Read storage. Reads are 2-bit packed and held in memory for the whole run,
// because every stage after counting (error correction, unitig mapping, paired
// end resolution, polishing) needs random access to them.
#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "kmer.h"

namespace ts {

// 3' quality trimming applied as reads are loaded, before anything is counted.
struct QualityTrim {
    bool enabled = true;
    int windowSize = 4;
    int meanQuality = 20;
    int phredOffset = 33;
};

// One input library.
struct Library {
    std::string r1;
    std::string r2;      // empty for single-end
    bool interleaved = false;
    // Filled in by insert-size estimation.
    double insertMean = 0;
    double insertStdDev = 0;
    int readLength = 0;
    bool oriented = true;   // true = FR (Illumina paired-end)
};

class SequenceStore {
public:
    SequenceStore() = default;

    // Loads every library. Paired reads are stored adjacently so the mate of
    // read i is always i^1, which keeps pair lookups branch-free.
    bool load(const std::vector<Library>& libs, int threads, std::string& error);

    size_t size() const { return offsets_.size() > 0 ? offsets_.size() - 1 : 0; }
    size_t totalBases() const { return totalBases_; }

    uint32_t length(size_t i) const {
        return static_cast<uint32_t>(offsets_[i + 1] - offsets_[i]);
    }

    // Base code 0..3, or -1 when the position was an ambiguous base.
    int baseAt(size_t read, uint32_t pos) const {
        uint64_t bit = offsets_[read] + pos;
        if (isAmbiguous(bit)) return -1;
        uint64_t word = data_[bit >> 5];
        return static_cast<int>((word >> ((bit & 31) * 2)) & 3);
    }

    void decode(size_t read, std::string& out) const;
    std::string decode(size_t read) const {
        std::string s;
        decode(read, s);
        return s;
    }

    bool paired() const { return paired_; }
    size_t mateOf(size_t i) const { return i ^ 1; }
    size_t pairCount() const { return paired_ ? size() / 2 : 0; }

    // Overwrite a base, used by the read error corrector.
    void setBase(size_t read, uint32_t pos, int code);

    // Marks [from, to) of a read ambiguous so every k-mer spanning it is
    // skipped. The corrector uses this to discard stretches it cannot vouch
    // for -- the k-mer-spectrum equivalent of quality trimming.
    void maskRange(size_t read, uint32_t from, uint32_t to);

    uint32_t maxReadLength() const { return maxLen_; }

    // Must be set before load() to have any effect.
    void setQualityTrim(const QualityTrim& qt) { qtrim_ = qt; }
    uint64_t trimmedBases() const { return trimmedBases_; }

private:
    bool isAmbiguous(uint64_t bit) const {
        return !ambiguous_.empty() && (ambiguous_[bit >> 6] >> (bit & 63)) & 1ULL;
    }

    std::vector<uint64_t> data_;        // 2 bits per base
    std::vector<uint64_t> ambiguous_;   // 1 bit per base, set for non-ACGT
    std::vector<uint64_t> offsets_;     // size()+1 entries, in bases
    uint64_t totalBases_ = 0;
    uint32_t maxLen_ = 0;
    bool paired_ = false;
    QualityTrim qtrim_;
    uint64_t trimmedBases_ = 0;
};

// Streams every valid k-mer of a read to `fn` as (canonicalKmer, position).
// Runs of ambiguous bases restart the rolling k-mer.
template <typename Fn>
inline void forEachKmer(const SequenceStore& store, size_t read, int k, Fn&& fn) {
    const uint32_t len = store.length(read);
    if (len < static_cast<uint32_t>(k)) return;
    Kmer fwd = 0, rc = 0;
    int valid = 0;
    for (uint32_t p = 0; p < len; ++p) {
        int c = store.baseAt(read, p);
        if (c < 0) { valid = 0; fwd = 0; rc = 0; continue; }
        fwd = pushBack(fwd, c, k);
        rc = pushFrontRc(rc, c, k);
        if (++valid >= k) {
            fn(fwd < rc ? fwd : rc, p + 1 - static_cast<uint32_t>(k));
        }
    }
}

// Writes contigs to a FASTA file with standard line wrapping.
bool writeFasta(const std::string& path, const std::vector<std::string>& seqs,
                const std::vector<std::string>& names, int lineWidth, std::string& error);

}  // namespace ts
