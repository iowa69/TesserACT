// Phase 2, emit-B: one pass over the RAW input read files at the end of a run.
//
// Three report features need the reads as sequenced, not the corrected, masked and trimmed copy the
// assembler keeps in memory:
//   * F2 self-QA (kqa port): hash-sampled (1 in 8) canonical 31-mer read counts;
//   * ends.tsv unsupported tips: exact read counts of every 31-mer near a contig end;
//   * F9 detection report: counts of the detection sketch's marker k-mers.
// They share one pass so the files are decompressed once. One reader thread per file decompresses and
// cuts records into batches; `threads` workers count k-mers. Every count is a sum over reads, so the
// result is the same at every thread count.
//
// Input formats: FASTQ (4-line records; the format of the first non-empty line decides) or FASTA
// (multi-line records joined), gzipped or not. Every base of every record is used: no quality
// trimming, no correction (the kqa prototype streamed the files the same way).
#pragma once

#include <atomic>
#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "p2_kmer.h"

namespace ts {
namespace p2 {

struct ReadPassSinks {
    ShardedCountTable* sampled = nullptr;   // kqa table (read counts of sampled k-mers)
    uint64_t sampleS = 8;
    FixedCountTable* exact = nullptr;       // exact counts of a fixed k-mer set
    // Detection markers: k-mer -> slot, per-slot counts, and one hash threshold per sketch density.
    // A k-mer is counted once per density whose threshold its markerHash (= mix64) does not exceed,
    // exactly as organism_detect.cpp's countSequence does.
    const std::unordered_map<uint64_t, uint32_t>* detectSlot = nullptr;
    std::vector<uint64_t> detectThresholds;
    std::atomic<uint32_t>* detectCounts = nullptr;
    // Locus reads (ends.tsv tips): a callback for the records whose index in the assembler's read store
    // has its bit set in `wantBits`. File f's record j is store read layout[f].first + j * layout[f].second
    // (SequenceStore's layout); a file whose layout entry has stride 0 is not collected. FASTQ records are
    // trimmed at the 3' end while the base quality is below `trimQ` (Phred+33; the sp2 prototype's qtrim)
    // before the callback sees them. Called from reader threads: it must be thread-safe.
    const std::vector<uint64_t>* wantBits = nullptr;
    std::vector<std::pair<uint64_t, uint64_t>> layout;
    int trimQ = 10;
    std::function<void(uint64_t, const char*, size_t)> onWanted;
    bool any() const { return sampled || exact || detectSlot || (wantBits && onWanted); }
};

struct ReadPassStats {
    size_t files = 0;
    uint64_t reads = 0, bases = 0, wanted = 0;
    std::vector<uint64_t> fileRecords;   // records read from each file
    double seconds = 0;
    std::string error;
};

// Runs the pass. False (with st.error) when a file cannot be opened or decompressed.
bool runReadPass(const std::vector<std::string>& files, int threads, ReadPassSinks& sinks, ReadPassStats& st);

// The same counting applied to in-memory sequences (tests and single-threaded callers).
void countSequences(const std::vector<std::string>& seqs, ReadPassSinks& sinks, ReadPassStats& st);

}  // namespace p2
}  // namespace ts
