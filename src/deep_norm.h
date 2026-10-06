// DEEP_NORM: thinning of ultra-deep read pairs before read error correction.
//
// Phase 2, workstream plasmid (EVAL_PLAN_P2 s5.2, R6a). Off unless TESSERACT_DEEP_NORM > 0;
// with the flag unset nothing below runs and the reads are untouched.
//
// Why. A small multicopy plasmid at 100-500x the chromosome's depth carries every sequencing
// error tens to hundreds of times: each clears any abundance cutoff, and the replicon
// assembles as a tangle of error branches instead of a circle. Thinning only those reads to
// about the chromosome's depth puts the errors back under the cutoff.
//
// What. Every read gets the median count, in the correction k-mer table, of the k-mers that
// end at every 4th read position. A pair (decided once, at its even mate) or a single-end
// read is a candidate when its depth reaches TESSERACT_DEEP_NORM x the coverage peak, and is
// kept with probability TESSERACT_DEEP_NORM_TARGET x peak / depth, decided by a hash of the
// pair (or read) index so the run is reproducible. A dropped read is masked whole, which
// every later stage skips; its mate is masked with it.
//
// Two mate rules, chosen by TESSERACT_DEEP_NORM_MATE:
//   * both  -- the rule of the dn7 binary (md5 96630dc0) that produced arm sg_dn30b: BOTH
//              mates' medians must reach the selection depth, and the keep probability uses
//              the LOWER median. A read in a short high-copy chromosomal repeat, whose mate
//              lies in unique flank, is never selected.
//   * F     -- (a number in [0, 1]; the default, 1/6) the rule of dev/salfix b3b6d41 (dn8,
//              md5 0517b8aa): one mate must reach the selection depth and the other F x it,
//              and the keep probability uses the HIGHER median.
// Both rules are exact ports: the decision and its arithmetic are those of the two binaries.
#pragma once

#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <vector>

namespace ts {

class SequenceStore;
class KmerTable;

struct DeepNormConfig {
    double sel = 0.0;            // TESSERACT_DEEP_NORM: select at >= sel x peak; 0 = off
    double target = 2.0;         // TESSERACT_DEEP_NORM_TARGET: keep about target x peak
    bool both = false;           // TESSERACT_DEEP_NORM_MATE=both (dn7)
    double mate = 1.0 / 6.0;     // TESSERACT_DEEP_NORM_MATE=<F> (b3b6d41); unused when `both`
    bool enabled() const { return sel > 0; }
};

// Reads the three flags (validated by the env:: readers).
DeepNormConfig deepNormConfigFromEnv();

struct DeepNormStats {
    bool ran = false;            // the selection loop ran (enabled, a positive peak, and a
                                 // read-correction stage to run it in)
    double peak = 0.0;
    size_t selectedPairs = 0, thinnedPairs = 0;     // pairs (decision units of 2 reads)
    size_t selectedSingle = 0, thinnedSingle = 0;   // single-end reads
    size_t candidates() const { return selectedPairs + selectedSingle; }   // b3b6d41's count
    size_t dropped() const { return thinnedPairs + thinnedSingle; }        // b3b6d41's count
    size_t maskedReads() const { return 2 * thinnedPairs + thinnedSingle; }
};

enum class DeepNormVerdict { Skip, Keep, Thin };

// One decision. A pair passes its two mates' medians (medA = even mate, medB = odd mate) and
// its pair index; a single-end read passes its median as medA (medB is ignored) and its read
// index. `peak` is the correction k-mer table's coverage peak.
DeepNormVerdict deepNormDecide(const DeepNormConfig& c, double peak, float medA, float medB,
                               bool paired, uint64_t unitIndex);

// Per read: the median count in `trusted` of the k-mers ending at every 4th position (0 for a
// read shorter than k or without a valid k-mer). Threads split the reads; the result does not
// depend on the thread count.
std::vector<float> deepNormMedians(const SequenceStore& reads, const KmerTable& trusted, int k,
                                   int threads);

// Selects and masks. Does nothing (ran = false) unless c.enabled() and peak > 0.
DeepNormStats deepNormApply(SequenceStore& reads, const KmerTable& trusted, int k, int threads,
                            double peak, const DeepNormConfig& c);

// The registered counter line, printed on every run (EVAL_PLAN_P2 s5.2, s5.3 IDW):
//   [deep_norm_counts] enabled=E ran=R rule=off|both|mate:F selected_pairs=.. thinned_pairs=..
//                      selected_single=.. thinned_single=.. masked_reads=..
// With the flag off it reads enabled=0 ran=0 rule=off and zeros.
void deepNormPrintCounts(std::FILE* out, const DeepNormConfig& c, const DeepNormStats& s);

}  // namespace ts
