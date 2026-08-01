// K-mer spectrum read error correction.
//
// A sequencing error creates a run of k-mers that appear nowhere else in the
// data. Anchoring on a stretch of trusted k-mers and extending outward, each
// base that breaks the run is replaced with the one that restores it. Removing
// those errors before the larger-k iterations both shrinks the error cloud the
// counter has to hold and recovers sequence that would otherwise fall below the
// abundance cutoff.
#pragma once

#include <cstdint>

#include "counter.h"
#include "seqio.h"

namespace ts {

struct CorrectionStats {
    size_t readsExamined = 0;
    size_t readsCorrected = 0;
    size_t basesCorrected = 0;
    size_t readsUncorrectable = 0;   // no trusted k-mer to anchor on
    size_t basesMasked = 0;          // stretches correction could not vouch for
};

// Rewrites `reads` in place. `solid` is the trusted k-mer set at `k`.
CorrectionStats correctReads(SequenceStore& reads, const KmerTable& solid, int k, int threads);

}  // namespace ts
