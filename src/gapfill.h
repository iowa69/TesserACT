// Scaffold gap closing.
//
// Scaffolding joins two chains whose adjacency the paired reads vouch for but
// which no graph path connects, and writes the space between them as a run of
// Ns. The sequence is usually not missing from the data -- it is missing from
// the *graph*, because the k-mers covering it fell below the abundance cutoff
// or sat in a branch the walker refused. Reads anchored on the two flanks
// still carry it.
//
// This stage recruits those reads, rebuilds a small de Bruijn graph from them
// alone -- where a locally sane abundance threshold replaces the global one --
// and looks for a unique path from the left flank into the right. When one
// exists the Ns are replaced by real sequence, which is what turns scaffold
// contiguity into contig contiguity.
#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "seqio.h"

namespace ts {

struct GapFillStats {
    size_t gapsSeen = 0;
    size_t gapsClosed = 0;         // replaced with real sequence
    size_t gapsAmbiguous = 0;      // more than one path fit; left as Ns
    size_t gapsNoPath = 0;         // nothing spanned it
    size_t gapsThinPool = 0;       // too few (or absurdly many) reads recruited
    size_t gapsOutOfBudget = 0;    // search gave up before exhausting the space
    // T09: the search that decided a gap stopped at the expansion budget while branches
    // were still unexplored, so its solution set is a prefix of the space, not the space.
    // Counted on every run; TESSERACT_FIX_GAPFILL_STRICT_BUDGET=1 refuses such closures.
    bool strictBudget = false;
    size_t gapsTruncated = 0;          // release-anchor deciding search truncated (any outcome)
    size_t gapsTruncatedAccepted = 0;  // closures accepted from a truncated search
    size_t gapsTruncatedRefused = 0;   // closures withheld by the strict flag
    size_t gapsCapped = 0;             // release-anchor deciding search stopped at kMaxSolutions
    size_t gapsCappedAccepted = 0;     // closures accepted from a capped search (measured only)
    size_t searches = 0;               // walks run, all floors and anchors (T38: no repeated floor)
    size_t seedBelowFloor = 0;     // left anchor absent from the local pool
    size_t targetBelowFloor = 0;   // right anchor absent from the local pool
    // T11: TESSERACT_FIX_GAPFILL_BACKOFF=<B>. Tried only where the release anchors left a
    // gap open and an anchor sat below the floor; every release closure is unchanged.
    int backoff = 0;                   // B in bases, 0 = off
    size_t seedBackedOff = 0;          // gaps where a moved left anchor was tried
    size_t targetBackedOff = 0;        // gaps where a moved right anchor was tried
    size_t backoffAttempts = 0;        // ladders run from moved anchors
    size_t closedAfterBackoff = 0;     // closures that needed a moved anchor
    size_t trimmedBp = 0;              // flank-tip bases re-spelled by those closures
    size_t backoffCancelled = 0;       // back-off closures dropped for overlapping a neighbour
    // T27: TESSERACT_FIX_GAPFILL_SKIP_INPUT_N=1. A read that carried a non-ACGT base in the
    // input is read with its N positions skipped instead of as 'A'.
    bool skipInputN = false;
    bool inputNProvenance = false;     // the caller supplied which reads had an input N
    size_t readsInputN = 0;            // reads flagged by that provenance
    size_t recruitedInputN = 0;        // of them, recruitments into some gap's pool
    size_t kmersSkippedInputN = 0;     // k-mer windows not counted because they span an input N
    double meanLocalDepth = 0;
    double meanFloor = 0;
    size_t readsRecruited = 0;
    size_t nBasesRemoved = 0;      // Ns that became sequence
    size_t basesInserted = 0;      // real bases written in their place
    double seconds = 0;
};

// Rewrites `contigs` in place, replacing closable N runs with sequence.
// `k` is the de Bruijn size used for the local reassembly; `flank` is how much
// sequence either side of a gap is used to recruit reads and anchor the walk.
//
// `inputAmbiguousReads`, when given, has one entry per read: nonzero for a read that held a
// non-ACGT base in the INPUT (before correction masked anything). It is consulted only under
// TESSERACT_FIX_GAPFILL_SKIP_INPUT_N; see readsWithInputAmbiguity().
GapFillStats closeGaps(std::vector<std::string>& contigs, const SequenceStore& reads,
                       int threads, int k, int flank,
                       const std::vector<uint8_t>* inputAmbiguousReads = nullptr);

// Organism Model 2.0 C1e: the same, but a gap whose entry in `allow` is 0 is never attempted
// (entries follow the order in which N-runs occur: contig by contig, left to right). `blocked`
// receives how many gaps were skipped. With `allow` null this is exactly the call above.
GapFillStats closeGaps(std::vector<std::string>& contigs, const SequenceStore& reads,
                       int threads, int k, int flank,
                       const std::vector<uint8_t>* inputAmbiguousReads,
                       const std::vector<uint8_t>* allow, size_t* blocked);

// One entry per read, nonzero where the read has any ambiguous position. Called right after
// loading and before read correction, the ambiguity bits are exactly the input Ns: the
// corrector only ever masks reads that were all ACGT, and it skips reads holding an N.
std::vector<uint8_t> readsWithInputAmbiguity(const SequenceStore& reads, int threads);

// The flag state a run would use, with every counter zero: for a run in which gap closing
// does not execute, so its counter lines are still printed.
GapFillStats gapFillIdleStats();

// The always-on counter lines of the G-emit gap-filler fixes ([gapfill-budget] T09,
// [gapfill-backoff] T11, [gapfill-inputn] T27), zeros included.
void logGapFillCounters(const GapFillStats& s);

// Internal to closeGaps (the T11 back-off), exposed for tests (build_v3: the campaign
// package's PKG-GAP test drives it). One gap's closure: the stretch [repStart, repStart +
// repLen) of `contig` it replaces, in gap order (contig, then position). Walking each contig
// left to right, a closure's stretch must end at or before the next closure's seed k-mer
// starts (repStart - k). A back-off closure that breaks this is dropped (closed = false) --
// the later one when both are -- and release closures never are. Returns the number dropped.
// The rule is the package's, kept verbatim by G-emit.
struct GapClosure {
    uint32_t contig = 0;
    uint32_t repStart = 0;
    uint32_t repLen = 0;
    bool backedOff = false;
    bool closed = false;
};
size_t cancelBackoffClashes(std::vector<GapClosure>& closures, int k);

}  // namespace ts
