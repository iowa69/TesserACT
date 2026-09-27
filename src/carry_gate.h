// T01: carry-only junctions (build_v3 G-graph, 2026-09-25).
//
// A rung counts the previous rung's contigs into the same table as the reads, at weight
// w = 4 per k-mer occurrence, and cuts at an abundance cutoff of about 2. So every carried
// k-mer is solid whether or not any read contains it. That is the point of the carry --
// it keeps sequence a smaller k resolved -- but it also keeps sequence a smaller k got
// WRONG. Two loci sharing an h-bp word W (h >= k0-1) are joined at k0 through W; if
// simplification at k0 then deletes the true continuations, the carried unitig spells
// A-left + W + B-right, and at every later rung k the k-h-1 k-mers that span the junction
// have zero read support and survive on the carry weight alone. Being interior to one
// unitig, no per-unitig rule ever sees them. (NODE_337 of A. baumannii GCF016919505v2.)
//
// The same zero-read signature -- an interior run of at most k-1 read-absent k-mers
// between read-supported k-mers -- is also what a genuine thin read overlap looks like,
// and keeping that bridge is what the carry is for. The two differ at the flanks: at a
// false junction the reads show a DIFFERENT continuation out of the flank (the true A-right
// or B-left), at a thin overlap they show none. So a run is gated only when a flank has a
// competing read-supported continuation.
#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "counter.h"

namespace ts {

struct CarryGateStats {
    size_t carriedKmers = 0;       // k-mer positions scanned in the carried contigs
    size_t readAbsent = 0;         // of those, solid only through the carry weight (read count 0)
    size_t runs = 0;               // maximal interior zero-read runs between read-supported k-mers
    size_t runsTooLong = 0;        // of those, longer than k-1 (not a junction signature)
    size_t runsNoCompetitor = 0;   // spared: no competing read-supported continuation at the flank(s)
    size_t candidates = 0;         // runs that meet the rule (withheld only when mode > 0)
    size_t runsGated = 0;          // runs actually withheld
    size_t kmersWithheld = 0;      // k-mers that dropped out of the solid table
};

// `carry`: the contigs that were counted into `solid` at `carryWeight` per k-mer occurrence
// (KmerCounter::count(reads, carry, carryWeight) followed by extractSolid). The read count of
// a solid k-mer x is then solid(x) - carryWeight * occurrences(x in carry), exactly.
//
// mode 0: dry run -- count everything, change nothing.
// mode 1: withhold a candidate run when EITHER flank has a competing continuation whose read
//         count is at least max(cutoff, 2).
// mode 2: require a competing continuation at BOTH flanks.
//
// Withholding removes the carry weight of each occurrence in the run; a k-mer whose count
// then falls below `cutoff` leaves the table (its read count is 0, so it always does unless
// another carried occurrence still holds it). Decisions are taken on the unmodified table
// and applied afterwards, so contig order does not matter.
CarryGateStats gateCarryOnlyJunctions(const std::vector<std::string>& carry, KmerTable& solid,
                                      int k, uint32_t carryWeight, uint32_t cutoff, int mode);

}  // namespace ts
