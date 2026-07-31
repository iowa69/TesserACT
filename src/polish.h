// Consensus polishing.
//
// Repeat resolution stitches unitigs together, and any base the graph got
// wrong survives into the contigs. Re-anchoring the reads onto the finished
// sequence and taking a per-position majority corrects those residual
// substitutions, which is what drives the mismatch rate toward zero.
#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "seqio.h"

namespace ts {

struct PolishStats {
    size_t readsUsed = 0;
    size_t basesChanged = 0;
    size_t positionsCovered = 0;
    size_t lowCoveragePositions = 0;
    double meanDepth = 0;
};

// Rewrites `contigs` in place. `minDepth` is the coverage a position needs
// before it may be changed, and `minFraction` the share of reads the winning
// base needs.
PolishStats polishContigs(std::vector<std::string>& contigs, const SequenceStore& reads,
                          int threads, int anchorK, int minDepth, double minFraction);

}  // namespace ts
