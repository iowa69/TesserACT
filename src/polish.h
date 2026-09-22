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

struct QualityPolishStats {
    bool enabled = false;
    bool provenanceAvailable = false;
    size_t originalBytes = 0, placementBytes = 0, candidateBytes = 0;
    size_t candidatePositions = 0, alternateProposals = 0, extraChanges = 0;
    size_t informativeFragmentSites = 0, uniformFragmentSites = 0, maxFragmentDepth = 0;
    size_t eligibleObservations = 0, excludedOriginalAmbiguous = 0, excludedMasked = 0;
    size_t excludedModified = 0, excludedMissingQuality = 0;
    size_t uniformObservations = 0, cappedObservations = 0;
    size_t belowDepth = 0, oneOrientation = 0, belowPosterior = 0, tied = 0, overflow = 0;
};

struct PolishStats {
    size_t readsUsed = 0;
    size_t basesChanged = 0;
    size_t positionsCovered = 0;
    size_t lowCoveragePositions = 0;
    double meanDepth = 0;
    QualityPolishStats quality;
};

// Rewrites `contigs` in place. `minDepth` is the coverage a position needs
// before it may be changed, and `minFraction` the share of reads the winning
// base needs.
PolishStats polishContigs(std::vector<std::string>& contigs, const SequenceStore& reads,
                          int threads, int anchorK, int minDepth, double minFraction);

}  // namespace ts
