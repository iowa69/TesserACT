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
    size_t seedBelowFloor = 0;     // left anchor absent from the local pool
    size_t targetBelowFloor = 0;   // right anchor absent from the local pool
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
GapFillStats closeGaps(std::vector<std::string>& contigs, const SequenceStore& reads,
                       int threads, int k, int flank);

}  // namespace ts
