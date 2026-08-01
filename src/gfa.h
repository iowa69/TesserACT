// GFA1 export of the assembly graph.
//
// Segments are unitigs, links carry the (k-1) overlap the de Bruijn graph
// implies, and paths record how each output contig traverses the graph -- which
// is what lets a viewer such as Bandage show where a contig stops and why.
#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "graph.h"

namespace ts {

// One resolved contig expressed as a walk over oriented unitigs.
struct GfaPath {
    std::string name;
    std::vector<uint64_t> oriented;   // (unitig << 1) | orientation
    std::vector<int> gaps;            // gaps[i] = N-gap inserted before element i
};

// Writes `path` in GFA1. Returns false and fills `error` on I/O failure.
// `segmentCount` and `linkCount` receive what was emitted.
bool writeGfa(const std::string& path, const UnitigGraph& g,
              const std::vector<GfaPath>& paths, size_t& segmentCount, size_t& linkCount,
              std::string& error);

}  // namespace ts
