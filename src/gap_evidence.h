#pragma once

#include <cstddef>
#include <cstdint>
#include <set>
#include <utility>

namespace ts {
class UnitigGraph;
class SequenceStore;

// Experimental FR paired-read evidence for gap closing. The graph is unchanged.
// Every accepted read has a consistent strand-aware projection to one terminal
// end using unique graph seeds. Competing eligible partners are not arbitrated
// by node order. An end is encoded as (unitig << 1) | end.
std::set<std::pair<uint64_t, uint64_t>> nominateOrientedGapJoins(
    const UnitigGraph& graph, const SequenceStore& reads, size_t maxDistToTip,
    int probeK, uint32_t minVotes, size_t* votedOut = nullptr);
}  // namespace ts
