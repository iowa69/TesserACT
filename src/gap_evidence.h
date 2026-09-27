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

// What the gap-close nominator did on one call (T08/T15 counters).
struct GapNominatorStats {
    int fix = 0;                 // TESSERACT_FIX_GAP_NOMINATOR level in force
    int oriented = 0;            // 1 when TESSERACT_GAP_ORIENTED=1 routed the call
    size_t tips = 0;             // dead ends seeded (legacy path only)
    size_t walked = 0;           // oriented node ids in the tip neighbourhoods
    size_t indexedKmers = 0;     // probe k-mers indexed
    size_t ambiguousKmers = 0;   // fix only: k-mers placed at two node/strand positions
    size_t strandRejects = 0;    // fix only: reads whose strand points away from every tip
};

// The paired-read nominator the gap closer calls (defined in assembler.cpp). Routes to
// nominateOrientedGapJoins when TESSERACT_GAP_ORIENTED=1; otherwise runs the legacy walk,
// corrected by TESSERACT_FIX_GAP_NOMINATOR=1 (T08; default off).
std::set<std::pair<uint64_t, uint64_t>> nominateGapJoins(
    const UnitigGraph& g, const SequenceStore& reads, size_t maxDistToTip, int probeK,
    uint32_t minVotes, size_t* votedOut, GapNominatorStats* statsOut = nullptr);
}  // namespace ts
