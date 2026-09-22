// Experimental measurement of final-graph k-mer multiplicity in the current
// corrected, mask-respecting read store. No graph mutation or copy inference.
#pragma once
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace ts {
class UnitigGraph;
class SequenceStore;

struct ReadCoverageMemoryPlan {
    size_t targetSlots = 0;
    size_t allocationBytes = 0; // table + returned vector + allocator allowance
};

// Public preflight for deterministic budget/overflow validation without an
// allocation. Position bound counts all len-k+1 windows, including invalid ones.
// The table never grows and never inserts a read-only/non-graph k-mer.
bool planReadCoverageMemory(size_t graphPositionBound, size_t nodeCount,
                            size_t byteBudget, ReadCoverageMemoryPlan& plan,
                            std::string& error);

struct ReadCoverageStats {
    uint64_t graphPositions = 0, distinctTargets = 0;
    uint64_t zeroTargets = 0, zeroPositions = 0;
    uint64_t readWindows = 0, matchedReadWindows = 0;
    size_t liveNodes = 0, zeroDepthNodes = 0, noValidKmerNodes = 0;
    size_t targetSlots = 0, allocationBytes = 0;
};
struct ReadCoverageResult {
    // One value per slot. Deleted slots preserve their previous coverage;
    // live nodes with no valid k-mer (including len<k) explicitly receive zero.
    std::vector<double> nodeDepths;
    ReadCoverageStats stats;
};

// Each valid corrected-read occurrence contributes one count; both mates and
// repeated occurrences count separately. Canonical graph duplicates share the
// same observed count, which contributes at each node position. This is read
// k-mer multiplicity, not fragment depth or a genomic copy-number estimate.
// Input and result are unchanged on failure. Counts and allocation arithmetic
// are overflow-checked; no synthetic carry count is added or subtracted.
bool measureObservedGraphCoverage(const UnitigGraph& graph,
                                   const SequenceStore& reads,
                                   size_t byteBudget,
                                   ReadCoverageResult& result,
                                   std::string& error);
} // namespace ts
