// Exact single-read graph paths. This module is not called by default.
#pragma once
#include <cstddef>
#include <cstdint>
#include <vector>

namespace ts {
class UnitigGraph;
class SequenceStore;
using ReadThreadPath = std::vector<uint64_t>;  // (node << 1) | traversal orientation
struct ReadThreadLimits {
    size_t maxSearchStates = 4096;
    size_t maxPathNodes = 64;
};
struct ReadThreadRoute {
    // Canonical between this full path and its reverse-conjugate. Includes
    // both flanks and every interior node; endpoints alone are insufficient.
    ReadThreadPath oriented;
    // Pair p has ID p. Unpaired read r has ID pairCount + r - pairedReads.
    // Each molecule occurs at most once in a route; IDs allow a downstream
    // caller to exclude molecules already used in paired-end evidence.
    std::vector<size_t> fragments;
};
struct ReadThreadStats {
    size_t readsExamined = 0, readsUnknown = 0, readsNoFlanks = 0;
    size_t readsNoExactPath = 0, readsAmbiguous = 0, readsSearchLimited = 0;
    size_t readsAccepted = 0, fragmentsAccepted = 0;
    size_t fragmentsConflicting = 0, matesDeduplicated = 0;
};
struct ReadThreadEvidence {
    std::vector<ReadThreadRoute> routes;
    ReadThreadStats stats;
};
ReadThreadPath reverseReadThreadPath(const ReadThreadPath& path);

// The caller supplies endpoint eligibility (e.g. its existing length/depth
// rules). Anchors must occur at least twice at coherent distinct read offsets
// on both flanks and be globally unique exact min(31, graph.k())-mers.
// Every observed read base and every traversed k-1 overlap must match exactly.
// Any ambiguous complete placement or exhausted search is rejected. The read
// must begin within its first eligible flank; no upstream prefix is inferred.
// If two mates yield different accepted routes, this molecule is rejected.
// No insert model, reference sequence, graph mutation or support decision is
// involved. Integration must separately reject conflicting routes and exclude
// fragments already supplying paired support; these are not insert spans.
ReadThreadEvidence collectReadThreadEvidence(
    const UnitigGraph& graph, const SequenceStore& reads,
    const std::vector<uint8_t>& eligibleEndpoint,
    const ReadThreadLimits& limits = ReadThreadLimits{});
}  // namespace ts
