// Internal resolver coverage calibration. Independent of graph simplification.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <utility>
#include <vector>

#include "graph.h"

namespace ts {
namespace detail {

inline double graphWeightedMedianCoverage(const UnitigGraph& graph, size_t minLength) {
    if (graph.k() <= 0) return 0.0;
    std::vector<std::pair<double, uint64_t>> depthMass;
    depthMass.reserve(graph.nodes.size());
    uint64_t total = 0;
    for (const Unitig& u : graph.nodes) {
        if (u.deleted || !std::isfinite(u.coverage) || u.coverage <= 0.0 ||
            u.seq.size() < minLength) continue;
        const uint64_t mass = u.seq.size() - static_cast<size_t>(graph.k()) + 1;
        depthMass.emplace_back(u.coverage, mass);
        total += mass;
    }
    if (total == 0) return 0.0;
    std::sort(depthMass.begin(), depthMass.end());
    const uint64_t middle = total / 2 + total % 2;
    uint64_t cumulative = 0;
    for (const auto& entry : depthMass) {
        cumulative += entry.second;
        if (cumulative >= middle) return entry.first;
    }
    return 0.0;
}

}  // namespace detail

// A unitig stores k-1 boundary bases also stored in its neighbours. Its
// independent graph length is therefore seq.size()-k+1, exactly the number
// of graph k-mers. Splitting a constant-depth unitig into arbitrary overlapping
// pieces preserves both this mass and the weighted coverage distribution.
//
// An unweighted unitig median instead gives every fragment one vote: many
// tiny low-depth pieces can classify an intact chromosome as a repeat. The
// weighted median remains robust to repeat sequence occupying a minority of
// graph length. Zero/nonfinite coverage, deleted nodes and invalid lengths
// have no calibrated mass. All valid lengths participate, including pieces
// below 2k, otherwise splitting a unitig could remove its weight entirely.
inline double graphLengthWeightedMedianCoverage(const UnitigGraph& graph) {
    if (graph.k() <= 0) return 0.0;
    return detail::graphWeightedMedianCoverage(graph, static_cast<size_t>(graph.k()));
}

// Preserve the legacy coverage estimator's >=2k population, changing only the
// weight of each valid unitig. Reference-absent short sequence can dominate
// whole-graph mass on pathological libraries; this conditional estimator does
// not admit those previously excluded pieces. Split invariance applies only
// when every resulting piece remains eligible. Empty valid populations return
// zero so the caller can keep its legacy fallback.
inline double graphEligibleLengthWeightedMedianCoverage(const UnitigGraph& graph) {
    if (graph.k() <= 0) return 0.0;
    return detail::graphWeightedMedianCoverage(graph, 2 * static_cast<size_t>(graph.k()));
}

}  // namespace ts
