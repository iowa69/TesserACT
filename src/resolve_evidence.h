// Internal scoring experiment for branch selection; no reference information.
#pragma once

#include <algorithm>
#include <cstdint>
#include <vector>

namespace ts {
namespace detail {

struct BranchEvidence {
    uint64_t source = 0;
    bool repeat = false;
    std::vector<double> scores;  // one entry per enumerated candidate route
};

// A collapsed-repeat source that supports every distinct destination carries no
// information about which copy's entry the current chain used. Its count
// imbalance can reflect copy number or sampling, so do not let it dominate
// distinctive source evidence. Routes sharing a destination are one alternative.
// Unique sources and repeat sources distinguishing alternatives retain all votes.
inline std::vector<bool> sharedRepeatSources(const std::vector<uint64_t>& destinations,
                                            const std::vector<BranchEvidence>& evidence) {
    std::vector<uint64_t> distinct = destinations;
    std::sort(distinct.begin(), distinct.end());
    distinct.erase(std::unique(distinct.begin(), distinct.end()), distinct.end());
    std::vector<bool> shared(evidence.size(), false);
    if (distinct.size() < 2) return shared;
    for (size_t m = 0; m < evidence.size(); ++m) {
        if (!evidence[m].repeat) continue;
        bool everyDestination = true;
        for (uint64_t destination : distinct) {
            bool supported = false;
            for (size_t i = 0; i < destinations.size(); ++i) {
                if (destinations[i] == destination && evidence[m].scores[i] > 0.0) {
                    supported = true;
                    break;
                }
            }
            if (!supported) { everyDestination = false; break; }
        }
        shared[m] = everyDestination;
    }
    return shared;
}

inline std::vector<double> distinctiveBranchScores(const std::vector<uint64_t>& destinations,
                                                   const std::vector<BranchEvidence>& evidence) {
    const std::vector<bool> shared = sharedRepeatSources(destinations, evidence);
    std::vector<double> scores(destinations.size(), 0.0);
    for (size_t m = 0; m < evidence.size(); ++m) {
        if (shared[m]) continue;
        for (size_t i = 0; i < scores.size(); ++i) scores[i] += evidence[m].scores[i];
    }
    return scores;
}

}  // namespace detail
}  // namespace ts
