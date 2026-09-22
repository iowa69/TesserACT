#include <iostream>
#include <stdexcept>
#include "resolve_evidence.h"

int main() {
    using ts::detail::BranchEvidence;
    using ts::detail::distinctiveBranchScores;
    int checks = 0;
    auto expect = [&](const std::vector<double>& actual,
                      const std::vector<double>& wanted, const char* name) {
        ++checks;
        if (actual != wanted) throw std::runtime_error(name);
    };
    // A depth imbalance from a shared repeat would rank destination 10 first
    // (100 versus 54). The distinct anchor selects 20 with four real votes.
    expect(distinctiveBranchScores({10, 20}, {{2, true, {100, 50}}, {4, false, {0, 4}}}),
           {0, 4}, "shared repeat overwhelms distinct source");
    expect(distinctiveBranchScores({10, 20}, {{2, true, {100, 50}}}),
           {0, 0}, "shared repeat alone selects destination");
    expect(distinctiveBranchScores({10, 20}, {{2, true, {0, 7}}, {4, false, {3, 0}}}),
           {3, 7}, "distinctive repeat evidence was lost");
    expect(distinctiveBranchScores({10, 20}, {{4, false, {100, 50}}}),
           {100, 50}, "unique anchor evidence was changed");
    expect(distinctiveBranchScores({10, 10}, {{2, true, {100, 50}}}),
           {100, 50}, "routes to same destination mistaken for a branch");
    expect(distinctiveBranchScores({10, 10, 20}, {{2, true, {100, 0, 50}}, {4, false, {2, 2, 4}}}),
           {2, 2, 4}, "alternative routes inflate destination count");
    expect(distinctiveBranchScores({10, 20, 30}, {{2, true, {100, 50, 0}}}),
           {100, 50, 0}, "source distinguishing one destination was lost");
    expect(distinctiveBranchScores({}, {}), {}, "empty branch changed");
    std::cout << checks << " branch-evidence checks passed\n";
}
