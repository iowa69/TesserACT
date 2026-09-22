#include "resolve_route_density.h"
#include <cmath>
#include <iostream>
#include <stdexcept>

int main() {
    using namespace ts::detail;
    int checks = 0;
    auto check = [&](bool good, const char* message) {
        ++checks;
        if (!good) throw std::runtime_error(message);
    };
    // The abundant short node cannot contain the long mode. Raw counts skew
    // 6601:1401 although the underlying two fragment modes have equal density.
    std::vector<RouteTrainingSequence> population{{550,100}, {2000,1}};
    std::vector<uint64_t> h(1001, 0);
    h[500] = 6601; h[600] = 1401;
    const auto bimodal = fitRouteInsertDensity(h, population, 300, 900);
    check(bimodal.usable, "censored model unusable");
    check(std::abs(bimodal.opportunity[500] - 6601) < 1e-9, "short-mode opportunity wrong");
    check(std::abs(bimodal.opportunity[600] - 1401) < 1e-9, "long-mode opportunity wrong");
    check(std::abs(bimodal.at(500) - bimodal.at(600)) < 1e-12, "depth/censor correction loses equal modes");
    auto votes = routeDistanceAllocation(bimodal, {400}, {100,200});
    check(std::abs(votes[0] + votes[1] - 1) < 1e-12, "pair counted more than once");
    check(decisiveRouteLength(votes, 1.02) == -1, "symmetric modes create a route winner");
    check(bimodal.at(550) > 0 && std::isfinite(bimodal.at(550)), "histogram hole has infinite/zero density");
    votes = routeDistanceAllocation(bimodal, {450}, {100,200});
    check(decisiveRouteLength(votes, 1.02) == -1, "equal-density holes create a route winner");
    check(routeDistanceAllocation(bimodal, {400}, {100}).front() == 0, "one distinct length is discriminated");
    check(decisiveRouteLength(routeDistanceAllocation(bimodal, {900}, {100,200}), 1.02) == -1,
          "unobserved/out-of-window alternatives create a winner");

    std::vector<uint64_t> flat(1001, 0);
    for (int l = 300; l <= 900; ++l) flat[l] = 2001 - l;
    const auto broad = fitRouteInsertDensity(flat, {{2000,1}}, 300, 900);
    check(broad.usable, "broad model unusable");
    check(decisiveRouteLength(routeDistanceAllocation(broad, {400}, {100,200}), 1.02) == -1,
          "flat corrected density creates a winner");
    check(decisiveRouteLength(routeDistanceAllocation(broad, {100}, {200,210}), 1.02) == -1,
          "fit boundary creates a spurious route winner from flat corrected density");
    check(decisiveRouteLength(routeDistanceAllocation(broad, {690}, {200,210}), 1.02) == -1,
          "upper fit boundary creates a spurious route winner from flat corrected density");
    auto low = h; low[500] = 2; low[600] = 1;
    check(!fitRouteInsertDensity(low, population, 300, 900).usable, "low training enables model");
    check(!fitRouteInsertDensity({}, population, 300, 900).usable, "no data enables model");
    check(!fitRouteInsertDensity(h, {{2000,0}}, 300, 900).usable, "zero exposure enables model");
    check(!fitRouteInsertDensity(h, {{2000,NAN}}, 300, 900).usable, "nonfinite exposure enables model");
    std::vector<uint64_t> censored(601,0); censored[500] = 5000;
    const auto shortOnly = fitRouteInsertDensity(censored, {{550,1}}, 300,900);
    check(shortOnly.usable, "single-mode training unexpectedly unusable");
    check(decisiveRouteLength(routeDistanceAllocation(shortOnly, {400}, {100,200}), 1.02) == -1,
          "training population unable to see long fragments favors short route");
    auto longFavored = h; longFavored[600] *= 4;
    const auto longModel = fitRouteInsertDensity(longFavored, population,300,900);
    check(decisiveRouteLength(routeDistanceAllocation(longModel,{400,400,400},{100,200}),1.02)==1,
          "long mode support ignored");
    auto shortFavored = h; shortFavored[500] *= 4;
    const auto shortModel = fitRouteInsertDensity(shortFavored,population,300,900);
    check(decisiveRouteLength(routeDistanceAllocation(shortModel,{400,400,400},{100,200}),1.02)==0,
          "short mode support ignored");
    const auto forward = routeDistanceAllocation(longModel,{390,400,410},{100,200});
    const auto reordered = routeDistanceAllocation(longModel,{410,390,400},{200,100});
    check(std::abs(forward[0]-reordered[1]) < 1e-12 && std::abs(forward[1]-reordered[0]) < 1e-12,
          "order changes evidence allocation");
    // Exactly mappable base intervals, including the k-1 terminal bases.
    check(uniqueSeedIntervals(std::vector<uint8_t>(520,0),31).empty(),"inaccessible short node has training exposure");
    const auto oneInterval=uniqueSeedIntervals(std::vector<uint8_t>(1970,1),31);
    check(oneInterval.size()==1 && oneInterval[0].begin==0 && oneInterval[0].end==2000,
          "unique start coordinates do not include k-1 bases");
    check(containedTrainingFragment(oneInterval,0,2000),"fully contained fragment rejected");
    check(!containedTrainingFragment(oneInterval,0,2001),"overhanging fragment admitted");
    const auto splitIntervals=uniqueSeedIntervals({1,1,0,1,1},3);
    check(splitIntervals.size()==2 && splitIntervals[0].end==4 && splitIntervals[1].begin==3 && splitIntervals[1].end==7,
          "ambiguous seed does not split calibration intervals");
    check(!containedTrainingFragment(splitIntervals,1,6),"fragment crossing ambiguous seed admitted");
    std::vector<uint64_t> accessible(1001,0);accessible[500]=1501;accessible[600]=1401;
    const auto accessibleModel=fitRouteInsertDensity(accessible,{{oneInterval[0].end,1}},300,900);
    check(accessibleModel.usable && std::abs(accessibleModel.at(500)-accessibleModel.at(600))<1e-12,
          "excluding inaccessible short nodes fails equal-mode recovery");

    // The former absolute-allocation floor was inflated by uninformative pairs.
    RouteInsertDensity simple;
    simple.usable=true; simple.probability.assign(1100,1);simple.opportunity.assign(1100,1);
    simple.localEffectiveObservations.assign(1100,1000);
    simple.probability[700]=10;
    std::vector<int> padded(40,400);padded.push_back(600);
    std::vector<double> contrast;
    auto allocated=routeDistanceAllocation(simple,padded,{100,200},&contrast,2);
    check(decisiveRouteLength(allocated,1.02)==0 && allocated[0]>2,
          "flat-padding counterexample no longer exercises the original weakness");
    check(supportedRouteLength(allocated,contrast,1.02,2)==-1,"flat pairs inflate discriminating support");
    allocated=routeDistanceAllocation(simple,{600,600},{100,200},&contrast,2);
    check(supportedRouteLength(allocated,contrast,1.02,2)==-1,"two fractional observations silently rounded to bar2");
    allocated=routeDistanceAllocation(simple,{600,600,600},{100,200},&contrast,2);
    check(supportedRouteLength(allocated,contrast,1.02,2)==0,"three strong observations do not clear bar2");

    RouteInsertDensity three;
    three.usable=true;three.probability.assign(1200,1);three.opportunity.assign(1200,1);
    three.localEffectiveObservations.assign(1200,1000);
    three.probability[700]=.001;three.probability[900]=10;three.probability[1100]=.001;
    std::vector<int> sharedAB(40,400);sharedAB.push_back(800);
    allocated=routeDistanceAllocation(three,sharedAB,{100,200,300},&contrast,2);
    check(decisiveRouteLength(allocated,1.02)==0 && allocated[0]>2,
          "three-route common-evidence fixture misses its intended score imbalance");
    check(supportedRouteLength(allocated,contrast,1.02,2)==-1,
          "A/B-common evidence becomes decisive merely because a third route is weaker");
    const auto reorderedThree=routeDistanceAllocation(three,sharedAB,{300,100,200},&contrast,2);
    check(supportedRouteLength(reorderedThree,contrast,1.02,2)==-1,
          "three-route alternative ordering changes insufficient-contrast fallback");
    std::vector<uint64_t> tail(951,0);tail[500]=1000000;tail[900]=1;
    const auto tailModel=fitRouteInsertDensity(tail,{{550,1000},{1000,1}},300,950);
    check(tailModel.usable && tailModel.effectiveObservations>100000,"tail fixture lacks large global ESS");
    check(std::abs(tailModel.localEffectiveObservations[900]-1)<1e-9,"isolated tail training count inflated");
    allocated=routeDistanceAllocation(tailModel,{700,700,700,700},{100,200},&contrast,2);
    check(supportedRouteLength(allocated,contrast,1.02,2)==-1,"one tail observation certifies route density");
    std::cout << checks << " route density checks passed\n";
}
