// Experimental conditional route evidence; independent of endpoint support bars.
#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <vector>

namespace ts { namespace detail {
struct RouteTrainingSequence { size_t length; double coverage; };
struct RouteTrainingInterval { size_t begin; size_t end; };
// A run of unique k-mer STARTS [a,b) covers bases [a,b+k-1).
// Every probe of a read wholly inside this interval has a unique graph locus.
inline std::vector<RouteTrainingInterval> uniqueSeedIntervals(const std::vector<uint8_t>& unique,
                                                            size_t k) {
    std::vector<RouteTrainingInterval> intervals;
    if (!k) return intervals;
    for (size_t start = 0; start < unique.size();) {
        if (!unique[start]) { ++start; continue; }
        size_t end = start + 1;
        while (end < unique.size() && unique[end]) ++end;
        intervals.push_back({start, end + k - 1});
        start = end;
    }
    return intervals;
}
inline bool containedTrainingFragment(const std::vector<RouteTrainingInterval>& intervals,
                                      size_t begin, size_t end) {
    auto it = std::upper_bound(intervals.begin(), intervals.end(), begin,
        [](size_t position, const auto& interval) { return position < interval.begin; });
    if (it == intervals.begin()) return false;
    --it;
    return begin >= it->begin && end <= it->end;
}
struct RouteInsertDensity {
    std::vector<double> probability;
    std::vector<double> opportunity;
    std::vector<double> localEffectiveObservations;
    uint64_t observations = 0;
    double effectiveObservations = 0;
    int bandwidth = 0;
    bool usable = false;

    double at(int length) const {
        if (length < 0 || size_t(length) >= probability.size()) return 0;
        return probability[size_t(length)];
    }
};

// H(L)/sum_u coverage(u)*(length(u)-L+1) corrects the length-dependent
// opportunity to observe a fully contained fragment on the training nodes.
// This assumes uniform fragment starts and representative unique-anchor mapping.
inline RouteInsertDensity fitRouteInsertDensity(const std::vector<uint64_t>& histogram,
                                                std::vector<RouteTrainingSequence> population,
                                                int minimum, int maximum) {
    RouteInsertDensity result;
    if (minimum < 0 || maximum < minimum || population.empty()) return result;
    population.erase(std::remove_if(population.begin(), population.end(), [](const auto& u) {
        return !u.length || !std::isfinite(u.coverage) || u.coverage <= 0;
    }), population.end());
    if (population.empty()) return result;
    std::sort(population.begin(), population.end(), [](const auto& a, const auto& b) {
        return a.length < b.length;
    });
    maximum = std::min<int64_t>(maximum, population.back().length);
    if (maximum < minimum) return result;
    result.probability.assign(size_t(maximum) + 1, 0);
    result.opportunity.assign(size_t(maximum) + 1, 0);
    result.localEffectiveObservations.assign(size_t(maximum) + 1, 0);
    std::vector<double> mass(size_t(maximum) + 1, 0);
    long double sumCoverage = 0, sumCoverageLength = 0;
    for (const auto& u : population) {
        sumCoverage += u.coverage;
        sumCoverageLength += static_cast<long double>(u.coverage) * (u.length + 1);
    }
    size_t expired = 0, exposed = 0;
    long double total = 0, squaredWeight = 0;
    for (int length = 0; length <= maximum; ++length) {
        while (expired < population.size() && population[expired].length < size_t(length)) {
            const auto& u = population[expired++];
            sumCoverage -= u.coverage;
            sumCoverageLength -= static_cast<long double>(u.coverage) * (u.length + 1);
        }
        const double opportunity = std::max(0.0, double(sumCoverageLength - length * sumCoverage));
        if (length < minimum || opportunity <= 0) continue;
        result.opportunity[size_t(length)] = opportunity;
        ++exposed;
        const uint64_t count = size_t(length) < histogram.size() ? histogram[size_t(length)] : 0;
        result.observations += count;
        mass[size_t(length)] = count / opportunity;
        total += mass[size_t(length)];
        squaredWeight += static_cast<long double>(count) / (opportunity * opportunity);
    }
    if (!exposed || total <= 0 || squaredWeight <= 0) return result;
    result.effectiveObservations = double(total * total / squaredWeight);
    // Same sample requirement as the pre-existing insert fit, applied also
    // after inverse-opportunity weighting rather than treating weights as reads.
    if (result.observations < 1000 || result.effectiveObservations < 1000) return result;
    double mean = 0, variance = 0, cumulative = 0;
    int q1 = -1, q3 = -1;
    for (int length = minimum; length <= maximum; ++length) {
        mass[size_t(length)] /= double(total);
        mean += length * mass[size_t(length)];
        cumulative += mass[size_t(length)];
        if (q1 < 0 && cumulative >= .25) q1 = length;
        if (q3 < 0 && cumulative >= .75) q3 = length;
    }
    for (int length = minimum; length <= maximum; ++length)
        variance += (length - mean) * (length - mean) * mass[size_t(length)];
    const double spread = std::min(std::sqrt(variance), std::max(0, q3 - q1) / 1.34);
    result.bandwidth = std::max(1, int(std::ceil(.9 * spread * std::pow(result.effectiveObservations, -.2))));
    const int h = result.bandwidth;
    // Local training support uses exactly the density kernel and inverse
    // opportunity weights, excluding the prior. A huge global sample cannot
    // make one isolated tail observation an adequately trained density peak.
    std::vector<long double> localWeight(mass.size(), 0), localSquared(mass.size(), 0);
    for (int observed = minimum; observed <= maximum; ++observed) {
        const uint64_t count = size_t(observed) < histogram.size() ? histogram[size_t(observed)] : 0;
        if (!count || result.opportunity[size_t(observed)] <= 0) continue;
        const long double atom = 1.0L / result.opportunity[size_t(observed)];
        for (int query = std::max(minimum, observed - h); query <= std::min(maximum, observed + h); ++query) {
            const long double weight = atom * (h + 1 - std::abs(query - observed));
            localWeight[size_t(query)] += count * weight;
            localSquared[size_t(query)] += count * weight * weight;
        }
    }
    for (size_t i = 0; i < mass.size(); ++i) if (localSquared[i] > 0)
        result.localEffectiveObservations[i] = double(localWeight[i] * localWeight[i] / localSquared[i]);
    // Linear-time triangular convolution from cumulative mass and first moment.
    std::vector<double> prefix(mass.size() + 1, 0), moment(mass.size() + 1, 0);
    for (size_t i = 0; i < mass.size(); ++i) {
        prefix[i + 1] = prefix[i] + mass[i];
        moment[i + 1] = moment[i] + i * mass[i];
    }
    for (int length = minimum; length <= maximum; ++length) {
        if (result.opportunity[size_t(length)] <= 0) continue;
        const size_t low = size_t(std::max(minimum, length - h));
        const size_t mid = size_t(length), high = size_t(std::min(maximum, length + h) + 1);
        const double left = (h + 1 - length) * (prefix[mid + 1] - prefix[low]) +
                            (moment[mid + 1] - moment[low]);
        const double right = (h + 1 + length) * (prefix[high] - prefix[mid + 1]) -
                             (moment[high] - moment[mid + 1]);
        // Normalize over the actually exposed part of the kernel. A flat
        // corrected density must remain flat at the fit/exposure boundaries.
        const double nl = double(mid - low), nr = double(high - mid - 1);
        const double kernelSum = (h + 1) * (nl + nr + 1) - nl * (nl + 1) / 2 - nr * (nr + 1) / 2;
        const double smooth = std::max(0.0, (left + right) / kernelSum);
        // One uncertain effective fragment, spread uniformly over exposed bins:
        // finite holes/tails, no infinite likelihood or invented hard cutoff.
        result.probability[mid] = (result.effectiveObservations * smooth + 1.0 / exposed) /
                                  (result.effectiveObservations + 1.0);
    }
    result.usable = true;
    return result;
}

// Every entry is one physical pair. Each contributes at most one unit across
// DISTINCT candidate lengths, never one whole vote per alternative route.
inline std::vector<double> routeDistanceAllocation(const RouteInsertDensity& density,
                                                   const std::vector<int>& spans,
                                                   const std::vector<int>& lengths,
                                                   std::vector<double>* contrast = nullptr,
                                                   double minimumLocalTraining = 0) {
    std::vector<double> totals(lengths.size(), 0), weights(lengths.size(), 0);
    if (contrast) contrast->assign(lengths.size(), 0);
    if (!density.usable || lengths.size() < 2) return totals;
    for (int span : spans) {
        double sum = 0;
        bool comparable = true;
        for (size_t j = 0; j < lengths.size(); ++j) {
            const int64_t implied = int64_t(span) + lengths[j];
            if (implied < 0 || implied >= int64_t(density.probability.size()) ||
                density.opportunity[size_t(implied)] <= 0) { comparable = false; break; }
            weights[j] = density.at(int(implied));
            sum += weights[j];
        }
        // No model comparison when any alternative lacks training opportunity
        // or leaves the existing plausible window. Legacy evidence still acts.
        if (!comparable || sum <= 0) continue;
        const size_t favored = size_t(std::max_element(weights.begin(), weights.end()) - weights.begin());
        const size_t favoredLength = size_t(int64_t(span) + lengths[favored]);
        if (minimumLocalTraining > 0 &&
            (favoredLength >= density.localEffectiveObservations.size() ||
             density.localEffectiveObservations[favoredLength] < minimumLocalTraining)) continue;
        for (size_t j = 0; j < lengths.size(); ++j) totals[j] += weights[j] / sum;
    }
    if (contrast) for (size_t j = 0; j < totals.size(); ++j) {
        double rival = 0;
        for (size_t other = 0; other < totals.size(); ++other)
            if (other != j) rival = std::max(rival, totals[other]);
        // Net contrast must distinguish this route from EVERY alternative.
        // A/B-common support cannot become distinguishing just because C is
        // weak: the strongest competing total, not the weakest, is subtracted.
        (*contrast)[j] = std::max(0.0, totals[j] - rival);
    }
    return totals;
}

inline int decisiveRouteLength(const std::vector<double>& totals, double tieRatio) {
    if (totals.size() < 2) return -1;
    size_t best = 0;
    double second = -1;
    for (size_t j = 1; j < totals.size(); ++j) {
        if (totals[j] > totals[best]) { second = totals[best]; best = j; }
        else second = std::max(second, totals[j]);
    }
    if (totals[best] <= 0 || totals[best] <= second || totals[best] < tieRatio * second) return -1;
    return int(best);
}
inline int supportedRouteLength(const std::vector<double>& totals,
                                const std::vector<double>& contrast,
                                double tieRatio, double linkBar) {
    const int winner = decisiveRouteLength(totals, tieRatio);
    if (winner < 0 || contrast.size() != totals.size() || contrast[size_t(winner)] < linkBar) return -1;
    return winner;
}
}} // namespace ts::detail
