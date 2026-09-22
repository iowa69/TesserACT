// Default-off original-call consensus experiment. These are uncalibrated
// likelihood scores: neither FASTQ Q nor unique seeds prove genomic origin.
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <limits>

namespace ts { namespace quality_consensus {

inline bool enabled() {
    const char* value = std::getenv("TESSERACT_POLISH_ORIGINAL_QUALITY");
    return value && std::strcmp(value, "1") == 0;
}

constexpr unsigned kMaxPhred = 40;  // Fixed before experiments, not calibrated.
constexpr uint64_t kLogScale = 1ULL << 20;
constexpr double kScoreScale = 2.0 * kLogScale;

// Constants common to all four alleles cancel. Store the nonnegative log odds
// for the observed base, rounding once before summation. Q0 (and below-chance
// Q1) is uniform evidence, rather than evidence against the called base.
inline uint64_t logOdds(unsigned q) {
    static const std::array<uint64_t, kMaxPhred + 1> table = [] {
        std::array<uint64_t, kMaxPhred + 1> t{};
        for (unsigned i = 1; i <= kMaxPhred; ++i) {
            const double error = std::min(0.75, std::pow(10.0, -double(i) / 10.0));
            t[i] = static_cast<uint64_t>(std::llround(kLogScale * std::log(3.0 * (1.0 - error) / error)));
        }
        return t;
    }();
    return table[std::min(q, kMaxPhred)];
}

struct Observation {
    uint8_t base = 0;        // original allele in contig-forward orientation
    uint8_t phred = 0;       // original numerical Q, before capping
    uint8_t orientation = 0; // native read placement, not target strand
};

struct Pileup {
    std::array<uint64_t, 4> scores{};
    uint64_t fragments = 0;  // distinct fragments with nonuniform evidence
    uint64_t observations = 0;
    std::array<uint8_t, 4> orientations{};
    bool overflow = false;
};

enum class AddResult { Informative, Uniform, Invalid, Overflow };

// Call exactly once for a physical fragment at a single contig copy/site.
// Single observations receive 2*w; two mates each receive w. Dividing scores
// by 2*kLogScale therefore gives the mean of the fragment's log likelihoods.
// A Q0 mate is retained in the mean, but cannot supply depth or strand support.
inline AddResult addFragment(Pileup& p, const Observation* observations, unsigned count) {
    if (!count || count > 2 || !observations) return AddResult::Invalid;
    if (p.overflow) return AddResult::Overflow;
    std::array<uint64_t, 4> increment{};
    std::array<uint8_t, 4> strands{};
    bool informative = false;
    for (unsigned i = 0; i < count; ++i) {
        const auto& o = observations[i];
        if (o.base > 3 || o.phred > 93 || o.orientation > 1) return AddResult::Invalid;
        const uint64_t weight = logOdds(o.phred);
        increment[o.base] += weight * (count == 1 ? 2 : 1);
        if (weight) {
            informative = true;
            strands[o.base] |= static_cast<uint8_t>(1U << o.orientation);
        }
    }
    const uint64_t maximum = std::numeric_limits<uint64_t>::max();
    if (p.observations > maximum - count || (informative && p.fragments == maximum)) {
        p.overflow = true;
        return AddResult::Overflow;
    }
    for (unsigned b = 0; b < 4; ++b) {
        if (p.scores[b] > maximum - increment[b]) {
            p.overflow = true;
            return AddResult::Overflow;
        }
    }
    p.observations += count;
    if (informative) ++p.fragments;
    for (unsigned b = 0; b < 4; ++b) {
        p.scores[b] += increment[b];
        p.orientations[b] |= strands[b];
    }
    return informative ? AddResult::Informative : AddResult::Uniform;
}

enum class DecisionReason { Accept, CurrentAllele, BelowDepth, OneOrientation,
                            BelowPosterior, Tie, Overflow, InvalidThreshold };
struct Decision {
    int best = -1;
    double posterior = 0; // normalized likelihood under a uniform allele prior
    DecisionReason reason = DecisionReason::Tie;
};

inline Decision decide(const Pileup& p, int current, int minFragments, double minPosterior) {
    Decision d;
    if (p.overflow) { d.reason = DecisionReason::Overflow; return d; }
    if (minFragments < 1 || !std::isfinite(minPosterior) || minPosterior < 0 || minPosterior > 1) {
        d.reason = DecisionReason::InvalidThreshold; return d;
    }
    int best = 0;
    for (int b = 1; b < 4; ++b) if (p.scores[b] > p.scores[best]) best = b;
    unsigned ties = 0;
    double denominator = 0;
    for (int b = 0; b < 4; ++b) {
        ties += p.scores[b] == p.scores[best];
        denominator += std::exp(-double(p.scores[best] - p.scores[b]) / kScoreScale);
    }
    d.posterior = 1.0 / denominator;
    if (ties > 1) { d.reason = DecisionReason::Tie; return d; }
    d.best = best;
    if (best == current) d.reason = DecisionReason::CurrentAllele;
    else if (p.fragments < static_cast<uint64_t>(minFragments)) d.reason = DecisionReason::BelowDepth;
    else if (p.orientations[best] != 3) d.reason = DecisionReason::OneOrientation;
    else if (d.posterior < minPosterior) d.reason = DecisionReason::BelowPosterior;
    else d.reason = DecisionReason::Accept;
    return d;
}

} } // namespace ts::quality_consensus
