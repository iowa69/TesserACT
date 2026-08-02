// Everything the run observed about itself, collected as it goes so the report
// can explain not just the final assembly but how each stage got there.
#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "correct.h"
#include "gapfill.h"
#include "mappolish.h"
#include "organism_join.h"
#include "graph.h"
#include "polish.h"
#include "resolve.h"

namespace ts {

// One rung of the multi-k ladder.
struct KIteration {
    int k = 0;
    uint64_t totalKmers = 0;
    uint64_t distinctKmers = 0;
    uint64_t solidKmers = 0;
    uint32_t cutoff = 0;
    double peakCoverage = 0;
    // count -> distinct k-mers with that count, truncated for reporting.
    std::vector<uint64_t> countHistogram;

    size_t unitigsBuilt = 0;
    size_t lengthBuilt = 0;
    size_t n50Built = 0;
    size_t unitigsFinal = 0;
    size_t lengthFinal = 0;
    size_t n50Final = 0;
    double medianCoverage = 0;

    std::vector<SimplifyRoundStats> rounds;
    double countSeconds = 0;
    double graphSeconds = 0;
    double simplifySeconds = 0;
    size_t carryOverContigs = 0;
};

struct ContigRecord {
    size_t length = 0;
    double coverage = 0;
    double gcPercent = 0;
    size_t gapBases = 0;   // N runs introduced by scaffolding
};

// The whole run, in the order it happened.
struct AssemblyReport {
    std::string version;
    std::string command;
    std::string startedAt;
    std::string mode;
    int threads = 0;
    double totalSeconds = 0;
    double peakMemoryBytes = 0;

    // input
    size_t reads = 0;
    uint64_t inputBases = 0;
    uint32_t maxReadLength = 0;
    uint64_t qualityTrimmedBases = 0;
    bool paired = false;
    std::vector<std::string> inputFiles;

    // stages
    bool correctionRun = false;
    CorrectionStats correction;
    double correctionSeconds = 0;
    int correctionK = 0;

    std::vector<KIteration> iterations;

    bool resolveRun = false;
    ResolveStats resolve;
    double resolveSeconds = 0;
    std::vector<uint64_t> insertHistogram;   // index = fragment length

    bool gapFillRun = false;
    GapFillStats gapFill;

    // The organism model stage: which panel it consulted, what it joined,
    // and which accessions the model was built without.
    MapPolishStats mapPolish;
    bool organismRun = false;
    OrganismJoinStats organism;
    std::string organismName;
    uint32_t organismGenomes = 0;
    uint32_t organismPlasmids = 0;
    std::vector<std::string> organismExcluded;

    bool polishRun = false;
    PolishStats polish;
    double polishSeconds = 0;

    // output
    std::vector<ContigRecord> contigs;   // sorted longest first
    size_t totalLength = 0;
    size_t n50 = 0, n90 = 0, l50 = 0;
    size_t largest = 0;
    double gcPercent = 0;
    double meanCoverage = 0;
    size_t gapBases = 0;
    size_t gfaSegments = 0;
    size_t gfaLinks = 0;

    // Recomputes the derived summary fields from `contigs`.
    void finalize();
};

// Self-contained HTML: all CSS and JS inlined, no external requests.
bool writeHtmlReport(const std::string& path, const AssemblyReport& rep, std::string& error);
bool writeJsonReport(const std::string& path, const AssemblyReport& rep, std::string& error);

}  // namespace ts
