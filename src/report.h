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
#include "organism_layout.h"
#include "replicon.h"
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
    size_t carryOverRescued = 0;   // k-mers admitted below the cutoff as trusted
    size_t gapsClosed = 0;         // dead-end pairs bridged by the gap closer
    size_t deadEndsBuilt = 0;      // dead ends in the raw graph, before any simplification
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
    size_t   ladderRescued = 0;        // contigs recovered from a smaller rung
    uint64_t ladderRescuedBases = 0;
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
    LayoutStats layout;
    RepliconAssignment replicons;
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

    // The contigs as written to contigs.fasta: split at every N, after the terminal-overlap trim
    // (and any split post-processing). Scaffolding raises n50/largest above these by asserting
    // an order across gaps; it adds no assembled sequence. Reporting only the scaffold figure is
    // how a layout change gets read as an assembly improvement, so both travel together
    // everywhere they are shown. (Release 1.3.0 counted the scaffolds split at runs of 10+ N,
    // before the trim; that described no written file -- G-emit T17.)
    size_t scaffoldGaps = 0;     // number of N-runs of any length (= scaffolds.agp N rows)
    size_t contigPieces = 0;     // contigs.fasta records
    size_t contigN50 = 0;
    size_t contigLargest = 0;
    size_t contigTotal = 0;      // bases in contigs.fasta (called bases; the records hold no N)
    size_t trimmedOverlaps = 0;      // terminal repeat overlaps trimmed from the records
    size_t trimmedOverlapBases = 0;  // bases they removed

    // Organism Model 2.0 (C3): the report.json "om2" object, rendered by om2_output.cpp. Empty
    // unless an om2 flag ran (TESSERACT_OM2_OUTPUT, TESSERACT_OM2_AGP_EVIDENCE, detection),
    // and then no key is written: a flags-off report is the release's byte for byte.
    std::string om2Json;

    // Phase 2 (EVAL_PLAN_P2 s5.2): the report.json "p2" object -- one block for every W1 feature,
    // rendered by p2::p2BlockJson (emit-B) with emit-A's r2/r3/f5 members (p2::reportJsonMembers).
    // Every assembler run sets it; with every TESSERACT_P2_* flag off it holds only enabled=0 and
    // zero counters. `p2` is what report.html shows of it.
    std::string p2Json;
    struct P2Html {
        bool selfQa = false;
        bool alarm = false;
        std::string verdict;
        double completeness = 0, qvRead0 = 0;
        uint64_t missingGt10x = 0, missing3to10 = 0, spikeSolid = 0, spikeRef = 0;
        bool ends = false;
        size_t endRows = 0, endsAudited = 0, tipEnds = 0, tipBases = 0, lowercasedBases = 0;
        bool detect = false;
        std::string detectCall, detectBest;
        double detectBestScore = 0;
        bool provenance = false;
    } p2;

    // Recomputes the derived summary fields from `contigs`.
    void finalize();
};

// Self-contained HTML: all CSS and JS inlined, no external requests.
bool writeHtmlReport(const std::string& path, const AssemblyReport& rep, std::string& error);
bool writeJsonReport(const std::string& path, const AssemblyReport& rep, std::string& error);

}  // namespace ts
