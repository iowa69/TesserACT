// Assembly orchestration: the multi-k loop and the stages that run around it.
#pragma once

#include <string>
#include <vector>

#include "graph.h"
#include "seqio.h"

namespace ts {

struct AssemblyOptions {
    std::vector<Library> libraries;
    std::string outDir = "tessera_out";
    std::vector<int> kValues;        // empty means "choose from read length"
    int threads = 0;                 // 0 = hardware concurrency
    uint32_t forcedCutoff = 0;       // 0 = auto
    size_t minContigLen = 0;         // 0 = 2*maxK
    bool verbose = true;
    bool correctReads = true;
    bool resolveRepeats = true;
    bool scaffold = true;
    bool polish = true;
    int minLinkSupport = 2;          // paired reads needed to trust a join
    double tieRatio = 1.15;           // winning branch must beat the runner-up by this
    double bubbleCoverageLimit = 0.35;   // see UnitigGraph::simplify
};

struct AssemblyStats {
    size_t contigs = 0;
    size_t totalLength = 0;
    size_t n50 = 0;
    size_t largest = 0;
    double meanCoverage = 0;
    double gcPercent = 0;
    double seconds = 0;
    std::vector<int> kUsed;
};

class Assembler {
public:
    explicit Assembler(AssemblyOptions opt);

    bool run(std::string& error);
    const AssemblyStats& stats() const { return stats_; }

private:
    // Picks the k ladder from the observed read length when none was given.
    std::vector<int> resolveKLadder() const;

    // One de Bruijn iteration: count, filter, build, simplify.
    bool iterate(int k, const std::vector<std::string>& carryOver,
                 UnitigGraph& graph, double& meanCoverage, std::string& error);

    AssemblyOptions opt_;
    SequenceStore reads_;
    AssemblyStats stats_;
};

// Estimates the fragment-length distribution of a paired library by mapping a
// sample of pairs onto long contigs.
bool estimateInsertSize(const SequenceStore& reads, const std::vector<std::string>& contigs,
                        int k, double& mean, double& stddev, size_t& mapped);

}  // namespace ts
