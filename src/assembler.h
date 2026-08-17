// Assembly orchestration: the multi-k loop and the stages that run around it.
#pragma once

#include <string>
#include <vector>

#include "graph.h"
#include "libqc.h"
#include "mappolish.h"
#include "organism_join.h"
#include "organism_layout.h"
#include "organism.h"
#include "report.h"
#include "seqio.h"

namespace ts {

// Presets trading speed against contiguity and caution. A mode only fills in
// knobs the user left alone; an explicit flag always wins.
enum class RunMode { Fast, Standard, Careful, Aggressive };

const char* runModeName(RunMode m);
bool parseRunMode(const std::string& s, RunMode& out);

struct AssemblyOptions {
    std::vector<Library> libraries;
    std::string outDir = "tesseract_out";
    std::vector<int> kValues;        // empty means "choose from mode and read length"
    int threads = 0;                 // 0 = hardware concurrency
    uint32_t forcedCutoff = 0;       // 0 = auto
    size_t minContigLen = 0;         // 0 = 2*maxK
    bool verbose = true;
    bool correctReads = true;
    bool resolveRepeats = true;
    bool scaffold = true;
    bool gapFill = true;
    bool polish = true;
    bool emitGfa = true;
    bool emitHtml = true;
    bool emitUnitigs = false;

    RunMode mode = RunMode::Standard;
    int minLinkSupport = 2;          // floor on paired reads needed to trust a join
    // Paired support required per unit of median unitig coverage before a
    // contested join is taken. A flat count cannot mean the same thing at every
    // depth, and the trade it controls is contiguity against chimeric joins --
    // measured across 20 closed-reference isolates:
    //     flat 2   median NGA50 303,624   12 misassemblies
    //     0.07                  291,265    7
    //     0.10                  277,799    3
    //
    // Retuned 0.10 -> 0.02 (with tieRatio 1.15 -> 1.02) against 50 closed-reference
    // strains chosen for being ones SPAdes wins. The 20-isolate measurement above is not
    // contradicted -- it swept ONE knob on a panel that was not selected for difficulty,
    // and it was right that loosening this alone costs misassemblies. What changed is that
    // the cost is now paid for elsewhere: the carry-weight boost that this release also
    // disables was itself generating misassemblies, and removing it buys back more than
    // the looser threshold spends.
    //
    // Measured together on those 50 strains, against SPAdes:
    //     setting                       NGA50 W/L   median deficit   misassemblies W/L
    //     0.10 / 1.15 (previous)          4/46        -25,414           16/13
    //     0.05 / 1.05                    14/36        -20,901           15/16
    //     0.02 / 1.02 (this)             18/32        -10,072           15/16  (p=0.91)
    //     0.00 / 1.00                    19/31        -10,072           14/19
    //
    // 0.00/1.00 buys one more contiguity win and gives back misassemblies for it, so the
    // knob is set one notch short of its limit rather than at it.
    double linkSupportPerX = 0.02;
    double tieRatio = 1.02;          // winning branch must beat the runner-up by this
    double bubbleCoverageLimit = 0.35;   // see UnitigGraph::simplify
    int simplifyRounds = 12;
    int polishPasses = 1;
    long long maxMemoryBytes = 0;   // 0 = 80% of physical RAM
    QualityTrim qtrim;              // 3' trimming applied as reads are loaded

    // The command line as invoked, recorded into the reports so a run can be
    // explained and reproduced from its own output.
    std::string commandLine;

    // A prior learned from closed genomes of the same organism, used only at
    // junctions no fragment can span. Empty means the assembly stops there
    // rather than guessing, which is the default.
    std::string organism;            // informational; names the expected model
    std::string organismModelPath;
    // Known insertion sequences. Contig ends inside one are left unjoined,
    // because that is where this isolate differs from the panel by construction.
    std::string isPanelPath;
    // Recurrent insertion sites: flank signatures that let an element be
    // placed on the chromosome rather than merely avoided.
    std::string isSitesPath;

    // A QC report over the same reads, produced by scepter (formerly fastplus). It
    // carries measurements the assembler would otherwise have to infer from its own
    // k-mer histograms -- read depth, genome size, and an error rate observed by
    // comparing overlapping mates -- and it carries them BEFORE the first rung is
    // counted, which is the only time several of the decisions below can still be made.
    // Empty means every such decision falls back to the self-derived estimate.
    std::string qcPath;
    // Lay the finished contigs out against the panel chromosome they most resemble, when
    // the model carries layout tracks. On by default because it is what closes the
    // chromosome; --no-layout leaves the assembly de novo in its ordering as well as its
    // sequence, which is the right choice when the ordering must not rest on a relative.
    bool layout = true;

    // Polishing against a full read alignment, after the k-mer polisher has
    // done what it can. Off unless a mapper is named, since it shells out.
    Mapper mapPolisher = Mapper::None;
    std::string mapperDir;

    // Set when the user named the knob explicitly, so applyMode() leaves it be.
    bool userSetK = false, userSetMinLink = false, userSetTie = false, userSetLinkPerX = false,
         userSetBubble = false, userSetRounds = false, userSetPolishPasses = false;

    void applyMode();
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
    const AssemblyReport& report() const { return report_; }

private:
    // Picks the k ladder from the mode and the observed read length, then lets the QC
    // report veto rungs the library is too thin to support.
    std::vector<int> resolveKLadder() const;
    std::vector<int> baseKLadder() const;
    std::vector<int> trimLadderToCoverage(std::vector<int> ladder) const;

    // One de Bruijn iteration: count, filter, build, simplify.
    bool iterate(int k, const std::vector<std::string>& carryOver, UnitigGraph& graph,
                 double& meanCoverage, KIteration& it, std::string& error,
                 double prevPeak = 0.0);

    // Coverage peak of the first k (no carry-over) and that k, used to project how
    // thin the later rungs will be without letting the boost feed back on itself.
    double basePeak_ = 0.0;
    int baseK_ = 0;
    int finalK_ = 0;   // largest k in the ladder, set before the loop

    AssemblyOptions opt_;
    LibraryQC qc_;
    SequenceStore reads_;
    OrganismModel organismModel_;
    IsPanel isPanel_;
    AssemblyStats stats_;
    AssemblyReport report_;
};

}  // namespace ts
