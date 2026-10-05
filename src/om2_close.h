// Organism Model 2.0, component C2 (CLOSE): stage [4c/7], between the model layout and scaffold
// gap closing. Default OFF (TESSERACT_OM2_CLOSE); with it unset the stage returns without touching the
// assembly and prints only its counter line, so every output file is unchanged.
//
// What it does when on, for every N-run the earlier stages wrote (and, with TESSERACT_OM2_CIRC,
// for the wrap of every record):
//
//   1. anchors both flanks in the final-k graph on 31-mers that occur once in the graph and lie in
//      single-copy unitigs (copy number from depth against the component's median);
//   2. finds every walk between the anchors up to 12 kb (bounded Dijkstra both ways, then the
//      sub-graph of nodes on some walk): its cut nodes and its bubble sites, whether it is
//      cyclic, whether it uses a unitig in both orientations (hairpin), and whether it needs a
//      single-copy unitig that is already placed elsewhere (double use: never bridged);
//   3. counts, in one pass over the reads, the fragments (single reads = threads, read pairs)
//      that link each junction's unique flanks to each other and to each bubble branch;
//   4. bridges junctions in order forced > pairs/threads > model tie-break > open, with a
//      budget (single-copy unitigs once, hard; repeat unitigs at most `hi` times, soft);
//   5. allocates one branch per bubble per locus with a basis label (om2_alloc.h);
//   6. writes a fill into the sequence (contig grade: all views) only when the walk is exhaustive,
//      unique or phased by pairs/threads at every site, hairpin-free, length-consistent, and the
//      adjacency is forced or spanned by fragments; every other admitted fill is genome-only and
//      stays in the ledger and in <out>/om2_close/ for the owner view;
//   7. closes a circle only on an exhaustive unique wrap walk with no pair contradiction, and then
//      records the dnaA rotation offset (TESSERACT_OM2_DNAA) -- never rotating an open record.
//
// Every base written comes from this isolate's graph. The sidecars (TESSERACT_OM2_RRN_PRIOR,
// TESSERACT_OM2_DNAA) carry hashes and counts only and are pinned to the loaded model's md5.
#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "om2_alloc.h"
#include "om2_types.h"

namespace ts {
class UnitigGraph;
class SequenceStore;
struct InsertModel;

namespace om2 {

enum class FillMode : uint8_t { None, Genome, Contig };

struct CloseOptions {
    bool enabled = false;              // TESSERACT_OM2_CLOSE
    FillMode fill = FillMode::Genome;  // TESSERACT_OM2_FILL=none|genome|contig
    AllocMode alloc = AllocMode::Phased;  // TESSERACT_OM2_ALLOC=off|consensus|phased|prior
    bool circ = false;                 // TESSERACT_OM2_CIRC
    bool flow = false;                 // TESSERACT_OM2_FLOW
    std::string rrnPrior;              // TESSERACT_OM2_RRN_PRIOR
    std::string dnaa;                  // TESSERACT_OM2_DNAA
    bool anyFlagSet = false;           // some TESSERACT_OM2_* of this stage is set

    // Fixed by the design (DESIGN.md section 6); not flags.
    int anchorWindow = 2000;           // bp searched for a unique anchor on each side
    int maxFill = 12000;               // longest walk considered (repeat traversal cap)
    int maxStates = 20000;             // per-junction search budget before abstaining
    int maxBranches = 16;              // paths enumerated per bubble site
    int circMinLen = 1000;             // shortest record whose wrap is tested
    double singleCopy = 1.5;           // copy <= this: single-copy for anchoring
    double hardSingleRatio = 1.6;      // hard c = 1: length >= 1 kb and depth <= 1.6 theta
    size_t hardSingleLen = 1000;
    uint32_t minSpan = 2;              // fragments linking two flanks: pairs/threads bridge
    uint32_t minContra = 3;            // fragments naming another record end: contradiction
    double tieMargin = 2.0;            // model tie-break: best >= 2x runner-up ...
    uint32_t tieGenomes = 20;          // ... and >= 20 genomes
    AllocParams allocParams;

    static CloseOptions fromEnv();
};

struct CloseStats {
    bool enabled = false;
    size_t junctions = 0, wraps = 0, anchored = 0, unanchored = 0;
    size_t clusters = 0, bridges = 0, forced = 0, pairThread = 0, tiebreak = 0, tiebreakUnverified = 0;
    size_t open = 0, budgetRefused = 0, hairpinRefused = 0, cyclic = 0, doubleUse = 0, noWalk = 0;
    size_t flankMismatch = 0, overlap = 0, lengthInconsistent = 0;
    size_t fillsContig = 0, fillsGenome = 0, fillBpContig = 0, fillBpGenome = 0;
    AllocCounts alloc;
    size_t variantSites = 0;
    size_t circlesTested = 0, circlesClosed = 0, circlesContig = 0, circlesOverlap = 0, rotated = 0;
    size_t abstainTangled = 0;
    // integration: N-runs matched to a C1 row, bridges C1's verdict refused, bridges C1's graph
    // verdict admitted (no pairs, not forced), wraps left open (writtenN = -1 in the ledger)
    size_t c1Attached = 0, c1Refused = 0, c1Gated = 0, wrapsOpen = 0, variantRows = 0;
    bool priorLoaded = false, dnaaLoaded = false;
    double seconds = 0;
};

struct CloseInputs {
    const UnitigGraph& graph;
    const SequenceStore& reads;
    const InsertModel& insert;
    std::string outDir;                // "" = write no side files
    std::string modelPath;             // the loaded organism model, for the sidecar md5 pins
    int threads = 1;
    bool verbose = false;
    // Integration (build_om2), all optional; null reproduces the standalone component.
    //   c1      C1's junction ledger. An N-run C1 recorded is not recorded again: its C1 row takes
    //           C2's fill (in the row's own orientation), allocation labels, contig-grade flag and
    //           closeGaps decision, and C1's verdict gates the bridge (see om2_close.cpp step 6).
    //   c1Runs  C1's row behind every N-run of `seqs` (SeamContext::runOwners, one vector per record).
    //   variants  receives one RepeatVariant per (site, locus) of every genome-view fill, in the
    //           orientation of the fill's ledger row, for genome/repeat_variants.tsv.
    Ledger* c1 = nullptr;
    const std::vector<std::vector<RunOwner>>* c1Runs = nullptr;
    std::vector<RepeatVariant>* variants = nullptr;
};

// Runs the stage. `seqs` gains contig-grade fills (N-runs replaced, circles completed); nothing
// else in the assembler's state is touched. `ledger` receives one Junction per N-run (unless
// in.c1 already holds a row for it, which is then updated instead) and per wrap, in that order,
// appended after any junctions already in it. The new rows' ids continue after in.c1's rows, so
// C1's rows followed by `ledger`'s form one ledger whose index equals the id (C3 reads it so).
// A wrap that is not closed is recorded with writtenN = -1; a wrap closed by an overlap of the
// record's two ends has claimedN = -overlap, writtenN = 0 and no fill. ledger.rotateOffset[w] for
// the w-th wrap junction (w < 64), in the stage's coordinates of the closed record: >= 0 = dnaA
// starts there on the forward strand; <= -2 = dnaA starts at (-2 - value) on the reverse strand;
// -1 = not located, or the wrap is not closed.
CloseStats runCloseStage(const CloseInputs& in, const CloseOptions& opt, std::vector<std::string>& seqs,
                         Ledger& ledger);

// "[om2-close] ..." -- printed on every run, zeros when the stage is off.
std::string formatCloseCounters(const CloseStats& s);
void logCloseCounters(const CloseStats& s);

}  // namespace om2
}  // namespace ts
