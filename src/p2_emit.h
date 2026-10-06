// Phase 2, wave W1, workstream emit-A (om2/design/EVAL_PLAN_P2.md; design inputs
// work/backward/FINAL.md R2, R3, F5 and work/backward/finishing/DESIGN.md). Every feature here
// is behind its own TESSERACT_P2_* flag, default off; with all of them unset the assembler's
// four output files are the release's byte for byte, report.json gains only the "p2" block
// (enabled=0, zero counters), and the counter lines print zeros.
//
//   R2  TESSERACT_P2_LIBGUARD=1        library orientation / insert-model guard (p2_libguard.h;
//                                      measured and applied inside the paired resolver).
//   R3  TESSERACT_P2_CIRC=close|verify exact circles. For every record the pair-based call
//                                      tagged `_circular`, the record's own graph walk is closed
//                                      the way backward/finishing/scripts/graph_close.py does
//                                      (self / via / repeat_end), giving the circle without the
//                                      duplicated (k-1) end or the missing connector bases.
//                                      `verify` also makes `_circular` a verified claim: the
//                                      circle is written circular only if a closure exists (G),
//                                      it is graph-isolated (I), reads span the closed join (J)
//                                      and few mates near the join leave the record (X) -- the
//                                      frozen rule of EVAL_PLAN_P2 s7.4. A candidate that fails
//                                      is written LINEAR and duplicate-free (owner decision D9).
//                                      `close` applies the closure but keeps the release tag rule
//                                      (EVAL_PLAN_P2 s12 row W1-4 fallback).
//   F5  TESSERACT_P2_SPIKEIN=1         PhiX174 spike-in screen: a record with >= 90 % of its
//                                      31-mers in the built-in phiX174 set is renamed with the
//                                      `_spikein` label in place of its replicon tag (never
//                                      `_plas`, never `_circular`). Name-only: no base changes.
//
// Side files, written only when R3 or F5 is on: replicons.tsv (one row per contigs.fasta record,
// with the R3 evidence) and p2_edits.tsv (every record edit, EVAL_PLAN_P2 s5.2).
#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "emit_post.h"
#include "gfa.h"
#include "graph.h"
#include "p2_libguard.h"
#include "seqio.h"

namespace ts {
namespace p2 {

enum class CircMode : uint8_t { Off = 0, Close = 1, Verify = 2 };

struct EmitConfig {
    bool libGuard = false;          // TESSERACT_P2_LIBGUARD
    CircMode circ = CircMode::Off;  // TESSERACT_P2_CIRC
    bool spikein = false;           // TESSERACT_P2_SPIKEIN
    bool recordStage() const { return circ != CircMode::Off || spikein; }
    bool any() const { return libGuard || recordStage(); }
};
// Read per call from the environment (never cached).
EmitConfig readConfig();
const char* circModeName(CircMode m);

// ---- R3 frozen thresholds (EVAL_PLAN_P2 s7.4, registered with the binary) ----------------------
struct VerifyThresholds {
    size_t minSpan = 3;          // J: spanning reads >= 3
    double minSpanCtrl = 0.3;    // J: spanning >= 0.3 x the control position's spanning count
    double maxClipFrac = 0.10;   // J: clip_frac <= 0.10
    double maxExitRatio = 0.02;  // X: exit_ratio <= 0.02
};
constexpr int kJoinAnchor = 20;      // a spanning read covers [J-20, J+20) (junction_verify.py --anchor)
constexpr int kClipWindow = 3;       // a clip point within +-3 bp of J
constexpr int kMinClip = 5;          // a soft clip of >= 5 bases
constexpr int kAnchorK = 31;         // read placement k (unique over every written record)
constexpr double kSpikeinMinFrac = 0.90;   // F5: >= 90 % of a record's 31-mers in phiX174

// ---- R3 graph closure --------------------------------------------------------------------------
enum class Closure : uint8_t { NoPath, EndsChanged, RepeatEnd, Self, Via, None, Gapped };
const char* closureName(Closure c);
inline bool isClosure(Closure c) { return c == Closure::RepeatEnd || c == Closure::Self || c == Closure::Via; }

struct ClosureResult {
    Closure mode = Closure::NoPath;
    bool closed = false;         // `seq` holds the closed circle
    std::string seq;
    size_t pathN = 0;
    int outLast = -1, inFirst = -1, isolated = -1;   // -1: no walk, not computed
    size_t alts = 0;             // distinct connectors closing the walk (a bubble at the join)
    std::string viaSeg;          // chosen connector, "<segment><+|->"
    size_t viaLen = 0;
    double viaDp = 0;
    std::string altDp;           // the other connectors, "<seg><sign>:<dp>" comma-separated
    long interior = 0;           // connector length - 2(k-1); negative = the ends overlap
    size_t viaOtherLinks = 0;
    std::string note;
};
// `s` is the record as it would be written, `walk` its gap-free walk ((unitig << 1) | orient,
// spelling `s` exactly), `g` the final graph. Mirrors graph_close.py decision for decision.
ClosureResult closeWalk(const std::string& s, const std::vector<uint64_t>& walk, const UnitigGraph& g);

// Largest o in [lo, min(hi, |s| - 1)] with prefix(o) == suffix(o); 0 if none (circ_census.py).
size_t terminalSelfOverlap(const std::string& s, size_t lo = 20, size_t hi = 3000);
// Largest (k - 1), k in `ladder`, with prefix(k-1) == suffix(k-1) and |s| > 2(k-1); 0 if none.
size_t kMinusOneSelfOverlap(const std::string& s, const std::vector<int>& ladder);

// ---- R3 junction verification ------------------------------------------------------------------
struct JunctionStats {
    bool computed = false;
    size_t circleLen = 0;        // length of the tested circle T
    size_t selfOvTrimmed = 0;    // exact self-overlap removed to form T when no closure exists
    size_t join = 0, ctrl = 0;   // positions on the half-rotated T
    size_t spanReads = 0, clipReads = 0, ctrlSpan = 0, ctrlClip = 0;
    double clipFrac = 0;
    size_t reads = 0, mateElsewhere = 0, nearJoin = 0, exitsNearJoin = 0;
    double exitRatio = 0;
};
struct VerifyBatch {
    // Every written record (as it stands before R3 edits).
    const std::vector<std::string>* records = nullptr;
    // Candidate c tests circle candCircle[c] in place of records[candRecord[c]].
    std::vector<size_t> candRecord;
    std::vector<std::string> candCircle;
};
struct VerifyRunStats {
    size_t readsPlaced = 0;
    size_t indexKmers = 0;
    int windowLo = 0, windowHi = 0;   // library window (p0.5 / p99.5 of FR template lengths)
    size_t windowPairs = 0;
};
// One placement pass of every read over the records, each candidate replaced by its
// half-rotated circle (junction_verify.py with minimap2 replaced by unique 31-mer placement and a
// gapless local score, match +2 / mismatch -8, standing in for the aligner's soft clipping).
std::vector<JunctionStats> verifyJunctions(const VerifyBatch& batch, const SequenceStore& reads,
                                           int threads, VerifyRunStats& run);
bool passJ(const JunctionStats& j, const VerifyThresholds& t);
bool passX(const JunctionStats& j, const VerifyThresholds& t);

// ---- F5 ----------------------------------------------------------------------------------------
const std::string& phix174();                 // NC_001422.1 (p2_phix.cpp)
// Fraction of the record's valid 31-mers that occur (either strand) in circular phiX174.
double phixKmerFraction(const std::string& s);

// ---- the record stage (runs on the written records, after the split and the terminal trim) -----
struct RecordRow {
    size_t piece = 0;
    std::string oldTag;
    bool spikein = false;
    double phixFrac = -1;        // -1: not computed
    bool candidate = false;      // carried `_circular` from the pair call
    ClosureResult closure;
    JunctionStats junction;
    bool passG = false, passI = false, passJ = false, passX = false;
    std::string verdict;         // verified | failed_<first missing part> | close | spikein | -
    std::string topology;        // circular | linear
    size_t lenBefore = 0, lenAfter = 0;
    size_t trimBp = 0, insertBp = 0, kMinusOneTrim = 0;
    std::string edit;            // none | closed | trimmed_duplicate | kminus1_trim
};
struct EditRow {
    size_t piece = 0;
    std::string feature;         // R3 | F5
    std::string operation;
    size_t start = 0, end = 0;   // 0-based, half-open, in the record before the edit
    size_t bases = 0;            // bases removed or inserted (0 for a rename)
};
struct RecordStageStats {
    size_t screened = 0, spikein = 0, spikeinBp = 0, spikeinCircular = 0;
    size_t candidates = 0;
    size_t byMode[7] = {0, 0, 0, 0, 0, 0, 0};    // indexed by Closure
    size_t isolated = 0, passJ = 0, passX = 0, verified = 0;
    size_t writtenCircular = 0, writtenLinear = 0;
    size_t closed = 0, trimRecords = 0, trimBp = 0, insertRecords = 0, insertBp = 0, kMinusOneTrims = 0;
    size_t scaffoldClaimsDropped = 0;   // verify: `_circular` removed from a gapped scaffold's name
    VerifyRunStats verify;
};
struct RecordStageIn {
    EmitConfig cfg;
    std::vector<EmitPiece>* pieces = nullptr;            // .drop != 0: not written
    std::vector<std::string>* scaffolds = nullptr;       // the scaffold records (outSeqs)
    std::vector<std::string>* scaffoldTags = nullptr;    // their replicon tags (outTags)
    const std::vector<GfaPath>* scaffoldPaths = nullptr; // their walks (outPathOf)
    const UnitigGraph* graph = nullptr;
    const SequenceStore* reads = nullptr;
    std::vector<int> ladder;
    int threads = 1;
    VerifyThresholds thresholds;
};
struct RecordStageResult {
    std::vector<RecordRow> rows;                // one per written piece, in output order
    std::vector<EditRow> edits;
    std::vector<size_t> scaffoldsChanged;       // scaffold sequence or tag changed
    RecordStageStats st;
};
// Applies F5 then R3 to the pieces (sequence and tag) and, where a piece is a whole untrimmed
// scaffold, to the scaffold record and tag as well. Names are assigned afterwards by the caller
// from piece.seq / piece.tag, so every name carries the final length and label.
RecordStageResult applyRecordStage(RecordStageIn& in);

// Replicon tag with the `_circular` suffix removed or added.
std::string dropCircular(const std::string& tag);

// Side files (after the caller named the pieces).
bool writeRepliconsTsv(const std::string& path, const RecordStageResult& r,
                       const std::vector<EmitPiece>& pieces, std::string& error);
bool writeEditsTsv(const std::string& path, const RecordStageResult& r,
                   const std::vector<EmitPiece>& pieces, std::string& error);
// The rows of p2_edits.tsv for this stage's edits ("record\tfeature\toperation\tstart\tend\tbases",
// 0-based half-open coordinates in the record before the edit), in edit order.
std::vector<std::string> formatEditRows(const RecordStageResult& r, const std::vector<EmitPiece>& pieces);
// Writes p2_edits.tsv: the header line, then `rows` as given. The integration build writes every W1
// feature's rows through this one writer (EVAL_PLAN_P2 s5.2).
bool writeEditLog(const std::string& path, const std::vector<std::string>& rows, std::string& error);

// Counter lines, printed on every run (zeros when the flag is off).
std::string formatLibGuardCounters(const LibGuardStats& s);
std::string formatCircCounters(const EmitConfig& c, const RecordStageStats& s);
std::string formatSpikeinCounters(const EmitConfig& c, const RecordStageStats& s);
// The report.json "p2" object (raw JSON text). With every flag off it holds only enabled=0 and
// zero counters.
std::string reportJson(const EmitConfig& c, const LibGuardStats& lg, const RecordStageResult* r,
                       const std::vector<EmitPiece>* pieces);
// The same object's members without the braces and without "enabled" ("r2": {...}, "r3": {...},
// "f5": {...}): the integration build merges them into the single report.json "p2" block that
// p2::p2BlockJson (p2_emitb.h) renders, so there is one "p2" key with one "enabled" flag.
std::string reportJsonMembers(const EmitConfig& c, const LibGuardStats& lg, const RecordStageResult* r,
                              const std::vector<EmitPiece>* pieces);

}  // namespace p2
}  // namespace ts
