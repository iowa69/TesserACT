// Organism Model 2.0, C1b-e -- the junction ledger and the seam gate.
//
// Flags (all unset by default; with every one unset `SeamContext::create` returns nullptr and
// no om2 code touches the assembly):
//
//   TESSERACT_OM2_EVIDENCE=1          build the evidence index and print its counters only
//   TESSERACT_OM2_SEAM=0|audit|act    audit: judge and record every junction, change nothing
//                                     (writes om2_seams.tsv); act: apply the verdicts
//   TESSERACT_OM2_UNSIZED=butt|break|sized   the 1-N butt, act mode (default break)
//   TESSERACT_OM2_CAP=keep|size       the layout's 2000-N cap, act mode (default size)
//   TESSERACT_OM2_PAIRS=0|1           read-pair seam test (default 1 when SEAM is set)
//   TESSERACT_OM2_PLASMID_RULE=0|1    plasmid joins need positive evidence, act (default 1)
//   TESSERACT_OM2_CLOSEGAPS_POLICY=0|1  the gap filler skips gaps the ledger does not allow,
//                                     act (default 1)
//   TESSERACT_OM2_ADMIT=<path>        class table (C3 calibration): p_misjoin per class, τ_S, τ_G
//   TESSERACT_OM2_TANGLE_MAX=<n>      abstain above n live graph nodes (default 5000)
//
// Counter lines, printed on every run (zeros when off): [om2-evidence] and [om2-seam].
#pragma once

#include <cstdint>
#include <cstdio>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "om2_evidence.h"
#include "om2_types.h"

namespace ts {
namespace om2 {

struct SeamConfig {
    enum class Mode { Off, Audit, Act };
    enum class Unsized { Butt, Break, Sized };
    enum class Cap { Keep, Size };
    bool evidence = false;
    Mode mode = Mode::Off;
    Unsized unsized = Unsized::Break;
    Cap cap = Cap::Size;
    bool pairs = true;
    bool plasmidRule = true;
    bool closeGapsPolicy = true;
    std::string admitPath;
    long tangleMax = 5000;
    bool policyWithoutAct = false;   // a policy flag was set outside act mode (ignored)
    // C4 round 3 (appended; read from TESSERACT_OM2_CLONAL=1 + TESSERACT_OM2_CLONAL_CAPWALK=1, never set otherwise): a
    // capped layout gap the cap rule would break is sized at the isolate's isolated shortest walk instead (cap_walk_c1)
    bool capWalkC1 = false;

    bool active() const { return evidence || mode != Mode::Off; }
    bool act() const { return mode == Mode::Act; }
    static SeamConfig fromEnv();
    static const char* modeName(Mode m);
    static const char* unsizedName(Unsized u);
    static const char* capName(Cap c);
};

struct SeamCounters {
    size_t judged = 0, resolverJudged = 0;
    size_t passExact = 0, passWalk = 0, resize = 0, breakDoubleUse = 0, contraPairsOk = 0,
           contraRefused = 0, pairsDeny = 0, silent = 0, unanchored = 0, abstain = 0;
    size_t buttSeen = 0, buttRetired = 0, capSeen = 0, capRetired = 0;
    size_t overlapSeen = 0, overlapRefused = 0, plasmidBlocked = 0;
    size_t pairTestable = 0, pairAbstain = 0, pairConfirm = 0, pairDeny = 0;
    size_t admitScaf = 0, admitGenomeOnly = 0, admitBreak = 0;
    size_t broken = 0, resized = 0;
    size_t closegapsSeen = 0, closegapsDisallowed = 0, closegapsBlocked = 0, closegapsClosed = 0;
    size_t plasmidDepthExcluded = 0;
    double judgeMs = 0, pairBuildMs = 0;
};

struct JudgeRequest {
    Source source = Source::Join;
    const char* stage = "";
    const std::string* left = nullptr;    // oriented so the junction is at its RIGHT end
    const std::string* right = nullptr;   // oriented so the junction is at its LEFT end
    int32_t claimedN = 0;                 // N the release writes; -L for an overlap merge of L
    int32_t panelGap = 0;                 // the model's / track's own gap before any clamp
    Port a, b;                            // pass-local ports, provenance only
    bool plasmidPass = false;
    bool chromosomalA = false, chromosomalB = false;
    bool plasmidicA = false, plasmidicB = false;
    uint32_t panelSupport = 0, panelPairs = 0;
    bool judgeOnly = false;               // recorded, never acted on (resolver N-runs)
    // Round 2 (clonal only, TESSERACT_OM2_CLONAL_VOUCH): the nearest relatives agree on this clonal join (all of them,
    // at the size the isolate's own walk or exact overlap gives) and no read pair names a third record end. The
    // pair test's DENY (too few crossing pairs) is then not a veto where an end is a repeat: pairs whose mate falls
    // in a collapsed repeat are not uniquely placed, so a repeat-bounded seam starves the crossing count (dev E. coli
    // clone GCF051549665v1: 7 of 7 PAIRS_DENY seams at repeat ends were true adjacencies).
    bool clonalVouch = false;
};

struct Decision {
    bool keep = true;
    int32_t writeN = INT32_MIN;           // INT32_MIN: write exactly what the release writes
    uint32_t junction = UINT32_MAX;
    bool vouched = false;                 // a clonal vouch kept a PAIRS_DENY seam at a repeat end
};

// Per-junction audit detail beyond the frozen Junction record.
struct JunctionAudit {
    std::string stage;
    Verdict graphVerdict = Verdict::NotJudged;   // before the pair test
    std::string action;                           // keep | break | resize:N | would_* (audit)
    std::string why;                              // the rule that decided the action
    std::string anchorA, anchorB;
    std::string anchorSeqA, anchorSeqB;   // the anchor 63-mers as they read in the oriented pieces
    int32_t slideA = -1, slideB = -1;     // their distance to the junction
    // Location signature: 63-mers at 0, 250, 1000, 3000, 8000 and 20000 bp from the junction on
    // each side, as they read in the oriented pieces. Repeat-bounded junctions share their
    // flanks and even their graph anchors (a collapsed repeat is copied into several contigs);
    // the farther probes tell them apart.
    std::vector<std::pair<int32_t, std::string>> probesL, probesR;
    uint32_t distinct = 0;
    bool capped = false;
    bool reachable = false;
    uint32_t innerUnique = 0, innerUniquePlaced = 0;
    uint32_t uOther = 0;
    bool pairTestable = false, pairAbstain = true;
    float pairLog10LR = 0, ctrlA = -1, ctrlB = -1;
    int pairCall = 0;
    bool overlapSingle = false;
    uint32_t originA = UINT32_MAX, originB = UINT32_MAX;
    Port passA, passB;
    int32_t panelGap = 0;
    bool judgeOnly = false;
    bool closedByGapfill = false;
    bool gapSeenAtClose = false;
    std::string finalRecord;
    long finalStart = -1, finalLen = -1;
    char finalOrient = '.';
    bool finalFuzzy = false;
};

// The class-level calibration table (C3). Missing file or class -> the v0 admission rule.
struct AdmitTable {
    bool loaded = false;
    double tauS = -1, tauG = -1;
    std::map<std::string, double> upper;
    bool load(const std::string& path, std::string& err);
};

class SeamContext {
public:
    // nullptr unless a TESSERACT_OM2_* flag asks for it. `graph` must outlive the context.
    static std::unique_ptr<SeamContext> create(const UnitigGraph& graph,
                                               const std::vector<std::string>& pieces,
                                               const SequenceStore& reads,
                                               const InsertModel& insert, int threads);
    // For tests: assemble a context from parts.
    SeamContext(const SeamConfig& cfg, std::unique_ptr<EvidenceIndex> ev,
                std::unique_ptr<PairIndex> pairs, std::unique_ptr<PortMap> ports);

    const SeamConfig& config() const { return cfg_; }
    const EvidenceIndex* evidence() const { return ev_.get(); }
    bool gating() const { return cfg_.mode != SeamConfig::Mode::Off && ev_ != nullptr; }
    bool plasmidRuleActive() const { return cfg_.act() && cfg_.plasmidRule && ev_ != nullptr; }
    // Depth says plasmid: >= 1.75 theta (C1d). Only meaningful when plasmidRuleActive().
    bool plasmidByDepth(double cov) const;
    void countPlasmidDepthExcluded() { ++counters.plasmidDepthExcluded; }

    void beginBatch(const std::vector<std::string>& pieces);
    Decision judge(const JudgeRequest& rq);
    // Judges (never acts on) every N-run already in `pieces` (the resolver's scaffold joins).
    void recordResolverGaps(const std::vector<std::string>& pieces);

    // C1e. Fills `allow` with one entry per N-run in closeGaps order and returns true when the
    // mask is to be applied (act mode and the policy on); always maps gaps to junctions.
    bool beforeCloseGaps(const std::vector<std::string>& seqs, std::vector<uint8_t>& allow);
    void afterCloseGaps(const std::vector<std::string>& seqs, size_t blocked);

    bool writeTsv(const std::string& path, const std::vector<std::string>& outSeqs,
                  const std::vector<std::string>& outNames, std::string& err);
    static void printCounters(const SeamContext* ctx, std::FILE* f);

    const Ledger& ledger() const { return ledger_; }
    // Integration (build_om2). C2 attaches its fills, allocation labels and closeGaps decisions to
    // the rows C1 recorded (never adding or removing a row); nothing else writes through this.
    Ledger& ledgerForUpdate() { return ledger_; }
    // Integration (build_om2): the C1 row behind every N-run of `seqs` (one vector per record, runs
    // in order), located exactly as for closeGaps and om2_seams.tsv. Rows that do not locate, or
    // lose a tie, own no run (id = -1).
    std::vector<std::vector<RunOwner>> runOwners(const std::vector<std::string>& seqs) const;
    const std::vector<JunctionAudit>& audit() const { return audit_; }

    // C4 (om2_clonal; used only when TESSERACT_OM2_CLONAL is set):
    //   pairIndex / portMap  read-only access to the pair index and the origin-port map
    //   clonalAdmit  a join the clonal stage proposed and C1 kept: genome view only, never handed to the
    //                gap filler, optionally never filled (noFill), and the clonal rule in om2_seams.tsv
    //   clonalBreak  a panel-only join the nearest relatives contradict: the clonal stage splits the
    //                record at its N-run, and the row reads as broken everywhere afterwards
    const PairIndex* pairIndex() const { return pairs_.get(); }
    const PortMap* portMap() const { return ports_.get(); }
    void clonalAdmit(uint32_t id, bool noFill, const std::string& why, int32_t mergeOverlap = 0);
    void clonalBreak(uint32_t id, const std::string& why);
    //   clonalSize   round 2: the relatives and a walk class agree on this junction's size; n > 0 resizes its N-run
    void clonalSize(uint32_t id, int32_t n);
    //   clonalNoFill round 3 (KC2): a non-positional junction whose walk passes an element stays a labelled gap:
    //                no stage writes bases into it (C2 honours Junction::noFill); n > 0 also resizes its N-run
    void clonalNoFill(uint32_t id, const std::string& why);
    SeamCounters counters;
    double createMs = 0;

private:
    struct Loc { long seq = -1; long run = -1; char orient = '.'; bool fuzzy = false; int score = 0; };
    int probeScore(const std::string& sq, long a, long b, uint32_t id, char orient) const;
    std::string classKey(const Junction& j, const JunctionAudit& x) const;
    std::vector<std::pair<long, long>> nRuns(const std::string& s) const;
    // Where each junction's N-run is in `seqs`: by its graph-unique anchor 63-mers and their
    // slides (both orientations), falling back to the 32 bp flanks. `owner` lists, per N-run in
    // closeGaps order, the junction it belongs to or -1; `runs` receives the N-runs per sequence.
    void locate(const std::vector<std::string>& seqs, std::vector<Loc>& loc, std::vector<long>& owner,
                std::vector<std::vector<std::pair<long, long>>>& runs) const;

    SeamConfig cfg_;
    std::unique_ptr<EvidenceIndex> ev_;
    std::unique_ptr<PairIndex> pairs_;
    std::unique_ptr<PortMap> ports_;
    AdmitTable admit_;
    Ledger ledger_;
    std::vector<JunctionAudit> audit_;
    std::map<std::string, std::vector<std::pair<uint32_t, char>>> flankAll_;   // ... -> every (id, orient)
    std::vector<long> gapJunction_;          // closeGaps order -> junction id, -1 none
};

}  // namespace om2
}  // namespace ts
