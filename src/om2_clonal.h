// Organism Model 2.0, component C4 (CLONAL): clonal-anchored closure. Default OFF.
//
// The layout (organism_layout.cpp) lays the assembly against the ONE panel chromosome it most
// resembles and asks nothing more of the panel. Bacterial populations are clonal: an isolate's
// nearest relatives in a dense panel (Klebsiella has thousands of closed genomes) usually share
// its gene order, its rRNA operon loci and most of its IS insertion sites. This stage asks the
// isolate's k NEAREST relatives, found from its own k-mers, to propose the order of the records
// the earlier stages could not join, and it asks them about every junction the genome view
// asserts -- and it keeps the isolate's own data as the judge:
//
//   1. nearest relatives   Jaccard of the isolate's sampled 31-mer markers (every live unitig of
//                          the final graph) against every layout track of the model, turned into a
//                          Mash-style distance; the k nearest are listed in genome/closure.txt
//   2. clonal layout       every record end is placed on each relative's chromosome by its nearest
//                          assembly-unique markers; a join is proposed only where at least k_min
//                          relatives name the same partner end and the proposal is mutual. Joins
//                          are sized from this isolate's graph walk when a walk class matches the
//                          relatives' gap, else from the relatives when there is no walk at all;
//                          a join the graph or read pairs contradict is never proposed. C1 judges
//                          every proposal (it may resize or break it) and C2 fills it from the
//                          isolate's own graph, with its allocation labels
//   3. junction classes    positional (the relatives' gap recurs at >= 90% of the whole panel:
//                          rRNA operons, tRNA arrays, conserved loci) vs non-positional (the gap
//                          varies across the panel: IS, transposons, prophages, islands). A
//                          non-positional element is written with bases only when k_min
//                          relatives carry it at the same size, the graph has flank->element and
//                          element->flank edges, and no read pair contradicts; otherwise the
//                          junction stays a labelled gap sized by the isolate's graph (noFill)
//   4. contradictions      a panel-only join (tier D/E) that k_min relatives contradict and none
//                          supports is broken (P2: a 5-nearest-donor prior flagged 29-70% of
//                          misjoin loci with no false alarm on correct seams)
//   5. labels              every junction of the genome view gets class, confidence, claim
//                          (confident | supported | provisional | gap), nearest-relative support,
//                          graph edges and contradicting pairs in genome/junctions.tsv (and a
//                          comment line in genome.agp)
//
// Every base written comes from the isolate's graph (C2); the panel contributes order, size and
// labels, never sequence. With TESSERACT_OM2_CLONAL unset nothing here runs: the stage prints its
// counter line with enabled=0 and zeros, and every output file is unchanged.
//
// Flags (all default OFF / release behaviour):
//   TESSERACT_OM2_CLONAL=1               the stage (needs --organism with layout tracks)
//   TESSERACT_OM2_CLONAL_K=<n>           relatives consulted (default 5)
//   TESSERACT_OM2_CLONAL_KMIN=<n>        relatives that must agree (default 3)
//   TESSERACT_OM2_CLONAL_LAYOUT=0|1      propose joins (default 1; 0 = labels only)
//   TESSERACT_OM2_CLONAL_BREAK=0|1       break contradicted panel-only joins (default 1)
//   TESSERACT_OM2_CLONAL_NONPOS=fill|gap non-positional elements with bases when agreed (default fill)
//   TESSERACT_OM2_CLONAL_DMAX=<d>        no clonal joins when the tool's d_near exceeds d (default 1)
//   TESSERACT_OM2_CLONAL_CONFIDENT=<p>   claim `confident` at p >= this (default 0.99)
//   TESSERACT_OM2_CLONAL_CONF_DMAX=<d>   no `confident` claim when d_near exceeds d (default 1)
//   TESSERACT_OM2_CLONAL_CAL=<path>      calibration table (class|tier|nr|group -> p), fitted on dev
//   TESSERACT_OM2_CLONAL_MINMARK=<n>     assembly-unique markers that must agree to place an end (3)
//   TESSERACT_OM2_CLONAL_WEIGHT=distance|count  relatives' votes weighted by distance (default) or counted
//   TESSERACT_OM2_CLONAL_DSCALE=<d>      distance scale of the vote weight exp(-(d - d_nearest)/scale) (0.0005)
//   TESSERACT_OM2_CLONAL_NRP=<path>      nearest-relative plasmid sidecar (<org>.om2nrp, devtools/om2_nrp_build):
//                                        plasmid records are then chained on the relatives' plasmids (2 agreeing,
//                                        depth within 1.6x, and C1's plasmid rule: pairs or an exact walk)
// Round 2 (dev round 1 attribution, om2/iterate/clonal_round1.md):
//   TESSERACT_OM2_CLONAL_OVERLAP=0|1     exact-overlap joins (default 1): the relatives place two record ends
//                                        overlapping and the isolate's two ends share exactly that many bases (>= 25,
//                                        within max(30, 10%) of the relatives' overlap): joined, the shared bases written
//                                        once in the genome view (Junction::mergeOverlap); a graph contradiction is
//                                        overruled only when every placed relative agrees (>= 2 placed)
//   TESSERACT_OM2_CLONAL_VOUCH=0|1       clonal vouch (default 1): a clonal join every placed relative agrees on (>= 2),
//                                        at the walk / exact-overlap size, with no pair to a third end, is not vetoed by
//                                        C1's pair DENY where an end is a repeat (JudgeRequest::clonalVouch)
//   TESSERACT_OM2_CLONAL_RESIZE=0|1      clonal size check of panel-sized joins (default 1): every placed relative
//                                        (>= 2) agrees on one size and a walk class of the isolate has it: the N-run is
//                                        resized to it when it differs by > max(500, 10%), and C2 refuses a fill of
//                                        another length (Junction::clonalSized)
//   TESSERACT_OM2_CLONAL_BRACKET=0|1     size bracket (default 0; diagnostic columns always): one walk class and the
//                                        relatives' gap differ by <= 2 kb: an N-run of their midpoint (never filled,
//                                        never `confident`)
// Round 3 (dev round 2 attribution, om2/iterate/clonal_round2.md); all default OFF:
//   TESSERACT_OM2_CLONAL_WALKSIZE=0|1    data-driven gap sizing from THIS isolate's graph: (a) a clonal link whose
//                                        relatives' gap matches no walk class is sized at the shortest walk class when
//                                        that class is isolated (exhaustive single class, or the next class >= 5 kb
//                                        longer: the longer classes loop through a repeat) and <= 20 kb; (b) a
//                                        panel-sized N-run that is no walk class of the isolate is resized to the walk
//                                        class the relatives' median gap matches, else (>= 1 relative placed, none at the
//                                        asserted size) to the isolated shortest class. Dev Kp round 2: the isolate's
//                                        shortest walk was the true gap where the panel size was wrong (5,567 N at a true
//                                        572; 24 of 26 truth-checked size_ambiguous links)
//   TESSERACT_OM2_CLONAL_KC2=0|1         kill criterion KC2 (IS-PREC-J failed on dev): a junction whose relatives' gap does
//                                        not recur at >= 90% of the panel (non-positional) and whose walk passes repeat
//                                        sequence (an element) is never filled: a labelled gap of its size (noFill), for
//                                        panel-sized N-runs (before C2) and clonal joins alike
//   TESSERACT_OM2_CLONAL_CONFRULE=0|1    `confident` only with >= 2 relatives placed and ALL agreeing by count, both
//                                        N-free sides >= 2 kb, the isolate's walk (if any) at the asserted size, a gap
//                                        <= 10 kb, and -- when the run shows discord (read pairs refusing >= DISCORD
//                                        clonal links or breaking contradicted joins) -- no repeat end; else capped at
//                                        0.985 (`supported`). Dev Kp round 2 counterfactual: 0.995 / 0.996 (LTO / LCO2)
//   TESSERACT_OM2_CLONAL_DISCORD=<n>     the discord threshold of CONFRULE (default 2; 0 = never)
//   TESSERACT_OM2_CLONAL_UNVERIFIED=0|1  a panel-sized gap (tier D/E) that the nearest relatives do not size (fewer than 2
//                                        placed, or not all at the asserted size) and that WALKSIZE did not resize is broken
//                                        when its size is no walk class of this isolate (join or layout gap), or when a
//                                        layout gap has no walk at all: the size is unverifiable, a wrong size a
//                                        misassembly (dev om2r1, 7 species, Rule-Q: join gaps no walk class 273 of 396
//                                        wrong, layout gaps without a walk 192 of 354 wrong)
//   TESSERACT_OM2_CLONAL_INDEL=0|1       the record-collinearity check of end placements tolerates an indel between the
//                                        isolate and a relative (span on the relative within max(10 kb, 5%) + 250 kb of the
//                                        record length, same direction): a record joined across a site where the relative
//                                        carries a 43-122 kb island this isolate lacks keeps its end placements (dev Kp
//                                        GCF001022035v1, CAPWALK on: 2 of 5 relatives lost both ends of a 1.24 Mb record)
//   TESSERACT_OM2_CLONAL_CAPWALK=0|1     C1's cap rule (a capped layout gap whose walks neither agree nor match the panel
//                                        gap is broken: cap_break) sizes the gap at the isolated shortest walk instead
//                                        (cap_walk_c1; read by SeamConfig, needs TESSERACT_OM2_CLONAL=1). Dev Kp round 2:
//                                        13 of 13 truth-checked cap_break rows with such a walk were true at that walk
// Round 3c (dev round 3 attribution, om2/iterate/clonal_round3.md); default OFF:
//   TESSERACT_OM2_CLONAL_STRICT=0|1      emission rule for every model-proposed junction of tier C/D/E (panel-sized
//                                        N-run or clonal link) at its final size n: it is written only when THIS
//                                        isolate's graph sizes it (n is its shortest walk class: the only class of an
//                                        exhaustive enumeration, or isolated by >= 10 kb, or -- when the distance-weighted
//                                        relatives do not contradict it more than they agree -- the shortest of several),
//                                        or when >= 2 relatives are placed, all read n and none contradicts (walk-less, or
//                                        n a walk class); a periodic walk set (>= 3 equally spaced classes: a tandem array
//                                        whose copy number the walks cannot fix) also needs weighted agreement >= 1.5 and
//                                        >= 3x the weighted contradiction. Otherwise the run is broken (labelled gap between
//                                        records) and the link refused. Overlaps and 1-N butts always pass. Dev round-3b
//                                        counterfactual (Rule-Q, genome view): FALSE junctions Kp 37 -> 26 (LTO) / 29
//                                        (LCO2), Ab 171 -> 102, no true join of a chr98 chromosome record lost
//   TESSERACT_OM2_CLONAL_STRICT_NR0=0|1  with STRICT: a junction no relative places (placed = 0) is broken too, even at
//                                        the isolate's shortest walk (dev counterfactual Kp LTO 26 -> 18 FALSE, one chr98
//                                        join lost; A/B only)
//
// Counter line, printed on every run (zeros when off): [om2-clonal].
#pragma once

#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <functional>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "om2_types.h"
#include "organism.h"

namespace ts {
class UnitigGraph;

namespace om2 {
class SeamContext;
class EvidenceIndex;
class PairIndex;

struct ClonalOptions {
    bool enabled = false;
    int k = 5;
    int kmin = 3;
    int report = 10;             // relatives listed in closure.txt (>= k)
    bool layout = true;
    bool breakContra = true;
    bool nonposFill = true;
    double dmax = 1.0;
    double confident = 0.99;
    double confDmax = 1.0;
    std::string calPath;
    int minMarkers = 3;
    // Relatives vote with weight exp(-(d - d_nearest) / dscale): the closest relatives propose (CLONAL CLOSURE
    // PRINCIPLE); a relative 10x farther than a clonemate barely counts, equally close relatives vote as a majority.
    // `count` = the unweighted rule (k_min of the k relatives).
    bool weighted = true;        // TESSERACT_OM2_CLONAL_WEIGHT=distance|count
    double dscale = 0.0005;      // TESSERACT_OM2_CLONAL_DSCALE
    double wmin = 0.9;           // weighted support a proposal needs (the nearest relative alone qualifies)
    std::string nrpPath;         // nearest-relative plasmid sidecar (<org>.om2nrp), md5-pinned to the model
    int kminPlasmid = 2;         // relatives whose plasmid must name the same partner end
    double plasmidDepthRatio = 1.6;   // plasmid records joined only within this depth ratio
    bool anyFlagSet = false;     // some TESSERACT_OM2_CLONAL_* is set
    bool overlapMerge = true;    // TESSERACT_OM2_CLONAL_OVERLAP (round 2)
    bool vouch = true;           // TESSERACT_OM2_CLONAL_VOUCH (round 2)
    bool bracket = false;        // TESSERACT_OM2_CLONAL_BRACKET (round 2)
    bool resize = true;          // TESSERACT_OM2_CLONAL_RESIZE (round 2)
    int32_t minOverlap = 25;     // exact overlap accepted from this length (fixed)
    int64_t bracketMax = 2000;   // walk vs relatives difference a bracket may span (fixed)
    bool unverified = false;     // TESSERACT_OM2_CLONAL_UNVERIFIED (round 3b)
    bool indel = false;          // TESSERACT_OM2_CLONAL_INDEL (round 3b)
    int64_t indelMax = 250000;   // INDEL: largest isolate / relative length difference tolerated (fixed)
    bool walksize = false;       // TESSERACT_OM2_CLONAL_WALKSIZE (round 3)
    bool kc2 = false;            // TESSERACT_OM2_CLONAL_KC2 (round 3)
    bool confRule = false;       // TESSERACT_OM2_CLONAL_CONFRULE (round 3)
    int discord = 2;             // TESSERACT_OM2_CLONAL_DISCORD (round 3)
    int64_t c1Gap = 5000;        // the next walk class at least this much longer: the shortest class is isolated (fixed)
    int64_t c1GapLayout = 10000; // round 3b: the same for a panel-sized N-run (the relatives' order only) (fixed)
    int64_t c1Max = 20000;       // an isolated shortest class sizes a gap only up to this length (fixed)
    int64_t confSide = 2000;     // CONFRULE: N-free bases on each side (fixed)
    int64_t confMaxGap = 10000;  // CONFRULE: largest asserted gap (fixed)
    // 1.5 (KC1, EVAL_PLAN_CLONAL s9: triggered under the binding scorer at round 3): a junction that reaches the
    // `confident` threshold is written `supported` unless TESSERACT_OM2_CLONAL_ALLOW_CONFIDENT=1 (dev A/B and the E0g
    // identity check against the round-3c candidate only). The numeric confidence column is unchanged: it is the
    // hand-set prior, never fitted (no TESSERACT_OM2_CLONAL_CAL table exists).
    bool allowConfident = false; // TESSERACT_OM2_CLONAL_ALLOW_CONFIDENT
    bool strict = false;         // TESSERACT_OM2_CLONAL_STRICT (round 3c)
    bool strictNr0 = false;      // TESSERACT_OM2_CLONAL_STRICT_NR0 (round 3c, A/B)
    double strictPeriodicWa = 1.5;     // STRICT: weighted agreement a periodic walk set needs (fixed)
    double strictPeriodicRatio = 3.0;  // STRICT: ... and at least this multiple of the weighted contradiction (fixed)

    // Fixed by this design (not flags).
    int window = 30000;          // context searched for markers at an end
    int useMarkers = 6;          // nearest present markers used to place one side
    int64_t maxAdjacent = 100000;// a relative's gap beyond this is not an adjacency
    int64_t maxOverlap = 1000;   // records may overlap this much on a relative
    double positional = 0.90;    // panel share at the relatives' gap: positional locus
    uint32_t panelMin = 10;      // tracks with both sides placed before the panel share is used
    uint32_t minContra = 3;      // read pairs naming another record end: contradiction
    double agreeShare = 0.66;    // share of placed relatives naming the modal partner

    static ClonalOptions fromEnv();
};

struct Relative {
    uint32_t track = 0;
    std::string name;            // panel accession (the track name)
    uint32_t shared = 0, trackMarkers = 0;
    double jaccard = 0, containment = 0, dist = 1;
};

struct NearestSet {
    bool ok = false;
    uint32_t isoMarkers = 0, isoChrMarkers = 0, tracks = 0;
    std::vector<Relative> ranked;   // nearest first, up to options.report
    double dNear() const { return ranked.empty() ? 1.0 : ranked.front().dist; }
    // The analysis group of EVAL_PLAN_CLONAL s3 by the tool's OWN estimate: C (<= 0.0005),
    // L (<= 0.002), F (beyond), '-' when no relative was found.
    char group() const;
    double ms = 0;
};

// One side of a junction placed on a panel chromosome: `x` is the first base beyond the side in
// the junction's direction (exit side: the base after the left context; entry side: the first
// base of the right context), `dir` +1 when the junction's direction runs along the chromosome.
struct SidePlace {
    bool ok = false;
    int64_t x = 0;
    int8_t dir = 0;
    uint16_t n = 0;              // markers agreeing on the placement
    int32_t nearest = 0;         // distance of the nearest agreeing marker from the junction
};

enum class TrackCall : uint8_t { Absent, Agree, GapDiff, Contra };

struct MarkerHit {
    uint32_t id = 0;
    int32_t pos = 0;             // k-mer start in the context
    uint8_t o = 0;               // 0 when the forward k-mer of the context is canonical
    int32_t dist = 0;            // bases between the k-mer and the junction
};

// A junction between two oriented contexts: the junction is at the END of `left` and at the START
// of `right`. `asserted` = bases asserted between them (N-run, fill length, or -overlap).
struct JunctionQuery {
    std::vector<MarkerHit> left, right;   // assembly-unique marker hits, nearest the junction first
    int64_t leftLen = 0;
    int64_t asserted = 0;
    bool sized = true;
};

struct TrackSupport {
    uint32_t tracks = 0, placed = 0, agree = 0, gapDiff = 0, contra = 0;
    double wPlaced = 0, wAgree = 0, wGapDiff = 0, wContra = 0;   // relatives only: distance-weighted
    std::vector<int64_t> gaps;     // implied gap of every track whose two sides agree in direction
    std::vector<int64_t> agreeGaps;
    std::vector<uint32_t> agreeRank;   // for relatives: rank (0 = nearest) of every agreeing one
};

// The clonal classes of EVAL_PLAN_CLONAL s6.1 / s6.7.
const char* clonalClassName(int c);
enum ClonalClass { kPositional = 0, kNonPositional, kLayout, kPlasmid, kWrap, kGraph, kUnknownClass };

// What the stage and the writer know about one junction (junctions.tsv columns).
struct ClonalAnnot {
    bool set = false;
    int cls = kUnknownClass;
    uint32_t nrK = 0, nrPlaced = 0, nrAgree = 0, nrGapDiff = 0, nrContra = 0;
    double wPlaced = -1, wAgree = 0, wContra = 0;   // distance-weighted (wPlaced < 0: not computed)
    int64_t nrGap = INT64_MIN;       // median gap of the agreeing relatives (else of the placed ones)
    double dNear = 1.0;
    char group = '-';
    uint32_t panelPlaced = 0, panelAgree = 0;
    int edgeA = -1, edgeB = -1, walk = -1, walkMatch = -1;   // -1 unknown, 0 no, 1 yes
    int element = -1;                // the shortest walk passes repeat sequence (an element) or not
    uint32_t pairContra = 0;
    bool sized = true;
    double confidence = -1;
    std::string claim = ".";
    std::string action = ".";        // join | break_contra | keep | refuse:<why> | .
    std::string basis = ".";         // graph+nr | graph | nr | pairs | ledger
    std::string nrBin() const;       // none | one | all | maj | min | contra
    double panelFrac() const { return panelPlaced ? static_cast<double>(panelAgree) / panelPlaced : -1; }
};

// The evidence behind one proposed clonal join and the rule that decides it (pure; unit-tested).
struct LinkEvidence {
    uint32_t k = 0, placed = 0, agree = 0;
    double wAgree = 0, wContra = 0;  // distance-weighted votes for the partner / for another partner
    int64_t nrGap = 0, nrSpread = 0;
    uint32_t panelPlaced = 0, panelAgree = 0;
    bool anchoredA = false, anchoredB = false, edgeA = false, edgeB = false;
    bool reachable = false, exhaustive = false, contraGraph = false;
    bool repeatA = false, repeatB = false, element = false;
    std::vector<int32_t> walkLens;   // walk-length classes (bases between the two ends)
    std::vector<uint16_t> walkCounts; // round 3: walks per class (parallel to walkLens; empty = unknown)
    uint32_t pairContra = 0;
    bool plasmid = false;
    double wPlaced = 0;              // round 2: distance-weighted relatives placing the end
    int32_t overlapExact = -1;       // round 2: exact suffix/prefix overlap of the two ends near the relatives' overlap
};
struct LinkDecision {
    bool join = false;
    int32_t n = 0;
    bool noFill = false;
    bool sized = true;
    int cls = kUnknownClass;
    bool walkMatch = false;
    std::string basis, why;
    bool overlap = false;            // round 2: joined on an exact overlap (n = -overlap)
    bool bracket = false;            // round 2: sized at the midpoint of walk and relatives
    int32_t bracketN = INT32_MIN;    // the bracket size, computed whenever it applies (diagnostic)
    bool c1 = false;                 // round 3: sized at the isolated shortest walk class (WALKSIZE)
};
// Round 3: the isolated shortest walk class of `lens` (counts parallel, may be empty), or INT32_MIN (unit-tested).
int32_t isolatedShortestWalk(const std::vector<int32_t>& lens, const std::vector<uint16_t>& counts, bool exhaustive,
                             int64_t gap, int64_t maxLen);
LinkDecision decideLink(const LinkEvidence& e, const ClonalOptions& o);
// Round 3c (STRICT): the evidence of one model-proposed junction at its final size, and the emission rule (pure;
// unit-tested). strictVerdict returns nullptr to write the junction, else the reason it is broken / refused.
struct StrictInput {
    int64_t n = 0;                    // the size to be written (bases between the two ends; <= 0 an overlap)
    bool overlap = false;             // an exact-overlap join
    std::vector<int32_t> walkLens;    // this isolate's walk classes between the two ends
    std::vector<uint16_t> walkCounts; // walks per class (parallel; may be empty)
    bool exhaustive = false;          // the walk enumeration was exhaustive
    uint32_t placed = 0, agree = 0, contra = 0;   // nearest relatives placing both ends / reading n / contradicting
    double wAgree = 0, wContra = 0;               // the same, distance-weighted
};
bool strictPeriodic(const std::vector<int32_t>& sortedLens);
const char* strictVerdict(const StrictInput& s, const ClonalOptions& o);
// Round 2: exact suffix/prefix overlap of L and R near `want` bases (-1: none); unit-tested.
int32_t exactEndOverlap(const std::string& L, const std::string& R, int64_t want, int32_t minOv, int64_t maxOv);

// Confidence model (PROPOSED defaults; a table fitted on dev replaces them per key).
struct ClonalCalibration {
    bool loaded = false;
    std::map<std::string, std::pair<double, uint32_t>> table;   // key -> (p, n)
    bool load(const std::string& path, std::string& err);
    // key: <class>|<tier>|<nr bin>|<group>
    static std::string key(const ClonalAnnot& a, char tier);
    double confidence(const ClonalAnnot& a, char tier, bool hasLedger) const;
};

class ClonalEngine {
public:
    ClonalEngine(const OrganismModel& model, const ClonalOptions& opt);

    const ClonalOptions& options() const { return opt_; }
    const NearestSet& nearest() const { return near_; }
    const OrganismModel& model() const { return model_; }

    // Nearest relatives from every live unitig of the graph (or, for tests, from sequences).
    void selectNearest(const UnitigGraph& g);
    void selectNearestFromSeqs(const std::vector<std::string>& seqs);

    // The assembly whose markers count as unique (a marker seen twice in it names a repeat).
    void setAssembly(const std::vector<std::string>& seqs);

    // Marker hits of a context: `exitSide` = the junction is at the END of s (hits from its last
    // `window` N-free bases), else at its START.
    std::vector<MarkerHit> hits(const std::string& s, bool exitSide) const;
    JunctionQuery query(const std::string& left, const std::string& right, int64_t asserted, bool sized) const;

    // One side on one relative (rank < k) / on any track given its marker index.
    SidePlace placeOnRelative(const std::vector<MarkerHit>& h, bool exitSide, int64_t ctxLen, size_t rank,
                              size_t need) const;
    static TrackCall call(const SidePlace& a, const SidePlace& b, int64_t trackLen, int64_t asserted, bool sized,
                          int64_t maxAdjacent, int64_t maxOverlap, int64_t& gap);

    // Support of each query among the k relatives, and among every track of the panel.
    TrackSupport relativeSupport(const JunctionQuery& q) const;
    void panelSupport(const std::vector<JunctionQuery>& qs, std::vector<TrackSupport>& out) const;

    size_t relatives() const { return nr_.size(); }
    // Vote weight of the relative of rank r: exp(-(d_r - d_0) / dscale), or 1 under the count rule.
    double weight(size_t rank) const;

    // Nearest-relative plasmids (after selectNearest): the plasmid records of the k nearest relatives, from the
    // sidecar. Refused unless its #tsm_md5 is `modelMd5`.
    bool loadPlasmids(const std::string& path, const std::string& modelMd5, std::string& err);
    size_t plasmidTracks() const { return pl_.size(); }
    size_t plasmidRank(size_t pi) const { return pl_[pi].rank; }
    int64_t plasmidLength(size_t pi) const { return pl_[pi].length; }
    const std::string& plasmidName(size_t pi) const { return pl_[pi].t.name; }
    SidePlace placeOnPlasmid(const std::vector<MarkerHit>& h, bool exitSide, size_t pi, size_t need) const;
    // What the stage decided per ledger row (junction id -> annotation); the writer reads it.
    std::map<uint32_t, ClonalAnnot> rows;
    const std::string& relativeName(size_t rank) const;
    int64_t relativeLength(size_t rank) const;

private:
    struct NrTrack {
        uint32_t track = 0;
        std::unordered_map<uint32_t, uint32_t> at;   // marker id -> index on the track
        int64_t length = 0;
    };
    SidePlace placeWith(const std::vector<MarkerHit>& h, bool exitSide, int64_t ctxLen, const LayoutTrack& t,
                        const std::function<bool(uint32_t, uint32_t&)>& find, size_t need) const;
    void buildRelatives();

    const OrganismModel& model_;
    ClonalOptions opt_;
    NearestSet near_;
    std::vector<NrTrack> nr_;
    struct PlTrack {
        size_t rank = 0;             // the relative carrying it
        LayoutTrack t;               // name = <genome>:<record>
        std::unordered_map<uint32_t, uint32_t> at;
        int64_t length = 0;          // exact record length (from the sidecar)
    };
    std::vector<PlTrack> pl_;
    std::unordered_map<uint32_t, uint8_t> asmCount_;   // marker id -> occurrences in the assembly (<= 2)
    bool asmSet_ = false;
};

int64_t trackLength(const LayoutTrack& t, uint64_t denom);

struct ClonalStats {
    bool enabled = false, ran = false;
    std::string skip = "-";          // why the join stage did not run
    int k = 0;
    double dNear = -1;
    char group = '-';
    size_t records = 0, candidates = 0, endsPlaced = 0, proposals = 0, mutual = 0;
    size_t refusedPairs = 0, refusedGraph = 0, refusedSize = 0, refusedC1 = 0, cycleCut = 0;
    size_t joins = 0, joinsPositional = 0, joinsNonPositional = 0, joinsLayout = 0, noFill = 0;
    size_t sizedGraph = 0, sizedNr = 0;
    size_t contraChecked = 0, contraBroken = 0, contraKeptProven = 0;
    size_t chains = 0, recordsOut = 0;
    bool nrpLoaded = false;
    size_t plasmidTracks = 0, plasmidCandidates = 0, plasmidMutual = 0, plasmidJoins = 0, refusedDepth = 0;
    // filled in by the writer
    size_t annotated = 0, confident = 0, supported = 0, provisional = 0, gapClaims = 0;
    double ms = 0;
    // round 2
    size_t overlapJoins = 0, overlapC1Conflict = 0, vouched = 0, bracketWould = 0, bracketJoins = 0, merged = 0;
    size_t sizeChecked = 0, sizeAgree = 0, resized = 0;
    // round 3
    size_t walksizeLinks = 0, walksizeChecked = 0, walksizeResized = 0, kc2Checked = 0, kc2NoFill = 0, kc2Links = 0;
    size_t confCapped = 0, discordEvents = 0, capWalkC1 = 0, unverifiedBroken = 0;
    bool discord = false;
    // round 3c
    bool strict = false, strictNr0 = false;
    size_t strictChecked = 0, strictBroken = 0, strictRefused = 0, strictPeriodic = 0, strictNoRelative = 0;
    // 1.5: `confident` claims written as `supported` (KC1). Not part of the [om2-clonal] counter line (whose format
    // the round-3c identity check E0g compares); printed on the [om2-layout] line and in report.json when non-zero.
    size_t confidentSuppressed = 0;
    bool kc1 = false;            // suppression active (allowConfident off) and the writer annotated
};

struct ClonalInputs {
    const OrganismModel& model;
    const UnitigGraph& graph;
    SeamContext* seam = nullptr;
    bool verbose = false;
    std::string outDir;          // "" = no side files; else <outDir>/om2_clonal/{nearest,ends,links}.tsv
};

// Stage [4b2/7], after the layout and before C2. Rewrites `seqs` (and `covs`, `layoutMembers`)
// only when it joins or breaks something; `changed` says so (gfaPaths no longer hold then).
ClonalStats runClonalStage(const ClonalInputs& in, ClonalEngine& eng, std::vector<std::string>& seqs,
                           std::vector<double>& covs, std::vector<char>& layoutMembers, bool& changed);

// The writer's view (om2_output.cpp): annotations for junctions of the finished genome records.
class ClonalSurface {
public:
    ClonalSurface(ClonalEngine& eng, const SeamContext* seam, ClonalStats& stats);
    struct Req {
        std::string left, right;     // N-free contexts: junction at the END of left, START of right
        int64_t asserted = 0;
        bool sized = true;
        const Junction* j = nullptr; // ledger row, or null (unrecorded seam)
        bool wrap = false;
        char replicon = 'u';         // c | p | u
    };
    void prepare(const std::vector<std::string>& records);   // assembly-unique markers of the view
    std::vector<ClonalAnnot> annotate(const std::vector<Req>& reqs);
    static std::string header();                              // junctions.tsv extra columns
    static std::string columns(const ClonalAnnot& a);
    std::string closureHeader() const;                        // closure.txt lines
    std::string agpComment(const std::string& id, const ClonalAnnot& a) const;
    std::string jsonBlock(int indent) const;
    const ClonalStats& stats() const { return stats_; }
    void noteMerge() { ++stats_.merged; }                     // the writer merged an exact-overlap join
    // The stage's own annotation of a ledger row (a clonal join, a refused or broken one), or null.
    const ClonalAnnot* stageRow(uint32_t id) const {
        auto it = eng_.rows.find(id);
        return it == eng_.rows.end() ? nullptr : &it->second;
    }

private:
    ClonalEngine& eng_;
    const SeamContext* seam_;
    ClonalStats& stats_;
    ClonalCalibration cal_;
};

std::string formatClonalCounters(const ClonalStats& s);
void logClonalCounters(const ClonalStats& s, std::FILE* f = stderr);

}  // namespace om2
}  // namespace ts
