// Organism Model 2.0 -- shared junction types (DESIGN.md §5.1, frozen interface).
//
// One Junction is recorded wherever the assembler writes an N-run or merges two pieces on
// the model's word: the resolver's scaffold joins, the adjacency join (organism_join.cpp)
// and the layout (organism_layout.cpp). C1 (SEAM) judges each one against this isolate's
// own final-k graph and read pairs; C2 (CLOSE) and C3 (SURFACE) read the same records.
//
// Nothing in this header does anything by itself. Every om2 code path is reached only when
// a TESSERACT_OM2_* flag is set, so with the flags unset the release output is unchanged.
//
// Interface change log: om2/design/INTERFACE_CHANGES.md (C1 appended Source::LayoutButt1,
// the layout's 1-N butt, which the frozen list had folded into Source::Layout).
#pragma once

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace ts {
namespace om2 {

enum class Source : uint8_t {
    Resolver,            // PairedResolver renderChain N-run
    ResolverUnknown100,  // resolver N-run of exactly 100 (unknown size)
    Join,                // organism_join: N-run of the panel gap
    JoinButt1,           // organism_join: panel overlap the sequence denied -> 1 N
    JoinOverlap,         // organism_join: exact overlap merge
    Layout,              // organism_layout: N-run of the track gap (<= 2000)
    LayoutCap2000,       // organism_layout: track gap > 2000 written as 2000 N
    LayoutOverlap,       // organism_layout: exact overlap merge
    Wrap,                // circular wrap junction (C2)
    LayoutButt1,         // organism_layout: overlap the sequence denied -> 1 N (appended by C1)
    Clonal               // om2_clonal: a join the isolate's nearest panel relatives propose (appended by C4)
};

enum class Verdict : uint8_t {
    NotJudged,
    PassExact,       // exactly one walk, |gmin - N| <= slack
    PassWalk,        // a walk exists, gmin <= N + slack
    Resize,          // every walk > N + slack; the shortest uses no placed single-copy node
    BreakDoubleUse,  // every walk > N + slack; the shortest needs a placed single-copy node
    ContraPairsOk,   // no walk, a flank's first single-copy continuation is elsewhere; pairs agree
    ContraRefused,   // ... and pairs do not confirm
    PairsDeny,       // spannable seam, LR >= 100 for an insertion
    Silent,          // no walk and no single-copy continuation on either side
    Unanchored,      // a flank has no graph-unique anchor
    Abstain          // tangled graph, or index not built
};

enum class Tier : uint8_t { A, B, C, D, E };
enum class Admit : uint8_t { Scaffold, GenomeOnly, Break };
enum class EndClass : uint8_t { Unknown, Unique, Repeat, DeadEnd };
enum class Basis : uint8_t {
    Observed, GraphWalk, PairPhased, ThreadPhased, IsolateRepeatExact, PriorAllocated,
    Multiplicity, Consensus
};

struct Port {
    uint32_t piece = UINT32_MAX;
    bool tail = false;
};

struct WalkClass {
    int32_t len = 0;   // gap length (bases between the two piece ends) of the class
    uint16_t n = 0;    // walks in the class
};

struct FillSpan {
    uint32_t off = 0, len = 0;
    Basis basis = Basis::Observed;
};

struct Junction {
    uint32_t id = 0;
    Source source = Source::Join;
    Port a, b;
    int32_t claimedN = 0, writtenN = 0;   // writtenN = -1 when broken
    Verdict verdict = Verdict::NotJudged;
    Tier tier = Tier::E;
    Admit admit = Admit::Scaffold;
    std::string cls;                      // class key for calibration
    EndClass endA = EndClass::Unknown, endB = EndClass::Unknown;
    float copyA = 0, copyB = 0;
    int32_t gmin = 0;
    std::vector<WalkClass> walks;
    bool exhaustive = false, hairpin = false;
    uint16_t uA = 0, uB = 0;
    float pairLambda = 0;
    uint32_t pairK = 0, pairContra = 0;
    uint32_t panelSupport = 0, panelGenomes = 0;
    float tieMargin = 0;
    float pMisjoin = -1;                  // calibrated class upper bound; -1 if none
    bool allowCloseGaps = true;
    bool contigFilled = false;
    std::string fillSeq;
    std::vector<FillSpan> fillSpans;
    std::string flankL32, flankR32;       // carried across closeGaps / polish
    // C4 (om2_clonal, appended): no stage may write bases into this junction. Set only by the clonal
    // stage, on a non-positional junction its nearest relatives do not agree on (a labelled gap).
    bool noFill = false;
    // C4 round 2 (appended): a clonal join whose two records overlap EXACTLY by this many bases (the nearest relatives
    // agree on the overlap): the genome view writes the shared bases once (the second component starts after them);
    // scaffolds.fasta keeps both records whole. 0 = no merge.
    int32_t mergeOverlap = 0;
    // C4 round 2 (appended): every placed nearest relative (>= 2) and a walk class of this isolate agree on this
    // junction's size (the clonal stage resized the N-run to it when it differed): C2 then refuses a fill of another
    // length (as at clonal joins), because another length is another route.
    bool clonalSized = false;
};

struct Ledger {
    std::vector<Junction> j;
    int64_t rotateOffset[64];
    Ledger() { for (int64_t& x : rotateOffset) x = -1; }
};

// ---- integration additions (build_om2; logged in om2/design/INTERFACE_CHANGES.md) --------------

// The ledger row behind one N-run of a set of records: C1 locates its junctions in any record set
// (graph-unique anchors, far probes, then the exact 32-bp flanks), and hands that location to C2
// (stage [4c/7], so C2 attaches its fills to C1's rows instead of recording the gap twice) and to
// C3 (the finished records, so repeat-bounded junctions with identical flanks are still told apart).
// `orient` is '+' when the record reads the junction as recorded (flankL32 left of the run), '-'
// when it reads it reverse-complemented; id = -1 when no row owns the run.
struct RunOwner {
    int64_t start = 0, end = 0;   // [start, end) of the N-run in its record
    int64_t id = -1;              // index into Ledger::j
    char orient = '.';
};

// Inter-copy variation at one site of a multi-copy repeat (C2 -> C3; proposed by C3 in
// om2/c3/README.txt, moved here from om2_output.h by the integration). One row per (site, locus):
// the allele written at that locus and its basis label, plus every allele the isolate carries at
// the site with its depth-estimated carrier copies. Alleles and fillOffset are in the orientation
// of the junction's ledger row (its fillSeq). A row with junction = -1 describes a site the isolate
// carries but did not place at any locus.
struct RepeatVariant {
    std::string family;                                   // rrn, IS, repeat, ...
    std::string cluster;                                  // repeat cluster / bubble id
    std::string site;                                     // site id inside the repeat unit
    int64_t junction = -1;                                // locus: the junction whose fill holds the site
    uint32_t fillOffset = 0;                              // 0-based offset of the site in that fill
    uint32_t siteLen = 1;                                 // bases the site spans in the fill
    std::string placed;                                   // allele written at the locus ("" if unplaced)
    Basis basis = Basis::Consensus;                       // label of the placed allele
    std::vector<std::pair<std::string, float>> alleles;   // allele -> estimated carrier copies
    std::string phasing = "none";                         // pair | thread | prior | none
    std::string loci = ".";                               // candidate loci of the minority alleles
};

const char* sourceName(Source s);
const char* verdictName(Verdict v);
const char* tierName(Tier t);
const char* admitName(Admit a);
const char* endClassName(EndClass e);

}  // namespace om2
}  // namespace ts
