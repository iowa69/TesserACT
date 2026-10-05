// Organism Model 2.0, C1a -- the evidence index.
//
// What this isolate's own data say about a junction the model proposes, asked of the FINAL-k
// unitig graph (the graph written to assembly_graph.gfa) and of the read pairs:
//
//   anchors      each piece end is placed in the graph by its first graph-unique 63-mer,
//                sliding inward from the end in 25 bp steps (up to 20 kb)
//   multiplicity copy number per unitig against the length-weighted median depth of its
//                component (nodes >= 2k); a unitig >= 1 kb at <= 1.6 theta is hard single-copy
//   walks        the shortest graph walk between two anchors (bounded Dijkstra, 60 kb) and the
//                set of walks up to 30 kb (<= 32 walks, 200k steps), grouped into +-50 bp
//                length classes, with exhaustive and hairpin flags
//   first-unique the first single-copy unitigs reachable from an end through repeat nodes only
//   pairs        crossing fragments at a seam against the crossing rate 1 kb inside each flank
//                (the locus_prior pilot rule, seam_pairs2.py), Poisson LR, abstain below
//                lambda = 10
//
// The graph-side rules mirror om2/diverge/graph_first/tools/gate_eval.py (the oracle), except
// that theta is per connected component (a multi-copy plasmid component is not a repeat).
// Everything here is read-only with respect to the assembly.
#pragma once

#include <cstdint>
#include <functional>
#include <memory>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

#include "graph.h"
#include "om2_types.h"
#include "resolve.h"   // InsertModel
#include "seqio.h"

namespace ts {
namespace om2 {

struct EvidenceOptions {
    // multiplicity
    int scMinLen = 1000;          // hard single-copy: length >= this ...
    double scMaxTheta = 1.6;      // ... and depth <= this x theta
    bool perComponentTheta = true;
    long compMinEligible = 20000; // a component needs this much >= 2k sequence for its own theta
    // anchors
    int anchorLen = 63;
    int anchorSlide = 20000;
    int anchorStep = 25;
    // walks
    int gminCap = 60000;          // Dijkstra bound (gate_eval BOUND)
    int walkFillMax = 30000;      // enumeration bound (census MAXFILL)
    int maxWalks = 32;            // census MAXPATHS
    long walkSteps = 200000;      // census MAXSTEPS
    int repeatTraverseMax = 12000;// length classes: walks through at most this much repeat
    int uBound = 30000;           // first-unique expansion bound (UBOUND)
    long uSteps = 200000;         // (USTEPS)
    int classTol = 50;            // walk-length class width
    // tangle abstention
    long tangleMax = 5000;        // live nodes
    double tangleScFrac = 0.80;   // share of graph bases in hard single-copy nodes
    // pairs (seam_pairs2.py)
    int pairMaxN = 100;           // pair-testable seams: N <= this
    int pairCtrl = 1000;          // control points this far inside each flank
    int pairAnch = 30;            // a crossing fragment reaches this far past the seam
    double pairMaxFragSd = 4.0;   // MAXF = mean + this x sd
    double pairLambdaMin = 10.0;  // abstain below
    double pairLog10LR = 2.0;     // LR >= 100
};

struct ONode {
    uint32_t id = UINT32_MAX;
    uint8_t o = 0;   // 0 = forward, 1 = reverse complement, in walking direction
};
inline bool operator==(ONode a, ONode b) { return a.id == b.id && a.o == b.o; }
inline bool operator!=(ONode a, ONode b) { return !(a == b); }
inline ONode flip(ONode n) { return ONode{n.id, static_cast<uint8_t>(n.o ^ 1)}; }

// A piece end placed in the graph. For an EXIT anchor (junction at the right end of the
// oriented piece) `off` is where the anchor 63-mer ENDS in the oriented node; for an ENTRY
// anchor (junction at the left end) it is where the anchor STARTS. `slide` is the number of
// piece bases between the anchor and the piece end.
struct EndAnchor {
    bool ok = false;
    const char* why = "none";   // edge | hit_N | slide_absent | slide_multi | ok
    ONode node;
    int32_t off = 0;
    int32_t slide = 0;
};

struct Multiplicity {
    float copy = 0;          // depth / theta of its component
    uint16_t lo = 0, hi = 0; // copy-number interval (depth only)
    bool hardSingle = false;
    bool uncertain = false;  // depth sits between two integers
    uint32_t component = 0;
};

struct WalkSet {
    bool reachable = false;        // a walk exists within gminCap
    int32_t fillMin = -1;          // shortest fill between the anchors
    int32_t gmin = 0;              // ... minus both slides: bases between the piece ends
    std::vector<WalkClass> classes;// enumerated walks by gap-length class (sorted by length)
    uint32_t walks = 0;            // enumerated walks
    uint32_t distinct = 0;         // distinct fill sequences among them
    bool capped = false;           // enumeration stopped at maxWalks or walkSteps
    bool repeatPruned = false;     // class enumeration cut a walk at repeatTraverseMax
    bool classesCapped = false;    // class enumeration stopped at maxWalks / walkSteps
    bool exhaustive = false;       // !capped (census semantics)
    bool hairpin = false;          // some walk uses a unitig in both orientations
    // the shortest walk (skip_path_check.spath), inner nodes only
    std::vector<ONode> inner;
    uint32_t innerUnique = 0, innerUniquePlaced = 0;
    uint64_t innerUniqueBp = 0, innerRepeatBp = 0;
    float innerMaxCopy = 0;
};

struct FirstUnique {
    std::vector<uint32_t> ids;   // sorted
    bool capped = false;
};

class EvidenceIndex {
public:
    // Builds the index over `g`, which must outlive it. Never modifies the graph.
    static std::unique_ptr<EvidenceIndex> build(const UnitigGraph& g, const EvidenceOptions& opt);

    bool abstained() const { return abstained_; }
    const char* abstainReason() const { return abstainWhy_; }
    const EvidenceOptions& options() const { return opt_; }
    int k() const { return k_; }
    double theta() const { return theta_; }
    const UnitigGraph& graph() const { return *g_; }

    EndAnchor anchorExit(const std::string& oriented, int32_t minSlide = 0) const;
    EndAnchor anchorEntry(const std::string& oriented, int32_t minSlide = 0) const;
    const Multiplicity& mult(uint32_t unitig) const { return mult_[unitig]; }
    WalkSet walks(const EndAnchor& exitA, const EndAnchor& entryB) const;
    FirstUnique firstUnique(ONode start) const;

    // Marks the hard single-copy unitigs whose middle 63-mer occurs in `pieces`.
    void setPlaced(const std::vector<std::string>& pieces);
    bool placed(uint32_t unitig) const { return unitig < placed_.size() && placed_[unitig]; }

    // The unitig holding the graph-unique 63-mer at seq[at, at+63), or UINT32_MAX.
    uint32_t locateUnique(const std::string& seq, size_t at) const;

    // Counters for the [om2-evidence] line.
    size_t liveNodes = 0, components = 0, singleCopy = 0, repeatNodes = 0, uncertainNodes = 0;
    double scFraction = 0;
    double buildMs = 0;
    mutable size_t anchorsOk = 0, anchorsFailed = 0;

private:
    struct Ent { uint64_t fp; uint32_t node; uint32_t posStrand; };
    struct Hit { ONode node; int32_t off; };
    // 0 = absent, 1 = unique (hit filled), 2 = multi
    int locate(const std::string& s, size_t from, Hit& hit) const;
    int64_t dmin(ONode sn, int32_t e, ONode tn, int32_t f) const;
    bool spath(ONode sn, int32_t e, ONode tn, int32_t f, std::vector<ONode>& path) const;
    struct Enumeration {
        uint32_t walks = 0;
        bool capped = false, pruned = false;
        std::set<std::string> fills;
        std::vector<std::vector<ONode>> paths;
    };
    void enumerate(ONode sn, int32_t e, ONode tn, int32_t f, bool prune, Enumeration& en) const;
    std::string oseq(ONode n) const;
    int32_t len(uint32_t id) const { return static_cast<int32_t>(g_->nodes[id].seq.size()); }

    const UnitigGraph* g_ = nullptr;
    EvidenceOptions opt_;
    int k_ = 0;
    double theta_ = 0;
    bool abstained_ = false;
    const char* abstainWhy_ = "";
    std::vector<Ent> index_;                     // sorted by fp
    std::vector<Multiplicity> mult_;
    std::vector<char> placed_;
    std::unordered_map<uint64_t, uint32_t> midIndex_;   // hard single-copy middle 63-mers
};

// Canonical 63-mer fingerprints (the index's own hash), shared with the ledger's locator.
bool fingerprint63(const std::string& s, size_t from, uint64_t& fp);
void forEachFingerprint63(const std::string& s, const std::function<void(size_t, uint64_t)>& fn);

// Pair evidence at a seam, computed on the ORIGIN pieces (the contigs the model stage starts
// from). Every later piece end is an origin piece end, found by its terminal 250 bases.
struct PairCall {
    bool testable = false;   // origin ports known, N <= pairMaxN, index usable
    bool abstain = true;     // lambda below the floor or no control
    float lambda = 0;        // mean crossing count at the two control points
    float ctrlA = -1, ctrlB = -1;
    uint32_t k = 0;          // crossing fragments at the seam
    uint32_t contra = 0;     // fragments leaving the exit end to some other end
    float log10LR = 0;       // insertion vs adjacency
    int call = 0;            // +1 deny (insertion), -1 confirm (adjacent), 0 none
};

class PairIndex {
public:
    PairIndex(const std::vector<std::string>& origin, const SequenceStore& reads,
              const InsertModel& insert, int threads, const EvidenceOptions& opt);

    bool usable() const { return usable_; }
    // Origin port of an oriented piece's junction end, UINT32_MAX when unknown or ambiguous.
    uint32_t exitPort(const std::string& oriented) const;
    uint32_t entryPort(const std::string& oriented) const;
    PairCall test(uint32_t exitPort, uint32_t entryPort, int32_t claimedN) const;
    size_t linkWeight(uint32_t pieceA, uint32_t pieceB) const;
    size_t supportFloor() const;

    // The Poisson rule alone (seam_pairs2 / pilot_summary.py), for tests and the oracle.
    static PairCall callFrom(uint32_t k, float lambda, const EvidenceOptions& opt);

    // C4 (om2_clonal, read-only additions): every port the crossing fragments of `port` reach, with
    // their count (each fragment has both reads within the fragment bound of their piece ends), and
    // the length of an origin piece. Used to count read pairs that contradict a proposed adjacency.
    std::vector<std::pair<uint32_t, uint32_t>> partners(uint32_t port) const;
    size_t pieceLength(uint32_t piece) const { return piece < len_.size() ? len_[piece] : 0; }
    size_t pieces() const { return len_.size(); }

    size_t readsAnchored = 0, pairsSpanning = 0, pairsControl = 0;
    int32_t maxFrag = 0;

private:
    bool usable_ = false;
    EvidenceOptions opt_;
    std::vector<size_t> len_;
    std::vector<int32_t> ctrl_;                          // per port; -1 = unavailable
    std::vector<uint32_t> outTotal_;                     // per port
    std::unordered_map<uint64_t, std::vector<int32_t>> spans_;
    std::unordered_map<uint64_t, size_t> links_;
    std::unordered_map<std::string, uint32_t> portOf_;   // exit-oriented terminal 63 -> port
    double meanDepth_ = 0;
};

// Terminal-250 port map without reads (ledger provenance when pairs are off).
class PortMap {
public:
    explicit PortMap(const std::vector<std::string>& origin);
    uint32_t exitPort(const std::string& oriented) const;
    uint32_t entryPort(const std::string& oriented) const;
private:
    std::unordered_map<std::string, uint32_t> portOf_;
};

}  // namespace om2
}  // namespace ts
