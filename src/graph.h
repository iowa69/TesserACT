// The compacted (unitig) de Bruijn graph and its simplification operations.
//
// The graph is bidirected: a unitig is a double-stranded sequence with two
// ends, and a link records which end of the neighbour it attaches to. Walking
// out of end 1 of u into a link {v, 0} means continuing through v forward;
// {v, 1} means continuing through v reverse-complemented. Every link is stored
// on both sides, so `u.ends[eu]` containing {v, ev} implies `v.ends[ev]`
// contains {u, eu}. Every operation below preserves that invariant.
#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "counter.h"
#include "kmer.h"

namespace ts {

struct Link {
    uint32_t to = 0;
    uint8_t toEnd = 0;
    bool operator==(const Link& o) const { return to == o.to && toEnd == o.toEnd; }
};

struct Unitig {
    std::string seq;
    double coverage = 0;      // mean k-mer multiplicity along the unitig
    bool deleted = false;
    std::vector<Link> ends[2];   // ends[0] = 5' side, ends[1] = 3' side

    size_t length() const { return seq.size(); }
    size_t degree(int end) const { return ends[end].size(); }
};

std::string reverseComplement(const std::string& s);

// What one pass of the simplification schedule changed, for the run report.
struct SimplifyRoundStats {
    int round = 0;
    size_t tipsRemoved = 0;
    size_t bubblesPopped = 0;
    size_t chimerasRemoved = 0;
    size_t isolatedRemoved = 0;
    size_t lowDepthRemoved = 0;
    size_t deadEnds = 0;
    size_t merged = 0;
    size_t unitigs = 0;
    size_t n50 = 0;
    size_t totalLength = 0;
};

class UnitigGraph {
public:
    UnitigGraph() = default;

    // Builds the compacted graph from the solid k-mer set.
    static UnitigGraph build(const KmerTable& solid, int k, int threads);

    int k() const { return k_; }
    void setK(int k) { k_ = k; }

    std::vector<Unitig> nodes;

    // ---- structural maintenance ----
    // Merges maximal non-branching chains into single unitigs.
    size_t compact();
    // Physically drops deleted nodes and renumbers links.
    void removeDeleted();
    // Drops a node and unlinks it from its neighbours.
    void deleteNode(uint32_t u);
    // Verifies the bidirected invariant; returns a description of the first
    // violation found, or an empty string. Used by the test suite.
    std::string validate() const;

    // ---- simplification ----
    // Removes short dead-end branches whose coverage is low relative to the
    // alternative at the branch point.
    size_t removeTips(size_t maxLen, double covRatio);
    // Collapses pairs of near-identical parallel paths, keeping the better
    // supported one. `maxLoserCoverage` bounds how well covered the discarded
    // side may be: an error bubble has a clearly weak side, whereas two
    // genuinely diverged repeat copies carry comparable depth and collapsing
    // them would fuse distinct loci into one misassembled contig.
    size_t popBubbles(size_t maxLen, double minIdentity, double maxLoserCoverage = 0);
    // Drops low-coverage unitigs that connect otherwise well-supported regions
    // (the classic erroneous-connection pattern).
    size_t removeErroneousConnections(double covThreshold, size_t maxLen);
    // Drops isolated short low-coverage fragments.
    size_t removeIsolated(double covThreshold, size_t maxLen);

    // How many of a unitig's two ends lead nowhere.
    int deadEnds(uint32_t u) const;
    // Net change in the graph's dead-end count if `u` were deleted: removing a
    // unitig closes its own loose ends but strands any neighbour that was
    // attached only through it. Negative means deleting tidies the graph.
    int deadEndChangeIfDeleted(uint32_t u) const;
    size_t totalDeadEnds() const;

    // Connected components, as lists of live unitig ids.
    std::vector<std::vector<uint32_t>> components() const;

    // Length-weighted median coverage of the longest `topN` unitigs in `ids`.
    // Anchoring on the longest rather than on every unitig is what stops a
    // swarm of short low-depth fragments from dragging the estimate down.
    double weightedMedianCoverage(const std::vector<uint32_t>& ids, size_t topN = 10) const;

    // Removes unitigs whose coverage falls below `fraction` of the local
    // median, judged both against the whole graph and against their own
    // connected component -- a low-copy plasmid is well below the global
    // median without being spurious, and only the component-relative test
    // tells the two apart. A unitig is kept regardless if deleting it would
    // strand its neighbours.
    size_t filterByReadDepth(double fraction);
    // Runs the full simplification schedule until it converges.
    // `bubbleCoverageLimit` is the fraction of mean coverage below which a
    // bubble side may be discarded; raising it collapses diverged repeat
    // copies too, which buys contiguity at the risk of misassembly.
    void simplify(double meanCoverage, int readLength, bool verbose,
                  double bubbleCoverageLimit = 0.35, int maxRounds = 12,
                  std::vector<SimplifyRoundStats>* rounds = nullptr);

    // ---- traversal helpers ----
    // Sequence of `u` in the given orientation (0 = forward, 1 = reverse).
    std::string oriented(uint32_t u, int orient) const;
    // Links leaving `u` when travelling in `orient`.
    const std::vector<Link>& exits(uint32_t u, int orient) const {
        return nodes[u].ends[orient == 0 ? 1 : 0];
    }
    // Orientation adopted after following a link.
    static int enterOrient(const Link& l) { return l.toEnd == 0 ? 0 : 1; }

    // ---- reporting ----
    size_t liveCount() const;
    size_t totalLength() const;
    size_t n50() const;
    double medianCoverage() const;
    void stats(const char* label, bool verbose) const;

    // Emits every unitig longer than `minLen` as a contig.
    void toContigs(size_t minLen, std::vector<std::string>& seqs,
                   std::vector<double>& covs) const;

private:
    void addLink(uint32_t u, int ue, uint32_t v, int ve);
    void unlink(uint32_t u, int ue, uint32_t v, int ve);
    bool mergeInto(uint32_t u, int ue);

    int k_ = 0;
};

// Sequence identity between two strings, computed with banded edit distance.
// Returns a value in [0,1]. Used to decide whether two paths form a bubble.
double sequenceIdentity(const std::string& a, const std::string& b, int maxBand);

}  // namespace ts
