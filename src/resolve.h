// Paired-end repeat resolution.
//
// The compacted graph stops wherever a repeat creates a branch, even when the
// correct way through is unambiguous given the read pairs. This stage anchors
// reads onto unitigs, learns the fragment-length distribution, and then walks
// paths through the graph choosing branches by paired-end support -- turning a
// graph of unitigs into a much smaller set of resolved contigs.
#pragma once

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

#include "graph.h"
#include "seqio.h"
#include "resolve_route_density.h"

namespace ts {

// Where a read sits on a unitig. `pos` is always the start of the read's
// footprint in the unitig's forward frame; `orient` says which way the read
// runs along it.
struct Anchor {
    uint32_t unitig = UINT32_MAX;
    // Signed, because a read overlapping a junction legitimately hangs off
    // either end of the unitig -- and those are exactly the informative reads.
    int32_t pos = 0;
    uint8_t orient = 0;   // 0 = read runs along the unitig's forward strand
    bool mapped() const { return unitig != UINT32_MAX; }
};

struct InsertModel {
    double mean = 0;
    double stddev = 0;
    int minPlausible = 0;
    int maxPlausible = 0;
    size_t observations = 0;
    bool usable = false;
};

// One output contig expressed as a walk over oriented unitigs, so the GFA and
// the report can show which graph nodes each contig actually traverses.
struct ResolvedPath {
    std::vector<uint64_t> oriented;   // (unitig << 1) | orientation
    std::vector<int> gaps;            // gaps[i] = N-gap inserted before element i
};

struct ResolveStats {
    size_t readsMapped = 0;
    size_t pairsLinking = 0;
    size_t distinctLinks = 0;
    size_t pathsBuilt = 0;
    size_t unitigsJoined = 0;
    size_t scaffoldJoins = 0;
    size_t gapBases = 0;
    InsertModel insert;
    // T35: the single-copy depth every repeat decision is taken against, the repeat
    // threshold derived from it, and how it was estimated. `legacyMedian` is the unweighted
    // median over live nodes >= 2k (UnitigGraph::medianCoverage, which is what report.json
    // iterations[].median_coverage shows); `thetaPopulation` counts the nodes the chosen
    // estimator weighed. Filled by the constructor; not yet written to report.json.
    double theta = 0;
    double repeatThreshold = 0;
    double legacyMedian = 0;
    const char* thetaEstimator = "unweighted";
    size_t thetaPopulation = 0;
};

// T03 (TESSERACT_FIX_GAP_FLANK): how one scaffold gap must look if it is still OPEN after
// gap filling. The resolver writes the release layout (`nWritten` Ns standing in for the
// first k-1 bases of the unitig after the gap), which is what the gap filler needs: its
// target k-mer then lies in full-depth sequence, and a gap it closes gets those bases back
// from the reads. A gap it leaves open gets them back from here: the N-run becomes `nOpen`
// (the estimated true gap; 1 for a verified overlap; kUnknownGapN for an unknown length)
// followed by `restore`. `left` / `right` are the flank bases around the N-run (in the
// orientation the resolver wrote it) that identify the gap after gap filling.
struct GapFlankRecord {
    std::string left, right;
    int nWritten = 0;
    int nOpen = 0;
    std::string restore;
};

// Applies the records to the scaffolds after gap filling: an N-run of exactly `nWritten`
// between `left` and `right` (either orientation) becomes `nOpen` Ns plus `restore`.
// A gap the filler closed, or one whose flanks changed, is left alone. Prints
// `[gapflank-restore]` on every call, zeros included.
void restoreGapFlanks(std::vector<std::string>& seqs, const std::vector<GapFlankRecord>& records);

// Prints every always-on resolver counter line ([resolver], [resolveflags], [enumtrunc],
// [revisit], [scafcycle], [gapest], [gapflank], [cov_contrib], [routeorder], [mirrorroute])
// with run=0 and
// zero counts. For runs in which the paired resolver is not constructed (single-end input,
// --no-resolve), so that each line appears on every run (OBJECTIVE amendment A2).
void printResolverCountersNotRun();

class PairedResolver {
public:
    // Markers the organism model says occur only on panel plasmids (1) or only on panel
    // chromosomes (2), keyed by canonical k-mer. Optional; when absent every decision
    // below behaves exactly as before.
    //
    // isRepeat() below is depth-only, and cannot tell a two-copy repeat from a plasmid at
    // two copies per cell -- the resolve.cpp comment on repeatThreshold says so, and puts
    // the cost at 4.8% of one assembly. Measured on 119 reference plasmids across 93
    // S. aureus isolates, the cost is concentrated exactly where that ambiguity bites:
    // below 1.6x chromosome depth TesserACT recovers 37% of plasmids whole against
    // SPAdes' 33%, and at or above 8x it recovers 7% against SPAdes' 60%. The model
    // carries the signal depth does not.
    void setExclusiveMarkers(const std::unordered_map<uint64_t, uint8_t>* m, uint32_t denom) {
        exclusiveMarkers_ = m;
        markerDenom_ = denom ? denom : 1;
    }
    PairedResolver(const UnitigGraph& graph, const SequenceStore& reads, int threads,
                   int minLinkSupport, double tieRatio, double linkSupportPerX = 0.10,
                   int minScaffoldSupport = 0);

    // Anchors reads and accumulates the oriented link support table.
    void buildSupport();

    // Fragment-length distribution, learned from pairs landing on one unitig.
    const InsertModel& insertModel() const { return insert_; }

    // Walks paths through the graph, resolving branches with paired support.
    // Returns the resolved contig sequences and their mean coverages.
    void resolve(std::vector<std::string>& contigs, std::vector<double>& covs);

    // When enabled, contig ends with paired support but no path through the
    // graph are joined across a gap of Ns sized from the fragment model.
    void setScaffolding(bool on) { scaffolding_ = on; }

    // Supply the plausible-length window from an outside measurement rather than fitting
    // it here. buildSupport() still learns the mean and spread from anchored pairs -- both
    // are reported and remain useful -- but these BOUNDS override the fitted ones.
    //
    // The bounds are the part worth replacing. Fitted here they are mean +/- 4 sd against
    // an sd taken from a sample whose tails survive only a 1% trim, which on a typical
    // library in this cohort works out to about [0, 3.3x mean]: every candidate distance
    // falls inside, and the plausibility term in scoreCandidate stops discriminating at
    // all. Percentiles of a histogram built on the raw reads cannot degenerate that way,
    // and they are measured without reference to the graph they are used to improve.
    void setInsertBounds(int minPlausible, int maxPlausible) {
        forcedMin_ = minPlausible;
        forcedMax_ = maxPlausible;
        haveForcedBounds_ = maxPlausible > minPlausible;
    }

    const ResolveStats& stats() const { return stats_; }

    // How each emitted contig walks the graph, in output order.
    const std::vector<ResolvedPath>& paths() const { return paths_; }

    // T03: one record per scaffold gap written by resolve() (empty unless
    // TESSERACT_FIX_GAP_FLANK is on). To be applied by restoreGapFlanks() after gap filling.
    const std::vector<GapFlankRecord>& gapFlankRecords() const { return gapFlank_; }

    // Fragment-length histogram observed from same-unitig pairs.
    const std::vector<uint64_t>& insertHistogram() const { return insertHistogram_; }

private:
    const std::unordered_map<uint64_t, uint8_t>* exclusiveMarkers_ = nullptr;
    uint32_t markerDenom_ = 1;

    // Oriented unitig identity, packed as (unitig << 1 | orientation).
    static uint64_t orientedId(uint32_t u, int d) {
        return (static_cast<uint64_t>(u) << 1) | static_cast<uint64_t>(d & 1);
    }
    static uint32_t unitigOf(uint64_t oid) { return static_cast<uint32_t>(oid >> 1); }
    static int orientOf(uint64_t oid) { return static_cast<int>(oid & 1); }

    void buildIndex();
    Anchor anchorRead(size_t read) const;

    // Paired support for continuing from `from` along a candidate path whose
    // intermediate nodes add `interLen` bases. A pair only counts when the
    // fragment length it implies for *this* path is plausible, which is what
    // separates a real continuation from a coincidental link.
    double scoreCandidate(uint64_t from, uint64_t terminal, int interLen) const;

    const UnitigGraph& g_;
    const SequenceStore& reads_;
    int threads_;
    int k_;
    // Anchoring uses a shorter k than the graph: at k=63 a single sequencing
    // error invalidates nearly every k-mer in a 150 bp read, so most reads --
    // including the ones spanning junctions -- would fail to anchor at all.
    int kMap_;

    // Anchor k-mer -> (unitig, position, strand flag), packed into 64 bits.
    // A k-mer occurring in more than one place is stored as `kAmbiguous` and
    // ignored, since it cannot identify a location.
    std::unordered_map<Kmer, uint64_t, KmerHasher> index_;

    // orientedUnitig -> orientedUnitig -> the fragment span each supporting
    // pair implies, excluding whatever sequence lies between the two unitigs.
    using SpanList = std::vector<int32_t>;
    std::unordered_map<uint64_t, std::unordered_map<uint64_t, SpanList>> support_;

    InsertModel insert_;
    // combo2/ends: every read's anchor, kept by buildSupport() only when
    // TESSERACT_PAIR_ANCHORED_PREFIX is set (12 bytes per read); empty otherwise.
    std::vector<Anchor> readAnchors_;
    std::vector<ResolvedPath> paths_;
    std::vector<GapFlankRecord> gapFlank_;
    std::vector<uint64_t> insertHistogram_;
    // Default-off: conditional route lengths only; endpoint count scores stay intact.
    bool routeDistance_ = false;
    detail::RouteInsertDensity routeDensity_;
    ResolveStats stats_;
    double medianCoverage_ = 0;
    int minLinkSupport_;
    double tieRatio_;
    double linkSupportPerX_;
    int minScaffoldSupport_ = 0;   // 0 = scale with observed depth
    bool scaffolding_ = false;
    int forcedMin_ = 0, forcedMax_ = 0;
    bool haveForcedBounds_ = false;
};

}  // namespace ts
