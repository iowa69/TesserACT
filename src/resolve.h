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
};

class PairedResolver {
public:
    PairedResolver(const UnitigGraph& graph, const SequenceStore& reads, int threads,
                   int minLinkSupport, double tieRatio);

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

    const ResolveStats& stats() const { return stats_; }

    // How each emitted contig walks the graph, in output order.
    const std::vector<ResolvedPath>& paths() const { return paths_; }

    // Fragment-length histogram observed from same-unitig pairs.
    const std::vector<uint64_t>& insertHistogram() const { return insertHistogram_; }

private:
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
    std::vector<ResolvedPath> paths_;
    std::vector<uint64_t> insertHistogram_;
    ResolveStats stats_;
    double medianCoverage_ = 0;
    int minLinkSupport_;
    double tieRatio_;
    bool scaffolding_ = false;
};

}  // namespace ts
