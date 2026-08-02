// A genus-specific prior, learned from closed genomes of the same organism.
//
// Paired reads settle a junction only when a fragment can span it. On the
// Klebsiella panel the fragments run 217 +/- 110 bp while the genome carries
// 177 repeat copies of 500 bp or more, so at those junctions there is no
// paired evidence to weigh -- not weak evidence, none. An assembler can either
// stop there, or guess from coverage and accept the misassemblies that brings.
//
// There is a third option. Klebsiella pneumoniae genomes are not arbitrary
// strings: the sequence flanking a repeat is largely conserved across the
// species, so a panel of closed genomes can say which flank normally follows
// which. This model records that. It samples canonical k-mers at a fixed
// density, keeps those that occur exactly once in a replicon set (so they name
// a locus rather than a repeat family), and stores how often each ordered pair
// of them is observed in the panel and at what distance.
//
// Chromosome and plasmid are kept apart, and that separation is the point
// rather than tidiness. Core gene order is strongly conserved along the
// chromosome, so adjacency there is close to a rule. Plasmids are mosaic:
// they recombine, carry the same transposons in different orders, and vary in
// copy number, so an adjacency observed on one plasmid is far weaker evidence
// about another. Pooling the two would let plasmid rearrangements vote on
// chromosomal junctions. Each table is therefore built, queried and thresholded
// on its own, in two passes -- the chromosome first, then the plasmids.
//
// The model must never be built from the genome being assembled. Every
// evaluation here excludes the isolate under test by accession; see
// OrganismModel::excluded.
#pragma once

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace ts {

// Sampling density: a marker is kept when its hash falls below this fraction
// of the 64-bit range. One in 512 positions gives ~500 bp spacing, so a 20 kb
// window near a contig end holds ~40 markers -- enough that a handful of
// missing ones does not silence the query.
constexpr uint64_t kMarkerSampleDenom = 512;
constexpr int kMarkerK = 31;

// How many downstream markers each marker links to. Sixteen at ~500 bp spacing
// reaches ~8 kb, which covers the rRNA operons and IS elements that defeat the
// paired evidence, without inflating the model with links no query will use.
constexpr int kMarkerNeighbours = 16;
constexpr int32_t kMaxMarkerDistance = 30000;

enum class Replicon { Chromosome = 0, Plasmid = 1 };

struct MarkerEdge {
    uint32_t support = 0;      // panel replicon sets observing this ordered pair
    int32_t medianDist = 0;    // median start-to-start distance, in bases
};

class OrganismModel {
public:
    bool load(const std::string& path, std::string& error);

    bool loaded() const { return loaded_; }
    uint32_t genomes(Replicon r) const {
        return r == Replicon::Chromosome ? genomesChr_ : genomesPls_;
    }
    size_t markerCount() const { return genomesPerMarker_[0].size(); }
    size_t edgeCount(Replicon r) const { return edges_[static_cast<int>(r)].size(); }
    const std::string& organism() const { return organism_; }
    // Accessions deliberately withheld when the model was built, recorded so a
    // run can prove it never saw the genome it is being scored against.
    const std::vector<std::string>& excluded() const { return excluded_; }

    // Marker id for a canonical k-mer, or UINT32_MAX.
    uint32_t markerOf(uint64_t canonicalKmer) const {
        auto it = index_.find(canonicalKmer);
        return it == index_.end() ? UINT32_MAX : it->second;
    }

    // How many panel replicon sets of this class carry the marker exactly once.
    uint32_t markerGenomes(uint32_t id, Replicon r) const {
        return genomesPerMarker_[static_cast<int>(r)][id];
    }

    // Support for walking from `from` to `to`, both oriented as (id << 1 | orient).
    const MarkerEdge* edge(uint64_t fromOriented, uint64_t toOriented, Replicon r) const {
        const auto& tbl = edges_[static_cast<int>(r)];
        auto it = tbl.find((fromOriented << 32) | toOriented);
        return it == tbl.end() ? nullptr : &it->second;
    }

    // ---- construction (used by the model builder) ----
    void beginBuild(const std::string& organism, int k);
    uint32_t internMarker(uint64_t canonicalKmer);
    void noteMarkerGenome(uint32_t id, Replicon r) { ++genomesPerMarker_[static_cast<int>(r)][id]; }
    void addObservation(uint64_t fromOriented, uint64_t toOriented, int32_t dist, Replicon r);
    void noteGenome(Replicon r) {
        if (r == Replicon::Chromosome) ++genomesChr_; else ++genomesPls_;
    }
    void noteExcluded(const std::string& acc) { excluded_.push_back(acc); }
    // Collapses the per-observation distance lists into medians and drops
    // pairs seen in too few panel members to mean anything.
    void finalise(uint32_t minSupportChr, uint32_t minSupportPls);
    bool save(const std::string& path, std::string& error) const;

private:
    bool loaded_ = false;
    std::string organism_;
    int k_ = kMarkerK;
    uint32_t genomesChr_ = 0, genomesPls_ = 0;
    std::vector<std::string> excluded_;

    std::unordered_map<uint64_t, uint32_t> index_;   // canonical k-mer -> marker id
    std::vector<uint32_t> genomesPerMarker_[2];
    std::unordered_map<uint64_t, MarkerEdge> edges_[2];

    // Only populated while building. A reservoir of a few distances rather
    // than every one: the median of eight samples is well inside the 500 bp
    // tolerance the query applies, and a heap-allocated vector per pair costs
    // more memory than the whole rest of the model at panel sizes in the
    // thousands.
    static constexpr int kDistSamples = 8;
    struct Pending {
        uint32_t count = 0;
        int32_t dist[kDistSamples] = {0, 0, 0, 0, 0, 0, 0, 0};
    };
    std::unordered_map<uint64_t, Pending> pending_[2];
};

// Hash used both to sample markers and to key the tables. Splitmix64.
inline uint64_t markerHash(uint64_t x) {
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}

inline int markerBaseCode(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return -1;
    }
}

// Rolling canonical 31-mer scan. `fn(canonicalKmer, position, orientation)`
// fires for every position whose k-mer is sampled as a marker; orientation is
// 0 when the forward k-mer is the canonical one. Sampling by hash rather than
// by position means the same loci are picked in every genome and in the
// assembly, which is what makes the panel's observations comparable.
template <typename F>
void forEachMarkerKmer(const std::string& seq, F&& fn) {
    const int k = kMarkerK;
    const uint64_t mask = (1ULL << (2 * k)) - 1;
    const uint64_t threshold = ~0ULL / kMarkerSampleDenom;
    uint64_t fwd = 0, rev = 0;
    int valid = 0;
    for (uint32_t i = 0; i < seq.size(); ++i) {
        const int c = markerBaseCode(seq[i]);
        if (c < 0) { valid = 0; fwd = rev = 0; continue; }
        fwd = ((fwd << 2) | static_cast<uint64_t>(c)) & mask;
        rev = (rev >> 2) | (static_cast<uint64_t>(3 - c) << (2 * (k - 1)));
        if (++valid < k) continue;
        const bool fwdCanon = fwd <= rev;
        const uint64_t canon = fwdCanon ? fwd : rev;
        if (markerHash(canon) > threshold) continue;
        fn(canon, i + 1 - static_cast<uint32_t>(k), fwdCanon ? 0 : 1);
    }
}

}  // namespace ts
