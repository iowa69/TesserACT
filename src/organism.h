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
// Overridable at compile time only, so the density can be measured. Sampling is by hash
// threshold, so a denser setting is a strict superset of a sparser one and the two are
// directly comparable. Build and query must agree: a model carries no record of the
// density it was built at, and a mismatched pair silently shares almost no markers.
#ifndef TS_MARKER_SAMPLE_DENOM
#define TS_MARKER_SAMPLE_DENOM 512
#endif
constexpr uint64_t kMarkerSampleDenom = TS_MARKER_SAMPLE_DENOM;
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

// One panel chromosome's marker order, kept whole.
//
// The adjacency table above answers "do these two flanks belong together". That is a
// local question, and closing a chromosome means answering about forty of them in a row
// without a single mistake -- one miss leaves the chromosome in two pieces. Measured on
// this cohort, per-junction success would have to reach about 99.7% for nine genomes in
// ten to close, against 83.9% actually achieved.
//
// A track answers a different, global question: where does this contig sit on a genome
// shaped like this one? Ordering every contig against one relative replaces forty
// independent decisions with a single alignment, and a locus the relative gets wrong
// costs that locus rather than the whole chromosome.
//
// The track holds markers, not sequence: positions along a genome the assembler will
// never see. So it needs no aligner and no panel FASTA at run time, and it cannot leak
// sequence from the panel into the output.
struct LayoutTrack {
    std::string name;                    // panel accession, for provenance
    std::vector<uint32_t> oriented;      // (marker id << 1) | orientation, in order
    std::vector<uint32_t> pos;           // start position on that chromosome
};

// One panel plasmid's marker CONTENT, without any order.
//
// The adjacency table answers "does marker A follow marker B", which is an order question,
// and order is exactly what plasmids do not conserve -- they recombine and carry the same
// transposons in different arrangements. Membership asks instead "do these markers occur
// together on one real plasmid", which mosaicism does not disturb.
//
// Measured before this was built, on 5,370 contig pairs with known truth: the share of a
// contig pair's markers explained by a single panel plasmid separates same-plasmid from
// different-plasmid pairs at AUC 0.857, and under a degree-preserving shuffle of the
// marker-plasmid graph that falls to 0.416 -- below chance. So the signal is co-residence
// and not the ubiquity of IS elements and backbone genes. It is the only one of the three
// candidates that works: depth leaves 43% of same-isolate plasmid pairs within 1.5x of
// each other, and read pairs cannot span the repeats that separate plasmid contigs.
struct PlasmidMembership {
    std::string name;                 // panel accession, for provenance
    uint32_t length = 0;
    std::vector<uint32_t> markers;    // marker ids, sorted, deduplicated
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

    // The sampling denominator this model was built at. A model carries it because build
    // and query must agree exactly: sampling is by hash threshold, so a query at 1/512
    // against a model built at 1/64 sees only the eighth of the markers the two happen to
    // share and silently reports almost nothing. Models written before version 5 predate
    // the field and were all built at 512.
    uint32_t markerDenom() const { return markerDenom_; }

    // ---- construction (used by the model builder) ----
    void beginBuild(const std::string& organism, int k);
    void setMarkerDenom(uint32_t d) { if (d >= 1) markerDenom_ = d; }
    uint32_t internMarker(uint64_t canonicalKmer);
    void noteMarkerGenome(uint32_t id, Replicon r) { ++genomesPerMarker_[static_cast<int>(r)][id]; }
    void addObservation(uint64_t fromOriented, uint64_t toOriented, int32_t dist, Replicon r);
    void noteGenome(Replicon r) {
        if (r == Replicon::Chromosome) ++genomesChr_; else ++genomesPls_;
    }
    void noteExcluded(const std::string& acc) { excluded_.push_back(acc); }
    void addTrack(LayoutTrack&& t) { tracks_.push_back(std::move(t)); }
    void addPlasmid(PlasmidMembership&& m) { plasmids_.push_back(std::move(m)); }

    const std::vector<PlasmidMembership>& plasmids() const { return plasmids_; }
    size_t plasmidCount() const { return plasmids_.size(); }

    const std::vector<LayoutTrack>& tracks() const { return tracks_; }
    size_t trackCount() const { return tracks_.size(); }
    // Collapses the per-observation distance lists into medians and drops
    // pairs seen in too few panel members to mean anything.
    void finalise(uint32_t minSupportChr, uint32_t minSupportPls);
    bool save(const std::string& path, std::string& error) const;

private:
    bool loaded_ = false;
    std::string organism_;
    int k_ = kMarkerK;
    uint32_t markerDenom_ = kMarkerSampleDenom;
    uint32_t genomesChr_ = 0, genomesPls_ = 0;
    std::vector<std::string> excluded_;

    std::unordered_map<uint64_t, uint32_t> index_;   // canonical k-mer -> marker id
    std::vector<uint32_t> genomesPerMarker_[2];
    std::unordered_map<uint64_t, MarkerEdge> edges_[2];
    std::vector<LayoutTrack> tracks_;
    std::vector<PlasmidMembership> plasmids_;

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
// `denom` is the sampling denominator; it must be the one the model being queried was built
// at, which the model reports through markerDenom(). It is a parameter rather than the
// constant because a model file records its own density, and a query at a different one
// shares almost no markers with it while failing silently.
template <typename F>
void forEachMarkerKmer(const std::string& seq, F&& fn, uint64_t denom = kMarkerSampleDenom) {
    const int k = kMarkerK;
    const uint64_t mask = (1ULL << (2 * k)) - 1;
    const uint64_t threshold = ~0ULL / (denom ? denom : kMarkerSampleDenom);
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
