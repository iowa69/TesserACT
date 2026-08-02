// Joining contigs on what the panel says about this organism.
//
// This stage runs after repeat resolution, on exactly the junctions paired
// reads could not reach: the ones where the repeat is longer than a fragment,
// so there is no pair to weigh and the resolver correctly stopped. Rather than
// guess from coverage, it asks whether closed genomes of the same organism
// agree that these two flanks belong together, and at what distance.
//
// The join is taken only when several independent marker pairs agree, a
// substantial share of the panel genomes carrying both markers support it, and
// the implied gap is consistent across them. Everything else is left broken.
#pragma once

#include <cstddef>
#include <string>
#include <unordered_set>
#include <vector>

#include "organism.h"

namespace ts {

struct OrganismJoinStats {
    size_t contigsIn = 0;
    size_t contigsOut = 0;
    size_t markersFound = 0;    // model markers located uniquely in the assembly
    size_t endsWithMarkers = 0;
    size_t candidates = 0;      // end pairs with any panel support at all
    size_t joins = 0;
    size_t chromosomeJoins = 0;
    size_t plasmidJoins = 0;
    size_t gapBases = 0;
    size_t rejectedWeak = 0;    // support below the panel threshold
    size_t rejectedNotMutual = 0;
    size_t rejectedInconsistent = 0;   // no agreeing gap cluster
    size_t rejectedInsertionSeq = 0;   // an end sitting inside a mobile element
};

// Canonical 31-mers of known insertion sequences. A contig that ends inside a
// mobile element is ending there *because* this isolate carries an IS the panel
// may not, and the panel's adjacency across that point describes a genome
// without the insertion. Joining on it deletes the element and misassembles.
class IsPanel {
public:
    bool load(const std::string& fastaPath, std::string& error);
    bool loaded() const { return !kmers_.empty(); }
    size_t size() const { return kmers_.size(); }
    size_t copies() const { return copies_; }
    // Fraction of sampled positions in `window` that fall in a known IS.
    double density(const std::string& seq, size_t from, size_t len) const;

private:
    std::unordered_set<uint64_t> kmers_;
    size_t copies_ = 0;
};

// Joins `contigs` in place using the model. `covs` is reordered to match.
// `k` is the assembly's k, used only to bound a negative gap. When `isPanel` is
// supplied, ends sitting inside a mobile element are left alone.
OrganismJoinStats joinByModel(const OrganismModel& model, std::vector<std::string>& contigs,
                              std::vector<double>& covs, int k, bool verbose,
                              const IsPanel* isPanel = nullptr);

}  // namespace ts
