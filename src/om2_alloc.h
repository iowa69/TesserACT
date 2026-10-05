// Organism Model 2.0, component C2: labelled allocation of repeat variants to loci, and the two
// leave-clone-out sidecars C2 reads (rRNA locus prior, dnaA locator sketch).
//
// A collapsed repeat (an rRNA operon, an IS family) appears in the graph once, with a bubble at
// every site where its copies differ. When the stage in om2_close.cpp bridges several loci
// through the same collapsed repeat, each locus needs one branch at each site. This module
// decides which, and says on what basis (DESIGN.md section 3.3):
//
//   PAIR_PHASED / THREAD_PHASED  flank-anchored read pairs / single reads link this locus's
//                                unique flank to the branch; posterior >= 0.99
//   PRIOR_ALLOCATED              ITS-class site only: the leave-clone-out locus prior names the ITS
//                                class (P >= 0.99, >= 10 genomes), a branch of that class is left in
//                                the budget, and the class chosen FLIPS when the prior is removed
//   MULTIPLICITY                 the branch budget (carrier counts from depth) drives the choice
//   CONSENSUS                    the isolate's majority branch
//
// Phased choices are made first and consume the budget; the prior is consulted next; every
// remaining locus takes the highest-marginal branch if its marginal is >= 0.8, else the majority
// branch. By construction the expected per-site accuracy is never below consensus paste.
//
// A prior never decides bases: it chooses among branches that exist in THIS isolate's graph.
#pragma once

#include <cstdint>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "om2_types.h"

namespace ts {
namespace om2 {

const char* c2BasisName(Basis b);

enum class AllocMode : uint8_t { Off, Consensus, Phased, Prior };

struct AllocParams {
    AllocMode mode = AllocMode::Phased;
    double phaseMin = 0.99;       // posterior for PAIR_/THREAD_PHASED
    double marginalMin = 0.8;     // budget marginal for the highest-marginal rule
    double eps = 0.02;            // probability that a linking fragment names the wrong branch
    double priorMin = 0.99;       // prior probability of the named branch
    uint32_t priorGenomesMin = 10;
};

struct SiteBranch {
    int carriers = 1;             // copies carrying this branch, from depth (the budget)
    double cov = 0;               // depth, for the majority rule's tie-break
};

struct SiteLocus {
    uint32_t junction = 0;
    std::vector<uint32_t> threadN, pairN;   // per branch: fragments linking this locus to it
    std::vector<double> prior;              // per ITS CLASS (see branchClass); empty = no usable prior
    uint32_t priorGenomes = 0;
};

struct SiteInput {
    bool itsClass = false;        // branch lengths differ by >= 50 bp: the prior may be consulted
    bool complex = false;         // too many paths to enumerate: one greedy branch, consensus only
    std::vector<SiteBranch> branches;
    // ITS class of each branch (-1 unknown). Real ITS bubbles carry several branches per class (the
    // ITS types times their own SNVs), so the prior names a CLASS and the budget picks the branch
    // within it. Empty: branch b is class b.
    std::vector<int> branchClass;
    std::vector<SiteLocus> loci;
};

struct SiteChoice {
    int branch = -1;              // -1: no allocation (mode off)
    Basis basis = Basis::Consensus;
    double posterior = 0;         // phased posterior, prior probability, or budget marginal
    bool priorFlip = false;       // PRIOR_ALLOCATED: the no-prior choice differs
    int noPriorBranch = -1;       // the counterfactual choice (for the report)
};

struct AllocCounts {
    size_t pair = 0, thread = 0, prior = 0, multiplicity = 0, consensus = 0;
    size_t priorNoFlip = 0;       // prior consulted, agreed with the no-prior choice
    size_t budgetExhausted = 0;   // every branch budget spent; majority used anyway
};

// The isolate's majority branch: most carriers, then most depth, then lowest index.
int majorityBranch(const SiteInput& s);

// Posterior of branch `b` given per-branch linking counts, under a uniform prior and a
// symmetric error model (a fragment names the true branch with probability 1 - eps).
double phasePosterior(const std::vector<uint32_t>& n, size_t b, double eps);

// Allocates one branch to every locus of one site. Deterministic: loci are visited in the
// order given (phased first, then prior by strength, then the rest in order).
std::vector<SiteChoice> allocateSite(const SiteInput& s, const AllocParams& p, AllocCounts& counts);

// ---- sidecars --------------------------------------------------------------------------------

// MD5 of a file as 32 lowercase hex digits; empty on a read error. Sidecars pin the .tsm they
// were built beside by its md5, and the loaders refuse a mismatch.
std::string md5File(const std::string& path);

// Canonical 31-mer hash shared with the sidecar builders: 2 bits per base (A0 C1 G2 T3), first
// base most significant, canonical = min(forward, reverse complement), then splitmix64.
uint64_t sidecarHash(uint64_t canonical31);

// <org>.om2rrn: per rRNA locus of the leave-clone-out panel, the flank marker hashes that name
// it and the ITS class counts over panel genomes; per ITS class, the hashes of 31-mers that
// discriminate it. Hashes and counts only: no sequence.
class RrnPrior {
public:
    bool load(const std::string& path, const std::string& modelMd5, std::string& error);
    bool loaded() const { return loaded_; }
    uint64_t flankDenom() const { return denom_; }
    int flankBp() const { return flankBp_; }
    // Best locus for a set of flank marker hashes: -1 unless >= 3 shared and >= 2x the runner-up.
    int matchLocus(const std::vector<uint64_t>& markers, int* score = nullptr) const;
    // ITS class of a branch sequence: -1 unless >= 3 discriminating hits and >= 2x the runner-up.
    int classOf(const std::string& seq) const;
    // Genomes at `locus` and the count of those carrying `cls`.
    uint32_t locusGenomes(int locus) const;
    uint32_t classCount(int locus, int cls) const;
    size_t loci() const { return loci_.size(); }
    size_t classes() const { return classes_.size(); }

private:
    struct Locus { uint32_t genomes = 0; std::unordered_map<int, uint32_t> counts; };
    bool loaded_ = false;
    uint64_t denom_ = 64;
    int flankBp_ = 4000;
    std::vector<std::unordered_set<uint64_t>> classes_;
    std::vector<Locus> loci_;
    std::unordered_map<uint64_t, std::vector<uint32_t>> markerLoci_;
};

// <org>.om2dnaa: hashes of the 31-mers of the first 1.5 kb of leave-clone-out panel chromosomes
// (present in >= 50% of them), with their offset from the chromosome start. A locator only:
// it says where dnaA begins in a closed record; it never contributes a base.
class DnaaSketch {
public:
    bool load(const std::string& path, const std::string& modelMd5, std::string& error);
    bool loaded() const { return loaded_; }
    struct Hit { bool found = false; int64_t offset = -1; bool reverse = false; uint32_t votes = 0, hits = 0; };
    Hit locate(const std::string& circularRecord) const;
    size_t size() const { return entries_.size(); }

private:
    struct Entry { int32_t offset = 0; uint8_t strand = 0; };
    bool loaded_ = false;
    std::unordered_map<uint64_t, Entry> entries_;
};

}  // namespace om2
}  // namespace ts
