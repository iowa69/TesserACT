// Deciding what molecule each contig belongs to, and saying so in the output.
//
// A pile of contigs is not what a clinical read-out needs. It needs "this is the
// chromosome, these are plasmid 1 and plasmid 2, and this is what did not fit anywhere" --
// with the evidence for each call attached, so a reader can tell an observation from an
// inference.
//
// Four independent signals, none of which is sufficient alone:
//
//   LAYOUT     a contig the panel layout placed on a chromosome track is chromosomal. This
//              is the strongest signal and covers most of the sequence, but only after the
//              chromosome closes.
//   MARKERS    the model already counts, per contig, how many of its markers are seen on
//              panel chromosomes versus panel plasmids. Decisive when a contig carries
//              markers at all; silent for the many short ones that do not.
//   DEPTH      the chromosome is single copy. A contig at several times the modal depth is
//              a multi-copy plasmid, and one far below it is contamination or a fragment.
//              Free, already computed, and the only signal that works on a contig with no
//              recognisable sequence at all.
//   PAIRS      contigs on one molecule are linked by read pairs; contigs on different
//              molecules should not be. This is what attaches the single-copy plasmid
//              contigs that depth cannot separate from the chromosome.
//
// Circularity is reported separately rather than used to classify: a circular contig is
// more likely a plasmid, but a closed chromosome is circular too, so it identifies a
// finished molecule rather than which molecule it is.
#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "organism.h"
#include "pairends.h"

namespace ts {

enum class RepliconClass : uint8_t {
    Chromosome = 0,
    Plasmid = 1,
    Unassigned = 2,   // no signal reached it; reported as-is rather than guessed at
};

struct RepliconCall {
    RepliconClass cls = RepliconClass::Unassigned;
    uint32_t group = 0;         // plasmid number, 1-based; 0 for chromosome/unassigned
    bool circular = false;      // its two ends are joined by read pairs
    double depthRatio = 0;      // coverage / modal coverage
    const char* basis = "none"; // which signal decided it, for the report
};

struct RepliconAssignment {
    std::vector<RepliconCall> calls;   // parallel to contigs
    size_t chromosomeContigs = 0;
    size_t plasmidContigs = 0;
    size_t unassignedContigs = 0;
    size_t plasmidGroups = 0;
    size_t circularContigs = 0;
    double modalDepth = 0;
};

// Classifies `contigs`. `layoutMembers` marks contigs the layout stage placed on the
// chromosome track (empty if layout did not run). `links` may be null.
RepliconAssignment assignReplicons(const std::vector<std::string>& contigs,
                                   const std::vector<double>& covs,
                                   const std::vector<char>& layoutMembers,
                                   const OrganismModel& model,
                                   const ContigEndLinks* links,
                                   bool verbose);

const char* repliconName(RepliconClass c);

}  // namespace ts
