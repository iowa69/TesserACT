// Laying an assembly out against the panel chromosome it most resembles.
//
// The adjacency join in organism_join.h answers a local question -- do these two flanks
// belong together -- and answers it well. Closing a chromosome means answering about forty
// such questions in a row without a single mistake, because one miss leaves the chromosome
// in two pieces. On this cohort that would require roughly 99.7% success per junction
// against 83.9% achieved, and no threshold reaches that.
//
// This stage asks a global question instead: given a panel chromosome's marker order,
// where does each contig sit on a genome shaped like that? Forty independent decisions
// collapse into one, and a locus the template gets wrong costs that locus rather than the
// whole chromosome.
//
// What that does and does not claim:
//   - the sequence is still assembled de novo; only ORDER and ORIENTATION come from the
//     panel, and they are reported separately for exactly that reason
//   - a gap is written as N and is not sequence
//   - a contig whose markers do not agree on one location is left out of the scaffold
//     rather than forced into it, because forcing it would invent a rearrangement
//
// Measured on 666 closed-reference isolates, laying out against a *randomly chosen* panel
// chromosome recovers a median 96.5% of the chromosome in one piece against 98.7% for the
// nearest relative. The method therefore rests mainly on species-level synteny, which in
// Klebsiella pneumoniae is strong, and only secondarily on picking a close relative. That
// is worth knowing before trusting it on an organism whose chromosome is less conserved.
#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "organism.h"

namespace ts {

struct LayoutStats {
    bool run = false;
    std::string track;          // panel chromosome used, for provenance
    size_t trackMarkers = 0;
    size_t sharedMarkers = 0;   // markers the assembly and the track have in common
    size_t contigsIn = 0;
    size_t contigsOut = 0;
    size_t placed = 0;          // contigs given a position on the track
    size_t unplaced = 0;        // left as they were: too few markers, or no agreement
    size_t incoherent = 0;      // markers agreed on no single location
    size_t dropped = 0;         // placement overlapped one already taken
    size_t scaffolds = 0;       // pieces the placed contigs formed
    size_t gapBases = 0;
    size_t overlapMerges = 0;   // seams closed on an exact sequence overlap
    size_t plasmidSkipped = 0;  // contigs the panel calls plasmid, kept off the track
};

// Orders and orients `contigs` against the best-matching layout track in `model`.
// `covs` is reordered to match. Contigs that cannot be placed are kept, untouched, after
// the scaffolds. Returns what happened; `run` is false when the model carries no tracks
// or too few markers could be located to choose one.
LayoutStats layoutByModel(const OrganismModel& model, std::vector<std::string>& contigs,
                          std::vector<double>& covs, int k, bool verbose);

}  // namespace ts
