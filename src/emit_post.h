// Emission-stage post-processing, factored out of Assembler::run so it can be tested:
// exact-containment dedup, the split of scaffolds into contigs at every N, the terminal
// dovetail trim, and the auxiliary AGP / GFA path records. With every TESSERACT_FIX_* switch
// off (TESSERACT_FIXES=0; since 1.4.0 the umbrella is on when unset) each function reproduces
// the release 1.3.0 behaviour exactly.
//
// G-emit fixes in here (combo3, 2026-09-25):
//   T37  dedup prefilter (output-neutral, unconditional)
//   T10  TESSERACT_FIX_BOUNDARY_SAFE_TRIM  -- one victim per physical overlap; keep k-1
//   T22  TESSERACT_FIX_TRIM_COPY_GUARD     -- copy-number-aware victim choice with evidence
//   T25  TESSERACT_FIX_SPLIT_POSTPROCESS   -- per-piece coverage/tag, 2k floor, dedup
//   T26  TESSERACT_FIX_AGP_GFA_V2          -- AGP/P-lines describe the written contigs
#pragma once

#include <cstddef>
#include <cstdint>
#include <functional>
#include <string>
#include <vector>

#include "gfa.h"
#include "graph.h"
#include "report.h"
#include "resolve.h"
#include "seqio.h"

namespace ts {

// ---- terminal dovetails (moved unchanged from assembler.cpp) -------------------------------
struct Dovetail {
    size_t a = 0, b = 0;        // a's suffix meets b's prefix, in orientation `bRc`
    size_t len = 0;
    bool bRc = false;
    bool aRc = false;           // the match used rc(a)'s suffix, i.e. a's 5' end
};

std::vector<Dovetail> findTerminalDovetails(const std::vector<std::string>& seqs,
                                            size_t minOverlap, size_t window);

// ---- T37: exact containment with an exact negative prefilter -------------------------------
struct DedupStats {
    size_t queried = 0;          // sequences tested against the kept text
    size_t exactSearches = 0;    // std::string::find calls actually made
    size_t dropped = 0;
    size_t droppedBases = 0;
};

// Longest first (ties by sequence), each sequence that occurs exactly, on either strand,
// inside the concatenation of the sequences kept before it is marked contained. Protected
// sequences are never dropped; `candidate`, when given, limits which may be dropped (the rest
// are kept and join the text). Identical decisions to the release loop: a window filter over
// every 32-byte window of the kept text only skips a find that cannot succeed.
std::vector<char> markContained(const std::vector<std::string>& seqs,
                                const std::vector<char>& protectedSeq,
                                const std::vector<char>* candidate, DedupStats* stats);

// ---- the split at every N ------------------------------------------------------------------
struct EmitPiece {
    size_t scaffold = 0;          // index of the scaffold record it comes from
    size_t ordinal = 0;           // 0-based piece number inside that scaffold
    size_t piecesInScaffold = 1;
    size_t objBegin = 0, objEnd = 0;   // [begin, end) inside the scaffold, before trimming
    std::string seq;              // the record as written (after trimming)
    std::string tag;
    double cov = 0;
    std::string cutFrontSeq, cutBackSeq;   // bases the terminal trim removed, in piece orientation
    uint8_t drop = 0;             // 0 written; 1 below the length floor, 2 exact duplicate (T25)
    std::string name;             // final contigs.fasta name, once assigned
};

// Release split: maximal non-N runs of every scaffold ('N' and 'n' both split).
std::vector<EmitPiece> splitAtGaps(const std::vector<std::string>& scaffolds,
                                   const std::vector<std::string>& tags,
                                   const std::vector<double>& covs);

// ---- scaffold paths cut at their gap steps -------------------------------------------------
struct PathSegment {
    std::vector<uint64_t> oriented;
    double covWeighted = 0;       // sum of unitig coverage x unitig length (the resolver's weighting)
    double covLength = 0;
};
std::vector<PathSegment> pathSegments(const GfaPath& path, const UnitigGraph& g);
// Spelling of a walk: first unitig whole, each next one after its (k-1) overlap.
std::string spellWalk(const std::vector<uint64_t>& oriented, const UnitigGraph& g);
// Which piece of `scaffold` each segment ended up in, found by locating each segment's
// spelled suffix (the part rendered whatever the gap-flank rule) in order. Empty on failure.
std::vector<size_t> mapSegmentsToPieces(const std::string& scaffold,
                                        const std::vector<PathSegment>& segs,
                                        const UnitigGraph& g);

// ---- T25 ------------------------------------------------------------------------------------
struct SplitPostStats {
    bool enabled = false;
    size_t pieces = 0, multiPieceScaffolds = 0;
    size_t droppedShort = 0, droppedShortBp = 0;
    size_t droppedDup = 0, droppedDupBp = 0;
    size_t circularTagsDropped = 0;
    size_t covRelabelled = 0, covFallback = 0;
};
// `paths[i]` is scaffold i's walk (empty when it has none); `graph` may be null.
void splitPostprocess(std::vector<EmitPiece>& pieces, const std::vector<std::string>& scaffolds,
                      const std::vector<GfaPath>& paths, const UnitigGraph* graph,
                      const std::vector<char>& protectedScaffold, size_t minLen,
                      SplitPostStats& st);
void logSplitPost(const SplitPostStats& st);

// ---- T10 / T22: the terminal dovetail trim -------------------------------------------------
// Evidence about one contig end holding a terminal segment (T22 stage 2).
struct EndEvidence {
    size_t consistent = 0;     // pairs anchored in the unique flank whose mate lies in the segment where predicted
    size_t inconsistent = 0;   // pairs predicted wholly inside the contig whose mate is not found there
    double copyRatio = 0;      // read k-mer depth of the segment / depth of the unique flank (0 = unknown)
    bool supported() const { return consistent >= 3 && inconsistent * 5 <= consistent + inconsistent; }
    bool contradicted() const { return inconsistent >= 3 && inconsistent * 2 > consistent + inconsistent; }
};
struct EndQuery {
    size_t seq = 0;
    bool front = false;        // the segment is at the sequence's 5' end
    size_t segLen = 0;
};
using EndEvidenceFn = std::function<std::vector<EndEvidence>(const std::vector<EndQuery>&)>;

struct TrimConfig {
    size_t minOverlap = 0;       // release: max(--trim-overlap, finalK)
    int k = 0;                   // final k
    bool boundarySafe = false;   // T10
    bool copyGuard = false;      // T22
};
struct TrimStats {
    size_t trimmed = 0, trimmedBases = 0;
    // T10
    size_t boundaryBpKept = 0, equalLenFlips = 0, skippedByBoundary = 0;
    // T22
    size_t components = 0, nonBipartite = 0, guardDiffers = 0, guardApplied = 0;
    size_t evidenceSwitched = 0, evidenceRefused = 0, endsRescued = 0, releaseCutsAvoided = 0;
};
// Trims `seqs` in place. cutFront/cutBack receive the bases removed at each end. With both
// switches off this is the release loop byte for byte. `evidence` is consulted only by the
// copy guard; without it (or with an unusable insert model) the guard keeps release decisions.
TrimStats trimTerminalOverlaps(std::vector<std::string>& seqs, const TrimConfig& cfg,
                               std::vector<size_t>& cutFront, std::vector<size_t>& cutBack,
                               const EndEvidenceFn* evidence);
void logTrimCounters(const TrimConfig& cfg, const TrimStats& st);

// Read-pair evidence for the copy guard (T22 stage 2): one pass over the reads.
std::vector<EndEvidence> gatherEndEvidence(const std::vector<std::string>& seqs,
                                           const std::vector<EndQuery>& queries,
                                           const SequenceStore& reads, const InsertModel& insert,
                                           int threads);

// ---- T26: AGP and GFA paths that describe contigs.fasta ------------------------------------
struct AgpGfaStats {
    bool enabled = false;
    size_t w = 0, wPartner = 0, wUnlisted = 0, nPaired = 0, nUnspecified = 0;
    size_t pEmitted = 0, pDropped = 0;
};
// Writes scaffolds.agp: every scaffold tiled by W rows naming contigs.fasta records (and, for
// bases the trim or the dedup removed from a piece, the record that still holds them) and by N
// rows whose evidence is `gapEvidence`. Returns false on an I/O error.
bool writeAgpV2(const std::string& path, const std::vector<std::string>& scaffolds,
                const std::vector<std::string>& scaffoldNames, const std::vector<EmitPiece>& pieces,
                const std::string& gapEvidence, AgpGfaStats& st, std::string& error);
// One P-line per written contig whose own walk spells it exactly; the rest are dropped.
std::vector<GfaPath> contigPathsV2(const std::vector<std::string>& scaffolds,
                                   const std::vector<GfaPath>& paths, const UnitigGraph& g,
                                   const std::vector<EmitPiece>& pieces, AgpGfaStats& st);
void logAgpGfa(const AgpGfaStats& st);

// ---- T17 -------------------------------------------------------------------------------------
// Fills rep.contigPieces/contigTotal/contigLargest/contigN50 from the records contigs.fasta
// holds, and rep.scaffoldGaps with the number of N-runs of any length in the scaffolds.
void computeContigStats(const std::vector<std::string>& records,
                        const std::vector<std::string>& scaffolds, AssemblyReport& rep);

}  // namespace ts
