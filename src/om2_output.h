// Organism Model 2.0, component C3 (SURFACE): the owner's genome view and its labels.
//
// Default OFF. With TESSERACT_OM2_OUTPUT unset nothing here writes a file, nothing reads the
// assembly beyond what is needed to print the counter line, and report.json carries no om2
// block, so every release output is byte-identical (EVAL_PLAN E0d). One [om2-out] counter
// line is printed on every run, with enabled=0 and zeros when off.
//
// With TESSERACT_OM2_OUTPUT=1 the finished records are written again, unchanged in sequence,
// as the output contract of DESIGN s5.3 / eval/invariants.py:
//
//   genome/genome.fasta      chromosome first, then plasmid groups, then unplaced records;
//                            circular and dnaA-first only when a closed wrap junction says so
//   genome/genome.agp        AGP 2.1; W components are scaffolds.fasta records (or contigs.fasta
//                            records in a gap-free run) or om2fill_<junction>; every gap line
//                            carries its linkage evidence
//   genome/fills.fasta       one record per om2fill_ component
//   genome/genome.mask.bed   basis label of every filled base; tiles every fill exactly
//   genome/junctions.tsv     one row per seam or ledger junction, with every ledger field
//   genome/repeat_variants.tsv  inter-copy variant sites of multi-copy repeats (from C2)
//   genome/closure.txt       per record: closed or open, and what would close each open gap
//   genome/README.txt        the three views and what each label means
//
// The writer never invents a join, a base or an order. Joins and fills come only from the
// ledger (C1 junctions, C2 fills); without a ledger every N-run is an `unrecorded` seam kept
// exactly as the release wrote it. Allocated bases (PRIOR_ALLOCATED, MULTIPLICITY, CONSENSUS)
// are written in lowercase and only into genome/, never into contigs.fasta or scaffolds.fasta.
#pragma once

#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <functional>
#include <string>
#include <utility>
#include <vector>

#include "om2_types.h"

namespace ts {
namespace om2 {

class ClonalSurface;   // om2_clonal.h (C4)

// RepeatVariant (one row of genome/repeat_variants.tsv) is declared in om2_types.h since the
// integration (build_om2): C2 produces the rows, the writer renders them.

struct SurfaceInput {
    std::string outDir;
    // The records the view is built from, in file order: the scaffolds.fasta records when the
    // run wrote one (scaffoldsFile), else the contigs.fasta records. Names are the headers.
    const std::vector<std::string>* seqs = nullptr;
    const std::vector<std::string>* names = nullptr;
    bool scaffoldsFile = false;
    const Ledger* ledger = nullptr;                       // C1/C2; null: every seam is unrecorded
    const std::vector<RepeatVariant>* variants = nullptr; // C2; null: header-only table
    bool modelRan = false;                                // a model joined or laid out contigs
    std::string modelPath, modelName;
    std::vector<std::pair<std::string, std::string>> sidecars;   // label -> path, hashed into report.json
    // Integration (build_om2), both optional:
    //   runOwner       per record of `seqs`, per N-run in order: the ledger row C1 located there
    //                  (SeamContext::runOwners). Used before the exact-flank match, which remains the
    //                  rule for runs it leaves unowned (C2's own rows, runs C1 could not place).
    //   originLocator  where the replication origin (dnaA) starts in a closed circular record, given
    //                  as the genome record's text: >= 0 forward strand, <= -2 reverse strand at
    //                  (-2 - value), -1 not found. When set it replaces Ledger::rotateOffset, whose
    //                  stage-[4c/7] coordinates no longer hold once gaps are closed or filled.
    const std::vector<std::vector<RunOwner>>* runOwner = nullptr;
    std::function<int64_t(const std::string&)> originLocator;
    // The other stages' counter lines ([om2-evidence], [om2-seam], [om2-close]) as printed, copied into
    // the report.json om2 block under "stage_counters" (label -> line).
    std::vector<std::pair<std::string, std::string>> stageCounters;
    // C4 (om2_clonal): when set, every junction of the genome view is labelled (class, confidence, claim,
    // nearest-relative, graph and pair support: extra junctions.tsv columns and genome.agp comment lines),
    // closure.txt starts with the nearest relatives, and report.json gets a clonal block. Null: unchanged.
    ClonalSurface* clonal = nullptr;
};

struct SurfaceStats {
    bool enabled = false;
    bool agpEvidence = false;
    bool ok = true;
    std::string error;
    size_t records = 0, chrRecords = 0, circular = 0, rotated = 0, rotateRefused = 0;
    size_t genomeOnlyJoins = 0, broken = 0, fills = 0, allocatedBp = 0, agpEvidenceGaps = 0;
    size_t seams = 0, unrecorded = 0, ambiguous = 0, scaffoldSplits = 0, variantRows = 0;
    size_t ownedByC1 = 0, overlapTrimmed = 0, reversed = 0;   // integration counters
    double detectScore = -1, detectSecond = -1;           // -1: detection did not run
};

bool outputEnabled();        // TESSERACT_OM2_OUTPUT=1
bool agpEvidenceEnabled();   // TESSERACT_OM2_AGP_EVIDENCE=1

// Writes genome/ (when enabled), rewrites scaffolds.fasta/.agp only where the ledger admits a
// junction to the genome view alone, adds per-gap evidence to scaffolds.agp (when enabled),
// and renders the report.json om2 block into `om2Json` (left empty when nothing om2 ran).
SurfaceStats writeSurface(const SurfaceInput& in, std::string& om2Json);
void logSurfaceCounters(const SurfaceStats& s, std::FILE* log = stderr);

// ---- pieces exposed for the component tests -------------------------------------------
std::string revcomp(const std::string& s);
const char* basisName(Basis b);
bool isAllocation(Basis b);   // PRIOR_ALLOCATED, MULTIPLICITY, CONSENSUS
const char* sourceName(Source s);
const char* verdictName(Verdict v);
const char* tierName(Tier t);
const char* admitName(Admit a);
const char* endClassName(EndClass e);
// Replicon tag at the end of a TesserACT record name: "_chr", "_plas_3_circular", "_unk", ...
struct TagInfo { char cls = 'u'; uint32_t group = 0; bool circular = false; };
TagInfo parseTag(const std::string& name);
// AGP 2.1 linkage evidence for a gap: from the junction when there is one, else from whether a
// model ran (the release's own rule: without a model every N comes from the pair scaffolder).
std::string agpEvidence(const Junction* j, bool modelRan);
// Checks an AGP 2.1 file for the structure the contract needs: 9 columns, contiguous object
// coordinates from 1, part numbers from 1, component spans equal object spans, U gaps 100 bp,
// linkage evidence present on every gap. Returns an empty string when valid.
std::string validateAgp(const std::string& path);
std::string sha256File(const std::string& path);
std::string md5File(const std::string& path);

}  // namespace om2
}  // namespace ts
