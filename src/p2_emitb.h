// Phase 2 (EVAL_PLAN_P2, wave W1, workstream emit, part B): report and label features, every one
// default OFF. With every flag unset no stage below runs, no file is written, contigs.fasta,
// scaffolds.fasta, scaffolds.agp and assembly_graph.gfa are byte-identical to TesserACT 1.5.0, and
// report.json gains only the registered block "p2": {"enabled": 0, "counters": {... all 0}}.
//
//   TESSERACT_P2_ENDS=1            ends.tsv: one row per end of every contigs.fasta record >= 500 bp
//                                  (copy number, non-unique tail, pair partners, unsupported tip bases);
//                                  report.json p2.ends. Footprint: report-only.
//   TESSERACT_P2_TIPS_LOWERCASE=1  the unsupported tip bases of ends.tsv are written in lower case in
//                                  contigs.fasta (never trimmed); implies ends.tsv; p2_edits.tsv lists
//                                  every edit. Footprint: case-only (scaffolds.fasta, AGP, GFA unchanged).
//   TESSERACT_P2_SELF_QA=1         k-mer self-QA (kqa port): report.json p2.self_qa, report.html content
//                                  check, p2_self_qa.{tsv,contigs.tsv,missing_hi.tsv,absent.bed}.
//                                  Footprint: report-only.
//   TESSERACT_P2_PROVENANCE=1      provenance manifest: report.json p2.provenance. Footprint: report-only.
//   TESSERACT_P2_DETECT_REPORT=1   organism detection as a report only: the reads are scored against a
//                                  detection sketch (TESSERACT_P2_DETECT_SKETCH, else
//                                  TESSERACT_OM2_DETECT_SKETCH, else $TESSERACT_MODEL_DIR/om2detect.sketch)
//                                  and the scores go to report.json p2.organism_detect. No model is ever
//                                  selected or applied. Footprint: report-only.
//
// Every run prints the counter lines [p2-ends] [p2-sqa] [p2-prov] [p2-detect] [p2-readpass], with
// enabled=0 and zeros when the flags are off.
#pragma once

#include <cstdio>
#include <string>
#include <vector>

#include "organism_detect.h"
#include "p2_ends.h"
#include "p2_provenance.h"
#include "p2_readpass.h"
#include "p2_selfqa.h"
#include "resolve.h"
#include "seqio.h"

namespace ts {
namespace p2 {

struct EmitBOptions {
    bool ends = false, lowercase = false, selfQa = false, provenance = false, detect = false;
    std::string detectSketch;   // resolved sketch path ("" = none found)
    bool any() const { return ends || lowercase || selfQa || provenance || detect; }
    bool endsAudit() const { return ends || lowercase; }
    bool needsReadPass() const { return endsAudit() || selfQa || detect; }
    static EmitBOptions fromEnv();
};

struct DetectReport {
    bool ran = false;
    std::string error, sketchPath, sketchMd5;
    std::vector<om2::DetectModelScore> scores;   // best first (stable)
    std::vector<std::string> tsmMd5, holdList;   // per entry of `scores`
    double best() const { return scores.empty() ? -1 : scores[0].score; }
    double second() const { return scores.size() < 2 ? -1 : scores[1].score; }
    bool accepted() const { return best() >= om2::kDetectAccept; }
    bool mixture() const { return second() >= om2::kDetectMixture; }
};

struct EmitBState {
    EmitBOptions opt;
    bool readPassRan = false;
    ReadPassStats rp;
    bool endsRan = false;
    EndsStats ends;
    std::string endsError;
    SelfQaResult sqa;
    DetectReport det;
    std::string provenanceJson;
    size_t provFiles = 0;
    uint64_t provBytes = 0;
    std::vector<std::string> edits;   // p2_edits.tsv rows
    std::string fileError;            // a side file that could not be written
};

struct EmitBRecords {
    std::string outDir;
    const std::vector<std::string>* seqs = nullptr;    // contigs.fasta records, final
    const std::vector<std::string>* names = nullptr;
    std::vector<double> covs;
    const SequenceStore* reads = nullptr;
    InsertModel insert;
    int threads = 1;
    std::vector<Library> libraries;                     // the run's libraries, as loaded into `reads`
};

// The input read files in library order, and for each one its (base, stride) in the read store's
// index layout when it can be matched record for record (a run with a single library), else (0, 0).
void readFileLayout(const std::vector<Library>& libs, std::vector<std::string>& files,
                    std::vector<std::pair<uint64_t, uint64_t>>& layout);

// Runs the record-level features just before contigs.fasta is written: the shared raw read pass,
// ends.tsv (+ lower-case tips), self-QA and the detection report; writes their side files. Returns
// true when `written` holds the records contigs.fasta must receive instead of *in.seqs (lower-case
// tips); false leaves contigs.fasta to be written from *in.seqs. Never fails the run: a feature
// that cannot run records why in `st`.
bool runEmitBRecords(const EmitBRecords& in, EmitBState& st, std::vector<std::string>& written);

// The report.json "p2" object (always written).
std::string p2BlockJson(const EmitBState& st);
// The counter lines (always printed).
void logEmitBCounters(const EmitBState& st, std::FILE* log);

}  // namespace p2
}  // namespace ts
