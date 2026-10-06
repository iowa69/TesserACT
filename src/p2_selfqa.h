// Phase 2, emit-B, feature F2: k-mer self-QA of the assembly against its own raw reads (`self_qa`).
//
// A port of the kqa prototype (backward/novel/tools/kqa.cpp; EVAL_PLAN_P2 s6.1 "reference prototypes"),
// with kqa's parameters: canonical 31-mers hash-sampled 1 in 8, the contigs.fasta records of at least
// 500 bp (kqa -m 500), every read base of every input file, PhiX174 (NC_001422.1) as the spike-in
// reference (kqa -x phix.fa). From one table of sampled k-mers with their read count and assembly
// multiplicity it reports:
//   * completeness: the share of solid read k-mers (count >= the spectrum valley) present in the assembly;
//   * QV from read-absent assembly k-mers (and from k-mers below the valley);
//   * missing sequence by depth class (count / spectrum peak: <0.1, 0.1-0.5, 0.5-1.5, 1.5-3, 3-10, >10x),
//     each as distinct sampled k-mers and as an estimate in bp (k-mers x 8);
//   * spike-in k-mers (PhiX): reference k-mers, solid in the reads, present in the assembly, copy number.
// The lost-replicon alarm is EVAL_PLAN_P2 L-SQA-b's frozen rule (kqa "flag A",
// backward/novel/kqa/lostrep_flagging.txt): missing_bp_est_gt10x >= 1000 bp.
//
// The collapsed-repeat and genome-size estimators of kqa are not ported: the novel track measured them
// as a negative result (genome-size error 2.2 % against 1.2 % for the plain assembly length).
#pragma once

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "p2_kmer.h"

namespace ts {
namespace p2 {

constexpr uint64_t kSqaSample = 8;
constexpr size_t kSqaMinRecord = 500;
constexpr uint64_t kSqaAlarmBp = 1000;      // frozen: EVAL_PLAN_P2 L-SQA-b (kqa flag A)
constexpr int kSqaBins = 6;
extern const char* const kSqaBinName[kSqaBins];   // lt0.1x 0.1-0.5x 0.5-1.5x 1.5-3x 3-10x gt10x

struct SqaContig {
    std::string name;
    size_t len = 0;
    uint64_t sampled = 0, read0 = 0, belowValley = 0;
    double medianCn = 0;
};

struct SelfQaResult {
    bool ran = false;
    std::string error;
    size_t asmRecords = 0;
    uint64_t asmBases = 0;
    uint64_t reads = 0, readBases = 0;
    uint32_t valley = 0, peak = 0;
    uint64_t solid = 0, solidFound = 0;
    double completeness = 0;
    uint64_t asmPos = 0, asmRead0 = 0, asmBelowValley = 0;
    double qvRead0 = 0, qvSolid = 0;
    uint64_t spikeRef = 0, spikeSolid = 0, spikeInAsm = 0;
    double spikeMedianCn = 0;
    uint64_t solidBin[kSqaBins] = {0, 0, 0, 0, 0, 0};
    uint64_t missBin[kSqaBins] = {0, 0, 0, 0, 0, 0};
    uint64_t missingBpEst(int b) const { return missBin[b] * kSqaSample; }
    bool alarm = false;
    std::vector<SqaContig> contigs;
    std::vector<std::pair<uint64_t, uint32_t>> missingHi;    // solid, absent, count/peak >= 3 (sorted)
    std::vector<std::pair<std::string, std::pair<long, long>>> absentWindows;  // read-absent windows
    std::vector<long> absentCounts;
    size_t tableDistinct = 0, tableBytes = 0;
    double seconds = 0;
};

class SelfQa {
public:
    // Records shorter than kSqaMinRecord are skipped (kqa -m 500). `spike` holds the spike-in
    // reference sequences (PhiX174 by default).
    void prepare(const std::vector<std::string>& records, const std::vector<std::string>& names,
                 const std::vector<std::string>& spike);
    ShardedCountTable& table() { return table_; }
    // Everything after the read pass. `reads`/`readBases` come from the pass statistics.
    SelfQaResult finish(uint64_t reads, uint64_t readBases) const;

private:
    ShardedCountTable table_;
    std::vector<const std::string*> recs_;
    std::vector<std::string> names_;
    std::vector<uint64_t> spike_;
    uint64_t asmPos_ = 0, asmBases_ = 0;
};

// The prototype's QV from a count of bad sampled positions out of `positions` (kqa.cpp `qv`).
double sqaQv(uint64_t bad, uint64_t positions);

// kqa's summary lines ("key\tvalue"), same keys as <prefix>.summary.tsv (minus the negative-result
// estimators), plus the alarm.
std::string selfQaSummaryTsv(const SelfQaResult& r);
// The report.json `self_qa` object (no trailing comma or newline).
std::string selfQaJson(const SelfQaResult& r);
// One-sentence content verdict for the HTML report and the counter line.
std::string selfQaVerdict(const SelfQaResult& r);
// Side files: <dir>/p2_self_qa.tsv, p2_self_qa.contigs.tsv, p2_self_qa.missing_hi.tsv,
// p2_self_qa.absent.bed. False with `error` on a write failure.
bool writeSelfQaFiles(const std::string& dir, const SelfQaResult& r, std::string& error);

std::string kmerString(uint64_t canon);

}  // namespace p2
}  // namespace ts
