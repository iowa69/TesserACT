// Phase 2, emit-B: the contig-end audit (ends.tsv) and lower-case unsupported tips.
//
// One row per end of every contigs.fasta record of at least 500 bp (EVAL_PLAN_P2 L-END-a): left end
// (L, record start) and right end (R, record end); a circular record (name tag `_circular`) gets the
// same two rows with topology=circular. Columns, with the second-pass prototype they port
// (backward/secondpass/proto/sp2v5.py; TRACK_RESULT s2):
//   * copy number: the record's k-mer coverage (the `cov` of its name) over the genome coverage, the
//     length-weighted median coverage of the records >= 5 kb (sp2: depth / genome depth); the end's raw
//     31-mer read depth (median over its last 300 bp) and that over the genome's (the median end depth of
//     the records >= 5 kb);
//   * non-unique tail: terminal bases covered by a 31-mer that occurs more than once in contigs.fasta
//     (both strands), within the last 741 bp (sp2 `nonunique_tail`, span = 700 + 31 + 10);
//   * pair-partner structure from read pairs anchored on the records' end windows (the corrected reads the
//     assembler holds, anchored like ContigEndLinks: 12 probes, >= 2 agreeing unique 31-mer votes):
//     outward reads within the insert reach, their mates, links to other ends, and the class
//     unique | branching | repeat_only | none | na (sp2: >= 3 links; a partner counts when its record is a
//     single-copy anchor, cn <= 1.6 at >= 250 bp or cn <= 4 at >= 5 kb; branching when the second
//     anchor partner has >= 1/3 of the best one's links);
//   * unsupported tip: terminal bases whose 31-mer the end's OWN reads do not vouch for (sp2's end audit
//     `tip_unvouched`): the locus pool is every read anchored (unique 31-mer votes, as above) on the last
//     400 bp of the end plus the mate of every read pointing out of the end within the insert reach; its
//     31-mers are counted on both strands from the RAW reads (3' trimmed below Q10, sp2 `qtrim`), fetched
//     by their index in the assembler's read store during the raw read pass; from the end inward, every
//     position whose k-mer has a pool count below max(2, 0.25 x the pool's median count over the end's
//     last 300 bp) is unsupported, stopping at the first vouched k-mer, at most 200 positions + 1. As in
//     sp2, only single-copy anchor records with a pool of >= 10 reads are audited (tip_status says why
//     another end is not). No mapper is run. When the raw records cannot be matched to the store (more
//     than one library, or a record count that differs), the pool's corrected reads from the store are
//     counted instead and pool_source says so.
// With TESSERACT_P2_TIPS_LOWERCASE=1 the unsupported tip bases are written in lower case in contigs.fasta
// (never trimmed; a run of exactly `unsupported_tip_bp` bases at that end), and every such edit is listed
// in p2_edits.tsv.
#pragma once

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "p2_kmer.h"
#include "resolve.h"
#include "seqio.h"

namespace ts {
namespace p2 {

constexpr size_t kEndsMinRecord = 500;
constexpr size_t kEndsSpan = 741;
constexpr size_t kEndsDepthWindow = 300;
constexpr size_t kEndsTipScan = 200;
constexpr size_t kEndsMinLinks = 3;
constexpr size_t kEndsPoolFlank = 400;
constexpr size_t kEndsMinPool = 10;

struct EndRow {
    size_t record = 0;
    char end = 'L';                 // 'L' record start, 'R' record end
    bool circular = false;
    double cn = 0;                  // record coverage / genome coverage
    bool anchor = false;
    double endKdepth = 0;           // raw 31-mer median count, last 300 bp
    double endCn = 0;               // endKdepth / genome end depth
    size_t nonuniqueTail = 0;
    size_t outward = 0, matesUnplaced = 0, matesOtherEnd = 0, selfWrap = 0;
    std::string partnerClass = "na";
    size_t nPartners = 0;           // partners with >= 3 links
    std::string bestPartner;        // name:end:links
    std::string partners;           // top 5, name:end:links;...
    size_t pool = 0;                // locus pool reads
    double poolKdepth = 0;          // pool 31-mer median count over the last 300 bp
    std::string tipStatus = "not_run";   // audited | not_anchor | small_pool | not_run
    size_t tip = 0;                 // unsupported tip bases
    size_t lowercased = 0;          // bases written in lower case at this end
};

struct EndsStats {
    size_t records = 0, ends = 0, tipEnds = 0, tipBases = 0, lowercasedBases = 0;
    size_t unique = 0, branching = 0, repeatOnly = 0, none = 0, na = 0;
    size_t readsPlaced = 0, pairsUsed = 0;
    size_t audited = 0, poolReads = 0;
    std::string poolSource = "none";
    double genomeCov = 0, genomeEndKdepth = 0;
    double seconds = 0;
};

struct EndsInput {
    const std::vector<std::string>* seqs = nullptr;    // contigs.fasta records as they will be written
    const std::vector<std::string>* names = nullptr;
    const std::vector<double>* covs = nullptr;         // per record k-mer coverage
    const SequenceStore* reads = nullptr;              // the assembler's reads (pair placement)
    InsertModel insert;
    int threads = 1;
};

class EndsAudit {
public:
    explicit EndsAudit(const EndsInput& in);
    // The k-mers whose raw read counts the tips need (every 31-mer of the last 300 bp of each end).
    void addRawKeys(FixedCountTable& t) const;
    // Copy number, non-unique tails, pair partners and the locus pools (no raw counts needed).
    void computeStructure();
    // The locus-pool read indices (into in.reads) as a bitset, for the raw read pass.
    const std::vector<uint64_t>& poolBits() const { return poolBits_; }
    // Counts one pool read's 31-mers into the tables of every end whose pool holds it. Thread-safe.
    void countPoolRead(uint64_t readIndex, const char* seq, size_t len);
    // Discards the raw pool counts and counts the pools' corrected reads from the store instead.
    void countPoolsFromStore();
    // Tips (and end depths) from the raw read counts; then rows() is complete.
    void finishTips(const FixedCountTable& raw, const std::string& poolSource);
    const std::vector<EndRow>& rows() const { return rows_; }
    EndsStats& stats() { return stats_; }
    const EndsStats& stats() const { return stats_; }

private:
    const EndsInput in_;
    std::vector<EndRow> rows_;
    std::vector<size_t> rowRecords_;    // records with rows (>= 500 bp)
    EndsStats stats_;
    std::vector<std::vector<uint64_t>> rowPool_;          // per row: pool read indices, sorted
    std::vector<std::pair<uint64_t, uint32_t>> readRows_; // (read, row) sorted by read
    std::vector<uint64_t> poolBits_;
    std::vector<FixedCountTable> local_;                  // per row: the end's last-300-bp 31-mers
    std::vector<char> audit_;                             // per row: audited
};

// contigs.fasta records with the tips of `rows` in lower case (records without rows unchanged), and the
// edit log rows "record\tfeature\toperation\tstart\tend\tbases" (1-based, inclusive).
std::vector<std::string> lowercaseTips(const std::vector<std::string>& seqs, const std::vector<std::string>& names,
                                       std::vector<EndRow>& rows, EndsStats& st, std::vector<std::string>& edits);

std::string endsTsv(const std::vector<EndRow>& rows, const std::vector<std::string>& names,
                    const std::vector<std::string>& seqs);
std::string endsJson(const EndsStats& st, bool lowercase);

// The prototype's tip rule on a list of k-mer counts ordered from the end inward (counts[j] = the count
// of the k-mer ending j bases before the end) with the end's median count `lda`: positions scanned
// before the first count >= max(2, 0.25 lda), at most counts.size().
size_t tipFromCounts(const std::vector<uint32_t>& counts, double lda);
// Python statistics.median (mean of the two middle values for an even count); 0 for none.
double medianOf(std::vector<uint32_t> v);

}  // namespace p2
}  // namespace ts
