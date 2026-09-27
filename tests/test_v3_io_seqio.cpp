// Component test for the read loader (combo3 G-io: T18, T33, and the library-level part of T32).
//
// T18  a gzip input that ends inside a member (no final block, or no CRC/ISIZE trailer) must
//      fail to load. Release 1.3.0 took zlib's "unexpected end of file" (Z_BUF_ERROR, returned
//      as a 0-byte read) for a clean end and loaded the records it happened to hold.
// T33  a blank line where a FASTQ header is due (trailing, or between records) is skipped, and
//      mates named SRR.spot.1 / SRR.spot.2 (fastq-dump --readids) pair up. Release 1.3.0 failed
//      both with exit 1.
// T32  (library level) the same file given twice -- as both mates, via a symlink, or as two
//      single-end libraries -- and an interleaved library that also names an r2 file are refused.
//
// Every tolerated input must load exactly the reads of its plain counterpart, and every input
// the release accepted must still load the same reads; both are checked by read digests.
// All fixtures are generated in-process (zlib for the gzip ones); nothing is read from disk
// outside a private temporary directory.
#include "seqio.h"

#include <zlib.h>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <random>
#include <string>
#include <vector>

#include <unistd.h>
#include "test_env.h"   // build_v3: T40, clear the ambient TESSERACT_* first

namespace {

int checks = 0;
int failures = 0;

void check(bool ok, const std::string& what) {
    ++checks;
    if (!ok) {
        ++failures;
        std::fprintf(stderr, "FAIL: %s\n", what.c_str());
    }
}

std::string dir;
std::string at(const std::string& name) { return dir + "/" + name; }

constexpr size_t kReads = 2000;
constexpr size_t kLen = 100;
std::vector<std::string> seqs;

void writeText(const std::string& path, const std::string& text) {
    std::ofstream f(path, std::ios::binary);
    f << text;
}

// Writes `text` gzipped. When `flushAt` is inside the text, the stream is sync-flushed there and
// the compressed offset of that point is returned: cutting the file at it leaves a stream that
// decodes to exactly text[0, flushAt) and then stops inside the member.
size_t writeGz(const std::string& path, const std::string& text, size_t flushAt = std::string::npos) {
    gzFile g = gzopen(path.c_str(), "wb");
    if (!g) { std::fprintf(stderr, "cannot write %s\n", path.c_str()); std::exit(2); }
    size_t cut = 0;
    if (flushAt < text.size()) {
        gzwrite(g, text.data(), static_cast<unsigned>(flushAt));
        gzflush(g, Z_SYNC_FLUSH);
        cut = static_cast<size_t>(gzoffset(g));
        gzwrite(g, text.data() + flushAt, static_cast<unsigned>(text.size() - flushAt));
    } else if (!text.empty()) {
        gzwrite(g, text.data(), static_cast<unsigned>(text.size()));
    }
    gzclose(g);
    return cut;
}

void truncateTo(const std::string& path, size_t size) { std::filesystem::resize_file(path, size); }
size_t sizeOf(const std::string& path) { return static_cast<size_t>(std::filesystem::file_size(path)); }

std::string fastqRecord(const std::string& name, const std::string& seq, const char* eol = "\n") {
    return "@" + name + eol + seq + eol + "+" + eol + std::string(seq.size(), 'I') + eol;
}

// name(i) gives the header text of read i.
template <typename Name>
std::string fastq(size_t from, size_t to, Name name, const char* eol = "\n") {
    std::string s;
    for (size_t i = from; i < to; ++i) s += fastqRecord(name(i), seqs[i], eol);
    return s;
}

std::string fasta(size_t from, size_t to) {
    std::string s;
    for (size_t i = from; i < to; ++i) s += ">f" + std::to_string(i) + "\n" + seqs[i] + "\n";
    return s;
}

// Interleaved: pair i is R1 = seqs[i], R2 = seqs[kReads-1-i].
std::string interleaved(size_t pairs, const char* s1, const char* s2) {
    std::string s;
    for (size_t i = 0; i < pairs; ++i) {
        s += fastqRecord("SRR1." + std::to_string(i + 1) + s1, seqs[i]);
        s += fastqRecord("SRR1." + std::to_string(i + 1) + s2, seqs[kReads - 1 - i]);
    }
    return s;
}

// The loader's counters. Built with -DTS_V3_IO_RELEASE_BASELINE, the test compiles against the
// release 1.3.0 header, which has no counters: they then read as zero and their checks are
// skipped, so the behavioural checks can be shown failing on the release loader.
struct Stats {
    size_t files = 0, records = 0, interRecordBlankLines = 0, dotSuffixMatePairs = 0;
};
#ifdef TS_V3_IO_RELEASE_BASELINE
constexpr bool kCounters = false;
Stats statsOf(const ts::SequenceStore&) { return Stats(); }
#else
constexpr bool kCounters = true;
Stats statsOf(const ts::SequenceStore& st) {
    const ts::LoadStats& l = st.loadStats();
    return Stats{l.files, l.records, l.interRecordBlankLines, l.dotSuffixMatePairs};
}
#endif

void counterCheck(bool ok, const std::string& what) {
    if (kCounters) check(ok, what);
}

struct Result {
    bool ok = false;
    std::string error;
    size_t reads = 0;
    uint64_t digest = 0;
    Stats stats;
};

uint64_t digestOf(const ts::SequenceStore& st) {
    uint64_t h = 1469598103934665603ULL;
    std::string s;
    for (size_t i = 0; i < st.size(); ++i) {
        st.decode(i, s);
        for (char c : s) { h ^= static_cast<unsigned char>(c); h *= 1099511628211ULL; }
        h ^= 0xff; h *= 1099511628211ULL;
    }
    return h;
}

Result load(const std::vector<ts::Library>& libs) {
    ts::SequenceStore st;
    Result r;
    r.ok = st.load(libs, 2, r.error);
    r.reads = st.size();
    r.stats = statsOf(st);
    if (r.ok) r.digest = digestOf(st);
    return r;
}

ts::Library single(const std::string& f) { ts::Library l; l.r1 = f; return l; }
ts::Library pair(const std::string& a, const std::string& b) { ts::Library l; l.r1 = a; l.r2 = b; return l; }
ts::Library inter(const std::string& f) { ts::Library l; l.r1 = f; l.interleaved = true; return l; }

bool has(const std::string& text, const char* needle) { return text.find(needle) != std::string::npos; }

void expectOk(const std::string& label, const Result& r, size_t reads) {
    check(r.ok, label + ": loads (error: " + r.error + ")");
    check(r.reads == reads, label + ": " + std::to_string(reads) + " reads (got " + std::to_string(r.reads) + ")");
}

void expectFail(const std::string& label, const Result& r, const char* needle) {
    check(!r.ok, label + ": refused (loaded " + std::to_string(r.reads) + " reads)");
    check(needle == nullptr || has(r.error, needle),
          label + ": error mentions '" + (needle ? needle : "") + "' (got: " + r.error + ")");
    counterCheck(r.stats.records == 0 && r.stats.files == 0, label + ": counters zero after a failed load");
}

auto slashName(int mate) {
    return [mate](size_t i) { return "r" + std::to_string(i) + "/" + std::to_string(mate); };
}

}  // namespace

int main() {
    testenv::clearTesseractEnv();
    char tmpl[] = "/tmp/tesseract_v3_io_XXXXXX";
    const char* tmpdir = std::getenv("TMPDIR");
    std::string pattern = std::string(tmpdir && *tmpdir ? tmpdir : "/tmp") + "/tesseract_v3_io_XXXXXX";
    std::vector<char> buf(pattern.begin(), pattern.end());
    buf.push_back('\0');
    char* made = mkdtemp(buf.data());
    if (!made) made = mkdtemp(tmpl);
    if (!made) { std::perror("mkdtemp"); return 2; }
    dir = made;

    std::mt19937 rng(20260925);
    seqs.resize(kReads);
    for (auto& s : seqs) {
        s.resize(kLen);
        for (char& c : s) c = "ACGT"[rng() % 4];
    }
    auto plainName = [](size_t i) { return "q" + std::to_string(i); };

    // ---------------------------------------------------------------- reference loads
    const std::string se = fastq(0, kReads, plainName);
    writeText(at("se.fq"), se);
    const Result ref = load({single(at("se.fq"))});
    expectOk("plain single-end FASTQ", ref, kReads);
    counterCheck(ref.stats.files == 1 && ref.stats.records == kReads && ref.stats.interRecordBlankLines == 0 &&
              ref.stats.dotSuffixMatePairs == 0,
          "plain single-end FASTQ: counters files=1 records=2000 blank=0 dot=0");

    writeText(at("p_1.fq"), fastq(0, kReads, slashName(1)));
    writeText(at("p_2.fq"), fastq(0, kReads, slashName(2)));
    const Result refPe = load({pair(at("p_1.fq"), at("p_2.fq"))});
    expectOk("paired /1 /2 FASTQ", refPe, 2 * kReads);
    counterCheck(refPe.stats.files == 2 && refPe.stats.dotSuffixMatePairs == 0, "paired /1 /2: files=2 dot=0");

    // ---------------------------------------------------------------- T18: truncated gzip
    {
        writeGz(at("se.fq.gz"), se);
        const Result r = load({single(at("se.fq.gz"))});
        expectOk("T18 control: intact gzip", r, kReads);
        check(r.digest == ref.digest, "T18 control: intact gzip reads == plain reads");
    }
    {
        // Cut exactly at a record boundary: every record decoded is complete, so only the
        // gzip status can tell that the file is incomplete.
        const std::string head = fastq(0, 1200, plainName);
        const size_t cut = writeGz(at("cut_se.fq.gz"), se, head.size());
        truncateTo(at("cut_se.fq.gz"), cut);
        const Result r = load({single(at("cut_se.fq.gz"))});
        expectFail("T18 gzip cut at a record boundary (single-end)", r, "unexpected end of file");
        check(has(r.error, "truncated or incomplete"), "T18 record-boundary cut: message names truncation");
    }
    {
        writeGz(at("notrailer.fq.gz"), se);
        truncateTo(at("notrailer.fq.gz"), sizeOf(at("notrailer.fq.gz")) - 8);
        expectFail("T18 gzip without its CRC/ISIZE trailer", load({single(at("notrailer.fq.gz"))}),
                   "unexpected end of file");
    }
    {
        writeGz(at("fa.fa.gz"), fasta(0, kReads));
        writeText(at("fa.fa"), fasta(0, kReads));
        const Result whole = load({single(at("fa.fa.gz"))});
        expectOk("T18 control: intact FASTA gzip", whole, kReads);
        check(whole.digest == load({single(at("fa.fa"))}).digest, "T18 control: FASTA gzip reads == plain FASTA");
        truncateTo(at("fa.fa.gz"), sizeOf(at("fa.fa.gz")) / 2);
        // FASTA had it worst: almost any cut silently shortened the last record.
        expectFail("T18 FASTA gzip cut mid-stream", load({single(at("fa.fa.gz"))}), "unexpected end of file");
    }
    {
        const std::string r1 = fastq(0, kReads, slashName(1)), r2 = fastq(0, kReads, slashName(2));
        const size_t c1 = writeGz(at("cut_1.fq.gz"), r1, fastq(0, 1500, slashName(1)).size());
        const size_t c2 = writeGz(at("cut_2.fq.gz"), r2, fastq(0, 1500, slashName(2)).size());
        truncateTo(at("cut_1.fq.gz"), c1);
        truncateTo(at("cut_2.fq.gz"), c2);
        // Both mates losing the same records keeps the pair counts equal, so the pairing
        // check alone cannot notice.
        expectFail("T18 both mates cut at the same record", load({pair(at("cut_1.fq.gz"), at("cut_2.fq.gz"))}),
                   "unexpected end of file");
    }
    {
        const std::string il = interleaved(kReads / 2, "/1", "/2");
        const size_t c = writeGz(at("cut_il.fq.gz"), il, interleaved(700, "/1", "/2").size());
        truncateTo(at("cut_il.fq.gz"), c);
        expectFail("T18 interleaved gzip cut at a pair boundary", load({inter(at("cut_il.fq.gz"))}),
                   "unexpected end of file");
    }
    {
        writeGz(at("mid.fq.gz"), se);
        truncateTo(at("mid.fq.gz"), sizeOf(at("mid.fq.gz")) / 2);
        // Refused by the release too, but as "truncated FASTQ record" or "quality length
        // differs"; the gzip status is the cause and is now what is reported.
        expectFail("T18 FASTQ gzip cut at an arbitrary byte", load({single(at("mid.fq.gz"))}),
                   "unexpected end of file");
    }
    {
        // Two concatenated members are one valid gzip file (gzip -t accepts it).
        const std::string first = fastq(0, 800, plainName), rest = fastq(800, kReads, plainName);
        writeGz(at("m1.gz"), first);
        writeGz(at("m2.gz"), rest);
        std::ifstream a(at("m1.gz"), std::ios::binary), b(at("m2.gz"), std::ios::binary);
        std::ofstream o(at("multi.fq.gz"), std::ios::binary);
        o << a.rdbuf() << b.rdbuf();
        o.close();
        const Result r = load({single(at("multi.fq.gz"))});
        expectOk("T18 control: two gzip members", r, kReads);
        check(r.digest == ref.digest, "T18 control: two gzip members == plain reads");
    }
    {
        // zlib ignores non-gzip bytes after a complete member with Z_OK (tape padding).
        writeGz(at("pad.fq.gz"), se);
        std::ofstream o(at("pad.fq.gz"), std::ios::binary | std::ios::app);
        o << std::string(512, '\0');
        o.close();
        const Result r = load({single(at("pad.fq.gz"))});
        expectOk("T18 control: gzip followed by 512 NUL bytes", r, kReads);
        check(r.digest == ref.digest, "T18 control: padded gzip == plain reads");
    }
    {
        writeGz(at("empty.fq.gz"), "");
        expectOk("T18 control: empty gzip member", load({single(at("empty.fq.gz"))}), 0);
        writeText(at("zero.fq"), "");
        expectOk("T18 control: zero-byte file", load({single(at("zero.fq"))}), 0);
    }

    // ---------------------------------------------------------------- T33: blank lines
    {
        writeText(at("trail.fq"), se + "\n");
        const Result r = load({single(at("trail.fq"))});
        expectOk("T33 trailing blank line", r, kReads);
        check(r.digest == ref.digest, "T33 trailing blank line: same reads as plain");
        counterCheck(r.stats.interRecordBlankLines == 1, "T33 trailing blank line: counter 1");
    }
    {
        writeText(at("trail_crlf.fq"), fastq(0, kReads, plainName, "\r\n") + "\r\n");
        const Result r = load({single(at("trail_crlf.fq"))});
        expectOk("T33 CRLF file with a trailing blank line", r, kReads);
        check(r.digest == ref.digest, "T33 CRLF trailing blank: same reads as plain");
    }
    {
        writeText(at("trail2.fq"), se + "\n\n");
        const Result r = load({single(at("trail2.fq"))});
        expectOk("T33 two trailing blank lines", r, kReads);
        counterCheck(r.stats.interRecordBlankLines == 2, "T33 two trailing blank lines: counter 2");
    }
    {
        writeText(at("mid_blank.fq"), fastq(0, 1000, plainName) + "\n" + fastq(1000, kReads, plainName));
        const Result r = load({single(at("mid_blank.fq"))});
        expectOk("T33 blank line between records", r, kReads);
        check(r.digest == ref.digest, "T33 blank line between records: same reads as plain");
        counterCheck(r.stats.interRecordBlankLines == 1, "T33 blank line between records: counter 1");
    }
    {
        writeGz(at("trail.fq.gz"), se + "\n");
        const Result r = load({single(at("trail.fq.gz"))});
        expectOk("T33 gzip with a trailing blank line", r, kReads);
        check(r.digest == ref.digest, "T33 gzip trailing blank: same reads as plain");
    }
    {
        writeText(at("tp_1.fq"), fastq(0, kReads, slashName(1)) + "\n");
        writeText(at("tp_2.fq"), fastq(0, kReads, slashName(2)) + "\n");
        const Result r = load({pair(at("tp_1.fq"), at("tp_2.fq"))});
        expectOk("T33 paired files with trailing blank lines", r, 2 * kReads);
        check(r.digest == refPe.digest, "T33 paired trailing blank: same reads as /1 /2 pair");
        counterCheck(r.stats.interRecordBlankLines == 2, "T33 paired trailing blank: counter 2");
    }
    {
        writeText(at("lead.fq"), "\n\n" + se);
        const Result r = load({single(at("lead.fq"))});
        expectOk("T33 control: leading blank lines (always accepted)", r, kReads);
        counterCheck(r.stats.interRecordBlankLines == 0, "T33 control: leading blank lines are not counted");
    }
    {
        std::string noNl = se;
        noNl.pop_back();
        writeText(at("nonl.fq"), noNl);
        const Result r = load({single(at("nonl.fq"))});
        expectOk("T33 control: no final newline", r, kReads);
        check(r.digest == ref.digest, "T33 control: no final newline: same reads");
    }
    {
        writeText(at("zlen.fq"), fastq(0, 10, plainName) + "@empty\n\n+\n\n" + fastq(10, 20, plainName));
        expectOk("T33 control: a zero-length record is still a record", load({single(at("zlen.fq"))}), 21);
    }

    // ---------------------------------------------------------------- T33: '.1' / '.2' mate names
    auto spotName = [](const char* mate, const char* desc) {
        return [mate, desc](size_t i) { return "SRR1." + std::to_string(i + 1) + mate + desc; };
    };
    {
        writeText(at("dot_1.fq"), fastq(0, kReads, spotName(".1", " 1 length=100")));
        writeText(at("dot_2.fq"), fastq(0, kReads, spotName(".2", " 1 length=100")));
        const Result r = load({pair(at("dot_1.fq"), at("dot_2.fq"))});
        expectOk("T33 fastq-dump --readids mates (SRR.spot.1 / .2)", r, 2 * kReads);
        check(r.digest == refPe.digest, "T33 --readids: same reads as the /1 /2 pair");
        counterCheck(r.stats.dotSuffixMatePairs == kReads, "T33 --readids: counter = every pair");
    }
    {
        writeText(at("dotn_1.fq"), fastq(0, kReads, spotName(".1", "")));
        writeText(at("dotn_2.fq"), fastq(0, kReads, spotName(".2", "")));
        expectOk("T33 --readids mates without a description", load({pair(at("dotn_1.fq"), at("dotn_2.fq"))}),
                 2 * kReads);
    }
    {
        writeText(at("dot_il.fq"), interleaved(kReads / 2, ".1", ".2"));
        const Result r = load({inter(at("dot_il.fq"))});
        expectOk("T33 --readids mates, interleaved", r, kReads);
        counterCheck(r.stats.dotSuffixMatePairs == kReads / 2, "T33 interleaved --readids: counter = every pair");
    }
    {
        // fastq-dump without --readids: both mates are SRR1.<spot>. Spots 1 and 2 end in
        // '.1' / '.2' themselves; identical names still match and are not counted.
        writeText(at("spot_1.fq"), fastq(0, kReads, spotName("", " 1 length=100")));
        writeText(at("spot_2.fq"), fastq(0, kReads, spotName("", " 1 length=100")));
        const Result r = load({pair(at("spot_1.fq"), at("spot_2.fq"))});
        expectOk("T33 control: fastq-dump names without --readids", r, 2 * kReads);
        counterCheck(r.stats.dotSuffixMatePairs == 0, "T33 control: identical names are not counted as dot pairs");
    }
    {
        auto casava = [](int mate) {
            return [mate](size_t i) { return "M01:1:FC:1:1:" + std::to_string(i) + " " + std::to_string(mate) + ":N:0:ACGT"; };
        };
        writeText(at("cas_1.fq"), fastq(0, kReads, casava(1)));
        writeText(at("cas_2.fq"), fastq(0, kReads, casava(2)));
        expectOk("T33 control: CASAVA 1.8 names", load({pair(at("cas_1.fq"), at("cas_2.fq"))}), 2 * kReads);
    }
    // Negative controls: these must keep failing.
    {
        writeText(at("off_2.fq"), fastq(0, kReads, [](size_t i) { return "SRR1." + std::to_string(i + 2) + ".2"; }));
        expectFail("T33 negative: --readids mates offset by one spot",
                   load({pair(at("dotn_1.fq"), at("off_2.fq"))}), "out of sync");
    }
    {
        writeText(at("x_1.fq"), fastq(0, kReads, plainName));
        writeText(at("x_2.fq"), fastq(0, kReads, [](size_t i) { return "q" + std::to_string(i) + ".1"; }));
        expectFail("T33 negative: 'X' against 'X.1' is not a mate pair",
                   load({pair(at("x_1.fq"), at("x_2.fq"))}), "out of sync");
    }
    {
        writeText(at("sw_2.fq"), fastq(0, kReads, [](size_t i) { return "r" + std::to_string(i + 1) + "/2"; }));
        expectFail("T33 negative: mates shifted by one record",
                   load({pair(at("p_1.fq"), at("sw_2.fq"))}), "out of sync");
    }
    {
        writeText(at("junk.fq"), fastq(0, 10, plainName) + "\nnot a header\n" + fastq(10, 20, plainName));
        expectFail("T33 negative: text that is not a header after a blank line",
                   load({single(at("junk.fq"))}), "expected a '@' header line");
    }
    {
        writeText(at("short.fq"), fastq(0, 10, plainName) + "@cut\nACGT\n+\n");
        expectFail("T33 negative: a record missing its quality line", load({single(at("short.fq"))}),
                   "truncated FASTQ record");
    }

    // ---------------------------------------------------------------- T32 (library level)
    expectFail("T32 the same file as both mates", load({pair(at("p_1.fq"), at("p_1.fq"))}), "same input file");
    {
        std::filesystem::create_symlink(at("p_1.fq"), at("link_2.fq"));
        expectFail("T32 a symlink to R1 given as R2", load({pair(at("p_1.fq"), at("link_2.fq"))}),
                   "same input file");
        expectFail("T32 another spelling of the path", load({pair(at("p_1.fq"), dir + "/./p_1.fq")}),
                   "same input file");
    }
    expectFail("T32 the same file in two single-end libraries",
               load({single(at("se.fq")), single(at("se.fq"))}), "same input file");
    expectFail("T32 a pair file repeated as a single-end library",
               load({pair(at("p_1.fq"), at("p_2.fq")), single(at("p_2.fq"))}), "same input file");
    {
        ts::Library l = pair(at("p_1.fq"), at("p_2.fq"));
        l.interleaved = true;
        expectFail("T32 interleaved library that also names r2", load({l}), "interleaved");
    }
    {
        const Result r = load({pair(at("p_1.fq"), at("p_2.fq")), single(at("se.fq"))});
        expectOk("T32 control: distinct files in a pair and a single-end library", r, 3 * kReads);
        counterCheck(r.stats.files == 3 && r.stats.records == 3 * kReads, "T32 control: counters files=3 records=6000");
    }

    std::filesystem::remove_all(dir);
    if (failures) {
        std::fprintf(stderr, "test_v3_io_seqio: %d of %d checks FAILED\n", failures, checks);
        return 1;
    }
    std::printf("test_v3_io_seqio: %d checks passed\n", checks);
    return 0;
}
