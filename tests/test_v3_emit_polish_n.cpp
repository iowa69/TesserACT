// G-emit T02: the consensus polisher must not write read bases over scaffold N-runs when
// TESSERACT_FIX_POLISH_SKIP_N=1 (or the umbrella TESSERACT_FIXES=1) is set.
//
// Case A (merge): truth L + X(50) + R, scaffold L + N + R (a gap estimate <= 1 becomes one N,
//   resolve.cpp joinGap = max(1, gap)). Reads overhanging the N by 1..7 bases vote X[0] from the
//   left and X[49] from the right; both are 'T'. Release writes 'T': one contiguous record with a
//   49-bp deletion instead of two records split at the N.
// Case B (island): scaffold L + N*40 + R, X[0] split 50/50 across reads. Release keeps X[0] as N
//   but overwrites X[1..3], so splitting leaves a 3-bp record between two N-runs.
//
// With the switch unset the release behaviour is asserted (so the flag-off path stays the
// release one); with it on, every N must survive. Against the release polish.cpp the flag-on
// assertions FAIL (build with -DT_RELEASE_API to compile against the release headers).
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <random>
#include <string>
#include <vector>
#include <unistd.h>

#include "polish.h"
#include "seqio.h"
#include "test_env.h"   // build_v3: T40, clear the ambient TESSERACT_* first

namespace {
int failures = 0, checks = 0;
void check(bool ok, const std::string& what) {
    ++checks;
    std::printf("%s %s\n", ok ? "ok  " : "FAIL", what.c_str());
    if (!ok) ++failures;
}
std::mt19937 rng(20260925);
std::string randomSeq(size_t n) {
    std::string s(n, 'A');
    for (char& c : s) c = "ACGT"[rng() % 4];
    return s;
}
ts::SequenceStore load(const std::vector<std::string>& reads, const char* tag) {
    const auto dir = std::filesystem::temp_directory_path() /
                     (std::string("v3emit-t02-") + tag + "-" + std::to_string(getpid()));
    std::filesystem::create_directories(dir);
    const auto p = dir / "r.fa";
    {
        std::ofstream f(p);
        for (size_t i = 0; i < reads.size(); ++i) f << ">" << i << "\n" << reads[i] << "\n";
    }
    ts::SequenceStore st;
    ts::Library lib;
    lib.r1 = p.string();
    std::string err;
    if (!st.load({lib}, 1, err)) { std::fprintf(stderr, "load failed: %s\n", err.c_str()); std::exit(2); }
    std::filesystem::remove_all(dir);
    return st;
}
size_t countN(const std::string& s) {
    size_t n = 0;
    for (char c : s) n += (c == 'N');
    return n;
}
std::vector<size_t> pieceLengths(const std::string& s) {   // the contigs.fasta split rule
    std::vector<size_t> out;
    size_t p = 0;
    while (p < s.size()) {
        while (p < s.size() && s[p] == 'N') ++p;
        if (p >= s.size()) break;
        size_t e = p;
        while (e < s.size() && s[e] != 'N') ++e;
        out.push_back(e - p);
        p = e;
    }
    return out;
}
const int kRead = 150, kCopies = 4, kMaxOverhang = 7;

void clearEnv() {
    unsetenv("TESSERACT_FIXES");
    unsetenv("TESSERACT_FIX_POLISH_SKIP_N");
    unsetenv("TESSERACT_POLISH_ORIGINAL_QUALITY");
    // 1.4.0: an unset umbrella is on (src/defaults.h); the release 1.3.0 state is TESSERACT_FIXES=0.
    setenv("TESSERACT_FIXES", "0", 1);
}

struct Fixture { std::string contig; ts::SequenceStore reads; };

Fixture caseA() {
    const std::string L = randomSeq(600), R = randomSeq(600);
    std::string X = randomSeq(50);
    X[0] = 'T'; X.back() = 'T';
    const std::string truth = L + X + R;
    std::vector<std::string> reads;
    for (int c = 0; c < kCopies; ++c)
        for (int o = 1; o <= kMaxOverhang; ++o) {
            reads.push_back(truth.substr(L.size() + o - kRead, kRead));
            reads.push_back(truth.substr(L.size() + X.size() - o, kRead));
        }
    return {L + "N" + R, load(reads, "a")};
}

Fixture caseB() {
    const std::string L = randomSeq(600), R = randomSeq(600);
    const std::string X = randomSeq(40);
    std::vector<std::string> reads;
    for (int c = 0; c < kCopies; ++c)
        for (int o = 1; o <= kMaxOverhang; ++o) {
            std::string xv = X;
            xv[0] = (c % 2) ? 'A' : 'C';
            const std::string truth = L + xv + R;
            reads.push_back(truth.substr(L.size() + o - kRead, kRead));
            reads.push_back(truth.substr(L.size() + X.size() - o, kRead));
        }
    return {L + std::string(40, 'N') + R, load(reads, "b")};
}

ts::PolishStats polishOnce(const Fixture& f, std::string& out) {
    std::vector<std::string> seqs{f.contig};
    const ts::PolishStats ps = ts::polishContigs(seqs, f.reads, 1, 31, 15, 0.90);
    out = seqs[0];
    return ps;
}
}  // namespace

int main() {
    testenv::clearTesseractEnv();
    const Fixture a = caseA();
    const Fixture b = caseB();
    std::string out;

    // ---- switch unset: the release behaviour, unchanged --------------------------------
    clearEnv();
    {
        const auto ps = polishOnce(a, out);
        check(countN(out) == 0 && pieceLengths(out).size() == 1 && out.size() == 1201,
              "A off: release overwrites the 1-N join (one 1201-bp record)");
#ifndef T_RELEASE_API
        check(!ps.nRuns.skipN && ps.nRuns.replaced == 1 && ps.nRuns.wouldReplace == 1 &&
                  ps.nRuns.runs == 1 && ps.nRuns.runsFullyReplaced == 1,
              "A off: counters replaced=1 runs_fully_replaced=1");
#endif
        (void)ps;
    }
    {
        const auto ps = polishOnce(b, out);
        const auto pl = pieceLengths(out);
        bool island = false;
        for (size_t l : pl) island |= (l <= 3);
        check(countN(out) < 40 && island, "B off: release leaves a <=3-bp island inside the 40-N run");
#ifndef T_RELEASE_API
        check(ps.nRuns.replaced > 0 && ps.nRuns.runsTouched == 1 && ps.nRuns.runsFullyReplaced == 0,
              "B off: counters replaced>0 runs_touched=1");
#endif
        (void)ps;
    }

    // ---- TESSERACT_FIX_POLISH_SKIP_N=1 -------------------------------------------------
    clearEnv();
    setenv("TESSERACT_FIX_POLISH_SKIP_N", "1", 1);
    {
        const auto ps = polishOnce(a, out);
        const auto pl = pieceLengths(out);
        check(countN(out) == 1 && pl.size() == 2 && pl[0] == 600 && pl[1] == 600,
              "A on: the N survives, pieces 600+600");
#ifndef T_RELEASE_API
        check(ps.nRuns.skipN && ps.nRuns.replaced == 0 && ps.nRuns.wouldReplace == 1 &&
                  ps.nRuns.runsFullyReplaced == 0,
              "A on: counters would_replace=1 replaced=0");
#endif
        (void)ps;
    }
    {
        const auto ps = polishOnce(b, out);
        const auto pl = pieceLengths(out);
        check(countN(out) == 40 && pl.size() == 2 && pl[0] == 600 && pl[1] == 600,
              "B on: all 40 N survive, pieces 600+600, no island");
#ifndef T_RELEASE_API
        check(ps.nRuns.replaced == 0 && ps.nRuns.wouldReplace > 0 && ps.nRuns.runsTouched == 0,
              "B on: counters replaced=0 runs_touched=0");
#endif
        (void)ps;
    }

    // ---- the umbrella turns it on; an explicit 0 wins over the umbrella ------------------
    clearEnv();
    setenv("TESSERACT_FIXES", "1", 1);
    polishOnce(a, out);
    check(countN(out) == 1, "umbrella TESSERACT_FIXES=1: the N survives");
    setenv("TESSERACT_FIX_POLISH_SKIP_N", "0", 1);
    polishOnce(a, out);
    check(countN(out) == 0, "TESSERACT_FIXES=1 with TESSERACT_FIX_POLISH_SKIP_N=0: release behaviour");
    clearEnv();

    std::printf("test_v3_emit_polish_n: %d checks, %d failures\n", checks, failures);
    return failures ? 1 : 0;
}
