// T30 (build_v3, G-graph): the correction fix cap and TESSERACT_FIX_EC_CAP_MASK.
//
// Ported from the T30 verification test (combo3/verify/T30). maxFixes = max(2, len/10) is
// one budget shared by the right and left walks, and both loop guards test it. When the
// budget ran out the right walk ended with nothing masked and the left walk never ran, so
// raw error bases on both sides were neither corrected nor masked -- contrary to the
// masking policy for a stalled walk. With the fix a spent walk keeps checking observed
// k-mers without substituting, and the first unsupported one is a stall.
//
// k=21, 150-bp read (maxFixes 15), solid set = the true genome.
//   A: 1 correctable left error + 18 isolated right errors (cap hit).
//   D: 18 isolated right errors only (cap hit).
//   B: control, 15 errors in total (fully corrected).
//   C: control, a real stall at 100/101 (masks 100..149).
// Uses only release APIs and the environment, so the same source fails on release 1.3.0.
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <random>
#include <string>
#include <unistd.h>
#include <vector>

#include "correct.h"

extern char** environ;

namespace {

int failures = 0;
void check(bool ok, const std::string& what) {
    std::printf("  %s  %s\n", ok ? "PASS" : "FAIL", what.c_str());
    if (!ok) ++failures;
}

void clearTesseractEnv() {
    std::vector<std::string> names;
    for (char** e = environ; e && *e; ++e) {
        if (std::strncmp(*e, "TESSERACT_", 10) == 0) {
            const char* eq = std::strchr(*e, '=');
            names.emplace_back(*e, eq ? static_cast<size_t>(eq - *e) : std::strlen(*e));
        }
    }
    for (const std::string& n : names) unsetenv(n.c_str());
}

std::string tmpDir;
std::mt19937 generator(20260925);
std::string randomSeq(int n) {
    std::string s(static_cast<size_t>(n), 'A');
    for (char& c : s) c = "ACGT"[generator() % 4];
    return s;
}
void addAll(ts::KmerTable& solid, const std::string& s, int k) {
    ts::Kmer f = 0;
    for (size_t i = 0; i < s.size(); ++i) {
        f = ts::pushBack(f, ts::baseCode(s[i]), k);
        if (static_cast<int>(i) >= k - 1) solid.put(ts::canonical(f, k), 30);
    }
}
ts::SequenceStore load(const std::string& s) {
    const std::string path = tmpDir + "/reads.fa";
    std::ofstream(path) << ">r\n" << s << "\n";
    ts::Library library;
    library.r1 = path;
    ts::SequenceStore reads;
    std::string error;
    if (!reads.load({library}, 1, error)) {
        std::fprintf(stderr, "load failed: %s\n", error.c_str());
        std::exit(2);
    }
    return reads;
}
int nonSolidKmers(const std::string& s, const ts::KmerTable& solid, int k) {
    int bad = 0, run = 0;
    ts::Kmer f = 0;
    for (size_t i = 0; i < s.size(); ++i) {
        const int c = ts::baseCode(s[i]);
        if (c < 0) { run = 0; f = 0; continue; }
        f = ts::pushBack(f, c, k);
        ++run;
        if (run >= k && !solid.contains(ts::canonical(f, k))) ++bad;
    }
    return bad;
}
struct Result { std::string out; ts::CorrectionStats st; };
Result run(const std::string& truth, const std::vector<int>& errs, const ts::KmerTable& solid, int k) {
    std::string raw = truth;
    for (int p : errs) raw[static_cast<size_t>(p)] = ts::codeBase((ts::baseCode(raw[static_cast<size_t>(p)]) + 1) & 3);
    auto reads = load(raw);
    Result r;
    r.st = ts::correctReads(reads, solid, k, 1, 8);
    r.out = reads.decode(0);
    return r;
}
std::vector<int> errsA() { std::vector<int> v = {8}; for (int p = 60; p <= 145; p += 5) v.push_back(p); return v; }
std::vector<int> errsB() { std::vector<int> v = {8}; for (int p = 60; p <= 125; p += 5) v.push_back(p); return v; }
std::vector<int> errsC() { return {8, 60, 65, 100, 101, 120}; }
std::vector<int> errsD() { std::vector<int> v; for (int p = 60; p <= 145; p += 5) v.push_back(p); return v; }

struct Four { Result a, b, c, d; };
Four runAll(const std::string& truth, const ts::KmerTable& solid, int k) {
    return {run(truth, errsA(), solid, k), run(truth, errsB(), solid, k), run(truth, errsC(), solid, k),
            run(truth, errsD(), solid, k)};
}

}  // namespace

int main() {
    clearTesseractEnv();
    // 1.4.0: the umbrella is on by default (src/defaults.h); "flags unset" below means the
    // release 1.3.0 state, so the umbrella is pinned off until a check turns it on.
    setenv("TESSERACT_FIXES", "0", 1);
    char tmpl[] = "/tmp/tess_v3_eccap.XXXXXX";
    const char* d = mkdtemp(tmpl);
    if (!d) return 2;
    tmpDir = d;

    const int k = 21;
    const std::string genome = randomSeq(400);
    ts::KmerTable solid;
    addAll(solid, genome, k);
    const std::string truth = genome.substr(50, 150);

    std::printf("flags unset (release behaviour):\n");
    const Four off = runAll(truth, solid, k);
    check(off.a.st.basesCorrected == 15 && off.a.out[8] != truth[8] && off.a.out[8] != 'N' &&
              off.a.st.basesMasked == 0,
          "A: cap hit leaves the left error and the right tail unmasked (release, pinned)");
    check(nonSolidKmers(off.a.out, solid, k) > 0 && nonSolidKmers(off.d.out, solid, k) > 0,
          "A, D: non-solid k-mers survive correction (the defect)");
    check(off.b.out == truth && off.c.st.basesMasked == 50, "controls B (fully corrected) and C (stall masks 50)");

    setenv("TESSERACT_FIX_EC_CAP_MASK", "1", 1);
    std::printf("TESSERACT_FIX_EC_CAP_MASK=1:\n");
    const Four on = runAll(truth, solid, k);
    check(nonSolidKmers(on.a.out, solid, k) == 0 && nonSolidKmers(on.d.out, solid, k) == 0,
          "A, D: no retained non-solid k-mer after the cap");
    check(on.a.st.basesCorrected == 15 && on.a.out.substr(0, 9) == std::string(9, 'N') &&
              on.a.out.substr(135) == std::string(15, 'N'),
          "A: same 15 fixes; 0..8 and 135..149 masked");
    check(on.d.out.substr(135) == std::string(15, 'N') && on.d.out.substr(0, 135) == truth.substr(0, 135),
          "D: 135..149 masked, the rest corrected");
    check(on.b.out == off.b.out && on.c.out == off.c.out && on.c.st.basesMasked == off.c.st.basesMasked,
          "controls B and C unchanged");

    unsetenv("TESSERACT_FIX_EC_CAP_MASK");
    setenv("TESSERACT_FIXES", "1", 1);
    const Four um = runAll(truth, solid, k);
    check(um.a.out == on.a.out && um.d.out == on.d.out, "umbrella TESSERACT_FIXES=1 enables the fix");
    setenv("TESSERACT_FIX_EC_CAP_MASK", "0", 1);
    const Four ov = runAll(truth, solid, k);
    check(ov.a.out == off.a.out && ov.d.out == off.d.out, "own flag 0 overrides the umbrella");

    std::printf("test_v3_graph_ec_cap: %s (%d failure%s)\n", failures ? "FAILED" : "ALL PASSED", failures,
                failures == 1 ? "" : "s");
    const std::string cmd = "rm -rf '" + tmpDir + "'";
    if (std::system(cmd.c_str()) != 0) std::printf("note: could not remove %s\n", tmpDir.c_str());
    return failures ? 1 : 0;
}
