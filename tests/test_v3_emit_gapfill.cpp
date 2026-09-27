// G-emit gap-filler fixes: T09 (truncated search accepted as unique), T11 (no anchor back-off),
// T27 (input N read as 'A'), T38 (floor-2 walk run twice).
//
// Every fixture is run with the switches unset -- where the RELEASE outcome is asserted, so
// the flag-off path is pinned to release -- and with the switch on, where the fixed outcome is
// asserted. Against the release gapfill.cpp the flag-on assertions FAIL; build with
// -DT_RELEASE_API to compile against the release headers (the counter checks drop out).
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <functional>
#include <random>
#include <string>
#include <utility>
#include <vector>
#include <unistd.h>

#include "gapfill.h"
#include "kmer.h"
#include "seqio.h"
#include "test_env.h"   // build_v3: T40, clear the ambient TESSERACT_* first

namespace {
int checks = 0, failures = 0;
void check(bool ok, const std::string& what) {
    ++checks;
    std::printf("%s %s\n", ok ? "ok  " : "FAIL", what.c_str());
    if (!ok) ++failures;
}
std::string dna(size_t n, unsigned seed) {
    std::mt19937 rng(seed);
    std::string s(n, 'A');
    for (char& b : s) b = "ACGT"[rng() % 4];
    return s;
}
std::string rc(const std::string& s) {
    std::string out;
    for (auto i = s.rbegin(); i != s.rend(); ++i) out += "TGCA"[ts::baseCode(*i)];
    return out;
}
char other(char c, char avoid = 0) {
    for (char x : std::string("ACGT")) if (x != c && x != avoid) return x;
    return 'A';
}
std::string mutate(std::string s, size_t pos) { s[pos] = other(s[pos]); return s; }

void clearEnv() {
    for (const char* v : {"TESSERACT_FIXES", "TESSERACT_FIX_GAPFILL_STRICT_BUDGET",
                          "TESSERACT_FIX_GAPFILL_BACKOFF", "TESSERACT_FIX_GAPFILL_SKIP_INPUT_N",
                          "TESSERACT_GF_DEBUG"})
        unsetenv(v);
}

std::string scratch() {
    static int n = 0;
    const auto dir = std::filesystem::temp_directory_path() /
                     ("v3emit-gf-" + std::to_string(getpid()) + "-" + std::to_string(n++));
    std::filesystem::create_directories(dir);
    return dir.string();
}
ts::SequenceStore loadSingle(const std::vector<std::string>& reads) {
    const std::string dir = scratch();
    {
        std::ofstream f(dir + "/r.fa");
        for (size_t i = 0; i < reads.size(); ++i) f << ">r" << i << "\n" << reads[i] << "\n";
    }
    ts::SequenceStore st;
    ts::Library lib;
    lib.r1 = dir + "/r.fa";
    std::string err;
    if (!st.load({lib}, 1, err)) { std::fprintf(stderr, "load: %s\n", err.c_str()); std::exit(2); }
    std::filesystem::remove_all(dir);
    return st;
}
using Pair = std::pair<std::string, std::string>;
void tile(std::vector<Pair>& out, const std::string& g, size_t from, size_t to, size_t step,
          const std::function<bool(size_t)>& keep = nullptr, size_t len = 100, size_t frag = 300) {
    for (size_t p = from; p + frag <= to && p + frag <= g.size(); p += step) {
        if (keep && !keep(p)) continue;
        out.push_back({g.substr(p, len), rc(g.substr(p + frag - len, len))});
    }
}
ts::SequenceStore loadPairs(const std::vector<Pair>& pairs) {
    const std::string dir = scratch();
    {
        std::ofstream f(dir + "/r1.fa"), s(dir + "/r2.fa");
        for (size_t i = 0; i < pairs.size(); ++i) {
            f << ">" << i << "\n" << pairs[i].first << "\n";
            s << ">" << i << "\n" << pairs[i].second << "\n";
        }
    }
    ts::SequenceStore st;
    ts::Library lib;
    lib.r1 = dir + "/r1.fa";
    lib.r2 = dir + "/r2.fa";
    std::string err;
    if (!st.load({lib}, 1, err)) { std::fprintf(stderr, "load: %s\n", err.c_str()); std::exit(2); }
    std::filesystem::remove_all(dir);
    return st;
}
struct Run { std::vector<std::string> c; ts::GapFillStats st; };
Run close(const std::vector<std::string>& contigs, const ts::SequenceStore& reads, int threads = 2,
          const std::vector<uint8_t>* prov = nullptr) {
    Run r{contigs, {}};
#ifdef T_RELEASE_API
    (void)prov;
    r.st = ts::closeGaps(r.c, reads, threads, 31, 300);
#else
    r.st = ts::closeGaps(r.c, reads, threads, 31, 300, prov);
#endif
    return r;
}

// ---------------------------------------------------------------------------------------
// T09: seed -A(50x)-> target, seed -B(20x)-> target, and a dead-end side branch leaving A
// halfway made of NB SNP bubbles (2^NB paths that never reach the target). Strongest-first
// DFS finds A, then spends the whole budget in the side branch and never sees B.
struct T09Case { std::string contig, withA, withB; ts::SequenceStore reads; };
T09Case t09(bool sideBranch, int NB) {
    std::mt19937 rng(90909);
    auto rnd = [&](size_t n) { std::string s(n, 'A'); for (char& c : s) c = "ACGT"[rng() % 4]; return s; };
    const std::string L = rnd(400), R = rnd(400);
    const std::string a1 = rnd(100);
    std::string b1 = rnd(100);
    if (b1[0] == a1[0]) b1[0] = other(a1[0]);
    const char e = other(a1[50]);
    std::string hap1 = rnd(static_cast<size_t>(NB) * 32 + 20), hap2 = hap1;
    for (int i = 0; i < NB; ++i) {
        const size_t p = 16 + static_cast<size_t>(i) * 32;
        hap2[p] = other(hap1[p]);
    }
    std::vector<std::string> reads;
    const std::string lf = L.substr(L.size() - 60), rf = R.substr(0, 60);
    for (int i = 0; i < 50; ++i) reads.push_back(lf + a1 + rf);
    for (int i = 0; i < 20; ++i) reads.push_back(lf + b1 + rf);
    if (sideBranch) {
        for (int i = 0; i < 15; ++i) reads.push_back(lf + a1.substr(0, 50) + e + hap1);
        for (int i = 0; i < 15; ++i) reads.push_back(lf + a1.substr(0, 50) + e + hap2);
    }
    return {L + std::string(100, 'N') + R, L + a1 + R, L + b1 + R, loadSingle(reads)};
}

// T11 paralog fixture (the efaecium GCF900634805v1 shape): the contig's left flank ends in a
// 10-bp tip = 5 bases of a paralog's continuation + 5 error bases, spelled by 2 reads. The
// first supported k-mer stepping left is on the paralog branch; the true locus branches off
// 5 bases further back.
struct ParalogCase { std::string contig, truthClosed; ts::SequenceStore reads; };
ParalogCase paralog() {
    const std::string U = dna(700, 101), S = dna(60, 102), V = dna(700, 103);
    std::string T = dna(900, 104), W = dna(600, 105);
    if (W[0] == T[0]) W[0] = other(T[0]);
    std::string err = dna(5, 106);
    if (err[0] == W[5]) err[0] = other(W[5]);
    const std::string trueLocus = U + S + T;            // ... gap is T[0,80) ...
    const std::string paraLocus = V + S + W;
    const size_t gapLen = 80;
    const std::string tip = W.substr(0, 5) + err;
    const std::string contig = U + S + tip + std::string(gapLen, 'N') + T.substr(gapLen);
    std::vector<Pair> pairs;
    tile(pairs, trueLocus, 0, trueLocus.size(), 6);
    tile(pairs, paraLocus, 0, paraLocus.size(), 3);
    // the tip reads: end exactly at the tip, mates upstream in U
    const std::string tipped = U + S + tip;
    for (int i = 0; i < 2; ++i)
        pairs.push_back({rc(tipped.substr(tipped.size() - 100, 100)), tipped.substr(tipped.size() - 400, 100)});
    return {contig, U + S + T, loadPairs(pairs)};
}
}  // namespace

int main() {
    testenv::clearTesseractEnv();
    clearEnv();

    // ================================ T09 ==============================================
    {
        const T09Case ctl = t09(false, 14), tr = t09(true, 14), small = t09(true, 8);
        clearEnv();
        Run c = close({ctl.contig}, ctl.reads, 1);
        check(c.st.gapsClosed == 0 && c.st.gapsAmbiguous == 1, "T09 off control: complete search -> ambiguous");
        Run t = close({tr.contig}, tr.reads, 1);
        check(t.st.gapsClosed == 1 && t.c[0] == tr.withA,
              "T09 off: truncated search closes with A (release behaviour)");
#ifndef T_RELEASE_API
        check(t.st.gapsTruncated == 1 && t.st.gapsTruncatedAccepted == 1 && t.st.gapsTruncatedRefused == 0 &&
                  !t.st.strictBudget, "T09 off: counters truncated=1 truncAccepted=1");
#endif
        Run s = close({small.contig}, small.reads, 1);
        check(s.st.gapsClosed == 0 && s.st.gapsAmbiguous == 1, "T09 off NB=8: search completes -> ambiguous");

        setenv("TESSERACT_FIX_GAPFILL_STRICT_BUDGET", "1", 1);
        Run t2 = close({tr.contig}, tr.reads, 1);
        check(t2.st.gapsClosed == 0 && t2.c[0] == tr.contig && t2.st.gapsAmbiguous == 1,
              "T09 on: truncated 'unique' closure refused, Ns kept, counted ambiguous");
#ifndef T_RELEASE_API
        check(t2.st.strictBudget && t2.st.gapsTruncatedRefused == 1 && t2.st.gapsTruncatedAccepted == 0,
              "T09 on: counters strict=1 truncRefused=1");
#endif
        Run c2 = close({ctl.contig}, ctl.reads, 1);
        Run s2 = close({small.contig}, small.reads, 1);
        check(c2.st.gapsAmbiguous == 1 && s2.st.gapsAmbiguous == 1 && c2.c[0] == ctl.contig,
              "T09 on: complete searches unchanged");
#ifndef T_RELEASE_API
        check(c2.st.gapsTruncated == 0 && s2.st.gapsTruncated == 0, "T09 on: complete searches not truncated");
#endif
        clearEnv();
    }

    // ================================ T38 ==============================================
#ifndef T_RELEASE_API
    {
        // Low depth (floor 2) and nothing spans the gap: the ladder must run floor 2 once.
        const std::string G = dna(1200, 38);
        std::vector<std::string> reads;
        for (int i = 0; i < 2; ++i) reads.push_back(G.substr(350, 150));   // ends at the gap
        for (int i = 0; i < 2; ++i) reads.push_back(G.substr(560, 150));   // starts after it
        const ts::SequenceStore st = loadSingle(reads);
        const std::string contig = G.substr(0, 500) + std::string(60, 'N') + G.substr(560);
        Run r = close({contig}, st, 1);
        check(r.st.gapsClosed == 0 && r.st.gapsNoPath == 1 && r.st.meanFloor == 2.0,
              "T38: floor-2 gap with no path stays open");
        check(r.st.searches == 1, "T38: floor 2 walked once (release: twice)");
        // A floor-3 gap still runs {3,2}.
        std::vector<std::string> deep;
        for (int i = 0; i < 60; ++i) deep.push_back(G.substr(350 + (i % 3), 150));
        for (int i = 0; i < 60; ++i) deep.push_back(G.substr(560 + (i % 3), 150));
        const ts::SequenceStore st3 = loadSingle(deep);
        Run r3 = close({contig}, st3, 1);
        check(r3.st.gapsNoPath == 1 && r3.st.meanFloor >= 3.0 && r3.st.searches == 2,
              "T38: floor-3 gap still walks floors 3 then 2");
    }
#endif

    // ================================ T11 ==============================================
    {
        const std::string G = dna(5000, 7);
        // S1: error tip (substitution 15 bp before the Ns) spelled by 2 reads.
        const std::string Gerr = mutate(G, 1985);
        std::vector<Pair> p1;
        tile(p1, G, 0, 5000, 5);
        for (int i = 0; i < 2; ++i) p1.push_back({Gerr.substr(1900, 100), rc(G.substr(2100, 100))});
        const ts::SequenceStore r1 = loadPairs(p1);
        const std::string s1 = Gerr.substr(0, 2000) + std::string(50, 'N') + G.substr(2050, 1950);
        // S5: the mirror on the target side.
        const std::string Gerr2 = mutate(G, 2065);
        std::vector<Pair> p5;
        tile(p5, G, 0, 5000, 5);
        for (int i = 0; i < 2; ++i) p5.push_back({G.substr(1850, 100), rc(Gerr2.substr(2050, 100))});
        const ts::SequenceStore r5 = loadPairs(p5);
        const std::string s5 = G.substr(0, 2000) + std::string(50, 'N') + Gerr2.substr(2050, 1950);
        // S3: correct thin tip crossed by 2 reads; S4: error tip whose reads run into the gap.
        auto clearTip = [](size_t p) { return (p >= 2060 || p + 100 <= 1960) && (p + 200 >= 2060 || p + 300 <= 1960); };
        std::vector<Pair> p3;
        tile(p3, G, 0, 5000, 5, clearTip);
        for (int i = 0; i < 2; ++i) p3.push_back({G.substr(1950, 150), rc(G.substr(2250, 150))});
        const ts::SequenceStore r3 = loadPairs(p3);
        const std::string s3 = G.substr(0, 2000) + std::string(50, 'N') + G.substr(2050, 1950);
        std::vector<Pair> p4;
        tile(p4, G, 0, 5000, 5);
        for (int i = 0; i < 2; ++i) p4.push_back({Gerr.substr(1930, 150), rc(G.substr(2250, 150))});
        const ts::SequenceStore r4 = loadPairs(p4);
        const std::string s4 = s1;
        const ParalogCase pc = paralog();

        clearEnv();
        Run a1 = close({s1}, r1), a5 = close({s5}, r5), a3 = close({s3}, r3), a4 = close({s4}, r4);
        Run ap = close({pc.contig}, pc.reads);
        check(a1.st.gapsClosed == 0 && a1.st.seedBelowFloor == 1, "T11 off S1: error tip blocks closure (release)");
        check(a5.st.gapsClosed == 0 && a5.st.targetBelowFloor == 1, "T11 off S5: target-side error tip blocks closure");
        check(a3.st.gapsClosed == 1 && a3.c[0] == G.substr(0, 4000), "T11 off S3: correct thin tip closes");
        check(a4.st.gapsClosed == 1 && a4.c[0] == Gerr.substr(0, 4000), "T11 off S4: closes keeping the tip error");
        check(ap.st.gapsClosed == 0 && ap.st.seedBelowFloor == 1, "T11 off paralog: gap stays open (release)");

        setenv("TESSERACT_FIX_GAPFILL_BACKOFF", "60", 1);
        Run b1 = close({s1}, r1), b5 = close({s5}, r5), b3 = close({s3}, r3), b4 = close({s4}, r4);
        Run bp = close({pc.contig}, pc.reads);
        check(b1.st.gapsClosed == 1 && b1.c[0] == G.substr(0, 4000), "T11 on S1: closed and equals the truth");
        check(b5.st.gapsClosed == 1 && b5.c[0] == G.substr(0, 4000), "T11 on S5: closed and equals the truth");
        check(b3.c == a3.c && b4.c == a4.c, "T11 on: release closures (S3, S4) byte-identical");
        check(bp.st.gapsClosed == 1 && bp.c[0] == pc.truthClosed,
              "T11 on paralog: closed through the branch point to the true locus");
#ifndef T_RELEASE_API
        check(b1.st.backoff == 60 && b1.st.closedAfterBackoff == 1 && b1.st.trimmedBp == 15 &&
                  b1.st.seedBackedOff == 1, "T11 on S1: counters closedAfterBackoff=1 trimmedBp=15");
        check(b5.st.closedAfterBackoff == 1 && b5.st.targetBackedOff == 1 && b5.st.trimmedBp == 16,
              "T11 on S5: counters targetBackedOff=1 trimmedBp=16");
        check(b3.st.closedAfterBackoff == 0 && b4.st.closedAfterBackoff == 0 && b3.st.backoffAttempts == 0,
              "T11 on: no back-off where release closed");
        check(bp.st.closedAfterBackoff == 1 && bp.st.trimmedBp == 10 && bp.st.backoffAttempts >= 2,
              "T11 on paralog: the paralog-branch anchor fails first, trimmedBp=10");
#endif
        // "1" means the default budget; the umbrella turns it on too.
        setenv("TESSERACT_FIX_GAPFILL_BACKOFF", "1", 1);
        Run d1 = close({s1}, r1);
        check(d1.c[0] == G.substr(0, 4000), "T11 BACKOFF=1 (default budget) closes S1");
        clearEnv();
        setenv("TESSERACT_FIXES", "1", 1);
        Run u1 = close({s1}, r1);
        check(u1.c[0] == G.substr(0, 4000), "T11 umbrella closes S1");
        setenv("TESSERACT_FIX_GAPFILL_BACKOFF", "0", 1);
        Run z1 = close({s1}, r1);
        check(z1.c[0] == s1, "T11 umbrella with BACKOFF=0: release behaviour");
        clearEnv();
    }

    // ================================ T27 ==============================================
    {
        std::mt19937 rng(27027);
        std::string G(1200, 'A');
        for (char& c : G) c = "ACGT"[rng() % 4];
        const size_t X = 530;
        G[X] = 'C';
        const std::string scaffold = G.substr(0, 500) + std::string(60, 'N') + G.substr(560);
        auto build = [&](int dupN, int good) {
            std::vector<std::string> reads;
            std::string dup = G.substr(455, 150);
            dup[X - 455] = 'N';
            for (int i = 0; i < dupN; ++i) reads.push_back(dup);
            for (int i = 0; i < good; ++i) reads.push_back(G.substr(450 + 3 * i, 150));
            reads.push_back(G.substr(300, 150));
            reads.push_back(G.substr(650, 150));
            return loadSingle(reads);
        };
        const ts::SequenceStore sa = build(2, 0), sb = build(2, 1), sc = build(2, 3), sk = build(0, 3);
#ifndef T_RELEASE_API
        const auto pa = ts::readsWithInputAmbiguity(sa, 2), pb = ts::readsWithInputAmbiguity(sb, 2);
        const auto pcv = ts::readsWithInputAmbiguity(sc, 2), pk = ts::readsWithInputAmbiguity(sk, 2);
        check(pa.size() == 4 && pa[0] && pa[1] && !pa[2] && !pa[3], "T27 provenance flags exactly the N reads");
        const auto *ppa = &pa, *ppb = &pb, *ppc = &pcv, *ppk = &pk;
#else
        const std::vector<uint8_t> *ppa = nullptr, *ppb = nullptr, *ppc = nullptr, *ppk = nullptr;
#endif
        auto base = [&](const Run& r) { return (r.st.gapsClosed == 1 && r.c[0].size() == G.size()) ? r.c[0][X] : 'N'; };
        clearEnv();
        Run a = close({scaffold}, sa, 1, ppa), b = close({scaffold}, sb, 1, ppb), c = close({scaffold}, sc, 1, ppc);
        Run k = close({scaffold}, sk, 1, ppk);
        check(base(a) == 'A', "T27 off A: closed with 'A' where every covering read has N (release)");
        check(base(b) == 'A', "T27 off B: closed with 'A' although the only observation is 'C' (release)");
        check(c.st.gapsClosed == 0 && c.st.gapsAmbiguous == 1, "T27 off C: phantom 'A' makes it ambiguous (release)");
        check(base(k) == 'C', "T27 off control: closed with 'C'");

        setenv("TESSERACT_FIX_GAPFILL_SKIP_INPUT_N", "1", 1);
        Run a2 = close({scaffold}, sa, 1, ppa), b2 = close({scaffold}, sb, 1, ppb), c2 = close({scaffold}, sc, 1, ppc);
        Run k2 = close({scaffold}, sk, 1, ppk);
        check(a2.st.gapsClosed == 0, "T27 on A: gap stays open");
        check(b2.st.gapsClosed == 0, "T27 on B: gap stays open");
        check(base(c2) == 'C' && c2.c[0] == G, "T27 on C: closed with the observed 'C', equals the truth");
        check(k2.c == k.c, "T27 on control: unchanged");
#ifndef T_RELEASE_API
        check(a2.st.skipInputN && a2.st.inputNProvenance && a2.st.readsInputN == 2 &&
                  a2.st.recruitedInputN == 2 && a2.st.kmersSkippedInputN == 2 * 31,
              "T27 on A: counters reads_input_n=2 kmers_skipped_input_n=62");
        Run np = close({scaffold}, sa, 1, nullptr);
        check(base(np) == 'A' && np.st.skipInputN && !np.st.inputNProvenance,
              "T27 on without provenance: release reading, provenance=0 reported");
#endif
        clearEnv();
    }

    std::printf("test_v3_emit_gapfill: %d checks, %d failures\n", checks, failures);
    return failures ? 1 : 0;
}
