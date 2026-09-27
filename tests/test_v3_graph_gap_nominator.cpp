// T08 (build_v3, G-graph): the legacy gap-close nominator and TESSERACT_FIX_GAP_NOMINATOR.
//
// Ported from the T08 verification probe (combo3/verify/T08). The legacy walk continues
// from the tip's end after a link that flips stored orientation, tests non-branching on
// the far end, and indexes k-mers strand-blind, first writer wins -- so which read pairs
// vote for a tip depends on how the unitigs happen to be STORED.
//
// TEST 1: compacted graph Y->Z, Z->{W,Q}, W->{U,S}, S->Y, Q->Y, V->Y; dead ends U 3' and
//   V 5'. For each of the 16 stored orientations of Z/W/S/U and each source X, 5 FR pairs
//   X<->V (minVotes 2). Only U's reads may vote for (U-tip, V-tip).
//   Release / flag off: 36/64 wrong-source nominations (pinned: flag-off is unchanged).
//   Fix: 0/64 wrong, 0/16 missed -- the same as the oriented nominator.
// TEST 2: isolated unitig X (both ends dead), pairs from X's 3' region to V's 5' end.
//   Release / flag off: the wrong end (X 5') in the 4/8 cases where X is stored reversed.
//   Fix: 8/8 correct.
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <random>
#include <string>
#include <unistd.h>
#include <vector>

#include "gap_evidence.h"
#include "graph.h"
#include "seqio.h"

extern char** environ;
using namespace ts;

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

std::mt19937 rng(7);
std::string rnd(size_t n) {
    std::string s(n, 'A');
    for (auto& c : s) c = "ACGT"[rng() & 3];
    return s;
}
uint64_t endId(uint32_t u, int e) { return (static_cast<uint64_t>(u) << 1) | static_cast<uint64_t>(e); }

struct Bio {  // biological node: sequence + stored flip
    std::string seq;
    bool flip = false;
    uint32_t id = 0;
    int stored(int bioEnd) const { return flip ? 1 - bioEnd : bioEnd; }
};
void blink(UnitigGraph& g, const Bio& a, int aBio, const Bio& b, int bBio) {
    const int ae = a.stored(aBio), be = b.stored(bBio);
    g.nodes[a.id].ends[ae].push_back(Link{b.id, static_cast<uint8_t>(be)});
    g.nodes[b.id].ends[be].push_back(Link{a.id, static_cast<uint8_t>(ae)});
}

std::string tmpDir;
bool pairsStore(const std::string& a, const std::string& b, int n, SequenceStore& reads) {
    const std::string f1n = tmpDir + "/r1.fq", f2n = tmpDir + "/r2.fq";
    {
        std::ofstream f1(f1n), f2(f2n);
        for (int i = 0; i < n; ++i) {
            // read 1 forward on a (upstream), read 2 reverse-complemented on b: an FR pair.
            const std::string r1 = a.substr(a.size() - 150 + 3 * i, 100);
            const std::string r2 = reverseComplement(b.substr(10 + 3 * i, 100));
            f1 << "@p" << i << "/1\n" << r1 << "\n+\n" << std::string(r1.size(), 'I') << "\n";
            f2 << "@p" << i << "/2\n" << r2 << "\n+\n" << std::string(r2.size(), 'I') << "\n";
        }
    }
    std::string err;
    Library lib;
    lib.r1 = f1n;
    lib.r2 = f2n;
    QualityTrim qt;
    qt.enabled = false;
    reads.setQualityTrim(qt);
    if (!reads.load({lib}, 1, err)) { std::printf("load: %s\n", err.c_str()); return false; }
    return true;
}

// Nominator configurations. Each sets exactly the variables it names.
struct Mode { const char* name; const char* fix; const char* umbrella; const char* oriented; };
const Mode kModes[] = {
    {"legacy (flags unset)", nullptr, nullptr, nullptr},
    {"FIX_GAP_NOMINATOR=1", "1", nullptr, nullptr},
    {"FIXES=1 (umbrella)", nullptr, "1", nullptr},
    {"FIXES=1, FIX_GAP_NOMINATOR=0", "0", "1", nullptr},
    {"GAP_ORIENTED=1", nullptr, nullptr, "1"},
};
constexpr int kNumModes = 5;
void setMode(const Mode& m) {
    unsetenv("TESSERACT_FIX_GAP_NOMINATOR");
    unsetenv("TESSERACT_FIXES");
    unsetenv("TESSERACT_GAP_ORIENTED");
    if (m.fix) setenv("TESSERACT_FIX_GAP_NOMINATOR", m.fix, 1);
    if (m.umbrella) setenv("TESSERACT_FIXES", m.umbrella, 1);
    if (m.oriented) setenv("TESSERACT_GAP_ORIENTED", m.oriented, 1);
}

int test1(int wrong[kNumModes], int miss[kNumModes]) {
    rng.seed(7);
    const int k = 33;
    Bio Y{rnd(300)}, Z{rnd(400)}, W{rnd(300)}, Q{rnd(300)}, U{rnd(300)}, S{rnd(300)}, V{rnd(400)};
    for (int mode = 0; mode < kNumModes; ++mode) wrong[mode] = miss[mode] = 0;
    for (int mask = 0; mask < 16; ++mask) {
        Z.flip = mask & 8; W.flip = mask & 4; S.flip = mask & 2; U.flip = mask & 1;
        UnitigGraph g;
        g.setK(k);
        uint32_t next = 0;
        for (Bio* b : {&Y, &Z, &W, &Q, &U, &S, &V}) {
            Unitig n;
            n.seq = b->flip ? reverseComplement(b->seq) : b->seq;
            n.coverage = 30;
            g.nodes.push_back(n);
            b->id = next++;
        }
        blink(g, Y, 0, Y, 1);
        blink(g, Y, 1, Z, 0);
        blink(g, Z, 1, W, 0);
        blink(g, Z, 1, Q, 0);
        blink(g, Q, 1, Y, 0);
        blink(g, W, 1, U, 0);
        blink(g, W, 1, S, 0);
        blink(g, S, 1, Y, 0);
        blink(g, V, 1, Y, 0);
        if (!g.validate().empty()) { std::printf("invalid fixture\n"); return 2; }
        const uint64_t uTip = endId(U.id, U.stored(1)), vTip = endId(V.id, V.stored(0));
        const auto want = std::make_pair(std::min(uTip, vTip), std::max(uTip, vTip));
        const Bio* srcs[] = {&Z, &W, &Q, &S, &U};
        for (int si = 0; si < 5; ++si) {
            SequenceStore reads;
            if (!pairsStore(srcs[si]->seq, V.seq, 5, reads)) return 2;
            for (int mode = 0; mode < kNumModes; ++mode) {
                setMode(kModes[mode]);
                size_t voted = 0;
                GapNominatorStats st;
                const auto nom = nominateGapJoins(g, reads, 5000, 31, 2, &voted, &st);
                const bool yes = nom.count(want) > 0;
                if (yes && srcs[si] != &U) ++wrong[mode];
                if (!yes && srcs[si] == &U) ++miss[mode];
            }
        }
    }
    return 0;
}

int test2(int right[kNumModes], int wrongEnd[kNumModes]) {
    rng.seed(9);
    const int k = 33;
    for (int mode = 0; mode < kNumModes; ++mode) right[mode] = wrongEnd[mode] = 0;
    const std::string xs = rnd(600), vs = rnd(400), ys = rnd(300);
    SequenceStore reads;
    if (!pairsStore(xs, vs, 5, reads)) return 2;
    for (int pad = 0; pad < 4; ++pad) {
        for (int flip = 0; flip < 2; ++flip) {
            UnitigGraph g;
            g.setK(k);
            Bio Y{ys}, X{xs}, V{vs};
            X.flip = flip;
            auto push = [&](Bio& b) {
                Unitig n;
                n.seq = b.flip ? reverseComplement(b.seq) : b.seq;
                n.coverage = 30;
                g.nodes.push_back(n);
                b.id = static_cast<uint32_t>(g.nodes.size() - 1);
            };
            push(Y);
            std::vector<Bio> P(pad);
            for (auto& p : P) { p.seq = rnd(200); push(p); }
            push(X);
            push(V);
            blink(g, Y, 0, Y, 1);
            for (auto& p : P) { blink(g, Y, 1, p, 0); blink(g, p, 1, Y, 0); }
            blink(g, V, 1, Y, 0);
            if (!g.validate().empty()) { std::printf("invalid fixture\n"); return 2; }
            const uint64_t vTip = endId(V.id, V.stored(0));
            const uint64_t xGood = endId(X.id, X.stored(1)), xBad = endId(X.id, X.stored(0));
            const auto good = std::make_pair(std::min(xGood, vTip), std::max(xGood, vTip));
            const auto bad = std::make_pair(std::min(xBad, vTip), std::max(xBad, vTip));
            for (int mode = 0; mode < kNumModes; ++mode) {
                setMode(kModes[mode]);
                size_t voted = 0;
                const auto nom = nominateGapJoins(g, reads, 5000, 31, 2, &voted);
                if (nom.count(good)) ++right[mode];
                if (nom.count(bad)) ++wrongEnd[mode];
            }
        }
    }
    return 0;
}

}  // namespace

int main() {
    clearTesseractEnv();
    char tmpl[] = "/tmp/tess_v3_gapnom.XXXXXX";
    const char* d = mkdtemp(tmpl);
    if (!d) return 2;
    tmpDir = d;

    int wrong[kNumModes], miss[kNumModes], right[kNumModes], wrongEnd[kNumModes];
    if (test1(wrong, miss) || test2(right, wrongEnd)) return 2;
    for (int m = 0; m < kNumModes; ++m) {
        std::printf("  %-30s TEST 1 wrong-source %2d/64 missed %2d/16 | TEST 2 correct %d/8 wrong-end %d/8\n",
                    kModes[m].name, wrong[m], miss[m], right[m], wrongEnd[m]);
    }
    check(wrong[0] == 36 && miss[0] == 0, "flags unset: legacy walk unchanged (release: 36/64 wrong-source)");
    check(right[0] == 4 && wrongEnd[0] == 4, "flags unset: legacy index unchanged (release: wrong end 4/8)");
    check(wrong[1] == 0 && miss[1] == 0, "FIX: no wrong-source nomination, no missed true pair");
    check(right[1] == 8 && wrongEnd[1] == 0, "FIX: isolated unitig credits the end the reads point to, 8/8");
    check(wrong[2] == wrong[1] && miss[2] == miss[1] && right[2] == right[1] && wrongEnd[2] == wrongEnd[1],
          "umbrella TESSERACT_FIXES=1 enables the fix");
    check(wrong[3] == wrong[0] && right[3] == right[0] && wrongEnd[3] == wrongEnd[0],
          "TESSERACT_FIX_GAP_NOMINATOR=0 overrides the umbrella");
    check(wrong[1] == wrong[4] && miss[1] == miss[4] && right[1] == right[4] && wrongEnd[1] == wrongEnd[4],
          "FIX agrees with the oriented nominator on both tests");

    std::printf("test_v3_graph_gap_nominator: %s (%d failure%s)\n", failures ? "FAILED" : "ALL PASSED",
                failures, failures == 1 ? "" : "s");
    const std::string cmd = "rm -rf '" + tmpDir + "'";
    if (std::system(cmd.c_str()) != 0) std::printf("note: could not remove %s\n", tmpDir.c_str());
    return failures ? 1 : 0;
}
