// T29 (build_v3, G-graph): simplify() convergence and TESSERACT_FIX_SIMPLIFY_COUNT_ALL.
//
// Ported from the T29 verification test (combo3/verify/T29). The stop rule counted tips,
// bubbles, erroneous connections and isolated unitigs, but not locally-weak removals, so
// a round whose only deletion came from removeLocallyWeak ended the loop and the tip that
// deletion exposed was never evaluated.
//
// Graph (k=21, readLength=100, meanCoverage=40):
//   X(40x,1000) --end1--> {Y(40x,1000), T(5x,50)}
//   T --end1--> W(2x,300); W end0 also receives S(40x,1000, 5' dead)
//   W --end1--> {Z(40x,1000), Z2(40x,1000)}
// Round 1 (ramp 0.40): nothing qualifies but W, which removeLocallyWeak deletes; T is left a
// tip beside Y. Release: the loop stops, T survives (6 live unitigs). Fixed: round 2 clips T
// and X+Y merge (4 live unitigs).
//
// Uses only release APIs and the environment, so the same source fails on release 1.3.0.
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <string>
#include <vector>

#include "graph.h"

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

void link(UnitigGraph& g, uint32_t u, int ue, uint32_t v, int ve) {
    g.nodes[u].ends[ue].push_back({v, static_cast<uint8_t>(ve)});
    g.nodes[v].ends[ve].push_back({u, static_cast<uint8_t>(ue)});
}

struct Outcome { size_t rounds = 0, live = 0; bool tipSurvives = false, valid = false; };

Outcome run() {
    std::mt19937 rng(29);
    auto rnd = [&](size_t n) {
        std::string s(n, 'A');
        for (auto& c : s) c = "ACGT"[rng() & 3];
        return s;
    };
    const int k = 21;
    const size_t K1 = k - 1;
    UnitigGraph g;
    g.setK(k);
    const std::string X = rnd(1000);
    const std::string sx = X.substr(X.size() - K1);
    const std::string Y = sx + rnd(980);
    const std::string T = sx + rnd(30);
    const std::string tt = T.substr(T.size() - K1);
    const std::string W = tt + rnd(280);
    const std::string ww = W.substr(W.size() - K1);
    const std::string Z = ww + "A" + rnd(979);
    const std::string Z2 = ww + "C" + rnd(979);
    const std::string S = rnd(980) + tt;
    auto add = [&](const std::string& s, double c) {
        Unitig u;
        u.seq = s;
        u.coverage = c;
        g.nodes.push_back(u);
        return static_cast<uint32_t>(g.nodes.size() - 1);
    };
    const uint32_t x = add(X, 40), y = add(Y, 40), t = add(T, 5), w = add(W, 2), z = add(Z, 40),
                   s = add(S, 40), z2 = add(Z2, 40);
    link(g, x, 1, y, 0);
    link(g, x, 1, t, 0);
    link(g, t, 1, w, 0);
    link(g, s, 1, w, 0);
    link(g, w, 1, z, 0);
    link(g, w, 1, z2, 0);
    Outcome o;
    if (!g.validate().empty()) return o;
    std::vector<SimplifyRoundStats> rounds;
    g.simplify(40.0, 100, false, 0.35, 12, &rounds, 0.0);
    o.valid = g.validate().empty();
    o.rounds = rounds.size();
    for (const Unitig& u : g.nodes) {
        if (u.deleted) continue;
        ++o.live;
        if (u.seq == T) o.tipSurvives = true;
    }
    return o;
}

void show(const char* tag, const Outcome& o) {
    std::printf("    [%s] rounds=%zu live=%zu tipT_survives=%d valid=%d\n", tag, o.rounds, o.live,
                o.tipSurvives ? 1 : 0, o.valid ? 1 : 0);
}

}  // namespace

int main() {
    clearTesseractEnv();
    const Outcome off = run();
    show("flags unset", off);
    check(off.valid && off.rounds == 1 && off.tipSurvives && off.live == 6,
          "flag off: legacy stop rule unchanged (1 round, T survives, 6 unitigs)");

    setenv("TESSERACT_FIX_SIMPLIFY_COUNT_ALL", "1", 1);
    const Outcome on = run();
    show("FIX_SIMPLIFY_COUNT_ALL=1", on);
    check(on.valid && on.rounds >= 2 && !on.tipSurvives && on.live == 4,
          "fix: the round after the locally-weak deletion runs and clips the exposed tip");

    unsetenv("TESSERACT_FIX_SIMPLIFY_COUNT_ALL");
    setenv("TESSERACT_FIXES", "1", 1);
    const Outcome um = run();
    show("FIXES=1", um);
    check(um.rounds == on.rounds && um.live == on.live && !um.tipSurvives, "umbrella enables the fix");
    setenv("TESSERACT_FIX_SIMPLIFY_COUNT_ALL", "0", 1);
    const Outcome ov = run();
    show("FIXES=1, FIX_SIMPLIFY_COUNT_ALL=0", ov);
    check(ov.rounds == off.rounds && ov.live == off.live && ov.tipSurvives, "own flag 0 overrides the umbrella");

    std::printf("test_v3_graph_simplify_stop: %s (%d failure%s)\n", failures ? "FAILED" : "ALL PASSED",
                failures, failures == 1 ? "" : "s");
    return failures ? 1 : 0;
}
