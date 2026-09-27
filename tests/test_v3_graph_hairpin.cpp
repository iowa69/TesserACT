// T28 (build_v3, G-graph): mergeInto and hairpin links; TESSERACT_FIX_HAIRPIN_KEEP.
//
// Ported from the T28 verification test (combo3/verify/T28). A unitig end whose terminal
// (k-1)-mer is a reverse-complement palindrome carries a hairpin link to itself. When
// compaction absorbs the node that carries it, the release drops the hairpin; when the
// carrier absorbs its neighbour instead, the hairpin survives. So node numbering decides
// the topology. The oracle is order invariance plus the graph built without the tip.
//
// Uses only release APIs and the environment, so the same source fails on release 1.3.0
// (8 failures) and passes with the fix; the flag-off half pins the release behaviour.
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <string>
#include <vector>

#include "graph.h"
#include "kmer.h"

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

size_t hairpins(const UnitigGraph& g) {
    size_t h = 0;
    for (uint32_t u = 0; u < g.nodes.size(); ++u) {
        if (g.nodes[u].deleted) continue;
        for (int e = 0; e < 2; ++e)
            for (const Link& l : g.nodes[u].ends[e])
                if (l.to == u && l.toEnd == e) ++h;
    }
    return h;
}
size_t live(const UnitigGraph& g) {
    size_t n = 0;
    for (const auto& x : g.nodes) n += !x.deleted;
    return n;
}
// Every live end whose terminal (k-1)-mer is a palindrome must carry its hairpin: build()
// links k-mer x to rc(x) there.
size_t palindromicEndsWithoutHairpin(const UnitigGraph& g) {
    const size_t K1 = static_cast<size_t>(g.k() - 1);
    size_t missing = 0;
    for (uint32_t u = 0; u < g.nodes.size(); ++u) {
        const Unitig& n = g.nodes[u];
        if (n.deleted || n.seq.size() < K1) continue;
        for (int e = 0; e < 2; ++e) {
            const std::string t = e == 1 ? n.seq.substr(n.seq.size() - K1) : n.seq.substr(0, K1);
            if (t != reverseComplement(t)) continue;
            bool has = false;
            for (const Link& l : n.ends[e]) has |= (l.to == u && l.toEnd == e);
            if (!has) ++missing;
        }
    }
    return missing;
}
void link(UnitigGraph& g, uint32_t u, int ue, uint32_t v, int ve) {   // addLink is private
    g.nodes[u].ends[ue].push_back(Link{v, static_cast<uint8_t>(ve)});
    if (!(u == v && ue == ve)) g.nodes[v].ends[ve].push_back(Link{u, static_cast<uint8_t>(ue)});
}
Unitig mk(const std::string& s) {
    Unitig u;
    u.seq = s;
    u.coverage = 30;
    return u;
}

// k=5, palindromic 4-mer ACGT. X=GGATTCA, Y=TTCAGACGT; hairpin Y.1 <-> Y.1.
// Returns (hairpins, dead ends) after compact.
std::pair<size_t, size_t> caseA(bool xFirst) {
    UnitigGraph g;
    g.setK(5);
    const std::string X = "GGATTCA", Y = "TTCAGACGT";
    uint32_t x, y;
    if (xFirst) { g.nodes.push_back(mk(X)); g.nodes.push_back(mk(Y)); x = 0; y = 1; }
    else        { g.nodes.push_back(mk(Y)); g.nodes.push_back(mk(X)); y = 0; x = 1; }
    link(g, x, 1, y, 0);
    link(g, y, 1, y, 1);
    const size_t merged = g.compact();
    const std::string s = g.nodes[0].seq;
    const bool seqOk = s == "GGATTCAGACGT" || s == reverseComplement(std::string("GGATTCAGACGT"));
    if (merged != 1 || live(g) != 1 || !seqOk || !g.validate().empty()) return {99, 99};
    return {hairpins(g), g.totalDeadEnds()};
}

// Y.1 -> {Y.1 hairpin, Z.0}; Z = ACGTAAC also carries its own hairpin at Z.0. Returns the
// number of links on the merged palindromic end and the total hairpins.
std::pair<size_t, size_t> caseB(bool xFirst) {
    UnitigGraph g;
    g.setK(5);
    const std::string X = "GGATTCA", Y = "TTCAGACGT", Z = "ACGTAAC";
    uint32_t x, y;
    const uint32_t z = 2;
    if (xFirst) { g.nodes.push_back(mk(X)); g.nodes.push_back(mk(Y)); x = 0; y = 1; }
    else        { g.nodes.push_back(mk(Y)); g.nodes.push_back(mk(X)); y = 0; x = 1; }
    g.nodes.push_back(mk(Z));
    link(g, x, 1, y, 0);
    link(g, y, 1, y, 1);
    link(g, y, 1, z, 0);
    link(g, z, 0, z, 0);
    g.compact();
    if (!g.validate().empty()) return {99, 99};
    uint32_t m = UINT32_MAX;
    for (uint32_t u = 0; u < g.nodes.size(); ++u)
        if (!g.nodes[u].deleted && g.nodes[u].seq != Z) m = u;
    if (m == UINT32_MAX) return {99, 99};
    const int pe = g.nodes[m].seq.substr(g.nodes[m].seq.size() - 4) == "ACGT" ? 1 : 0;
    return {g.nodes[m].ends[pe].size(), hairpins(g)};
}

std::mt19937 rng;
std::string rnd(size_t n) {
    std::string s(n, 'A');
    for (auto& c : s) c = "ACGT"[rng() & 3];
    return s;
}
void add(KmerTable& t, const std::string& s, int k, uint32_t c) {
    Kmer f = 0, r = 0;
    int v = 0;
    for (char ch : s) {
        const int b = baseCode(ch);
        f = pushBack(f, b, k);
        r = pushFrontRc(r, b, k);
        if (++v >= k) {
            const Kmer x = f < r ? f : r;
            t.put(x, t.get(x) + c);
        }
    }
}
// k-mer built: S = R(300) + P or P + R, P a 20-bp palindrome (k=21), a 2x tip off R.
// After removeTips + compact the graph must equal the one built from S alone.
struct CResult { int usable = 0, lost = 0, mismatch = 0, invalid = 0, unmarked = 0; };
CResult caseC(bool palAtStart, int tipPos) {
    const int k = 21;
    CResult r;
    for (unsigned seed = 1; seed <= 200; ++seed) {
        rng.seed(seed * 7919u + (palAtStart ? 1u : 0u) + static_cast<unsigned>(tipPos));
        const std::string h = rnd(10), P = h + reverseComplement(h), R = rnd(300);
        const std::string S = palAtStart ? P + R : R + P;
        const std::string tip = S.substr(static_cast<size_t>(tipPos), 21) + rnd(25);
        KmerTable t;
        add(t, S, k, 30);
        add(t, tip, k, 2);
        UnitigGraph g = UnitigGraph::build(t, k, 1);
        g.compact();
        KmerTable t2;
        add(t2, S, k, 30);
        UnitigGraph c = UnitigGraph::build(t2, k, 1);
        c.compact();
        if (hairpins(c) != 1 || live(c) != 1) continue;
        ++r.usable;
        g.removeTips(100, 0.5);
        g.compact();
        if (hairpins(g) != hairpins(c)) ++r.lost;
        if (hairpins(g) != hairpins(c) || g.totalDeadEnds() != c.totalDeadEnds() || live(g) != live(c))
            ++r.mismatch;
        if (!g.validate().empty()) ++r.invalid;
        if (palindromicEndsWithoutHairpin(g) != 0) ++r.unmarked;
    }
    return r;
}

struct Summary { std::pair<size_t, size_t> a0, a1, b0, b1; int cLost = 0, cMismatch = 0, cInvalid = 0, cUnmarked = 0, cUsable = 0; };
Summary runAll() {
    Summary s;
    s.a0 = caseA(true);
    s.a1 = caseA(false);
    s.b0 = caseB(true);
    s.b1 = caseB(false);
    const struct { bool start; int tip; } vars[] = {{false, 130}, {false, 60}, {true, 150}, {true, 250}};
    for (const auto& v : vars) {
        const CResult r = caseC(v.start, v.tip);
        s.cUsable += r.usable;
        s.cLost += r.lost;
        s.cMismatch += r.mismatch;
        s.cInvalid += r.invalid;
        s.cUnmarked += r.unmarked;
    }
    std::printf("    A: X first -> hairpins=%zu deadEnds=%zu; Y first -> hairpins=%zu deadEnds=%zu\n",
                s.a0.first, s.a0.second, s.a1.first, s.a1.second);
    std::printf("    B: X first -> palindrome-end links=%zu hairpins=%zu; Y first -> %zu / %zu\n",
                s.b0.first, s.b0.second, s.b1.first, s.b1.second);
    std::printf("    C: usable=%d hairpin_lost=%d topology_mismatch=%d palindromic_end_unmarked=%d validate_fail=%d (of 800 seeds)\n",
                s.cUsable, s.cLost, s.cMismatch, s.cUnmarked, s.cInvalid);
    return s;
}

}  // namespace

int main() {
    clearTesseractEnv();

    std::printf("flags unset (release behaviour):\n");
    const Summary off = runAll();
    check(off.a1.first == 1 && off.a1.second == 1, "flag off: the carrier absorbing keeps its hairpin (release)");
    check(off.a0.first == 0 && off.a0.second == 2, "flag off: the carrier being absorbed loses it (release defect, pinned)");
    check(off.cLost > 0 && off.cInvalid == 0, "flag off: case C loses hairpins, validate() cannot see it (release)");

    setenv("TESSERACT_FIX_HAIRPIN_KEEP", "1", 1);
    std::printf("TESSERACT_FIX_HAIRPIN_KEEP=1:\n");
    const Summary on = runAll();
    check(on.a0 == on.a1 && on.a0.first == 1 && on.a0.second == 1,
          "A: same topology whichever node absorbs (1 hairpin, 1 dead end)");
    check(on.b0 == on.b1 && on.b0.first == 2 && on.b0.second == 2,
          "B: merged palindrome end keeps both the hairpin and Z, in either order");
    check(on.cUsable == off.cUsable && on.cLost == 0 && on.cMismatch == 0 && on.cUnmarked == 0 &&
              on.cInvalid == 0,
          "C: 0 of 800 differ from the graph built without the tip");

    unsetenv("TESSERACT_FIX_HAIRPIN_KEEP");
    setenv("TESSERACT_FIXES", "1", 1);
    std::printf("TESSERACT_FIXES=1 (umbrella):\n");
    const Summary um = runAll();
    check(um.a0 == on.a0 && um.b0 == on.b0 && um.cLost == 0 && um.cMismatch == 0,
          "umbrella enables the fix");
    setenv("TESSERACT_FIX_HAIRPIN_KEEP", "0", 1);
    std::printf("TESSERACT_FIXES=1, TESSERACT_FIX_HAIRPIN_KEEP=0:\n");
    const Summary ov = runAll();
    check(ov.a0 == off.a0 && ov.cLost == off.cLost && ov.cMismatch == off.cMismatch,
          "own flag 0 overrides the umbrella");

    std::printf("test_v3_graph_hairpin: %s (%d failure%s)\n", failures ? "FAILED" : "ALL PASSED", failures,
                failures == 1 ? "" : "s");
    return failures ? 1 : 0;
}
