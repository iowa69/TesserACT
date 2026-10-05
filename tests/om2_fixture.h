// Shared fixtures for the Organism Model 2.0 C1 tests (header only; not a test of its own).
//
// Graphs are built the way the assembler builds them: every k-mer of a set of replicons is
// counted with a chosen depth, UnitigGraph::build + compact make the unitigs. Pieces (the
// contigs the gate is asked about) are unitig sequences oriented along the replicon, which is
// what the resolver emits at a repeat it cannot cross.
#pragma once

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <string>
#include <vector>

#include "counter.h"
#include "graph.h"
#include "kmer.h"
#include "om2_evidence.h"
#include "om2_ledger.h"
#include "test_env.h"

namespace om2fx {

using namespace ts;

inline int& failures() { static int f = 0; return f; }
inline void check(bool ok, const std::string& what) {
    std::printf("  %s  %s\n", ok ? "PASS" : "FAIL", what.c_str());
    if (!ok) ++failures();
}

struct Rng {
    std::mt19937 g;
    explicit Rng(unsigned seed) : g(seed) {}
    std::string seq(size_t n) {
        static const char* b = "ACGT";
        std::string s(n, 'A');
        for (char& c : s) c = b[g() & 3];
        return s;
    }
    // A copy of `s` with `n` substitutions at the given positions.
    static std::string mutateAt(std::string s, const std::vector<size_t>& at) {
        for (size_t p : at) s[p] = s[p] == 'A' ? 'C' : s[p] == 'C' ? 'G' : s[p] == 'G' ? 'T' : 'A';
        return s;
    }
};

struct Replicon {
    std::string seq;
    uint32_t depth;
    bool circular;
};

inline void addKmers(KmerTable& t, const std::string& s, int k, uint32_t depth, bool circular) {
    const std::string x = circular ? s + s.substr(0, static_cast<size_t>(k - 1)) : s;
    bool ok = true;
    for (size_t i = 0; i + static_cast<size_t>(k) <= x.size(); ++i) {
        const Kmer km = canonical(stringToKmer(x.substr(i, static_cast<size_t>(k)), k, ok), k);
        t.put(km, t.get(km) + depth);
    }
}

inline UnitigGraph buildGraph(const std::vector<Replicon>& reps, int k) {
    KmerTable t;
    for (const Replicon& r : reps) addKmers(t, r.seq, k, r.depth, r.circular);
    UnitigGraph g = UnitigGraph::build(t, k, 1);
    g.compact();
    g.removeDeleted();
    return g;
}

// The unitig holding `probe` (forward or reverse complemented), oriented so that `probe`
// reads forward in it. Empty when not exactly one unitig holds it.
inline std::string orientedNodeWith(const UnitigGraph& g, const std::string& probe, uint32_t* idOut = nullptr) {
    std::string found;
    int hits = 0;
    const std::string rp = reverseComplement(probe);
    for (uint32_t u = 0; u < g.nodes.size(); ++u) {
        if (g.nodes[u].deleted) continue;
        const std::string& s = g.nodes[u].seq;
        if (s.find(probe) != std::string::npos) { found = s; ++hits; if (idOut) *idOut = u; }
        else if (s.find(rp) != std::string::npos) { found = reverseComplement(s); ++hits; if (idOut) *idOut = u; }
    }
    return hits == 1 ? found : std::string();
}

inline std::vector<std::string> nodeSeqs(const UnitigGraph& g) {
    std::vector<std::string> v;
    for (const Unitig& u : g.nodes) if (!u.deleted) v.push_back(u.seq);
    return v;
}

// Gap between the end of `left` and the start of `right` in `genome` (both must occur once).
inline long trueGap(const std::string& genome, const std::string& left, const std::string& right) {
    const size_t a = genome.find(left);
    const size_t b = genome.find(right);
    if (a == std::string::npos || b == std::string::npos) return -999999;
    return static_cast<long>(b) - static_cast<long>(a + left.size());
}

inline om2::SeamConfig cfgMode(om2::SeamConfig::Mode m) {
    om2::SeamConfig c;
    c.mode = m;
    c.evidence = true;
    return c;
}

inline om2::Decision judgeJoin(om2::SeamContext& ctx, const std::string& left, const std::string& right,
                               int32_t claimedN, om2::Source src = om2::Source::Join) {
    om2::JudgeRequest rq;
    rq.stage = "test";
    rq.source = src;
    rq.left = &left;
    rq.right = &right;
    rq.claimedN = claimedN;
    rq.panelGap = claimedN;
    return ctx.judge(rq);
}


}  // namespace om2fx
