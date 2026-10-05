#include "om2_close.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <functional>
#include <map>
#include <mutex>
#include <queue>
#include <set>
#include <thread>
#include <unordered_map>
#include <unordered_set>

#include "envflags.h"
#include "graph.h"
#include "kmer.h"
#include "resolve.h"
#include "seqio.h"
#include "util.h"

namespace ts {
namespace om2 {

// ---- options ---------------------------------------------------------------------------------

CloseOptions CloseOptions::fromEnv() {
    CloseOptions o;
    static const char* const kNames[] = {"TESSERACT_OM2_CLOSE", "TESSERACT_OM2_FILL", "TESSERACT_OM2_ALLOC",
                                         "TESSERACT_OM2_CIRC",  "TESSERACT_OM2_RRN_PRIOR", "TESSERACT_OM2_DNAA",
                                         "TESSERACT_OM2_FLOW"};
    for (const char* n : kNames) if (env::isSet(n)) o.anyFlagSet = true;
    o.enabled = env::on("TESSERACT_OM2_CLOSE", false);
    if (const char* f = env::text("TESSERACT_OM2_FILL")) {
        o.fill = std::strcmp(f, "none") == 0 ? FillMode::None
               : std::strcmp(f, "contig") == 0 ? FillMode::Contig : FillMode::Genome;
    }
    if (const char* a = env::text("TESSERACT_OM2_ALLOC")) {
        o.alloc = std::strcmp(a, "off") == 0 ? AllocMode::Off
                : std::strcmp(a, "consensus") == 0 ? AllocMode::Consensus
                : std::strcmp(a, "prior") == 0 ? AllocMode::Prior : AllocMode::Phased;
    }
    o.circ = env::on("TESSERACT_OM2_CIRC", false);
    o.flow = env::on("TESSERACT_OM2_FLOW", false);
    if (const char* p = env::text("TESSERACT_OM2_RRN_PRIOR")) o.rrnPrior = p;
    if (const char* p = env::text("TESSERACT_OM2_DNAA")) o.dnaa = p;
    o.allocParams.mode = o.alloc;
    return o;
}

std::string formatCloseCounters(const CloseStats& s) {
    char buf[1600];
    std::snprintf(buf, sizeof buf,
                  "[om2-close] enabled=%d clusters=%zu bridges=%zu forced=%zu pair_thread=%zu tiebreak=%zu open=%zu "
                  "budget_refused=%zu hairpin_refused=%zu fills_contig=%zu fills_genome=%zu alloc_pair=%zu "
                  "alloc_thread=%zu alloc_prior=%zu alloc_multiplicity=%zu alloc_consensus=%zu variant_sites=%zu "
                  "circles_tested=%zu circles_closed=%zu rotated=%zu abstain_tangled=%zu | junctions=%zu wraps=%zu "
                  "anchored=%zu unanchored=%zu tiebreak_unverified=%zu cyclic=%zu double_use=%zu no_walk=%zu "
                  "flank_mismatch=%zu overlap=%zu length_inconsistent=%zu fill_bp_contig=%zu fill_bp_genome=%zu "
                  "prior_noflip=%zu budget_exhausted=%zu circles_contig=%zu circles_overlap=%zu c1_attached=%zu "
                  "c1_refused=%zu c1_gated=%zu wraps_open=%zu variant_rows=%zu prior_loaded=%d dnaa_loaded=%d "
                  "seconds=%.1f",
                  s.enabled ? 1 : 0, s.clusters, s.bridges, s.forced, s.pairThread, s.tiebreak, s.open,
                  s.budgetRefused, s.hairpinRefused, s.fillsContig, s.fillsGenome, s.alloc.pair, s.alloc.thread,
                  s.alloc.prior, s.alloc.multiplicity, s.alloc.consensus, s.variantSites, s.circlesTested,
                  s.circlesClosed, s.rotated, s.abstainTangled, s.junctions, s.wraps, s.anchored, s.unanchored,
                  s.tiebreakUnverified, s.cyclic, s.doubleUse, s.noWalk, s.flankMismatch, s.overlap,
                  s.lengthInconsistent, s.fillBpContig, s.fillBpGenome, s.alloc.priorNoFlip, s.alloc.budgetExhausted,
                  s.circlesContig, s.circlesOverlap, s.c1Attached, s.c1Refused, s.c1Gated, s.wrapsOpen,
                  s.variantRows, s.priorLoaded ? 1 : 0, s.dnaaLoaded ? 1 : 0, s.seconds);
    return buf;
}

void logCloseCounters(const CloseStats& s) { std::fprintf(stderr, "%s\n", formatCloseCounters(s).c_str()); }

namespace {

constexpr int K31 = 31;
constexpr uint64_t kMask31 = (1ULL << 62) - 1;
constexpr uint32_t kTarget = 0xFFFFFFFEu;   // the target instance in a route search
constexpr int64_t kInf = INT64_MAX / 4;

template <typename F>
void forEach31(const std::string& s, size_t from, size_t to, F&& fn) {
    uint64_t fwd = 0, rev = 0;
    int valid = 0;
    to = std::min(to, s.size());
    for (size_t i = from; i < to; ++i) {
        const int c = baseCode(s[i]);
        if (c < 0) { valid = 0; fwd = rev = 0; continue; }
        fwd = ((fwd << 2) | static_cast<uint64_t>(c)) & kMask31;
        rev = (rev >> 2) | (static_cast<uint64_t>(3 - c) << 60);
        if (++valid < K31) continue;
        fn(fwd <= rev ? fwd : rev, i + 1 - static_cast<size_t>(K31), fwd <= rev);
    }
}

inline uint32_t onode(uint32_t n, int o) { return (n << 1) | static_cast<uint32_t>(o & 1); }

// ---- the graph, with copy numbers -------------------------------------------------------------

struct GraphView {
    const UnitigGraph& g;
    int k;
    int64_t kk;   // k - 1: the overlap between linked unitigs
    std::vector<double> theta;
    std::vector<float> copy;
    std::vector<uint8_t> single, hardSingle;
    std::vector<int> hi;

    GraphView(const UnitigGraph& gr, const CloseOptions& opt) : g(gr), k(gr.k()), kk(gr.k() - 1) {
        const size_t n = g.nodes.size();
        theta.assign(n, 0);
        copy.assign(n, 0);
        single.assign(n, 0);
        hardSingle.assign(n, 0);
        hi.assign(n, 1);
        const std::vector<std::vector<uint32_t>> comps = g.components();
        for (const auto& c : comps) {
            const double th = g.weightedMedianCoverage(c, 10);
            for (uint32_t v : c) {
                if (v >= n) continue;
                const double cov = g.nodes[v].coverage;
                const double cv = th > 0 ? cov / th : 1.0;
                theta[v] = th;
                copy[v] = static_cast<float>(cv);
                single[v] = cv <= opt.singleCopy ? 1 : 0;
                hardSingle[v] = (g.nodes[v].seq.size() >= opt.hardSingleLen && cov <= opt.hardSingleRatio * th) ? 1 : 0;
                const long r = std::lround(cv);
                hi[v] = std::max(1, static_cast<int>(cv > 4 ? r + 2 : r + 1));
            }
        }
    }
    bool live(uint32_t node) const { return node < g.nodes.size() && !g.nodes[node].deleted; }
    int64_t len(uint32_t o) const { return static_cast<int64_t>(g.nodes[o >> 1].seq.size()); }
    int64_t lenp(uint32_t o) const { return len(o) - kk; }
    std::string oseq(uint32_t o) const { return g.oriented(o >> 1, static_cast<int>(o & 1)); }
    void succ(uint32_t o, std::vector<uint32_t>& out) const {
        out.clear();
        for (const Link& l : g.exits(o >> 1, static_cast<int>(o & 1))) {
            if (!live(l.to)) continue;
            out.push_back(onode(l.to, UnitigGraph::enterOrient(l)));
        }
    }
    void pred(uint32_t o, std::vector<uint32_t>& out) const {
        succ(o ^ 1u, out);
        for (uint32_t& x : out) x ^= 1u;
    }
};

// ---- 31-mer census over the graph --------------------------------------------------------------

struct Occ {
    uint32_t count = 0;
    uint32_t node = UINT32_MAX;
    uint32_t pos = 0;
    uint8_t fwd = 1;          // the node's forward 31-mer is the canonical one
    uint8_t admissible = 0;   // not wholly inside an overlap shared with a neighbour
};

struct Census {
    std::unordered_map<uint64_t, uint32_t> idx;
    std::vector<Occ> occ;
    std::vector<uint64_t> key;

    uint32_t add(uint64_t canon) {
        const auto it = idx.find(canon);
        if (it != idx.end()) return it->second;
        const uint32_t id = static_cast<uint32_t>(occ.size());
        idx.emplace(canon, id);
        occ.emplace_back();
        key.push_back(canon);
        return id;
    }
    const Occ* find(uint64_t canon) const {
        const auto it = idx.find(canon);
        return it == idx.end() ? nullptr : &occ[it->second];
    }
    // Counts every occurrence of every registered 31-mer in the graph.
    void run(const GraphView& gv) {
        if (occ.empty()) return;
        const size_t ov = static_cast<size_t>(gv.kk);
        for (uint32_t v = 0; v < gv.g.nodes.size(); ++v) {
            const Unitig& u = gv.g.nodes[v];
            if (u.deleted) continue;
            const bool hasPred = !u.ends[0].empty(), hasSucc = !u.ends[1].empty();
            const size_t L = u.seq.size();
            forEach31(u.seq, 0, L, [&](uint64_t canon, size_t p, bool fwd) {
                const auto it = idx.find(canon);
                if (it == idx.end()) return;
                Occ& o = occ[it->second];
                if (++o.count == 1) {
                    o.node = v;
                    o.pos = static_cast<uint32_t>(p);
                    o.fwd = fwd ? 1 : 0;
                    const bool inHead = hasPred && p + K31 <= ov;
                    const bool inTail = hasSucc && L >= ov && p >= L - ov;
                    o.admissible = (inHead || inTail) ? 0 : 1;
                }
            });
        }
    }
    bool unique(const Occ& o) const { return o.count == 1 && o.admissible; }
};

// ---- routes ------------------------------------------------------------------------------------

struct Site {
    std::vector<std::vector<uint32_t>> branches;   // interior oriented nodes of each path (may be empty)
    bool complex = false;
    std::string key;                               // orientation-free identity of the bubble
    std::vector<std::string> branchKey;
};

struct Route {
    bool found = false, zero = false, cyclic = false, hairpin = false, abstain = false, exhaustive = true;
    bool usesRepeat = false, usesHardSingle = false;
    std::vector<uint32_t> spine;     // cut nodes, spine[0] = source, spine.back() = kTarget
    std::vector<int> siteAt;         // per spine gap: site index or -1
    std::vector<Site> sites;
    int64_t fmin = 0, fmax = 0;      // fill length range before the flank trims
    double walks = 0;
    std::vector<int64_t> classes;    // walk-length classes (+-50 bp)
    std::vector<uint16_t> classN;
};

std::string nodeSetKey(std::vector<uint32_t> nodes) {
    std::sort(nodes.begin(), nodes.end());
    nodes.erase(std::unique(nodes.begin(), nodes.end()), nodes.end());
    std::string s;
    for (uint32_t n : nodes) { s += std::to_string(n); s += ','; }
    return s.empty() ? std::string("-") : s;
}

struct RouteQuery {
    uint32_t oa = 0;
    int64_t posEndA = 0;
    uint32_t ob = 0;
    int64_t posStartB = 0;
    int64_t maxFill = 12000;
};

// All walks from anchor A (end of its 31-mer, in oriented node oa) to anchor B (start of its
// 31-mer, in ob) whose fill is at most maxFill. `blocked` nodes may not be used as interior
// nodes; neither may either anchor's unitig.
Route searchRoute(const GraphView& gv, const RouteQuery& q, const std::vector<uint8_t>& blocked,
                  const CloseOptions& opt) {
    Route r;
    const int64_t kk = gv.kk;
    if (q.oa == q.ob && q.posStartB >= q.posEndA) {
        r.found = r.zero = true;
        r.fmin = r.fmax = q.posStartB - q.posEndA;
        r.walks = 1;
        r.classes = {r.fmin};
        r.classN = {1};
        return r;
    }
    const uint32_t na = q.oa >> 1, nb = q.ob >> 1;
    auto interiorOk = [&](uint32_t o) {
        const uint32_t n = o >> 1;
        return n != na && n != nb && !(n < blocked.size() && blocked[n]);
    };
    std::vector<uint32_t> nb_;
    // forward
    std::unordered_map<uint32_t, int64_t> D;
    int64_t DT = kInf;
    using QE = std::pair<int64_t, uint32_t>;
    std::priority_queue<QE, std::vector<QE>, std::greater<QE>> pq;
    const int64_t d0 = gv.len(q.oa) - q.posEndA;
    const int64_t bound = q.maxFill + kk;
    auto relaxF = [&](uint32_t s, int64_t d) {
        if (d > bound) return;
        auto it = D.find(s);
        if (it == D.end() || d < it->second) { D[s] = d; pq.emplace(d, s); }
    };
    gv.succ(q.oa, nb_);
    for (uint32_t s : nb_) {
        if (s == q.ob) DT = std::min(DT, d0);
        else if (interiorOk(s)) relaxF(s, d0);
    }
    int pops = 0;
    while (!pq.empty()) {
        const QE e = pq.top();
        pq.pop();
        if (D[e.second] != e.first) continue;
        if (++pops > opt.maxStates) { r.abstain = true; break; }
        const int64_t nd = e.first + gv.lenp(e.second);
        gv.succ(e.second, nb_);
        for (uint32_t s : nb_) {
            if (s == q.ob) DT = std::min(DT, nd);
            else if (interiorOk(s)) relaxF(s, nd);
        }
    }
    if (r.abstain) return r;
    const int64_t ET = q.posStartB - kk;
    if (DT >= kInf || DT + ET > q.maxFill) return r;   // no walk within the cap
    // backward
    std::unordered_map<uint32_t, int64_t> E;
    std::priority_queue<QE, std::vector<QE>, std::greater<QE>> pb;
    auto relaxB = [&](uint32_t p, int64_t e) {
        if (e > q.maxFill) return;
        auto it = E.find(p);
        if (it == E.end() || e < it->second) { E[p] = e; pb.emplace(e, p); }
    };
    gv.pred(q.ob, nb_);
    for (uint32_t p : nb_) if (p != q.oa && interiorOk(p)) relaxB(p, gv.lenp(p) + ET);
    pops = 0;
    while (!pb.empty()) {
        const QE e = pb.top();
        pb.pop();
        if (E[e.second] != e.first) continue;
        if (++pops > opt.maxStates) { r.abstain = true; break; }
        gv.pred(e.second, nb_);
        for (uint32_t p : nb_) if (p != q.oa && interiorOk(p)) relaxB(p, gv.lenp(p) + e.first);
    }
    if (r.abstain) return r;
    // nodes on some walk within the cap
    std::vector<uint32_t> S;
    for (const auto& kv : D) {
        const auto it = E.find(kv.first);
        if (it != E.end() && kv.second + it->second <= q.maxFill) S.push_back(kv.first);
    }
    std::sort(S.begin(), S.end());
    std::unordered_map<uint32_t, uint32_t> id;   // oriented node -> local index (0 source, 1 target)
    id[q.oa] = 0;
    for (size_t i = 0; i < S.size(); ++i) id[S[i]] = static_cast<uint32_t>(i + 2);
    const size_t V = S.size() + 2;
    std::vector<uint32_t> nodeOf(V);
    nodeOf[0] = q.oa;
    nodeOf[1] = kTarget;
    for (size_t i = 0; i < S.size(); ++i) nodeOf[i + 2] = S[i];
    for (uint32_t o : S) if (std::binary_search(S.begin(), S.end(), o ^ 1u)) r.hairpin = true;
    std::vector<std::vector<uint32_t>> adj(V);
    std::vector<int> indeg(V, 0);
    auto addEdges = [&](uint32_t from, uint32_t o, int64_t dFrom) {
        gv.succ(o, nb_);
        const int64_t base = dFrom + (from == 0 ? 0 : gv.lenp(o));
        for (uint32_t s : nb_) {
            if (s == q.ob) {
                if (base + ET <= q.maxFill) { adj[from].push_back(1); ++indeg[1]; }
                continue;
            }
            const auto it = id.find(s);
            if (it == id.end() || it->second < 2) continue;
            if (base + E[s] <= q.maxFill) { adj[from].push_back(it->second); ++indeg[it->second]; }
        }
    };
    addEdges(0, q.oa, d0);
    for (size_t i = 0; i < S.size(); ++i) addEdges(static_cast<uint32_t>(i + 2), S[i], D[S[i]]);
    // Kahn, deterministic
    std::vector<int> deg = indeg;
    std::priority_queue<std::pair<int64_t, uint32_t>, std::vector<std::pair<int64_t, uint32_t>>,
                        std::greater<std::pair<int64_t, uint32_t>>> ready;
    auto dOf = [&](uint32_t v) -> int64_t { return v == 0 ? -1 : v == 1 ? kInf : D[nodeOf[v]]; };
    ready.emplace(dOf(0), 0);
    std::vector<uint32_t> order;
    while (!ready.empty()) {
        const uint32_t v = ready.top().second;
        ready.pop();
        order.push_back(v);
        for (uint32_t w : adj[v]) if (--deg[w] == 0) ready.emplace(dOf(w), w);
    }
    r.found = true;
    for (uint32_t o : S) {
        if (!gv.single[o >> 1]) r.usesRepeat = true;
        if (gv.hardSingle[o >> 1]) r.usesHardSingle = true;
    }
    if (order.size() != V || order.back() != 1) {
        r.cyclic = true;
        r.exhaustive = false;
        r.fmin = DT + ET;
        r.fmax = q.maxFill;
        r.walks = 2;   // at least: a cycle within the cap
        return r;
    }
    // cut nodes: every edge leaving the processed set enters v
    std::vector<int> outdeg(V, 0);
    for (size_t v = 0; v < V; ++v) outdeg[v] = static_cast<int>(adj[v].size());
    std::vector<size_t> posInOrder(V);
    for (size_t i = 0; i < order.size(); ++i) posInOrder[order[i]] = i;
    std::vector<size_t> cutPos;
    long open = 0;
    for (size_t i = 0; i < order.size(); ++i) {
        const uint32_t v = order[i];
        if (open == indeg[v]) cutPos.push_back(i);
        open -= indeg[v];
        open += outdeg[v];
    }
    // spine and sites
    int64_t fixedLen = d0 + ET;   // contributions of the source tail and target head
    std::vector<std::pair<int64_t, int64_t>> siteRange;   // per site: min / max branch contribution
    std::vector<std::vector<int64_t>> siteLens;
    for (size_t c = 0; c < cutPos.size(); ++c) {
        const uint32_t v = order[cutPos[c]];
        r.spine.push_back(nodeOf[v]);
        if (v >= 2) fixedLen += gv.lenp(nodeOf[v]);
        if (c + 1 == cutPos.size()) break;
        const size_t lo = cutPos[c], hiPos = cutPos[c + 1];
        const uint32_t w = order[hiPos];
        if (hiPos == lo + 1) {   // adjacent in order: only the direct edge
            r.siteAt.push_back(-1);
            continue;
        }
        // enumerate the paths v -> w through the nodes between them
        Site site;
        std::vector<uint32_t> stack;
        bool overflow = false;
        std::function<void(uint32_t)> dfs = [&](uint32_t x) {
            if (overflow) return;
            for (uint32_t y : adj[x]) {
                if (y == w) {
                    if (site.branches.size() >= static_cast<size_t>(opt.maxBranches)) { overflow = true; return; }
                    std::vector<uint32_t> b;
                    for (uint32_t s : stack) b.push_back(nodeOf[s]);
                    site.branches.push_back(b);
                    continue;
                }
                if (posInOrder[y] <= lo || posInOrder[y] >= hiPos) continue;
                stack.push_back(y);
                dfs(y);
                stack.pop_back();
                if (overflow) return;
            }
        };
        dfs(v);
        if (overflow) {
            // too many paths: one greedy path by depth, consensus only
            site.branches.clear();
            site.complex = true;
            r.exhaustive = false;
            std::vector<uint32_t> b;
            uint32_t x = v;
            int guard = 0;
            while (x != w && ++guard < 100000) {
                uint32_t best = UINT32_MAX;
                double bestCov = -1;
                for (uint32_t y : adj[x]) {
                    if (y == w) { best = w; break; }
                    if (posInOrder[y] <= lo || posInOrder[y] >= hiPos) continue;
                    const double cv = gv.g.nodes[nodeOf[y] >> 1].coverage;
                    if (cv > bestCov) { bestCov = cv; best = y; }
                }
                if (best == UINT32_MAX) break;
                if (best != w) b.push_back(nodeOf[best]);
                x = best;
            }
            site.branches.push_back(b);
        }
        std::vector<uint32_t> all;
        all.push_back(nodeOf[v] == kTarget ? q.ob >> 1 : nodeOf[v] >> 1);
        all.push_back(nodeOf[w] == kTarget ? q.ob >> 1 : nodeOf[w] >> 1);
        std::vector<int64_t> lens;
        for (const auto& b : site.branches) {
            std::vector<uint32_t> ns;
            int64_t L = 0;
            for (uint32_t o : b) { ns.push_back(o >> 1); all.push_back(o >> 1); L += gv.lenp(o); }
            site.branchKey.push_back(nodeSetKey(ns));
            lens.push_back(L);
        }
        site.key = nodeSetKey(all);
        siteLens.push_back(lens);
        r.siteAt.push_back(static_cast<int>(r.sites.size()));
        r.sites.push_back(site);
    }
    // walk count and length classes
    r.walks = 1;
    std::vector<int64_t> sums = {fixedLen};
    for (const auto& lens : siteLens) {
        r.walks *= static_cast<double>(lens.size());
        std::vector<int64_t> next;
        for (int64_t s : sums) for (int64_t l : lens) next.push_back(s + l);
        std::sort(next.begin(), next.end());
        next.erase(std::unique(next.begin(), next.end()), next.end());
        if (next.size() > 256) next = {next.front(), next.back()};
        sums.swap(next);
    }
    r.fmin = sums.front();
    r.fmax = sums.back();
    for (int64_t s : sums) {
        if (!r.classes.empty() && s - r.classes.back() <= 50) { ++r.classN.back(); continue; }
        r.classes.push_back(s);
        r.classN.push_back(1);
    }
    return r;
}

// The single-copy unitigs first reached from an oriented node through repeat unitigs only.
std::set<uint32_t> firstUnique(const GraphView& gv, uint32_t from, bool forward, int64_t maxLen, int maxStates,
                               bool& capped) {
    std::set<uint32_t> U;
    capped = false;
    std::unordered_map<uint32_t, int64_t> seen;
    std::queue<std::pair<uint32_t, int64_t>> q;
    std::vector<uint32_t> nb;
    if (forward) gv.succ(from, nb); else gv.pred(from, nb);
    for (uint32_t s : nb) q.emplace(s, 0);
    int states = 0;
    while (!q.empty()) {
        const auto e = q.front();
        q.pop();
        if (seen.count(e.first)) continue;
        seen[e.first] = e.second;
        if (++states > maxStates) { capped = true; break; }
        if (gv.single[e.first >> 1]) { U.insert(e.first >> 1); continue; }
        const int64_t nd = e.second + gv.lenp(e.first);
        if (nd > maxLen) { capped = true; continue; }
        if (forward) gv.succ(e.first, nb); else gv.pred(e.first, nb);
        for (uint32_t s : nb) q.emplace(s, nd);
    }
    return U;
}

// ---- junction working state ---------------------------------------------------------------------

// Gated (integration): no fragments span and the graph does not force it, but C1's gate judged the
// adjacency graph-backed (PASS_EXACT, PASS_WALK, RESIZE, CONTRA_PAIRS_OK). Ranked with the model
// tie-break it stands in for; its fills stay genome-only unless forced or spanned, as before.
enum class Bridge : uint8_t { None, Forced, PairThread, Tiebreak, Gated, TiebreakUnverified, Open };

const char* bridgeName(Bridge b) {
    switch (b) {
        case Bridge::None: return "none";
        case Bridge::Forced: return "forced";
        case Bridge::PairThread: return "pair_thread";
        case Bridge::Tiebreak: return "tiebreak";
        case Bridge::Gated: return "c1_gated";
        case Bridge::TiebreakUnverified: return "tiebreak_unverified";
        case Bridge::Open: return "open";
    }
    return "?";
}

const char* verdictStr(Verdict v) {
    switch (v) {
        case Verdict::NotJudged: return "NOT_JUDGED";
        case Verdict::PassExact: return "PASS_EXACT";
        case Verdict::PassWalk: return "PASS_WALK";
        case Verdict::Resize: return "RESIZE";
        case Verdict::BreakDoubleUse: return "BREAK_DOUBLE_USE";
        case Verdict::ContraPairsOk: return "CONTRA_PAIRS_OK";
        case Verdict::ContraRefused: return "CONTRA_REFUSED";
        case Verdict::PairsDeny: return "PAIRS_DENY";
        case Verdict::Silent: return "SILENT";
        case Verdict::Unanchored: return "UNANCHORED";
        case Verdict::Abstain: return "ABSTAIN";
    }
    return "?";
}

const char* endStr(EndClass e) {
    switch (e) {
        case EndClass::Unknown: return "unknown";
        case EndClass::Unique: return "unique";
        case EndClass::Repeat: return "repeat";
        case EndClass::DeadEnd: return "dead_end";
    }
    return "?";
}

struct CJ {
    uint32_t ledgerIdx = 0;
    const Junction* c1 = nullptr;      // C1's record of this junction, when C1 made one
    int64_t c1Id = -1;                 // its row in in.c1 (integration), -1 none
    char c1Orient = '.';               // '+': the record reads C1's row as recorded; '-': reversed
    std::vector<RepeatVariant> vrows;  // repeat_variants rows of this locus, in record orientation
    bool wrap = false;
    uint32_t rec = 0;
    size_t nStart = 0, nEnd = 0;       // the N-run [nStart, nEnd); a wrap has nStart = nEnd = record length
    int32_t claimedN = 0;
    // anchors
    bool anchoredA = false, anchoredB = false;
    uint32_t oa = 0, ob = 0;
    int64_t posEndA = 0, posStartB = 0, dA = 0, dB = 0;
    size_t recAEnd = 0, recBStart = 0;
    EndClass endA = EndClass::Unknown, endB = EndClass::Unknown;
    float copyA = 0, copyB = 0;
    std::vector<uint64_t> tagsL, tagsR;   // canonical 31-mers
    // route
    Route route;
    Verdict verdict = Verdict::NotJudged;
    bool forced = false;
    size_t uA = 0, uB = 0;
    // evidence
    uint32_t spanThread = 0, spanPair = 0;
    // decisions
    Bridge bridge = Bridge::None;
    bool accepted = false, contig = false, genome = false;
    std::string reason;
    std::vector<int> choice;               // per route site: canonical branch index
    std::vector<Basis> choiceBasis;
    std::vector<int> localBranch;          // per route site: local branch index chosen
    std::string fill;
    std::vector<FillSpan> spans;
    int64_t overlapBp = 0;
    int priorLocus = -1;
};

bool sized(int32_t n) { return n > 1 && n != 100 && n != 2000; }

// ---- link counting over the reads ----------------------------------------------------------------

struct LinkCount { uint32_t thread = 0, pair = 0; };

inline uint64_t pairKey(uint32_t a, uint32_t b) {
    if (a > b) std::swap(a, b);
    return (static_cast<uint64_t>(a) << 32) | b;
}

class TagIndex {
public:
    void add(uint64_t canon, uint32_t group) {
        auto& v = map_[canon];
        if (std::find(v.begin(), v.end(), group) == v.end() && v.size() < 4) v.push_back(group);
        const uint64_t h = mix64(canon);
        bits_[(h >> 38) >> 6] |= 1ULL << ((h >> 38) & 63);
    }
    TagIndex() : bits_((1ULL << 26) / 64, 0) {}
    bool empty() const { return map_.empty(); }
    const std::vector<uint32_t>* find(uint64_t canon) const {
        const uint64_t h = mix64(canon);
        if (!((bits_[(h >> 38) >> 6] >> ((h >> 38) & 63)) & 1ULL)) return nullptr;
        const auto it = map_.find(canon);
        return it == map_.end() ? nullptr : &it->second;
    }

private:
    std::unordered_map<uint64_t, std::vector<uint32_t>> map_;
    std::vector<uint64_t> bits_;
};

void readGroups(const SequenceStore& rs, size_t r, const TagIndex& tags, std::vector<uint32_t>& out) {
    out.clear();
    const uint32_t len = rs.length(r);
    if (len < static_cast<uint32_t>(K31)) return;
    uint64_t fwd = 0, rev = 0;
    int valid = 0;
    for (uint32_t p = 0; p < len; ++p) {
        const int c = rs.baseAt(r, p);
        if (c < 0) { valid = 0; fwd = rev = 0; continue; }
        fwd = ((fwd << 2) | static_cast<uint64_t>(c)) & kMask31;
        rev = (rev >> 2) | (static_cast<uint64_t>(3 - c) << 60);
        if (++valid < K31) continue;
        const std::vector<uint32_t>* g = tags.find(fwd <= rev ? fwd : rev);
        if (g) out.insert(out.end(), g->begin(), g->end());
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    if (out.size() > 16) out.resize(16);
}

std::unordered_map<uint64_t, LinkCount> countLinks(const SequenceStore& rs, const TagIndex& tags, int threads) {
    std::unordered_map<uint64_t, LinkCount> total;
    if (tags.empty() || rs.size() == 0) return total;
    const size_t pairs = rs.pairCount();
    const size_t singles = rs.size() - rs.pairedReads();
    const size_t units = pairs + singles;
    const int T = std::max(1, std::min(threads, 64));
    std::mutex mu;
    auto worker = [&](int t) {
        std::unordered_map<uint64_t, LinkCount> local;
        std::vector<uint32_t> g1, g2;
        const size_t from = units * static_cast<size_t>(t) / static_cast<size_t>(T);
        const size_t to = units * static_cast<size_t>(t + 1) / static_cast<size_t>(T);
        for (size_t u = from; u < to; ++u) {
            if (u < pairs) {
                readGroups(rs, 2 * u, tags, g1);
                readGroups(rs, 2 * u + 1, tags, g2);
            } else {
                readGroups(rs, rs.pairedReads() + (u - pairs), tags, g1);
                g2.clear();
            }
            if (g1.size() + g2.size() < 2) continue;
            std::unordered_set<uint64_t> doneKeys;
            for (int side = 0; side < 2; ++side) {
                const std::vector<uint32_t>& g = side == 0 ? g1 : g2;
                for (size_t i = 0; i < g.size(); ++i) {
                    for (size_t j = i + 1; j < g.size(); ++j) {
                        const uint64_t key = pairKey(g[i], g[j]);
                        if (doneKeys.insert(key).second) ++local[key].thread;
                    }
                }
            }
            for (uint32_t a : g1) {
                for (uint32_t b : g2) {
                    if (a == b) continue;
                    const uint64_t key = pairKey(a, b);
                    if (doneKeys.insert(key).second) ++local[key].pair;
                }
            }
        }
        std::lock_guard<std::mutex> lk(mu);
        for (const auto& kv : local) {
            LinkCount& c = total[kv.first];
            c.thread += kv.second.thread;
            c.pair += kv.second.pair;
        }
    };
    std::vector<std::thread> pool;
    for (int t = 0; t < T; ++t) pool.emplace_back(worker, t);
    for (auto& th : pool) th.join();
    return total;
}

// ---- small helpers -------------------------------------------------------------------------------

std::vector<std::pair<size_t, size_t>> nRuns(const std::string& s) {
    std::vector<std::pair<size_t, size_t>> out;
    size_t i = 0;
    while (i < s.size()) {
        if (s[i] == 'N' || s[i] == 'n') {
            size_t j = i;
            while (j < s.size() && (s[j] == 'N' || s[j] == 'n')) ++j;
            out.emplace_back(i, j);
            i = j;
        } else {
            ++i;
        }
    }
    return out;
}

// [from, to) of the flank window left of `pos` (stops at an N) and right of `pos`.
std::pair<size_t, size_t> leftWindow(const std::string& s, size_t pos, size_t w) {
    size_t from = pos > w ? pos - w : 0;
    for (size_t i = pos; i > from; --i) {
        const char c = s[i - 1];
        if (c == 'N' || c == 'n') { from = i; break; }
    }
    return {from, pos};
}
std::pair<size_t, size_t> rightWindow(const std::string& s, size_t pos, size_t w) {
    size_t to = std::min(s.size(), pos + w);
    for (size_t i = pos; i < to; ++i) {
        const char c = s[i];
        if (c == 'N' || c == 'n') { to = i; break; }
    }
    return {pos, to};
}

// Orientation and oriented position of a record 31-mer found in a unitig.
void placeAnchor(const GraphView& gv, const Occ& o, bool recFwdCanon, uint32_t& onodeOut, int64_t& posOut) {
    const bool sameStrand = (recFwdCanon ? 1 : 0) == o.fwd;
    const int64_t L = static_cast<int64_t>(gv.g.nodes[o.node].seq.size());
    onodeOut = onode(o.node, sameStrand ? 0 : 1);
    posOut = sameStrand ? static_cast<int64_t>(o.pos) : L - static_cast<int64_t>(o.pos) - K31;
}

std::string fmtAnchor(const GraphView& gv, bool ok, uint32_t o, int64_t pos) {
    if (!ok) return "-";
    char b[64];
    std::snprintf(b, sizeof b, "%u%c:%lld", o >> 1, (o & 1) ? '-' : '+', static_cast<long long>(pos));
    (void)gv;
    return b;
}

// Largest-remainder integer split of `total` in proportion to `w`.
std::vector<int> largestRemainder(const std::vector<double>& w, int total) {
    std::vector<int> out(w.size(), 0);
    double sum = 0;
    for (double x : w) sum += std::max(0.0, x);
    if (sum <= 0 || total <= 0) return out;
    std::vector<std::pair<double, size_t>> rem;
    int given = 0;
    for (size_t i = 0; i < w.size(); ++i) {
        const double exact = std::max(0.0, w[i]) / sum * total;
        out[i] = static_cast<int>(std::floor(exact));
        given += out[i];
        rem.emplace_back(-(exact - out[i]), i);
    }
    std::sort(rem.begin(), rem.end());
    for (size_t i = 0; given < total && i < rem.size(); ++i, ++given) ++out[rem[i].second];
    return out;
}

}  // namespace

// ---- the stage ------------------------------------------------------------------------------------

CloseStats runCloseStage(const CloseInputs& in, const CloseOptions& opt, std::vector<std::string>& seqs,
                         Ledger& ledger) {
    CloseStats st;
    if (!opt.enabled) {
        if (opt.anyFlagSet) {
            std::fprintf(stderr, "[om2-close] note: a TESSERACT_OM2_* option of this stage is set but "
                                 "TESSERACT_OM2_CLOSE is not; the stage is off\n");
        }
        return st;
    }
    st.enabled = true;
    const auto t0 = std::chrono::steady_clock::now();
    const UnitigGraph& g = in.graph;
    if (g.k() < K31 || seqs.empty()) {
        st.seconds = 0;
        return st;
    }
    GraphView gv(g, opt);
    const int64_t kk = gv.kk;

    // sidecars
    RrnPrior prior;
    DnaaSketch dnaa;
    std::string modelMd5;
    if ((!opt.rrnPrior.empty() || !opt.dnaa.empty()) && !in.modelPath.empty()) modelMd5 = md5File(in.modelPath);
    if (!opt.rrnPrior.empty()) {
        std::string err;
        st.priorLoaded = prior.load(opt.rrnPrior, modelMd5, err);
        if (!st.priorLoaded) std::fprintf(stderr, "[om2-close] rRNA prior refused: %s\n", err.c_str());
    }
    if (!opt.dnaa.empty()) {
        std::string err;
        st.dnaaLoaded = dnaa.load(opt.dnaa, modelMd5, err);
        if (!st.dnaaLoaded) std::fprintf(stderr, "[om2-close] dnaA sketch refused: %s\n", err.c_str());
    }
    AllocParams ap = opt.allocParams;
    if (ap.mode == AllocMode::Prior && !st.priorLoaded) ap.mode = AllocMode::Phased;

    // tag window: how far a read pair reaches from a flank
    int64_t reach = 600;
    if (in.insert.usable && in.insert.mean > 0) {
        reach = static_cast<int64_t>(in.insert.mean + 3.0 * in.insert.stddev);
        reach = std::max<int64_t>(300, std::min<int64_t>(1500, reach));
    }

    // ---- 1. junctions: every N-run, and the wrap of each record long enough ----
    std::vector<CJ> J;
    // New rows' ids continue after C1's rows (integration) and after any rows already in `ledger`.
    const size_t c1Rows = in.c1 ? in.c1->j.size() : 0;
    size_t newRows = 0;
    auto nextId = [&]() { return static_cast<uint32_t>(c1Rows + ledger.j.size() + newRows++); };
    for (size_t r = 0; r < seqs.size(); ++r) {
        for (const auto& nr : nRuns(seqs[r])) {
            if (nr.first == 0 || nr.second >= seqs[r].size()) continue;
            CJ j;
            j.rec = static_cast<uint32_t>(r);
            j.nStart = nr.first;
            j.nEnd = nr.second;
            j.claimedN = static_cast<int32_t>(nr.second - nr.first);
            if (in.c1 && in.c1Runs && r < in.c1Runs->size()) {
                for (const RunOwner& ro : (*in.c1Runs)[r]) {
                    if (ro.start == static_cast<int64_t>(nr.first) && ro.end == static_cast<int64_t>(nr.second) &&
                        ro.id >= 0 && static_cast<size_t>(ro.id) < c1Rows) {
                        j.c1Id = ro.id;
                        j.c1Orient = ro.orient;
                        j.c1 = &in.c1->j[static_cast<size_t>(ro.id)];
                    }
                }
            }
            if (j.c1Id >= 0) { j.ledgerIdx = static_cast<uint32_t>(j.c1Id); ++st.c1Attached; }
            else j.ledgerIdx = nextId();
            J.push_back(j);
        }
    }
    st.junctions = J.size();
    if (opt.circ) {
        std::vector<std::pair<size_t, size_t>> bySize;
        for (size_t r = 0; r < seqs.size(); ++r) {
            if (seqs[r].size() >= static_cast<size_t>(opt.circMinLen)) bySize.emplace_back(seqs[r].size(), r);
        }
        std::sort(bySize.begin(), bySize.end(), [](const std::pair<size_t, size_t>& a,
                                                   const std::pair<size_t, size_t>& b) {
            return a.first != b.first ? a.first > b.first : a.second < b.second;
        });
        for (size_t i = 0; i < bySize.size() && i < 64; ++i) {
            CJ j;
            j.ledgerIdx = nextId();
            j.wrap = true;
            j.rec = static_cast<uint32_t>(bySize[i].second);
            j.nStart = j.nEnd = seqs[j.rec].size();
            j.claimedN = 0;
            J.push_back(j);
            ++st.wraps;
        }
    }

    // w-th wrap junction, in J order (= ledger order of the wrap rows): the rotateOffset index
    std::vector<int> wrapOrd(J.size(), -1);
    {
        int w = 0;
        for (size_t i = 0; i < J.size(); ++i) if (J[i].wrap) wrapOrd[i] = w++;
    }

    auto leftWin = [&](const CJ& j, size_t w) {
        return leftWindow(seqs[j.rec], j.nStart, w);
    };
    auto rightWin = [&](const CJ& j, size_t w) {
        return j.wrap ? rightWindow(seqs[j.rec], 0, w) : rightWindow(seqs[j.rec], j.nEnd, w);
    };

    // ---- 2. census 1: flank 31-mers ----
    Census c1;
    const size_t aw = static_cast<size_t>(opt.anchorWindow);
    for (const CJ& j : J) {
        const auto lw = leftWin(j, aw), rw = rightWin(j, aw);
        forEach31(seqs[j.rec], lw.first, lw.second, [&](uint64_t c, size_t, bool) { c1.add(c); });
        forEach31(seqs[j.rec], rw.first, rw.second, [&](uint64_t c, size_t, bool) { c1.add(c); });
    }
    c1.run(gv);

    // anchors and flank tags
    for (CJ& j : J) {
        const std::string& s = seqs[j.rec];
        const auto lw = leftWin(j, aw), rw = rightWin(j, aw);
        const size_t leftTagFrom = j.nStart > static_cast<size_t>(reach) ? j.nStart - static_cast<size_t>(reach) : 0;
        const size_t rightTagEnd = (j.wrap ? 0 : j.nEnd) + static_cast<size_t>(reach);
        bool sawRepeatL = false, sawRepeatR = false;
        forEach31(s, lw.first, lw.second, [&](uint64_t c, size_t p, bool fwd) {
            const Occ* o = c1.find(c);
            if (!o || !c1.unique(*o)) return;
            if (!gv.single[o->node]) { sawRepeatL = true; return; }
            if (p >= leftTagFrom) j.tagsL.push_back(c);
            j.anchoredA = true;   // rightmost wins: later positions overwrite
            placeAnchor(gv, *o, fwd, j.oa, j.posEndA);
            j.posEndA += K31;
            j.recAEnd = p + K31;
            j.copyA = gv.copy[o->node];
        });
        bool firstB = true;
        forEach31(s, rw.first, rw.second, [&](uint64_t c, size_t p, bool fwd) {
            const Occ* o = c1.find(c);
            if (!o || !c1.unique(*o)) return;
            if (!gv.single[o->node]) { sawRepeatR = true; return; }
            if (p + K31 <= rightTagEnd) j.tagsR.push_back(c);
            if (!firstB) return;
            firstB = false;
            j.anchoredB = true;
            placeAnchor(gv, *o, fwd, j.ob, j.posStartB);
            j.recBStart = p;
            j.copyB = gv.copy[o->node];
        });
        j.endA = j.anchoredA ? EndClass::Unique : sawRepeatL ? EndClass::Repeat : EndClass::Unknown;
        j.endB = j.anchoredB ? EndClass::Unique : sawRepeatR ? EndClass::Repeat : EndClass::Unknown;
        if (j.anchoredA) j.dA = static_cast<int64_t>(j.nStart) - static_cast<int64_t>(j.recAEnd);
        if (j.anchoredB) j.dB = static_cast<int64_t>(j.recBStart) - static_cast<int64_t>(j.wrap ? 0 : j.nEnd);
        if (j.anchoredA && j.anchoredB) ++st.anchored; else ++st.unanchored;
    }

    // ---- 3. placed single-copy unitigs (may not be walked through again) ----
    std::vector<uint8_t> placed(g.nodes.size(), 0);
    {
        Census c2;
        std::unordered_map<uint64_t, uint32_t> probeNode;
        for (uint32_t v = 0; v < g.nodes.size(); ++v) {
            if (!gv.live(v) || !gv.hardSingle[v]) continue;
            const std::string& sq = g.nodes[v].seq;
            const size_t L = sq.size();
            for (size_t p : {static_cast<size_t>(kk) + 1, L / 2, L > static_cast<size_t>(kk) + 2 * K31
                                                                ? L - static_cast<size_t>(kk) - K31 - 1 : L / 2}) {
                if (p + K31 > L) continue;
                forEach31(sq, p, p + K31, [&](uint64_t c, size_t, bool) { probeNode[c] = v; });
            }
        }
        for (const std::string& s : seqs) {
            forEach31(s, 0, s.size(), [&](uint64_t c, size_t, bool) {
                const auto it = probeNode.find(c);
                if (it != probeNode.end()) placed[it->second] = 1;
            });
        }
    }

    // ---- 4. routes ----
    for (CJ& j : J) {
        if (!j.anchoredA || !j.anchoredB) { j.verdict = Verdict::Unanchored; continue; }
        RouteQuery q;
        q.oa = j.oa;
        q.posEndA = j.posEndA;
        q.ob = j.ob;
        q.posStartB = j.posStartB;
        q.maxFill = opt.maxFill;
        j.route = searchRoute(gv, q, placed, opt);
        if (j.route.abstain) { ++st.abstainTangled; j.verdict = Verdict::Abstain; j.reason = "search_cap"; continue; }
        if (!j.route.found) {
            const std::vector<uint8_t> none;
            const Route any = searchRoute(gv, q, none, opt);
            if (any.found) { j.verdict = Verdict::BreakDoubleUse; ++st.doubleUse; j.reason = "needs_placed_single_copy"; }
            else { j.verdict = Verdict::Silent; ++st.noWalk; j.reason = "no_walk"; }
            continue;
        }
        const int64_t fAdj = j.route.fmin - j.dA - j.dB;
        const bool unique = !j.route.cyclic && j.route.sites.empty();
        if (!j.wrap && sized(j.claimedN)) {
            const int64_t slack = std::max<int64_t>(1000, j.claimedN / 4);
            if (fAdj > j.claimedN + slack) j.verdict = Verdict::Resize;
            else if (unique && std::llabs(fAdj - j.claimedN) <= 1000) j.verdict = Verdict::PassExact;
            else j.verdict = Verdict::PassWalk;
        } else {
            j.verdict = unique ? Verdict::PassExact : Verdict::PassWalk;
        }
        bool cA = false, cB = false;
        const std::set<uint32_t> UA = firstUnique(gv, j.oa, true, opt.maxFill, opt.maxStates, cA);
        const std::set<uint32_t> UB = firstUnique(gv, j.ob, false, opt.maxFill, opt.maxStates, cB);
        j.uA = UA.size();
        j.uB = UB.size();
        j.forced = j.route.zero || (!cA && !cB && UA.size() == 1 && UB.size() == 1 && *UA.begin() == (j.ob >> 1) &&
                                    *UB.begin() == (j.oa >> 1));
    }

    // ---- 5. branch tags, then one pass over the reads ----
    // groups: 2i = left flank of junction i, 2i+1 = right flank; then one per (site, branch)
    std::map<std::string, uint32_t> siteId;             // site key -> index
    std::vector<std::vector<std::string>> siteBranchKeys;   // canonical (sorted) branch keys
    std::vector<std::vector<std::vector<uint32_t>>> siteBranchNodes;   // one oriented instance per branch
    for (const CJ& j : J) {
        if (!j.route.found || j.route.cyclic) continue;
        for (const Site& s : j.route.sites) {
            auto it = siteId.find(s.key);
            if (it == siteId.end()) {
                it = siteId.emplace(s.key, static_cast<uint32_t>(siteBranchKeys.size())).first;
                std::vector<std::pair<std::string, size_t>> order;
                for (size_t b = 0; b < s.branchKey.size(); ++b) order.emplace_back(s.branchKey[b], b);
                std::sort(order.begin(), order.end());
                std::vector<std::string> keys;
                std::vector<std::vector<uint32_t>> nodes;
                for (const auto& ob : order) { keys.push_back(ob.first); nodes.push_back(s.branches[ob.second]); }
                siteBranchKeys.push_back(keys);
                siteBranchNodes.push_back(nodes);
            }
        }
    }
    std::vector<uint32_t> branchGroupBase(siteBranchKeys.size(), 0);
    uint32_t nextGroup = static_cast<uint32_t>(2 * J.size());
    for (size_t s = 0; s < siteBranchKeys.size(); ++s) {
        branchGroupBase[s] = nextGroup;
        nextGroup += static_cast<uint32_t>(siteBranchKeys[s].size());
    }
    TagIndex tags;
    for (size_t i = 0; i < J.size(); ++i) {
        for (uint64_t c : J[i].tagsL) tags.add(c, static_cast<uint32_t>(2 * i));
        for (uint64_t c : J[i].tagsR) tags.add(c, static_cast<uint32_t>(2 * i + 1));
    }
    {
        Census c3;
        for (const auto& nodesOfSite : siteBranchNodes) {
            for (const auto& br : nodesOfSite) {
                for (uint32_t o : br) {
                    const std::string& sq = g.nodes[o >> 1].seq;
                    forEach31(sq, 0, sq.size(), [&](uint64_t c, size_t, bool) { c3.add(c); });
                }
            }
        }
        c3.run(gv);
        for (size_t s = 0; s < siteBranchNodes.size(); ++s) {
            for (size_t b = 0; b < siteBranchNodes[s].size(); ++b) {
                for (uint32_t o : siteBranchNodes[s][b]) {
                    const std::string& sq = g.nodes[o >> 1].seq;
                    forEach31(sq, 0, sq.size(), [&](uint64_t c, size_t, bool) {
                        const Occ* oc = c3.find(c);
                        if (oc && c3.unique(*oc) && oc->node == (o >> 1)) {
                            tags.add(c, branchGroupBase[s] + static_cast<uint32_t>(b));
                        }
                    });
                }
            }
        }
    }
    const std::unordered_map<uint64_t, LinkCount> links = countLinks(in.reads, tags, in.threads);
    auto link = [&](uint32_t a, uint32_t b) -> LinkCount {
        const auto it = links.find(pairKey(a, b));
        return it == links.end() ? LinkCount() : it->second;
    };
    for (size_t i = 0; i < J.size(); ++i) {
        const LinkCount lc = link(static_cast<uint32_t>(2 * i), static_cast<uint32_t>(2 * i + 1));
        J[i].spanThread = lc.thread;
        J[i].spanPair = lc.pair;
    }

    // ---- 6. bridges, in order forced > pairs/threads > tie-break > open, under the budget ----
    std::vector<size_t> order;
    for (size_t i = 0; i < J.size(); ++i) {
        CJ& j = J[i];
        if (!j.route.found || j.route.abstain) continue;
        if (j.route.hairpin) { ++st.hairpinRefused; j.bridge = Bridge::Open; j.reason = "hairpin"; continue; }
        if (j.route.cyclic) { ++st.cyclic; j.bridge = Bridge::Open; j.reason = "cyclic"; continue; }
        ++st.bridges;
        // C1's judgement of this junction, when C1 recorded it (the integrated build); standalone
        // scans have none, and a model-proposed join then counts as an unverified tie-break.
        const Junction* led = j.c1;
        // Integration: C1 judged this adjacency. A junction C1 would break (a skipped single-copy
        // node, a contradiction the pairs do not confirm, a pair deny) or does not admit is never
        // bridged here, whatever the walk: a fill would only give a misjoin bases.
        const bool c1Refuses = led && (led->admit == Admit::Break || led->writtenN < 0 ||
                                       led->verdict == Verdict::BreakDoubleUse ||
                                       led->verdict == Verdict::ContraRefused || led->verdict == Verdict::PairsDeny);
        const bool c1Graph = led && (led->verdict == Verdict::PassExact || led->verdict == Verdict::PassWalk ||
                                     led->verdict == Verdict::Resize || led->verdict == Verdict::ContraPairsOk);
        if (c1Refuses) { j.bridge = Bridge::Open; ++st.open; ++st.c1Refused; j.reason = std::string("c1_") + verdictName(led->verdict); continue; }
        // C4 (om2_clonal): a non-positional junction the nearest relatives do not agree on stays a
        // labelled gap; no walk may write bases into it. Only the clonal stage sets noFill.
        if (led && led->noFill) { j.bridge = Bridge::Open; ++st.open; j.reason = "clonal_nofill"; continue; }
        if (j.forced) j.bridge = Bridge::Forced;
        else if (j.spanThread + j.spanPair >= opt.minSpan) j.bridge = Bridge::PairThread;
        else if (led && led->panelGenomes >= opt.tieGenomes && led->tieMargin >= opt.tieMargin) j.bridge = Bridge::Tiebreak;
        else if (c1Graph) j.bridge = Bridge::Gated;
        else if (!led || led->panelGenomes == 0) j.bridge = Bridge::TiebreakUnverified;
        else j.bridge = Bridge::Open;
        switch (j.bridge) {
            case Bridge::Forced: ++st.forced; break;
            case Bridge::PairThread: ++st.pairThread; break;
            case Bridge::Tiebreak: ++st.tiebreak; break;
            case Bridge::Gated: ++st.c1Gated; break;
            case Bridge::TiebreakUnverified: ++st.tiebreakUnverified; break;
            default: ++st.open; j.reason = "tiebreak_margin"; break;
        }
        if (j.bridge != Bridge::Open) order.push_back(i);
    }
    std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) {
        return static_cast<int>(J[a].bridge) < static_cast<int>(J[b].bridge);
    });
    std::vector<int> used(g.nodes.size(), 0);
    for (size_t i : order) {
        CJ& j = J[i];
        if (opt.fill == FillMode::None) { j.reason = "fill_none"; continue; }
        bool ok = true;
        const std::vector<uint32_t>& sp = j.route.spine;
        for (size_t s = 1; s + 1 < sp.size(); ++s) {
            const uint32_t n = sp[s] >> 1;
            if (gv.hardSingle[n] ? used[n] >= 1 : used[n] >= gv.hi[n]) { ok = false; break; }
        }
        if (!ok) { ++st.budgetRefused; j.reason = "budget"; continue; }
        if (!j.route.sites.empty() && ap.mode == AllocMode::Off) { j.reason = "alloc_off"; continue; }
        const int64_t fAdjMin = j.route.fmin - j.dA - j.dB;
        if (!j.wrap && sized(j.claimedN) && j.verdict == Verdict::Resize && j.route.usesHardSingle) {
            ++st.lengthInconsistent;
            j.reason = "length_inconsistent";
            continue;
        }
        if (j.wrap && (!j.route.sites.empty() || !j.route.exhaustive)) { j.reason = "wrap_not_unique"; continue; }
        // A unique wrap walk through a repeat unitig cannot tell a circle from a segment the same
        // repeat flanks on both sides (IS-X-IS; MEASURED on S. aureus GCF000010465v1: a 2.5 kb
        // chromosomal record closed through its flanking 965 bp repeat). The adjacency must be
        // forced by the graph, spanned by fragments, or the walk repeat-free.
        if (j.wrap && !(j.forced || j.spanThread + j.spanPair >= opt.minSpan || !j.route.usesRepeat)) {
            j.reason = "wrap_through_repeat_unverified";
            continue;
        }
        // Even the longest walk far shorter than a sized gap: a shortcut through a repeat, not the gap.
        if (!j.wrap && sized(j.claimedN) &&
            j.route.fmax - j.dA - j.dB < j.claimedN - std::max<int64_t>(1000, j.claimedN / 4)) {
            ++st.lengthInconsistent;
            j.reason = "walk_shorter_than_gap";
            continue;
        }
        (void)fAdjMin;
        for (size_t s = 1; s + 1 < sp.size(); ++s) ++used[sp[s] >> 1];
        j.accepted = true;
    }

    // repeat clusters touched by accepted bridges
    {
        std::vector<uint32_t> parent(g.nodes.size());
        for (uint32_t v = 0; v < parent.size(); ++v) parent[v] = v;
        std::function<uint32_t(uint32_t)> findp = [&](uint32_t x) {
            while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; }
            return x;
        };
        std::set<uint32_t> rep;
        for (const CJ& j : J) {
            if (!j.accepted) continue;
            auto addRep = [&](uint32_t o) { if (o != kTarget && !gv.single[o >> 1]) rep.insert(o >> 1); };
            for (uint32_t o : j.route.spine) addRep(o);
            for (const Site& s : j.route.sites) for (const auto& b : s.branches) for (uint32_t o : b) addRep(o);
        }
        std::vector<uint32_t> nb;
        for (uint32_t v : rep) {
            for (int o = 0; o < 2; ++o) {
                gv.succ(onode(v, o), nb);
                for (uint32_t w : nb) if (rep.count(w >> 1)) parent[findp(v)] = findp(w >> 1);
            }
        }
        std::set<uint32_t> roots;
        for (uint32_t v : rep) roots.insert(findp(v));
        st.clusters = roots.size();
    }

    // ---- 7. allocation, per site, over every accepted locus that walks it ----
    // junction flank markers for the locus prior
    if (st.priorLoaded) {
        const uint64_t thr = ~0ULL / prior.flankDenom();
        for (CJ& j : J) {
            if (!j.accepted || j.route.sites.empty()) continue;
            std::vector<uint64_t> markers;
            const size_t fb = static_cast<size_t>(std::max(500, prior.flankBp()));
            const auto lw = leftWin(j, fb), rw = rightWin(j, fb);
            auto take = [&](uint64_t c, size_t, bool) { const uint64_t h = sidecarHash(c); if (h <= thr) markers.push_back(h); };
            forEach31(seqs[j.rec], lw.first, lw.second, take);
            forEach31(seqs[j.rec], rw.first, rw.second, take);
            j.priorLocus = prior.matchLocus(markers);
        }
    }
    struct LocusRef { size_t junction; size_t site; std::vector<int> localToCanon; };
    std::vector<std::vector<LocusRef>> siteLoci(siteBranchKeys.size());
    for (size_t i = 0; i < J.size(); ++i) {
        CJ& j = J[i];
        if (!j.accepted) continue;
        j.choice.assign(j.route.sites.size(), -1);
        j.choiceBasis.assign(j.route.sites.size(), Basis::Consensus);
        j.localBranch.assign(j.route.sites.size(), -1);
        for (size_t s = 0; s < j.route.sites.size(); ++s) {
            const Site& site = j.route.sites[s];
            const uint32_t sid = siteId[site.key];
            LocusRef ref;
            ref.junction = i;
            ref.site = s;
            for (const std::string& bk : site.branchKey) {
                const auto& keys = siteBranchKeys[sid];
                ref.localToCanon.push_back(static_cast<int>(std::lower_bound(keys.begin(), keys.end(), bk) - keys.begin()));
            }
            siteLoci[sid].push_back(ref);
        }
    }
    struct SiteReport {
        std::string key; bool its = false, complex = false; std::vector<int64_t> blen; std::vector<double> bcov;
        std::vector<int> carriers, flow, bclass; int majority = -1; std::vector<std::string> allele;
        std::vector<std::string> assigned;
    };
    std::vector<SiteReport> reports(siteBranchKeys.size());
    for (size_t sid = 0; sid < siteBranchKeys.size(); ++sid) {
        if (siteLoci[sid].empty()) continue;
        const size_t m = siteBranchKeys[sid].size();
        SiteInput in2;
        in2.branches.resize(m);
        SiteReport& rep = reports[sid];
        rep.blen.resize(m);
        rep.bcov.resize(m);
        rep.allele.resize(m);
        // cut-node depth around the site, from the first locus that walks it
        const LocusRef& first = siteLoci[sid].front();
        const Route& fr = J[first.junction].route;
        size_t gap = 0;
        for (size_t x = 0; x < fr.siteAt.size(); ++x) if (fr.siteAt[x] == static_cast<int>(first.site)) gap = x;
        auto nodeCov = [&](uint32_t o) {
            const uint32_t n = o == kTarget ? (J[first.junction].ob >> 1) : (o >> 1);
            return g.nodes[n].coverage;
        };
        const uint32_t leftCut = fr.spine[gap], rightCut = fr.spine[gap + 1];
        const double covCut = std::min(nodeCov(leftCut), nodeCov(rightCut));
        const double th = gv.theta[leftCut == kTarget ? (J[first.junction].ob >> 1) : (leftCut >> 1)];
        const Site& fsite = fr.sites[first.site];
        const bool complexSite = fsite.complex;
        // branch b (canonical order) as the first locus walks it
        std::vector<std::vector<uint32_t>> bnodes(m);
        for (size_t lb = 0; lb < first.localToCanon.size(); ++lb) {
            const int cb = first.localToCanon[lb];
            if (cb >= 0 && static_cast<size_t>(cb) < m) bnodes[static_cast<size_t>(cb)] = fsite.branches[lb];
        }
        int emptyIdx = -1;
        int64_t maxVar = 0;
        for (size_t b = 0; b < m; ++b) {
            const auto& nodes = bnodes[b];
            double cmin = -1;
            int64_t L = 0;
            std::string bseq;
            for (uint32_t o : nodes) {
                const double c = g.nodes[o >> 1].coverage;
                cmin = cmin < 0 ? c : std::min(cmin, c);
                L += gv.lenp(o);
                const std::string os = gv.oseq(o);
                bseq += bseq.empty() ? os : os.substr(static_cast<size_t>(kk));
            }
            rep.blen[b] = L;
            maxVar = std::max<int64_t>(maxVar, L - kk);
            if (nodes.empty()) { emptyIdx = static_cast<int>(b); rep.allele[b] = "-"; continue; }
            in2.branches[b].cov = cmin;
            in2.branches[b].carriers = std::max(1, static_cast<int>(std::lround(th > 0 ? cmin / th : 1.0)));
            rep.bcov[b] = cmin;
            if (bseq.size() <= 60) {
                rep.allele[b] = bseq;
            } else {
                char h[64];
                std::snprintf(h, sizeof h, "len=%zu;h=%016llx", bseq.size(),
                              static_cast<unsigned long long>(mix64(std::hash<std::string>()(bseq))));
                rep.allele[b] = h;
            }
        }
        if (emptyIdx >= 0) {
            int others = 0;
            for (size_t b = 0; b < m; ++b) if (static_cast<int>(b) != emptyIdx) others += in2.branches[b].carriers;
            const int total = static_cast<int>(std::lround(th > 0 ? covCut / th : 1.0));
            in2.branches[static_cast<size_t>(emptyIdx)].carriers = std::max(1, total - others);
            in2.branches[static_cast<size_t>(emptyIdx)].cov = std::max(0.0, covCut);
            rep.bcov[static_cast<size_t>(emptyIdx)] = std::max(0.0, covCut);
        }
        int64_t minLen = INT64_MAX, maxLen = 0;
        for (size_t b = 0; b < m; ++b) { minLen = std::min(minLen, rep.blen[b]); maxLen = std::max(maxLen, rep.blen[b]); }
        (void)maxVar;
        in2.itsClass = m >= 2 && maxLen - minLen >= 50;   // ITS types differ in length; SNV clusters do not
        in2.complex = complexSite;
        rep.its = in2.itsClass;
        rep.complex = complexSite;
        // branch classes for the prior
        std::vector<int> bclass(m, -1);
        bool classesOk = false;
        if (st.priorLoaded && in2.itsClass && !complexSite) {
            for (size_t b = 0; b < m; ++b) {
                std::string ctx;
                const std::string ls = leftCut == kTarget ? gv.oseq(J[first.junction].ob) : gv.oseq(leftCut);
                ctx = ls.substr(ls.size() > 300 ? ls.size() - 300 : 0);
                for (uint32_t o : bnodes[b]) ctx += gv.oseq(o).substr(static_cast<size_t>(kk));
                const std::string rs = rightCut == kTarget ? gv.oseq(J[first.junction].ob) : gv.oseq(rightCut);
                ctx += rs.substr(static_cast<size_t>(kk), 300);
                bclass[b] = prior.classOf(ctx);
                if (bclass[b] >= 0) classesOk = true;
            }
            in2.branchClass = bclass;
        }
        rep.bclass = bclass;
        for (const LocusRef& ref : siteLoci[sid]) {
            SiteLocus L;
            L.junction = static_cast<uint32_t>(ref.junction);
            L.threadN.assign(m, 0);
            L.pairN.assign(m, 0);
            for (size_t b = 0; b < m; ++b) {
                const uint32_t gB = branchGroupBase[sid] + static_cast<uint32_t>(b);
                for (uint32_t side = 0; side < 2; ++side) {
                    const LinkCount lc = link(static_cast<uint32_t>(2 * ref.junction + side), gB);
                    L.threadN[b] += lc.thread;
                    L.pairN[b] += lc.pair;
                }
            }
            const int locus = J[ref.junction].priorLocus;
            if (classesOk && locus >= 0) {
                const uint32_t n = prior.locusGenomes(locus);
                const size_t C = prior.classes();
                L.priorGenomes = n;
                L.prior.resize(C);
                for (size_t c = 0; c < C; ++c) {
                    L.prior[c] = (prior.classCount(locus, static_cast<int>(c)) + 0.5) / (n + 0.5 * static_cast<double>(C));
                }
            }
            in2.loci.push_back(L);
        }
        const std::vector<SiteChoice> ch = allocateSite(in2, ap, st.alloc);
        ++st.variantSites;
        for (size_t li = 0; li < siteLoci[sid].size(); ++li) {
            const LocusRef& ref = siteLoci[sid][li];
            CJ& j = J[ref.junction];
            j.choice[ref.site] = ch[li].branch;
            j.choiceBasis[ref.site] = ch[li].basis;
            for (size_t lb = 0; lb < ref.localToCanon.size(); ++lb) {
                if (ref.localToCanon[lb] == ch[li].branch) j.localBranch[ref.site] = static_cast<int>(lb);
            }
            char a[96];
            std::snprintf(a, sizeof a, "%u:%d:%s%s", j.ledgerIdx, ch[li].branch, c2BasisName(ch[li].basis),
                          ch[li].priorFlip ? "(flip)" : "");
            rep.assigned.push_back(a);
        }
        rep.majority = majorityBranch(in2);
        for (size_t b = 0; b < m; ++b) rep.carriers.push_back(in2.branches[b].carriers);
        if (opt.flow && !complexSite) {
            std::vector<double> w;
            for (size_t b = 0; b < m; ++b) w.push_back(rep.bcov[b]);
            rep.flow = largestRemainder(w, static_cast<int>(siteLoci[sid].size()));
        }
    }

    // ---- 8. fills ----
    for (CJ& j : J) {
        if (!j.accepted) continue;
        const Route& r = j.route;
        bool allPhased = true;
        for (size_t s = 0; s < r.sites.size(); ++s) {
            if (j.localBranch[s] < 0) { allPhased = false; j.reason = "no_branch"; break; }
            if (j.choiceBasis[s] != Basis::PairPhased && j.choiceBasis[s] != Basis::ThreadPhased) allPhased = false;
        }
        bool anyUnchosen = false;
        for (size_t s = 0; s < r.sites.size(); ++s) if (j.localBranch[s] < 0) anyUnchosen = true;
        if (anyUnchosen) { j.accepted = false; continue; }
        // path string P and element ranges
        std::string P = gv.oseq(j.oa);
        struct Elem { size_t from, to; int site; };
        std::vector<Elem> el;
        el.push_back({0, P.size(), -1});
        size_t tStart = 0;
        if (r.zero) {
            tStart = 0;
        } else {
            for (size_t x = 0; x + 1 < r.spine.size(); ++x) {
                const int s = r.siteAt[x];
                if (s >= 0) {
                    const size_t from = P.size();
                    for (uint32_t o : r.sites[static_cast<size_t>(s)].branches[static_cast<size_t>(j.localBranch[static_cast<size_t>(s)])]) {
                        P += gv.oseq(o).substr(static_cast<size_t>(kk));
                    }
                    el.push_back({from, P.size(), s});
                }
                const uint32_t nx = r.spine[x + 1];
                const size_t from = P.size();
                if (nx == kTarget) {
                    P += gv.oseq(j.ob).substr(static_cast<size_t>(kk));
                    tStart = P.size() - static_cast<size_t>(gv.len(j.ob));
                } else {
                    P += gv.oseq(nx).substr(static_cast<size_t>(kk));
                }
                el.push_back({from, P.size(), -1});
            }
        }
        const int64_t fs = j.posEndA + j.dA;
        const int64_t fe = static_cast<int64_t>(tStart) + j.posStartB - j.dB;
        // the record's own bases between each anchor and the gap must lie on the walk
        const std::string& rec = seqs[j.rec];
        const bool okA = j.dA >= 0 && static_cast<size_t>(j.posEndA + j.dA) <= P.size() &&
                         P.compare(static_cast<size_t>(j.posEndA), static_cast<size_t>(j.dA), rec, j.recAEnd,
                                   static_cast<size_t>(j.dA)) == 0;
        const size_t bFrom = j.wrap ? 0 : j.nEnd;
        const int64_t pbStart = static_cast<int64_t>(tStart) + j.posStartB - j.dB;
        const bool okB = j.dB >= 0 && pbStart >= 0 &&
                         P.compare(static_cast<size_t>(pbStart), static_cast<size_t>(j.dB), rec, bFrom,
                                   static_cast<size_t>(j.dB)) == 0;
        if (!okA || !okB) { ++st.flankMismatch; j.accepted = false; j.reason = "flank_mismatch"; continue; }
        if (fe < fs) {
            j.overlapBp = fs - fe;
            if (!j.wrap) { ++st.overlap; j.accepted = false; j.reason = "overlap"; continue; }
        } else {
            j.fill = P.substr(static_cast<size_t>(fs), static_cast<size_t>(fe - fs));
        }
        const int64_t F = fe - fs;
        // C4 (om2_clonal): a clonal join is sized by the isolate's own walk at the relatives' distance; a fill of
        // another length is another route (typically through a different copy of the element, the exact walk being
        // blocked by a placed node) and would assert the wrong sequence (MEASURED, dev E. coli GCF051549665v1:
        // N 2069 = true gap, filled with a 3134 bp route). Such a fill is refused; the sized gap stays.
        if (j.c1 && (j.c1->source == Source::Clonal || j.c1->clonalSized) && !j.wrap &&
            std::llabs(F - j.claimedN) > std::max<int64_t>(300, j.claimedN / 20)) {
            j.accepted = false;
            j.reason = "clonal_fill_length";
            continue;
        }
        const bool lengthOk = j.wrap || !sized(j.claimedN) ||
                              std::llabs(F - j.claimedN) <= std::max<int64_t>(1000, j.claimedN / 4);
        const bool spanned = j.spanThread + j.spanPair >= opt.minSpan;
        j.contig = opt.fill == FillMode::Contig && r.exhaustive && !r.cyclic && !r.hairpin && lengthOk &&
                   allPhased && (j.forced || spanned) && F >= 0;
        j.genome = !j.contig;
        const Basis base = j.contig ? Basis::GraphWalk : (r.usesRepeat ? Basis::IsolateRepeatExact : Basis::GraphWalk);
        if (F > 0) {
            for (const Elem& e : el) {
                const int64_t a = std::max<int64_t>(fs, static_cast<int64_t>(e.from));
                int64_t b = std::min<int64_t>(fe, static_cast<int64_t>(e.to));
                if (e.site >= 0 && e.to == e.from && a <= fe - 1) b = a + 1;   // empty branch: label one base
                if (b <= a) continue;
                const Basis bs = e.site >= 0 ? j.choiceBasis[static_cast<size_t>(e.site)] : base;
                FillSpan sp;
                sp.off = static_cast<uint32_t>(a - fs);
                sp.len = static_cast<uint32_t>(b - a);
                sp.basis = bs;
                if (!j.spans.empty() && j.spans.back().basis == bs &&
                    j.spans.back().off + j.spans.back().len == sp.off) {
                    j.spans.back().len += sp.len;
                } else if (!j.spans.empty() && j.spans.back().off + j.spans.back().len > sp.off) {
                    // an empty-branch label overlapping the next element: shorten the next
                    const uint32_t end = sp.off + sp.len;
                    const uint32_t prevEnd = j.spans.back().off + j.spans.back().len;
                    if (end > prevEnd) { sp.off = prevEnd; sp.len = end - prevEnd; j.spans.push_back(sp); }
                } else {
                    j.spans.push_back(sp);
                }
            }
        }
        // Integration: one repeat_variants row per bubble site this locus walks, with every branch
        // the isolate carries there (its depth-estimated carrier copies) and the branch placed here.
        // Raw sequences in the record's orientation; oriented and formatted when the row is written.
        if (in.variants && F > 0) {
            const size_t ji = static_cast<size_t>(&j - &J[0]);
            for (const Elem& e : el) {
                if (e.site < 0) continue;
                const size_t sidx = static_cast<size_t>(e.site);
                const Site& site = r.sites[sidx];
                const auto sit = siteId.find(site.key);
                if (sit == siteId.end()) continue;
                const uint32_t sid = sit->second;
                const SiteReport& rp = reports[sid];
                const LocusRef* ref = nullptr;
                for (const LocusRef& lr : siteLoci[sid]) if (lr.junction == ji && lr.site == sidx) ref = &lr;
                if (!ref) continue;
                const int64_t a = std::max<int64_t>(fs, static_cast<int64_t>(e.from));
                const int64_t b = std::min<int64_t>(fe, static_cast<int64_t>(e.to));
                if (a >= fe) continue;   // the site lies outside the bases this fill writes
                // a site cut by the fill's ends (it overlaps the record's own flank bases) is not listed
                if (a != static_cast<int64_t>(e.from) || b != static_cast<int64_t>(e.to)) continue;
                RepeatVariant v;
                v.family = "repeat";
                v.cluster = "bubble" + std::to_string(sid);
                v.site = rp.complex ? "COMPLEX" : rp.its ? "ITS" : "SNV";
                v.fillOffset = static_cast<uint32_t>(a - fs);
                v.siteLen = static_cast<uint32_t>(std::max<int64_t>(1, b - a));
                v.basis = j.choiceBasis[sidx];
                v.phasing = v.basis == Basis::PairPhased ? "pair" : v.basis == Basis::ThreadPhased ? "thread"
                          : v.basis == Basis::PriorAllocated ? "prior" : "none";
                for (size_t lb = 0; lb < site.branches.size(); ++lb) {
                    std::string al;
                    for (uint32_t o : site.branches[lb]) al += gv.oseq(o).substr(static_cast<size_t>(kk));
                    const int cb = lb < ref->localToCanon.size() ? ref->localToCanon[lb] : -1;
                    const float carriers = cb >= 0 && static_cast<size_t>(cb) < rp.carriers.size()
                                               ? static_cast<float>(rp.carriers[static_cast<size_t>(cb)]) : 0.0f;
                    if (static_cast<int>(lb) == j.localBranch[sidx]) v.placed = al;
                    v.alleles.emplace_back(std::move(al), carriers);
                }
                // Report the variant, not the whole bubble branch: drop the bases every allele shares at
                // both ends (an SNV bubble's branches are ~k bases long and differ in one), and move the
                // site to where the alleles differ. The placed allele is still exactly the fill's bases there.
                {
                    size_t minLen = SIZE_MAX;
                    for (const auto& al : v.alleles) minLen = std::min(minLen, al.first.size());
                    size_t P = 0, S = 0;
                    auto sameAt = [&](auto pick) {
                        const char c = pick(v.alleles.front().first);
                        for (const auto& al : v.alleles) if (pick(al.first) != c) return false;
                        return true;
                    };
                    while (minLen != SIZE_MAX && P < minLen && sameAt([&](const std::string& x) { return x[P]; })) ++P;
                    while (minLen != SIZE_MAX && P + S < minLen &&
                           sameAt([&](const std::string& x) { return x[x.size() - 1 - S]; })) ++S;
                    if (P + S > 0 && !v.alleles.empty()) {
                        for (auto& al : v.alleles) al.first = al.first.substr(P, al.first.size() - P - S);
                        v.placed = v.placed.substr(P, v.placed.size() - P - S);
                        v.fillOffset += static_cast<uint32_t>(P);
                        v.siteLen = static_cast<uint32_t>(std::max<size_t>(1, v.placed.size()));
                        if (v.placed.empty() && static_cast<int64_t>(v.fillOffset) >= fe - fs) continue;
                    }
                }
                std::string loci;
                for (const LocusRef& lr : siteLoci[sid]) {
                    if (lr.junction == ji) continue;
                    if (!loci.empty()) loci += ',';
                    loci += std::to_string(J[lr.junction].ledgerIdx);
                }
                v.loci = loci.empty() ? "." : loci;
                j.vrows.push_back(std::move(v));
            }
        }
    }

    // ---- 9. circles and dnaA ----
    for (CJ& j : J) {
        if (!j.wrap) continue;
        ++st.circlesTested;
        if (!j.accepted) continue;
        // no pair contradiction: neither record end links more strongly to another record's ends
        bool contra = false;
        const uint32_t gT = static_cast<uint32_t>(2 * (&j - &J[0])), gH = gT + 1;
        const LinkCount self = link(gT, gH);
        const uint32_t selfN = self.thread + self.pair;
        for (size_t i2 = 0; i2 < J.size() && !contra; ++i2) {
            if (J[i2].rec == j.rec) continue;
            for (uint32_t side = 0; side < 2; ++side) {
                const uint32_t go = static_cast<uint32_t>(2 * i2 + side);
                for (uint32_t mine : {gT, gH}) {
                    const LinkCount lc = link(mine, go);
                    const uint32_t n = lc.thread + lc.pair;
                    if (n >= opt.minContra && n >= selfN) contra = true;
                }
            }
        }
        if (contra) { j.accepted = false; j.contig = j.genome = false; j.reason = "pair_contradiction"; continue; }
        ++st.circlesClosed;
        if (j.overlapBp > 0) { ++st.circlesOverlap; j.contig = false; j.genome = true; }
        std::string circ = seqs[j.rec];
        if (j.overlapBp > 0) circ.resize(circ.size() - static_cast<size_t>(std::min<int64_t>(j.overlapBp, circ.size())));
        else circ += j.fill;
        if (st.dnaaLoaded) {
            const DnaaSketch::Hit h = dnaa.locate(circ);
            const int w = wrapOrd[static_cast<size_t>(&j - &J[0])];
            if (h.found && w >= 0 && w < 64) {
                ledger.rotateOffset[w] = h.reverse ? -2 - h.offset : h.offset;
                ++st.rotated;
                char b[96];
                std::snprintf(b, sizeof b, "dnaA=%lld%s votes=%u/%u", static_cast<long long>(h.offset),
                              h.reverse ? "(-)" : "(+)", h.votes, h.hits);
                j.reason = b;
            }
        }
    }

    // ledger
    auto fmtAllele = [](const std::string& al) -> std::string {
        if (al.empty()) return "-";
        if (al.size() <= 60) return al;
        uint64_t h = 1469598103934665603ULL;   // FNV-1a 64 of the oriented allele
        for (char c : al) { h ^= static_cast<uint8_t>(c); h *= 1099511628211ULL; }
        char b[64];
        std::snprintf(b, sizeof b, "len=%zu;h=%016llx", al.size(), static_cast<unsigned long long>(h));
        return b;
    };
    auto emitVariants = [&](const CJ& j, bool rc) {
        if (!in.variants || !j.accepted) return;
        const size_t F = j.fill.size();
        for (RepeatVariant v : j.vrows) {
            if (rc) {
                v.fillOffset = F >= static_cast<size_t>(v.fillOffset) + v.siteLen
                                   ? static_cast<uint32_t>(F - v.fillOffset - v.siteLen) : 0;
                v.placed = reverseComplement(v.placed);
                for (auto& al : v.alleles) al.first = reverseComplement(al.first);
            }
            v.placed = fmtAllele(v.placed);
            for (auto& al : v.alleles) al.first = fmtAllele(al.first);
            v.junction = j.ledgerIdx;
            in.variants->push_back(std::move(v));
            ++st.variantRows;
        }
    };
    ledger.j.reserve(ledger.j.size() + J.size());
    std::vector<std::pair<size_t, size_t>> wrapRows;   // (row in ledger.j, index in J): flanks after step 10
    for (size_t ci = 0; ci < J.size(); ++ci) {
        CJ& j = J[ci];
        if (j.c1Id >= 0) {
            // Integration: C1 recorded this junction. Its row takes the fill in its own orientation,
            // the contig-grade flag and the closeGaps decision; C1's verdict, tier (A once a
            // contig-grade fill replaces the gap), admission, class and flanks stand.
            Junction& R = in.c1->j[static_cast<size_t>(j.c1Id)];
            const bool rc = j.c1Orient == '-';
            if (j.bridge == Bridge::Open) R.allowCloseGaps = false;
            if (j.accepted && j.contig) { R.contigFilled = true; R.writtenN = 0; R.tier = Tier::A; }
            if (j.accepted && j.genome) {
                R.fillSeq = rc ? reverseComplement(j.fill) : j.fill;
                R.fillSpans.clear();
                if (rc) {
                    const uint32_t F = static_cast<uint32_t>(j.fill.size());
                    for (auto it = j.spans.rbegin(); it != j.spans.rend(); ++it) {
                        FillSpan sp = *it;
                        sp.off = F - it->off - it->len;
                        R.fillSpans.push_back(sp);
                    }
                } else {
                    R.fillSpans = j.spans;
                }
            }
            emitVariants(j, rc);
            continue;
        }
        Junction L;
        L.id = j.ledgerIdx;   // == c1Rows + ledger.j.size() here: assigned in creation order
        L.source = j.wrap ? Source::Wrap
                 : j.claimedN == 1 ? Source::JoinButt1
                 : j.claimedN == 2000 ? Source::LayoutCap2000
                 : j.claimedN == 100 ? Source::ResolverUnknown100 : Source::Join;
        L.a = {j.rec, true};
        L.b = {j.rec, false};
        L.claimedN = j.claimedN;
        L.writtenN = j.contig ? 0 : j.claimedN;
        // A wrap is a junction only once it is closed: an open wrap is recorded as not made (C3
        // reads writtenN >= 0 with no fill as "closed exactly"). An overlap wrap is closed by
        // removing the record's repeated end: claimedN = -overlap (C1's overlap convention).
        if (j.wrap && !j.accepted) L.writtenN = -1;
        if (j.wrap && j.accepted && j.overlapBp > 0) L.claimedN = -static_cast<int32_t>(j.overlapBp);
        L.verdict = j.verdict;
        L.tier = j.contig ? Tier::A : j.bridge == Bridge::PairThread ? Tier::B : j.bridge == Bridge::Forced ? Tier::C
               : j.bridge == Bridge::Tiebreak ? Tier::D : Tier::E;
        L.admit = Admit::Scaffold;
        L.cls = std::string("c2scan:") + (j.wrap ? "wrap" : "N" + std::to_string(j.claimedN));
        L.endA = j.endA;
        L.endB = j.endB;
        L.copyA = j.copyA;
        L.copyB = j.copyB;
        L.gmin = j.route.found ? static_cast<int32_t>(j.route.fmin - j.dA - j.dB) : -1;
        for (size_t c = 0; c < j.route.classes.size(); ++c) {
            WalkClass w;
            w.len = static_cast<int32_t>(j.route.classes[c] - j.dA - j.dB);
            w.n = j.route.classN[c];
            L.walks.push_back(w);
        }
        L.exhaustive = j.route.found && j.route.exhaustive;
        L.hairpin = j.route.hairpin;
        L.uA = static_cast<uint16_t>(std::min<size_t>(j.uA, 65535));
        L.uB = static_cast<uint16_t>(std::min<size_t>(j.uB, 65535));
        L.pairK = j.spanThread + j.spanPair;
        L.allowCloseGaps = !(j.bridge == Bridge::Open);
        L.contigFilled = j.accepted && j.contig;
        if (j.accepted && j.genome) { L.fillSeq = j.fill; L.fillSpans = j.spans; }
        const std::string& s = seqs[j.rec];
        if (!j.wrap) {   // captured before any contig-grade fill moves the coordinates
            L.flankL32 = s.substr(j.nStart >= 32 ? j.nStart - 32 : 0, std::min<size_t>(32, j.nStart));
            L.flankR32 = s.substr(j.nEnd, 32);
        }
        ledger.j.push_back(L);
        if (j.wrap) wrapRows.emplace_back(ledger.j.size() - 1, ci);
        emitVariants(j, false);
        if (j.wrap && !j.accepted) ++st.wrapsOpen;
    }

    // the records as the stage found them, for the side files' flank context
    const std::vector<std::string> preFill = in.outDir.empty() ? std::vector<std::string>() : seqs;

    // ---- 10. write: contig-grade fills into the sequence, genome-only fills to the ledger ----
    std::vector<size_t> contigOrder;
    for (size_t i = 0; i < J.size(); ++i) if (J[i].accepted && J[i].contig) contigOrder.push_back(i);
    std::sort(contigOrder.begin(), contigOrder.end(), [&](size_t a, size_t b) {
        if (J[a].rec != J[b].rec) return J[a].rec < J[b].rec;
        return J[a].nStart > J[b].nStart;   // right to left keeps earlier offsets valid
    });
    for (size_t i : contigOrder) {
        CJ& j = J[i];
        std::string& s = seqs[j.rec];
        if (j.wrap) { s += j.fill; ++st.circlesContig; }
        else s.replace(j.nStart, j.nEnd - j.nStart, j.fill);
        ++st.fillsContig;
        st.fillBpContig += j.fill.size();
    }
    for (const CJ& j : J) {
        if (j.accepted && j.genome) {
            ++st.fillsGenome;
            st.fillBpGenome += j.fill.size();
        }
    }
    // Wrap rows carry the record's two ends as their flanks (flankL32 = the last 32 bases, flankR32 =
    // the first 32), taken after the contig-grade fills above (a closed contig-grade wrap has its
    // fill appended): the writer finds a wrap's record by them.
    for (const auto& wr : wrapRows) {
        const std::string& s = seqs[J[wr.second].rec];
        if (s.size() < 32) continue;
        ledger.j[wr.first].flankL32 = s.substr(s.size() - 32);
        ledger.j[wr.first].flankR32 = s.substr(0, 32);
    }

    // ---- 11. side files ----
    if (!in.outDir.empty()) {
        const std::string dir = in.outDir + "/om2_close";
        if (util::makeDirs(dir)) {
            std::ofstream jt(dir + "/junctions.tsv");
            jt << "#om2-close junctions (C2). Coordinates are record coordinates at stage [4c/7], before gap "
                  "closing and polishing; the ledger carries 32-bp flanks to relocate them.\n";
            jt << "id\trecord\twrap\tn_start\tclaimed_n\tanchor_a\tanchor_b\tend_a\tend_b\tcopy_a\tcopy_b\tverdict"
                  "\twalks\tsites\texhaustive\thairpin\tcyclic\tfill_min\tfill_max\tu_a\tu_b\tforced\tspan_thread"
                  "\tspan_pair\tbridge\tgrade\tfill_len\toverlap\tbasis\treason\tprior_locus\n";
            for (const CJ& j : J) {
                std::string basis;
                for (const FillSpan& sp : j.spans) {
                    if (!basis.empty()) basis += ',';
                    basis += std::string(c2BasisName(sp.basis)) + ":" + std::to_string(sp.len);
                }
                const bool rf = j.route.found;
                jt << j.ledgerIdx << '\t' << j.rec << '\t' << (j.wrap ? 1 : 0) << '\t' << j.nStart << '\t' << j.claimedN
                   << '\t' << fmtAnchor(gv, j.anchoredA, j.oa, j.posEndA) << '\t'
                   << fmtAnchor(gv, j.anchoredB, j.ob, j.posStartB) << '\t' << endStr(j.endA) << '\t'
                   << endStr(j.endB) << '\t' << j.copyA << '\t' << j.copyB << '\t' << verdictStr(j.verdict) << '\t'
                   << (rf ? j.route.walks : 0) << '\t' << j.route.sites.size() << '\t'
                   << (rf && j.route.exhaustive ? 1 : 0) << '\t' << (j.route.hairpin ? 1 : 0) << '\t'
                   << (j.route.cyclic ? 1 : 0) << '\t' << (rf ? j.route.fmin - j.dA - j.dB : -1) << '\t'
                   << (rf ? j.route.fmax - j.dA - j.dB : -1) << '\t' << j.uA << '\t' << j.uB << '\t'
                   << (j.forced ? 1 : 0) << '\t' << j.spanThread << '\t' << j.spanPair << '\t' << bridgeName(j.bridge)
                   << '\t' << (!j.accepted ? "none" : j.contig ? "contig" : "genome") << '\t' << j.fill.size() << '\t'
                   << j.overlapBp << '\t' << (basis.empty() ? "-" : basis) << '\t'
                   << (j.reason.empty() ? "-" : j.reason) << '\t' << j.priorLocus << '\n';
            }
            std::ofstream fa(dir + "/fills.fasta");
            std::ofstream bed(dir + "/fills.mask.bed");
            std::ofstream ctx(dir + "/fills.context.fasta");   // fill with 1 kb of the record's own flanks
            for (const CJ& j : J) {
                if (!j.accepted || j.fill.empty()) continue;
                {
                    const std::string& s = preFill[j.rec];
                    const size_t lFrom = j.nStart > 1000 ? j.nStart - 1000 : 0;
                    const std::string left = s.substr(lFrom, j.nStart - lFrom);
                    const size_t rFrom = j.wrap ? 0 : j.nEnd;
                    const std::string right = s.substr(rFrom, std::min<size_t>(1000, s.size() - rFrom));
                    ctx << ">om2fill_" << j.ledgerIdx << "_ctx left=" << left.size() << " fill=" << j.fill.size()
                        << " right=" << right.size() << " grade=" << (j.contig ? "contig" : "genome") << '\n';
                    const std::string all = left + j.fill + right;
                    for (size_t p = 0; p < all.size(); p += 80) ctx << all.substr(p, 80) << '\n';
                }
                fa << ">om2fill_" << j.ledgerIdx << " record=" << j.rec << " n_start=" << j.nStart << " claimed_n="
                   << j.claimedN << " grade=" << (j.contig ? "contig" : "genome") << " bridge=" << bridgeName(j.bridge)
                   << (j.wrap ? " wrap=1" : "") << '\n';
                for (size_t p = 0; p < j.fill.size(); p += 80) fa << j.fill.substr(p, 80) << '\n';
                for (const FillSpan& sp : j.spans) {
                    bed << "om2fill_" << j.ledgerIdx << '\t' << sp.off << '\t' << sp.off + sp.len << '\t'
                        << c2BasisName(sp.basis) << '\t' << j.ledgerIdx << '\t' << bridgeName(j.bridge) << '\n';
                }
            }
            std::ofstream rv(dir + "/repeat_variants.tsv");
            rv << "#om2-close repeat variants: every branch of every bubble site on a bridged walk, the depth "
                  "estimate of the copies carrying it, and the loci it was allocated to (junction:branch:basis)\n";
            rv << "site\tclass\tbranch\tbranch_len\tdepth\tdepth_copies\tflow_copies\tmajority\tits_class\tallele\tloci\n";
            for (size_t sid = 0; sid < reports.size(); ++sid) {
                const SiteReport& rp = reports[sid];
                if (rp.carriers.empty()) continue;
                for (size_t b = 0; b < rp.carriers.size(); ++b) {
                    std::string loci;
                    for (const std::string& a : rp.assigned) {
                        const size_t c1p = a.find(':');
                        const size_t c2p = a.find(':', c1p + 1);
                        if (std::atoi(a.substr(c1p + 1, c2p - c1p - 1).c_str()) == static_cast<int>(b)) {
                            if (!loci.empty()) loci += ',';
                            loci += a;
                        }
                    }
                    rv << sid << '\t' << (rp.complex ? "COMPLEX" : rp.its ? "ITS" : "SNV") << '\t' << b << '\t'
                       << rp.blen[b] << '\t' << rp.bcov[b] << '\t' << rp.carriers[b] << '\t'
                       << (rp.flow.empty() ? std::string("NA") : std::to_string(rp.flow[b])) << '\t'
                       << (static_cast<int>(b) == rp.majority ? 1 : 0) << '\t'
                       << (rp.bclass.empty() || rp.bclass[b] < 0 ? std::string("-") : std::to_string(rp.bclass[b])) << '\t'
                       << rp.allele[b] << '\t'
                       << (loci.empty() ? "-" : loci) << '\n';
                }
            }
            std::ofstream cc(dir + "/circles.tsv");
            cc << "record\tlength\twalks\tclosed\tgrade\tfill_len\toverlap\trotate_offset\treason\n";
            for (const CJ& j : J) {
                if (!j.wrap) continue;
                cc << j.rec << '\t' << j.nStart << '\t' << (j.route.found ? j.route.walks : 0) << '\t'
                   << (j.accepted ? 1 : 0) << '\t' << (!j.accepted ? "none" : j.contig ? "contig" : "genome") << '\t'
                   << j.fill.size() << '\t' << j.overlapBp << '\t'
                   << (wrapOrd[static_cast<size_t>(&j - &J[0])] >= 0 && wrapOrd[static_cast<size_t>(&j - &J[0])] < 64
                           ? ledger.rotateOffset[wrapOrd[static_cast<size_t>(&j - &J[0])]] : -1)
                   << '\t' << (j.reason.empty() ? "-" : j.reason)
                   << '\n';
            }
        }
    }
    st.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    if (in.verbose) {
        std::fprintf(stderr, "[4c/7] om2 close: %zu junctions (%zu anchored), %zu bridges, %zu contig fills (%zu bp), "
                             "%zu genome-only fills (%zu bp), %zu variant sites, %zu/%zu circles closed, %.1fs\n",
                     st.junctions + st.wraps, st.anchored, st.bridges, st.fillsContig, st.fillBpContig, st.fillsGenome,
                     st.fillBpGenome, st.variantSites, st.circlesClosed, st.circlesTested, st.seconds);
    }
    return st;
}

}  // namespace om2
}  // namespace ts
