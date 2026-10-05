#include "om2_evidence.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <functional>
#include <map>
#include <queue>
#include <set>
#include <thread>
#include <unordered_set>

#include "kmer.h"

namespace ts {
namespace om2 {

namespace {

using u128 = unsigned __int128;
constexpr int kA = 63;   // anchor length (census FL)
const u128 kMask126 = (static_cast<u128>(1) << 126) - 1;

inline int code(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return -1;
    }
}

inline uint64_t fp128(u128 x) {
    const uint64_t lo = static_cast<uint64_t>(x);
    const uint64_t hi = static_cast<uint64_t>(x >> 64);
    return mix64(lo ^ mix64(hi + 0x632BE59BD9B4E019ULL));
}

// Canonical fingerprint of the 63-mer at s[from, from+63). False on a non-ACGT base.
bool fp63(const std::string& s, size_t from, uint64_t& fp) {
    if (from + kA > s.size()) return false;
    u128 f = 0, r = 0;
    for (int i = 0; i < kA; ++i) {
        const int b = code(s[from + static_cast<size_t>(i)]);
        if (b < 0) return false;
        f = (f << 2) | static_cast<u128>(b);
        r |= static_cast<u128>(3 - b) << (2 * i);
    }
    fp = fp128(f < r ? f : r);
    return true;
}

// Every 63-mer of s: fn(start, canonical fingerprint, strand) with strand 0 when the forward
// 63-mer is the canonical one.
template <typename F>
void forEach63(const std::string& s, F&& fn) {
    if (s.size() < static_cast<size_t>(kA)) return;
    u128 f = 0, r = 0;
    int valid = 0;
    for (size_t i = 0; i < s.size(); ++i) {
        const int b = code(s[i]);
        if (b < 0) { valid = 0; f = r = 0; continue; }
        f = ((f << 2) | static_cast<u128>(b)) & kMask126;
        r = (r >> 2) | (static_cast<u128>(3 - b) << 124);
        if (++valid < kA) continue;
        fn(i + 1 - static_cast<size_t>(kA), fp128(f < r ? f : r), f < r ? 0 : 1);
    }
}

std::string rcString(const std::string& s) { return reverseComplement(s); }

double nowMs() {
    return std::chrono::duration<double, std::milli>(
               std::chrono::steady_clock::now().time_since_epoch()).count();
}

// Length-weighted median depth over nodes >= 2k (gate_eval.theta_of).
double weightedMedian(std::vector<std::pair<double, long>>& xs) {
    if (xs.empty()) return 0;
    std::sort(xs.begin(), xs.end());
    long tot = 0;
    for (const auto& x : xs) tot += x.second;
    long acc = 0;
    for (const auto& x : xs) {
        acc += x.second;
        if (acc >= (tot + 1) / 2) return x.first;
    }
    return xs.back().first;
}

}  // namespace

const char* sourceName(Source s) {
    switch (s) {
        case Source::Resolver: return "resolver";
        case Source::ResolverUnknown100: return "resolver_u100";
        case Source::Join: return "join";
        case Source::JoinButt1: return "join_butt1";
        case Source::JoinOverlap: return "join_overlap";
        case Source::Layout: return "layout";
        case Source::LayoutCap2000: return "layout_cap2000";
        case Source::LayoutOverlap: return "layout_overlap";
        case Source::Wrap: return "wrap";
        case Source::LayoutButt1: return "layout_butt1";
        case Source::Clonal: return "clonal";
    }
    return "?";
}
const char* verdictName(Verdict v) {
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
const char* tierName(Tier t) {
    static const char* n[] = {"A", "B", "C", "D", "E"};
    return n[static_cast<int>(t)];
}
const char* admitName(Admit a) {
    switch (a) {
        case Admit::Scaffold: return "scaffold";
        case Admit::GenomeOnly: return "genome_only";
        case Admit::Break: return "break";
    }
    return "?";
}
const char* endClassName(EndClass e) {
    switch (e) {
        case EndClass::Unknown: return "unknown";
        case EndClass::Unique: return "unique";
        case EndClass::Repeat: return "repeat";
        case EndClass::DeadEnd: return "dead_end";
    }
    return "?";
}

bool fingerprint63(const std::string& s, size_t from, uint64_t& fp) { return fp63(s, from, fp); }
void forEachFingerprint63(const std::string& s, const std::function<void(size_t, uint64_t)>& fn) {
    forEach63(s, [&](size_t pos, uint64_t fp, int) { fn(pos, fp); });
}

// ============================ EvidenceIndex ==========================================

std::unique_ptr<EvidenceIndex> EvidenceIndex::build(const UnitigGraph& g,
                                                    const EvidenceOptions& opt) {
    const double t0 = nowMs();
    std::unique_ptr<EvidenceIndex> ix(new EvidenceIndex());
    ix->g_ = &g;
    ix->opt_ = opt;
    ix->k_ = g.k();
    const size_t n = g.nodes.size();
    ix->mult_.assign(n, Multiplicity());
    ix->placed_.assign(n, 0);
    const long k = ix->k_;

    // ---- theta: global and per component ------------------------------------
    std::vector<std::pair<double, long>> xs;
    size_t totalBases = 0;
    for (uint32_t u = 0; u < n; ++u) {
        const Unitig& x = g.nodes[u];
        if (x.deleted) continue;
        ++ix->liveNodes;
        totalBases += x.seq.size();
        const long L = static_cast<long>(x.seq.size());
        if (L >= 2 * k && x.coverage > 0) xs.emplace_back(x.coverage, L - k + 1);
    }
    ix->theta_ = weightedMedian(xs);

    const std::vector<std::vector<uint32_t>> comps = g.components();
    ix->components = comps.size();
    size_t scBases = 0;
    for (size_t ci = 0; ci < comps.size(); ++ci) {
        std::vector<std::pair<double, long>> cx;
        long eligible = 0;
        for (uint32_t u : comps[ci]) {
            const Unitig& x = g.nodes[u];
            const long L = static_cast<long>(x.seq.size());
            if (L >= 2 * k && x.coverage > 0) { cx.emplace_back(x.coverage, L - k + 1); eligible += L - k + 1; }
        }
        double th = ix->theta_;
        if (opt.perComponentTheta && eligible >= opt.compMinEligible) {
            const double tc = weightedMedian(cx);
            if (tc > 0) th = tc;
        }
        for (uint32_t u : comps[ci]) {
            const Unitig& x = g.nodes[u];
            Multiplicity& m = ix->mult_[u];
            m.component = static_cast<uint32_t>(ci);
            m.copy = th > 0 ? static_cast<float>(x.coverage / th) : 0.0f;
            m.hardSingle = static_cast<long>(x.seq.size()) >= opt.scMinLen && th > 0 &&
                           x.coverage <= opt.scMaxTheta * th;
            if (m.hardSingle) {
                m.lo = m.hi = 1;
                ++ix->singleCopy;
                scBases += x.seq.size();
            } else {
                const int est = std::max(1, static_cast<int>(std::lround(m.copy)));
                const int spread = est > 4 ? 2 : 1;
                m.lo = static_cast<uint16_t>(std::max(1, est - spread));
                m.hi = static_cast<uint16_t>(std::min(65535, est + spread));
                m.uncertain = std::fabs(m.copy - static_cast<float>(est)) > 0.35f;
                if (m.copy >= 1.5f) ++ix->repeatNodes;
                if (m.uncertain) ++ix->uncertainNodes;
            }
        }
    }
    ix->scFraction = totalBases ? static_cast<double>(scBases) / static_cast<double>(totalBases) : 0;

    if (static_cast<long>(ix->liveNodes) > opt.tangleMax) {
        ix->abstained_ = true;
        ix->abstainWhy_ = "nodes";
    } else if (ix->scFraction < opt.tangleScFrac) {
        ix->abstained_ = true;
        ix->abstainWhy_ = "single_copy_fraction";
    }
    if (ix->abstained_) {
        ix->buildMs = nowMs() - t0;
        return ix;
    }

    // ---- 63-mer index over every live unitig --------------------------------
    ix->index_.reserve(totalBases);
    for (uint32_t u = 0; u < n; ++u) {
        const Unitig& x = g.nodes[u];
        if (x.deleted) continue;
        forEach63(x.seq, [&](size_t pos, uint64_t fp, int strand) {
            ix->index_.push_back(Ent{fp, u, static_cast<uint32_t>((pos << 1) | static_cast<size_t>(strand))});
        });
    }
    std::sort(ix->index_.begin(), ix->index_.end(), [](const Ent& a, const Ent& b) {
        if (a.fp != b.fp) return a.fp < b.fp;
        if (a.node != b.node) return a.node < b.node;
        return a.posStrand < b.posStrand;
    });

    // ---- middle 63-mers of the hard single-copy unitigs (placement test) ----
    for (uint32_t u = 0; u < n; ++u) {
        const Unitig& x = g.nodes[u];
        if (x.deleted || !ix->mult_[u].hardSingle || x.seq.size() < 2 * static_cast<size_t>(kA)) continue;
        uint64_t fp;
        if (!fp63(x.seq, x.seq.size() / 2 - 31, fp)) continue;
        auto it = ix->midIndex_.find(fp);
        if (it == ix->midIndex_.end()) ix->midIndex_.emplace(fp, u);
        else it->second = UINT32_MAX;
    }
    ix->buildMs = nowMs() - t0;
    return ix;
}

std::string EvidenceIndex::oseq(ONode nd) const {
    const std::string& s = g_->nodes[nd.id].seq;
    return nd.o == 0 ? s : rcString(s);
}

int EvidenceIndex::locate(const std::string& s, size_t from, Hit& hit) const {
    uint64_t fp;
    if (!fp63(s, from, fp)) return -1;
    auto lo = std::lower_bound(index_.begin(), index_.end(), fp,
                               [](const Ent& e, uint64_t v) { return e.fp < v; });
    int found = 0;
    std::string q, qrc;
    for (auto it = lo; it != index_.end() && it->fp == fp; ++it) {
        if (q.empty()) {
            q = s.substr(from, kA);
            for (char& c : q) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
            qrc = rcString(q);
        }
        const std::string& ns = g_->nodes[it->node].seq;
        const size_t pos = it->posStrand >> 1;
        if (ns.compare(pos, kA, q) == 0) {
            if (++found == 1) hit = Hit{ONode{it->node, 0}, static_cast<int32_t>(pos)};
        } else if (ns.compare(pos, kA, qrc) == 0) {
            if (++found == 1)
                hit = Hit{ONode{it->node, 1},
                          static_cast<int32_t>(ns.size() - (pos + static_cast<size_t>(kA)))};
        }
        if (found > 1) return 2;
    }
    return found;
}

EndAnchor EvidenceIndex::anchorExit(const std::string& s, int32_t minSlide) const {
    EndAnchor a;
    if (abstained_) { a.why = "abstain"; return a; }
    const long L = static_cast<long>(s.size());
    const int step = std::max(1, opt_.anchorStep);
    const long d0 = minSlide > 0 ? ((minSlide + step - 1) / step) * step : 0;
    const char* last = "none";
    for (long d = d0; d < opt_.anchorSlide; d += step) {
        const long lo = L - d - kA;
        if (lo < 0) { a.why = "edge"; ++anchorsFailed; return a; }
        Hit h;
        const int r = locate(s, static_cast<size_t>(lo), h);
        if (r < 0) { a.why = "hit_N"; ++anchorsFailed; return a; }
        if (r == 1) {
            a.ok = true;
            a.why = "ok";
            a.node = h.node;
            a.off = h.off + kA;
            a.slide = static_cast<int32_t>(d);
            ++anchorsOk;
            return a;
        }
        last = r == 0 ? "slide_absent" : "slide_multi";
    }
    a.why = last;
    ++anchorsFailed;
    return a;
}

EndAnchor EvidenceIndex::anchorEntry(const std::string& s, int32_t minSlide) const {
    EndAnchor a;
    if (abstained_) { a.why = "abstain"; return a; }
    const long L = static_cast<long>(s.size());
    const int step = std::max(1, opt_.anchorStep);
    const long d0 = minSlide > 0 ? ((minSlide + step - 1) / step) * step : 0;
    const char* last = "none";
    for (long d = d0; d < opt_.anchorSlide; d += step) {
        if (d + kA > L) { a.why = "edge"; ++anchorsFailed; return a; }
        Hit h;
        const int r = locate(s, static_cast<size_t>(d), h);
        if (r < 0) { a.why = "hit_N"; ++anchorsFailed; return a; }
        if (r == 1) {
            a.ok = true;
            a.why = "ok";
            a.node = h.node;
            a.off = h.off;
            a.slide = static_cast<int32_t>(d);
            ++anchorsOk;
            return a;
        }
        last = r == 0 ? "slide_absent" : "slide_multi";
    }
    a.why = last;
    ++anchorsFailed;
    return a;
}

uint32_t EvidenceIndex::locateUnique(const std::string& seq, size_t at) const {
    if (abstained_) return UINT32_MAX;
    Hit h;
    return locate(seq, at, h) == 1 ? h.node.id : UINT32_MAX;
}

// gate_eval.dmin: shortest fill from offset e of `sn` to offset f of `tn`; -1 when none.
int64_t EvidenceIndex::dmin(ONode sn, int32_t e, ONode tn, int32_t f) const {
    const int64_t k1 = k_ - 1;
    int64_t best = -1;
    if (sn == tn && f >= e) best = f - e;
    const size_t n2 = g_->nodes.size() * 2;
    std::vector<int64_t> seen(n2, INT64_MAX);
    using QE = std::pair<int64_t, uint32_t>;
    std::priority_queue<QE, std::vector<QE>, std::greater<QE>> pq;
    pq.push({len(sn.id), sn.id * 2 + sn.o});
    while (!pq.empty()) {
        const QE top = pq.top();
        pq.pop();
        const int64_t wlen = top.first;
        const ONode node{top.second >> 1, static_cast<uint8_t>(top.second & 1)};
        if (best >= 0 && wlen - k1 - e >= best) break;
        if (seen[top.second] < wlen) continue;
        for (const Link& l : g_->exits(node.id, node.o)) {
            if (g_->nodes[l.to].deleted) continue;
            const ONode nxt{l.to, static_cast<uint8_t>(UnitigGraph::enterOrient(l))};
            const int64_t posT = wlen - k1;
            if (nxt == tn) {
                const int64_t fl = posT + f - e;
                if (fl >= 0 && (best < 0 || fl < best)) best = fl;
            }
            const int64_t nw = wlen + len(nxt.id) - k1;
            if (nw - e > opt_.gminCap) continue;
            const uint32_t si = nxt.id * 2 + nxt.o;
            if (nw < seen[si]) {
                seen[si] = nw;
                pq.push({nw, si});
            }
        }
    }
    return best;
}

// skip_path_check.spath: the shortest walk itself, as oriented nodes from sn to tn.
bool EvidenceIndex::spath(ONode sn, int32_t e, ONode tn, int32_t f, std::vector<ONode>& path) const {
    path.clear();
    const int64_t k1 = k_ - 1;
    if (sn == tn && f >= e) { path.push_back(sn); return true; }
    const size_t n2 = g_->nodes.size() * 2;
    std::vector<int64_t> seen(n2, INT64_MAX);
    std::vector<uint32_t> prev(n2, UINT32_MAX);
    using QE = std::pair<int64_t, uint32_t>;
    std::priority_queue<QE, std::vector<QE>, std::greater<QE>> pq;
    const uint32_t s0 = sn.id * 2 + sn.o;
    seen[s0] = len(sn.id);
    pq.push({len(sn.id), s0});
    int64_t best = -1;
    uint32_t bestPrev = UINT32_MAX;
    while (!pq.empty()) {
        const QE top = pq.top();
        pq.pop();
        const int64_t wlen = top.first;
        if (best >= 0 && wlen - k1 - e >= best) break;
        if (seen[top.second] < wlen) continue;
        const ONode node{top.second >> 1, static_cast<uint8_t>(top.second & 1)};
        for (const Link& l : g_->exits(node.id, node.o)) {
            if (g_->nodes[l.to].deleted) continue;
            const ONode nxt{l.to, static_cast<uint8_t>(UnitigGraph::enterOrient(l))};
            if (nxt == tn) {
                const int64_t fl = wlen - k1 + f - e;
                if (fl >= 0 && (best < 0 || fl < best)) { best = fl; bestPrev = top.second; }
            }
            const int64_t nw = wlen + len(nxt.id) - k1;
            if (nw - e > opt_.gminCap) continue;
            const uint32_t si = nxt.id * 2 + nxt.o;
            if (nw < seen[si]) {
                seen[si] = nw;
                prev[si] = top.second;
                pq.push({nw, si});
            }
        }
    }
    if (best < 0) return false;
    std::vector<uint32_t> rev{tn.id * 2 + tn.o, bestPrev};
    for (size_t guard = 0; rev.back() != s0 && guard < n2 + 2; ++guard) {
        const uint32_t p = prev[rev.back()];
        if (p == UINT32_MAX) return false;
        rev.push_back(p);
    }
    if (rev.back() != s0) return false;
    for (auto it = rev.rbegin(); it != rev.rend(); ++it)
        path.push_back(ONode{*it >> 1, static_cast<uint8_t>(*it & 1)});
    return true;
}

// gap_graph_census_lib.Graph.walks: every walk up to walkFillMax, at most maxWalks, with the
// step cap; fills are compared as sequences. With `prune`, a walk may pass through at most
// repeatTraverseMax bases of repeat nodes (the walk-length classes of DESIGN C1a).
void EvidenceIndex::enumerate(ONode sn, int32_t e, ONode tn, int32_t f, bool prune, Enumeration& en) const {
    const int64_t k1 = k_ - 1;
    struct T { ONode n; int32_t parent; int64_t wlen; int64_t rep; };
    std::vector<T> arena;
    std::vector<int32_t> stack;
    long steps = 0;
    auto addWalk = [&](const std::string& fill, const std::vector<ONode>& p) {
        ++en.walks;
        if (en.fills.insert(fill).second) en.paths.push_back(p);
    };
    if (sn == tn && f >= e) addWalk(oseq(sn).substr(static_cast<size_t>(e), static_cast<size_t>(f - e)), {sn});
    arena.push_back(T{sn, -1, len(sn.id), 0});
    stack.push_back(0);
    bool stop = false;
    while (!stack.empty() && !stop) {
        const int32_t idx = stack.back();
        stack.pop_back();
        const T t = arena[static_cast<size_t>(idx)];
        for (const Link& l : g_->exits(t.n.id, t.n.o)) {
            if (g_->nodes[l.to].deleted) continue;
            const ONode nxt{l.to, static_cast<uint8_t>(UnitigGraph::enterOrient(l))};
            if (++steps > opt_.walkSteps) { en.capped = true; stop = true; break; }
            const int64_t nl = len(nxt.id);
            const int64_t posT = t.wlen - k1;
            if (nxt == tn) {
                const int64_t fillEnd = posT + f;
                if (fillEnd >= e && fillEnd - e <= opt_.walkFillMax) {
                    std::vector<ONode> p;
                    for (int32_t c = idx; c >= 0; c = arena[static_cast<size_t>(c)].parent)
                        p.push_back(arena[static_cast<size_t>(c)].n);
                    std::reverse(p.begin(), p.end());
                    p.push_back(nxt);
                    std::string spelled = oseq(p[0]);
                    for (size_t i = 1; i < p.size(); ++i) spelled += oseq(p[i]).substr(static_cast<size_t>(k1));
                    addWalk(spelled.substr(static_cast<size_t>(e), static_cast<size_t>(fillEnd - e)), p);
                    if (static_cast<int>(en.walks) >= opt_.maxWalks) { en.capped = true; stop = true; break; }
                }
            }
            if (t.wlen + nl - k1 - e <= opt_.walkFillMax) {
                const int64_t rep = t.rep + (mult_[nxt.id].hardSingle ? 0 : nl - k1);
                if (prune && rep > opt_.repeatTraverseMax) { en.pruned = true; continue; }
                arena.push_back(T{nxt, idx, t.wlen + nl - k1, rep});
                stack.push_back(static_cast<int32_t>(arena.size() - 1));
            }
        }
    }
}

namespace {
bool pathHairpin(const std::vector<ONode>& p) {
    std::unordered_map<uint32_t, uint8_t> seenO;
    for (const ONode& x : p) {
        auto it = seenO.find(x.id);
        if (it == seenO.end()) seenO.emplace(x.id, static_cast<uint8_t>(1u << x.o));
        else it->second |= static_cast<uint8_t>(1u << x.o);
    }
    for (const auto& kv : seenO) if (kv.second == 3) return true;
    return false;
}
}  // namespace

WalkSet EvidenceIndex::walks(const EndAnchor& A, const EndAnchor& B) const {
    WalkSet ws;
    if (abstained_ || !A.ok || !B.ok) return ws;
    const int64_t dm = dmin(A.node, A.off, B.node, B.off);
    if (dm < 0) return ws;
    ws.reachable = true;
    ws.fillMin = static_cast<int32_t>(dm);
    ws.gmin = static_cast<int32_t>(dm) - A.slide - B.slide;
    std::vector<ONode> p;
    if (spath(A.node, A.off, B.node, B.off, p) && p.size() >= 2) {
        if (pathHairpin(p)) ws.hairpin = true;
        for (size_t i = 1; i + 1 < p.size(); ++i) {
            const ONode& x = p[i];
            ws.inner.push_back(x);
            const Multiplicity& m = mult_[x.id];
            const uint64_t L = g_->nodes[x.id].seq.size();
            if (m.hardSingle) {
                ++ws.innerUnique;
                ws.innerUniqueBp += L;
                if (placed(x.id)) ++ws.innerUniquePlaced;
            } else {
                ws.innerRepeatBp += L > static_cast<uint64_t>(k_) ? L - static_cast<uint64_t>(k_) + 1 : 0;
                ws.innerMaxCopy = std::max(ws.innerMaxCopy, m.copy);
            }
        }
    }
    // The census enumeration decides the verdict (distinct fills, capped) exactly as the
    // gate_eval oracle does; the repeat-bounded one gives the reported walk-length classes.
    Enumeration census;
    enumerate(A.node, A.off, B.node, B.off, false, census);
    ws.walks = census.walks;
    ws.distinct = static_cast<uint32_t>(census.fills.size());
    ws.capped = census.capped;
    ws.exhaustive = !census.capped;
    for (const auto& path : census.paths) if (pathHairpin(path)) ws.hairpin = true;
    Enumeration bounded;
    enumerate(A.node, A.off, B.node, B.off, true, bounded);
    ws.repeatPruned = bounded.pruned;
    ws.classesCapped = bounded.capped;
    for (const auto& path : bounded.paths) if (pathHairpin(path)) ws.hairpin = true;
    std::vector<int32_t> lens;
    for (const std::string& fl : bounded.fills) lens.push_back(static_cast<int32_t>(fl.size()));
    std::sort(lens.begin(), lens.end());
    for (int32_t L : lens) {
        const int32_t gap = L - A.slide - B.slide;
        if (!ws.classes.empty() && gap - ws.classes.back().len <= opt_.classTol) {
            ++ws.classes.back().n;
        } else {
            ws.classes.push_back(WalkClass{gap, 1});
        }
    }
    // The shortest walk exists (Dijkstra) even when a depth-first enumeration spent its budget
    // elsewhere in a tangle; its class is always reported.
    if (ws.classes.empty() || ws.classes.front().len > ws.gmin + opt_.classTol)
        ws.classes.insert(ws.classes.begin(), WalkClass{ws.gmin, 1});
    return ws;
}

FirstUnique EvidenceIndex::firstUnique(ONode start) const {
    FirstUnique fu;
    if (abstained_ || start.id == UINT32_MAX) return fu;
    const int64_t k1 = k_ - 1;
    std::set<uint32_t> out;
    std::vector<std::pair<ONode, int64_t>> stack{{start, 0}};
    std::unordered_set<uint32_t> seen;
    long steps = 0;
    while (!stack.empty()) {
        const auto cur = stack.back();
        stack.pop_back();
        for (const Link& l : g_->exits(cur.first.id, cur.first.o)) {
            if (g_->nodes[l.to].deleted) continue;
            const ONode nxt{l.to, static_cast<uint8_t>(UnitigGraph::enterOrient(l))};
            if (++steps > opt_.uSteps) {
                fu.capped = true;
                fu.ids.assign(out.begin(), out.end());
                return fu;
            }
            if (mult_[nxt.id].hardSingle && nxt.id != start.id) { out.insert(nxt.id); continue; }
            const int64_t nd = cur.second + len(nxt.id) - k1;
            const uint32_t si = nxt.id * 2 + nxt.o;
            if (nd > opt_.uBound || seen.count(si)) continue;
            seen.insert(si);
            stack.push_back({nxt, nd});
        }
    }
    fu.ids.assign(out.begin(), out.end());
    return fu;
}

void EvidenceIndex::setPlaced(const std::vector<std::string>& pieces) {
    placed_.assign(g_->nodes.size(), 0);
    if (abstained_ || midIndex_.empty()) return;
    for (const std::string& s : pieces) {
        forEach63(s, [&](size_t pos, uint64_t fp, int) {
            auto it = midIndex_.find(fp);
            if (it == midIndex_.end() || it->second == UINT32_MAX) return;
            const std::string& ns = g_->nodes[it->second].seq;
            const size_t mid = ns.size() / 2 - 31;
            const std::string q = s.substr(pos, kA);
            if (ns.compare(mid, kA, q) == 0 || ns.compare(mid, kA, rcString(q)) == 0)
                placed_[it->second] = 1;
        });
    }
}

// ============================ pair evidence ===========================================

namespace {

// Open addressing on 62-bit canonical 31-mers; key stored +1 so 0 marks an empty slot.
class U64Table {
public:
    explicit U64Table(size_t expected) {
        size_t cap = 16;
        while (cap < expected * 2) cap <<= 1;
        keys_.assign(cap, 0);
        vals_.assign(cap, 0);
        mask_ = cap - 1;
    }
    uint64_t* find(uint64_t key) {
        size_t i = mix64(key) & mask_;
        while (keys_[i] != 0) {
            if (keys_[i] == key + 1) return &vals_[i];
            i = (i + 1) & mask_;
        }
        return nullptr;
    }
    const uint64_t* find(uint64_t key) const { return const_cast<U64Table*>(this)->find(key); }
    // Returns true when newly inserted.
    bool insert(uint64_t key, uint64_t val, uint64_t*& slot) {
        size_t i = mix64(key) & mask_;
        while (keys_[i] != 0) {
            if (keys_[i] == key + 1) { slot = &vals_[i]; return false; }
            i = (i + 1) & mask_;
        }
        keys_[i] = key + 1;
        vals_[i] = val;
        slot = &vals_[i];
        return true;
    }
private:
    std::vector<uint64_t> keys_, vals_;
    size_t mask_ = 0;
};

constexpr uint64_t kAmbig = UINT64_MAX;
constexpr int kPK = 31;
constexpr int kMaxProbes = 12;
constexpr int kMinVotes = 2;

inline uint64_t packPos(uint32_t piece, uint32_t pos, int strand) {
    return (static_cast<uint64_t>(piece) << 33) | (static_cast<uint64_t>(pos) << 1) |
           static_cast<uint64_t>(strand & 1);
}
inline uint64_t pairKey(uint32_t a, uint32_t b) {
    return (static_cast<uint64_t>(a) << 32) | static_cast<uint64_t>(b);
}

// Ends are keyed by their terminal kPortKey bases. A contig that stops at a repeat carries the
// repeat's first k-1 (<= 127) bases, shared by every contig stopping there, so the key must
// reach past them into the contig's own sequence.
constexpr size_t kPortKey = 250;

// Exit-oriented terminal bases of an origin end.
std::string exitKeyOf(const std::string& s, int end) {
    if (s.size() < kPortKey) return std::string();
    return end == 1 ? s.substr(s.size() - kPortKey) : rcString(s.substr(0, kPortKey));
}

void buildPortMap(const std::vector<std::string>& origin,
                  std::unordered_map<std::string, uint32_t>& portOf) {
    for (uint32_t c = 0; c < origin.size(); ++c) {
        for (int e = 0; e < 2; ++e) {
            const std::string key = exitKeyOf(origin[c], e);
            if (key.empty()) continue;
            auto it = portOf.find(key);
            if (it == portOf.end()) portOf.emplace(key, c * 2 + static_cast<uint32_t>(e));
            else it->second = UINT32_MAX;
        }
    }
}

uint32_t exitLookup(const std::unordered_map<std::string, uint32_t>& portOf, const std::string& s) {
    if (s.size() < kPortKey) return UINT32_MAX;
    auto it = portOf.find(s.substr(s.size() - kPortKey));
    return it == portOf.end() ? UINT32_MAX : it->second;
}
uint32_t entryLookup(const std::unordered_map<std::string, uint32_t>& portOf, const std::string& s) {
    if (s.size() < kPortKey) return UINT32_MAX;
    auto it = portOf.find(rcString(s.substr(0, kPortKey)));
    return it == portOf.end() ? UINT32_MAX : it->second;
}

}  // namespace

PortMap::PortMap(const std::vector<std::string>& origin) { buildPortMap(origin, portOf_); }
uint32_t PortMap::exitPort(const std::string& s) const { return exitLookup(portOf_, s); }
uint32_t PortMap::entryPort(const std::string& s) const { return entryLookup(portOf_, s); }

PairIndex::PairIndex(const std::vector<std::string>& origin, const SequenceStore& reads,
                     const InsertModel& insert, int threads, const EvidenceOptions& opt)
    : opt_(opt) {
    const uint32_t n = static_cast<uint32_t>(origin.size());
    len_.resize(n);
    for (uint32_t c = 0; c < n; ++c) len_[c] = origin[c].size();
    ctrl_.assign(static_cast<size_t>(n) * 2, -1);
    outTotal_.assign(static_cast<size_t>(n) * 2, 0);
    buildPortMap(origin, portOf_);
    if (!insert.usable || reads.pairCount() == 0 || origin.empty()) return;
    if (threads <= 0) threads = 1;
    maxFrag = static_cast<int32_t>(insert.mean + opt.pairMaxFragSd * insert.stddev);
    if (maxFrag <= 0) return;
    const long W = opt.pairCtrl + maxFrag + static_cast<long>(reads.maxReadLength()) + 64;

    // ---- index the 31-mers near every end; uniqueness judged over all pieces ----
    auto inWindow = [&](size_t L, size_t p) {
        return static_cast<long>(p) < W || static_cast<long>(p + kPK) > static_cast<long>(L) - W;
    };
    size_t windowPositions = 0;
    for (uint32_t c = 0; c < n; ++c) {
        const size_t L = origin[c].size();
        if (L < static_cast<size_t>(kPK)) continue;
        windowPositions += std::min<size_t>(L - kPK + 1, static_cast<size_t>(2 * W));
    }
    U64Table table(windowPositions + 16);
    const uint64_t mask62 = (1ULL << (2 * kPK)) - 1;
    for (int pass = 0; pass < 2; ++pass) {
        for (uint32_t c = 0; c < n; ++c) {
            const std::string& s = origin[c];
            uint64_t f = 0, r = 0;
            int valid = 0;
            for (size_t i = 0; i < s.size(); ++i) {
                const int b = code(s[i]);
                if (b < 0) { valid = 0; f = r = 0; continue; }
                f = ((f << 2) | static_cast<uint64_t>(b)) & mask62;
                r = (r >> 2) | (static_cast<uint64_t>(3 - b) << (2 * (kPK - 1)));
                if (++valid < kPK) continue;
                const size_t p = i + 1 - kPK;
                const uint64_t canon = f < r ? f : r;
                const bool w = inWindow(s.size(), p);
                if (pass == 0 && w) {
                    uint64_t* slot = nullptr;
                    if (!table.insert(canon, packPos(c, static_cast<uint32_t>(p), f == canon ? 0 : 1), slot))
                        *slot = kAmbig;
                } else if (pass == 1 && !w) {
                    uint64_t* slot = table.find(canon);
                    if (slot) *slot = kAmbig;
                }
            }
        }
    }

    // ---- anchor reads pairwise; accumulate crossings, controls and links ----------
    const int32_t anch = opt.pairAnch;
    const int32_t ctrlD = opt.pairCtrl;
    struct Local {
        std::vector<uint32_t> ctrl;
        std::vector<uint32_t> out;
        std::unordered_map<uint64_t, std::vector<int32_t>> spans;
        std::unordered_map<uint64_t, size_t> links;
        size_t anchored = 0, spanning = 0, control = 0;
        double bases = 0;
    };
    const size_t pairedReads = reads.pairedReads();
    const size_t nthreads = static_cast<size_t>(threads);
    std::vector<Local> locals(nthreads);
    for (Local& l : locals) { l.ctrl.assign(static_cast<size_t>(n) * 2, 0); l.out.assign(static_cast<size_t>(n) * 2, 0); }
    std::atomic<size_t> next{0};
    const size_t block = 4096;   // even: a pair never straddles two blocks
    struct Placed { uint32_t piece = UINT32_MAX; int32_t pos = 0; uint8_t orient = 0; uint32_t len = 0; };

    auto place = [&](size_t r) {
        Placed pl;
        const uint32_t len = reads.length(r);
        if (len < static_cast<uint32_t>(kPK)) return pl;
        struct Vote { uint32_t piece; int32_t pos; uint8_t orient; int count; };
        Vote votes[kMaxProbes];
        int distinct = 0;
        const uint32_t probes = std::min<uint32_t>(kMaxProbes, len - kPK + 1);
        const uint32_t step = std::max<uint32_t>(1, (len - kPK + 1) / probes);
        for (uint32_t rp = 0; rp + kPK <= len; rp += step) {
            uint64_t f = 0, rc = 0;
            bool bad = false;
            for (int q = 0; q < kPK; ++q) {
                const int b = reads.baseAt(r, rp + static_cast<uint32_t>(q));
                if (b < 0) { bad = true; break; }
                f = ((f << 2) | static_cast<uint64_t>(b)) & mask62;
                rc = (rc >> 2) | (static_cast<uint64_t>(3 - b) << (2 * (kPK - 1)));
            }
            if (bad) continue;
            const uint64_t canon = f < rc ? f : rc;
            const uint64_t* v = table.find(canon);
            if (!v || *v == kAmbig) continue;
            const uint32_t c = static_cast<uint32_t>(*v >> 33);
            const int up = static_cast<int>((*v >> 1) & 0xFFFFFFFFULL);
            const int sflag = static_cast<int>(*v & 1);
            const int orient = ((f == canon) ? 0 : 1) ^ sflag;
            const int32_t startPos = orient == 0
                ? static_cast<int32_t>(up - static_cast<int>(rp))
                : static_cast<int32_t>(up - static_cast<int>(len - rp - static_cast<uint32_t>(kPK)));
            int found = -1;
            for (int q = 0; q < distinct; ++q) {
                if (votes[q].piece == c && votes[q].pos == startPos &&
                    votes[q].orient == static_cast<uint8_t>(orient)) { found = q; break; }
            }
            if (found >= 0) ++votes[found].count;
            else if (distinct < kMaxProbes) votes[distinct++] = {c, startPos, static_cast<uint8_t>(orient), 1};
        }
        int best = -1, bestVotes = 0;
        for (int q = 0; q < distinct; ++q) {
            if (votes[q].count > bestVotes) { bestVotes = votes[q].count; best = q; }
        }
        if (best < 0 || bestVotes < kMinVotes) return pl;
        pl.piece = votes[best].piece;
        pl.pos = votes[best].pos;
        pl.orient = votes[best].orient;
        pl.len = len;
        return pl;
    };

    auto worker = [&](size_t tid) {
        Local& L = locals[tid];
        while (true) {
            const size_t begin = next.fetch_add(block);
            if (begin >= pairedReads) break;
            const size_t end = std::min(begin + block, pairedReads);
            for (size_t r = begin; r + 1 < end; r += 2) {
                const Placed a = place(r);
                const Placed b = place(r + 1);
                if (a.piece != UINT32_MAX) { ++L.anchored; L.bases += a.len; }
                if (b.piece != UINT32_MAX) { ++L.anchored; L.bases += b.len; }
                if (a.piece == UINT32_MAX || b.piece == UINT32_MAX) continue;
                const int32_t la = static_cast<int32_t>(len_[a.piece]);
                const int32_t lb = static_cast<int32_t>(len_[b.piece]);
                const int aEnd = a.orient == 0 ? 1 : 0;
                const int bEnd = b.orient == 0 ? 1 : 0;
                if (a.piece == b.piece) {
                    if (a.orient == b.orient) continue;
                    const Placed& fw = a.orient == 0 ? a : b;
                    const Placed& rv = a.orient == 0 ? b : a;
                    if (fw.pos <= rv.pos) {
                        // an ordinary inward pair: a control crossing when it straddles a
                        // control point
                        const int32_t s0 = fw.pos;
                        const int32_t s1 = rv.pos + static_cast<int32_t>(rv.len);
                        if (s1 - s0 > maxFrag) continue;
                        const int32_t cL = ctrlD, cR = la - ctrlD;
                        if (s0 <= cL - anch && s1 >= cL + anch) { ++L.ctrl[a.piece * 2 + 0]; ++L.control; }
                        if (s0 <= cR - anch && s1 >= cR + anch) { ++L.ctrl[a.piece * 2 + 1]; ++L.control; }
                        continue;
                    }
                    if (aEnd == bEnd) continue;
                }
                const int32_t da = aEnd == 1 ? la - a.pos : a.pos + static_cast<int32_t>(a.len);
                const int32_t db = bEnd == 1 ? lb - b.pos : b.pos + static_cast<int32_t>(b.len);
                if (a.piece != b.piece)
                    ++L.links[pairKey(std::min(a.piece, b.piece), std::max(a.piece, b.piece))];
                if (da < anch || db < anch || da + db > maxFrag) continue;
                const uint32_t pa = a.piece * 2 + static_cast<uint32_t>(aEnd);
                const uint32_t pb = b.piece * 2 + static_cast<uint32_t>(bEnd);
                L.spans[pairKey(pa, pb)].push_back(da + db);
                L.spans[pairKey(pb, pa)].push_back(da + db);
                ++L.out[pa];
                ++L.out[pb];
                ++L.spanning;
            }
        }
    };
    {
        std::vector<std::thread> pool;
        for (size_t t = 0; t < nthreads; ++t) pool.emplace_back(worker, t);
        for (std::thread& t : pool) t.join();
    }
    std::vector<uint32_t> ctrlSum(static_cast<size_t>(n) * 2, 0);
    double bases = 0;
    for (Local& L : locals) {
        for (size_t i = 0; i < ctrlSum.size(); ++i) { ctrlSum[i] += L.ctrl[i]; outTotal_[i] += L.out[i]; }
        for (auto& kv : L.spans) {
            auto& dst = spans_[kv.first];
            dst.insert(dst.end(), kv.second.begin(), kv.second.end());
        }
        for (auto& kv : L.links) links_[kv.first] += kv.second;
        readsAnchored += L.anchored;
        pairsSpanning += L.spanning;
        pairsControl += L.control;
        bases += L.bases;
    }
    // A control exists only where the pilot could take one: the point lies inside the piece
    // with a clean (N-free) stretch between it and the end.
    for (uint32_t c = 0; c < n; ++c) {
        const std::string& s = origin[c];
        const long L = static_cast<long>(s.size());
        if (L <= ctrlD + anch) continue;
        const bool leftClean = s.find('N', 0) == std::string::npos ||
                               static_cast<long>(s.find('N', 0)) > ctrlD;
        const size_t lastN = s.rfind('N');
        const bool rightClean = lastN == std::string::npos ||
                                static_cast<long>(lastN) < L - ctrlD - 1;
        if (leftClean) ctrl_[c * 2 + 0] = static_cast<int32_t>(ctrlSum[c * 2 + 0]);
        if (rightClean) ctrl_[c * 2 + 1] = static_cast<int32_t>(ctrlSum[c * 2 + 1]);
    }
    size_t windowBases = 0;
    for (uint32_t c = 0; c < n; ++c) windowBases += std::min<size_t>(len_[c], static_cast<size_t>(2 * W));
    meanDepth_ = windowBases ? bases / static_cast<double>(windowBases) : 0;
    usable_ = true;
}

uint32_t PairIndex::exitPort(const std::string& s) const { return exitLookup(portOf_, s); }
uint32_t PairIndex::entryPort(const std::string& s) const { return entryLookup(portOf_, s); }

size_t PairIndex::linkWeight(uint32_t a, uint32_t b) const {
    auto it = links_.find(pairKey(std::min(a, b), std::max(a, b)));
    return it == links_.end() ? 0 : it->second;
}

std::vector<std::pair<uint32_t, uint32_t>> PairIndex::partners(uint32_t port) const {
    std::vector<std::pair<uint32_t, uint32_t>> out;
    if (!usable_ || port == UINT32_MAX) return out;
    // spans_ holds both orders of every crossing pair, so the keys that start with `port` are all
    // of its partners. A linear pass: the clonal stage asks for a few hundred ports per run.
    for (const auto& kv : spans_) {
        if (static_cast<uint32_t>(kv.first >> 32) != port) continue;
        out.emplace_back(static_cast<uint32_t>(kv.first & 0xFFFFFFFFULL), static_cast<uint32_t>(kv.second.size()));
    }
    std::sort(out.begin(), out.end());
    return out;
}

size_t PairIndex::supportFloor() const {
    return static_cast<size_t>(std::min(10.0, std::max(3.0, meanDepth_ * 0.06)));
}

PairCall PairIndex::callFrom(uint32_t k, float lambda, const EvidenceOptions& opt) {
    PairCall pc;
    pc.testable = true;
    pc.k = k;
    pc.lambda = lambda;
    if (!(lambda >= opt.pairLambdaMin)) { pc.abstain = true; return pc; }
    pc.abstain = false;
    // pilot_summary.py: under adjacency expect 0.5 * lambda crossing pairs (seam efficiency),
    // under an insertion 0.5.
    const double le = 0.5 * static_cast<double>(lambda), li = 0.5;
    const double kk = static_cast<double>(k);
    pc.log10LR = static_cast<float>((kk * std::log(li) - li - (kk * std::log(le) - le)) / std::log(10.0));
    if (pc.log10LR >= opt.pairLog10LR) pc.call = +1;
    else if (pc.log10LR <= -opt.pairLog10LR) pc.call = -1;
    return pc;
}

PairCall PairIndex::test(uint32_t P, uint32_t Q, int32_t claimedN) const {
    PairCall pc;
    if (!usable_ || P == UINT32_MAX || Q == UINT32_MAX || claimedN > opt_.pairMaxN) return pc;
    if (P >= ctrl_.size() || Q >= ctrl_.size()) return pc;
    uint32_t k = 0;
    auto it = spans_.find(pairKey(P, Q));
    if (it != spans_.end()) k = static_cast<uint32_t>(it->second.size());
    const int32_t cA = ctrl_[P], cB = ctrl_[Q];
    double sum = 0;
    int nc = 0;
    if (cA >= 0) { sum += cA; ++nc; }
    if (cB >= 0) { sum += cB; ++nc; }
    if (nc == 0) {
        pc.testable = true;
        pc.k = k;
        pc.abstain = true;
        return pc;
    }
    pc = callFrom(k, static_cast<float>(sum / nc), opt_);
    pc.ctrlA = static_cast<float>(cA);
    pc.ctrlB = static_cast<float>(cB);
    pc.contra = outTotal_[P] >= k ? outTotal_[P] - k : 0;
    return pc;
}

}  // namespace om2
}  // namespace ts
