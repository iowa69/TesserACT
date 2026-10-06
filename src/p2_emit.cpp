// Phase 2 W1 emit-A: R3 exact / verified circles and F5 spike-in labels on the written records,
// plus the counter lines and the report.json "p2" block of R2, R3 and F5. See p2_emit.h.
#include "p2_emit.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <set>
#include <thread>
#include <unordered_set>
#include <utility>

#include "envflags.h"
#include "kmer.h"

namespace ts {
namespace p2 {

EmitConfig readConfig() {
    EmitConfig c;
    c.libGuard = env::on("TESSERACT_P2_LIBGUARD", false);
    if (const char* m = env::text("TESSERACT_P2_CIRC")) {
        if (std::strcmp(m, "close") == 0) c.circ = CircMode::Close;
        else if (std::strcmp(m, "verify") == 0) c.circ = CircMode::Verify;
    }
    c.spikein = env::on("TESSERACT_P2_SPIKEIN", false);
    return c;
}

const char* circModeName(CircMode m) {
    switch (m) {
        case CircMode::Off: return "off";
        case CircMode::Close: return "close";
        case CircMode::Verify: return "verify";
    }
    return "?";
}

const char* closureName(Closure c) {
    switch (c) {
        case Closure::NoPath: return "nopath";
        case Closure::EndsChanged: return "ends_changed";
        case Closure::RepeatEnd: return "repeat_end";
        case Closure::Self: return "self";
        case Closure::Via: return "via";
        case Closure::None: return "none";
        case Closure::Gapped: return "gapped";
    }
    return "?";
}

std::string dropCircular(const std::string& tag) {
    std::string t = tag;
    const size_t c = t.find("_circular");
    if (c != std::string::npos) t.erase(c, std::strlen("_circular"));
    return t;
}

// ---- graph helpers ---------------------------------------------------------------------------
namespace {

inline bool liveNode(const UnitigGraph& g, uint64_t o) {
    const uint64_t u = o >> 1;
    return u < g.nodes.size() && !g.nodes[u].deleted;
}

std::string orientedSeq(const UnitigGraph& g, uint64_t o) {
    return g.oriented(static_cast<uint32_t>(o >> 1), static_cast<int>(o & 1));
}

// Oriented successors of `o` ((unitig << 1) | orient, 0 = forward), as the GFA L lines list them:
// leaving a forward unitig through its 3' end, entering {v, 0} forward and {v, 1} reversed.
void successors(const UnitigGraph& g, uint64_t o, std::vector<uint64_t>& out) {
    out.clear();
    if (!liveNode(g, o)) return;
    const uint32_t u = static_cast<uint32_t>(o >> 1);
    for (const Link& l : g.exits(u, static_cast<int>(o & 1))) {
        if (l.to >= g.nodes.size() || g.nodes[l.to].deleted) continue;
        out.push_back((static_cast<uint64_t>(l.to) << 1) | static_cast<uint64_t>(l.toEnd & 1));
    }
}

void predecessors(const UnitigGraph& g, uint64_t o, std::vector<uint64_t>& out) {
    successors(g, o ^ 1ULL, out);
    for (uint64_t& x : out) x ^= 1ULL;
}

std::string segName(uint64_t o) { return std::to_string(o >> 1); }
char segSign(uint64_t o) { return (o & 1) ? '-' : '+'; }

// The coverage as the GFA prints it (dp:f:%.4f), which is what the prototype compared.
double dp4(const UnitigGraph& g, uint64_t o) {
    char buf[64];
    std::snprintf(buf, sizeof buf, "%.4f", g.nodes[o >> 1].coverage);
    return std::strtod(buf, nullptr);
}

}  // namespace

ClosureResult closeWalk(const std::string& s, const std::vector<uint64_t>& walk, const UnitigGraph& g) {
    ClosureResult r;
    r.pathN = walk.size();
    if (walk.empty()) { r.mode = Closure::NoPath; return r; }
    const uint64_t first = walk.front(), last = walk.back();
    if (!liveNode(g, first) || !liveNode(g, last) || g.k() < 2) {
        r.mode = Closure::NoPath;
        r.note = "walk_off_graph";
        return r;
    }
    const size_t k1 = static_cast<size_t>(g.k() - 1);
    const std::string fseq = orientedSeq(g, first), lseq = orientedSeq(g, last);
    const bool endsOk = s.size() >= k1 && fseq.size() >= k1 && lseq.size() >= k1 &&
                        s.compare(0, k1, fseq, 0, k1) == 0 &&
                        s.compare(s.size() - k1, k1, lseq, lseq.size() - k1, k1) == 0;

    std::vector<uint64_t> outs, tmp;
    successors(g, last, outs);
    bool direct = false;
    std::vector<uint64_t> vias;
    for (uint64_t b : outs) {
        if (b == first) { direct = true; continue; }
        successors(g, b, tmp);
        for (uint64_t c : tmp)
            if (c == first) vias.push_back(b);
    }
    const std::set<uint64_t> outSet(outs.begin(), outs.end());
    std::vector<uint64_t> ins;
    predecessors(g, first, ins);
    const std::set<uint64_t> inSet(ins.begin(), ins.end());
    r.outLast = static_cast<int>(outSet.size());
    r.inFirst = static_cast<int>(inSet.size());
    // Isolated closure: every exit of the last segment re-enters the first one, directly or
    // through one connector, and every entry of the first segment comes from the last one the
    // same way (Unicycler's completed-circle rule extended to a bubble at the join).
    auto leadsToFirst = [&](uint64_t b) {
        if (b == first) return true;
        std::vector<uint64_t> nx;
        successors(g, b, nx);
        return std::find(nx.begin(), nx.end(), first) != nx.end();
    };
    auto comesFromLast = [&](uint64_t a) { return a == last || outSet.count(a) > 0; };
    bool iso = !outSet.empty() && !inSet.empty();
    if (iso) for (uint64_t b : outSet) if (!leadsToFirst(b)) { iso = false; break; }
    if (iso) for (uint64_t a : inSet) if (!comesFromLast(a)) { iso = false; break; }
    r.isolated = iso ? 1 : 0;

    if (!endsOk) {
        r.mode = Closure::EndsChanged;
    } else if (walk.size() > 1 && first == last) {
        // The walk re-enters its first segment: the circle is the walk minus that repeated
        // segment (its new bases plus the closing overlap).
        r.mode = Closure::RepeatEnd;
        const size_t cut = fseq.size();
        if (s.size() > cut && s.compare(s.size() - cut, cut, s, 0, cut) == 0) {
            r.closed = true;
            r.seq = s.substr(0, s.size() - cut);
        } else {
            r.note = "repeat_end_mismatch";
        }
    } else if (direct) {
        r.mode = Closure::Self;
        if (s.size() > k1 && s.compare(s.size() - k1, k1, s, 0, k1) == 0) {
            r.closed = true;
            r.seq = s.substr(0, s.size() - k1);
        } else {
            r.note = "overlap_mismatch";
        }
    } else if (!vias.empty()) {
        r.mode = Closure::Via;
        const std::set<uint64_t> vs(vias.begin(), vias.end());
        r.alts = vs.size();
        // The deepest connector; ties as the prototype breaks them (segment name as text, sign).
        uint64_t y = *vs.begin();
        auto key = [&](uint64_t o) { return std::make_tuple(dp4(g, o), segName(o), segSign(o)); };
        for (uint64_t o : vs)
            if (key(o) > key(y)) y = o;
        const std::string ys = orientedSeq(g, y);
        const long interior = static_cast<long>(ys.size()) - 2 * static_cast<long>(k1);
        r.viaSeg = segName(y) + segSign(y);
        r.viaLen = ys.size();
        r.viaDp = dp4(g, y);
        std::vector<uint64_t> others(vs.begin(), vs.end());
        std::sort(others.begin(), others.end(), [](uint64_t a, uint64_t b) {
            return std::make_pair(segName(a), segSign(a)) < std::make_pair(segName(b), segSign(b));
        });
        for (uint64_t o : others) {
            if (o == y) continue;
            char buf[96];
            std::snprintf(buf, sizeof buf, "%s%s%c:%.1f", r.altDp.empty() ? "" : ",", segName(o).c_str(),
                          segSign(o), dp4(g, o));
            r.altDp += buf;
        }
        r.interior = interior;
        if (interior >= 0) {
            r.closed = true;
            r.seq = s + ys.substr(k1, ys.size() - 2 * k1);
        } else {
            const size_t o = static_cast<size_t>(-interior);
            if (s.size() > o && s.compare(s.size() - o, o, s, 0, o) == 0) {
                r.closed = true;
                r.seq = s.substr(0, s.size() - o);
            } else {
                r.note = "short_via_mismatch";
            }
        }
        std::vector<uint64_t> yo, yi;
        successors(g, y, yo);
        predecessors(g, y, yi);
        std::set<uint64_t> yOut(yo.begin(), yo.end()), yIn(yi.begin(), yi.end());
        yOut.erase(first);
        yIn.erase(last);
        r.viaOtherLinks = yOut.size() + yIn.size();
    } else {
        r.mode = Closure::None;
    }
    return r;
}

size_t terminalSelfOverlap(const std::string& s, size_t lo, size_t hi) {
    if (s.size() < 2) return 0;
    size_t best = 0;
    const size_t top = std::min(hi, s.size() - 1);
    for (size_t o = lo; o <= top; ++o)
        if (s.compare(0, o, s, s.size() - o, o) == 0) best = o;
    return best;
}

size_t kMinusOneSelfOverlap(const std::string& s, const std::vector<int>& ladder) {
    size_t best = 0;
    for (int k : ladder) {
        if (k < 2) continue;
        const size_t o = static_cast<size_t>(k - 1);
        if (s.size() <= 2 * o) continue;
        if (o > best && s.compare(0, o, s, s.size() - o, o) == 0) best = o;
    }
    return best;
}

// ---- 31-mers in 64 bits ------------------------------------------------------------------------
namespace {

constexpr uint64_t kKmerMask = (uint64_t(1) << (2 * kAnchorK)) - 1;
constexpr int kTopShift = 2 * (kAnchorK - 1);

// fn(canonical, position of the first base, strand: 0 when the forward k-mer is canonical)
template <class F>
void forEach31(const std::string& s, F fn) {
    uint64_t f = 0, r = 0;
    int valid = 0;
    for (size_t i = 0; i < s.size(); ++i) {
        const int b = baseCode(s[i]);
        if (b < 0) { valid = 0; f = r = 0; continue; }
        f = ((f << 2) | static_cast<uint64_t>(b)) & kKmerMask;
        r = (r >> 2) | (static_cast<uint64_t>(3 - b) << kTopShift);
        if (++valid >= kAnchorK) fn(f <= r ? f : r, i + 1 - static_cast<size_t>(kAnchorK), f <= r ? 0 : 1);
    }
}

const std::unordered_set<uint64_t>& phixSet() {
    static const std::unordered_set<uint64_t> kSet = [] {
        std::unordered_set<uint64_t> s;
        const std::string& p = phix174();
        const std::string circ = p + p.substr(0, static_cast<size_t>(kAnchorK - 1));
        forEach31(circ, [&](uint64_t c, size_t, int) { s.insert(c); });
        return s;
    }();
    return kSet;
}

}  // namespace

double phixKmerFraction(const std::string& s) {
    const auto& set = phixSet();
    size_t n = 0, hit = 0;
    forEach31(s, [&](uint64_t c, size_t, int) {
        ++n;
        if (set.count(c)) ++hit;
    });
    return n ? static_cast<double>(hit) / static_cast<double>(n) : 0.0;
}

// ---- R3 junction verification -------------------------------------------------------------------
namespace {

struct Placed {
    uint32_t rec = UINT32_MAX;
    int32_t pos = 0;
    uint8_t orient = 0;
    bool ok() const { return rec != UINT32_MAX; }
};

// Every occurrence of every 31-mer of the records, sorted. A k-mer present more than kMaxOcc
// times is ignored (it says nothing about where a read belongs).
constexpr size_t kMaxOcc = 16;
struct Index31 {
    std::vector<uint64_t> keys, vals;
    std::pair<size_t, size_t> find(uint64_t k) const {
        const auto lo = std::lower_bound(keys.begin(), keys.end(), k);
        if (lo == keys.end() || *lo != k) return {0, 0};
        const auto hi = std::upper_bound(lo, keys.end(), k);
        return {static_cast<size_t>(lo - keys.begin()), static_cast<size_t>(hi - keys.begin())};
    }
};

Index31 buildIndex(const std::vector<const std::string*>& recs) {
    std::vector<std::pair<uint64_t, uint64_t>> all;
    size_t total = 0;
    for (const std::string* s : recs) total += s->size();
    all.reserve(total);
    for (size_t r = 0; r < recs.size(); ++r) {
        forEach31(*recs[r], [&](uint64_t c, size_t pos, int strand) {
            all.emplace_back(c, (static_cast<uint64_t>(r) << 33) | (static_cast<uint64_t>(pos) << 1) |
                                    static_cast<uint64_t>(strand));
        });
    }
    std::sort(all.begin(), all.end());
    Index31 ix;
    ix.keys.reserve(all.size());
    ix.vals.reserve(all.size());
    for (const auto& kv : all) {
        ix.keys.push_back(kv.first);
        ix.vals.push_back(kv.second);
    }
    return ix;
}

uint64_t mix64(uint64_t x) {   // splitmix64 finaliser: a fixed pseudo-random choice per read
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}

// Placement by k-mer votes: up to 12 probes, each voting for every occurrence of its k-mer; the
// (record, start, strand) placements with the most votes (at least 2) are returned, at most
// kMaxTies of them. One placement = a unique placement; several = a read inside sequence the
// records share, which an aligner places on one copy at random (MAPQ 0) and so does the caller.
constexpr size_t kMaxTies = 8;
std::vector<Placed> placeCandidates(const SequenceStore& reads, size_t r, const Index31& ix) {
    std::vector<Placed> out;
    const uint32_t len = reads.length(r);
    if (len < static_cast<uint32_t>(kAnchorK)) return out;
    struct Vote { uint32_t rec; int32_t pos; uint8_t orient; int count; };
    std::vector<Vote> votes;
    const uint32_t span = len - static_cast<uint32_t>(kAnchorK) + 1;
    const uint32_t probes = std::min<uint32_t>(12, span);
    const uint32_t step = std::max<uint32_t>(1, span / probes);
    for (uint32_t rp = 0; rp + static_cast<uint32_t>(kAnchorK) <= len; rp += step) {
        uint64_t f = 0, rc = 0;
        for (int q = 0; q < kAnchorK; ++q) {
            const int b = reads.rawBaseAt(r, rp + static_cast<uint32_t>(q));
            f = ((f << 2) | static_cast<uint64_t>(b)) & kKmerMask;
            rc = (rc >> 2) | (static_cast<uint64_t>(3 - b) << kTopShift);
        }
        const bool fwdCanon = f <= rc;
        const std::pair<size_t, size_t> rg = ix.find(fwdCanon ? f : rc);
        if (rg.second == rg.first || rg.second - rg.first > kMaxOcc) continue;
        for (size_t x = rg.first; x < rg.second; ++x) {
            const uint64_t v = ix.vals[x];
            const uint32_t rec = static_cast<uint32_t>(v >> 33);
            const int up = static_cast<int>((v >> 1) & 0xFFFFFFFFULL);
            const int orient = (fwdCanon ? 0 : 1) ^ static_cast<int>(v & 1);
            const int32_t start = orient == 0 ? up - static_cast<int>(rp)
                                              : up - static_cast<int>(len - rp - static_cast<uint32_t>(kAnchorK));
            bool found = false;
            for (Vote& w : votes)
                if (w.rec == rec && w.pos == start && w.orient == orient) { ++w.count; found = true; break; }
            if (!found) votes.push_back({rec, start, static_cast<uint8_t>(orient), 1});
        }
    }
    int bv = 0;
    for (const Vote& w : votes) bv = std::max(bv, w.count);
    if (bv < 2) return out;
    for (const Vote& w : votes) {
        if (w.count != bv) continue;
        Placed p;
        p.rec = w.rec;
        p.pos = w.pos;
        p.orient = w.orient;
        out.push_back(p);
        if (out.size() >= kMaxTies) break;
    }
    return out;
}

// One placement per read. A read with tied placements goes where its mate makes a proper pair
// (same record, opposite strands, within 5 kb) when it can, else to a fixed pseudo-random one.
void placeAll(const SequenceStore& reads, const Index31& ix, int threads, std::vector<Placed>& placed,
              size_t& nplacedOut) {
    const size_t n = reads.size();
    placed.assign(n, Placed());
    const size_t paired = reads.pairedReads();
    std::atomic<size_t> next{0}, nplaced{0};
    auto proper = [](const Placed& a, const Placed& b) {
        return a.rec == b.rec && a.orient != b.orient && std::labs(static_cast<long>(a.pos) - b.pos) <= 5000;
    };
    auto worker = [&]() {
        size_t local = 0;
        while (true) {
            const size_t b0 = next.fetch_add(4096);
            if (b0 >= n) break;
            const size_t e0 = std::min(n, b0 + 4096);
            for (size_t r = b0; r < e0; ++r) {
                if (r < paired && (r & 1)) continue;          // handled with its first mate
                if (r < paired) {
                    const std::vector<Placed> A = placeCandidates(reads, r, ix);
                    const std::vector<Placed> B = placeCandidates(reads, r + 1, ix);
                    Placed pa, pb;
                    if (A.size() > 1 || B.size() > 1) {
                        std::vector<std::pair<size_t, size_t>> good;
                        for (size_t i = 0; i < A.size(); ++i)
                            for (size_t j = 0; j < B.size(); ++j)
                                if (proper(A[i], B[j])) good.emplace_back(i, j);
                        if (!good.empty()) {
                            const auto& g = good[mix64(r) % good.size()];
                            pa = A[g.first];
                            pb = B[g.second];
                        } else {
                            if (!A.empty()) pa = A[mix64(r) % A.size()];
                            if (!B.empty()) pb = B[mix64(r + 1) % B.size()];
                        }
                    } else {
                        if (!A.empty()) pa = A[0];
                        if (!B.empty()) pb = B[0];
                    }
                    placed[r] = pa;
                    placed[r + 1] = pb;
                    local += (pa.ok() ? 1 : 0) + (pb.ok() ? 1 : 0);
                } else {
                    const std::vector<Placed> A = placeCandidates(reads, r, ix);
                    if (!A.empty()) { placed[r] = A[mix64(r) % A.size()]; ++local; }
                }
            }
        }
        nplaced += local;
    };
    std::vector<std::thread> pool;
    for (int t = 0; t < std::max(1, threads); ++t) pool.emplace_back(worker);
    for (auto& t : pool) t.join();
    nplacedOut = nplaced.load();
}

// Gapless local alignment of a placed read against its record: match +2, mismatch -8 (the
// minimap2 short-read scores), an ambiguous record base -1, a base beyond the record's ends a
// mismatch. [a, b) is the best-scoring stretch of the read; the rest is "soft-clipped".
struct Aln { int a = 0, b = 0, score = 0; };
Aln alignGapless(const SequenceStore& reads, size_t r, const Placed& p, const std::string& rec) {
    const int len = static_cast<int>(reads.length(r));
    Aln best;
    int cur = 0, curStart = 0;
    best.score = 0;
    for (int i = 0; i < len; ++i) {
        const long rp = static_cast<long>(p.pos) + i;
        int sc;
        if (rp < 0 || rp >= static_cast<long>(rec.size())) {
            sc = -8;
        } else {
            const int rb = baseCode(rec[static_cast<size_t>(rp)]);
            const int qb = p.orient == 0 ? reads.rawBaseAt(r, static_cast<uint32_t>(i))
                                         : 3 - reads.rawBaseAt(r, static_cast<uint32_t>(len - 1 - i));
            sc = rb < 0 ? -1 : (rb == qb ? 2 : -8);
        }
        if (cur <= 0) { cur = 0; curStart = i; }
        cur += sc;
        if (cur > best.score) { best.score = cur; best.a = curStart; best.b = i + 1; }
    }
    return best;
}

constexpr int kMinAlnScore = 40;   // minimap2 -x sr minimal peak DP score (-s40)

}  // namespace

bool passJ(const JunctionStats& j, const VerifyThresholds& t) {
    return j.computed && j.spanReads >= t.minSpan && j.clipFrac <= t.maxClipFrac &&
           static_cast<double>(j.spanReads) >= t.minSpanCtrl * static_cast<double>(std::max<size_t>(1, j.ctrlSpan));
}

bool passX(const JunctionStats& j, const VerifyThresholds& t) {
    return j.computed && j.exitRatio <= t.maxExitRatio;
}

std::vector<JunctionStats> verifyJunctions(const VerifyBatch& batch, const SequenceStore& reads,
                                           int threads, VerifyRunStats& run) {
    std::vector<JunctionStats> out(batch.candRecord.size());
    run = VerifyRunStats();
    if (!batch.records || batch.candRecord.empty() || reads.size() == 0) return out;
    const std::vector<std::string>& records = *batch.records;
    std::vector<const std::string*> recs(records.size());
    for (size_t i = 0; i < records.size(); ++i) recs[i] = &records[i];
    std::vector<int> candOf(records.size(), -1);
    for (size_t c = 0; c < batch.candRecord.size(); ++c) {
        recs[batch.candRecord[c]] = &batch.candCircle[c];
        candOf[batch.candRecord[c]] = static_cast<int>(c);
    }
    const Index31 ix = buildIndex(recs);
    run.indexKmers = ix.keys.size();

    // ---- place every read --------------------------------------------------------------------
    const size_t n = reads.size();
    std::vector<Placed> placed;
    placeAll(reads, ix, threads, placed, run.readsPlaced);

    // ---- library window: FR template lengths on the longest record (junction_verify.py) -------
    size_t longest = 0;
    for (size_t i = 1; i < recs.size(); ++i)
        if (recs[i]->size() > recs[longest]->size()) longest = i;
    std::vector<int> tl;
    for (size_t r = 0; r + 1 < reads.pairedReads() && tl.size() < 200000; r += 2) {
        const Placed& x = placed[r];
        const Placed& y = placed[r + 1];
        if (!x.ok() || !y.ok() || x.rec != longest || y.rec != longest || x.orient == y.orient) continue;
        const Placed& f = x.orient == 0 ? x : y;
        const Placed& v = x.orient == 0 ? y : x;
        const int vlen = static_cast<int>(reads.length(x.orient == 0 ? r + 1 : r));
        if (f.pos > v.pos) continue;
        const int t = v.pos + vlen - f.pos;
        if (t > 0) tl.push_back(t);
    }
    int lo = 50, hi = 1000;
    if (!tl.empty()) {
        std::sort(tl.begin(), tl.end());
        const int med = tl[tl.size() / 2];
        std::vector<int> kept;
        for (int t : tl) if (t <= 10 * med) kept.push_back(t);
        if (!kept.empty()) {
            lo = kept[static_cast<size_t>(static_cast<double>(kept.size()) * 0.005)];
            hi = kept[std::min(kept.size() - 1, static_cast<size_t>(static_cast<double>(kept.size()) * 0.995))];
        }
        run.windowPairs = kept.size();
    }
    run.windowLo = lo;
    run.windowHi = hi;

    // ---- reads on each candidate circle ------------------------------------------------------
    std::vector<std::vector<size_t>> readsOf(batch.candRecord.size());
    for (size_t r = 0; r < n; ++r)
        if (placed[r].ok() && candOf[placed[r].rec] >= 0) readsOf[static_cast<size_t>(candOf[placed[r].rec])].push_back(r);

    for (size_t c = 0; c < batch.candRecord.size(); ++c) {
        JunctionStats& js = out[c];
        const std::string& R = batch.candCircle[c];
        const long T = static_cast<long>(R.size());
        js.computed = true;
        js.circleLen = R.size();
        const long h = T / 2;
        const long J = T - h;
        const long C = T / 4;
        js.join = static_cast<size_t>(J);
        js.ctrl = static_cast<size_t>(C);
        const long A = kJoinAnchor;
        const uint32_t rec = static_cast<uint32_t>(batch.candRecord[c]);
        for (size_t r : readsOf[c]) {
            const Placed& p = placed[r];
            const Aln al = alignGapless(reads, r, p, R);
            if (al.score < kMinAlnScore) continue;
            const long len = static_cast<long>(reads.length(r));
            const long rs = static_cast<long>(p.pos) + al.a, re = static_cast<long>(p.pos) + al.b;
            const long leftClip = al.a, rightClip = len - al.b;
            for (int which = 0; which < 2; ++which) {
                const long X = which == 0 ? J : C;
                const bool clipped = (leftClip >= kMinClip && std::labs(rs - X) <= kClipWindow) ||
                                     (rightClip >= kMinClip && std::labs(re - X) <= kClipWindow);
                size_t& span = which == 0 ? js.spanReads : js.ctrlSpan;
                size_t& clip = which == 0 ? js.clipReads : js.ctrlClip;
                if (clipped) ++clip;
                else if (rs <= X - A && re >= X + A) ++span;
            }
            ++js.reads;
            const bool near = std::labs(rs - J) <= hi || std::labs(re - J) <= hi;
            if (near) ++js.nearJoin;
            if (reads.hasMate(r)) {
                const Placed& m = placed[reads.mateOf(r)];
                if (m.ok() && m.rec != rec) {
                    ++js.mateElsewhere;
                    if (near) ++js.exitsNearJoin;
                }
            }
        }
        js.clipFrac = (js.spanReads + js.clipReads)
                          ? static_cast<double>(js.clipReads) / static_cast<double>(js.spanReads + js.clipReads)
                          : 0.0;
        js.exitRatio = js.nearJoin ? static_cast<double>(js.exitsNearJoin) / static_cast<double>(js.nearJoin) : 0.0;
    }
    return out;
}

// ---- the record stage ---------------------------------------------------------------------------
namespace {

bool hasCircular(const std::string& tag) { return tag.find("_circular") != std::string::npos; }

// The piece's own gap-free walk, exactly as contigPathsV2 finds it (empty when none spells it).
std::vector<uint64_t> walkOfPiece(const EmitPiece& p, const std::vector<EmitPiece>& pieces,
                                  const std::vector<std::string>& scaffolds, const std::vector<GfaPath>& paths,
                                  const UnitigGraph& g) {
    const size_t s = p.scaffold;
    if (s >= paths.size() || s >= scaffolds.size() || paths[s].oriented.empty()) return {};
    size_t inScaffold = 0;
    for (const EmitPiece& q : pieces) inScaffold += q.scaffold == s ? 1 : 0;
    const std::vector<PathSegment> segs = pathSegments(paths[s], g);
    const std::vector<size_t> pieceOf = mapSegmentsToPieces(scaffolds[s], segs, g);
    if (pieceOf.empty() || pieceOf.back() + 1 != inScaffold) return {};
    size_t nSeg = 0, seg = 0;
    for (size_t x = 0; x < pieceOf.size(); ++x)
        if (pieceOf[x] == p.ordinal) { ++nSeg; seg = x; }
    if (nSeg != 1 || spellWalk(segs[seg].oriented, g) != p.seq) return {};
    return segs[seg].oriented;
}

void recordEdit(RecordStageResult& res, size_t piece, const char* feature, const std::string& op,
                size_t start, size_t end, size_t bases) {
    EditRow e;
    e.piece = piece;
    e.feature = feature;
    e.operation = op;
    e.start = start;
    e.end = end;
    e.bases = bases;
    res.edits.push_back(std::move(e));
}

// Applies `next` to piece i (sequence edit at the 3' end only) and records it.
void applySeqEdit(RecordStageIn& in, RecordStageResult& res, RecordRow& row, size_t i, const std::string& next,
                  const std::string& why) {
    EmitPiece& p = (*in.pieces)[i];
    const std::string prev = p.seq;
    if (next == prev) return;
    if (next.size() < prev.size() && prev.compare(0, next.size(), next) == 0) {
        const size_t cut = prev.size() - next.size();
        recordEdit(res, i, "R3", "trim_end:" + why, next.size(), prev.size(), cut);
        row.trimBp += cut;
        // The removed bases are a duplicate of the record's own start. Where the piece is not a
        // whole untrimmed scaffold the scaffold keeps them and the AGP tiles them from the record
        // that holds them (writeAgpV2), like the bases the terminal trim removes.
        p.cutBackSeq = prev.substr(next.size()) + p.cutBackSeq;
    } else if (next.size() > prev.size() && next.compare(0, prev.size(), prev) == 0) {
        const size_t add = next.size() - prev.size();
        recordEdit(res, i, "R3", "insert_end:" + why, prev.size(), prev.size(), add);
        row.insertBp += add;
    } else {
        recordEdit(res, i, "R3", "replace:" + why, 0, prev.size(), next.size());
    }
    p.seq = next;
}

}  // namespace

RecordStageResult applyRecordStage(RecordStageIn& in) {
    RecordStageResult res;
    std::vector<EmitPiece>& P = *in.pieces;
    std::vector<size_t> live;
    for (size_t i = 0; i < P.size(); ++i)
        if (!P[i].drop) live.push_back(i);
    res.rows.resize(live.size());
    for (size_t j = 0; j < live.size(); ++j) {
        RecordRow& row = res.rows[j];
        row.piece = live[j];
        row.oldTag = P[live[j]].tag;
        row.lenBefore = P[live[j]].seq.size();
        row.topology = hasCircular(row.oldTag) ? "circular" : "linear";
        row.verdict = "-";
        row.edit = "none";
    }
    // Snapshot for the scaffold propagation below: which pieces are whole, untrimmed scaffolds.
    std::vector<char> wholeScaffold(P.size(), 0);
    for (size_t i = 0; i < P.size(); ++i) {
        const EmitPiece& p = P[i];
        wholeScaffold[i] = p.piecesInScaffold == 1 && p.cutFrontSeq.empty() && p.cutBackSeq.empty() &&
                           p.scaffold < in.scaffolds->size() && (*in.scaffolds)[p.scaffold] == p.seq;
    }

    // ---- F5: phiX174 spike-in screen ---------------------------------------------------------
    if (in.cfg.spikein) {
        for (RecordRow& row : res.rows) {
            EmitPiece& p = P[row.piece];
            row.phixFrac = phixKmerFraction(p.seq);
            ++res.st.screened;
            if (p.seq.size() < static_cast<size_t>(kAnchorK) || row.phixFrac < kSpikeinMinFrac) continue;
            row.spikein = true;
            ++res.st.spikein;
            res.st.spikeinBp += p.seq.size();
            if (hasCircular(p.tag)) ++res.st.spikeinCircular;
            recordEdit(res, row.piece, "F5", "rename:" + p.tag + ">_spikein", 0, 0, 0);
            p.tag = "_spikein";
            row.topology = "linear";
            row.verdict = "spikein";
        }
    }

    // ---- R3: closure, verification, duplicate-free linear fallback ------------------------------
    if (in.cfg.circ != CircMode::Off && in.graph) {
        const UnitigGraph& g = *in.graph;
        VerifyBatch batch;
        std::vector<std::string> records;
        records.reserve(live.size());
        for (size_t i : live) records.push_back(P[i].seq);
        std::vector<size_t> candRows;
        for (size_t j = 0; j < res.rows.size(); ++j) {
            RecordRow& row = res.rows[j];
            const EmitPiece& p = P[row.piece];
            if (row.spikein || !hasCircular(p.tag)) continue;
            row.candidate = true;
            ++res.st.candidates;
            if (row.phixFrac < 0) row.phixFrac = phixKmerFraction(p.seq);
            if (p.piecesInScaffold > 1) {
                row.closure.mode = Closure::Gapped;
            } else {
                const std::vector<uint64_t> walk = walkOfPiece(p, P, *in.scaffolds, *in.scaffoldPaths, g);
                row.closure = closeWalk(p.seq, walk, g);
            }
            ++res.st.byMode[static_cast<int>(row.closure.mode)];
            if (row.closure.isolated == 1) ++res.st.isolated;
            // The circle the reads are asked about: the graph closure, else the exact terminal
            // self-overlap removed, else the record as emitted (prep_jv_inputs.py).
            std::string T;
            if (row.closure.closed) {
                T = row.closure.seq;
            } else {
                const size_t ov = p.seq.size() >= 40 ? terminalSelfOverlap(p.seq) : 0;
                row.junction.selfOvTrimmed = ov;
                T = p.seq.substr(0, p.seq.size() - ov);
            }
            const size_t h = T.size() / 2;
            batch.candRecord.push_back(j);
            batch.candCircle.push_back(T.substr(h) + T.substr(0, h));
            candRows.push_back(j);
        }
        batch.records = &records;
        if (!candRows.empty() && in.reads) {
            const std::vector<JunctionStats> js = verifyJunctions(batch, *in.reads, in.threads, res.st.verify);
            for (size_t c = 0; c < candRows.size(); ++c) {
                const size_t keepOv = res.rows[candRows[c]].junction.selfOvTrimmed;
                res.rows[candRows[c]].junction = js[c];
                res.rows[candRows[c]].junction.selfOvTrimmed = keepOv;
            }
        }
        for (size_t j : candRows) {
            RecordRow& row = res.rows[j];
            EmitPiece& p = P[row.piece];
            const std::string before = p.seq;
            row.passG = isClosure(row.closure.mode) && row.closure.closed;
            row.passI = row.closure.isolated == 1;
            row.passJ = passJ(row.junction, in.thresholds);
            row.passX = passX(row.junction, in.thresholds);
            if (row.passJ) ++res.st.passJ;
            if (row.passX) ++res.st.passX;
            std::string next = before;
            std::string why;
            bool keepCircular = true;   // close mode keeps the release tag rule
            if (in.cfg.circ == CircMode::Verify) {
                const bool verified = row.passG && row.passI && row.passJ && row.passX;
                if (verified) {
                    ++res.st.verified;
                    next = row.closure.seq;
                    why = std::string("closure_") + closureName(row.closure.mode);
                    row.verdict = "verified";
                    row.edit = "closed";
                } else {
                    keepCircular = false;
                    std::string failed;
                    if (!row.passG) failed += "G";
                    if (!row.passI) failed += failed.empty() ? "I" : ",I";
                    if (!row.passJ) failed += failed.empty() ? "J" : ",J";
                    if (!row.passX) failed += failed.empty() ? "X" : ",X";
                    row.verdict = "linear_fails_" + failed;
                    // D9: written linear and duplicate-free. The graph proves a self / repeat_end
                    // closure's removed bases are written twice; otherwise only an exact (k-1)
                    // terminal self-overlap of a ladder rung is the de Bruijn duplicate.
                    if (row.closure.closed &&
                        (row.closure.mode == Closure::Self || row.closure.mode == Closure::RepeatEnd)) {
                        next = row.closure.seq;
                        why = std::string("duplicate_") + closureName(row.closure.mode);
                        row.edit = "trimmed_duplicate";
                    }
                }
            } else {   // close: closure and duplicate removal, the release tag rule
                if (row.passG) {
                    next = row.closure.seq;
                    why = std::string("closure_") + closureName(row.closure.mode);
                    row.verdict = "closed";
                    row.edit = "closed";
                } else {
                    row.verdict = "no_closure";
                }
            }
            // A record written linear (verify, failed), or one the graph could not close (close),
            // keeps no exact (k-1) terminal self-overlap of any rung of the ladder (L-R3a): each is
            // the de Bruijn duplicate of the record's own start.
            const bool linearOut = in.cfg.circ == CircMode::Verify ? !keepCircular : !row.passG;
            if (linearOut) {
                for (size_t guard = 0; guard < in.ladder.size() + 1; ++guard) {
                    const size_t kov = kMinusOneSelfOverlap(next, in.ladder);
                    if (!kov) break;
                    next = next.substr(0, next.size() - kov);
                    row.kMinusOneTrim += kov;
                }
                if (row.kMinusOneTrim) {
                    why = why.empty() ? std::string("kminus1_self_overlap") : why + "+kminus1_self_overlap";
                    if (row.edit == "none") row.edit = "kminus1_trim";
                    ++res.st.kMinusOneTrims;
                }
            }
            applySeqEdit(in, res, row, row.piece, next, why);
            if (!keepCircular) {
                recordEdit(res, row.piece, "R3", "rename:" + p.tag + ">" + dropCircular(p.tag), 0, 0, 0);
                p.tag = dropCircular(p.tag);
                row.topology = "linear";
            } else {
                row.topology = "circular";
            }
            if (row.topology == "circular") ++res.st.writtenCircular; else ++res.st.writtenLinear;
            if (row.edit == "closed") ++res.st.closed;
            if (row.trimBp) { ++res.st.trimRecords; res.st.trimBp += row.trimBp; }
            if (row.insertBp) { ++res.st.insertRecords; res.st.insertBp += row.insertBp; }
        }
    }

    // ---- propagate to the scaffold records and their tags ---------------------------------------
    std::set<size_t> changed;
    for (RecordRow& row : res.rows) {
        EmitPiece& p = P[row.piece];
        row.lenAfter = p.seq.size();
        const size_t s = p.scaffold;
        if (s >= in.scaffolds->size()) continue;
        const bool seqChanged = row.lenAfter != row.lenBefore || row.edit != "none";
        const bool tagChanged = p.tag != row.oldTag;
        if (!seqChanged && !tagChanged) continue;
        if (wholeScaffold[row.piece]) {
            // A whole untrimmed scaffold is the record: the scaffold takes the record's bases.
            if ((*in.scaffolds)[s] != p.seq) {
                (*in.scaffolds)[s] = p.seq;
                p.cutBackSeq.clear();
                p.objBegin = 0;
                p.objEnd = p.seq.size();
            }
            if ((*in.scaffoldTags)[s] != p.tag) (*in.scaffoldTags)[s] = p.tag;
            changed.insert(s);
        } else if (tagChanged && hasCircular(row.oldTag) && !hasCircular(p.tag) &&
                   hasCircular((*in.scaffoldTags)[s])) {
            // A scaffold whose `_circular` claim did not survive verification is no claim either.
            (*in.scaffoldTags)[s] = dropCircular((*in.scaffoldTags)[s]);
            changed.insert(s);
        } else if (seqChanged) {
            changed.insert(s);   // only the AGP tiling of the trimmed bases changes
        }
    }
    // verify: a scaffold split at a gap carries the pair call's `_circular` in its scaffolds.fasta
    // name although none of its pieces does (the T25 split rule). Such a claim cannot be verified,
    // so under `verify` no record anywhere keeps an unverified `_circular`.
    if (in.cfg.circ == CircMode::Verify) {
        std::vector<size_t> pieceCount(in.scaffolds->size(), 0);
        for (const EmitPiece& p : P)
            if (p.scaffold < pieceCount.size()) ++pieceCount[p.scaffold];
        for (size_t s = 0; s < in.scaffolds->size(); ++s) {
            if (pieceCount[s] < 2 || !hasCircular((*in.scaffoldTags)[s])) continue;
            (*in.scaffoldTags)[s] = dropCircular((*in.scaffoldTags)[s]);
            ++res.st.scaffoldClaimsDropped;
            changed.insert(s);
        }
    }
    res.scaffoldsChanged.assign(changed.begin(), changed.end());
    return res;
}

// ---- side files ----------------------------------------------------------------------------------
namespace {

std::string fmtD(double v, const char* f = "%.4f") {
    char b[64];
    std::snprintf(b, sizeof b, f, v);
    return b;
}

std::string classOf(const std::string& tag) {
    std::string t = dropCircular(tag);
    if (!t.empty() && t[0] == '_') t.erase(0, 1);
    return t.empty() ? "-" : t;
}

}  // namespace

bool writeRepliconsTsv(const std::string& path, const RecordStageResult& r,
                       const std::vector<EmitPiece>& pieces, std::string& error) {
    std::FILE* f = std::fopen(path.c_str(), "w");
    if (!f) { error = "cannot write " + path; return false; }
    std::fprintf(f, "# TesserACT phase-2 replicons.tsv (EVAL_PLAN_P2 W1 R3/F5): one row per contigs.fasta record.\n"
                    "# topology=circular only for a record written with the `_circular` claim. candidate=1: the pair call\n"
                    "# tagged the record circular; closure/junction columns are the R3 evidence (G: graph closure,\n"
                    "# I: graph-isolated, J: reads span the closed join, X: exit ratio). '-' = not computed.\n");
    std::fprintf(f, "record\tlength\tclass\ttopology\tcandidate\tverdict\tedit\tlength_before\ttrim_bp\tinsert_bp\t"
                    "kminus1_trim_bp\tclosure\tclosed_len\tpath_n\tgraph_isolated\tout_last\tin_first\talts\tvia_seg\t"
                    "via_len\tvia_dp\talt_dp\tinterior\tvia_other_links\tcircle_len\tself_ov_trimmed\tspan_reads\t"
                    "clip_reads\tclip_frac\tctrl_span\tctrl_clip\treads\tmate_elsewhere\tnear_join\texits_near_join\t"
                    "exit_ratio\tpass_g\tpass_i\tpass_j\tpass_x\tphix_kmer_frac\tspikein\told_tag\n");
    for (const RecordRow& row : r.rows) {
        const EmitPiece& p = pieces[row.piece];
        const ClosureResult& c = row.closure;
        const JunctionStats& j = row.junction;
        const bool cand = row.candidate;
        auto num = [&](bool show, size_t v) { return show ? std::to_string(v) : std::string("-"); };
        std::fprintf(f, "%s\t%zu\t%s\t%s\t%d\t%s\t%s\t%zu\t%zu\t%zu\t%zu\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t"
                        "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%s\n",
                     p.name.c_str(), p.seq.size(), row.spikein ? "spikein" : classOf(p.tag).c_str(),
                     row.topology.c_str(), cand ? 1 : 0, row.verdict.c_str(), row.edit.c_str(), row.lenBefore,
                     row.trimBp, row.insertBp, row.kMinusOneTrim,
                     cand ? closureName(c.mode) : "-",
                     cand && c.closed ? std::to_string(c.seq.size()).c_str() : "-",
                     num(cand, c.pathN).c_str(),
                     cand && c.isolated >= 0 ? std::to_string(c.isolated).c_str() : "-",
                     cand && c.outLast >= 0 ? std::to_string(c.outLast).c_str() : "-",
                     cand && c.inFirst >= 0 ? std::to_string(c.inFirst).c_str() : "-",
                     num(cand, c.alts).c_str(),
                     cand && !c.viaSeg.empty() ? c.viaSeg.c_str() : "-",
                     cand && !c.viaSeg.empty() ? std::to_string(c.viaLen).c_str() : "-",
                     cand && !c.viaSeg.empty() ? fmtD(c.viaDp).c_str() : "-",
                     cand && !c.altDp.empty() ? c.altDp.c_str() : "-",
                     cand && c.mode == Closure::Via ? std::to_string(c.interior).c_str() : "-",
                     cand && c.mode == Closure::Via ? std::to_string(c.viaOtherLinks).c_str() : "-",
                     num(cand && j.computed, j.circleLen).c_str(), num(cand, j.selfOvTrimmed).c_str(),
                     num(cand && j.computed, j.spanReads).c_str(), num(cand && j.computed, j.clipReads).c_str(),
                     cand && j.computed ? fmtD(j.clipFrac).c_str() : "-",
                     num(cand && j.computed, j.ctrlSpan).c_str(), num(cand && j.computed, j.ctrlClip).c_str(),
                     num(cand && j.computed, j.reads).c_str(), num(cand && j.computed, j.mateElsewhere).c_str(),
                     num(cand && j.computed, j.nearJoin).c_str(), num(cand && j.computed, j.exitsNearJoin).c_str(),
                     cand && j.computed ? fmtD(j.exitRatio).c_str() : "-",
                     cand ? (row.passG ? "1" : "0") : "-", cand ? (row.passI ? "1" : "0") : "-",
                     cand ? (row.passJ ? "1" : "0") : "-", cand ? (row.passX ? "1" : "0") : "-",
                     row.phixFrac >= 0 ? fmtD(row.phixFrac).c_str() : "-", row.spikein ? 1 : 0,
                     row.oldTag.empty() ? "-" : row.oldTag.c_str());
    }
    const bool ok = std::ferror(f) == 0;
    const bool closed = std::fclose(f) == 0;
    if (!ok || !closed) { error = "write failed on " + path; return false; }
    return true;
}

std::vector<std::string> formatEditRows(const RecordStageResult& r, const std::vector<EmitPiece>& pieces) {
    // EVAL_PLAN_P2 s5.2: record, feature, operation, start, end, bases. start/end are 0-based,
    // half-open, in the record as assembled before the edit; bases = bases removed or inserted.
    std::vector<std::string> rows;
    rows.reserve(r.edits.size());
    for (const EditRow& e : r.edits) {
        rows.push_back(pieces[e.piece].name + "\t" + e.feature + "\t" + e.operation + "\t" + std::to_string(e.start) +
                       "\t" + std::to_string(e.end) + "\t" + std::to_string(e.bases));
    }
    return rows;
}

bool writeEditLog(const std::string& path, const std::vector<std::string>& rows, std::string& error) {
    std::FILE* f = std::fopen(path.c_str(), "w");
    if (!f) { error = "cannot write " + path; return false; }
    std::fprintf(f, "record\tfeature\toperation\tstart\tend\tbases\n");
    for (const std::string& row : rows) std::fprintf(f, "%s\n", row.c_str());
    const bool ok = std::ferror(f) == 0;
    const bool closed = std::fclose(f) == 0;
    if (!ok || !closed) { error = "write failed on " + path; return false; }
    return true;
}

bool writeEditsTsv(const std::string& path, const RecordStageResult& r, const std::vector<EmitPiece>& pieces,
                   std::string& error) {
    return writeEditLog(path, formatEditRows(r, pieces), error);
}

// ---- counter lines and report.json ---------------------------------------------------------------
std::string formatLibGuardCounters(const LibGuardStats& s) {
    char b[512];
    std::snprintf(b, sizeof b,
                  "[p2-r2] enabled=%d evaluated=%d fired=%d reason=%d pairs=%zu same_unitig=%zu inward=%zu "
                  "outward=%zu dovetail=%zu same_strand=%zu fr_observations=%zu outward_frac=%.4f fr_fit=%.4f",
                  s.enabled ? 1 : 0, s.evaluated ? 1 : 0, s.fired ? 1 : 0, s.reason, s.pairs, s.sameUnitig, s.inward,
                  s.outward, s.dovetail, s.sameStrand, s.frObservations, s.outwardFrac, s.frFit);
    return b;
}

std::string formatCircCounters(const EmitConfig& c, const RecordStageStats& s) {
    char b[1024];
    std::snprintf(b, sizeof b,
                  "[p2-r3] enabled=%d mode=%d candidates=%zu closure_self=%zu closure_via=%zu closure_repeat_end=%zu "
                  "closure_none=%zu closure_nopath=%zu closure_ends_changed=%zu closure_gapped=%zu isolated=%zu "
                  "pass_j=%zu pass_x=%zu verified=%zu written_circular=%zu written_linear=%zu closed=%zu "
                  "trim_records=%zu trim_bp=%zu insert_records=%zu insert_bp=%zu kminus1_trims=%zu "
                  "scaffold_claims_dropped=%zu reads_placed=%zu window_hi=%d",
                  c.circ != CircMode::Off ? 1 : 0, static_cast<int>(c.circ), s.candidates,
                  s.byMode[static_cast<int>(Closure::Self)], s.byMode[static_cast<int>(Closure::Via)],
                  s.byMode[static_cast<int>(Closure::RepeatEnd)], s.byMode[static_cast<int>(Closure::None)],
                  s.byMode[static_cast<int>(Closure::NoPath)], s.byMode[static_cast<int>(Closure::EndsChanged)],
                  s.byMode[static_cast<int>(Closure::Gapped)], s.isolated, s.passJ, s.passX, s.verified,
                  s.writtenCircular, s.writtenLinear, s.closed, s.trimRecords, s.trimBp, s.insertRecords, s.insertBp,
                  s.kMinusOneTrims, s.scaffoldClaimsDropped, s.verify.readsPlaced, s.verify.windowHi);
    return b;
}

std::string formatSpikeinCounters(const EmitConfig& c, const RecordStageStats& s) {
    char b[256];
    std::snprintf(b, sizeof b, "[p2-f5] enabled=%d screened=%zu spikein=%zu spikein_bp=%zu spikein_circular=%zu",
                  c.spikein ? 1 : 0, s.screened, s.spikein, s.spikeinBp, s.spikeinCircular);
    return b;
}

namespace {

std::string jstr(const std::string& v) {
    std::string o = "\"";
    for (char ch : v) {
        if (ch == '"' || ch == '\\') { o += '\\'; o += ch; }
        else if (static_cast<unsigned char>(ch) < 0x20) { char b[8]; std::snprintf(b, sizeof b, "\\u%04x", ch); o += b; }
        else o += ch;
    }
    return o + "\"";
}

std::string jnum(double v) {
    char b[64];
    std::snprintf(b, sizeof b, "%.6g", v);
    return b;
}

}  // namespace

std::string reportJson(const EmitConfig& c, const LibGuardStats& lg, const RecordStageResult* r,
                       const std::vector<EmitPiece>* pieces) {
    return "{\"enabled\": " + std::to_string(c.any() ? 1 : 0) + ", " + reportJsonMembers(c, lg, r, pieces) + "}";
}

std::string reportJsonMembers(const EmitConfig& c, const LibGuardStats& lg, const RecordStageResult* r,
                              const std::vector<EmitPiece>* pieces) {
    static const RecordStageResult kEmpty;
    const RecordStageResult& rr = r ? *r : kEmpty;
    const RecordStageStats& s = rr.st;
    std::string o;
    // R2
    o += "\"r2\": {\"enabled\": " + std::to_string(lg.enabled ? 1 : 0) +
         ", \"evaluated\": " + std::to_string(lg.evaluated ? 1 : 0) + ", \"fired\": " + std::to_string(lg.fired ? 1 : 0) +
         ", \"reason\": " + std::to_string(lg.reason) + ", \"pairs\": " + std::to_string(lg.pairs) +
         ", \"same_unitig\": " + std::to_string(lg.sameUnitig) + ", \"inward\": " + std::to_string(lg.inward) +
         ", \"outward\": " + std::to_string(lg.outward) + ", \"dovetail\": " + std::to_string(lg.dovetail) +
         ", \"same_strand\": " + std::to_string(lg.sameStrand) +
         ", \"fr_observations\": " + std::to_string(lg.frObservations) + ", \"outward_frac\": " + jnum(lg.outwardFrac) +
         ", \"fr_fit\": " + jnum(lg.frFit);
    if (lg.enabled) {
        o += ", \"verdict\": " + jstr(libguard::verdict(lg)) + ", \"max_outward_frac\": " + jnum(libguard::kMaxOutwardFrac) +
             ", \"min_fr_fit\": " + jnum(libguard::kMinFrFit);
    }
    o += "}";
    // R3
    const bool r3 = c.circ != CircMode::Off;
    o += ", \"r3\": {\"enabled\": " + std::to_string(r3 ? 1 : 0) + ", \"mode\": " + std::to_string(static_cast<int>(c.circ)) +
         ", \"candidates\": " + std::to_string(s.candidates) +
         ", \"closure_self\": " + std::to_string(s.byMode[static_cast<int>(Closure::Self)]) +
         ", \"closure_via\": " + std::to_string(s.byMode[static_cast<int>(Closure::Via)]) +
         ", \"closure_repeat_end\": " + std::to_string(s.byMode[static_cast<int>(Closure::RepeatEnd)]) +
         ", \"closure_none\": " + std::to_string(s.byMode[static_cast<int>(Closure::None)]) +
         ", \"closure_nopath\": " + std::to_string(s.byMode[static_cast<int>(Closure::NoPath)]) +
         ", \"closure_ends_changed\": " + std::to_string(s.byMode[static_cast<int>(Closure::EndsChanged)]) +
         ", \"closure_gapped\": " + std::to_string(s.byMode[static_cast<int>(Closure::Gapped)]) +
         ", \"isolated\": " + std::to_string(s.isolated) + ", \"pass_j\": " + std::to_string(s.passJ) +
         ", \"pass_x\": " + std::to_string(s.passX) + ", \"verified\": " + std::to_string(s.verified) +
         ", \"written_circular\": " + std::to_string(s.writtenCircular) +
         ", \"written_linear\": " + std::to_string(s.writtenLinear) + ", \"closed\": " + std::to_string(s.closed) +
         ", \"trim_records\": " + std::to_string(s.trimRecords) + ", \"trim_bp\": " + std::to_string(s.trimBp) +
         ", \"insert_records\": " + std::to_string(s.insertRecords) + ", \"insert_bp\": " + std::to_string(s.insertBp) +
         ", \"kminus1_trims\": " + std::to_string(s.kMinusOneTrims) +
         ", \"scaffold_claims_dropped\": " + std::to_string(s.scaffoldClaimsDropped) +
         ", \"reads_placed\": " + std::to_string(s.verify.readsPlaced);
    if (r3) {
        const VerifyThresholds t;
        o += ", \"mode_name\": " + jstr(circModeName(c.circ)) + ", \"window\": [" + std::to_string(s.verify.windowLo) +
             ", " + std::to_string(s.verify.windowHi) + "]" + ", \"thresholds\": {\"min_span\": " + std::to_string(t.minSpan) +
             ", \"min_span_ctrl\": " + jnum(t.minSpanCtrl) + ", \"max_clip_frac\": " + jnum(t.maxClipFrac) +
             ", \"max_exit_ratio\": " + jnum(t.maxExitRatio) + ", \"anchor_k\": " + std::to_string(kAnchorK) + "}";
        o += ", \"circles\": [";
        bool firstRow = true;
        for (const RecordRow& row : rr.rows) {
            if (!row.candidate || !pieces) continue;
            const EmitPiece& p = (*pieces)[row.piece];
            const JunctionStats& j = row.junction;
            o += firstRow ? "" : ", ";
            firstRow = false;
            o += "{\"record\": " + jstr(p.name) + ", \"topology\": " + jstr(row.topology) +
                 ", \"verdict\": " + jstr(row.verdict) + ", \"edit\": " + jstr(row.edit) +
                 ", \"closure\": " + jstr(closureName(row.closure.mode)) +
                 ", \"length_before\": " + std::to_string(row.lenBefore) + ", \"length\": " + std::to_string(row.lenAfter) +
                 ", \"trim_bp\": " + std::to_string(row.trimBp) + ", \"insert_bp\": " + std::to_string(row.insertBp) +
                 ", \"kminus1_trim_bp\": " + std::to_string(row.kMinusOneTrim) +
                 ", \"graph_isolated\": " + std::to_string(row.closure.isolated) +
                 ", \"alts\": " + std::to_string(row.closure.alts) +
                 ", \"circle_len\": " + std::to_string(j.circleLen) + ", \"span_reads\": " + std::to_string(j.spanReads) +
                 ", \"clip_reads\": " + std::to_string(j.clipReads) + ", \"clip_frac\": " + jnum(j.clipFrac) +
                 ", \"ctrl_span\": " + std::to_string(j.ctrlSpan) + ", \"near_join\": " + std::to_string(j.nearJoin) +
                 ", \"exits_near_join\": " + std::to_string(j.exitsNearJoin) + ", \"exit_ratio\": " + jnum(j.exitRatio) +
                 ", \"pass_g\": " + std::to_string(row.passG ? 1 : 0) + ", \"pass_i\": " + std::to_string(row.passI ? 1 : 0) +
                 ", \"pass_j\": " + std::to_string(row.passJ ? 1 : 0) + ", \"pass_x\": " + std::to_string(row.passX ? 1 : 0) +
                 ", \"phix_kmer_frac\": " + jnum(row.phixFrac) + "}";
        }
        o += "]";
    }
    o += "}";
    // F5
    o += ", \"f5\": {\"enabled\": " + std::to_string(c.spikein ? 1 : 0) + ", \"screened\": " + std::to_string(s.screened) +
         ", \"spikein\": " + std::to_string(s.spikein) + ", \"spikein_bp\": " + std::to_string(s.spikeinBp) +
         ", \"spikein_circular\": " + std::to_string(s.spikeinCircular);
    if (c.spikein) {
        o += ", \"set\": \"phiX174 NC_001422.1\", \"min_kmer_frac\": " + jnum(kSpikeinMinFrac) + ", \"records\": [";
        bool firstRow = true;
        for (const RecordRow& row : rr.rows) {
            if (!row.spikein || !pieces) continue;
            o += (firstRow ? "" : ", ") + jstr((*pieces)[row.piece].name);
            firstRow = false;
        }
        o += "]";
    }
    o += "}";
    return o;
}

}  // namespace p2
}  // namespace ts
