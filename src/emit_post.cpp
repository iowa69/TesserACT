#include "emit_post.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <memory>
#include <thread>
#include <unordered_map>
#include <utility>

#include "envflags.h"
#include "kmer.h"

namespace ts {

// ---- terminal dovetails ----------------------------------------------------------------------
// Moved unchanged from assembler.cpp (release 1.3.0 lines 391-452), so the trim can be tested.
// NOTE for G-xcut (T16): the TESSERACT_RC_DOVETAIL function-local static moved with it.
std::vector<Dovetail> findTerminalDovetails(const std::vector<std::string>& seqs,
                                            size_t minOverlap, size_t window) {
    constexpr size_t kProbe = 32;
    std::vector<Dovetail> out;
    if (seqs.size() < 2 || minOverlap < kProbe) return out;
    auto hash = [](const char* p) {
        uint64_t h = 1469598103934665603ULL;
        for (size_t i = 0; i < kProbe; ++i) { h ^= static_cast<unsigned char>(p[i]); h *= 1099511628211ULL; }
        return h;
    };
    // index: hash of a 32-mer inside the PREFIX window of each oriented contig
    struct Ent { uint32_t idx; uint32_t pos; bool rc; };
    std::unordered_map<uint64_t, std::vector<Ent>> index;
    std::vector<std::string> rcs(seqs.size());
    for (size_t j = 0; j < seqs.size(); ++j) {
        if (seqs[j].size() < minOverlap) continue;
        rcs[j] = reverseComplement(seqs[j]);
        for (int o = 0; o < 2; ++o) {
            const std::string& t = o ? rcs[j] : seqs[j];
            const size_t lim = std::min(window, t.size()) - kProbe + 1;
            for (size_t p = 0; p < lim; ++p)
                index[hash(t.data() + p)].push_back({static_cast<uint32_t>(j),
                                                     static_cast<uint32_t>(p), o != 0});
        }
    }
    // The index holds both orientations of every contig's PREFIX window, but the query
    // below probed only the forward 3' end, so a head-to-head overlap (a's 5' end meeting
    // b's 5' end) was structurally invisible. Probing rc(a)'s 3' end as well closes that
    // class.
    // ON by default since 1.3.0 (duplication LOSS -> tie on 146 isolates).
    // TESSERACT_RC_DOVETAIL=0 restores the forward-only probe.
    // build_v3 (T16): read per call through the envflags table, as G-xcut's retrofit did at
    // the function's old home in assembler.cpp.
    const bool rcQuery = env::on("TESSERACT_RC_DOVETAIL", true);
    for (size_t i = 0; i < seqs.size(); ++i) {
        if (seqs[i].size() < minOverlap) continue;
        for (int qo = 0; qo < (rcQuery ? 2 : 1); ++qo) {
            const std::string& a = qo ? rcs[i] : seqs[i];
            if (a.size() < minOverlap) continue;
            auto it = index.find(hash(a.data() + a.size() - kProbe));
            if (it == index.end()) continue;
            for (const Ent& e : it->second) {
                if (e.idx == i) continue;
                const size_t L = static_cast<size_t>(e.pos) + kProbe;
                if (L <= minOverlap) continue;
                const std::string& b = e.rc ? rcs[e.idx] : seqs[e.idx];
                if (L > a.size() || L > b.size()) continue;
                if (std::memcmp(a.data() + a.size() - L, b.data(), L) != 0) continue;
                out.push_back({i, e.idx, L, e.rc, qo != 0});
            }
        }
    }
    std::sort(out.begin(), out.end(), [](const Dovetail& x, const Dovetail& y) {
        if (x.len != y.len) return x.len > y.len;      // longest first: it decides the end
        if (x.a != y.a) return x.a < y.a;
        if (x.b != y.b) return x.b < y.b;
        if (x.aRc != y.aRc) return !x.aRc;
        return !x.bRc && y.bRc;                         // total order, so runs reproduce
    });
    return out;
}

// ---- T37: exact containment with an exact negative prefilter -------------------------------
namespace {

// A one-hash bit filter over every 32-byte window of every string appended to the text.
// mayOccur(x) is false only if some window of x was never added -- and then x cannot occur in
// the text, because an occurrence never spans the '\x01' separator (x holds no such byte), so
// it lies inside one appended string whose every window was added. A false positive only runs
// the exact find the release always ran. Strings shorter than one window always fall through.
struct WindowFilter {
    static constexpr size_t kWin = 32;
    static constexpr uint64_t kBase = 1099511628211ULL;
    std::vector<uint64_t> bits;
    uint64_t mask = 0;
    uint64_t basePowW = 1;
    explicit WindowFilter(size_t totalBases) {
        size_t nbits = size_t(1) << 16;
        while (nbits < 8 * totalBases) nbits <<= 1;
        bits.assign(nbits / 64, 0);
        mask = nbits - 1;
        for (size_t i = 0; i < kWin; ++i) basePowW *= kBase;
    }
    uint64_t slot(uint64_t h) const { return ((h * 0x9E3779B97F4A7C15ULL) >> 17) & mask; }
    template <class F> bool forEachWindow(const std::string& s, F f) const {
        if (s.size() < kWin) return true;
        uint64_t h = 0;
        for (size_t i = 0; i < kWin; ++i) h = h * kBase + static_cast<unsigned char>(s[i]);
        if (!f(h)) return false;
        for (size_t i = kWin; i < s.size(); ++i) {
            h = h * kBase + static_cast<unsigned char>(s[i]) -
                basePowW * static_cast<unsigned char>(s[i - kWin]);
            if (!f(h)) return false;
        }
        return true;
    }
    void add(const std::string& s) {
        forEachWindow(s, [&](uint64_t h) {
            const uint64_t b = slot(h);
            bits[b >> 6] |= uint64_t(1) << (b & 63);
            return true;
        });
    }
    bool mayOccur(const std::string& x) const {
        return forEachWindow(x, [&](uint64_t h) {
            const uint64_t b = slot(h);
            return ((bits[b >> 6] >> (b & 63)) & 1) != 0;
        });
    }
};

}  // namespace

std::vector<char> markContained(const std::vector<std::string>& seqs,
                                const std::vector<char>& protectedSeq,
                                const std::vector<char>* candidate, DedupStats* stats) {
    std::vector<char> contained(seqs.size(), 0);
    DedupStats local;
    DedupStats& st = stats ? *stats : local;
    if (seqs.size() <= 1) return contained;
    std::vector<size_t> byLen(seqs.size());
    for (size_t i = 0; i < byLen.size(); ++i) byLen[i] = i;
    std::sort(byLen.begin(), byLen.end(), [&](size_t a, size_t b) {
        if (seqs[a].size() != seqs[b].size()) return seqs[a].size() > seqs[b].size();
        return seqs[a] < seqs[b];          // total order, so the result is reproducible
    });
    size_t total = 0;
    for (const std::string& s : seqs) total += s.size();
    std::string text;                       // the sequences kept so far, longest first
    text.reserve(total + seqs.size() + 1);
    WindowFilter filter(total);
    for (size_t idx : byLen) {
        const std::string& s = seqs[idx];
        const bool mayDrop = (idx >= protectedSeq.size() || !protectedSeq[idx]) &&
                             (!candidate || (*candidate)[idx]);
        bool hit = false;
        if (!s.empty() && !text.empty() && mayDrop) {
            ++st.queried;
            if (filter.mayOccur(s)) { ++st.exactSearches; hit = text.find(s) != std::string::npos; }
            if (!hit) {
                const std::string r = reverseComplement(s);
                if (filter.mayOccur(r)) { ++st.exactSearches; hit = text.find(r) != std::string::npos; }
            }
        }
        if (hit) {
            contained[idx] = 1;
            ++st.dropped;
            st.droppedBases += s.size();
        } else {
            text.push_back('\x01');        // a separator no base can match across
            text.append(s);
            filter.add(s);
        }
    }
    return contained;
}

// ---- the split at every N ------------------------------------------------------------------
std::vector<EmitPiece> splitAtGaps(const std::vector<std::string>& scaffolds,
                                   const std::vector<std::string>& tags,
                                   const std::vector<double>& covs) {
    std::vector<EmitPiece> out;
    out.reserve(scaffolds.size());
    for (size_t i = 0; i < scaffolds.size(); ++i) {
        const std::string& sq = scaffolds[i];
        const std::string tag = i < tags.size() ? tags[i] : std::string();
        const double cov = i < covs.size() ? covs[i] : 0.0;
        const size_t first = out.size();
        size_t pos = 0;
        while (pos < sq.size()) {
            while (pos < sq.size() && (sq[pos] == 'N' || sq[pos] == 'n')) ++pos;
            if (pos >= sq.size()) break;
            size_t e = pos;
            while (e < sq.size() && sq[e] != 'N' && sq[e] != 'n') ++e;
            EmitPiece p;
            p.scaffold = i;
            p.ordinal = out.size() - first;
            p.objBegin = pos;
            p.objEnd = e;
            p.seq = sq.substr(pos, e - pos);
            p.tag = tag;
            p.cov = cov;
            out.push_back(std::move(p));
            pos = e;
        }
        for (size_t j = first; j < out.size(); ++j) out[j].piecesInScaffold = out.size() - first;
    }
    return out;
}

// ---- scaffold paths cut at their gap steps -------------------------------------------------
std::vector<PathSegment> pathSegments(const GfaPath& path, const UnitigGraph& g) {
    std::vector<PathSegment> segs;
    for (size_t i = 0; i < path.oriented.size(); ++i) {
        const int gap = i < path.gaps.size() ? path.gaps[i] : 0;
        if (segs.empty() || (i > 0 && gap > 0)) segs.emplace_back();
        const uint32_t u = static_cast<uint32_t>(path.oriented[i] >> 1);
        segs.back().oriented.push_back(path.oriented[i]);
        if (u < g.nodes.size()) {
            const double len = static_cast<double>(g.nodes[u].seq.size());
            segs.back().covWeighted += g.nodes[u].coverage * len;
            segs.back().covLength += len;
        }
    }
    return segs;
}

std::string spellWalk(const std::vector<uint64_t>& oriented, const UnitigGraph& g) {
    std::string s;
    const size_t ov = static_cast<size_t>(std::max(0, g.k() - 1));
    for (uint64_t o : oriented) {
        const uint32_t u = static_cast<uint32_t>(o >> 1);
        if (u >= g.nodes.size()) return std::string();
        const std::string piece = g.oriented(u, static_cast<int>(o & 1));
        if (s.empty()) s = piece;
        else s += piece.size() > ov ? piece.substr(ov) : std::string();
    }
    return s;
}

std::vector<size_t> mapSegmentsToPieces(const std::string& scaffold,
                                        const std::vector<PathSegment>& segs,
                                        const UnitigGraph& g) {
    constexpr size_t kMinProbe = 24;       // shorter spellings could land anywhere
    constexpr size_t kTipAllowance = 200;  // bases a gap-filler back-off may re-spell at a segment end
    std::vector<size_t> runStart;          // starts of the N-runs, in order
    for (size_t p = 0; p < scaffold.size();) {
        if (scaffold[p] != 'N' && scaffold[p] != 'n') { ++p; continue; }
        runStart.push_back(p);
        while (p < scaffold.size() && (scaffold[p] == 'N' || scaffold[p] == 'n')) ++p;
    }
    const size_t nPieces = [&] {
        size_t n = 0;
        for (size_t p = 0; p < scaffold.size();) {
            while (p < scaffold.size() && (scaffold[p] == 'N' || scaffold[p] == 'n')) ++p;
            if (p >= scaffold.size()) break;
            ++n;
            while (p < scaffold.size() && scaffold[p] != 'N' && scaffold[p] != 'n') ++p;
        }
        return n;
    }();
    const size_t ov = static_cast<size_t>(std::max(0, g.k() - 1));
    std::vector<size_t> pieceOf;
    size_t from = 0;
    for (size_t j = 0; j < segs.size(); ++j) {
        const std::string spelled = spellWalk(segs[j].oriented, g);
        // After a gap the release renderer drops the first k-1 bases (the gap-flank fix keeps
        // them); the part after them is rendered either way.
        const std::string sfx = j == 0 ? spelled : (spelled.size() > ov ? spelled.substr(ov) : std::string());
        size_t hitEnd = std::string::npos;
        if (sfx.size() >= kMinProbe) {
            const size_t q = scaffold.find(sfx, from);
            if (q != std::string::npos) hitEnd = q + sfx.size();
        }
        if (hitEnd == std::string::npos && sfx.size() >= 2 * kTipAllowance + kMinProbe) {
            // A closure after back-off re-spells a few flank-tip bases at either end.
            const std::string core = sfx.substr(kTipAllowance, sfx.size() - 2 * kTipAllowance);
            const size_t q = scaffold.find(core, from);
            if (q != std::string::npos) hitEnd = q + core.size();
        }
        if (hitEnd == std::string::npos) return {};
        // The piece holding base hitEnd-1 is the number of N-runs that start before it.
        const size_t piece = static_cast<size_t>(
            std::lower_bound(runStart.begin(), runStart.end(), hitEnd - 1) - runStart.begin());
        if (!pieceOf.empty() && piece < pieceOf.back()) return {};
        pieceOf.push_back(piece);
        from = hitEnd;
    }
    // Every piece must be accounted for by a segment, and the walk must end in the last one.
    if (pieceOf.empty() || pieceOf.back() + 1 != nPieces) return {};
    for (size_t j = 1; j < pieceOf.size(); ++j)
        if (pieceOf[j] > pieceOf[j - 1] + 1) return {};
    if (pieceOf.front() != 0) return {};
    return pieceOf;
}

// ---- T25 ------------------------------------------------------------------------------------
void splitPostprocess(std::vector<EmitPiece>& pieces, const std::vector<std::string>& scaffolds,
                      const std::vector<GfaPath>& paths, const UnitigGraph* graph,
                      const std::vector<char>& protectedScaffold, size_t minLen,
                      SplitPostStats& st) {
    st.pieces = pieces.size();
    std::vector<std::vector<size_t>> byScaffold(scaffolds.size());
    for (size_t i = 0; i < pieces.size(); ++i) byScaffold[pieces[i].scaffold].push_back(i);
    for (const auto& v : byScaffold) st.multiPieceScaffolds += v.size() > 1 ? 1 : 0;
    if (!st.enabled) return;

    for (size_t s = 0; s < scaffolds.size(); ++s) {
        const std::vector<size_t>& mine = byScaffold[s];
        if (mine.size() < 2) continue;
        // (a) Coverage from the piece's own stretch of the walk. The scaffold label is the
        //     length-weighted mean over every unitig in the chain, so a 60x chromosome piece
        //     and a 180x repeat piece were both written with the chain's mean.
        bool relabelled = false;
        if (graph && s < paths.size() && !paths[s].oriented.empty()) {
            const std::vector<PathSegment> segs = pathSegments(paths[s], *graph);
            const std::vector<size_t> pieceOf = mapSegmentsToPieces(scaffolds[s], segs, *graph);
            if (!pieceOf.empty() && pieceOf.back() + 1 == mine.size()) {
                std::vector<double> w(mine.size(), 0), l(mine.size(), 0);
                for (size_t j = 0; j < segs.size(); ++j) {
                    w[pieceOf[j]] += segs[j].covWeighted;
                    l[pieceOf[j]] += segs[j].covLength;
                }
                bool ok = true;
                for (size_t j = 0; j < mine.size(); ++j) ok = ok && l[j] > 0;
                if (ok) {
                    for (size_t j = 0; j < mine.size(); ++j) pieces[mine[j]].cov = w[j] / l[j];
                    st.covRelabelled += mine.size();
                    relabelled = true;
                }
            }
        }
        if (!relabelled) st.covFallback += mine.size();
        // (b) A scaffold whose own two ends the pairs joined is circular; none of its pieces is.
        for (size_t idx : mine) {
            const size_t c = pieces[idx].tag.find("_circular");
            if (c != std::string::npos) {
                pieces[idx].tag.erase(c, std::strlen("_circular"));
                ++st.circularTagsDropped;
            }
        }
        // (c) The same length floor every whole record met before the split (--min-contig / 2k).
        for (size_t idx : mine) {
            if (pieces[idx].seq.size() < minLen) {
                pieces[idx].drop = 1;
                ++st.droppedShort;
                st.droppedShortBp += pieces[idx].seq.size();
            }
        }
    }
    // (d) Exact containment again, now that pieces exist: the pre-split test saw scaffolds with
    //     their Ns and could not see a piece inside another record. Only pieces of multi-piece
    //     scaffolds may go (every other record already passed the test), never a piece of a
    //     plasmid or circular call.
    std::vector<size_t> live;
    for (size_t i = 0; i < pieces.size(); ++i) if (!pieces[i].drop) live.push_back(i);
    std::vector<std::string> seqs;
    std::vector<char> prot, cand;
    seqs.reserve(live.size());
    for (size_t i : live) {
        seqs.push_back(pieces[i].seq);
        const size_t s = pieces[i].scaffold;
        prot.push_back(s < protectedScaffold.size() && protectedScaffold[s] ? 1 : 0);
        cand.push_back(pieces[i].piecesInScaffold > 1 ? 1 : 0);
    }
    const std::vector<char> dup = markContained(seqs, prot, &cand, nullptr);
    for (size_t j = 0; j < live.size(); ++j) {
        if (!dup[j]) continue;
        pieces[live[j]].drop = 2;
        ++st.droppedDup;
        st.droppedDupBp += pieces[live[j]].seq.size();
    }
}

void logSplitPost(const SplitPostStats& st) {
    std::fprintf(stderr,
                 "[splitpost] enabled=%d pieces=%zu multi_piece_scaffolds=%zu dropped_short=%zu "
                 "dropped_short_bp=%zu dropped_dup=%zu dropped_dup_bp=%zu circular_tags_dropped=%zu "
                 "cov_relabelled=%zu cov_fallback=%zu\n",
                 st.enabled ? 1 : 0, st.pieces, st.multiPieceScaffolds, st.droppedShort,
                 st.droppedShortBp, st.droppedDup, st.droppedDupBp, st.circularTagsDropped,
                 st.covRelabelled, st.covFallback);
}

// ---- T10 / T22: the terminal dovetail trim -------------------------------------------------
namespace {

// Contig end ids: 2*seq for the 5' end, 2*seq+1 for the 3' end.
inline size_t endOfA(const Dovetail& d) { return 2 * d.a + (d.aRc ? 0 : 1); }
inline size_t endOfB(const Dovetail& d) { return 2 * d.b + (d.bRc ? 1 : 0); }

// Canonical k-mers of s[from, to) (N resets the window), with the start position.
template <class F> void forEachCanonicalKmer(const std::string& s, size_t from, size_t to, int k, F fn) {
    Kmer fwd = 0, rc = 0;
    int valid = 0;
    for (size_t p = from; p < to; ++p) {
        const int b = baseCode(s[p]);
        if (b < 0) { valid = 0; continue; }
        fwd = pushBack(fwd, b, k);
        rc = pushFrontRc(rc, b, k);
        if (++valid < k) continue;
        fn(fwd < rc ? fwd : rc, p + 1 - static_cast<size_t>(k));
    }
}

// The non-guard victim end of one record: the release rule, or T10's strict order.
inline size_t baselineVictim(const Dovetail& d, const std::vector<std::string>& seqs, bool oneVictim,
                             bool* flipped = nullptr) {
    bool aLonger = seqs[d.a].size() >= seqs[d.b].size();
    if (oneVictim && seqs[d.a].size() == seqs[d.b].size() && d.a > d.b) {
        aLonger = false;              // equal lengths: the higher index is always the victim
        if (flipped) *flipped = true;
    }
    return aLonger ? endOfB(d) : endOfA(d);
}

}  // namespace

TrimStats trimTerminalOverlaps(std::vector<std::string>& seqs, const TrimConfig& cfg,
                               std::vector<size_t>& cutFront, std::vector<size_t>& cutBack,
                               const EndEvidenceFn* evidence) {
    TrimStats st;
    cutFront.assign(seqs.size(), 0);
    cutBack.assign(seqs.size(), 0);
    const std::vector<Dovetail> dv = findTerminalDovetails(seqs, cfg.minOverlap, 8000);
    const bool oneVictim = cfg.boundarySafe;

    // ---- T22: which ends each overlap component keeps ------------------------------------
    // Ends that share one terminal segment form a bipartite graph: the ends that run INTO the
    // segment and the ends that run OUT of it. The release rule cuts the shorter partner of
    // each record and never looks at the keeper, so an end keeps its copy only if it is longer
    // than every partner, and a segment present at two loci keeps one copy. Keeping one whole
    // side keeps max(entering, exiting) copies -- the lower bound on the copy number when the
    // walks are right -- and is string-lossless by construction. Because a longer copy set is
    // exactly what re-admits a wrongly walked terminal copy, a component departs from the
    // non-guard decision only when read pairs support every end it rescues and the segment's
    // read depth carries that many copies.
    std::unordered_map<size_t, uint8_t> guardCut;   // end id -> 1 cut / 0 kept, applied components only
    if (cfg.copyGuard && !dv.empty()) {
        std::unordered_map<size_t, std::vector<size_t>> adj;   // end -> neighbouring ends
        for (const Dovetail& d : dv) {
            adj[endOfA(d)].push_back(endOfB(d));
            adj[endOfB(d)].push_back(endOfA(d));
        }
        // Non-guard cut set (victim of some record, stub guard aside).
        std::unordered_map<size_t, uint8_t> baseCut;
        for (const Dovetail& d : dv) baseCut[baselineVictim(d, seqs, oneVictim)] = 1;
        std::vector<size_t> ends;
        for (const auto& kv : adj) ends.push_back(kv.first);
        std::sort(ends.begin(), ends.end());
        std::unordered_map<size_t, int> colour;
        struct Comp { std::vector<size_t> side[2]; bool bip = true; };
        std::vector<Comp> comps;
        for (size_t e0 : ends) {
            if (colour.count(e0)) continue;
            Comp c;
            std::vector<size_t> stack{e0};
            colour[e0] = 0;
            while (!stack.empty()) {
                const size_t e = stack.back();
                stack.pop_back();
                c.side[colour[e]].push_back(e);
                for (size_t f : adj[e]) {
                    auto it = colour.find(f);
                    if (it == colour.end()) { colour[f] = 1 - colour[e]; stack.push_back(f); }
                    else if (it->second == colour[e]) c.bip = false;
                }
            }
            std::sort(c.side[0].begin(), c.side[0].end());
            std::sort(c.side[1].begin(), c.side[1].end());
            comps.push_back(std::move(c));
        }
        st.components = comps.size();
        // Longest record per end: the terminal segment the end holds.
        std::unordered_map<size_t, size_t> segLen;
        for (const Dovetail& d : dv) {
            size_t& a = segLen[endOfA(d)];
            size_t& b = segLen[endOfB(d)];
            a = std::max(a, d.len);
            b = std::max(b, d.len);
        }
        std::vector<EndQuery> queries;
        std::unordered_map<size_t, size_t> queryOf;
        for (const Comp& c : comps) {
            if (!c.bip) continue;
            for (int s = 0; s < 2; ++s)
                for (size_t e : c.side[s]) {
                    queryOf[e] = queries.size();
                    queries.push_back({e / 2, (e % 2) == 0, segLen[e]});
                }
        }
        std::vector<EndEvidence> ev(queries.size());
        if (evidence && *evidence && !queries.empty()) {
            ev = (*evidence)(queries);
            if (ev.size() != queries.size()) ev.assign(queries.size(), EndEvidence());
        }
        auto longest = [&](const std::vector<size_t>& side) {
            size_t best = 0;
            for (size_t e : side) best = std::max(best, seqs[e / 2].size());
            return best;
        };
        // Pass 1: the side each component would keep, and the ends that keeping it would rescue.
        struct Plan { int keep = 0; bool switched = false, differs = false; };
        std::vector<Plan> plans(comps.size());
        std::vector<size_t> rescueEnds;
        for (size_t ci = 0; ci < comps.size(); ++ci) {
            const Comp& c = comps[ci];
            if (!c.bip) continue;
            Plan& pl = plans[ci];
            pl.keep = c.side[0].size() != c.side[1].size()
                          ? (c.side[0].size() > c.side[1].size() ? 0 : 1)
                          : (longest(c.side[0]) != longest(c.side[1])
                                 ? (longest(c.side[0]) > longest(c.side[1]) ? 0 : 1)
                                 : (c.side[0].front() < c.side[1].front() ? 0 : 1));
            auto anyOf = [&](int s, bool (EndEvidence::*pred)() const) {
                for (size_t e : c.side[s]) if ((ev[queryOf[e]].*pred)()) return true;
                return false;
            };
            // Placement evidence picks the victim: a kept end whose mates say its segment is not
            // what follows its flank hands the copy to the other side when that side is clean.
            if (anyOf(pl.keep, &EndEvidence::contradicted) && !anyOf(1 - pl.keep, &EndEvidence::contradicted) &&
                anyOf(1 - pl.keep, &EndEvidence::supported)) {
                pl.keep = 1 - pl.keep;
                pl.switched = true;
            }
            // Same decision as without the guard? Then nothing to prove.
            bool same = true;
            for (size_t e : c.side[1 - pl.keep]) same = same && baseCut.count(e);
            for (size_t e : c.side[pl.keep]) same = same && !baseCut.count(e);
            pl.differs = !same;
            if (pl.differs)
                for (size_t e : c.side[pl.keep]) if (baseCut.count(e)) rescueEnds.push_back(e);
        }
        // Copies of each rescued segment the assembly already holds away from every contig end
        // taking part in a dovetail (for instance a tandem copy inside the rescued end's own
        // flank): mean over the segment's 31-mers of their occurrences in contig interiors.
        std::unordered_map<size_t, double> interior;
        if (!rescueEnds.empty()) {
            constexpr int kc = 31;
            std::unordered_map<Kmer, uint32_t, KmerHasher> count;
            std::unordered_map<size_t, std::vector<Kmer>> segKm;
            for (size_t e : rescueEnds) {
                const std::string& s = seqs[e / 2];
                const size_t L = std::min(segLen[e], s.size());
                const size_t from = (e % 2) == 0 ? 0 : s.size() - L;
                forEachCanonicalKmer(s, from, from + L, kc, [&](const Kmer& km, size_t) {
                    segKm[e].push_back(km);
                    count.emplace(km, 0);
                });
            }
            std::vector<size_t> maskFront(seqs.size(), 0), maskBack(seqs.size(), 0);
            for (const auto& kv : segLen) {
                if (kv.first % 2 == 0) maskFront[kv.first / 2] = std::max(maskFront[kv.first / 2], kv.second);
                else maskBack[kv.first / 2] = std::max(maskBack[kv.first / 2], kv.second);
            }
            for (size_t i = 0; i < seqs.size(); ++i) {
                const size_t n = seqs[i].size();
                if (maskFront[i] + maskBack[i] >= n) continue;
                forEachCanonicalKmer(seqs[i], maskFront[i], n - maskBack[i], kc, [&](const Kmer& km, size_t) {
                    auto it = count.find(km);
                    if (it != count.end()) ++it->second;
                });
            }
            for (const auto& kv : segKm) {
                double sum = 0;
                for (const Kmer& km : kv.second) sum += count[km];
                interior[kv.first] = kv.second.empty() ? 0 : sum / static_cast<double>(kv.second.size());
            }
        }
        // Pass 2: a component departs from the non-guard decision only when every end it keeps
        // is clean, every end it rescues is supported by pairs, and the reads carry as many
        // copies as the assembly would then hold (kept ends, cut-side ends too short to cut,
        // and interior copies).
        for (size_t ci = 0; ci < comps.size(); ++ci) {
            const Comp& c = comps[ci];
            if (!c.bip) { ++st.nonBipartite; continue; }
            const Plan& pl = plans[ci];
            if (!pl.differs) continue;
            const int keep = pl.keep;
            ++st.guardDiffers;
            size_t rescued = 0;
            bool ok = true;
            // Cut-side ends the stub guard will spare still hold a copy.
            size_t spared = 0;
            for (size_t e : c.side[1 - keep])
                if (segLen[e] + 200 >= seqs[e / 2].size()) ++spared;
            const double kept = static_cast<double>(c.side[keep].size() + spared);
            for (size_t e : c.side[keep]) {
                const EndEvidence& x = ev[queryOf[e]];
                if (x.contradicted()) ok = false;
                if (!baseCut.count(e)) continue;
                ++rescued;
                if (!x.supported()) ok = false;
                if (!(x.copyRatio > 0) || x.copyRatio < static_cast<double>(c.side[keep].size()) - 0.5) ok = false;
                // The depth ratio is a noisy count (a true two-copy segment reads 1.6-2.0), hence
                // the half-copy margin -- but not when the assembly already holds the segment
                // inside a contig: then the reads must carry every copy outright, because the
                // copy being rescued competes with one whose placement is already settled.
                const double margin = interior[e] < 0.5 ? 0.5 : 0.0;
                if (kept + interior[e] > x.copyRatio + margin) ok = false;
            }
            // Audit trail: one line per end of a component the guard would change (flag on only).
            for (int sd = 0; sd < 2; ++sd)
                for (size_t e : c.side[sd]) {
                    const EndEvidence& x = ev[queryOf[e]];
                    const auto it = interior.find(e);
                    std::fprintf(stderr,
                                 "[trimguard-end] seq=%zu end=%s seg=%zu len=%zu consistent=%zu inconsistent=%zu "
                                 "copy_ratio=%.2f interior_copies=%.2f side=%s release=%s decision=%s\n",
                                 e / 2, e % 2 ? "3'" : "5'", segLen[e], seqs[e / 2].size(), x.consistent,
                                 x.inconsistent, x.copyRatio, it == interior.end() ? 0.0 : it->second,
                                 sd == keep ? "keep" : "cut", baseCut.count(e) ? "cut" : "keep",
                                 ok ? "guard" : "release");
                }
            if (!ok) { ++st.evidenceRefused; continue; }
            ++st.guardApplied;
            if (pl.switched) ++st.evidenceSwitched;
            st.endsRescued += rescued;
            st.releaseCutsAvoided += rescued;
            for (size_t e : c.side[keep]) guardCut[e] = 0;
            for (size_t e : c.side[1 - keep]) guardCut[e] = 1;
        }
    }

    // ---- the release loop (victim overridden where the guard applies) ----------------------
    std::vector<char> frontSet(seqs.size(), 0), backSet(seqs.size(), 0);
    for (const Dovetail& d : dv) {
        // Longest dovetails first, and an end already decided is not revisited: two
        // overlaps sharing an end would otherwise cut the same bases twice.
        size_t victimEnd;
        auto g = guardCut.find(endOfA(d));
        if (g != guardCut.end()) {
            victimEnd = g->second ? endOfA(d) : endOfB(d);
        } else {
            bool flipped = false;
            victimEnd = baselineVictim(d, seqs, oneVictim, &flipped);
            if (flipped) ++st.equalLenFlips;
        }
        const size_t victim = victimEnd / 2;
        // a's SUFFIX meets b's PREFIX in b's own orientation; if b was matched
        // reverse-complemented, its prefix there is its suffix here.
        const bool cutAtFront = (victimEnd % 2) == 0;
        std::vector<char>& set = cutAtFront ? frontSet : backSet;
        std::vector<size_t>& cut = cutAtFront ? cutFront : cutBack;
        if (set[victim]) continue;
        // Never leave a stub: a contig that is mostly overlap is the wholly-contained
        // case, which --dedup-contained handles by removing it outright.
        const size_t other = cutAtFront ? cutBack[victim] : cutFront[victim];
        if (d.len + other + 200 >= seqs[victim].size()) continue;
        set[victim] = 1;
        cut[victim] = d.len;
    }
    // T10: stop each cut k-1 bases short, so the victim keeps every k-mer spanning its
    // flank-to-overlap boundary; the keeper still holds the overlap's interior.
    const size_t keep = cfg.boundarySafe && cfg.k > 1 ? static_cast<size_t>(cfg.k - 1) : 0;
    for (size_t i = 0; i < seqs.size(); ++i) {
        for (size_t* c : {&cutFront[i], &cutBack[i]}) {
            if (*c == 0) continue;
            if (keep) {
                if (*c <= keep) { *c = 0; ++st.skippedByBoundary; continue; }
                *c -= keep;
                st.boundaryBpKept += keep;
            }
            ++st.trimmed;
            st.trimmedBases += *c;
        }
        if (!cutFront[i] && !cutBack[i]) continue;
        seqs[i] = seqs[i].substr(cutFront[i], seqs[i].size() - cutFront[i] - cutBack[i]);
    }
    return st;
}

void logTrimCounters(const TrimConfig& cfg, const TrimStats& st) {
    std::fprintf(stderr,
                 "[boundarytrim] enabled=%d trims=%zu bp_trimmed=%zu boundary_bp_kept=%zu "
                 "equal_len_flips=%zu skipped_short=%zu\n",
                 cfg.boundarySafe ? 1 : 0, st.trimmed, st.trimmedBases, st.boundaryBpKept,
                 st.equalLenFlips, st.skippedByBoundary);
    std::fprintf(stderr,
                 "[trimguard] enabled=%d components=%zu nonbipartite=%zu guard_differs=%zu "
                 "guard_applied=%zu evidence_refused=%zu evidence_switched=%zu ends_rescued=%zu "
                 "release_cuts_avoided=%zu\n",
                 cfg.copyGuard ? 1 : 0, st.components, st.nonBipartite, st.guardDiffers,
                 st.guardApplied, st.evidenceRefused, st.evidenceSwitched, st.endsRescued,
                 st.releaseCutsAvoided);
}

// ---- T22 stage 2: read-pair evidence at contig ends ------------------------------------------
std::vector<EndEvidence> gatherEndEvidence(const std::vector<std::string>& seqs,
                                           const std::vector<EndQuery>& queries,
                                           const SequenceStore& reads, const InsertModel& insert,
                                           int threads) {
    std::vector<EndEvidence> out(queries.size());
    if (queries.empty() || !insert.usable || reads.pairCount() == 0 || insert.mean <= 0) return out;
    if (threads <= 0) threads = 1;
    constexpr int k = 31;
    constexpr int kStride = 3;             // every 3rd read k-mer is probed
    constexpr int kSlack = 30;
    const double sd = std::max(1.0, insert.stddev);
    const long insMin = std::max<long>(insert.minPlausible, static_cast<long>(insert.mean - 3 * sd));
    const long insMax = std::max<long>(insMin + 1, std::min<long>(insert.maxPlausible > 0 ? insert.maxPlausible
                                                                                          : static_cast<long>(insert.mean + 3 * sd),
                                                                  static_cast<long>(insert.mean + 3 * sd)));
    const long reach = insMax + static_cast<long>(reads.maxReadLength());

    // Each query as an oriented sequence C with its segment at the END: C = F S.
    struct Q { std::string c; size_t n = 0, L = 0; };
    std::vector<Q> qs(queries.size());
    for (size_t i = 0; i < queries.size(); ++i) {
        const std::string& s = seqs[queries[i].seq];
        qs[i].c = queries[i].front ? reverseComplement(s) : s;
        qs[i].n = qs[i].c.size();
        qs[i].L = std::min(queries[i].segLen, qs[i].n);
    }
    // Key table: flank-window k-mers (anchors, must be unique in the assembly) and segment
    // k-mers (for depth and for placing mates).
    struct Key { uint32_t query = 0; int32_t pos = 0; uint8_t strand = 0; uint8_t kind = 0; uint32_t occ = 0; };
    // kind bit 1 = flank anchor, bit 2 = segment k-mer
    std::unordered_map<Kmer, uint32_t, KmerHasher> keyId;
    std::vector<Key> keys;
    auto forKmers = [&](const std::string& s, size_t from, size_t to, auto&& fn) {
        Kmer fwd = 0, rc = 0;
        int valid = 0;
        for (size_t p = from; p < to; ++p) {
            const int b = baseCode(s[p]);
            if (b < 0) { valid = 0; continue; }
            fwd = pushBack(fwd, b, k);
            rc = pushFrontRc(rc, b, k);
            if (++valid < k) continue;
            fn(fwd < rc ? fwd : rc, static_cast<int32_t>(p + 1 - k), fwd < rc ? 0 : 1);
        }
    };
    std::vector<std::vector<std::pair<Kmer, int32_t>>> segKmers(qs.size());   // forward k-mers of S (+ margin)
    for (size_t i = 0; i < qs.size(); ++i) {
        const Q& q = qs[i];
        if (q.n < static_cast<size_t>(k) || q.L == 0 || q.L >= q.n) continue;
        const size_t fEnd = q.n - q.L;
        const size_t fBeg = fEnd > static_cast<size_t>(reach) ? fEnd - static_cast<size_t>(reach) : 0;
        forKmers(q.c, fBeg, fEnd, [&](const Kmer& km, int32_t pos, int strand) {
            auto it = keyId.find(km);
            if (it == keyId.end()) {
                keyId.emplace(km, static_cast<uint32_t>(keys.size()));
                keys.push_back({static_cast<uint32_t>(i), pos, static_cast<uint8_t>(strand), 1, 0});
            } else {
                keys[it->second].kind |= 1;
                keys[it->second].query = UINT32_MAX;       // two anchor sites: not an anchor
            }
        });
        const size_t margin = reads.maxReadLength() + static_cast<size_t>(std::max(0.0, insert.mean));
        const size_t sBeg = fEnd > margin ? fEnd - margin : 0;
        Kmer fwd = 0;
        int valid = 0;
        for (size_t p = sBeg; p < q.n; ++p) {
            const int b = baseCode(q.c[p]);
            if (b < 0) { valid = 0; continue; }
            fwd = pushBack(fwd, b, k);
            if (++valid < k) continue;
            const int32_t pos = static_cast<int32_t>(p + 1 - k);
            segKmers[i].push_back({fwd, pos});
            if (static_cast<size_t>(pos) < fEnd) continue;
            const Kmer can = canonical(fwd, k);
            auto it = keyId.find(can);
            if (it == keyId.end()) {
                keyId.emplace(can, static_cast<uint32_t>(keys.size()));
                keys.push_back({UINT32_MAX, 0, 0, 2, 0});
            } else {
                keys[it->second].kind |= 2;
            }
        }
    }
    if (keys.empty()) return out;
    // Uniqueness in the assembly: count every occurrence of every key over all sequences.
    for (const std::string& s : seqs) {
        forKmers(s, 0, s.size(), [&](const Kmer& km, int32_t, int) {
            auto it = keyId.find(km);
            if (it != keyId.end()) ++keys[it->second].occ;
        });
    }
    // Per query: forward-strand position index of the segment region.
    std::vector<std::unordered_multimap<Kmer, int32_t, KmerHasher>> segIndex(qs.size());
    for (size_t i = 0; i < qs.size(); ++i)
        for (const auto& kp : segKmers[i]) segIndex[i].emplace(kp.first, kp.second);

    struct Acc { std::vector<size_t> consistent, inconsistent; std::vector<uint64_t> hits; };
    std::vector<Acc> acc(static_cast<size_t>(threads));
    for (Acc& a : acc) {
        a.consistent.assign(qs.size(), 0);
        a.inconsistent.assign(qs.size(), 0);
        a.hits.assign(keys.size(), 0);
    }
    auto worker = [&](int tid) {
        Acc& a = acc[static_cast<size_t>(tid)];
        std::string mate;
        for (size_t r = static_cast<size_t>(tid); r < reads.size(); r += static_cast<size_t>(threads)) {
            const uint32_t len = reads.length(r);
            if (len < static_cast<uint32_t>(k)) continue;
            // Probe every kStride-th k-mer: depth for segment/flank keys, votes for anchors.
            struct Vote { uint32_t q; int32_t start; uint8_t fwd; int n; };
            Vote votes[8];
            int nv = 0;
            {
                Kmer fwd = 0, rcv = 0;
                int run = 0;
                for (uint32_t p = 0; p < len; ++p) {
                    // Read through correction masks (see rawBaseAt): the corrector masks exactly
                    // the stretches where a read stops agreeing with the trusted k-mers, which is
                    // what a read from the true locus does at a wrongly walked junction. Masked,
                    // that disagreement would vanish from the evidence instead of counting.
                    const int b = reads.rawBaseAt(r, p);
                    fwd = pushBack(fwd, b, k);
                    rcv = pushFrontRc(rcv, b, k);
                    if (++run < k) continue;
                    const uint32_t kp = p + 1 - static_cast<uint32_t>(k);
                    if (kp % kStride) continue;
                    const bool readIsCanon = fwd < rcv;
                    const Kmer can = readIsCanon ? fwd : rcv;
                    auto it = keyId.find(can);
                    if (it == keyId.end()) continue;
                    const Key& key = keys[it->second];
                    ++a.hits[it->second];
                    if (!(key.kind & 1) || key.query == UINT32_MAX || key.occ != 1) continue;
                    // The read's own k-mer orientation relative to C decides its strand on C.
                    const bool onCfwd = (readIsCanon ? 0 : 1) == key.strand;
                    const int32_t start = onCfwd ? key.pos - static_cast<int32_t>(kp)
                                                 : key.pos - static_cast<int32_t>(len - static_cast<uint32_t>(k) - kp);
                    bool merged = false;
                    for (int v = 0; v < nv && !merged; ++v)
                        if (votes[v].q == key.query && votes[v].start == start && votes[v].fwd == onCfwd) {
                            ++votes[v].n;
                            merged = true;
                        }
                    if (!merged && nv < 8) votes[nv++] = {key.query, start, static_cast<uint8_t>(onCfwd), 1};
                }
            }
            if (!reads.hasMate(r) || nv == 0) continue;
            int best = -1;
            for (int v = 0; v < nv; ++v)
                if (votes[v].n >= 2 && (best < 0 || votes[v].n > votes[best].n)) best = v;
            if (best < 0 || !votes[best].fwd) continue;           // must point towards the segment
            const uint32_t qi = votes[best].q;
            const Q& q = qs[qi];
            const long p0 = votes[best].start;
            const size_t m = reads.mateOf(r);
            const long lm = reads.length(m);
            if (lm < k) continue;
            const long winStart = p0 + insMin - lm;               // earliest start of rc(mate) on C
            const long winEnd = p0 + insMax;                      // latest end
            // Pairs whose mate is expected (at the mean insert) to start inside the segment, or
            // within one insert before it, speak to the end: a walk that left the right locus a
            // little before the exact overlap begins (a chimeric junction inside the flank, which
            // the release trim happened to leave too short for QUAST to see) shows up as mates
            // that do not match the contig there.
            if (p0 + static_cast<long>(insert.mean) - lm <
                static_cast<long>(q.n - q.L) - static_cast<long>(insert.mean)) continue;
            // Place rc(mate) on C: vote over its k-mers' positions in the segment index.
            mate.resize(static_cast<size_t>(lm));
            for (uint32_t p = 0; p < static_cast<uint32_t>(lm); ++p) mate[p] = codeBase(reads.rawBaseAt(m, p));
            const std::string mr = reverseComplement(mate);
            std::map<long, int> place;
            int valid = 0;
            Kmer fwd = 0;
            int run = 0;
            for (size_t p = 0; p < mr.size(); ++p) {
                const int b = baseCode(mr[p]);
                if (b < 0) { run = 0; continue; }
                fwd = pushBack(fwd, b, k);
                if (++run < k) continue;
                if (((p + 1 - k) % kStride) != 0) continue;
                ++valid;
                auto range = segIndex[qi].equal_range(fwd);
                for (auto it = range.first; it != range.second; ++it)
                    ++place[static_cast<long>(it->second) - static_cast<long>(p + 1 - k)];
            }
            if (valid == 0) continue;
            bool consistent = false;
            for (const auto& kv : place) {
                if (kv.second * 2 < valid) continue;
                if (kv.first >= winStart - kSlack && kv.first + lm <= winEnd + kSlack) consistent = true;
            }
            if (consistent) ++a.consistent[qi];
            else if (winEnd <= static_cast<long>(q.n)) ++a.inconsistent[qi];
        }
    };
    std::vector<std::thread> pool;
    for (int t = 0; t < threads; ++t) pool.emplace_back(worker, t);
    for (auto& th : pool) th.join();

    std::vector<uint64_t> hits(keys.size(), 0);
    for (const Acc& a : acc)
        for (size_t i = 0; i < keys.size(); ++i) hits[i] += a.hits[i];
    for (size_t i = 0; i < qs.size(); ++i) {
        for (const Acc& a : acc) {
            out[i].consistent += a.consistent[i];
            out[i].inconsistent += a.inconsistent[i];
        }
    }
    // Copy ratio: mean probe hits per segment k-mer over mean hits per unique flank anchor.
    std::vector<double> sSum(qs.size(), 0), sN(qs.size(), 0), fSum(qs.size(), 0), fN(qs.size(), 0);
    for (size_t i = 0; i < keys.size(); ++i) {
        const Key& key = keys[i];
        if ((key.kind & 1) && key.query != UINT32_MAX && key.occ == 1) {
            fSum[key.query] += static_cast<double>(hits[i]);
            fN[key.query] += 1;
        }
    }
    for (size_t qi = 0; qi < qs.size(); ++qi) {
        const size_t fEnd = qs[qi].n - std::min(qs[qi].n, qs[qi].L);
        for (const auto& kp : segKmers[qi]) {
            if (static_cast<size_t>(kp.second) < fEnd) continue;
            auto it = keyId.find(canonical(kp.first, k));
            if (it == keyId.end()) continue;
            sSum[qi] += static_cast<double>(hits[it->second]);
            sN[qi] += 1;
        }
        if (sN[qi] > 0 && fN[qi] > 0 && fSum[qi] > 0)
            out[qi].copyRatio = (sSum[qi] / sN[qi]) / (fSum[qi] / fN[qi]);
    }
    return out;
}

// ---- T26: AGP and GFA paths that describe contigs.fasta ------------------------------------
bool writeAgpV2(const std::string& path, const std::vector<std::string>& scaffolds,
                const std::vector<std::string>& scaffoldNames, const std::vector<EmitPiece>& pieces,
                const std::string& gapEvidence, AgpGfaStats& st, std::string& error) {
    // Every written record, so bases a piece lost to the trim or the dedup can be pointed at
    // in the record that still holds them.
    std::string text;
    std::vector<size_t> recStart, recPiece;
    for (size_t i = 0; i < pieces.size(); ++i) {
        if (pieces[i].drop) continue;
        text.push_back('\x01');
        recStart.push_back(text.size());
        recPiece.push_back(i);
        text += pieces[i].seq;
    }
    struct Hit { bool found = false; size_t piece = 0, off = 0; bool minus = false; };
    auto locate = [&](const std::string& x) {
        Hit h;
        if (x.empty() || text.empty()) return h;
        size_t q = text.find(x);
        if (q == std::string::npos) {
            q = text.find(reverseComplement(x));
            if (q == std::string::npos) return h;
            h.minus = true;
        }
        const size_t r = static_cast<size_t>(std::upper_bound(recStart.begin(), recStart.end(), q) -
                                             recStart.begin()) - 1;
        h.found = true;
        h.piece = recPiece[r];
        h.off = q - recStart[r];
        return h;
    };
    std::FILE* agp = std::fopen(path.c_str(), "w");
    if (!agp) { error = "cannot write " + path; return false; }
    std::fprintf(agp, "##agp-version\t2.1\n");
    std::fprintf(agp, "# TesserACT -- order and orientation of contigs.fasta records within scaffolds.fasta.\n");
    std::fprintf(agp, "# An N row is an adjacency the assembly asserts without\n"
                      "# sequence; it is not evidence of the intervening bases.\n");
    std::vector<std::vector<size_t>> byScaffold(scaffolds.size());
    for (size_t i = 0; i < pieces.size(); ++i) byScaffold[pieces[i].scaffold].push_back(i);
    for (size_t s = 0; s < scaffolds.size(); ++s) {
        const std::string& sq = scaffolds[s];
        const char* obj = scaffoldNames[s].c_str();
        size_t part = 0, pos = 0, unlisted = 0;
        auto nRow = [&](size_t b, size_t e) {
            std::fprintf(agp, "%s\t%zu\t%zu\t%zu\tN\t%zu\tscaffold\tyes\t%s\n", obj, b + 1, e, ++part,
                         e - b, gapEvidence.c_str());
            if (gapEvidence == "paired-ends") ++st.nPaired;
            else ++st.nUnspecified;
        };
        auto wRegion = [&](size_t b, size_t e) {   // bases no written record holds as-is here
            if (e <= b) return;
            const Hit h = locate(sq.substr(b, e - b));
            if (h.found) {
                std::fprintf(agp, "%s\t%zu\t%zu\t%zu\tW\t%s\t%zu\t%zu\t%c\n", obj, b + 1, e, ++part,
                             pieces[h.piece].name.c_str(), h.off + 1, h.off + (e - b), h.minus ? '-' : '+');
                ++st.w;
                ++st.wPartner;
            } else {
                std::fprintf(agp, "%s\t%zu\t%zu\t%zu\tW\t%s_unlisted_%zu\t1\t%zu\t+\n", obj, b + 1, e, ++part,
                             obj, ++unlisted, e - b);
                ++st.w;
                ++st.wUnlisted;
            }
        };
        for (size_t idx : byScaffold[s]) {
            const EmitPiece& p = pieces[idx];
            if (p.objBegin > pos) nRow(pos, p.objBegin);
            if (p.drop) {
                wRegion(p.objBegin, p.objEnd);
            } else {
                const size_t f = p.cutFrontSeq.size(), b = p.cutBackSeq.size();
                wRegion(p.objBegin, p.objBegin + f);
                std::fprintf(agp, "%s\t%zu\t%zu\t%zu\tW\t%s\t1\t%zu\t+\n", obj, p.objBegin + f + 1,
                             p.objEnd - b, ++part, p.name.c_str(), p.seq.size());
                ++st.w;
                wRegion(p.objEnd - b, p.objEnd);
            }
            pos = p.objEnd;
        }
        if (!byScaffold[s].empty() && pos < sq.size()) nRow(pos, sq.size());
    }
    const bool ok = std::ferror(agp) == 0;
    const bool closed = std::fclose(agp) == 0;
    if (!ok || !closed) { error = "write failed on " + path; return false; }
    return true;
}

std::vector<GfaPath> contigPathsV2(const std::vector<std::string>& scaffolds,
                                   const std::vector<GfaPath>& paths, const UnitigGraph& g,
                                   const std::vector<EmitPiece>& pieces, AgpGfaStats& st) {
    std::vector<GfaPath> out;
    std::vector<std::vector<size_t>> byScaffold(scaffolds.size());
    for (size_t i = 0; i < pieces.size(); ++i) byScaffold[pieces[i].scaffold].push_back(i);
    for (size_t s = 0; s < scaffolds.size(); ++s) {
        if (s >= paths.size() || paths[s].oriented.empty()) continue;
        size_t written = 0;
        for (size_t idx : byScaffold[s]) written += pieces[idx].drop ? 0 : 1;
        const std::vector<PathSegment> segs = pathSegments(paths[s], g);
        const std::vector<size_t> pieceOf = mapSegmentsToPieces(scaffolds[s], segs, g);
        if (pieceOf.empty() || pieceOf.back() + 1 != byScaffold[s].size()) {
            st.pDropped += written;
            continue;
        }
        for (size_t j = 0; j < byScaffold[s].size(); ++j) {
            const EmitPiece& p = pieces[byScaffold[s][j]];
            if (p.drop) continue;
            // Only a piece that is exactly one gap-free stretch of the walk, spelled base for base,
            // gets a path. A closed gap has no graph step (its bases came from the reads), and the
            // polisher, the trim and the release gap-flank rule all change the spelling.
            size_t nSeg = 0, seg = 0;
            for (size_t x = 0; x < pieceOf.size(); ++x)
                if (pieceOf[x] == j) { ++nSeg; seg = x; }
            if (nSeg != 1 || spellWalk(segs[seg].oriented, g) != p.seq) { ++st.pDropped; continue; }
            GfaPath gp;
            gp.name = p.name;
            gp.oriented = segs[seg].oriented;
            gp.gaps.assign(gp.oriented.size(), 0);
            out.push_back(std::move(gp));
            ++st.pEmitted;
        }
    }
    return out;
}

void logAgpGfa(const AgpGfaStats& st) {
    std::fprintf(stderr,
                 "[agpgfa] enabled=%d W=%zu W_partner=%zu W_unlisted=%zu N_paired=%zu N_unspecified=%zu "
                 "P_emitted=%zu P_dropped=%zu\n",
                 st.enabled ? 1 : 0, st.w, st.wPartner, st.wUnlisted, st.nPaired, st.nUnspecified,
                 st.pEmitted, st.pDropped);
}

// Contig-level statistics (T17): the records contigs.fasta holds -- split at every N, after the
// terminal-overlap trim and any split post-processing -- and the number of N-runs, of any
// length, in the scaffolds (= the N rows of scaffolds.agp). The release counted runs of 10+ N
// in the scaffolds before the trim, which described no file the run wrote.
void computeContigStats(const std::vector<std::string>& records,
                               const std::vector<std::string>& scaffolds, AssemblyReport& rep) {
    auto isN = [](char c) { return c == 'N' || c == 'n'; };
    rep.scaffoldGaps = 0;
    for (const std::string& s : scaffolds) {
        for (size_t i = 0; i < s.size();) {
            if (!isN(s[i])) { ++i; continue; }
            ++rep.scaffoldGaps;
            while (i < s.size() && isN(s[i])) ++i;
        }
    }
    std::vector<size_t> lens;
    lens.reserve(records.size());
    for (const std::string& r : records) lens.push_back(r.size());
    rep.contigPieces = lens.size();
    rep.contigTotal = 0;
    rep.contigLargest = 0;
    for (size_t l : lens) {
        rep.contigTotal += l;
        rep.contigLargest = std::max(rep.contigLargest, l);
    }
    std::sort(lens.begin(), lens.end(), std::greater<size_t>());
    size_t acc = 0;
    rep.contigN50 = 0;
    for (size_t l : lens) {
        acc += l;
        if (acc * 2 >= rep.contigTotal) { rep.contigN50 = l; break; }
    }
}

}  // namespace ts
