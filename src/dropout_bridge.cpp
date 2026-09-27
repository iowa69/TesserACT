// Dropout bridging; see dropout_bridge.h. Default OFF (TESSERACT_DROPOUT_BRIDGE=1).
#include "dropout_bridge.h"

#include <algorithm>
#include <cstdio>
#include <map>
#include <set>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

#include "kmer.h"
#include "util.h"

namespace ts {

namespace {

constexpr int kA = 31;                 // anchor length; one 64-bit word
constexpr int kW = 60;                 // terminal window
constexpr int kMaxMismatches = 2;      // read vs both windows, in total
constexpr uint64_t kAMask = (uint64_t{1} << (2 * kA)) - 1;
constexpr uint64_t kEmptyKey = ~uint64_t{0};   // never a 62-bit A-mer

// Verbatim copy of the anonymous-namespace lowComplexity() in graph.cpp, which this
// file cannot reach; the design names that exact test.
bool lowComplexity(const std::string& s) {
    if (s.empty()) return true;
    size_t counts[4] = {0, 0, 0, 0};
    size_t total = 0;
    for (char c : s) {
        switch (c) {
            case 'A': ++counts[0]; ++total; break;
            case 'C': ++counts[1]; ++total; break;
            case 'G': ++counts[2]; ++total; break;
            case 'T': ++counts[3]; ++total; break;
            default: break;
        }
    }
    if (total == 0) return true;
    size_t best = 0;
    for (size_t c : counts) best = std::max(best, c);
    const double limit = total >= 40 ? 0.90 : (total >= 20 ? 0.82 : 0.75);
    return static_cast<double>(best) > limit * static_cast<double>(total);
}

// Open-addressing map from a canonical A-mer to a small index, with a bit prefilter so
// the ~all read A-mers that miss cost one cache-resident probe.
class AmerMap {
public:
    void build(const std::vector<std::pair<uint64_t, uint32_t>>& items) {
        size_t cap = 1024;
        while (cap < items.size() * 4) cap <<= 1;
        keys_.assign(cap, kEmptyKey);
        vals_.assign(cap, 0);
        mask_ = cap - 1;
        size_t bits = size_t{1} << 16;
        while (bits < items.size() * 64) bits <<= 1;
        filter_.assign(bits / 64, 0);
        fmask_ = bits - 1;
        for (const auto& kv : items) {
            const uint64_t h = mix64(kv.first);
            const uint64_t fb = (h >> 24) & fmask_;
            filter_[fb >> 6] |= uint64_t{1} << (fb & 63);
            size_t i = h & mask_;
            while (keys_[i] != kEmptyKey && keys_[i] != kv.first) i = (i + 1) & mask_;
            keys_[i] = kv.first;
            vals_[i] = kv.second;
        }
    }
    const uint32_t* find(uint64_t key) const {
        if (keys_.empty()) return nullptr;
        const uint64_t h = mix64(key);
        const uint64_t fb = (h >> 24) & fmask_;
        if (!((filter_[fb >> 6] >> (fb & 63)) & 1)) return nullptr;
        for (size_t i = h & mask_;; i = (i + 1) & mask_) {
            if (keys_[i] == key) return &vals_[i];
            if (keys_[i] == kEmptyKey) return nullptr;
        }
    }

private:
    std::vector<uint64_t> keys_;
    std::vector<uint32_t> vals_;
    std::vector<uint64_t> filter_;
    size_t mask_ = 0;
    uint64_t fmask_ = 0;
};

// Rolls canonical A-mers over a string; fn(fwd, rc, startPos). Non-ACGT resets.
template <typename Fn>
void forEachAmer(const std::string& s, Fn&& fn) {
    uint64_t fwd = 0, rev = 0;
    int valid = 0;
    for (size_t p = 0; p < s.size(); ++p) {
        const int c = baseCode(s[p]);
        if (c < 0) { valid = 0; fwd = rev = 0; continue; }
        fwd = ((fwd << 2) | static_cast<uint64_t>(c)) & kAMask;
        rev = (rev >> 2) | (static_cast<uint64_t>(3 - c) << (2 * (kA - 1)));
        if (++valid >= kA) fn(fwd, rev, p + 1 - static_cast<size_t>(kA));
    }
}

// One eligible dead end (u, e). OUT is the sequence running out of the end (u.seq for
// e = 1, rc(u.seq) for e = 0); IN = rc(OUT) is what a traversal entering the end reads.
struct End {
    uint32_t u = 0;
    uint8_t e = 0;
    double cov = 0;
    size_t len = 0;
    int window = 0;          // min(W, |tail|)
    std::string tail;        // OUT.last(T), T = min(|u|, max(W, K-1))
    std::string head;        // IN.first(T) = rc(tail)
};

// A window A-mer: which end, the distance from its first base to the tip in OUT
// coordinates (A..W), and whether its canonical form is the OUT-strand sequence.
struct Anchor {
    uint32_t end = 0;
    uint16_t off = 0;
    bool canonIsOut = false;
};

struct Hit {
    uint32_t end;
    int32_t implied;     // exit: read position of OUT's tip; enter: read position of IN's start
    int32_t pos;         // A-mer start in the read
    uint8_t enter;
};

// One read's vote for "end a runs into end b with `g` bases between them". Direction-
// normalised so a <= b: a read that natively ran b -> a is recorded reverse-complemented.
struct Nom {
    uint32_t a = 0, b = 0;
    int32_t g = 0;
    uint32_t read = 0;
    uint64_t pairId = 0;
    int32_t key = 0;         // native 5' end's distance to the native exit tip
    uint8_t strand = 0;      // 0 = native a -> b, 1 = native b -> a
    uint16_t rawVotes = 0;   // gap bases that the corrector had masked (read raw)
    std::string gap;         // normalised gap bases, g > 0 only
};

// Same effect as UnitigGraph::addLink (private): record the link on both sides.
void linkEnds(UnitigGraph& g, uint32_t u, int ue, uint32_t v, int ve) {
    const Link a{v, static_cast<uint8_t>(ve)};
    auto& lu = g.nodes[u].ends[ue];
    if (std::find(lu.begin(), lu.end(), a) == lu.end()) lu.push_back(a);
    const Link b{u, static_cast<uint8_t>(ue)};
    auto& lv = g.nodes[v].ends[ve];
    if (std::find(lv.begin(), lv.end(), b) == lv.end()) lv.push_back(b);
}

}  // namespace

DropoutBridgeStats bridgeDropouts(UnitigGraph& g, const SequenceStore& reads, int threads,
                                  uint32_t minReads, bool trace) {
    util::Timer timer;
    const long long rss0 = util::currentMemoryBytes();
    long long rssPeak = rss0;
    auto sampleRss = [&] { rssPeak = std::max(rssPeak, util::currentMemoryBytes()); };

    DropoutBridgeStats st;
    st.enabled = true;
    const int K = g.k();
    st.k = K;
    if (threads <= 0) threads = 1;
    if (minReads < 2) minReads = 2;
    st.minReads = minReads;
    const double med = g.medianCoverage();
    st.medianCoverage = med;
    const int gMin = -(K - 2);
    const int gMax = static_cast<int>(reads.maxReadLength()) - 2 * kW;

    // ---- 1. eligible dead ends ------------------------------------------------------
    std::vector<End> ends;
    const size_t tailWant = static_cast<size_t>(std::max(kW, K - 1));
    for (uint32_t u = 0; u < g.nodes.size(); ++u) {
        const Unitig& U = g.nodes[u];
        if (U.deleted) continue;
        for (int e = 0; e < 2; ++e) {
            if (!U.ends[e].empty()) continue;
            ++st.deadEnds;
            if (U.seq.size() < static_cast<size_t>(2 * K)) continue;
            if (!(U.coverage >= 0.4 * med && U.coverage <= 1.6 * med)) continue;
            End en;
            en.u = u;
            en.e = static_cast<uint8_t>(e);
            en.cov = U.coverage;
            en.len = U.seq.size();
            const size_t T = std::min(U.seq.size(), tailWant);
            en.tail = e == 1 ? U.seq.substr(U.seq.size() - T)
                             : reverseComplement(U.seq.substr(0, T));
            en.head = reverseComplement(en.tail);
            en.window = static_cast<int>(std::min<size_t>(kW, T));
            ends.push_back(std::move(en));
        }
    }

    // ---- 2. anchor index --------------------------------------------------------------
    // Every canonical A-mer of every window, then two uniqueness filters: exactly one
    // occurrence over all windows, and exactly one occurrence in the live graph.
    std::vector<Anchor> anchors;
    std::vector<uint64_t> anchorKey;
    std::vector<uint32_t> windowCount;
    std::vector<uint32_t> graphCount;
    AmerMap table;
    {
        std::unordered_map<uint64_t, uint32_t> idx;
        for (uint32_t i = 0; i < ends.size(); ++i) {
            const End& en = ends[i];
            const std::string& t = en.tail;
            const size_t T = t.size();
            if (en.window < kA) continue;
            const std::string win = t.substr(T - static_cast<size_t>(en.window));
            forEachAmer(win, [&](uint64_t fwd, uint64_t rev, size_t p) {
                const uint64_t canon = fwd < rev ? fwd : rev;
                auto it = idx.find(canon);
                if (it != idx.end()) { ++windowCount[it->second]; return; }
                const uint32_t id = static_cast<uint32_t>(anchors.size());
                idx.emplace(canon, id);
                Anchor an;
                an.end = i;
                an.off = static_cast<uint16_t>(static_cast<size_t>(en.window) - p);
                an.canonIsOut = fwd < rev;
                anchors.push_back(an);
                anchorKey.push_back(canon);
                windowCount.push_back(1);
                graphCount.push_back(0);
            });
        }
        std::vector<std::pair<uint64_t, uint32_t>> items;
        items.reserve(anchors.size());
        for (uint32_t id = 0; id < anchors.size(); ++id) items.push_back({anchorKey[id], id});
        AmerMap all;
        all.build(items);
        for (const Unitig& U : g.nodes) {
            if (U.deleted) continue;
            forEachAmer(U.seq, [&](uint64_t fwd, uint64_t rev, size_t) {
                const uint32_t* id = all.find(fwd < rev ? fwd : rev);
                if (id && graphCount[*id] < 0xffffffffu) ++graphCount[*id];
            });
        }
        std::vector<uint32_t> keptPerEnd(ends.size(), 0);
        std::vector<Anchor> kept;
        items.clear();
        for (uint32_t id = 0; id < anchors.size(); ++id) {
            if (windowCount[id] != 1 || graphCount[id] != 1) { ++st.anchorsDropped; continue; }
            items.push_back({anchorKey[id], static_cast<uint32_t>(kept.size())});
            kept.push_back(anchors[id]);
            ++keptPerEnd[anchors[id].end];
        }
        anchors.swap(kept);
        table.build(items);
        st.anchors = anchors.size();
        for (uint32_t i = 0; i < ends.size(); ++i) {
            if (keptPerEnd[i] > 0) ++st.eligible;
            else ++st.endsNoAnchor;
        }
    }
    sampleRss();

    // ---- 3. read scan -------------------------------------------------------------------
    std::vector<Nom> noms;
    if (!anchors.empty() && reads.size() > 0) {
        const size_t pairedReads = reads.pairedReads();
        std::vector<std::vector<Nom>> local(static_cast<size_t>(threads));
        std::vector<size_t> localHit(static_cast<size_t>(threads), 0);
        std::vector<size_t> localMaskedHit(static_cast<size_t>(threads), 0);
        auto worker = [&](int tid) {
            auto& out = local[static_cast<size_t>(tid)];
            std::vector<Hit> hits;
            struct Grp { uint32_t end; bool enter; int32_t implied, first, last; };
            std::vector<Grp> exits, enters;
            for (size_t r = static_cast<size_t>(tid); r < reads.size();
                 r += static_cast<size_t>(threads)) {
                const uint32_t len = reads.length(r);
                if (len < static_cast<uint32_t>(kA)) continue;
                hits.clear();
                // Rolled THROUGH the corrector's mask (rawBaseAt), as forEachKmerRaw does for the
                // gap filler. correctReads masks from the point where correction stalls to the
                // read's end, so a read that spans a dropout has the far flank masked: with a
                // masked roll the 2x301 sentinel gave 3 nominations from 25,804 anchor-hitting
                // reads. An anchor is an exact 31-mer match to sequence unique in the graph, so
                // reading raw bases here admits no base the graph does not already carry.
                uint64_t fwd = 0, rev = 0;
                for (uint32_t p = 0; p < len; ++p) {
                    const int c = reads.rawBaseAt(r, p);
                    fwd = ((fwd << 2) | static_cast<uint64_t>(c)) & kAMask;
                    rev = (rev >> 2) | (static_cast<uint64_t>(3 - c) << (2 * (kA - 1)));
                    if (p + 1 < static_cast<uint32_t>(kA)) continue;
                    const bool fwdCanon = fwd < rev;
                    const uint32_t* id = table.find(fwdCanon ? fwd : rev);
                    if (!id) continue;
                    const Anchor& an = anchors[*id];
                    const int32_t start = static_cast<int32_t>(p) + 1 - kA;
                    for (int32_t x = start; x <= static_cast<int32_t>(p); ++x) {
                        if (reads.baseAt(r, static_cast<uint32_t>(x)) < 0) {
                            ++localMaskedHit[static_cast<size_t>(tid)];
                            break;
                        }
                    }
                    if (fwdCanon == an.canonIsOut) {
                        hits.push_back({an.end, start + an.off, start, 0});           // exits
                    } else {
                        hits.push_back({an.end, start - (an.off - kA), start, 1});    // enters
                    }
                }
                if (hits.empty()) continue;
                ++localHit[static_cast<size_t>(tid)];
                if (hits.size() < 2) continue;
                std::sort(hits.begin(), hits.end(), [](const Hit& x, const Hit& y) {
                    if (x.end != y.end) return x.end < y.end;
                    if (x.enter != y.enter) return x.enter < y.enter;
                    return x.pos < y.pos;
                });
                exits.clear();
                enters.clear();
                for (size_t lo = 0; lo < hits.size();) {
                    size_t hi = lo;
                    while (hi < hits.size() && hits[hi].end == hits[lo].end) ++hi;
                    // hits[lo..hi) share an end; sorted exits first, then enters.
                    const bool anyExit = hits[lo].enter == 0;
                    const bool anyEnter = hits[hi - 1].enter == 1;
                    bool consistent = !(anyExit && anyEnter);
                    for (size_t x = lo + 1; consistent && x < hi; ++x) {
                        if (hits[x].implied != hits[lo].implied) consistent = false;
                    }
                    if (consistent) {
                        Grp gr{hits[lo].end, anyEnter, hits[lo].implied, hits[lo].pos,
                               hits[hi - 1].pos};
                        (anyEnter ? enters : exits).push_back(gr);
                    }
                    lo = hi;
                }
                for (const Grp& X : exits) {
                    for (const Grp& Y : enters) {
                        const int32_t tipPos = X.implied;
                        const int32_t inStart = Y.implied;
                        const int32_t gap = inStart - tipPos;
                        if (gap < gMin || gap > gMax) continue;
                        if (tipPos <= 0 || tipPos > static_cast<int32_t>(len) || inStart < 0 ||
                            inStart >= static_cast<int32_t>(len)) continue;
                        const End& EI = ends[X.end];
                        const End& EJ = ends[Y.end];
                        const int32_t Ti = static_cast<int32_t>(EI.tail.size());
                        const int32_t Tj = static_cast<int32_t>(EJ.head.size());
                        int mm = 0;
                        // OUT_i from the first exit hit to the tip.
                        for (int32_t x = X.first; x < tipPos && mm <= kMaxMismatches; ++x) {
                            const int32_t t = Ti - (tipPos - x);
                            if (t < 0) { mm = kMaxMismatches + 1; break; }
                            if (codeBase(reads.rawBaseAt(r, static_cast<uint32_t>(x))) !=
                                EI.tail[static_cast<size_t>(t)]) ++mm;
                        }
                        // IN_j from its start (or the tip, when the flanks overlap and
                        // OUT_i already covered those read bases) to the end of the last hit.
                        const int32_t yEnd = std::min<int32_t>(Y.last + kA, static_cast<int32_t>(len));
                        for (int32_t x = std::max(inStart, tipPos); x < yEnd && mm <= kMaxMismatches; ++x) {
                            const int32_t t = x - inStart;
                            if (t >= Tj) { mm = kMaxMismatches + 1; break; }
                            if (codeBase(reads.rawBaseAt(r, static_cast<uint32_t>(x))) !=
                                EJ.head[static_cast<size_t>(t)]) ++mm;
                        }
                        if (mm > kMaxMismatches) continue;
                        Nom n;
                        n.g = gap;
                        n.read = static_cast<uint32_t>(r);
                        n.pairId = r < pairedReads ? static_cast<uint64_t>(r >> 1)
                                                   : static_cast<uint64_t>(pairedReads / 2 + (r - pairedReads));
                        n.key = tipPos;
                        std::string gs;
                        if (gap > 0) {
                            gs.resize(static_cast<size_t>(gap));
                            for (int32_t x = tipPos; x < inStart; ++x) {
                                gs[static_cast<size_t>(x - tipPos)] =
                                    codeBase(reads.rawBaseAt(r, static_cast<uint32_t>(x)));
                                if (reads.baseAt(r, static_cast<uint32_t>(x)) < 0 &&
                                    n.rawVotes < 0xffff) ++n.rawVotes;
                            }
                        }
                        if (X.end <= Y.end) {
                            n.a = X.end; n.b = Y.end; n.strand = 0; n.gap = std::move(gs);
                        } else {
                            n.a = Y.end; n.b = X.end; n.strand = 1; n.gap = reverseComplement(gs);
                        }
                        out.push_back(std::move(n));
                    }
                }
            }
        };
        std::vector<std::thread> pool;
        pool.reserve(static_cast<size_t>(threads));
        for (int t = 0; t < threads; ++t) pool.emplace_back(worker, t);
        for (auto& th : pool) th.join();
        for (int t = 0; t < threads; ++t) {
            st.readsHit += localHit[static_cast<size_t>(t)];
            st.maskedHits += localMaskedHit[static_cast<size_t>(t)];
            for (Nom& n : local[static_cast<size_t>(t)]) noms.push_back(std::move(n));
            std::vector<Nom>().swap(local[static_cast<size_t>(t)]);
        }
    }
    sampleRss();
    // Deterministic whatever the thread count.
    std::sort(noms.begin(), noms.end(), [](const Nom& x, const Nom& y) {
        if (x.a != y.a) return x.a < y.a;
        if (x.b != y.b) return x.b < y.b;
        return x.read < y.read;
    });
    st.nominations = noms.size();

    // ---- 4. acceptance ------------------------------------------------------------------
    struct Group { size_t lo, hi; std::vector<size_t> byPair; };
    std::vector<Group> groups;
    for (size_t lo = 0; lo < noms.size();) {
        size_t hi = lo;
        while (hi < noms.size() && noms[hi].a == noms[lo].a && noms[hi].b == noms[lo].b) ++hi;
        Group gr{lo, hi, {}};
        // Sorted by read id, hence by pair id: keep the first read of each pair.
        for (size_t x = lo; x < hi; ++x) {
            if (gr.byPair.empty() || noms[gr.byPair.back()].pairId != noms[x].pairId)
                gr.byPair.push_back(x);
        }
        groups.push_back(std::move(gr));
        lo = hi;
    }
    st.pairsNominated = groups.size();
    // A partner with >= 2 distinct pairs is a competitor, whatever else it fails.
    std::vector<uint32_t> strong(ends.size(), 0);
    for (const Group& gr : groups) {
        if (gr.byPair.size() < 2) continue;
        const Nom& n = noms[gr.lo];
        ++strong[n.a];
        if (n.b != n.a) ++strong[n.b];
    }

    struct Accepted { uint32_t a, b; int32_t g; std::string consensus; size_t reads, pairs, frags, modal; };
    std::vector<Accepted> accepted;
    for (const Group& gr : groups) {
        const Nom& n0 = noms[gr.lo];
        const End& EA = ends[n0.a];
        const End& EB = ends[n0.b];
        const size_t nReads = gr.hi - gr.lo;
        const size_t nPairs = gr.byPair.size();
        std::vector<size_t> frags;
        {
            std::set<std::pair<int, int32_t>> seen;
            for (size_t x : gr.byPair) {
                if (seen.insert({noms[x].strand, noms[x].key}).second) frags.push_back(x);
            }
        }
        std::map<int32_t, size_t> gCount;
        for (size_t x : frags) ++gCount[noms[x].g];
        int32_t mode = 0;
        size_t modal = 0;
        for (const auto& kv : gCount) {
            if (kv.second > modal) { modal = kv.second; mode = kv.first; }
        }
        auto refuse = [&](size_t& counter, const char* why) {
            ++counter;
            if (trace && nPairs >= minReads) {
                std::fprintf(stderr,
                             "[dropoutbridge-refused] reason=%s a=%u:%d b=%u:%d gap=%d reads=%zu "
                             "pairs=%zu frags=%zu modal=%zu cov_a=%.2f cov_b=%.2f len_a=%zu len_b=%zu\n",
                             why, EA.u, EA.e, EB.u, EB.e, static_cast<int>(mode), nReads, nPairs,
                             frags.size(), modal, EA.cov, EB.cov, EA.len, EB.len);
            }
        };
        if (nPairs < minReads) { refuse(st.refMinReads, "minreads"); continue; }
        if (frags.size() < minReads) { refuse(st.refDup, "dup"); continue; }
        if (modal < minReads || 3 * modal < 2 * frags.size()) { refuse(st.refModal, "modal"); continue; }
        std::string consensus;
        size_t rawVotes = 0;
        if (mode > 0) {
            std::vector<const Nom*> modalNoms;
            for (size_t x : frags) {
                if (noms[x].g == mode) { modalNoms.push_back(&noms[x]); rawVotes += noms[x].rawVotes; }
            }
            bool ok = true;
            consensus.resize(static_cast<size_t>(mode));
            for (size_t p = 0; ok && p < static_cast<size_t>(mode); ++p) {
                size_t cnt[4] = {0, 0, 0, 0};
                for (const Nom* m : modalNoms) {
                    const int c = baseCode(m->gap[p]);
                    if (c >= 0) ++cnt[c];
                }
                int best = 0;
                for (int c = 1; c < 4; ++c) if (cnt[c] > cnt[best]) best = c;
                if (3 * cnt[best] < 2 * modalNoms.size()) ok = false;
                consensus[p] = codeBase(best);
            }
            if (!ok) { refuse(st.refConsensus, "consensus"); continue; }
        }
        if (strong[n0.a] > 1 || strong[n0.b] > 1) { refuse(st.refCompeting, "competing"); continue; }
        const size_t o = mode < 0 ? static_cast<size_t>(-mode) : 0;
        if (mode < 0) {
            if (o >= EA.tail.size() || o >= EB.head.size() ||
                EA.tail.compare(EA.tail.size() - o, o, EB.head, 0, o) != 0) {
                refuse(st.refOverlap, "overlap");
                continue;
            }
        }
        {
            const std::string winA = EA.tail.substr(EA.tail.size() - static_cast<size_t>(EA.window));
            const size_t wb = static_cast<size_t>(EB.window);
            const std::string winB = o < wb ? EB.head.substr(o, wb - o) : std::string();
            if (lowComplexity(winA + consensus + winB)) { refuse(st.refLowComplexity, "lowcx"); continue; }
        }
        if (EA.u == EB.u) { refuse(st.refSelf, "self"); continue; }
        st.rawGapVotes += rawVotes;
        accepted.push_back({n0.a, n0.b, mode, std::move(consensus), nReads, nPairs, frags.size(), modal});
    }
    sampleRss();

    // ---- 5. insertion (closeGapsByOverlap::insertBridge convention) -----------------------
    const size_t K1 = static_cast<size_t>(K - 1);
    for (const Accepted& ac : accepted) {
        const End& EA = ends[ac.a];
        const End& EB = ends[ac.b];
        if (g.nodes[EA.u].deleted || g.nodes[EB.u].deleted) continue;
        if (!g.nodes[EA.u].ends[EA.e].empty() || !g.nodes[EB.u].ends[EB.e].empty()) continue;
        std::string bseq = EA.tail.substr(EA.tail.size() - K1);
        if (ac.g >= 0) {
            bseq += ac.consensus;
            bseq += EB.head.substr(0, K1);
        } else {
            const size_t o = static_cast<size_t>(-ac.g);
            if (o > K1) continue;
            bseq += EB.head.substr(o, K1 - o);
        }
        if (bseq.size() < static_cast<size_t>(K)) continue;   // must be a legal unitig
        if (trace) {
            const size_t fl = 30;
            std::string junction = EA.tail.substr(EA.tail.size() - fl);
            if (ac.g >= 0) junction += ac.consensus + EB.head.substr(0, fl);
            else junction += EB.head.substr(static_cast<size_t>(-ac.g), fl);
            std::fprintf(stderr,
                         "[dropoutbridge-join] a=%u:%d b=%u:%d gap=%d reads=%zu frags=%zu cov_a=%.2f "
                         "cov_b=%.2f len_a=%zu len_b=%zu pairs=%zu modal=%zu bridge_len=%zu junction=%s\n",
                         EA.u, EA.e, EB.u, EB.e, static_cast<int>(ac.g), ac.reads, ac.frags, EA.cov,
                         EB.cov, EA.len, EB.len, ac.pairs, ac.modal, bseq.size(), junction.c_str());
        }
        Unitig nb;
        nb.seq = std::move(bseq);
        nb.coverage = std::min(EA.cov, EB.cov);
        const uint32_t b = static_cast<uint32_t>(g.nodes.size());
        g.nodes.push_back(std::move(nb));
        linkEnds(g, EA.u, EA.e, b, 0);
        linkEnds(g, b, 1, EB.u, EB.e);
        ++st.bridged;
        if (ac.g < 0) ++st.overlapBridged;
        else st.bpAdded += static_cast<size_t>(ac.g);
    }
    if (st.bridged > 0) {
        st.merged = g.compact();
        st.error = g.validate();
        if (!st.error.empty()) st.error = "dropout bridge left an invalid graph: " + st.error;
    }
    sampleRss();
    st.rssDeltaMb = static_cast<double>(rssPeak - rss0) / (1024.0 * 1024.0);
    st.seconds = timer.elapsed();
    return st;
}

std::string formatDropoutBridgeStats(const DropoutBridgeStats& s) {
    char buf[1024];
    std::snprintf(buf, sizeof(buf),
                  "[dropoutbridge] K=%d med=%.2f deadends=%zu eligible=%zu anchors=%zu anchors_dropped=%zu "
                  "reads_hit=%zu pairs_nominated=%zu ref_minreads=%zu ref_dup=%zu ref_modal=%zu "
                  "ref_consensus=%zu ref_competing=%zu ref_overlap=%zu ref_lowcx=%zu ref_self=%zu "
                  "bridged=%zu overlap_bridged=%zu bp_added=%zu seconds=%.2f min_reads=%u nominations=%zu "
                  "ends_noanchor=%zu merged=%zu raw_gap_votes=%zu masked_hits=%zu rss_delta_mb=%.1f enabled=%d",
                  s.k, s.medianCoverage, s.deadEnds, s.eligible, s.anchors, s.anchorsDropped, s.readsHit,
                  s.pairsNominated, s.refMinReads, s.refDup, s.refModal, s.refConsensus, s.refCompeting,
                  s.refOverlap, s.refLowComplexity, s.refSelf, s.bridged, s.overlapBridged, s.bpAdded,
                  s.seconds, s.minReads, s.nominations, s.endsNoAnchor, s.merged, s.rawGapVotes,
                  s.maskedHits, s.rssDeltaMb, s.enabled ? 1 : 0);
    return buf;
}

}  // namespace ts
