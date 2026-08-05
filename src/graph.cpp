#include "graph.h"

#include "util.h"

#include <algorithm>
#include <thread>
#include <cstdio>
#include <cstdlib>
#include <unordered_map>
#include <unordered_set>

namespace ts {

std::string reverseComplement(const std::string& s) {
    std::string out(s.size(), 'N');
    for (size_t i = 0; i < s.size(); ++i) {
        char c = s[s.size() - 1 - i];
        switch (c) {
            case 'A': out[i] = 'T'; break;
            case 'C': out[i] = 'G'; break;
            case 'G': out[i] = 'C'; break;
            case 'T': out[i] = 'A'; break;
            default:  out[i] = 'N'; break;
        }
    }
    return out;
}

double sequenceIdentity(const std::string& a, const std::string& b, int maxBand) {
    const int n = static_cast<int>(a.size());
    const int m = static_cast<int>(b.size());
    if (n == 0 && m == 0) return 1.0;
    if (n == 0 || m == 0) return 0.0;
    if (std::abs(n - m) > maxBand) return 0.0;

    const int band = std::max(maxBand, std::abs(n - m)) + 1;
    // Banded Levenshtein: only |i-j| <= band is ever reachable within the
    // distance budget, so two rows of width 2*band+1 are enough.
    const int width = 2 * band + 1;
    std::vector<int> prev(static_cast<size_t>(width), INT32_MAX / 2);
    std::vector<int> cur(static_cast<size_t>(width), INT32_MAX / 2);

    auto at = [band](int j, int i) { return j - i + band; };

    for (int j = 0; j <= std::min(m, band); ++j) prev[static_cast<size_t>(at(j, 0))] = j;

    for (int i = 1; i <= n; ++i) {
        std::fill(cur.begin(), cur.end(), INT32_MAX / 2);
        const int lo = std::max(0, i - band);
        const int hi = std::min(m, i + band);
        for (int j = lo; j <= hi; ++j) {
            int best = INT32_MAX / 2;
            if (j > 0) {
                int idx = at(j - 1, i);
                if (idx >= 0 && idx < width) best = std::min(best, cur[static_cast<size_t>(idx)] + 1);
            }
            int idxUp = at(j, i - 1);
            if (idxUp >= 0 && idxUp < width) best = std::min(best, prev[static_cast<size_t>(idxUp)] + 1);
            if (j > 0) {
                int idxDiag = at(j - 1, i - 1);
                if (idxDiag >= 0 && idxDiag < width) {
                    int cost = (a[static_cast<size_t>(i - 1)] == b[static_cast<size_t>(j - 1)]) ? 0 : 1;
                    best = std::min(best, prev[static_cast<size_t>(idxDiag)] + cost);
                }
            }
            if (j == 0) best = std::min(best, i);
            cur[static_cast<size_t>(at(j, i))] = best;
        }
        prev.swap(cur);
    }
    int idx = at(m, n);
    if (idx < 0 || idx >= width) return 0.0;
    int dist = prev[static_cast<size_t>(idx)];
    int denom = std::max(n, m);
    if (denom == 0) return 1.0;
    double id = 1.0 - static_cast<double>(dist) / static_cast<double>(denom);
    return id < 0 ? 0 : id;
}

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------

UnitigGraph UnitigGraph::build(const KmerTable& solid, int k, int threads) {
    UnitigGraph g;
    g.k_ = k;
    if (solid.size() == 0) return g;

    auto isSolid = [&](Kmer x) { return solid.contains(canonical(x, k)); };
    auto covOf = [&](Kmer x) { return static_cast<double>(solid.get(canonical(x, k))); };

    auto outNbrs = [&](Kmer x, Kmer* dst) {
        int n = 0;
        for (int b = 0; b < 4; ++b) {
            Kmer y = pushBack(x, b, k);
            if (isSolid(y)) dst[n++] = y;
        }
        return n;
    };
    auto inNbrs = [&](Kmer x, Kmer* dst) {
        int n = 0;
        for (int b = 0; b < 4; ++b) {
            Kmer y = pushFront(x, b, k);
            if (isSolid(y)) dst[n++] = y;
        }
        return n;
    };

    std::vector<Kmer> keys;
    keys.reserve(solid.size());
    solid.forEach([&](Kmer key, uint32_t) { keys.push_back(key); });
    // forEach yields open-addressing slot order, and which of two colliding keys
    // takes the earlier slot depends on insertion order, which the sharded
    // counter does not fix. Every later decision -- node ids, which orientation
    // of a unitig is stored, the sweep order of the greedy simplifications --
    // follows from this order, so it is sorted rather than left to chance.
    std::sort(keys.begin(), keys.end());

    KmerTable visited;
    size_t placedKmers = 0;

    Kmer nb[4], nb2[4];

    // Walks right from `start`, consuming k-mers into a new unitig.
    auto extend = [&](Kmer start) {
        std::string seq = kmerToString(start, k);
        Kmer cur = start;
        visited.put(canonical(cur, k), 1);
        double covSum = covOf(cur);
        size_t nk = 1;

        while (true) {
            if (outNbrs(cur, nb) != 1) break;
            Kmer nxt = nb[0];
            if (inNbrs(nxt, nb2) != 1) break;
            if (visited.contains(canonical(nxt, k))) break;   // closed a cycle
            seq += codeBase(lowBase(nxt));
            cur = nxt;
            visited.put(canonical(cur, k), 1);
            covSum += covOf(cur);
            ++nk;
        }
        Unitig u;
        u.seq = std::move(seq);
        u.coverage = nk ? covSum / static_cast<double>(nk) : 0;
        g.nodes.push_back(std::move(u));
    };

    // Deciding which k-mers start a unitig is a pure function of `solid`: it
    // touches no shared state, and at eight table probes per k-mer per
    // orientation it is the bulk of construction. Doing it up front across
    // threads leaves the walk itself sequential, so unitigs come out in exactly
    // the same order and the assembly stays byte-identical.
    util::Timer phaseTimer;
    const bool phaseDebug = std::getenv("TESSERA_GRAPH_PHASES") != nullptr;
    const double tKeys = phaseTimer.elapsed();

    std::vector<uint8_t> startFlag(keys.size() * 2, 0);
    {
        int nt = threads > 0 ? threads : 1;
        if (nt > static_cast<int>(keys.size())) nt = static_cast<int>(keys.size());
        if (nt < 1) nt = 1;
        std::vector<std::thread> pool;
        pool.reserve(static_cast<size_t>(nt));
        for (int t = 0; t < nt; ++t) {
            pool.emplace_back([&, t]() {
                // Counts solid successors; only the count is needed here.
                auto outDegree = [&](Kmer x) {
                    int n = 0;
                    for (int b = 0; b < 4; ++b) {
                        if (solid.contains(canonical(pushBack(x, b, k), k))) ++n;
                    }
                    return n;
                };
                for (size_t i = static_cast<size_t>(t); i < keys.size();
                     i += static_cast<size_t>(nt)) {
                    // keys are canonical, so one reverseComplement per key
                    // serves both orientations.
                    const Kmer fwd = keys[i];
                    const Kmer rev = reverseComplement(fwd, k);
                    for (int orient = 0; orient < 2; ++orient) {
                        const Kmer x = orient == 0 ? fwd : rev;
                        // A k-mer starts a unitig unless it is entered from a
                        // unique predecessor that itself continues uniquely
                        // into it.
                        Kmer pred = 0;
                        int nIn = 0;
                        for (int b = 0; b < 4; ++b) {
                            const Kmer y = pushFront(x, b, k);
                            if (solid.contains(canonical(y, k))) { pred = y; ++nIn; }
                        }
                        const bool start = (nIn != 1) || (outDegree(pred) != 1);
                        startFlag[i * 2 + static_cast<size_t>(orient)] = start ? 1 : 0;
                    }
                }
            });
        }
        for (auto& th : pool) th.join();
    }

    const double tFlags = phaseTimer.elapsed();

    // The walk is the bulk of construction -- 78-84% of build time once the
    // start classification above was threaded -- and it parallelises exactly,
    // because the walks are disjoint by construction. A step is only taken from
    // `cur` to `nxt` when `cur` has one successor and `nxt` one predecessor,
    // which makes `nxt` not a start; so no other walk begins there, and any walk
    // reaching it must have come through `cur` and is therefore this one. The
    // shared `visited` table was never arbitrating between walks. It was doing
    // two other jobs, and each is handled separately below.
    //
    // Job one: stopping a walk that closes on itself. That collision can only
    // be with the current walk, so a per-walk set does it.
    // Job two: skipping a unitig's second discovery, from its other end. Both
    // orientations are walked here and one is dropped afterwards -- the second
    // walk retraces the first's nodes backwards and so yields exactly its
    // reverse complement.
    struct StartRef { uint32_t key; uint8_t orient; };
    std::vector<StartRef> starts;
    starts.reserve(keys.size());
    for (size_t i = 0; i < keys.size(); ++i) {
        for (int orient = 0; orient < 2; ++orient) {
            if (startFlag[i * 2 + static_cast<size_t>(orient)]) {
                starts.push_back({static_cast<uint32_t>(i), static_cast<uint8_t>(orient)});
            }
        }
    }

    std::vector<Unitig> walked(starts.size());
    std::vector<uint32_t> walkedKmers(starts.size(), 0);
    {
        int nt = threads > 0 ? threads : 1;
        if (nt > static_cast<int>(starts.size())) nt = static_cast<int>(starts.size());
        if (nt < 1) nt = 1;
        std::vector<std::thread> pool;
        pool.reserve(static_cast<size_t>(nt));
        for (int t = 0; t < nt; ++t) {
            pool.emplace_back([&, t]() {
                std::unordered_set<Kmer, KmerHasher> seen;
                Kmer lnb[4], lnb2[4];
                auto lOut = [&](Kmer x, Kmer* dst) {
                    int n = 0;
                    for (int b = 0; b < 4; ++b) {
                        Kmer y = pushBack(x, b, k);
                        if (solid.contains(canonical(y, k))) dst[n++] = y;
                    }
                    return n;
                };
                auto lIn = [&](Kmer x, Kmer* dst) {
                    int n = 0;
                    for (int b = 0; b < 4; ++b) {
                        Kmer y = pushFront(x, b, k);
                        if (solid.contains(canonical(y, k))) dst[n++] = y;
                    }
                    return n;
                };
                for (size_t s = static_cast<size_t>(t); s < starts.size();
                     s += static_cast<size_t>(nt)) {
                    Kmer cur = starts[s].orient == 0
                                   ? keys[starts[s].key]
                                   : reverseComplement(keys[starts[s].key], k);
                    seen.clear();
                    std::string seq = kmerToString(cur, k);
                    seen.insert(canonical(cur, k));
                    double covSum = static_cast<double>(solid.get(canonical(cur, k)));
                    uint32_t nk = 1;
                    while (true) {
                        if (lOut(cur, lnb) != 1) break;
                        const Kmer nxt = lnb[0];
                        if (lIn(nxt, lnb2) != 1) break;
                        const Kmer cn = canonical(nxt, k);
                        if (seen.count(cn)) break;   // closed a cycle
                        seq += codeBase(lowBase(nxt));
                        cur = nxt;
                        seen.insert(cn);
                        covSum += static_cast<double>(solid.get(cn));
                        ++nk;
                    }
                    walked[s].seq = std::move(seq);
                    walked[s].coverage = covSum / static_cast<double>(nk);
                    walkedKmers[s] = nk;
                }
            });
        }
        for (auto& th : pool) th.join();
    }

    // Keep the first orientation of each unitig, in the order the sequential
    // walk would have reached it.
    {
        std::unordered_set<std::string> emitted;
        emitted.reserve(starts.size());
        std::string rc;
        for (size_t s = 0; s < starts.size(); ++s) {
            const std::string& seq = walked[s].seq;
            rc.assign(seq.rbegin(), seq.rend());
            for (char& c : rc) {
                switch (c) {
                    case 'A': c = 'T'; break;
                    case 'C': c = 'G'; break;
                    case 'G': c = 'C'; break;
                    case 'T': c = 'A'; break;
                    default: break;
                }
            }
            if (emitted.count(seq < rc ? seq : rc)) continue;
            emitted.insert(seq < rc ? seq : rc);
            placedKmers += walkedKmers[s];
            g.nodes.push_back(std::move(walked[s]));
        }
    }

    // Whatever remains lies on a cycle where every k-mer has in-degree and
    // out-degree 1, so no start exists; break each such cycle arbitrarily.
    // Emitted unitigs partition the solid set, so anything left over shows up
    // as a shortfall in the k-mer count -- and when there is none, which is the
    // usual case, the marking pass this needs is skipped entirely.
    if (placedKmers < solid.size()) {
        for (const Unitig& u : g.nodes) {
            Kmer fwd = 0;
            int valid = 0;
            for (size_t p = 0; p < u.seq.size(); ++p) {
                const int b = baseCode(u.seq[p]);
                if (b < 0) { valid = 0; continue; }
                fwd = pushBack(fwd, b, k);
                if (++valid >= k) visited.put(canonical(fwd, k), 1);
            }
        }
        for (Kmer key : keys) {
            if (visited.contains(key)) continue;
            extend(key);
        }
    }

    const double tWalk = phaseTimer.elapsed();

    // ---- link the unitigs ----
    // (timing for this phase is reported at the end of build)
    std::unordered_map<Kmer, std::pair<uint32_t, uint8_t>, KmerHasher> entry;
    entry.reserve(g.nodes.size() * 2);
    for (uint32_t i = 0; i < g.nodes.size(); ++i) {
        const std::string& s = g.nodes[i].seq;
        bool ok = false;
        Kmer first = stringToKmer(s, k, ok);
        if (ok) entry.emplace(first, std::make_pair(i, static_cast<uint8_t>(0)));
        Kmer last = stringToKmer(s.substr(s.size() - static_cast<size_t>(k)), k, ok);
        if (ok) entry.emplace(reverseComplement(last, k), std::make_pair(i, static_cast<uint8_t>(1)));
    }

    for (uint32_t i = 0; i < g.nodes.size(); ++i) {
        const std::string& s = g.nodes[i].seq;
        bool ok = false;
        Kmer last = stringToKmer(s.substr(s.size() - static_cast<size_t>(k)), k, ok);
        if (ok) {
            for (int b = 0; b < 4; ++b) {
                Kmer y = pushBack(last, b, k);
                if (!isSolid(y)) continue;
                auto it = entry.find(y);
                if (it != entry.end()) g.addLink(i, 1, it->second.first, it->second.second);
            }
        }
        Kmer firstRc = 0;
        Kmer first = stringToKmer(s, k, ok);
        if (ok) {
            firstRc = reverseComplement(first, k);
            for (int b = 0; b < 4; ++b) {
                Kmer y = pushBack(firstRc, b, k);
                if (!isSolid(y)) continue;
                auto it = entry.find(y);
                if (it != entry.end()) g.addLink(i, 0, it->second.first, it->second.second);
            }
        }
    }
    if (phaseDebug) {
        const double tLink = phaseTimer.elapsed();
        std::fprintf(stderr,
                     "      [graph k=%d] keys %.2fs  startFlags %.2fs  walk %.2fs  link %.2fs\n",
                     k, tKeys, tFlags - tKeys, tWalk - tFlags, tLink - tWalk);
    }
    return g;
}

// ---------------------------------------------------------------------------
// Link maintenance
// ---------------------------------------------------------------------------

void UnitigGraph::addLink(uint32_t u, int ue, uint32_t v, int ve) {
    Link a{v, static_cast<uint8_t>(ve)};
    auto& lu = nodes[u].ends[ue];
    if (std::find(lu.begin(), lu.end(), a) == lu.end()) lu.push_back(a);

    Link b{u, static_cast<uint8_t>(ue)};
    auto& lv = nodes[v].ends[ve];
    if (std::find(lv.begin(), lv.end(), b) == lv.end()) lv.push_back(b);
}

void UnitigGraph::unlink(uint32_t u, int ue, uint32_t v, int ve) {
    auto rm = [](std::vector<Link>& vec, Link target) {
        vec.erase(std::remove(vec.begin(), vec.end(), target), vec.end());
    };
    rm(nodes[u].ends[ue], Link{v, static_cast<uint8_t>(ve)});
    rm(nodes[v].ends[ve], Link{u, static_cast<uint8_t>(ue)});
}

void UnitigGraph::deleteNode(uint32_t u) {
    if (nodes[u].deleted) return;
    for (int e = 0; e < 2; ++e) {
        // Copy: unlink mutates the vector being iterated.
        std::vector<Link> links = nodes[u].ends[e];
        for (const Link& l : links) unlink(u, e, l.to, l.toEnd);
    }
    nodes[u].ends[0].clear();
    nodes[u].ends[1].clear();
    nodes[u].deleted = true;
}

void UnitigGraph::removeDeleted() {
    std::vector<uint32_t> remap(nodes.size(), UINT32_MAX);
    uint32_t next = 0;
    for (uint32_t i = 0; i < nodes.size(); ++i) {
        if (!nodes[i].deleted) remap[i] = next++;
    }
    std::vector<Unitig> kept;
    kept.reserve(next);
    for (uint32_t i = 0; i < nodes.size(); ++i) {
        if (nodes[i].deleted) continue;
        Unitig u = std::move(nodes[i]);
        for (int e = 0; e < 2; ++e) {
            std::vector<Link> fixed;
            fixed.reserve(u.ends[e].size());
            for (const Link& l : u.ends[e]) {
                if (remap[l.to] != UINT32_MAX) fixed.push_back(Link{remap[l.to], l.toEnd});
            }
            u.ends[e].swap(fixed);
        }
        kept.push_back(std::move(u));
    }
    nodes.swap(kept);
}

std::string UnitigGraph::oriented(uint32_t u, int orient) const {
    return orient == 0 ? nodes[u].seq : reverseComplement(nodes[u].seq);
}

// ---------------------------------------------------------------------------
// Compaction
// ---------------------------------------------------------------------------

bool UnitigGraph::mergeInto(uint32_t u, int ue) {
    Unitig& U = nodes[u];
    if (U.deleted || U.ends[ue].size() != 1) return false;
    const Link l = U.ends[ue][0];
    const uint32_t v = l.to;
    const int ve = l.toEnd;
    if (v == u) return false;                       // self loop
    if (nodes[v].deleted) return false;
    if (nodes[v].ends[ve].size() != 1) return false;   // v is a branch point

    Unitig& V = nodes[v];
    const size_t ov = static_cast<size_t>(k_ - 1);
    // Orient v so its attaching end faces u, then splice off the shared (k-1)
    // overlap.
    const std::string vs = (ve == 0) ? V.seq : reverseComplement(V.seq);

    const size_t nkU = U.seq.size() >= static_cast<size_t>(k_) ? U.seq.size() - static_cast<size_t>(k_) + 1 : 1;
    const size_t nkV = V.seq.size() >= static_cast<size_t>(k_) ? V.seq.size() - static_cast<size_t>(k_) + 1 : 1;
    const double newCov = (U.coverage * static_cast<double>(nkU) + V.coverage * static_cast<double>(nkV)) /
                          static_cast<double>(nkU + nkV);

    // The far end of v becomes the new far end of the merged unitig.
    const int vFar = 1 - ve;
    std::vector<Link> vFarLinks = V.ends[vFar];

    if (ue == 1) {
        U.seq += vs.substr(ov);
    } else {
        // Entering u's 5' end: orient v so its 3' end meets u, then prepend.
        const std::string head = (ve == 1) ? V.seq : reverseComplement(V.seq);
        U.seq = head.substr(0, head.size() - ov) + U.seq;
    }
    U.coverage = newCov;

    // Detach v entirely, then reattach its far-end links to u's merged end.
    for (const Link& fl : vFarLinks) unlink(v, vFar, fl.to, fl.toEnd);
    unlink(u, ue, v, ve);
    V.deleted = true;
    V.ends[0].clear();
    V.ends[1].clear();

    for (const Link& fl : vFarLinks) {
        if (fl.to == v) continue;   // v had a self-loop; it dies with v
        addLink(u, ue, fl.to, fl.toEnd);
    }
    return true;
}

size_t UnitigGraph::compact() {
    size_t merged = 0;
    bool progress = true;
    while (progress) {
        progress = false;
        for (uint32_t u = 0; u < nodes.size(); ++u) {
            if (nodes[u].deleted) continue;
            for (int e = 0; e < 2; ++e) {
                while (mergeInto(u, e)) { ++merged; progress = true; }
            }
        }
    }
    if (merged) removeDeleted();
    return merged;
}

// ---------------------------------------------------------------------------
// Simplification
// ---------------------------------------------------------------------------

size_t UnitigGraph::removeTips(size_t maxLen, double covRatio) {
    size_t removed = 0;
    for (uint32_t u = 0; u < nodes.size(); ++u) {
        Unitig& U = nodes[u];
        if (U.deleted || U.seq.size() > maxLen) continue;

        int dead;
        if (U.ends[0].empty() && !U.ends[1].empty()) dead = 0;
        else if (U.ends[1].empty() && !U.ends[0].empty()) dead = 1;
        else continue;

        const int att = 1 - dead;
        if (U.ends[att].size() != 1) continue;
        const Link l = U.ends[att][0];
        if (l.to == u) continue;

        // The tip only matters if the junction offers a better alternative.
        double bestAlt = 0;
        size_t alts = 0;
        for (const Link& sib : nodes[l.to].ends[l.toEnd]) {
            if (sib.to == u) continue;
            ++alts;
            bestAlt = std::max(bestAlt, nodes[sib.to].coverage);
        }
        if (alts == 0) continue;
        if (U.coverage < covRatio * bestAlt) {
            deleteNode(u);
            ++removed;
        }
    }
    return removed;
}

namespace {

struct SimplePath {
    uint32_t endNode = 0;
    int endEnd = 0;
    std::string seq;
    double minCov = 0;
    std::vector<uint32_t> nodes;
    bool valid = false;
};

}  // namespace

size_t UnitigGraph::popBubbles(size_t maxLen, double minIdentity, double maxLoserCoverage) {
    size_t popped = 0;
    const size_t ov = static_cast<size_t>(k_ - 1);

    // Follows a link through unbranching nodes, collecting the traversed
    // sequence until a branch point or the length budget is reached.
    auto walk = [&](const Link& start, size_t budget) {
        SimplePath p;
        uint32_t cur = start.to;
        int entry = start.toEnd;
        for (int step = 0; step < 12; ++step) {
            if (nodes[cur].deleted) return p;
            const int exit = 1 - entry;
            std::string s = (entry == 0) ? nodes[cur].seq : reverseComplement(nodes[cur].seq);
            if (p.seq.empty()) p.seq = s;
            else p.seq += s.substr(ov);
            p.nodes.push_back(cur);
            p.minCov = p.nodes.size() == 1 ? nodes[cur].coverage
                                           : std::min(p.minCov, nodes[cur].coverage);
            if (p.seq.size() > budget) return p;

            if (nodes[cur].ends[entry].size() != 1) return p;   // re-entrant branch
            if (nodes[cur].ends[exit].size() != 1) return p;    // diverges again
            const Link nxt = nodes[cur].ends[exit][0];
            if (nodes[nxt.to].ends[nxt.toEnd].size() > 1) {
                // The node the path flows *into* is the convergence point. Two
                // sides of a bubble end on different last nodes but share this
                // one, so it is what the two walks must agree on.
                p.endNode = nxt.to;
                p.endEnd = nxt.toEnd;
                p.valid = true;
                return p;
            }
            cur = nxt.to;
            entry = nxt.toEnd;
        }
        return p;
    };

    for (uint32_t u = 0; u < nodes.size(); ++u) {
        if (nodes[u].deleted) continue;
        for (int e = 0; e < 2; ++e) {
            if (nodes[u].ends[e].size() < 2) continue;
            std::vector<Link> links = nodes[u].ends[e];
            bool changed = false;
            for (size_t i = 0; i < links.size() && !changed; ++i) {
                for (size_t j = i + 1; j < links.size() && !changed; ++j) {
                    if (nodes[links[i].to].deleted || nodes[links[j].to].deleted) continue;
                    SimplePath a = walk(links[i], maxLen);
                    SimplePath b = walk(links[j], maxLen);
                    if (!a.valid || !b.valid) continue;
                    if (a.endNode != b.endNode || a.endEnd != b.endEnd) continue;
                    if (a.nodes.empty() || b.nodes.empty()) continue;
                    // Same node on both sides means this is not a real bubble.
                    if (a.nodes.front() == b.nodes.front()) continue;

                    const int band = static_cast<int>(std::max<size_t>(8, maxLen / 10));
                    if (sequenceIdentity(a.seq, b.seq, band) < minIdentity) continue;

                    const SimplePath& loser = (a.minCov <= b.minCov) ? a : b;
                    const SimplePath& winner = (a.minCov <= b.minCov) ? b : a;
                    if (maxLoserCoverage > 0 && loser.minCov > maxLoserCoverage) continue;
                    bool shared = false;
                    for (uint32_t ln : loser.nodes) {
                        if (std::find(winner.nodes.begin(), winner.nodes.end(), ln) != winner.nodes.end()) {
                            shared = true;
                            break;
                        }
                    }
                    if (shared) continue;
                    for (uint32_t ln : loser.nodes) {
                        if (ln != u && ln != a.endNode) deleteNode(ln);
                    }
                    ++popped;
                    changed = true;
                }
            }
        }
    }
    return popped;
}


int UnitigGraph::deadEnds(uint32_t u) const {
    if (nodes[u].deleted) return 0;
    return (nodes[u].ends[0].empty() ? 1 : 0) + (nodes[u].ends[1].empty() ? 1 : 0);
}

int UnitigGraph::deadEndChangeIfDeleted(uint32_t u) const {
    if (nodes[u].deleted) return 0;
    const int before = deadEnds(u);
    int after = 0;
    for (int e = 0; e < 2; ++e) {
        for (const Link& l : nodes[u].ends[e]) {
            if (nodes[l.to].deleted) continue;
            // The neighbour keeps a connection only if something other than
            // `u` attaches to that end of it.
            size_t live = 0;
            for (const Link& back : nodes[l.to].ends[l.toEnd]) {
                if (!nodes[back.to].deleted && back.to != u) ++live;
            }
            if (live == 0) ++after;
        }
    }
    return after - before;
}

size_t UnitigGraph::totalDeadEnds() const {
    size_t n = 0;
    for (uint32_t u = 0; u < nodes.size(); ++u) {
        if (!nodes[u].deleted) n += static_cast<size_t>(deadEnds(u));
    }
    return n;
}

std::vector<std::vector<uint32_t>> UnitigGraph::components() const {
    std::vector<std::vector<uint32_t>> out;
    std::vector<char> seen(nodes.size(), 0);
    std::vector<uint32_t> stack;
    for (uint32_t start = 0; start < nodes.size(); ++start) {
        if (nodes[start].deleted || seen[start]) continue;
        out.emplace_back();
        stack.push_back(start);
        seen[start] = 1;
        while (!stack.empty()) {
            const uint32_t u = stack.back();
            stack.pop_back();
            out.back().push_back(u);
            for (int e = 0; e < 2; ++e) {
                for (const Link& l : nodes[u].ends[e]) {
                    if (nodes[l.to].deleted || seen[l.to]) continue;
                    seen[l.to] = 1;
                    stack.push_back(l.to);
                }
            }
        }
    }
    return out;
}

double UnitigGraph::weightedMedianCoverage(const std::vector<uint32_t>& ids, size_t topN) const {
    std::vector<uint32_t> live;
    live.reserve(ids.size());
    for (uint32_t u : ids) {
        if (!nodes[u].deleted) live.push_back(u);
    }
    if (live.empty()) return 0.0;
    std::sort(live.begin(), live.end(), [&](uint32_t a, uint32_t b) {
        return nodes[a].seq.size() > nodes[b].seq.size();
    });
    if (live.size() > topN) live.resize(topN);
    std::sort(live.begin(), live.end(), [&](uint32_t a, uint32_t b) {
        return nodes[a].coverage < nodes[b].coverage;
    });
    size_t total = 0;
    for (uint32_t u : live) total += nodes[u].seq.size();
    if (total == 0) return 0.0;
    size_t acc = 0;
    for (uint32_t u : live) {
        acc += nodes[u].seq.size();
        if (acc * 2 >= total) return nodes[u].coverage;
    }
    return nodes[live.back()].coverage;
}

size_t UnitigGraph::filterByReadDepth(double fraction) {
    if (fraction <= 0) return 0;

    std::vector<uint32_t> all;
    all.reserve(nodes.size());
    for (uint32_t u = 0; u < nodes.size(); ++u) {
        if (!nodes[u].deleted) all.push_back(u);
    }
    if (all.size() < 2) return 0;

    const double globalCut = fraction * weightedMedianCoverage(all);
    if (globalCut <= 0) return 0;

    size_t removed = 0;
    for (const std::vector<uint32_t>& comp : components()) {
        const double compMedian = weightedMedianCoverage(comp);
        const double compCut = fraction * compMedian;
        // A whole component sitting under the global bar is contamination or
        // an error cloud rather than a low-copy replicon, and may go wholesale.
        const bool wholeComponentLow = compMedian < globalCut;

        for (uint32_t u : comp) {
            if (nodes[u].deleted) continue;
            const double cov = nodes[u].coverage;
            if (cov >= globalCut && cov >= compCut) continue;
            // Keep anything whose removal would strand a neighbour, unless the
            // component was already condemned.
            if (!wholeComponentLow && deadEnds(u) == 0 && deadEndChangeIfDeleted(u) > 0) continue;
            deleteNode(u);
            ++removed;
        }
    }
    return removed;
}


size_t UnitigGraph::removeLocallyWeak(double fraction, double absoluteFloor) {
    if (fraction <= 0) return 0;
    // Decided against the graph as it stands, then applied, so that removing
    // one weak unitig cannot make its equally weak neighbour look strong.
    std::vector<char> doomed(nodes.size(), 0);
    for (uint32_t u = 0; u < nodes.size(); ++u) {
        const Unitig& U = nodes[u];
        if (U.deleted) continue;
        std::vector<double> nbr;
        for (int e = 0; e < 2; ++e) {
            for (const Link& l : U.ends[e]) {
                if (!nodes[l.to].deleted) nbr.push_back(nodes[l.to].coverage);
            }
        }
        // An isolated unitig has no local scale to be judged against.
        if (nbr.size() < 2) continue;
        std::sort(nbr.begin(), nbr.end());
        const double localMedian = nbr[nbr.size() / 2];
        if (localMedian <= 0) continue;
        // Relative weakness alone is not enough. A plasmid segment sitting at
        // thirty-fold next to a repeat at three hundred is a tenth of its
        // neighbours and still perfectly real -- deleting it cost 12.8 kb of a
        // 234 kb plasmid on one panel isolate. A sequencing error is weak in
        // both senses at once: a small fraction of its neighbours *and* far
        // below the depth at which the genome itself is assembled.
        if (U.coverage >= absoluteFloor) continue;
        if (U.coverage < fraction * localMedian) doomed[u] = 1;
    }
    size_t removed = 0;
    for (uint32_t u = 0; u < nodes.size(); ++u) {
        if (!doomed[u] || nodes[u].deleted) continue;
        deleteNode(u);
        ++removed;
    }
    return removed;
}

size_t UnitigGraph::removeErroneousConnections(double covThreshold, size_t maxLen) {
    size_t removed = 0;
    for (uint32_t u = 0; u < nodes.size(); ++u) {
        Unitig& U = nodes[u];
        if (U.deleted) continue;
        if (U.seq.size() > maxLen) continue;
        if (U.ends[0].empty() || U.ends[1].empty()) continue;   // not a connector
        if (U.coverage >= covThreshold) continue;

        // The previous guard here read
        //     if (U.ends[0].size() + U.ends[1].size() < 2) continue;
        // and could never fire: the two lines above already established that BOTH ends
        // carry at least one link, so the sum is always >= 2. It documented an intent
        // ("only cut when both junctions survive") that it did not implement, which is
        // why raising the length ceiling previously cost genome fraction -- long, real
        // connectors were being cut with no alternative path to carry the sequence.
        //
        // SPAdes gates the same operation on outdeg(start) > 1 && indeg(end) > 1: both
        // junctions must retain another edge once this one is gone. Implement that.
        bool alternatives = true;
        for (int e = 0; e < 2 && alternatives; ++e) {
            const Link& l = U.ends[e][0];
            size_t live = 0;
            for (const Link& sib : nodes[l.to].ends[l.toEnd]) {
                if (!nodes[sib.to].deleted) ++live;
            }
            if (live < 2) alternatives = false;   // only this edge -- cutting orphans it
        }
        if (!alternatives) continue;

        deleteNode(u);
        ++removed;
    }
    return removed;
}

size_t UnitigGraph::removeIsolated(double covThreshold, size_t maxLen) {
    size_t removed = 0;
    for (uint32_t u = 0; u < nodes.size(); ++u) {
        Unitig& U = nodes[u];
        if (U.deleted) continue;
        if (!U.ends[0].empty() || !U.ends[1].empty()) continue;
        if (U.seq.size() > maxLen) continue;
        if (U.coverage >= covThreshold) continue;
        deleteNode(u);
        ++removed;
    }
    return removed;
}


namespace {

// An overlap made mostly of one base carries no information: it matches
// anywhere, so accepting it would wire the graph together arbitrarily.
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
    // A longer overlap is allowed to be more skewed, since it is that much
    // less likely to occur by chance.
    const double limit = total >= 40 ? 0.90 : (total >= 20 ? 0.82 : 0.75);
    return static_cast<double>(best) > limit * static_cast<double>(total);
}

}  // namespace

size_t UnitigGraph::joinDeadEnds(size_t minOverlap) {
    if (minOverlap < 8) minOverlap = 8;
    const size_t maxOverlap = static_cast<size_t>(k_ - 1);
    if (maxOverlap < minOverlap) return 0;
    // Only an overlap of exactly k-1 can be taken. A link carries no length of
    // its own: mergeInto splices k-1 unconditionally and the GFA writes every
    // link as (k-1)M, so a join accepted at a shorter overlap would delete the
    // difference in real bases at the junction. Allowing shorter overlaps back
    // means giving Link its own overlap field first.
    minOverlap = maxOverlap;

    // Collect every loose end as the oriented sequence leaving it, so a tail
    // and a head can be compared directly.
    struct End { uint32_t u; int end; };
    std::vector<End> ends;
    for (uint32_t u = 0; u < nodes.size(); ++u) {
        if (nodes[u].deleted) continue;
        if (nodes[u].seq.size() < maxOverlap) continue;
        for (int e = 0; e < 2; ++e) {
            if (nodes[u].ends[e].empty()) ends.push_back({u, e});
        }
    }
    if (ends.size() < 2) return 0;

    // Index candidate heads by their first `minOverlap` bases, so each tail
    // only compares against ends that could possibly match.
    std::unordered_map<std::string, std::vector<size_t>> byPrefix;
    std::vector<std::string> headSeq(ends.size());
    for (size_t i = 0; i < ends.size(); ++i) {
        const Unitig& U = nodes[ends[i].u];
        // Leaving end 0 means the unitig is read reverse-complemented, so the
        // sequence that *arrives* at this end is the reverse complement.
        headSeq[i] = ends[i].end == 0 ? U.seq : reverseComplement(U.seq);
        byPrefix[headSeq[i].substr(0, minOverlap)].push_back(i);
    }

    size_t added = 0;
    std::vector<char> used(ends.size(), 0);
    for (size_t i = 0; i < ends.size(); ++i) {
        if (used[i]) continue;
        const Unitig& U = nodes[ends[i].u];
        if (U.deleted) continue;
        // The sequence running *out* of this end.
        const std::string tail = ends[i].end == 1 ? U.seq : reverseComplement(U.seq);

        bool joined = false;
        for (size_t ov = maxOverlap; ov >= minOverlap && !joined; --ov) {
            if (tail.size() < ov) continue;
            const std::string suffix = tail.substr(tail.size() - ov);
            auto it = byPrefix.find(suffix.substr(0, minOverlap));
            if (it == byPrefix.end()) continue;
            for (size_t j : it->second) {
                if (j == i || used[j] || used[i]) continue;
                if (ends[j].u == ends[i].u) continue;          // no self-join
                if (nodes[ends[j].u].deleted) continue;
                if (headSeq[j].size() < ov) continue;
                // Exact overlap only: a mismatch here is a different locus.
                if (headSeq[j].compare(0, ov, suffix) != 0) continue;
                if (lowComplexity(suffix)) continue;
                addLink(ends[i].u, ends[i].end, ends[j].u, ends[j].end);
                used[i] = used[j] = 1;
                ++added;
                joined = true;
                break;
            }
            if (ov == minOverlap) break;
        }
    }
    return added;
}

void UnitigGraph::simplify(double meanCoverage, int readLength, bool verbose,
                           double bubbleCoverageLimit, int maxRounds,
                           std::vector<SimplifyRoundStats>* rounds) {
    // Tip clipping is where TesserACT and SPAdes differ most. SPAdes clips tips up to
    //     max(min(k, read_length/2) * 3.5, read_length)
    // and removes any whose coverage is below 2.0x the best alternative (rctc 2.0).
    // TesserACT clips up to max(2k, read_length) at 0.19-0.40x: 1.7x shorter and 5-10x
    // more conservative. Unclipped tips are exactly what stops a graph resolving into
    // long unitigs, and they proliferate at low coverage -- where the measured gap lives.
    // Both bounds are overridable so the trade-off against genome fraction can be
    // measured rather than assumed.
    static const double kTipLenMult = [] {
        // Swept against SPAdes' tc_lb 3.5 on low-coverage strains. TesserACT's legacy
        // ceiling of max(2k, readLength) is ~1.7x shorter than SPAdes', and the LENGTH
        // bound -- not the coverage ratio -- is what dominates: on ERR11578724 the ratio
        // change alone gave 24,171 while ratio+length gave 353,308 against a legacy
        // 8,989, a 39x gain, with genome fraction RISING 95.48 -> 95.78.
        const char* e = std::getenv("TESSERA_TIP_LEN_MULT");
        return e ? std::atof(e) : 3.5;
    }();
    static const double kTipRatio = [] {
        // SPAdes uses rctc 2.0. Swept here, 1.0 beat 2.0 on the strain with the largest
        // gain (353,308 vs 173,770), so copying their constant would have left half the
        // improvement behind. Different strains prefer different values, which is the
        // argument for measuring rather than adopting.
        const char* e = std::getenv("TESSERA_TIP_RATIO");
        return e ? std::atof(e) : 1.0;
    }();
    size_t tipLen = static_cast<size_t>(std::max(2 * k_, readLength));
    if (kTipLenMult > 0) {
        const double kk = std::min<double>(k_, readLength / 2.0);
        tipLen = static_cast<size_t>(std::max(kk * kTipLenMult,
                                              static_cast<double>(readLength)));
    }
    const size_t bubbleLen = static_cast<size_t>(std::max(3 * k_, 2 * readLength));

    for (int round = 0; round < maxRounds; ++round) {
        // Aggressiveness ramps up across rounds: early passes only remove the
        // obviously spurious, so genuine low-coverage sequence survives long
        // enough to be joined into something defensible.
        const double ramp = 0.4 + 0.6 * static_cast<double>(round) /
                                  static_cast<double>(std::max(1, maxRounds - 1));
        SimplifyRoundStats st;
        st.round = round + 1;

        const double tipRatio = kTipRatio > 0 ? kTipRatio * ramp : 0.35 * ramp + 0.05;
        st.tipsRemoved = removeTips(tipLen, tipRatio);
        st.merged += compact();
        // A true error bubble sits well below even the half-depth a two-copy
        // split would give, which is what separates it from real divergence.
        st.bubblesPopped = popBubbles(bubbleLen, 0.95, meanCoverage * bubbleCoverageLimit);
        st.merged += compact();
        static const double chimeraFactor = [] {
            const char* e = std::getenv("TESSERA_CHIMERA_FACTOR");
            return e ? std::atof(e) : 0.12;
        }();
        // Cutting long low-coverage connectors risks severing real sequence,
        // so the operation stays confined to short ones. Raising this ceiling
        // was measured and cost genome fraction without buying contiguity.
        //
        // The THRESHOLD, though, was provably inert wherever it mattered. A unitig's
        // coverage can never fall below the k-mer abundance cutoff, which is 2 in 99% of
        // strain x k pairs; meanCoverage * 0.12 * ramp only clears 2 once meanCoverage
        // passes ~17 at full ramp and ~41 at the opening ramp. Measured across the
        // benchmark: 0 removals in 38/38 strains where peak*0.12 < 2, and those strains
        // have a median NGA50 ratio of 0.471 against SPAdes. Chimera removal is also the
        // single operation that most separates winners from losers at the top rung.
        //
        // So give it a floor just above the abundance cutoff. Below that the call is a
        // no-op dressed up as a simplification step.
        static const double chimeraFloor = [] {
            const char* e = std::getenv("TESSERA_CHIMERA_FLOOR");
            return e ? std::atof(e) : 3.0;
        }();
        const double chimeraCut = std::max(meanCoverage * chimeraFactor * ramp,
                                           chimeraFloor * ramp);
        // Length ceiling. SPAdes uses 2*max(5*min(k, RL/2), RL) - 1 -- at k=55/RL=150
        // that is 549 bp against TesserACT's 2k = 110, a 5x difference, and the largest
        // single numerical divergence between the two. Real chimeras and IS-boundary
        // mis-junctions run 200-500 bp, so a 110 bp ceiling cannot see them, and each
        // survivor is a contig break. Klebsiella carries dozens of IS copies and seven
        // rRNA operons, so those breaks are not rare.
        //
        // Raising this was measured before and cost genome fraction -- but that was with
        // the broken gate above, which cut connectors having no alternative path. With a
        // real outdeg/indeg gate the operation can only remove edges whose junctions both
        // survive, which is what makes the longer reach safe.
        static const double kEcLenMult = [] {
            const char* e = std::getenv("TESSERA_EC_LEN_MULT");
            return e ? std::atof(e) : 5.0;
        }();
        const double ecK = std::min<double>(k_, readLength / 2.0);
        const size_t ecLen = static_cast<size_t>(
            2.0 * std::max(kEcLenMult * ecK, static_cast<double>(readLength)) - 1.0);
        st.chimerasRemoved = removeErroneousConnections(chimeraCut, ecLen);
        st.merged += compact();
        st.isolatedRemoved = removeIsolated(meanCoverage * 0.25 * ramp, tipLen);

        // Local weakness, for the replicons a global cutoff cannot serve.
        // On by default. Over the whole 40-isolate panel: genome fraction
        // 98.83 -> 99.06, past SPAdes' 98.95, and of 78 plasmids those
        // recovered at 99% or better rise from 41 to 62 against SPAdes' 68,
        // with none left below half recovery and 78 of 78 past the 90% mark,
        // matching SPAdes. Misassemblies go 13 -> 16, still well under SPAdes'
        // 26, and mismatches 0.14 -> 0.21 against its 0.44 -- the extra
        // sequence recovered is repeat-adjacent and carries more of them.
        // Set TESSERA_LOCAL_WEAK=0 to turn it off.
        static const double localWeak = [] {
            const char* e = std::getenv("TESSERA_LOCAL_WEAK");
            return e ? std::atof(e) : 0.10;
        }();
        if (localWeak > 0) {
            st.lowDepthRemoved += removeLocallyWeak(localWeak, meanCoverage * 0.25);
            st.merged += compact();
        }

        // A missing k-mer breaks the graph where the sequence itself does not.
        // Rejoining loose ends that overlap exactly recovers edges no decision
        // rule could, because there was no edge to decide between.
        // Default was 0, so this never executed in any of 666 benchmark assemblies --
        // an implemented, tested recovery step that no user has ever run. Dead ends do
        // not clear across simplification (3,660 -> 1,895 over twelve rounds on a typical
        // strain), and an exact overlap between two loose ends is evidence no decision
        // rule can supply, because there was no edge to decide between. Default it on at
        // a conservative overlap; set 0 to restore the old behaviour.
        static const size_t joinOverlap = [] {
            const char* e = std::getenv("TESSERA_JOIN_DEADENDS");
            return e ? static_cast<size_t>(std::atoi(e)) : 31u;
        }();
        if (joinOverlap > 0) {
            st.deadEndsJoined = joinDeadEnds(joinOverlap);
            st.merged += compact();
        }
        st.merged += compact();

        st.unitigs = liveCount();
        st.n50 = n50();
        // A rising dead-end count means cleaning is fragmenting the graph
        // rather than tidying it, which is worth seeing in the report.
        st.deadEnds = totalDeadEnds();
        st.totalLength = totalLength();
        const size_t changes = st.tipsRemoved + st.bubblesPopped +
                               st.chimerasRemoved + st.isolatedRemoved;
        if (rounds) rounds->push_back(st);

        if (verbose) {
            std::fprintf(stderr,
                         "    round %-2d tips=%-6zu bubbles=%-5zu chimeras=%-5zu isolated=%-5zu"
                         "lowdepth=%-5zu deadends=%-5zu "
                         "  unitigs=%-8zu N50=%zu\n",
                         round + 1, st.tipsRemoved, st.bubblesPopped, st.chimerasRemoved,
                         st.isolatedRemoved, st.lowDepthRemoved, st.deadEnds,
                         st.unitigs, st.n50);
        }
        if (changes == 0) break;
    }
}

// ---------------------------------------------------------------------------
// Reporting
// ---------------------------------------------------------------------------

size_t UnitigGraph::liveCount() const {
    size_t n = 0;
    for (const Unitig& u : nodes) if (!u.deleted) ++n;
    return n;
}

size_t UnitigGraph::totalLength() const {
    size_t n = 0;
    for (const Unitig& u : nodes) if (!u.deleted) n += u.seq.size();
    return n;
}

size_t UnitigGraph::n50() const {
    std::vector<size_t> lens;
    lens.reserve(nodes.size());
    size_t total = 0;
    for (const Unitig& u : nodes) {
        if (u.deleted) continue;
        lens.push_back(u.seq.size());
        total += u.seq.size();
    }
    if (lens.empty()) return 0;
    std::sort(lens.begin(), lens.end(), std::greater<size_t>());
    size_t acc = 0;
    for (size_t l : lens) {
        acc += l;
        if (acc * 2 >= total) return l;
    }
    return lens.back();
}

double UnitigGraph::medianCoverage() const {
    std::vector<double> cov;
    for (const Unitig& u : nodes) {
        if (!u.deleted && u.seq.size() >= static_cast<size_t>(2 * k_)) cov.push_back(u.coverage);
    }
    if (cov.empty()) return 0;
    std::nth_element(cov.begin(), cov.begin() + static_cast<long>(cov.size() / 2), cov.end());
    return cov[cov.size() / 2];
}

void UnitigGraph::stats(const char* label, bool verbose) const {
    if (!verbose) return;
    std::fprintf(stderr, "  %-22s unitigs=%-8zu total=%-12zu N50=%-10zu medCov=%.1f\n",
                 label, liveCount(), totalLength(), n50(), medianCoverage());
}

void UnitigGraph::toContigs(size_t minLen, std::vector<std::string>& seqs,
                            std::vector<double>& covs) const {
    for (const Unitig& u : nodes) {
        if (u.deleted || u.seq.size() < minLen) continue;
        seqs.push_back(u.seq);
        covs.push_back(u.coverage);
    }
}

std::string UnitigGraph::validate() const {
    for (uint32_t u = 0; u < nodes.size(); ++u) {
        if (nodes[u].deleted) {
            if (!nodes[u].ends[0].empty() || !nodes[u].ends[1].empty()) {
                return "deleted node " + std::to_string(u) + " still has links";
            }
            continue;
        }
        if (nodes[u].seq.size() < static_cast<size_t>(k_)) {
            return "node " + std::to_string(u) + " is shorter than k";
        }
        for (int e = 0; e < 2; ++e) {
            for (const Link& l : nodes[u].ends[e]) {
                if (l.to >= nodes.size()) return "node " + std::to_string(u) + " links out of range";
                if (nodes[l.to].deleted) {
                    return "node " + std::to_string(u) + " links to deleted node " +
                           std::to_string(l.to);
                }
                const auto& back = nodes[l.to].ends[l.toEnd];
                if (std::find(back.begin(), back.end(), Link{u, static_cast<uint8_t>(e)}) == back.end()) {
                    return "asymmetric link " + std::to_string(u) + ":" + std::to_string(e) +
                           " -> " + std::to_string(l.to) + ":" + std::to_string(l.toEnd);
                }
            }
        }
    }
    return "";
}

}  // namespace ts
