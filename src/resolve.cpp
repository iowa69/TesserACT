#include "resolve.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <mutex>
#include <cstdlib>
#include <thread>

namespace ts {

namespace {

// The distance enumeration may search when it is only establishing which
// continuations exist. Repeats in an enterobacterial genome run to a few
// kilobases (IS elements ~0.8-2.5 kb, rRNA operons ~5 kb); this covers them
// without letting the search wander across the chromosome.
constexpr double kTopologyReach = 6000.0;

constexpr int kMaxProbes = 12;        // k-mer lookups per read before deciding
constexpr int kMinVotes = 2;
constexpr int kAnchorK = 31;
constexpr uint64_t kAmbiguous = UINT64_MAX;

inline uint64_t packIndex(uint32_t unitig, uint32_t pos, int strand) {
    return (static_cast<uint64_t>(unitig) << 33) | (static_cast<uint64_t>(pos) << 1) |
           static_cast<uint64_t>(strand & 1);
}
inline uint32_t idxUnitig(uint64_t v) { return static_cast<uint32_t>(v >> 33); }
inline uint32_t idxPos(uint64_t v) { return static_cast<uint32_t>((v >> 1) & 0xFFFFFFFFULL); }
inline int idxStrand(uint64_t v) { return static_cast<int>(v & 1); }

}  // namespace

PairedResolver::PairedResolver(const UnitigGraph& graph, const SequenceStore& reads, int threads,
                               int minLinkSupport, double tieRatio, int minScaffoldSupport)
    : g_(graph), reads_(reads), threads_(threads > 0 ? threads : 1), k_(graph.k()),
      kMap_(std::min(kAnchorK, graph.k())),
      minLinkSupport_(minLinkSupport), tieRatio_(tieRatio),
      minScaffoldSupport_(minScaffoldSupport) {
    medianCoverage_ = graph.medianCoverage();
}

void PairedResolver::buildIndex() {
    size_t totalKmers = 0;
    for (const Unitig& u : g_.nodes) {
        if (!u.deleted && u.seq.size() >= static_cast<size_t>(kMap_)) {
            totalKmers += u.seq.size() - static_cast<size_t>(kMap_) + 1;
        }
    }
    index_.reserve(totalKmers * 2);

    for (uint32_t i = 0; i < g_.nodes.size(); ++i) {
        const Unitig& u = g_.nodes[i];
        if (u.deleted || u.seq.size() < static_cast<size_t>(kMap_)) continue;
        Kmer fwd = 0, rc = 0;
        int valid = 0;
        for (uint32_t p = 0; p < u.seq.size(); ++p) {
            int c = baseCode(u.seq[p]);
            if (c < 0) { valid = 0; continue; }
            fwd = pushBack(fwd, c, kMap_);
            rc = pushFrontRc(rc, c, kMap_);
            if (++valid < kMap_) continue;
            const uint32_t start = p + 1 - static_cast<uint32_t>(kMap_);
            const Kmer canon = fwd < rc ? fwd : rc;
            // strand flag records whether the unitig's forward k-mer is the
            // canonical one, so a read hit can be resolved to an orientation.
            const int strand = (fwd == canon) ? 0 : 1;
            auto it = index_.find(canon);
            if (it == index_.end()) index_.emplace(canon, packIndex(i, start, strand));
            else it->second = kAmbiguous;   // occurs elsewhere; cannot locate a read
        }
    }
}

Anchor PairedResolver::anchorRead(size_t read) const {
    Anchor best;
    const int len = static_cast<int>(reads_.length(read));
    if (len < kMap_) return best;

    const int span = len - kMap_;
    const int probes = std::min(kMaxProbes, span + 1);

    struct Vote { uint32_t unitig; int32_t pos; uint8_t orient; int count; };
    Vote votes[kMaxProbes];
    int distinct = 0;

    // Probes are spread across the read because a sequencing error invalidates
    // every k-mer overlapping it; one error must not veto the whole anchor.
    for (int t = 0; t < probes; ++t) {
        const int rp = probes == 1 ? 0 : span * t / (probes - 1);
        Kmer fwd = 0, rcv = 0;
        bool ok = true;
        for (int j = 0; j < kMap_; ++j) {
            int c = reads_.baseAt(read, static_cast<uint32_t>(rp + j));
            if (c < 0) { ok = false; break; }
            fwd = pushBack(fwd, c, kMap_);
            rcv = pushFrontRc(rcv, c, kMap_);
        }
        if (!ok) continue;

        const Kmer canon = fwd < rcv ? fwd : rcv;
        auto it = index_.find(canon);
        if (it == index_.end() || it->second == kAmbiguous) continue;

        const uint32_t u = idxUnitig(it->second);
        const int up = static_cast<int>(idxPos(it->second));
        const int sflag = idxStrand(it->second);
        const int rflag = (fwd == canon) ? 0 : 1;
        const int orient = rflag ^ sflag;

        // A read spanning a junction legitimately hangs off the unitig, so the
        // position is deliberately allowed to fall outside [0, unitigLen).
        const int32_t startPos = (orient == 0)
                                     ? static_cast<int32_t>(up - rp)
                                     : static_cast<int32_t>(up - (len - rp - kMap_));

        int found = -1;
        for (int q = 0; q < distinct; ++q) {
            if (votes[q].unitig == u && votes[q].pos == startPos &&
                votes[q].orient == static_cast<uint8_t>(orient)) { found = q; break; }
        }
        if (found >= 0) ++votes[found].count;
        else if (distinct < kMaxProbes) {
            votes[distinct++] = {u, startPos, static_cast<uint8_t>(orient), 1};
        }
    }

    int bestIdx = -1, bestVotes = 0;
    for (int q = 0; q < distinct; ++q) {
        if (votes[q].count > bestVotes) { bestVotes = votes[q].count; bestIdx = q; }
    }
    if (bestIdx < 0 || bestVotes < kMinVotes) return best;

    best.unitig = votes[bestIdx].unitig;
    best.pos = votes[bestIdx].pos;
    best.orient = votes[bestIdx].orient;
    return best;
}

void PairedResolver::buildSupport() {
    buildIndex();
    if (!reads_.paired()) return;

    const size_t pairs = reads_.pairCount();
    std::vector<std::vector<int>> insertSamples(static_cast<size_t>(threads_));
    std::vector<std::unordered_map<uint64_t, std::unordered_map<uint64_t, std::vector<int32_t>>>>
        localSupport(static_cast<size_t>(threads_));
    std::atomic<size_t> mappedCount{0}, linkCount{0};

    auto worker = [&](int tid) {
        auto& samples = insertSamples[static_cast<size_t>(tid)];
        auto& sup = localSupport[static_cast<size_t>(tid)];
        size_t localMapped = 0, localLinks = 0;

        for (size_t p = static_cast<size_t>(tid); p < pairs; p += static_cast<size_t>(threads_)) {
            const size_t i1 = p * 2, i2 = p * 2 + 1;
            const Anchor a1 = anchorRead(i1);
            if (!a1.mapped()) continue;
            const Anchor a2 = anchorRead(i2);
            if (!a2.mapped()) continue;
            localMapped += 2;

            const int len1 = static_cast<int>(reads_.length(i1));
            const int len2 = static_cast<int>(reads_.length(i2));

            if (a1.unitig == a2.unitig) {
                // Same unitig: a properly oriented pair measures the fragment
                // length directly, which is how the insert model is learned.
                if (a1.orient == 0 && a2.orient == 1 && a2.pos + len2 > a1.pos) {
                    samples.push_back(a2.pos + len2 - a1.pos);
                } else if (a2.orient == 0 && a1.orient == 1 && a1.pos + len1 > a2.pos) {
                    samples.push_back(a1.pos + len1 - a2.pos);
                }
                continue;
            }

            const int eu = (a1.orient == 0) ? 1 : 0;
            const int ev = (a2.orient == 1) ? 0 : 1;
            // Orientation of each unitig as traversed by the fragment.
            const int dU = (eu == 1) ? 0 : 1;
            const int dV = (ev == 0) ? 0 : 1;

            const int lenU = static_cast<int>(g_.nodes[a1.unitig].seq.size());
            const int lenV = static_cast<int>(g_.nodes[a2.unitig].seq.size());
            // Bases of the fragment that fall inside each unitig, measured from
            // the read's 5' end to the boundary the fragment crosses.
            const int consumedU = (a1.orient == 0) ? (lenU - a1.pos) : (a1.pos + len1);
            const int consumedV = (dV == 0) ? (a2.pos + len2) : (lenV - a2.pos);
            // Concatenating two unitigs merges their shared (k-1) overlap.
            const int32_t span = static_cast<int32_t>(consumedU + consumedV - (k_ - 1));

            sup[orientedId(a1.unitig, dU)][orientedId(a2.unitig, dV)].push_back(span);
            // The mirrored traversal must be recorded too, so a path arriving
            // from the other side sees the same evidence.
            sup[orientedId(a2.unitig, 1 - dV)][orientedId(a1.unitig, 1 - dU)].push_back(span);
            ++localLinks;
        }
        mappedCount += localMapped;
        linkCount += localLinks;
    };

    std::vector<std::thread> pool;
    for (int t = 0; t < threads_; ++t) pool.emplace_back(worker, t);
    for (auto& th : pool) th.join();

    for (auto& sup : localSupport) {
        for (auto& kv : sup) {
            auto& dst = support_[kv.first];
            for (auto& kv2 : kv.second) {
                auto& v = dst[kv2.first];
                v.insert(v.end(), kv2.second.begin(), kv2.second.end());
            }
        }
    }

    std::vector<int> all;
    for (auto& s : insertSamples) all.insert(all.end(), s.begin(), s.end());
    for (int v : all) {
        if (v < 0) continue;
        const size_t bin = static_cast<size_t>(v);
        if (bin >= insertHistogram_.size()) insertHistogram_.resize(bin + 1, 0);
        ++insertHistogram_[bin];
    }
    if (all.size() >= 1000) {
        std::sort(all.begin(), all.end());
        // Trim the tails before fitting: chimeric pairs and mismapped reads sit
        // far out and would inflate the standard deviation badly.
        const size_t lo = all.size() / 100;
        const size_t hi = all.size() - all.size() / 100;
        double sum = 0;
        for (size_t i = lo; i < hi; ++i) sum += all[i];
        const double n = static_cast<double>(hi - lo);
        insert_.mean = sum / n;
        double var = 0;
        for (size_t i = lo; i < hi; ++i) {
            const double d = all[i] - insert_.mean;
            var += d * d;
        }
        insert_.stddev = std::sqrt(var / n);
        insert_.observations = all.size();
        insert_.minPlausible = std::max(0, static_cast<int>(insert_.mean - 4 * insert_.stddev));
        insert_.maxPlausible = static_cast<int>(insert_.mean + 4 * insert_.stddev);
        insert_.usable = true;
    }

    size_t distinct = 0;
    for (auto& kv : support_) distinct += kv.second.size();
    stats_.readsMapped = mappedCount.load();
    stats_.pairsLinking = linkCount.load();
    stats_.distinctLinks = distinct;
    stats_.insert = insert_;
}

double PairedResolver::scoreCandidate(uint64_t from, uint64_t terminal, int interLen) const {
    auto it = support_.find(from);
    if (it == support_.end()) return 0;
    auto jt = it->second.find(terminal);
    if (jt == it->second.end()) return 0;
    if (!insert_.usable) return static_cast<double>(jt->second.size());

    // A pair supports this particular path only if the fragment length it
    // implies once the path's intermediate sequence is inserted is one the
    // library could actually have produced.
    double score = 0;
    for (int32_t span : jt->second) {
        const int implied = span + interLen;
        if (implied >= insert_.minPlausible && implied <= insert_.maxPlausible) score += 1;
    }
    return score;
}

void PairedResolver::resolve(std::vector<std::string>& contigs, std::vector<double>& covs) {
    const size_t n = g_.nodes.size();
    const size_t ov = static_cast<size_t>(k_ - 1);

    // A collapsed repeat may be traversed more than once and carries no unique
    // paired evidence, so it cannot seed a chain, terminate a path, or be
    // traversed -- it is emitted verbatim and walls off both neighbours.
    //
    // Depth alone does not identify one. A plasmid at two copies per cell sits
    // at the same ~2x median as a two-copy repeat while its k-mers are
    // perfectly unique, and a bacterial isolate usually carries several
    // plasmids; on one panel isolate 4.8% of the assembly was excluded this
    // way, including a 68 kb unitig at 1.78x median that was the whole of that
    // assembly's NGA50. Raising the depth cutoff only trades that error for the
    // opposite one, and stops two-copy repeats resolving.
    //
    // Two refinements were tried and measured on the closed-reference panel,
    // and both were rejected. Raising the multiplier to 2.4 frees the plasmids
    // but stops genuine two-copy repeats resolving. Additionally requiring a
    // branching end -- on the theory that distinct copies of a repeat are
    // flanked by different sequence -- left contig NGA50 unchanged to the base
    // pair on both datasets tested while adding misassemblies (3 -> 5 and
    // 2 -> 3), because a collapsed tandem array has a single link at each end
    // and was then traversed once instead of many times.
    //
    // So this stays depth-only. The plasmid cost is real but is the cheaper
    // error, and the contiguity it was meant to buy came from the scaffolding
    // and tie-break changes instead.
    const double repeatThreshold = medianCoverage_ * 1.6;
    auto isRepeat = [&](uint32_t u) {
        return medianCoverage_ > 0 && g_.nodes[u].coverage > repeatThreshold;
    };

    // Reads only anchor uniquely inside a unitig long enough to own k-mers its
    // neighbours do not share across the (k-1) junction overlaps. Shorter
    // unitigs carry no paired evidence, so extension has to look *through* them
    // to the next unitig that does.
    const size_t anchorableLen = static_cast<size_t>(2 * k_);
    auto anchorable = [&](uint32_t u) {
        return !g_.nodes[u].deleted && g_.nodes[u].seq.size() >= anchorableLen && !isRepeat(u);
    };

    const double reach = insert_.usable ? insert_.maxPlausible : 1000.0;
    constexpr size_t kMaxNodes = 10;
    constexpr double kForcedReachFactor = 1.5;
    constexpr size_t kMaxCandidates = 64;

    // How far enumeration may look when it is only asking *what continuations
    // exist*, as opposed to asking the read pairs to choose between them.
    //
    // Tying the search itself to fragment reach conflates the two. On libraries
    // whose mates overlap -- most of the closed-reference panel, where the
    // fragment is ~230 bp against ~190 bp reads -- reach is a few hundred bases,
    // so any repeat longer than that ends the search before a single candidate
    // is generated: on one isolate 2,620 of 3,288 chain ends returned no
    // candidate at all, against 269 successful joins. The paired evidence never
    // got the chance to be insufficient; it was never consulted.
    //
    // Searching further does not weaken any decision. A candidate beyond
    // fragment reach collects no paired support and so cannot win a contested
    // choice -- it can only be taken when it is the sole way through, which is
    // a statement about the graph's topology that needs no reads to back it.
    const double searchBudget = std::max(reach, kTopologyReach);

    auto flip = [](uint64_t oid) { return orientedId(unitigOf(oid), 1 - orientOf(oid)); };
    auto addedLen = [&](uint64_t oid) { return g_.nodes[unitigOf(oid)].seq.size() - ov; };

    // Every way out of `from` that terminates on an anchorable unitig within
    // fragment reach; intermediate nodes are unanchorable repeats.
    auto enumerate = [&](uint64_t from, std::vector<std::vector<uint64_t>>& out) {
        // `forced` marks a path that has never had an alternative. Paired reads
        // cannot vouch for anything beyond one fragment length, but where the
        // graph offers no other way through, length is irrelevant -- so a forced
        // path may run past `reach` and cross repeats no pair could span.
        struct Frame { std::vector<uint64_t> nodes; size_t len; bool forced; };
        std::vector<Frame> work;
        const auto& firstExits = g_.exits(unitigOf(from), orientOf(from));
        for (const Link& l : firstExits) {
            if (g_.nodes[l.to].deleted) continue;
            work.push_back({{orientedId(l.to, UnitigGraph::enterOrient(l))},
                            g_.nodes[l.to].seq.size() - ov, firstExits.size() == 1});
        }
        while (!work.empty() && out.size() < kMaxCandidates) {
            Frame f = std::move(work.back());
            work.pop_back();
            const uint64_t tail = f.nodes.back();
            const uint32_t tu = unitigOf(tail);
            if (anchorable(tu)) { out.push_back(std::move(f.nodes)); continue; }
            if (f.nodes.size() >= kMaxNodes) continue;
            // A forced path (one the graph offered no alternative to) may run a
            // little past what a fragment can span, but not far: measured on the
            // benchmark panel, letting it run unbounded bought only a lower
            // contig count and cost a misassembly on K. pneumoniae.
            const double budget = f.forced ? searchBudget * kForcedReachFactor : searchBudget;
            if (static_cast<double>(f.len) > budget) continue;
            const auto& exits = g_.exits(tu, orientOf(tail));
            for (const Link& l : exits) {
                if (g_.nodes[l.to].deleted) continue;
                const uint64_t nxt = orientedId(l.to, UnitigGraph::enterOrient(l));
                bool seen = false;
                for (uint64_t o : f.nodes) if (o == nxt) { seen = true; break; }
                if (seen) continue;
                Frame g2 = f;
                g2.nodes.push_back(nxt);
                g2.len += g_.nodes[l.to].seq.size() - ov;
                g2.forced = f.forced && exits.size() == 1;
                work.push_back(std::move(g2));
            }
        }
    };

    // Start with every anchorable unitig as a chain of one, then repeatedly
    // join chains whose paired evidence mutually prefers each other. Growing
    // chains lets support accumulate over a whole contig tail rather than just
    // the last unitig, which is what carries a walk past a repeat.
    std::vector<std::vector<uint64_t>> chains;
    for (uint32_t u = 0; u < n; ++u) {
        if (anchorable(u)) chains.push_back({orientedId(u, 0)});
    }

    std::vector<uint32_t> ownerChain(n, UINT32_MAX);
    std::vector<uint32_t> ownerPos(n, 0);
    auto reindex = [&]() {
        std::fill(ownerChain.begin(), ownerChain.end(), UINT32_MAX);
        for (uint32_t c = 0; c < chains.size(); ++c) {
            for (uint32_t j = 0; j < chains[c].size(); ++j) {
                const uint32_t u = unitigOf(chains[c][j]);
                ownerChain[u] = c;
                ownerPos[u] = j;
            }
        }
    };

    struct Cont {
        uint32_t chainB = UINT32_MAX;
        int endB = 0;
        std::vector<uint64_t> connector;
        double score = 0;
        bool ok = false;
    };

    // Best continuation off one end of a chain, scored with every member of
    // that chain that still lies within fragment reach of the boundary.
    long long dbgNoCand = 0, dbgLowSupport = 0, dbgTie = 0, dbgMidChain = 0, dbgOk = 0;
    long long dbgTieBest = 0, dbgTieSecond = 0;

    auto bestContinuation = [&](uint32_t c, int end) {
        Cont result;
        const std::vector<uint64_t>& ch = chains[c];
        if (ch.empty()) return result;

        // Orient the chain so the end under consideration is the tail.
        std::vector<uint64_t> tailFirst;
        if (end == 1) tailFirst.assign(ch.begin(), ch.end());
        else {
            tailFirst.reserve(ch.size());
            for (size_t i = ch.size(); i-- > 0;) tailFirst.push_back(flip(ch[i]));
        }

        // Distance from the far edge of each member to the chain boundary.
        std::vector<size_t> distToEnd(tailFirst.size(), 0);
        size_t acc = 0;
        for (size_t i = tailFirst.size(); i-- > 0;) {
            distToEnd[i] = acc;
            acc += addedLen(tailFirst[i]);
        }

        std::vector<std::vector<uint64_t>> cands;
        enumerate(tailFirst.back(), cands);
        if (cands.empty()) { ++dbgNoCand; return result; }

        double best = -1, second = -1;
        size_t pick = 0;
        std::vector<double> scores(cands.size(), 0.0);
        for (size_t i = 0; i < cands.size(); ++i) {
            int interLen = 0;
            for (size_t j = 0; j + 1 < cands[i].size(); ++j) {
                interLen += static_cast<int>(g_.nodes[unitigOf(cands[i][j])].seq.size()) - (k_ - 1);
            }
            const uint64_t terminal = cands[i].back();
            double sc = 0;
            for (size_t m = tailFirst.size(); m-- > 0;) {
                if (static_cast<double>(distToEnd[m]) > reach) break;
                sc += scoreCandidate(tailFirst[m], terminal,
                                     interLen + static_cast<int>(distToEnd[m]));
            }
            scores[i] = sc;
            if (sc > best) { second = best; best = sc; pick = i; }
            else if (sc > second) second = sc;
        }

        if (cands.size() > 1) {
            if (best < minLinkSupport_) { ++dbgLowSupport; return result; }
            // A near tie means the repeat is genuinely unresolved; guessing
            // would manufacture a misassembly.
            if (second > 0 && best < tieRatio_ * second) {
                // Unless every near-tied path ends on the same unitig in the
                // same orientation. Then the destination is not in doubt at
                // all -- the paired evidence cannot separate them precisely
                // because it supports the identical join -- and only the route
                // between is ambiguous. Refusing here throws away a join that
                // is certain; take it, and settle the interior by coverage,
                // which is the evidence that does distinguish the arms.
                const uint64_t bestTerm = cands[pick].back();
                bool sameDestination = true;
                for (size_t i = 0; i < cands.size(); ++i) {
                    if (scores[i] * tieRatio_ < best) continue;   // not near-tied
                    if (cands[i].back() != bestTerm) { sameDestination = false; break; }
                }
                if (!sameDestination) {
                    ++dbgTie;
                    dbgTieBest += static_cast<long long>(best);
                    dbgTieSecond += static_cast<long long>(second);
                    return result;
                }
                double bestInterior = -1.0;
                for (size_t i = 0; i < cands.size(); ++i) {
                    if (scores[i] * tieRatio_ < best) continue;
                    double minCov = std::numeric_limits<double>::max();
                    for (size_t j = 0; j + 1 < cands[i].size(); ++j) {
                        minCov = std::min(minCov, g_.nodes[unitigOf(cands[i][j])].coverage);
                    }
                    // A path with no interior is the most direct route.
                    if (cands[i].size() < 2) minCov = std::numeric_limits<double>::max();
                    if (minCov > bestInterior) { bestInterior = minCov; pick = i; }
                }
            }
        }

        // The terminal has to be an *end* of its chain; landing in the middle
        // would mean cutting an already-supported contig in half.
        const uint64_t terminal = cands[pick].back();
        const uint32_t tu = unitigOf(terminal);
        const uint32_t cb = ownerChain[tu];
        if (cb == UINT32_MAX || cb == c) return result;
        const std::vector<uint64_t>& other = chains[cb];
        const uint32_t j = ownerPos[tu];
        const int storedOrient = orientOf(other[j]);

        if (storedOrient == orientOf(terminal) && j == 0) result.endB = 0;
        else if (storedOrient != orientOf(terminal) && j + 1 == other.size()) result.endB = 1;
        else { ++dbgMidChain; return result; }
        ++dbgOk;

        result.chainB = cb;
        result.connector.assign(cands[pick].begin(), cands[pick].end() - 1);
        result.score = best;
        result.ok = true;
        return result;
    };

    for (int round = 0; round < 24; ++round) {
        reindex();
        const size_t nc = chains.size();
        std::vector<Cont> cont(nc * 2);
        for (uint32_t c = 0; c < nc; ++c) {
            if (chains[c].empty()) continue;
            cont[c * 2 + 0] = bestContinuation(c, 0);
            cont[c * 2 + 1] = bestContinuation(c, 1);
        }

        std::vector<char> merged(nc, 0);
        size_t joins = 0;
        for (uint32_t c = 0; c < nc; ++c) {
            if (merged[c] || chains[c].empty()) continue;
            for (int e = 0; e < 2; ++e) {
                const Cont& f = cont[c * 2 + e];
                if (!f.ok || merged[f.chainB] || merged[c]) continue;
                // Only join when both ends independently chose each other.
                const Cont& back = cont[f.chainB * 2 + f.endB];
                if (!back.ok || back.chainB != c || back.endB != e) continue;

                std::vector<uint64_t> a;
                if (e == 1) a = chains[c];
                else {
                    a.reserve(chains[c].size());
                    for (size_t i = chains[c].size(); i-- > 0;) a.push_back(flip(chains[c][i]));
                }
                std::vector<uint64_t> b;
                if (f.endB == 0) b = chains[f.chainB];
                else {
                    b.reserve(chains[f.chainB].size());
                    for (size_t i = chains[f.chainB].size(); i-- > 0;) {
                        b.push_back(flip(chains[f.chainB][i]));
                    }
                }

                a.insert(a.end(), f.connector.begin(), f.connector.end());
                a.insert(a.end(), b.begin(), b.end());
                chains[c].swap(a);
                chains[f.chainB].clear();
                merged[c] = 1;
                merged[f.chainB] = 1;
                ++joins;
                break;
            }
        }
        if (joins == 0) break;
    }

    if (getenv("TESSERA_DEBUG_RESOLVE")) {
        std::fprintf(stderr,
                     "      [debug] continuation outcomes: ok=%lld no-candidate=%lld "
                     "low-support=%lld tie=%lld mid-chain=%lld  (tie mean best=%.1f second=%.1f)\n",
                     dbgOk, dbgNoCand, dbgLowSupport, dbgTie, dbgMidChain,
                     dbgTie ? static_cast<double>(dbgTieBest) / static_cast<double>(dbgTie) : 0.0,
                     dbgTie ? static_cast<double>(dbgTieSecond) / static_cast<double>(dbgTie) : 0.0);
    }

    // ---- scaffolding -----------------------------------------------------
    // Chain ends with strong paired support but no path through the graph are
    // joined across a gap of Ns whose length comes from the fragment model.
    // Each chain has two ports; a join consumes one port at each end, so the
    // accepted joins form simple paths that are then walked in order.
    const size_t nc = chains.size();
    std::vector<uint32_t> joinTo(nc * 2, UINT32_MAX);   // port -> port
    std::vector<int> joinGap(nc * 2, 0);

    if (scaffolding_ && insert_.usable) {
        reindex();
        struct End { uint64_t oriented; std::vector<uint64_t> members; std::vector<size_t> dist; };
        std::vector<End> ends(nc * 2);
        for (uint32_t c = 0; c < nc; ++c) {
            if (chains[c].empty()) continue;
            for (int e = 0; e < 2; ++e) {
                std::vector<uint64_t> tailFirst;
                if (e == 1) tailFirst.assign(chains[c].begin(), chains[c].end());
                else {
                    for (size_t i = chains[c].size(); i-- > 0;) tailFirst.push_back(flip(chains[c][i]));
                }
                End& E = ends[c * 2 + e];
                E.oriented = tailFirst.back();
                size_t acc = 0;
                for (size_t i = tailFirst.size(); i-- > 0;) {
                    if (static_cast<double>(acc) > reach) break;
                    E.members.push_back(tailFirst[i]);
                    E.dist.push_back(acc);
                    acc += addedLen(tailFirst[i]);
                }
            }
        }

        // The port an oriented unitig represents when a fragment enters it.
        std::unordered_map<uint64_t, uint32_t> entryOf;
        for (uint32_t c = 0; c < nc; ++c) {
            if (chains[c].empty()) continue;
            entryOf[chains[c].front()] = c * 2 + 0;
            entryOf[flip(chains[c].back())] = c * 2 + 1;
        }

        struct Best { uint32_t partner = UINT32_MAX; double score = 0; int gap = 0; };
        std::vector<Best> best(ends.size());
        for (size_t idx = 0; idx < ends.size(); ++idx) {
            const End& E = ends[idx];
            if (E.members.empty()) continue;
            // Anything reachable through the graph was already handled by chain
            // extension; scaffolding only spans true gaps.
            std::vector<std::vector<uint64_t>> reachable;
            enumerate(E.oriented, reachable);
            std::unordered_map<uint64_t, char> viaGraph;
            for (const auto& r : reachable) viaGraph[r.back()] = 1;

            std::unordered_map<uint32_t, std::vector<int>> gaps;
            for (size_t m = 0; m < E.members.size(); ++m) {
                auto it = support_.find(E.members[m]);
                if (it == support_.end()) continue;
                for (const auto& kv : it->second) {
                    if (viaGraph.count(kv.first)) continue;
                    auto pt = entryOf.find(kv.first);
                    if (pt == entryOf.end()) continue;
                    if (pt->second / 2 == idx / 2) continue;   // same chain
                    for (int32_t span : kv.second) {
                        const int gap = static_cast<int>(insert_.mean) - span -
                                        static_cast<int>(E.dist[m]);
                        if (gap < -static_cast<int>(ov) || gap > insert_.maxPlausible) continue;
                        gaps[pt->second].push_back(gap);
                    }
                }
            }
            for (auto& kv : gaps) {
                if (static_cast<double>(kv.second.size()) > best[idx].score) {
                    std::sort(kv.second.begin(), kv.second.end());
                    best[idx].score = static_cast<double>(kv.second.size());
                    best[idx].partner = kv.first;
                    best[idx].gap = kv.second[kv.second.size() / 2];
                }
            }
        }

        // How many spanning pairs a join needs. A fixed count cannot be right
        // at both ends of the depth range: the number of pairs crossing a
        // junction scales with coverage, so a flat 5 is nearly unreachable on a
        // 40x library and trivially met on a 200x one. On one panel isolate 153
        // of the 173 dead ends that found a partner were rejected on this
        // threshold alone -- the mutual-best rule turned away none of them.
        // Scaling with the observed depth keeps the evidence bar constant in
        // the units that matter, with a floor of 3 so a shallow library still
        // needs corroboration.
        const double minScaffoldSupport =
            minScaffoldSupport_ > 0
                ? static_cast<double>(minScaffoldSupport_)
                : std::min(10.0, std::max(3.0, medianCoverage_ * 0.06));
        for (uint32_t idx = 0; idx < best.size(); ++idx) {
            const Best& b = best[idx];
            if (b.partner == UINT32_MAX || b.score < minScaffoldSupport) continue;
            if (best[b.partner].partner != idx) continue;      // must be mutual
            if (joinTo[idx] != UINT32_MAX || joinTo[b.partner] != UINT32_MAX) continue;
            joinTo[idx] = b.partner;
            joinTo[b.partner] = idx;
            joinGap[idx] = std::max(1, b.gap);
            joinGap[b.partner] = joinGap[idx];
        }
    }

    // ---- emit ------------------------------------------------------------
    std::vector<char> placed(n, 0);
    std::vector<char> done(nc, 0);

    ResolvedPath curPath;
    auto renderChain = [&](uint32_t c, bool flipped, std::string& seq,
                           double& covWeighted, size_t& covLen) {
        const size_t sz = chains[c].size();
        for (size_t i = 0; i < sz; ++i) {
            const uint64_t oid = flipped ? flip(chains[c][sz - 1 - i]) : chains[c][i];
            curPath.oriented.push_back(oid);
            curPath.gaps.push_back(0);
            const uint32_t u = unitigOf(oid);
            const std::string piece = g_.oriented(u, orientOf(oid));
            if (seq.empty()) seq = piece;
            else seq += piece.substr(ov);
            covWeighted += g_.nodes[u].coverage * static_cast<double>(piece.size());
            covLen += piece.size();
            placed[u] = 1;
        }
        if (sz > 1) stats_.unitigsJoined += sz - 1;
    };

    for (uint32_t c = 0; c < nc; ++c) {
        if (chains[c].empty() || done[c]) continue;
        // Start scaffolds at a free port so each path is walked from one end.
        if (joinTo[c * 2 + 0] != UINT32_MAX && joinTo[c * 2 + 1] != UINT32_MAX) continue;

        std::string seq;
        double covWeighted = 0;
        size_t covLen = 0;
        curPath = ResolvedPath();
        uint32_t cur = c;
        // Entering through the port that has no join leaves the other free to
        // continue; a chain joined only at port 1 is therefore used forward.
        int inPort = (joinTo[c * 2 + 0] == UINT32_MAX) ? 0 : 1;
        size_t steps = 0;

        int pendingGap = 0;
        while (true) {
            done[cur] = 1;
            const size_t before = curPath.oriented.size();
            renderChain(cur, inPort == 1, seq, covWeighted, covLen);
            if (pendingGap > 0 && before < curPath.gaps.size()) curPath.gaps[before] = pendingGap;
            pendingGap = 0;
            const uint32_t outPort = cur * 2 + (inPort == 0 ? 1 : 0);
            const uint32_t partner = joinTo[outPort];
            if (partner == UINT32_MAX || ++steps > nc) break;
            const uint32_t nextC = partner / 2;
            if (done[nextC] || chains[nextC].empty()) break;
            seq.append(static_cast<size_t>(joinGap[outPort]), 'N');
            stats_.gapBases += static_cast<size_t>(joinGap[outPort]);
            ++stats_.scaffoldJoins;
            pendingGap = joinGap[outPort];
            cur = nextC;
            inPort = static_cast<int>(partner % 2);
        }

        if (seq.empty()) continue;
        contigs.push_back(std::move(seq));
        covs.push_back(covLen ? covWeighted / static_cast<double>(covLen) : 0);
        paths_.push_back(curPath);
        ++stats_.pathsBuilt;
    }

    // Chains inside a scaffold cycle have no free port; break them arbitrarily.
    for (uint32_t c = 0; c < nc; ++c) {
        if (chains[c].empty() || done[c]) continue;
        std::string seq;
        double covWeighted = 0;
        size_t covLen = 0;
        done[c] = 1;
        curPath = ResolvedPath();
        renderChain(c, false, seq, covWeighted, covLen);
        if (seq.empty()) continue;
        contigs.push_back(std::move(seq));
        covs.push_back(covLen ? covWeighted / static_cast<double>(covLen) : 0);
        paths_.push_back(curPath);
        ++stats_.pathsBuilt;
    }

    // Anything no chain walked through would otherwise be lost.
    for (uint32_t u = 0; u < n; ++u) {
        if (g_.nodes[u].deleted || placed[u]) continue;
        contigs.push_back(g_.nodes[u].seq);
        covs.push_back(g_.nodes[u].coverage);
        ResolvedPath solo;
        solo.oriented.push_back(orientedId(u, 0));
        solo.gaps.push_back(0);
        paths_.push_back(solo);
    }
}

}  // namespace ts
