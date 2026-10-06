#include "p2_ends.h"

#include <algorithm>
#include <atomic>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <map>
#include <thread>
#include <unordered_map>

#include "p2_json.h"

namespace ts {
namespace p2 {

namespace {

constexpr int kProbes = 12;        // ContigEndLinks / polisher anchoring constants
constexpr int kMinVotes = 2;
constexpr uint64_t kAmbig = UINT64_MAX;

inline uint64_t packPlace(uint32_t rec, uint32_t pos, int strand) {
    return (static_cast<uint64_t>(rec) << 33) | (static_cast<uint64_t>(pos) << 1) | static_cast<uint64_t>(strand & 1);
}

struct Placed {
    uint32_t rec = UINT32_MAX;
    int32_t pos = 0;
    uint8_t orient = 0;     // 0: read runs along the record forward (points at the right end)
    uint32_t len = 0;
    bool ok() const { return rec != UINT32_MAX; }
};

// Canonical 31-mer of a stored read at [p, p+31), or false when a base is masked/ambiguous.
inline bool readKmer(const SequenceStore& rs, size_t r, uint32_t p, uint64_t& fwOut, uint64_t& canon) {
    uint64_t fw = 0, rv = 0;
    constexpr int shift = 2 * (kK - 1);
    for (int q = 0; q < kK; ++q) {
        const int c = rs.baseAt(r, p + static_cast<uint32_t>(q));
        if (c < 0) return false;
        fw = ((fw << 2) | static_cast<uint64_t>(c)) & kKMask;
        rv = (rv >> 2) | (static_cast<uint64_t>(3 - c) << shift);
    }
    fwOut = fw;
    canon = fw < rv ? fw : rv;
    return true;
}

double lengthWeightedMedian(std::vector<std::pair<double, size_t>> v) {
    if (v.empty()) return 0;
    std::sort(v.begin(), v.end());
    double tot = 0;
    for (const auto& x : v) tot += static_cast<double>(x.second);
    double acc = 0;
    for (const auto& x : v) {
        acc += static_cast<double>(x.second);
        if (acc >= tot / 2) return x.first;
    }
    return v.back().first;
}

std::string endName(const std::vector<std::string>& names, size_t rec, int end, size_t n) {
    return names[rec] + ":" + (end ? "R" : "L") + ":" + std::to_string(n);
}

}  // namespace

double medianOf(std::vector<uint32_t> v) {
    if (v.empty()) return 0;
    std::sort(v.begin(), v.end());
    const size_t n = v.size();
    if (n % 2) return v[n / 2];
    return (static_cast<double>(v[n / 2 - 1]) + static_cast<double>(v[n / 2])) / 2.0;
}

size_t tipFromCounts(const std::vector<uint32_t>& counts, double lda) {
    const double thr = std::max(2.0, 0.25 * lda);
    size_t tip = 0;
    for (uint32_t c : counts) {
        if (static_cast<double>(c) >= thr) break;
        ++tip;
    }
    return tip;
}

EndsAudit::EndsAudit(const EndsInput& in) : in_(in) {
    const std::vector<std::string>& s = *in_.seqs;
    for (size_t i = 0; i < s.size(); ++i) {
        if (s[i].size() < kEndsMinRecord) continue;
        rowRecords_.push_back(i);
        const bool circ = (*in_.names)[i].find("_circular") != std::string::npos;
        for (int e = 0; e < 2; ++e) {
            EndRow r;
            r.record = i;
            r.end = e ? 'R' : 'L';
            r.circular = circ;
            rows_.push_back(r);
        }
    }
    stats_.records = rowRecords_.size();
    stats_.ends = rows_.size();
}

void EndsAudit::addRawKeys(FixedCountTable& t) const {
    const std::vector<std::string>& s = *in_.seqs;
    for (size_t i : rowRecords_) {
        const std::string& q = s[i];
        const size_t L = q.size();
        const size_t w = std::min(L, kEndsDepthWindow);
        for (size_t p = 0; p + kK <= w; ++p) {
            const uint64_t c = canon31At(q, p);
            if (c != UINT64_MAX) t.addKey(c);
        }
        for (size_t p = L - w; p + kK <= L; ++p) {
            const uint64_t c = canon31At(q, p);
            if (c != UINT64_MAX) t.addKey(c);
        }
    }
}

void EndsAudit::computeStructure() {
    const auto t0 = std::chrono::steady_clock::now();
    const std::vector<std::string>& seqs = *in_.seqs;
    const size_t nrec = seqs.size();

    // ---- copy number --------------------------------------------------------------------------
    std::vector<std::pair<double, size_t>> longCov, rowCov, allCov;
    for (size_t i = 0; i < nrec; ++i) {
        const double c = i < in_.covs->size() ? (*in_.covs)[i] : 0;
        allCov.push_back({c, seqs[i].size()});
        if (seqs[i].size() >= kEndsMinRecord) rowCov.push_back({c, seqs[i].size()});
        if (seqs[i].size() >= 5000) longCov.push_back({c, seqs[i].size()});
    }
    const double gcov = !longCov.empty() ? lengthWeightedMedian(longCov)
                      : !rowCov.empty() ? lengthWeightedMedian(rowCov) : lengthWeightedMedian(allCov);
    stats_.genomeCov = gcov;
    std::vector<double> cn(nrec, 0);
    std::vector<char> anchor(nrec, 0);
    for (size_t i = 0; i < nrec; ++i) {
        const double c = i < in_.covs->size() ? (*in_.covs)[i] : 0;
        cn[i] = gcov > 0 ? c / gcov : 0;
        const size_t L = seqs[i].size();
        anchor[i] = (L >= 250 && cn[i] <= 1.6) || (L >= 5000 && cn[i] <= 4.0);
    }
    for (EndRow& r : rows_) {
        r.cn = cn[r.record];
        r.anchor = anchor[r.record] != 0;
    }

    // ---- non-unique tails (31-mer multiplicity over every record) ------------------------------
    {
        FixedCountTable mult;
        for (size_t i : rowRecords_) {
            const std::string& q = seqs[i];
            const size_t L = q.size(), w = std::min(L, kEndsSpan);
            for (size_t p = 0; p + kK <= w; ++p) { const uint64_t c = canon31At(q, p); if (c != UINT64_MAX) mult.addKey(c); }
            for (size_t p = L - w; p + kK <= L; ++p) { const uint64_t c = canon31At(q, p); if (c != UINT64_MAX) mult.addKey(c); }
        }
        mult.freeze();
        for (const std::string& q : seqs)
            forEachCanon31(q.data(), q.size(), [&](size_t, uint64_t c) { mult.hitSerial(c); });
        for (EndRow& r : rows_) {
            const std::string& q = seqs[r.record];
            const size_t L = q.size(), w = std::min(L, kEndsSpan);
            const size_t base = r.end == 'L' ? 0 : L - w;
            std::vector<char> masked(w, 0);
            for (size_t p = 0; p + kK <= w; ++p) {
                const uint64_t c = canon31At(q, base + p);
                if (c != UINT64_MAX && mult.count(c) > 1)
                    for (size_t x = p; x < p + kK; ++x) masked[x] = 1;
            }
            size_t t = 0;
            if (r.end == 'L') { while (t < w && masked[t]) ++t; }
            else { while (t < w && masked[w - 1 - t]) ++t; }
            r.nonuniqueTail = t;
        }
    }

    // ---- pair partners ---------------------------------------------------------------------------
    const SequenceStore* rs = in_.reads;
    const bool pairsUsable = rs && in_.insert.usable && rs->size() > 0;
    if (!pairsUsable) {
        for (EndRow& r : rows_) r.partnerClass = "na";
    } else {
        const int32_t W = in_.insert.maxPlausible > 0 ? in_.insert.maxPlausible
                                                      : static_cast<int32_t>(in_.insert.mean * 3);
        // index the end windows of every record; a k-mer seen twice anywhere in contigs.fasta is ambiguous
        std::unordered_map<uint64_t, uint64_t> index;
        std::vector<std::pair<size_t, size_t>> win(nrec, {0, 0});    // [0, a) and [b, L)
        // the window must hold the reads that reach the last 400 bp (locus pools) and the reads that point
        // out of the end within the insert reach (partners)
        const size_t winLen = std::max(static_cast<size_t>(std::max(W, 0)),
                                       kEndsPoolFlank + static_cast<size_t>(rs->maxReadLength()));
        for (size_t i = 0; i < nrec; ++i) {
            const size_t L = seqs[i].size();
            const size_t w = std::min(L, winLen);
            win[i] = {w, L - w};
        }
        for (uint32_t i = 0; i < nrec; ++i) {
            const std::string& q = seqs[i];
            const size_t L = q.size();
            if (L < static_cast<size_t>(kK)) continue;
            uint64_t fw = 0, rv = 0;
            int valid = 0;
            constexpr int shift = 2 * (kK - 1);
            for (size_t p = 0; p < L; ++p) {
                const int c = baseCode2(q[p]);
                if (c < 0) { valid = 0; fw = rv = 0; continue; }
                fw = ((fw << 2) | static_cast<uint64_t>(c)) & kKMask;
                rv = (rv >> 2) | (static_cast<uint64_t>(3 - c) << shift);
                if (++valid < kK) continue;
                const size_t start = p + 1 - kK;
                const uint64_t canon = fw < rv ? fw : rv;
                const bool inWin = start + kK <= win[i].first || start >= win[i].second;
                if (inWin) {
                    auto it = index.find(canon);
                    if (it == index.end())
                        index.emplace(canon, packPlace(i, static_cast<uint32_t>(start), fw == canon ? 0 : 1));
                    else
                        it->second = kAmbig;
                }
            }
        }
        // a window k-mer that also occurs outside every window
        for (uint32_t i = 0; i < nrec; ++i) {
            const std::string& q = seqs[i];
            forEachCanon31(q.data(), q.size(), [&](size_t start, uint64_t canon) {
                const bool inWin = start + kK <= win[i].first || start >= win[i].second;
                if (inWin) return;
                auto it = index.find(canon);
                if (it != index.end()) it->second = kAmbig;
            });
        }

        const size_t nEnds = nrec * 2;
        struct Acc {
            std::vector<size_t> outward, unplaced, otherEnd, selfWrap;
            std::map<std::pair<uint32_t, uint32_t>, size_t> links;   // (port, port) -> pairs
            std::vector<std::pair<uint64_t, uint32_t>> pool;         // (read, row)
            size_t placed = 0, pairs = 0;
        };
        std::vector<int32_t> rowOf(nEnds, -1);
        for (size_t ri = 0; ri < rows_.size(); ++ri)
            rowOf[rows_[ri].record * 2 + (rows_[ri].end == 'R' ? 1 : 0)] = static_cast<int32_t>(ri);
        // a placed read overlapping the last 400 bp of an end with a row belongs to that end's pool
        auto addOwn = [&](Acc& a, const Placed& x, uint64_t rx) {
            const int64_t L = static_cast<int64_t>(seqs[x.rec].size());
            const int64_t s0 = x.pos, e0 = static_cast<int64_t>(x.pos) + static_cast<int64_t>(x.len);
            const int32_t rr = rowOf[x.rec * 2 + 1], rl = rowOf[x.rec * 2];
            if (rr >= 0 && e0 > L - static_cast<int64_t>(kEndsPoolFlank)) a.pool.push_back({rx, static_cast<uint32_t>(rr)});
            if (rl >= 0 && s0 < static_cast<int64_t>(kEndsPoolFlank)) a.pool.push_back({rx, static_cast<uint32_t>(rl)});
        };
        const int nt = std::max(1, in_.threads);
        std::vector<Acc> acc(static_cast<size_t>(nt));
        for (Acc& a : acc) {
            a.outward.assign(nEnds, 0);
            a.unplaced.assign(nEnds, 0);
            a.otherEnd.assign(nEnds, 0);
            a.selfWrap.assign(nEnds, 0);
        }
        auto place = [&](size_t r) -> Placed {
            Placed out;
            const uint32_t len = rs->length(r);
            if (len < static_cast<uint32_t>(kK)) return out;
            struct Vote { uint32_t rec; int32_t pos; uint8_t orient; int count; };
            Vote votes[kProbes];
            int distinct = 0;
            const uint32_t probes = std::min<uint32_t>(kProbes, len - kK + 1);
            const uint32_t step = std::max<uint32_t>(1, (len - kK + 1) / probes);
            for (uint32_t rp = 0; rp + kK <= len; rp += step) {
                uint64_t fw = 0, canon = 0;
                if (!readKmer(*rs, r, rp, fw, canon)) continue;
                auto it = index.find(canon);
                if (it == index.end() || it->second == kAmbig) continue;
                const uint32_t c = static_cast<uint32_t>(it->second >> 33);
                const int up = static_cast<int>((it->second >> 1) & 0xFFFFFFFFULL);
                const int sflag = static_cast<int>(it->second & 1);
                const int orient = ((fw == canon) ? 0 : 1) ^ sflag;
                const int32_t startPos = orient == 0 ? static_cast<int32_t>(up - static_cast<int>(rp))
                                                     : static_cast<int32_t>(up - static_cast<int>(len - rp - kK));
                int found = -1;
                for (int q = 0; q < distinct; ++q)
                    if (votes[q].rec == c && votes[q].pos == startPos && votes[q].orient == orient) { found = q; break; }
                if (found >= 0) ++votes[found].count;
                else if (distinct < kProbes) votes[distinct++] = {c, startPos, static_cast<uint8_t>(orient), 1};
            }
            int best = -1, bestVotes = 0;
            for (int q = 0; q < distinct; ++q)
                if (votes[q].count > bestVotes) { bestVotes = votes[q].count; best = q; }
            if (best < 0 || bestVotes < kMinVotes) return out;
            out.rec = votes[best].rec;
            out.pos = votes[best].pos;
            out.orient = votes[best].orient;
            out.len = len;
            return out;
        };
        // distance from a placed read to the end it points at, and that end (0 left, 1 right)
        auto pointing = [&](const Placed& x, int& end) -> int32_t {
            const int32_t L = static_cast<int32_t>(seqs[x.rec].size());
            end = x.orient == 0 ? 1 : 0;
            return x.orient == 0 ? L - x.pos : x.pos + static_cast<int32_t>(x.len);
        };
        const size_t nreads = rs->size();
        const size_t paired = rs->pairedReads();
        std::atomic<size_t> next{0};
        constexpr size_t kBlock = 8192;   // even: blocks never split a pair
        auto worker = [&](Acc& a) {
            for (;;) {
                const size_t begin = next.fetch_add(kBlock);
                if (begin >= nreads) break;
                const size_t end = std::min(begin + kBlock, nreads);
                for (size_t r = begin; r < end;) {
                    // pairs are (r, r + 1) for even r below pairedReads(); a block starts even
                    const bool isPair = r < paired && (r % 2 == 0) && r + 1 < end;
                    if (isPair) {
                        const Placed pa = place(r), pb = place(r + 1);
                        a.placed += pa.ok() + pb.ok();
                        ++a.pairs;
                        const Placed* pr[2] = {&pa, &pb};
                        for (int s = 0; s < 2; ++s) {
                            const Placed& x = *pr[s];
                            const Placed& y = *pr[1 - s];
                            if (!x.ok()) continue;
                            addOwn(a, x, r + static_cast<size_t>(s));
                            int xe = 0;
                            const int32_t dx = pointing(x, xe);
                            if (dx < 0 || dx > W) continue;
                            const uint32_t E = x.rec * 2 + static_cast<uint32_t>(xe);
                            ++a.outward[E];
                            if (rowOf[E] >= 0) a.pool.push_back({r + static_cast<size_t>(1 - s), static_cast<uint32_t>(rowOf[E])});
                            if (!y.ok()) { ++a.unplaced[E]; continue; }
                            int ye = 0;
                            const int32_t dy = pointing(y, ye);
                            if (y.rec == x.rec) {
                                // a pair leaving the record through opposite ends, facing away (circle)
                                if (ye != xe && y.orient != x.orient && dy >= 0 && dy <= W) {
                                    const Placed& fwd = x.orient == 0 ? x : y;
                                    const Placed& rev = x.orient == 0 ? y : x;
                                    if (fwd.pos > rev.pos) ++a.selfWrap[E];
                                }
                                continue;
                            }
                            ++a.otherEnd[E];
                            if (dy >= 0 && dy <= W) ++a.links[{E, y.rec * 2 + static_cast<uint32_t>(ye)}];
                        }
                        r += 2;
                    } else {
                        const Placed px = place(r);
                        if (px.ok()) {
                            ++a.placed;
                            addOwn(a, px, r);
                            int xe = 0;
                            const int32_t dx = pointing(px, xe);
                            if (dx >= 0 && dx <= W) ++a.outward[px.rec * 2 + static_cast<uint32_t>(xe)];
                        }
                        r += 1;
                    }
                }
            }
        };
        {
            std::vector<std::thread> pool;
            for (int t = 0; t < nt; ++t) pool.emplace_back(worker, std::ref(acc[static_cast<size_t>(t)]));
            for (std::thread& t : pool) t.join();
        }
        Acc tot;
        tot.outward.assign(nEnds, 0);
        tot.unplaced.assign(nEnds, 0);
        tot.otherEnd.assign(nEnds, 0);
        tot.selfWrap.assign(nEnds, 0);
        for (const Acc& a : acc) {
            for (size_t e = 0; e < nEnds; ++e) {
                tot.outward[e] += a.outward[e];
                tot.unplaced[e] += a.unplaced[e];
                tot.otherEnd[e] += a.otherEnd[e];
                tot.selfWrap[e] += a.selfWrap[e];
            }
            for (const auto& kv : a.links) tot.links[kv.first] += kv.second;
            tot.placed += a.placed;
            tot.pairs += a.pairs;
            readRows_.insert(readRows_.end(), a.pool.begin(), a.pool.end());
        }
        std::sort(readRows_.begin(), readRows_.end());
        readRows_.erase(std::unique(readRows_.begin(), readRows_.end()), readRows_.end());
        rowPool_.assign(rows_.size(), std::vector<uint64_t>());
        poolBits_.assign((nreads + 63) / 64, 0);
        for (const auto& x : readRows_) {
            rowPool_[x.second].push_back(x.first);
            poolBits_[x.first >> 6] |= 1ULL << (x.first & 63);
        }
        stats_.readsPlaced = tot.placed;
        stats_.pairsUsed = tot.pairs;
        // partners per end
        std::vector<std::vector<std::pair<size_t, uint32_t>>> part(nEnds);   // (links, port)
        for (const auto& kv : tot.links) part[kv.first.first].push_back({kv.second, kv.first.second});
        for (EndRow& r : rows_) {
            const uint32_t E = static_cast<uint32_t>(r.record * 2 + (r.end == 'R' ? 1 : 0));
            r.outward = tot.outward[E];
            r.matesUnplaced = tot.unplaced[E];
            r.matesOtherEnd = tot.otherEnd[E];
            r.selfWrap = tot.selfWrap[E];
            std::vector<std::pair<size_t, uint32_t>> p = part[E];
            std::sort(p.begin(), p.end(), [&](const auto& a, const auto& b) {
                if (a.first != b.first) return a.first > b.first;
                return a.second < b.second;
            });
            std::string top;
            for (size_t j = 0; j < p.size() && j < 5; ++j) {
                if (j) top += ";";
                top += endName(*in_.names, p[j].second / 2, static_cast<int>(p[j].second & 1), p[j].first);
            }
            r.partners = top;
            std::vector<std::pair<size_t, uint32_t>> good, any3;
            for (const auto& x : p) {
                if (x.first < kEndsMinLinks) continue;
                any3.push_back(x);
                if (anchor[x.second / 2]) good.push_back(x);
            }
            r.nPartners = any3.size();
            if (good.empty()) {
                r.partnerClass = any3.empty() ? "none" : "repeat_only";
                if (!any3.empty())
                    r.bestPartner = endName(*in_.names, any3[0].second / 2, static_cast<int>(any3[0].second & 1), any3[0].first);
            } else {
                r.partnerClass = (good.size() >= 2 && good[1].first * 3 >= good[0].first) ? "branching" : "unique";
                r.bestPartner = endName(*in_.names, good[0].second / 2, static_cast<int>(good[0].second & 1), good[0].first);
            }
        }
    }
    // the ends audited for tips (sp2: anchor records with >= 10 pool reads) and their k-mer tables
    rowPool_.resize(rows_.size());
    audit_.assign(rows_.size(), 0);
    local_.clear();
    local_.resize(rows_.size());
    for (size_t ri = 0; ri < rows_.size(); ++ri) {
        EndRow& r = rows_[ri];
        r.pool = rowPool_[ri].size();
        if (!r.anchor) { r.tipStatus = "not_anchor"; continue; }
        if (r.pool < kEndsMinPool) { r.tipStatus = "small_pool"; continue; }
        r.tipStatus = "audited";
        audit_[ri] = 1;
        const std::string& q = seqs[r.record];
        const size_t L = q.size(), w = std::min(L, kEndsDepthWindow);
        const size_t from = r.end == 'L' ? 0 : L - w;
        for (size_t p = from; p + kK <= from + w; ++p) {
            const uint64_t c = canon31At(q, p);
            if (c != UINT64_MAX) local_[ri].addKey(c);
        }
        local_[ri].freeze(false);
    }
    for (const EndRow& r : rows_) {
        if (r.partnerClass == "unique") ++stats_.unique;
        else if (r.partnerClass == "branching") ++stats_.branching;
        else if (r.partnerClass == "repeat_only") ++stats_.repeatOnly;
        else if (r.partnerClass == "none") ++stats_.none;
        else ++stats_.na;
    }
    stats_.seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

void EndsAudit::countPoolRead(uint64_t readIndex, const char* seq, size_t len) {
    auto it = std::lower_bound(readRows_.begin(), readRows_.end(), std::make_pair(readIndex, uint32_t(0)));
    for (; it != readRows_.end() && it->first == readIndex; ++it) {
        if (!audit_[it->second]) continue;
        FixedCountTable& t = local_[it->second];
        forEachCanon31(seq, len, [&](size_t, uint64_t c) { t.hit(c); });
    }
}

void EndsAudit::countPoolsFromStore() {
    for (FixedCountTable& t : local_) if (t.size()) t.resetCounts();
    if (!in_.reads) return;
    std::string s;
    uint64_t last = UINT64_MAX;
    for (const auto& x : readRows_) {
        if (!audit_[x.second]) continue;
        if (x.first != last) { in_.reads->decode(x.first, s); last = x.first; }
        FixedCountTable& t = local_[x.second];
        forEachCanon31(s.data(), s.size(), [&](size_t, uint64_t c) { t.hit(c); });
    }
}

void EndsAudit::finishTips(const FixedCountTable& raw, const std::string& poolSource) {
    const auto t0 = std::chrono::steady_clock::now();
    const std::vector<std::string>& seqs = *in_.seqs;
    stats_.poolSource = poolSource;
    stats_.audited = 0;
    stats_.poolReads = 0;
    for (size_t ri = 0; ri < rows_.size(); ++ri) {
        EndRow& r = rows_[ri];
        const std::string& q = seqs[r.record];
        const size_t L = q.size();
        const size_t w = std::min(L, kEndsDepthWindow);
        // k-mers of the last 300 bp (window) and from the end inward (counts[j]: k-mer ending j bases
        // before the end), in the end's own coordinates
        std::vector<uint64_t> window, inward;
        const size_t jmax = std::min(kEndsTipScan, L - kK);
        if (r.end == 'R') {
            for (size_t p = L - w; p + kK <= L; ++p) window.push_back(canon31At(q, p));
            for (size_t j = 0; j <= jmax; ++j) inward.push_back(canon31At(q, L - j - kK));
        } else {
            for (size_t p = 0; p + kK <= w; ++p) window.push_back(canon31At(q, p));
            for (size_t j = 0; j <= jmax; ++j) inward.push_back(canon31At(q, j));
        }
        std::vector<uint32_t> global;
        for (uint64_t c : window) global.push_back(c == UINT64_MAX ? 0 : raw.count(c));
        r.endKdepth = medianOf(global);
        r.tip = 0;
        stats_.poolReads += r.pool;
        if (!audit_.empty() && audit_[ri]) {
            ++stats_.audited;
            const FixedCountTable& t = local_[ri];
            std::vector<uint32_t> wl, il;
            for (uint64_t c : window) wl.push_back(c == UINT64_MAX ? 0 : t.count(c));
            for (uint64_t c : inward) il.push_back(c == UINT64_MAX ? 0 : t.count(c));
            r.poolKdepth = medianOf(wl);
            r.tip = tipFromCounts(il, r.poolKdepth);
        }
    }
    // genome end depth: median end depth over the ends of records >= 5 kb (else over every end)
    std::vector<double> d5, dall;
    for (const EndRow& r : rows_) {
        dall.push_back(r.endKdepth);
        if (seqs[r.record].size() >= 5000) d5.push_back(r.endKdepth);
    }
    std::vector<double>& use = d5.empty() ? dall : d5;
    double g = 0;
    if (!use.empty()) {
        std::sort(use.begin(), use.end());
        const size_t n = use.size();
        g = n % 2 ? use[n / 2] : (use[n / 2 - 1] + use[n / 2]) / 2.0;
    }
    stats_.genomeEndKdepth = g;
    stats_.tipEnds = stats_.tipBases = 0;
    for (EndRow& r : rows_) {
        r.endCn = g > 0 ? r.endKdepth / g : 0;
        if (r.tip) { ++stats_.tipEnds; stats_.tipBases += r.tip; }
    }
    stats_.seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

std::vector<std::string> lowercaseTips(const std::vector<std::string>& seqs, const std::vector<std::string>& names,
                                       std::vector<EndRow>& rows, EndsStats& st, std::vector<std::string>& edits) {
    std::vector<std::string> out = seqs;
    st.lowercasedBases = 0;
    for (EndRow& r : rows) {
        r.lowercased = 0;
        if (!r.tip) continue;
        std::string& q = out[r.record];
        const size_t L = q.size();
        const size_t n = std::min(r.tip, L);
        const size_t from = r.end == 'L' ? 0 : L - n;
        for (size_t p = from; p < from + n; ++p)
            q[p] = static_cast<char>(std::tolower(static_cast<unsigned char>(q[p])));
        r.lowercased = n;
        st.lowercasedBases += n;
        edits.push_back(names[r.record] + "\ttips_lowercase\tlowercase_" + (r.end == 'L' ? "left" : "right") +
                        "\t" + std::to_string(from) + "\t" + std::to_string(from + n) + "\t" + std::to_string(n));
    }
    return out;
}

std::string endsTsv(const std::vector<EndRow>& rows, const std::vector<std::string>& names,
                    const std::vector<std::string>& seqs) {
    std::string o =
        "record\tend\tlength\ttopology\tcn\tanchor\tend_kdepth\tend_cn\tnonunique_tail_bp\toutward_reads\t"
        "mates_unplaced\tmates_other_end\tself_wrap_pairs\tpartner_class\tn_partners\tbest_partner\tpartners\t"
        "pool_reads\tpool_kdepth\ttip_status\tunsupported_tip_bp\tlowercased_bp\n";
    char b[512];
    for (const EndRow& r : rows) {
        std::snprintf(b, sizeof b, "\t%c\t%zu\t%s\t%.3f\t%d\t%.1f\t%.3f\t%zu\t%zu\t%zu\t%zu\t%zu\t%s\t%zu\t",
                      r.end, seqs[r.record].size(), r.circular ? "circular" : "linear", r.cn, r.anchor ? 1 : 0,
                      r.endKdepth, r.endCn, r.nonuniqueTail, r.outward, r.matesUnplaced, r.matesOtherEnd,
                      r.selfWrap, r.partnerClass.c_str(), r.nPartners);
        o += names[r.record];
        o += b;
        o += r.bestPartner.empty() ? "-" : r.bestPartner;
        o += "\t";
        o += r.partners.empty() ? "-" : r.partners;
        std::snprintf(b, sizeof b, "\t%zu\t%.1f\t%s\t%zu\t%zu\n", r.pool, r.poolKdepth, r.tipStatus.c_str(), r.tip,
                      r.lowercased);
        o += b;
    }
    return o;
}

std::string endsJson(const EndsStats& st, bool lowercase) {
    JsonObj o;
    o.str("file", "ends.tsv");
    o.u("min_record_length", kEndsMinRecord);
    o.u("records", st.records).u("ends", st.ends);
    o.u("tip_ends", st.tipEnds).u("tip_bases", st.tipBases);
    o.b("lowercase", lowercase);
    o.u("lowercased_bases", st.lowercasedBases);
    JsonObj pc;
    pc.u("unique", st.unique).u("branching", st.branching).u("repeat_only", st.repeatOnly).u("none", st.none).u("na", st.na);
    o.raw("partner_class", pc.done());
    o.num("genome_cov", st.genomeCov, 4);
    o.num("genome_end_kdepth", st.genomeEndKdepth, 1);
    o.u("reads_placed", st.readsPlaced).u("pairs_examined", st.pairsUsed);
    o.u("ends_audited", st.audited).u("pool_reads_total", st.poolReads);
    o.str("pool_source", st.poolSource);
    o.str("tip_rule", "anchor ends with >= 10 locus-pool reads: terminal bases before the first 31-mer with pool "
                      "count >= max(2, 0.25 x the pool's median count over the end's last 300 bp), at most 201 "
                      "(sp2 end audit)");
    o.str("lowercase_class", "unsupported_tip: ingested as masked bases (EVAL_PLAN_P2 s3.3 item 3)");
    o.num("seconds", st.seconds, 3);
    return o.done();
}

}  // namespace p2
}  // namespace ts
