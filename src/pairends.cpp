#include "pairends.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <mutex>
#include <thread>

#include "kmer.h"

namespace ts {

namespace {

// Same anchoring constants as the polisher and the resolver. They agree because they are
// solving the same problem -- put a read where it belongs, or nowhere -- and a read placed
// by one rule and rejected by another would make the stages disagree about the assembly.
constexpr int kMaxProbes = 12;
constexpr int kMinVotes = 2;
constexpr uint64_t kAmbiguous = UINT64_MAX;

inline uint64_t packIndex(uint32_t contig, uint32_t pos, int strand) {
    return (static_cast<uint64_t>(contig) << 33) | (static_cast<uint64_t>(pos) << 1) |
           static_cast<uint64_t>(strand & 1);
}

inline uint64_t pairKey(uint32_t a, uint32_t b) {
    return (static_cast<uint64_t>(a) << 32) | static_cast<uint64_t>(b);
}

struct Placed {
    uint32_t contig = UINT32_MAX;
    int32_t pos = 0;      // start of the read on the contig's forward strand
    uint8_t orient = 0;   // 0 = read runs along the contig forward
    uint32_t len = 0;
    bool ok() const { return contig != UINT32_MAX; }
};

}  // namespace

ContigEndLinks::ContigEndLinks(const std::vector<std::string>& contigs,
                               const SequenceStore& reads, const InsertModel& insert,
                               int threads, int anchorK) {
    insert_ = insert;
    contigLen_.resize(contigs.size());
    for (size_t i = 0; i < contigs.size(); ++i) contigLen_[i] = contigs[i].size();

    // Without a fragment distribution there is no way to turn a pair into a distance, and
    // a link with no distance cannot confirm or deny anything. Better to be unusable and
    // say so than to invent a scale.
    if (!insert.usable || contigs.empty() || reads.pairCount() == 0) return;
    if (threads <= 0) threads = 1;
    const int k = std::max(15, std::min(anchorK, 31));

    // ---- index the contigs -------------------------------------------------
    // Uniqueness is judged over ALL contigs, not just their end windows: a k-mer unique
    // among the ends but repeated in some contig's interior would place reads at an end
    // they do not come from, which is exactly the mistake that would manufacture a link.
    std::unordered_map<Kmer, uint64_t, KmerHasher> index;
    size_t totalBases = 0;
    for (const std::string& s : contigs) totalBases += s.size();
    index.reserve(totalBases * 2);
    for (uint32_t c = 0; c < contigs.size(); ++c) {
        const std::string& s = contigs[c];
        if (s.size() < static_cast<size_t>(k)) continue;
        Kmer fwd = 0, rc = 0;
        int valid = 0;
        for (uint32_t p = 0; p < s.size(); ++p) {
            const int b = baseCode(s[p]);
            if (b < 0) { valid = 0; continue; }
            fwd = pushBack(fwd, b, k);
            rc = pushFrontRc(rc, b, k);
            if (++valid < k) continue;
            const Kmer canon = fwd < rc ? fwd : rc;
            auto it = index.find(canon);
            if (it == index.end()) {
                index.emplace(canon, packIndex(c, p + 1 - static_cast<uint32_t>(k),
                                               (fwd == canon) ? 0 : 1));
            } else {
                it->second = kAmbiguous;
            }
        }
    }

    // ---- anchor every read -------------------------------------------------
    std::vector<Placed> placed(reads.size());
    const size_t nreads = reads.size();
    std::atomic<size_t> next{0};
    std::atomic<size_t> anchored{0};
    const size_t block = 4096;

    auto worker = [&]() {
        struct Vote { uint32_t contig; int32_t pos; uint8_t orient; int count; };
        Vote votes[kMaxProbes];
        while (true) {
            const size_t begin = next.fetch_add(block);
            if (begin >= nreads) break;
            const size_t end = std::min(begin + block, nreads);
            for (size_t r = begin; r < end; ++r) {
                const uint32_t len = reads.length(r);
                if (len < static_cast<uint32_t>(k)) continue;
                int distinct = 0;
                const uint32_t probes = std::min<uint32_t>(kMaxProbes, len - k + 1);
                const uint32_t step = std::max<uint32_t>(1, (len - k + 1) / probes);
                for (uint32_t rp = 0; rp + k <= len; rp += step) {
                    Kmer fwd = 0, rcm = 0;
                    bool bad = false;
                    for (int q = 0; q < k; ++q) {
                        const int b = reads.baseAt(r, rp + static_cast<uint32_t>(q));
                        if (b < 0) { bad = true; break; }
                        fwd = pushBack(fwd, b, k);
                        rcm = pushFrontRc(rcm, b, k);
                    }
                    if (bad) continue;
                    const Kmer canon = fwd < rcm ? fwd : rcm;
                    auto it = index.find(canon);
                    if (it == index.end() || it->second == kAmbiguous) continue;

                    const uint32_t c = static_cast<uint32_t>(it->second >> 33);
                    const int up = static_cast<int>((it->second >> 1) & 0xFFFFFFFFULL);
                    const int sflag = static_cast<int>(it->second & 1);
                    const int orient = ((fwd == canon) ? 0 : 1) ^ sflag;
                    const int32_t startPos =
                        (orient == 0) ? static_cast<int32_t>(up - static_cast<int>(rp))
                                      : static_cast<int32_t>(
                                            up - static_cast<int>(len - rp - static_cast<uint32_t>(k)));
                    int found = -1;
                    for (int q = 0; q < distinct; ++q) {
                        if (votes[q].contig == c && votes[q].pos == startPos &&
                            votes[q].orient == static_cast<uint8_t>(orient)) { found = q; break; }
                    }
                    if (found >= 0) ++votes[found].count;
                    else if (distinct < kMaxProbes) {
                        votes[distinct++] = {c, startPos, static_cast<uint8_t>(orient), 1};
                    }
                }
                int bestIdx = -1, bestVotes = 0;
                for (int q = 0; q < distinct; ++q) {
                    if (votes[q].count > bestVotes) { bestVotes = votes[q].count; bestIdx = q; }
                }
                if (bestIdx < 0 || bestVotes < kMinVotes) continue;
                placed[r] = Placed{votes[bestIdx].contig, votes[bestIdx].pos,
                                   votes[bestIdx].orient, len};
                anchored.fetch_add(1, std::memory_order_relaxed);
            }
        }
    };

    {
        std::vector<std::thread> pool;
        pool.reserve(static_cast<size_t>(threads));
        for (int t = 0; t < threads; ++t) pool.emplace_back(worker);
        for (std::thread& t : pool) t.join();
    }

    if (totalBases > 0) {
        // Mean depth over anchored reads: the floor below scales with it so a 20x and a
        // 200x library are not held to the same absolute count of supporting pairs.
        double bases = 0;
        for (size_t r = 0; r < nreads; ++r) if (placed[r].ok()) bases += placed[r].len;
        meanDepth_ = bases / static_cast<double>(totalBases);
    }

    // ---- turn placed pairs into end links -----------------------------------
    // Only pairs whose two reads point OUTWARD through different contig ends say anything
    // about what follows those ends. A pair pointing inward is inside one contig and
    // describes nothing beyond it.
    const int32_t reach = insert_.maxPlausible > 0 ? insert_.maxPlausible
                                                   : static_cast<int32_t>(insert_.mean * 3);
    for (size_t r = 0; r < nreads; r += 2) {
        if (!reads.hasMate(r)) break;
        const Placed& a = placed[r];
        const Placed& b = placed[r + 1];
        if (!a.ok() || !b.ok()) continue;

        if (a.contig != b.contig) {
            ++contigLinks_[pairKey(std::min(a.contig, b.contig), std::max(a.contig, b.contig))];
        }

        // Distance from each read to the end it points at. orient 0 means the read runs
        // along the contig forward, so it points at the RIGHT end.
        const int32_t la = static_cast<int32_t>(contigLen_[a.contig]);
        const int32_t lb = static_cast<int32_t>(contigLen_[b.contig]);
        const int aEnd = a.orient == 0 ? 1 : 0;
        const int bEnd = b.orient == 0 ? 1 : 0;

        // Same contig needs more care than "opposite ends", and getting it wrong makes the
        // circularity call meaningless: an ordinary inward-facing pair in the middle of a
        // contig ALSO has its forward read exiting the right end and its reverse read
        // exiting the left, so an end-only test counts every normal fragment as evidence
        // of a circle. Measured before this check existed: a median of 50 contigs per
        // isolate were called circular, against the one to six replicons a Klebsiella
        // isolate actually has.
        //
        // What separates them is position, not orientation. In a normal pair the two reads
        // face each other, so the forward-strand read lies UPSTREAM of the reverse-strand
        // one. In a pair spanning a circular junction they face away from each other, and
        // the forward read lies downstream of the reverse read -- the fragment leaves the
        // right end and re-enters at the left.
        if (a.contig == b.contig) {
            if (aEnd == bEnd) continue;                 // both off one end: no molecule
            if (a.orient == b.orient) continue;         // same strand: not a proper pair
            const Placed& fwd = a.orient == 0 ? a : b;
            const Placed& rev = a.orient == 0 ? b : a;
            if (fwd.pos <= rev.pos) continue;           // ordinary inward pair, says nothing
        }

        const int32_t da = aEnd == 1 ? la - a.pos : a.pos + static_cast<int32_t>(a.len);
        const int32_t db = bEnd == 1 ? lb - b.pos : b.pos + static_cast<int32_t>(b.len);
        if (da < 0 || db < 0 || da > reach || db > reach) continue;

        // Span consumed inside the two contigs; whatever lies between them is the gap.
        const int32_t span = da + db;
        ends_[pairKey(endPort(a.contig, aEnd), endPort(b.contig, bEnd))].spans.push_back(span);
        ends_[pairKey(endPort(b.contig, bEnd), endPort(a.contig, aEnd))].spans.push_back(span);
    }

    usable_ = true;
}

size_t ContigEndLinks::supportFloor() const {
    // Same shape as the scaffolder's: at least 3, never more than 10, scaling with depth.
    const double scaled = meanDepth_ * 0.06;
    const double f = std::min(10.0, std::max(3.0, scaled));
    return static_cast<size_t>(f);
}

PairEvidence ContigEndLinks::query(uint32_t exitPort, uint32_t entryPort,
                                   int32_t predictedGap) const {
    PairEvidence ev;
    if (!usable_) return ev;

    // Could a fragment have spanned this at all? If not, the absence of pairs is a fact
    // about the library and not about the join, and the caller must not read it as a veto.
    const int32_t maxSpan = insert_.maxPlausible > 0 ? insert_.maxPlausible
                                                     : static_cast<int32_t>(insert_.mean * 3);
    ev.spannable = predictedGap < maxSpan;

    auto it = ends_.find(pairKey(exitPort, entryPort));
    if (it == ends_.end()) return ev;
    const std::vector<int32_t>& spans = it->second.spans;
    ev.total = spans.size();
    if (spans.empty()) return ev;

    // Each pair implies a gap of mean - span. Tolerance is the fragment spread itself:
    // a pair cannot locate a junction more precisely than the distribution it came from.
    std::vector<int32_t> gaps;
    gaps.reserve(spans.size());
    for (int32_t s : spans) gaps.push_back(static_cast<int32_t>(insert_.mean) - s);
    std::sort(gaps.begin(), gaps.end());
    ev.medianGap = gaps[gaps.size() / 2];

    const int32_t tol = static_cast<int32_t>(std::max(200.0, 3.0 * insert_.stddev));
    for (int32_t g : gaps) {
        if (std::abs(g - predictedGap) <= tol) ++ev.consistent;
    }
    return ev;
}

PairEvidence ContigEndLinks::circular(uint32_t contig) const {
    PairEvidence ev;
    if (!usable_ || contig >= contigLen_.size()) return ev;
    // Below one fragment length the question cannot be asked. A single ordinary fragment
    // then reaches from one end of the contig to the other, so its two reads sit at
    // opposite ends on opposite strands whatever the topology, and the position test that
    // separates a junction-spanning pair from an inward-facing one has nothing left to
    // separate. The three contigs that survived the position fix on KP1000 were 336, 388
    // and 489 bp -- all of them well inside a fragment, none of them a plausible replicon.
    const int64_t minLen = static_cast<int64_t>(insert_.mean + 3.0 * insert_.stddev);
    if (static_cast<int64_t>(contigLen_[contig]) < minLen) return ev;
    // A circle joins its own right end to its own left end with no sequence between them,
    // so the pairs that prove it are the ones leaving the contig through opposite ends.
    // Those are recorded like any other end link; only the predicted gap differs, being
    // zero rather than something the panel proposed.
    auto it = ends_.find(pairKey(endPort(contig, 1), endPort(contig, 0)));
    if (it == ends_.end()) return ev;
    ev.spannable = true;
    ev.total = it->second.spans.size();
    std::vector<int32_t> gaps;
    gaps.reserve(ev.total);
    for (int32_t s : it->second.spans) gaps.push_back(static_cast<int32_t>(insert_.mean) - s);
    if (gaps.empty()) return ev;
    std::sort(gaps.begin(), gaps.end());
    ev.medianGap = gaps[gaps.size() / 2];
    const int32_t tol = static_cast<int32_t>(std::max(200.0, 3.0 * insert_.stddev));
    for (int32_t g : gaps) {
        if (std::abs(g) <= tol) ++ev.consistent;
    }
    return ev;
}

size_t ContigEndLinks::linkWeight(uint32_t a, uint32_t b) const {
    auto it = contigLinks_.find(pairKey(std::min(a, b), std::max(a, b)));
    return it == contigLinks_.end() ? 0 : it->second;
}

}  // namespace ts
