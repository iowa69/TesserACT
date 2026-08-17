#include "organism_layout.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <unordered_map>

#include "graph.h"

namespace ts {

namespace {

// A contig needs this many markers agreeing on one place before it may be placed.
// Two markers can agree by chance in a genome carrying thousands; four cannot.
constexpr size_t kMinMarkersToPlace = 4;

// Share of a contig's markers that must agree on one location before it may be placed.
//
// Held low deliberately. By the time this stage runs the contigs already carry N gaps that
// the adjacency join inserted, and those gaps exist in this assembly and not in the panel
// genome the track came from. Every gap shifts the start position that each downstream
// marker implies, so a long contig's votes arrive in several clusters separated by the
// cumulative gap length -- about 49 kb on a typical isolate here. Demanding that 60% of
// votes fall in ONE cluster rejected the largest contig in the assembly for a reason that
// has nothing to do with whether it belongs where its markers say.
//
// A third of the markers agreeing on one place, out of thousands, is still far beyond
// coincidence, and the tolerance below absorbs the rest.
constexpr double kMinCoherent = 0.35;

// Two placements overlapping by more than this share of the shorter one cannot both be
// right; the one resting on more markers is kept.
constexpr double kMaxPlacementOverlap = 0.5;

// Beyond this distance on the track, two neighbouring contigs are not asserted to be one
// scaffold. The template is a different genome and at some separation it stops being
// evidence about this one.
constexpr int64_t kMaxTrackGap = 400000;

// Longest N run written for a gap.
//
// The track's distance is what the *template* has between these two loci, which is not
// the same as what this isolate is missing -- an isolate whose sequence is already fully
// placed and still shows a large track gap is one whose relative carries an insertion it
// lacks. Writing the track distance as Ns would assert absent sequence that is not
// absent. Capping asserts the order and orientation, which the track does support, and
// declines to assert a quantity, which it does not.
constexpr int64_t kGapCap = 2000;

constexpr size_t kMinOverlapMerge = 50;
constexpr int64_t kOverlapSlack = 600;

// Largest exact suffix/prefix overlap near `predicted`; 0 if none. Same contract as the
// join stage: the template proposes, the sequence confirms, and an overlap the sequence
// denies is never trimmed.
size_t exactOverlap(const std::string& left, const std::string& right, int64_t predicted) {
    if (predicted <= 0) return 0;
    const size_t hi = std::min<size_t>({static_cast<size_t>(predicted + kOverlapSlack),
                                        left.size(), right.size()});
    if (hi < kMinOverlapMerge) return 0;
    const size_t lo = std::max<size_t>(
        kMinOverlapMerge, predicted > kOverlapSlack
                              ? static_cast<size_t>(predicted - kOverlapSlack)
                              : kMinOverlapMerge);
    const std::string seed = right.substr(0, kMinOverlapMerge);
    size_t from = left.size() - hi;
    while (true) {
        const size_t at = left.find(seed, from);
        if (at == std::string::npos || at > left.size() - kMinOverlapMerge) break;
        const size_t L = left.size() - at;
        if (L >= lo && L <= hi &&
            std::equal(right.begin(), right.begin() + static_cast<long>(L),
                       left.begin() + static_cast<long>(at))) {
            return L;
        }
        from = at + 1;
    }
    return 0;
}

struct Placement {
    size_t contig;
    int64_t start;       // where the contig's first base falls on the track
    int64_t end;
    int orient;          // 1 = reverse complement of the contig as written
    size_t markers;      // markers supporting this placement
};

// Where a marker of the assembly sits: which contig, at what offset, on which strand.
struct Located {
    size_t contig;
    uint32_t pos;
    int orient;
};

}  // namespace

LayoutStats layoutByModel(const OrganismModel& model, std::vector<std::string>& contigs,
                          std::vector<double>& covs, int k, bool verbose) {
    LayoutStats st;
    st.contigsIn = contigs.size();
    st.contigsOut = contigs.size();
    if (!model.loaded() || model.trackCount() == 0 || contigs.size() < 2) return st;

    // ---- locate the model's markers in this assembly -------------------------
    // A marker seen more than once names a repeat here and cannot identify a locus, so it
    // is dropped -- the same rule the adjacency join uses, for the same reason.
    std::unordered_map<uint32_t, Located> located;
    std::unordered_map<uint32_t, uint8_t> seen;
    for (size_t c = 0; c < contigs.size(); ++c) {
        forEachMarkerKmer(contigs[c], [&](uint64_t km, uint32_t pos, int orient) {
            const uint32_t id = model.markerOf(km);
            if (id == UINT32_MAX) return;
            if (seen[id] < 2) ++seen[id];
            located[id] = Located{c, pos, orient};
        }, model.markerDenom());
    }
    for (auto it = located.begin(); it != located.end();) {
        it = seen[it->first] > 1 ? located.erase(it) : std::next(it);
    }
    if (located.empty()) return st;

    // ---- choose the track this assembly most resembles ------------------------
    // Scored by containment in both directions, not by the raw count of shared markers.
    // A raw count rewards whichever panel genome simply carries the most markers -- the
    // longest chromosome, or the one with the least repeat content -- rather than the one
    // this isolate resembles, and picking on that basis put a poor template in front of a
    // good one. Dividing by the larger of the two marker sets penalises a track that is
    // large but mostly unmatched just as it penalises one that matches only a fragment.
    size_t bestTrack = SIZE_MAX, bestShared = 0;
    double bestScore = 0;
    for (size_t t = 0; t < model.tracks().size(); ++t) {
        const LayoutTrack& tr = model.tracks()[t];
        size_t shared = 0;
        for (uint32_t o : tr.oriented) {
            if (located.count(o >> 1)) ++shared;
        }
        if (shared == 0) continue;
        const double denom =
            static_cast<double>(std::max(tr.oriented.size(), located.size()));
        const double score = static_cast<double>(shared) / denom;
        if (score > bestScore) {
            bestScore = score;
            bestShared = shared;
            bestTrack = t;
        }
    }
    if (bestTrack == SIZE_MAX || bestShared < 200) return st;

    const LayoutTrack& track = model.tracks()[bestTrack];
    st.run = true;
    st.track = track.name;
    st.trackMarkers = track.oriented.size();
    st.sharedMarkers = bestShared;

    // ---- place each contig ----------------------------------------------------
    // Every shared marker votes for where its contig starts on the track. A marker at
    // track position P, found at offset `pos` in the contig, implies the contig begins at
    // P - pos when the two agree on strand, and at P + pos + k - len when they do not.
    struct Vote { int64_t start; int orient; };
    std::vector<std::vector<Vote>> votes(contigs.size());
    for (size_t i = 0; i < track.oriented.size(); ++i) {
        const uint32_t id = track.oriented[i] >> 1;
        const int trackOrient = static_cast<int>(track.oriented[i] & 1u);
        auto it = located.find(id);
        if (it == located.end()) continue;
        const Located& L = it->second;
        const int64_t tp = static_cast<int64_t>(track.pos[i]);
        const int64_t len = static_cast<int64_t>(contigs[L.contig].size());
        if (L.orient == trackOrient) {
            votes[L.contig].push_back({tp - static_cast<int64_t>(L.pos), 0});
        } else {
            votes[L.contig].push_back(
                {tp + static_cast<int64_t>(L.pos) + kMarkerK - len, 1});
        }
    }

    // A contig the panel recognises as plasmid must not be laid onto a chromosome track.
    //
    // The track is a chromosome's marker order, so placing a plasmid contig on it splices
    // plasmid sequence into the chromosome scaffold -- and plasmids share IS elements and
    // whole genes with the chromosome, so they have markers that look placeable. Measured
    // on the 666: closure lifts plasmid reconstruction from 37.8% to 73.7% of a plasmid in
    // one contig, and at the same time takes chimeric plasmid contigs from 16 to 58. A
    // chimera is wrong sequence rather than missing sequence, which is the worse failure,
    // so it is worth refusing some placements to avoid it.
    //
    // The test is deliberately one-sided: a contig is skipped only when the panel is
    // lopsidedly sure it is plasmid. Contigs with no markers either way are still placed,
    // because "unrecognised" is not evidence of anything.
    std::vector<char> plasmidContig(contigs.size(), 0);
    {
        std::unordered_map<uint32_t, uint32_t> seenAll;
        std::vector<std::pair<size_t, uint32_t>> hits;
        for (size_t c = 0; c < contigs.size(); ++c) {
            forEachMarkerKmer(contigs[c], [&](uint64_t km, uint32_t, int) {
                const uint32_t id = model.markerOf(km);
                if (id == UINT32_MAX) return;
                hits.emplace_back(c, id);
                ++seenAll[id];
            }, model.markerDenom());
        }
        std::vector<uint32_t> chrV(contigs.size(), 0), plsV(contigs.size(), 0);
        for (const auto& h : hits) {
            if (seenAll[h.second] != 1) continue;
            const uint32_t gc = model.markerGenomes(h.second, Replicon::Chromosome);
            const uint32_t gp = model.markerGenomes(h.second, Replicon::Plasmid);
            if (gc > 0 && gp == 0) ++chrV[h.first];
            else if (gp > 0 && gc == 0) ++plsV[h.first];
        }
        // Whether to act on this is a measured trade, not a fact, so it is a switch with
        // a default rather than a silent choice. Over 1,355 plasmids of the 666:
        //
        //   skip ON  (plasmid-safe)     59.0% of a plasmid in one contig, 97.9% collinear,
        //                               43 chimeric
        //   skip OFF (plasmid-complete) 73.7% in one contig, 94.2% collinear, 58 chimeric
        //
        // Off buys roughly 88 more plasmids reaching 90% and pays about 50 more that come
        // back in the wrong order plus 15 more carrying foreign sequence. The vanilla arm
        // is 99.9% collinear, so every one of those mis-orderings is introduced here.
        //
        // ON by default: a wrong sequence is worse than a missing one for anything that
        // gets reported clinically, and the residual is meant to be a work list rather
        // than a guess. TESSERACT_PLASMID_COMPLETE=1 takes the other side.
        const bool keepPlasmidsOff = [] {
            const char* e = std::getenv("TESSERACT_PLASMID_COMPLETE");
            return !(e && *e && *e != '0');
        }();
        if (keepPlasmidsOff) {
            for (size_t c = 0; c < contigs.size(); ++c) {
                if (plsV[c] >= 4 && plsV[c] >= chrV[c] * 3) {
                    plasmidContig[c] = 1;
                    ++st.plasmidSkipped;
                }
            }
        }
    }

    std::vector<Placement> placed;
    for (size_t c = 0; c < contigs.size(); ++c) {
        if (plasmidContig[c]) {
            ++st.unplaced;
            continue;
        }
        if (votes[c].size() < kMinMarkersToPlace) {
            ++st.unplaced;
            continue;
        }
        // The modal start, taken as the largest cluster of votes agreeing within a
        // tolerance. Markers scattered over a repeat family produce a wide spread, and a
        // mean over that spread would place the contig where none of its markers say.
        //
        // The tolerance has to scale with the contig. Every indel between this isolate and
        // the template shifts the start each downstream marker implies, and those shifts
        // accumulate along the contig: markers at the two ends of a 2 Mb contig can imply
        // starts tens of kilobases apart while agreeing perfectly about where the contig
        // goes. A fixed 2 kb window placed three contigs out of 142 for exactly that
        // reason. Two per cent of the contig's length tracks the accumulation; the floor
        // keeps short contigs from being judged on a window narrower than one repeat.
        std::sort(votes[c].begin(), votes[c].end(),
                  [](const Vote& a, const Vote& b) { return a.start < b.start; });
        const int64_t kTol = std::max<int64_t>(
            10000, static_cast<int64_t>(contigs[c].size()) / 25);
        size_t bestBegin = 0, bestLen = 0;
        for (size_t i = 0, j = 0; i < votes[c].size(); ++i) {
            if (j < i) j = i;
            while (j < votes[c].size() && votes[c][j].start - votes[c][i].start <= kTol) ++j;
            if (j - i > bestLen) { bestLen = j - i; bestBegin = i; }
        }
        if (bestLen < kMinMarkersToPlace ||
            static_cast<double>(bestLen) < kMinCoherent * static_cast<double>(votes[c].size())) {
            ++st.incoherent;
            ++st.unplaced;
            continue;
        }
        int fwd = 0;
        for (size_t i = bestBegin; i < bestBegin + bestLen; ++i) {
            if (votes[c][i].orient == 0) ++fwd;
        }
        const int orient = fwd * 2 >= static_cast<int>(bestLen) ? 0 : 1;
        const int64_t start = votes[c][bestBegin + bestLen / 2].start;
        placed.push_back({c, start, start + static_cast<int64_t>(contigs[c].size()),
                          orient, bestLen});
    }

    if (placed.empty()) return st;

    // Two contigs cannot occupy the same stretch of the layout. Taking the better
    // supported one first means a repeat copy never displaces the sequence that belongs.
    std::sort(placed.begin(), placed.end(),
              [](const Placement& a, const Placement& b) { return a.markers > b.markers; });
    std::vector<Placement> chosen;
    for (const Placement& p : placed) {
        bool clash = false;
        for (const Placement& q : chosen) {
            const int64_t ov = std::min(p.end, q.end) - std::max(p.start, q.start);
            if (ov > 0 && static_cast<double>(ov) >
                              kMaxPlacementOverlap * static_cast<double>(
                                  std::min(p.end - p.start, q.end - q.start))) {
                clash = true;
                break;
            }
        }
        if (clash) { ++st.dropped; ++st.unplaced; continue; }
        chosen.push_back(p);
    }
    std::sort(chosen.begin(), chosen.end(),
              [](const Placement& a, const Placement& b) { return a.start < b.start; });
    st.placed = chosen.size();

    if (std::getenv("TESSERACT_DEBUG_LAYOUT")) {
        for (size_t i = 0; i < chosen.size(); ++i) {
            const Placement& p = chosen[i];
            std::fprintf(stderr,
                         "      [layout] contig %zu len %zu -> track %lld..%lld orient %d "
                         "markers %zu%s\n",
                         p.contig, contigs[p.contig].size(),
                         static_cast<long long>(p.start), static_cast<long long>(p.end),
                         p.orient, p.markers,
                         i ? (chosen[i].start - chosen[i - 1].end > kMaxTrackGap
                                  ? "  <-- SPLIT" : "") : "");
        }
    }

    // ---- emit ------------------------------------------------------------------
    std::vector<char> used(contigs.size(), 0);
    std::vector<std::string> out;
    std::vector<double> outCov;

    std::string cur;
    double covW = 0;
    size_t covLen = 0;
    auto flush = [&]() {
        if (cur.empty()) return;
        out.push_back(std::move(cur));
        outCov.push_back(covLen ? covW / static_cast<double>(covLen) : 0.0);
        cur.clear();
        covW = 0;
        covLen = 0;
        ++st.scaffolds;
    };

    for (size_t i = 0; i < chosen.size(); ++i) {
        const Placement& p = chosen[i];
        used[p.contig] = 1;
        std::string piece =
            p.orient == 0 ? contigs[p.contig] : reverseComplement(contigs[p.contig]);

        if (!cur.empty()) {
            const int64_t gap = p.start - chosen[i - 1].end;
            if (gap > kMaxTrackGap) {
                flush();
            } else if (gap < 0) {
                const size_t trim = exactOverlap(cur, piece, -gap);
                if (trim) {
                    piece.erase(0, trim);
                    ++st.overlapMerges;
                } else {
                    cur.append(1, 'N');
                    ++st.gapBases;
                }
            } else {
                const size_t n = static_cast<size_t>(std::min<int64_t>(std::max<int64_t>(gap, 1), kGapCap));
                cur.append(n, 'N');
                st.gapBases += n;
            }
        }
        cur += piece;
        covW += covs[p.contig] * static_cast<double>(piece.size());
        covLen += piece.size();
    }
    flush();

    // Everything the track could not place is kept exactly as it was. These are the
    // contigs worth looking at: what escapes a layout is disproportionately mobile.
    for (size_t c = 0; c < contigs.size(); ++c) {
        if (used[c]) continue;
        out.push_back(contigs[c]);
        outCov.push_back(covs[c]);
    }

    contigs.swap(out);
    covs.swap(outCov);
    st.contigsOut = contigs.size();

    if (verbose) {
        std::fprintf(stderr,
                     "      laid out against %s (%zu of %zu track markers shared): "
                     "%zu contigs placed into %zu scaffolds, %zu left unplaced "
                     "(%zu incoherent, %zu displaced), %zu gap bases, %zu overlap merges\n",
                     st.track.c_str(), st.sharedMarkers, st.trackMarkers, st.placed,
                     st.scaffolds, st.unplaced, st.incoherent, st.dropped, st.gapBases,
                     st.overlapMerges);
    }
    return st;
}

}  // namespace ts
