#include "gap_evidence.h"

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <map>
#include <unordered_map>
#include <vector>

#include "graph.h"
#include "kmer.h"
#include "seqio.h"

namespace ts {
namespace {
constexpr uint64_t noTip = std::numeric_limits<uint64_t>::max();

struct TipProjection {
    uint64_t tip = noTip;
    size_t extension = 0;  // bases beyond this oriented unitig's terminal base
    bool ambiguous = false;
};

struct SeedPlacement {
    uint32_t unitig = 0;
    uint32_t position = 0;  // forward unitig coordinates
    bool canonicalReverse = false;
    bool ambiguous = false;
};

template <typename Fn>
void scanSeeds(const std::string& seq, int k, Fn&& fn) {
    Kmer fwd = 0, rc = 0;
    int valid = 0;
    for (size_t p = 0; p < seq.size(); ++p) {
        const int base = baseCode(seq[p]);
        if (base < 0) {
            fwd = rc = 0;
            valid = 0;
            continue;
        }
        fwd = pushBack(fwd, base, k);
        rc = pushFrontRc(rc, base, k);
        if (++valid < k || fwd == rc) continue;
        const bool reverse = rc < fwd;
        fn(reverse ? rc : fwd, p + 1 - static_cast<size_t>(k), reverse);
    }
}
}  // namespace

std::set<std::pair<uint64_t, uint64_t>> nominateOrientedGapJoins(
    const UnitigGraph& graph, const SequenceStore& reads, size_t maxDistToTip,
    int probeK, uint32_t minVotes, size_t* votedOut) {
    std::set<std::pair<uint64_t, uint64_t>> out;
    if (votedOut) *votedOut = 0;
    if (!reads.paired() || probeK < 1 || probeK > kMaxK || graph.k() < 1 || minVotes == 0)
        return out;
    const size_t n = graph.nodes.size();
    if (n > std::numeric_limits<uint32_t>::max()) return out;
    const size_t overlap = static_cast<size_t>(graph.k() - 1);
    std::vector<TipProjection> tips(2 * n);

    // Carry the exit end with every stack entry. Upon traversing a reverse link,
    // the next node's exit is toEnd, not the end of the original terminal node.
    // The added distance belongs to the downstream node already traversed.
    for (uint32_t u = 0; u < n; ++u) {
        if (graph.nodes[u].deleted) continue;
        for (unsigned end = 0; end < 2; ++end) {
            if (!graph.nodes[u].ends[end].empty()) continue;
            const uint64_t tip = (static_cast<uint64_t>(u) << 1) | end;
            struct Walk { uint64_t oriented; size_t extension; };
            std::vector<Walk> stack{{tip, 0}};
            while (!stack.empty()) {
                const Walk walk = stack.back();
                stack.pop_back();
                TipProjection& placement = tips[walk.oriented];
                if (placement.tip != noTip) {
                    if (placement.tip != tip || placement.extension != walk.extension)
                        placement.ambiguous = true;
                    continue;
                }
                placement.tip = tip;
                placement.extension = walk.extension;
                const uint32_t cur = static_cast<uint32_t>(walk.oriented >> 1);
                const unsigned back = 1 - (walk.oriented & 1);
                const Unitig& node = graph.nodes[cur];
                const size_t increment = node.seq.size() > overlap ? node.seq.size() - overlap : 0;
                if (increment > maxDistToTip || walk.extension > maxDistToTip - increment)
                    continue;
                for (const Link& link : node.ends[back]) {
                    if (link.to >= n || link.toEnd > 1 || graph.nodes[link.to].deleted) continue;
                    // Multiple predecessors may converge on this path. A predecessor
                    // with more than one onward choice does not identify one terminal tip.
                    if (graph.nodes[link.to].ends[link.toEnd].size() != 1) continue;
                    stack.push_back({(static_cast<uint64_t>(link.to) << 1) | link.toEnd,
                                     walk.extension + increment});
                }
            }
        }
    }

    std::unordered_map<Kmer, SeedPlacement, KmerHasher> index;
    for (uint32_t u = 0; u < n; ++u) {
        const Unitig& node = graph.nodes[u];
        if (node.deleted) continue;
        const TipProjection& f = tips[static_cast<size_t>(u) * 2 + 1];
        const TipProjection& r = tips[static_cast<size_t>(u) * 2];
        if ((f.tip == noTip || f.ambiguous) && (r.tip == noTip || r.ambiguous)) continue;
        scanSeeds(node.seq, probeK, [&](const Kmer& key, size_t pos, bool reverse) {
            // Index only positions in the bounded terminal neighborhood. The read's
            // actual projected start is checked below, after its strand is known.
            const bool forwardNear = f.tip != noTip && !f.ambiguous && f.extension <= maxDistToTip &&
                                     node.seq.size() - pos <= maxDistToTip - f.extension;
            const bool reverseNear = r.tip != noTip && !r.ambiguous && r.extension <= maxDistToTip &&
                                     pos + static_cast<size_t>(probeK) <= maxDistToTip - r.extension;
            if (!forwardNear && !reverseNear) return;
            auto inserted = index.emplace(key, SeedPlacement{u, static_cast<uint32_t>(pos), reverse, false});
            if (!inserted.second) inserted.first->second.ambiguous = true;
        });
    }
    if (index.empty()) return out;

    // A seed unique inside a terminal neighborhood may repeat elsewhere in the
    // graph. Reject those copies too: sequence order must never assign a repeat
    // to whichever tip happened to be indexed first.
    for (uint32_t u = 0; u < n; ++u) {
        if (graph.nodes[u].deleted) continue;
        scanSeeds(graph.nodes[u].seq, probeK, [&](const Kmer& key, size_t pos, bool reverse) {
            auto found = index.find(key);
            if (found == index.end() || found->second.ambiguous) return;
            SeedPlacement& hit = found->second;
            if (hit.unitig != u || hit.position != pos || hit.canonicalReverse != reverse)
                hit.ambiguous = true;
        });
    }

    auto readTip = [&](size_t read) -> uint64_t {
        Kmer fwd = 0, rc = 0;
        int valid = 0;
        uint64_t chosen = noTip;
        int64_t projectedStart = 0;
        unsigned anchors = 0;
        for (uint32_t p = 0; p < reads.length(read); ++p) {
            const int base = reads.baseAt(read, p);
            if (base < 0) { fwd = rc = 0; valid = 0; continue; }
            fwd = pushBack(fwd, base, probeK);
            rc = pushFrontRc(rc, base, probeK);
            if (++valid < probeK || fwd == rc) continue;
            const bool canonicalReverse = rc < fwd;
            auto found = index.find(canonicalReverse ? rc : fwd);
            if (found == index.end() || found->second.ambiguous) continue;
            const SeedPlacement& hit = found->second;
            const bool readReverse = canonicalReverse != hit.canonicalReverse;
            const TipProjection& target = tips[static_cast<size_t>(hit.unitig) * 2 + (readReverse ? 0 : 1)];
            if (target.tip == noTip || target.ambiguous) return noTip;
            const int64_t length = static_cast<int64_t>(graph.nodes[hit.unitig].seq.size());
            const int64_t orientedPosition = readReverse ? length - probeK - hit.position : hit.position;
            const int64_t readPosition = p + 1 - probeK;
            const int64_t distance = static_cast<int64_t>(target.extension) + length - orientedPosition + readPosition;
            if (distance < 0 || static_cast<uint64_t>(distance) > maxDistToTip) return noTip;
            if (chosen != noTip && (chosen != target.tip || projectedStart != distance)) return noTip;
            chosen = target.tip;
            projectedStart = distance;
            ++anchors;
        }
        // Two consistent positions protect against a lone accidental exact seed.
        return anchors >= 2 ? chosen : noTip;
    };

    std::map<std::pair<uint64_t, uint64_t>, uint32_t> votes;
    for (size_t p = 0; p < reads.pairCount(); ++p) {
        const uint64_t a = readTip(2 * p);
        if (a == noTip) continue;
        const uint64_t b = readTip(2 * p + 1);
        if (b == noTip || (a >> 1) == (b >> 1)) continue;
        uint32_t& count = votes[{std::min(a, b), std::max(a, b)}];
        if (count != std::numeric_limits<uint32_t>::max()) ++count;
    }
    if (votedOut) *votedOut = votes.size();
    // In the sequence-competition experiment, graph.closeGapsByOverlap will
    // arbitrate only after exact overlap and complexity eligibility are known.
    // The same flag enables that two-phase path in the assembler; it must not
    // send this expanded set into the legacy greedy closer.
    const char* sequenceCompetition = std::getenv("TESSERACT_GAP_SEQUENCE_COMPETITION");
    if (sequenceCompetition && std::strcmp(sequenceCompetition, "1") == 0) {
        for (const auto& vote : votes)
            if (vote.second >= minVotes) out.insert(vote.first);
        return out;
    }
    std::unordered_map<uint64_t, unsigned> partners;
    for (const auto& vote : votes) {
        if (vote.second < minVotes) continue;
        ++partners[vote.first.first];
        ++partners[vote.first.second];
    }
    for (const auto& vote : votes) {
        if (vote.second >= minVotes && partners[vote.first.first] == 1 && partners[vote.first.second] == 1)
            out.insert(vote.first);
    }
    return out;
}
}  // namespace ts
