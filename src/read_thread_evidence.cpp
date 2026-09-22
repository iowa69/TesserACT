#include "read_thread_evidence.h"
#include "graph.h"
#include "seqio.h"

#include <algorithm>
#include <functional>
#include <limits>
#include <map>
#include <set>
#include <tuple>
#include <unordered_map>

namespace ts {
ReadThreadPath reverseReadThreadPath(const ReadThreadPath& path) {
    ReadThreadPath result;
    result.reserve(path.size());
    for (auto i = path.rbegin(); i != path.rend(); ++i) result.push_back(*i ^ 1);
    return result;
}
namespace {
constexpr size_t kMinimumAnchors = 2;
struct Seed { uint32_t node = 0, pos = 0; bool reverse = false, ambiguous = false; };
using SeedIndex = std::unordered_map<uint64_t, Seed>;
using AnchorKey = std::tuple<uint32_t, uint8_t, int64_t>;
using AnchorCounts = std::map<AnchorKey, size_t>;
uint64_t maskFor(int k) { return (uint64_t(1) << (2 * k)) - 1; }
int orientedBase(const UnitigGraph& graph, uint64_t node, size_t pos) {
    const std::string& seq = graph.nodes[node >> 1].seq;
    if (pos >= seq.size()) return -1;
    const int code = baseCode(seq[(node & 1) ? seq.size() - pos - 1 : pos]);
    return code < 0 ? -1 : ((node & 1) ? 3 - code : code);
}
SeedIndex buildSeedIndex(const UnitigGraph& graph, int k) {
    SeedIndex index;
    const uint64_t mask = maskFor(k);
    for (size_t u = 0; u < graph.nodes.size(); ++u) {
        const auto& node = graph.nodes[u];
        if (node.deleted) continue;
        uint64_t f = 0, r = 0;
        int valid = 0;
        for (size_t p = 0; p < node.seq.size(); ++p) {
            const int c = baseCode(node.seq[p]);
            if (c < 0) { f = r = 0; valid = 0; continue; }
            f = ((f << 2) | uint64_t(c)) & mask;
            r = (r >> 2) | (uint64_t(3 - c) << (2 * (k - 1)));
            if (++valid < k) continue;
            // Refuse a nonrepresentable graph coordinate rather than wrap it.
            if (u > UINT32_MAX || p + 1 - size_t(k) > UINT32_MAX) return {};
            const uint64_t canonical = std::min(f, r);
            // Even graph k below 31 can produce reverse-complement palindromes.
            // Such a seed has no strand and cannot anchor an oriented route.
            if (f == r) continue;
            const auto inserted = index.emplace(canonical,
                Seed{uint32_t(u), uint32_t(p + 1 - size_t(k)), f != canonical, false});
            if (!inserted.second) inserted.first->second.ambiguous = true;
        }
    }
    return index;
}
AnchorCounts anchors(const UnitigGraph& graph, const std::string& read,
                     int k, const SeedIndex& index) {
    AnchorCounts counts;
    uint64_t f = 0, r = 0;
    const uint64_t mask = maskFor(k);
    for (size_t p = 0; p < read.size(); ++p) {
        const int c = baseCode(read[p]);  // unknown reads are rejected before this call
        f = ((f << 2) | uint64_t(c)) & mask;
        r = (r >> 2) | (uint64_t(3 - c) << (2 * (k - 1)));
        if (p + 1 < size_t(k)) continue;
        const auto found = index.find(std::min(f, r));
        if (found == index.end() || found->second.ambiguous) continue;
        const Seed& seed = found->second;
        const uint8_t orient = uint8_t(seed.reverse ^ (f > r));
        const int64_t position = orient
            ? int64_t(graph.nodes[seed.node].seq.size()) - seed.pos - k : seed.pos;
        const int64_t offset = position - int64_t(p + 1 - size_t(k));
        ++counts[AnchorKey{seed.node, orient, offset}];
    }
    return counts;
}
bool validTransition(const UnitigGraph& graph, uint64_t from, const Link& link,
                     size_t overlap) {
    if (link.to >= graph.nodes.size() || link.toEnd > 1 || graph.nodes[link.to].deleted)
        return false;
    const uint64_t to = (uint64_t(link.to) << 1) | UnitigGraph::enterOrient(link);
    const auto& a = graph.nodes[from >> 1];
    const auto& b = graph.nodes[link.to];
    if (a.seq.size() <= overlap || b.seq.size() <= overlap) return false;
    const int fromEnd = (from & 1) ? 0 : 1;
    if (std::find(b.ends[link.toEnd].begin(), b.ends[link.toEnd].end(),
                  Link{uint32_t(from >> 1), uint8_t(fromEnd)}) == b.ends[link.toEnd].end())
        return false;
    for (size_t i = 0; i < overlap; ++i) {
        const int x = orientedBase(graph, from, a.seq.size() - overlap + i);
        const int y = orientedBase(graph, to, i);
        if (x < 0 || x != y) return false;
    }
    return true;
}
ReadThreadPath mapRead(const UnitigGraph& graph, const SequenceStore& reads, size_t readId,
                       const std::vector<uint8_t>& eligible, int k, const SeedIndex& index,
                       const ReadThreadLimits& limits, ReadThreadStats& stats) {
    ++stats.readsExamined;
    const std::string read = reads.decode(readId);
    if (std::any_of(read.begin(), read.end(), [](char c) { return baseCode(c) < 0; })) {
        ++stats.readsUnknown; return {};
    }
    if (read.size() < size_t(k)) { ++stats.readsNoFlanks; return {}; }
    const AnchorCounts hits = anchors(graph, read, k, index);
    std::vector<std::pair<uint64_t, size_t>> starts;
    for (const auto& hit : hits) {
        const uint32_t node = std::get<0>(hit.first);
        const uint8_t orient = std::get<1>(hit.first);
        const int64_t offset = std::get<2>(hit.first);
        if (eligible[node] && hit.second >= kMinimumAnchors && offset >= 0 &&
            uint64_t(offset) < graph.nodes[node].seq.size())
            starts.emplace_back((uint64_t(node) << 1) | orient, size_t(offset));
    }
    if (starts.empty()) { ++stats.readsNoFlanks; return {}; }
    // Retain the starting coordinate too: two placements along the same walk
    // are still ambiguous. Duplicated adjacency entries are not two routes.
    using Placement = std::pair<ReadThreadPath, size_t>;
    std::set<Placement> complete;
    bool exhausted = false;
    size_t states = 0;
    const size_t overlap = size_t(graph.k() - 1);
    std::function<void(uint64_t, size_t, size_t, size_t, ReadThreadPath&)> walk;
    walk = [&](uint64_t node, size_t pos, size_t rp, size_t firstPos, ReadThreadPath& path) {
        if (exhausted || complete.size() > 1) return;
        if (++states > limits.maxSearchStates || path.size() >= limits.maxPathNodes) {
            exhausted = true; return;
        }
        const size_t nodeLength = graph.nodes[node >> 1].seq.size();
        if (pos >= nodeLength) return;
        const size_t count = std::min(nodeLength - pos, read.size() - rp);
        for (size_t i = 0; i < count; ++i)
            if (orientedBase(graph, node, pos + i) != baseCode(read[rp + i])) return;
        path.push_back(node);
        rp += count;
        if (rp == read.size()) complete.insert({path, firstPos});
        else {
            for (const Link& link : graph.exits(uint32_t(node >> 1), int(node & 1))) {
                if (!validTransition(graph, node, link, overlap)) continue;
                const uint64_t next = (uint64_t(link.to) << 1) | UnitigGraph::enterOrient(link);
                walk(next, overlap, rp, firstPos, path);
                if (exhausted || complete.size() > 1) break;
            }
        }
        path.pop_back();
    };
    for (const auto& start : starts) {
        ReadThreadPath path;
        walk(start.first, start.second, 0, start.second, path);
        if (exhausted || complete.size() > 1) break;
    }
    if (exhausted) { ++stats.readsSearchLimited; return {}; }
    if (complete.size() > 1) { ++stats.readsAmbiguous; return {}; }
    if (complete.empty()) { ++stats.readsNoExactPath; return {}; }
    ReadThreadPath path = complete.begin()->first;
    if (path.size() < 2 || (path.front() >> 1) == (path.back() >> 1) || !eligible[path.back() >> 1]) {
        ++stats.readsNoFlanks; return {};
    }
    int64_t lastOffset = int64_t(complete.begin()->second);
    for (size_t i = 0; i + 1 < path.size(); ++i)
        lastOffset -= int64_t(graph.nodes[path[i] >> 1].seq.size() - overlap);
    const AnchorKey endKey{uint32_t(path.back() >> 1), uint8_t(path.back() & 1), lastOffset};
    const auto endHits = hits.find(endKey);
    if (endHits == hits.end() || endHits->second < kMinimumAnchors) {
        ++stats.readsNoFlanks; return {};
    }
    ++stats.readsAccepted;
    ReadThreadPath reversed = reverseReadThreadPath(path);
    return std::min(path, reversed);
}
}  // namespace
ReadThreadEvidence collectReadThreadEvidence(const UnitigGraph& graph, const SequenceStore& reads,
                                             const std::vector<uint8_t>& eligibleEndpoint,
                                             const ReadThreadLimits& limits) {
    ReadThreadEvidence evidence;
    if (graph.k() < 2 || eligibleEndpoint.size() != graph.nodes.size()) return evidence;
    const int k = std::min(31, graph.k());
    const SeedIndex index = buildSeedIndex(graph, k);
    std::map<ReadThreadPath, std::vector<size_t>> byRoute;
    for (size_t read = 0; read < reads.size();) {
        const bool paired = reads.hasMate(read);
        const size_t fragment = paired ? read / 2 : reads.pairCount() + read - reads.pairedReads();
        auto first = mapRead(graph, reads, read, eligibleEndpoint, k, index, limits, evidence.stats);
        ReadThreadPath second;
        if (paired)
            second = mapRead(graph, reads, read + 1, eligibleEndpoint, k, index, limits, evidence.stats);
        read += paired ? 2 : 1;
        if (!first.empty() && !second.empty() && first != second) {
            ++evidence.stats.fragmentsConflicting; continue;
        }
        if (!first.empty() && !second.empty()) ++evidence.stats.matesDeduplicated;
        const ReadThreadPath& route = first.empty() ? second : first;
        if (route.empty()) continue;
        byRoute[route].push_back(fragment);
        ++evidence.stats.fragmentsAccepted;
    }
    for (auto& entry : byRoute) evidence.routes.push_back({entry.first, std::move(entry.second)});
    return evidence;
}
}  // namespace ts
