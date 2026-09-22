#include "read_coverage.h"
#include "graph.h"
#include "seqio.h"

#include <algorithm>
#include <limits>
#include <new>
#include <stdexcept>
#include <utility>

namespace ts { namespace {
struct Target {
    Kmer key;
    uint64_t observations = 0;
    bool present = false; // zero observations must not mean absent target
};
constexpr size_t allocatorAllowance = 8192;

template<class Fn> void graphKmers(const std::string& seq, int k, Fn fn) {
    Kmer forward = 0, reverse = 0;
    int valid = 0;
    for (char letter : seq) {
        const int code = baseCode(letter);
        if (code < 0) { forward = reverse = 0; valid = 0; continue; }
        forward = pushBack(forward, code, k);
        reverse = pushFrontRc(reverse, code, k);
        if (valid < k) ++valid;
        if (valid >= k) fn(forward < reverse ? forward : reverse);
    }
}
bool add(uint64_t& value, uint64_t amount) {
    if (amount > std::numeric_limits<uint64_t>::max() - value) return false;
    value += amount;
    return true;
}
} // namespace

bool planReadCoverageMemory(size_t positions, size_t nodes, size_t byteBudget,
                            ReadCoverageMemoryPlan& plan, std::string& error) {
    ReadCoverageMemoryPlan pending;
    const size_t maximum = std::numeric_limits<size_t>::max();
    if (positions) {
        pending.targetSlots = 4;
        while (positions > pending.targetSlots - pending.targetSlots / 4) {
            if (pending.targetSlots > maximum / 2) {
                error = "observed coverage target capacity overflow"; return false;
            }
            pending.targetSlots *= 2;
        }
    }
    if (pending.targetSlots > maximum / sizeof(Target) || nodes > maximum / sizeof(double)) {
        error = "observed coverage allocation size overflow"; return false;
    }
    const size_t tableBytes = pending.targetSlots * sizeof(Target);
    const size_t outputBytes = nodes * sizeof(double);
    if (outputBytes > maximum - allocatorAllowance ||
        tableBytes > maximum - allocatorAllowance - outputBytes) {
        error = "observed coverage combined allocation size overflow"; return false;
    }
    pending.allocationBytes = tableBytes + outputBytes + allocatorAllowance;
    if (pending.allocationBytes > byteBudget) {
        error = "observed coverage needs " + std::to_string(pending.allocationBytes) +
            " additional bytes but budget permits " + std::to_string(byteBudget);
        return false;
    }
    plan = pending;
    return true;
}

bool measureObservedGraphCoverage(const UnitigGraph& graph, const SequenceStore& reads,
                                   size_t byteBudget, ReadCoverageResult& result,
                                   std::string& error) {
    const int k = graph.k();
    if (k < 1 || k > kMaxK) { error = "observed coverage requires a valid graph k"; return false; }
    size_t positionBound = 0;
    for (const auto& node : graph.nodes) {
        if (node.deleted || node.seq.size() < static_cast<size_t>(k)) continue;
        const size_t positions = node.seq.size() - static_cast<size_t>(k) + 1;
        if (positions > std::numeric_limits<size_t>::max() - positionBound) {
            error = "observed coverage graph position overflow"; return false;
        }
        positionBound += positions;
    }
    ReadCoverageMemoryPlan plan;
    if (!planReadCoverageMemory(positionBound, graph.nodes.size(), byteBudget, plan, error)) return false;
    try {
        std::vector<Target> targets(plan.targetSlots);
        ReadCoverageResult pending;
        pending.nodeDepths.resize(graph.nodes.size());
        auto& stats = pending.stats;
        stats.targetSlots = plan.targetSlots;
        stats.allocationBytes = plan.allocationBytes;
        const size_t mask = plan.targetSlots ? plan.targetSlots - 1 : 0;
        auto locate = [&](Kmer key) -> Target* {
            if (targets.empty()) return nullptr;
            size_t slot = kmerHash(key) & mask;
            while (targets[slot].present && targets[slot].key != key) slot = (slot + 1) & mask;
            return &targets[slot];
        };
        for (const auto& node : graph.nodes) {
            if (node.deleted) continue;
            ++stats.liveNodes;
            graphKmers(node.seq, k, [&](Kmer key) {
                Target* target = locate(key);
                // Every valid key is included in the preflight position bound.
                if (!target->present) { target->present = true; target->key = key; ++stats.distinctTargets; }
                ++stats.graphPositions;
            });
        }
        bool valid = true;
        for (size_t read = 0; read < reads.size() && valid; ++read) {
            forEachKmer(reads, read, k, [&](Kmer key, uint32_t) {
                if (!valid) return;
                valid = add(stats.readWindows, 1);
                Target* target = locate(key);
                if (valid && target && target->present) {
                    valid = add(target->observations, 1) && add(stats.matchedReadWindows, 1);
                }
            });
        }
        if (!valid) { error = "observed coverage read count overflow"; return false; }
        for (const auto& target : targets) if (target.present && !target.observations) ++stats.zeroTargets;
        for (size_t i = 0; i < graph.nodes.size(); ++i) {
            const auto& node = graph.nodes[i];
            if (node.deleted) { pending.nodeDepths[i] = node.coverage; continue; }
            uint64_t sum = 0, positions = 0;
            graphKmers(node.seq, k, [&](Kmer key) {
                const auto* target = locate(key);
                valid = valid && add(sum, target->observations);
                ++positions;
                if (!target->observations) ++stats.zeroPositions;
            });
            if (!valid) { error = "observed coverage node sum overflow"; return false; }
            pending.nodeDepths[i] = positions ? static_cast<double>(sum) / static_cast<double>(positions) : 0;
            if (!positions) ++stats.noValidKmerNodes;
            if (!sum) ++stats.zeroDepthNodes;
        }
        result = std::move(pending); // Publish only the complete measurement.
        return true;
    } catch (const std::bad_alloc&) {
        error = "observed coverage target/output allocation failed";
    } catch (const std::length_error&) {
        error = "observed coverage target/output allocation exceeds vector limit";
    }
    return false;
}
} // namespace ts
