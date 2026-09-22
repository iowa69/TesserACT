#include "read_coverage.h"
#include "graph.h"
#include "seqio.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <random>
#include <string>
#include <unistd.h>

namespace {
size_t checks = 0;
void check(bool ok, const char* what) {
    ++checks;
    if (!ok) { std::fprintf(stderr, "FAIL: %s\n", what); std::exit(1); }
}
std::string directory;
std::mt19937 rng(907193);
std::string dna(size_t n) {
    std::string s(n, 'A');
    for (char& c : s) c = "ACGT"[rng() % 4];
    return s;
}
// String oracle deliberately uses none of the packed k-mer/canonical/hash API.
std::string reverse(std::string s) {
    std::reverse(s.begin(), s.end());
    for (char& c : s) {
        switch (c) {
            case 'A': case 'a': c = 'T'; break;
            case 'C': case 'c': c = 'G'; break;
            case 'G': case 'g': c = 'C'; break;
            case 'T': case 't': c = 'A'; break;
            default: c = 'N';
        }
    }
    return s;
}
template<class Fn> void strings(std::string s, size_t k, Fn fn) {
    for (char& c : s) if (c >= 'a' && c <= 'z') c -= 'a' - 'A';
    for (size_t pos = 0; pos + k <= s.size(); ++pos) {
        const std::string word = s.substr(pos, k);
        if (word.find_first_not_of("ACGT") == std::string::npos) fn(std::min(word, reverse(word)));
    }
}
ts::SequenceStore load(const std::vector<std::string>& seqs, bool paired = false) {
    std::ofstream first(directory + "/a.fa"), second(directory + "/b.fa");
    for (size_t i = 0; i < seqs.size(); ++i) {
        auto& file = paired && i % 2 ? second : first;
        file << ">read" << (paired ? i / 2 : i) << (paired ? (i % 2 ? "/2" : "/1") : "") << '\n' << seqs[i] << '\n';
    }
    first.close(); second.close();
    ts::Library lib; lib.r1 = directory + "/a.fa";
    if (paired) lib.r2 = directory + "/b.fa";
    ts::SequenceStore reads; std::string error;
    check(reads.load({lib}, 1, error), "load synthetic reads");
    return reads;
}
void node(ts::UnitigGraph& graph, const std::string& seq, bool deleted = false) {
    ts::Unitig n; n.seq = seq; n.coverage = 49.125 + graph.nodes.size(); n.deleted = deleted;
    graph.nodes.push_back(std::move(n));
}
bool sameGraph(const ts::UnitigGraph& a, const ts::UnitigGraph& b) {
    if (a.k() != b.k() || a.nodes.size() != b.nodes.size()) return false;
    for (size_t i = 0; i < a.nodes.size(); ++i) {
        const auto& x = a.nodes[i]; const auto& y = b.nodes[i];
        if (x.seq != y.seq || x.coverage != y.coverage || x.deleted != y.deleted ||
            x.ends[0] != y.ends[0] || x.ends[1] != y.ends[1]) return false;
    }
    return true;
}
ts::ReadCoverageResult compare(const ts::UnitigGraph& graph, const ts::SequenceStore& reads) {
    const auto before = graph;
    std::map<std::string, uint64_t> counts, targets;
    uint64_t readWindows = 0, matched = 0, positions = 0, zeroPositions = 0;
    size_t live = 0, noKmers = 0, zeroNodes = 0;
    for (size_t i = 0; i < reads.size(); ++i) strings(reads.decode(i), graph.k(), [&](const std::string& s) { ++counts[s]; ++readWindows; });
    std::vector<double> expected;
    for (const auto& n : graph.nodes) {
        if (n.deleted) { expected.push_back(n.coverage); continue; }
        ++live; uint64_t sum = 0, npos = 0;
        strings(n.seq, graph.k(), [&](const std::string& s) {
            const uint64_t count = counts[s]; targets[s] = count; sum += count; ++npos;
            if (!count) ++zeroPositions;
        });
        positions += npos;
        expected.push_back(npos ? double(sum) / npos : 0);
        if (!npos) ++noKmers;
        if (!sum) ++zeroNodes;
    }
    uint64_t zeroTargets = 0;
    for (const auto& t : targets) { matched += t.second; zeroTargets += t.second == 0; }
    ts::ReadCoverageResult actual; std::string error;
    check(ts::measureObservedGraphCoverage(graph, reads, 32U << 20, actual, error), "measurement succeeds");
    check(sameGraph(graph, before), "successful measurement never mutates graph");
    check(actual.nodeDepths.size() == expected.size(), "all graph slots returned");
    for (size_t i = 0; i < expected.size(); ++i) check(actual.nodeDepths[i] == expected[i], "node depth agrees with independent string oracle");
    check(actual.stats.graphPositions == positions && actual.stats.distinctTargets == targets.size(), "target multiplicity matches oracle");
    check(actual.stats.readWindows == readWindows && actual.stats.matchedReadWindows == matched, "read and matched counts agree with oracle");
    check(actual.stats.zeroTargets == zeroTargets && actual.stats.zeroPositions == zeroPositions, "zero observations explicit per key and position");
    check(actual.stats.liveNodes == live && actual.stats.noValidKmerNodes == noKmers && actual.stats.zeroDepthNodes == zeroNodes, "live/zero/short node statistics");
    check(actual.stats.allocationBytes <= (32U << 20), "allocation preflight bounds working payload");
    return actual;
}
} // namespace

int main() {
    char pattern[] = "/tmp/tesseract-read-coverage-XXXXXX";
    const char* tmp = mkdtemp(pattern); check(tmp != nullptr, "temporary directory"); directory = tmp;
    for (int k : {21, 33, 127}) {
        for (int attempt = 0; attempt < 12; ++attempt) {
            const std::string genome = dna(1200);
            std::vector<std::string> sequences;
            for (size_t i = 0; i < 24; ++i) {
                auto s = genome.substr(rng() % 700, k + 130);
                if (i % 2) s = reverse(s);
                if (i % 4 == 0) s[k / 2] = 'N';
                sequences.push_back(s);
            }
            sequences.push_back(std::string(k + 8, 'A'));
            sequences.push_back(std::string(k + 8, 'T'));
            auto reads = load(sequences, true);
            // EC edits and masks must be read from the current store, not raw bases.
            reads.setBase(3, 40, (reads.baseAt(3, 40) + 2) % 4);
            reads.maskRange(4, 19, 26);
            reads.maskRange(7, 0, reads.length(7));
            const auto masked = reads.decode(4);
            ts::UnitigGraph graph; graph.setK(k);
            for (size_t i = 0; i < 15; ++i) node(graph, genome.substr(i * 47, k + 70), i % 7 == 6);
            node(graph, graph.nodes[1].seq); // same canonical targets, no count inflation
            node(graph, reverse(graph.nodes[1].seq));
            node(graph, std::string(k + 3, 'A')); // within-node repeated target
            node(graph, std::string(k + 1, 'C')); // absent target
            node(graph, std::string(k - 1, 'A'));
            node(graph, ""); node(graph, std::string(k + 12, 'N'));
            node(graph, genome.substr(0, k + 30) + "NN" + genome.substr(500, k + 30));
            auto lower = graph.nodes[0].seq; for (char& c : lower) c += 'a' - 'A'; node(graph, lower);
            graph.nodes[0].ends[1].push_back({1, 0}); graph.nodes[1].ends[0].push_back({0, 1});
            const auto result = compare(graph, reads);
            check(reads.decode(4) == masked, "measurement preserves corrected read masks");
            const auto unchanged = graph;
            ts::ReadCoverageResult failed; failed.nodeDepths = {123.5}; failed.stats.readWindows = 77;
            std::string error;
            check(!ts::measureObservedGraphCoverage(graph, reads, 1, failed, error), "insufficient memory rejects before allocation");
            check(failed.nodeDepths == std::vector<double>{123.5} && failed.stats.readWindows == 77, "failed measurement leaves previous result untouched");
            check(sameGraph(graph, unchanged), "failed measurement leaves graph untouched");
            // RC every graph/read preserves means and target counts. New reads
            // come from decoded state so existing masks remain Ns.
            std::vector<std::string> flippedReads;
            for (size_t i = 0; i < reads.size(); ++i) flippedReads.push_back(reverse(reads.decode(i)));
            auto flippedStore = load(flippedReads, true);
            auto flippedGraph = graph; for (auto& n : flippedGraph.nodes) n.seq = reverse(n.seq);
            const auto flipped = compare(flippedGraph, flippedStore);
            check(result.nodeDepths == flipped.nodeDepths, "reverse-complement invariance");
            std::reverse(flippedGraph.nodes.begin(), flippedGraph.nodes.end());
            auto reordered = compare(flippedGraph, flippedStore);
            std::reverse(reordered.nodeDepths.begin(), reordered.nodeDepths.end());
            check(reordered.nodeDepths == result.nodeDepths, "graph node-order invariance");
        }
        auto reads = load({std::string(k + 2, 'A'), std::string(k + 2, 'T')}, true);
        ts::UnitigGraph graph; graph.setK(k);
        node(graph, std::string(k + 1, 'A')); node(graph, std::string(k + 1, 'T'));
        node(graph, std::string(k, 'C')); node(graph, std::string(k, 'G'), true);
        auto result = compare(graph, reads);
        check(result.nodeDepths[0] == 6 && result.nodeDepths[1] == 6, "both mates and repeated read windows count once each, never per graph copy");
        check(result.nodeDepths[2] == 0, "zero-support carried target explicitly zero");
        auto noReads = load({"NNNNN", "N"}); compare(graph, noReads);
        ts::UnitigGraph empty; empty.setK(k); compare(empty, reads);
        ts::UnitigGraph shortOnly; shortOnly.setK(k); node(shortOnly, std::string(k - 1, 'A')); compare(shortOnly, reads);
        graph.setK(0); const auto before = graph;
        ts::ReadCoverageResult failure; failure.nodeDepths = {12}; std::string error;
        check(!ts::measureObservedGraphCoverage(graph, reads, 32U << 20, failure, error), "invalid k rejected");
        check(failure.nodeDepths == std::vector<double>{12} && sameGraph(graph, before), "invalid-k failure transactional");
        graph.setK(ts::kMaxK + 1);
        check(!ts::measureObservedGraphCoverage(graph, reads, 32U << 20, failure, error), "oversize k rejected");
    }
    ts::ReadCoverageMemoryPlan plan; std::string error;
    const size_t maximum = std::numeric_limits<size_t>::max();
    check(!ts::planReadCoverageMemory(maximum, 1, maximum, plan, error), "target capacity overflow rejected without allocation");
    check(!ts::planReadCoverageMemory(0, maximum, maximum, plan, error), "output allocation overflow rejected");
    check(!ts::planReadCoverageMemory(maximum / 32, 0, maximum, plan, error), "target byte-size overflow rejected");
    check(ts::planReadCoverageMemory(100, 10, maximum, plan, error), "ordinary memory plan succeeds");
    const size_t exact = plan.allocationBytes;
    check(ts::planReadCoverageMemory(100, 10, exact, plan, error), "exact allocation budget succeeds");
    check(!ts::planReadCoverageMemory(100, 10, exact - 1, plan, error), "one-byte short budget rejects");
    std::filesystem::remove_all(directory);
    std::printf("observed read coverage: %zu checks passed\n", checks);
}
