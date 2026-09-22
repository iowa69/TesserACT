// Focused regression for terminal repeat completion. Run separately from the
// broad unit suite; no reads or reference sequence are used by this fixture.
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <random>
#include <string>
#include <vector>
#include "resolve.h"

namespace {
std::string dna(size_t n, unsigned seed) {
    std::mt19937 rng(seed);
    std::string s(n, 'A');
    for (char& b : s) b = "ACGT"[rng() % 4];
    return s;
}
void link(ts::UnitigGraph& g, unsigned u, int eu, unsigned v, int ev) {
    g.nodes[u].ends[eu].push_back({v, static_cast<uint8_t>(ev)});
    g.nodes[v].ends[ev].push_back({u, static_cast<uint8_t>(eu)});
}
struct Fixture {
    ts::UnitigGraph graph;
    std::string anchor, baseline, completed;
};
Fixture fixture(bool head = false, bool divergent = false, bool twoSnps = false) {
    Fixture f;
    auto& g = f.graph;
    g.setK(5);
    g.nodes.resize(9);
    f.anchor = dna(200, 1);
    g.nodes[0].seq = f.anchor;
    g.nodes[0].coverage = 10;
    g.nodes[1].seq = f.anchor.substr(196) + dna(56, 2);
    g.nodes[1].coverage = 70;
    g.nodes[2].seq = g.nodes[1].seq.substr(56) + "A" + "CTGA";
    g.nodes[3].seq = g.nodes[1].seq.substr(56) + "G" + "CTGA";
    if (twoSnps) g.nodes[3].seq[3] = g.nodes[3].seq[3] == 'A' ? 'C' : 'A';
    g.nodes[2].coverage = 30;
    g.nodes[3].coverage = 40;
    g.nodes[4].seq = "CTGA" + dna(46, 3);
    g.nodes[4].coverage = 70;
    for (size_t i = 5; i < 9; ++i) {
        g.nodes[i].seq = dna(100, static_cast<unsigned>(i + 1));
        g.nodes[i].coverage = 10;
    }
    link(g, 0, 1, 1, 0);
    link(g, 1, 1, 2, 0);
    link(g, 1, 1, 3, 0);
    link(g, 2, 1, 4, 0);
    if (!divergent) link(g, 3, 1, 4, 0);
    f.baseline = f.anchor + g.nodes[1].seq.substr(4);
    f.completed = f.baseline + g.nodes[3].seq.substr(4) + g.nodes[4].seq.substr(4);
    if (head) {
        for (auto& node : g.nodes) {
            node.seq = ts::reverseComplement(node.seq);
            std::swap(node.ends[0], node.ends[1]);
            for (auto& end : node.ends) for (auto& edge : end) edge.toEnd = 1 - edge.toEnd;
        }
        f.anchor = ts::reverseComplement(f.anchor);
        f.baseline = ts::reverseComplement(f.baseline);
        f.completed = ts::reverseComplement(f.completed);
    }
    return f;
}
std::string resolveAnchor(const Fixture& f, bool enabled) {
    setenv("TESSERACT_PREFIX_SNP_BUBBLES", enabled ? "1" : "0", 1);
    ts::SequenceStore reads;
    ts::PairedResolver resolver(f.graph, reads, 1, 2, 1.02, 0.02);
    std::vector<std::string> seqs;
    std::vector<double> covs;
    resolver.resolve(seqs, covs);
    for (size_t i = 0; i < seqs.size(); ++i) {
        if (seqs[i].find(f.anchor) == std::string::npos) continue;
        // The exposed path must spell the exact sequence at both orientations.
        std::string spelled;
        for (uint64_t oid : resolver.paths()[i].oriented) {
            const std::string piece = f.graph.oriented(static_cast<uint32_t>(oid >> 1), oid & 1);
            spelled += spelled.empty() ? piece : piece.substr(4);
        }
        if (spelled != seqs[i]) throw std::runtime_error("path/sequence mismatch");
        return seqs[i];
    }
    throw std::runtime_error("anchor missing");
}
}
int main() {
    int checks = 0;
    auto check = [&](bool ok, const char* message) {
        ++checks;
        if (!ok) throw std::runtime_error(message);
    };
    for (bool head : {false, true}) {
        const auto f = fixture(head);
        check(f.graph.validate().empty(), "invalid graph");
        check(resolveAnchor(f, false) == f.baseline, "disabled experiment changed baseline");
        check(resolveAnchor(f, true) == f.completed, "SNP repeat not completed");
        check(f.completed.size() > f.baseline.size(), "completion must add sequence");
        const auto d = fixture(head, true);
        check(resolveAnchor(d, true) == d.baseline, "divergent destinations selected");
        const auto m = fixture(head, false, true);
        check(resolveAnchor(m, true) == m.baseline, "multi-variant route selected");
    }
    // Exactly one accepted join must be reported by the actual arbitration
    // point when TESSERACT_JOIN_TRACE=1 is used to run this executable.
    ts::UnitigGraph linear;
    linear.setK(5);
    linear.nodes.resize(2);
    linear.nodes[0].seq = dna(100, 40);
    linear.nodes[1].seq = linear.nodes[0].seq.substr(96) + dna(96, 41);
    linear.nodes[0].coverage = linear.nodes[1].coverage = 10;
    link(linear, 0, 1, 1, 0);
    ts::SequenceStore reads;
    ts::PairedResolver resolver(linear, reads, 1, 2, 1.02, 0.02);
    std::vector<std::string> contigs;
    std::vector<double> covs;
    resolver.resolve(contigs, covs);
    check(contigs.size() == 1, "linear graph must join");
    check(contigs.front() == linear.nodes[0].seq + linear.nodes[1].seq.substr(4),
          "accepted join spells wrong sequence");
    check(resolver.stats().unitigsJoined == 1, "accepted join count mismatch");
    unsetenv("TESSERACT_PREFIX_SNP_BUBBLES");
    std::cout << checks << " resolver prefix checks passed\n";
}
