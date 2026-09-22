// Standalone scientific regression tests for the default-off gap-evidence arm.
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <random>
#include <string>
#include <utility>
#include <vector>
#include <unistd.h>

#include "gap_evidence.h"
#include "graph.h"
#include "seqio.h"

namespace {
int checks = 0;
int failures = 0;
void check(bool ok, const char* what) {
    ++checks;
    if (!ok) { ++failures; std::fprintf(stderr, "FAIL: %s\n", what); }
}
std::mt19937 randomBases(948023);
std::string sequence(size_t length) {
    std::string seq(length, 'A');
    for (char& c : seq) c = "ACGT"[randomBases() % 4];
    return seq;
}
std::string rc(const std::string& seq) {
    std::string out;
    for (auto i = seq.rbegin(); i != seq.rend(); ++i) out += "TGCA"[ts::baseCode(*i)];
    return out;
}
ts::UnitigGraph graph(const std::vector<std::string>& seqs) {
    ts::UnitigGraph g;
    g.setK(5);
    for (const std::string& s : seqs) {
        ts::Unitig u;
        u.seq = s;
        g.nodes.push_back(u);
    }
    return g;
}
void link(ts::UnitigGraph& g, unsigned a, unsigned ae, unsigned b, unsigned be) {
    g.nodes[a].ends[ae].push_back({b, static_cast<uint8_t>(be)});
    g.nodes[b].ends[be].push_back({a, static_cast<uint8_t>(ae)});
}
ts::SequenceStore reads(const std::vector<std::pair<std::string, std::string>>& pairs) {
    const auto dir = std::filesystem::temp_directory_path() / ("tesseract-gap-test-" + std::to_string(getpid()));
    std::filesystem::create_directories(dir);
    const auto first = dir / "r1.fa", second = dir / "r2.fa";
    {
        std::ofstream f(first), s(second);
        for (size_t i = 0; i < pairs.size(); ++i) {
            f << ">" << i << "\n" << pairs[i].first << "\n";
            s << ">" << i << "\n" << pairs[i].second << "\n";
        }
    }
    ts::SequenceStore store;
    ts::Library lib;
    lib.r1 = first.string(); lib.r2 = second.string();
    std::string error;
    if (!store.load({lib}, 1, error)) {
        std::fprintf(stderr, "load failed: %s\n", error.c_str());
        std::exit(2);
    }
    std::filesystem::remove_all(dir);
    return store;
}
auto nominate(const ts::UnitigGraph& g, const ts::SequenceStore& r, size_t distance = 500) {
    return ts::nominateOrientedGapJoins(g, r, distance, 15, 2);
}
}

int main() {
    const std::string a = sequence(180), b = sequence(180), c = sequence(180);
    const auto g = graph({a, b});
    const std::pair<std::string, std::string> fr{a.substr(100, 60), rc(b.substr(20, 60))};
    auto paired = reads({fr, fr});
    const auto baseline = nominate(g, paired);
    check(baseline.size() == 1 && baseline.count({1, 2}) == 1, "FR mates nominate actual facing ends");
    check(nominate(graph({rc(a), b}), paired).count({0, 2}) == 1, "unitig reverse-complement representation preserves biological pair");
    check(nominate(graph({b, a}), paired).count({0, 3}) == 1, "node reordering does not change biological pair");
    check(nominate(g, reads({fr})).empty(), "one read pair cannot meet two-vote threshold");
    check(ts::nominateOrientedGapJoins(g, paired, 500, 0, 2).empty(), "invalid k does not shift unchecked");
    const std::pair<std::string, std::string> ff{fr.first, b.substr(20, 60)};
    const auto opposite = nominate(g, reads({ff, ff}));
    check(opposite.count({1, 2}) == 0 && opposite.count({1, 3}) == 1, "reverse-complementing a mate changes nominated end");

    auto duplicate = graph({a, b, a});
    check(nominate(duplicate, paired).empty(), "identical graph copies never use first-writer placement");
    // The duplicate is internal to a distant connected path, outside all indexed
    // terminal neighborhoods. The whole-graph duplicate screen must still reject it.
    auto external = graph({a, b, sequence(900), a, sequence(900)});
    link(external, 2, 1, 3, 0); link(external, 3, 1, 4, 0);
    check(nominate(external, paired, 500).empty(), "repeat outside tip neighborhood is not falsely unique");

    const std::pair<std::string, std::string> conflict{a.substr(90, 30) + c.substr(90, 30), fr.second};
    check(nominate(graph({a, b, c}), reads({conflict, conflict})).empty(), "conflicting placements within one read are rejected");
    const std::pair<std::string, std::string> toC{fr.first, rc(c.substr(20, 60))};
    check(nominate(graph({a, b, c}), reads({fr, fr, toC, toC})).empty(), "competing eligible gap partners remain unresolved");

    const std::string upstream = sequence(180);
    const std::string middleForward = upstream.substr(upstream.size() - 4) + sequence(176);
    const std::string terminal = middleForward.substr(middleForward.size() - 4) + sequence(176);
    auto inverted = graph({terminal, rc(middleForward), upstream, b});
    link(inverted, 2, 1, 1, 1); link(inverted, 1, 0, 0, 0);
    const std::pair<std::string, std::string> through{upstream.substr(100, 60), fr.second};
    const auto throughPairs = reads({through, through});
    const auto walked = nominate(inverted, throughPairs, 500);
    check(walked.size() == 1 && walked.count({1, 6}) == 1, "walk updates orientation through reversing links");
    check(nominate(inverted, throughPairs, 400).empty(), "projected read-to-tip distance is enforced");

    auto branched = inverted;
    ts::Unitig branch; branch.seq = sequence(180); branched.nodes.push_back(branch);
    link(branched, 1, 0, 4, 0);
    check(nominate(branched, throughPairs).empty(), "branch in onward direction prevents terminal projection");

    auto far = graph({sequence(900), upstream, b});
    link(far, 1, 1, 0, 0);
    const auto farVotes = nominate(far, throughPairs, 500);
    check(farVotes.empty(), "neighborhood cost uses long downstream node, not short predecessor");
    std::printf("gap evidence: %d checks, %d failures\n", checks, failures);
    return failures ? 1 : 0;
}
