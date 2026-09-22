// End-to-end resolver scoring regression using graph-consistent synthetic pairs.
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <vector>
#include "resolve.h"

namespace {
std::string dna(size_t length, unsigned seed) {
    std::mt19937 rng(seed);
    std::string s(length, 'A');
    for (char& c : s) c = "ACGT"[rng() % 4];
    return s;
}
void edge(ts::UnitigGraph& g, unsigned a, unsigned b) {
    g.nodes[a].ends[1].push_back({b, 0});
    g.nodes[b].ends[0].push_back({a, 1});
}
struct TempReads {
    std::filesystem::path path;
    TempReads() : path(std::filesystem::temp_directory_path() /
                      ("tesseract_shared_evidence_" + std::to_string(getpid()))) {
        std::filesystem::create_directory(path);
    }
    ~TempReads() { std::filesystem::remove_all(path); }
};
}

int main() {
    // Isolate chain arbitration from the independently tested terminal emitter.
    setenv("TESSERACT_COMMON_PREFIX", "0", 1);
    unsetenv("TESSERACT_WEIGHTED_RESOLVER_COVERAGE");
    unsetenv("TESSERACT_PREFIX_SNP_BUBBLES");
    ts::UnitigGraph g;
    g.setK(31);
    g.nodes.resize(9);
    g.nodes[0].seq = dna(300, 1);  // unique source A
    g.nodes[1].seq = g.nodes[0].seq.substr(270) + dna(160, 2); // shared repeat R
    g.nodes[2].seq = g.nodes[1].seq.substr(160) + dna(40, 3);  // short anchor B
    g.nodes[3].seq = g.nodes[2].seq.substr(40) + dna(270, 4);  // destination C
    g.nodes[4].seq = g.nodes[2].seq.substr(40) + dna(270, 5);  // destination D
    for (size_t i = 5; i < g.nodes.size(); ++i) g.nodes[i].seq = dna(150, i + 10);
    for (auto& node : g.nodes) node.coverage = 10;
    g.nodes[1].coverage = 50;
    edge(g, 0, 1);
    edge(g, 1, 2);
    edge(g, 2, 3);
    edge(g, 2, 4);
    if (!g.validate().empty()) throw std::runtime_error("invalid fixture graph");

    TempReads temp;
    const auto r1 = temp.path / "r1.fastq";
    const auto r2 = temp.path / "r2.fastq";
    std::ofstream first(r1), second(r2);
    size_t count = 0;
    auto pairs = [&](unsigned source, size_t pos, unsigned target, int number) {
        const std::string f = g.nodes[source].seq.substr(pos, 70);
        const std::string r = ts::reverseComplement(g.nodes[target].seq.substr(0, 70));
        for (int i = 0; i < number; ++i) {
            const std::string name = "@pair" + std::to_string(count++);
            first << name << "/1\n" << f << "\n+\n" << std::string(f.size(), 'I') << '\n';
            second << name << "/2\n" << r << "\n+\n" << std::string(r.size(), 'I') << '\n';
        }
    };
    pairs(0, 180, 2, 10);  // establishes A-R-B in the first merge round
    pairs(1, 50, 3, 100);  // repeat evidence favors C due to count imbalance
    pairs(1, 50, 4, 50);   // but the repeat also supports D
    pairs(2, 0, 4, 4);     // distinctive B evidence supports only D
    first.close();
    second.close();
    ts::Library lib;
    lib.r1 = r1.string();
    lib.r2 = r2.string();
    ts::SequenceStore reads;
    std::string error;
    if (!reads.load({lib}, 1, error)) throw std::runtime_error(error);

    auto run = [&](bool exclusion, bool audit) {
        setenv("TESSERACT_EXCLUDE_SHARED_REPEAT_SUPPORT", exclusion ? "1" : "0", 1);
        setenv("TESSERACT_SHARED_SUPPORT_AUDIT", audit ? "1" : "0", 1);
        ts::PairedResolver resolver(g, reads, 1, 2, 1.02, 0.02);
        resolver.setInsertBounds(0, 2000);
        resolver.buildSupport();
        if (resolver.stats().pairsLinking != count) throw std::runtime_error("fixture anchors lost");
        std::vector<std::string> contigs;
        std::vector<double> covs;
        resolver.resolve(contigs, covs);
        return contigs;
    };
    const auto baseline = run(false, false);
    const auto audited = run(false, true);
    const auto candidate = run(true, false);
    int checks = 0;
    auto check = [&](bool ok, const char* why) {
        ++checks;
        if (!ok) throw std::runtime_error(why);
    };
    check(baseline == audited, "audit mode changes output");
    auto sourcePath = [&](const std::vector<std::string>& contigs) -> std::string {
        for (const auto& c : contigs) if (c.find(g.nodes[0].seq) != std::string::npos) return c;
        throw std::runtime_error("source anchor missing");
    };
    const std::string before = sourcePath(baseline), after = sourcePath(candidate);
    check(before.find(g.nodes[3].seq) != std::string::npos, "fixture baseline does not pick C");
    check(before.find(g.nodes[4].seq) == std::string::npos, "baseline also contains D");
    check(after.find(g.nodes[4].seq) != std::string::npos, "distinctive B support did not select D");
    check(after.find(g.nodes[3].seq) == std::string::npos, "shared repeat still selected C");
    check(before.size() == after.size(), "fixture loses contiguity");
    check(baseline.size() == candidate.size(), "fixture increases contig count");
    std::cout << checks << " full resolver shared-evidence checks passed (" << count << " pairs)\n";
}
