// Ownership-aware terminal extension: actual resolver graph fixtures.
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
#include "resolve.h"

namespace {
int checks = 0;
void check(bool ok, const char* message) { ++checks; if (!ok) throw std::runtime_error(message); }
std::string dna(size_t n, unsigned seed) {
    std::mt19937 rng(seed); std::string s(n, 'A');
    for (char& c : s) c = "ACGT"[rng() % 4]; return s;
}
void link(ts::UnitigGraph& g, uint32_t a, uint32_t b) {
    g.nodes[a].ends[1].push_back({b, 0}); g.nodes[b].ends[0].push_back({a, 1});
}
ts::UnitigGraph graph(bool repeated, bool continueChain, bool reverse) {
    ts::UnitigGraph g; g.setK(31); g.nodes.resize(7);
    const std::string overlap = dna(30, 987);
    g.nodes[0].seq = dna(370, 1) + overlap;
    g.nodes[1].seq = dna(370, 2) + overlap;
    g.nodes[2].seq = overlap + dna(370, 3);
    g.nodes[3].seq = g.nodes[2].seq.substr(370) + dna(370, 4);
    for (size_t i=4;i<g.nodes.size();++i) g.nodes[i].seq = dna(400, unsigned(i+1));
    for (auto& n:g.nodes) n.coverage=10;
    if (repeated) g.nodes[2].coverage=100;
    link(g,0,2); link(g,1,2); if(continueChain) link(g,2,3);
    if(reverse) for(auto& n:g.nodes) {
        n.seq=ts::reverseComplement(n.seq); std::swap(n.ends[0],n.ends[1]);
        for(auto& e:n.ends) for(auto& l:e) l.toEnd=1-l.toEnd;
    }
    return g;
}
struct Result { std::vector<std::string> seqs; std::vector<ts::ResolvedPath> paths; size_t joins; };
Result resolve(const ts::UnitigGraph& g, bool enabled) {
    setenv("TESSERACT_OWNED_ANCHOR_PREFIX",enabled?"1":"0",1);
    ts::SequenceStore reads; ts::PairedResolver resolver(g,reads,1,2,1.02,0.02);
    Result r; std::vector<double> covs; resolver.resolve(r.seqs,covs);
    r.paths=resolver.paths(); r.joins=resolver.stats().unitigsJoined;
    check(r.seqs.size()==r.paths.size(),"sequence/path count differs");
    for(size_t i=0;i<r.seqs.size();++i) {
        std::string spell;
        for(uint64_t u:r.paths[i].oriented) {
            auto part=g.oriented(uint32_t(u>>1),u&1);
            spell+=spell.empty()?part:part.substr(size_t(g.k()-1));
        }
        check(spell==r.seqs[i],"path does not spell sequence");
    }
    return r;
}
size_t copies(const Result& r,uint32_t node) {
    size_t count=0; for(const auto& p:r.paths) for(auto u:p.oriented) count+=(u>>1)==node; return count;
}
std::set<std::string> kmers(const Result& r, size_t k = 31) {
    std::set<std::string> result;
    for(const auto& s:r.seqs) for(size_t i=0;i+k<=s.size();++i) {
        auto word=s.substr(i,k); result.insert(std::min(word,ts::reverseComplement(word)));
    } return result;
}
bool containsSequence(const Result& r, const std::string& sequence, size_t minimumLength = 0) {
    const std::string reversed = ts::reverseComplement(sequence);
    for (const auto& contig : r.seqs) {
        if (contig.size() < minimumLength) continue;
        if (contig.find(sequence) != std::string::npos || contig.find(reversed) != std::string::npos)
            return true;
    }
    return false;
}
struct BiasedRepeat {
    ts::UnitigGraph graph;
    std::string source, crossingContext, uniqueSourceSeed;
};
BiasedRepeat biasedRepeat(bool reverse) {
    BiasedRepeat fixture;
    auto& g = fixture.graph;
    g.setK(127); g.nodes.resize(5);
    const std::string repeat = dna(1000, 903);
    g.nodes[0].seq = dna(320, 901) + repeat.substr(0, 126);  // 446 bases
    g.nodes[1].seq = dna(474, 902) + repeat.substr(0, 126);
    g.nodes[2].seq = repeat;
    g.nodes[3].seq = repeat.substr(874) + dna(474, 904);
    g.nodes[4].seq = repeat.substr(874) + dna(474, 905);
    for (auto& node : g.nodes) node.coverage = 20;
    // Ground truth for this adverse fixture has two loci A-R-C and B-R-D.
    // Uneven local coverage makes their shared repeat depth 30, below 1.6*20.
    // This truth label is an assertion about the synthetic fixture and is
    // never passed to the resolver as reference or copy-number information.
    g.nodes[2].coverage = 30;
    link(g,0,2); link(g,1,2); link(g,2,3); link(g,2,4);
    fixture.source = g.nodes[0].seq;
    fixture.crossingContext = fixture.source.substr(200) + repeat.substr(126, 200);
    fixture.uniqueSourceSeed = fixture.source.substr(20, 31);
    if (reverse) for (auto& node : g.nodes) {
        node.seq = ts::reverseComplement(node.seq);
        std::swap(node.ends[0], node.ends[1]);
        for (auto& end : node.ends) for (auto& edge : end) edge.toEnd = 1 - edge.toEnd;
    }
    return fixture;
}
}
int main() {
    unsetenv("TESSERACT_PREFIX_SNP_BUBBLES"); unsetenv("TESSERACT_EXACT_READ_THREADS");
    unsetenv("TESSERACT_ROUTE_DISTANCE"); unsetenv("TESSERACT_WEIGHTED_ELIGIBLE_COVERAGE");
    unsetenv("TESSERACT_WEIGHTED_RESOLVER_COVERAGE");
    for(bool reverse:{false,true}) for(bool chain:{false,true}) {
        auto g=graph(false,chain,reverse); check(g.validate().empty(),"invalid fixture");
        auto off=resolve(g,false),on=resolve(g,true);
        check(copies(off,2)==3,"fixture failed to reproduce three copies");
        check(copies(on,2)==1,"owned anchor duplicated by terminal extension");
        check(on.joins==off.joins,"existing chain joins changed");
        check(on.joins==(chain?1u:0u),"unexpected chain arbitration");
        check(kmers(off)==kmers(on),"graph k-mer content lost");
        size_t before=0,after=0; for(const auto&s:off.seqs)before+=s.size();for(const auto&s:on.seqs)after+=s.size();
        check(after<before,"duplicate sequence did not decrease");
        // This fixture also makes the limitation explicit: unsupported terminal
        // context becomes shorter; graph k-mer retention is not an NGA50 proof.
        check(on.seqs.size()==off.seqs.size(),"contig unexpectedly removed");
    }
    for(bool reverse:{false,true}) {
        auto g=graph(true,false,reverse); check(g.validate().empty(),"invalid repeat fixture");
        auto off=resolve(g,false),on=resolve(g,true);
        check(off.seqs==on.seqs,"genuine high-depth repeat extension changed");
        check(copies(on,2)==3,"unowned repeat not retained");
    }
    for (bool reverse : {false, true}) {
        const auto fixture = biasedRepeat(reverse);
        const auto& g = fixture.graph;
        check(g.validate().empty(), "invalid depth-biased genuine-repeat fixture");
        check(g.medianCoverage() == 20 && g.nodes[2].coverage <= 1.6 * g.medianCoverage(),
              "genuine repeat must be classified as an anchor in this adverse fixture");
        const auto off = resolve(g, false), on = resolve(g, true);
        // These assertions document a real limitation; they do not ask this
        // coverage-only rule to infer the fixture's hidden true copy count.
        check(copies(off, 2) > 2, "baseline should over-represent the genuine two-copy repeat");
        check(copies(on, 2) == 1, "ownership rule can under-represent a genuine two-copy repeat");
        check(containsSequence(off, fixture.crossingContext), "baseline lacks the correct local repeat context");
        check(!containsSequence(on, fixture.crossingContext), "adverse fixture must expose lost contiguous locus context");
        check(kmers(off, 127) == kmers(on, 127), "global graph k-mer set should still be retained");
        check(off.joins == 0 && on.joins == 0, "limitation must arise in terminal emission, not chain arbitration");
        check(fixture.source.size() == 446 && containsSequence(on, fixture.source),
              "short source sequence must remain in the unfiltered output");
        check(containsSequence(off, fixture.uniqueSourceSeed, 500),
              "baseline extended source should pass the cohort's 500-base evaluation cutoff");
        check(containsSequence(on, fixture.uniqueSourceSeed), "unique source bases disappeared before length filtering");
        check(!containsSequence(on, fixture.uniqueSourceSeed, 500),
              "446-base source loses unique evaluated sequence at the 500-base cutoff");
        // No QUAST result or alignment gain is asserted. A genuine reference
        // locus using these unique source bases can lose genome fraction even
        // while all unfiltered graph k-mers remain present in some contig.
    }
    unsetenv("TESSERACT_OWNED_ANCHOR_PREFIX");
    std::cout<<checks<<" prefix ownership checks passed\n";
}
