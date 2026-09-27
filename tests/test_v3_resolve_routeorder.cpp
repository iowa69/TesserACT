// T31 (build_v3 G-resolve): with TESSERACT_ROUTE_DISTANCE=1 the route-distance allocation
// sums floating-point weights in the order support_ holds the spans, and that order depends
// on the thread count (pair p goes to thread p % T; thread-local lists are appended in
// thread order). The SAME multiset of spans therefore gives totals that differ in the last
// bits with -t, and a decision sitting on a bar can flip.
//
// Fixture: the shipped route scenario (tests/test_resolver_route.cpp, "correctLong"): two
// arms of 100 and 200 bp between unique flanks, 41 comparable pairs, internal training.
// The resolver runs with threads 1, 2, 3, 4 and 8, JOIN_TRACE on; [routeorder-trace] prints
// the allocation totals and contrasts in %a (exact bits).
//   release (flag unset): at least one thread count gives different bits from -t 1.
//   TESSERACT_FIX_ROUTE_ORDER=1 (or TESSERACT_FIXES=1): identical bits at every -t.
// Also a direct check of the mechanism on the header: stride-T striped permutations of one
// span multiset change the bits unless the spans are put in canonical (sorted) order.
#include "test_v3_resolve_util.h"

#include <cstring>
#include <map>
#include "test_env.h"   // build_v3: T40, clear the ambient TESSERACT_* first

using namespace v3r;

namespace {
void edge(ts::UnitigGraph& g, unsigned a, unsigned b) {
    g.nodes[a].ends[1].push_back({b, 0});
    g.nodes[b].ends[0].push_back({a, 1});
}

struct Fixture { ts::UnitigGraph g; ts::SequenceStore reads; };
Fixture build() {
    Fixture fx;
    fx.g.setK(31);
    fx.g.nodes.resize(9);
    auto& n = fx.g.nodes;
    const int endpointLength = 1000;
    n[0].seq = dna(endpointLength, 100);
    const std::string suffix = dna(30, 101);
    n[1].seq = n[0].seq.substr(endpointLength - 30) + dna(70, 102) + suffix;
    n[2].seq = n[0].seq.substr(endpointLength - 30) + dna(170, 103) + suffix;
    n[3].seq = suffix + dna(endpointLength - 30, 104);
    n[4].seq = dna(4000, 105);
    for (size_t i = 5; i < n.size(); ++i) n[i].seq = dna(200, 200 + static_cast<unsigned>(i));
    for (auto& node : n) node.coverage = 10;
    n[1].coverage = 60;
    n[2].coverage = 40;
    edge(fx.g, 0, 1); edge(fx.g, 1, 3); edge(fx.g, 0, 2); edge(fx.g, 2, 3);
    std::vector<Pair> pairs;
    size_t trained = 0;
    auto pair = [&](unsigned a, int pos, unsigned b, int pos2) {
        pairs.push_back({n[a].seq.substr(static_cast<size_t>(pos), 70), rc(n[b].seq.substr(static_cast<size_t>(pos2), 70))});
    };
    for (int z = -10; z <= 10; ++z) {
        const int fragment = 600 + z * 10;
        for (int rep = 0; rep < 20 * (11 - std::abs(z)); ++rep) {
            const int start = 100 + static_cast<int>((trained * 37) % 2000);
            pair(4, start, 4, start + fragment - 70);
            ++trained;
        }
    }
    for (int i = 0; i < 41; ++i) pair(0, endpointLength - 200 + i - 20, 3, 160);
    fx.reads = store(pairs, "t31");
    return fx;
}

std::vector<std::string> traces(const Fixture& fx, int threads) {
    StderrCapture cap;
    cap.start();
    ts::PairedResolver r(fx.g, fx.reads, threads, 2, 1.02, .02);
    r.buildSupport();
    std::vector<std::string> seqs;
    std::vector<double> covs;
    r.resolve(seqs, covs);
    const std::string log = cap.stop();
    auto t = lines(log, "[routeorder-trace]");
    std::sort(t.begin(), t.end());
    const auto ro = lines(log, "[routeorder]");
    std::printf("  -t %d: %zu decision(s) %s\n", threads, t.size(), ro.empty() ? "" : ro[0].c_str());
    for (const auto& x : t) std::printf("    %s\n", x.c_str());
    return t;
}
}  // namespace

int main() {
    testenv::clearTesseractEnv();
    releaseDefaults();
    // JOIN_TRACE is cached by the first resolve() of the process; set before any run.
    setenv("TESSERACT_JOIN_TRACE", "1", 1);
    setenv("TESSERACT_ROUTE_DISTANCE", "1", 1);
    const Fixture fx = build();
    const int threadCounts[] = {1, 2, 3, 4, 8};
    for (int mode = 0; mode < 3; ++mode) {
        unsetenv("TESSERACT_FIX_ROUTE_ORDER");
        unsetenv("TESSERACT_FIXES");
        if (mode == 1) setenv("TESSERACT_FIX_ROUTE_ORDER", "1", 1);
        if (mode == 2) setenv("TESSERACT_FIXES", "1", 1);
        std::printf("[%s]\n", mode == 0 ? "flag unset" : mode == 1 ? "TESSERACT_FIX_ROUTE_ORDER=1" : "TESSERACT_FIXES=1");
        std::vector<std::string> ref;
        size_t differ = 0;
        bool anyDecision = true;
        for (int t : threadCounts) {
            const auto tr = traces(fx, t);
            if (tr.empty()) anyDecision = false;
            if (t == 1) { ref = tr; continue; }
            // Compare the totals/contrasts, ignoring the sorted= field.
            auto strip = [](std::vector<std::string> v) {
                for (auto& s : v) { const size_t p = s.find(" sorted="); if (p != std::string::npos) s.erase(p, 9); }
                return v;
            };
            if (strip(tr) != strip(ref)) ++differ;
        }
        expect(anyDecision, "the fixture reaches the route-distance allocation at every -t");
        if (mode == 0) expect(differ > 0, "release: allocation bits depend on the thread count (the defect)");
        else expect(differ == 0, "fixed: allocation bits identical at -t 1, 2, 3, 4, 8");
    }
    // Mechanism on the header alone.
    {
        std::vector<ts::detail::RouteTrainingSequence> population;
        std::vector<uint64_t> histogram(1200, 0);
        std::mt19937 rng(5);
        std::normal_distribution<double> nd(600, 60);
        for (int i = 0; i < 200000; ++i) {
            const int f = static_cast<int>(std::lround(nd(rng)));
            if (f > 0 && f < 1200) ++histogram[static_cast<size_t>(f)];
        }
        population.push_back({100000, 10});
        const auto density = ts::detail::fitRouteInsertDensity(histogram, population, 32, 1100);
        size_t changedRaw = 0, changedSorted = 0, trials = 0;
        for (int d = 0; d < 400 && density.usable; ++d) {
            std::vector<int> spans;
            const int n = 5 + static_cast<int>(rng() % 60);
            for (int i = 0; i < n; ++i) spans.push_back(420 + static_cast<int>(rng() % 200));
            const std::vector<int> lengths{30, 30 + 1 + static_cast<int>(rng() % 120)};
            std::vector<int> sorted = spans;
            std::sort(sorted.begin(), sorted.end());
            const auto base = ts::detail::routeDistanceAllocation(density, sorted, lengths);
            for (int T : {2, 3, 4, 8}) {
                std::vector<int> striped;
                for (int t = 0; t < T; ++t)
                    for (size_t i = 0; i < spans.size(); ++i)
                        if (static_cast<int>(i % static_cast<size_t>(T)) == t) striped.push_back(spans[i]);
                const auto raw = ts::detail::routeDistanceAllocation(density, striped, lengths);
                std::vector<int> canon = striped;
                std::sort(canon.begin(), canon.end());
                const auto fixedAlloc = ts::detail::routeDistanceAllocation(density, canon, lengths);
                const auto rawSpans = ts::detail::routeDistanceAllocation(density, spans, lengths);
                ++trials;
                if (std::memcmp(raw.data(), rawSpans.data(), raw.size() * sizeof(double)) != 0) ++changedRaw;
                if (std::memcmp(fixedAlloc.data(), base.data(), base.size() * sizeof(double)) != 0) ++changedSorted;
            }
        }
        std::printf("[header] %zu permuted allocations: bits changed %zu unsorted, %zu sorted\n", trials, changedRaw,
                    changedSorted);
        expect(trials > 0 && changedRaw > 0, "header: a striped permutation changes the bits (order dependence)");
        expect(trials > 0 && changedSorted == 0, "header: canonical (sorted) order gives identical bits");
    }
    unsetenv("TESSERACT_ROUTE_DISTANCE");
    releaseDefaults();
    return finish("test_v3_resolve_routeorder (T31)");
}
