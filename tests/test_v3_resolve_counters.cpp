// T14 + T35 (build_v3 G-resolve), output-neutral.
//
// T35: the resolver's single-copy depth theta (length-weighted median since 1.3.0), the
// repeat threshold (1.6 theta) and the estimator were never reported; report.json shows
// the legacy unweighted node median, which is not what the resolver used. Now in
// ResolveStats (theta, repeatThreshold, legacyMedian, thetaEstimator, thetaPopulation) and
// in a [resolver] line printed on every construction.
//
// T14: the release resolver flags K2 runs with (REQUIRE_SUPPORT_SINGLE, MIN_FALLBACK_DEST,
// EXCLUDE_SHARED_REPEAT_SUPPORT) printed no decision counters unless DEBUG/JOIN_TRACE was
// set, and MIN_FALLBACK_DEST never. Now one [resolveflags] line per run with the values in
// force and the decision counts; every resolver counter line also appears (run=0) on runs
// where the resolver is not constructed, via printResolverCountersNotRun().
#include "test_v3_resolve_util.h"
#include "test_env.h"   // build_v3: T40, clear the ambient TESSERACT_* first

using namespace v3r;

namespace {
void add(ts::UnitigGraph& graph, size_t mass, double coverage) {
    ts::Unitig u;
    u.seq.assign(mass + static_cast<size_t>(graph.k()) - 1, 'A');
    u.coverage = coverage;
    graph.nodes.push_back(std::move(u));
}
struct Theta { ts::ResolveStats st; std::string line; };
Theta theta(const ts::UnitigGraph& graph, const char* all, const char* eligible) {
    setFlag("TESSERACT_WEIGHTED_RESOLVER_COVERAGE", all);
    setFlag("TESSERACT_WEIGHTED_ELIGIBLE_COVERAGE", eligible);
    ts::SequenceStore reads;
    StderrCapture cap;
    cap.start();
    ts::PairedResolver resolver(graph, reads, 1, 3, 1.5);
    const std::string log = cap.stop();
    Theta t;
    t.st = resolver.stats();
    const auto l = lines(log, "[resolver]");
    t.line = l.size() == 1 ? l[0] : "";
    std::printf("  %s\n", t.line.c_str());
    return t;
}
const char* const kTags[] = {"[resolveflags]", "[enumtrunc]", "[revisit]", "[scafcycle]", "[gapest]",
                             "[gapflank]", "[cov_contrib]", "[routeorder]", "[mirrorroute]"};
}  // namespace

int main() {
    testenv::clearTesseractEnv();
    releaseDefaults();
    // ---- T35 ----
    {
        std::printf("[T35] theta reporting\n");
        ts::UnitigGraph biased;
        biased.setK(31);
        add(biased, 5000, 30);
        add(biased, 500, 90);
        for (int i = 0; i < 50; ++i) add(biased, 32, 4);
        const Theta def = theta(biased, nullptr, nullptr);
        expect(def.st.theta == 30 && def.st.legacyMedian == 4 && def.st.repeatThreshold == 30 * 1.6,
               "default: theta=30 (weighted), legacy median 4, repeat threshold 48");
        expect(std::string(def.st.thetaEstimator) == "weighted" && def.st.thetaPopulation == 52,
               "default: estimator=weighted, population=52 nodes");
        expect(field(def.line, "run") == 1 && def.line.find("theta=30.0000") != std::string::npos &&
               def.line.find("legacyMedian=4.0000") != std::string::npos, "[resolver] line carries theta and the legacy median");
        const Theta legacy = theta(biased, "0", nullptr);
        expect(legacy.st.theta == 4 && std::string(legacy.st.thetaEstimator) == "unweighted",
               "TESSERACT_WEIGHTED_RESOLVER_COVERAGE=0: theta=4, estimator=unweighted");
        const Theta elig = theta(biased, nullptr, "1");
        expect(elig.st.theta == 30 && std::string(elig.st.thetaEstimator) == "eligible_weighted",
               "TESSERACT_WEIGHTED_ELIGIBLE_COVERAGE=1: estimator=eligible_weighted");
        ts::UnitigGraph empty;
        empty.setK(99);
        add(empty, 10, 2);
        const Theta fb = theta(empty, nullptr, "1");
        expect(std::string(fb.st.thetaEstimator) == "unweighted_fallback" && fb.st.theta == 0,
               "no eligible population: estimator=unweighted_fallback, theta = the legacy value");
        setFlag("TESSERACT_WEIGHTED_RESOLVER_COVERAGE", nullptr);
        setFlag("TESSERACT_WEIGHTED_ELIGIBLE_COVERAGE", nullptr);
    }
    // ---- T14: a real resolve() in the K2 configuration ----
    {
        std::printf("[T14] K2 flags, lone candidate through a repeat (audit fixture C)\n");
        // Cached by the first resolve() of this process, so set before any run.
        setenv("TESSERACT_REQUIRE_SUPPORT_SINGLE", "1", 1);
        setenv("TESSERACT_MIN_FALLBACK_DEST", "1000000000", 1);
        setenv("TESSERACT_EXCLUDE_SHARED_REPEAT_SUPPORT", "1", 1);
        const int k = 31;
        const size_t ov = k - 1;
        const std::string U = dna(1500, 31), R = dna(400, 32), L = dna(1, 33), X = dna(1500, 34);
        const std::string genome = U + R + L + R + X;
        ts::UnitigGraph g;
        g.setK(k);
        g.nodes.resize(4);
        g.nodes[0].seq = U + R.substr(0, ov);
        g.nodes[1].seq = R;
        g.nodes[2].seq = R.substr(R.size() - ov) + L + R.substr(0, ov);
        g.nodes[3].seq = R.substr(R.size() - ov) + X;
        for (auto& nd : g.nodes) nd.coverage = 40;
        g.nodes[1].coverage = 80;
        auto edge = [&](unsigned a, unsigned b) { g.nodes[a].ends[1].push_back({b, 0}); g.nodes[b].ends[0].push_back({a, 1}); };
        edge(0, 1); edge(1, 3); edge(1, 2); edge(2, 1);
        std::vector<Pair> pairs;
        for (size_t p = 0; p + 300 <= genome.size(); p += 3) pairs.push_back({genome.substr(p, 100), rc(genome.substr(p + 200, 100))});
        const auto reads = store(pairs, "t14");
        StderrCapture cap;
        cap.start();
        ts::PairedResolver r(g, reads, 1, 2, 2.0, 0.10);
        r.setScaffolding(true);
        r.buildSupport();
        std::vector<std::string> seqs;
        std::vector<double> covs;
        r.resolve(seqs, covs);
        const std::string log = cap.stop();
        for (const char* tag : kTags) {
            const auto l = lines(log, tag);
            expect(l.size() == 1 && field(l[0], "run") == 1, std::string(tag) + " printed exactly once, run=1");
        }
        const auto rf = lines(log, "[resolveflags]");
        const std::string line = rf.empty() ? "" : rf[0];
        std::printf("  %s\n", line.c_str());
        expect(field(line, "requireSupportSingle") == 1 && field(line, "minFallbackDest") == 1000000000 &&
               field(line, "excludeSharedRepeat") == 1, "[resolveflags] reports the values in force");
        expect(field(line, "loneRepeat") >= 1 && field(line, "evals") >= 2,
               "[resolveflags] counts the REQUIRE_SUPPORT_SINGLE refusal (loneRepeat >= 1)");
        expect(lines(log, "[resolver]").size() == 1, "[resolver] printed once per construction");
    }
    // ---- T14: a run without the resolver prints the same lines, run=0 ----
    {
        std::printf("[T14] printResolverCountersNotRun\n");
        StderrCapture cap;
        cap.start();
        ts::printResolverCountersNotRun();
        const std::string log = cap.stop();
        for (const char* tag : kTags) {
            const auto l = lines(log, tag);
            expect(l.size() == 1 && field(l[0], "run") == 0, std::string(tag) + " printed exactly once, run=0");
        }
        const auto rv = lines(log, "[resolver]");
        expect(rv.size() == 1 && field(rv[0], "run") == 0, "[resolver] run=0 printed");
        const auto rf = lines(log, "[resolveflags]");
        expect(!rf.empty() && field(rf[0], "minFallbackDest") == 1000000000 && field(rf[0], "requireSupportSingle") == 1,
               "run=0 line still reports the flag values in force");
        setenv("TESSERACT_FIXES", "1", 1);
        cap.start();
        ts::printResolverCountersNotRun();
        const std::string on = cap.stop();
        expect(field(lines(on, "[gapflank]")[0], "enabled") == 1 && field(lines(on, "[revisit]")[0], "guard") == 1,
               "run=0 lines report the resolved fix flags (umbrella on)");
        unsetenv("TESSERACT_FIXES");
    }
    releaseDefaults();
    return finish("test_v3_resolve_counters (T14, T35)");
}
