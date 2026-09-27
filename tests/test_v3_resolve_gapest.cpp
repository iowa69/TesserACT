// T20 (build_v3 G-resolve): the scaffold gap of a join was the estimate of whichever port has
// the lower index (the median over ITS members' pairs to the partner's ENTRY unitig), so the
// N-run changed with node numbering; forced insert bounds without a fitted mean
// (TESSERACT_QC_INSERT, fewer than 1000 same-unitig pairs) estimated every gap from mean = 0;
// and (with T03, test_v3_resolve_gapflank) an estimate <= 0 without a verified overlap was
// written as a 1-N adjacency claim.
//
// P: A | gap | B1 (a short entry unitig) -> B2, the same reads under two node numberings.
//    release: the N-runs differ with the numbering; TESSERACT_FIX_GAP_ESTIMATE=1: identical.
// Q: forced bounds with fewer than 1000 same-unitig pairs: the fix skips scaffolding (counted).
// S: two long unitigs, reads 150 bp, fragments Normal(400, sigma). Both ports see the same
//    pairs, so the fixed estimate must EQUAL the release one: the fix does not change the
//    estimator. The short bias of that estimator at wide inserts (printed here) is NOT
//    corrected: an analytic correction regressed on a real library (README).
#include "test_v3_resolve_util.h"

#include <algorithm>
#include <cmath>
#include "test_env.h"   // build_v3: T40, clear the ambient TESSERACT_* first

using namespace v3r;

namespace {

std::vector<Pair> simulate(const std::string& g, size_t nPairs, int L, double mu, double sigma, unsigned seed,
                           int minFrag, int maxFrag) {
    std::mt19937 rng(seed);
    std::normal_distribution<double> nd(mu, sigma > 0 ? sigma : 1e-9);
    std::vector<Pair> out;
    while (out.size() < nPairs) {
        const int f = static_cast<int>(std::lround(nd(rng)));
        if (f < minFrag || f > maxFrag || static_cast<size_t>(f) > g.size()) continue;
        std::uniform_int_distribution<size_t> ud(0, g.size() - static_cast<size_t>(f));
        const size_t s = ud(rng);
        std::string r1 = g.substr(s, static_cast<size_t>(L));
        std::string r2 = rc(g.substr(s + static_cast<size_t>(f) - static_cast<size_t>(L), static_cast<size_t>(L)));
        if (rng() & 1) std::swap(r1, r2);
        out.push_back({r1, r2});
    }
    return out;
}

struct Out { size_t joins = 0, gapBases = 0; std::string log; };
Out run(const ts::UnitigGraph& g, const ts::SequenceStore& reads, bool fixed, int forcedMin = 0, int forcedMax = 0) {
    releaseDefaults();
    if (fixed) setFlag("TESSERACT_FIX_GAP_ESTIMATE", "1");
    StderrCapture cap;
    cap.start();
    ts::PairedResolver r(g, reads, 1, 2, 2.0, 0.10);
    r.setScaffolding(true);
    if (forcedMax > forcedMin) r.setInsertBounds(forcedMin, forcedMax);
    r.buildSupport();
    std::vector<std::string> seqs;
    std::vector<double> covs;
    r.resolve(seqs, covs);
    Out o;
    o.log = cap.stop();
    o.joins = r.stats().scaffoldJoins;
    o.gapBases = r.stats().gapBases;
    return o;
}
double median(std::vector<double> v) {
    if (v.empty()) return NAN;
    std::sort(v.begin(), v.end());
    return v.size() % 2 ? v[v.size() / 2] : 0.5 * (v[v.size() / 2 - 1] + v[v.size() / 2]);
}

// Median (N - truth) over `reps` replicates, release and fixed.
std::pair<double, double> biasCase(int k, double sigma, long gap, int reps) {
    std::vector<double> errRel, errFix;
    for (int rep = 0; rep < reps; ++rep) {
        const std::string G = dna(8000, 1000 + static_cast<unsigned>(rep));
        const size_t LA = 3000, LB = 3000;
        const size_t bStart = static_cast<size_t>(static_cast<long>(LA) + gap);
        const std::string locus = G.substr(0, bStart + LB);
        ts::UnitigGraph g;
        g.setK(k);
        g.nodes.resize(2);
        g.nodes[0].seq = locus.substr(0, LA);
        g.nodes[1].seq = locus.substr(bStart, LB);
        g.nodes[0].coverage = g.nodes[1].coverage = 40;
        const auto reads = store(simulate(locus, 5000, 150, 400, sigma, 77 + static_cast<unsigned>(rep), 150, 2000), "S");
        const long truth = gap + k - 1;
        const Out rel = run(g, reads, false), fix = run(g, reads, true);
        if (rel.joins == 1) errRel.push_back(static_cast<double>(static_cast<long>(rel.gapBases) - truth));
        if (fix.joins == 1) errFix.push_back(static_cast<double>(static_cast<long>(fix.gapBases) - truth));
        if (rep == 0) {
            const auto ge = lines(fix.log, "[gapest]");
            expect(ge.size() == 1 && field(ge[0], "enabled") == 1 && field(ge[0], "estimated") == 1 &&
                   field(ge[0], "changed") == 0, "[gapest] enabled=1 estimated=1 changed=0 with the fix");
        }
    }
    const double mr = median(errRel), mf = median(errFix);
    std::printf("[S k=%d sigma=%3.0f gap=%4ld] release median(N-truth)=%+.1f (n=%zu) | fixed %+.1f (n=%zu)\n", k, sigma,
                gap, mr, errRel.size(), mf, errFix.size());
    return {mr, mf};
}

// Replicates in which the two node numberings give different N-runs.
void portCase(int k, double sigma, long gap, size_t b1core, int reps) {
    size_t differRel = 0, differFix = 0, n = 0;
    std::vector<double> errFix;
    for (int rep = 0; rep < reps; ++rep) {
        const std::string G = dna(9000, 5000 + static_cast<unsigned>(rep));
        const size_t LA = 3000, LB = 3000, ov = static_cast<size_t>(k - 1);
        const size_t bStart = static_cast<size_t>(static_cast<long>(LA) + gap);
        const std::string locus = G.substr(0, bStart + b1core + LB);
        const std::string A = locus.substr(0, LA);
        const std::string B1 = locus.substr(bStart, b1core + ov);
        const std::string B2 = locus.substr(bStart + b1core, LB);
        const auto reads = store(simulate(locus, 5000, 150, 400, sigma, 99 + static_cast<unsigned>(rep), 150, 2000), "P");
        size_t rel[2] = {0, 0}, fix[2] = {0, 0};
        bool ok = true;
        for (int order = 0; order < 2; ++order) {
            ts::UnitigGraph g;
            g.setK(k);
            g.nodes.resize(3);
            const unsigned iA = order ? 2 : 0, iB1 = order ? 0 : 1, iB2 = order ? 1 : 2;
            g.nodes[iA].seq = A;
            g.nodes[iB1].seq = B1;
            g.nodes[iB2].seq = B2;
            for (auto& nd : g.nodes) nd.coverage = 40;
            g.nodes[iB1].ends[1].push_back({iB2, 0});
            g.nodes[iB2].ends[0].push_back({iB1, 1});
            const Out r = run(g, reads, false), f = run(g, reads, true);
            ok = ok && r.joins == 1 && f.joins == 1;
            rel[order] = r.gapBases;
            fix[order] = f.gapBases;
        }
        if (!ok) continue;
        ++n;
        differRel += rel[0] != rel[1];
        differFix += fix[0] != fix[1];
        errFix.push_back(static_cast<double>(static_cast<long>(fix[0]) - (gap + k - 1)));
    }
    std::printf("[P k=%d sigma=%.0f gap=%ld entryCore=%zu] numberings disagree: release %zu/%zu, fixed %zu/%zu; "
                "fixed median(N-truth)=%+.1f\n", k, sigma, gap, b1core, differRel, n, differFix, n, median(errFix));
    expect(n > 0 && differFix == 0, "P: with the fix the N-run does not depend on node numbering");
    expect(differRel > 0, "P: release N-run depends on node numbering (the defect)");
}
}  // namespace

int main() {
    testenv::clearTesseractEnv();
    const int reps = 6;
    size_t sameAll = 0, cellsAll = 0;
    for (int k : {31, 99}) {
        for (double sigma : {50.0, 200.0}) {
            for (long gap : {0L, 150L}) {
                const auto m = biasCase(k, sigma, gap, reps);
                ++cellsAll;
                if (m.first == m.second) ++sameAll;
            }
        }
    }
    std::printf("S cells where the fixed median error equals the release one: %zu/%zu\n", sameAll, cellsAll);
    expect(sameAll == cellsAll, "S: symmetric ports, the fixed estimate equals the release estimate");
    for (double sigma : {50.0, 150.0})
        for (size_t core : {60u, 150u}) portCase(31, sigma, 60, core, 6);
    {
        // Q: forced bounds, 900 pairs (< 1000 same-unitig pairs, no fitted mean).
        const std::string G = dna(8000, 4242);
        const long gap = 100;
        const std::string locus = G.substr(0, static_cast<size_t>(3000 + gap + 3000));
        ts::UnitigGraph g;
        g.setK(127);
        g.nodes.resize(2);
        g.nodes[0].seq = locus.substr(0, 3000);
        g.nodes[1].seq = locus.substr(static_cast<size_t>(3000 + gap), 3000);
        g.nodes[0].coverage = g.nodes[1].coverage = 40;
        const auto reads = store(simulate(locus, 900, 150, 300, 100, 4243, 150, 2000), "Q");
        const Out rel = run(g, reads, false, 250, 600), fix = run(g, reads, true, 250, 600);
        const auto ge = lines(fix.log, "[gapest]");
        std::printf("[Q k=127 forced [250,600], 900 pairs] release joins=%zu Ns=%zu | fixed joins=%zu Ns=%zu | %s\n",
                    rel.joins, rel.gapBases, fix.joins, fix.gapBases, ge.empty() ? "" : ge[0].c_str());
        expect(rel.joins == 1, "Q release: a join sized from mean = 0 (the defect)");
        expect(fix.joins == 0 && !ge.empty() && field(ge[0], "skippedNoMean") == 1,
               "Q fixed: forced bounds without a fitted mean -> scaffolding skipped and counted");
    }
    releaseDefaults();
    return finish("test_v3_resolve_gapest (T20)");
}
