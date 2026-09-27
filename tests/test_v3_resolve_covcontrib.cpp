// T24 (build_v3 G-resolve): contig coverage weighted every unitig by its FULL length, so a
// k-long connector that adds one base carried k bases of weight and a contig of long
// unique unitigs joined by short deep connectors reported an inflated depth (cov_ in
// names; the replicon depth call).
//
// Fixture (after combo3/verify/T24): A (3000 bp, 40x) -> C (one 31-mer, 1000x) -> B (40x)
// resolves to one 6,001 bp contig.
//   release: 44.9101x (full-length weighting)
//   TESSERACT_FIX_COV_CONTRIB=1 (or TESSERACT_FIXES=1): 40.1600x = the contributed-base mean,
//   within 0.05x of the graph's own k-mer-count merge rule (40.1608x).
// Plus the scaffold case with T03 on: the Ns never count, and the piece after the gap is
// weighted by the bases the resolver appended (its first k-1 are restored only after gap
// filling, outside the resolver).
#include "test_v3_resolve_util.h"

#include <cmath>
#include "test_env.h"   // build_v3: T40, clear the ambient TESSERACT_* first

using namespace v3r;

namespace {
struct Out { std::vector<std::string> seqs; std::vector<double> covs; std::string log; };
Out run(const ts::UnitigGraph& g, const ts::SequenceStore& reads) {
    StderrCapture cap;
    cap.start();
    ts::PairedResolver r(g, reads, 1, 2, 2.0, 0.10);
    r.setScaffolding(true);
    r.buildSupport();
    Out o;
    r.resolve(o.seqs, o.covs);
    o.log = cap.stop();
    return o;
}
}  // namespace

int main() {
    testenv::clearTesseractEnv();
    const int k = 31;
    const size_t ov = static_cast<size_t>(k - 1);
    const std::string G = dna(6001, 7);
    ts::UnitigGraph g;
    g.setK(k);
    g.nodes.resize(3);
    g.nodes[0].seq = G.substr(0, 3000);
    g.nodes[1].seq = G.substr(3000 - ov, k);
    g.nodes[2].seq = G.substr(3001 - ov);
    g.nodes[0].coverage = 40;
    g.nodes[1].coverage = 1000;
    g.nodes[2].coverage = 40;
    g.nodes[0].ends[1].push_back({1, 0}); g.nodes[1].ends[0].push_back({0, 1});
    g.nodes[1].ends[1].push_back({2, 0}); g.nodes[2].ends[0].push_back({1, 1});
    std::vector<Pair> pairs;
    for (size_t p = 0; p + 300 <= G.size(); p += 5) pairs.push_back({G.substr(p, 100), rc(G.substr(p + 200, 100))});
    const auto reads = store(pairs, "t24");
    const double lA = 3000, lC = k, lB = static_cast<double>(g.nodes[2].seq.size());
    const double released = (40 * lA + 1000 * lC + 40 * lB) / (lA + lC + lB);
    const double contrib = (40 * lA + 1000 * (lC - ov) + 40 * (lB - ov)) / (lA + (lC - ov) + (lB - ov));
    const double kmers = (40 * (lA - ov) + 1000 * (lC - ov) + 40 * (lB - ov)) / ((lA - ov) + (lC - ov) + (lB - ov));
    struct Mode { const char* fix; const char* umbrella; bool on; };
    for (const Mode m : {Mode{nullptr, nullptr, false}, Mode{"1", nullptr, true}, Mode{nullptr, "1", true},
                         Mode{"0", "1", false}}) {
        releaseDefaults();
        setFlag("TESSERACT_FIX_COV_CONTRIB", m.fix);
        setFlag("TESSERACT_FIXES", m.umbrella);
        const Out o = run(g, reads);
        size_t idx = o.seqs.size();
        for (size_t i = 0; i < o.seqs.size(); ++i) if (o.seqs[i] == G || o.seqs[i] == rc(G)) idx = i;
        const auto cl = lines(o.log, "[cov_contrib]");
        std::printf("[fix=%s fixes=%s] cov=%.4f (full %.4f, contributed %.4f, k-mer %.4f) | %s\n", m.fix ? m.fix : "-",
                    m.umbrella ? m.umbrella : "-", idx < o.seqs.size() ? o.covs[idx] : -1.0, released, contrib, kmers,
                    cl.empty() ? "" : cl[0].c_str());
        expect(idx < o.seqs.size(), "A-C-B emitted as one contig spelling the locus");
        if (idx == o.seqs.size()) continue;
        expect(cl.size() == 1 && field(cl[0], "enabled") == (m.on ? 1 : 0) && field(cl[0], "full_vs_contributed_off10") == 1,
               "[cov_contrib] reports the flag and 1 path >10% apart");
        if (!m.on) expect(std::fabs(o.covs[idx] - released) < 1e-9, "release: full-length weighting (44.91x, the defect)");
        else {
            expect(std::fabs(o.covs[idx] - contrib) < 1e-9, "fixed: contributed-base weighted mean (40.16x)");
            expect(std::fabs(o.covs[idx] - kmers) < 0.05, "fixed: within 0.05x of the k-mer-count merge rule");
        }
    }
    // Scaffold with T03: two 3 kb unitigs at 40x and 60x, 100 bp apart. Contributed weighting
    // counts A whole and B less its first k-1 (never the Ns): (40*3000 + 60*2970)/5970.
    {
        const std::string T = dna(6100, 21);
        ts::UnitigGraph s;
        s.setK(k);
        s.nodes.resize(2);
        s.nodes[0].seq = T.substr(0, 3000);
        s.nodes[1].seq = T.substr(3100, 3000);
        s.nodes[0].coverage = 40;
        s.nodes[1].coverage = 60;
        std::vector<Pair> sp;
        for (size_t p = 0; p + 300 <= T.size(); p += 5) sp.push_back({T.substr(p, 100), rc(T.substr(p + 200, 100))});
        const auto sreads = store(sp, "t24s");
        releaseDefaults();
        setFlag("TESSERACT_FIX_COV_CONTRIB", "1");
        setFlag("TESSERACT_FIX_GAP_FLANK", "1");
        const Out o = run(s, sreads);
        size_t idx = o.seqs.size();
        for (size_t i = 0; i < o.seqs.size(); ++i) if (o.seqs[i].find('N') != std::string::npos) idx = i;
        expect(idx < o.seqs.size(), "scaffold emitted");
        if (idx < o.seqs.size()) {
            const double expectCov = (40.0 * 3000 + 60.0 * 2970) / 5970.0;
            std::printf("[scaffold T03+T24] cov=%.4f expected %.4f\n", o.covs[idx], expectCov);
            expect(std::fabs(o.covs[idx] - expectCov) < 1e-9, "T03+T24: the bases the resolver appended, Ns never weighted");
        }
    }
    releaseDefaults();
    return finish("test_v3_resolve_covcontrib (T24)");
}
