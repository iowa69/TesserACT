// R6 housekeeping (build_v3 G-resolve; triage: "fix if cheap, flag-gated").
//
// R6: a legacy mutual join takes the interior route of whichever end is visited first (the
// lower chain index), so when the two ends chose different interiors between the same
// endpoints the emitted arm depends on node numbering. Fixture: A -> {X | Y} -> B, two arms
// of equal length and equal depth (a SNP bubble), no paired evidence can separate them,
// and the two ends meet the arms in opposite orders.
//   release: the emitted arm changes when the node ids are permuted (order dependence).
//   TESSERACT_FIX_MIRROR_ROUTE=1 (or TESSERACT_FIXES=1): the same arm under every numbering;
//   [mirrorroute] counts the disagreement.
// (R7, a repeated oriented source counted twice, is NOT fixed: 0 misassemblies in 26,508 K2
// joins, and no fixture short of a multi-round tandem chain exercises it; see README.)
#include "test_v3_resolve_util.h"

#include <map>
#include "test_env.h"   // build_v3: T40, clear the ambient TESSERACT_* first

using namespace v3r;

namespace {
void edge(ts::UnitigGraph& g, unsigned a, unsigned b) {
    g.nodes[a].ends[1].push_back({b, 0});
    g.nodes[b].ends[0].push_back({a, 1});
}
struct Out { std::vector<std::string> seqs; std::string log; };
Out run(const ts::UnitigGraph& g, const ts::SequenceStore& reads) {
    StderrCapture cap;
    cap.start();
    ts::PairedResolver r(g, reads, 1, 2, 2.0, 0.10);
    r.setScaffolding(false);
    r.buildSupport();
    Out o;
    std::vector<double> covs;
    r.resolve(o.seqs, covs);
    o.log = cap.stop();
    return o;
}
// Permute node ids of a graph.
ts::UnitigGraph permute(const ts::UnitigGraph& g, const std::vector<unsigned>& perm) {
    ts::UnitigGraph h;
    h.setK(g.k());
    h.nodes.resize(g.nodes.size());
    for (size_t i = 0; i < g.nodes.size(); ++i) {
        h.nodes[perm[i]] = g.nodes[i];
        for (auto& links : h.nodes[perm[i]].ends)
            for (auto& l : links) l.to = perm[l.to];
    }
    return h;
}
// Which arm the output joining A..B carries: 'X', 'Y', '-' (not joined).
char armOf(const Out& o, const std::string& A, const std::string& B, const std::string& X, const std::string& Y) {
    for (const auto& s0 : o.seqs) {
        for (const std::string& s : {s0, rc(s0)}) {
            if (s.find(A.substr(0, 40)) == std::string::npos || s.find(B.substr(B.size() - 40)) == std::string::npos) continue;
            if (s.find(X) != std::string::npos) return 'X';
            if (s.find(Y) != std::string::npos) return 'Y';
        }
    }
    return '-';
}
}  // namespace

int main() {
    testenv::clearTesseractEnv();
    // ---- R6 ----
    {
        const int k = 31;
        const size_t ov = k - 1;
        const std::string Aseq = dna(3000, 61), Bseq = dna(3000, 62);
        const std::string arm = "A";                    // a SNP bubble: the arms differ at one base
        const std::string armY = "C";
        // A -> X -> B and A -> Y -> B; X, Y = A[-30:] + base + B[:30] (61 bp < 2k, not anchorable).
        ts::UnitigGraph g;
        g.setK(k);
        g.nodes.resize(4);
        g.nodes[0].seq = Aseq;
        g.nodes[1].seq = Aseq.substr(Aseq.size() - ov) + arm + Bseq.substr(0, ov);
        g.nodes[2].seq = Aseq.substr(Aseq.size() - ov) + armY + Bseq.substr(0, ov);
        g.nodes[3].seq = Bseq;
        for (auto& nd : g.nodes) nd.coverage = 40;
        g.nodes[1].coverage = g.nodes[2].coverage = 20;
        // A lists its exits X, Y; B lists its entrances Y, X: each end's search meets the
        // two equal arms in the opposite order, so an exact tie is settled differently.
        edge(g, 0, 1); edge(g, 0, 2); edge(g, 2, 3); edge(g, 1, 3);
        const std::string genome = Aseq + arm + Bseq;
        std::vector<Pair> pairs;
        for (size_t p = 0; p + 400 <= genome.size(); p += 3) pairs.push_back({genome.substr(p, 100), rc(genome.substr(p + 300, 100))});
        const auto reads = store(pairs, "r6");
        const std::string X = Aseq.substr(Aseq.size() - 10) + arm + Bseq.substr(0, 10);
        const std::string Y = Aseq.substr(Aseq.size() - 10) + armY + Bseq.substr(0, 10);
        const std::vector<std::vector<unsigned>> perms = {{0, 1, 2, 3}, {3, 1, 2, 0}, {0, 2, 1, 3}, {3, 2, 1, 0}};
        for (int mode = 0; mode < 3; ++mode) {
            releaseDefaults();
            if (mode == 1) setFlag("TESSERACT_FIX_MIRROR_ROUTE", "1");
            if (mode == 2) setFlag("TESSERACT_FIXES", "1");
            std::string arms;
            long long disagreements = 0;
            for (const auto& perm : perms) {
                const Out o = run(permute(g, perm), reads);
                if (std::getenv("V3R_DUMP_LOG")) std::printf("%s", o.log.c_str());
                arms += armOf(o, Aseq, Bseq, X, Y);
                const auto mr = lines(o.log, "[mirrorroute]");
                if (mr.size() == 1) disagreements += field(mr[0], "disagreements");
            }
            std::printf("[R6 %s] arms under 4 numberings: %s, disagreements=%lld\n",
                        mode == 0 ? "release" : mode == 1 ? "FIX_MIRROR_ROUTE=1" : "FIXES=1", arms.c_str(), disagreements);
            expect(arms.find('-') == std::string::npos, "R6: A and B joined through one arm under every numbering");
            if (mode == 0) {
                expect(arms.find('X') != std::string::npos && arms.find('Y') != std::string::npos,
                       "R6 release: the emitted arm depends on node numbering (the defect)");
                expect(disagreements > 0, "R6 release: [mirrorroute] counts the disagreement");
            } else {
                expect(arms == std::string(arms.size(), arms[0]), "R6 fixed: the same arm under every numbering");
            }
        }
    }
    releaseDefaults();
    return finish("test_v3_resolve_housekeeping (R6)");
}
