// T21 (build_v3 G-resolve): enumerate() discards frames at the search budget, at kMaxNodes and
// at the 64-candidate cap without recording it, so a side whose rival was discarded looks
// like "the sole way through" and a lone candidate with no paired support is joined.
//
// Budget fixture (after combo3/verify/T21): U -> A, and U -> R -> D with R a depth-repeat
// longer than the search budget (6 kb); no read pair crosses any junction.
//   R = 5,900 bp: U's side has 2 candidates, best 0 -> refused as low support in every mode.
//   R = 7,000 bp: the R route is discarded at the budget; release joins U-A with score 0.
//   TESSERACT_FIX_TRUNC_GUARD=1 (or TESSERACT_FIXES=1): refused; [enumtrunc] refused=... > 0.
// [enumtrunc] is printed on every run, zeros included.
#include "test_v3_resolve_util.h"
#include "test_env.h"   // build_v3: T40, clear the ambient TESSERACT_* first

using namespace v3r;

namespace {
void tile(std::vector<Pair>& out, const std::string& g, size_t step, size_t len = 100, size_t frag = 300) {
    if (g.size() < frag) return;
    for (size_t p = 0; p + frag <= g.size(); p += step) out.push_back({g.substr(p, len), rc(g.substr(p + frag - len, len))});
}
void edge(ts::UnitigGraph& g, unsigned a, unsigned b) {
    g.nodes[a].ends[1].push_back({b, 0});
    g.nodes[b].ends[0].push_back({a, 1});
}
struct Out { bool joined = false; std::string line; };
Out budgetCase(size_t Lr, const char* guard, const char* umbrella) {
    releaseDefaults();
    setFlag("TESSERACT_FIX_TRUNC_GUARD", guard);
    setFlag("TESSERACT_FIXES", umbrella);
    const int k = 31;
    const size_t ov = k - 1;
    const std::string U = dna(5000, 1), Anew = dna(5000, 2), Rnew = dna(Lr, 3), Dnew = dna(5000, 4);
    ts::UnitigGraph g;
    g.setK(k);
    g.nodes.resize(4);
    g.nodes[0].seq = U;
    g.nodes[1].seq = U.substr(U.size() - ov) + Anew;
    g.nodes[2].seq = U.substr(U.size() - ov) + Rnew;
    g.nodes[3].seq = Rnew.substr(Rnew.size() - ov) + Dnew;
    for (auto& nd : g.nodes) nd.coverage = 40;
    g.nodes[2].coverage = 200;
    edge(g, 0, 1); edge(g, 0, 2); edge(g, 2, 3);
    std::vector<Pair> pairs;
    for (const auto& nd : g.nodes) tile(pairs, nd.seq, 5);
    const auto reads = store(pairs, "t21");
    StderrCapture cap;
    cap.start();
    ts::PairedResolver r(g, reads, 1, 2, 2.0, 0.10);
    r.setScaffolding(false);
    r.buildSupport();
    std::vector<std::string> seqs;
    std::vector<double> covs;
    r.resolve(seqs, covs);
    const std::string log = cap.stop();
    Out o;
    const std::string junction = U.substr(U.size() - 60) + Anew.substr(0, 60);
    for (const auto& s : seqs)
        if (s.find(junction) != std::string::npos || rc(s).find(junction) != std::string::npos) o.joined = true;
    const auto et = lines(log, "[enumtrunc]");
    o.line = et.size() == 1 ? et[0] : "";
    std::printf("[budget R=%zu guard=%s fixes=%s] JOINED_UA=%d | %s\n", Lr, guard ? guard : "-",
                umbrella ? umbrella : "-", o.joined ? 1 : 0, o.line.c_str());
    return o;
}
}  // namespace

int main() {
    testenv::clearTesseractEnv();
    {
        const Out rel = budgetCase(5900, nullptr, nullptr), fix = budgetCase(5900, "1", nullptr);
        expect(!rel.joined && !fix.joined, "R=5,900: two candidates, no support -> refused in both modes");
        expect(field(rel.line, "sides") == 0 && field(fix.line, "refused") == 0, "R=5,900: nothing truncated, nothing refused");
    }
    {
        const Out rel = budgetCase(7000, nullptr, nullptr);
        expect(rel.joined, "R=7,000 release: the budget-discarded rival is ignored and U-A joined unsupported (the defect)");
        expect(field(rel.line, "guard") == 0 && field(rel.line, "budget") > 0 && field(rel.line, "loneBelowBar") > 0 &&
               field(rel.line, "refused") == 0, "R=7,000 release: [enumtrunc] counts budget>0 loneBelowBar>0 refused=0");
        for (int umbrella = 0; umbrella < 2; ++umbrella) {
            const Out fix = umbrella ? budgetCase(7000, nullptr, "1") : budgetCase(7000, "1", nullptr);
            expect(!fix.joined, std::string("R=7,000 ") + (umbrella ? "TESSERACT_FIXES=1" : "GUARD=1") +
                                    ": the lone side with a discarded rival is refused");
            expect(field(fix.line, "guard") == 1 && field(fix.line, "refused") > 0, "R=7,000 fixed: [enumtrunc] guard=1 refused>0");
        }
        const Out optOut = budgetCase(7000, "0", "1");
        expect(optOut.joined, "R=7,000 TESSERACT_FIXES=1 with GUARD=0: release behaviour");
    }
    releaseDefaults();
    return finish("test_v3_resolve_trunc (T21)");
}
