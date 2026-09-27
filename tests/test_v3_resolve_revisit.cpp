// T07 (build_v3 G-resolve): enumerate() never extends a frame to a node it already holds,
// so a tandem R -> L -> R is traversed ONCE whatever the genome has, and the join deletes
// whole tandem units (QUAST: local misassemblies of +1/2/3x the loop length, or indels).
//
// Fixture (after combo3/verify/T07): genome U R (L R)^(c-1) X, graph U->R, R->L, L->R, R->X,
// R = 100 bp (collapsed, depth c x), L = one base (a 61 bp unitig at k=31, not anchorable),
// U and X 6 kb, FR pairs with 100 bp reads and fragments uniform in 400-600 (pairs span the
// tandem, so the collapsed route is accepted WITH paired support). Cases: c0 = control (no
// tandem, no loop edge), c1 = the audit fixture (R=400, fixed 300 bp fragments; no pair spans;
// flank 1.5 kb), c2 = 2 copies, c3 = 3 copies.
//   release:                                         c1/c2/c3 joined with ONE copy (wrong)
//   TESSERACT_FIX_REVISIT_GUARD=1 (or TESSERACT_FIXES=1): c1/c2/c3 refused (not joined);
//                                                    the c0 control is unchanged
#include "test_v3_resolve_util.h"
#include "test_env.h"   // build_v3: T40, clear the ambient TESSERACT_* first

using namespace v3r;

namespace {

void tile(std::vector<Pair>& out, const std::string& g, size_t step, size_t len, size_t fmin, size_t fmax,
          unsigned seed) {
    std::mt19937 rng(seed);
    for (size_t p = 0; p + fmax <= g.size(); p += step) {
        const size_t frag = fmin + (fmax > fmin ? rng() % (fmax - fmin + 1) : 0);
        out.push_back({g.substr(p, len), rc(g.substr(p + frag - len, len))});
    }
}

enum class Verdict { Collapsed, Correct, NotJoined, Other };
const char* verdictName(Verdict v) {
    switch (v) {
        case Verdict::Collapsed: return "joined with the WRONG copy number";
        case Verdict::Correct: return "joined with the TRUE copy number";
        case Verdict::NotJoined: return "not joined";
        default: return "other";
    }
}

struct Result { Verdict v = Verdict::Other; std::string log; };

Result runCase(const std::string& which, const char* guard, const char* umbrella) {
    releaseDefaults();
    setFlag("TESSERACT_FIX_REVISIT_GUARD", guard);
    setFlag("TESSERACT_FIXES", umbrella ? umbrella : "0");   // 1.4.0: unset would follow the default umbrella (on)
    const int k = 31;
    const size_t ov = k - 1;
    size_t rlen = 100, fmin = 400, fmax = 600, copies = 2, flank = 6000;
    if (which == "c1") { rlen = 400; fmin = fmax = 300; copies = 2; flank = 1500; }
    if (which == "c3") copies = 3;
    if (which == "c0") copies = 1;
    const std::string U = dna(flank, 31), R = dna(rlen, 32), L = dna(1, 33), X = dna(flank, 34);
    std::string genome = U + R;
    for (size_t c = 1; c < copies; ++c) genome += L + R;
    genome += X;
    ts::UnitigGraph g;
    g.setK(k);
    g.nodes.resize(4);
    g.nodes[0].seq = U + R.substr(0, ov);
    g.nodes[1].seq = R;
    g.nodes[2].seq = R.substr(R.size() - ov) + L + R.substr(0, ov);
    g.nodes[3].seq = R.substr(R.size() - ov) + X;
    for (auto& nd : g.nodes) nd.coverage = 40;
    g.nodes[1].coverage = 40.0 * static_cast<double>(copies);
    g.nodes[2].coverage = copies > 1 ? 40.0 * static_cast<double>(copies - 1) : 40.0;
    auto edge = [&](unsigned a, unsigned b) {
        g.nodes[a].ends[1].push_back({b, 0});
        g.nodes[b].ends[0].push_back({a, 1});
    };
    edge(0, 1); edge(1, 3);
    if (copies > 1) { edge(1, 2); edge(2, 1); }
    std::vector<Pair> pairs;
    tile(pairs, genome, flank > 1500 ? 2 : 3, 100, fmin, fmax, 7);
    const ts::SequenceStore reads = store(pairs, which);
    StderrCapture cap;
    cap.start();
    ts::PairedResolver r(g, reads, 1, 2, 2.0, 0.10);
    r.setScaffolding(true);
    r.buildSupport();
    std::vector<std::string> seqs;
    std::vector<double> covs;
    r.resolve(seqs, covs);
    Result out;
    out.log = cap.stop();
    const std::string c1s = U + R + X;
    std::string truth = U + R;
    for (size_t c = 1; c < copies; ++c) truth += L + R;
    truth += X;
    bool joined = false, correct = false, collapsed = false;
    for (const auto& s0 : seqs) {
        const std::string s = s0.find(U.substr(0, 50)) != std::string::npos ? s0 : rc(s0);
        const bool u = s.find(U.substr(0, 50)) != std::string::npos;
        const bool x = s.find(X.substr(X.size() - 50)) != std::string::npos;
        if (u && x) {
            joined = true;
            if (s.find(truth) != std::string::npos) correct = true;
            else if (s.find(c1s) != std::string::npos) collapsed = true;
        }
    }
    out.v = !joined ? Verdict::NotJoined : correct ? Verdict::Correct : collapsed ? Verdict::Collapsed : Verdict::Other;
    const auto rv = lines(out.log, "[revisit]");
    std::printf("[%s guard=%s fixes=%s] %s | %s\n", which.c_str(), guard ? guard : "-",
                umbrella ? umbrella : "-", verdictName(out.v), rv.empty() ? "(no [revisit] line)" : rv[0].c_str());
    if (std::getenv("V3R_DUMP_LOG")) std::printf("%s", out.log.c_str());
    return out;
}

void checkCounter(const Result& r, int guard) {
    const auto rv = lines(r.log, "[revisit]");
    expect(rv.size() == 1, "exactly one [revisit] counter line");
    if (rv.empty()) return;
    expect(field(rv[0], "guard") == guard, "[revisit] reports the resolved guard flag");
}
}  // namespace

int main() {
    testenv::clearTesseractEnv();
    // Release: the tandem is collapsed (defect); the control joins correctly.
    for (const char* c : {"c1", "c2", "c3"}) {
        const Result r = runCase(c, nullptr, nullptr);
        expect(r.v == Verdict::Collapsed, std::string(c) + " release: joined with ONE tandem copy (the defect)");
        checkCounter(r, 0);
        const auto rv = lines(r.log, "[revisit]");
        expect(!rv.empty() && field(rv[0], "pickedSkipped") > 0 && field(rv[0], "refused") == 0,
               std::string(c) + " release: pickedSkipped > 0 is counted, nothing refused");
    }
    {
        const Result r = runCase("c0", nullptr, nullptr);
        expect(r.v == Verdict::Correct, "c0 release control: joined correctly");
    }
    // Guard, and the umbrella: every collapsed pick is refused, the control is untouched.
    for (int umbrella = 0; umbrella < 2; ++umbrella) {
        const char* g = umbrella ? nullptr : "1";
        const char* u = umbrella ? "1" : nullptr;
        const std::string mode = umbrella ? " TESSERACT_FIXES=1" : " GUARD=1";
        for (const char* c : {"c1", "c2", "c3"}) {
            const Result r = runCase(c, g, u);
            expect(r.v == Verdict::NotJoined, std::string(c) + mode + ": the collapsed join is refused");
            checkCounter(r, 1);
            const auto rv = lines(r.log, "[revisit]");
            expect(!rv.empty() && field(rv[0], "refused") > 0, std::string(c) + mode + ": refused > 0");
        }
        const Result r = runCase("c0", g, u);
        expect(r.v == Verdict::Correct, "c0" + mode + " control: unchanged (no loop, joined correctly)");
        const auto rv = lines(r.log, "[revisit]");
        expect(!rv.empty() && field(rv[0], "pickedSkipped") == 0, "c0" + mode + ": pickedSkipped=0");
    }
    // Umbrella with an explicit opt-out.
    {
        const Result r = runCase("c2", "0", "1");
        expect(r.v == Verdict::Collapsed, "c2 TESSERACT_FIXES=1 with GUARD=0: release behaviour");
    }
    releaseDefaults();
    return finish("test_v3_resolve_revisit (T07)");
}
