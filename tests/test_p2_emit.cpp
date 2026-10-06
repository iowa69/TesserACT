// Phase 2 W1 emit-A (src/p2_emit.*, src/p2_libguard.h), unit level:
//   R2  the orientation classes and the guard decision at its thresholds;
//   R3  graph closure in every mode the prototype (graph_close.py) knows, isolation, the
//       self-overlap helpers, junction verification on simulated reads (a true circle passes,
//       a chromosome fragment called circular fails), and the record stage in verify and close
//       mode (closed / linear duplicate-free, tags, scaffold propagation, edit log);
//   F5  the phiX174 k-mer screen and the `_spikein` label;
//   the report.json "p2" block with every flag off (enabled=0, zero counters only).
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <random>
#include <string>
#include <vector>
#include <unistd.h>

#include "emit_post.h"
#include "graph.h"
#include "p2_emit.h"
#include "p2_libguard.h"
#include "seqio.h"
#include "test_env.h"

namespace {
int checks = 0, failures = 0;
void check(bool ok, const std::string& what) {
    ++checks;
    std::printf("%s %s\n", ok ? "ok  " : "FAIL", what.c_str());
    if (!ok) ++failures;
}
std::mt19937 rng(20261006);
std::string dna(size_t n) {
    std::string s(n, 'A');
    for (char& c : s) c = "ACGT"[rng() % 4];
    return s;
}
std::string rc(const std::string& s) { return ts::reverseComplement(s); }

const int K = 31;
const size_t K1 = K - 1;

ts::UnitigGraph graphOf(const std::vector<std::pair<std::string, double>>& nodes) {
    ts::UnitigGraph g;
    g.setK(K);
    for (const auto& n : nodes) {
        ts::Unitig u;
        u.seq = n.first;
        u.coverage = n.second;
        g.nodes.push_back(u);
    }
    return g;
}
// u's end eu meets v's end ev (stored on both sides, as the graph keeps it)
void link(ts::UnitigGraph& g, uint32_t u, int eu, uint32_t v, int ev) {
    g.nodes[u].ends[eu].push_back({v, static_cast<uint8_t>(ev)});
    if (!(u == v && eu == ev)) g.nodes[v].ends[ev].push_back({u, static_cast<uint8_t>(eu)});
}
uint64_t fwd(uint32_t u) { return static_cast<uint64_t>(u) << 1; }

bool isRotation(const std::string& a, const std::string& b) {
    return a.size() == b.size() && ((b + b).find(a) != std::string::npos || (b + b).find(rc(a)) != std::string::npos);
}

// ---- R2 -------------------------------------------------------------------------------------
void testLibGuard() {
    using namespace ts::libguard;
    check(classifyOpposite(0, 150, 200, 150) == 0, "R2: forward mate upstream is inward (FR)");
    check(classifyOpposite(0, 150, 0, 150) == 0, "R2: fully overlapping mates are inward");
    check(classifyOpposite(6000, 150, 0, 150) == 1, "R2: forward mate past the reverse mate's end is outward (RF)");
    check(classifyOpposite(150, 150, 0, 150) == 1, "R2: abutting outward mates are outward");
    check(classifyOpposite(100, 150, 0, 150) == 2, "R2: overlapping offset mates are dovetail (neither)");
    check(classifyOpposite(0, 150, -30, 150) == 2, "R2: adapter read-through geometry is dovetail, never outward");

    ts::LibGuardStats s;
    s.enabled = true;
    s.pairs = 100000; s.inward = 60000; s.outward = 900; s.frObservations = 60000;
    decide(s);
    check(!s.fired && s.outwardEvaluable && std::string(verdict(s)) == "pass", "R2: a normal FR library passes");
    s = ts::LibGuardStats(); s.enabled = true;
    s.pairs = 100000; s.inward = 399; s.outward = 600; s.frObservations = 30000;
    decide(s);
    check(!s.fired && !s.outwardEvaluable, "R2: outward branch needs >= 1000 oriented pairs");
    s = ts::LibGuardStats(); s.enabled = true;
    s.pairs = 100000; s.inward = 5000; s.outward = 5001; s.frObservations = 30000;
    decide(s);
    check(s.fired && s.reason == 1, "R2: > 50 % outward fires (outward branch)");
    s = ts::LibGuardStats(); s.enabled = true;
    s.pairs = 100000; s.inward = 5000; s.outward = 5000; s.frObservations = 30000;
    decide(s);
    check(!s.fired, "R2: exactly 50 % outward does not fire");
    s = ts::LibGuardStats(); s.enabled = true;
    s.pairs = 100000; s.inward = 20000; s.outward = 10; s.frObservations = 24999;
    decide(s);
    check(s.fired && s.reason == 2, "R2: FR model on < 25 % of pairs fires (fit branch)");
    s = ts::LibGuardStats(); s.enabled = true;
    s.pairs = 100000; s.inward = 20000; s.outward = 10; s.frObservations = 25000;
    decide(s);
    check(!s.fired, "R2: FR fit of exactly 25 % does not fire");
    s = ts::LibGuardStats(); s.enabled = false;
    s.pairs = 100000; s.inward = 10; s.outward = 99990; s.frObservations = 10;
    decide(s);
    check(!s.fired && s.reason == 3, "R2: with the flag off the guard never fires (reason still measured)");
}

// ---- R3 closure ----------------------------------------------------------------------------
void testClosure() {
    const std::string P = dna(K1), Q = dna(K1);
    {   // self: one unitig whose 3' end links to its own 5' end
        const std::string A = P + dna(500) + P;
        ts::UnitigGraph g = graphOf({{A, 50}});
        link(g, 0, 1, 0, 0);
        const ts::p2::ClosureResult r = ts::p2::closeWalk(A, {fwd(0)}, g);
        check(r.mode == ts::p2::Closure::Self && r.closed && r.seq == A.substr(0, A.size() - K1),
              "R3: self closure removes the (k-1) duplicated end");
        check(r.isolated == 1 && r.outLast == 1 && r.inFirst == 1, "R3: a lone self-loop is graph-isolated");
        // an extra exit from the circle's end into another unitig: not isolated
        ts::UnitigGraph g2 = graphOf({{A, 50}, {P + dna(300), 50}});
        link(g2, 0, 1, 0, 0);
        link(g2, 0, 1, 1, 0);
        const ts::p2::ClosureResult r2 = ts::p2::closeWalk(A, {fwd(0)}, g2);
        check(r2.mode == ts::p2::Closure::Self && r2.isolated == 0, "R3: a second exit makes the closure non-isolated");
    }
    {   // via: A ... Q -> Y = Q + 3 bases + P -> A
        const std::string A = P + dna(600) + Q;
        const std::string Y = Q + "ACG" + P;
        ts::UnitigGraph g = graphOf({{A, 50}, {Y, 40}});
        link(g, 0, 1, 1, 0);
        link(g, 1, 1, 0, 0);
        const ts::p2::ClosureResult r = ts::p2::closeWalk(A, {fwd(0)}, g);
        check(r.mode == ts::p2::Closure::Via && r.closed && r.seq == A + "ACG" && r.interior == 3 && r.alts == 1,
              "R3: via closure inserts the connector's interior bases");
        check(r.isolated == 1 && r.viaSeg == "1+", "R3: a single connector at the join is isolated");
        // the same with a bubble: two connectors, the deeper one is taken
        const std::string Y2 = Q + "TTT" + P;
        ts::UnitigGraph gb = graphOf({{A, 50}, {Y, 10}, {Y2, 20}});
        link(gb, 0, 1, 1, 0); link(gb, 1, 1, 0, 0);
        link(gb, 0, 1, 2, 0); link(gb, 2, 1, 0, 0);
        const ts::p2::ClosureResult rb = ts::p2::closeWalk(A, {fwd(0)}, gb);
        check(rb.mode == ts::p2::Closure::Via && rb.seq == A + "TTT" && rb.alts == 2 && rb.isolated == 1 &&
                  rb.altDp == "1+:10.0",
              "R3: a bubble at the join takes the deeper connector and reports the other");
    }
    {   // via with a negative interior: the record's ends already overlap by 2 bases
        const std::string Qn = dna(K1 - 2) + P.substr(0, 2);
        const std::string A = P + dna(400) + Qn;
        const std::string Y = Qn + P.substr(2);
        ts::UnitigGraph g = graphOf({{A, 50}, {Y, 40}});
        link(g, 0, 1, 1, 0);
        link(g, 1, 1, 0, 0);
        const ts::p2::ClosureResult r = ts::p2::closeWalk(A, {fwd(0)}, g);
        check(r.mode == ts::p2::Closure::Via && r.closed && r.interior == -2 && r.seq == A.substr(0, A.size() - 2),
              "R3: a connector shorter than two overlaps trims the overlapping ends");
    }
    {   // repeat_end: walk A B A
        const std::string A = P + dna(300) + Q, B = Q + dna(200) + P;
        ts::UnitigGraph g = graphOf({{A, 50}, {B, 50}});
        link(g, 0, 1, 1, 0);
        link(g, 1, 1, 0, 0);
        const std::string s = A + B.substr(K1) + A.substr(K1);
        const ts::p2::ClosureResult r = ts::p2::closeWalk(s, {fwd(0), fwd(1), fwd(0)}, g);
        check(r.mode == ts::p2::Closure::RepeatEnd && r.closed && r.seq.size() == A.size() + B.size() - 2 * K1 &&
                  isRotation(r.seq, A + B.substr(K1, B.size() - 2 * K1)),
              "R3: repeat_end drops the repeated segment");
    }
    {   // none, no path, ends changed
        const std::string A = P + dna(300) + Q;
        ts::UnitigGraph g = graphOf({{A, 50}});
        check(ts::p2::closeWalk(A, {fwd(0)}, g).mode == ts::p2::Closure::None, "R3: no link closes the walk -> none");
        check(ts::p2::closeWalk(A, {}, g).mode == ts::p2::Closure::NoPath, "R3: no walk -> nopath");
        std::string polished = A;
        polished[0] = polished[0] == 'A' ? 'C' : 'A';
        check(ts::p2::closeWalk(polished, {fwd(0)}, g).mode == ts::p2::Closure::EndsChanged,
              "R3: record ends that differ from the walk -> ends_changed");
    }
    {   // reverse orientation in the walk
        const std::string A = P + dna(500) + P;
        ts::UnitigGraph g = graphOf({{A, 50}});
        link(g, 0, 1, 0, 0);
        const ts::p2::ClosureResult r = ts::p2::closeWalk(rc(A), {fwd(0) | 1}, g);
        check(r.mode == ts::p2::Closure::Self && r.seq == rc(A).substr(0, A.size() - K1),
              "R3: a reverse-oriented walk closes the same way");
    }
}

void testOverlaps() {
    const std::string core = dna(1000);
    const std::string s = core + core.substr(0, 126);
    check(ts::p2::terminalSelfOverlap(s) == 126, "R3: largest exact terminal self-overlap 20..3000");
    check(ts::p2::kMinusOneSelfOverlap(s, {21, 33, 55, 77, 99, 127}) == 126, "R3: a (k-1) self-overlap of the ladder");
    check(ts::p2::kMinusOneSelfOverlap(s, {21, 33, 55, 77, 99}) == 0, "R3: an overlap that is no rung's k-1 is not one");
    check(ts::p2::terminalSelfOverlap(core) == 0, "R3: no self-overlap in random sequence");
}

// ---- F5 -------------------------------------------------------------------------------------
void testSpikein() {
    const std::string& phix = ts::p2::phix174();
    check(phix.size() == 5386, "F5: the built-in phiX174 is 5,386 bp");
    check(std::fabs(ts::p2::phixKmerFraction(phix) - 1.0) < 1e-12, "F5: phiX174 itself scores 1.0");
    const std::string wrap = phix.substr(3000) + phix.substr(0, 3000);
    check(std::fabs(ts::p2::phixKmerFraction(rc(wrap)) - 1.0) < 1e-12, "F5: any rotation, either strand, scores 1.0");
    check(ts::p2::phixKmerFraction(phix + phix.substr(0, 98)) == 1.0, "F5: phiX with its duplicated k-1 end scores 1.0");
    check(ts::p2::phixKmerFraction(dna(5000)) < 0.01, "F5: random sequence scores ~0");
    const std::string mixed = phix.substr(0, 4000) + dna(1500);
    check(ts::p2::phixKmerFraction(mixed) < 0.9, "F5: a record with 27 % foreign sequence is not a spike-in");
}

// ---- reads for the verification and record-stage tests --------------------------------------
struct Temp {
    std::filesystem::path dir = std::filesystem::temp_directory_path() / ("tesseract-p2emit-" + std::to_string(getpid()));
    Temp() { std::filesystem::create_directories(dir); }
    ~Temp() { std::filesystem::remove_all(dir); }
};
// FR pairs, 100 bp mates, fragments 260..340, sampled uniformly from each molecule.
ts::SequenceStore pairsFrom(const std::vector<std::pair<std::string, bool>>& molecules, size_t perKb) {
    Temp t;
    const auto p1 = t.dir / "r1.fa", p2 = t.dir / "r2.fa";
    std::ofstream a(p1), b(p2);
    size_t id = 0;
    for (const auto& m : molecules) {
        const std::string& s = m.first;
        const bool circ = m.second;
        const std::string ext = circ ? s + s.substr(0, 1000) : s;
        const size_t n = s.size() * perKb / 1000;
        for (size_t i = 0; i < n; ++i) {
            const size_t f = 260 + rng() % 81;
            const size_t maxStart = circ ? s.size() - 1 : s.size() - f;
            const size_t st = rng() % (maxStart + 1);
            std::string frag = ext.substr(st, f);
            if (frag.size() < f) continue;
            if (rng() % 2) frag = rc(frag);
            a << '>' << id << "/1\n" << frag.substr(0, 100) << '\n';
            b << '>' << id << "/2\n" << rc(frag).substr(0, 100) << '\n';
            ++id;
        }
    }
    a.close();
    b.close();
    ts::Library l;
    l.r1 = p1.string();
    l.r2 = p2.string();
    ts::SequenceStore out;
    std::string err;
    const bool ok = out.load({l}, 1, err);
    check(ok, "reads load: " + err);
    return out;
}

void testVerify() {
    const ts::p2::VerifyThresholds th;
    // a true 3 kb circle, a chromosome X + F + Y whose fragment F the pairs called circular
    const std::string M = dna(3000);
    const std::string X = dna(20000), F = dna(2500), Y = dna(20000);
    const ts::SequenceStore reads = pairsFrom({{M, true}, {X + F + Y, false}}, 300);
    std::vector<std::string> records = {M + M.substr(0, K1), X, F + F.substr(0, K1), Y};
    ts::p2::VerifyBatch b;
    b.records = &records;
    b.candRecord = {0, 2};
    b.candCircle = {M.substr(1500) + M.substr(0, 1500), F.substr(1250) + F.substr(0, 1250)};
    ts::p2::VerifyRunStats run;
    const auto js = ts::p2::verifyJunctions(b, reads, 2, run);
    check(js.size() == 2 && js[0].computed && js[1].computed, "R3: verification computed for both candidates");
    std::printf("      true circle: span=%zu clip=%zu ctrl=%zu exit=%.4f near=%zu | fragment: span=%zu clip=%zu "
                "ctrl=%zu exit=%.4f near=%zu | window %d-%d\n",
                js[0].spanReads, js[0].clipReads, js[0].ctrlSpan, js[0].exitRatio, js[0].nearJoin, js[1].spanReads,
                js[1].clipReads, js[1].ctrlSpan, js[1].exitRatio, js[1].nearJoin, run.windowLo, run.windowHi);
    check(ts::p2::passJ(js[0], th) && ts::p2::passX(js[0], th), "R3: the true circle's join is spanned, no exits");
    check(!ts::p2::passJ(js[1], th) && !ts::p2::passX(js[1], th),
          "R3: the chromosome fragment's join is clipped and its mates exit");
    // thread-count invariance of the evidence
    const auto js4 = ts::p2::verifyJunctions(b, reads, 4, run);
    check(js4[0].spanReads == js[0].spanReads && js4[1].clipReads == js[1].clipReads &&
              js4[0].exitsNearJoin == js[0].exitsNearJoin && js4[1].exitsNearJoin == js[1].exitsNearJoin,
          "R3: verification is identical across thread counts");
}

// ---- the record stage -------------------------------------------------------------------------
void testRecordStage() {
    const std::string P = dna(K1);
    const std::string A = P + dna(2940) + P;   // a 2,970-bp circle written with its (k-1) duplicate
    const std::string chr = dna(30000);
    const std::string phixRec = ts::p2::phix174() + ts::p2::phix174().substr(0, K1);
    ts::UnitigGraph g = graphOf({{A, 300}, {chr, 30}, {phixRec, 40}});
    link(g, 0, 1, 0, 0);
    const ts::SequenceStore reads = pairsFrom({{A.substr(0, A.size() - K1), true}, {chr, false}}, 200);

    auto setup = [&](std::vector<ts::EmitPiece>& pieces, std::vector<std::string>& scaf, std::vector<std::string>& tags,
                     std::vector<ts::GfaPath>& paths) {
        scaf = {chr, A, phixRec};
        tags = {"_chr", "_plas_1_circular", "_unk_circular"};
        pieces = ts::splitAtGaps(scaf, tags, {30, 300, 40});
        paths.assign(3, ts::GfaPath());
        paths[0].oriented = {fwd(1)}; paths[0].gaps = {0};
        paths[1].oriented = {fwd(0)}; paths[1].gaps = {0};
        paths[2].oriented = {fwd(2)}; paths[2].gaps = {0};
    };
    for (int mode = 0; mode < 3; ++mode) {
        std::vector<ts::EmitPiece> pieces;
        std::vector<std::string> scaf, tags;
        std::vector<ts::GfaPath> paths;
        setup(pieces, scaf, tags, paths);
        ts::p2::RecordStageIn in;
        in.cfg.circ = mode == 2 ? ts::p2::CircMode::Close : ts::p2::CircMode::Verify;
        in.cfg.spikein = true;
        in.pieces = &pieces;
        in.scaffolds = &scaf;
        in.scaffoldTags = &tags;
        in.scaffoldPaths = &paths;
        in.graph = &g;
        const ts::SequenceStore none;
        in.reads = mode == 1 ? &none : &reads;   // mode 1: no reads, so J cannot pass
        in.ladder = {21, 31};
        in.threads = 2;
        const ts::p2::RecordStageResult r = ts::p2::applyRecordStage(in);
        const std::string closed = A.substr(0, A.size() - K1);
        const std::string label = mode == 0 ? "verify" : mode == 1 ? "verify, no reads" : "close";
        check(pieces[2].tag == "_spikein" && pieces[2].seq == phixRec && tags[2] == "_spikein",
              "F5 (" + label + "): the phiX record is labelled _spikein, bases untouched, never _plas/_circular");
        check(r.st.candidates == 1, "R3 (" + label + "): the spike-in is no circle candidate");
        check(pieces[1].seq == closed && scaf[1] == closed && pieces[1].objEnd == closed.size(),
              "R3 (" + label + "): the record and its scaffold lose the duplicated end");
        if (mode == 1) {
            check(pieces[1].tag == "_plas_1" && tags[1] == "_plas_1" && r.st.writtenLinear == 1,
                  "R3 (verify, no reads): unverified -> written linear, duplicate-free (D9)");
        } else {
            check(pieces[1].tag == "_plas_1_circular" && tags[1] == "_plas_1_circular" && r.st.writtenCircular == 1,
                  "R3 (" + label + "): the circle keeps its _circular claim");
        }
        size_t trims = 0, renames = 0;
        for (const auto& e : r.edits) {
            if (e.operation.rfind("trim_end", 0) == 0 && e.bases == K1 && e.start == closed.size()) ++trims;
            if (e.operation.rfind("rename", 0) == 0) ++renames;
        }
        check(trims == 1 && renames == (mode == 1 ? 2u : 1u), "R3/F5 (" + label + "): every edit is logged");
        check(pieces[0].seq == chr && pieces[0].tag == "_chr", "R3/F5 (" + label + "): other records untouched");
        const std::vector<size_t> want = {1, 2};
        check(r.scaffoldsChanged == want, "R3/F5 (" + label + "): the changed scaffolds are listed");
    }
}

void testReportFlagsOff() {
    const ts::p2::EmitConfig off;
    const ts::LibGuardStats lg;
    const std::string j = ts::p2::reportJson(off, lg, nullptr, nullptr);
    // every value is a number, and every number is 0
    bool onlyZeros = true;
    for (size_t i = 0; i < j.size(); ++i) {
        if (j[i] == ':') {
            size_t v = i + 1;
            while (v < j.size() && j[v] == ' ') ++v;
            if (j[v] == '{') continue;
            if (j[v] != '0' || (v + 1 < j.size() && j[v + 1] != ',' && j[v + 1] != '}')) onlyZeros = false;
        }
    }
    check(onlyZeros && j.find('[') == std::string::npos, "report.json p2 block with every flag off: enabled=0, zero counters only");
    check(ts::p2::formatLibGuardCounters(lg).find("[p2-r2] enabled=0") == 0 &&
              ts::p2::formatCircCounters(off, ts::p2::RecordStageStats()).find("[p2-r3] enabled=0") == 0 &&
              ts::p2::formatSpikeinCounters(off, ts::p2::RecordStageStats()).find("[p2-f5] enabled=0") == 0,
          "counter lines print with enabled=0");
}

void testConfig() {
    unsetenv("TESSERACT_P2_LIBGUARD"); unsetenv("TESSERACT_P2_CIRC"); unsetenv("TESSERACT_P2_SPIKEIN");
    check(!ts::p2::readConfig().any(), "config: every emit-A flag is off when unset");
    setenv("TESSERACT_P2_CIRC", "verify", 1);
    setenv("TESSERACT_P2_SPIKEIN", "1", 1);
    setenv("TESSERACT_P2_LIBGUARD", "1", 1);
    const ts::p2::EmitConfig c = ts::p2::readConfig();
    check(c.circ == ts::p2::CircMode::Verify && c.spikein && c.libGuard, "config: verify, spike-in and guard read per call");
    setenv("TESSERACT_P2_CIRC", "close", 1);
    check(ts::p2::readConfig().circ == ts::p2::CircMode::Close, "config: close mode");
    unsetenv("TESSERACT_P2_LIBGUARD"); unsetenv("TESSERACT_P2_CIRC"); unsetenv("TESSERACT_P2_SPIKEIN");
}

}  // namespace

int main() {
    testenv::clearTesseractEnv();
    testLibGuard();
    testClosure();
    testOverlaps();
    testSpikein();
    testVerify();
    testRecordStage();
    testReportFlagsOff();
    testConfig();
    std::printf("test_p2_emit: %d checks, %d failures\n", checks, failures);
    return failures ? 1 : 0;
}
