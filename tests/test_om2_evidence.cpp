// Organism Model 2.0, C1a: the evidence index and the graph side of the seam gate, on
// synthetic graphs built from k-mers exactly as the assembler builds them.
//
// Cases (DESIGN.md C1 unit tests), each with the truth known from the genome string:
//   1  an IS copy skipped at a 1-N butt is resized to the IS walk length
//   2  an island skip whose walk needs an already-placed single-copy node is broken
//   3  a hairpin walk (inverted repeat) is flagged
//   4  tandem operons with a 210 bp spacer give two walk classes
//   5  dead ends on both sides give SILENT
//   6  CONTRA: no walk, but the flank's first single-copy continuation is another piece
//   7  a multi-copy plasmid component is not called a repeat
//   8  a tangled graph makes the index abstain
//   9  multi-copy rRNA operons with inter-copy variants (a 5'-end SNP, ITS length types):
//      every allele combination is a walk, ITS types are separate length classes, the
//      correct join passes, a butt is resized to the walk, the exit is not identifiable
#include <algorithm>
#include <cstdio>
#include <string>
#include <vector>

#include "om2_fixture.h"

using namespace ts;
using namespace om2fx;

namespace {

const int K = 31;

std::unique_ptr<om2::SeamContext> ctxFor(const UnitigGraph& g, om2::SeamConfig::Mode mode,
                                         const om2::EvidenceOptions& opt = om2::EvidenceOptions()) {
    auto ev = om2::EvidenceIndex::build(g, opt);
    return std::unique_ptr<om2::SeamContext>(
        new om2::SeamContext(cfgMode(mode), std::move(ev), nullptr, nullptr));
}

void caseIsResize() {
    std::printf("case 1: IS skipped at a 1-N butt -> RESIZE to the walk\n");
    Rng r(11);
    const std::string A = r.seq(3000), B = r.seq(3000), C = r.seq(3000), IS = r.seq(1200);
    const std::string G = A + IS + B + IS + C;
    UnitigGraph g = buildGraph({{G, 30, false}}, K);
    const std::string X = orientedNodeWith(g, A.substr(1000, 100));
    const std::string Y = orientedNodeWith(g, B.substr(1000, 100));
    check(!X.empty() && !Y.empty(), "flank unitigs found");
    const long truth = trueGap(G, X, Y);
    auto ctx = ctxFor(g, om2::SeamConfig::Mode::Act);
    ctx->beginBatch({X, Y});
    const om2::Decision d = judgeJoin(*ctx, X, Y, 1, om2::Source::JoinButt1);
    const om2::Junction& j = ctx->ledger().j.back();
    check(j.verdict == om2::Verdict::Resize, std::string("verdict RESIZE (got ") + om2::verdictName(j.verdict) + ")");
    check(j.gmin == truth, "gmin " + std::to_string(j.gmin) + " == true gap " + std::to_string(truth));
    check(d.keep && d.writeN == truth, "act: kept, N sized to the IS walk (" + std::to_string(d.writeN) + ")");
    check(j.tier == om2::Tier::C, "tier C (graph)");
    const auto* ev = ctx->evidence();
    uint32_t isId = UINT32_MAX;
    orientedNodeWith(g, IS.substr(500, 100), &isId);
    check(isId != UINT32_MAX && !ev->mult(isId).hardSingle && ev->mult(isId).lo <= 2 && ev->mult(isId).hi >= 2,
          "IS unitig is a repeat with copy interval holding 2");
}

void caseIslandBreak() {
    std::printf("case 2: island skip through a placed single-copy node -> BREAK_DOUBLE_USE\n");
    Rng r(12);
    const std::string A = r.seq(3000), B = r.seq(3000), C = r.seq(3000), D = r.seq(3000),
                      E = r.seq(3000), F = r.seq(3000), I = r.seq(2000), R1 = r.seq(600), R2 = r.seq(600);
    const std::string G = A + R1 + I + R2 + B + C + R1 + D + E + R2 + F;
    UnitigGraph g = buildGraph({{G, 30, false}}, K);
    const std::string X = orientedNodeWith(g, A.substr(1000, 100));
    const std::string Y = orientedNodeWith(g, B.substr(1000, 100));
    const std::string Inode = orientedNodeWith(g, I.substr(900, 100));
    check(!X.empty() && !Y.empty() && !Inode.empty(), "unitigs found");
    {
        auto ctx = ctxFor(g, om2::SeamConfig::Mode::Act);
        ctx->beginBatch({X, Y, Inode});
        const om2::Decision d = judgeJoin(*ctx, X, Y, 100);
        const om2::Junction& j = ctx->ledger().j.back();
        check(j.verdict == om2::Verdict::BreakDoubleUse,
              std::string("island placed: BREAK_DOUBLE_USE (got ") + om2::verdictName(j.verdict) + ")");
        check(!d.keep, "act: the join is broken");
        check(ctx->audit().back().innerUniquePlaced == 1, "one placed single-copy node on the shortest walk");
    }
    {
        auto ctx = ctxFor(g, om2::SeamConfig::Mode::Act);
        ctx->beginBatch({X, Y});   // the island is in no piece
        const om2::Decision d = judgeJoin(*ctx, X, Y, 100);
        const om2::Junction& j = ctx->ledger().j.back();
        check(j.verdict == om2::Verdict::Resize,
              std::string("island unplaced: RESIZE (got ") + om2::verdictName(j.verdict) + ")");
        check(d.keep && d.writeN == trueGap(G, X, Y), "act: sized to the island walk");
    }
    {
        auto ctx = ctxFor(g, om2::SeamConfig::Mode::Audit);
        ctx->beginBatch({X, Y, Inode});
        const om2::Decision d = judgeJoin(*ctx, X, Y, 100);
        check(d.keep && d.writeN == INT32_MIN && ctx->audit().back().action == "would_break",
              "audit: nothing changes, action recorded as would_break");
    }
}

void caseHairpin() {
    std::printf("case 3: inverted repeat -> hairpin flag\n");
    Rng r(13);
    const std::string A = r.seq(3000), B = r.seq(3000), M = r.seq(1500), R = r.seq(800);
    const std::string G = A + R + M + reverseComplement(R) + B;
    UnitigGraph g = buildGraph({{G, 30, false}}, K);
    const std::string X = orientedNodeWith(g, A.substr(1000, 100));
    const std::string Y = orientedNodeWith(g, B.substr(1000, 100));
    auto ev = om2::EvidenceIndex::build(g, om2::EvidenceOptions());
    const om2::EndAnchor a = ev->anchorExit(X), b = ev->anchorEntry(Y);
    check(a.ok && b.ok, "both ends anchored");
    const om2::WalkSet ws = ev->walks(a, b);
    check(ws.reachable && ws.hairpin, "walk exists and uses a unitig in both orientations");
    check(ws.gmin == trueGap(G, X, Y), "gmin equals the true gap");
    // A plain unique-to-unique adjacency is not a hairpin.
    const std::string G2 = r.seq(3000) + R + r.seq(3000) + R + r.seq(3000);
    UnitigGraph g2 = buildGraph({{G2, 30, false}}, K);
    auto ev2 = om2::EvidenceIndex::build(g2, om2::EvidenceOptions());
    const std::string X2 = orientedNodeWith(g2, G2.substr(1000, 100));
    const std::string Y2 = orientedNodeWith(g2, G2.substr(3800 + 1000, 100));
    const om2::WalkSet ws2 = ev2->walks(ev2->anchorExit(X2), ev2->anchorEntry(Y2));
    check(ws2.reachable && !ws2.hairpin, "direct repeat: no hairpin flag");
}

void caseTandem() {
    std::printf("case 4: tandem operons, 210 bp spacer -> two walk classes\n");
    Rng r(14);
    // Flanks long enough that single-copy sequence dominates, as in a real genome (the index
    // abstains below 80% single-copy bases).
    // C is long enough that a walk back round through B+C exceeds the 30 kb enumeration bound.
    const std::string A = r.seq(8000), B = r.seq(8000), C = r.seq(25000), D = r.seq(8000),
                      OP = r.seq(5000), S = r.seq(210);
    const std::string G = A + OP + S + OP + B + C + OP + D;
    UnitigGraph g = buildGraph({{G, 30, false}}, K);
    const std::string X = orientedNodeWith(g, A.substr(1000, 100));
    const std::string Y = orientedNodeWith(g, B.substr(1000, 100));
    auto ev = om2::EvidenceIndex::build(g, om2::EvidenceOptions());
    check(!ev->abstained(), "index built (single-copy fraction " + std::to_string(ev->scFraction) + ")");
    const om2::WalkSet ws = ev->walks(ev->anchorExit(X), ev->anchorEntry(Y));
    check(ws.classes.size() == 2, "two walk-length classes (got " + std::to_string(ws.classes.size()) + ")");
    if (ws.classes.size() == 2) {
        check(ws.classes[1].len - ws.classes[0].len == 5210, "classes differ by operon + spacer (5210)");
        check(ws.classes[1].len == trueGap(G, X, Y), "the longer class is the true tandem gap");
    }
    check(ws.repeatPruned && ws.distinct > 2,
          "longer tandem walks (> 12 kb of repeat) seen by the census, left out of the classes");
}

void caseSilent() {
    std::printf("case 5: dead ends on both sides -> SILENT\n");
    Rng r(15);
    const std::string A = r.seq(4000), B = r.seq(4000);
    UnitigGraph g = buildGraph({{A, 30, false}, {B, 30, false}}, K);
    auto ctx = ctxFor(g, om2::SeamConfig::Mode::Act);
    ctx->beginBatch({A, B});
    const om2::Decision d = judgeJoin(*ctx, A, B, 100);
    const om2::Junction& j = ctx->ledger().j.back();
    check(j.verdict == om2::Verdict::Silent, std::string("SILENT (got ") + om2::verdictName(j.verdict) + ")");
    check(j.endA == om2::EndClass::Unique && j.endB == om2::EndClass::Unique, "both ends single-copy");
    check(d.keep && j.admit == om2::Admit::Scaffold && j.tier == om2::Tier::E,
          "kept, tier E, admitted to scaffolds at unique ends (v0)");
    check(!j.allowCloseGaps, "a SILENT model gap is not handed to the gap filler");
}

void caseContra() {
    std::printf("case 6: first single-copy continuation is another piece -> CONTRA (no pairs: refused)\n");
    Rng r(16);
    const std::string A = r.seq(3000), C = r.seq(3000), D = r.seq(3000), E = r.seq(3000),
                      R = r.seq(700), B = r.seq(3000);
    const std::string G = A + R + C + D + R + E;
    UnitigGraph g = buildGraph({{G, 30, false}, {B, 30, false}}, K);
    const std::string X = orientedNodeWith(g, A.substr(1000, 100));
    auto ctx = ctxFor(g, om2::SeamConfig::Mode::Act);
    ctx->beginBatch({X, B});
    const om2::Decision d = judgeJoin(*ctx, X, B, 100);
    const om2::Junction& j = ctx->ledger().j.back();
    check(j.verdict == om2::Verdict::ContraRefused,
          std::string("CONTRA without pair support -> refused (got ") + om2::verdictName(j.verdict) + ")");
    check(ctx->audit().back().uOther >= 1, "the other continuation is counted");
    check(!d.keep, "act: broken");
}

void casePlasmid() {
    std::printf("case 7: multi-copy plasmid component is not a repeat\n");
    Rng r(17);
    const std::string A = r.seq(20000), B = r.seq(20000), C = r.seq(20000), R = r.seq(1500),
                      P = r.seq(25000);
    const std::string G = A + R + B + R + C;
    UnitigGraph g = buildGraph({{G, 30, false}, {P, 150, true}}, K);
    uint32_t pid = UINT32_MAX, rid = UINT32_MAX, aid = UINT32_MAX;
    orientedNodeWith(g, P.substr(5000, 100), &pid);
    orientedNodeWith(g, R.substr(500, 100), &rid);
    orientedNodeWith(g, A.substr(5000, 100), &aid);
    auto ev = om2::EvidenceIndex::build(g, om2::EvidenceOptions());
    check(pid != UINT32_MAX && ev->mult(pid).hardSingle, "plasmid unitig (5x depth, own component) is single-copy");
    check(rid != UINT32_MAX && !ev->mult(rid).hardSingle, "chromosomal 2-copy repeat is a repeat");
    check(aid != UINT32_MAX && ev->mult(aid).hardSingle, "chromosomal unique unitig is single-copy");
    om2::EvidenceOptions glob;
    glob.perComponentTheta = false;
    auto evg = om2::EvidenceIndex::build(g, glob);
    check(!evg->mult(pid).hardSingle, "control: a global theta would call the plasmid a repeat");
}

void caseTangle() {
    std::printf("case 8: tangled graph -> abstain\n");
    Rng r(18);
    const std::string A = r.seq(3000), B = r.seq(3000), IS = r.seq(1200);
    UnitigGraph g = buildGraph({{A + IS + B + IS + A.substr(0, 500), 30, false}}, K);
    om2::EvidenceOptions opt;
    opt.tangleMax = 1;
    auto ctx = ctxFor(g, om2::SeamConfig::Mode::Act, opt);
    check(ctx->evidence()->abstained(), "index abstains above the node limit");
    const om2::Decision d = judgeJoin(*ctx, A, B, 100);
    check(ctx->ledger().j.back().verdict == om2::Verdict::Abstain && d.keep,
          "the junction is judged ABSTAIN and kept as today");
}

void caseRrna() {
    std::printf("case 9: 4 rRNA operon copies with inter-copy variants\n");
    Rng r(19);
    const std::string OP = r.seq(6500);
    const std::string itsX = OP.substr(2500, 300);
    const std::string itsY = r.seq(450);
    // allele b: a substitution 100 bp into the operon (within pair reach of the flank);
    // ITS type X (300 bp, the base sequence) or Y (450 bp): copies aX, aY, bX, bY.
    auto copyOf = [&](bool b, bool y) {
        std::string s = OP.substr(0, 2500) + (y ? itsY : itsX) + OP.substr(2800);
        if (b) s = Rng::mutateAt(s, {100});
        return s;
    };
    std::vector<std::string> U, Dn;
    for (int i = 0; i < 4; ++i) { U.push_back(r.seq(6000)); Dn.push_back(r.seq(6000)); }
    const std::string G = U[0] + copyOf(false, false) + Dn[0] + U[1] + copyOf(false, true) + Dn[1] +
                          U[2] + copyOf(true, false) + Dn[2] + U[3] + copyOf(true, true) + Dn[3];
    UnitigGraph g = buildGraph({{G, 30, false}}, K);
    std::vector<std::string> X(4), Y(4);
    for (int i = 0; i < 4; ++i) {
        X[i] = orientedNodeWith(g, U[i].substr(1000, 100));
        Y[i] = orientedNodeWith(g, Dn[i].substr(1000, 100));
    }
    bool found = true;
    for (int i = 0; i < 4; ++i) found = found && !X[i].empty() && !Y[i].empty();
    check(found, "all 8 flank unitigs found");
    if (!found) return;
    auto ev = om2::EvidenceIndex::build(g, om2::EvidenceOptions());
    check(!ev->abstained(), "index built (single-copy fraction " + std::to_string(ev->scFraction) + ")");
    uint32_t sharedId = UINT32_MAX, bubbleId = UINT32_MAX;
    orientedNodeWith(g, OP.substr(4000, 100), &sharedId);
    orientedNodeWith(g, itsY.substr(200, 60), &bubbleId);
    check(sharedId != UINT32_MAX && ev->mult(sharedId).lo <= 4 && ev->mult(sharedId).hi >= 4 &&
          !ev->mult(sharedId).hardSingle, "shared operon segment: copy interval holds 4");
    check(bubbleId != UINT32_MAX && ev->mult(bubbleId).lo <= 2 && ev->mult(bubbleId).hi >= 2,
          "ITS-Y branch: copy interval holds 2");

    const om2::EndAnchor a0 = ev->anchorExit(X[0]), b0 = ev->anchorEntry(Y[0]);
    const om2::WalkSet ws = ev->walks(a0, b0);
    std::string cl;
    for (const auto& c : ws.classes) cl += " " + std::to_string(c.len) + "x" + std::to_string(c.n);
    std::printf("    walks=%u distinct=%u capped=%d classes:%s\n", ws.walks, ws.distinct, ws.capped ? 1 : 0, cl.c_str());
    check(ws.classes.size() == 2 && ws.classes[0].n + ws.classes[1].n == 4,
          "4 walks U1->D1 within one operon traversal, one per allele combination");
    check(ws.distinct > 4 && ws.repeatPruned,
          "walks round other operon copies exist (census) but pass > 12 kb of repeat (not in the classes)");
    check(ws.classes.size() == 2 && ws.classes[0].n == 2 && ws.classes[1].n == 2,
          "two length classes (ITS X, ITS Y), two walks each");
    const long t0 = trueGap(G, X[0], Y[0]), t1 = trueGap(G, X[1], Y[1]);
    check(ws.classes.size() == 2 && ws.classes[0].len == t0 && ws.classes[1].len == t1,
          "class lengths equal the true ITS-X and ITS-Y operon gaps");
    const om2::FirstUnique fu = ev->firstUnique(a0.node);
    check(fu.ids.size() == 4, "first single-copy continuations from U1: the 4 downstream flanks");

    auto ctx = ctxFor(g, om2::SeamConfig::Mode::Act);
    ctx->beginBatch({X[0], X[1], X[2], X[3], Y[0], Y[1], Y[2], Y[3]});
    const om2::Decision dOk = judgeJoin(*ctx, X[0], Y[0], static_cast<int32_t>(t0), om2::Source::Layout);
    const om2::Junction jOk = ctx->ledger().j.back();
    check(jOk.verdict == om2::Verdict::PassWalk && dOk.keep && dOk.writeN == INT32_MIN,
          std::string("correct join at the panel gap: PASS_WALK, kept as written (got ") + om2::verdictName(jOk.verdict) + ")");
    check(jOk.tier == om2::Tier::D, "tier D: panel order with a graph-consistent walk");
    const om2::Decision dButt = judgeJoin(*ctx, X[1], Y[1], 1, om2::Source::LayoutButt1);
    const om2::Junction jButt = ctx->ledger().j.back();
    check(jButt.verdict == om2::Verdict::Resize && dButt.keep && dButt.writeN == static_cast<int32_t>(t0),
          "1-N butt across copy 2: RESIZE to the shortest walk class (" + std::to_string(dButt.writeN) + ")");
    const om2::Decision dX = judgeJoin(*ctx, X[0], Y[2], static_cast<int32_t>(t0), om2::Source::Layout);
    check(ctx->ledger().j.back().verdict == om2::Verdict::PassWalk && dX.keep,
          "which exit follows which entry is not identifiable from the graph: U1->D3 also PASS_WALK");
    check(!ctx->audit().back().pairTestable, "a 6.5 kb operon gap is not pair-testable (N > 100)");
}

}  // namespace

int main() {
    testenv::clearTesseractEnv();
    caseIsResize();
    caseIslandBreak();
    caseHairpin();
    caseTandem();
    caseSilent();
    caseContra();
    casePlasmid();
    caseTangle();
    caseRrna();
    if (failures()) {
        std::printf("test_om2_evidence: %d FAILED\n", failures());
        return 1;
    }
    std::printf("test_om2_evidence: all passed\n");
    return 0;
}
