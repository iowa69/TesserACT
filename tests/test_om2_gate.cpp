// Organism Model 2.0, C1b-e: the seam gate's pair test, overlap rule, act-mode policies,
// plasmid rule, closeGaps allow-mask, ledger TSV and the flag reader.
//
// Pair evidence is exercised end to end: read pairs are simulated from a known genome,
// written as FASTQ, loaded into a SequenceStore and anchored by PairIndex, and the graph is
// built independently -- so the two sources can be made to agree or to contradict.
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <unistd.h>
#include <vector>

#include "om2_fixture.h"
#include "seqio.h"

using namespace ts;
using namespace om2fx;

namespace {

const int K = 31;

// Pairs from `genome` at `depth`, fragments N(400, 40) clipped to [320, 480], 100 bp reads.
std::string writeReads(const std::string& dir, const std::string& genome, double depth, unsigned seed,
                       Library& lib) {
    std::mt19937 g(seed);
    std::normal_distribution<double> fl(400.0, 40.0);
    const int rl = 100;
    const size_t pairs = static_cast<size_t>(depth * static_cast<double>(genome.size()) / (2.0 * rl));
    lib.r1 = dir + "/r1.fq";
    lib.r2 = dir + "/r2.fq";
    std::ofstream f1(lib.r1), f2(lib.r2);
    const std::string q(rl, 'I');
    for (size_t i = 0; i < pairs; ++i) {
        int L = static_cast<int>(fl(g));
        L = std::max(320, std::min(480, L));
        const size_t p = g() % (genome.size() - static_cast<size_t>(L));
        const std::string frag = genome.substr(p, static_cast<size_t>(L));
        f1 << "@p" << i << "/1\n" << frag.substr(0, rl) << "\n+\n" << q << "\n";
        f2 << "@p" << i << "/2\n" << reverseComplement(frag.substr(frag.size() - rl)) << "\n+\n" << q << "\n";
    }
    return lib.r1;
}

InsertModel insert400() {
    InsertModel m;
    m.mean = 400;
    m.stddev = 40;
    m.minPlausible = 200;
    m.maxPlausible = 600;
    m.observations = 10000;
    m.usable = true;
    return m;
}

std::string tmpDir() {
    const char* t = std::getenv("TMPDIR");
    std::string tpl = std::string(t && *t ? t : "/tmp") + "/om2gate.XXXXXX";
    std::vector<char> buf(tpl.begin(), tpl.end());
    buf.push_back('\0');
    const char* d = mkdtemp(buf.data());
    return d ? std::string(d) : std::string();
}

struct PairFixture {
    SequenceStore reads;
    bool ok = false;
    PairFixture(const std::string& genome, unsigned seed) {
        const std::string d = tmpDir();
        Library lib;
        writeReads(d, genome, 30.0, seed, lib);
        QualityTrim qt;
        qt.enabled = false;
        reads.setQualityTrim(qt);
        std::string err;
        ok = reads.load({lib}, 2, err);
        if (!ok) std::printf("  read load failed: %s\n", err.c_str());
        std::remove(lib.r1.c_str());
        std::remove(lib.r2.c_str());
        rmdir(d.c_str());
    }
};

std::unique_ptr<om2::SeamContext> ctxWithPairs(const UnitigGraph& g, om2::SeamConfig::Mode mode,
                                               const std::vector<std::string>& origin,
                                               const SequenceStore* reads) {
    om2::EvidenceOptions opt;
    auto ev = om2::EvidenceIndex::build(g, opt);
    std::unique_ptr<om2::PairIndex> pi;
    if (reads) pi.reset(new om2::PairIndex(origin, *reads, insert400(), 2, opt));
    std::unique_ptr<om2::PortMap> pm(new om2::PortMap(origin));
    auto ctx = std::unique_ptr<om2::SeamContext>(
        new om2::SeamContext(cfgMode(mode), std::move(ev), std::move(pi), std::move(pm)));
    ctx->beginBatch(origin);
    return ctx;
}

void casePairRule() {
    std::printf("pair rule (seam_pairs2 / pilot_summary LR)\n");
    om2::EvidenceOptions o;
    const om2::PairCall a = om2::PairIndex::callFrom(3, 5.0f, o);
    check(a.abstain && a.call == 0, "lambda 5 < 10: abstain");
    const om2::PairCall b = om2::PairIndex::callFrom(18, 20.0f, o);
    check(!b.abstain && b.call == -1, "k 18 at lambda 20: adjacency confirmed");
    const om2::PairCall c = om2::PairIndex::callFrom(0, 20.0f, o);
    check(!c.abstain && c.call == +1 && c.log10LR >= 2.0f, "k 0 at lambda 20: insertion (LR >= 100)");
    const om2::PairCall d = om2::PairIndex::callFrom(0, 9.9f, o);
    check(d.abstain, "k 0 at lambda 9.9: abstain, silence is not a veto");
}

void casePairsLive() {
    std::printf("pairs from reads: confirm a SILENT seam, deny an insertion, rescue a CONTRA\n");
    Rng r(21);
    const std::string A = r.seq(6000), B = r.seq(6000), IS = r.seq(1500);
    UnitigGraph g = buildGraph({{A, 30, false}, {B, 30, false}}, K);   // graph: no walk A->B
    {
        PairFixture pf(A + B, 101);   // the reads say A and B are adjacent
        check(pf.ok && pf.reads.pairCount() > 0, "reads loaded");
        auto ctx = ctxWithPairs(g, om2::SeamConfig::Mode::Act, {A, B}, &pf.reads);
        const om2::Decision d = judgeJoin(*ctx, A, B, 1, om2::Source::JoinButt1);
        const om2::Junction& j = ctx->ledger().j.back();
        const om2::JunctionAudit& x = ctx->audit().back();
        check(x.pairTestable && !x.pairAbstain && x.pairCall == -1,
              "adjacent: pairs confirm (k=" + std::to_string(j.pairK) + " lambda=" + std::to_string(j.pairLambda) + ")");
        check(j.verdict == om2::Verdict::Silent && j.tier == om2::Tier::B, "SILENT in the graph, tier B by pairs");
        check(j.pairK >= 0.4 * j.pairLambda, "crossing count near the control rate");
        (void)d;
    }
    {
        // An insertion the join would delete, where the graph connects the flanks through a
        // short repeat (PASS_WALK): the pairs' k = 0 at lambda >= 10 is a veto.
        const std::string P = r.seq(6000), Q = r.seq(6000), R2 = r.seq(300), W = r.seq(6000), Z = r.seq(6000);
        UnitigGraph gw = buildGraph({{P + R2 + Q + W + R2 + Z, 30, false}}, K);
        const std::string X = orientedNodeWith(gw, P.substr(1000, 100));
        const std::string Y = orientedNodeWith(gw, Q.substr(1000, 100));
        PairFixture pf(X + IS + Y, 102);
        auto ctx = ctxWithPairs(gw, om2::SeamConfig::Mode::Act, {X, Y}, &pf.reads);
        const om2::Decision d = judgeJoin(*ctx, X, Y, 1, om2::Source::JoinButt1);
        const om2::Junction& j = ctx->ledger().j.back();
        check(ctx->audit().back().graphVerdict == om2::Verdict::PassWalk && j.verdict == om2::Verdict::PairsDeny &&
              j.pairK == 0 && j.pairLambda >= 10,
              std::string("graph-connected insertion: PAIRS_DENY with k=0 (got ") + om2::verdictName(j.verdict) + ")");
        check(!d.keep || d.writeN != INT32_MIN, "act: the butt is not kept as written (broken or graph-sized)");
    }
    {
        // The same reads against dead-end flanks (SILENT): no crossing pairs is what a coverage
        // dropout looks like too, so the pairs record a deny but do not veto.
        PairFixture pf(A + IS + B, 104);
        auto ctx = ctxWithPairs(g, om2::SeamConfig::Mode::Audit, {A, B}, &pf.reads);
        judgeJoin(*ctx, A, B, 20);
        const om2::Junction& j = ctx->ledger().j.back();
        check(j.verdict == om2::Verdict::Silent && ctx->audit().back().pairCall == 1,
              "SILENT junction: pair deny recorded (pair_call=+1), verdict stays SILENT");
    }
    {
        // CONTRA: A's first single-copy continuation is C in the graph, yet the reads join A to B.
        const std::string C = r.seq(6000), D = r.seq(6000), E = r.seq(6000), R = r.seq(700);
        UnitigGraph g2 = buildGraph({{A + R + C + D + R + E, 30, false}, {B, 30, false}}, K);
        const std::string X = orientedNodeWith(g2, A.substr(1000, 100));
        PairFixture pf(X + B, 103);
        auto ctx = ctxWithPairs(g2, om2::SeamConfig::Mode::Act, {X, B}, &pf.reads);
        const om2::Decision d = judgeJoin(*ctx, X, B, 1, om2::Source::JoinButt1);
        const om2::Junction& j = ctx->ledger().j.back();
        check(j.verdict == om2::Verdict::ContraPairsOk,
              std::string("CONTRA with pair support: kept (got ") + om2::verdictName(j.verdict) + ")");
        check(d.keep, "act: kept (a pair-confirmed butt survives UNSIZED=break)");
        auto ctx2 = ctxWithPairs(g2, om2::SeamConfig::Mode::Act, {X, B}, nullptr);
        const om2::Decision d2 = judgeJoin(*ctx2, X, B, 1, om2::Source::JoinButt1);
        check(ctx2->ledger().j.back().verdict == om2::Verdict::ContraRefused && !d2.keep,
              "CONTRA without pairs: refused");
    }
}

void caseOverlap() {
    std::printf("overlap rule: merge only through single-copy sequence\n");
    Rng r(22);
    const std::string A = r.seq(10000), B = r.seq(10000), C = r.seq(10000), IS = r.seq(1200);
    const std::string G = A + IS + B + IS + C;
    UnitigGraph g = buildGraph({{G, 30, false}}, K);
    auto ev = om2::EvidenceIndex::build(g, om2::EvidenceOptions());
    auto ctx = std::unique_ptr<om2::SeamContext>(
        new om2::SeamContext(cfgMode(om2::SeamConfig::Mode::Act), std::move(ev), nullptr, nullptr));
    // A+IS merged onto IS+C on their shared IS: the collapsed-repeat merge.
    const std::string left = A + IS, right = IS + C;
    ctx->beginBatch({left, right});
    const om2::Decision d = judgeJoin(*ctx, left, right, -1200, om2::Source::JoinOverlap);
    check(!ctx->audit().back().overlapSingle && !d.keep,
          std::string("overlap on an IS copy: refused (verdict ") + om2::verdictName(ctx->ledger().j.back().verdict) + ")");
    // A+IS+B[:1300] merged onto B+IS+C on B's first 1300 bases: single-copy overlap.
    const std::string left2 = A + IS + B.substr(0, 1300), right2 = B + IS + C;
    ctx->beginBatch({left2, right2});
    const om2::Decision d2 = judgeJoin(*ctx, left2, right2, -1300, om2::Source::JoinOverlap);
    check(ctx->audit().back().overlapSingle && d2.keep && ctx->ledger().j.back().tier == om2::Tier::A,
          "overlap inside a single-copy unitig: kept, tier A");
    // Two single-copy unitigs that are graph neighbours overlap by exactly k-1 bases.
    const std::string P = r.seq(8000), Q = r.seq(8000), Z = r.seq(8000);
    UnitigGraph g2 = buildGraph({{P + Q, 30, false}, {P.substr(P.size() - 40) + Z, 30, false}}, K);
    const std::string Xn = orientedNodeWith(g2, P.substr(1000, 100));
    const std::string Yn = orientedNodeWith(g2, Q.substr(1000, 100));
    auto ev2 = om2::EvidenceIndex::build(g2, om2::EvidenceOptions());
    auto ctx2 = std::unique_ptr<om2::SeamContext>(
        new om2::SeamContext(cfgMode(om2::SeamConfig::Mode::Act), std::move(ev2), nullptr, nullptr));
    ctx2->beginBatch({Xn, Yn});
    const om2::Decision d3 = judgeJoin(*ctx2, Xn, Yn, -(K - 1), om2::Source::JoinOverlap);
    check(d3.keep && ctx2->ledger().j.back().gmin == -(K - 1),
          "k-1 overlap of single-copy graph neighbours: kept (gmin " + std::to_string(ctx2->ledger().j.back().gmin) + ")");
}

void casePolicies() {
    std::printf("act-mode policies: UNSIZED, CAP, plasmid rule; audit changes nothing\n");
    Rng r(23);
    const std::string A = r.seq(6000), B = r.seq(6000);
    UnitigGraph gs = buildGraph({{A, 30, false}, {B, 30, false}}, K);   // SILENT pair
    const std::string P = r.seq(6000), Bq = r.seq(16000), Dq = r.seq(16000), E = r.seq(6000), R = r.seq(600);
    const std::string G = P + R + Bq + Dq + R + E;                          // one short repeat
    UnitigGraph gr = buildGraph({{G, 30, false}}, K);
    const std::string X = orientedNodeWith(gr, P.substr(1000, 100));
    const std::string Y = orientedNodeWith(gr, Bq.substr(1000, 100));
    const long t = trueGap(G, X, Y);

    auto mk = [&](const UnitigGraph& g, om2::SeamConfig c) {
        auto ev = om2::EvidenceIndex::build(g, om2::EvidenceOptions());
        auto ctx = std::unique_ptr<om2::SeamContext>(new om2::SeamContext(c, std::move(ev), nullptr, nullptr));
        return ctx;
    };
    om2::SeamConfig c = cfgMode(om2::SeamConfig::Mode::Act);
    c.unsized = om2::SeamConfig::Unsized::Break;
    {
        auto ctx = mk(gs, c);
        ctx->beginBatch({A, B});
        check(!judgeJoin(*ctx, A, B, 1, om2::Source::JoinButt1).keep, "UNSIZED=break: a SILENT butt is broken");
        check(ctx->counters.buttRetired == 1, "butt_retired counted");
    }
    c.unsized = om2::SeamConfig::Unsized::Butt;
    {
        auto ctx = mk(gs, c);
        ctx->beginBatch({A, B});
        const om2::Decision d = judgeJoin(*ctx, A, B, 1, om2::Source::JoinButt1);
        check(d.keep && d.writeN == INT32_MIN, "UNSIZED=butt: the 1-N butt stays as released");
    }
    c.unsized = om2::SeamConfig::Unsized::Sized;
    {
        auto ctx = mk(gr, c);
        ctx->beginBatch({X, Y});
        const om2::Decision d = judgeJoin(*ctx, X, Y, 1, om2::Source::JoinButt1);
        check(d.keep && d.writeN == t, "UNSIZED=sized: butt sized to the graph walk (" + std::to_string(d.writeN) +
              " vs truth " + std::to_string(t) + ")");
        auto ctx2 = mk(gs, c);
        ctx2->beginBatch({A, B});
        check(!judgeJoin(*ctx2, A, B, 1, om2::Source::JoinButt1).keep, "UNSIZED=sized: no walk -> broken");
    }
    c.unsized = om2::SeamConfig::Unsized::Break;
    c.cap = om2::SeamConfig::Cap::Size;
    {
        auto ctx = mk(gr, c);
        ctx->beginBatch({X, Y});
        om2::JudgeRequest rq;
        rq.source = om2::Source::LayoutCap2000;
        rq.left = &X;
        rq.right = &Y;
        rq.claimedN = 2000;
        rq.panelGap = 4000;
        const om2::Decision d = ctx->judge(rq);
        check(d.keep && d.writeN == t, "CAP=size: the 2000-N cap sized from the single walk class (" +
              std::to_string(d.writeN) + ")");
        check(ctx->counters.capRetired == 1, "cap_retired counted");
    }
    c.cap = om2::SeamConfig::Cap::Keep;
    {
        auto ctx = mk(gr, c);
        ctx->beginBatch({X, Y});
        om2::JudgeRequest rq;
        rq.source = om2::Source::LayoutCap2000;
        rq.left = &X;
        rq.right = &Y;
        rq.claimedN = 2000;
        rq.panelGap = 4000;
        const om2::Decision d = ctx->judge(rq);
        check(d.keep && d.writeN == INT32_MIN, "CAP=keep: 2000 N as released");
    }
    {
        auto ctx = mk(gr, cfgMode(om2::SeamConfig::Mode::Act));
        ctx->beginBatch({X, Y});
        om2::JudgeRequest rq;
        rq.source = om2::Source::Join;
        rq.left = &X;
        rq.right = &Y;
        rq.claimedN = static_cast<int32_t>(t);
        rq.plasmidPass = true;
        rq.chromosomalA = true;
        check(!ctx->judge(rq).keep && ctx->counters.plasmidBlocked == 1,
              "plasmid rule: a chromosome-voted end cannot join in the plasmid pass");
        rq.chromosomalA = false;
        check(ctx->judge(rq).keep == (ctx->ledger().j.back().verdict == om2::Verdict::PassExact),
              "plasmid rule: a plasmid join is kept only on positive evidence (PASS_EXACT here)");
    }
    {
        auto ctx = mk(gr, cfgMode(om2::SeamConfig::Mode::Audit));
        ctx->beginBatch({X, Y});
        om2::JudgeRequest rq;
        rq.source = om2::Source::Join;
        rq.left = &X;
        rq.right = &Y;
        rq.claimedN = static_cast<int32_t>(t);
        rq.plasmidPass = true;
        rq.chromosomalA = true;
        const om2::Decision d = ctx->judge(rq);
        check(d.keep && d.writeN == INT32_MIN && ctx->audit().back().action == "would_break",
              "audit: the plasmid rule is recorded, not applied");
    }
}

void caseCloseGapsAndTsv() {
    std::printf("closeGaps allow-mask, closure tracking and the ledger TSV\n");
    Rng r(24);
    const std::string A = r.seq(6000), B = r.seq(6000);
    const std::string P = r.seq(6000), Bq = r.seq(16000), Dq = r.seq(16000), E = r.seq(6000), R = r.seq(600);
    UnitigGraph g = buildGraph({{A, 30, false}, {B, 30, false}, {P + R + Bq + Dq + R + E, 30, false}}, K);
    const std::string X = orientedNodeWith(g, P.substr(1000, 100));
    const std::string Y = orientedNodeWith(g, Bq.substr(1000, 100));
    const std::string G = P + R + Bq + Dq + R + E;
    const int32_t t = static_cast<int32_t>(trueGap(G, X, Y));
    for (int mode = 0; mode < 2; ++mode) {
        auto ev = om2::EvidenceIndex::build(g, om2::EvidenceOptions());
        auto ctx = std::unique_ptr<om2::SeamContext>(new om2::SeamContext(
            cfgMode(mode ? om2::SeamConfig::Mode::Act : om2::SeamConfig::Mode::Audit), std::move(ev), nullptr, nullptr));
        ctx->beginBatch({A, B, X, Y});
        judgeJoin(*ctx, A, B, 100);                    // SILENT: not graph-backed
        judgeJoin(*ctx, X, Y, t, om2::Source::Layout); // graph-backed, unique ends
        const std::vector<om2::Junction>& J = ctx->ledger().j;
        check(!J[0].allowCloseGaps && J[1].allowCloseGaps, "allowCloseGaps: SILENT no, graph-backed yes");
        std::vector<std::string> seqs = {A + std::string(100, 'N') + B,
                                         reverseComplement(X + std::string(static_cast<size_t>(t), 'N') + Y)};
        std::vector<uint8_t> allow;
        const bool apply = ctx->beforeCloseGaps(seqs, allow);
        check(allow.size() == 2 && allow[0] == 0 && allow[1] == 1,
              "mask follows the ledger, including a reverse-complemented scaffold");
        check(apply == (mode == 1), mode ? "act: the mask is applied" : "audit: the mask is only recorded");
        // the gap filler closes the second gap
        std::vector<std::string> after = {seqs[0], reverseComplement(X + G.substr(G.find(X) + X.size(), static_cast<size_t>(t)) + Y)};
        ctx->afterCloseGaps(after, mode ? 1 : 0);
        check(ctx->audit()[1].closedByGapfill && !ctx->audit()[0].closedByGapfill, "closure attributed to junction 1");
        if (mode == 1) {
            char tpl[] = "/tmp/om2tsvXXXXXX";
            const int fd = mkstemp(tpl);
            if (fd >= 0) close(fd);
            std::string err;
            const bool w = ctx->writeTsv(tpl, {reverseComplement(seqs[0])}, {"NODE_1_length_x"}, err);
            std::ifstream in(tpl);
            std::stringstream ss;
            ss << in.rdbuf();
            const std::string txt = ss.str();
            std::remove(tpl);
            size_t lines = 0;
            for (char ch : txt) lines += ch == '\n';
            check(w && lines == 4, "TSV: comment, header and one row per junction");
            check(txt.find("NODE_1_length_x\t6000\t100\t-") != std::string::npos,
                  "TSV: junction 0 located in the written scaffold (reverse strand)");
            check(txt.find("SILENT") != std::string::npos && txt.find("\tE\t") != std::string::npos,
                  "TSV: verdict and tier columns");
        }
    }
}

void caseRepeatBoundedLocate() {
    std::printf("locating repeat-bounded junctions that share their flanks and graph anchors\n");
    Rng r(25);
    const std::string U1 = r.seq(6000), U2 = r.seq(6000), V1 = r.seq(6000), V2 = r.seq(6000),
                      M1 = r.seq(2000), M2 = r.seq(2000), R = r.seq(1500), R2 = r.seq(1500);
    UnitigGraph g = buildGraph({{U1 + R + M1 + R2 + V1 + U2 + R + M2 + R2 + V2, 30, false}}, K);
    auto ev = om2::EvidenceIndex::build(g, om2::EvidenceOptions());
    auto ctx = std::unique_ptr<om2::SeamContext>(
        new om2::SeamContext(cfgMode(om2::SeamConfig::Mode::Act), std::move(ev), nullptr, nullptr));
    // Each piece ends in (starts with) a whole copy of the collapsed repeat: identical 32 bp
    // flanks and identical graph anchors for the two junctions.
    const std::string X1 = U1 + R, Y1 = R2 + V1, X2 = U2 + R, Y2 = R2 + V2;
    ctx->beginBatch({X1, Y1, X2, Y2});
    judgeJoin(*ctx, X1, Y1, 2000, om2::Source::Layout);
    judgeJoin(*ctx, X2, Y2, 2000, om2::Source::Layout);
    const auto& A = ctx->audit();
    check(A[0].anchorSeqA == A[1].anchorSeqA && ctx->ledger().j[0].flankL32 == ctx->ledger().j[1].flankL32,
          "fixture: same flanks and same graph anchor");
    const std::vector<std::string> seqs = {X2 + std::string(2000, 'N') + Y2,
                                           reverseComplement(X1 + std::string(2000, 'N') + Y1)};
    std::vector<uint8_t> allow;
    ctx->beforeCloseGaps(seqs, allow);
    check(ctx->counters.closegapsSeen == 2, "both N-runs mapped to their own junction by the far probes");
    char tpl[] = "/tmp/om2locXXXXXX";
    const int fd = mkstemp(tpl);
    if (fd >= 0) close(fd);
    std::string err;
    ctx->writeTsv(tpl, seqs, {"REC_A", "REC_B"}, err);
    std::ifstream in(tpl);
    std::stringstream ss;
    ss << in.rdbuf();
    std::remove(tpl);
    const std::string txt = ss.str();
    check(txt.find("REC_B\t" + std::to_string(Y1.size()) + "\t2000\t-") != std::string::npos &&
          txt.find("REC_A\t" + std::to_string(X2.size()) + "\t2000\t+") != std::string::npos,
          "TSV: junction 0 in REC_B (reverse), junction 1 in REC_A (forward)");
}

void caseEnv() {
    std::printf("flags: all unset = no context, zero counters\n");
    testenv::clearTesseractEnv();
    om2::SeamConfig c = om2::SeamConfig::fromEnv();
    check(!c.active() && c.mode == om2::SeamConfig::Mode::Off, "no TESSERACT_OM2_*: inactive");
    UnitigGraph empty;
    SequenceStore noReads;
    check(om2::SeamContext::create(empty, {}, noReads, InsertModel(), 1) == nullptr, "create() returns null");
    std::FILE* f = std::tmpfile();
    om2::SeamContext::printCounters(nullptr, f);
    std::rewind(f);
    char buf[4096];
    std::string out;
    while (std::fgets(buf, sizeof buf, f)) out += buf;
    std::fclose(f);
    check(out.find("[om2-evidence] enabled=0 nodes=0") == 0 && out.find("[om2-seam] enabled=0 mode=off") != std::string::npos &&
              out.find(" judged=0 ") != std::string::npos && out.find(" broken=0 ") != std::string::npos,
          "counter lines print with zeros");
    setenv("TESSERACT_OM2_SEAM", "audit", 1);
    c = om2::SeamConfig::fromEnv();
    check(c.active() && c.mode == om2::SeamConfig::Mode::Audit && c.pairs, "SEAM=audit: audit mode, pairs on");
    unsetenv("TESSERACT_OM2_SEAM");
    setenv("TESSERACT_OM2_UNSIZED", "sized", 1);
    c = om2::SeamConfig::fromEnv();
    check(c.evidence && c.mode == om2::SeamConfig::Mode::Off && c.policyWithoutAct,
          "a policy flag without SEAM=act builds evidence only and is reported as ignored");
    unsetenv("TESSERACT_OM2_UNSIZED");
    setenv("TESSERACT_OM2_SEAM", "act", 1);
    setenv("TESSERACT_OM2_CAP", "keep", 1);
    setenv("TESSERACT_OM2_PAIRS", "0", 1);
    c = om2::SeamConfig::fromEnv();
    check(c.act() && c.cap == om2::SeamConfig::Cap::Keep && !c.pairs &&
          c.unsized == om2::SeamConfig::Unsized::Break, "SEAM=act, CAP=keep, PAIRS=0; UNSIZED defaults to break");
    testenv::clearTesseractEnv();
}

}  // namespace

int main() {
    testenv::clearTesseractEnv();
    casePairRule();
    casePairsLive();
    caseOverlap();
    casePolicies();
    caseCloseGapsAndTsv();
    caseRepeatBoundedLocate();
    caseEnv();
    if (failures()) {
        std::printf("test_om2_gate: %d FAILED\n", failures());
        return 1;
    }
    std::printf("test_om2_gate: all passed\n");
    return 0;
}
