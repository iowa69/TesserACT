// Component test for the dropout bridge (TESSERACT_DROPOUT_BRIDGE, combo2 PKG-BRIDGE).
// T-B1..T-B12 are the design's cases; T-B13.. cover orientation, zero gap, min_reads,
// modal refusal, low complexity, thread-count determinism and two bridges in one graph.
#include "dropout_bridge.h"
#include "graph.h"
#include "seqio.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <random>
#include <string>
#include <unistd.h>
#include <utility>
#include <vector>
#include "test_env.h"   // build_v3: T40, clear the ambient TESSERACT_* first

namespace {

int checks = 0;
void check(bool ok, const std::string& what) {
    ++checks;
    if (!ok) { std::fprintf(stderr, "FAIL: %s\n", what.c_str()); std::exit(1); }
}

std::mt19937 rng(20260925);
std::string randomSeq(size_t n) {
    std::string s(n, 'A');
    for (char& c : s) c = "ACGT"[rng() % 4];
    return s;
}
std::string rc(const std::string& s) { return ts::reverseComplement(s); }
char otherBase(char c) { return c == 'A' ? 'C' : 'A'; }

std::string directory;
int fixtureNo = 0;
constexpr int K = 61;

struct Fixture {
    ts::UnitigGraph g;
    std::vector<std::pair<std::string, std::string>> pairs;

    Fixture() { g.setK(K); }
    uint32_t node(const std::string& seq, double cov = 30.0) {
        ts::Unitig u;
        u.seq = seq;
        u.coverage = cov;
        g.nodes.push_back(u);
        return static_cast<uint32_t>(g.nodes.size() - 1);
    }
    // A pair whose mate lands nowhere near any window.
    void read(const std::string& r1) { pairs.push_back({r1, randomSeq(200)}); }
    void pair(const std::string& r1, const std::string& r2) { pairs.push_back({r1, r2}); }

    ts::SequenceStore load() {
        const std::string base = directory + "/f" + std::to_string(++fixtureNo);
        {
            std::ofstream a(base + "_1.fa"), b(base + "_2.fa");
            for (size_t i = 0; i < pairs.size(); ++i) {
                a << ">p" << i << "/1\n" << pairs[i].first << '\n';
                b << ">p" << i << "/2\n" << pairs[i].second << '\n';
            }
        }
        ts::Library lib;
        lib.r1 = base + "_1.fa";
        lib.r2 = base + "_2.fa";
        ts::SequenceStore reads;
        std::string error;
        check(reads.load({lib}, 2, error), "load fixture: " + error);
        return reads;
    }
    ts::DropoutBridgeStats run(uint32_t minReads = 2, int threads = 2,
                               const std::function<void(ts::SequenceStore&)>& mutate = nullptr) {
        ts::SequenceStore reads = load();
        if (mutate) mutate(reads);
        return ts::bridgeDropouts(g, reads, threads, minReads, false);
    }
    bool hasLive(const std::string& s) const {
        for (const ts::Unitig& u : g.nodes) {
            if (!u.deleted && (u.seq == s || u.seq == rc(s))) return true;
        }
        return false;
    }
};

bool sameGraph(const ts::UnitigGraph& a, const ts::UnitigGraph& b) {
    if (a.k() != b.k() || a.nodes.size() != b.nodes.size()) return false;
    for (size_t i = 0; i < a.nodes.size(); ++i) {
        const ts::Unitig& x = a.nodes[i];
        const ts::Unitig& y = b.nodes[i];
        if (x.seq != y.seq || x.coverage != y.coverage || x.deleted != y.deleted) return false;
        for (int e = 0; e < 2; ++e) if (!(x.ends[e] == y.ends[e])) return false;
    }
    return true;
}

}  // namespace

int main() {
    testenv::clearTesseractEnv();
    char templ[] = "/tmp/tesseract-dropout-bridge-XXXXXX";
    const char* tmp = mkdtemp(templ);
    check(tmp != nullptr, "temporary directory");
    directory = tmp;

    // Truth: U = G[0,500), a 25 bp dropout G[500,525), V = G[525,1025).
    const std::string G = randomSeq(1200);
    const std::string U = G.substr(0, 500), V = G.substr(525, 500), truth = G.substr(0, 1025);

    // T-B1: 3 reads span a 25 bp gap -> bridged; joined sequence = truth; validate passes.
    {
        Fixture f;
        f.node(U); f.node(V);
        for (int s : {380, 390, 400}) f.read(G.substr(s, 200));
        const auto st = f.run();
        check(st.error.empty(), "T-B1 no error");
        check(st.bridged == 1 && st.overlapBridged == 0 && st.bpAdded == 25, "T-B1 bridged once, 25 bp");
        check(st.eligible == 4 && st.deadEnds == 4 && st.anchors == 120 && st.anchorsDropped == 0,
              "T-B1 4 eligible ends, 120 anchors");
        check(st.nominations == 3 && st.pairsNominated == 1, "T-B1 3 nominations, 1 pair");
        check(st.merged == 2, "T-B1 compaction merged exactly the bridge chain");
        check(f.g.validate().empty(), "T-B1 validate");
        check(f.g.liveCount() == 1 && f.hasLive(truth), "T-B1 joined sequence equals the truth");
        const std::string line = ts::formatDropoutBridgeStats(st);
        const char* keys[] = {"[dropoutbridge] K=", " med=", " deadends=", " eligible=", " anchors=",
                              " anchors_dropped=", " reads_hit=", " pairs_nominated=", " ref_minreads=",
                              " ref_dup=", " ref_modal=", " ref_consensus=", " ref_competing=",
                              " ref_overlap=", " ref_lowcx=", " ref_self=", " bridged=",
                              " overlap_bridged=", " bp_added=", " seconds="};
        size_t at = 0;
        for (const char* key : keys) {
            const size_t p = line.find(key, at);
            check(p != std::string::npos, std::string("T-B1 counter line has ") + key + " in order");
            at = p + 1;
        }
        check(line.find(" bridged=1 ") != std::string::npos, "T-B1 counter line says bridged=1");
    }

    // T-B2: 1 read -> not bridged.
    {
        Fixture f;
        f.node(U); f.node(V);
        f.read(G.substr(390, 200));
        const ts::UnitigGraph before = f.g;
        const auto st = f.run();
        check(st.bridged == 0 && st.refMinReads == 1 && st.nominations == 1, "T-B2 one read refused (minreads)");
        check(sameGraph(before, f.g), "T-B2 graph untouched");
    }

    // T-B3: 2 PCR-duplicate reads (same start, different pairs) -> not bridged.
    {
        Fixture f;
        f.node(U); f.node(V);
        f.read(G.substr(390, 200));
        f.read(G.substr(390, 200));
        const auto st = f.run();
        check(st.bridged == 0 && st.refDup == 1 && st.refMinReads == 0, "T-B3 duplicates refused (dup)");
    }

    // T-B4: competing partners, U->V and U->W at 2 reads each -> neither bridged.
    {
        Fixture f;
        const std::string W = randomSeq(500), gapW = randomSeq(20);
        f.node(U); f.node(V); f.node(W);
        const std::string GW = U + gapW + W;
        for (int s : {380, 400}) f.read(G.substr(s, 200));
        for (int s : {385, 405}) f.read(GW.substr(s, 200));
        const auto st = f.run();
        check(st.bridged == 0 && st.refCompeting == 2 && st.pairsNominated == 2, "T-B4 both competitors refused");
    }

    // T-B5: non-mutual: U->V at 3, X->V at 2 -> refused. X->V at 1 does not compete.
    {
        const std::string X = randomSeq(500), gapX = randomSeq(18);
        const std::string GX = X + gapX + V;
        {
            Fixture f;
            f.node(U); f.node(V); f.node(X);
            for (int s : {380, 390, 400}) f.read(G.substr(s, 200));
            for (int s : {383, 403}) f.read(GX.substr(s, 200));
            const auto st = f.run();
            check(st.bridged == 0 && st.refCompeting == 2, "T-B5 V nominated by U and X: refused");
        }
        {
            Fixture f;
            f.node(U); f.node(V); f.node(X);
            for (int s : {380, 390, 400}) f.read(G.substr(s, 200));
            f.read(GX.substr(393, 200));
            const auto st = f.run();
            check(st.bridged == 1 && st.refCompeting == 0 && st.refMinReads == 1,
                  "T-B5 single-read rival does not block");
            check(f.hasLive(truth) && f.hasLive(X), "T-B5 U+V joined, X untouched");
        }
    }

    // T-B6: g = -40 (flanks overlap by 40): consistent -> bridged at |b| = 2K-42; inconsistent -> refused.
    {
        const std::string V2 = G.substr(460, 500);   // U's last 40 bases == V2's first 40
        const std::string truth2 = G.substr(0, 960);
        {
            Fixture f;
            f.node(U); f.node(V2);
            for (int s : {380, 390, 400}) f.read(G.substr(s, 200));
            const auto st = f.run();
            check(st.bridged == 1 && st.overlapBridged == 1 && st.bpAdded == 0, "T-B6 overlap bridged");
            check(f.g.validate().empty() && f.g.liveCount() == 1 && f.hasLive(truth2), "T-B6 joined = truth");
            size_t merged = 0;
            for (const auto& u : f.g.nodes) if (!u.deleted) merged = u.seq.size();
            const size_t bridgeLen = merged - U.size() - V2.size() + 2 * (K - 1);
            check(bridgeLen == static_cast<size_t>(2 * K - 42), "T-B6 |b| = 2K-42");
        }
        {
            std::string V3 = V2;
            V3[10] = otherBase(V3[10]);                // inside the 40 bp overlap
            Fixture f;
            f.node(U); f.node(V3);
            for (int s : {380, 390, 400}) f.read(G.substr(s, 200));
            const auto st = f.run();
            check(st.bridged == 0 && st.refOverlap == 1, "T-B6 inconsistent overlap refused");
        }
    }

    // T-B7: a window A-mer repeated elsewhere in the graph is dropped (bridge still made).
    {
        std::string Z = randomSeq(500);
        Z.replace(200, 31, U.substr(450, 31));       // one of U's end-1 window A-mers
        // Make the flanking bases differ from U's, or the neighbouring window A-mers repeat too.
        if (Z[199] == U[449]) Z[199] = otherBase(U[449]);
        if (Z[231] == U[481]) Z[231] = otherBase(U[481]);
        Fixture f;
        f.node(U); f.node(V); f.node(Z);
        for (int s : {380, 390, 400}) f.read(G.substr(s, 200));
        const auto st = f.run();
        check(st.anchorsDropped == 1 && st.anchors == 179, "T-B7 repeated anchor dropped");
        check(st.bridged == 1 && f.hasLive(truth), "T-B7 bridge made from the other anchors");
    }

    // T-B8: a flank at 3x the median depth is ineligible.
    {
        Fixture f;
        f.node(U, 90.0); f.node(V);
        for (int i = 0; i < 3; ++i) f.node(randomSeq(400));
        for (int s : {380, 390, 400}) f.read(G.substr(s, 200));
        const auto st = f.run();
        check(st.medianCoverage == 30.0 && st.deadEnds == 10 && st.eligible == 8, "T-B8 U's ends ineligible");
        check(st.bridged == 0 && st.nominations == 0, "T-B8 nothing nominated");
    }

    // T-B9: one-base consensus disagreement with 2 reads -> refused.
    {
        Fixture f;
        f.node(U); f.node(V);
        f.read(G.substr(380, 200));
        std::string r = G.substr(400, 200);
        r[512 - 400] = otherBase(r[512 - 400]);      // gap position 12
        f.read(r);
        const auto st = f.run();
        check(st.bridged == 0 && st.refConsensus == 1, "T-B9 consensus disagreement refused");
    }

    // T-B10: reverse-strand mates of one pair count once.
    {
        Fixture f;
        f.node(U); f.node(V);
        f.pair(G.substr(380, 200), rc(G.substr(400, 200)));
        const auto st = f.run();
        check(st.nominations == 2 && st.pairsNominated == 1, "T-B10 both mates nominate the same end pair");
        check(st.bridged == 0 && st.refMinReads == 1, "T-B10 ... and count as one pair");
        Fixture f2;
        f2.node(U); f2.node(V);
        f2.pair(G.substr(380, 200), rc(G.substr(400, 200)));
        f2.pair(G.substr(390, 200), rc(G.substr(405, 200)));
        const auto st2 = f2.run();
        check(st2.bridged == 1 && f2.hasLive(truth), "T-B10 two such pairs bridge");
    }

    // T-B11: a same-unitig nomination (circularising) is refused and counted.
    {
        const std::string gapC = randomSeq(25);
        const std::string C = U + gapC;              // circular replicon: U, then gapC, then U again
        const std::string CC = C + C;
        Fixture f;
        f.node(U);
        for (int s : {380, 390, 400}) f.read(CC.substr(s, 200));
        const auto st = f.run();
        check(st.bridged == 0 && st.refSelf == 1, "T-B11 self pair refused and counted");
    }

    // T-B12: with nothing to bridge the graph is left bit-for-bit alone (no compaction, no
    // renumbering). The flag-off identity of the binary is the A1 assembly check.
    {
        Fixture f;
        f.node(U); f.node(V);
        f.node(randomSeq(404));                      // an unrelated unitig, never touched
        f.read(randomSeq(200));
        const ts::UnitigGraph before = f.g;
        const auto st = f.run();
        check(st.bridged == 0 && st.merged == 0 && sameGraph(before, f.g), "T-B12 graph untouched");
    }

    // T-B13: orientation. V stored reverse-complemented (its facing end is end 1) and every
    // spanning read on the reverse strand.
    {
        Fixture f;
        f.node(U); f.node(rc(V));
        for (int s : {380, 390, 400}) f.read(rc(G.substr(s, 200)));
        const auto st = f.run();
        check(st.bridged == 1 && st.bpAdded == 25 && f.hasLive(truth) && f.g.validate().empty(),
              "T-B13 reverse-stored flank, reverse-strand reads");
    }

    // T-B14: zero gap (reads join the flanks directly): |b| = 2K-2.
    {
        const std::string V0 = G.substr(500, 500);
        Fixture f;
        f.node(U); f.node(V0);
        for (int s : {380, 390, 400}) f.read(G.substr(s, 200));
        const auto st = f.run();
        check(st.bridged == 1 && st.bpAdded == 0 && st.overlapBridged == 0 && f.hasLive(G.substr(0, 1000)),
              "T-B14 zero gap bridged");
    }

    // T-B15: min_reads = 3.
    {
        Fixture f;
        f.node(U); f.node(V);
        for (int s : {380, 400}) f.read(G.substr(s, 200));
        const auto st = f.run(3);
        check(st.minReads == 3 && st.bridged == 0 && st.refMinReads == 1, "T-B15 2 reads < min_reads 3");
        Fixture f2;
        f2.node(U); f2.node(V);
        for (int s : {380, 390, 400}) f2.read(G.substr(s, 200));
        check(f2.run(3).bridged == 1, "T-B15 3 reads >= min_reads 3");
    }

    // T-B16: no dominant gap length (2 reads say 25, 2 say 26) -> refused (modal).
    {
        Fixture f;
        f.node(U); f.node(V);
        const std::string Gi = G.substr(0, 512) + "A" + G.substr(512);   // one-base insertion in the gap
        for (int s : {380, 400}) f.read(G.substr(s, 200));
        for (int s : {385, 405}) f.read(Gi.substr(s, 200));
        const auto st = f.run();
        check(st.bridged == 0 && st.refModal == 1, "T-B16 split gap lengths refused");
    }

    // T-B17: a low-complexity junction (windows and gap nearly all A) is refused, although
    // every anchor is unique.
    {
        const std::string Ulc = randomSeq(440) + std::string(30, 'A') + "C" + std::string(29, 'A');
        const std::string Vlc = std::string(29, 'A') + "G" + std::string(30, 'A') + randomSeq(440);
        const std::string Glc = Ulc + std::string(25, 'A') + Vlc;
        Fixture f;
        f.node(Ulc); f.node(Vlc);
        for (int s : {380, 390, 400}) f.read(Glc.substr(s, 200));
        const auto st = f.run();
        check(st.bridged == 0 && st.refLowComplexity == 1, "T-B17 low-complexity junction refused");
    }

    // T-B18: the result does not depend on the thread count.
    {
        Fixture a, b;
        const std::string R = randomSeq(300);
        for (Fixture* f : {&a, &b}) {
            f->node(U); f->node(V); f->node(R);
        }
        rng.seed(7);
        for (int s : {380, 390, 400}) a.read(G.substr(s, 200));
        rng.seed(7);
        for (int s : {380, 390, 400}) b.read(G.substr(s, 200));
        const auto sa = a.run(2, 1), sb = b.run(2, 5);
        check(sa.bridged == 1 && sb.bridged == 1 && sameGraph(a.g, b.g), "T-B18 1 vs 5 threads identical");
    }

    // T-B19: two independent bridges in one graph.
    {
        const std::string H = randomSeq(1100);
        const std::string P = H.substr(0, 450), Q = H.substr(470, 500);   // 20 bp dropout
        Fixture f;
        f.node(Q); f.node(U); f.node(P); f.node(V);
        for (int s : {380, 390, 400}) f.read(G.substr(s, 200));
        for (int s : {330, 345, 360}) f.read(rc(H.substr(s, 200)));
        const auto st = f.run();
        check(st.bridged == 2 && st.merged == 4 && st.bpAdded == 45, "T-B19 two bridges");
        check(f.g.validate().empty() && f.hasLive(truth) && f.hasLive(H.substr(0, 970)), "T-B19 both joined");
    }

    // T-B20: correctReads masks from where correction stalls to the read's END, so a read
    // spanning a dropout has its far flank masked. Anchors are rolled through the mask
    // (raw bases): the bridge is still made, and the counters say the mask was read through.
    {
        Fixture f;
        f.node(U); f.node(V);
        for (int s : {380, 390, 400}) f.read(G.substr(s, 200));
        const auto st = f.run(2, 2, [](ts::SequenceStore& reads) {
            // Reads 0, 2, 4 are the spanning mates; mask from gap position 5 to the end.
            for (size_t r : {0u, 2u, 4u}) {
                const uint32_t from = static_cast<uint32_t>(505 - (380 + 10 * (r / 2)));
                reads.maskRange(r, from, reads.length(r));
            }
        });
        check(st.bridged == 1 && f.hasLive(truth), "T-B20 far flank masked: still bridged, sequence exact");
        check(st.maskedHits > 0 && st.rawGapVotes == 3 * 20, "T-B20 masked hits and raw gap votes counted");
    }

    std::printf("test_dropout_bridge: %d checks passed\n", checks);
    return 0;
}
