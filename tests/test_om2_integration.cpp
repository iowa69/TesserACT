// Organism Model 2.0 integration (build_om2): the hand-over between the three components.
//
//   c1_to_c2   C2 attaches to the rows C1 recorded instead of recording the gaps again. A row C1
//              located reverse-complemented receives the fill reverse-complemented (spans mirrored);
//              a row C1 would break is never bridged and never cemented by closeGaps; a gap C1 did
//              not record gets a new row whose id continues after C1's rows; the repeat-variant
//              rows name the placed allele exactly where it sits in the row's fill.
//   c2_to_c3   the writer takes C1's location of each junction before the exact-flank rule, so two
//              repeat-bounded junctions with identical flanks each get their own fill (one of them
//              recorded reverse-complemented); a wrap closed by an overlap of the record's ends is
//              written with the repeated end once; the origin locator turns a record whose dnaA is
//              on the reverse strand and rotates it to start there; a wrap C2 left open (writtenN
//              -1) is not called circular.
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <vector>

#include "counter.h"
#include "graph.h"
#include "kmer.h"
#include "om2_close.h"
#include "om2_output.h"
#include "resolve.h"
#include "seqio.h"
#include "test_env.h"

using namespace ts;
using namespace ts::om2;

namespace {
int checks = 0;
void check(bool ok, const std::string& why) {
    ++checks;
    if (!ok) throw std::runtime_error(why);
}

std::mt19937 rng(20260929);
std::string dna(size_t n) {
    std::string s(n, 'A');
    for (char& c : s) c = "ACGT"[rng() & 3];
    return s;
}
std::string rc(const std::string& s) { return reverseComplement(s); }
char other(char c) { return c == 'A' ? 'C' : c == 'C' ? 'G' : c == 'G' ? 'T' : 'A'; }

struct Mol { std::string seq; bool circular; };
constexpr int kK = 41;

UnitigGraph buildGraph(const std::vector<Mol>& mols) {
    KmerTable t(1 << 16);
    for (const Mol& m : mols) {
        const std::string s = m.circular ? m.seq + m.seq.substr(0, kK - 1) : m.seq;
        for (size_t i = 0; i + kK <= s.size(); ++i) {
            bool ok = false;
            const Kmer km = stringToKmer(s.substr(i, kK), kK, ok);
            const Kmer c = canonical(km, kK);
            t.put(c, t.get(c) + 1);
        }
    }
    return UnitigGraph::build(t, kK, 1);
}

struct Temp {
    std::filesystem::path dir;
    explicit Temp(const std::string& tag)
        : dir(std::filesystem::temp_directory_path() / ("tesseract-om2int-" + tag + "-" + std::to_string(getpid()))) {
        std::filesystem::create_directories(dir);
    }
    ~Temp() { std::filesystem::remove_all(dir); }
};

SequenceStore simulate(const std::vector<Mol>& mols, const Temp& tmp) {
    std::normal_distribution<double> ins(350, 30);
    std::ofstream a(tmp.dir / "r1.fa"), b(tmp.dir / "r2.fa");
    size_t id = 0;
    for (const Mol& m : mols) {
        const size_t L = m.seq.size();
        const size_t pairs = L * 40 / 300;
        for (size_t p = 0; p < pairs; ++p) {
            const size_t f = static_cast<size_t>(std::max(250.0, std::min(450.0, ins(rng))));
            if (L <= f) continue;
            const size_t start = rng() % (L - f);
            const std::string frag = m.seq.substr(start, f);
            a << '>' << id << "/1\n" << frag.substr(0, 150) << '\n';
            b << '>' << id << "/2\n" << rc(frag.substr(f - 150)) << '\n';
            ++id;
        }
    }
    a.close();
    b.close();
    Library lib;
    lib.r1 = (tmp.dir / "r1.fa").string();
    lib.r2 = (tmp.dir / "r2.fa").string();
    SequenceStore s;
    std::string err;
    check(s.load({lib}, 1, err), "reads load: " + err);
    return s;
}

bool tiles(const Junction& j) {
    uint32_t at = 0;
    for (const FillSpan& s : j.fillSpans) {
        if (s.off != at || s.len == 0) return false;
        at += s.len;
    }
    return at == j.fillSeq.size() && at > 0;
}

// ---- C1 -> C2 ---------------------------------------------------------------------------------
void testC1ToC2() {
    std::vector<std::string> U;
    for (int i = 0; i < 4; ++i) U.push_back(dna(3000));
    const std::string IS = dna(1200);
    std::string IS1 = IS;
    IS1[600] = other(IS1[600]);   // copy 2 differs mid-element: a bubble no pair reaches
    const std::string genome = U[0] + IS + U[1] + IS1 + U[2] + IS + U[3];
    const std::string N(1200, 'N');
    std::vector<std::string> seqs = {U[0] + N + U[1] + N + U[2] + N + U[3]};
    const std::vector<Mol> mols = {{genome, false}};
    UnitigGraph g = buildGraph(mols);
    Temp tmp("c1c2");
    const SequenceStore reads = simulate(mols, tmp);
    InsertModel ins;
    ins.mean = 350;
    ins.stddev = 30;
    ins.usable = true;

    // C1's ledger: row 0 graph-backed and recorded reverse-complemented; row 1 a junction C1 breaks.
    Ledger c1;
    Junction r0;
    r0.id = 0;
    r0.source = Source::Layout;
    r0.claimedN = r0.writtenN = 1200;
    r0.verdict = Verdict::PassWalk;
    r0.tier = Tier::D;
    r0.admit = Admit::Scaffold;
    r0.panelGenomes = 40;
    r0.flankL32 = rc(U[1].substr(0, 32));
    r0.flankR32 = rc(U[0].substr(U[0].size() - 32));
    Junction r1 = r0;
    r1.id = 1;
    r1.verdict = Verdict::BreakDoubleUse;
    r1.admit = Admit::Break;
    r1.flankL32 = U[1].substr(U[1].size() - 32);
    r1.flankR32 = U[2].substr(0, 32);
    c1.j = {r0, r1};
    const std::vector<std::vector<RunOwner>> runs = {{{3000, 4200, 0, '-'}, {7200, 8400, 1, '+'}, {11400, 12600, -1, '.'}}};
    std::vector<RepeatVariant> variants;

    CloseOptions o;
    o.enabled = true;
    o.fill = FillMode::Genome;
    o.alloc = AllocMode::Phased;
    o.allocParams.mode = AllocMode::Phased;
    CloseInputs in{g, reads, ins, "", "", 1, false};
    in.c1 = &c1;
    in.c1Runs = &runs;
    in.variants = &variants;
    Ledger led;
    const std::string before = seqs[0];
    const CloseStats st = runCloseStage(in, o, seqs, led);
    check(seqs[0] == before, "genome-grade fills never touch the record");
    check(st.c1Attached == 2, "two N-runs attach to C1's rows, got " + std::to_string(st.c1Attached));
    check(st.c1Refused == 1, "the row C1 breaks is refused, got " + std::to_string(st.c1Refused));
    check(st.c1Gated == 1, "the graph-backed row is bridged on C1's verdict, got " + std::to_string(st.c1Gated));
    check(c1.j.size() == 2, "C1's ledger keeps its two rows");
    check(c1.j[0].fillSeq == rc(IS), "row 0 receives the IS fill in its own (reverse) orientation");
    check(tiles(c1.j[0]), "row 0's spans tile its fill");
    check(c1.j[0].allowCloseGaps, "a bridged row stays open to closeGaps");
    check(c1.j[0].verdict == Verdict::PassWalk && c1.j[0].tier == Tier::D && c1.j[0].admit == Admit::Scaffold,
          "C1's verdict, tier and admission stand");
    check(c1.j[0].flankL32 == r0.flankL32, "C1's flanks stand");
    check(c1.j[1].fillSeq.empty(), "the row C1 breaks gets no fill");
    check(!c1.j[1].allowCloseGaps, "and closeGaps may not cement it");
    check(led.j.size() == 1, "one new row, for the gap C1 did not record; got " + std::to_string(led.j.size()));
    check(led.j[0].id == 2, "its id continues after C1's rows");
    check(led.j[0].fillSeq == IS, "and its fill is the third copy");
    check(led.j[0].flankL32 == U[2].substr(U[2].size() - 32) && led.j[0].flankR32 == U[3].substr(0, 32),
          "the new row carries the record's flanks");
    // repeat_variants rows: the placed allele is exactly where it sits in the row's fill
    size_t rows0 = 0, rows2 = 0;
    for (const RepeatVariant& v : variants) {
        check(v.junction == 0 || v.junction == 2, "variant rows only for the bridged loci");
        const std::string& fill = v.junction == 0 ? c1.j[0].fillSeq : led.j[0].fillSeq;
        check(static_cast<size_t>(v.fillOffset) + v.siteLen <= fill.size(), "site inside its fill");
        check(fill.substr(v.fillOffset, v.siteLen) == v.placed,
              "placed allele equals the fill at the site (junction " + std::to_string(v.junction) + ")");
        check(v.alleles.size() == 2, "both alleles of the bubble are listed");
        check(v.siteLen == 1 && v.placed.size() == 1 && v.alleles[0].first.size() == 1 && v.alleles[1].first.size() == 1,
              "the SNV is reported as the one differing base, not the whole bubble branch");
        check(v.site == "SNV", "an equal-length bubble is an SNV site");
        check(v.fillOffset == (v.junction == 2 ? 600u : 599u),
              "the site sits at the SNV (IS offset 600; mirrored to 599 in the reversed row), got " +
                  std::to_string(v.fillOffset));
        bool listed = false;
        for (const auto& al : v.alleles) listed = listed || al.first == v.placed;
        check(listed, "the placed allele is one of them");
        check(v.basis != Basis::PairPhased && v.basis != Basis::ThreadPhased, "a mid-IS site is not phased");
        (v.junction == 0 ? rows0 : rows2)++;
    }
    check(rows0 == 1 && rows2 == 1, "one variant row per locus");
    check(st.variantRows == 2, "counted");
    // the true allele at both bridged loci is the majority one (copies 1 and 3 carry it)
    check(c1.j[0].fillSeq.find(rc(IS.substr(580, 41))) != std::string::npos, "locus 1 carries the majority allele");
}

// ---- C2 -> C3 ---------------------------------------------------------------------------------
std::vector<std::pair<std::string, std::string>> readFasta(const std::string& p) {
    std::vector<std::pair<std::string, std::string>> out;
    std::ifstream in(p);
    std::string l;
    while (std::getline(in, l)) {
        if (!l.empty() && l[0] == '>') out.push_back({l.substr(1), ""});
        else if (!out.empty()) out.back().second += l;
    }
    return out;
}

Junction row(uint32_t id, Source src) {
    Junction j;
    j.id = id;
    j.source = src;
    j.verdict = Verdict::PassWalk;
    j.tier = Tier::C;
    j.admit = Admit::Scaffold;
    return j;
}

void testC2ToC3() {
    Temp tmp("c2c3");
    const std::string dir = tmp.dir.string();
    // record 1: two repeat-bounded junctions with identical flanks; the second recorded reversed
    const std::string FL = dna(32), FR = dna(32);
    const std::string P1 = dna(1000) + FL, Q1 = FR + dna(1000) + FL, P2 = FR + dna(1000);
    const std::string N(200, 'N');
    const std::string fill0 = dna(200), fill1 = dna(200);
    // record 2: a plasmid whose two ends overlap by 50 bp, dnaA-like marker on the reverse strand
    const std::string M = dna(40);
    std::string X = dna(3000);
    X.replace(1000, 40, rc(M));
    // record 3: a plasmid whose wrap C2 left open
    const std::string Y = dna(2000);
    std::vector<std::string> seqs = {P1 + N + Q1 + N + P2, X + X.substr(0, 50), Y};
    std::vector<std::string> names = {"NODE_1_length_3264_cov_30.0000_chr",
                                      "NODE_2_length_3050_cov_60.0000_plas_1",
                                      "NODE_3_length_2000_cov_60.0000_plas_2"};
    Ledger L;
    Junction j0 = row(0, Source::Layout);
    j0.claimedN = j0.writtenN = 200;
    j0.flankL32 = FL;
    j0.flankR32 = FR;
    j0.fillSeq = fill0;
    j0.fillSpans = {{0, 200, Basis::GraphWalk}};
    Junction j1 = row(1, Source::Layout);
    j1.claimedN = j1.writtenN = 200;
    j1.flankL32 = rc(FR);
    j1.flankR32 = rc(FL);
    j1.fillSeq = rc(fill1);
    j1.fillSpans = {{0, 150, Basis::GraphWalk}, {150, 50, Basis::Consensus}};
    Junction w0 = row(2, Source::Wrap);
    w0.claimedN = -50;
    w0.writtenN = 0;
    w0.flankL32 = seqs[1].substr(seqs[1].size() - 32);
    w0.flankR32 = seqs[1].substr(0, 32);
    Junction w1 = row(3, Source::Wrap);
    w1.writtenN = -1;
    w1.flankL32 = Y.substr(Y.size() - 32);
    w1.flankR32 = Y.substr(0, 32);
    L.j = {j0, j1, w0, w1};
    const std::vector<std::vector<RunOwner>> owner = {{{1032, 1232, 0, '+'}, {2296, 2496, 1, '-'}}, {}, {}};

    setenv("TESSERACT_OM2_OUTPUT", "1", 1);
    SurfaceInput si;
    si.outDir = dir;
    si.seqs = &seqs;
    si.names = &names;
    si.scaffoldsFile = true;
    si.ledger = &L;
    si.runOwner = &owner;
    si.originLocator = [&M](const std::string& t) -> int64_t {
        const size_t f = t.find(M);
        if (f != std::string::npos) return static_cast<int64_t>(f);
        const size_t r = t.find(rc(M));
        if (r != std::string::npos) return -2 - static_cast<int64_t>(r + M.size() - 1);
        return -1;
    };
    si.stageCounters = {{"om2_close", "[om2-close] enabled=1 bridges=2"}};
    std::string js;
    const SurfaceStats st = writeSurface(si, js);
    check(js.find("\"stage_counters\"") != std::string::npos &&
              js.find("\"om2_close\": \"[om2-close] enabled=1 bridges=2\"") != std::string::npos,
          "the other stages' counter lines are carried into the report.json om2 block");
    unsetenv("TESSERACT_OM2_OUTPUT");
    check(st.ok, "writer ok: " + st.error);
    check(st.ownedByC1 == 2 && st.unrecorded == 0 && st.ambiguous == 0,
          "both identical-flank junctions located by C1's owner map (owned " + std::to_string(st.ownedByC1) + ")");
    const auto G = readFasta(dir + "/genome/genome.fasta");
    check(G.size() == 3, "three genome records, got " + std::to_string(G.size()));
    std::string chr, pl1, pl2, h1, h2;
    for (const auto& r : G) {
        const std::string nm = r.first.substr(0, r.first.find(' '));
        std::string up = r.second;
        for (char& c : up) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
        if (nm == "chromosome_1") chr = up;
        if (nm == "plasmid_1") { pl1 = up; h1 = r.first; }
        if (nm == "plasmid_2") { pl2 = up; h2 = r.first; }
    }
    check(chr == P1 + fill0 + Q1 + fill1 + P2, "each fill at its own junction, the reversed row's turned back");
    check(st.overlapTrimmed == 1, "the overlap wrap is closed by writing the repeated end once");
    check(pl1.size() == X.size(), "plasmid_1 is the circle's length, got " + std::to_string(pl1.size()));
    check(st.reversed == 1 && pl1.compare(0, M.size(), M) == 0, "reverse-strand origin: record turned and starts at it");
    check((pl1 + pl1).find(rc(X)) != std::string::npos, "and is a rotation of the reverse-complemented circle");
    check(h1.find("topology=circular") != std::string::npos, "plasmid_1 is circular");
    check(h2.find("topology=circular") == std::string::npos && pl2 == Y, "the open wrap stays linear and untouched");
    const std::string bed = [&] {
        std::ifstream in(dir + "/genome/genome.mask.bed");
        std::stringstream ss;
        ss << in.rdbuf();
        return ss.str();
    }();
    check(bed.find("CONSENSUS") != std::string::npos, "the reversed row's allocation span is in the mask");
    // the allocated span is lowercase in genome.fasta: the last 50 bases of the reversed fill as
    // stored are the first 50 of fill1 in record orientation
    std::string raw;
    for (const auto& r : G) if (r.first.rfind("chromosome_1", 0) == 0) raw = r.second;
    const size_t at = P1.size() + fill0.size() + Q1.size();
    bool lower = true, upperRest = true;
    for (size_t i = 0; i < 50; ++i) lower = lower && std::islower(static_cast<unsigned char>(raw[at + i]));
    for (size_t i = 50; i < 200; ++i) upperRest = upperRest && std::isupper(static_cast<unsigned char>(raw[at + i]));
    check(lower && upperRest, "allocation lowercase exactly on its span, mirrored for the reversed row");
}

}  // namespace

int main() {
    testenv::clearTesseractEnv();
    try {
        testC1ToC2();
        testC2ToC3();
    } catch (const std::exception& e) {
        std::fprintf(stderr, "FAIL after %d checks: %s\n", checks, e.what());
        return 1;
    }
    std::printf("test_om2_integration: %d checks passed\n", checks);
    return 0;
}
