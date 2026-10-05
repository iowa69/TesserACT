// Component tests for the 1.5 release additions:
//
// 1. The layout-only genome view (TESSERACT_OM2_LAYOUT_ONLY, src/om2_output.cpp). One scaffold
//    record carries eight junctions of every kind the confirmation rule separates:
//      a resolver N-run (read_pairs_scaffolder), a panel-only join (cut), a join that read pairs
//      cross (read_pairs), a layout gap with C1 verdict PASS_WALK (cut under `isolate`, kept under
//      `c1pass`), a gap filled from the graph (graph_fill), a gap filled with a CONSENSUS span
//      (allocation: cut), an exact-overlap clonal merge nothing confirms (cut: abutting pieces),
//      and a seam no ledger row owns (unrecorded: kept).
//    Checked: the pieces of genome.fasta, their names and header fields; genome/layout.agp (valid
//    AGP 2.1, one object, the pieces in order with a gap row at every cut junction of the size the
//    view asserted); scaffolds.fasta cut at the same junctions; junctions.tsv statuses and the
//    extra columns; that no base of a cut fill is written; that the sequence of the pieces is the
//    sequence of the uncut view with the cut junctions removed; the report.json layout block and
//    the scaffold record table handed back to the caller. Then `c1pass`, and the flag unset (no cut,
//    no layout.agp, no extra columns: the round-3c writer).
// 2. The `salmonella` alias of --organism (src/organism.cpp).
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <random>
#include <sstream>
#include <string>
#include <vector>
#include <unistd.h>

#include "test_env.h"
#include "om2_output.h"
#include "organism.h"
#include "seqio.h"

using namespace ts;
using namespace ts::om2;

namespace {

int checks = 0, failures = 0;
void check(bool ok, const std::string& what) {
    ++checks;
    if (!ok) { ++failures; std::fprintf(stderr, "FAIL: %s\n", what.c_str()); }
}

std::mt19937 rng(20260930);
std::string randomSeq(size_t n) {
    std::string s(n, 'A');
    for (char& c : s) c = "ACGT"[rng() % 4];
    return s;
}
std::string slurp(const std::string& p) {
    std::ifstream in(p, std::ios::binary);
    std::ostringstream o;
    o << in.rdbuf();
    return o.str();
}
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
std::string firstToken(const std::string& h) { return h.substr(0, h.find(' ')); }
std::vector<std::string> tabs(const std::string& l) {
    std::vector<std::string> v;
    std::stringstream ss(l);
    std::string x;
    while (std::getline(ss, x, '\t')) v.push_back(x);
    return v;
}
std::vector<std::vector<std::string>> readTable(const std::string& p) {
    std::vector<std::vector<std::string>> rows;
    std::ifstream in(p);
    std::string l;
    while (std::getline(in, l)) if (!l.empty() && l[0] != '#') rows.push_back(tabs(l));
    return rows;
}

Junction blank(uint32_t id, Source src) {
    Junction j{};
    j.id = id;
    j.source = src;
    j.a = {0, true};
    j.b = {0, false};
    j.verdict = Verdict::Silent;
    j.tier = Tier::E;
    j.admit = Admit::Scaffold;
    j.endA = j.endB = EndClass::Unique;
    j.copyA = j.copyB = 1;
    j.gmin = -1;
    j.pMisjoin = -1;
    return j;
}

// The fixture: pieces P0..P8 and the runs between them.
struct Fixture {
    std::vector<std::string> P;         // the N-free pieces, as the uncut genome view writes them
    std::string record;                 // the scaffold record (N-runs; merges as a 1-N butt)
    Ledger ledger;
    std::string fillGraph, fillCons;
};

Fixture build() {
    Fixture f;
    for (int i = 0; i < 9; ++i) f.P.push_back(randomSeq(1500 + 100 * static_cast<size_t>(i)));
    // P6 ends with 40 bases that P7 starts with: the exact-overlap merge (written once in the view).
    const std::string ov = f.P[6].substr(f.P[6].size() - 40);
    f.P[7] = ov + f.P[7];
    f.fillGraph = randomSeq(150);
    f.fillCons = randomSeq(240);
    const size_t runs[8] = {100, 500, 300, 800, 200, 250, 1, 60};
    std::string r = f.P[0];
    for (int k = 0; k < 8; ++k) {
        r += std::string(runs[k], 'N');
        r += f.P[k + 1];
    }
    f.record = r;
    // ledger rows, flanks from the record
    size_t at = f.P[0].size();
    const Source src[8] = {Source::Resolver, Source::Join, Source::Join, Source::Layout,
                           Source::Join, Source::Join, Source::Clonal, Source::Join};
    for (int k = 0; k < 8; ++k) {
        Junction j = blank(static_cast<uint32_t>(k), src[k]);
        j.claimedN = j.writtenN = static_cast<int32_t>(runs[k]);
        j.flankL32 = r.substr(at - 32, 32);
        j.flankR32 = r.substr(at + runs[k], 32);
        if (k == 0) { j.tier = Tier::B; j.verdict = Verdict::Silent; }
        if (k == 2) { j.pairK = 3; j.pairContra = 0; j.tier = Tier::D; }
        if (k == 3) { j.verdict = Verdict::PassWalk; j.tier = Tier::D; }
        if (k == 4) {
            j.fillSeq = f.fillGraph;
            j.fillSpans = {{0, 150, Basis::GraphWalk}};
        }
        if (k == 5) {
            j.fillSeq = f.fillCons;
            j.fillSpans = {{0, 100, Basis::PairPhased}, {100, 140, Basis::Consensus}};
            j.pairK = 1;   // one pair only: not read-pair confirmed
        }
        if (k == 6) { j.mergeOverlap = 40; j.claimedN = -40; j.writtenN = 1; }
        at += runs[k] + f.P[k + 1].size();
        if (k == 7) continue;   // the last run is left without a ledger row: `unrecorded`
        f.ledger.j.push_back(j);
    }
    return f;
}

struct RunOut {
    SurfaceStats st;
    std::string json;
};
RunOut runWriter(const std::string& dir, const Fixture& f) {
    std::filesystem::remove_all(dir);
    std::filesystem::create_directories(dir);
    std::vector<std::string> seqs = {f.record};
    std::vector<std::string> names = {"NODE_1_length_" + std::to_string(f.record.size()) + "_cov_30.0000_chr"};
    std::string err;
    writeFasta(dir + "/scaffolds.fasta", seqs, names, 80, err);
    SurfaceInput in;
    in.outDir = dir;
    in.seqs = &seqs;
    in.names = &names;
    in.scaffoldsFile = true;
    in.ledger = &f.ledger;
    in.modelRan = true;
    RunOut o;
    o.st = writeSurface(in, o.json);
    return o;
}

void testLayoutOnly(const std::string& base) {
    const Fixture f = build();

    // ---- the flag unset: the round-3c writer (reference for the uncut sequence) ----------------
    testenv::clearTesseractEnv();
    setenv("TESSERACT_OM2_OUTPUT", "1", 1);
    const RunOut ref = runWriter(base + "/off", f);
    check(ref.st.ok && !ref.st.layoutOnly, "flag unset: writer ran, not layout-only");
    const auto gOff = readFasta(base + "/off/genome/genome.fasta");
    check(gOff.size() == 1 && firstToken(gOff[0].first) == "chromosome_1", "flag unset: one record chromosome_1");
    check(!std::filesystem::exists(base + "/off/genome/layout.agp"), "flag unset: no layout.agp");
    check(gOff[0].first.find(" layout=") == std::string::npos, "flag unset: no layout header field");
    const std::string hdrOff = slurp(base + "/off/genome/junctions.tsv");
    check(hdrOff.find("layout_view") == std::string::npos, "flag unset: no layout columns in junctions.tsv");
    check(ref.json.find("layout_only") == std::string::npos, "flag unset: no layout block in report.json");
    check(readFasta(base + "/off/scaffolds.fasta").size() == 1, "flag unset: scaffolds.fasta not cut");

    // ---- layout-only, isolate rule --------------------------------------------------------------
    setenv("TESSERACT_OM2_LAYOUT_ONLY", "1", 1);
    const RunOut lo = runWriter(base + "/on", f);
    check(lo.st.ok && lo.st.layoutOnly && lo.st.layoutConfirm == "isolate", "layout-only ran (isolate)");
    check(lo.st.layoutCuts == 4, "4 junctions cut (panel-only join, PASS_WALK layout, allocated fill, merge): " +
                                     std::to_string(lo.st.layoutCuts));
    check(lo.st.confirmedPairsScaffolder == 1 && lo.st.confirmedReadPairs == 1 && lo.st.confirmedGraphFill == 1 &&
              lo.st.confirmedUnrecorded == 1,
          "confirmed: 1 scaffolder, 1 read pairs, 1 graph fill, 1 unrecorded");
    const auto g = readFasta(base + "/on/genome/genome.fasta");
    check(g.size() == 5, "5 genome pieces: " + std::to_string(g.size()));
    const std::vector<std::string> expectNames = {"chromosome_1.1", "chromosome_1.2", "chromosome_1.3",
                                                  "chromosome_1.4", "chromosome_1.5"};
    for (size_t i = 0; i < g.size() && i < expectNames.size(); ++i) {
        check(firstToken(g[i].first) == expectNames[i], "piece name " + expectNames[i]);
        check(g[i].first.find(" layout=chromosome_1 layout_part=" + std::to_string(i + 1) + "/5") != std::string::npos,
              "piece header carries its layout part: " + g[i].first);
    }
    if (g.size() == 5) {
        const std::string N = "N";
        check(g[0].second == f.P[0] + std::string(100, 'N') + f.P[1], "piece 1 = P0 + resolver gap + P1");
        check(g[1].second == f.P[2] + std::string(300, 'N') + f.P[3], "piece 2 = P2 + pair-confirmed gap + P3");
        check(g[2].second == f.P[4] + f.fillGraph + f.P[5], "piece 3 = P4 + graph fill + P5");
        check(g[3].second == f.P[6], "piece 4 = P6 (ends at the unconfirmed merge)");
        check(g[4].second == f.P[7].substr(40) + std::string(60, 'N') + f.P[8],
              "piece 5 = P7 after the shared 40 bases + unrecorded gap + P8");
        check(g[3].second.find(f.fillCons.substr(0, 50)) == std::string::npos &&
                  g[2].second.find(f.fillCons.substr(0, 50)) == std::string::npos,
              "no base of the cut (allocated) fill is written");
        // The pieces, joined at the cut junctions, are the uncut view's record with those junctions removed.
        std::string offNoCut = gOff.empty() ? std::string() : gOff[0].second;
        std::string joined;
        for (const auto& p : g) joined += p.second;
        const size_t removed = 500 + 800 + 240;   // the cut N-runs and the cut fill
        check(offNoCut.size() == joined.size() + removed, "pieces = uncut view minus the cut junctions (length)");
    }
    // layout.agp
    const std::string lagp = base + "/on/genome/layout.agp";
    check(std::filesystem::exists(lagp), "layout.agp written");
    check(validateAgp(lagp).empty(), "layout.agp is valid AGP 2.1: " + validateAgp(lagp));
    const auto rows = readTable(lagp);
    std::string kinds;
    for (const auto& r : rows) kinds += r.size() > 4 ? r[4] : "?";
    check(kinds == "WNWNWNWW", "layout.agp rows W N W N W N W W (the merge abuts): " + kinds);
    if (rows.size() == 8) {
        check(rows[1][5] == "500" && rows[3][5] == "800" && rows[5][5] == "240", "layout gap rows keep the asserted sizes");
        check(rows[1][8].find("align_genus") != std::string::npos, "a panel-ordered layout gap says align_genus");
        check(rows[0][5] == "chromosome_1.1" && rows[7][5] == "chromosome_1.5", "layout components are the pieces");
        check(rows[0][0] == "chromosome_1", "the layout object is the uncut record's name");
    }
    check(validateAgp(base + "/on/genome/genome.agp").empty(), "genome.agp still valid: " +
                                                                   validateAgp(base + "/on/genome/genome.agp"));
    // scaffolds.fasta cut at the same junctions (the merge's butt included); fills are genome-only
    const auto sc = readFasta(base + "/on/scaffolds.fasta");
    check(sc.size() == 5, "scaffolds.fasta cut into 5 records: " + std::to_string(sc.size()));
    check(lo.st.scaffoldParts.size() == 5, "scaffold parts handed back for report.json: " +
                                               std::to_string(lo.st.scaffoldParts.size()));
    // junctions.tsv
    const auto jt = readTable(base + "/on/genome/junctions.tsv");
    check(!jt.empty() && jt[0].size() >= 5 && jt[0][jt[0].size() - 5] == "layout_view", "junctions.tsv layout columns");
    size_t nLayout = 0, nJoined = 0;
    for (size_t i = 1; i < jt.size(); ++i) {
        const auto& r = jt[i];
        if (r.size() < 5) continue;
        const std::string view = r[r.size() - 5];
        if (view == "layout_only") {
            ++nLayout;
            check(r[1] == "layout_gap", "a cut junction's status is layout_gap");
            check(r[r.size() - 4] == "none", "a cut junction has no confirm basis");
        }
        if (view == "joined") ++nJoined;
    }
    check(nLayout == 4 && nJoined == 4, "4 layout_gap rows and 4 joined rows: " + std::to_string(nLayout) + "/" +
                                           std::to_string(nJoined));
    check(lo.json.find("\"layout_only\"") != std::string::npos && lo.json.find("\"layout_gaps\": 4") != std::string::npos,
          "report.json layout block");
    const std::string closure = slurp(base + "/on/genome/closure.txt");
    check(closure.find("#layout\tchromosome_1\tpieces=5\tlayout_gaps=4") != std::string::npos, "closure.txt layout line");

    // ---- c1pass: the PASS_WALK layout gap is kept; the pair-confirmed join and the fills are not -----
    setenv("TESSERACT_OM2_LAYOUT_CONFIRM", "c1pass", 1);
    const RunOut c1 = runWriter(base + "/c1", f);
    check(c1.st.layoutConfirm == "c1pass" && c1.st.confirmedC1Pass == 1, "c1pass keeps the PASS_WALK junction");
    check(c1.st.layoutCuts == 5, "c1pass cuts 5 (panel join, pairs join, two fills, merge): " +
                                     std::to_string(c1.st.layoutCuts));
    unsetenv("TESSERACT_OM2_LAYOUT_CONFIRM");

    // ---- LAYOUT_ONLY alone implies the view ------------------------------------------------------
    unsetenv("TESSERACT_OM2_OUTPUT");
    check(outputEnabled(), "TESSERACT_OM2_LAYOUT_ONLY=1 alone enables the genome view");
    testenv::clearTesseractEnv();
    check(!outputEnabled() && !layoutOnlyEnabled(), "flags cleared: view off");
}

void testConfirmRule() {
    Junction j = blank(1, Source::Join);
    check(confirmBasis(nullptr, false, true, true) == "unrecorded", "no ledger row, ledger present: unrecorded");
    check(confirmBasis(nullptr, false, false, true).empty(), "no ledger after a model ran: not confirmed");
    check(confirmBasis(nullptr, false, false, false) == "unrecorded", "no model, no ledger: the release's own run");
    check(confirmBasis(&j, false, true, true).empty(), "a silent panel join is not confirmed");
    j.pairK = 2; j.pairContra = 1;
    check(confirmBasis(&j, false, true, true).empty(), "pairs crossing but one leaving elsewhere: not confirmed");
    j.pairContra = 0;
    check(confirmBasis(&j, false, true, true) == "read_pairs", ">= 2 crossing pairs, none elsewhere: read_pairs");
    Junction r = blank(2, Source::ResolverUnknown100);
    check(confirmBasis(&r, false, true, true) == "read_pairs_scaffolder", "resolver U100: scaffolder pairs");
    Junction w = blank(3, Source::Clonal);
    w.verdict = Verdict::PassExact;
    check(confirmBasis(&w, false, true, true).empty(), "PASS_EXACT alone does not confirm under `isolate`");
    check(confirmBasis(&w, true, true, true) == "c1_pass_exact", "PASS_EXACT confirms under `c1pass`");
    w.fillSeq = "ACGTACGTAC";
    w.fillSpans = {{0, 10, Basis::ThreadPhased}};
    check(confirmBasis(&w, false, true, true) == "graph_fill", "a thread-phased fill confirms");
    w.fillSpans = {{0, 5, Basis::ThreadPhased}, {5, 5, Basis::PriorAllocated}};
    check(confirmBasis(&w, false, true, true).empty(), "a fill with a prior-allocated span does not");
    w.contigFilled = true;
    w.fillSpans = {{0, 10, Basis::GraphWalk}};
    check(confirmBasis(&w, false, true, true).empty(), "a contig-filled row has no genome fill (no seam either)");
}

void testAlias() {
    check(canonicalOrganism("salmonella") == "senterica", "salmonella -> senterica");
    check(canonicalOrganism("Salmonella") == "senterica", "Salmonella -> senterica");
    check(canonicalOrganism("SALMONELLA") == "senterica", "SALMONELLA -> senterica");
    check(canonicalOrganism("senterica") == "senterica", "senterica unchanged");
    check(canonicalOrganism("klebsiella") == "kpneumoniae", "klebsiella -> kpneumoniae (N22, unchanged)");
    check(canonicalOrganism("ecoli") == "ecoli", "other names unchanged");
    check(canonicalOrganism("salmonella_enterica") == "salmonella_enterica", "only the exact alias is mapped");
}

}  // namespace

int main(int argc, char** argv) {
    testenv::clearTesseractEnv();
    const std::string base = argc > 1 ? std::string(argv[1])
                                      : (std::filesystem::temp_directory_path() /
                                         ("rel150_layout_" + std::to_string(::getpid()))).string();
    testLayoutOnly(base);
    testConfirmRule();
    testAlias();
    if (argc <= 1) std::filesystem::remove_all(base);
    std::printf("test_rel150_layout: %d checks, %d failures\n", checks, failures);
    return failures ? 1 : 0;
}
