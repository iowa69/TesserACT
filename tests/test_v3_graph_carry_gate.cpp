// T01 (build_v3, G-graph): carry-only junctions and TESSERACT_FIX_CARRY_READ_GATE.
//
// Uses the real loci of the K2-only misassembly NODE_337 (A. baumannii GCF016919505v2,
// NZ_CP070362.2 near 2219123 and 1850346): both contain the 20 bp word W at 300..319,
// and the carried chimera A-left + W + B-right is the shape of NODE_337. Counting runs
// through the release path the assembler uses at every rung after the first --
// KmerCounter::count(reads, carry, 4) then extractSolid(2) -- and the gate is the one
// Assembler::iterate() calls.
//
// Before the fix (mode 0, which is also exactly what the release does) the 12 junction
// 33-mers have zero read support, are solid through the carry weight alone, and the k=33
// graph spells the chimeric adjacency. After it (mode 1/2) the junction is withheld, the
// two true junctions survive, and a genuine thin read overlap -- the same zero-read
// signature, but with no competing continuation -- is spared.
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>
#include <unistd.h>
#include <vector>

#include "carry_gate.h"
#include "counter.h"
#include "graph.h"
#include "kmer.h"
#include "seqio.h"

extern char** environ;

namespace {

int failures = 0;
void check(bool ok, const char* what) {
    std::printf("  %s  %s\n", ok ? "PASS" : "FAIL", what);
    if (!ok) ++failures;
}

void clearTesseractEnv() {
    std::vector<std::string> names;
    for (char** e = environ; e && *e; ++e) {
        if (std::strncmp(*e, "TESSERACT_", 10) == 0) {
            const char* eq = std::strchr(*e, '=');
            names.emplace_back(*e, eq ? static_cast<size_t>(eq - *e) : std::strlen(*e));
        }
    }
    for (const std::string& n : names) unsetenv(n.c_str());
}

const char* const kLocusA =
    "GTTCAAACTTAAGTACAATTTTATAAATTTTGAAGAATGCTTAAATCAATTTAATAAAGAAAAATCAGTACATATTGATT"
    "TAAGTCTTATGCCCAAGTTTAAAAATGAAGATGAATTCATTTTGTGGCTTGCAGGCTTTATTGAGCGAATGACAATTGGA"
    "GGGAAAGAGAAGCTTCCTCCAATTTCAAACTATATTCCAAAAGGCTTTAAAATTGATCAAGCCGAAATTCCTGTGGCTGC"
    "TGAACCAAGCAAAGAAGAAAATGCGGAGATGATTATTAACTACTTTAAATCTGAGGATTTCATTAAAAATAATAAAATAT"
    "CCTCTTAAATATAGAAGTATTTGTCCACCTTACTTCTATATCTGGCTGATAATATTCAGCTACCGCCTTCTGGGGCGGTT"
    "TTTTTATTTAACAAGCTTAAATCATCAGCATATAAAACTTATCAAAAATGTCAGTAAATATGGAAATGTTTTATAAAGAA"
    "TTAATAAAGACATAAAAAATAATAAGTTGGAGATTATATTGTTATTTTTATTAAAATATTTTTCACATAAAAATAATCTT"
    "AAATTACGTATAGAAACATTTACTGGTACATAATTTCCAAAAAATAGGATGATGAAGAGA";
const char* const kLocusB =
    "TTGCTATATTGTTAATTATCACAAGCGATAGGGATTCTAGCTTATGAGCACTTCACGAGTAAAGCTAACTAAGTCATTTA"
    "TCGATCAACTTGAACTCACACCGGCTATTTATCGTGACAGTGAGATTATTGGCTTTGCTATTCGTGTAAATAACTCTTAT"
    "AAAACTTACATTGTGGAGAAGAAAGTAAAAGGGAAATCAATTCGTTGTAAATTAGGGGATTATGAAAAAATTACTGTTGA"
    "AGACGCCCGTATACTTGCACAGCAGAAGTTAAAAGAGTTAACTGACTCAAACCCTTTAAGCATTAAAAATAATAAAATAT"
    "TAAAAAATAGTTTAAATGAGAAAATTGATCAGCCTTATCTCAAAGAAGCTTTCCAAGTATATATAAATCACCACGAGTTA"
    "AAAGAAAGGACATTAGCCGACTATAGAGAAGTTATAGAAAAATATTTAATAGACTTAAGTGAACTTAAACTGATTGATAT"
    "TACTGAACAGAGGATTGAGGAAAGGTATACTCAGTTATCTCAGTATAGTTATGCAAAAGCTAATTTATCTATGCGTGTAC"
    "TTAGAGCGGTCTATCGTTTTTCAATTAAGTATTATCAAAATAAAAATTGTGAAGTAATTA";

std::string rc(const std::string& s) { return ts::reverseComplement(s); }

ts::Kmer canon(const std::string& s, int k) {
    bool ok = false;
    return ts::canonical(ts::stringToKmer(s, k, ok), k);
}

bool spells(const ts::UnitigGraph& g, const std::string& q) {
    const std::string r = rc(q);
    for (const ts::Unitig& u : g.nodes) {
        if (u.deleted) continue;
        if (u.seq.find(q) != std::string::npos || u.seq.find(r) != std::string::npos) return true;
    }
    return false;
}

void tile(const std::string& s, size_t from, size_t to, size_t len, size_t step,
          std::vector<std::string>& reads) {
    for (size_t p = from; p + len <= to; p += step) {
        reads.push_back(s.substr(p, len));
        reads.push_back(rc(s.substr(p, len)));
    }
}

std::string tmpDir;

bool load(const std::vector<std::string>& reads, const std::string& name, ts::SequenceStore& store) {
    const std::string path = tmpDir + "/" + name;
    {
        std::ofstream o(path);
        for (size_t i = 0; i < reads.size(); ++i) o << ">r" << i << "\n" << reads[i] << "\n";
    }
    ts::Library lib;
    lib.r1 = path;
    std::string err;
    ts::QualityTrim qt;
    qt.enabled = false;
    store.setQualityTrim(qt);
    if (!store.load({lib}, 1, err)) {
        std::printf("load failed: %s\n", err.c_str());
        return false;
    }
    return true;
}

// The release counting path for a rung with carry-over.
void countRung(const ts::SequenceStore& store, const std::vector<std::string>& carry, int k,
               uint32_t weight, uint32_t cutoff, ts::KmerTable& solid) {
    ts::KmerCounter c(k, 2);
    c.count(store, carry, weight);
    c.extractSolid(cutoff, solid);
}

ts::UnitigGraph graphOf(const ts::KmerTable& solid, int k) {
    ts::UnitigGraph g = ts::UnitigGraph::build(solid, k, 1);
    g.compact();
    return g;
}

void show(const char* tag, const ts::CarryGateStats& s) {
    std::printf("    [%s] carriedKmers=%zu readAbsent=%zu runs=%zu runsTooLong=%zu "
                "runsNoCompetitor=%zu candidates=%zu runsGated=%zu kmersWithheld=%zu\n",
                tag, s.carriedKmers, s.readAbsent, s.runs, s.runsTooLong, s.runsNoCompetitor,
                s.candidates, s.runsGated, s.kmersWithheld);
}

}  // namespace

int main() {
    clearTesseractEnv();
    char tmpl[] = "/tmp/tess_v3_carrygate.XXXXXX";
    const char* d = mkdtemp(tmpl);
    if (!d) { std::printf("mkdtemp failed\n"); return 2; }
    tmpDir = d;

    const std::string A = kLocusA, B = kLocusB;
    const std::string W = "CATTAAAAATAATAAAATAT";
    const int k = 33;
    const uint32_t w = 4, cutoff = 2;
    check(A.substr(300, 20) == W && B.substr(300, 20) == W, "fixture: W at 300..319 in both loci");
    const std::string chimera = A.substr(0, 320) + B.substr(320);
    const std::string chimJunction = chimera.substr(299, 22);   // the 22-mer no read holds

    // ---------------- Case 1: the NODE_337 chimera ----------------
    std::printf("case 1: carried chimera A-left+W+B-right, k=%d, weight %u, cutoff %u\n", k, w, cutoff);
    std::vector<std::string> reads;
    tile(A, 0, A.size(), 150, 10, reads);
    tile(B, 0, B.size(), 150, 10, reads);
    ts::SequenceStore store;
    if (!load(reads, "case1.fa", store)) return 2;
    ts::KmerTable base;
    countRung(store, {chimera}, k, w, cutoff, base);
    size_t junctionSolid = 0;
    for (size_t s = 320 - k + 1; s < 300; ++s)
        if (base.get(canon(chimera.substr(s, k), k)) != 0) ++junctionSolid;
    check(junctionSolid == 12, "release path: all 12 zero-read junction 33-mers are solid via the carry");
    {
        ts::KmerTable solid = base;
        const ts::CarryGateStats s = ts::gateCarryOnlyJunctions({chimera}, solid, k, w, cutoff, 0);
        show("mode 0", s);
        check(s.readAbsent == 12 && s.runs == 1 && s.candidates == 1, "dry run finds exactly the junction run");
        check(s.runsGated == 0 && s.kmersWithheld == 0 && solid.size() == base.size(),
              "dry run changes nothing");
        check(spells(graphOf(solid, k), chimJunction),
              "mode 0 (= release): k=33 graph spells the chimeric A|W|B adjacency (the defect)");
    }
    for (int mode = 1; mode <= 2; ++mode) {
        ts::KmerTable solid = base;
        const ts::CarryGateStats s = ts::gateCarryOnlyJunctions({chimera}, solid, k, w, cutoff, mode);
        show(mode == 1 ? "mode 1" : "mode 2", s);
        check(s.runsGated == 1 && s.kmersWithheld == 12, "junction run withheld (12 k-mers)");
        check(solid.size() + 12 == base.size(), "nothing else leaves the table");
        const ts::UnitigGraph g = graphOf(solid, k);
        check(!spells(g, chimJunction), "chimeric adjacency gone from the k=33 graph");
        check(spells(g, A.substr(280, 60)) && spells(g, B.substr(280, 60)),
              "both true junctions still spelled");
        check(g.validate().empty(), "graph valid");
    }

    // ---------------- Case 2: genuine thin overlap (control) ----------------
    std::printf("case 2: control, one locus whose reads overlap by only 25 bp\n");
    std::vector<std::string> reads2;
    tile(A, 0, 330, 150, 10, reads2);
    tile(A, 305, A.size(), 150, 10, reads2);
    ts::SequenceStore store2;
    if (!load(reads2, "case2.fa", store2)) return 2;
    ts::KmerTable base2;
    countRung(store2, {A}, k, w, cutoff, base2);
    for (int mode = 0; mode <= 2; ++mode) {
        ts::KmerTable solid = base2;
        const ts::CarryGateStats s = ts::gateCarryOnlyJunctions({A}, solid, k, w, cutoff, mode);
        show(mode == 0 ? "mode 0" : (mode == 1 ? "mode 1" : "mode 2"), s);
        check(s.runs == 1 && s.readAbsent >= 7, "same zero-read signature as the chimera (a flanked run)");
        check(s.runsNoCompetitor == 1 && s.candidates == 0 && s.kmersWithheld == 0,
              "no competing continuation: the genuine bridge is spared");
        check(spells(graphOf(solid, k), A.substr(280, 80)), "thin bridge still spelled");
    }

    // ---------------- Case 3: exact multiplicity arithmetic ----------------
    std::printf("case 3: a contig carried twice\n");
    {
        ts::KmerTable twice;
        countRung(store, {chimera, chimera}, k, w, cutoff, twice);
        ts::KmerTable solid = twice;
        const ts::CarryGateStats s =
            ts::gateCarryOnlyJunctions({chimera, chimera}, solid, k, w, cutoff, 1);
        show("mode 1", s);
        check(s.readAbsent == 24 && s.runs == 2 && s.runsGated == 2,
              "read counts subtract the carry weight of every occurrence");
        check(s.kmersWithheld == 12 && !spells(graphOf(solid, k), chimJunction),
              "both occurrences withheld: the junction leaves the table");
        ts::KmerTable aTwice;
        countRung(store, {A, A}, k, w, cutoff, aTwice);
        const ts::CarryGateStats sa = ts::gateCarryOnlyJunctions({A, A}, aTwice, k, w, cutoff, 1);
        check(sa.readAbsent == 0 && sa.runs == 0, "a fully read-covered contig carried twice has no zero-read k-mer");
    }

    // ---------------- Case 4: competitor on one flank only ----------------
    std::printf("case 4: only the left flank shows a competing continuation\n");
    {
        std::vector<std::string> reads4;
        tile(A, 0, A.size(), 150, 10, reads4);
        tile(B, 300, B.size(), 150, 10, reads4);   // B-left never sequenced
        ts::SequenceStore store4;
        if (!load(reads4, "case4.fa", store4)) return 2;
        ts::KmerTable base4;
        countRung(store4, {chimera}, k, w, cutoff, base4);
        ts::KmerTable s1 = base4, s2 = base4;
        const ts::CarryGateStats g1 = ts::gateCarryOnlyJunctions({chimera}, s1, k, w, cutoff, 1);
        const ts::CarryGateStats g2 = ts::gateCarryOnlyJunctions({chimera}, s2, k, w, cutoff, 2);
        show("mode 1", g1);
        show("mode 2", g2);
        check(g1.runsGated == 1 && !spells(graphOf(s1, k), chimJunction), "mode 1: either flank suffices");
        check(g2.runsGated == 0 && g2.runsNoCompetitor == 1 && spells(graphOf(s2, k), chimJunction),
              "mode 2: both flanks required");
    }

    // ---------------- Case 5: nothing carried ----------------
    {
        ts::KmerTable solid = base;
        const ts::CarryGateStats s = ts::gateCarryOnlyJunctions({}, solid, k, w, cutoff, 1);
        check(s.carriedKmers == 0 && s.kmersWithheld == 0 && solid.size() == base.size(),
              "first rung (no carry): no-op");
        const ts::CarryGateStats z = ts::gateCarryOnlyJunctions({chimera}, solid, k, 0, cutoff, 1);
        check(z.carriedKmers == 0 && solid.size() == base.size(), "carry weight 0: no-op");
    }

    std::printf("test_v3_graph_carry_gate: %s (%d failure%s)\n", failures ? "FAILED" : "ALL PASSED",
                failures, failures == 1 ? "" : "s");
    std::string cmd = "rm -rf '" + tmpDir + "'";
    if (std::system(cmd.c_str()) != 0) std::printf("note: could not remove %s\n", tmpDir.c_str());
    return failures ? 1 : 0;
}
