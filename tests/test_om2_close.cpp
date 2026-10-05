// Organism Model 2.0, C2: the close stage end to end, on exact de Bruijn graphs of synthetic
// genomes and error-free paired reads simulated from them.
//
//   rrna8     an 8-copy rRNA-like operon (1.8-2 kb, longer than any fragment) with three ITS types
//             (4/3/1 copies), a variant 60 bp inside each operon start and one 100 bp before each end
//             (reachable by reads anchored in the unique flank), and a mid-operon variant nobody can
//             phase. The layout's scaffold holds 8 N-runs. Checked: every locus gets the TRUE end
//             alleles (phased), ITS and mid variants are allocated and labelled, nothing is written
//             into the sequence (genome-only), the variant table lists every branch; with the
//             leave-clone-out locus prior every ITS type is right and the flips are PRIOR_ALLOCATED.
//   shortrep  a 150 bp repeat that read pairs span: contig-grade fills, sequence restored exactly.
//   is3       three identical 1.2 kb IS copies: ISOLATE_REPEAT_EXACT genome-only fills; the same with
//             one copy carrying a SNV (a bubble, allocated); and with depth saying one copy (the
//             single-copy budget refuses the second and third use).
//   hairpin   a walk that needs an inverted repeat in both orientations: refused.
//   circle    a circular plasmid closes on its unique wrap walk and dnaA is located on it; a linear
//             record is never rotated; a plasmid whose wrap crosses a bubble does not close.
//   off       the stage off: nothing read, nothing changed, counters zero.
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

std::mt19937 rng(20260928);
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
        : dir(std::filesystem::temp_directory_path() / ("tesseract-om2close-" + tag + "-" + std::to_string(getpid()))) {
        std::filesystem::create_directories(dir);
    }
    ~Temp() { std::filesystem::remove_all(dir); }
};

// Error-free pairs, 2x150, fragments 350 +- 30, ~40x.
SequenceStore simulate(const std::vector<Mol>& mols, const Temp& tmp) {
    std::normal_distribution<double> ins(350, 30);
    std::ofstream a(tmp.dir / "r1.fa"), b(tmp.dir / "r2.fa");
    size_t id = 0;
    for (const Mol& m : mols) {
        const size_t L = m.seq.size();
        const size_t pairs = L * 40 / 300;
        const std::string ext = m.circular ? m.seq + m.seq.substr(0, 600) : m.seq;
        for (size_t p = 0; p < pairs; ++p) {
            const size_t f = static_cast<size_t>(std::max(250.0, std::min(450.0, ins(rng))));
            if (!m.circular && L <= f) continue;
            const size_t start = m.circular ? rng() % L : rng() % (L - f);
            const std::string frag = ext.substr(start, f);
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

InsertModel insert350() {
    InsertModel m;
    m.mean = 350;
    m.stddev = 30;
    m.usable = true;
    return m;
}

CloseOptions onOpts(FillMode fill, AllocMode alloc) {
    CloseOptions o;
    o.enabled = true;
    o.fill = fill;
    o.alloc = alloc;
    o.allocParams.mode = alloc;
    return o;
}

std::string readFile(const std::filesystem::path& p) {
    std::ifstream in(p);
    std::stringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

size_t countBasis(const Junction& j, Basis b) {
    size_t n = 0;
    for (const FillSpan& s : j.fillSpans) n += s.basis == b;
    return n;
}

bool spansTile(const Junction& j) {
    uint32_t at = 0;
    for (const FillSpan& s : j.fillSpans) {
        if (s.off != at || s.len == 0) return false;
        at += s.len;
    }
    return at == j.fillSeq.size();
}

// ---- rrna8 ------------------------------------------------------------------------------------

struct Rrna {
    std::vector<std::string> U;       // 9 unique segments
    std::vector<std::string> O;       // 8 operon copies
    std::string core1, core2, its[3];
    const int itsOf[8] = {0, 0, 0, 0, 1, 1, 1, 2};
    const int e1Of[8] = {0, 1, 0, 0, 1, 0, 0, 1};
    const int e2Of[8] = {1, 0, 0, 1, 0, 0, 1, 0};
    const int mOf[8] = {0, 0, 1, 0, 0, 0, 1, 0};
    std::string genome, record;
    Rrna() {
        for (int i = 0; i < 9; ++i) U.push_back(dna(3000));
        core1 = dna(600);
        core2 = dna(1000);
        its[0] = dna(150);
        its[1] = dna(260);
        its[2] = dna(210);
        for (int i = 0; i < 8; ++i) {
            std::string c1 = core1, c2 = core2;
            if (e1Of[i]) c1[60] = other(c1[60]);
            if (mOf[i]) c2[500] = other(c2[500]);
            if (e2Of[i]) c2[900] = other(c2[900]);
            O.push_back(c1 + its[itsOf[i]] + c2);
        }
        genome = U[0];
        record = U[0];
        for (int i = 0; i < 8; ++i) {
            genome += O[static_cast<size_t>(i)] + U[static_cast<size_t>(i + 1)];
            record += std::string(2000, 'N') + U[static_cast<size_t>(i + 1)];
        }
    }
    // The 41 bp around a variant of copy i, as it must read in a fill of locus i.
    std::string e1Window(int i) const { return O[static_cast<size_t>(i)].substr(40, 41); }
    std::string e2Window(int i) const {
        const size_t p = core1.size() + its[itsOf[i]].size() + 900;
        return O[static_cast<size_t>(i)].substr(p - 20, 41);
    }
};

void testRrna8() {
    Rrna R;
    const std::vector<Mol> mols = {{R.genome, false}};
    UnitigGraph g = buildGraph(mols);
    check(g.validate().empty(), "rrna8 graph is valid");
    Temp tmp("rrna8");
    const SequenceStore reads = simulate(mols, tmp);
    const InsertModel ins = insert350();

    // -- phased allocation, genome-only fills --
    {
        std::vector<std::string> seqs = {R.record};
        Ledger led;
        CloseInputs in{g, reads, ins, (tmp.dir / "out1").string(), "", 1, false};
        const CloseStats st = runCloseStage(in, onOpts(FillMode::Genome, AllocMode::Phased), seqs, led);
        check(st.enabled && st.junctions == 8 && st.anchored == 8, "8 junctions, all anchored");
        check(st.bridges == 8 && st.tiebreakUnverified == 8, "8 bridges, model tie-break (standalone, unverified)");
        check(st.fillsGenome == 8 && st.fillsContig == 0, "8 genome-only fills, none written into the sequence");
        check(seqs[0] == R.record, "genome-only: the record is unchanged");
        check(st.variantSites == 4, "4 variant sites (E1, ITS, M, E2), got " + std::to_string(st.variantSites));
        check(st.alloc.pair + st.alloc.thread == 16, "both end variants phased at all 8 loci, got " +
                                                         std::to_string(st.alloc.pair + st.alloc.thread));
        check(st.alloc.prior == 0, "no prior in phased mode");
        check(st.clusters == 1, "one repeat cluster");
        check(led.j.size() == 8, "8 ledger rows");
        int itsRight = 0;
        for (int i = 0; i < 8; ++i) {
            const Junction& j = led.j[static_cast<size_t>(i)];
            check(!j.contigFilled && !j.fillSeq.empty(), "locus " + std::to_string(i) + " has a genome-only fill");
            check(spansTile(j), "mask spans tile the fill exactly at locus " + std::to_string(i));
            check(j.fillSeq.find(R.e1Window(i)) != std::string::npos, "true 5' end allele at locus " + std::to_string(i));
            check(j.fillSeq.find(R.e2Window(i)) != std::string::npos, "true 3' end allele at locus " + std::to_string(i));
            check(countBasis(j, Basis::ThreadPhased) + countBasis(j, Basis::PairPhased) >= 2, "phased spans");
            check(countBasis(j, Basis::IsolateRepeatExact) >= 1, "shared operon body labelled ISOLATE_REPEAT_EXACT");
            check(countBasis(j, Basis::PriorAllocated) == 0, "no prior label");
            check(j.fillSeq.size() >= 1700 && j.fillSeq.size() <= 1900, "fill is operon-sized");
            itsRight += j.fillSeq.find(R.its[R.itsOf[i]]) != std::string::npos;
            check(j.exhaustive && !j.hairpin && j.walks.size() >= 2, "exhaustive, 3 ITS length classes");
        }
        check(itsRight == 4, "without a prior the ITS is the consensus type: right at 4/8 loci");
        const std::string rv = readFile(tmp.dir / "out1" / "om2_close" / "repeat_variants.tsv");
        check(std::count(rv.begin(), rv.end(), '\n') == 2 + 2 + 3 + 2 + 2, "repeat_variants.tsv: 9 branch rows");
        check(rv.find("\tITS\t") != std::string::npos, "the ITS site is reported as ITS-class");
        check(rv.find("THREAD_PHASED") != std::string::npos || rv.find("PAIR_PHASED") != std::string::npos,
              "phased loci listed");
        const std::string bed = readFile(tmp.dir / "out1" / "om2_close" / "fills.mask.bed");
        check(std::count(bed.begin(), bed.end(), '\n') >= 8 * 5, "mask BED rows for every labelled span");
    }

    // -- the leave-clone-out locus prior --
    {
        const std::string model = (tmp.dir / "model.tsm").string();
        std::ofstream(model) << "stand-in model";
        const std::string md5 = md5File(model);
        auto canon = [](const std::string& s) {
            uint64_t f = 0, r = 0;
            for (size_t i = 0; i < 31; ++i) {
                const int c = baseCode(s[i]);
                f = (f << 2) | static_cast<uint64_t>(c);
                r = (r >> 2) | (static_cast<uint64_t>(3 - c) << 60);
            }
            return f <= r ? f : r;
        };
        auto hashes = [&](const std::string& s, uint64_t denom) {
            std::string out;
            for (size_t p = 0; p + 31 <= s.size(); ++p) {
                const uint64_t h = sidecarHash(canon(s.substr(p, 31)));
                if (denom && h > ~0ULL / denom) continue;
                char b[24];
                std::snprintf(b, sizeof b, "%s%016llx", out.empty() ? "" : ",", static_cast<unsigned long long>(h));
                out += b;
            }
            return out;
        };
        std::string pr = "#om2rrn v1\n#organism synthetic\n#tsm_md5 " + md5 + "\n#genomes 200\n#flank_denom 32\n"
                         "#flank_bp 1500\n";
        for (int c = 0; c < 3; ++c) pr += "C " + std::to_string(c) + " 0 " + hashes(R.its[c], 0) + "\n";
        for (int i = 0; i < 8; ++i) {
            const std::string flank = R.U[static_cast<size_t>(i)].substr(1500) + R.U[static_cast<size_t>(i + 1)].substr(0, 1500);
            pr += "L " + std::to_string(i) + " 200 " + std::to_string(R.itsOf[i]) + ":200 " + hashes(flank, 32) + "\n";
        }
        const std::string prPath = (tmp.dir / "syn.om2rrn").string();
        std::ofstream(prPath) << pr;

        std::vector<std::string> seqs = {R.record};
        Ledger led;
        CloseOptions o = onOpts(FillMode::Genome, AllocMode::Prior);
        o.rrnPrior = prPath;
        CloseInputs in{g, reads, ins, "", model, 1, false};
        const CloseStats st = runCloseStage(in, o, seqs, led);
        check(st.priorLoaded, "prior loaded under its md5 pin");
        check(st.alloc.prior == 3, "the three type-B loci flip to PRIOR_ALLOCATED, got " + std::to_string(st.alloc.prior));
        for (int i = 0; i < 8; ++i) {
            const Junction& j = led.j[static_cast<size_t>(i)];
            check(j.fillSeq.find(R.its[R.itsOf[i]]) != std::string::npos,
                  "with the prior the ITS type is right at locus " + std::to_string(i));
            check((countBasis(j, Basis::PriorAllocated) > 0) == (R.itsOf[i] == 1),
                  "PRIOR_ALLOCATED exactly where the prior flips the choice (locus " + std::to_string(i) + ")");
            check(j.fillSeq.find(R.e1Window(i)) != std::string::npos, "end alleles still phased with the prior");
        }
        // pinned to another model: refused, and the stage falls back to phased
        std::vector<std::string> seqs2 = {R.record};
        Ledger led2;
        const std::string other = (tmp.dir / "other.tsm").string();
        std::ofstream(other) << "another model";
        CloseInputs in2{g, reads, ins, "", other, 1, false};
        const CloseStats st2 = runCloseStage(in2, o, seqs2, led2);
        check(!st2.priorLoaded && st2.alloc.prior == 0, "a prior pinned to another model is refused");
    }

    // -- allocation off: sites cannot be filled; consensus: majority everywhere --
    {
        std::vector<std::string> seqs = {R.record};
        Ledger led;
        CloseInputs in{g, reads, ins, "", "", 1, false};
        const CloseStats st = runCloseStage(in, onOpts(FillMode::Genome, AllocMode::Off), seqs, led);
        check(st.fillsGenome == 0 && st.fillsContig == 0, "ALLOC=off: no fill crosses a bubble");
        std::vector<std::string> seqs2 = {R.record};
        Ledger led2;
        const CloseStats st2 = runCloseStage(in, onOpts(FillMode::Genome, AllocMode::Consensus), seqs2, led2);
        check(st2.fillsGenome == 8 && st2.alloc.consensus == 32, "ALLOC=consensus: 32 consensus allocations");
        // FILL=contig does not write an allocated operon into the sequence
        std::vector<std::string> seqs3 = {R.record};
        Ledger led3;
        const CloseStats st3 = runCloseStage(in, onOpts(FillMode::Contig, AllocMode::Phased), seqs3, led3);
        check(st3.fillsContig == 0 && seqs3[0] == R.record, "allocation never reaches contigs.fasta (I4)");
    }
}

// ---- shortrep: a repeat read pairs span --------------------------------------------------------

void testShortRepeat() {
    const std::string V0 = dna(2000), V1 = dna(2000), V2 = dna(2000), R = dna(150);
    const std::string genome = V0 + R + V1 + R + V2;
    const std::vector<Mol> mols = {{genome, false}};
    UnitigGraph g = buildGraph(mols);
    Temp tmp("short");
    const SequenceStore reads = simulate(mols, tmp);
    const InsertModel ins = insert350();
    const std::string record = V0 + std::string(150, 'N') + V1 + std::string(1, 'N') + V2;
    std::vector<std::string> seqs = {record};
    Ledger led;
    CloseInputs in{g, reads, ins, "", "", 1, false};
    const CloseStats st = runCloseStage(in, onOpts(FillMode::Contig, AllocMode::Phased), seqs, led);
    check(st.pairThread == 2, "both seams spanned by read pairs, got " + std::to_string(st.pairThread));
    check(st.fillsContig == 2, "two contig-grade fills");
    check(seqs[0] == genome, "the record now reads exactly as the genome");
    check(led.j[0].contigFilled && led.j[0].tier == Tier::A, "ledger: contig-filled, tier A");
    // genome mode: the same fills, kept out of the sequence
    std::vector<std::string> seqs2 = {record};
    Ledger led2;
    const CloseStats st2 = runCloseStage(in, onOpts(FillMode::Genome, AllocMode::Phased), seqs2, led2);
    check(st2.fillsGenome == 2 && seqs2[0] == record, "FILL=genome writes nothing into the sequence");
    check(led2.j[0].fillSeq == R && led2.j[1].fillSeq == R, "genome-only fills carry the repeat");
}

// ---- is3: IS copies -----------------------------------------------------------------------------

void testIs3() {
    const std::string W0 = dna(2500), W1 = dna(2500), W2 = dna(2500), W3 = dna(2500), IS = dna(1200);
    std::string IS2 = IS;
    IS2[600] = other(IS2[600]);
    const std::string rec = W0 + std::string(1, 'N') + W1 + std::string(1, 'N') + W2 + std::string(1, 'N') + W3;
    Temp tmp("is3");
    const InsertModel ins = insert350();
    // uniform cluster
    {
        const std::vector<Mol> mols = {{W0 + IS + W1 + IS + W2 + IS + W3, false}};
        UnitigGraph g = buildGraph(mols);
        const SequenceStore reads = simulate(mols, tmp);
        std::vector<std::string> seqs = {rec};
        Ledger led;
        CloseInputs in{g, reads, ins, "", "", 1, false};
        const CloseStats st = runCloseStage(in, onOpts(FillMode::Contig, AllocMode::Phased), seqs, led);
        check(st.bridges == 3 && st.fillsGenome == 3 && st.fillsContig == 0,
              "uniform IS: 3 genome-only fills (1.2 kb is beyond read-pair reach)");
        check(st.variantSites == 0, "uniform IS: no variant site");
        for (int i = 0; i < 3; ++i) {
            const Junction& j = led.j[static_cast<size_t>(i)];
            check(j.fillSeq == IS, "the IS copy is filled exactly");
            check(j.fillSpans.size() == 1 && j.fillSpans[0].basis == Basis::IsolateRepeatExact,
                  "labelled ISOLATE_REPEAT_EXACT");
            check(j.verdict == Verdict::PassExact, "a unique walk: PASS_EXACT");
        }
        // depth says one copy: the hard single-copy budget allows one use
        for (Unitig& u : g.nodes) {
            if (u.seq.size() >= 1000 && u.coverage > 2.5) u.coverage = 1.0;
        }
        std::vector<std::string> seqs2 = {rec};
        Ledger led2;
        const CloseStats st2 = runCloseStage(in, onOpts(FillMode::Genome, AllocMode::Phased), seqs2, led2);
        check(st2.budgetRefused == 2 && st2.fillsGenome == 1, "single-copy budget: one use, two refused");
        // a sized gap far longer than every walk: the walk is a shortcut, refused
        const UnitigGraph g3 = buildGraph(mols);
        std::vector<std::string> seqs3 = {W0 + std::string(4500, 'N') + W1, W2, W3};   // W2, W3 placed
        Ledger led3;
        CloseInputs in3{g3, reads, ins, "", "", 1, false};
        const CloseStats st3 = runCloseStage(in3, onOpts(FillMode::Genome, AllocMode::Phased), seqs3, led3);
        check(st3.lengthInconsistent == 1 && st3.fillsGenome == 0, "a 1.2 kb walk does not fill a 4.5 kb gap");
    }
    // one copy carries a SNV: a bubble
    {
        const std::vector<Mol> mols = {{W0 + IS + W1 + IS2 + W2 + IS + W3, false}};
        UnitigGraph g = buildGraph(mols);
        const SequenceStore reads = simulate(mols, tmp);
        std::vector<std::string> seqs = {rec};
        Ledger led;
        CloseInputs in{g, reads, ins, "", "", 1, false};
        const CloseStats st = runCloseStage(in, onOpts(FillMode::Genome, AllocMode::Phased), seqs, led);
        check(st.fillsGenome == 3 && st.variantSites == 1, "bubbled IS: 3 fills, 1 variant site");
        check(st.alloc.consensus + st.alloc.multiplicity == 3, "the mid-IS SNV is beyond reach: allocated, labelled");
        size_t labelled = 0;
        for (const Junction& j : led.j) labelled += countBasis(j, Basis::Consensus) + countBasis(j, Basis::Multiplicity);
        check(labelled == 3, "one allocation span per locus");
    }
}

// ---- hairpin ----------------------------------------------------------------------------------

void testHairpin() {
    const std::string A = dna(2500), B = dna(300), C = dna(2500), R = dna(500);
    const std::vector<Mol> mols = {{A + R + B + rc(R) + C, false}};
    UnitigGraph g = buildGraph(mols);
    Temp tmp("hairpin");
    const SequenceStore reads = simulate(mols, tmp);
    const InsertModel ins = insert350();
    std::vector<std::string> seqs = {A + std::string(1500, 'N') + C};
    const std::string before = seqs[0];
    Ledger led;
    CloseInputs in{g, reads, ins, "", "", 1, false};
    const CloseStats st = runCloseStage(in, onOpts(FillMode::Contig, AllocMode::Phased), seqs, led);
    check(st.hairpinRefused == 1, "a walk through both orientations of a unitig is refused");
    check(st.fillsContig == 0 && st.fillsGenome == 0 && seqs[0] == before, "nothing filled");
    check(led.j[0].hairpin, "ledger: hairpin flag");
}

// ---- circles and dnaA -------------------------------------------------------------------------

void testCircle() {
    const std::string P = dna(8000);
    const std::string A = dna(2000), B = dna(2000), R1 = dna(400);
    std::string R2 = R1;
    R2[200] = other(R2[200]);
    const std::string Q = A + R1 + B + R2;
    const std::string lin = dna(6000);
    // IS-X-IS: a chromosomal segment X flanked by the same 965 bp repeat on both sides. Its wrap
    // walk (X tail -> IS -> X head) is unique, and it is not a circle.
    const std::string IS = dna(965), X = dna(2539), A2 = dna(3000), B2 = dna(3000);
    const std::vector<Mol> mols = {{P, true}, {Q, true}, {lin, false}, {A2 + IS + X + IS + B2, false}};
    UnitigGraph g = buildGraph(mols);
    Temp tmp("circle");
    const SequenceStore reads = simulate(mols, tmp);
    const InsertModel ins = insert350();
    // records: P from 3000; Q missing its R2 copy at the wrap; a linear molecule
    const std::string recP = P.substr(3000) + P.substr(0, 3000);
    const std::string recQ = A + R1 + B;
    std::vector<std::string> seqs = {recP, recQ, lin, X};
    // a dnaA sketch: P's "chromosome start" at P[5000], and the linear record's first 1.5 kb
    const std::string model = (tmp.dir / "model.tsm").string();
    std::ofstream(model) << "stand-in model";
    std::string sk = "#om2dnaa v1\n#organism synthetic\n#tsm_md5 " + md5File(model) + "\n#genomes 100\n#window 1500\n";
    auto addSketch = [&](const std::string& start) {
        for (size_t p = 0; p + 31 <= 1500; ++p) {
            const std::string km = start.substr(p, 31);
            uint64_t f = 0, r = 0;
            for (size_t i = 0; i < 31; ++i) {
                const int c = baseCode(km[i]);
                f = (f << 2) | static_cast<uint64_t>(c);
                r = (r >> 2) | (static_cast<uint64_t>(3 - c) << 60);
            }
            char line[80];
            std::snprintf(line, sizeof line, "%016llx %zu %d 90\n",
                          static_cast<unsigned long long>(sidecarHash(f <= r ? f : r)), p, f <= r ? 0 : 1);
            sk += line;
        }
    };
    addSketch(P.substr(5000, 1500));
    addSketch(lin.substr(0, 1500));
    const std::string skPath = (tmp.dir / "syn.om2dnaa").string();
    std::ofstream(skPath) << sk;
    Ledger led;
    CloseOptions o = onOpts(FillMode::Contig, AllocMode::Phased);
    o.circ = true;
    o.dnaa = skPath;
    CloseInputs in{g, reads, ins, (tmp.dir / "out").string(), model, 1, false};
    const CloseStats st = runCloseStage(in, o, seqs, led);
    check(st.dnaaLoaded, "sketch loaded");
    check(st.circlesTested == 4, "four wraps tested");
    check(st.circlesClosed == 1, "only the plasmid with a unique wrap walk closes, got " + std::to_string(st.circlesClosed));
    check(seqs[0] == recP, "P closes with no bases to add (its ends meet exactly)");
    check(led.rotateOffset[0] == 2000, "dnaA located at 2000 on the closed record, got " +
                                           std::to_string(led.rotateOffset[0]));
    check(led.rotateOffset[1] == -1, "the plasmid whose wrap crosses a bubble is not closed and not rotated");
    check(led.rotateOffset[2] == -1, "a linear record is never rotated, although dnaA k-mers are on it");
    check(st.rotated == 1, "one rotation");
    check(led.rotateOffset[3] == -1 && seqs[3] == X, "IS-X-IS: the unique walk through the repeat does not close X");
    const std::string circTsv = readFile(tmp.dir / "out" / "om2_close" / "circles.tsv");
    check(circTsv.find("wrap_through_repeat_unverified") != std::string::npos, "and circles.tsv says why");
    const std::string circ = readFile(tmp.dir / "out" / "om2_close" / "circles.tsv");
    check(circ.find("wrap_not_unique") != std::string::npos, "circles.tsv names why Q stays open");
}

// ---- off --------------------------------------------------------------------------------------

void testOff() {
    const std::vector<Mol> mols = {{dna(3000), false}};
    UnitigGraph g = buildGraph(mols);
    Temp tmp("off");
    const SequenceStore reads = simulate(mols, tmp);
    const InsertModel ins = insert350();
    std::vector<std::string> seqs = {mols[0].seq.substr(0, 1000) + "NNNNN" + mols[0].seq.substr(1005)};
    const std::vector<std::string> before = seqs;
    Ledger led;
    CloseInputs in{g, reads, ins, (tmp.dir / "out").string(), "", 1, false};
    const CloseStats st = runCloseStage(in, CloseOptions(), seqs, led);
    check(!st.enabled && seqs == before && led.j.empty(), "off: nothing changed, no ledger");
    check(!std::filesystem::exists(tmp.dir / "out"), "off: no side files");
    const std::string line = formatCloseCounters(st);
    check(line.rfind("[om2-close] enabled=0 clusters=0 bridges=0 forced=0 pair_thread=0 tiebreak=0 open=0 "
                     "budget_refused=0 hairpin_refused=0 fills_contig=0 fills_genome=0 alloc_pair=0 alloc_thread=0 "
                     "alloc_prior=0 alloc_multiplicity=0 alloc_consensus=0 variant_sites=0 circles_tested=0 "
                     "circles_closed=0 rotated=0 abstain_tangled=0",
                     0) == 0,
          "counter line: the design's fields, in order, zeros when off");
    // environment parsing
    setenv("TESSERACT_OM2_CLOSE", "1", 1);
    setenv("TESSERACT_OM2_FILL", "contig", 1);
    setenv("TESSERACT_OM2_ALLOC", "prior", 1);
    setenv("TESSERACT_OM2_CIRC", "1", 1);
    CloseOptions e = CloseOptions::fromEnv();
    check(e.enabled && e.fill == FillMode::Contig && e.alloc == AllocMode::Prior && e.circ && !e.flow,
          "flags parse");
    testenv::clearTesseractEnv();
    setenv("TESSERACT_OM2_FILL", "genome", 1);
    CloseOptions e2 = CloseOptions::fromEnv();
    check(!e2.enabled && e2.anyFlagSet, "a sub-flag without TESSERACT_OM2_CLOSE leaves the stage off");
    testenv::clearTesseractEnv();
    CloseOptions e3 = CloseOptions::fromEnv();
    check(!e3.enabled && !e3.anyFlagSet && e3.fill == FillMode::Genome && e3.alloc == AllocMode::Phased,
          "defaults");
}

}  // namespace

int main() {
    testenv::clearTesseractEnv();
    try {
        testOff();
        testShortRepeat();
        testIs3();
        testHairpin();
        testCircle();
        testRrna8();
    } catch (const std::exception& e) {
        std::fprintf(stderr, "test_om2_close FAILED: %s\n", e.what());
        return 1;
    }
    std::printf("test_om2_close: %d checks passed\n", checks);
    return 0;
}
