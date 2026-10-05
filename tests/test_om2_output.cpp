// Component tests for the Organism Model 2.0 C3 output writer (src/om2_output.*).
//
// 1. Flags off: nothing is written, report.json gets no om2 block, the counter line is zeros.
// 2. A synthetic multi-copy rRNA case with inter-copy variants, written from a synthetic
//    junction ledger the way C1 and C2 will hand it over: seven operon copies (three on the
//    minus strand) whose copies differ at a left-end SNP and a right-end SNP (phased by pairs),
//    an ITS region of three types carried 4/2/1 (allocated by the population prior, or by
//    consensus), and a mid-operon SNP no pair reaches (consensus). Also: a junction admitted
//    to the genome view only, a junction the admission breaks, a seam without a ledger entry,
//    a pairs-only circular plasmid, a two-record plasmid group, and a chromosome wrap closed
//    by a graph walk and rotated to a dnaA offset. The output is checked here against the
//    contract of om2/design/eval/invariants.py (I2, I4 re-implemented; I3 needs the frozen
//    script, which tests/../om2 drivers run on the directory this test can keep).
// 3. A fill whose basis spans do not tile it is never written.
// 4. AGP evidence, tag parsing, the AGP validator, SHA-256.
//
// Usage: test_om2_output [KEEP_DIR]   (KEEP_DIR receives the synthetic run for invariants.py)
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <map>
#include <random>
#include <sstream>
#include <string>
#include <vector>
#include <unistd.h>

#include "test_env.h"
#include "om2_output.h"
#include "seqio.h"

using namespace ts;
using namespace ts::om2;

namespace {

int checks = 0, failures = 0;
void check(bool ok, const std::string& what) {
    ++checks;
    if (!ok) { ++failures; std::fprintf(stderr, "FAIL: %s\n", what.c_str()); }
}

std::mt19937 rng(20260928);
std::string randomSeq(size_t n) {
    std::string s(n, 'A');
    for (char& c : s) c = "ACGT"[rng() % 4];
    return s;
}
char otherBase(char c) { return c == 'A' ? 'C' : c == 'C' ? 'G' : c == 'G' ? 'T' : 'A'; }
std::string upper(std::string s) { for (char& c : s) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c))); return s; }

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

Junction blank(uint32_t id) {
    Junction j{};
    j.id = id;
    j.source = Source::Layout;
    j.a = {0, true};
    j.b = {0, false};
    j.claimedN = j.writtenN = 0;
    j.verdict = Verdict::NotJudged;
    j.tier = Tier::E;
    j.admit = Admit::Scaffold;
    j.endA = j.endB = EndClass::Unknown;
    j.copyA = j.copyB = 1;
    j.gmin = -1;
    j.exhaustive = j.hairpin = false;
    j.uA = j.uB = 0;
    j.pairLambda = 0;
    j.pairK = j.pairContra = 0;
    j.panelSupport = j.panelGenomes = 0;
    j.tieMargin = 0;
    j.pMisjoin = -1;
    j.allowCloseGaps = true;
    j.contigFilled = false;
    return j;
}

std::string counterLine(const SurfaceStats& s) {
    char path[] = "/tmp/om2_counterXXXXXX";
    const int fd = mkstemp(path);
    std::FILE* f = fdopen(fd, "w+");
    logSurfaceCounters(s, f);
    std::fflush(f);
    std::rewind(f);
    char buf[2048] = {0};
    const size_t n = std::fread(buf, 1, sizeof buf - 1, f);
    std::fclose(f);
    std::remove(path);
    return std::string(buf, n);
}

// ---- 1. flags off -----------------------------------------------------------------------
void testFlagsOff(const std::string& dir) {
    std::filesystem::create_directories(dir);
    std::vector<std::string> seqs{randomSeq(5000) + std::string(100, 'N') + randomSeq(4000)};
    std::vector<std::string> names{"NODE_1_length_9100_cov_30.0000_chr"};
    SurfaceInput in;
    in.outDir = dir;
    in.seqs = &seqs;
    in.names = &names;
    in.scaffoldsFile = true;
    std::string js = "sentinel";
    const SurfaceStats st = writeSurface(in, js);
    check(!st.enabled && st.ok, "flags off: disabled and ok");
    check(js.empty(), "flags off: no report.json om2 block");
    check(!std::filesystem::exists(dir + "/genome"), "flags off: no genome/ directory");
    check(!std::filesystem::exists(dir + "/scaffolds.fasta"), "flags off: scaffolds.fasta untouched (none written)");
    const std::string line = counterLine(st);
    check(line.rfind("[om2-out] enabled=0 records=0 chr_records=0 circular=0 rotated=0 genome_only_joins=0 fills=0 "
                     "allocated_bp=0 agp_evidence_gaps=0 detect_score=NA detect_second=NA", 0) == 0,
          "flags off: counter line is zeros: " + line);
}

// ---- 2. the rRNA case ---------------------------------------------------------------------
struct Site { const char* name; uint32_t off, len; };

void testRrnCase(const std::string& dir) {
    std::filesystem::create_directories(dir);
    const size_t OP = 5000;
    const std::string unit = randomSeq(OP);
    // Sites of the repeat unit (forward orientation); >= 40 bp apart so no 31-mer spans two.
    const Site endL{"endL:60", 60, 1}, its{"ITS:1700-2000", 1700, 300}, mid{"mid:3000", 3000, 1},
        endR{"endR:4940", 4940, 1};
    // ITS types: B and C differ from A at six positions each, 40-50 bp apart.
    std::map<char, std::string> itsType;
    itsType['A'] = unit.substr(its.off, its.len);
    itsType['B'] = itsType['A'];
    itsType['C'] = itsType['A'];
    for (uint32_t p : {10u, 50u, 100u, 150u, 200u, 250u}) itsType['B'][p] = otherBase(itsType['B'][p]);
    for (uint32_t p : {20u, 80u, 130u, 180u, 230u, 290u}) itsType['C'][p] = otherBase(otherBase(itsType['C'][p]));
    const char refL = unit[endL.off], altL = otherBase(refL);
    const char refM = unit[mid.off], altM = otherBase(refM);
    const char refR = unit[endR.off], altR = otherBase(refR);
    struct Copy { bool l; char t; bool m; bool r; bool minus; };
    const Copy truth[7] = {{true, 'A', false, false, false}, {false, 'A', false, true, false},
                           {false, 'B', true, false, false}, {true, 'A', false, false, false},
                           {false, 'B', false, false, true}, {false, 'A', false, true, true},
                           {false, 'C', true, false, true}};
    auto build = [&](bool l, char t, bool m, bool r) {
        std::string s = unit;
        s[endL.off] = l ? altL : refL;
        s.replace(its.off, its.len, itsType[t]);
        s[mid.off] = m ? altM : refM;
        s[endR.off] = r ? altR : refR;
        return s;
    };
    std::vector<std::string> trueCopy;
    for (const Copy& c : truth) trueCopy.push_back(build(c.l, c.t, c.m, c.r));
    // What C2 writes at each locus: ends phased by pairs (right), ITS by the prior at loci 0-5
    // (wrong at locus 4: truth B, prior A) and by consensus at locus 6, the mid SNP by
    // consensus everywhere (wrong at loci 2 and 6).
    const char allocIts[7] = {'A', 'A', 'B', 'A', 'A', 'A', 'A'};
    std::vector<std::string> fill;
    for (int i = 0; i < 7; ++i) fill.push_back(build(truth[i].l, allocIts[i], false, truth[i].r));

    // Chromosome: S0 O0 S1 O1 S2 O2 S3a X S3b O3 S4 O4' S5 O5' S6a Y S6b O6' S7, closed by W.
    std::vector<std::string> S;
    for (int i = 0; i < 8; ++i) S.push_back(randomSeq(3000));
    const std::string X = randomSeq(400), Y = randomSeq(50), W = randomSeq(300);
    const std::string S3a = S[3].substr(0, 1500), S3b = S[3].substr(1500), S6a = S[6].substr(0, 1400),
                      S6b = S[6].substr(1400);
    const std::string N5k(OP, 'N');
    std::string chr;
    std::vector<size_t> opPos;   // record position of each operon N-run
    auto add = [&](const std::string& s) { chr += s; };
    add(S[0]); opPos.push_back(chr.size()); add(N5k); add(S[1]); opPos.push_back(chr.size()); add(N5k);
    add(S[2]); opPos.push_back(chr.size()); add(N5k); add(S3a);
    const size_t xPos = chr.size();
    add(std::string(400, 'N')); add(S3b); opPos.push_back(chr.size()); add(N5k); add(S[4]);
    opPos.push_back(chr.size()); add(N5k); add(S[5]); opPos.push_back(chr.size()); add(N5k); add(S6a);
    const size_t yPos = chr.size();
    add(std::string(50, 'N')); add(S6b); opPos.push_back(chr.size()); add(N5k); add(S[7]);
    (void)xPos; (void)yPos;

    const std::string P1 = randomSeq(6000), P2a = randomSeq(2500), P2b = randomSeq(1800);
    const std::string U1 = randomSeq(1200), U2 = randomSeq(900);
    std::vector<std::string> seqs{chr, P1, P2a, P2b, U1 + std::string(100, 'N') + U2};
    std::vector<std::string> names;
    const char* tags[5] = {"_chr", "_plas_1_circular", "_plas_2", "_plas_2", "_unk"};
    for (int i = 0; i < 5; ++i)
        names.push_back("NODE_" + std::to_string(i + 1) + "_length_" + std::to_string(seqs[i].size()) +
                        "_cov_30.0000" + tags[i]);

    // ---- the ledger ----
    Ledger L{};
    for (int64_t& o : L.rotateOffset) o = -1;
    auto flanks = [&](Junction& j, const std::string& rec, size_t pos, size_t len, bool minus) {
        const std::string fl = rec.substr(pos - 32, 32), fr = rec.substr(pos + len, 32);
        if (!minus) { j.flankL32 = fl; j.flankR32 = fr; }
        else { j.flankL32 = revcomp(fr); j.flankR32 = revcomp(fl); }
    };
    const Site sites[4] = {endL, its, mid, endR};
    std::vector<RepeatVariant> variants;
    for (int i = 0; i < 7; ++i) {
        Junction j = blank(100 + i);
        j.source = Source::Layout;
        j.claimedN = j.writtenN = static_cast<int32_t>(OP);
        j.verdict = Verdict::Resize;
        j.tier = Tier::C;
        j.admit = Admit::Scaffold;
        j.endA = j.endB = EndClass::Repeat;
        j.walks = {{static_cast<int32_t>(OP), 3}};
        j.gmin = static_cast<int32_t>(OP);
        j.exhaustive = true;
        j.pMisjoin = 0.02f;
        j.cls = "kpneumoniae|layout|RESIZE|C|N>3000|repeat-repeat|pairs_abstain";
        flanks(j, chr, opPos[i], OP, truth[i].minus);
        j.fillSeq = fill[i];
        const Basis itsBasis = i == 6 ? Basis::Consensus : Basis::PriorAllocated;
        j.fillSpans = {{0, 40, Basis::IsolateRepeatExact},   {40, 41, Basis::PairPhased},
                       {81, 1619, Basis::IsolateRepeatExact}, {1700, 300, itsBasis},
                       {2000, 990, Basis::IsolateRepeatExact}, {2990, 21, Basis::Consensus},
                       {3011, 1909, Basis::IsolateRepeatExact}, {4920, 41, Basis::PairPhased},
                       {4961, 39, Basis::IsolateRepeatExact}};
        L.j.push_back(j);
        for (const Site& s : sites) {
            RepeatVariant v;
            v.family = "rrn";
            v.cluster = "rrn_1";
            v.site = s.name;
            v.junction = 100 + i;
            v.fillOffset = s.off;
            v.siteLen = s.len;
            if (s.off == its.off) {
                v.placed = std::string(1, allocIts[i]);
                v.basis = itsBasis;
                v.alleles = {{"A", 4.f}, {"B", 2.f}, {"C", 1.f}};
                v.phasing = i == 6 ? "none" : "prior";
            } else if (s.off == mid.off) {
                v.placed = std::string(1, refM);
                v.basis = Basis::Consensus;
                v.alleles = {{std::string(1, refM), 5.f}, {std::string(1, altM), 2.f}};
                v.loci = "102;103;104;106;107;108;109";
            } else {
                const char ref = s.off == endL.off ? refL : refR, alt = s.off == endL.off ? altL : altR;
                const bool isAlt = s.off == endL.off ? truth[i].l : truth[i].r;
                v.placed = std::string(1, isAlt ? alt : ref);
                v.basis = Basis::PairPhased;
                v.alleles = {{std::string(1, ref), 5.f}, {std::string(1, alt), 2.f}};
                v.phasing = "pair";
            }
            variants.push_back(v);
        }
    }
    RepeatVariant unplacedSite;   // a minority allele the isolate carries at no placed locus
    unplacedSite.family = "rrn";
    unplacedSite.cluster = "rrn_1";
    unplacedSite.site = "23S:4100";
    unplacedSite.alleles = {{"G", 6.f}, {"T", 1.f}};
    variants.push_back(unplacedSite);
    {   // X: a panel-ordered silent junction, admitted to the genome view only
        Junction j = blank(200);
        j.source = Source::Join;
        j.claimedN = j.writtenN = 400;
        j.verdict = Verdict::Silent;
        j.tier = Tier::E;
        j.admit = Admit::GenomeOnly;
        j.endA = j.endB = EndClass::DeadEnd;
        j.pMisjoin = 0.04f;
        flanks(j, chr, xPos, 400, false);
        L.j.push_back(j);
    }
    {   // the unplaced record's seam: the admission breaks it
        Junction j = blank(201);
        j.source = Source::Join;
        j.claimedN = j.writtenN = 100;
        j.verdict = Verdict::BreakDoubleUse;
        j.tier = Tier::E;
        j.admit = Admit::Break;
        flanks(j, seqs[4], U1.size(), 100, false);
        L.j.push_back(j);
    }
    {   // a junction C1 recorded and the gap filler closed later: not an N-run any more
        Junction j = blank(202);
        j.source = Source::Resolver;
        j.claimedN = j.writtenN = 20;
        j.verdict = Verdict::PassExact;
        j.tier = Tier::B;
        j.flankL32 = randomSeq(32);
        j.flankR32 = randomSeq(32);
        L.j.push_back(j);
    }
    const int64_t dnaA = 16100;   // record + fills coordinates: 100 bp into S2
    {   // the chromosome wrap, closed by a unique graph walk W
        Junction j = blank(300);
        j.source = Source::Wrap;
        j.claimedN = j.writtenN = 300;
        j.verdict = Verdict::PassExact;
        j.tier = Tier::C;
        j.endA = j.endB = EndClass::Unique;
        j.exhaustive = true;
        j.pMisjoin = 0.01f;
        j.flankL32 = chr.substr(chr.size() - 32);
        j.flankR32 = chr.substr(0, 32);
        j.fillSeq = W;
        j.fillSpans = {{0, 300, Basis::GraphWalk}};
        L.j.push_back(j);
        L.rotateOffset[0] = dnaA;
    }

    // ---- the files the release run would have written before the writer runs ----
    std::string err;
    check(writeFasta(dir + "/scaffolds.fasta", seqs, names, 80, err), "write scaffolds.fasta");
    {
        std::vector<std::string> cs, cn;
        std::ofstream agp(dir + "/scaffolds.agp");
        agp << "##agp-version\t2.1\n# release-format AGP\n";
        for (size_t i = 0; i < seqs.size(); ++i) {
            const std::string& q = seqs[i];
            size_t pos = 0, part = 0, piece = 0;
            while (pos < q.size()) {
                size_t e = pos;
                while (e < q.size() && q[e] != 'N') ++e;
                if (e > pos) {
                    ++part; ++piece;
                    agp << names[i] << '\t' << pos + 1 << '\t' << e << '\t' << part << "\tW\t" << names[i] << '_'
                        << piece << "\t1\t" << e - pos << "\t+\n";
                    cs.push_back(q.substr(pos, e - pos));
                    cn.push_back("NODE_c" + std::to_string(cs.size()));
                }
                size_t g = e;
                while (g < q.size() && q[g] == 'N') ++g;
                if (g > e) {
                    ++part;
                    agp << names[i] << '\t' << e + 1 << '\t' << g << '\t' << part << "\tN\t" << g - e
                        << "\tscaffold\tyes\talign_genus\n";
                }
                pos = g > e ? g : e;
            }
        }
        check(writeFasta(dir + "/contigs.fasta", cs, cn, 80, err), "write contigs.fasta");
        std::ofstream gfa(dir + "/assembly_graph.gfa");
        gfa << "H\tVN:Z:1.0\n";
        int n = 0;
        for (const std::string& s : S) gfa << "S\t" << ++n << '\t' << s << '\n';
        for (const std::string& s : trueCopy) gfa << "S\t" << ++n << '\t' << s << '\n';   // operon copies as graph paths
        gfa << "S\t" << ++n << '\t' << S[7].substr(S[7].size() - 40) + W + S[0].substr(0, 40) << '\n';
        for (const std::string& s : {P1, P2a, P2b, U1, U2}) gfa << "S\t" << ++n << '\t' << s << '\n';
    }

    setenv("TESSERACT_OM2_OUTPUT", "1", 1);
    setenv("TESSERACT_OM2_AGP_EVIDENCE", "1", 1);
    SurfaceInput in;
    in.outDir = dir;
    in.seqs = &seqs;
    in.names = &names;
    in.scaffoldsFile = true;
    in.ledger = &L;
    in.variants = &variants;
    in.modelRan = true;
    std::string js;
    const SurfaceStats st = writeSurface(in, js);
    unsetenv("TESSERACT_OM2_OUTPUT");
    unsetenv("TESSERACT_OM2_AGP_EVIDENCE");
    check(st.ok, "rRNA case: writer ok: " + st.error);
    check(st.fills == 8, "8 fills (7 operons + the wrap), got " + std::to_string(st.fills));
    check(st.genomeOnlyJoins == 1 && st.broken == 1 && st.unrecorded == 1, "one genome-only join, one break, one unrecorded seam");
    check(st.circular == 1 && st.rotated == 1 && st.rotateRefused == 0, "the chromosome is closed and rotated");
    check(st.scaffoldSplits == 2, "scaffolds.fasta split at the genome-only and the broken seam");
    check(st.variantRows == 29, "29 repeat-variant rows (7 loci x 4 sites + 1 unplaced site)");
    check(st.agpEvidenceGaps == 8, "scaffolds.agp: evidence on every remaining gap (8), got " +
                                        std::to_string(st.agpEvidenceGaps));
    const size_t allocPerCopy = 300 + 21;   // ITS (prior or consensus) + the mid-SNP window
    check(st.allocatedBp == 7 * allocPerCopy, "allocated bases = 7 x (ITS + mid window)");

    // ---- genome.fasta: order, names, headers, the closed and rotated chromosome ----
    const auto G = readFasta(dir + "/genome/genome.fasta");
    std::vector<std::string> gnames;
    for (const auto& r : G) gnames.push_back(firstToken(r.first));
    const std::vector<std::string> wantNames{"chromosome_1", "plasmid_1", "plasmid_2_1", "plasmid_2_2", "unplaced_1",
                                             "unplaced_2"};
    check(gnames == wantNames, "genome records in order chromosome, plasmid groups, unplaced");
    check(G.size() == 6 && G[0].first.find("topology=circular") != std::string::npos &&
              G[0].first.find("circular_basis=graph") != std::string::npos &&
              G[0].first.find("replicon=chr") != std::string::npos,
          "chromosome header: circular by a graph walk");
    check(G.size() == 6 && G[0].first.find("gaps=2 ") != std::string::npos, "chromosome keeps two open gaps (X, Y)");
    check(G.size() == 6 && G[0].first.find("allocated_bp=" + std::to_string(7 * allocPerCopy)) != std::string::npos,
          "chromosome header allocated_bp");
    check(G.size() == 6 && G[0].first.find("exp_misjoins=NA") != std::string::npos,
          "exp_misjoins is NA when a gap (the unrecorded one) has no calibrated class");
    check(G.size() == 6 && G[1].first.find("topology=linear") != std::string::npos &&
              G[1].first.find("circular_basis=pairs") != std::string::npos,
          "a pairs-only circular plasmid stays linear, circular_basis=pairs");
    // The expected chromosome: the record with fills in place of the operon N-runs, the X and Y
    // gaps kept as N, the wrap fill appended, rotated to 100 bp into S2.
    std::string expect;
    {
        std::string noRot;
        size_t pos = 0;
        for (int i = 0; i < 7; ++i) {
            noRot += chr.substr(pos, opPos[i] - pos);
            noRot += truth[i].minus ? revcomp(fill[i]) : fill[i];
            pos = opPos[i] + OP;
        }
        noRot += chr.substr(pos) + W;
        expect = noRot.substr(dnaA) + noRot.substr(0, dnaA);
    }
    check(G.size() == 6 && upper(G[0].second) == upper(expect), "chromosome sequence = fills in place + wrap, rotated to dnaA");
    check(G.size() == 6 && upper(G[0].second).compare(0, 200, S[2].substr(100, 200)) == 0, "chromosome starts at the dnaA offset");
    // lowercase exactly on allocation
    size_t lower = 0;
    for (char c : G[0].second) lower += std::islower(static_cast<unsigned char>(c)) ? 1 : 0;
    check(lower == 7 * allocPerCopy, "lowercase bases = allocated bases");

    // ---- I2: AGP reconstruction and BED tiling (the invariants.py contract) ----
    check(validateAgp(dir + "/genome/genome.agp").empty(), "genome.agp is valid AGP 2.1: " + validateAgp(dir + "/genome/genome.agp"));
    check(validateAgp(dir + "/scaffolds.agp").empty(), "rewritten scaffolds.agp is valid AGP 2.1: " + validateAgp(dir + "/scaffolds.agp"));
    std::map<std::string, std::string> scaf, fills, genome;
    for (const auto& r : readFasta(dir + "/scaffolds.fasta")) scaf[firstToken(r.first)] = r.second;
    for (const auto& r : readFasta(dir + "/genome/fills.fasta")) fills[firstToken(r.first)] = r.second;
    for (const auto& r : G) genome[firstToken(r.first)] = r.second;
    check(fills.size() == 8, "fills.fasta holds 8 records");
    std::map<std::string, std::string> recon;
    std::map<std::string, std::vector<std::pair<size_t, size_t>>> fillIv;
    bool compsFound = true;
    for (const auto& x : readTable(dir + "/genome/genome.agp")) {
        const size_t ob = std::stoul(x[1]), oe = std::stoul(x[2]);
        if (x[4] == "N" || x[4] == "U") { recon[x[0]] += std::string(oe - ob + 1, 'N'); continue; }
        const bool isFill = x[5].rfind("om2fill_", 0) == 0;
        const auto& src = isFill ? fills : scaf;
        auto it = src.find(x[5]);
        if (it == src.end()) { compsFound = false; continue; }
        std::string s = it->second.substr(std::stoul(x[6]) - 1, std::stoul(x[7]) - std::stoul(x[6]) + 1);
        if (x[8] == "-") s = revcomp(s);
        recon[x[0]] += s;
        if (isFill) fillIv[x[0]].push_back({ob - 1, oe});
    }
    check(compsFound, "every genome.agp component is a scaffolds.fasta or fills.fasta record");
    bool reconOk = recon.size() == genome.size();
    for (const auto& kv : genome) reconOk = reconOk && upper(recon[kv.first]) == upper(kv.second);
    check(reconOk, "I2: genome.fasta == reconstruction from AGP + scaffolds.fasta + fills.fasta");
    std::map<std::string, std::vector<std::vector<std::string>>> bed;
    for (const auto& x : readTable(dir + "/genome/genome.mask.bed")) bed[x[0]].push_back(x);
    bool tiled = true;
    for (const auto& kv : fillIv)
        for (const auto& iv : kv.second) {
            size_t cov = 0;
            for (const auto& b : bed[kv.first]) {
                const size_t s0 = std::stoul(b[1]), e = std::stoul(b[2]);
                if (s0 >= iv.first && e <= iv.second) cov += e - s0;
            }
            tiled = tiled && cov == iv.second - iv.first;
        }
    check(tiled, "I2: every fill component is tiled exactly by basis intervals");
    bool lowerOnlyAlloc = true;
    for (const auto& kv : genome) {
        std::vector<char> alloc(kv.second.size(), 0);
        for (const auto& b : bed[kv.first])
            if (b[3] == "PRIOR_ALLOCATED" || b[3] == "MULTIPLICITY" || b[3] == "CONSENSUS")
                std::fill(alloc.begin() + std::stol(b[1]), alloc.begin() + std::stol(b[2]), 1);
        for (size_t i = 0; i < kv.second.size(); ++i)
            if (std::islower(static_cast<unsigned char>(kv.second[i])) && !alloc[i]) lowerOnlyAlloc = false;
    }
    check(lowerOnlyAlloc, "I2: lowercase only inside allocation intervals");
    // I4: no allocated interval (with 50 bp of context) in contigs.fasta or scaffolds.fasta
    std::vector<std::string> training;
    for (const auto& r : readFasta(dir + "/contigs.fasta")) training.push_back(upper(r.second));
    for (const auto& kv : scaf) training.push_back(upper(kv.second));
    bool hygiene = true;
    size_t allocIntervals = 0;
    for (const auto& kv : bed)
        for (const auto& b : kv.second) {
            if (b[3] != "PRIOR_ALLOCATED" && b[3] != "MULTIPLICITY" && b[3] != "CONSENSUS") continue;
            ++allocIntervals;
            const std::string& g = genome[kv.first];
            const size_t s0 = std::stoul(b[1]), e = std::stoul(b[2]);
            const std::string ctx = upper(g.substr(s0 >= 50 ? s0 - 50 : 0, std::min(g.size(), e + 50) - (s0 >= 50 ? s0 - 50 : 0)));
            for (const std::string& t : training)
                if (t.find(ctx) != std::string::npos || t.find(revcomp(ctx)) != std::string::npos) hygiene = false;
        }
    check(allocIntervals == 14 && hygiene, "I4: 14 allocation intervals, none in contigs.fasta or scaffolds.fasta");
    // basis intervals of the minus-strand copies are mirrored
    bool mirrored = false;
    for (const auto& b : bed["chromosome_1"])
        if (b[4] == "104" && b[3] == "PRIOR_ALLOCATED" && b[5].find(":rc") != std::string::npos) mirrored = true;
    check(mirrored, "minus-strand fill written reverse-complemented with its basis intervals mirrored");

    // ---- scaffolds.fasta: split, never filled ----
    check(scaf.size() == 7, "scaffolds.fasta: 5 records + 2 splits = 7 records");
    check(scaf.count("NODE_1.1_length_" + std::to_string(xPos) + "_cov_30.0000_chr") == 1, "split part named NODE_1.1_length_<part>_...");
    bool noFillInScaf = true;
    for (const auto& kv : scaf)
        for (int i = 0; i < 7; ++i)
            if (upper(kv.second).find(fill[i].substr(1000, 500)) != std::string::npos) noFillInScaf = false;
    check(noFillInScaf, "no genome-view fill reaches scaffolds.fasta");
    {
        // every scaffolds.agp object is a scaffolds.fasta record of the same length
        std::map<std::string, size_t> objEnd;
        for (const auto& x : readTable(dir + "/scaffolds.agp")) objEnd[x[0]] = std::stoul(x[2]);
        bool same = objEnd.size() == scaf.size();
        for (const auto& kv : objEnd) same = same && scaf.count(kv.first) && scaf[kv.first].size() == kv.second;
        check(same, "scaffolds.agp objects are exactly the split scaffolds.fasta records");
        const std::string agpText = slurp(dir + "/scaffolds.agp");
        check(agpText.find("\trepeat\tyes\talign_genus;map\n") != std::string::npos,
              "scaffolds.agp: an operon gap carries gap type repeat and evidence align_genus;map");
        check(agpText.find("\tscaffold\tyes\tunspecified\n") != std::string::npos,
              "scaffolds.agp: the unrecorded seam keeps `unspecified` (a model ran)");
    }

    // ---- junctions.tsv ----
    std::map<std::string, int> status;
    const auto jrows = readTable(dir + "/genome/junctions.tsv");
    for (size_t i = 1; i < jrows.size(); ++i) ++status[jrows[i][1]];
    check(status["filled_genome"] == 7 && status["wrap_closed"] == 1 && status["open_gap"] == 1 &&
              status["unrecorded_gap"] == 1 && status["broken"] == 1 && status["closed_or_absent"] == 1,
          "junctions.tsv statuses: 7 filled, wrap closed, 1 open, 1 unrecorded, 1 broken, 1 closed elsewhere");
    bool widths = true;
    for (const auto& r : jrows) widths = widths && r.size() == jrows[0].size();
    check(widths && jrows[0].size() == 46, "junctions.tsv: 46 columns on every row");

    // ---- repeat_variants.tsv: the inter-copy diversity, and where each placed allele sits ----
    const auto rv = readTable(dir + "/genome/repeat_variants.tsv");
    check(rv.size() == 30, "repeat_variants.tsv: header + 29 rows");
    bool placedRight = true, itsCarriers = true;
    size_t snpRows = 0;
    for (size_t i = 1; i < rv.size(); ++i) {
        const auto& x = rv[i];
        if (x[2] == "ITS:1700-2000") itsCarriers = itsCarriers && x[10] == "A:4.00;B:2.00;C:1.00" && x[11] == "3";
        if (x[5] == "." || x[6] != "1") continue;
        ++snpRows;
        const std::string& g = genome[x[4]];
        const char base = static_cast<char>(std::toupper(static_cast<unsigned char>(g[std::stoul(x[5])])));
        const char want = x[7] == "+" ? x[8][0] : revcomp(x[8])[0];
        placedRight = placedRight && base == want;
    }
    check(snpRows == 21 && placedRight, "every placed SNP allele is at its genome_pos0, on either strand (21 rows)");
    check(itsCarriers, "ITS rows carry all three types with carriers 4/2/1");
    const std::string rvText = slurp(dir + "/genome/repeat_variants.tsv");
    check(rvText.find("rrn\trrn_1\t23S:4100\t.\t.\t.\t1\t.\t.\t.\tG:6.00;T:1.00\t2\t7.00\tnone\t.") != std::string::npos,
          "a site placed at no locus is still reported with its alleles");

    // ---- closure.txt, README, report block ----
    const std::string closure = slurp(dir + "/genome/closure.txt");
    check(closure.find("chromosome_1\tlength=") == 0 && closure.find("topology=circular\tstatus=open") != std::string::npos,
          "closure.txt: chromosome circular but open (two gaps)");
    check(closure.find("1 gap dead_end") != std::string::npos && closure.find("1 gap unrecorded") != std::string::npos,
          "closure.txt: gaps listed by class");
    check(closure.find("joined by read pairs") != std::string::npos, "closure.txt: the pairs-only plasmid wrap is explained");
    check(slurp(dir + "/genome/README.txt").find("PRIOR_ALLOCATED") != std::string::npos, "README.txt explains the labels");
    check(js.find("\"ledger\": true") != std::string::npos && js.find("\"circular\": 1") != std::string::npos &&
              js.find("\"allocated_bp\": " + std::to_string(7 * allocPerCopy)) != std::string::npos,
          "report.json om2 block carries the ledger flag and counters");
    check(js.front() == '{' && js.substr(js.size() - 5) == "    }", "om2 block is one JSON object at the report's indentation");
    int depth = 0;
    bool balanced = true;
    for (char c : js) { if (c == '{' || c == '[') ++depth; if (c == '}' || c == ']') --depth; balanced = balanced && depth >= 0; }
    check(balanced && depth == 0, "om2 block brackets balance");
}

// ---- 3. a fill whose spans do not tile it is never written ----------------------------------
void testBadFill(const std::string& dir) {
    std::filesystem::create_directories(dir);
    const std::string a = randomSeq(2000), b = randomSeq(2000), f = randomSeq(500);
    std::vector<std::string> seqs{a + std::string(500, 'N') + b};
    std::vector<std::string> names{"NODE_1_length_4500_cov_20.0000_chr"};
    std::string err;
    writeFasta(dir + "/scaffolds.fasta", seqs, names, 80, err);
    Ledger L{};
    for (int64_t& o : L.rotateOffset) o = -1;
    Junction j = blank(7);
    j.claimedN = j.writtenN = 500;
    j.flankL32 = a.substr(a.size() - 32);
    j.flankR32 = b.substr(0, 32);
    j.fillSeq = f;
    j.fillSpans = {{0, 200, Basis::GraphWalk}, {250, 250, Basis::Consensus}};   // a hole at 200-250
    L.j.push_back(j);
    setenv("TESSERACT_OM2_OUTPUT", "1", 1);
    SurfaceInput in;
    in.outDir = dir;
    in.seqs = &seqs;
    in.names = &names;
    in.scaffoldsFile = true;
    in.ledger = &L;
    std::string js;
    const SurfaceStats st = writeSurface(in, js);
    unsetenv("TESSERACT_OM2_OUTPUT");
    check(st.ok && st.fills == 0, "a fill with a hole in its basis tiling is not written");
    const auto G = readFasta(dir + "/genome/genome.fasta");
    check(G.size() == 1 && G[0].second == seqs[0], "its gap stays N in the genome view");
}

// ---- 4. pieces ---------------------------------------------------------------------------------
void testPieces(const std::string& dir) {
    Junction j = blank(1);
    j.source = Source::Resolver;
    check(agpEvidence(&j, true) == "paired-ends", "resolver gap: paired-ends");
    j.source = Source::JoinButt1;
    j.tier = Tier::E;
    check(agpEvidence(&j, true) == "align_genus", "model join, tier E: align_genus");
    j.tier = Tier::B;
    check(agpEvidence(&j, true) == "paired-ends;align_genus", "model join confirmed by pairs");
    j.tier = Tier::C;
    check(agpEvidence(&j, true) == "align_genus;map", "model join, graph tier: align_genus;map");
    check(agpEvidence(nullptr, false) == "paired-ends" && agpEvidence(nullptr, true) == "unspecified",
          "unrecorded gap: the release's rule");
    TagInfo t = parseTag("NODE_3_length_100_cov_12.5000_plas_4_circular");
    check(t.cls == 'p' && t.group == 4 && t.circular, "tag _plas_4_circular");
    t = parseTag("NODE_1_length_100_cov_12.5000_chr");
    check(t.cls == 'c' && !t.circular, "tag _chr");
    t = parseTag("NODE_9_length_100_cov_1.0000_plas_circular");
    check(t.cls == 'p' && t.group == 0 && t.circular, "tag _plas_circular (ungrouped)");
    t = parseTag("NODE_9_length_100_cov_1.0000_unk");
    check(t.cls == 'u', "tag _unk");
    std::filesystem::create_directories(dir);
    {
        std::ofstream o(dir + "/bad.agp");
        o << "##agp-version\t2.1\nr1\t1\t100\t1\tW\tc1\t1\t100\t+\nr1\t102\t200\t2\tW\tc2\t1\t99\t+\n";
    }
    check(!validateAgp(dir + "/bad.agp").empty(), "AGP validator rejects non-contiguous coordinates");
    {
        std::ofstream o(dir + "/bad2.agp");
        o << "##agp-version\t2.1\nr1\t1\t100\t1\tW\tc1\t1\t100\t+\nr1\t101\t150\t2\tU\t50\tscaffold\tyes\tpaired-ends\n"
             "r1\t151\t200\t3\tW\tc2\t1\t50\t+\n";
    }
    check(!validateAgp(dir + "/bad2.agp").empty(), "AGP validator rejects a U gap that is not 100 bp");
    {
        std::ofstream o(dir + "/abc.txt");
        o << "abc";
    }
    check(sha256File(dir + "/abc.txt") == "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad", "SHA-256 of abc");
    check(md5File(dir + "/abc.txt") == "900150983cd24fb0d6963f7d28e17f72", "MD5 of abc");
    check(revcomp("AACGTn") == "nACGTT", "revcomp keeps case and N");
}

}  // namespace

int main(int argc, char** argv) {
    testenv::clearTesseractEnv();
    const std::string root = argc > 1 ? std::string(argv[1])
                                      : (std::filesystem::temp_directory_path() /
                                         ("tesseract-om2-output-" + std::to_string(getpid()))).string();
    std::filesystem::remove_all(root);
    testFlagsOff(root + "/off");
    testRrnCase(root + "/rrn");
    testBadFill(root + "/badfill");
    testPieces(root + "/pieces");
    if (argc <= 1) std::filesystem::remove_all(root);
    std::printf("test_om2_output: %d checks, %d failures%s\n", checks, failures,
                argc > 1 ? (" (kept " + root + ")").c_str() : "");
    return failures ? 1 : 0;
}
