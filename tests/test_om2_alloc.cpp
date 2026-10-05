// Organism Model 2.0, C2: allocation rules, the sidecar loaders and their md5 pin.
//
// The allocation tests are the rules of DESIGN.md section 3.3 on hand-built sites: an 8-copy rRNA
// operon with three ITS types (4/3/1 carriers) and a pair-reachable end variant. What must hold:
//   - flank-anchored evidence phases a locus only at posterior >= 0.99 (two agreeing fragments);
//   - the locus prior is used only on ITS-class sites, and is labelled PRIOR_ALLOCATED only where
//     the choice FLIPS without it (the counterfactual is recorded);
//   - the budget drives MULTIPLICITY; an exhausted budget falls back to the majority branch;
//   - consensus mode pastes the majority everywhere; off allocates nothing.
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <vector>

#include "om2_alloc.h"
#include "test_env.h"

using namespace ts::om2;

namespace {
int checks = 0;
void check(bool ok, const std::string& why) {
    ++checks;
    if (!ok) throw std::runtime_error(why);
}

SiteInput itsSite() {
    SiteInput s;
    s.itsClass = true;
    s.branches.resize(3);
    s.branches[0].carriers = 4; s.branches[0].cov = 160;   // ITS type A
    s.branches[1].carriers = 3; s.branches[1].cov = 120;   // ITS type B
    s.branches[2].carriers = 1; s.branches[2].cov = 40;    // ITS type C
    for (uint32_t j = 0; j < 8; ++j) {
        SiteLocus L;
        L.junction = j;
        L.threadN.assign(3, 0);
        L.pairN.assign(3, 0);
        s.loci.push_back(L);
    }
    return s;
}

// Truth for the 8 loci: A A A A B B B C, in junction order.
const int kTruth[8] = {0, 0, 0, 0, 1, 1, 1, 2};

void testPosterior() {
    check(std::fabs(phasePosterior({1, 0}, 0, 0.02) - 0.98) < 1e-9, "one fragment is not enough (0.98)");
    check(phasePosterior({2, 0}, 0, 0.02) >= 0.99, "two agreeing fragments phase");
    check(phasePosterior({5, 1}, 0, 0.02) >= 0.99, "5 vs 1 phases");
    check(phasePosterior({3, 3}, 0, 0.02) < 0.6, "a tie does not phase");
    check(phasePosterior({2, 0, 0}, 0, 0.02) >= 0.99, "three branches, two agreeing");
}

void testConsensusFloor() {
    SiteInput s = itsSite();
    AllocParams p;
    p.mode = AllocMode::Phased;
    AllocCounts c;
    const auto ch = allocateSite(s, p, c);
    int correct = 0;
    for (size_t i = 0; i < 8; ++i) {
        check(ch[i].branch == 0, "no evidence, no prior: the majority branch at every locus");
        check(ch[i].basis == Basis::Consensus, "labelled CONSENSUS");
        correct += ch[i].branch == kTruth[i];
    }
    check(correct == 4, "consensus paste is right at 4/8 ITS loci");
    check(c.consensus == 8 && c.prior == 0, "counts");
}

void testPriorFlip() {
    SiteInput s = itsSite();
    for (size_t i = 0; i < 8; ++i) {
        s.loci[i].prior.assign(3, 0.0025);
        s.loci[i].prior[static_cast<size_t>(kTruth[i])] = 0.995;
        s.loci[i].priorGenomes = 200;
    }
    AllocParams p;
    p.mode = AllocMode::Prior;
    AllocCounts c;
    const auto ch = allocateSite(s, p, c);
    for (size_t i = 0; i < 8; ++i) {
        check(ch[i].branch == kTruth[i], "the prior allocates the true ITS type at locus " + std::to_string(i));
        const bool flipped = ch[i].noPriorBranch != ch[i].branch;
        check(ch[i].priorFlip == flipped, "priorFlip records exactly the counterfactual flip");
        check((ch[i].basis == Basis::PriorAllocated) == flipped,
              "PRIOR_ALLOCATED only where the no-prior choice differs (locus " + std::to_string(i) + ")");
    }
    // Order 0..7 (equal strength). The A loci agree with the no-prior choice (majority); the three B
    // loci flip; by the C locus only C is left in budget, so multiplicity already names it: no flip.
    check(c.prior == 3, "the three B loci flip, got " + std::to_string(c.prior));
    check(c.priorNoFlip == 5, "A loci and the budget-forced C locus agree without the prior");
    check(ch[7].basis == Basis::Multiplicity, "the C locus is MULTIPLICITY, not PRIOR_ALLOCATED");

    // Too few genomes, or a sub-0.99 prior: never used.
    SiteInput weak = s;
    for (auto& L : weak.loci) L.priorGenomes = 9;
    AllocCounts c2;
    const auto ch2 = allocateSite(weak, p, c2);
    for (const auto& x : ch2) check(x.basis != Basis::PriorAllocated, "fewer than 10 genomes: no prior");
    SiteInput soft = s;
    for (size_t i = 0; i < 8; ++i) soft.loci[i].prior[static_cast<size_t>(kTruth[i])] = 0.98;
    AllocCounts c3;
    for (const auto& x : allocateSite(soft, p, c3)) check(x.basis != Basis::PriorAllocated, "P < 0.99: no prior");

    // Not an ITS-class site: the prior is never consulted.
    SiteInput snv = s;
    snv.itsClass = false;
    AllocCounts c4;
    for (const auto& x : allocateSite(snv, p, c4)) check(x.basis != Basis::PriorAllocated, "SNV site: no prior");

    // Phased mode ignores the prior.
    AllocParams ph;
    ph.mode = AllocMode::Phased;
    AllocCounts c5;
    for (const auto& x : allocateSite(s, ph, c5)) check(x.basis != Basis::PriorAllocated, "phased mode: no prior");
}

void testClassPrior() {
    // The real shape: an ITS bubble with 5 branches, two ITS classes carrying SNV sub-variants.
    // Classes: branches 0,1 = class 0 (5 copies); 2,3 = class 1 (2 copies); 4 = unknown (1 copy).
    SiteInput s;
    s.itsClass = true;
    s.branches.resize(5);
    const int carriers[5] = {3, 2, 1, 1, 1};
    for (size_t b = 0; b < 5; ++b) { s.branches[b].carriers = carriers[b]; s.branches[b].cov = 40.0 * carriers[b]; }
    s.branchClass = {0, 0, 1, 1, -1};
    for (uint32_t j = 0; j < 4; ++j) {
        SiteLocus L;
        L.junction = j;
        L.threadN.assign(5, 0);
        L.pairN.assign(5, 0);
        L.priorGenomes = 1600;
        L.prior = j == 3 ? std::vector<double>{0.004, 0.994} : std::vector<double>{0.994, 0.004};
        s.loci.push_back(L);
    }
    AllocParams p;
    p.mode = AllocMode::Prior;
    AllocCounts c;
    const auto ch = allocateSite(s, p, c);
    for (size_t i = 0; i < 3; ++i) {
        check(s.branchClass[static_cast<size_t>(ch[i].branch)] == 0, "class-0 loci get a class-0 branch");
        check(ch[i].basis != Basis::PriorAllocated, "the majority is already class 0: no flip");
    }
    check(s.branchClass[static_cast<size_t>(ch[3].branch)] == 1, "the class-1 locus gets a class-1 branch");
    check(ch[3].basis == Basis::PriorAllocated && ch[3].priorFlip, "and it is a flip: PRIOR_ALLOCATED");
    check(s.branchClass[static_cast<size_t>(ch[3].noPriorBranch)] != 1, "the counterfactual names another class");
    check(c.prior == 1 && c.priorNoFlip == 3, "one flip, three agreements");
    // a class whose branches are all spent in the budget: the prior abstains
    SiteInput t = s;
    t.branches[2].carriers = 0;
    t.branches[3].carriers = 0;
    AllocCounts c2;
    const auto ch2 = allocateSite(t, p, c2);
    check(ch2[3].basis != Basis::PriorAllocated, "no class-1 branch in budget: the prior abstains");
}

void testPhasedEndVariant() {
    // A pair-reachable end variant, 2 alleles, 5/3 carriers; loci 0-5 have linking fragments.
    SiteInput s;
    s.branches.resize(2);
    s.branches[0].carriers = 5; s.branches[0].cov = 200;
    s.branches[1].carriers = 3; s.branches[1].cov = 120;
    const int truth[8] = {0, 1, 0, 0, 1, 0, 0, 1};
    for (uint32_t j = 0; j < 8; ++j) {
        SiteLocus L;
        L.junction = j;
        L.threadN.assign(2, 0);
        L.pairN.assign(2, 0);
        if (j < 3) L.threadN[static_cast<size_t>(truth[j])] = 4;          // reads through the flank
        else if (j < 6) L.pairN[static_cast<size_t>(truth[j])] = 3;       // mates in the flank
        s.loci.push_back(L);
    }
    s.loci[5].pairN[1] = 1;   // one stray fragment: 3 vs 1 still phases
    AllocParams p;
    AllocCounts c;
    const auto ch = allocateSite(s, p, c);
    for (size_t i = 0; i < 3; ++i) {
        check(ch[i].branch == truth[i] && ch[i].basis == Basis::ThreadPhased, "thread-phased locus");
    }
    for (size_t i = 3; i < 6; ++i) {
        check(ch[i].branch == truth[i] && ch[i].basis == Basis::PairPhased, "pair-phased locus");
    }
    // Phased: 0,1,0,0,1,0 -> budget left 5-4=1 of allele 0, 3-2=1 of allele 1: marginal 0.5 -> majority.
    check(ch[6].basis == Basis::Consensus && ch[6].branch == 0, "unphased locus 6: majority");
    // Locus 7: budget now 0 of allele 0 and 1 of allele 1 -> forced by multiplicity.
    check(ch[7].branch == 1 && ch[7].basis == Basis::Multiplicity, "unphased locus 7: MULTIPLICITY");
    check(c.thread == 3 && c.pair == 3 && c.multiplicity == 1 && c.consensus == 1, "counts");
}

void testBudgetExhaustion() {
    SiteInput s;
    s.branches.resize(2);
    s.branches[0].carriers = 1; s.branches[0].cov = 40;
    s.branches[1].carriers = 1; s.branches[1].cov = 38;
    for (uint32_t j = 0; j < 4; ++j) {
        SiteLocus L;
        L.junction = j;
        L.threadN.assign(2, 0);
        L.pairN.assign(2, 0);
        s.loci.push_back(L);
    }
    AllocParams p;
    AllocCounts c;
    const auto ch = allocateSite(s, p, c);
    check(ch[0].branch == 0 && ch[0].basis == Basis::Consensus, "first: majority (0.5 marginal)");
    check(ch[1].branch == 1 && ch[1].basis == Basis::Multiplicity, "second: the only branch left in budget");
    check(ch[2].branch == 0 && ch[3].branch == 0, "budget spent: majority");
    check(c.budgetExhausted == 2, "two loci past the budget");
}

void testModes() {
    SiteInput s = itsSite();
    s.loci[5].threadN[1] = 9;
    AllocParams p;
    p.mode = AllocMode::Consensus;
    AllocCounts c;
    for (const auto& x : allocateSite(s, p, c)) check(x.branch == 0 && x.basis == Basis::Consensus, "consensus mode");
    p.mode = AllocMode::Off;
    AllocCounts c2;
    for (const auto& x : allocateSite(s, p, c2)) check(x.branch == -1, "off allocates nothing");
    SiteInput cx = itsSite();
    cx.complex = true;
    AllocParams ph;
    AllocCounts c3;
    for (const auto& x : allocateSite(cx, ph, c3)) check(x.basis == Basis::Consensus, "complex site: consensus only");
}

struct Temp {
    std::filesystem::path dir =
        std::filesystem::temp_directory_path() / ("tesseract-om2alloc-" + std::to_string(getpid()));
    Temp() { std::filesystem::create_directories(dir); }
    ~Temp() { std::filesystem::remove_all(dir); }
    std::string write(const std::string& name, const std::string& text) const {
        const std::string p = (dir / name).string();
        std::ofstream(p, std::ios::binary) << text;
        return p;
    }
};

void testMd5() {
    Temp t;
    check(md5File(t.write("e", "")) == "d41d8cd98f00b204e9800998ecf8427e", "md5 of empty");
    check(md5File(t.write("a", "abc")) == "900150983cd24fb0d6963f7d28e17f72", "md5 of abc");
    check(md5File(t.write("f", "The quick brown fox jumps over the lazy dog")) == "9e107d9d372bb6826bd81d3542a419d6",
          "md5 of the fox");
    check(md5File(t.write("n", "12345678901234567890123456789012345678901234567890123456789012345678901234567890")) ==
              "57edf4a22be3c955ac49da2e2107b67a",
          "md5 over two blocks");
    check(md5File((t.dir / "missing").string()).empty(), "md5 of a missing file is empty");
}

uint64_t canon31(const std::string& s) {
    uint64_t f = 0, r = 0;
    for (size_t i = 0; i < 31; ++i) {
        const int c = std::string("ACGT").find(s[i]);
        f = (f << 2) | static_cast<uint64_t>(c);
        r = (r >> 2) | (static_cast<uint64_t>(3 - c) << 60);
    }
    return f <= r ? f : r;
}

std::string rc(const std::string& s) {
    std::string o(s.rbegin(), s.rend());
    for (char& c : o) c = c == 'A' ? 'T' : c == 'C' ? 'G' : c == 'G' ? 'C' : 'A';
    return o;
}

void testSidecars() {
    Temp t;
    const std::string model = t.write("model.tsm", "not really a model");
    const std::string md5 = md5File(model);
    // a dnaA sketch from a 1.5 kb "chromosome start"
    std::string chr;
    unsigned x = 12345;
    for (int i = 0; i < 6000; ++i) { x = x * 1103515245u + 12345u; chr += "ACGT"[(x >> 16) & 3]; }
    std::string sk = "#om2dnaa v1\n#organism test\n#tsm_md5 " + md5 + "\n#genomes 50\n#window 1500\n";
    for (int p = 0; p + 31 <= 1500; ++p) {
        const std::string km = chr.substr(static_cast<size_t>(p), 31);
        const uint64_t c = canon31(km);
        const uint64_t f = canon31(km) == [&] {
            uint64_t v = 0;
            for (char ch : km) v = (v << 2) | static_cast<uint64_t>(std::string("ACGT").find(ch));
            return v;
        }() ? 0 : 1;
        char line[80];
        std::snprintf(line, sizeof line, "%016llx %d %llu 40\n", static_cast<unsigned long long>(sidecarHash(c)), p,
                      static_cast<unsigned long long>(f));
        sk += line;
    }
    const std::string skPath = t.write("s.om2dnaa", sk);
    DnaaSketch d;
    std::string err;
    check(d.load(skPath, md5, err), "sketch loads with the right pin: " + err);
    // record = the chromosome rotated so that its start sits at 2000
    const std::string rec = chr.substr(4000) + chr.substr(0, 4000);
    const DnaaSketch::Hit h = d.locate(rec);
    check(h.found && !h.reverse && h.offset == 2000, "dnaA located at 2000 (+), got " + std::to_string(h.offset));
    const DnaaSketch::Hit hr = d.locate(rc(rec));
    // on the reverse record, dnaA's first base sits at L-1-2000
    check(hr.found && hr.reverse && hr.offset == static_cast<int64_t>(rec.size()) - 1 - 2000,
          "dnaA located on the reverse strand");
    DnaaSketch bad;
    check(!bad.load(skPath, "0123456789abcdef0123456789abcdef", err), "sketch refused for another model");
    check(err.find("md5") != std::string::npos, "the refusal names the md5 pin");
    check(!bad.load(skPath, "", err), "no model: the pin cannot be checked, refused");
    check(!bad.load(t.write("nopin", "#om2dnaa v1\n0000000000000001 0 0 1\n"), md5, err), "no pin: refused");

    // rRNA prior: two classes, one locus
    std::string its1, its2, fl;
    for (int i = 0; i < 200; ++i) { x = x * 1103515245u + 12345u; its1 += "ACGT"[(x >> 16) & 3]; }
    for (int i = 0; i < 300; ++i) { x = x * 1103515245u + 12345u; its2 += "ACGT"[(x >> 16) & 3]; }
    for (int i = 0; i < 3000; ++i) { x = x * 1103515245u + 12345u; fl += "ACGT"[(x >> 16) & 3]; }
    auto hashes = [&](const std::string& s, uint64_t denom) {
        std::string out;
        for (size_t p = 0; p + 31 <= s.size(); ++p) {
            const uint64_t h = sidecarHash(canon31(s.substr(p, 31)));
            if (denom && h > ~0ULL / denom) continue;
            char b[24];
            std::snprintf(b, sizeof b, "%s%016llx", out.empty() ? "" : ",", static_cast<unsigned long long>(h));
            out += b;
        }
        return out;
    };
    std::string pr = "#om2rrn v1\n#organism test\n#tsm_md5 " + md5 + "\n#genomes 200\n#flank_denom 16\n#flank_bp 1500\n";
    pr += "C 0 0 " + hashes(its1, 0) + "\nC 1 0 " + hashes(its2, 0) + "\n";
    pr += "L 0 200 0:199;1:1 " + hashes(fl, 16) + "\n";
    const std::string prPath = t.write("p.om2rrn", pr);
    RrnPrior rp;
    check(rp.load(prPath, md5, err), "prior loads: " + err);
    check(rp.classOf("GG" + its1 + "TT") == 0 && rp.classOf(its2) == 1, "ITS classes by discriminating hashes");
    check(rp.classOf(fl.substr(0, 400)) == -1, "no class for unrelated sequence");
    std::vector<uint64_t> mk;
    for (size_t p = 0; p + 31 <= 1200; ++p) {
        const uint64_t h = sidecarHash(canon31(fl.substr(p, 31)));
        if (h <= ~0ULL / 16) mk.push_back(h);
    }
    check(rp.matchLocus(mk) == 0, "locus found from its flank markers");
    check(rp.locusGenomes(0) == 200 && rp.classCount(0, 0) == 199 && rp.classCount(0, 1) == 1, "locus counts");
    RrnPrior rp2;
    check(!rp2.load(prPath, "ffffffffffffffffffffffffffffffff", err), "prior refused for another model");
}

}  // namespace

int main() {
    testenv::clearTesseractEnv();
    try {
        testPosterior();
        testConsensusFloor();
        testPriorFlip();
        testClassPrior();
        testPhasedEndVariant();
        testBudgetExhaustion();
        testModes();
        testMd5();
        testSidecars();
    } catch (const std::exception& e) {
        std::fprintf(stderr, "test_om2_alloc FAILED: %s\n", e.what());
        return 1;
    }
    std::printf("test_om2_alloc: %d checks passed\n", checks);
    return 0;
}
