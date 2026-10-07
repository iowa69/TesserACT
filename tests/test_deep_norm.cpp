// p2/plasmid: DEEP_NORM read thinning (src/deep_norm.{h,cpp}).
//
// 1. Flags: unset = off with the b3b6d41 defaults; TESSERACT_DEEP_NORM_MATE=both selects the
//    dn7 rule, a number the b3b6d41 rule with that fraction.
// 2. Decision: the dn7 rule needs BOTH mates at the bar and keeps with probability from the
//    LOWER median; the b3b6d41 rule needs one mate at the bar and the other at F x it, and
//    keeps from the HIGHER median; the draw is the splitmix64 hash of the pair index.
// 3. Apply on real reads: off touches nothing; on, only deep pairs are selected, thinned
//    reads are masked whole with their mates, and the medians do not depend on the threads.
// 4. The counter line, idle and active.
#include "deep_norm.h"
#include "counter.h"
#include "kmer.h"
#include "seqio.h"
#include "test_env.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <random>
#include <string>
#include <unistd.h>
#include <vector>

namespace {

int checks = 0;
void check(bool ok, const std::string& what) {
    ++checks;
    if (!ok) { std::fprintf(stderr, "FAIL: %s\n", what.c_str()); std::exit(1); }
}

double drawOf(uint64_t unit) {   // the keep draw of the two source binaries
    return static_cast<double>(ts::mix64(unit ^ 0x5a17f1c3ULL) >> 11) / 9007199254740992.0;
}

void flags() {
    ts::DeepNormConfig c = ts::deepNormConfigFromEnv();
    check(!c.enabled() && c.sel == 0.0, "unset: off");
    check(c.target == 2.0 && !c.both && c.mate == 1.0 / 6.0, "unset: b3b6d41 defaults");
    setenv("TESSERACT_DEEP_NORM", "30", 1);
    setenv("TESSERACT_DEEP_NORM_MATE", "both", 1);
    c = ts::deepNormConfigFromEnv();
    check(c.enabled() && c.sel == 30.0 && c.both && c.target == 2.0, "DEEP_NORM=30 MATE=both: dn7");
    setenv("TESSERACT_DEEP_NORM_MATE", "0.16666666666666666", 1);
    c = ts::deepNormConfigFromEnv();
    check(!c.both && c.mate == 1.0 / 6.0, "MATE=0.16666666666666666 is exactly the 1/6 default");
    setenv("TESSERACT_DEEP_NORM_MATE", "0.25", 1);
    setenv("TESSERACT_DEEP_NORM_TARGET", "3", 1);
    c = ts::deepNormConfigFromEnv();
    check(!c.both && c.mate == 0.25 && c.target == 3.0, "MATE=0.25 TARGET=3");
    unsetenv("TESSERACT_DEEP_NORM");
    unsetenv("TESSERACT_DEEP_NORM_MATE");
    unsetenv("TESSERACT_DEEP_NORM_TARGET");
}

void decision() {
    using V = ts::DeepNormVerdict;
    ts::DeepNormConfig dn7; dn7.sel = 30; dn7.both = true;
    ts::DeepNormConfig dn8; dn8.sel = 30;                       // mate 1/6
    const double peak = 100;                                   // bar = 3000
    // One mate below the bar: dn7 never selects; dn8 selects when the other is >= 500.
    for (uint64_t u = 0; u < 200; ++u) {
        check(ts::deepNormDecide(dn7, peak, 9000.0f, 2999.0f, true, u) == V::Skip, "dn7: low mate skips");
        check(ts::deepNormDecide(dn7, peak, 2999.0f, 9000.0f, true, u) == V::Skip, "dn7: low mate skips (odd)");
        check(ts::deepNormDecide(dn8, peak, 9000.0f, 2999.0f, true, u) != V::Skip, "dn8: other mate >= 1/6 selects");
        check(ts::deepNormDecide(dn8, peak, 9000.0f, 499.0f, true, u) == V::Skip, "dn8: other mate < 1/6 skips");
        check(ts::deepNormDecide(dn8, peak, 2999.0f, 2999.0f, true, u) == V::Skip, "dn8: no mate at the bar");
    }
    // At the bar is selected (the test is `< bar`).
    bool sawSelect = false;
    for (uint64_t u = 0; u < 50; ++u) sawSelect |= ts::deepNormDecide(dn7, peak, 3000.0f, 3000.0f, true, u) != V::Skip;
    check(sawSelect, "dn7: both mates exactly at the bar are selected");
    // The keep probability: dn7 from the lower median (200/3000), dn8 from the higher
    // (200/60000). Every unit whose draw lies between the two is kept by dn7, thinned by dn8.
    const double keep7 = 2.0 * peak / 3000.0, keep8 = 2.0 * peak / 60000.0;
    size_t between = 0, agree = 0;
    for (uint64_t u = 0; u < 20000; ++u) {
        const double d = drawOf(u);
        const V v7 = ts::deepNormDecide(dn7, peak, 3000.0f, 60000.0f, true, u);
        const V v8 = ts::deepNormDecide(dn8, peak, 60000.0f, 3000.0f, true, u);
        check(v7 == (d < keep7 ? V::Keep : V::Thin), "dn7 verdict is draw < target x peak / lower");
        check(v8 == (d < keep8 ? V::Keep : V::Thin), "dn8 verdict is draw < target x peak / higher");
        if (d >= keep8 && d < keep7) { ++between; check(v7 == V::Keep && v8 == V::Thin, "lower vs higher"); }
        agree += v7 == v8;
    }
    check(between > 1000, "the lower-median rule keeps more pairs (" + std::to_string(between) + ")");
    check(agree > 0, "the rules agree where the draw is outside both");
    // A single-end read: medA alone, the read index as the unit.
    for (uint64_t u = 0; u < 200; ++u) {
        const V v = ts::deepNormDecide(dn7, peak, 5000.0f, 0.0f, false, u);
        check(v == (drawOf(u) < 2.0 * peak / 5000.0 ? V::Keep : V::Thin), "single read: own median");
        check(ts::deepNormDecide(dn7, peak, 2000.0f, 9e9f, false, u) == V::Skip, "single read below the bar");
        check(ts::deepNormDecide(dn8, peak, 5000.0f, 0.0f, false, u) == v, "single read: both rules agree");
    }
    // The target: 0 thins every selected unit.
    ts::DeepNormConfig zero = dn7; zero.target = 0;
    for (uint64_t u = 0; u < 200; ++u)
        check(ts::deepNormDecide(zero, peak, 5000.0f, 5000.0f, true, u) == V::Thin, "target 0 thins all");
}

std::mt19937 rng(61);
std::string randomSeq(size_t n) {
    std::string s(n, 'A');
    for (char& c : s) c = "ACGT"[rng() % 4];
    return s;
}
std::string rc(std::string s) {
    std::reverse(s.begin(), s.end());
    for (char& c : s) c = ts::codeBase(3 - ts::baseCode(c));
    return s;
}
void putKmers(ts::KmerTable& t, const std::string& s, int k, uint32_t count) {
    ts::Kmer f = 0;
    for (size_t i = 0; i < s.size(); ++i) {
        f = ts::pushBack(f, ts::baseCode(s[i]), k);
        if (i + 1 >= static_cast<size_t>(k)) t.put(ts::canonical(f, k), count);
    }
}
bool fullyMasked(const ts::SequenceStore& r, size_t i) {
    for (uint32_t p = 0; p < r.length(i); ++p) if (r.baseAt(i, p) >= 0) return false;
    return true;
}
bool untouched(const ts::SequenceStore& r, size_t i) {
    for (uint32_t p = 0; p < r.length(i); ++p) if (r.baseAt(i, p) < 0) return false;
    return true;
}

void apply() {
    char templ[] = "/tmp/tesseract-deep-norm-XXXXXX";
    const char* tmp = mkdtemp(templ);
    check(tmp != nullptr, "temporary directory");
    const std::string dir = tmp;
    const int k = 21;
    const std::string plasmid = randomSeq(3000), chrom = randomSeq(20000), repeat = randomSeq(400);
    // Pairs, in this order: 300 plasmid-plasmid, 200 chromosome-chromosome, 200 repeat-flank
    // (one mate inside a short high-copy repeat, the other in unique flank); then 50 single
    // plasmid reads.
    std::ofstream a(dir + "/r1.fa"), b(dir + "/r2.fa"), s(dir + "/s.fa");
    auto pair = [&](const std::string& x, const std::string& y, int i) {
        a << ">p" << i << "/1\n" << x << '\n';
        b << ">p" << i << "/2\n" << rc(y) << '\n';
    };
    int id = 0;
    for (int i = 0; i < 300; ++i) { const size_t p = rng() % (plasmid.size() - 300); pair(plasmid.substr(p, 100), plasmid.substr(p + 200, 100), id++); }
    for (int i = 0; i < 200; ++i) { const size_t p = rng() % (chrom.size() - 300); pair(chrom.substr(p, 100), chrom.substr(p + 200, 100), id++); }
    for (int i = 0; i < 200; ++i) { const size_t p = rng() % (repeat.size() - 100); pair(repeat.substr(p, 100), chrom.substr(rng() % (chrom.size() - 100), 100), id++); }
    for (int i = 0; i < 50; ++i) s << ">s" << i << '\n' << plasmid.substr(rng() % (plasmid.size() - 100), 100) << '\n';
    a.close(); b.close(); s.close();
    ts::KmerTable t(200000);
    putKmers(t, chrom, k, 50);        // peak 50 -> bar 1500 (sel 30), 1/6 bar 250
    putKmers(t, repeat, k, 3000);     // deep, but its pairs' other mate is at 50
    putKmers(t, plasmid, k, 5000);
    const double peak = 50;

    auto load = [&]() {
        ts::Library pl; pl.r1 = dir + "/r1.fa"; pl.r2 = dir + "/r2.fa";
        ts::Library sl; sl.r1 = dir + "/s.fa";
        ts::SequenceStore r;
        std::string error;
        check(r.load({pl, sl}, 2, error), "load fixture " + error);
        return r;
    };
    {   // medians: thread-count independent, and as built
        ts::SequenceStore r = load();
        check(r.pairedReads() == 1400 && r.size() == 1450, "fixture shape");
        const std::vector<float> m1 = ts::deepNormMedians(r, t, k, 1);
        const std::vector<float> m3 = ts::deepNormMedians(r, t, k, 3);
        check(m1 == m3, "medians do not depend on the thread count");
        check(m1[0] == 5000.0f && m1[1] == 5000.0f && m1[600] == 50.0f && m1[1000] == 3000.0f &&
              m1[1001] == 50.0f && m1[1400] == 5000.0f, "medians");
    }
    {   // off: nothing runs, nothing is masked
        ts::SequenceStore r = load();
        const ts::DeepNormStats st = ts::deepNormApply(r, t, k, 2, peak, ts::DeepNormConfig());
        check(!st.ran && st.candidates() == 0 && st.maskedReads() == 0, "off: no work");
        bool all = true;
        for (size_t i = 0; i < r.size(); ++i) all = all && untouched(r, i);
        check(all, "off: every read untouched");
        ts::DeepNormConfig on; on.sel = 30;
        check(!ts::deepNormApply(r, t, k, 2, 0.0, on).ran, "no peak: no work");
    }
    for (int rule = 0; rule < 2; ++rule) {
        ts::SequenceStore r = load();
        ts::DeepNormConfig c; c.sel = 30; c.both = rule == 0;
        const ts::DeepNormStats st = ts::deepNormApply(r, t, k, 3, peak, c);
        const std::string tag = rule == 0 ? "dn7: " : "dn8: ";
        check(st.ran, tag + "ran");
        // dn7 selects the 300 plasmid pairs only; dn8 also the 200 repeat-flank pairs (the
        // flank mate, 50, is below 1/6 of the bar, 250 -- so not those either). Both select
        // the 50 single plasmid reads.
        check(st.selectedPairs == 300, tag + "selected pairs " + std::to_string(st.selectedPairs));
        check(st.selectedSingle == 50, tag + "selected singles");
        check(st.thinnedPairs > 250 && st.thinnedPairs < 300, tag + "most deep pairs thinned " + std::to_string(st.thinnedPairs));
        size_t masked = 0;
        for (size_t i = 0; i < r.size(); ++i) {
            if (fullyMasked(r, i)) { ++masked; continue; }
            check(untouched(r, i), tag + "a read is masked whole or not at all");
        }
        check(masked == st.maskedReads(), tag + "masked reads = 2 x thinned pairs + thinned singles");
        for (size_t i = 0; i < 1400; i += 2) check(fullyMasked(r, i) == fullyMasked(r, i + 1), tag + "mates masked together");
        for (size_t i = 600; i < 1400; ++i) check(untouched(r, i), tag + "chromosome and repeat-flank pairs untouched");
        // Each thinned pair is exactly the one whose draw is >= keep.
        for (size_t i = 0; i < 600; i += 2)
            check(fullyMasked(r, i) == !(drawOf(i / 2) < 2.0 * peak / 5000.0), tag + "pair verdict by its index");
    }
    {   // the dn8 rule selects the repeat-flank pairs once the flank mate reaches 1/6 of the bar
        ts::KmerTable t2(200000);
        putKmers(t2, chrom, k, 300);      // flank now at 300 >= 250
        putKmers(t2, repeat, k, 3000);
        putKmers(t2, plasmid, k, 5000);
        ts::SequenceStore r7 = load(), r8 = load();
        ts::DeepNormConfig c7; c7.sel = 30; c7.both = true;
        ts::DeepNormConfig c8; c8.sel = 30;
        const ts::DeepNormStats s7 = ts::deepNormApply(r7, t2, k, 2, peak, c7);
        const ts::DeepNormStats s8 = ts::deepNormApply(r8, t2, k, 2, peak, c8);
        check(s7.selectedPairs == 300, "dn7 ignores repeat-flank pairs");
        check(s8.selectedPairs == 500, "dn8 selects repeat-flank pairs (" + std::to_string(s8.selectedPairs) + ")");
    }
    std::remove((dir + "/r1.fa").c_str());
    std::remove((dir + "/r2.fa").c_str());
    std::remove((dir + "/s.fa").c_str());
    rmdir(dir.c_str());
}

std::string printed(const ts::DeepNormConfig& c, const ts::DeepNormStats& s) {
    char* buf = nullptr;
    size_t len = 0;
    std::FILE* f = open_memstream(&buf, &len);
    ts::deepNormPrintCounts(f, c, s);
    std::fclose(f);
    std::string out(buf, len);
    std::free(buf);
    return out;
}

void counterLine() {
    check(printed(ts::DeepNormConfig(), ts::DeepNormStats()) ==
              "  [deep_norm_counts] enabled=0 ran=0 rule=off selected_pairs=0 thinned_pairs=0 "
              "selected_single=0 thinned_single=0 masked_reads=0\n",
          "idle counter line");
    ts::DeepNormConfig c; c.sel = 30; c.both = true;
    ts::DeepNormStats s; s.ran = true; s.selectedPairs = 7; s.thinnedPairs = 5; s.selectedSingle = 3; s.thinnedSingle = 2;
    check(printed(c, s) ==
              "  [deep_norm_counts] enabled=1 ran=1 rule=both selected_pairs=7 thinned_pairs=5 "
              "selected_single=3 thinned_single=2 masked_reads=12\n",
          "active counter line");
    c.both = false;
    check(printed(c, ts::DeepNormStats()).find("rule=mate:0.16666666666666666 ") != std::string::npos,
          "b3b6d41 rule named with its fraction");
}

}  // namespace

int main() {
    testenv::clearTesseractEnv();
    flags();
    decision();
    apply();
    counterLine();
    std::printf("test_deep_norm: %d checks passed\n", checks);
    return 0;
}
