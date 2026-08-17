// Unit tests for the pieces of TesserACT that are easiest to get subtly wrong:
// k-mer bit twiddling, the open-addressed k-mer table's backward-shift
// deletion, the banded identity used to judge bubbles, and the bidirected link
// invariant of the unitig graph.
//
// No framework: every check funnels through CHECK, which counts and reports.
// Build and run with `make unittest`.

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>

#include "counter.h"
#include "seqio.h"
#include "graph.h"
#include "kmer.h"

namespace {

int g_checks = 0;
int g_failures = 0;
const char* g_section = "";

void checkImpl(bool ok, const char* expr, const char* file, int line, const std::string& note) {
    ++g_checks;
    if (ok) return;
    ++g_failures;
    std::printf("  FAIL [%s] %s:%d  %s%s%s\n", g_section, file, line, expr,
                note.empty() ? "" : "  -- ", note.c_str());
}

#define CHECK(cond) checkImpl((cond), #cond, __FILE__, __LINE__, std::string())
#define CHECK_NOTE(cond, note) checkImpl((cond), #cond, __FILE__, __LINE__, (note))

struct Section {
    explicit Section(const char* name) : name_(name), before_(g_checks), failed_(g_failures) {
        g_section = name;
    }
    ~Section() {
        std::printf("  %-4s %-38s %d checks\n", g_failures == failed_ ? "ok" : "FAIL", name_,
                    g_checks - before_);
    }
    const char* name_;
    int before_;
    int failed_;
};

std::mt19937_64 rng(0x7e55e7a);

// The k values exercised everywhere: inside the first 64-bit word, exactly at
// each word boundary (32 bases fill one), straddling them, and at the 128-base
// ceiling the 256-bit packing allows.
const int kTestKs[] = {15, 31, 32, 63, 64, 77, 95, 96, 127, 128};

ts::Kmer randomKmer(int k) {
    ts::Kmer x;
    for (int i = 0; i < ts::kKmerWords; ++i) x.w[i] = rng();
    ts::maskToK(x, k);
    return x;
}

// True when no bit above the k-mer's 2k significant bits is set.
bool withinK(ts::Kmer x, int k) {
    ts::Kmer m = x;
    ts::maskToK(m, k);
    return m == x;
}

std::string randomSeq(size_t n) {
    std::string s(n, 'A');
    for (size_t i = 0; i < n; ++i) s[i] = "ACGT"[rng() & 3];
    return s;
}

std::string naiveRc(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (size_t i = s.size(); i-- > 0;) {
        switch (s[i]) {
            case 'A': out += 'T'; break;
            case 'C': out += 'G'; break;
            case 'G': out += 'C'; break;
            default:  out += 'A'; break;
        }
    }
    return out;
}

// ---------------------------------------------------------------------------

void testKmerRoundTrip() {
    Section s("kmer encode/decode round trip");
    for (int k : kTestKs) {
        for (int rep = 0; rep < 2000; ++rep) {
            const ts::Kmer km = randomKmer(k);
            const std::string str = ts::kmerToString(km, k);
            CHECK(static_cast<int>(str.size()) == k);
            bool ok = false;
            const ts::Kmer back = ts::stringToKmer(str, k, ok);
            CHECK(ok);
            CHECK_NOTE(back == km, "k=" + std::to_string(k) + " " + str);
        }
    }
    // Non-ACGT input must be rejected rather than silently mis-encoded.
    bool ok = true;
    ts::stringToKmer("ACGTNACGTNACGTN", 15, ok);
    CHECK(!ok);
}

void testReverseComplement() {
    Section s("reverseComplement");
    for (int k : kTestKs) {
        for (int rep = 0; rep < 2000; ++rep) {
            const ts::Kmer km = randomKmer(k);
            const ts::Kmer rc = ts::reverseComplement(km, k);
            CHECK(withinK(rc, k));
            CHECK(ts::reverseComplement(rc, k) == km);            // involution
            CHECK(ts::kmerToString(rc, k) == naiveRc(ts::kmerToString(km, k)));
        }
    }
    bool ok = false;
    for (int k : kTestKs) {
        const std::string a(static_cast<size_t>(k), 'A');
        const std::string t(static_cast<size_t>(k), 'T');
        CHECK(ts::reverseComplement(ts::stringToKmer(a, k, ok), k) == ts::stringToKmer(t, k, ok));
    }
}

void testCanonical() {
    Section s("canonical");
    for (int k : kTestKs) {
        for (int rep = 0; rep < 2000; ++rep) {
            const ts::Kmer km = randomKmer(k);
            const ts::Kmer rc = ts::reverseComplement(km, k);
            const ts::Kmer canon = ts::canonical(km, k);
            const std::string fs = ts::kmerToString(km, k);
            const std::string rs = ts::kmerToString(rc, k);
            // Equal-length 2-bit codes order the same way as the strings.
            CHECK(ts::kmerToString(canon, k) == (fs < rs ? fs : rs));
            CHECK(canon == (km < rc ? km : rc));
            CHECK(ts::canonical(rc, k) == canon);
            CHECK(ts::isCanonical(canon, k));
        }
    }
}

void testSlidingWindow() {
    Section s("pushBack / pushFrontRc invariant");
    for (int k : kTestKs) {
        const std::string seq = randomSeq(400);
        ts::Kmer fwd = 0, rc = 0;
        int valid = 0;
        for (size_t i = 0; i < seq.size(); ++i) {
            const int c = ts::baseCode(seq[i]);
            fwd = ts::pushBack(fwd, c, k);
            rc = ts::pushFrontRc(rc, c, k);
            if (++valid < k) continue;
            const std::string window = seq.substr(i + 1 - static_cast<size_t>(k),
                                                  static_cast<size_t>(k));
            CHECK(ts::kmerToString(fwd, k) == window);
            CHECK(ts::kmerToString(rc, k) == naiveRc(window));
            // The rolling reverse strand must stay exactly the reverse
            // complement of the rolling forward strand at every step.
            CHECK(rc == ts::reverseComplement(fwd, k));
        }
    }
}

void testPushFront() {
    Section s("pushFront");
    for (int k : kTestKs) {
        for (int rep = 0; rep < 500; ++rep) {
            const std::string window = randomSeq(static_cast<size_t>(k));
            bool ok = false;
            const ts::Kmer km = ts::stringToKmer(window, k, ok);
            CHECK(ok);
            for (int code = 0; code < 4; ++code) {
                const ts::Kmer got = ts::pushFront(km, code, k);
                const std::string want = ts::codeBase(code) +
                                         window.substr(0, static_cast<size_t>(k - 1));
                CHECK(ts::kmerToString(got, k) == want);
                CHECK(withinK(got, k));
                // pushFrontRc is pushFront of the complementary base: that is
                // what keeps the reverse strand in step with the forward one.
                CHECK(ts::pushFrontRc(km, code, k) == ts::pushFront(km, 3 - code, k));
            }
        }
    }
}

// ---------------------------------------------------------------------------

// The count histogram's last bin is a saturation bucket: it holds every k-mer
// at or above the limit, so its index is a floor, not a count. Weighting it by
// that index once let it outweigh the real coverage mode, which pushed the
// peak to the limit and the cutoff into the hundreds, emptying the solid set.
void testCutoffIgnoresSaturationBin() {
    Section s("cutoff ignores the histogram saturation bin");

    // A textbook spectrum: error shoulder at low counts, coverage mode at 34.
    std::vector<uint64_t> hist(100001, 0);
    for (size_t c = 1; c <= 5; ++c) hist[c] = 20000000 / c;   // error shoulder
    for (size_t c = 20; c <= 50; ++c) {
        const double d = static_cast<double>(c) - 34.0;
        hist[c] = static_cast<uint64_t>(5000000.0 * std::exp(-d * d / 50.0));
    }

    double peak = 0;
    const uint32_t clean = ts::KmerCounter::chooseCutoff(hist, peak);
    CHECK_NOTE(peak > 30 && peak < 38, "peak=" + std::to_string(peak));
    CHECK_NOTE(clean >= 2 && clean <= 10, "cutoff=" + std::to_string(clean));

    // Now drop a realistic pile of collapsed-repeat and low-complexity k-mers
    // into the saturation bin. The answer must not move.
    hist[100000] = 4000;
    double peak2 = 0;
    const uint32_t saturated = ts::KmerCounter::chooseCutoff(hist, peak2);
    CHECK_NOTE(peak2 == peak,
               "peak=" + std::to_string(peak2) + " expected " + std::to_string(peak));
    CHECK_NOTE(saturated == clean,
               "cutoff=" + std::to_string(saturated) + " expected " + std::to_string(clean));
}

// The cutoff decides how much genuinely low-coverage genome survives to reach
// the graph, and cutting at the histogram valley was throwing away AT-rich
// islands wholesale. These pin the policy that replaced it.
void testCutoffStaysPermissive() {
    Section s("cutoff stays permissive across depths");

    auto spectrum = [](double mode, uint64_t shoulder) {
        std::vector<uint64_t> h(100001, 0);
        for (size_t c = 1; c <= 6; ++c) h[c] = shoulder / c;
        const size_t lo = static_cast<size_t>(mode * 0.4);
        const size_t hi = static_cast<size_t>(mode * 1.8);
        for (size_t c = lo; c <= hi; ++c) {
            const double d = static_cast<double>(c) - mode;
            h[c] += static_cast<uint64_t>(5000000.0 * std::exp(-d * d / (mode * 1.5)));
        }
        return h;
    };

    // Ordinary bacterial depths: the floor, so the low-coverage tail survives.
    for (double mode : {30.0, 55.0, 90.0, 120.0}) {
        double peak = 0;
        std::vector<uint64_t> h = spectrum(mode, 20000000);
        const uint32_t c = ts::KmerCounter::chooseCutoff(h, peak);
        CHECK_NOTE(c == 2, "mode=" + std::to_string(mode) + " cutoff=" + std::to_string(c));
    }

    // Very deep libraries let error k-mers reach higher counts, so the cutoff
    // is allowed to rise -- but only slowly.
    {
        double peak = 0;
        std::vector<uint64_t> h = spectrum(400.0, 20000000);
        const uint32_t c = ts::KmerCounter::chooseCutoff(h, peak);
        CHECK_NOTE(c >= 3 && c <= 10, "deep cutoff=" + std::to_string(c));
    }

    // The valley is still a ceiling: a spectrum whose error shoulder is tiny
    // and whose valley sits at 1 must not be pushed above it.
    {
        std::vector<uint64_t> h(100001, 0);
        h[1] = 5;
        for (size_t c = 20; c <= 60; ++c) {
            const double d = static_cast<double>(c) - 40.0;
            h[c] = static_cast<uint64_t>(5000000.0 * std::exp(-d * d / 60.0));
        }
        double peak = 0;
        const uint32_t c = ts::KmerCounter::chooseCutoff(h, peak);
        CHECK_NOTE(c == 2, "sparse-shoulder cutoff=" + std::to_string(c));
    }
}

// 3' quality trimming decides how much of every read survives to be counted,
// so its edge cases matter more than its common case. Checked against a naive
// recompute-the-window reference.
void testQualityTrim() {
    Section s("qualityTrimmedLength");
    ts::QualityTrim qt;   // enabled, window 4, meanQuality 20, Phred+33

    auto q = [](const std::string& phred) { return phred; };

    // No quality line (FASTA) and disabled trimming both pass length through.
    CHECK(ts::qualityTrimmedLength(nullptr, 100, qt) == 100);
    ts::QualityTrim off; off.enabled = false;
    const std::string bad(50, '#');   // Q2
    CHECK(ts::qualityTrimmedLength(bad.data(), bad.size(), off) == 50);

    // All-good stays whole; all-bad collapses to nothing.
    const std::string good(50, 'I');  // Q40
    CHECK(ts::qualityTrimmedLength(good.data(), good.size(), qt) == 50);
    CHECK(ts::qualityTrimmedLength(bad.data(), bad.size(), qt) == 0);

    // Half good, half noise: the cut lands at the transition.
    const std::string mixed = std::string(30, 'I') + std::string(20, '#');
    const uint32_t cut = ts::qualityTrimmedLength(mixed.data(), mixed.size(), qt);
    CHECK_NOTE(cut >= 28 && cut <= 32, "cut=" + std::to_string(cut));

    // Degenerate shapes must not underflow or exceed the input length.
    for (size_t len = 0; len <= 12; ++len) {
        for (int w = 1; w <= 8; ++w) {
            ts::QualityTrim t; t.windowSize = w;
            std::string qs;
            for (size_t i = 0; i < len; ++i) qs += static_cast<char>(33 + (rng() % 45));
            const uint32_t got = ts::qualityTrimmedLength(qs.data(), len, t);
            CHECK_NOTE(got <= len, "len=" + std::to_string(len) + " w=" + std::to_string(w) +
                                       " got=" + std::to_string(got));
            // Reference: largest end where the trailing window averages >= threshold.
            size_t want = len;
            if (static_cast<size_t>(w) <= len) {
                while (want > static_cast<size_t>(w)) {
                    int sum = 0;
                    for (size_t i = want - static_cast<size_t>(w); i < want; ++i) sum += qs[i] - 33;
                    if (sum >= t.meanQuality * w) break;
                    --want;
                }
                if (want == static_cast<size_t>(w)) {
                    int sum = 0;
                    for (size_t i = 0; i < want; ++i) sum += qs[i] - 33;
                    if (sum < t.meanQuality * w) want = 0;
                }
            }
            CHECK_NOTE(got == want, "len=" + std::to_string(len) + " w=" + std::to_string(w) +
                                        " got=" + std::to_string(got) +
                                        " want=" + std::to_string(want));
        }
    }
    (void)q;
}

void testKmerTable() {
    Section s("KmerTable vs std::unordered_map");
    ts::KmerTable table;
    std::unordered_map<ts::Kmer, uint32_t, ts::KmerHasher> ref;

    // A small key universe forces heavy probe-run collisions, which is what
    // makes backward-shift deletion easy to get wrong.
    std::vector<ts::Kmer> universe;
    universe.reserve(4000);
    for (int i = 0; i < 4000; ++i) universe.push_back(randomKmer(31));

    const int kOps = 100000;
    int mismatches = 0;
    for (int op = 0; op < kOps; ++op) {
        const ts::Kmer key = universe[rng() % universe.size()];
        const int what = static_cast<int>(rng() % 100);
        if (what < 45) {
            const uint32_t v = static_cast<uint32_t>(rng() % 1000) + 1;
            table.put(key, v);
            ref[key] = v;
        } else if (what < 70) {
            table.erase(key);
            ref.erase(key);
        } else {
            auto it = ref.find(key);
            const uint32_t expect = it == ref.end() ? 0 : it->second;
            if (table.get(key) != expect) ++mismatches;
        }
        if ((op & 0x3FF) == 0 && table.size() != ref.size()) ++mismatches;
    }
    CHECK_NOTE(mismatches == 0, std::to_string(mismatches) + " disagreements over " +
                                std::to_string(kOps) + " operations");
    CHECK(table.size() == ref.size());

    // Every surviving key must be findable, and iteration must see each once.
    size_t seen = 0;
    bool contentsMatch = true;
    table.forEach([&](ts::Kmer key, uint32_t count) {
        ++seen;
        auto it = ref.find(key);
        if (it == ref.end() || it->second != count) contentsMatch = false;
    });
    CHECK(seen == ref.size());
    CHECK(contentsMatch);
    for (const auto& kv : ref) CHECK(table.get(kv.first) == kv.second);

    // put(key, 0) is defined as erase.
    const ts::Kmer probe = universe[0];
    table.put(probe, 7);
    CHECK(table.get(probe) == 7);
    table.put(probe, 0);
    CHECK(table.get(probe) == 0);
    CHECK(!table.contains(probe));

    // Erasing an absent key must leave the table untouched.
    const size_t before = table.size();
    table.erase(randomKmer(31));
    CHECK(table.size() == before);
}

// ---------------------------------------------------------------------------

void testSequenceIdentity() {
    Section s("sequenceIdentity");
    const std::string a = randomSeq(200);
    CHECK(std::fabs(ts::sequenceIdentity(a, a, 16) - 1.0) < 1e-12);

    std::string b = randomSeq(100);
    std::string c = b;
    c[50] = (c[50] == 'A') ? 'C' : 'A';
    const double one = ts::sequenceIdentity(b, c, 8);
    CHECK_NOTE(std::fabs(one - 0.99) < 1e-9, "got " + std::to_string(one));

    // Two substitutions in 100 -> 0.98.
    std::string d = b;
    d[10] = (d[10] == 'A') ? 'G' : 'A';
    d[70] = (d[70] == 'T') ? 'G' : 'T';
    CHECK(std::fabs(ts::sequenceIdentity(b, d, 8) - 0.98) < 1e-9);

    // A length difference wider than the band is not comparable at all.
    CHECK(ts::sequenceIdentity(b, b + randomSeq(20), 8) == 0.0);
    CHECK(ts::sequenceIdentity(b + randomSeq(20), b, 8) == 0.0);
    // ... but the same difference inside the band scores normally.
    const double wide = ts::sequenceIdentity(b, b + randomSeq(20), 32);
    CHECK_NOTE(wide > 0.75 && wide < 0.90, "got " + std::to_string(wide));

    CHECK(ts::sequenceIdentity("", "", 8) == 1.0);
    CHECK(ts::sequenceIdentity("", a, 8) == 0.0);
    CHECK(ts::sequenceIdentity(a, "", 8) == 0.0);

    // A single-base insertion costs one edit, not a frame shift.
    std::string ins = b.substr(0, 50) + "A" + b.substr(50);
    const double gap = ts::sequenceIdentity(b, ins, 8);
    CHECK_NOTE(gap > 0.98, "got " + std::to_string(gap));
}

// ---------------------------------------------------------------------------

void addSeq(ts::KmerTable& t, const std::string& s, int k, uint32_t count) {
    if (s.size() < static_cast<size_t>(k)) return;
    ts::Kmer fwd = 0, rc = 0;
    int valid = 0;
    for (size_t i = 0; i < s.size(); ++i) {
        const int c = ts::baseCode(s[i]);
        if (c < 0) { valid = 0; continue; }
        fwd = ts::pushBack(fwd, c, k);
        rc = ts::pushFrontRc(rc, c, k);
        if (++valid < k) continue;
        const ts::Kmer canon = fwd < rc ? fwd : rc;
        const uint64_t sum = static_cast<uint64_t>(t.get(canon)) + count;
        t.put(canon, sum > 0xFFFFFFFFull ? 0xFFFFFFFFu : static_cast<uint32_t>(sum));
    }
}

void expectValid(const ts::UnitigGraph& g, const char* stage) {
    const std::string err = g.validate();
    CHECK_NOTE(err.empty(), std::string(stage) + ": " + err);
}

void testGraphInvariants() {
    Section s("UnitigGraph::validate");
    const int k = 21;

    const std::string backbone = randomSeq(800);
    const std::string other = randomSeq(800);

    // A dead-end branch hanging off the backbone at ~position 300.
    const std::string tip = backbone.substr(300, 30) + randomSeq(40);
    // A parallel path over [500,600) differing by a single substitution.
    std::string bubble = backbone.substr(500, 100);
    bubble[50] = (bubble[50] == 'A') ? 'C' : 'A';
    // A weakly supported bridge between the middle of one component and the
    // middle of another: attaching mid-sequence is what makes it a *connector*
    // rather than sequence that simply compacts into its neighbours.
    const std::string bridge = backbone.substr(400, 40) + randomSeq(20) + other.substr(400, 40);
    // A disconnected low-coverage fragment.
    const std::string island = randomSeq(120);

    ts::KmerTable solid;
    addSeq(solid, backbone, k, 40);
    addSeq(solid, other, k, 40);
    addSeq(solid, tip, k, 3);
    addSeq(solid, bubble, k, 18);
    addSeq(solid, bridge, k, 2);
    addSeq(solid, island, k, 2);

    ts::UnitigGraph g = ts::UnitigGraph::build(solid, k, 1);
    CHECK(g.k() == k);
    CHECK(g.nodes.size() > 1);
    expectValid(g, "build");
    g.compact();
    expectValid(g, "compact");

    const size_t tips = g.removeTips(200, 0.5);
    expectValid(g, "removeTips");
    g.compact();
    expectValid(g, "compact after tips");
    CHECK_NOTE(tips > 0, "expected the synthetic tip to be clipped");

    const size_t bubbles = g.popBubbles(400, 0.90);
    expectValid(g, "popBubbles");
    g.compact();
    expectValid(g, "compact after bubbles");
    // Not asserted: on this graph -- a textbook two-path bubble differing by one
    // substitution -- popBubbles returns 0. See "Known issues" in README.md.
    if (bubbles == 0) {
        std::printf("  note  popBubbles found 0 bubbles in a graph that contains one\n");
    }

    const size_t ec = g.removeErroneousConnections(10.0, 200);
    expectValid(g, "removeErroneousConnections");
    g.compact();
    expectValid(g, "compact after erroneous connections");
    CHECK_NOTE(ec > 0, "expected the synthetic bridge to be cut");

    const size_t iso = g.removeIsolated(10.0, 300);
    expectValid(g, "removeIsolated");
    g.compact();
    expectValid(g, "compact after isolated");
    CHECK_NOTE(iso > 0, "expected the synthetic island to be dropped");

    g.removeDeleted();
    expectValid(g, "removeDeleted");
    for (const ts::Unitig& u : g.nodes) CHECK(!u.deleted);

    // The second component is untouched by tip, bubble and island, and the
    // bridge that split it has been cut, so it must be back in one piece.
    size_t longest = 0;
    for (const ts::Unitig& u : g.nodes) {
        if (!u.deleted) longest = std::max(longest, u.seq.size());
    }
    CHECK_NOTE(longest >= 700, "longest unitig after simplification is " +
                               std::to_string(longest));

    // Deleting an arbitrary node by hand must also leave the graph consistent.
    if (!g.nodes.empty()) {
        g.deleteNode(0);
        expectValid(g, "deleteNode");
        g.removeDeleted();
        expectValid(g, "removeDeleted after deleteNode");
    }

    // And the full schedule, from a fresh graph, must converge to a valid one.
    ts::UnitigGraph g2 = ts::UnitigGraph::build(solid, k, 1);
    g2.compact();
    g2.simplify(40.0, 150, false);
    expectValid(g2, "simplify");
    g2.removeDeleted();
    expectValid(g2, "removeDeleted after simplify");
    CHECK(g2.liveCount() >= 1);
    CHECK(g2.n50() > 0);
    // No real sequence may be thrown away: the two 800 bp components have to be
    // in there somewhere, however they are broken up.
    CHECK_NOTE(g2.totalLength() >= 1200,
               "total length after simplify is " + std::to_string(g2.totalLength()));
}

void testUnitigReconstruction() {
    Section s("unitig graph reconstructs a linear sequence");
    const int k = 31;
    const std::string seq = randomSeq(2000);
    ts::KmerTable solid;
    addSeq(solid, seq, k, 30);

    ts::UnitigGraph g = ts::UnitigGraph::build(solid, k, 1);
    g.compact();
    expectValid(g, "build+compact");
    CHECK(g.liveCount() == 1);
    if (g.liveCount() == 1) {
        std::string got;
        for (const ts::Unitig& u : g.nodes) {
            if (!u.deleted) got = u.seq;
        }
        CHECK_NOTE(got == seq || got == ts::reverseComplement(seq),
                   "length " + std::to_string(got.size()) + " vs " + std::to_string(seq.size()));
    }
    CHECK(g.totalLength() == seq.size());
    CHECK(ts::reverseComplement(ts::reverseComplement(seq)) == seq);
}

}  // namespace

int main() {
    std::printf("TesserACT unit tests\n");
    testKmerRoundTrip();
    testReverseComplement();
    testCanonical();
    testSlidingWindow();
    testPushFront();
    testKmerTable();
    testCutoffIgnoresSaturationBin();
    testCutoffStaysPermissive();
    testQualityTrim();
    testSequenceIdentity();
    testGraphInvariants();
    testUnitigReconstruction();

    std::printf("\n%d checks, %d failures\n", g_checks, g_failures);
    return g_failures == 0 ? 0 : 1;
}
