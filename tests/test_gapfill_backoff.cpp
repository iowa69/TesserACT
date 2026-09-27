// Component tests for PKG-GAP (combo2 DESIGN.md section 4), both default off:
//   TESSERACT_GAP_KEEP_FLANK=1   the first unitig after a scaffold gap is emitted whole;
//   TESSERACT_GAPFILL_BACKOFF=B  the gap filler steps a below-floor anchor back by <= B
//                                bases and re-spells the stepped-over tip from the reads.
// T-G0 keep-flank on a scaffolded pair of unitigs; T-G1..T-G4 the back-off on
// synthetic genomes with deterministic read tilings; parsing of both flags.
// T-G5/T-G6/T-G7 are run on real isolates (see the tree's README.md).
//
// build_v3 (integration): build_v3 keeps ONE implementation of each fix these flags carried.
//   TESSERACT_GAP_KEEP_FLANK=1   is an alias of TESSERACT_FIX_GAP_FLANK (G-resolve T03, the
//                                record-and-restore design; test_v3_resolve_gapflank holds its
//                                own fixtures). T-G0 now checks the alias wiring and T03's
//                                open-gap result on the package fixtures; the package's
//                                render-whole expectations (N count unchanged, "=all", the
//                                one-mismatch overlap) describe the retired implementation.
//   TESSERACT_GAPFILL_BACKOFF=B  is an alias of T11's back-off (G-emit) with B bases; the
//                                T-G1..T-G4 fixtures run unchanged against it.
// Malformed values are rejected by the envflags table (exit 2 at startup, T16) instead of
// reading as "off"; the parsing checks test the table directly.
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <functional>
#include <random>
#include <string>
#include <utility>
#include <vector>
#include <unistd.h>

#include "envflags.h"
#include "gapfill.h"
#include "graph.h"
#include "resolve.h"
#include "seqio.h"
#include "test_env.h"

namespace {
int checks = 0;
int failures = 0;
void check(bool ok, const std::string& what) {
    ++checks;
    if (!ok) { ++failures; std::fprintf(stderr, "FAIL: %s\n", what.c_str()); }
}

std::string dna(size_t n, unsigned seed) {
    std::mt19937 rng(seed);
    std::string s(n, 'A');
    for (char& b : s) b = "ACGT"[rng() % 4];
    return s;
}
std::string rc(const std::string& s) {
    std::string out;
    for (auto i = s.rbegin(); i != s.rend(); ++i) out += "TGCA"[ts::baseCode(*i)];
    return out;
}
// A substitution that is never the base it replaces.
std::string mutate(std::string s, size_t pos) {
    s[pos] = s[pos] == 'A' ? 'C' : 'A';
    return s;
}

using Pair = std::pair<std::string, std::string>;

// FR pairs tiled along `g`: a `frag`-base fragment starts every `step` bases in
// [from, to - frag], read 1 its first `len` bases, read 2 the reverse
// complement of its last `len`. `keep` may veto a fragment by its start.
void tile(std::vector<Pair>& out, const std::string& g, size_t from, size_t to, size_t step,
          const std::function<bool(size_t)>& keep = nullptr, size_t len = 100, size_t frag = 300) {
    for (size_t p = from; p + frag <= to && p + frag <= g.size(); p += step) {
        if (keep && !keep(p)) continue;
        out.push_back({g.substr(p, len), rc(g.substr(p + frag - len, len))});
    }
}

ts::SequenceStore store(const std::vector<Pair>& pairs) {
    const auto dir = std::filesystem::temp_directory_path() /
                     ("tesseract-gapfill-backoff-" + std::to_string(getpid()));
    std::filesystem::create_directories(dir);
    const auto first = dir / "r1.fa", second = dir / "r2.fa";
    {
        std::ofstream f(first), s(second);
        for (size_t i = 0; i < pairs.size(); ++i) {
            f << ">" << i << "\n" << pairs[i].first << "\n";
            s << ">" << i << "\n" << pairs[i].second << "\n";
        }
    }
    ts::SequenceStore st;
    ts::Library lib;
    lib.r1 = first.string();
    lib.r2 = second.string();
    std::string error;
    if (!st.load({lib}, 1, error)) {
        std::fprintf(stderr, "load failed: %s\n", error.c_str());
        std::exit(2);
    }
    std::filesystem::remove_all(dir);
    return st;
}

struct Closed {
    std::vector<std::string> contigs;
    ts::GapFillStats st;
};
// backoff == nullptr: TESSERACT_GAPFILL_BACKOFF unset (release); otherwise set to it.
Closed closeWith(const std::vector<std::string>& contigs, const ts::SequenceStore& reads,
                 const char* backoff, int flank = 300) {
    if (backoff) setenv("TESSERACT_GAPFILL_BACKOFF", backoff, 1);
    else unsetenv("TESSERACT_GAPFILL_BACKOFF");
    Closed r{contigs, {}};
    r.st = ts::closeGaps(r.contigs, reads, 2, 31, flank);
    unsetenv("TESSERACT_GAPFILL_BACKOFF");
    return r;
}

// Every N-free stretch of `s` occurs in `g`: no splice wrote a base the genome
// does not have there, and none duplicated or dropped sequence inside a stretch.
bool piecesInGenome(const std::string& s, const std::string& g) {
    size_t pos = 0;
    while (pos < s.size()) {
        while (pos < s.size() && s[pos] == 'N') ++pos;
        size_t e = pos;
        while (e < s.size() && s[e] != 'N') ++e;
        if (e > pos && g.find(s.substr(pos, e - pos)) == std::string::npos) return false;
        pos = e;
    }
    return true;
}

// ---- T-G0: keep-flank emission (build_v3: alias of T03) ---------------------------
// Two unitigs with no graph path between them, which the scaffolder joins across
// Ns: B starts `gapOrOverlap` bases after A ends (negative: B's first bases
// repeat A's last ones -- two dead ends overlapping by less than k-1, so the
// graph has no edge). `mismatch` >= 0 plants a substitution in B's copy of the
// overlap. Release drops B's first k-1 bases. build_v3: GAP_KEEP_FLANK=1 must be exactly
// TESSERACT_FIX_GAP_FLANK=1 (resolve() keeps the release layout; restoreGapFlanks() after gap
// filling gives B's k-1 bases back), and the fix's own variable wins over the alias.
void keepFlankCase(int k, long gapOrOverlap, long mismatch, const std::string& tag) {
    const std::string G = dna(6400, 11);
    const size_t ov = static_cast<size_t>(k - 1);
    ts::UnitigGraph g;
    g.setK(k);
    g.nodes.resize(2);
    g.nodes[0].seq = G.substr(0, 3000);
    g.nodes[1].seq = G.substr(static_cast<size_t>(3000 + gapOrOverlap), 3000);
    if (mismatch >= 0) g.nodes[1].seq = mutate(g.nodes[1].seq, static_cast<size_t>(mismatch));
    g.nodes[0].coverage = g.nodes[1].coverage = 40;
    std::vector<Pair> pairs;
    tile(pairs, G, 0, 6100, 5);
    const ts::SequenceStore reads = store(pairs);

    struct Out {
        std::vector<std::string> seqs, restored;
        std::vector<ts::ResolvedPath> paths;
        ts::ResolveStats st;
    };
    auto run = [&](const char* alias, const char* own) {
        if (alias) setenv("TESSERACT_GAP_KEEP_FLANK", alias, 1);
        else unsetenv("TESSERACT_GAP_KEEP_FLANK");
        if (own) setenv("TESSERACT_FIX_GAP_FLANK", own, 1);
        else unsetenv("TESSERACT_FIX_GAP_FLANK");
        ts::PairedResolver r(g, reads, 1, 2, 2.0, 0.10);
        r.setScaffolding(true);
        r.buildSupport();
        Out o;
        std::vector<double> covs;
        r.resolve(o.seqs, covs);
        o.paths = r.paths();
        o.st = r.stats();
        o.restored = o.seqs;
        ts::restoreGapFlanks(o.restored, r.gapFlankRecords());   // no gap filling here: every gap open
        unsetenv("TESSERACT_GAP_KEEP_FLANK");
        unsetenv("TESSERACT_FIX_GAP_FLANK");
        return o;
    };
    const Out off = run(nullptr, nullptr), alias = run("1", nullptr), own = run(nullptr, "1");
    const Out aliasOff = run("0", nullptr), ownWins = run("1", "0");
    const std::string A = g.nodes[0].seq, B = g.nodes[1].seq;
    auto scaffold = [&](const std::vector<std::string>& seqs) -> std::string {
        for (const std::string& s : seqs) {
            if (s.find('N') == std::string::npos) continue;
            return s.compare(0, 60, A, 0, 60) == 0 ? s : rc(s);   // either orientation; A first
        }
        return std::string();
    };
    const std::string sOff = scaffold(off.seqs);
    check(off.st.scaffoldJoins == 1 && !sOff.empty(), tag + " fixture: the two unitigs are scaffolded");
    const size_t n = off.st.gapBases;
    check(sOff == A + std::string(n, 'N') + B.substr(ov), tag + " flag off: release drops k-1 bases after the gap");
    check(off.restored == off.seqs, tag + " flag off: nothing to restore");
    check(alias.seqs == off.seqs, tag + " alias: resolve() keeps the release layout (T03)");
    check(alias.restored == own.restored && alias.seqs == own.seqs,
          tag + " GAP_KEEP_FLANK=1 is byte-identical to TESSERACT_FIX_GAP_FLANK=1");
    check(aliasOff.restored == off.restored, tag + " GAP_KEEP_FLANK=0 is release");
    check(ownWins.restored == off.restored, tag + " TESSERACT_FIX_GAP_FLANK=0 wins over the alias");
    const std::string sOn = scaffold(alias.restored);
    check(sOn != sOff, tag + " alias: the open gap is restored");
    // Whatever the overlap rule decided, B's sequence after the Ns is B less at most k-1 bases
    // (a verified overlap), and nothing before the Ns changed.
    const size_t lastN = sOn.find_last_of('N');
    const size_t firstN = sOn.find('N');
    check(firstN == A.size() && sOn.compare(0, A.size(), A) == 0, tag + " alias: A intact before the Ns");
    const std::string after = lastN == std::string::npos ? std::string() : sOn.substr(lastN + 1);
    check(after.size() >= B.size() - ov && after.size() <= B.size() &&
              B.compare(B.size() - after.size(), after.size(), after) == 0,
          tag + " alias: B follows the Ns, whole or less a verified overlap");
    if (gapOrOverlap > static_cast<long>(ov))
        check(sOn == A + std::string(n - ov, 'N') + B,
              tag + " alias: a true gap restores B whole and writes the estimated gap (N = gap - (k-1))");
    for (const Out* o : {&alias, &own}) {
        check(o->st.scaffoldJoins == off.st.scaffoldJoins && o->seqs.size() == off.seqs.size(),
              tag + " join count and contig count unchanged");
    }
}
void testKeepFlank() {
    keepFlankCase(31, 100, -1, "T-G0 gap 100");
    keepFlankCase(31, -20, -1, "T-G0 exact overlap 20");
    keepFlankCase(61, -45, 10, "T-G0 overlap 45, one mismatch, k=61");
    keepFlankCase(61, -25, 10, "T-G0 overlap 25 with a mismatch, k=61");
    // The package accepted "all" (keep every k-1 base, the literal rule); build_v3 retires it,
    // and the table refuses it at startup like any other malformed value.
    const ts::env::Spec* spec = ts::env::find("TESSERACT_GAP_KEEP_FLANK");
    ts::env::Parsed p;
    std::string why;
    check(spec && !ts::env::parse(*spec, "all", p, why), "T-G0 GAP_KEEP_FLANK=all is refused (retired)");
    check(spec && !ts::env::parse(*spec, "yes", p, why), "T-G0 GAP_KEEP_FLANK=yes is refused");
}

// ---- back-off fixtures -----------------------------------------------------------
const std::string& genome() {
    static const std::string G = dna(5000, 7);
    return G;
}

// T-G1: the last 30 bp of the left flank are spelled by two reads only; they carry
// a substitution 15 bp before the gap that the 40x tiling of the locus does not.
void testTaperedErrorTip() {
    const std::string& G = genome();
    const std::string Gerr = mutate(G, 1985);
    std::vector<Pair> pairs;
    tile(pairs, G, 0, 5000, 5);
    for (int i = 0; i < 2; ++i) pairs.push_back({Gerr.substr(1900, 100), rc(G.substr(2100, 100))});
    const ts::SequenceStore reads = store(pairs);
    const std::vector<std::string> contig{Gerr.substr(0, 2000) + std::string(50, 'N') + G.substr(2050, 1950)};

    const Closed off = closeWith(contig, reads, nullptr);
    check(off.st.gapsClosed == 0 && off.contigs == contig, "T-G1 B=0: not closed");
    check(off.st.seedBelowFloor == 1, "T-G1 B=0: seedgone");
    check(off.st.seedBackedOff == 0 && off.st.closedAfterBackoff == 0, "T-G1 B=0: back-off counters zero");

    const Closed on = closeWith(contig, reads, "60");
    check(on.st.gapsClosed == 1, "T-G1 B=60: closed");
    check(on.contigs.size() == 1 && on.contigs[0] == G.substr(0, 4000), "T-G1 B=60: contig equals the truth");
    check(on.st.seedBackedOff == 1 && on.st.targetBackedOff == 0, "T-G1 B=60: only the seed moved");
    check(on.st.closedAfterBackoff == 1 && on.st.backoffCancelled == 0, "T-G1 B=60: closed after back-off");
    check(on.st.trimmedBp == 15, "T-G1 B=60: the 15 tip bases after the moved seed were re-spelled");
    check(on.st.nBasesRemoved == 50 && on.st.basesInserted == 65, "T-G1 B=60: 50 Ns + 15 tip -> 65 bases");

    const Closed shortB = closeWith(contig, reads, "10");
    check(shortB.st.gapsClosed == 0 && shortB.contigs == contig && shortB.st.seedBackedOff == 0,
          "T-G1 B=10: no supported seed within reach, anchor kept, not closed");
}

// T-G1b/c: a tip whose bases are right but whose coverage genuinely thins. The
// moved seed's walk must run through the old one, so the back-off changes nothing.
void testTrueTaper() {
    const std::string& G = genome();
    // Fragments with a read overlapping [1960, 2060) are dropped from the tiling.
    auto clear = [](size_t p) { return (p >= 2060 || p + 100 <= 1960) && (p + 200 >= 2060 || p + 300 <= 1960); };
    const std::vector<std::string> contig{G.substr(0, 2000) + std::string(50, 'N') + G.substr(2050, 1950)};
    {
        // Two 150 bp reads cross the whole dropout: the release closes it at floor 2.
        std::vector<Pair> pairs;
        tile(pairs, G, 0, 5000, 5, clear);
        for (int i = 0; i < 2; ++i) pairs.push_back({G.substr(1950, 150), rc(G.substr(2250, 150))});
        const ts::SequenceStore reads = store(pairs);
        const Closed off = closeWith(contig, reads, nullptr), on = closeWith(contig, reads, "60");
        check(off.st.gapsClosed == 1 && off.contigs[0] == G.substr(0, 4000), "T-G1b release closes the true taper");
        check(on.contigs == off.contigs && on.st.gapsClosed == off.st.gapsClosed &&
                  on.st.closedAfterBackoff == 0 && on.st.seedBackedOff == 0,
              "T-G1b B=60: release closure kept byte-identical, back-off not tried");
    }
    {
        // Nothing crosses the dropout: open with and without the back-off.
        std::vector<Pair> pairs;
        tile(pairs, G, 0, 5000, 5, clear);
        const ts::SequenceStore reads = store(pairs);
        const Closed off = closeWith(contig, reads, nullptr), on = closeWith(contig, reads, "60");
        check(off.st.gapsClosed == 0 && on.st.gapsClosed == 0 && on.contigs == contig,
              "T-G1c unspanned dropout stays open with B=60");
        check(on.st.seedBackedOff + on.st.targetBackedOff > 0 && on.st.closedAfterBackoff == 0,
              "T-G1c back-off tried and closed nothing");
    }
}

// T-G2: two haplotypes differ inside the gap and are equally deep: two paths of
// equal weakest-link support, so the gap stays open with or without back-off --
// including when the back-off is what reaches the gap (the T-G1 error tip).
void testTwoEqualPaths() {
    const std::string& G = genome();
    const std::string H = mutate(G, 2025);
    for (int errorTip = 0; errorTip < 2; ++errorTip) {
        std::vector<Pair> pairs;
        tile(pairs, G, 0, 5000, 10);
        tile(pairs, H, 0, 5000, 10);
        std::string left = G.substr(0, 2000);
        if (errorTip) {
            const std::string Gerr = mutate(G, 1985);
            for (int i = 0; i < 2; ++i) pairs.push_back({Gerr.substr(1900, 100), rc(G.substr(2100, 100))});
            left = Gerr.substr(0, 2000);
        }
        const ts::SequenceStore reads = store(pairs);
        const std::vector<std::string> contig{left + std::string(50, 'N') + G.substr(2050, 1950)};
        const Closed off = closeWith(contig, reads, nullptr), on = closeWith(contig, reads, "60");
        const std::string tag = errorTip ? "T-G2 (error tip)" : "T-G2";
        check(off.st.gapsClosed == 0 && off.contigs == contig, tag + " B=0: open");
        check(on.st.gapsClosed == 0 && on.contigs == contig, tag + " B=60: open");
        if (errorTip) {
            // Counted as the release attempt left it (unspanned from the error seed).
            check(on.st.seedBackedOff == 1 && on.st.closedAfterBackoff == 0 && on.st.gapsNoPath == 1,
                  tag + " B=60: back-off tried, two paths, nothing closed");
        } else {
            check(off.st.gapsAmbiguous == 1 && on.st.seedBackedOff == 0 && on.st.targetBackedOff == 0,
                  tag + " anchors at the floor: ambiguous, back-off not needed");
        }
    }
}

// T-G3: the flanks overlap on the genome by 50 bp. No walk goes backwards, so
// the gap stays open and nothing is duplicated or removed -- also when the
// left tip is an error tip and the seed is moved.
void testOverlappingFlanks() {
    const std::string& G = genome();
    for (int errorTip = 0; errorTip < 2; ++errorTip) {
        std::vector<Pair> pairs;
        tile(pairs, G, 0, 5000, 5);
        std::string left = G.substr(0, 2000);
        if (errorTip) {
            const std::string Gerr = mutate(G, 1985);
            for (int i = 0; i < 2; ++i) pairs.push_back({Gerr.substr(1900, 100), rc(G.substr(2100, 100))});
            left = Gerr.substr(0, 2000);
        }
        const ts::SequenceStore reads = store(pairs);
        const std::vector<std::string> contig{left + std::string(50, 'N') + G.substr(1950, 2050)};
        const Closed off = closeWith(contig, reads, nullptr), on = closeWith(contig, reads, "60");
        const std::string tag = errorTip ? "T-G3 (error tip)" : "T-G3";
        check(off.contigs == contig && on.contigs == contig && on.st.gapsClosed == 0,
              tag + ": output unchanged by the back-off, gap open, nothing duplicated");
        check(on.st.seedBackedOff == static_cast<size_t>(errorTip), tag + ": seed moved only for the error tip");
    }
}

// T-G4: two gaps 70 bp apart whose facing tips both carry a two-read error, so
// the first gap's target and the second gap's seed both move into the 70 bp
// between them. Both close, the splices do not meet, and the contig is the truth.
// flank 40 keeps the two gaps' recruitment windows apart; at the production
// 300 they share every k-mer between the gaps, so the invariant is checked there.
void testNearbyGaps() {
    const std::string& G = genome();
    const std::string mid = mutate(mutate(G, 2055), 2115).substr(2050, 70);
    std::vector<Pair> pairs;
    tile(pairs, G, 0, 5000, 2);
    for (int i = 0; i < 2; ++i) {
        pairs.push_back({mutate(G, 2055).substr(2050, 100), rc(G.substr(2250, 100))});
        pairs.push_back({mutate(G, 2115).substr(2020, 100), rc(G.substr(2220, 100))});
    }
    const ts::SequenceStore reads = store(pairs);
    const std::vector<std::string> contig{G.substr(0, 2000) + std::string(50, 'N') + mid +
                                          std::string(50, 'N') + G.substr(2170, 1830)};
    const Closed off = closeWith(contig, reads, nullptr, 40), on = closeWith(contig, reads, "60", 40);
    check(off.st.gapsClosed == 0 && off.contigs == contig, "T-G4 B=0: both gaps open");
    check(on.st.gapsClosed == 2 && on.contigs[0] == G.substr(0, 4000), "T-G4 B=60: both closed, contig is the truth");
    check(on.st.targetBackedOff == 1 && on.st.seedBackedOff == 1 && on.st.backoffCancelled == 0 &&
              on.st.trimmedBp == 11,
          "T-G4 B=60: target moved 6, seed moved 5, splices apart, none cancelled");
    for (int flank : {40, 300}) {
        const Closed r = closeWith(contig, reads, "200", flank);
        check(r.contigs.size() == 1 && piecesInGenome(r.contigs[0], G),
              "T-G4 B=200 flank " + std::to_string(flank) + ": every N-free stretch is genome sequence");
    }
}

// T-G4, the rule itself: a back-off closure whose replaced stretch reaches the
// next closure's seed k-mer is dropped (the later one when both backed off) and
// counted; release closures are never dropped.
void testClashRule() {
    const int k = 31;
    using C = ts::GapClosure;
    {
        std::vector<C> c{{0, 2000, 50, false, true}, {0, 2080, 60, true, true}};   // 2080 < 2050 + 31
        check(ts::cancelBackoffClashes(c, k) == 1 && c[0].closed && !c[1].closed,
              "T-G4 rule: back-off closure reaching a release closure's footprint is dropped");
    }
    {
        std::vector<C> c{{0, 2000, 100, true, true}, {0, 2120, 50, false, true}};  // 2120 < 2100 + 31
        check(ts::cancelBackoffClashes(c, k) == 1 && !c[0].closed && c[1].closed,
              "T-G4 rule: the release closure survives, the earlier back-off one goes");
    }
    {
        std::vector<C> c{{0, 2000, 100, true, true}, {0, 2110, 50, true, true}, {0, 2131, 40, true, true}};
        check(ts::cancelBackoffClashes(c, k) == 1 && c[0].closed && !c[1].closed && c[2].closed,
              "T-G4 rule: of two back-off closures the later is dropped; the next is compared with the survivor");
    }
    {
        std::vector<C> c{{0, 2000, 50, false, true}, {0, 2081, 50, false, true}, {1, 10, 50, true, true},
                         {1, 60, 20, true, false}, {1, 91, 5, true, true}};
        check(ts::cancelBackoffClashes(c, k) == 0 && c[0].closed && c[1].closed && c[2].closed && c[4].closed,
              "T-G4 rule: exactly k apart, other contigs and open gaps never clash");
    }
}

// Flag parsing: only an integer 1..200 turns the back-off on (build_v3: the table refuses the
// malformed values at startup instead of reading them as off; see the header).
void testBackoffParsing() {
    const std::string& G = genome();
    const std::string Gerr = mutate(G, 1985);
    std::vector<Pair> pairs;
    tile(pairs, G, 0, 5000, 5);
    for (int i = 0; i < 2; ++i) pairs.push_back({Gerr.substr(1900, 100), rc(G.substr(2100, 100))});
    const ts::SequenceStore reads = store(pairs);
    const std::vector<std::string> contig{Gerr.substr(0, 2000) + std::string(50, 'N') + G.substr(2050, 1950)};
    {
        const Closed r = closeWith(contig, reads, "0");
        check(r.contigs == contig && r.st.gapsClosed == 0 && r.st.seedBackedOff == 0,
              "TESSERACT_GAPFILL_BACKOFF=\"0\" is off");
    }
    const ts::env::Spec* spec = ts::env::find("TESSERACT_GAPFILL_BACKOFF");
    for (const char* v : {"201", "-5", "60x", "", "abc"}) {
        ts::env::Parsed p;
        std::string why;
        check(spec && !ts::env::parse(*spec, v, p, why),
              std::string("TESSERACT_GAPFILL_BACKOFF=\"") + v + "\" is refused by the flag table");
    }
    const Closed one = closeWith(contig, reads, "200");
    check(one.st.gapsClosed == 1 && one.contigs[0] == G.substr(0, 4000), "TESSERACT_GAPFILL_BACKOFF=200 is on");
    // A contig with no gaps and an empty set both return untouched.
    std::vector<std::string> none{G.substr(0, 3000)};
    const Closed plain = closeWith(none, reads, "60");
    check(plain.contigs == none && plain.st.gapsSeen == 0, "no gaps: untouched");
    const Closed empty = closeWith({}, reads, "60");
    check(empty.contigs.empty(), "no contigs: untouched");
}
}  // namespace

int main() {
    testenv::clearTesseractEnv();
    // 1.4.0: an unset umbrella is on (src/defaults.h); the release 1.3.0 state is TESSERACT_FIXES=0.
    setenv("TESSERACT_FIXES", "0", 1);
    testKeepFlank();
    testTaperedErrorTip();
    testTrueTaper();
    testTwoEqualPaths();
    testOverlappingFlanks();
    testNearbyGaps();
    testClashRule();
    testBackoffParsing();
    std::printf("test_gapfill_backoff: %d checks, %d failures\n", checks, failures);
    return failures ? 1 : 0;
}
