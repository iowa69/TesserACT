// T-M4 (combo2 PKG-MERGE): the terminal dovetail trim, release vs TESSERACT_BOUNDARY_SAFE_TRIM.
//
// Astra Q40's counterexample: k = 99, S1 = C^505 A^100, S2 = A^100 G^600. Release trims the
// whole 100-bp overlap from the shorter S1, which leaves 505 bases (both the 200-bp stub guard
// and QUAST's 500-bp floor are respected) and still holds every removed base in S2. Yet the 98
// flank-spanning 99-mers C^a A^(99-a), a = 1..98, vanish from the assembly. Keeping k-1
// overlap bases on the victim keeps all of them, while the keeper holds the interior k-mers.
//
// Also covered: the mirrored record (input order; equal lengths, where release cuts both
// copies), reverse-complement variants, coupled front and back cuts on one victim, the 200-bp
// stub guard, and 400 random fixtures with planted dovetails in both orientations.
// Standalone component test: `make componenttest`.
//
// build_v3 (integration): the package's own trim (ts::trimTerminalDovetails in assembler.cpp)
// and G-emit's T10 (ts::trimTerminalOverlaps in emit_post.cpp) were two implementations of the
// same fix; build_v3 keeps ONE, T10, and TESSERACT_BOUNDARY_SAFE_TRIM is an alias of
// TESSERACT_FIX_BOUNDARY_SAFE_TRIM. This test now drives T10 through a local adapter with the
// package's signature, so every package expectation below checks the kept implementation.
// T10 couples the package's two knobs (one victim per overlap AND keep k-1), so the adapter
// maps (0, false) to release and (keep, true) to boundarySafe with k = keep + 1; the one
// package call with one-victim alone (case 5, keep 0) has no T10 form and was removed.
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <set>
#include <string>
#include <vector>

#include "assembler.h"
#include "emit_post.h"
#include "seqio.h"
#include "test_env.h"

namespace {

int g_checks = 0, g_failures = 0;
#define CHECK(cond)                                                                   \
    do {                                                                              \
        ++g_checks;                                                                   \
        if (!(cond)) {                                                                \
            ++g_failures;                                                             \
            std::printf("  FAIL %s:%d  %s\n", __FILE__, __LINE__, #cond);             \
        }                                                                             \
    } while (0)

constexpr size_t K = 99;
constexpr size_t KEEP = K - 1;

}  // namespace

namespace ts {
// The package's stats and signature, over the single build_v3 implementation (T10).
struct TerminalTrimStats {
    size_t trimmed = 0, trimmedBases = 0, keptBoundaryBp = 0, equalLenFlips = 0;
};
TerminalTrimStats trimTerminalDovetails(std::vector<std::string>& seqs, size_t minOv,
                                        size_t keepBoundary, bool oneVictimPerOverlap) {
    if ((keepBoundary == 0) == oneVictimPerOverlap) {
        std::printf("  FAIL adapter: keepBoundary=%zu oneVictim=%d has no T10 form\n", keepBoundary,
                    oneVictimPerOverlap ? 1 : 0);
        std::exit(1);
    }
    TrimConfig cfg;
    cfg.minOverlap = minOv;
    cfg.k = static_cast<int>(keepBoundary + 1);
    cfg.boundarySafe = oneVictimPerOverlap;
    cfg.copyGuard = false;
    std::vector<size_t> cutFront, cutBack;
    const TrimStats st = trimTerminalOverlaps(seqs, cfg, cutFront, cutBack, nullptr);
    TerminalTrimStats out;
    out.trimmed = st.trimmed;
    out.trimmedBases = st.trimmedBases;
    out.keptBoundaryBp = st.boundaryBpKept;
    out.equalLenFlips = st.equalLenFlips;
    return out;
}
}  // namespace ts

namespace {

std::string rep(char c, size_t n) { return std::string(n, c); }

std::string canon(const std::string& s) {
    const std::string r = ts::reverseComplement(s);
    return std::min(s, r);
}

std::set<std::string> kmerClasses(const std::vector<std::string>& seqs, size_t k) {
    std::set<std::string> out;
    for (const std::string& s : seqs)
        for (size_t i = 0; i + k <= s.size(); ++i) out.insert(canon(s.substr(i, k)));
    return out;
}

// Classes present before and absent after.
std::set<std::string> lost(const std::vector<std::string>& before,
                           const std::vector<std::string>& after, size_t k) {
    const std::set<std::string> a = kmerClasses(before, k), b = kmerClasses(after, k);
    std::set<std::string> out;
    for (const std::string& x : a) if (!b.count(x)) out.insert(x);
    return out;
}

// The 98 flank-spanning classes C^a A^(99-a), a = 1..98.
std::set<std::string> flankClasses() {
    std::set<std::string> out;
    for (size_t a = 1; a <= 98; ++a) out.insert(canon(rep('C', a) + rep('A', K - a)));
    return out;
}

// Random contigs with planted exact dovetails, some longer than k and some not, in both
// orientations, and now and then two partners of equal length. (That keepBoundary = 0 is the
// release consumer is checked on real data by combo2/merge/tools/trim_replay, which
// reproduces the release binary's contigs.fasta from its scaffolds.fasta.)
struct Planted { std::vector<std::string> seqs; };

Planted randomFixture(std::mt19937_64& rng) {
    auto rnd = [&](size_t n) {
        std::string s(n, 'A');
        for (char& c : s) c = "ACGT"[rng() & 3];
        return s;
    };
    Planted p;
    const size_t n = 2 + rng() % 4;
    for (size_t i = 0; i < n; ++i) p.seqs.push_back(rnd(300 + rng() % 900));
    const size_t plants = 1 + rng() % 4;
    for (size_t t = 0; t < plants; ++t) {
        const size_t a = rng() % n, b = rng() % n;
        if (a == b) continue;
        const size_t len = 60 + rng() % 300;
        if (len + 10 >= p.seqs[a].size() || len + 10 >= p.seqs[b].size()) continue;
        const std::string ov = p.seqs[a].substr(p.seqs[a].size() - len);
        const bool equalLen = rng() % 4 == 0;
        std::string& sb = p.seqs[b];
        const size_t L = p.seqs[a].size();
        if (rng() & 1) {                                   // a's 3' end meets b's 5' end
            sb.replace(0, len, ov);
            if (equalLen) { if (sb.size() > L) sb.resize(L); else sb += rnd(L - sb.size()); }
        } else {                                           // a's 3' end meets rc(b)'s 5' end
            sb.replace(sb.size() - len, len, ts::reverseComplement(ov));
            if (equalLen) {
                if (sb.size() > L) sb = sb.substr(sb.size() - L);
                else sb = rnd(L - sb.size()) + sb;
            }
        }
    }
    return p;
}

}  // namespace

int main() {
    testenv::clearTesseractEnv();
    std::printf("test_boundary_safe_trim\n");
    const std::string S1 = rep('C', 505) + rep('A', 100);
    const std::string S2 = rep('A', 100) + rep('G', 600);

    // ---- 1. Astra Q40, release: the whole overlap is cut and 98 classes vanish ----------
    {
        std::vector<std::string> v = {S1, S2};
        const ts::TerminalTrimStats st = ts::trimTerminalDovetails(v, K, 0, false);
        CHECK(st.trimmed == 1);
        CHECK(st.trimmedBases == 100);
        CHECK(st.keptBoundaryBp == 0);
        CHECK(v[0] == rep('C', 505));
        CHECK(v[1] == S2);
        CHECK(v[0].size() >= 500);                          // QUAST floor respected
        const std::set<std::string> gone = lost({S1, S2}, v, K);
        CHECK(gone.size() == 98);
        CHECK(gone == flankClasses());
        // every removed base still has a holder: A^100 is in S2
        CHECK(v[1].find(rep('A', 100)) != std::string::npos);
        std::printf("  info release loses %zu flank-spanning 99-mer classes\n", gone.size());
    }
    // ---- 2. Astra Q40, boundary-safe: all 98 survive -------------------------------------
    {
        std::vector<std::string> v = {S1, S2};
        const ts::TerminalTrimStats st = ts::trimTerminalDovetails(v, K, KEEP, true);
        CHECK(st.trimmed == 1);
        CHECK(st.trimmedBases == 2);
        CHECK(st.keptBoundaryBp == 98);
        CHECK(v[0] == rep('C', 505) + rep('A', 98));
        CHECK(v[1] == S2);
        CHECK(lost({S1, S2}, v, K).empty());
        const std::set<std::string> after = kmerClasses(v, K);
        for (const std::string& f : flankClasses()) CHECK(after.count(f) == 1);
        // no class is created either: the output is made of substrings of the input
        const std::set<std::string> before = kmerClasses({S1, S2}, K);
        CHECK(after == before);
        std::printf("  info boundary-safe keeps all 98 (cut %zu bp, kept %zu bp)\n",
                    st.trimmedBases, st.keptBoundaryBp);
    }
    // ---- 3. mirrored record: input order swapped (the record reaching the victim is the
    //         mirror probed from rc(S2)) gives the same single cut ----------------------------
    for (int mode = 0; mode < 2; ++mode) {
        std::vector<std::string> v = {S2, S1};
        const ts::TerminalTrimStats st =
            ts::trimTerminalDovetails(v, K, mode ? KEEP : 0, mode == 1);
        CHECK(st.trimmed == 1);
        CHECK(v[0] == S2);
        CHECK(v[1] == (mode ? rep('C', 505) + rep('A', 98) : rep('C', 505)));
        CHECK(lost({S2, S1}, v, K).size() == (mode ? 0u : 98u));
    }
    // ---- 4. reverse-complement variants: S1 reversed, S2 reversed, both reversed -----------
    for (int variant = 1; variant < 4; ++variant) {
        const std::string a = (variant & 1) ? ts::reverseComplement(S1) : S1;
        const std::string b = (variant & 2) ? ts::reverseComplement(S2) : S2;
        for (int mode = 0; mode < 2; ++mode) {
            std::vector<std::string> v = {a, b};
            const ts::TerminalTrimStats st =
                ts::trimTerminalDovetails(v, K, mode ? KEEP : 0, mode == 1);
            CHECK(st.trimmed == 1);
            CHECK(v[1] == b);                                // the longer contig is the keeper
            CHECK(v[0].size() == (mode ? 505u + 98u : 505u));
            const std::set<std::string> gone = lost({a, b}, v, K);
            if (mode == 0) CHECK(gone == flankClasses());
            else CHECK(gone.empty());
        }
    }
    // ---- 5. mirrored records of EQUAL length: release cuts BOTH copies of the overlap -----
    {
        const std::string E1 = rep('C', 600) + rep('A', 100);    // 700
        const std::string E2 = rep('A', 100) + rep('G', 600);    // 700
        std::vector<std::string> v = {E1, E2};
        const ts::TerminalTrimStats rel = ts::trimTerminalDovetails(v, K, 0, false);
        CHECK(rel.trimmed == 2);                                 // release defect, reproduced
        CHECK(v[0] == rep('C', 600));
        CHECK(v[1] == rep('G', 600));
        CHECK(v[0].find('A') == std::string::npos && v[1].find('A') == std::string::npos);

        std::vector<std::string> w = {E1, E2};
        const ts::TerminalTrimStats safe = ts::trimTerminalDovetails(w, K, KEEP, true);
        CHECK(safe.trimmed == 1);                                // one victim per overlap
        CHECK(safe.equalLenFlips == 1);
        CHECK(w[0] == E1);                                       // lower index keeps
        CHECK(w[1] == rep('A', 98) + rep('G', 600));
        CHECK(lost({E1, E2}, w, K).empty());
        // (build_v3: the package's "one-victim alone, keepBoundary 0" call has no T10 form;
        // see the header.)
        std::printf("  info equal lengths: release cuts %zu copies, safe cuts %zu\n",
                    rel.trimmed, safe.trimmed);
    }
    // ---- 6. coupled front and back cuts on one victim ------------------------------------
    // Random flanks (homopolymer fixtures would add spurious overlaps between the keepers).
    std::mt19937_64 frng(0x5eed40);
    auto rnd = [&](size_t n) {
        std::string s(n, 'A');
        for (char& c : s) c = "ACGT"[frng() & 3];
        return s;
    };
    const std::string X = rnd(1000), Y = rnd(1100), core = rnd(400), core2 = rnd(160);
    {
        const std::string V = X.substr(1000 - 150) + core + Y.substr(0, 120);   // 670: shorter
        std::vector<std::string> v = {X, V, Y};
        const ts::TerminalTrimStats rel = ts::trimTerminalDovetails(v, K, 0, false);
        CHECK(rel.trimmed == 2);
        CHECK(v[1] == core);
        std::vector<std::string> w = {X, V, Y};
        const ts::TerminalTrimStats safe = ts::trimTerminalDovetails(w, K, KEEP, true);
        CHECK(safe.trimmed == 2);
        CHECK(safe.keptBoundaryBp == 2 * KEEP);
        CHECK(safe.trimmedBases == (150 - KEEP) + (120 - KEEP));
        CHECK(w[1] == X.substr(1000 - KEEP) + core + Y.substr(0, KEEP));
        CHECK(w[0] == X && w[2] == Y);
        CHECK(lost({X, V, Y}, w, K).empty());
        CHECK(lost({X, V, Y}, v, K).size() == 2 * KEEP);   // 98 flank classes at each end
        // the same victim entered reverse-complemented
        const std::string Vr = ts::reverseComplement(V);
        std::vector<std::string> r = {X, Vr, Y};
        CHECK(ts::trimTerminalDovetails(r, K, KEEP, true).trimmed == 2);
        CHECK(r[1] == ts::reverseComplement(w[1]));
        CHECK(lost({X, Vr, Y}, r, K).empty());
    }
    // ---- 7. the 200-bp stub guard decides on the FULL overlaps in both modes ---------------
    {
        // 100-bp overlap on a 300-bp victim: 100 + 200 >= 300 -> never cut, even though the
        // boundary-safe cut (2 bp) would leave 298 bp.
        const std::string A = rep('C', 200) + rep('A', 100);      // 300
        std::vector<std::string> v = {A, S2};
        CHECK(ts::trimTerminalDovetails(v, K, 0, false).trimmed == 0);
        CHECK(v[0] == A);
        std::vector<std::string> w = {A, S2};
        CHECK(ts::trimTerminalDovetails(w, K, KEEP, true).trimmed == 0);
        CHECK(w[0] == A);
        // 301 bp: 100 + 200 < 301 -> cut in both modes
        const std::string B = rep('C', 201) + rep('A', 100);
        std::vector<std::string> x = {B, S2};
        CHECK(ts::trimTerminalDovetails(x, K, 0, false).trimmed == 1);
        CHECK(x[0] == rep('C', 201));
        std::vector<std::string> y = {B, S2};
        CHECK(ts::trimTerminalDovetails(y, K, KEEP, true).trimmed == 1);
        CHECK(y[0] == rep('C', 201) + rep('A', 98));
        // coupled guard: the second cut counts the first one's full length, in both modes.
        // V is 430 bp: front 150 is decided first (longer), then back 120 + 150 + 200 >= 430
        // is skipped, even though the boundary-safe cuts (52 + 22) would leave 356 bp.
        const std::string V = X.substr(1000 - 150) + core2 + Y.substr(0, 120);
        std::vector<std::string> p = {X, V, Y}, q = {X, V, Y};
        CHECK(ts::trimTerminalDovetails(p, K, 0, false).trimmed == 1);
        CHECK(ts::trimTerminalDovetails(q, K, KEEP, true).trimmed == 1);
        CHECK(p[1] == core2 + Y.substr(0, 120));
        CHECK(q[1] == X.substr(1000 - KEEP) + core2 + Y.substr(0, 120));
        CHECK(lost({X, V, Y}, q, K).empty());
    }
    // ---- 8. overlaps not longer than k are never touched; a cut <= keep is skipped --------
    {
        const std::string A = rep('C', 505) + rep('A', 99);
        const std::string B = rep('A', 99) + rep('G', 600);
        std::vector<std::string> v = {A, B};
        CHECK(ts::trimTerminalDovetails(v, K, 0, false).trimmed == 0);   // L must exceed minOv
        // keepBoundary >= overlap: the decided cut becomes <= 0 and is skipped (T10 with k = 101)
        std::vector<std::string> w = {S1, S2};
        const ts::TerminalTrimStats st = ts::trimTerminalDovetails(w, K, 100, true);
        CHECK(st.trimmed == 0);
        CHECK(w[0] == S1 && w[1] == S2);
    }
    // ---- 9. randomized: keep 0 / one-victim off is the release consumer; safe mode only ever
    //         cuts a suffix of what release cuts, never more, and is idempotent on its input ----
    {
        std::mt19937_64 rng(0xb0d5a7e);
        size_t fixtures = 0, cutFixtures = 0, lostRel = 0, lostSafe = 0;
        for (int it = 0; it < 400; ++it) {
            const Planted p = randomFixture(rng);
            std::vector<std::string> rel = p.seqs, safe = p.seqs;
            const ts::TerminalTrimStats a = ts::trimTerminalDovetails(rel, K, 0, false);
            const ts::TerminalTrimStats b = ts::trimTerminalDovetails(safe, K, KEEP, true);
            ++fixtures;
            if (a.trimmed) ++cutFixtures;
            CHECK(b.trimmed <= a.trimmed);
            for (size_t i = 0; i < p.seqs.size(); ++i) {
                // each output is a substring of its input and the safe one contains release's
                CHECK(p.seqs[i].find(rel[i]) != std::string::npos);
                CHECK(p.seqs[i].find(safe[i]) != std::string::npos);
                CHECK(safe[i].size() >= rel[i].size() || b.equalLenFlips > 0);
                if (safe[i].size() < p.seqs[i].size()) CHECK(safe[i].size() > 200);
            }
            lostRel += lost(p.seqs, rel, K).size();
            lostSafe += lost(p.seqs, safe, K).size();
        }
        CHECK(cutFixtures > 50);
        CHECK(lostSafe == 0);
        CHECK(lostRel > 0);
        std::printf("  info randomized: %zu fixtures, %zu with cuts; 99-mer classes lost: "
                    "release %zu, boundary-safe %zu\n", fixtures, cutFixtures, lostRel, lostSafe);
    }

    std::printf("%s test_boundary_safe_trim: %d checks, %d failures\n",
                g_failures ? "FAIL" : "ok", g_checks, g_failures);
    return g_failures ? 1 : 0;
}
