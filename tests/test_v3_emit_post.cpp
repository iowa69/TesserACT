// G-emit emission-stage fixes, unit level (src/emit_post.*):
//   T37 dedup prefilter          -- decisions identical to the release loop, far fewer finds
//   T10 boundary-safe trim       -- no lost k-mer classes; equal-length partners cut once
//   T22 copy guard               -- a segment at two loci keeps two copies, with evidence
//   T25 split post-processing    -- per-piece coverage, circular tag, 2k floor, dedup
//   T26 AGP / P-lines v2         -- rows name written records and tile every base
//   T17 contig statistics        -- describe the written records
// Every fix is exercised with its switch off (release behaviour asserted) and on.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <map>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>
#include <unistd.h>

#include "emit_post.h"
#include "graph.h"
#include "report.h"
#include "test_env.h"   // build_v3: T40, clear the ambient TESSERACT_* first

namespace {
int checks = 0, failures = 0;
void check(bool ok, const std::string& what) {
    ++checks;
    std::printf("%s %s\n", ok ? "ok  " : "FAIL", what.c_str());
    if (!ok) ++failures;
}
std::mt19937 rng(424242);
std::string dna(size_t n) {
    std::string s(n, 'A');
    for (char& c : s) c = "ACGT"[rng() % 4];
    return s;
}
std::string rc(const std::string& s) { return ts::reverseComplement(s); }

// ---- T37: the release loop, verbatim, as the oracle ---------------------------------------
std::vector<char> releaseDedup(const std::vector<std::string>& seqs, const std::vector<char>& prot) {
    std::vector<char> contained(seqs.size(), 0);
    std::vector<size_t> byLen(seqs.size());
    for (size_t i = 0; i < byLen.size(); ++i) byLen[i] = i;
    std::sort(byLen.begin(), byLen.end(), [&](size_t a, size_t b) {
        if (seqs[a].size() != seqs[b].size()) return seqs[a].size() > seqs[b].size();
        return seqs[a] < seqs[b];
    });
    std::string text;
    for (size_t idx : byLen) {
        const std::string& s = seqs[idx];
        bool hit = false;
        if (!s.empty() && !text.empty() && !prot[idx]) {
            hit = text.find(s) != std::string::npos;
            if (!hit) hit = text.find(rc(s)) != std::string::npos;
        }
        if (hit) contained[idx] = 1;
        else { text.push_back('\x01'); text.append(s); }
    }
    return contained;
}

void testDedup() {
    size_t diffs = 0, queries = 0;
    for (int round = 0; round < 6; ++round) {
        std::vector<std::string> seqs;
        for (int i = 0; i < 40; ++i) seqs.push_back(dna(200 + rng() % 3000));
        const size_t base = seqs.size();
        for (int i = 0; i < 160; ++i) {
            const std::string& h = seqs[rng() % base];
            const size_t len = 1 + rng() % std::min<size_t>(h.size(), 600);
            const size_t st = rng() % (h.size() - len + 1);
            std::string x = h.substr(st, len);
            switch (rng() % 9) {
                case 0: break;                                   // forward substring
                case 1: x = rc(x); break;                        // reverse-complement substring
                case 2: x[rng() % x.size()] = 'T'; break;        // near miss (or not)
                case 3: x = h; break;                            // forward twin
                case 4: x = rc(h); break;                        // rc twin
                case 5: x = x.substr(0, std::min<size_t>(x.size(), 20)); break;   // < one window
                case 6: x.insert(x.size() / 2, "NNNNN"); break;  // N-run
                case 7: for (char& c : x) if (c == 'A') c = 'n'; break;          // lowercase n
                case 8: x = seqs[rng() % base].substr(0, 40) + seqs[rng() % base].substr(0, 40); break;
            }
            seqs.push_back(x);
        }
        seqs.push_back("");
        seqs.push_back(std::string(50, 'N'));
        std::vector<char> prot(seqs.size(), 0);
        for (size_t i = 0; i < seqs.size(); ++i) prot[i] = (rng() % 10) == 0;
        ts::DedupStats st;
        const auto a = ts::markContained(seqs, prot, nullptr, &st);
        const auto b = releaseDedup(seqs, prot);
        for (size_t i = 0; i < seqs.size(); ++i) diffs += a[i] != b[i];
        queries += seqs.size();
    }
    check(diffs == 0, "T37: prefiltered dedup decisions identical to the release loop (" +
                          std::to_string(queries) + " sequences)");
    // Performance: thousands of non-contained contigs against one large text.
    std::vector<std::string> many{dna(2000000)};
    for (int i = 0; i < 3000; ++i) many.push_back(dna(250));
    ts::DedupStats st;
    const auto c = ts::markContained(many, std::vector<char>(many.size(), 0), nullptr, &st);
    size_t dropped = 0;
    for (char x : c) dropped += x;
    check(st.queried == 3000 && st.exactSearches < 60 && dropped == 0,
          "T37: 3000 non-contained queries run " + std::to_string(st.exactSearches) +
              " exact finds (release: 6000)");
}

// ---- k-mer class bookkeeping ---------------------------------------------------------------
std::unordered_set<std::string> classes(const std::vector<std::string>& seqs, int k) {
    std::unordered_set<std::string> out;
    for (const std::string& s : seqs)
        for (size_t i = 0; i + static_cast<size_t>(k) <= s.size(); ++i) {
            const std::string f = s.substr(i, static_cast<size_t>(k)), r = rc(f);
            out.insert(std::min(f, r));
        }
    return out;
}
size_t lost(const std::vector<std::string>& before, const std::vector<std::string>& after, int k) {
    const auto a = classes(before, k), b = classes(after, k);
    size_t n = 0;
    for (const auto& x : a) n += b.count(x) ? 0 : 1;
    return n;
}
size_t occurrences(const std::vector<std::string>& seqs, const std::string& x) {
    size_t n = 0;
    const std::string r = rc(x);
    for (const std::string& s : seqs) {
        for (size_t p = s.find(x); p != std::string::npos; p = s.find(x, p + 1)) ++n;
        if (r != x) for (size_t p = s.find(r); p != std::string::npos; p = s.find(r, p + 1)) ++n;
    }
    return n;
}
ts::TrimStats trim(std::vector<std::string>& seqs, int k, bool safe, bool guard,
                   const ts::EndEvidenceFn* ev = nullptr) {
    ts::TrimConfig cfg;
    cfg.minOverlap = static_cast<size_t>(k);
    cfg.k = k;
    cfg.boundarySafe = safe;
    cfg.copyGuard = guard;
    std::vector<size_t> f, b;
    return ts::trimTerminalOverlaps(seqs, cfg, f, b, ev);
}

void testBoundarySafe() {
    const int k = 99;
    // Case 1: X(505)+R(100) against R+Y(600): the release cut drops the 98 classes spanning X|R.
    {
        const std::string X = dna(505), R = dna(100), Y = dna(600);
        const std::vector<std::string> in{X + R, R + Y};
        std::vector<std::string> rel = in, fix = in;
        trim(rel, k, false, false);
        trim(fix, k, true, false);
        check(rel[0] == X && lost(in, rel, k) == 98, "T10 off: release cuts the whole overlap, 98 classes lost");
        check(fix[0].size() == X.size() + 98 && lost(in, fix, k) == 0,
              "T10 on: the victim keeps k-1 overlap bases, 0 classes lost");
    }
    // Case 2: equal lengths -- release cuts the overlap from both partners.
    {
        const std::string X = dna(600), R = dna(150), Y = dna(600);
        const std::vector<std::string> in{X + R, R + Y};
        std::vector<std::string> rel = in, fix = in;
        trim(rel, k, false, false);
        const ts::TrimStats st = trim(fix, k, true, false);
        check(rel[0] == X && rel[1] == Y && occurrences(rel, R) == 0,
              "T10 off: equal-length partners both cut, the overlap vanishes");
        check(occurrences(fix, R) == 1 && lost(in, fix, k) == 0 && st.trimmed == 1,
              "T10 on: one victim per physical overlap, overlap kept once, 0 classes lost");
    }
    // Random fixtures, 4 orientations each.
    size_t relLost = 0, fixLost = 0, fixAbsent = 0;
    for (int i = 0; i < 50; ++i) {
        const std::string X = dna(300 + rng() % 700), R = dna(100 + rng() % 300), Y = dna(300 + rng() % 700);
        for (int o = 0; o < 4; ++o) {
            std::vector<std::string> in{X + R, R + Y};
            if (o & 1) in[0] = rc(in[0]);
            if (o & 2) in[1] = rc(in[1]);
            std::vector<std::string> rel = in, fix = in;
            trim(rel, k, false, false);
            trim(fix, k, true, false);
            relLost += lost(in, rel, k);
            fixLost += lost(in, fix, k);
            fixAbsent += occurrences(fix, R) == 0;
        }
    }
    check(relLost > 0 && fixLost == 0 && fixAbsent == 0,
          "T10: 200 random fixtures, release loses " + std::to_string(relLost) + " classes, fix 0");
}

// The copy guard writes one audit line per end to stderr; keep the test output readable.
struct QuietStderr {
    int saved = -1;
    QuietStderr() { std::fflush(stderr); saved = dup(2); std::FILE* f = std::freopen("/dev/null", "w", stderr); (void)f; }
    ~QuietStderr() { std::fflush(stderr); dup2(saved, 2); close(saved); clearerr(stderr); }
};

void testCopyGuard() {
    QuietStderr quiet;
    const int k = 127;
    const std::string R = dna(1048);
    std::vector<std::string> flank{dna(1500), dna(1800), dna(2100), dna(2400)};
    ts::EndEvidenceFn supported = [](const std::vector<ts::EndQuery>& q) {
        std::vector<ts::EndEvidence> out(q.size());
        for (auto& e : out) { e.consistent = 12; e.copyRatio = 2.0; }
        return out;
    };
    ts::EndEvidenceFn oneCopy = [](const std::vector<ts::EndQuery>& q) {
        std::vector<ts::EndEvidence> out(q.size());
        for (auto& e : out) { e.consistent = 12; e.copyRatio = 1.0; }
        return out;
    };
    size_t cases = 0, relOnce = 0, guardTwice = 0, noEvidenceSame = 0, lowRatioSame = 0, lossless = 0, lossChecked = 0;
    std::vector<int> perm{0, 1, 2, 3};
    int permIdx = 0;
    do {
        if (permIdx++ % 4) continue;                 // 6 of the 24 length orders
        // A = X R, B = R Y, C = Z R, D = R W: R at two loci, all four flanks walked into it.
        const std::vector<std::string> base{flank[perm[0]] + R, R + flank[perm[1]], flank[perm[2]] + R,
                                            R + flank[perm[3]]};
        std::vector<int> order{0, 1, 2, 3};
        int orderIdx = 0;
        do {
            if (orderIdx++ % 3) continue;            // 8 of the 24 output orders
            for (int o = 0; o < 16; o += 3) {
                std::vector<std::string> in;
                for (int j = 0; j < 4; ++j) in.push_back((o >> j) & 1 ? rc(base[order[j]]) : base[order[j]]);
                std::vector<std::string> rel = in, g = in, none = in, low = in;
                trim(rel, k, false, false);
                trim(g, k, false, true, &supported);
                trim(none, k, false, true, nullptr);
                trim(low, k, false, true, &oneCopy);
                ++cases;
                relOnce += occurrences(rel, R) == 1;
                guardTwice += occurrences(g, R) == 2;
                noEvidenceSame += none == rel;
                lowRatioSame += low == rel;
                if (cases % 32 == 1) { ++lossChecked; lossless += lost(in, g, 31) <= lost(in, rel, 31); }
            }
        } while (std::next_permutation(order.begin(), order.end()));
    } while (std::next_permutation(perm.begin(), perm.end()));
    check(relOnce > 0, "T22 off: release keeps a 2-locus segment once in " + std::to_string(relOnce) + "/" +
                           std::to_string(cases) + " cases");
    check(guardTwice == cases, "T22 on, supported: two copies kept in " + std::to_string(guardTwice) + "/" +
                                   std::to_string(cases));
    check(noEvidenceSame == cases, "T22 on without evidence: release decisions in every case");
    check(lowRatioSame == cases, "T22 on, depth says one copy: release decisions in every case");
    check(lossChecked > 0 && lossless == lossChecked, "T22 on: never loses more k-mer classes than release (" + std::to_string(lossChecked) + " sampled cases)");
    // Copy-number-1 control: identical to release even with supporting evidence.
    {
        const std::string X = dna(1600), Y = dna(1900);
        const std::vector<std::string> in{X + R, R + Y};
        std::vector<std::string> rel = in, g = in;
        trim(rel, k, false, false);
        trim(g, k, false, true, &supported);
        check(g == rel, "T22 on: a one-copy overlap is trimmed exactly as release");
    }
    // Wrong walk longer than the right one: X R W[:600] against R W. Release cuts the correct
    // copy (the shorter). Evidence against the chimeric end hands the copy to the other side.
    {
        const std::string X = dna(1500), W = dna(1500);
        const std::vector<std::string> in{X + R + W.substr(0, 600), R + W};
        ts::EndEvidenceFn caseEv = [](const std::vector<ts::EndQuery>& q) {
            std::vector<ts::EndEvidence> out(q.size());
            for (size_t i = 0; i < q.size(); ++i) {
                if (q[i].seq == 0) { out[i].inconsistent = 9; out[i].consistent = 1; }
                else { out[i].consistent = 11; }
                out[i].copyRatio = 1.0;
            }
            return out;
        };
        std::vector<std::string> rel = in, g = in;
        trim(rel, k, false, false);
        const ts::TrimStats st = trim(g, k, false, true, &caseEv);
        check(rel[0] == in[0] && rel[1] == W.substr(600), "T22 off: the correct (shorter) copy is cut");
        check(g[0] == X && g[1] == R + W && st.evidenceSwitched == 1,
              "T22 on: contradicted chimeric end is cut, the supported copy kept");
    }
    // A tandem part-copy of R inside one flank: the assembly already holds more of R than the
    // reads' two copies allow once both entering ends keep theirs -- no rescue.
    {
        const std::string Zt = flank[2] + R.substr(448);
        const std::vector<std::string> in{flank[0] + R, R + flank[1], Zt + R, R + flank[3]};
        std::vector<std::string> rel = in, g = in;
        trim(rel, k, false, false);
        const ts::TrimStats st = trim(g, k, false, true, &supported);
        check(g == rel && st.guardApplied == 0 && st.evidenceRefused == 1,
              "T22 on: interior copies push the kept count past the reads' copy number -> release decision");
    }
    // Composition with T10: boundary-safe cuts under the guard.
    {
        const std::vector<std::string> in{flank[0] + R, R + flank[1], flank[2] + R, R + flank[3]};
        std::vector<std::string> g = in;
        trim(g, k, true, true, &supported);
        check(occurrences(g, R) == 2 && lost(in, g, k) == 0, "T10+T22: two copies, no k-mer class lost");
    }
}

// ---- T25 / T26: a tiny graph ----------------------------------------------------------------
ts::UnitigGraph tinyGraph(const std::vector<std::pair<std::string, double>>& unitigs, int k) {
    ts::UnitigGraph g;
    g.setK(k);
    for (const auto& u : unitigs) {
        ts::Unitig n;
        n.seq = u.first;
        n.coverage = u.second;
        g.nodes.push_back(n);
    }
    return g;
}

void testSplitPost() {
    const int k = 31;
    const std::string u0 = dna(200), u1 = dna(150);
    const ts::UnitigGraph g = tinyGraph({{u0, 60.0}, {u1, 180.0}}, k);
    ts::GfaPath path;
    path.oriented = {0 << 1, 1 << 1};
    path.gaps = {0, 5};
    const std::string scaffold = u0 + "NNNNN" + u1.substr(k - 1);   // the release gap-flank rule
    const double scafCov = (60.0 * 200 + 180.0 * 150) / 350.0;
    // Case A: coverage, circular tag, short piece (minLen 150 drops the 120-bp piece).
    {
        auto pieces = ts::splitAtGaps({scaffold}, {"_unk_circular"}, {scafCov});
        ts::SplitPostStats off;
        auto p0 = pieces;
        ts::splitPostprocess(p0, {scaffold}, {path}, &g, {0}, 150, off);
        check(p0.size() == 2 && p0[0].cov == scafCov && p0[1].cov == scafCov && p0[1].tag == "_unk_circular" &&
                  !p0[0].drop && !p0[1].drop && off.pieces == 2 && off.multiPieceScaffolds == 1,
              "T25 off: both pieces inherit the scaffold cov and _circular, none dropped");
        ts::SplitPostStats on;
        on.enabled = true;
        ts::splitPostprocess(pieces, {scaffold}, {path}, &g, {0}, 150, on);
        check(pieces[0].cov == 60.0 && pieces[1].cov == 180.0 && on.covRelabelled == 2,
              "T25 on: each piece gets its own segment's coverage (60, 180)");
        check(pieces[0].tag == "_unk" && pieces[1].tag == "_unk" && on.circularTagsDropped == 2,
              "T25 on: _circular dropped from the pieces of a multi-piece scaffold");
        check(!pieces[0].drop && pieces[1].drop == 1 && on.droppedShort == 1 && on.droppedShortBp == 120,
              "T25 on: the 120-bp piece is below the 150-bp floor and dropped");
        // The gap-flank fix renders u1 whole; the mapping must still find it.
        const std::string fixedScaffold = u0 + "NN" + u1;
        auto q = ts::splitAtGaps({fixedScaffold}, {"_chr"}, {scafCov});
        ts::SplitPostStats on2;
        on2.enabled = true;
        ts::splitPostprocess(q, {fixedScaffold}, {path}, &g, {0}, 62, on2);
        check(q[0].cov == 60.0 && q[1].cov == 180.0, "T25 on: coverage mapping also under the whole-flank render");
        // A closed gap merges two segments into one piece: a single piece keeps its label.
        const std::string closedScaffold = u0 + "ACGTACG" + u1.substr(k - 1);
        auto c = ts::splitAtGaps({closedScaffold}, {"_chr"}, {scafCov});
        ts::SplitPostStats on3;
        on3.enabled = true;
        ts::splitPostprocess(c, {closedScaffold}, {path}, &g, {0}, 62, on3);
        check(c.size() == 1 && c[0].cov == scafCov && on3.covRelabelled == 0, "T25 on: a single piece keeps its label");
        // Two gaps, the first closed: pieces {seg0+seg1, seg2}.
        const std::string u2 = dna(180);
        const ts::UnitigGraph g3 = tinyGraph({{u0, 60.0}, {u1, 180.0}, {u2, 30.0}}, k);
        ts::GfaPath p3;
        p3.oriented = {0 << 1, 1 << 1, 2 << 1};
        p3.gaps = {0, 5, 8};
        const std::string sc3 = u0 + "ACGTACG" + u1.substr(k - 1) + "NNNNNNNN" + u2.substr(k - 1);
        auto d = ts::splitAtGaps({sc3}, {"_chr"}, {1.0});
        ts::SplitPostStats on4;
        on4.enabled = true;
        ts::splitPostprocess(d, {sc3}, {p3}, &g3, {0}, 62, on4);
        check(d.size() == 2 && std::abs(d[0].cov - scafCov) < 1e-9 && d[1].cov == 30.0,
              "T25 on: a closed gap's two segments are pooled into one piece's coverage");
    }
    // Case B: a piece that is an exact substring (rc) of a standalone record is dropped,
    // unless its scaffold is a protected call.
    {
        const std::string X = dna(1000), P1 = dna(700);
        const std::string P2 = rc(X.substr(200, 300));
        const std::vector<std::string> sc{X, P1 + "NNNNNNNNNN" + P2};
        for (int prot = 0; prot < 2; ++prot) {
            auto pieces = ts::splitAtGaps(sc, {"_chr", "_chr"}, {10.0, 10.0});
            ts::SplitPostStats on;
            on.enabled = true;
            ts::splitPostprocess(pieces, sc, {ts::GfaPath(), ts::GfaPath()}, nullptr,
                                 {0, static_cast<char>(prot)}, 62, on);
            if (!prot)
                check(pieces.size() == 3 && pieces[2].drop == 2 && on.droppedDup == 1 && on.droppedDupBp == 300,
                      "T25 on: a 300-bp piece contained (rc) in a standalone record is dropped");
            else
                check(pieces[2].drop == 0 && on.droppedDup == 0, "T25 on: a protected scaffold's piece is kept");
        }
    }
    // Case C: the T02 shape -- a partially overwritten N-run leaves 1- and 3-bp pieces.
    {
        const std::string A = dna(400), B = dna(400);
        const std::string sc = A + "N" + "C" + "NNNN" + "TTT" + "N" + B;
        auto pieces = ts::splitAtGaps({sc}, {"_chr"}, {10.0});
        ts::SplitPostStats off, on;
        auto p0 = pieces;
        ts::splitPostprocess(p0, {sc}, {ts::GfaPath()}, nullptr, {0}, 62, off);
        on.enabled = true;
        ts::splitPostprocess(pieces, {sc}, {ts::GfaPath()}, nullptr, {0}, 62, on);
        size_t tinyOff = 0, tinyOn = 0;
        for (const auto& p : p0) tinyOff += (!p.drop && p.seq.size() <= 3);
        for (const auto& p : pieces) tinyOn += (!p.drop && p.seq.size() <= 3);
        check(tinyOff == 2 && tinyOn == 0 && on.droppedShort == 2, "T25 on: the 1-bp and 3-bp records are gone");
    }
}

// Parses an AGP and checks every row against the scaffolds and the written records.
struct AgpCheck { size_t w = 0, n = 0, badSeq = 0, badTile = 0, unknownComp = 0, compInRecords = 0; };
AgpCheck checkAgp(const std::string& path, const std::vector<std::string>& scaffolds,
                  const std::vector<std::string>& scafNames, const std::map<std::string, std::string>& records) {
    AgpCheck c;
    std::map<std::string, std::string> obj;
    for (size_t i = 0; i < scaffolds.size(); ++i) obj[scafNames[i]] = scaffolds[i];
    std::map<std::string, size_t> nextPos;
    std::ifstream f(path);
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream is(line);
        std::string o, type, comp;
        size_t b, e, part;
        is >> o >> b >> e >> part >> type;
        const std::string& s = obj[o];
        if (b != nextPos[o] + 1 || e > s.size()) ++c.badTile;
        nextPos[o] = e;
        const std::string region = s.substr(b - 1, e - b + 1);
        if (type == "N") {
            ++c.n;
            if (region.find_first_not_of("Nn") != std::string::npos) ++c.badSeq;
            continue;
        }
        ++c.w;
        size_t cb, ce;
        char orient;
        is >> comp >> cb >> ce >> orient;
        auto it = records.find(comp);
        if (it == records.end()) { ++c.unknownComp; continue; }
        ++c.compInRecords;
        std::string piece = it->second.substr(cb - 1, ce - cb + 1);
        if (orient == '-') piece = rc(piece);
        if (piece != region) ++c.badSeq;
    }
    for (size_t i = 0; i < scaffolds.size(); ++i) if (nextPos[scafNames[i]] != scaffolds[i].size()) ++c.badTile;
    return c;
}

void testAgpV2() {
    const int k = 99;
    // Scaffold 0: P0 N P1 where P1's tail is an exact copy of scaffold 1's head (a dovetail
    // the trim removes from the shorter partner); scaffold 1: one piece.
    const std::string R = dna(400), P0 = dna(900), X = dna(600), Y = dna(2500);
    const std::string s0 = P0 + std::string(20, 'N') + X + R;
    const std::string s1 = R + Y;
    const std::vector<std::string> sc{s0, s1};
    const std::vector<std::string> names{"S0", "S1"};
    for (int safe = 0; safe < 2; ++safe) {
        auto pieces = ts::splitAtGaps(sc, {"_chr", "_chr"}, {10.0, 10.0});
        std::vector<std::string> seqs;
        for (const auto& p : pieces) seqs.push_back(p.seq);
        ts::TrimConfig cfg;
        cfg.minOverlap = k;
        cfg.k = k;
        cfg.boundarySafe = safe;
        std::vector<size_t> f, b;
        ts::trimTerminalOverlaps(seqs, cfg, f, b, nullptr);
        std::map<std::string, std::string> records;
        for (size_t j = 0; j < pieces.size(); ++j) {
            pieces[j].cutFrontSeq = pieces[j].seq.substr(0, f[j]);
            pieces[j].cutBackSeq = pieces[j].seq.substr(pieces[j].seq.size() - b[j]);
            pieces[j].seq = seqs[j];
            pieces[j].name = "NODE_" + std::to_string(j + 1) + "_length_" + std::to_string(seqs[j].size());
            records[pieces[j].name] = pieces[j].seq;
        }
        const std::string dir = (std::filesystem::temp_directory_path() /
                                 ("v3emit-agp-" + std::to_string(getpid()))).string();
        std::filesystem::create_directories(dir);
        ts::AgpGfaStats st;
        std::string err;
        const bool ok = ts::writeAgpV2(dir + "/s.agp", sc, names, pieces, "paired-ends", st, err);
        const AgpCheck c = checkAgp(dir + "/s.agp", sc, names, records);
        if (getenv("V3EMIT_DEBUG")) {
            std::ifstream in(dir + "/s.agp");
            std::string l;
            while (std::getline(in, l)) std::printf("   agp| %s\n", l.c_str());
            std::printf("   w=%zu n=%zu badSeq=%zu badTile=%zu unknown=%zu partner=%zu unlisted=%zu\n", c.w, c.n,
                        c.badSeq, c.badTile, c.unknownComp, st.wPartner, st.wUnlisted);
        }
        std::filesystem::remove_all(dir);
        check(ok && c.badSeq == 0 && c.badTile == 0 && c.unknownComp == 0 && c.w == 4 && c.n == 1 &&
                  st.wPartner == 1 && st.wUnlisted == 0 && st.nPaired == 1,
              std::string("T26 on") + (safe ? " (+T10)" : "") +
                  ": every W row names a written record and matches it base for base; the trimmed "
                  "overlap points at the record still holding it");
    }
    // P-lines: a piece that is exactly its own gap-free walk gets one, named as the record.
    {
        const int kp = 31;
        const std::string u0 = dna(200), u1 = dna(150);
        const ts::UnitigGraph g = tinyGraph({{u0, 60.0}, {u1, 180.0}}, kp);
        ts::GfaPath path;
        path.oriented = {0 << 1, 1 << 1};
        path.gaps = {0, 5};
        for (int whole = 0; whole < 2; ++whole) {
            const std::string scaf = u0 + "NNNNN" + (whole ? u1 : u1.substr(kp - 1));
            auto pieces = ts::splitAtGaps({scaf}, {"_chr"}, {1.0});
            for (size_t j = 0; j < pieces.size(); ++j) pieces[j].name = "NODE_" + std::to_string(j + 1);
            ts::AgpGfaStats st;
            const auto paths = ts::contigPathsV2({scaf}, {path}, g, pieces, st);
            bool spelled = true;
            for (const auto& p : paths) {
                const size_t j = static_cast<size_t>(std::stoul(p.name.substr(5))) - 1;
                spelled = spelled && ts::spellWalk(p.oriented, g) == pieces[j].seq;
            }
            if (!whole)
                check(paths.size() == 1 && paths[0].name == "NODE_1" && st.pEmitted == 1 && st.pDropped == 1 && spelled,
                      "T26 on: P-line for the piece its walk spells; the k-1-short piece is dropped");
            else
                check(paths.size() == 2 && st.pEmitted == 2 && st.pDropped == 0 && spelled,
                      "T26 on: with the whole-flank render both pieces get exact P-lines");
        }
    }
}

void testStats() {
    ts::AssemblyReport rep;
    ts::computeContigStats({"ACGT", "ACGT", "AAAA", "CCC"}, {"ACGTNACGT", "AAAANNNNNNNNNNNCCC"}, rep);
    check(rep.contigPieces == 4 && rep.scaffoldGaps == 2 && rep.contigTotal == 15 && rep.contigLargest == 4 &&
              rep.contigN50 == 4,
          "T17: stats from the written records; a 1-N run counts as a gap");
}

}  // namespace

int main() {
    testenv::clearTesseractEnv();
    setenv("TESSERACT_FIXES", "0", 1);   // 1.4.0: unset would follow the default umbrella (on)
    testDedup();
    testBoundarySafe();
    testCopyGuard();
    testSplitPost();
    testAgpV2();
    testStats();
    std::printf("test_v3_emit_post: %d checks, %d failures\n", checks, failures);
    return failures ? 1 : 0;
}
