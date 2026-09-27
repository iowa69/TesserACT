// T03 (build_v3 G-resolve): renderChain dropped the first k-1 bases of the first unitig
// after every open scaffold N-gap; the N-run (an estimate of true gap + k-1) stood in for
// them and the GFA P-line (gN + the WHOLE segment) spelled k-1 more than the FASTA.
//
// Fixture: two unitigs A = T[0,3000) and B = T[3000+gap, 6000+gap) of one genome T, FR pairs
// with fixed 300 bp fragments (so the gap estimate is exact), scaffolding on. The scaffold
// and the P-line written by the real ts::writeGfa are checked.
//   release (flags unset):  A + N*(gap+k-1) + B[k-1:]   and the P-line spells k-1 more.
//   TESSERACT_FIX_GAP_FLANK=1 (or TESSERACT_FIXES=1): resolve() still writes the release
//   layout (that is what the gap filler must see: its target k-mer then lies in full-depth
//   sequence) and records each gap; restoreGapFlanks(), run after gap filling, turns every
//   gap still open into
//     gap >= 1, no verified overlap:   A + N*gap + B          (P-line == scaffold)
//     exact overlap L >= 4 (gap = -L): A + 'N' + B[L:]         (P-line spells +L, counted)
//     gap <= 0 and no overlap:         A + 'N' + B  (T20 on: 100 Ns, unknown length)
//   in either orientation, and leaves a gap the filler closed alone.
// Variants: k 31/61, B stored reverse-complemented, node order swapped (the walk then
// starts at B and the piece after the gap is A).
#include "test_v3_resolve_util.h"
#include "gfa.h"

#include <map>
#include <sstream>
#include "test_env.h"   // build_v3: T40, clear the ambient TESSERACT_* first

using namespace v3r;

namespace {

// Spell a GFA P-line with TesserACT's convention (gfa.cpp): `xM` = overlap x, `gN` = g Ns
// then the WHOLE next segment.
std::string spellPLine(const std::string& gfaPath) {
    std::ifstream in(gfaPath);
    std::map<std::string, std::string> seg;
    std::string line, P;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::vector<std::string> f;
        std::stringstream ss(line);
        std::string x;
        while (std::getline(ss, x, '\t')) f.push_back(x);
        if (f[0] == "S") seg[f[1]] = f[2];
        if (f[0] == "P") P = line;
    }
    std::vector<std::string> f, steps, ovs;
    { std::stringstream ss(P); std::string x; while (std::getline(ss, x, '\t')) f.push_back(x); }
    { std::stringstream ss(f[2]); std::string x; while (std::getline(ss, x, ',')) steps.push_back(x); }
    if (f[3] != "*") { std::stringstream ss(f[3]); std::string x; while (std::getline(ss, x, ',')) ovs.push_back(x); }
    std::string out;
    for (size_t i = 0; i < steps.size(); ++i) {
        std::string s = seg[steps[i].substr(0, steps[i].size() - 1)];
        if (steps[i].back() == '-') s = rc(s);
        if (i == 0) { out = s; continue; }
        const std::string& o = ovs[i - 1];
        const long n = std::atol(o.substr(0, o.size() - 1).c_str());
        if (o.back() == 'N') out += std::string(static_cast<size_t>(n), 'N') + s;
        else out += s.substr(static_cast<size_t>(n));
    }
    return out;
}

enum class Mode { Release, Fix, Umbrella, UmbrellaOptOut, FixWithT20 };
const char* modeName(Mode m) {
    switch (m) {
        case Mode::Release: return "flags unset";
        case Mode::Fix: return "TESSERACT_FIX_GAP_FLANK=1";
        case Mode::Umbrella: return "TESSERACT_FIXES=1";
        case Mode::UmbrellaOptOut: return "TESSERACT_FIXES=1 TESSERACT_FIX_GAP_FLANK=0";
        case Mode::FixWithT20: return "TESSERACT_FIX_GAP_FLANK=1 TESSERACT_FIX_GAP_ESTIMATE=1";
    }
    return "?";
}

void gapCase(int k, long gap, bool bRc, bool swapNodes, Mode mode) {
    releaseDefaults();
    if (mode == Mode::Fix) setFlag("TESSERACT_FIX_GAP_FLANK", "1");
    if (mode == Mode::Umbrella) setFlag("TESSERACT_FIXES", "1");
    if (mode == Mode::UmbrellaOptOut) {
        setFlag("TESSERACT_FIXES", "1");
        setFlag("TESSERACT_FIX_GAP_FLANK", "0");
        // everything else the umbrella turns on must not touch this fixture's scaffold
        setFlag("TESSERACT_FIX_GAP_ESTIMATE", "0");
    }
    if (mode == Mode::FixWithT20) {
        setFlag("TESSERACT_FIX_GAP_FLANK", "1");
        setFlag("TESSERACT_FIX_GAP_ESTIMATE", "1");
    }
    const bool fixed = mode == Mode::Fix || mode == Mode::Umbrella || mode == Mode::FixWithT20;
    const bool t20 = mode == Mode::Umbrella || mode == Mode::FixWithT20;

    const std::string G = dna(6400, 11);
    const std::string A = G.substr(0, 3000);
    const std::string B = G.substr(static_cast<size_t>(3000 + gap), 3000);
    const std::string T = G.substr(0, static_cast<size_t>(6000 + gap));
    ts::UnitigGraph g;
    g.setK(k);
    g.nodes.resize(2);
    const int ia = swapNodes ? 1 : 0, ib = swapNodes ? 0 : 1;
    g.nodes[ia].seq = A;
    g.nodes[ib].seq = bRc ? rc(B) : B;
    g.nodes[0].coverage = g.nodes[1].coverage = 40;
    std::vector<Pair> pairs;
    for (size_t p = 0; p + 300 <= T.size(); p += 5) pairs.push_back({T.substr(p, 100), rc(T.substr(p + 200, 100))});
    const ts::SequenceStore reads = store(pairs, "t03");
    StderrCapture cap;
    cap.start();
    ts::PairedResolver r(g, reads, 1, 2, 2.0, 0.10);
    r.setScaffolding(true);
    r.buildSupport();
    std::vector<std::string> seqs;
    std::vector<double> covs;
    r.resolve(seqs, covs);
    const std::vector<std::string> resolverOut = seqs;
    // What the assembler does after gap filling (no filler here: every gap is open).
    ts::restoreGapFlanks(seqs, r.gapFlankRecords());
    const std::string log = cap.stop();
    std::printf("[k=%d gap=%ld bRc=%d swap=%d | %s] contigs=%zu joins=%zu\n", k, gap, bRc, swapNodes,
                modeName(mode), seqs.size(), r.stats().scaffoldJoins);
    size_t pi = SIZE_MAX;
    for (size_t i = 0; i < seqs.size(); ++i) if (seqs[i].find('N') != std::string::npos) pi = i;
    if (pi == SIZE_MAX) { expect(false, "a scaffold with an N-run is emitted"); return; }
    std::string s = seqs[pi];
    const bool flippedOut = s.compare(0, 60, T, 0, 60) != 0;
    if (flippedOut) s = rc(s);

    std::vector<ts::GfaPath> gp(1);
    gp[0].name = "scaf";
    gp[0].oriented = r.paths()[pi].oriented;
    gp[0].gaps = r.paths()[pi].gaps;
    const std::string gfaFile =
        (std::filesystem::temp_directory_path() / ("v3r-t03-" + std::to_string(getpid()) + ".gfa")).string();
    size_t ns = 0, nl = 0;
    std::string err;
    ts::writeGfa(gfaFile, g, gp, ns, nl, err);
    std::string spelled = spellPLine(gfaFile);
    std::filesystem::remove(gfaFile);
    if (flippedOut) spelled = rc(spelled);

    const size_t ov = static_cast<size_t>(k - 1);
    const bool afterIsB = static_cast<int>(r.paths()[pi].oriented[0] >> 1) == ia;
    const auto gl = lines(log, "[gapflank]");
    expect(gl.size() == 1, "exactly one [gapflank] counter line");
    const std::string line = gl.empty() ? "" : gl[0];
    expect(field(line, "enabled") == (fixed ? 1 : 0), "[gapflank] enabled= reports the resolved flag");
    const auto rl = lines(log, "[gapflank-restore]");
    expect(rl.size() == 1 && field(rl[0], "restored") == (fixed ? 1 : 0), "[gapflank-restore] restored the open gap iff fixed");
    // The resolver's own output is the release layout either way.
    {
        std::string ro = resolverOut[pi];
        if (flippedOut) ro = rc(ro);
        const std::string Ns(static_cast<size_t>(gap + k - 1), 'N');
        const std::string expectRel = afterIsB ? T.substr(0, 3000) + Ns + T.substr(static_cast<size_t>(3000 + gap) + ov)
                                               : T.substr(0, 3000 - ov) + Ns + T.substr(static_cast<size_t>(3000 + gap));
        expect(ro == expectRel, "resolve() writes the release layout (what the gap filler sees)");
    }
    if (!fixed) {
        const std::string Ns(static_cast<size_t>(gap + k - 1), 'N');
        const std::string expectRel = afterIsB ? T.substr(0, 3000) + Ns + T.substr(static_cast<size_t>(3000 + gap) + ov)
                                               : T.substr(0, 3000 - ov) + Ns + T.substr(static_cast<size_t>(3000 + gap));
        expect(s == expectRel, "release: k-1 real bases after the Ns are replaced by N");
        expect(spelled.size() == s.size() + ov, "release: the GFA P-line spells k-1 more bases than the scaffold");
        expect(field(line, "gaps") == 0, "release: [gapflank] gaps=0");
        return;
    }
    // The fixed rule, computed independently from the genome.
    size_t L = 0;
    for (size_t l = ov; l > 0; --l) if (A.compare(A.size() - l, l, B, 0, l) == 0) { L = l; break; }
    const bool trim = L >= 4 && gap + static_cast<long>(L) <= 20;   // exact estimate: sigma 0
    long nN = 1;
    if (trim) nN = 1;
    else if (gap >= 1) nN = gap;
    else nN = t20 ? 100 : 1;
    const size_t from = trim ? L : 0;
    const std::string Ns(static_cast<size_t>(nN), 'N');
    const std::string expectFix = afterIsB ? A + Ns + B.substr(from) : A.substr(0, A.size() - from) + Ns + B;
    std::printf("  exact suffix/prefix overlap L=%zu trim=%d -> N-run %ld\n", L, trim ? 1 : 0, nN);
    expect(s == expectFix, "fix: the unitig after the Ns is kept whole (less a verified overlap); N-run = estimate");
    if (gap >= 1 && !trim)
        expect(s == T.substr(0, 3000) + std::string(static_cast<size_t>(gap), 'N') + T.substr(static_cast<size_t>(3000 + gap)),
               "fix: every non-N base equals the genome at its position and the N-run equals the true gap");
    if (trim) {
        expect(s.find(std::string(2, 'N')) == std::string::npos, "fix: a true overlap is dropped once and marked by one N");
        expect(spelled.size() == s.size() + L, "fix: overlap case, the P-line spells +L (GFA1 cannot express gap+trim)");
        expect(field(line, "pLineOverspell") == static_cast<long long>(L), "fix: [gapflank] pLineOverspell=L");
    } else {
        expect(spelled == s, "fix: the GFA P-line spells exactly the scaffold");
    }
    expect(field(line, "gaps") == 1 && field(line, "restorableBases") == static_cast<long long>(ov - from),
           "fix: [gapflank] gaps=1 restorableBases=k-1-L");
    expect(r.stats().gapBases == static_cast<size_t>(gap + k - 1), "fix: stats gapBases = Ns the resolver wrote (release layout)");
    // The same records applied to the reverse complement of the scaffold.
    {
        std::vector<std::string> rcv{rc(resolverOut[pi])};
        ts::restoreGapFlanks(rcv, r.gapFlankRecords());
        expect(rcv[0] == rc(seqs[pi]), "fix: restoration works on the reverse-complemented scaffold");
    }
    // A gap the filler closed (the true locus, no N) is left alone.
    {
        std::vector<std::string> closed{flippedOut ? rc(T) : T};
        StderrCapture c2;
        c2.start();
        ts::restoreGapFlanks(closed, r.gapFlankRecords());
        const std::string l2 = c2.stop();
        const auto rr = lines(l2, "[gapflank-restore]");
        expect(closed[0] == (flippedOut ? rc(T) : T) && !rr.empty() && field(rr[0], "restored") == 0 &&
               field(rr[0], "unmatchedRecords") == 1, "fix: a closed gap is left alone (unmatchedRecords=1)");
    }
}
}  // namespace

int main() {
    testenv::clearTesseractEnv();
    struct C { int k; long gap; bool rcB, swap; };
    const C cases[] = {C{31, 100, false, false}, C{31, 0, false, false}, C{31, -20, false, false},
                       C{61, 50, false, false}, C{31, 100, true, false}, C{31, 100, false, true},
                       C{31, -20, true, true}, C{61, 7, true, false}, C{61, -40, false, false}};
    for (const C& c : cases)
        for (Mode m : {Mode::Release, Mode::Fix})
            gapCase(c.k, c.gap, c.rcB, c.swap, m);
    // Umbrella resolution: TESSERACT_FIXES=1 turns the fix on; TESSERACT_FIX_GAP_FLANK=0 opts out.
    gapCase(31, 100, false, false, Mode::Umbrella);
    gapCase(31, 100, false, false, Mode::UmbrellaOptOut);
    // With T20: an estimate <= 0 without a verified overlap is an unknown-length gap.
    gapCase(31, 0, false, false, Mode::FixWithT20);
    gapCase(61, 50, false, false, Mode::FixWithT20);
    releaseDefaults();
    return finish("test_v3_resolve_gapflank (T03)");
}
