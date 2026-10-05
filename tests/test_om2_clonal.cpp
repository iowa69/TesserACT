// Organism Model 2.0, C4 (CLONAL): clonal-anchored closure.
//
//   rule        decideLink on the cases the design names: a positional repeat (rRNA-like) sized and
//               filled from the isolate's walk; a non-positional element the nearest relatives carry
//               (filled) and one they do not (a NOVEL insertion: labelled gap sized by the isolate's
//               walk, no fill); a clone MISSING an element its relatives carry (the isolate's short
//               walk wins over the relatives' gap); a rearranged clone (the graph's first single-copy
//               continuation names another end: refused); read-pair contradiction; too few agreeing
//               relatives; ambiguous sizes; the KC2 switch (non-positional elements as gaps only)
//   panels      a synthetic clonal panel (an ancestor, a clonal lineage, far lineages) built into a
//               model with layout tracks: the nearest relatives are ranked by true divergence; a
//               junction cut out of the isolate is placed on every relative at its true gap, on
//               either strand and across a relative's origin; a rearranged clone is contradicted;
//               a lineage IS the far lineages lack is non-positional across the panel, a conserved
//               operon locus positional; a novel insertion reads as a size difference
//   confidence  unsized joins are never confident; contradicting relatives and pairs lower it; a
//               dev table (n >= 20) replaces the prior
//   off         with TESSERACT_OM2_CLONAL unset the options are off and the counter line reads
//               enabled=0 with zeros
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <random>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <vector>

#include "graph.h"
#include "om2_clonal.h"
#include "organism.h"
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
char other(char c) { return c == 'A' ? 'C' : c == 'C' ? 'G' : c == 'G' ? 'T' : 'A'; }
std::string mutate(const std::string& s, double rate) {
    std::string t = s;
    std::uniform_real_distribution<double> u(0, 1);
    for (char& c : t) if (u(rng) < rate) c = other(c);
    return t;
}
std::string rc(const std::string& s) { return reverseComplement(s); }

// ---------------------------------------------------------------------------------- rule

LinkEvidence base() {
    LinkEvidence e;
    e.k = 5; e.placed = 5; e.agree = 5;
    e.wAgree = 5; e.wContra = 0;   // five equally close relatives, all for this partner
    e.nrGap = 5080; e.nrSpread = 40;
    e.panelPlaced = 100; e.panelAgree = 96;
    e.anchoredA = e.anchoredB = e.edgeA = e.edgeB = true;
    e.reachable = true; e.exhaustive = true;
    e.repeatA = false; e.repeatB = false; e.element = true;
    e.walkLens = {5100};
    return e;
}

void testRule() {
    ClonalOptions o;   // defaults: k 5, k_min 3, fill non-positional elements when agreed
    // positional repeat (an rRNA operon: the gap recurs across the panel), filled from the walk
    {
        const LinkDecision d = decideLink(base(), o);
        check(d.join && d.n == 5100 && d.walkMatch && d.cls == kPositional && !d.noFill && d.basis == "graph+nr",
              "positional operon junction: sized by the matching walk and filled");
    }
    // a unique-flanked gap the whole panel shares: layout class
    {
        LinkEvidence e = base();
        e.element = false; e.walkLens = {5070};
        const LinkDecision d = decideLink(e, o);
        check(d.join && d.cls == kLayout && !d.noFill, "conserved unique gap: layout class");
    }
    // non-positional element (IS) the nearest relatives carry at the same flanks and size
    LinkEvidence is = base();
    is.panelAgree = 35; is.nrGap = 1320; is.walkLens = {1332};
    {
        const LinkDecision d = decideLink(is, o);
        check(d.join && d.cls == kNonPositional && d.walkMatch && !d.noFill && d.n == 1332,
              "IS the relatives carry: filled");
    }
    // ... but without the element->flank edge it stays a labelled gap
    {
        LinkEvidence e = is;
        e.edgeB = false;
        const LinkDecision d = decideLink(e, o);
        check(d.join && d.noFill, "IS without both flank edges: labelled gap");
    }
    // ... and with only 2 of 5 relatives agreeing (k_min 3) there is no join at all
    {
        LinkEvidence e = is;
        e.agree = 2; e.wAgree = 2; e.wContra = 3;
        check(!decideLink(e, o).join && decideLink(e, o).why == "nr_agree", "2 of 5 equally close relatives: refused");
        ClonalOptions oc = o;
        oc.weighted = false;
        check(!decideLink(e, oc).join, "count rule: 2/5 refused");
        e.agree = 3; e.wAgree = 3; e.wContra = 2;
        check(!decideLink(e, oc).join && !decideLink(e, o).join, "3 of 5 (share 0.6 < 0.66): refused under both rules");
        e.agree = 4; e.wAgree = 4; e.wContra = 1;
        check(decideLink(e, oc).join && decideLink(e, o).join, "4 of 5: joined under both rules");
    }
    // distance-weighted votes: the nearest clonemate (weight 1) carries a join the lineage-level relatives
    // (weight ~0.05 each) either do not see or contradict; equally close contradicting relatives do not let it
    {
        LinkEvidence e = is;
        e.agree = 1; e.placed = 5; e.wAgree = 1.0; e.wContra = 4 * 0.05;
        check(decideLink(e, o).join, "weighted: the nearest clonemate outvotes 4 distant relatives");
        ClonalOptions oc = o;
        oc.weighted = false;
        check(!decideLink(e, oc).join, "count rule: 1/5 refused");
        e.wAgree = 1.0; e.wContra = 1.9;
        check(!decideLink(e, o).join, "weighted: two relatives as close as the nearest disagree: refused");
        e.wAgree = 0.6; e.wContra = 0;
        check(!decideLink(e, o).join, "weighted: support below one clonemate: refused");
    }
    // KC2: non-positional elements as labelled gaps only
    {
        ClonalOptions o2 = o;
        o2.nonposFill = false;
        const LinkDecision d = decideLink(is, o2);
        check(d.join && d.noFill, "NONPOS=gap: an agreed IS is a labelled gap");
    }
    // NOVEL insertion: the relatives are adjacent (gap ~0), the isolate's only walk carries an IS
    {
        LinkEvidence e = is;
        e.nrGap = 4; e.nrSpread = 2; e.walkLens = {1335};
        const LinkDecision d = decideLink(e, o);
        check(d.join && d.cls == kNonPositional && !d.walkMatch && d.noFill && d.n == 1335 && d.basis == "graph",
              "novel IS: labelled gap sized by the isolate's walk, no bases");
    }
    // clone MISSING an IS its relatives carry: relatives 1320 apart, the isolate's walk is 6 bp of unique
    // sequence (no element) -> sized by the isolate, plain graph walk allowed
    {
        LinkEvidence e = is;
        e.walkLens = {6}; e.element = false;
        const LinkDecision d = decideLink(e, o);
        check(d.join && d.n == 6 && !d.walkMatch && !d.noFill && d.cls == kNonPositional,
              "missing IS: the isolate's own short walk sizes the join");
    }
    // rearranged clone: no walk, and the first single-copy continuation names another record end
    {
        LinkEvidence e = base();
        e.reachable = false; e.walkLens.clear(); e.contraGraph = true;
        const LinkDecision d = decideLink(e, o);
        check(!d.join && d.why == "graph_contra", "rearranged clone (graph contradiction): refused");
    }
    // read pairs name a third record end
    {
        LinkEvidence e = base();
        e.pairContra = 4;
        check(!decideLink(e, o).join && decideLink(e, o).why == "pair_contra", "pair contradiction: refused");
    }
    // dead ends (no walk): the relatives size it when they agree, else no join
    {
        LinkEvidence e = base();
        e.reachable = false; e.walkLens.clear();
        const LinkDecision d = decideLink(e, o);
        check(d.join && d.n == 5080 && d.basis == "nr", "dead-end gap: sized by the relatives");
        e.nrSpread = 4000;
        check(!decideLink(e, o).join && decideLink(e, o).why == "nr_spread", "relatives disagree on the size: refused");
    }
    // several walk classes, none at the relatives' size
    {
        LinkEvidence e = base();
        e.walkLens = {1200, 9000};
        check(!decideLink(e, o).join && decideLink(e, o).why == "size_ambiguous", "ambiguous walk sizes: refused");
    }
    // a plasmid join needs the isolate's own walk at the relatives' plasmid size
    {
        LinkEvidence e = base();
        e.plasmid = true; e.k = 2; e.placed = 3; e.agree = 3; e.wAgree = 3;
        const LinkDecision d = decideLink(e, o);
        check(d.join && d.cls == kPlasmid && d.walkMatch, "plasmid join with a matching walk");
        e.reachable = false; e.walkLens.clear();
        check(!decideLink(e, o).join && decideLink(e, o).why == "plasmid_needs_walk",
              "plasmid join sized by the relatives alone: refused");
    }
    // a large overlap on the relatives is not a join
    {
        LinkEvidence e = base();
        e.reachable = false; e.walkLens.clear(); e.nrGap = -5000; e.nrSpread = 10;
        check(!decideLink(e, o).join, "records overlapping 5 kb on the relatives: refused");
    }
}

// ---------------------------------------------------------------------------------- round 2 rules

void testRound2() {
    ClonalOptions o;
    // exact end overlap: found at the relatives' size, refused when absent or far from it
    {
        const std::string core = dna(300), a = dna(2000), b = dna(2000);
        const std::string L = a + core, R = core + b;
        check(exactEndOverlap(L, R, 300, 25, 1000) == 300, "exact overlap of 300 found at 300");
        check(exactEndOverlap(L, R, 280, 25, 1000) == 300, "exact overlap found within 10% of the relatives' 280");
        check(exactEndOverlap(L, R, 150, 25, 1000) == -1, "no exact overlap near 150");
        check(exactEndOverlap(a, b, 300, 25, 1000) == -1, "unrelated ends: none");
        std::string R2 = R;
        R2[150] = other(R2[150]);
        check(exactEndOverlap(L, R2, 300, 25, 1000) == -1, "one mismatch inside the overlap: none");
        std::string Ln = L;
        Ln[Ln.size() - 10] = 'N';
        std::string Rn = R;
        Rn[290] = 'N';
        check(exactEndOverlap(Ln, Rn, 300, 25, 1000) == -1, "N inside the overlap: none");
    }
    // decideLink: an exact overlap joins even when the walks disagree (size_ambiguous before round 2)
    {
        LinkEvidence e = base();
        e.nrGap = -351; e.nrSpread = 0; e.walkLens = {23552}; e.exhaustive = false; e.overlapExact = 351;
        e.wPlaced = 5;
        LinkDecision d = decideLink(e, o);
        check(d.join && d.overlap && d.n == -351 && d.basis == "overlap_exact", "exact overlap: joined at -351");
        ClonalOptions off = o;
        off.overlapMerge = false;
        check(!decideLink(e, off).join && decideLink(e, off).why == "size_ambiguous", "OVERLAP=0: round-1 refusal");
        // a graph contradiction is overruled only when every placed relative agrees
        e.reachable = false; e.walkLens.clear(); e.contraGraph = true;
        check(decideLink(e, o).join, "exact overlap + all relatives: graph contradiction overruled");
        e.wAgree = 3.5; e.wContra = 1.5;
        check(!decideLink(e, o).join && decideLink(e, o).why == "graph_contra", "exact overlap, relatives split: graph contradiction stands");
        e = base(); e.nrGap = -300; e.overlapExact = 300; e.wPlaced = 5; e.pairContra = 5;
        check(!decideLink(e, o).join && decideLink(e, o).why == "pair_contra", "exact overlap never overrules pairs to a third end");
        e = base(); e.nrGap = -300; e.overlapExact = 20; e.wPlaced = 5; e.walkLens = {8000}; e.exhaustive = false;
        check(!decideLink(e, o).join, "overlap below 25 bp: not an overlap join");
    }
    // bracket: one walk class within 2 kb of the relatives: midpoint, only with the flag; never filled
    {
        LinkEvidence e = base();
        e.nrGap = -390; e.nrSpread = 0; e.walkLens = {889}; e.exhaustive = false;
        LinkDecision d = decideLink(e, o);
        check(!d.join && d.why == "size_ambiguous" && d.bracketN == 249, "bracket computed (249) but off by default");
        ClonalOptions ob = o;
        ob.bracket = true;
        d = decideLink(e, ob);
        check(d.join && d.bracket && d.n == 249 && d.noFill && d.basis == "bracket", "BRACKET=1: midpoint 249, no fill");
        e.walkLens = {6979}; e.nrGap = 10676;
        d = decideLink(e, ob);
        check(!d.join && d.bracketN == INT32_MIN, "walk and relatives 3.7 kb apart: no bracket");
        e.walkLens = {1200, 9000}; e.nrGap = 2000;
        check(!decideLink(e, ob).join, "two walk classes: no bracket");
        e = base(); e.nrGap = 3134; e.walkLens = {2069}; e.exhaustive = false; e.plasmid = true; e.k = 2;
        check(!decideLink(e, ob).join, "plasmid: no bracket");
    }
    // flags
    {
        testenv::clearTesseractEnv();
        ClonalOptions d = ClonalOptions::fromEnv();
        check(d.overlapMerge && d.vouch && !d.bracket && !d.anyFlagSet, "round-2 defaults: overlap 1, vouch 1, bracket 0");
        setenv("TESSERACT_OM2_CLONAL_OVERLAP", "0", 1);
        setenv("TESSERACT_OM2_CLONAL_VOUCH", "0", 1);
        setenv("TESSERACT_OM2_CLONAL_BRACKET", "1", 1);
        d = ClonalOptions::fromEnv();
        check(!d.overlapMerge && !d.vouch && d.bracket && d.anyFlagSet, "round-2 flags read");
        testenv::clearTesseractEnv();
    }
}

// ---------------------------------------------------------------------------------- round 3 rules

void testRound3() {
    // isolatedShortestWalk: the shortest class, isolated by >= gap from the next, walked once, <= maxLen
    check(isolatedShortestWalk({572}, {1}, true, 5000, 20000) == 572, "single exhaustive class: 572");
    check(isolatedShortestWalk({572}, {1}, false, 5000, 20000) == INT32_MIN, "single class, enumeration cut: none");
    check(isolatedShortestWalk({572, 28273}, {1, 12}, false, 5000, 20000) == 572, "572 then 28,273 (loops): 572");
    check(isolatedShortestWalk({28273, 572}, {12, 1}, false, 5000, 20000) == 572, "order of classes does not matter");
    check(isolatedShortestWalk({4028, 13805, 15006}, {1, 1, 1}, false, 5000, 20000) == 4028, "4,028 then 13,805: 4,028");
    check(isolatedShortestWalk({465, 3519, 6573}, {1, 1, 1}, false, 5000, 20000) == INT32_MIN,
          "tandem period 3,054 (< 5 kb): ambiguous");
    check(isolatedShortestWalk({-159, 8, 175}, {1, 1, 1}, false, 5000, 20000) == INT32_MIN, "short tandem: ambiguous");
    check(isolatedShortestWalk({1003}, {32}, true, 5000, 20000) == INT32_MIN, "shortest class walked 32 times: none");
    check(isolatedShortestWalk({26635}, {1}, true, 5000, 20000) == INT32_MIN, "shortest class above 20 kb: none");
    check(isolatedShortestWalk({}, {}, true, 5000, 20000) == INT32_MIN, "no walk: none");
    check(isolatedShortestWalk({572, 28273}, {}, false, 5000, 20000) == 572, "unknown counts: lengths decide");
    // decideLink with WALKSIZE: relatives agree on the order, their gap is no walk, the shortest walk is isolated
    ClonalOptions o;
    {
        LinkEvidence e = base();
        e.nrGap = -200; e.walkLens = {1000, 25968}; e.walkCounts = {1, 1}; e.exhaustive = false; e.element = false;
        check(!decideLink(e, o).join && decideLink(e, o).why == "size_ambiguous", "WALKSIZE off: size_ambiguous");
        ClonalOptions ow = o;
        ow.walksize = true;
        LinkDecision d = decideLink(e, ow);
        check(d.join && d.c1 && d.n == 1000 && d.basis == "graph_c1" && d.cls == kNonPositional && !d.walkMatch,
              "WALKSIZE: joined at the isolated shortest walk 1,000");
        e.walkLens = {465, 3519, 6573}; e.walkCounts = {1, 1, 1};
        check(!decideLink(e, ow).join, "WALKSIZE: a tandem-periodic walk set stays size_ambiguous");
        e.walkLens = {1000, 25968}; e.plasmid = true; e.k = 2;
        check(!decideLink(e, ow).join, "WALKSIZE never sizes a plasmid link");
        e = base(); e.nrGap = -200; e.walkLens = {1000, 25968}; e.walkCounts = {1, 1}; e.exhaustive = false;
        e.pairContra = 3;
        check(!decideLink(e, ow).join && decideLink(e, ow).why == "pair_contra", "WALKSIZE never overrules pairs");
        e.pairContra = 0; e.contraGraph = true;
        check(!decideLink(e, ow).join, "WALKSIZE never overrules a graph contradiction");
    }
    // KC2: a non-positional ELEMENT is never filled, even when the relatives carry it at the walk size
    {
        LinkEvidence e = base();
        e.panelAgree = 40;   // the relatives' gap recurs at 40% of the panel: non-positional
        e.element = true; e.walkLens = {5100}; e.nrGap = 5080;
        LinkDecision d = decideLink(e, o);
        check(d.join && d.cls == kNonPositional && !d.noFill, "KC2 off: a carried element is filled");
        ClonalOptions ok = o;
        ok.kc2 = true;
        d = decideLink(e, ok);
        check(d.join && d.cls == kNonPositional && d.noFill, "KC2: the element is a labelled gap (noFill)");
        e.element = false;
        check(!decideLink(e, ok).noFill, "KC2: a repeat-free non-positional walk may be filled");
        e.element = true; e.panelAgree = 96;
        d = decideLink(e, ok);
        check(d.cls == kPositional && !d.noFill, "KC2: a positional locus (96% of the panel) is filled");
    }
    // flags
    {
        testenv::clearTesseractEnv();
        ClonalOptions d = ClonalOptions::fromEnv();
        check(!d.walksize && !d.kc2 && !d.confRule && d.discord == 2 && !d.anyFlagSet, "round-3 defaults: all off");
        setenv("TESSERACT_OM2_CLONAL_WALKSIZE", "1", 1);
        setenv("TESSERACT_OM2_CLONAL_KC2", "1", 1);
        setenv("TESSERACT_OM2_CLONAL_CONFRULE", "1", 1);
        setenv("TESSERACT_OM2_CLONAL_DISCORD", "3", 1);
        d = ClonalOptions::fromEnv();
        check(d.walksize && d.kc2 && d.confRule && d.discord == 3 && d.anyFlagSet, "round-3 flags read");
        testenv::clearTesseractEnv();
        d = ClonalOptions::fromEnv();
        check(!d.unverified && !d.indel && d.c1GapLayout == 10000 && d.indelMax == 250000, "round-3b defaults: off");
        setenv("TESSERACT_OM2_CLONAL_UNVERIFIED", "1", 1);
        setenv("TESSERACT_OM2_CLONAL_INDEL", "1", 1);
        d = ClonalOptions::fromEnv();
        check(d.unverified && d.indel && d.anyFlagSet, "round-3b flags read");
        testenv::clearTesseractEnv();
    }
    // counter line carries the round-3 fields (zeros when off)
    {
        ClonalStats st;
        const std::string l = formatClonalCounters(st);
        check(l.find("walksize_resized=0") != std::string::npos && l.find("kc2_nofill=0") != std::string::npos &&
              l.find("conf_capped=0") != std::string::npos && l.find("discord=0") != std::string::npos &&
              l.find("unverified_broken=0") != std::string::npos && l.find("cap_walk_c1=0") != std::string::npos,
              "counter line: round-3 / 3b fields");
        // the layout rule isolates the shortest class by 10 kb: 1,940 then 10,835 (a true 4,217 in dev Kp) is ambiguous
        check(isolatedShortestWalk({1940, 10835, 11045}, {1, 1, 2}, false, 10000, 20000) == INT32_MIN,
              "layout rule: 1,940 / 10,835 not isolated at 10 kb");
        check(isolatedShortestWalk({572, 28273}, {1, 12}, false, 10000, 20000) == 572, "layout rule: 572 / 28,273 isolated");
    }
}

// ---------------------------------------------------------------------------------- panels

constexpr uint64_t kDenom = 32;   // dense sampling: small synthetic genomes still carry many markers

void addGenome(OrganismModel& m, const std::string& name, const std::string& chr) {
    std::vector<std::pair<uint64_t, std::pair<uint32_t, int>>> hits;
    forEachMarkerKmer(chr, [&](uint64_t km, uint32_t pos, int o) { hits.push_back({km, {pos, o}}); }, kDenom);
    std::unordered_map<uint64_t, int> occ;
    for (const auto& h : hits) ++occ[h.first];
    LayoutTrack t;
    t.name = name;
    for (const auto& h : hits) {
        if (occ[h.first] != 1) continue;
        const uint32_t id = m.internMarker(h.first);
        m.noteMarkerGenome(id, Replicon::Chromosome);
        t.oriented.push_back((id << 1) | static_cast<uint32_t>(h.second.second));
        t.pos.push_back(h.second.first);
    }
    m.noteGenome(Replicon::Chromosome);
    m.addTrack(std::move(t));
}

// Synthetic panel. Ancestor chromosome with 3 copies of an operon-like repeat R at fixed loci (every
// genome has them: positional) and IS copies. The isolate's lineage (8 clonal relatives) carries an
// extra IS at `lineageIS`; 12 far genomes do not (non-positional across the panel).
struct Panel {
    std::string anc, R, IS;
    size_t op[3] = {60000, 150000, 240000};
    size_t lineageIS = 110000;
    size_t novelSite = 190000;
    std::vector<std::string> names;
    std::vector<std::string> chrs;
};

std::string withInsert(const std::string& s, size_t at, const std::string& ins) {
    return s.substr(0, at) + ins + s.substr(at);
}

Panel makePanel(OrganismModel& m) {
    Panel p;
    p.anc = dna(300000);
    p.R = dna(5000);
    p.IS = dna(1300);
    // operons at the same loci everywhere (inserted right to left keeps offsets valid)
    std::string ancOps = p.anc;
    for (int i = 2; i >= 0; --i) ancOps = withInsert(ancOps, p.op[i], mutate(p.R, 0.002));
    // lineage ancestor: the lineage IS (coordinates before the operons are inserted refer to p.anc)
    const std::string lineage = mutate(ancOps, 0.0005);
    // positions shift by the operons inserted before them
    auto shifted = [&](size_t pos) {
        size_t s = pos;
        for (size_t o : p.op) if (o < pos) s += 5000;
        return s;
    };
    const std::string lineageWithIS = withInsert(lineage, shifted(p.lineageIS), p.IS);
    m.beginBuild("synthetica", kMarkerK);
    m.setMarkerDenom(kDenom);
    for (int i = 0; i < 8; ++i) {
        p.names.push_back("CLONE" + std::to_string(i + 1));
        p.chrs.push_back(mutate(lineageWithIS, 0.0002 * (i + 1)));   // increasing distance
    }
    for (int i = 0; i < 12; ++i) {
        p.names.push_back("FAR" + std::to_string(i + 1));
        p.chrs.push_back(mutate(ancOps, 0.012));
    }
    for (size_t i = 0; i < p.chrs.size(); ++i) addGenome(m, p.names[i], p.chrs[i]);
    // the isolate: another member of the lineage, closest to CLONE1
    p.anc = lineageWithIS;   // keep the lineage genome for the isolate below
    return p;
}

size_t shift(const Panel& p, size_t pos, bool withLineageIS) {
    size_t s = pos;
    for (size_t o : p.op) if (o < pos) s += 5000;
    if (withLineageIS && p.lineageIS < pos) s += 1300;
    return s;
}

void testRound3c() {
    // STRICT: the emission rule of a model-proposed junction at its final size (dev round-3b cases, clonal_round3.md)
    ClonalOptions o;
    o.enabled = true;
    o.strict = true;
    auto S = [](int64_t n, std::vector<int32_t> w, bool ex, uint32_t p, uint32_t a, uint32_t c, double wa, double wc) {
        StrictInput s;
        s.n = n; s.walkLens = std::move(w); s.exhaustive = ex; s.placed = p; s.agree = a; s.contra = c;
        s.wAgree = wa; s.wContra = wc;
        return s;
    };
    auto keep = [&](const StrictInput& s) { return strictVerdict(s, o) == nullptr; };
    auto why = [&](const StrictInput& s) { const char* w = strictVerdict(s, o); return std::string(w ? w : "keep"); };
    // overlaps and butts always pass
    StrictInput ov = S(-60, {}, false, 0, 0, 0, 0, 0); ov.overlap = true;
    check(keep(ov), "exact overlap: kept");
    check(keep(S(1, {}, false, 0, 0, 0, 0, 0)), "1-N butt: kept");
    // the isolate's own graph sizes the gap (Kp chr98 joins: GCF001022035v1 1,004 / 1,563; GCF030272075v1 715 / 900)
    check(keep(S(1004, {1004, 25972}, true, 5, 0, 1, 0.0, 0.241)), "isolated shortest walk 1,004, relatives other size: kept");
    check(keep(S(1563, {1563, 11482, 26598, 27877}, false, 0, 0, 0, 0, 0)), "shortest of several, no contradiction: kept");
    check(keep(S(900, {900}, false, 5, 0, 0, 0, 0)), "single class 900, relatives other size: kept");
    check(keep(S(1049, {1049}, true, 3, 1, 1, 0.542, 1.0)), "unique exhaustive walk, relatives contradict: kept");
    // shortest of several walk classes that the weighted relatives contradict: broken
    check(why(S(870, {946, 10377, 19808, 29235}, false, 3, 1, 2, 0.4, 1.9)) == "strict_walk_contradicted",
          "shortest-not-isolated, weighted relatives contradict: broken");
    check(keep(S(1398, {1398}, false, 4, 2, 1, 1.436, 0.852)), "single class, contradiction below agreement: kept");
    // a longer walk class the relatives do not all read: broken
    check(why(S(5068, {687, 5068, 9449, 13830}, false, 5, 2, 0, 1.2, 0)) == "strict_not_shortest", "not the shortest class: broken");
    // no walk class at the size, relatives not unanimous: broken
    check(why(S(4061, {1940, 10835, 11045}, false, 5, 1, 3, 0.5, 1.5)) == "strict_no_walk_class", "no walk class: broken");
    // walk-less: kept only when >= 2 relatives are placed and all read the size
    check(keep(S(2232, {}, false, 5, 5, 0, 4.0, 0)), "walk-less, 5 of 5 relatives at the size: kept");
    check(why(S(1118, {}, false, 5, 4, 0, 3.3, 0)) == "strict_walkless", "walk-less, 4 of 5 relatives: broken");
    check(why(S(1915, {}, false, 1, 1, 0, 1.0, 0)) == "strict_walkless", "walk-less, one relative: broken");
    // unanimous relatives at a walk class that is not the shortest: kept
    check(keep(S(5208, {4979, 5208, 5303}, false, 5, 5, 0, 2.3, 0)), "unanimous relatives at a walk class: kept");
    // unanimous relatives, no walk class at the size: broken
    check(why(S(3000, {900}, false, 3, 3, 0, 2.0, 0)) == "strict_no_walk_class", "unanimous relatives, no walk at the size: broken");
    // periodic walk set (tandem array): the shortest class needs near relatives agreeing
    check(strictPeriodic({465, 3519, 6573, 9627}), "465 / 3,519 / 6,573: periodic");
    check(strictPeriodic({-159, 8, 175, 342}), "-159 / 8 / 175: periodic");
    check(!strictPeriodic({1563, 11482, 26598}), "1,563 / 11,482 / 26,598: not periodic");
    check(!strictPeriodic({1004, 25972}), "two classes: not periodic");
    check(why(S(465, {465, 3519, 6573, 9627}, false, 5, 0, 0, 0, 0)) == "strict_periodic", "periodic, no relative agrees: broken");
    check(why(S(120, {-155, -37, 81, 199}, false, 5, 2, 3, 1.04, 0.13)) == "strict_periodic",
          "periodic 118-bp unit, weighted agreement 1.04 (GCF039653445v1): broken");
    check(keep(S(518, {518, 3183, 5848, 8513}, false, 5, 4, 1, 2.65, 0.245)), "periodic, near relatives agree (GCF904865805v1): kept");
    // NR0 (A/B): no relative places the ends -> broken even at the isolate's walk
    check(keep(S(1292, {1292}, true, 0, 0, 0, 0, 0)), "no relative, unique walk, NR0 off: kept");
    o.strictNr0 = true;
    check(why(S(1292, {1292}, true, 0, 0, 0, 0, 0)) == "strict_no_relative", "no relative, NR0 on: broken");
    check(keep(S(1004, {1004, 25972}, true, 5, 0, 1, 0.0, 0.241)), "relatives placed, NR0 on: kept");
    // flags: STRICT is read from the environment and is OFF by default
    testenv::clearTesseractEnv();
    check(!ClonalOptions::fromEnv().strict && !ClonalOptions::fromEnv().strictNr0, "STRICT / STRICT_NR0 off by default");
    setenv("TESSERACT_OM2_CLONAL_STRICT", "1", 1);
    setenv("TESSERACT_OM2_CLONAL_STRICT_NR0", "1", 1);
    check(ClonalOptions::fromEnv().strict && ClonalOptions::fromEnv().strictNr0 && ClonalOptions::fromEnv().anyFlagSet,
          "STRICT flags read");
    testenv::clearTesseractEnv();
    ClonalStats cs;
    cs.enabled = true; cs.strict = true; cs.strictChecked = 7; cs.strictBroken = 2; cs.strictRefused = 1;
    const std::string line = formatClonalCounters(cs);
    check(line.find("strict=1 strict_nr0=0 strict_checked=7 strict_broken=2 strict_refused=1 strict_periodic=0 "
                    "strict_no_relative=0") != std::string::npos, "counter line carries the STRICT counters");
    check(formatClonalCounters(ClonalStats{}).find("strict=0 strict_nr0=0 strict_checked=0") != std::string::npos,
          "STRICT counters are zero when off");
}

void testPanels() {
    OrganismModel m;
    Panel p = makePanel(m);
    const std::string iso = mutate(p.chrs[0], 0.0001);   // a clonemate of CLONE1
    ClonalOptions o;
    o.enabled = true;
    ClonalEngine eng(m, o);
    eng.selectNearestFromSeqs({iso});
    const NearestSet& ns = eng.nearest();
    check(ns.ok && ns.ranked.size() >= 5, "nearest relatives found");
    check(ns.ranked[0].name == "CLONE1", "the nearest relative is the closest clone (" + ns.ranked[0].name + ")");
    bool ordered = true;
    for (size_t i = 0; i < 5; ++i) if (ns.ranked[i].name.compare(0, 5, "CLONE") != 0) ordered = false;
    check(ordered, "the 5 nearest are the clonal lineage");
    check(ns.ranked[0].dist < ns.ranked[4].dist && ns.ranked[4].dist < ns.ranked[8].dist,
          "distances increase with divergence");
    check(ns.ranked[0].dist < 0.001, "clonemate distance below 0.001");
    check(eng.relatives() == 5, "k = 5 relatives indexed");
    eng.setAssembly({iso});

    // a junction cut out of the isolate at a unique locus: every relative places it at the true gap
    auto cut = [&](size_t a, size_t gap) {
        return std::make_pair(iso.substr(a - 20000, 20000), iso.substr(a + gap, 20000));
    };
    {
        const auto lr = cut(90000, 700);
        const TrackSupport ts = eng.relativeSupport(eng.query(lr.first, lr.second, 700, true));
        check(ts.placed == 5 && ts.agree == 5, "unique junction: 5/5 relatives agree");
        bool near = true;
        for (int64_t g : ts.agreeGaps) if (std::llabs(g - 700) > 20) near = false;
        check(near, "relatives' gap equals the true gap (+-20)");
        // the same junction read on the other strand
        const TrackSupport tr = eng.relativeSupport(eng.query(rc(lr.second), rc(lr.first), 700, true));
        check(tr.agree == 5, "reverse strand: 5/5 agree");
        // asserted wrongly by 5 kb: a size difference, not a contradiction
        const TrackSupport tw = eng.relativeSupport(eng.query(lr.first, lr.second, 5700, true));
        check(tw.gapDiff == 5 && tw.agree == 0 && tw.contra == 0, "wrong size: gap difference");
    }
    // across a relative's origin (the isolate's record ends and restarts there)
    {
        const TrackSupport ts = eng.relativeSupport(eng.query(iso.substr(iso.size() - 20000), iso.substr(0, 20000), 0, true));
        check(ts.agree == 5, "junction across the chromosome origin: agree");
    }
    // REARRANGED clone: an inversion of the segment between two unique points makes the new junction
    // (left of the segment -> reverse of its far end) contradict every relative
    {
        const size_t a = 70000, b = 100000;
        const std::string left = iso.substr(a - 20000, 20000);
        const std::string right = rc(iso.substr(b - 20000, 20000));   // the inverted segment starts with rc of its end
        const TrackSupport ts = eng.relativeSupport(eng.query(left, right, 0, true));
        check(ts.contra >= 3 && ts.agree == 0, "rearranged clone: relatives contradict the inversion junction");
    }
    // the lineage IS: junction between its two flanks, asserted with the element's size
    const size_t isAt = shift(p, p.lineageIS, false);   // left flank ends here in the lineage genome
    {
        const std::string left = iso.substr(isAt - 20000, 20000);
        const std::string right = iso.substr(isAt + 1300, 20000);
        const JunctionQuery q = eng.query(left, right, 1300, true);
        const TrackSupport ts = eng.relativeSupport(q);
        check(ts.agree == 5, "lineage IS: the 5 clonal relatives carry it at the same size");
        std::vector<TrackSupport> panel;
        eng.panelSupport({q}, panel);
        check(panel[0].placed >= 18, "lineage IS: most panel genomes place both flanks");
        check(panel[0].agree == 8 && panel[0].gapDiff >= 10,
              "lineage IS: 8 carriers agree, far lineages differ (non-positional)");
        check(panel[0].agree < 0.9 * panel[0].placed, "lineage IS locus is not positional");
        // a clone MISSING this IS: its flanks are adjacent; the relatives say 1300 -> gap difference
        const TrackSupport tm = eng.relativeSupport(eng.query(left, right, 0, true));
        check(tm.gapDiff == 5 && tm.contra == 0, "clone missing the IS: size difference, not contradiction");
    }
    // an operon locus: the gap recurs across the whole panel (positional)
    {
        const size_t o1 = shift(p, p.op[1], true);   // start of operon 2 in the lineage genome
        const std::string left = iso.substr(o1 - 20000, 20000);
        const std::string right = iso.substr(o1 + 5000, 20000);
        const JunctionQuery q = eng.query(left, right, 5000, true);
        std::vector<TrackSupport> panel;
        eng.panelSupport({q}, panel);
        check(panel[0].placed >= 18 && panel[0].agree >= 0.9 * panel[0].placed, "operon locus is positional");
    }
    // a NOVEL insertion in the isolate (no panel genome has it): relatives see adjacent flanks
    {
        const size_t at = shift(p, p.novelSite, true);
        const std::string left = iso.substr(at - 20000, 20000);
        const std::string right = iso.substr(at, 20000);
        const TrackSupport ts = eng.relativeSupport(eng.query(left, right, 1300, true));
        check(ts.gapDiff == 5, "novel IS: relatives adjacent, the asserted element is a size difference");
        bool zero = true;
        for (int64_t g : ts.gaps) if (std::llabs(g) > 20) zero = false;
        check(zero, "novel IS: relatives' gap ~0");
    }
    // side placement needs markers: an empty context places nothing
    {
        const SidePlace s = eng.placeOnRelative({}, true, 0, 0, 2);
        check(!s.ok, "no markers: not placed");
        int64_t g;
        check(ClonalEngine::call(s, s, 1000, 0, true, 100000, 1000, g) == TrackCall::Absent, "absent side: Absent");
    }
}

// ---------------------------------------------------------------------------------- plasmid sidecar

// A sidecar line for one relative's plasmid: its single-copy markers in order (the builder's format).
std::string sidecarLine(const OrganismModel& m, const std::string& genome, const std::string& rec, const std::string& seq) {
    std::vector<std::pair<uint64_t, std::pair<uint32_t, int>>> hits;
    forEachMarkerKmer(seq, [&](uint64_t km, uint32_t pos, int o) { hits.push_back({km, {pos, o}}); }, kDenom);
    std::unordered_map<uint64_t, int> occ;
    for (const auto& h : hits) ++occ[h.first];
    std::string ids, pos;
    for (const auto& h : hits) {
        if (occ[h.first] != 1) continue;
        const uint32_t id = m.markerOf(h.first);
        if (id == UINT32_MAX) continue;
        if (!ids.empty()) { ids += ','; pos += ','; }
        ids += std::to_string((id << 1) | static_cast<uint32_t>(h.second.second));
        pos += std::to_string(h.second.first);
    }
    return "P\t" + genome + "\t" + rec + "\t" + std::to_string(seq.size()) + "\t1\t" + ids + "\t" + pos + "\n";
}

void testPlasmidSidecar() {
    OrganismModel m;
    Panel p = makePanel(m);
    // the lineage plasmid, interned into the model as the builder's --plasmids would
    const std::string pl = dna(40000);
    std::vector<std::string> copies;
    for (int i = 0; i < 8; ++i) {
        copies.push_back(mutate(pl, 0.0002 * (i + 1)));
        forEachMarkerKmer(copies.back(), [&](uint64_t km, uint32_t, int) { m.internMarker(km); }, kDenom);
    }
    const std::string iso = mutate(p.chrs[0], 0.0001);
    const std::string isoPl = mutate(copies[0], 0.0001);
    ClonalOptions o;
    o.enabled = true;
    ClonalEngine eng(m, o);
    eng.selectNearestFromSeqs({iso});
    char tmpl[] = "/tmp/om2clnrpXXXXXX";
    const int fd = mkstemp(tmpl);
    check(fd >= 0, "temp file");
    close(fd);
    auto write = [&](const std::string& md5, uint64_t denom) {
        std::ofstream f(tmpl);
        f << "#om2nrp\t1\n#tsm_md5\t" << md5 << "\n#tracks_sha256\tx\n#hold_md5\tNA\n#scenario\tt\n#denom\t" << denom
          << "\n#genomes\t8\n#records\t8\n";
        for (int i = 0; i < 8; ++i) f << sidecarLine(m, "CLONE" + std::to_string(i + 1), "PL" + std::to_string(i + 1), copies[i]);
    };
    std::string err;
    write("0123456789abcdef0123456789abcdef", kDenom);
    check(!eng.loadPlasmids(tmpl, "ffffffffffffffffffffffffffffffff", err) && eng.plasmidTracks() == 0,
          "sidecar pinned to another model: refused");
    write("0123456789abcdef0123456789abcdef", 512);
    check(!eng.loadPlasmids(tmpl, "0123456789abcdef0123456789abcdef", err), "sidecar at another marker density: refused");
    write("0123456789abcdef0123456789abcdef", kDenom);
    check(eng.loadPlasmids(tmpl, "0123456789abcdef0123456789abcdef", err), "sidecar loads (" + err + ")");
    check(eng.plasmidTracks() == 5, "only the k = 5 nearest relatives' plasmids are kept");
    eng.setAssembly({iso, isoPl});
    // a plasmid junction cut out of the isolate's plasmid is placed on the relatives' plasmids at its true gap,
    // including across the plasmid's own origin (circular)
    for (const size_t at : {20000u, 39000u}) {
        const std::string twice = isoPl + isoPl + isoPl;
        const std::string left = twice.substr(at + isoPl.size() - 8000, 8000);
        const std::string right = twice.substr(at + isoPl.size() + 300, 8000);
        const std::vector<MarkerHit> hl = eng.hits(left, true), hr = eng.hits(right, false);
        size_t agree = 0;
        for (size_t pi = 0; pi < eng.plasmidTracks(); ++pi) {
            const SidePlace a = eng.placeOnPlasmid(hl, true, pi, 2), b = eng.placeOnPlasmid(hr, false, pi, 2);
            int64_t g;
            if (ClonalEngine::call(a, b, eng.plasmidLength(pi), 300, true, 100000, 1000, g) == TrackCall::Agree &&
                std::llabs(g - 300) <= 20)
                ++agree;
        }
        check(agree == 5, "plasmid junction at " + std::to_string(at) + ": 5/5 relatives' plasmids place it");
    }
    std::remove(tmpl);
}

// ---------------------------------------------------------------------------------- confidence

void testConfidence() {
    ClonalCalibration cal;
    ClonalAnnot a;
    a.set = true;
    a.cls = kPositional;
    a.nrPlaced = 5; a.nrAgree = 5; a.group = 'C';
    a.sized = true;
    const double pAll = cal.confidence(a, 'E', true);
    check(pAll > 0.99, "panel-only join, all 5 clonal relatives agree: above 0.99 under the prior");
    ClonalAnnot u = a;
    u.sized = false;
    check(cal.confidence(u, 'E', true) <= 0.98, "unsized: never confident");
    ClonalAnnot c = a;
    c.nrAgree = 0; c.nrContra = 4;
    check(cal.confidence(c, 'E', true) < 0.6, "contradicting relatives: low");
    ClonalAnnot pc = a;
    pc.pairContra = 5;
    check(cal.confidence(pc, 'E', true) < pAll, "contradicting pairs lower it");
    ClonalAnnot far = a;
    far.group = 'F';
    check(cal.confidence(far, 'E', true) < pAll, "far relatives: lower than clonal");
    check(a.nrBin() == "all" && c.nrBin() == "contra", "nr bins");
    // weighted bins: one relative placed is its own bin, and a panel-sized join it agrees with is not confident
    ClonalAnnot one = a;
    one.nrPlaced = 1; one.nrAgree = 1; one.wPlaced = 1; one.wAgree = 1; one.wContra = 0;
    check(one.nrBin() == "one" && cal.confidence(one, 'E', true) < 0.99, "one relative: bin 'one', not confident");
    ClonalAnnot two = a;
    two.nrPlaced = 2; two.nrAgree = 1; two.wPlaced = 1.05; two.wAgree = 1; two.wContra = 0;
    check(two.nrBin() == "all", "weighted: nearest agrees, a far one placed silent: 'all'");
    // a dev table replaces the prior (n >= 20 only)
    char tmpl[] = "/tmp/om2clcalXXXXXX";
    const int fd = mkstemp(tmpl);
    check(fd >= 0, "temp file");
    close(fd);
    {
        std::ofstream f(tmpl);
        f << "key\tp\tn\n" << ClonalCalibration::key(a, 'E') << "\t0.9123\t40\n"
          << ClonalCalibration::key(far, 'E') << "\t0.5\t3\n";
    }
    std::string err;
    check(cal.load(tmpl, err), "calibration table loads");
    check(std::fabs(cal.confidence(a, 'E', true) - 0.9123) < 1e-9, "table value used (n >= 20)");
    check(cal.confidence(far, 'E', true) != 0.5, "table row with n < 20 ignored");
    std::remove(tmpl);
}

void testOff() {
    const ClonalOptions o = ClonalOptions::fromEnv();
    check(!o.enabled && !o.anyFlagSet, "flags unset: off");
    ClonalStats s;
    const std::string line = formatClonalCounters(s);
    check(line.rfind("[om2-clonal] enabled=0 ran=0", 0) == 0, "counter line enabled=0");
    check(line.find("joins=0") != std::string::npos && line.find("annotated=0") != std::string::npos, "zeros");
    setenv("TESSERACT_OM2_CLONAL", "1", 1);
    setenv("TESSERACT_OM2_CLONAL_KMIN", "9", 1);
    const ClonalOptions on = ClonalOptions::fromEnv();
    check(on.enabled && on.anyFlagSet && on.kmin == on.k, "k_min is capped at k");
    unsetenv("TESSERACT_OM2_CLONAL");
    unsetenv("TESSERACT_OM2_CLONAL_KMIN");
}

}  // namespace

int main() {
    testenv::clearTesseractEnv();
    try {
        testRule();
        testRound2();
        testRound3();
        testRound3c();
        testPanels();
        testPlasmidSidecar();
        testConfidence();
        testOff();
    } catch (const std::exception& e) {
        std::fprintf(stderr, "FAIL after %d checks: %s\n", checks, e.what());
        return 1;
    }
    std::printf("test_om2_clonal: %d checks passed\n", checks);
    return 0;
}
