// Organism Model 2.0, component C4 (CLONAL): clonal-anchored closure. See om2_clonal.h.
#include "om2_clonal.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <sstream>

#include "envflags.h"
#include "graph.h"
#include "om2_evidence.h"
#include "om2_ledger.h"
#include "organism.h"
#include "util.h"

namespace ts {
namespace om2 {

namespace {

double nowMs() {
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

std::string fmtD(double v, const char* f = "%.4g") {
    char b[64];
    std::snprintf(b, sizeof b, f, v);
    return b;
}

// Canonical base check without the table: markerBaseCode.
bool isBase(char c) { return markerBaseCode(c) >= 0; }

int64_t median(std::vector<int64_t> v) {
    if (v.empty()) return INT64_MIN;
    std::sort(v.begin(), v.end());
    return v[(v.size() - 1) / 2];
}

// The N-free stretch of `s` touching the junction: its end (exit side) or its start (entry side),
// at most `w` bases.
std::string nfreeTail(const std::string& s, size_t w) {
    size_t from = s.size() > w ? s.size() - w : 0;
    for (size_t i = s.size(); i > from; --i) {
        if (!isBase(s[i - 1])) { from = i; break; }
    }
    return s.substr(from);
}
std::string nfreeHead(const std::string& s, size_t w) {
    size_t to = std::min(w, s.size());
    for (size_t i = 0; i < to; ++i) {
        if (!isBase(s[i])) { to = i; break; }
    }
    return s.substr(0, to);
}

char tierChar(Tier t) {
    switch (t) {
        case Tier::A: return 'A';
        case Tier::B: return 'B';
        case Tier::C: return 'C';
        case Tier::D: return 'D';
        case Tier::E: return 'E';
    }
    return 'E';
}

}  // namespace

// ------------------------------------------------------------------------------ options

ClonalOptions ClonalOptions::fromEnv() {
    ClonalOptions o;
    o.enabled = env::on("TESSERACT_OM2_CLONAL", false);
    o.k = static_cast<int>(env::integer("TESSERACT_OM2_CLONAL_K", 5));
    o.kmin = static_cast<int>(env::integer("TESSERACT_OM2_CLONAL_KMIN", 3));
    o.layout = env::on("TESSERACT_OM2_CLONAL_LAYOUT", true);
    o.breakContra = env::on("TESSERACT_OM2_CLONAL_BREAK", true);
    if (const char* np = env::text("TESSERACT_OM2_CLONAL_NONPOS")) o.nonposFill = std::strcmp(np, "gap") != 0;
    o.dmax = env::real("TESSERACT_OM2_CLONAL_DMAX", 1.0);
    o.confident = env::real("TESSERACT_OM2_CLONAL_CONFIDENT", 0.99);
    o.confDmax = env::real("TESSERACT_OM2_CLONAL_CONF_DMAX", 1.0);
    if (const char* p = env::text("TESSERACT_OM2_CLONAL_CAL")) o.calPath = p;
    o.minMarkers = static_cast<int>(env::integer("TESSERACT_OM2_CLONAL_MINMARK", 3));
    if (const char* w = env::text("TESSERACT_OM2_CLONAL_WEIGHT")) o.weighted = std::strcmp(w, "count") != 0;
    o.dscale = env::real("TESSERACT_OM2_CLONAL_DSCALE", 0.0005);
    if (const char* p = env::text("TESSERACT_OM2_CLONAL_NRP")) o.nrpPath = p;
    o.overlapMerge = env::on("TESSERACT_OM2_CLONAL_OVERLAP", true);
    o.vouch = env::on("TESSERACT_OM2_CLONAL_VOUCH", true);
    o.bracket = env::on("TESSERACT_OM2_CLONAL_BRACKET", false);
    o.resize = env::on("TESSERACT_OM2_CLONAL_RESIZE", true);
    o.unverified = env::on("TESSERACT_OM2_CLONAL_UNVERIFIED", false);
    o.indel = env::on("TESSERACT_OM2_CLONAL_INDEL", false);
    o.walksize = env::on("TESSERACT_OM2_CLONAL_WALKSIZE", false);
    o.kc2 = env::on("TESSERACT_OM2_CLONAL_KC2", false);
    o.confRule = env::on("TESSERACT_OM2_CLONAL_CONFRULE", false);
    o.discord = static_cast<int>(env::integer("TESSERACT_OM2_CLONAL_DISCORD", 2));
    o.strict = env::on("TESSERACT_OM2_CLONAL_STRICT", false);
    o.strictNr0 = env::on("TESSERACT_OM2_CLONAL_STRICT_NR0", false);
    o.allowConfident = env::on("TESSERACT_OM2_CLONAL_ALLOW_CONFIDENT", false);
    for (const char* f : {"TESSERACT_OM2_CLONAL_STRICT", "TESSERACT_OM2_CLONAL_STRICT_NR0", "TESSERACT_OM2_CLONAL_INDEL", "TESSERACT_OM2_CLONAL_UNVERIFIED", "TESSERACT_OM2_CLONAL_CAPWALK", "TESSERACT_OM2_CLONAL_WALKSIZE", "TESSERACT_OM2_CLONAL_KC2", "TESSERACT_OM2_CLONAL_CONFRULE",
                          "TESSERACT_OM2_CLONAL_DISCORD", "TESSERACT_OM2_CLONAL_RESIZE", "TESSERACT_OM2_CLONAL_OVERLAP", "TESSERACT_OM2_CLONAL_VOUCH", "TESSERACT_OM2_CLONAL_BRACKET",
                          "TESSERACT_OM2_CLONAL_NRP", "TESSERACT_OM2_CLONAL_WEIGHT", "TESSERACT_OM2_CLONAL_DSCALE", "TESSERACT_OM2_CLONAL_K", "TESSERACT_OM2_CLONAL_KMIN", "TESSERACT_OM2_CLONAL_LAYOUT",
                          "TESSERACT_OM2_CLONAL_BREAK", "TESSERACT_OM2_CLONAL_NONPOS", "TESSERACT_OM2_CLONAL_DMAX",
                          "TESSERACT_OM2_CLONAL_CONFIDENT", "TESSERACT_OM2_CLONAL_CONF_DMAX",
                          "TESSERACT_OM2_CLONAL_CAL", "TESSERACT_OM2_CLONAL_MINMARK"}) {
        if (env::isSet(f)) o.anyFlagSet = true;
    }
    if (o.kmin > o.k) o.kmin = o.k;
    o.report = std::max(10, o.k);
    return o;
}

const char* clonalClassName(int c) {
    switch (c) {
        case kPositional: return "positional";
        case kNonPositional: return "non_positional";
        case kLayout: return "layout";
        case kPlasmid: return "plasmid";
        case kWrap: return "wrap";
        case kGraph: return "graph";
        default: return "unknown";
    }
}

char NearestSet::group() const {
    if (ranked.empty()) return '-';
    const double d = dNear();
    return d <= 0.0005 ? 'C' : d <= 0.002 ? 'L' : 'F';
}

std::string ClonalAnnot::nrBin() const {
    if (wPlaced >= 0) {   // distance-weighted votes
        if (wPlaced < 0.5 || nrPlaced == 0) return "none";
        if (wContra > wAgree) return "contra";
        // One relative alone is not independent of a panel-sized join (the layout's template is usually the nearest
        // relative itself): its own bin, with a small odds factor (smoke4 #168: a 1/1 'agreement' on a FALSE size).
        if (nrPlaced < 2) return nrAgree ? "one" : "min";
        if (wAgree >= 0.9 * wPlaced) return "all";
        if (wAgree >= 0.6 * wPlaced) return "maj";
        return "min";
    }
    if (nrPlaced < 2) return "none";
    if (nrContra > 0 && nrAgree * 2 < nrPlaced) return "contra";
    if (nrAgree == nrPlaced) return "all";
    if (nrAgree * 10 >= nrPlaced * 6) return "maj";
    return "min";
}

int64_t trackLength(const LayoutTrack& t, uint64_t denom) {
    if (t.pos.empty()) return 0;
    // The model keeps marker positions, not the chromosome length. The last sampled marker sits on
    // average one sampling interval before the end (geometric spacing, mean `denom`).
    return static_cast<int64_t>(t.pos.back()) + static_cast<int64_t>(denom ? denom : kMarkerSampleDenom);
}

// ------------------------------------------------------------------------------ the rule

// Round 3 (WALKSIZE): the shortest walk class is the isolate's own answer when it is isolated -- the only class of an
// exhaustive enumeration, or at least `gap` shorter than every other class (the longer classes loop through a repeat
// the flanks share) -- walked once, and at most `maxLen` long. Else INT32_MIN.
int32_t isolatedShortestWalk(const std::vector<int32_t>& lens, const std::vector<uint16_t>& counts, bool exhaustive,
                             int64_t gap, int64_t maxLen) {
    if (lens.empty()) return INT32_MIN;
    size_t i0 = 0;
    for (size_t i = 1; i < lens.size(); ++i) if (lens[i] < lens[i0]) i0 = i;
    const int32_t w0 = lens[i0];
    if (static_cast<int64_t>(w0) > maxLen) return INT32_MIN;
    if (!counts.empty() && i0 < counts.size() && counts[i0] > 1) return INT32_MIN;
    int64_t next = INT64_MAX;
    for (size_t i = 0; i < lens.size(); ++i) {
        if (i == i0) continue;
        next = std::min<int64_t>(next, lens[i]);
    }
    if (next == INT64_MAX) return exhaustive ? w0 : INT32_MIN;
    return next - static_cast<int64_t>(w0) >= gap ? w0 : INT32_MIN;
}

// Round 3c (STRICT): >= 3 walk classes whose first three spacings are equal within max(30, 5%): a tandem array (or a
// repeat the flanks share several times) whose copy number the walk lengths cannot fix. `sortedLens` ascending.
bool strictPeriodic(const std::vector<int32_t>& sortedLens) {
    if (sortedLens.size() < 3) return false;
    const size_t m = std::min<size_t>(4, sortedLens.size());
    const int64_t d0 = static_cast<int64_t>(sortedLens[1]) - sortedLens[0];
    if (d0 <= 0) return false;
    const int64_t tolP = std::max<int64_t>(30, d0 / 20);
    for (size_t i = 2; i < m; ++i)
        if (std::llabs(static_cast<int64_t>(sortedLens[i]) - sortedLens[i - 1] - d0) > tolP) return false;
    return true;
}

// Round 3c (STRICT): the emission rule of one model-proposed junction at its final size (see om2_clonal.h). The
// isolate's own graph sizes the gap, or every placed relative (>= 2) reads that size; otherwise the reason is returned.
// Dev round-3b counterfactual (vrule.py V8 / V10, Rule-Q): Kp FALSE 37 -> 26, Ab 171 -> 102, no chr98 join lost.
const char* strictVerdict(const StrictInput& s, const ClonalOptions& o) {
    if (s.overlap || s.n <= 1) return nullptr;
    std::vector<int32_t> W = s.walkLens;
    std::sort(W.begin(), W.end());
    auto tolW = [](int64_t x) { return std::max<int64_t>(500, std::llabs(x) / 10); };
    bool match = false;
    for (int32_t w : W) if (std::llabs(static_cast<int64_t>(w) - s.n) <= tolW(w)) match = true;
    const bool shortest = !W.empty() && std::llabs(static_cast<int64_t>(W[0]) - s.n) <= tolW(W[0]);
    const bool iso = !W.empty() && (W.size() == 1 || static_cast<int64_t>(W[1]) - W[0] >= o.c1GapLayout);
    const bool unique = W.size() == 1 && s.exhaustive;
    const bool nrAll = s.placed >= 2 && s.agree == s.placed && s.contra == 0;
    bool keep = false;
    if (match && shortest && (unique || !(s.wContra > s.wAgree && !iso))) keep = true;
    if (!keep && nrAll && (W.empty() || match)) keep = true;
    if (!keep) {
        if (W.empty()) return "strict_walkless";
        if (!match) return "strict_no_walk_class";
        return shortest ? "strict_walk_contradicted" : "strict_not_shortest";
    }
    if (!nrAll && !unique && match && shortest && strictPeriodic(W) &&
        !(s.wAgree >= o.strictPeriodicWa && s.wAgree >= o.strictPeriodicRatio * s.wContra))
        return "strict_periodic";
    if (o.strictNr0 && !nrAll && s.placed == 0) return "strict_no_relative";
    return nullptr;
}

LinkDecision decideLink(const LinkEvidence& e, const ClonalOptions& o) {
    LinkDecision d;
    const uint32_t kminEff = static_cast<uint32_t>(std::max(1, std::min<int>(o.kmin, static_cast<int>(e.k))));
    const bool panelKnown = e.panelPlaced >= o.panelMin;
    const bool positional = panelKnown && static_cast<double>(e.panelAgree) >= o.positional * e.panelPlaced;
    const bool uniqueEnds = !e.repeatA && !e.repeatB && !e.element;
    d.cls = e.plasmid ? kPlasmid : positional ? (uniqueEnds ? kLayout : kPositional) : kNonPositional;
    const bool agreed = o.weighted ? (e.wAgree >= o.wmin && e.wAgree >= o.agreeShare * (e.wAgree + e.wContra))
                                   : (e.agree >= kminEff && static_cast<double>(e.agree) >= o.agreeShare * e.placed);
    if (!agreed) { d.why = "nr_agree"; return d; }
    if (e.pairContra >= o.minContra) { d.why = "pair_contra"; return d; }
    // Round 2: the two ends share their bases EXACTLY at (about) the overlap the relatives give. The isolate's own
    // sequence then asserts the adjacency's size; a graph contradiction (a first-unique successor elsewhere, usually a
    // repeat inside the overlap) is overruled only when every placed relative agrees.
    if (o.overlapMerge && e.overlapExact >= o.minOverlap) {
        const bool strong = e.placed >= 2 && (o.weighted ? (e.wPlaced > 0 && e.wAgree >= 0.9 * e.wPlaced)
                                                          : e.agree >= e.placed);
        if (!e.contraGraph || strong) {
            d.n = -e.overlapExact;
            d.overlap = true;
            d.walkMatch = true;
            d.basis = "overlap_exact";
            d.join = true;
            d.why = "join";
            return d;
        }
    }
    if (e.contraGraph) { d.why = "graph_contra"; return d; }
    const int64_t tol = std::max<int64_t>(500, std::llabs(e.nrGap) / 10);
    int64_t n = 0;
    if (e.reachable && !e.walkLens.empty()) {
        int64_t best = INT64_MAX;
        for (int32_t w : e.walkLens) {
            const int64_t dd = std::llabs(static_cast<int64_t>(w) - e.nrGap);
            if (dd <= tol && dd < best) { best = dd; n = w; }
        }
        if (best != INT64_MAX) {
            d.walkMatch = true;
            d.basis = "graph+nr";
        } else {
            bool oneClass = true;
            for (int32_t w : e.walkLens) if (std::abs(w - e.walkLens.front()) > 100) oneClass = false;
            // Round 2: one walk class and the relatives' gap within bracketMax: their midpoint is within bracketMax / 2
            // of both candidates (computed always, used only with TESSERACT_OM2_CLONAL_BRACKET=1).
            if (oneClass && std::llabs(static_cast<int64_t>(e.walkLens.front()) - e.nrGap) <= o.bracketMax &&
                e.nrSpread <= tol)
                d.bracketN = static_cast<int32_t>(std::max<int64_t>(1, (static_cast<int64_t>(e.walkLens.front()) + e.nrGap) / 2));
            const int32_t c1 = (o.walksize && !e.plasmid)
                ? isolatedShortestWalk(e.walkLens, e.walkCounts, e.exhaustive, o.c1Gap, o.c1Max) : INT32_MIN;
            if (oneClass && e.exhaustive) {
                // The isolate's own walk disagrees with the relatives' gap: one side carries an element
                // the other lacks (a novel insertion, or a relative's element this isolate does not have).
                n = e.walkLens.front();
                d.basis = "graph";
                if (d.cls != kPlasmid) d.cls = kNonPositional;
            } else if (c1 != INT32_MIN) {
                // Round 3 (WALKSIZE): the relatives agree on the ORDER, their gap is no walk of this isolate, and the
                // isolate's shortest walk is isolated (the longer classes loop through a shared repeat): the isolate's
                // own walk sizes the join (dev Kp round 2: 24 of 26 truth-checked size_ambiguous refusals were true
                // adjacencies at exactly that walk).
                n = c1;
                d.c1 = true;
                d.basis = "graph_c1";
                if (d.cls != kPlasmid) d.cls = kNonPositional;
            } else if (o.bracket && !e.plasmid && d.bracketN != INT32_MIN) {
                n = d.bracketN;
                d.bracket = true;
                d.basis = "bracket";
            } else {
                d.why = "size_ambiguous";
                return d;
            }
        }
    } else {
        if (e.nrSpread > tol) { d.why = "nr_spread"; return d; }
        n = e.nrGap;
        d.basis = "nr";
    }
    if (n < -o.maxOverlap) { d.why = "overlap"; return d; }
    // A plasmid join needs positive evidence from this isolate (C1d's principle): its own graph walk at the size the
    // relatives' plasmids give. The relatives alone never size a plasmid join.
    if (e.plasmid && !d.walkMatch) { d.why = "plasmid_needs_walk"; return d; }
    d.n = static_cast<int32_t>(std::max<int64_t>(1, std::min<int64_t>(n, 2000000000)));
    // A non-positional ELEMENT (the walk passes repeat sequence) is written with bases only when the
    // nearest relatives carry it at the same size, both flank edges exist and no pair contradicts.
    if (d.cls == kNonPositional && e.element) {
        const bool carried = o.weighted ? e.wAgree >= o.wmin : e.agree >= kminEff;
        const bool ok = o.nonposFill && d.walkMatch && carried && e.edgeA && e.edgeB && e.pairContra == 0;
        d.noFill = !ok || o.kc2;   // round 3, KC2: a non-positional element is never written with bases
    }
    if (d.bracket) d.noFill = true;   // an uncertain size is never filled
    d.join = true;
    d.why = "join";
    return d;
}

// Round 2: the exact suffix/prefix overlap of `L` (junction at its end) and `R` (junction at its start) nearest to
// `want` bases, within max(30, want / 10) of it and at least minOv long; -1 when there is none.
int32_t exactEndOverlap(const std::string& L, const std::string& R, int64_t want, int32_t minOv, int64_t maxOv) {
    if (want <= 0) return -1;
    const int64_t tol = std::max<int64_t>(30, want / 10);
    const int64_t lo = std::max<int64_t>(minOv, want - tol);
    const int64_t hi = std::min<int64_t>({maxOv, want + tol, static_cast<int64_t>(L.size()), static_cast<int64_t>(R.size())});
    int32_t best = -1;
    int64_t bestD = INT64_MAX;
    for (int64_t ov = lo; ov <= hi; ++ov) {
        const size_t a = L.size() - static_cast<size_t>(ov);
        bool ok = true;
        for (int64_t i = 0; i < ov && ok; ++i) {
            const char c = L[a + static_cast<size_t>(i)];
            ok = c == R[static_cast<size_t>(i)] && c != 'N' && c != 'n';
        }
        if (!ok) continue;
        const int64_t dd = std::llabs(ov - want);
        if (dd < bestD) { bestD = dd; best = static_cast<int32_t>(ov); }
    }
    return best;
}

// ------------------------------------------------------------------------------ calibration

bool ClonalCalibration::load(const std::string& path, std::string& err) {
    std::ifstream in(path);
    if (!in) { err = "cannot open " + path; return false; }
    std::string l;
    while (std::getline(in, l)) {
        if (l.empty() || l[0] == '#') continue;
        std::vector<std::string> c;
        std::stringstream ss(l);
        std::string f;
        while (std::getline(ss, f, '\t')) c.push_back(f);
        // key <TAB> p <TAB> n
        if (c.size() < 3 || c[0] == "key") continue;
        table[c[0]] = {std::atof(c[1].c_str()), static_cast<uint32_t>(std::atoi(c[2].c_str()))};
    }
    loaded = !table.empty();
    if (!loaded) err = path + ": no rows";
    return loaded;
}

std::string ClonalCalibration::key(const ClonalAnnot& a, char tier) {
    return std::string(clonalClassName(a.cls)) + "|" + tier + "|" + a.nrBin() + "|" + a.group;
}

double ClonalCalibration::confidence(const ClonalAnnot& a, char tier, bool hasLedger) const {
    double p = -1;
    if (loaded) {
        auto it = table.find(key(a, tier));
        if (it != table.end() && it->second.second >= 20) p = it->second.first;
    }
    if (p < 0) {
        // PROPOSED prior, to be replaced by the dev-fitted table (EVAL_PLAN_CLONAL s6.7). Base rates by
        // evidence tier (DESIGN s1: panel-only SILENT joins right 96.6-98.7%, tie-breaks 93-96%), then
        // odds multipliers for the nearest relatives, the pairs and the graph.
        double p0;
        if (a.cls == kGraph) p0 = 0.99;                // read-pair scaffold of the resolver
        else if (!hasLedger) p0 = 0.95;
        else switch (tier) {
            case 'A': p0 = 0.999; break;
            case 'B': p0 = 0.995; break;
            case 'C': p0 = 0.99; break;
            case 'D': p0 = 0.97; break;
            default: p0 = 0.95; break;
        }
        double odds = p0 / (1 - p0);
        const std::string bin = a.nrBin();
        if (bin == "contra") odds *= 0.05;
        else if (bin == "all") odds *= a.group == 'C' ? 20 : a.group == 'L' ? 5 : 2;
        else if (bin == "maj") odds *= 1.5;
        else if (bin == "one") odds *= a.group == 'C' ? 2 : 1;
        else if (bin == "min") odds *= 0.5;
        if (a.pairContra >= 3) odds *= 0.05;
        if (a.sized && a.walk == 1 && a.walkMatch == 0) odds *= 0.2;
        p = odds / (1 + odds);
    }
    if (!a.sized) p = std::min(p, 0.98);
    return p;
}

// ------------------------------------------------------------------------------ engine

ClonalEngine::ClonalEngine(const OrganismModel& model, const ClonalOptions& opt) : model_(model), opt_(opt) {}

double ClonalEngine::weight(size_t rank) const {
    if (!opt_.weighted || rank >= near_.ranked.size() || near_.ranked.empty()) return 1.0;
    const double dd = near_.ranked[rank].dist - near_.ranked.front().dist;
    return std::exp(-std::max(0.0, dd) / std::max(1e-9, opt_.dscale));
}

const std::string& ClonalEngine::relativeName(size_t rank) const {
    static const std::string none;
    return rank < nr_.size() ? model_.tracks()[nr_[rank].track].name : none;
}
int64_t ClonalEngine::relativeLength(size_t rank) const { return rank < nr_.size() ? nr_[rank].length : 0; }

void ClonalEngine::selectNearest(const UnitigGraph& g) {
    // Streamed, not copied: every live unitig is scanned once for sampled markers.
    const double t0 = nowMs();
    std::vector<uint8_t> present(model_.markerCount(), 0);
    const uint64_t denom = model_.markerDenom();
    for (const Unitig& u : g.nodes) {
        if (u.deleted) continue;
        forEachMarkerKmer(u.seq, [&](uint64_t km, uint32_t, int) {
            const uint32_t id = model_.markerOf(km);
            if (id != UINT32_MAX && id < present.size()) present[id] = 1;
        }, denom);
    }
    near_ = NearestSet();
    near_.tracks = static_cast<uint32_t>(model_.trackCount());
    for (size_t id = 0; id < present.size(); ++id) {
        if (!present[id]) continue;
        ++near_.isoMarkers;
        if (model_.markerGenomes(static_cast<uint32_t>(id), Replicon::Chromosome) > 0) ++near_.isoChrMarkers;
    }
    std::vector<Relative> all;
    all.reserve(model_.trackCount());
    for (size_t t = 0; t < model_.tracks().size(); ++t) {
        const LayoutTrack& tr = model_.tracks()[t];
        Relative r;
        r.track = static_cast<uint32_t>(t);
        r.name = tr.name;
        r.trackMarkers = static_cast<uint32_t>(tr.oriented.size());
        for (uint32_t o : tr.oriented) if ((o >> 1) < present.size() && present[o >> 1]) ++r.shared;
        const double uni = static_cast<double>(r.trackMarkers) + near_.isoChrMarkers - r.shared;
        r.jaccard = uni > 0 ? r.shared / uni : 0;
        r.containment = r.trackMarkers ? static_cast<double>(r.shared) / r.trackMarkers : 0;
        // Distance from the CONTAINMENT of the relative's chromosome markers in the isolate: the Jaccard form counts
        // the isolate's plasmid-borne and accessory markers against every chromosome and ran ~6x the reference
        // mash distance on dev Kp (GCF055901135v2: Jaccard 0.00121, containment 0.00021, panel mash 0.00018).
        r.dist = r.containment > 0 ? -std::log(r.containment) / kMarkerK : 1.0;
        if (r.dist < 0) r.dist = 0;
        all.push_back(std::move(r));
    }
    std::sort(all.begin(), all.end(), [](const Relative& a, const Relative& b) {
        if (a.dist != b.dist) return a.dist < b.dist;
        if (a.shared != b.shared) return a.shared > b.shared;
        return a.name < b.name;
    });
    if (all.size() > static_cast<size_t>(opt_.report)) all.resize(static_cast<size_t>(opt_.report));
    near_.ranked = std::move(all);
    near_.ok = !near_.ranked.empty() && near_.ranked.front().shared > 0;
    near_.ms = nowMs() - t0;
    buildRelatives();
}

void ClonalEngine::selectNearestFromSeqs(const std::vector<std::string>& seqs) {
    UnitigGraph g;
    for (const std::string& s : seqs) {
        Unitig u;
        u.seq = s;
        g.nodes.push_back(std::move(u));
    }
    selectNearest(g);
}

void ClonalEngine::buildRelatives() {
    nr_.clear();
    if (!near_.ok) return;
    const size_t k = std::min<size_t>(static_cast<size_t>(std::max(1, opt_.k)), near_.ranked.size());
    for (size_t r = 0; r < k; ++r) {
        NrTrack t;
        t.track = near_.ranked[r].track;
        const LayoutTrack& tr = model_.tracks()[t.track];
        t.at.reserve(tr.oriented.size() * 2);
        for (uint32_t i = 0; i < tr.oriented.size(); ++i) t.at.emplace(tr.oriented[i] >> 1, i);
        t.length = trackLength(tr, model_.markerDenom());
        nr_.push_back(std::move(t));
    }
}

void ClonalEngine::setAssembly(const std::vector<std::string>& seqs) {
    asmCount_.clear();
    for (const std::string& s : seqs) {
        forEachMarkerKmer(s, [&](uint64_t km, uint32_t, int) {
            const uint32_t id = model_.markerOf(km);
            if (id == UINT32_MAX) return;
            uint8_t& c = asmCount_[id];
            if (c < 2) ++c;
        }, model_.markerDenom());
    }
    asmSet_ = true;
}

std::vector<MarkerHit> ClonalEngine::hits(const std::string& s, bool exitSide) const {
    const std::string ctx = exitSide ? nfreeTail(s, static_cast<size_t>(opt_.window))
                                     : nfreeHead(s, static_cast<size_t>(opt_.window));
    std::vector<MarkerHit> out;
    std::unordered_map<uint32_t, uint8_t> local;
    const int64_t L = static_cast<int64_t>(ctx.size());
    forEachMarkerKmer(ctx, [&](uint64_t km, uint32_t pos, int orient) {
        const uint32_t id = model_.markerOf(km);
        if (id == UINT32_MAX) return;
        if (asmSet_) {
            auto it = asmCount_.find(id);
            if (it == asmCount_.end() || it->second != 1) return;
        }
        ++local[id];
        MarkerHit h;
        h.id = id;
        h.pos = static_cast<int32_t>(pos);
        h.o = static_cast<uint8_t>(orient);
        h.dist = exitSide ? static_cast<int32_t>(L - (static_cast<int64_t>(pos) + kMarkerK)) : static_cast<int32_t>(pos);
        out.push_back(h);
    }, model_.markerDenom());
    // A marker twice in one context names a local repeat.
    out.erase(std::remove_if(out.begin(), out.end(), [&](const MarkerHit& h) { return local[h.id] != 1; }), out.end());
    std::sort(out.begin(), out.end(), [](const MarkerHit& a, const MarkerHit& b) {
        return a.dist != b.dist ? a.dist < b.dist : a.id < b.id;
    });
    return out;
}

JunctionQuery ClonalEngine::query(const std::string& left, const std::string& right, int64_t asserted, bool sized) const {
    JunctionQuery q;
    q.left = hits(left, true);
    q.right = hits(right, false);
    q.leftLen = static_cast<int64_t>(left.size());
    q.asserted = asserted;
    q.sized = sized;
    return q;
}

SidePlace ClonalEngine::placeWith(const std::vector<MarkerHit>& h, bool exitSide, int64_t trackLen,
                                  const LayoutTrack& t, const std::function<bool(uint32_t, uint32_t&)>& find,
                                  size_t need) const {
    struct C { int64_t x; int8_t dir; int32_t dist; };
    std::vector<C> cs;
    const size_t use = static_cast<size_t>(std::max(1, opt_.useMarkers));
    for (const MarkerHit& m : h) {
        uint32_t i;
        if (!find(m.id, i)) continue;
        const int64_t P = static_cast<int64_t>(t.pos[i]);
        const uint8_t oT = static_cast<uint8_t>(t.oriented[i] & 1u);
        const bool same = m.o == oT;
        C c;
        c.dist = m.dist;
        if (exitSide) {
            if (same) { c.x = P + m.dist + kMarkerK; c.dir = 1; }
            else { c.x = P - m.dist - 1; c.dir = -1; }
        } else {
            if (same) { c.x = P - m.dist; c.dir = 1; }
            else { c.x = P + (kMarkerK - 1) + m.dist; c.dir = -1; }
        }
        cs.push_back(c);
        if (cs.size() >= use) break;
    }
    SidePlace sp;
    if (cs.size() < need) return sp;
    auto circ = [&](int64_t d) {
        if (trackLen > 0) {
            while (d > trackLen / 2) d -= trackLen;
            while (d < -trackLen / 2) d += trackLen;
        }
        return d;
    };
    size_t bestI = SIZE_MAX, bestN = 0;
    for (size_t i = 0; i < cs.size(); ++i) {
        size_t n = 0;
        for (size_t j = 0; j < cs.size(); ++j) {
            if (cs[j].dir != cs[i].dir) continue;
            const int64_t tol = 200 + static_cast<int64_t>(0.03 * std::max(cs[i].dist, cs[j].dist));
            if (std::llabs(circ(cs[j].x - cs[i].x)) <= tol) ++n;
        }
        if (n > bestN) { bestN = n; bestI = i; }   // first (nearest) wins ties
    }
    if (bestI == SIZE_MAX || bestN < need) return sp;
    sp.ok = true;
    sp.x = cs[bestI].x;
    if (trackLen > 0) { sp.x %= trackLen; if (sp.x < 0) sp.x += trackLen; }
    sp.dir = cs[bestI].dir;
    sp.n = static_cast<uint16_t>(std::min<size_t>(bestN, 65535));
    sp.nearest = cs[bestI].dist;
    return sp;
}

SidePlace ClonalEngine::placeOnRelative(const std::vector<MarkerHit>& h, bool exitSide, int64_t, size_t rank,
                                        size_t need) const {
    if (rank >= nr_.size()) return SidePlace();
    const NrTrack& nt = nr_[rank];
    const LayoutTrack& t = model_.tracks()[nt.track];
    return placeWith(h, exitSide, nt.length, t, [&](uint32_t id, uint32_t& i) {
        auto it = nt.at.find(id);
        if (it == nt.at.end()) return false;
        i = it->second;
        return true;
    }, need);
}

bool ClonalEngine::loadPlasmids(const std::string& path, const std::string& modelMd5, std::string& err) {
    pl_.clear();
    std::ifstream in(path);
    if (!in) { err = "cannot open " + path; return false; }
    std::map<std::string, size_t> rankOf;
    for (size_t r = 0; r < nr_.size(); ++r) rankOf[model_.tracks()[nr_[r].track].name] = r;
    std::string l, md5;
    bool magic = false;
    while (std::getline(in, l)) {
        if (l.empty()) continue;
        if (l[0] == '#') {
            if (l.compare(0, 8, "#om2nrp\t") == 0) magic = true;
            else if (l.compare(0, 9, "#tsm_md5\t") == 0) md5 = l.substr(9);
            else if (l.compare(0, 7, "#denom\t") == 0 &&
                     static_cast<uint64_t>(std::atoll(l.c_str() + 7)) != model_.markerDenom()) {
                err = path + ": marker density differs from the model's";
                pl_.clear();
                return false;
            }
            continue;
        }
        if (!magic || md5 != modelMd5) {
            err = path + (magic ? ": #tsm_md5 " + md5 + " is not the loaded model's " + modelMd5 : ": not an om2nrp sidecar");
            pl_.clear();
            return false;
        }
        // P genome record length n ids pos
        std::vector<std::string> c;
        std::stringstream ss(l);
        std::string f;
        while (std::getline(ss, f, '\t')) c.push_back(f);
        if (c.size() < 7 || c[0] != "P") continue;
        auto it = rankOf.find(c[1]);
        if (it == rankOf.end() || c[5] == ".") continue;
        PlTrack t;
        t.rank = it->second;
        t.t.name = c[1] + ":" + c[2];
        t.length = std::atoll(c[3].c_str());
        std::stringstream si(c[5]), sp(c[6]);
        std::string a, b;
        while (std::getline(si, a, ',') && std::getline(sp, b, ',')) {
            const uint32_t o = static_cast<uint32_t>(std::strtoul(a.c_str(), nullptr, 10));
            if ((o >> 1) >= model_.markerCount()) continue;
            t.at.emplace(o >> 1, static_cast<uint32_t>(t.t.oriented.size()));
            t.t.oriented.push_back(o);
            t.t.pos.push_back(static_cast<uint32_t>(std::strtoul(b.c_str(), nullptr, 10)));
        }
        if (!t.t.oriented.empty()) pl_.push_back(std::move(t));
    }
    if (!magic || md5 != modelMd5) {
        err = path + (magic ? ": #tsm_md5 " + md5 + " is not the loaded model's " + modelMd5 : ": not an om2nrp sidecar");
        pl_.clear();
        return false;
    }
    std::sort(pl_.begin(), pl_.end(), [](const PlTrack& a, const PlTrack& b) {
        return a.rank != b.rank ? a.rank < b.rank : a.t.name < b.t.name;
    });
    return true;
}

SidePlace ClonalEngine::placeOnPlasmid(const std::vector<MarkerHit>& h, bool exitSide, size_t pi, size_t need) const {
    if (pi >= pl_.size()) return SidePlace();
    const PlTrack& p = pl_[pi];
    return placeWith(h, exitSide, p.length, p.t, [&](uint32_t id, uint32_t& i) {
        auto it = p.at.find(id);
        if (it == p.at.end()) return false;
        i = it->second;
        return true;
    }, need);
}

TrackCall ClonalEngine::call(const SidePlace& a, const SidePlace& b, int64_t trackLen, int64_t asserted, bool sized,
                             int64_t maxAdjacent, int64_t maxOverlap, int64_t& gap) {
    gap = INT64_MIN;
    if (!a.ok || !b.ok) return TrackCall::Absent;
    if (a.dir != b.dir) return TrackCall::Contra;
    int64_t g = (b.x - a.x) * a.dir;
    if (trackLen > 0) {
        while (g < -trackLen / 2) g += trackLen;
        while (g > trackLen / 2) g -= trackLen;
    }
    gap = g;
    if (g < -maxOverlap || g > maxAdjacent) return TrackCall::Contra;
    if (!sized) return TrackCall::Agree;
    const int64_t tol = std::max<int64_t>(1000, std::max(std::llabs(g), std::llabs(asserted)) / 10);
    return std::llabs(g - asserted) <= tol ? TrackCall::Agree : TrackCall::GapDiff;
}

TrackSupport ClonalEngine::relativeSupport(const JunctionQuery& q) const {
    TrackSupport s;
    for (size_t r = 0; r < nr_.size(); ++r) {
        ++s.tracks;
        const SidePlace a = placeOnRelative(q.left, true, q.leftLen, r, 2);
        const SidePlace b = placeOnRelative(q.right, false, 0, r, 2);
        int64_t g;
        const TrackCall c = call(a, b, nr_[r].length, q.asserted, q.sized, opt_.maxAdjacent, opt_.maxOverlap, g);
        if (c == TrackCall::Absent) continue;
        const double w = weight(r);
        ++s.placed;
        s.wPlaced += w;
        if (g != INT64_MIN && c != TrackCall::Contra) s.gaps.push_back(g);
        if (c == TrackCall::Agree) { ++s.agree; s.wAgree += w; s.agreeGaps.push_back(g); s.agreeRank.push_back(static_cast<uint32_t>(r)); }
        else if (c == TrackCall::GapDiff) { ++s.gapDiff; s.wGapDiff += w; }
        else { ++s.contra; s.wContra += w; }
    }
    return s;
}

void ClonalEngine::panelSupport(const std::vector<JunctionQuery>& qs, std::vector<TrackSupport>& out) const {
    out.assign(qs.size(), TrackSupport());
    struct Ref { uint32_t q; uint8_t side; };
    std::unordered_map<uint32_t, std::vector<Ref>> want;
    const size_t use = static_cast<size_t>(std::max(1, opt_.useMarkers)) * 2;   // spare markers for absent ones
    for (uint32_t q = 0; q < qs.size(); ++q) {
        for (uint8_t side = 0; side < 2; ++side) {
            const auto& h = side == 0 ? qs[q].left : qs[q].right;
            for (size_t i = 0; i < h.size() && i < use; ++i) want[h[i].id].push_back(Ref{q, side});
        }
    }
    if (want.empty()) return;
    struct Found { uint32_t key; uint32_t id; uint32_t at; };
    std::vector<Found> found;
    for (const LayoutTrack& t : model_.tracks()) {
        found.clear();
        for (uint32_t i = 0; i < t.oriented.size(); ++i) {
            auto it = want.find(t.oriented[i] >> 1);
            if (it == want.end()) continue;
            for (const Ref& r : it->second) found.push_back(Found{r.q * 2u + r.side, t.oriented[i] >> 1, i});
        }
        if (found.empty()) continue;
        std::sort(found.begin(), found.end(), [](const Found& a, const Found& b) {
            return a.key != b.key ? a.key < b.key : a.id < b.id;
        });
        const int64_t L = trackLength(t, model_.markerDenom());
        size_t a = 0;
        while (a < found.size()) {
            const uint32_t q = found[a].key / 2;
            size_t b = a;
            while (b < found.size() && found[b].key / 2 == q) ++b;
            // both sides present on this track?
            size_t mid = a;
            while (mid < b && found[mid].key % 2 == 0) ++mid;
            if (mid > a && mid < b) {
                auto finder = [&](size_t from, size_t to) {
                    return [&found, from, to](uint32_t id, uint32_t& i) {
                        for (size_t x = from; x < to; ++x) if (found[x].id == id) { i = found[x].at; return true; }
                        return false;
                    };
                };
                const SidePlace pl = placeWith(qs[q].left, true, L, t, finder(a, mid), 2);
                const SidePlace pr = placeWith(qs[q].right, false, L, t, finder(mid, b), 2);
                int64_t g;
                const TrackCall c = call(pl, pr, L, qs[q].asserted, qs[q].sized, opt_.maxAdjacent, opt_.maxOverlap, g);
                TrackSupport& s = out[q];
                ++s.tracks;
                if (c != TrackCall::Absent) {
                    ++s.placed;
                    if (c == TrackCall::Agree) ++s.agree;
                    else if (c == TrackCall::GapDiff) ++s.gapDiff;
                    else ++s.contra;
                }
            }
            a = b;
        }
    }
}

// ------------------------------------------------------------------------------ the stage

namespace {

struct EndCtx {
    std::vector<MarkerHit> exitHits, entryHits;
};

struct Proposal {
    int64_t partner = -1;
    int64_t gap = 0;
};

std::string revcompWindowHead(const std::string& s, size_t w) {
    return reverseComplement(s.substr(0, std::min(w, s.size())));
}
std::string revcompWindowTail(const std::string& s, size_t w) {
    return reverseComplement(s.size() > w ? s.substr(s.size() - w) : s);
}

}  // namespace

ClonalStats runClonalStage(const ClonalInputs& in, ClonalEngine& eng, std::vector<std::string>& seqs,
                           std::vector<double>& covs, std::vector<char>& layoutMembers, bool& changed) {
    ClonalStats st;
    changed = false;
    const ClonalOptions& o = eng.options();
    st.enabled = o.enabled;
    if (!o.enabled) return st;
    st.strict = o.strict;
    st.strictNr0 = o.strict && o.strictNr0;
    const double t0 = nowMs();
    const NearestSet& ns = eng.nearest();
    st.k = static_cast<int>(eng.relatives());
    st.dNear = ns.ok ? ns.dNear() : -1;
    st.group = ns.group();
    st.records = seqs.size();
    // side files (diagnostics for the dev rounds): <outDir>/om2_clonal/
    std::string sideDir;
    if (!in.outDir.empty() && util::makeDirs(in.outDir + "/om2_clonal")) sideDir = in.outDir + "/om2_clonal";
    if (!sideDir.empty()) {
        std::ofstream f(sideDir + "/nearest.tsv");
        f << "rank\taccession\tdist\tjaccard\tcontainment\tshared\ttrack_markers\tused\n";
        for (size_t i = 0; i < ns.ranked.size(); ++i) {
            const Relative& r = ns.ranked[i];
            f << i + 1 << '\t' << r.name << '\t' << fmtD(r.dist, "%.6f") << '\t' << fmtD(r.jaccard, "%.5f") << '\t'
              << fmtD(r.containment, "%.5f") << '\t' << r.shared << '\t' << r.trackMarkers << '\t'
              << (i < eng.relatives() ? 1 : 0) << '\n';
        }
    }
    auto done = [&](const char* why) { st.skip = why; st.ms = nowMs() - t0; st.recordsOut = seqs.size(); return st; };
    if (!ns.ok || eng.relatives() == 0) return done("no_relatives");
    if (!o.layout) return done("layout_off");
    if (ns.dNear() > o.dmax) return done("dmax");
    if (!in.seam || !in.seam->gating() || !in.seam->config().act() || !in.seam->evidence() ||
        in.seam->evidence()->abstained())
        return done("needs_seam_act");
    SeamContext& seam = *in.seam;
    const EvidenceIndex& ev = *seam.evidence();
    const PairIndex* pairs = seam.pairIndex();
    const UnitigGraph& g = ev.graph();
    const size_t W = static_cast<size_t>(o.window);
    const size_t kW = 25000;   // the window handed to the gate and to the pair ports
    st.ran = true;
    for (const JunctionAudit& X : seam.audit()) if (X.why.find("cap_walk_c1") != std::string::npos) ++st.capWalkC1;
    const uint32_t kminEff = static_cast<uint32_t>(std::max(1, std::min<int>(o.kmin, static_cast<int>(eng.relatives()))));

    // ---- 1. break panel-only joins the nearest relatives contradict ------------------------------
    //         (round 2: and check the size of every panel-sized join the relatives agree on)
    if (o.breakContra || o.resize) {
        eng.setAssembly(seqs);
        const std::vector<std::vector<RunOwner>> owners = seam.runOwners(seqs);
        const Ledger& led = seam.ledger();
        std::vector<std::vector<std::pair<size_t, size_t>>> cuts(seqs.size());
        std::vector<std::vector<std::pair<std::pair<size_t, size_t>, int64_t>>> resizes(seqs.size());   // round 2
        bool any = false;
        struct Kc2Cand { uint32_t id; std::string left, right; int64_t n; };   // round 3: KC2 class of every kept run
        std::vector<Kc2Cand> kc2c;
        for (size_t r = 0; r < seqs.size() && r < owners.size(); ++r) {
            for (const RunOwner& ow : owners[r]) {
                if (ow.id < 0 || static_cast<size_t>(ow.id) >= led.j.size()) continue;
                const Junction& J = led.j[static_cast<size_t>(ow.id)];
                if (J.writtenN <= 0 || J.source == Source::Resolver || J.source == Source::ResolverUnknown100 ||
                    J.source == Source::Clonal || J.source == Source::Wrap)
                    continue;
                const std::string& s = seqs[r];
                const size_t a = static_cast<size_t>(ow.start), b = static_cast<size_t>(ow.end);
                std::string left = s.substr(a > W ? a - W : 0, a > W ? W : a);
                std::string right = s.substr(b, std::min(W, s.size() - b));
                if (ow.orient == '-') { std::string t = reverseComplement(right); right = reverseComplement(left); left = t; }
                const JunctionQuery q = eng.query(left, right, static_cast<int64_t>(b - a), true);
                const TrackSupport ts = eng.relativeSupport(q);
                ++st.contraChecked;
                const bool contradicted = o.weighted
                    ? (ts.wContra >= o.wmin && ts.wAgree + ts.wGapDiff < 0.1 * ts.wContra)
                    : (ts.contra >= kminEff && ts.agree == 0 && ts.gapDiff == 0);
                const size_t nResBefore = resizes[r].size();
                if (!contradicted && o.resize && ts.placed >= 2 && ts.contra == 0) {
                    // Round 2: the clonal size check. Every placed relative reads the same order; they agree on one size
                    // S (tight); a walk class of this isolate has S. S == the asserted size: the size is fixed (C2 fills
                    // only at it). S != asserted by > max(500, 10%): the N-run is resized to S (dev E. coli
                    // GCF051549665v1: a 4,446-N layout gap at a true 2,146, all 5 relatives and the shortest walk at 2,146).
                    ++st.sizeChecked;
                    const int64_t asserted = static_cast<int64_t>(b - a);
                    const bool allAgree = ts.agree == ts.placed;
                    const bool allDiff = ts.gapDiff == ts.placed && !ts.gaps.empty();
                    bool nrun = true;
                    for (size_t x = a; x < b && nrun; ++x) nrun = s[x] == 'N' || s[x] == 'n';
                    if (allAgree) {
                        seam.clonalSize(static_cast<uint32_t>(ow.id), 0);
                        ++st.sizeAgree;
                    } else if (allDiff && nrun) {
                        std::vector<int64_t> gs = ts.gaps;
                        std::sort(gs.begin(), gs.end());
                        const int64_t S = gs[(gs.size() - 1) / 2];
                        const int64_t tol = std::max<int64_t>(500, std::llabs(S) / 10);
                        bool walk = false;
                        for (const WalkClass& w : J.walks) if (std::llabs(static_cast<int64_t>(w.len) - S) <= tol) walk = true;
                        if (J.walks.empty() && J.gmin > 0 && std::llabs(static_cast<int64_t>(J.gmin) - S) <= tol) walk = true;
                        // no walk at all (a SILENT junction: the graph does not connect the flanks, typically a coverage
                        // dropout): the relatives alone size it, as for a clonal join (dev Kp GCF055901135v2 j54: 2,232 N at a
                        // true 181, both placed relatives at 186)
                        const bool noWalk = J.walks.empty() && J.gmin <= 0 && J.verdict == Verdict::Silent;
                        // an overlap on every relative (S <= 0, e.g. the isolate's own -60 walk) is written as a 1-N butt
                        const int64_t newN = std::max<int64_t>(1, S);
                        if (S >= -o.maxOverlap && gs.back() - gs.front() <= tol && (walk || noWalk) &&
                            std::llabs(S - asserted) > tol && newN != asserted) {
                            resizes[r].push_back({{a, b}, newN});
                            seam.clonalSize(static_cast<uint32_t>(ow.id), static_cast<int32_t>(newN));
                            ClonalAnnot an;
                            an.set = true;
                            an.nrK = ts.tracks; an.nrPlaced = ts.placed; an.nrAgree = ts.agree; an.nrGapDiff = ts.gapDiff;
                            an.nrContra = ts.contra; an.nrGap = S;
                            an.wPlaced = ts.wPlaced; an.wAgree = ts.wAgree; an.wContra = ts.wContra;
                            an.dNear = ns.dNear(); an.group = ns.group();
                            an.action = "clonal_resize";
                            an.basis = "graph+nr";
                            eng.rows[static_cast<uint32_t>(ow.id)] = an;
                            ++st.resized;
                            any = true;
                        }
                    }
                }
                const bool resizedNow = resizes[r].size() > nResBefore;
                int64_t finalN = resizedNow ? resizes[r].back().second : static_cast<int64_t>(b - a);
                if (!contradicted && o.walksize && !resizedNow && !J.walks.empty()) {
                    // Round 3 (WALKSIZE): the panel-sized N-run is no walk of this isolate. The isolate's own graph sizes
                    // it: the walk class the relatives' median gap matches, else (the relatives place the ends, none at
                    // the asserted size) the isolated shortest walk class. Dev Kp round 2: 12 of 16 changed runs right
                    // (5,567 N at a true 572 = the only walk; 4,505 N at a true 572, classes 572:1 and 28,273:12).
                    ++st.walksizeChecked;
                    const int64_t asserted = static_cast<int64_t>(b - a);
                    bool nrun = true;
                    for (size_t x = a; x < b && nrun; ++x) nrun = s[x] == 'N' || s[x] == 'n';
                    auto tolW = [](int64_t x) { return std::max<int64_t>(500, std::llabs(x) / 10); };
                    bool assertedIsWalk = false;
                    for (const WalkClass& w : J.walks)
                        if (std::llabs(static_cast<int64_t>(w.len) - asserted) <= tolW(w.len)) assertedIsWalk = true;
                    int64_t pick = INT64_MIN;
                    const char* basis = "graph";
                    if (nrun && !assertedIsWalk) {
                        if (ts.placed >= 2 && ts.contra == 0 && !ts.gaps.empty()) {
                            std::vector<int64_t> gs = ts.gaps;
                            std::sort(gs.begin(), gs.end());
                            const int64_t S = gs[(gs.size() - 1) / 2];
                            int64_t bestD = INT64_MAX;
                            for (const WalkClass& w : J.walks) {
                                const int64_t dd = std::llabs(static_cast<int64_t>(w.len) - S);
                                if (dd <= tolW(S) && dd < bestD) { bestD = dd; pick = w.len; basis = "graph+nr"; }
                            }
                        }
                        // round 3b: any placed relative (the relatives confirm the order; they may carry the element
                        // this isolate lacks, so their size is not asked), unless every placed relative (>= 2) reads the
                        // asserted size; the shortest class isolated by >= 10 kb (dev Kp round 2: 10 of 10 changed runs
                        // right on Rule-Q isolates, against 3 with the round-3 'no relative at the asserted size' rule)
                        const bool allAtAsserted = ts.placed >= 2 && ts.agree == ts.placed;
                        if (pick == INT64_MIN && ts.placed >= 1 && !allAtAsserted) {
                            std::vector<int32_t> lens;
                            std::vector<uint16_t> cnts;
                            for (const WalkClass& w : J.walks) { lens.push_back(w.len); cnts.push_back(w.n); }
                            const int32_t c1 = isolatedShortestWalk(lens, cnts, J.exhaustive, o.c1GapLayout, o.c1Max);
                            if (c1 != INT32_MIN) { pick = c1; basis = "graph_c1"; }
                        }
                    }
                    if (pick != INT64_MIN && pick >= -o.maxOverlap) {
                        const int64_t newN = std::max<int64_t>(1, pick);
                        if (newN != asserted && std::llabs(pick - asserted) > tolW(pick)) {
                            resizes[r].push_back({{a, b}, newN});
                            seam.clonalSize(static_cast<uint32_t>(ow.id), static_cast<int32_t>(newN));
                            ClonalAnnot an;
                            an.set = true;
                            an.nrK = ts.tracks; an.nrPlaced = ts.placed; an.nrAgree = ts.agree; an.nrGapDiff = ts.gapDiff;
                            an.nrContra = ts.contra; an.nrGap = ts.gaps.empty() ? INT64_MIN : median(ts.gaps);
                            an.wPlaced = ts.wPlaced; an.wAgree = ts.wAgree; an.wContra = ts.wContra;
                            an.dNear = ns.dNear(); an.group = ns.group();
                            an.action = "walk_resize";
                            an.basis = basis;
                            eng.rows[static_cast<uint32_t>(ow.id)] = an;
                            ++st.walksizeResized;
                            any = true;
                            finalN = newN;
                        }
                    }
                }
                bool unverifiedCut = false;
                const bool layoutSrc = J.source == Source::Layout || J.source == Source::LayoutCap2000;
                const bool walkless = J.walks.empty() && J.gmin <= 0;
                bool assertedWalk = false;
                for (const WalkClass& w : J.walks)
                    if (std::llabs(static_cast<int64_t>(w.len) - finalN) <= std::max<int64_t>(500, std::llabs(w.len) / 10))
                        assertedWalk = true;
                if (!contradicted && o.unverified && !resizedNow && finalN == static_cast<int64_t>(b - a) &&
                    (J.tier == Tier::D || J.tier == Tier::E) &&
                    ((walkless && layoutSrc) || (!walkless && !assertedWalk && (layoutSrc || J.source == Source::Join)))) {
                    // Round 3b (UNVERIFIED): a panel-sized gap whose size no walk of this isolate has (or, for a layout gap,
                    // with no walk at all) and that the relatives do not all read: the size is unverifiable and a wrong
                    // size is a misassembly. The run is broken. Dev om2r1 (all 7 species, Rule-Q, genome view): join gaps
                    // that are no walk class 273 of 396 wrong, layout gaps without a walk 192 of 354 wrong.
                    const bool agreed = ts.placed >= 2 && ts.agree == ts.placed && ts.contra == 0;
                    if (!agreed) {
                        cuts[r].emplace_back(a, b);
                        any = true;
                        unverifiedCut = true;
                        ClonalAnnot an;
                        an.set = true;
                        an.nrK = ts.tracks; an.nrPlaced = ts.placed; an.nrAgree = ts.agree; an.nrGapDiff = ts.gapDiff;
                        an.nrContra = ts.contra; an.nrGap = ts.gaps.empty() ? INT64_MIN : median(ts.gaps);
                        an.wPlaced = ts.wPlaced; an.wAgree = ts.wAgree; an.wContra = ts.wContra;
                        an.dNear = ns.dNear(); an.group = ns.group();
                        an.action = "break_unverified";
                        an.basis = "nr";
                        eng.rows[static_cast<uint32_t>(ow.id)] = an;
                        seam.clonalBreak(static_cast<uint32_t>(ow.id), "clonal_unverified");
                        ++st.unverifiedBroken;
                    }
                }
                bool strictCut = false;
                if (o.strict && !unverifiedCut && !(contradicted && !(J.tier == Tier::A || J.tier == Tier::B || J.tier == Tier::C)) &&
                    (J.tier == Tier::C || J.tier == Tier::D || J.tier == Tier::E)) {
                    // Round 3c (STRICT): the run is written only when this isolate's graph or every placed relative sizes
                    // it at its final size (strictVerdict). A resized run is re-read by the relatives at its new size.
                    ++st.strictChecked;
                    StrictInput si;
                    si.n = finalN;
                    si.exhaustive = J.exhaustive;
                    for (const WalkClass& w : J.walks) { si.walkLens.push_back(w.len); si.walkCounts.push_back(w.n); }
                    TrackSupport tf = ts;
                    if (finalN != static_cast<int64_t>(b - a)) tf = eng.relativeSupport(eng.query(left, right, finalN, true));
                    si.placed = tf.placed; si.agree = tf.agree; si.contra = tf.contra;
                    si.wAgree = tf.wAgree; si.wContra = tf.wContra;
                    if (const char* why = strictVerdict(si, o)) {
                        cuts[r].emplace_back(a, b);
                        any = true;
                        strictCut = true;
                        auto it = eng.rows.find(static_cast<uint32_t>(ow.id));
                        ClonalAnnot an = it != eng.rows.end() ? it->second : ClonalAnnot{};
                        an.set = true;
                        an.nrK = tf.tracks; an.nrPlaced = tf.placed; an.nrAgree = tf.agree; an.nrGapDiff = tf.gapDiff;
                        an.nrContra = tf.contra; an.nrGap = tf.gaps.empty() ? INT64_MIN : median(tf.gaps);
                        an.wPlaced = tf.wPlaced; an.wAgree = tf.wAgree; an.wContra = tf.wContra;
                        an.dNear = ns.dNear(); an.group = ns.group();
                        an.action = an.action == "." || an.action.empty() ? std::string("break_") + why
                                                                           : an.action + "+break_" + why;
                        an.basis = "strict";
                        eng.rows[static_cast<uint32_t>(ow.id)] = an;
                        seam.clonalBreak(static_cast<uint32_t>(ow.id), std::string("clonal_") + why);
                        ++st.strictBroken;
                        if (std::strcmp(why, "strict_periodic") == 0) ++st.strictPeriodic;
                        if (std::strcmp(why, "strict_no_relative") == 0) ++st.strictNoRelative;
                    }
                }
                const bool willCut = unverifiedCut || strictCut || (contradicted && !(J.tier == Tier::A || J.tier == Tier::B || J.tier == Tier::C));
                if (o.kc2 && !willCut) kc2c.push_back(Kc2Cand{static_cast<uint32_t>(ow.id), left, right, finalN});
                if (unverifiedCut || strictCut) continue;
                if (!contradicted) continue;
                if (J.tier == Tier::A || J.tier == Tier::B || J.tier == Tier::C) { ++st.contraKeptProven; continue; }
                cuts[r].emplace_back(a, b);
                any = true;
                ClonalAnnot an;
                an.set = true;
                an.nrK = ts.tracks; an.nrPlaced = ts.placed; an.nrAgree = ts.agree; an.nrGapDiff = ts.gapDiff;
                an.nrContra = ts.contra; an.nrGap = median(ts.gaps);
                an.wPlaced = ts.wPlaced; an.wAgree = ts.wAgree; an.wContra = ts.wContra;
                an.dNear = ns.dNear(); an.group = ns.group();
                an.action = "break_contra";
                an.basis = "nr";
                eng.rows[static_cast<uint32_t>(ow.id)] = an;
                seam.clonalBreak(static_cast<uint32_t>(ow.id), "clonal_contra");
                ++st.contraBroken;
            }
        }
        if (o.kc2 && !kc2c.empty()) {
            // Round 3 (KC2): the class of every kept panel-sized run at its final size, as the writer labels it: the
            // relatives' gap recurring at >= 90% of the panel is a positional locus; otherwise the run is non-positional
            // and, unless the isolate's shortest walk between its flanks is repeat-free, a labelled gap C2 never fills.
            std::vector<JunctionQuery> qs;
            qs.reserve(kc2c.size());
            for (const Kc2Cand& c : kc2c) qs.push_back(eng.query(c.left, c.right, c.n, true));
            std::vector<TrackSupport> panel;
            eng.panelSupport(qs, panel);
            for (size_t i = 0; i < kc2c.size() && i < panel.size(); ++i) {
                const Kc2Cand& c = kc2c[i];
                ++st.kc2Checked;
                const Junction& J = led.j[c.id];
                const bool panelKnown = panel[i].placed >= o.panelMin;
                const bool positional = panelKnown && panel[i].agree >= o.positional * panel[i].placed;
                const bool repeatEnds = J.endA == EndClass::Repeat || J.endB == EndClass::Repeat;
                int element = -1;
                const EndAnchor A = ev.anchorExit(c.left);
                const EndAnchor B = ev.anchorEntry(c.right);
                if (A.ok && B.ok) {
                    const WalkSet ws = ev.walks(A, B);
                    element = ws.reachable ? ((ws.innerRepeatBp > 0 || ws.innerMaxCopy > 1.5f) ? 1 : 0) : -1;
                }
                const int cls = positional ? ((repeatEnds || element == 1) ? kPositional : kLayout) : kNonPositional;
                auto it = eng.rows.find(c.id);
                if (it != eng.rows.end()) { it->second.cls = cls; if (it->second.element < 0) it->second.element = element; }
                if (cls != kNonPositional || element == 0 || J.noFill) continue;
                seam.clonalNoFill(c.id, "kc2_nonpositional_gap");
                ++st.kc2NoFill;
                if (it == eng.rows.end()) {
                    ClonalAnnot an;
                    an.set = true;
                    an.cls = kNonPositional;
                    an.element = element;
                    an.dNear = ns.dNear();
                    an.group = ns.group();
                    an.panelPlaced = panel[i].placed;
                    an.panelAgree = panel[i].agree;
                    an.action = "kc2_gap";
                    an.basis = "ledger";
                    eng.rows[c.id] = an;
                } else {
                    it->second.action += "+kc2_gap";
                }
            }
        }
        if (any) {
            std::vector<std::string> out;
            std::vector<double> outCov;
            std::vector<char> outLm;
            const bool hadLm = !layoutMembers.empty();
            for (size_t r = 0; r < seqs.size(); ++r) {
                const double cv = r < covs.size() ? covs[r] : 0.0;
                const char lm = r < layoutMembers.size() ? layoutMembers[r] : 0;
                // round 2: resized N-runs first (a run is either cut or resized, never both)
                if (!resizes[r].empty()) {
                    std::sort(resizes[r].begin(), resizes[r].end());
                    std::string t;
                    size_t from = 0;
                    for (const auto& z : resizes[r]) {
                        t += seqs[r].substr(from, z.first.first - from);
                        t.append(static_cast<size_t>(z.second), 'N');
                        from = z.first.second;
                    }
                    t += seqs[r].substr(from);
                    // re-map the cut runs of this record onto the resized text
                    std::vector<std::pair<size_t, size_t>> nc;
                    for (const auto& c : cuts[r]) {
                        int64_t d = 0;
                        for (const auto& z : resizes[r])
                            if (z.first.second <= c.first) d += z.second - static_cast<int64_t>(z.first.second - z.first.first);
                        nc.emplace_back(static_cast<size_t>(static_cast<int64_t>(c.first) + d), static_cast<size_t>(static_cast<int64_t>(c.second) + d));
                    }
                    cuts[r].swap(nc);
                    seqs[r].swap(t);
                }
                std::sort(cuts[r].begin(), cuts[r].end());
                size_t from = 0;
                for (const auto& c : cuts[r]) {
                    if (c.first > from) { out.push_back(seqs[r].substr(from, c.first - from)); outCov.push_back(cv); outLm.push_back(lm); }
                    from = c.second;
                }
                if (from < seqs[r].size()) { out.push_back(seqs[r].substr(from)); outCov.push_back(cv); outLm.push_back(lm); }
            }
            seqs.swap(out);
            covs.swap(outCov);
            if (hadLm) layoutMembers.swap(outLm);
            changed = true;
        }
    }

    // ---- 2. candidate records: chromosomal, long enough, placed ends -----------------------------
    eng.setAssembly(seqs);
    const OrganismModel& model = eng.model();
    const size_t R = seqs.size();
    std::vector<char> cand(R, 0);
    for (size_t r = 0; r < R; ++r) {
        if (seqs[r].size() < 500) continue;
        uint32_t chrV = 0, plsV = 0;
        forEachMarkerKmer(seqs[r], [&](uint64_t km, uint32_t, int) {
            const uint32_t id = model.markerOf(km);
            if (id == UINT32_MAX) return;
            const uint32_t gc = model.markerGenomes(id, Replicon::Chromosome);
            const uint32_t gp = model.markerGenomes(id, Replicon::Plasmid);
            if (gc > 0 && gp == 0) ++chrV;
            else if (gp > 0 && gc == 0) ++plsV;
        }, model.markerDenom());
        if ((plsV > 0 && plsV >= chrV * 3) || (r < covs.size() && seam.plasmidByDepth(covs[r]))) {
            cand[r] = 2;   // plasmid: chained only on the relatives' plasmids (4b)
            continue;
        }
        cand[r] = 1;
        ++st.candidates;
    }
    const size_t K = eng.relatives();
    const size_t need = static_cast<size_t>(std::max(1, o.minMarkers));
    // end u = 2r (head) / 2r+1 (tail); exit via u leaves the record through u, entry via u enters it
    std::vector<EndCtx> ends(2 * R);
    std::vector<std::vector<SidePlace>> exitPl(2 * R, std::vector<SidePlace>(K)), entryPl(2 * R, std::vector<SidePlace>(K));
    auto endHits = [&](size_t r) {
        const std::string& s = seqs[r];
        ends[2 * r + 1].exitHits = eng.hits(s, true);                                  // tail, record forward
        ends[2 * r].exitHits = eng.hits(revcompWindowHead(s, W), true);                 // head, record reversed
        ends[2 * r].entryHits = eng.hits(s, false);                                     // head, record forward
        ends[2 * r + 1].entryHits = eng.hits(revcompWindowTail(s, W), false);           // tail, record reversed
    };
    for (size_t r = 0; r < R; ++r) {
        if (cand[r] != 1) continue;
        const std::string& s = seqs[r];
        ends[2 * r + 1].exitHits = eng.hits(s, true);                                  // tail, record forward
        ends[2 * r].exitHits = eng.hits(revcompWindowHead(s, W), true);                 // head, record reversed
        ends[2 * r].entryHits = eng.hits(s, false);                                     // head, record forward
        ends[2 * r + 1].entryHits = eng.hits(revcompWindowTail(s, W), false);           // tail, record reversed
        for (size_t t = 0; t < K; ++t) {
            for (size_t u = 2 * r; u <= 2 * r + 1; ++u) {
                exitPl[u][t] = eng.placeOnRelative(ends[u].exitHits, true, 0, t, need);
                entryPl[u][t] = eng.placeOnRelative(ends[u].entryHits, false, 0, t, need);
            }
            // The record must read collinearly on the relative: its two ends the record's length apart.
            const SidePlace& h = entryPl[2 * r][t];
            const SidePlace& tl = exitPl[2 * r + 1][t];
            if (h.ok && tl.ok) {
                const int64_t Lt = eng.relativeLength(t);
                int64_t span = (tl.x - h.x) * h.dir;
                if (Lt > 0) { while (span < 0) span += Lt; while (span > Lt) span -= Lt; }
                const int64_t len = static_cast<int64_t>(s.size());
                const int64_t slack = std::max<int64_t>(10000, len / 20) + (o.indel ? o.indelMax : 0);
                if (h.dir != tl.dir || std::llabs(span - len) > slack) {
                    for (size_t u = 2 * r; u <= 2 * r + 1; ++u) { exitPl[u][t] = SidePlace(); entryPl[u][t] = SidePlace(); }
                }
            }
        }
        for (size_t u = 2 * r; u <= 2 * r + 1; ++u)
            for (size_t t = 0; t < K; ++t) if (exitPl[u][t].ok) { ++st.endsPlaced; break; }
    }

    // ---- 3. per relative: the next entry beyond every exit ---------------------------------------
    std::vector<std::vector<Proposal>> prop(2 * R, std::vector<Proposal>(K));
    for (size_t t = 0; t < K; ++t) {
        const int64_t Lt = eng.relativeLength(t);
        auto fwd = [&](int64_t from, int64_t to, int8_t dir) {
            int64_t gg = (to - from) * dir;
            if (Lt > 0) { while (gg < -Lt / 2) gg += Lt; while (gg > Lt / 2) gg -= Lt; }
            return gg;
        };
        for (size_t u = 0; u < 2 * R; ++u) {
            const SidePlace& ex = exitPl[u][t];
            if (!ex.ok) continue;
            int64_t best = INT64_MAX;
            int64_t bestV = -1;
            for (size_t v = 0; v < 2 * R; ++v) {
                if (v / 2 == u / 2) continue;
                const SidePlace& en = entryPl[v][t];
                if (!en.ok || en.dir != ex.dir) continue;
                const int64_t gg = fwd(ex.x, en.x, ex.dir);
                if (gg < -o.maxOverlap || gg > o.maxAdjacent) continue;
                if (gg < best || (gg == best && static_cast<int64_t>(v) < bestV)) { best = gg; bestV = static_cast<int64_t>(v); }
            }
            if (bestV < 0) continue;
            // Another record's exit inside the gap: that record overlaps this end on the relative.
            bool ambiguous = false;
            for (size_t w = 0; w < 2 * R && !ambiguous; ++w) {
                if (w / 2 == u / 2 || w / 2 == static_cast<size_t>(bestV) / 2) continue;
                const SidePlace& e2 = exitPl[w][t];
                if (!e2.ok || e2.dir != ex.dir) continue;
                const int64_t gw = fwd(ex.x, e2.x, ex.dir);
                if (gw >= -o.maxOverlap && gw < best) ambiguous = true;
            }
            if (ambiguous) continue;
            prop[u][t] = Proposal{bestV, best};
        }
    }

    if (!sideDir.empty()) {
        std::ofstream f(sideDir + "/ends.tsv");
        f << "record\tlength\tcandidate\tend\trelative\texit_x\texit_dir\texit_n\tentry_x\tentry_dir\tentry_n\t"
             "partner\tgap\n";
        for (size_t u = 0; u < 2 * R; ++u) {
            if (cand[u / 2] != 1) continue;
            for (size_t t = 0; t < K; ++t) {
                const SidePlace& a = exitPl[u][t];
                const SidePlace& b = entryPl[u][t];
                f << u / 2 << '\t' << seqs[u / 2].size() << "\tchr\t" << ((u & 1) ? "tail" : "head") << '\t'
                  << eng.relativeName(t) << '\t';
                if (a.ok) f << a.x << '\t' << int(a.dir) << '\t' << a.n; else f << ".\t.\t0";
                f << '\t';
                if (b.ok) f << b.x << '\t' << int(b.dir) << '\t' << b.n; else f << ".\t.\t0";
                const Proposal& pr = prop[u][t];
                if (pr.partner >= 0) f << '\t' << pr.partner / 2 << ((pr.partner & 1) ? ":tail" : ":head") << '\t' << pr.gap;
                else f << "\t.\t.";
                f << '\n';
            }
        }
    }

    // ---- 4. consensus over the relatives, mutual pairs --------------------------------------------
    struct Cons {
        int64_t partner = -1;
        uint32_t placed = 0, agree = 0;
        double wPlaced = 0, wAgree = 0, wContra = 0;
        int64_t gap = 0, spread = 0;
    };
    // The relatives' vote on one end: each relative on which the end is placed votes for its proposed partner with
    // weight eng.weight(rank); the partner with the largest weight (count, under the count rule) wins. The gap is the
    // median of the closest agreeing relatives (weight >= half the top agreeing weight), and so is the spread.
    auto vote = [&](const std::function<bool(size_t)>& placedOn, const std::function<Proposal(size_t)>& propOn) {
        Cons c;
        struct Acc { uint32_t n = 0; double w = 0; std::vector<std::pair<double, int64_t>> g; };
        std::map<int64_t, Acc> by;
        for (size_t t = 0; t < K; ++t) {
            if (!placedOn(t)) continue;
            const double w = eng.weight(t);
            ++c.placed;
            c.wPlaced += w;
            const Proposal p = propOn(t);
            if (p.partner < 0) continue;
            Acc& a = by[p.partner];
            ++a.n;
            a.w += w;
            a.g.emplace_back(w, p.gap);
        }
        double best = -1;
        for (const auto& kv : by) {
            const double key = o.weighted ? kv.second.w : static_cast<double>(kv.second.n);
            if (key > best + 1e-12) { best = key; c.partner = kv.first; }
        }
        if (c.partner >= 0) {
            const Acc& a = by[c.partner];
            c.agree = a.n;
            c.wAgree = a.w;
            for (const auto& kv : by) if (kv.first != c.partner) c.wContra += kv.second.w;
            double wmax = 0;
            for (const auto& x : a.g) wmax = std::max(wmax, x.first);
            std::vector<int64_t> gs;
            for (const auto& x : a.g) if (x.first >= 0.5 * wmax) gs.push_back(x.second);
            c.gap = median(gs);
            c.spread = *std::max_element(gs.begin(), gs.end()) - *std::min_element(gs.begin(), gs.end());
        }
        return c;
    };
    auto accepted = [&](const Cons& c, uint32_t kmin) {
        if (c.partner < 0) return false;
        if (o.weighted) return c.wAgree >= o.wmin && c.wAgree >= o.agreeShare * (c.wAgree + c.wContra);
        return c.agree >= kmin && c.agree >= o.agreeShare * c.placed;
    };
    std::vector<Cons> cons(2 * R);
    for (size_t u = 0; u < 2 * R; ++u) {
        Cons c = vote([&](size_t t) { return exitPl[u][t].ok; }, [&](size_t t) { return prop[u][t]; });
        if (accepted(c, kminEff)) ++st.proposals;
        else c.partner = -1;
        cons[u] = c;
    }
    struct Link {
        size_t u = 0, v = 0;
        bool plasmid = false;
        Cons c;
        LinkEvidence e;
        LinkDecision d;
        std::string L, Rs;   // gate windows
    };
    std::vector<Link> links;
    for (size_t u = 0; u < 2 * R; ++u) {
        const int64_t v = cons[u].partner;
        if (v < 0 || static_cast<size_t>(v) <= u) continue;
        if (cons[static_cast<size_t>(v)].partner != static_cast<int64_t>(u)) continue;
        Link lk;
        lk.u = u;
        lk.v = static_cast<size_t>(v);
        lk.c = cons[u];
        // the two directions measure the same gap; keep the more conservative spread
        lk.c.spread = std::max(cons[u].spread, cons[static_cast<size_t>(v)].spread);
        links.push_back(std::move(lk));
    }
    st.mutual = links.size();

    // ---- 4b. plasmids: the same chaining on the nearest relatives' own plasmids (sidecar) -----------
    const size_t PT = eng.plasmidTracks();
    st.plasmidTracks = PT;
    std::vector<char> pcand(R, 0), pplaced(R, 0);
    if (PT > 0) {
        for (size_t r = 0; r < R; ++r) {
            // plasmid-voted or plasmid-deep records only: a chromosome record that no relative's chromosome places
            // is not thereby a plasmid (smoke4: one such record joined on a relative's plasmid was FALSE)
            if (cand[r] == 2) pcand[r] = 1;
            if (pcand[r]) ++st.plasmidCandidates;
        }
        const size_t needP = 2;
        struct PP { int64_t pi = -1; SidePlace s; bool tie = false; };
        std::vector<std::vector<PP>> pex(2 * R, std::vector<PP>(K)), pen(2 * R, std::vector<PP>(K));
        auto better = [](PP& cur, int64_t pi, const SidePlace& sp) {
            if (!sp.ok) return;
            if (cur.pi < 0 || sp.n > cur.s.n) { cur.pi = pi; cur.s = sp; cur.tie = false; }
            else if (sp.n == cur.s.n) cur.tie = true;
        };
        for (size_t r = 0; r < R; ++r) {
            if (!pcand[r]) continue;
            if (cand[r] != 1) endHits(r);
            for (size_t pi = 0; pi < PT; ++pi) {
                const size_t t = eng.plasmidRank(pi);
                if (t >= K) continue;
                for (size_t u = 2 * r; u <= 2 * r + 1; ++u) {
                    better(pex[u][t], static_cast<int64_t>(pi), eng.placeOnPlasmid(ends[u].exitHits, true, pi, needP));
                    better(pen[u][t], static_cast<int64_t>(pi), eng.placeOnPlasmid(ends[u].entryHits, false, pi, needP));
                }
            }
            for (size_t t = 0; t < K; ++t) {
                for (size_t u = 2 * r; u <= 2 * r + 1; ++u) {
                    if (pex[u][t].tie) pex[u][t].pi = -1;   // two of the relative's plasmids place it equally
                    if (pen[u][t].tie) pen[u][t].pi = -1;
                    if (pex[u][t].pi >= 0 || pen[u][t].pi >= 0) pplaced[r] = 1;
                }
                // collinear on the plasmid: the record's two ends its length apart
                const PP& h = pen[2 * r][t];
                const PP& tl = pex[2 * r + 1][t];
                if (h.pi >= 0 && tl.pi >= 0) {
                    bool bad = h.pi != tl.pi || h.s.dir != tl.s.dir;
                    if (!bad) {
                        const int64_t Lp = eng.plasmidLength(static_cast<size_t>(h.pi));
                        int64_t span = (tl.s.x - h.s.x) * h.s.dir;
                        if (Lp > 0) { while (span < 0) span += Lp; while (span > Lp) span -= Lp; }
                        const int64_t len = static_cast<int64_t>(seqs[r].size());
                        bad = std::llabs(span - len) > std::max<int64_t>(3000, len / 20);
                    }
                    if (bad) for (size_t u = 2 * r; u <= 2 * r + 1; ++u) { pex[u][t].pi = -1; pen[u][t].pi = -1; }
                }
            }
        }
        std::vector<std::vector<Proposal>> pprop(2 * R, std::vector<Proposal>(K));
        for (size_t t = 0; t < K; ++t) {
            for (size_t u = 0; u < 2 * R; ++u) {
                const PP& ex = pex[u][t];
                if (ex.pi < 0) continue;
                const int64_t Lp = eng.plasmidLength(static_cast<size_t>(ex.pi));
                auto fwd = [&](int64_t from, int64_t to, int8_t dir) {
                    int64_t gg = (to - from) * dir;
                    if (Lp > 0) { while (gg < -o.maxOverlap) gg += Lp; while (gg > Lp) gg -= Lp; }
                    return gg;
                };
                int64_t best = INT64_MAX, bestV = -1;
                for (size_t v = 0; v < 2 * R; ++v) {
                    if (v / 2 == u / 2) continue;
                    const PP& en = pen[v][t];
                    if (en.pi != ex.pi || en.s.dir != ex.s.dir) continue;
                    const int64_t gg = fwd(ex.s.x, en.s.x, ex.s.dir);
                    if (gg < -o.maxOverlap || gg > std::min<int64_t>(Lp, o.maxAdjacent)) continue;
                    if (gg < best || (gg == best && static_cast<int64_t>(v) < bestV)) { best = gg; bestV = static_cast<int64_t>(v); }
                }
                if (bestV < 0) continue;
                bool ambiguous = false;
                for (size_t w = 0; w < 2 * R && !ambiguous; ++w) {
                    if (w / 2 == u / 2 || w / 2 == static_cast<size_t>(bestV) / 2) continue;
                    const PP& e2 = pex[w][t];
                    if (e2.pi != ex.pi || e2.s.dir != ex.s.dir) continue;
                    const int64_t gw = fwd(ex.s.x, e2.s.x, ex.s.dir);
                    if (gw >= -o.maxOverlap && gw < best) ambiguous = true;
                }
                if (!ambiguous) pprop[u][t] = Proposal{bestV, best};
            }
        }
        const uint32_t kminP = static_cast<uint32_t>(std::max(1, o.kminPlasmid));
        std::vector<Cons> pcons(2 * R);
        for (size_t u = 0; u < 2 * R; ++u) {
            Cons c = vote([&](size_t t) { return pex[u][t].pi >= 0; }, [&](size_t t) { return pprop[u][t]; });
            // Plasmids are mosaic: one relative's plasmid is never enough, whatever its weight (smoke4: 2 of 8
            // single-relative plasmid joins through an IS walk were FALSE). >= k_min_plasmid relatives must agree.
            if (!accepted(c, kminP) || c.agree < kminP) c.partner = -1;
            pcons[u] = c;
        }
        for (size_t u = 0; u < 2 * R; ++u) {
            const int64_t v = pcons[u].partner;
            if (v < 0 || static_cast<size_t>(v) <= u || pcons[static_cast<size_t>(v)].partner != static_cast<int64_t>(u)) continue;
            // an end already linked on the chromosome relatives keeps that link
            bool taken = false;
            for (const Link& lk : links) if (lk.u == u || lk.v == u || lk.u == static_cast<size_t>(v) || lk.v == static_cast<size_t>(v)) taken = true;
            if (taken) continue;
            Link lk;
            lk.u = u;
            lk.v = static_cast<size_t>(v);
            lk.plasmid = true;
            lk.c = pcons[u];
            lk.c.spread = std::max(pcons[u].spread, pcons[static_cast<size_t>(v)].spread);
            links.push_back(std::move(lk));
            ++st.plasmidMutual;
        }
    }

    // ---- 5. evidence for each mutual link -------------------------------------------------------
    auto orientedExit = [&](size_t u) {   // the record oriented so end u is its right end
        const std::string& s = seqs[u / 2];
        return (u & 1) ? (s.size() > kW ? s.substr(s.size() - kW) : s) : revcompWindowHead(s, kW);
    };
    auto orientedEntry = [&](size_t v) {  // the record oriented so end v is its left end
        const std::string& s = seqs[v / 2];
        return (v & 1) ? revcompWindowTail(s, kW) : s.substr(0, std::min(kW, s.size()));
    };
    // origin port of every candidate record end (the pair index is keyed on the origin pieces): chromosome candidates
    // for chromosome links; plasmid records placed on a relative's plasmid for plasmid links. Repeat contigs (collapsed
    // rRNA / IS copies, deep like a plasmid) are neither: a flank's pairs into its own repeat are no contradiction.
    std::unordered_map<uint32_t, size_t> endOfPortChr, endOfPortPls;
    if (pairs && pairs->usable()) {
        for (size_t r = 0; r < R; ++r) {
            const bool c = cand[r] == 1, pp = pcand[r] && pplaced[r];
            if (!c && !pp) continue;
            for (size_t u = 2 * r; u <= 2 * r + 1; ++u) {
                const uint32_t port = pairs->exitPort(orientedExit(u));
                if (port == UINT32_MAX) continue;
                if (c) endOfPortChr.emplace(port, u);
                if (pp) endOfPortPls.emplace(port, u);
            }
        }
    }
    std::vector<JunctionQuery> pq;
    for (Link& lk : links) {
        lk.L = orientedExit(lk.u);
        lk.Rs = orientedEntry(lk.v);
        LinkEvidence& e = lk.e;
        e.k = lk.plasmid ? static_cast<uint32_t>(std::max(1, o.kminPlasmid)) : static_cast<uint32_t>(K);
        e.plasmid = lk.plasmid;
        e.placed = lk.c.placed;
        e.agree = lk.c.agree;
        e.wAgree = lk.c.wAgree;
        e.wContra = lk.c.wContra;
        e.nrGap = lk.c.gap;
        e.nrSpread = lk.c.spread;
        e.wPlaced = lk.c.wPlaced;
        if (o.overlapMerge && lk.c.gap < 0)
            e.overlapExact = exactEndOverlap(lk.L, lk.Rs, -lk.c.gap, o.minOverlap, o.maxOverlap);
        const EndAnchor A = ev.anchorExit(lk.L);
        const EndAnchor B = ev.anchorEntry(lk.Rs);
        e.anchoredA = A.ok;
        e.anchoredB = B.ok;
        if (A.ok) {
            e.edgeA = !g.exits(A.node.id, A.node.o).empty();
            e.repeatA = !ev.mult(A.node.id).hardSingle;
        }
        if (B.ok) {
            const ONode bb = flip(B.node);
            e.edgeB = !g.exits(bb.id, bb.o).empty();
            e.repeatB = !ev.mult(B.node.id).hardSingle;
        }
        if (A.ok && B.ok) {
            const WalkSet ws = ev.walks(A, B);
            e.reachable = ws.reachable;
            e.exhaustive = ws.exhaustive;
            for (const WalkClass& wc : ws.classes) { e.walkLens.push_back(wc.len); e.walkCounts.push_back(wc.n); }
            if (e.walkLens.empty() && ws.reachable) { e.walkLens.push_back(ws.gmin); e.walkCounts.clear(); }
            e.element = ws.reachable && (ws.innerRepeatBp > 0 || ws.innerMaxCopy > 1.5f);
            if (!ws.reachable) {
                const FirstUnique uL = ev.firstUnique(A.node);
                const FirstUnique uR = ev.firstUnique(flip(B.node));
                uint32_t other = 0;
                for (uint32_t id : uL.ids) if (id != B.node.id) ++other;
                for (uint32_t id : uR.ids) if (id != A.node.id) ++other;
                e.contraGraph = other > 0;
            }
        }
        // read pairs from either end to a THIRD candidate record's end. Pairs into repeat contigs, plasmid
        // pieces or short unplaced pieces are not contradictions: an element this isolate carries inside
        // the gap (a novel insertion) is exactly such a piece.
        if (pairs && pairs->usable()) {
            const uint32_t P = pairs->exitPort(lk.L), Q = pairs->entryPort(lk.Rs);
            const std::unordered_map<uint32_t, size_t>& endOfPort = lk.plasmid ? endOfPortPls : endOfPortChr;
            auto contra = [&](uint32_t port, uint32_t other, size_t selfEnd, size_t otherEnd) -> uint32_t {
                if (port == UINT32_MAX) return 0;
                uint32_t n = 0, toOther = 0;
                for (const auto& pr : pairs->partners(port)) {
                    if (other != UINT32_MAX && pr.first == other) { toOther += pr.second; continue; }
                    auto it = endOfPort.find(pr.first);
                    if (it == endOfPort.end()) continue;
                    if (it->second / 2 == selfEnd / 2 || it->second == otherEnd) continue;
                    n += pr.second;
                }
                return n > toOther ? n : 0;
            };
            e.pairContra = std::max(contra(P, Q, lk.u, lk.v), contra(Q, P, lk.v, lk.u));
        }
        pq.push_back(eng.query(lk.L, lk.Rs, lk.c.gap, true));
    }
    std::vector<TrackSupport> panel;
    eng.panelSupport(pq, panel);
    for (size_t i = 0; i < links.size(); ++i) {
        links[i].e.panelPlaced = panel[i].placed;
        links[i].e.panelAgree = panel[i].agree;
        links[i].d = decideLink(links[i].e, o);
        if (o.strict && links[i].d.join && !links[i].plasmid) {
            // Round 3c (STRICT): the same emission rule for a clonal link at its decided size; the relatives that place
            // the end and name another partner count against it.
            ++st.strictChecked;
            StrictInput si;
            si.n = links[i].d.n;
            si.overlap = links[i].d.overlap;
            si.walkLens = links[i].e.walkLens;
            si.walkCounts = links[i].e.walkCounts;
            si.exhaustive = links[i].e.exhaustive;
            si.placed = links[i].c.placed;
            si.agree = links[i].c.agree;
            si.contra = links[i].c.placed > links[i].c.agree ? links[i].c.placed - links[i].c.agree : 0;
            si.wAgree = links[i].c.wAgree;
            si.wContra = links[i].c.wContra;
            if (const char* why = strictVerdict(si, o)) {
                links[i].d.join = false;
                links[i].d.why = why;
                ++st.strictRefused;
                if (std::strcmp(why, "strict_periodic") == 0) ++st.strictPeriodic;
                if (std::strcmp(why, "strict_no_relative") == 0) ++st.strictNoRelative;
                continue;
            }
        }
        if (links[i].plasmid && links[i].d.join) {
            // plasmid records of one molecule share its copy number
            const double ca = links[i].u / 2 < covs.size() ? covs[links[i].u / 2] : 0.0;
            const double cb = links[i].v / 2 < covs.size() ? covs[links[i].v / 2] : 0.0;
            if (ca <= 0 || cb <= 0 || std::max(ca, cb) > o.plasmidDepthRatio * std::min(ca, cb)) {
                links[i].d.join = false;
                links[i].d.why = "depth";
                ++st.refusedDepth;
                continue;
            }
        }
        const std::string& why = links[i].d.why;
        if (!links[i].d.join && why == "size_ambiguous" && links[i].d.bracketN != INT32_MIN) ++st.bracketWould;
        if (!links[i].d.join) {
            if (why == "pair_contra") ++st.refusedPairs;
            else if (why == "graph_contra") ++st.refusedGraph;
            else ++st.refusedSize;
        }
    }

    if (!sideDir.empty()) {
        std::ofstream f(sideDir + "/links.tsv");
        f << "u\tv\tplasmid\tplaced\tagree\tw_agree\tw_contra\tnr_gap\tnr_spread\tanchored_a\tanchored_b\tedge_a\tedge_b\treachable\t"
             "exhaustive\twalks\telement\tcontra_graph\tpair_contra\tpanel_placed\tpanel_agree\tjoin\tn\tno_fill\t"
             "class\tbasis\twhy\tw_placed\toverlap_exact\tbracket_n\tflank_a32\tflank_b32\n";
        for (const Link& lk : links) {
            std::string w;
            for (int32_t x : lk.e.walkLens) w += (w.empty() ? "" : ";") + std::to_string(x);
            f << lk.u / 2 << ((lk.u & 1) ? ":tail" : ":head") << '\t' << lk.v / 2 << ((lk.v & 1) ? ":tail" : ":head")
              << '\t' << (lk.plasmid ? 1 : 0) << '\t' << lk.c.placed << '\t' << lk.c.agree << '\t'
              << fmtD(lk.c.wAgree, "%.3f") << '\t' << fmtD(lk.c.wContra, "%.3f") << '\t' << lk.c.gap << '\t'
              << lk.c.spread << '\t' << lk.e.anchoredA << '\t' << lk.e.anchoredB << '\t' << lk.e.edgeA << '\t'
              << lk.e.edgeB << '\t' << lk.e.reachable << '\t' << lk.e.exhaustive << '\t' << (w.empty() ? "." : w) << '\t'
              << lk.e.element << '\t' << lk.e.contraGraph << '\t' << lk.e.pairContra << '\t' << lk.e.panelPlaced
              << '\t' << lk.e.panelAgree << '\t' << lk.d.join << '\t' << lk.d.n << '\t' << lk.d.noFill << '\t'
              << clonalClassName(lk.d.cls) << '\t' << (lk.d.basis.empty() ? "." : lk.d.basis) << '\t'
              << (lk.d.why.empty() ? "." : lk.d.why) << '\t' << fmtD(lk.c.wPlaced, "%.3f") << '\t'
              << lk.e.overlapExact << '\t'
              << (lk.d.bracketN == INT32_MIN ? std::string(".") : std::to_string(lk.d.bracketN)) << '\t'
              << (lk.L.size() >= 32 ? lk.L.substr(lk.L.size() - 32) : lk.L) << '\t'
              << lk.Rs.substr(0, std::min<size_t>(32, lk.Rs.size())) << '\n';
        }
    }

    // ---- 6. chains -------------------------------------------------------------------------------
    std::vector<int64_t> partner(2 * R, -1);
    std::vector<int64_t> linkOf(2 * R, -1);
    for (size_t i = 0; i < links.size(); ++i) {
        if (!links[i].d.join) continue;
        partner[links[i].u] = static_cast<int64_t>(links[i].v);
        partner[links[i].v] = static_cast<int64_t>(links[i].u);
        linkOf[links[i].u] = linkOf[links[i].v] = static_cast<int64_t>(i);
    }
    // A cycle is cut at its weakest link (fewest agreeing relatives, then the largest gap); the wrap
    // itself is C2's to close, from the graph.
    {
        std::vector<char> seen(R, 0);
        for (size_t r0 = 0; r0 < R; ++r0) {
            if (seen[r0] || partner[2 * r0] < 0 || partner[2 * r0 + 1] < 0) continue;
            // walk from r0's tail; a cycle returns to r0's head
            std::vector<size_t> cyc;
            size_t r = r0, leave = 2 * r0 + 1;
            bool cycle = false;
            for (size_t steps = 0; steps <= R; ++steps) {
                seen[r] = 1;
                const int64_t lk = linkOf[leave];
                if (lk < 0) break;
                cyc.push_back(static_cast<size_t>(lk));
                const size_t enter = static_cast<size_t>(partner[leave]);
                r = enter / 2;
                if (r == r0) { cycle = true; break; }
                leave = enter ^ 1u;
            }
            if (!cycle || cyc.empty()) continue;
            size_t cut = cyc.front();
            for (size_t i : cyc) {
                const Link& a = links[i];
                const Link& b = links[cut];
                if (a.c.wAgree < b.c.wAgree - 1e-9 || (std::fabs(a.c.wAgree - b.c.wAgree) <= 1e-9 && a.c.gap > b.c.gap)) cut = i;
            }
            partner[links[cut].u] = partner[links[cut].v] = -1;
            linkOf[links[cut].u] = linkOf[links[cut].v] = -1;
            links[cut].d.join = false;
            links[cut].d.why = "cycle_cut";
            ++st.cycleCut;
        }
    }
    bool anyJoin = false;
    for (const Link& lk : links) if (lk.d.join) anyJoin = true;
    if (!anyJoin) return done(changed ? "no_join_after_breaks" : "no_join");

    // ---- 7. write the chains, every join judged by C1 --------------------------------------------
    seam.beginBatch(seqs);
    std::vector<char> used(R, 0);
    std::vector<std::string> out;
    std::vector<double> outCov;
    std::vector<char> outLm;
    const bool hadLm = !layoutMembers.empty();
    auto startsChain = [&](size_t r) { return partner[2 * r] < 0 || partner[2 * r + 1] < 0; };
    for (size_t r0 = 0; r0 < R; ++r0) {
        if (used[r0]) continue;
        if (!startsChain(r0)) continue;   // interior of a chain; reached from its start
        // orient the first record so the chain continues from its right end
        bool fwd = partner[2 * r0] < 0;   // head free: start forward
        size_t r = r0;
        std::string cur;
        double covW = 0;
        size_t covLen = 0;
        char lm = 0;
        auto add = [&](size_t rr, bool f) {
            std::string piece = f ? seqs[rr] : reverseComplement(seqs[rr]);
            covW += (rr < covs.size() ? covs[rr] : 0.0) * static_cast<double>(piece.size());
            covLen += piece.size();
            if (hadLm && rr < layoutMembers.size() && layoutMembers[rr]) lm = 1;
            used[rr] = 1;
            return piece;
        };
        auto flush = [&]() {
            if (cur.empty()) return;
            out.push_back(std::move(cur));
            outCov.push_back(covLen ? covW / static_cast<double>(covLen) : 0.0);
            outLm.push_back(lm);
            cur.clear();
            covW = 0;
            covLen = 0;
            lm = 0;
            ++st.chains;
        };
        cur = add(r, fwd);
        while (true) {
            const size_t leave = fwd ? 2 * r + 1 : 2 * r;
            const int64_t li = linkOf[leave];
            if (li < 0) break;
            const size_t enter = static_cast<size_t>(partner[leave]);
            const size_t nr = enter / 2;
            if (used[nr]) break;
            const bool nfwd = (enter & 1u) == 0;   // entered through its head: forward
            Link& lk = links[static_cast<size_t>(li)];
            std::string piece = nfwd ? seqs[nr] : reverseComplement(seqs[nr]);
            JudgeRequest rq;
            rq.source = Source::Clonal;
            rq.stage = "clonal";
            rq.left = &cur;
            rq.right = &piece;
            rq.claimedN = lk.d.n;
            rq.panelGap = static_cast<int32_t>(std::max<int64_t>(INT32_MIN / 2, std::min<int64_t>(lk.c.gap, INT32_MAX / 2)));
            rq.a = Port{static_cast<uint32_t>(r), fwd};
            rq.b = Port{static_cast<uint32_t>(nr), !nfwd};
            rq.panelSupport = lk.c.agree;
            rq.panelPairs = lk.c.placed;
            // Plasmid links carry their own positive evidence (decideLink: the isolate's walk at the relatives'
            // plasmid size, >= 2 relatives' plasmids agreeing, depth within 1.6x, no contradicting pairs), so C1's
            // plasmid rule (pairs confirm / PASS_EXACT / pair link) is not asked again: with the other record of the
            // plasmid already placed, a longer walk through it keeps the verdict from PASS_EXACT.
            rq.plasmidPass = false;
            // Round 2: the relatives' vouch (every placed relative agrees, >= 2 placed, at the walk / exact-overlap size,
            // no pair to a third end): C1's pair DENY is then no veto at a repeat end.
            {
                const bool strong = lk.c.placed >= 2 && (o.weighted ? (lk.c.wPlaced > 0 && lk.c.wAgree >= 0.9 * lk.c.wPlaced)
                                                                    : lk.c.agree >= lk.c.placed);
                rq.clonalVouch = o.vouch && strong && lk.e.pairContra == 0 && (lk.d.walkMatch || lk.d.overlap);
            }
            // an exact-overlap join is judged at its overlap (claimedN = -overlap); the record keeps a 1-N butt
            // (C1 records writtenN = 1) and the genome view merges the overlap
            const Decision dec = seam.judge(rq);
            if (dec.vouched) ++st.vouched;
            ClonalAnnot an;
            an.set = true;
            an.cls = lk.d.cls;
            an.nrK = static_cast<uint32_t>(K);
            an.nrPlaced = lk.c.placed;
            an.nrAgree = lk.c.agree;
            an.wPlaced = lk.c.wPlaced;
            an.wAgree = lk.c.wAgree;
            an.wContra = lk.c.wContra;
            an.nrGap = lk.c.gap;
            an.dNear = ns.dNear();
            an.group = ns.group();
            an.panelPlaced = lk.e.panelPlaced;
            an.panelAgree = lk.e.panelAgree;
            an.edgeA = lk.e.anchoredA ? (lk.e.edgeA ? 1 : 0) : -1;
            an.edgeB = lk.e.anchoredB ? (lk.e.edgeB ? 1 : 0) : -1;
            an.walk = lk.e.anchoredA && lk.e.anchoredB ? (lk.e.reachable ? 1 : 0) : -1;
            an.walkMatch = lk.d.walkMatch ? 1 : 0;
            an.element = lk.e.reachable ? (lk.e.element ? 1 : 0) : -1;
            an.pairContra = lk.e.pairContra;
            an.basis = lk.d.basis;
            // An exact-overlap join that C1 would resize (a longer walk between the ends) is a conflict: refused.
            const bool overlapConflict = dec.keep && lk.d.overlap && dec.writeN != INT32_MIN;
            if (overlapConflict) {
                ++st.overlapC1Conflict;
                seam.clonalBreak(dec.junction, "clonal_overlap_resize_conflict");
            }
            if (!dec.keep || overlapConflict) {
                ++st.refusedC1;
                an.action = overlapConflict ? "refuse:c1_overlap_conflict" : "refuse:c1";
                eng.rows[dec.junction] = an;
                flush();
                fwd = nfwd;
                r = nr;
                cur = add(r, fwd);
                continue;
            }
            const size_t nn = static_cast<size_t>(std::max<int32_t>(1, dec.writeN != INT32_MIN ? dec.writeN : lk.d.n));
            seam.clonalAdmit(dec.junction, lk.d.noFill, std::string("clonal_") + clonalClassName(lk.d.cls),
                             lk.d.overlap ? -lk.d.n : 0);
            if (lk.d.overlap) ++st.overlapJoins;
            if (lk.d.bracket) ++st.bracketJoins;
            if (lk.d.c1) ++st.walksizeLinks;
            if (o.kc2 && lk.d.noFill && lk.d.cls == kNonPositional) ++st.kc2Links;
            an.action = "join";
            eng.rows[dec.junction] = an;
            ++st.joins;
            if (lk.plasmid) ++st.plasmidJoins;
            if (lk.d.cls == kPositional) ++st.joinsPositional;
            else if (lk.d.cls == kLayout) ++st.joinsLayout;
            else if (lk.d.cls != kPlasmid) ++st.joinsNonPositional;
            if (lk.d.noFill) ++st.noFill;
            if (lk.d.basis == "nr") ++st.sizedNr; else ++st.sizedGraph;
            cur.append(nn, 'N');
            piece = add(nr, nfwd);
            cur += piece;
            fwd = nfwd;
            r = nr;
        }
        flush();
    }
    // records never reached from a chain start (cycles cut above leave none) and the rest, in order
    for (size_t r = 0; r < R; ++r) {
        if (used[r]) continue;
        out.push_back(seqs[r]);
        outCov.push_back(r < covs.size() ? covs[r] : 0.0);
        outLm.push_back(hadLm && r < layoutMembers.size() ? layoutMembers[r] : 0);
    }
    seqs.swap(out);
    covs.swap(outCov);
    if (hadLm) layoutMembers.swap(outLm);
    changed = true;
    st.recordsOut = seqs.size();
    st.ms = nowMs() - t0;
    if (in.verbose) {
        std::fprintf(stderr, "[4b2/7] om2 clonal: %zu relatives (d_near %.5f, group %c), %zu mutual proposals, "
                             "%zu joins (%zu positional, %zu layout, %zu non-positional, %zu no-fill), "
                             "%zu contradicted joins broken, %zu -> %zu records, %.1fs\n",
                     eng.relatives(), ns.dNear(), ns.group(), st.mutual, st.joins, st.joinsPositional,
                     st.joinsLayout, st.joinsNonPositional, st.noFill, st.contraBroken, st.records,
                     st.recordsOut, st.ms / 1000.0);
    }
    return st;
}

// ------------------------------------------------------------------------------ the writer's view

ClonalSurface::ClonalSurface(ClonalEngine& eng, const SeamContext* seam, ClonalStats& stats)
    : eng_(eng), seam_(seam), stats_(stats) {
    if (!eng_.options().calPath.empty()) {
        std::string err;
        if (!cal_.load(eng_.options().calPath, err))
            std::fprintf(stderr, "[om2-clonal] warning: calibration table not used: %s\n", err.c_str());
    }
}

void ClonalSurface::prepare(const std::vector<std::string>& records) {
    eng_.setAssembly(records);
    // round 3 (CONFRULE): the run's own read pairs refused the relatives' order at >= DISCORD links / joins
    stats_.discordEvents = stats_.refusedPairs + stats_.contraBroken;
    stats_.discord = eng_.options().confRule && eng_.options().discord > 0 &&
                     stats_.discordEvents >= static_cast<size_t>(eng_.options().discord);
}

std::vector<ClonalAnnot> ClonalSurface::annotate(const std::vector<Req>& reqs) {
    std::vector<ClonalAnnot> out(reqs.size());
    std::vector<JunctionQuery> qs;
    qs.reserve(reqs.size());
    for (const Req& r : reqs) qs.push_back(eng_.query(r.left, r.right, r.asserted, r.sized));
    std::vector<TrackSupport> panel;
    eng_.panelSupport(qs, panel);
    const EvidenceIndex* ev = seam_ ? seam_->evidence() : nullptr;
    const bool evOk = ev && !ev->abstained();
    const NearestSet& ns = eng_.nearest();
    const ClonalOptions& o = eng_.options();
    for (size_t i = 0; i < reqs.size(); ++i) {
        const Req& rq = reqs[i];
        ClonalAnnot a;
        const Junction* j = rq.j;
        const ClonalAnnot* stage = nullptr;
        if (j) {
            auto it = eng_.rows.find(j->id);
            if (it != eng_.rows.end()) stage = &it->second;
        }
        const TrackSupport ts = eng_.relativeSupport(qs[i]);
        a.set = true;
        a.nrK = static_cast<uint32_t>(eng_.relatives());
        a.nrPlaced = ts.placed;
        a.nrAgree = ts.agree;
        a.wPlaced = ts.wPlaced;
        a.wAgree = ts.wAgree;
        a.wContra = ts.wContra;
        a.nrGapDiff = ts.gapDiff;
        a.nrContra = ts.contra;
        a.nrGap = ts.agreeGaps.empty() ? median(ts.gaps) : median(ts.agreeGaps);
        a.dNear = ns.dNear();
        a.group = ns.group();
        a.panelPlaced = panel[i].placed;
        a.panelAgree = panel[i].agree;
        a.sized = rq.sized;
        // graph: the flank->element and element->flank edges, and a walk of the asserted size
        if (evOk) {
            const EndAnchor A = ev->anchorExit(rq.left);
            const EndAnchor B = ev->anchorEntry(rq.right);
            const UnitigGraph& g = ev->graph();
            if (A.ok) a.edgeA = g.exits(A.node.id, A.node.o).empty() ? 0 : 1;
            if (B.ok) { const ONode bb = flip(B.node); a.edgeB = g.exits(bb.id, bb.o).empty() ? 0 : 1; }
            std::vector<int32_t> lens;
            if (j && (!j->walks.empty() || j->gmin > 0)) {
                a.walk = 1;
                for (const WalkClass& w : j->walks) lens.push_back(w.len);
                if (lens.empty()) lens.push_back(j->gmin);
            } else if (A.ok && B.ok) {
                const WalkSet ws = ev->walks(A, B);
                a.walk = ws.reachable ? 1 : 0;
                for (const WalkClass& w : ws.classes) lens.push_back(w.len);
                if (lens.empty() && ws.reachable) lens.push_back(ws.gmin);
                a.element = ws.reachable ? ((ws.innerRepeatBp > 0 || ws.innerMaxCopy > 1.5f) ? 1 : 0) : -1;
            }
            if (a.walk == 1 && rq.sized) {
                a.walkMatch = 0;
                const int64_t tol = std::max<int64_t>(500, std::llabs(rq.asserted) / 10);
                for (int32_t w : lens) if (std::llabs(static_cast<int64_t>(w) - rq.asserted) <= tol) a.walkMatch = 1;
                // a filled junction: the fill IS a walk of the graph
                if (j && !j->fillSeq.empty()) a.walkMatch = 1;
            }
        }
        a.pairContra = j ? j->pairContra : 0;
        // class
        const bool repeatEnds = j ? (j->endA == EndClass::Repeat || j->endB == EndClass::Repeat) : false;
        const bool panelKnown = a.panelPlaced >= o.panelMin;
        const bool positional = panelKnown && a.panelAgree >= o.positional * a.panelPlaced;
        if (rq.wrap || (j && j->source == Source::Wrap)) a.cls = kWrap;
        else if (rq.replicon == 'p') a.cls = kPlasmid;
        else if (j && (j->source == Source::Resolver || j->source == Source::ResolverUnknown100)) a.cls = kGraph;
        else if (stage) a.cls = stage->cls;
        else if (positional) a.cls = (repeatEnds || a.element == 1) ? kPositional : kLayout;
        else a.cls = kNonPositional;
        if (stage) {
            a.action = stage->action;
            a.basis = stage->basis;
            a.pairContra = std::max(a.pairContra, stage->pairContra);
            if (a.element < 0) a.element = stage->element;
        } else {
            a.action = "keep";
            a.basis = j ? "ledger" : "unrecorded";
        }
        const char tier = j ? tierChar(j->tier) : 'E';
        a.confidence = cal_.confidence(a, tier, j != nullptr);
        // Round 2: a junction side with fewer than 2 assembly-unique markers in its N-free context (a short piece between
        // two gaps, or repeat-only sequence) can be checked neither by the relatives nor against any reference: its
        // probability is capped below the `confident` claim (dev E. coli GCF051549665v1: both UNSCORABLE `confident`
        // junctions had a 2.4 / 3.9 kb repeat-bounded side).
        if (std::min(qs[i].left.size(), qs[i].right.size()) < 2) a.confidence = std::min(a.confidence, 0.985);
        if (o.confRule && a.confidence >= o.confident) {
            // Round 3 (CONFRULE): `confident` needs every placed relative (>= 2, by count) at this size, both N-free
            // sides >= 2 kb (the join can be checked at all), the isolate's walk (if any) at the asserted size, a gap
            // <= 10 kb, and no repeat end in a run whose own read pairs refused the relatives' order (discord). Dev Kp
            // round 2 counterfactual: 0.995 (LTO) / 0.996 (LCO2) against 0.973 / 0.970 for the round-2 labels.
            size_t lf = 0, rf = 0;
            while (lf < rq.left.size() && rq.left[rq.left.size() - 1 - lf] != 'N' && rq.left[rq.left.size() - 1 - lf] != 'n') ++lf;
            while (rf < rq.right.size() && rq.right[rf] != 'N' && rq.right[rf] != 'n') ++rf;
            const bool repeatEnd = j && (j->endA == EndClass::Repeat || j->endB == EndClass::Repeat);
            const bool ok = a.sized && a.nrPlaced >= 2 && a.nrAgree == a.nrPlaced &&
                            static_cast<int64_t>(std::min(lf, rf)) >= o.confSide && !(a.walk == 1 && a.walkMatch == 0) &&
                            std::llabs(rq.asserted) <= o.confMaxGap && !(stats_.discord && repeatEnd);
            if (!ok) { a.confidence = std::min(a.confidence, 0.985); ++stats_.confCapped; }
        }
        if (!a.sized) a.claim = "gap";
        else if (a.confidence >= o.confident && a.dNear <= o.confDmax) a.claim = "confident";
        else if (a.confidence >= 0.9) a.claim = "supported";
        else a.claim = "provisional";
        // 1.5, KC1: no `confident` label ships (no threshold of the hand-set prior reached >= 0.99 with Wilson LB >=
        // 0.98 under the binding scorer on dev). The junction keeps its evidence columns and its numeric prior.
        if (!o.allowConfident) {
            stats_.kc1 = true;
            if (a.claim == "confident") { a.claim = "supported"; ++stats_.confidentSuppressed; }
        }
        ++stats_.annotated;
        if (a.claim == "confident") ++stats_.confident;
        else if (a.claim == "supported") ++stats_.supported;
        else if (a.claim == "gap") ++stats_.gapClaims;
        else ++stats_.provisional;
        out[i] = a;
    }
    return out;
}

std::string ClonalSurface::header() {
    return "clonal_class\tconfidence\tclaim\tnr_k\tnr_placed\tnr_agree\tnr_gapdiff\tnr_contra\tnr_gap\td_near_tool\t"
           "d_group\tnr_w_placed\tnr_w_agree\tnr_w_contra\tpanel_placed\tpanel_agree_frac\tgraph_edge_a\tgraph_edge_b\tgraph_walk\twalk_matches\telement\t"
           "pair_contra\tclonal_action\tclonal_basis";
}

std::string ClonalSurface::columns(const ClonalAnnot& a) {
    if (!a.set) {
        std::string s;
        for (int c = 0; c < 24; ++c) s += c ? "\t." : ".";
        return s;
    }
    auto tri = [](int v) { return v < 0 ? std::string(".") : std::to_string(v); };
    std::ostringstream o;
    o << clonalClassName(a.cls) << '\t' << (a.confidence < 0 ? std::string(".") : fmtD(a.confidence, "%.5f")) << '\t'
      << a.claim << '\t' << a.nrK << '\t' << a.nrPlaced << '\t' << a.nrAgree << '\t' << a.nrGapDiff << '\t' << a.nrContra
      << '\t' << (a.nrGap == INT64_MIN ? std::string(".") : std::to_string(a.nrGap)) << '\t' << fmtD(a.dNear, "%.6f")
      << '\t' << a.group << '\t' << (a.wPlaced < 0 ? std::string(".") : fmtD(a.wPlaced, "%.3f")) << '\t'
      << (a.wPlaced < 0 ? std::string(".") : fmtD(a.wAgree, "%.3f")) << '\t'
      << (a.wPlaced < 0 ? std::string(".") : fmtD(a.wContra, "%.3f")) << '\t' << a.panelPlaced << '\t'
      << (a.panelPlaced ? fmtD(a.panelFrac(), "%.3f") : std::string(".")) << '\t' << tri(a.edgeA) << '\t'
      << tri(a.edgeB) << '\t' << tri(a.walk) << '\t' << tri(a.walkMatch) << '\t' << tri(a.element) << '\t'
      << a.pairContra << '\t' << a.action << '\t' << a.basis;
    return o.str();
}

std::string ClonalSurface::closureHeader() const {
    const NearestSet& ns = eng_.nearest();
    const ClonalOptions& o = eng_.options();
    std::string s = "#clonal\tenabled=1\tk=" + std::to_string(eng_.relatives()) + "\tk_min=" + std::to_string(o.kmin) +
                    "\td_near_tool=" + fmtD(ns.ok ? ns.dNear() : -1, "%.6f") + "\tgroup=" + std::string(1, ns.group()) +
                    "\tiso_markers=" + std::to_string(ns.isoMarkers) + "\tiso_chr_markers=" +
                    std::to_string(ns.isoChrMarkers) + "\ttracks=" + std::to_string(ns.tracks) +
                    "\tjoins=" + std::to_string(stats_.joins) + "\tbroken_contra=" + std::to_string(stats_.contraBroken) +
                    "\tstage=" + (stats_.ran ? std::string("ran") : "skipped:" + stats_.skip) + "\n";
    for (size_t i = 0; i < ns.ranked.size(); ++i) {
        const Relative& r = ns.ranked[i];
        s += "#nearest\t" + std::to_string(i + 1) + "\t" + r.name + "\tdist=" + fmtD(r.dist, "%.6f") + "\tjaccard=" +
             fmtD(r.jaccard, "%.5f") + "\tcontainment=" + fmtD(r.containment, "%.5f") + "\tshared=" +
             std::to_string(r.shared) + "\ttrack_markers=" + std::to_string(r.trackMarkers) +
             (i < eng_.relatives() ? "\tused=1" : "\tused=0") + "\n";
    }
    return s;
}

std::string ClonalSurface::agpComment(const std::string& id, const ClonalAnnot& a) const {
    if (!a.set) return std::string();
    return "# om2-clonal junction=" + id + " class=" + clonalClassName(a.cls) + " claim=" + a.claim +
           " confidence=" + (a.confidence < 0 ? std::string("NA") : fmtD(a.confidence, "%.5f")) + " nr_agree=" +
           std::to_string(a.nrAgree) + "/" + std::to_string(a.nrPlaced) + " basis=" + a.basis + "\n";
}

std::string ClonalSurface::jsonBlock(int indent) const {
    const NearestSet& ns = eng_.nearest();
    const std::string pad(static_cast<size_t>(indent) * 4, ' '), pad2(static_cast<size_t>(indent + 1) * 4, ' '),
        pad3(static_cast<size_t>(indent + 2) * 4, ' ');
    std::string s = pad + "\"clonal\": {\n";
    s += pad2 + "\"enabled\": true,\n";
    s += pad2 + "\"k\": " + std::to_string(eng_.relatives()) + ",\n";
    s += pad2 + "\"d_near_tool\": " + fmtD(ns.ok ? ns.dNear() : -1, "%.6f") + ",\n";
    s += pad2 + "\"group\": \"" + std::string(1, ns.group()) + "\",\n";
    s += pad2 + "\"counters\": \"" + formatClonalCounters(stats_) + "\",\n";
    s += pad2 + "\"nearest\": [\n";
    for (size_t i = 0; i < ns.ranked.size(); ++i) {
        const Relative& r = ns.ranked[i];
        s += pad3 + "{\"rank\": " + std::to_string(i + 1) + ", \"accession\": \"" + r.name + "\", \"dist\": " +
             fmtD(r.dist, "%.6f") + ", \"jaccard\": " + fmtD(r.jaccard, "%.5f") + ", \"shared\": " +
             std::to_string(r.shared) + ", \"used\": " + (i < eng_.relatives() ? "true" : "false") + "}" +
             (i + 1 < ns.ranked.size() ? "," : "") + "\n";
    }
    s += pad2 + "]\n";
    s += pad + "},\n";
    return s;
}

std::string formatClonalCounters(const ClonalStats& s) {
    char b[4096];
    std::snprintf(b, sizeof b,
                  "[om2-clonal] enabled=%d ran=%d skip=%s k=%d d_near=%s group=%c records=%zu candidates=%zu "
                  "ends_placed=%zu proposals=%zu mutual=%zu refused_pairs=%zu refused_graph=%zu refused_size=%zu "
                  "refused_c1=%zu cycle_cut=%zu joins=%zu joins_positional=%zu joins_layout=%zu "
                  "joins_non_positional=%zu no_fill=%zu sized_graph=%zu sized_nr=%zu contra_checked=%zu "
                  "contra_broken=%zu contra_kept_proven=%zu chains=%zu records_out=%zu nrp=%d plasmid_tracks=%zu "
                  "plasmid_candidates=%zu plasmid_mutual=%zu plasmid_joins=%zu refused_depth=%zu annotated=%zu "
                  "confident=%zu supported=%zu provisional=%zu gap_claims=%zu ms=%.0f overlap_joins=%zu "
                  "overlap_c1_conflict=%zu vouched=%zu bracket_would=%zu bracket_joins=%zu merged=%zu size_checked=%zu "
                  "size_agree=%zu resized=%zu walksize_checked=%zu walksize_resized=%zu walksize_links=%zu "
                  "kc2_checked=%zu kc2_nofill=%zu kc2_links=%zu conf_capped=%zu discord_events=%zu discord=%d cap_walk_c1=%zu unverified_broken=%zu "
                  "strict=%d strict_nr0=%d strict_checked=%zu strict_broken=%zu strict_refused=%zu strict_periodic=%zu strict_no_relative=%zu",
                  s.enabled ? 1 : 0, s.ran ? 1 : 0, s.skip.c_str(), s.k,
                  s.dNear < 0 ? "NA" : fmtD(s.dNear, "%.6f").c_str(), s.group, s.records, s.candidates, s.endsPlaced,
                  s.proposals, s.mutual, s.refusedPairs, s.refusedGraph, s.refusedSize, s.refusedC1, s.cycleCut,
                  s.joins, s.joinsPositional, s.joinsLayout, s.joinsNonPositional, s.noFill, s.sizedGraph, s.sizedNr,
                  s.contraChecked, s.contraBroken, s.contraKeptProven, s.chains, s.recordsOut, s.nrpLoaded ? 1 : 0,
                  s.plasmidTracks, s.plasmidCandidates, s.plasmidMutual, s.plasmidJoins, s.refusedDepth, s.annotated,
                  s.confident, s.supported, s.provisional, s.gapClaims, s.ms, s.overlapJoins, s.overlapC1Conflict,
                  s.vouched, s.bracketWould, s.bracketJoins, s.merged, s.sizeChecked, s.sizeAgree, s.resized,
                  s.walksizeChecked, s.walksizeResized, s.walksizeLinks, s.kc2Checked, s.kc2NoFill, s.kc2Links,
                  s.confCapped, s.discordEvents, s.discord ? 1 : 0, s.capWalkC1, s.unverifiedBroken, s.strict ? 1 : 0,
                  s.strictNr0 ? 1 : 0, s.strictChecked, s.strictBroken, s.strictRefused, s.strictPeriodic, s.strictNoRelative);
    return b;
}

void logClonalCounters(const ClonalStats& s, std::FILE* f) { std::fprintf(f, "%s\n", formatClonalCounters(s).c_str()); }

}  // namespace om2
}  // namespace ts
