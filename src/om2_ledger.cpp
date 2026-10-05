#include "om2_ledger.h"
#include "om2_clonal.h"   // round 3: isolatedShortestWalk (cap_walk_c1)

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <fstream>
#include <sstream>
#include <unordered_map>

#include "envflags.h"
#include "graph.h"

namespace ts {
namespace om2 {

namespace {

constexpr int32_t kSlack = 1000;       // gate_eval: |gmin - n| and skip slack
constexpr double kPlasmidDepth = 1.75; // C1d: depth >= 1.75 theta votes plasmid
constexpr size_t kFlank = 32;
constexpr size_t kWin = 25000;         // piece-end window handed to the judge

double nowMs() {
    return std::chrono::duration<double, std::milli>(
               std::chrono::steady_clock::now().time_since_epoch()).count();
}

bool eq(const char* a, const char* b) { return a && b && std::strcmp(a, b) == 0; }

std::string tail(const std::string& s, size_t n) { return s.size() <= n ? s : s.substr(s.size() - n); }
std::string head(const std::string& s, size_t n) { return s.size() <= n ? s : s.substr(0, n); }

std::string anchorStr(const EndAnchor& a) {
    if (!a.ok) return a.why;
    char b[96];
    std::snprintf(b, sizeof(b), "%u%c:%d:%d", a.node.id, a.node.o ? '-' : '+', a.off, a.slide);
    return b;
}

const char* nClass(Source s, int32_t n) {
    if (s == Source::JoinOverlap || s == Source::LayoutOverlap) return "overlap";
    if (s == Source::JoinButt1 || s == Source::LayoutButt1) return "butt1";
    if (s == Source::LayoutCap2000) return "cap2000";
    if (n <= 100) return "le100";
    if (n <= 3000) return "le3000";
    return "gt3000";
}

bool isModelSource(Source s) {
    return s != Source::Resolver && s != Source::ResolverUnknown100;
}

char endChar(EndClass e) {
    switch (e) {
        case EndClass::Unique: return 'U';
        case EndClass::Repeat: return 'R';
        case EndClass::DeadEnd: return 'D';
        default: return '?';
    }
}

}  // namespace

// ------------------------------------------------------------------- config

const char* SeamConfig::modeName(Mode m) {
    switch (m) {
        case Mode::Off: return "off";
        case Mode::Audit: return "audit";
        case Mode::Act: return "act";
    }
    return "?";
}
const char* SeamConfig::unsizedName(Unsized u) {
    switch (u) {
        case Unsized::Butt: return "butt";
        case Unsized::Break: return "break";
        case Unsized::Sized: return "sized";
    }
    return "?";
}
const char* SeamConfig::capName(Cap c) { return c == Cap::Keep ? "keep" : "size"; }

SeamConfig SeamConfig::fromEnv() {
    SeamConfig c;
    const char* seam = env::text("TESSERACT_OM2_SEAM");
    if (eq(seam, "audit")) c.mode = Mode::Audit;
    else if (eq(seam, "act")) c.mode = Mode::Act;
    const char* uns = env::text("TESSERACT_OM2_UNSIZED");
    if (eq(uns, "butt")) c.unsized = Unsized::Butt;
    else if (eq(uns, "sized")) c.unsized = Unsized::Sized;
    const char* cap = env::text("TESSERACT_OM2_CAP");
    if (eq(cap, "keep")) c.cap = Cap::Keep;
    c.pairs = env::on("TESSERACT_OM2_PAIRS", true);
    c.plasmidRule = env::on("TESSERACT_OM2_PLASMID_RULE", true);
    c.closeGapsPolicy = env::on("TESSERACT_OM2_CLOSEGAPS_POLICY", true);
    if (const char* p = env::text("TESSERACT_OM2_ADMIT")) c.admitPath = p;
    c.tangleMax = static_cast<long>(env::integer("TESSERACT_OM2_TANGLE_MAX", 5000));
    const bool anyPolicy = uns || cap || env::isSet("TESSERACT_OM2_PLASMID_RULE") ||
                           env::isSet("TESSERACT_OM2_CLOSEGAPS_POLICY") || !c.admitPath.empty();
    c.evidence = env::on("TESSERACT_OM2_EVIDENCE", false) || anyPolicy ||
                 env::isSet("TESSERACT_OM2_PAIRS") || env::isSet("TESSERACT_OM2_TANGLE_MAX") ||
                 c.mode != Mode::Off;
    c.policyWithoutAct = anyPolicy && c.mode != Mode::Act;
    c.capWalkC1 = env::on("TESSERACT_OM2_CLONAL", false) && env::on("TESSERACT_OM2_CLONAL_CAPWALK", false);
    return c;
}

bool AdmitTable::load(const std::string& path, std::string& err) {
    std::ifstream in(path);
    if (!in) { err = "cannot open " + path; return false; }
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '#') {
            if (line.compare(0, 7, "#tau_S=") == 0) tauS = std::atof(line.c_str() + 7);
            else if (line.compare(0, 7, "#tau_G=") == 0) tauG = std::atof(line.c_str() + 7);
            continue;
        }
        // cls <TAB> n <TAB> misjoins <TAB> upper
        std::vector<std::string> col;
        std::stringstream ss(line);
        std::string f;
        while (std::getline(ss, f, '\t')) col.push_back(f);
        if (col.size() < 4 || col[0] == "cls") continue;
        upper[col[0]] = std::atof(col[3].c_str());
    }
    if (tauS < 0 || tauG < 0) { err = path + ": missing #tau_S= / #tau_G= header"; return false; }
    loaded = true;
    return true;
}

// ------------------------------------------------------------------ context

SeamContext::SeamContext(const SeamConfig& cfg, std::unique_ptr<EvidenceIndex> ev,
                         std::unique_ptr<PairIndex> pairs, std::unique_ptr<PortMap> ports)
    : cfg_(cfg), ev_(std::move(ev)), pairs_(std::move(pairs)), ports_(std::move(ports)) {
    if (!cfg_.admitPath.empty()) {
        std::string err;
        if (!admit_.load(cfg_.admitPath, err))
            std::fprintf(stderr, "[om2-seam] warning: admission table not used: %s\n", err.c_str());
    }
}

std::unique_ptr<SeamContext> SeamContext::create(const UnitigGraph& graph,
                                                 const std::vector<std::string>& pieces,
                                                 const SequenceStore& reads,
                                                 const InsertModel& insert, int threads) {
    const SeamConfig cfg = SeamConfig::fromEnv();
    if (!cfg.active()) return nullptr;
    const double t0 = nowMs();
    EvidenceOptions opt;
    opt.tangleMax = cfg.tangleMax;
    std::unique_ptr<EvidenceIndex> ev = EvidenceIndex::build(graph, opt);
    std::unique_ptr<PairIndex> pairs;
    double pairMs = 0;
    if (cfg.mode != SeamConfig::Mode::Off && cfg.pairs) {
        const double tp = nowMs();
        pairs.reset(new PairIndex(pieces, reads, insert, threads, opt));
        pairMs = nowMs() - tp;
    }
    std::unique_ptr<PortMap> ports(new PortMap(pieces));
    std::unique_ptr<SeamContext> ctx(new SeamContext(cfg, std::move(ev), std::move(pairs), std::move(ports)));
    ctx->counters.pairBuildMs = pairMs;
    if (cfg.policyWithoutAct) {
        std::fprintf(stderr, "[om2-seam] note: policy flags act only with TESSERACT_OM2_SEAM=act; "
                             "recorded as would_* in om2_seams.tsv\n");
    }
    if (ctx->gating()) ctx->recordResolverGaps(pieces);
    ctx->createMs = nowMs() - t0;
    return ctx;
}

bool SeamContext::plasmidByDepth(double cov) const {
    return ev_ && ev_->theta() > 0 && cov >= kPlasmidDepth * ev_->theta();
}

void SeamContext::beginBatch(const std::vector<std::string>& pieces) {
    if (ev_) ev_->setPlaced(pieces);
}

std::string SeamContext::classKey(const Junction& j, const JunctionAudit& x) const {
    const char* pair = !x.pairTestable ? "na" : x.pairAbstain ? "abstain"
                     : x.pairCall > 0 ? "deny" : x.pairCall < 0 ? "confirm" : "none";
    char b[256];
    std::snprintf(b, sizeof(b), "src=%s|v=%s|t=%s|n=%s|e=%c%c|p=%s", sourceName(j.source),
                  verdictName(j.verdict), tierName(j.tier), nClass(j.source, j.claimedN),
                  endChar(j.endA), endChar(j.endB), pair);
    return b;
}

Decision SeamContext::judge(const JudgeRequest& rq) {
    const double t0 = nowMs();
    Decision dec;
    Junction J;
    JunctionAudit X;
    J.id = static_cast<uint32_t>(ledger_.j.size());
    dec.junction = J.id;
    J.source = rq.source;
    J.claimedN = rq.claimedN;
    J.panelSupport = rq.panelSupport;
    J.panelGenomes = rq.panelPairs;
    X.stage = rq.stage;
    X.passA = rq.a;
    X.passB = rq.b;
    X.panelGap = rq.panelGap;
    X.judgeOnly = rq.judgeOnly;
    const std::string& L = *rq.left;
    const std::string& R = *rq.right;
    J.flankL32 = tail(L, kFlank);
    J.flankR32 = head(R, kFlank);
    for (int32_t d : {0, 250, 1000, 3000, 8000, 20000}) {
        if (L.size() >= static_cast<size_t>(d) + 63) {
            std::string q = L.substr(L.size() - static_cast<size_t>(d) - 63, 63);
            if (q.find('N') == std::string::npos) X.probesL.emplace_back(d, std::move(q));
        }
        if (R.size() >= static_cast<size_t>(d) + 63) {
            std::string q = R.substr(static_cast<size_t>(d), 63);
            if (q.find('N') == std::string::npos) X.probesR.emplace_back(d, std::move(q));
        }
    }
    if (ports_) {
        X.originA = ports_->exitPort(L);
        X.originB = ports_->entryPort(R);
    }
    J.a = X.originA != UINT32_MAX ? Port{X.originA >> 1, (X.originA & 1u) != 0} : rq.a;
    J.b = X.originB != UINT32_MAX ? Port{X.originB >> 1, (X.originB & 1u) != 0} : rq.b;

    const bool isOverlap = rq.source == Source::JoinOverlap || rq.source == Source::LayoutOverlap;
    const int32_t n = rq.claimedN;
    bool contra = false;
    WalkSet ws;
    EndAnchor A, B;
    if (!ev_ || ev_->abstained()) {
        J.verdict = Verdict::Abstain;
    } else {
        A = ev_->anchorExit(L, n < 0 ? -n : 0);
        B = ev_->anchorEntry(R);
        X.anchorA = anchorStr(A);
        X.anchorB = anchorStr(B);
        if (A.ok) { X.anchorSeqA = L.substr(L.size() - static_cast<size_t>(A.slide) - 63, 63); X.slideA = A.slide; }
        if (B.ok) { X.anchorSeqB = R.substr(static_cast<size_t>(B.slide), 63); X.slideB = B.slide; }
        if (!A.ok || !B.ok) {
            J.verdict = Verdict::Unanchored;
        } else {
            const Multiplicity& ma = ev_->mult(A.node.id);
            const Multiplicity& mb = ev_->mult(B.node.id);
            J.copyA = ma.copy;
            J.copyB = mb.copy;
            const UnitigGraph& g = ev_->graph();
            J.endA = ma.hardSingle ? EndClass::Unique
                   : g.exits(A.node.id, A.node.o).empty() ? EndClass::DeadEnd : EndClass::Repeat;
            const ONode bb = flip(B.node);
            J.endB = mb.hardSingle ? EndClass::Unique
                   : g.exits(bb.id, bb.o).empty() ? EndClass::DeadEnd : EndClass::Repeat;
            ws = ev_->walks(A, B);
            const FirstUnique uL = ev_->firstUnique(A.node);
            const FirstUnique uR = ev_->firstUnique(bb);
            J.uA = static_cast<uint16_t>(std::min<size_t>(uL.ids.size(), 65535));
            J.uB = static_cast<uint16_t>(std::min<size_t>(uR.ids.size(), 65535));
            for (uint32_t id : uL.ids) if (id != B.node.id) ++X.uOther;
            for (uint32_t id : uR.ids) if (id != A.node.id) ++X.uOther;
            J.walks = ws.classes;
            J.exhaustive = ws.exhaustive;
            J.hairpin = ws.hairpin;
            X.distinct = ws.distinct;
            X.capped = ws.capped;
            X.innerUnique = ws.innerUnique;
            X.innerUniquePlaced = ws.innerUniquePlaced;
            X.reachable = ws.reachable;
            if (ws.reachable) {
                J.gmin = ws.gmin;
                if (ws.gmin > n + kSlack) {
                    J.verdict = ws.innerUniquePlaced > 0 ? Verdict::BreakDoubleUse : Verdict::Resize;
                } else if (std::abs(ws.gmin - n) <= kSlack && ws.distinct == 1 && !ws.capped) {
                    J.verdict = Verdict::PassExact;
                } else {
                    J.verdict = Verdict::PassWalk;
                }
            } else if (X.uOther > 0) {
                contra = true;
                J.verdict = Verdict::ContraRefused;   // until pairs confirm
            } else {
                J.verdict = Verdict::Silent;
            }
        }
    }
    X.graphVerdict = contra ? Verdict::ContraRefused : J.verdict;

    // ---- pairs ---------------------------------------------------------------------
    if (pairs_ && pairs_->usable() && !rq.judgeOnly) {
        const PairCall pc = pairs_->test(X.originA, X.originB, n);
        X.pairTestable = pc.testable;
        X.pairAbstain = pc.abstain;
        X.pairLog10LR = pc.log10LR;
        X.pairCall = pc.testable && !pc.abstain ? pc.call : 0;
        X.ctrlA = pc.ctrlA;
        X.ctrlB = pc.ctrlB;
        J.pairLambda = pc.lambda;
        J.pairK = pc.k;
        J.pairContra = pc.contra;
        if (pc.testable) {
            ++counters.pairTestable;
            if (pc.abstain) ++counters.pairAbstain;
            else if (X.pairCall < 0) ++counters.pairConfirm;
            else if (X.pairCall > 0) ++counters.pairDeny;
        }
    }
    if (contra && X.pairCall < 0) J.verdict = Verdict::ContraPairsOk;
    // Pairs veto a junction only where the graph connects its flanks. At a SILENT junction the
    // flanks are graph dead ends, most often a coverage dropout, and a dropout starves the
    // crossing count too: on dev E. coli GCF013376575v1 (insert 186 +- 115) the pilot rule
    // denied 55 seams, 53 SILENT and all 53 truly adjacent, while the 6 true deletions were
    // all graph-connected (PASS_WALK / RESIZE / BREAK). Silence is not a veto.
    if (!contra && X.pairCall > 0 &&
        (J.verdict == Verdict::PassExact || J.verdict == Verdict::PassWalk)) {
        if (rq.clonalVouch && rq.source == Source::Clonal &&
            (J.endA == EndClass::Repeat || J.endB == EndClass::Repeat)) {
            dec.vouched = true;          // round 2: see JudgeRequest::clonalVouch
        } else {
            J.verdict = Verdict::PairsDeny;
        }
    }

    // ---- overlap rule: an exact overlap confirms a merge only through single-copy sequence
    // Allowed when some graph-unique 63-mer of the overlap lies in a single-copy unitig, or when
    // the graph joins the two single-copy ends at exactly the overlap (unitig neighbours: their
    // shared k-1 bases are never graph-unique themselves), or on PASS_EXACT.
    bool overlapRefused = false;
    if (isOverlap) {
        const int32_t ol = -n;
        if (ev_ && !ev_->abstained() && ol >= 63 && R.size() >= static_cast<size_t>(ol)) {
            for (int32_t at = 0; at + 63 <= ol && !X.overlapSingle; at += 25) {
                const uint32_t u = ev_->locateUnique(R, static_cast<size_t>(at));
                if (u != UINT32_MAX && ev_->mult(u).hardSingle) X.overlapSingle = true;
            }
        }
        const bool exactNeighbours = X.reachable && std::abs(ws.gmin - n) <= 10 &&
                                     J.endA == EndClass::Unique && J.endB == EndClass::Unique;
        overlapRefused = !X.overlapSingle && !exactNeighbours && J.verdict != Verdict::PassExact;
    }

    // ---- tier -----------------------------------------------------------------------
    if (isOverlap && X.overlapSingle) J.tier = Tier::A;
    else if (X.pairCall < 0 && J.verdict != Verdict::PairsDeny) J.tier = Tier::B;
    else if (J.verdict == Verdict::PassExact || J.verdict == Verdict::Resize) J.tier = Tier::C;
    else if (J.verdict == Verdict::PassWalk || J.verdict == Verdict::ContraPairsOk) J.tier = Tier::D;
    else J.tier = Tier::E;

    // ---- plasmid rule (C1d) ------------------------------------------------------------
    bool plasmidBlocked = false;
    if (rq.plasmidPass && !rq.judgeOnly) {
        const bool chrEnd = rq.chromosomalA || rq.chromosomalB;
        const bool pairLink = pairs_ && pairs_->usable() && X.originA != UINT32_MAX &&
                              X.originB != UINT32_MAX &&
                              pairs_->linkWeight(X.originA >> 1, X.originB >> 1) >= pairs_->supportFloor();
        const bool confirmed = X.pairCall < 0 || J.verdict == Verdict::PassExact || pairLink;
        plasmidBlocked = chrEnd || !confirmed;
    }

    // ---- admission ----------------------------------------------------------------------
    J.cls = classKey(J, X);
    const bool breakVerdict = J.verdict == Verdict::BreakDoubleUse ||
                              J.verdict == Verdict::ContraRefused || J.verdict == Verdict::PairsDeny;
    bool pairsDenyResize = false;
    if (J.verdict == Verdict::PairsDeny && ws.reachable && ws.gmin > n + 100) pairsDenyResize = true;
    if (admit_.loaded && admit_.upper.count(J.cls)) {
        J.pMisjoin = static_cast<float>(admit_.upper[J.cls]);
        J.admit = J.pMisjoin <= admit_.tauS ? Admit::Scaffold
                : J.pMisjoin <= admit_.tauG ? Admit::GenomeOnly : Admit::Break;
    } else if ((breakVerdict && !pairsDenyResize) || overlapRefused || plasmidBlocked) {
        J.admit = Admit::Break;
    } else if (J.verdict == Verdict::Silent) {
        J.admit = (J.endA == EndClass::Unique && J.endB == EndClass::Unique) ? Admit::Scaffold
                                                                              : Admit::GenomeOnly;
    } else {
        J.admit = Admit::Scaffold;
    }

    // ---- action ---------------------------------------------------------------------------
    bool keep = J.admit != Admit::Break;
    int32_t writeN = INT32_MIN;
    std::string why = keep ? "admit" : overlapRefused ? "overlap_not_single_copy"
                                     : plasmidBlocked ? "plasmid_rule" : "verdict";
    if (keep && (J.verdict == Verdict::Resize || pairsDenyResize)) {
        writeN = std::max(1, ws.gmin);
        why = "graph_sized";
    }
    const bool butt = rq.source == Source::JoinButt1 || rq.source == Source::LayoutButt1;
    const bool capJ = rq.source == Source::LayoutCap2000;
    if (butt) ++counters.buttSeen;
    if (capJ) ++counters.capSeen;
    if (isOverlap) ++counters.overlapSeen;
    // The 1-N butt is retired where nothing but the panel stands behind it. A butt the pairs
    // confirm (tier B) keeps its single N; one the graph pins to a single walk (PASS_EXACT) is
    // sized from that walk; the rest are broken (break) or sized from any walk (sized).
    if (keep && butt && writeN == INT32_MIN && cfg_.unsized != SeamConfig::Unsized::Butt) {
        if (X.pairCall < 0) {
            why = "butt_pairs_confirm";
        } else if (J.verdict == Verdict::PassExact && ws.gmin >= 1) {
            writeN = ws.gmin;
            why = "unsized_graph_exact";
        } else if (cfg_.unsized == SeamConfig::Unsized::Break) {
            keep = false;
            why = "unsized_break";
        } else if (ws.reachable && ws.gmin >= 1) {
            writeN = ws.gmin;
            why = "unsized_sized";
        } else if (!ws.reachable) {
            keep = false;
            why = "unsized_no_walk";
        }
    }
    if (keep && capJ && writeN == INT32_MIN && cfg_.cap == SeamConfig::Cap::Size) {
        bool agree = false;
        if (ws.reachable && !J.walks.empty() && !ws.classesCapped) {
            const int32_t lo = J.walks.front().len, hi = J.walks.back().len;
            agree = hi - lo <= std::max<int32_t>(100, lo / 20);
        }
        bool panelConsistent = false;
        for (const WalkClass& c : J.walks) {
            if (std::abs(c.len - rq.panelGap) <= std::max<int32_t>(500, rq.panelGap / 10)) panelConsistent = true;
        }
        const bool uniqueDead = J.verdict == Verdict::Silent &&
                                J.endA != EndClass::Repeat && J.endB != EndClass::Repeat;
        // C4 round 3 (WALKSIZE): the isolate's shortest walk class, isolated from the others (they loop through a repeat
        // the flanks share), sizes the capped gap instead of breaking it (dev Kp round 2: all 13 truth-checked cap_break
        // rows with such a walk were true adjacencies at exactly that walk).
        int32_t c1 = INT32_MIN;
        if (!agree && !(panelConsistent || uniqueDead) && cfg_.capWalkC1 && ws.reachable && !J.walks.empty()) {
            std::vector<int32_t> lens;
            std::vector<uint16_t> cnts;
            for (const WalkClass& c : J.walks) { lens.push_back(c.len); cnts.push_back(c.n); }
            c1 = isolatedShortestWalk(lens, cnts, ws.exhaustive, 5000, 20000);
        }
        if (agree) { writeN = std::max(1, J.walks.front().len); why = "cap_walk"; }
        else if (panelConsistent || uniqueDead) { writeN = std::max(1, rq.panelGap); why = "cap_panel"; }
        else if (c1 != INT32_MIN) { writeN = std::max(1, c1); why = "cap_walk_c1"; }
        else { keep = false; why = "cap_break"; }
    }

    // allowCloseGaps (C1e): never cement a model gap that is repeat-bounded, not graph-backed
    // or genome-only.
    if (isModelSource(rq.source)) {
        const bool graphBacked = J.verdict == Verdict::PassExact || J.verdict == Verdict::PassWalk ||
                                 J.verdict == Verdict::Resize;
        J.allowCloseGaps = graphBacked && J.endA != EndClass::Repeat && J.endB != EndClass::Repeat &&
                           J.admit != Admit::GenomeOnly;
    }

    // C4 round 2: a clonal exact-overlap join is judged at its overlap (claimedN = -overlap: the exit anchor slides back
    // past the shared bases) and written as a 1-N butt in the record (the genome view merges the overlap).
    const int32_t releaseN = isOverlap ? 0 : (rq.source == Source::Clonal && n < 0) ? 1 : n;
    if (rq.judgeOnly) {
        X.action = "record";
        J.writtenN = releaseN;
        ++counters.resolverJudged;
    } else {
        ++counters.judged;
        switch (J.verdict) {
            case Verdict::PassExact: ++counters.passExact; break;
            case Verdict::PassWalk: ++counters.passWalk; break;
            case Verdict::Resize: ++counters.resize; break;
            case Verdict::BreakDoubleUse: ++counters.breakDoubleUse; break;
            case Verdict::ContraPairsOk: ++counters.contraPairsOk; break;
            case Verdict::ContraRefused: ++counters.contraRefused; break;
            case Verdict::PairsDeny: ++counters.pairsDeny; break;
            case Verdict::Silent: ++counters.silent; break;
            case Verdict::Unanchored: ++counters.unanchored; break;
            case Verdict::Abstain: ++counters.abstain; break;
            default: break;
        }
        if (J.admit == Admit::Scaffold) ++counters.admitScaf;
        else if (J.admit == Admit::GenomeOnly) ++counters.admitGenomeOnly;
        else ++counters.admitBreak;
        if (overlapRefused) ++counters.overlapRefused;
        if (plasmidBlocked) ++counters.plasmidBlocked;
        std::string act = !keep ? "break" : writeN != INT32_MIN ? "resize:" + std::to_string(writeN) : "keep";
        if (cfg_.act()) {
            X.action = act;
            if (!keep) ++counters.broken;
            else if (writeN != INT32_MIN) ++counters.resized;
            if (butt && (!keep || writeN != INT32_MIN)) ++counters.buttRetired;
            if (capJ && (!keep || writeN != INT32_MIN)) ++counters.capRetired;
            dec.keep = keep;
            dec.writeN = writeN;
            J.writtenN = !keep ? -1 : writeN != INT32_MIN ? writeN : releaseN;
        } else {
            X.action = act == "keep" ? "keep" : "would_" + act;
            J.writtenN = releaseN;
        }
    }
    X.why = dec.vouched ? why + "+clonal_vouch" : why;
    ledger_.j.push_back(std::move(J));
    audit_.push_back(std::move(X));
    const Junction& jj = ledger_.j.back();
    if (!jj.flankL32.empty() && !jj.flankR32.empty() && jj.writtenN > 0) {
        const std::string k1 = jj.flankL32 + "|" + jj.flankR32;
        const std::string k2 = reverseComplement(jj.flankR32) + "|" + reverseComplement(jj.flankL32);
        flankAll_[k1].push_back({jj.id, '+'});
        if (k2 != k1) flankAll_[k2].push_back({jj.id, '-'});
    }
    counters.judgeMs += nowMs() - t0;
    return dec;
}

void SeamContext::clonalAdmit(uint32_t id, bool noFill, const std::string& why, int32_t mergeOverlap) {
    if (id >= ledger_.j.size()) return;
    Junction& J = ledger_.j[id];
    if (J.admit == Admit::Scaffold) J.admit = Admit::GenomeOnly;
    J.allowCloseGaps = false;
    J.noFill = noFill || mergeOverlap > 0;
    J.mergeOverlap = mergeOverlap > 0 ? mergeOverlap : 0;
    JunctionAudit& X = audit_[id];
    X.why = X.why.empty() ? why : X.why + "+" + why;
}

void SeamContext::clonalSize(uint32_t id, int32_t n) {
    if (id >= ledger_.j.size()) return;
    Junction& J = ledger_.j[id];
    J.clonalSized = true;
    JunctionAudit& X = audit_[id];
    if (n > 0 && n != J.writtenN) {
        X.why = X.why.empty() ? "clonal_resize" : X.why + "+clonal_resize";
        X.action = "resize:" + std::to_string(n);
        J.writtenN = n;
    } else {
        X.why = X.why.empty() ? "clonal_size_agree" : X.why + "+clonal_size_agree";
    }
}

void SeamContext::clonalNoFill(uint32_t id, const std::string& why) {
    if (id >= ledger_.j.size()) return;
    Junction& J = ledger_.j[id];
    J.noFill = true;
    J.allowCloseGaps = false;   // the legacy gap closer may not cement bases either
    JunctionAudit& X = audit_[id];
    X.why = X.why.empty() ? why : X.why + "+" + why;
}

void SeamContext::clonalBreak(uint32_t id, const std::string& why) {
    if (id >= ledger_.j.size()) return;
    Junction& J = ledger_.j[id];
    // The run leaves the records, so the exact-flank candidates of the locator must not name it.
    if (!J.flankL32.empty() && !J.flankR32.empty()) {
        for (const std::string& k : {J.flankL32 + "|" + J.flankR32,
                                     reverseComplement(J.flankR32) + "|" + reverseComplement(J.flankL32)}) {
            auto it = flankAll_.find(k);
            if (it == flankAll_.end()) continue;
            auto& v = it->second;
            v.erase(std::remove_if(v.begin(), v.end(),
                                   [&](const std::pair<uint32_t, char>& x) { return x.first == id; }), v.end());
            if (v.empty()) flankAll_.erase(it);
        }
    }
    J.writtenN = -1;
    J.admit = Admit::Break;
    J.allowCloseGaps = false;
    JunctionAudit& X = audit_[id];
    X.action = "break";
    X.why = why;
}

std::vector<std::pair<long, long>> SeamContext::nRuns(const std::string& s) const {
    std::vector<std::pair<long, long>> out;
    size_t i = 0;
    while (i < s.size()) {
        if (s[i] != 'N' && s[i] != 'n') { ++i; continue; }
        size_t j = i;
        while (j < s.size() && (s[j] == 'N' || s[j] == 'n')) ++j;
        out.emplace_back(static_cast<long>(i), static_cast<long>(j));
        i = j;
    }
    return out;
}

void SeamContext::recordResolverGaps(const std::vector<std::string>& pieces) {
    if (!gating()) return;
    beginBatch(pieces);
    for (size_t c = 0; c < pieces.size(); ++c) {
        const std::string& s = pieces[c];
        for (const auto& r : nRuns(s)) {
            const size_t a = static_cast<size_t>(r.first), b = static_cast<size_t>(r.second);
            const std::string left = s.substr(a > kWin ? a - kWin : 0, a > kWin ? kWin : a);
            const std::string right = s.substr(b, std::min(kWin, s.size() - b));
            JudgeRequest rq;
            rq.stage = "resolver";
            rq.left = &left;
            rq.right = &right;
            rq.claimedN = static_cast<int32_t>(b - a);
            rq.panelGap = rq.claimedN;
            rq.source = rq.claimedN == 100 ? Source::ResolverUnknown100 : Source::Resolver;
            rq.a = Port{static_cast<uint32_t>(c), true};
            rq.b = Port{static_cast<uint32_t>(c), false};
            rq.judgeOnly = true;
            judge(rq);
        }
    }
}

int SeamContext::probeScore(const std::string& sq, long a, long b, uint32_t id, char orient) const {
    // A probe counts when it occurs within +-150 bp of where the junction puts it (T03 restores
    // up to k-1 bases next to a resolver gap after closeGaps, shifting everything by that much).
    const long kTol = 150;
    const long n = static_cast<long>(sq.size());
    auto near = [&](long pos, const std::string& q) {
        const long lo = std::max(0L, pos - kTol), hi = std::min(n, pos + 63 + kTol);
        return hi - lo >= 63 &&
               sq.substr(static_cast<size_t>(lo), static_cast<size_t>(hi - lo)).find(q) != std::string::npos;
    };
    const JunctionAudit& X = audit_[id];
    int sc = 0;
    if (orient == '+') {
        for (const auto& p : X.probesL) sc += near(a - p.first - 63, p.second);
        for (const auto& p : X.probesR) sc += near(b + p.first, p.second);
    } else {
        for (const auto& p : X.probesL) sc += near(b + p.first, reverseComplement(p.second));
        for (const auto& p : X.probesR) sc += near(a - p.first - 63, reverseComplement(p.second));
    }
    return sc;
}

void SeamContext::locate(const std::vector<std::string>& seqs, std::vector<Loc>& loc,
                         std::vector<long>& owner,
                         std::vector<std::vector<std::pair<long, long>>>& runs) const {
    const size_t nj = ledger_.j.size();
    loc.assign(nj, Loc());
    runs.assign(seqs.size(), {});
    for (size_t i = 0; i < seqs.size(); ++i) runs[i] = nRuns(seqs[i]);
    const long kTol = 150;
    // ---- candidates: exact 32 bp flanks, and graph anchors beside an N-run (+-150 bp) ------
    std::map<std::pair<long, long>, std::map<std::pair<uint32_t, char>, bool>> cand;   // run -> (id, orient) -> exact
    for (size_t i = 0; i < seqs.size(); ++i) {
        const std::string& sq = seqs[i];
        for (size_t r = 0; r < runs[i].size(); ++r) {
            const long a = runs[i][r].first, b = runs[i][r].second;
            if (a < static_cast<long>(kFlank) || b + static_cast<long>(kFlank) > static_cast<long>(sq.size())) continue;
            auto it = flankAll_.find(sq.substr(static_cast<size_t>(a) - kFlank, kFlank) + "|" +
                                     sq.substr(static_cast<size_t>(b), kFlank));
            if (it == flankAll_.end()) continue;
            for (const auto& c : it->second) cand[{static_cast<long>(i), static_cast<long>(r)}][c] = true;
        }
    }
    std::unordered_map<uint64_t, std::vector<std::pair<uint32_t, char>>> amap;
    for (uint32_t id = 0; id < nj; ++id) {
        const JunctionAudit& X = audit_[id];
        if (ledger_.j[id].writtenN <= 0) continue;
        uint64_t fp;
        if (X.slideA >= 0 && fingerprint63(X.anchorSeqA, 0, fp)) amap[fp].push_back({id, 'A'});
        if (X.slideB >= 0 && fingerprint63(X.anchorSeqB, 0, fp)) amap[fp].push_back({id, 'B'});
    }
    if (!amap.empty()) {
        for (size_t i = 0; i < seqs.size(); ++i) {
            const std::string& sq = seqs[i];
            const auto& rr = runs[i];
            if (rr.empty()) continue;
            auto startNear = [&](long pred) -> long {
                auto it = std::lower_bound(rr.begin(), rr.end(), std::make_pair(pred, 0L));
                return (it != rr.end() && it->first <= pred + kTol) ? static_cast<long>(it - rr.begin()) : -1;
            };
            auto endNear = [&](long pred) -> long {
                auto it = std::upper_bound(rr.begin(), rr.end(), pred,
                                           [](long v, const std::pair<long, long>& x) { return v < x.second; });
                if (it == rr.begin()) return -1;
                --it;
                return (it->second >= pred - kTol) ? static_cast<long>(it - rr.begin()) : -1;
            };
            forEachFingerprint63(sq, [&](size_t pos, uint64_t fp) {
                auto it = amap.find(fp);
                if (it == amap.end()) return;
                const std::string here = sq.substr(pos, 63);
                for (const auto& hit : it->second) {
                    const JunctionAudit& X = audit_[hit.first];
                    const std::string& an = hit.second == 'A' ? X.anchorSeqA : X.anchorSeqB;
                    const long slide = hit.second == 'A' ? X.slideA : X.slideB;
                    const bool fwd = here == an;
                    if (!fwd && here != reverseComplement(an)) continue;
                    const long p = static_cast<long>(pos);
                    long r;
                    if (hit.second == 'A') r = fwd ? startNear(p + 63 + slide) : endNear(p - slide);
                    else r = fwd ? endNear(p - slide) : startNear(p + 63 + slide);
                    if (r < 0) continue;
                    auto& m = cand[{static_cast<long>(i), r}];
                    m.emplace(std::make_pair(hit.first, fwd ? '+' : '-'), false);
                }
            });
        }
    }
    // ---- score every candidate by its probes; keep mutual best pairs ------------------------
    struct Best { long seq = -1, run = -1; char orient = '.'; int score = 0; bool tie = false; bool exact = false; };
    std::vector<Best> bestOfJ(nj);
    std::map<std::pair<long, long>, Best> bestOfRun;   // seq,run -> junction as seq=id
    for (const auto& kv : cand) {
        const long i = kv.first.first, r = kv.first.second;
        const std::string& sq = seqs[static_cast<size_t>(i)];
        const long a = runs[static_cast<size_t>(i)][static_cast<size_t>(r)].first;
        const long b = runs[static_cast<size_t>(i)][static_cast<size_t>(r)].second;
        Best& br = bestOfRun[kv.first];
        for (const auto& c : kv.second) {
            const uint32_t id = c.first.first;
            const int sc = probeScore(sq, a, b, id, c.first.second);
            if (sc <= 0) continue;
            Best& bj = bestOfJ[id];
            if (sc > bj.score) bj = Best{i, r, c.first.second, sc, false, c.second};
            else if (sc == bj.score && (bj.seq != i || bj.run != r)) bj.tie = true;
            if (sc > br.score) br = Best{static_cast<long>(id), 0, c.first.second, sc, false, c.second};
            else if (sc == br.score && br.seq != static_cast<long>(id)) br.tie = true;
        }
    }
    for (uint32_t id = 0; id < nj; ++id) {
        const Best& bj = bestOfJ[id];
        if (bj.score == 0 || bj.tie) continue;
        auto it = bestOfRun.find({bj.seq, bj.run});
        if (it == bestOfRun.end() || it->second.tie || it->second.seq != static_cast<long>(id)) continue;
        const long a = runs[static_cast<size_t>(bj.seq)][static_cast<size_t>(bj.run)].first;
        const std::string& sq = seqs[static_cast<size_t>(bj.seq)];
        // exact: the 32 bp flanks still sit right against the N-run
        const bool exact = bj.exact && a >= static_cast<long>(kFlank) &&
                           (bj.orient == '+' ? sq.compare(static_cast<size_t>(a) - kFlank, kFlank, ledger_.j[id].flankL32) == 0
                                             : sq.compare(static_cast<size_t>(a) - kFlank, kFlank, reverseComplement(ledger_.j[id].flankR32)) == 0);
        loc[id] = Loc{bj.seq, bj.run, bj.orient, !exact, bj.score};
    }
    owner.clear();
    std::map<std::pair<long, long>, long> claim;
    for (uint32_t id = 0; id < nj; ++id) if (loc[id].seq >= 0) claim[{loc[id].seq, loc[id].run}] = id;
    for (size_t i = 0; i < seqs.size(); ++i) {
        for (size_t r = 0; r < runs[i].size(); ++r) {
            auto it = claim.find({static_cast<long>(i), static_cast<long>(r)});
            owner.push_back(it == claim.end() ? -1 : it->second);
        }
    }
}

std::vector<std::vector<RunOwner>> SeamContext::runOwners(const std::vector<std::string>& seqs) const {
    std::vector<Loc> loc;
    std::vector<long> owner;
    std::vector<std::vector<std::pair<long, long>>> runs;
    locate(seqs, loc, owner, runs);
    std::vector<std::vector<RunOwner>> out(seqs.size());
    size_t g = 0;
    for (size_t i = 0; i < seqs.size(); ++i) {
        for (const auto& r : runs[i]) {
            RunOwner o;
            o.start = r.first;
            o.end = r.second;
            const long id = g < owner.size() ? owner[g] : -1;
            ++g;
            if (id >= 0) {
                o.id = id;
                o.orient = loc[static_cast<size_t>(id)].orient;
            }
            out[i].push_back(o);
        }
    }
    return out;
}

bool SeamContext::beforeCloseGaps(const std::vector<std::string>& seqs, std::vector<uint8_t>& allow) {
    std::vector<Loc> loc;
    std::vector<std::vector<std::pair<long, long>>> runs;
    locate(seqs, loc, gapJunction_, runs);
    allow.assign(gapJunction_.size(), 1);
    for (size_t g = 0; g < gapJunction_.size(); ++g) {
        const long id = gapJunction_[g];
        if (id < 0) continue;
        ++counters.closegapsSeen;
        audit_[static_cast<size_t>(id)].gapSeenAtClose = true;
        if (!ledger_.j[static_cast<size_t>(id)].allowCloseGaps) { allow[g] = 0; ++counters.closegapsDisallowed; }
    }
    return cfg_.act() && cfg_.closeGapsPolicy;
}

void SeamContext::afterCloseGaps(const std::vector<std::string>& seqs, size_t blocked) {
    counters.closegapsBlocked = blocked;
    std::vector<Loc> loc;
    std::vector<long> owner;
    std::vector<std::vector<std::pair<long, long>>> runs;
    locate(seqs, loc, owner, runs);
    std::vector<char> present(ledger_.j.size(), 0);
    for (long id : owner) if (id >= 0) present[static_cast<size_t>(id)] = 1;
    for (long id : gapJunction_) {
        if (id >= 0 && !present[static_cast<size_t>(id)]) {
            audit_[static_cast<size_t>(id)].closedByGapfill = true;
            ++counters.closegapsClosed;
        }
    }
}

bool SeamContext::writeTsv(const std::string& path, const std::vector<std::string>& outSeqs,
                           const std::vector<std::string>& outNames, std::string& err) {
    // ---- where did each written N-run end up? anchors, then flanks -----------------------
    {
        std::vector<Loc> loc;
        std::vector<long> owner;
        std::vector<std::vector<std::pair<long, long>>> runs;
        locate(outSeqs, loc, owner, runs);
        for (size_t id = 0; id < ledger_.j.size(); ++id) {
            JunctionAudit& X = audit_[id];
            if (loc[id].seq < 0 || X.closedByGapfill) continue;
            const std::string& nm = outNames[static_cast<size_t>(loc[id].seq)];
            const size_t sp = nm.find_first_of(" \t");
            const auto& r = runs[static_cast<size_t>(loc[id].seq)][static_cast<size_t>(loc[id].run)];
            X.finalRecord = sp == std::string::npos ? nm : nm.substr(0, sp);
            X.finalStart = r.first;
            X.finalLen = r.second - r.first;
            X.finalOrient = loc[id].orient;
            X.finalFuzzy = loc[id].fuzzy;
        }
    }

    std::FILE* f = std::fopen(path.c_str(), "w");
    if (!f) { err = "cannot write " + path; return false; }
    std::fprintf(f, "# om2_seams.tsv -- Organism Model 2.0 C1 junction ledger (mode=%s unsized=%s cap=%s "
                    "pairs=%d plasmid_rule=%d closegaps_policy=%d admit=%s theta=%.3f k=%d)\n",
                 SeamConfig::modeName(cfg_.mode), SeamConfig::unsizedName(cfg_.unsized),
                 SeamConfig::capName(cfg_.cap), cfg_.pairs ? 1 : 0, cfg_.plasmidRule ? 1 : 0,
                 cfg_.closeGapsPolicy ? 1 : 0, admit_.loaded ? cfg_.admitPath.c_str() : "v0",
                 ev_ ? ev_->theta() : 0.0, ev_ ? ev_->k() : 0);
    std::fprintf(f, "id\tstage\tsource\tverdict\tgraph_verdict\ttier\tadmit\taction\twhy\tclaimed_n\t"
                    "written_n\tpanel_gap\tgmin\twalk_classes\tdistinct\tcapped\thairpin\t"
                    "inner_unique\tinner_unique_placed\tuA\tuB\tu_other\tanchorA\tanchorB\tendA\tendB\t"
                    "copyA\tcopyB\tpair_testable\tpair_abstain\tpair_k\tpair_lambda\tpair_ctrlA\t"
                    "pair_ctrlB\tpair_log10LR\tpair_call\tpair_contra\toverlap_single\torigin_a\t"
                    "origin_b\tpanel_support\tpanel_pairs\tallow_closegaps\tseen_at_closegaps\t"
                    "closed_by_gapfill\tfinal_record\tfinal_start\tfinal_len\tfinal_orient\t"
                    "final_fuzzy\tflankL32\tflankR32\tcls\n");
    for (size_t id = 0; id < ledger_.j.size(); ++id) {
        const Junction& J = ledger_.j[id];
        const JunctionAudit& X = audit_[id];
        std::string wc;
        for (const WalkClass& c : J.walks) {
            if (!wc.empty()) wc += ";";
            wc += std::to_string(c.len) + "x" + std::to_string(c.n);
        }
        if (wc.empty()) wc = ".";
        const std::string gminStr = X.reachable ? std::to_string(J.gmin) : std::string(".");
        auto port = [](uint32_t p) { return p == UINT32_MAX ? std::string(".") : std::to_string(p >> 1) + ((p & 1u) ? "R" : "L"); };
        std::fprintf(f, "%zu\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%s\t%s\t%u\t%d\t%d\t%u\t%u\t%u\t%u\t%u\t%s\t%s\t"
                        "%s\t%s\t%.3f\t%.3f\t%d\t%d\t%u\t%.2f\t%.0f\t%.0f\t%.2f\t%d\t%u\t%d\t%s\t%s\t%u\t%u\t%d\t%d\t%d\t"
                        "%s\t%ld\t%ld\t%c\t%d\t%s\t%s\t%s\n",
                     id, X.stage.c_str(), sourceName(J.source), verdictName(J.verdict),
                     verdictName(X.graphVerdict), tierName(J.tier), admitName(J.admit), X.action.c_str(),
                     X.why.c_str(), J.claimedN, J.writtenN, X.panelGap,
                     gminStr.c_str(),
                     wc.c_str(), X.distinct, X.capped ? 1 : 0, J.hairpin ? 1 : 0, X.innerUnique,
                     X.innerUniquePlaced, J.uA, J.uB, X.uOther,
                     X.anchorA.empty() ? "." : X.anchorA.c_str(), X.anchorB.empty() ? "." : X.anchorB.c_str(),
                     endClassName(J.endA), endClassName(J.endB), J.copyA, J.copyB,
                     X.pairTestable ? 1 : 0, X.pairAbstain ? 1 : 0, J.pairK, J.pairLambda, X.ctrlA,
                     X.ctrlB, X.pairLog10LR, X.pairCall, J.pairContra, X.overlapSingle ? 1 : 0,
                     port(X.originA).c_str(), port(X.originB).c_str(), J.panelSupport, J.panelGenomes,
                     J.allowCloseGaps ? 1 : 0, X.gapSeenAtClose ? 1 : 0, X.closedByGapfill ? 1 : 0,
                     X.finalRecord.empty() ? "." : X.finalRecord.c_str(), X.finalStart, X.finalLen,
                     X.finalOrient, X.finalFuzzy ? 1 : 0, J.flankL32.c_str(), J.flankR32.c_str(),
                     J.cls.c_str());
    }
    std::fclose(f);
    return true;
}

void SeamContext::printCounters(const SeamContext* ctx, std::FILE* f) {
    const EvidenceIndex* ev = ctx ? ctx->ev_.get() : nullptr;
    const PairIndex* pi = ctx ? ctx->pairs_.get() : nullptr;
    std::fprintf(f,
                 "[om2-evidence] enabled=%d nodes=%zu components=%zu single_copy=%zu repeat=%zu "
                 "uncertain=%zu sc_fraction=%.3f theta=%.2f anchored=%zu unanchored=%zu "
                 "abstain_tangled=%d pair_index=%d reads_anchored=%zu pairs_spanning=%zu "
                 "pairs_control=%zu build_ms=%.0f pair_ms=%.0f\n",
                 ev ? 1 : 0, ev ? ev->liveNodes : 0, ev ? ev->components : 0, ev ? ev->singleCopy : 0,
                 ev ? ev->repeatNodes : 0, ev ? ev->uncertainNodes : 0, ev ? ev->scFraction : 0.0,
                 ev ? ev->theta() : 0.0, ev ? ev->anchorsOk : 0, ev ? ev->anchorsFailed : 0,
                 ev && ev->abstained() ? 1 : 0, pi && pi->usable() ? 1 : 0, pi ? pi->readsAnchored : 0,
                 pi ? pi->pairsSpanning : 0, pi ? pi->pairsControl : 0, ev ? ev->buildMs : 0.0,
                 ctx ? ctx->counters.pairBuildMs : 0.0);
    const SeamCounters z;
    const SeamCounters& c = ctx ? ctx->counters : z;
    std::fprintf(f,
                 "[om2-seam] enabled=%d mode=%s unsized=%s cap=%s pairs=%d plasmid_rule=%d "
                 "closegaps_policy=%d admit=%s judged=%zu resolver_judged=%zu pass_exact=%zu "
                 "pass_walk=%zu resize=%zu break_double_use=%zu contra_pairs_ok=%zu "
                 "contra_refused=%zu pairs_deny=%zu silent=%zu unanchored=%zu abstain=%zu "
                 "butt_seen=%zu butt_retired=%zu cap_seen=%zu cap_retired=%zu overlap_seen=%zu "
                 "overlap_refused=%zu plasmid_blocked=%zu plasmid_depth_excluded=%zu "
                 "pair_testable=%zu pair_abstain=%zu pair_confirm=%zu pair_deny=%zu admit_scaf=%zu "
                 "admit_genome_only=%zu admit_break=%zu broken=%zu resized=%zu closegaps_seen=%zu "
                 "closegaps_disallowed=%zu closegaps_blocked=%zu closegaps_closed=%zu judge_ms=%.0f\n",
                 ctx && ctx->gating() ? 1 : 0, ctx ? SeamConfig::modeName(ctx->cfg_.mode) : "off",
                 ctx ? SeamConfig::unsizedName(ctx->cfg_.unsized) : "-",
                 ctx ? SeamConfig::capName(ctx->cfg_.cap) : "-", ctx && ctx->cfg_.pairs ? 1 : 0,
                 ctx && ctx->cfg_.plasmidRule ? 1 : 0, ctx && ctx->cfg_.closeGapsPolicy ? 1 : 0,
                 ctx && ctx->admit_.loaded ? "table" : "v0", c.judged, c.resolverJudged, c.passExact,
                 c.passWalk, c.resize, c.breakDoubleUse, c.contraPairsOk, c.contraRefused, c.pairsDeny,
                 c.silent, c.unanchored, c.abstain, c.buttSeen, c.buttRetired, c.capSeen, c.capRetired,
                 c.overlapSeen, c.overlapRefused, c.plasmidBlocked, c.plasmidDepthExcluded,
                 c.pairTestable, c.pairAbstain, c.pairConfirm, c.pairDeny, c.admitScaf,
                 c.admitGenomeOnly, c.admitBreak, c.broken, c.resized, c.closegapsSeen,
                 c.closegapsDisallowed, c.closegapsBlocked, c.closegapsClosed, c.judgeMs);
}

}  // namespace om2
}  // namespace ts
