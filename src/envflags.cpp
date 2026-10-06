#include "envflags.h"

#include <algorithm>
#include <cctype>
#include <cerrno>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <initializer_list>
#include <utility>
#include <vector>
#include <unistd.h>

namespace ts { namespace env {
namespace {

constexpr double kAny = 1.7976931348623157e308;   // DBL_MAX: any finite value
constexpr long long kIntMin = INT_MIN;
constexpr long long kIntMax = INT_MAX;

constexpr Spec I(const char* n, long long lo, long long hi) { return {n, Kind::Integer, lo, hi, 0.0, 0.0, nullptr}; }
constexpr Spec R(const char* n, double lo, double hi) { return {n, Kind::Real, 0, 0, lo, hi, nullptr}; }
constexpr Spec K(const char* n, Kind k) { return {n, k, 0, 0, 0.0, 0.0, nullptr}; }
constexpr Spec C(const char* n, const char* words) { return {n, Kind::Choice, 0, 0, 0.0, 0.0, words}; }
// build_v3: a Choice whose word list holds the token "#" also accepts a base-10 integer in
// [lo, hi] (TESSERACT_PREFIX_MIN_BODY=<bp>|reach, TESSERACT_PAIR_ANCHORED_PREFIX=<m>|bar).
constexpr Spec CI(const char* n, const char* words, long long lo, long long hi) {
    return {n, Kind::Choice, lo, hi, 0.0, 0.0, words};
}

// Ranges. An integer flag the release read with atoi accepts exactly the values atoi was
// defined on (int), narrowed to >= 0 where the result was cast to an unsigned type, which
// is where -1 used to wrap to SIZE_MAX or 4294967295. A real flag accepts any finite
// number, narrowed where the result is later converted to an integer type and a larger
// value would make that conversion undefined. Defaults live at the call sites.
const Spec kTable[] = {
    // k-mer ladder, counting and carry-over (assembler.cpp)
    R("TESSERACT_QC_MIN_KCOV", -kAny, kAny),                 // <= 0 disables
    I("TESSERACT_QC_MIN_RUNGS", 0, kIntMax),                 // size_t
    K("TESSERACT_RC_DOVETAIL", Kind::Switch),                // default on
    I("TESSERACT_DENSE_LADDER", kIntMin, kIntMax),           // int; <= 0 is the shortest ladder
    K("TESSERACT_LEGACY_LADDER", Kind::Switch),
    R("TESSERACT_MIN_KMER_MASS", -kAny, kAny),
    I("TESSERACT_CARRY_WEIGHT", 0, kIntMax),                 // uint32_t
    R("TESSERACT_CARRY_BOOST_BELOW", -kAny, kAny),
    I("TESSERACT_CARRY_WEIGHT_LOW", 0, kIntMax),             // uint32_t
    K("TESSERACT_QC_CARRY", Kind::Switch),
    K("TESSERACT_TRUSTED_CARRY", Kind::Switch),
    K("TESSERACT_QC_INSERT", Kind::Switch),
    // gap closing (assembler.cpp, gap_evidence.cpp)
    I("TESSERACT_GAPCLOSE", 0, kIntMax),                     // size_t
    I("TESSERACT_GAPCLOSE_VOTES", 0, kIntMax),               // uint32_t
    K("TESSERACT_GAP_ORIENTED", Kind::Binary),
    K("TESSERACT_GAP_SEQUENCE_COMPETITION", Kind::Binary),
    // read correction, coverage and polishing
    K("TESSERACT_EC_REQUIRE_UNIQUE_BEST", Kind::Binary),
    K("TESSERACT_MATE_RESCUE", Kind::Binary),
    K("TESSERACT_OBSERVED_GRAPH_COVERAGE", Kind::Binary),
    R("TESSERACT_POLISH_FRACTION", -kAny, kAny),
    K("TESSERACT_POLISH_ORIGINAL_QUALITY", Kind::Binary),
    R("TESSERACT_MAPPOLISH_FRACTION", -kAny, kAny),
    I("TESSERACT_MAPPOLISH_DEPTH", kIntMin, kIntMax),        // int
    // graph simplification (graph.cpp)
    R("TESSERACT_TIP_LEN_MULT", -1e15, 1e15),                // times k, then cast to size_t
    R("TESSERACT_TIP_RATIO", -kAny, kAny),
    R("TESSERACT_TIP_ABS_MULT", -1e15, 1e15),                // times k, then cast to size_t
    K("TESSERACT_FITTED_EC", Kind::Switch),
    R("TESSERACT_CHIMERA_FACTOR", -kAny, kAny),
    R("TESSERACT_CHIMERA_FLOOR", -kAny, kAny),
    R("TESSERACT_EC_LEN_MULT", -1e15, 1e15),                 // times k, then cast to size_t
    R("TESSERACT_LOCAL_WEAK", -kAny, kAny),
    I("TESSERACT_JOIN_DEADENDS", 0, kIntMax),                // size_t
    K("TESSERACT_GRAPH_PHASES", Kind::Presence),             // timing log only
    K("TESSERACT_GF_DEBUG", Kind::Presence),                 // gap-fill log only
    // paired resolver (resolve.cpp)
    K("TESSERACT_WEIGHTED_RESOLVER_COVERAGE", Kind::Binary), // default on
    K("TESSERACT_WEIGHTED_ELIGIBLE_COVERAGE", Kind::Binary),
    K("TESSERACT_ROUTE_DISTANCE", Kind::Binary),
    K("TESSERACT_EXACT_READ_THREADS", Kind::Binary),         // default on
    R("TESSERACT_MATCH_DOMINANCE", -kAny, kAny),
    K("TESSERACT_PREFIX_SNP_BUBBLES", Kind::Switch),
    K("TESSERACT_REQUIRE_SUPPORT_SINGLE", Kind::Switch),
    K("TESSERACT_JOIN_TRACE", Kind::Switch),
    K("TESSERACT_NO_UNSPANNED_FALLBACK", Kind::Switch),
    I("TESSERACT_MIN_FALLBACK_DEST", 0, kIntMax),            // size_t; 0 disables
    K("TESSERACT_EXCLUDE_SHARED_REPEAT_SUPPORT", Kind::Switch),
    K("TESSERACT_SHARED_SUPPORT_AUDIT", Kind::Switch),
    R("TESSERACT_SCAF_SUPPORT", -kAny, kAny),                // <= 0 keeps the depth-scaled bar
    R("TESSERACT_COMMON_PREFIX", -1e18, 1e18),               // cast to size_t; <= 0 disables
    K("TESSERACT_OWNED_ANCHOR_PREFIX", Kind::Binary),
    K("TESSERACT_NO_PLASMID_VOUCH", Kind::One),
    K("TESSERACT_DEBUG_RESOLVE", Kind::Presence),            // log only
    // build_v3 campaign packages (combo2, build_v2_merge package.patch), all default off.
    // joinrule L4/XO1/L7 + FB-M (resolve.cpp)
    K("TESSERACT_ONE_SIDED_TIE", Kind::Switch),
    I("TESSERACT_ONE_SIDED_TIE_LONG", 0, kIntMax),           // size_t; default 2000
    K("TESSERACT_TIE_SINGLE_EXIT", Kind::Switch),
    I("TESSERACT_TIE_SINGLE_EXIT_LONG", 0, kIntMax),         // size_t; unset = no gate
    K("TESSERACT_FALLBACK_MATCHING_ONLY", Kind::Switch),
    // ENDS terminal extension (resolve.cpp)
    CI("TESSERACT_PREFIX_MIN_BODY", "reach|#", 0, kIntMax),   // bp, or the insert reach; 0 = off
    K("TESSERACT_CERTIFIED_PREFIX", Kind::Switch),
    R("TESSERACT_CERTIFIED_PREFIX_FLOW", -kAny, kAny),       // default 0.8
    K("TESSERACT_PREFIX_NO_REVISIT", Kind::Switch),
    CI("TESSERACT_PAIR_ANCHORED_PREFIX", "bar|#", 0, kIntMax),// mates, or the resolver's linkBar; 0 = off
    C("TESSERACT_PAIR_ANCHORED_MODE", "node|pos"),           // node is the default
    // GAP: aliases of the single implementations kept in build_v3 (see BUILD_NOTES.md)
    K("TESSERACT_GAP_KEEP_FLANK", Kind::Binary),             // = TESSERACT_FIX_GAP_FLANK (T03)
    I("TESSERACT_GAPFILL_BACKOFF", 0, 200),                  // bases; T11's back-off, 0 = off
    // BRIDGE (assembler.cpp, dropout_bridge.cpp)
    K("TESSERACT_DROPOUT_BRIDGE", Kind::Binary),
    I("TESSERACT_DROPOUT_BRIDGE_MIN_READS", 2, 10),
    // Q40 boundary-safe trim: alias of TESSERACT_FIX_BOUNDARY_SAFE_TRIM (T10)
    K("TESSERACT_BOUNDARY_SAFE_TRIM", Kind::Binary),
    // build_v3 defect fixes (combo3), all default off; TESSERACT_FIXES=1 turns on every fix
    // whose own variable is unset, and TESSERACT_FIX_<NAME>=0 keeps one fix off under it.
    K("TESSERACT_FIXES", Kind::Binary),
    I("TESSERACT_FIX_CARRY_READ_GATE", 0, 2),                // T01: 0 dry run, 1 either flank, 2 both
    K("TESSERACT_FIX_GAP_NOMINATOR", Kind::Binary),          // T08
    K("TESSERACT_FIX_HAIRPIN_KEEP", Kind::Binary),           // T28
    K("TESSERACT_FIX_SIMPLIFY_COUNT_ALL", Kind::Binary),     // T29
    K("TESSERACT_FIX_EC_CAP_MASK", Kind::Binary),            // T30
    K("TESSERACT_FIX_POLISH_SKIP_N", Kind::Binary),          // T02
    K("TESSERACT_FIX_GAPFILL_STRICT_BUDGET", Kind::Binary),  // T09
    I("TESSERACT_FIX_GAPFILL_BACKOFF", 0, 200),              // T11: 1 = 60 bp, 2..200 = bp
    K("TESSERACT_FIX_GAPFILL_SKIP_INPUT_N", Kind::Binary),   // T27
    K("TESSERACT_FIX_BOUNDARY_SAFE_TRIM", Kind::Binary),     // T10
    K("TESSERACT_FIX_TRIM_COPY_GUARD", Kind::Binary),        // T22
    K("TESSERACT_FIX_SPLIT_POSTPROCESS", Kind::Binary),      // T25
    K("TESSERACT_FIX_AGP_GFA_V2", Kind::Binary),             // T26
    K("TESSERACT_FIX_GAP_FLANK", Kind::Binary),              // T03
    K("TESSERACT_FIX_REVISIT_GUARD", Kind::Binary),          // T07
    K("TESSERACT_FIX_SCAFFOLD_CYCLE", Kind::Binary),         // T13
    K("TESSERACT_FIX_GAP_ESTIMATE", Kind::Binary),           // T20
    K("TESSERACT_FIX_TRUNC_GUARD", Kind::Binary),            // T21
    K("TESSERACT_FIX_COV_CONTRIB", Kind::Binary),            // T24
    K("TESSERACT_FIX_ROUTE_ORDER", Kind::Binary),            // T31
    K("TESSERACT_FIX_MIRROR_ROUTE", Kind::Binary),           // R6
    // organism model, layout and replicon calls
    R("TESSERACT_MODEL_MIN_PANEL", 0.0, 4294967295.0),       // cast to uint32_t
    R("TESSERACT_MODEL_MIN_FRACTION", -kAny, kAny),
    R("TESSERACT_MODEL_MIN_PAIRS", 0.0, 1e18),               // cast to size_t
    R("TESSERACT_MODEL_MAX_OVERLAP", -2147483648.0, 2147483647.0),  // cast to int32_t
    R("TESSERACT_MODEL_MIN_OVERLAP", 0.0, 1e18),             // cast to size_t
    I("TESSERACT_MODEL_ROUNDS", kIntMin, kIntMax),           // clamped to 1..32 where used
    C("TESSERACT_MODEL_MATCH", "greedy|mutual"),             // mutual is the default
    K("TESSERACT_IS_VETO_CHR", Kind::Binary),                // default on
    K("TESSERACT_PLASMID_COMPLETE", Kind::Binary),
    I("TESSERACT_COMEMBER_WEAK", LLONG_MIN, LLONG_MAX),      // <= 0 keeps the default
    I("TESSERACT_COMEMBER_STRONG", LLONG_MIN, LLONG_MAX),    // <= 0 keeps the default
    R("TESSERACT_COMEMBER_MIN", -kAny, kAny),                // outside (0, 1] keeps the default
    K("TESSERACT_HUB_GROUP", Kind::Switch),
    R("TESSERACT_HUB_MIN", -kAny, kAny),                     // outside (0, 1] keeps the default
    K("TESSERACT_DEBUG_ORGANISM", Kind::Presence),           // log only
    K("TESSERACT_DEBUG_LAYOUT", Kind::Presence),             // log only
    K("TESSERACT_MODEL_AUTHOR", Kind::Presence),             // unlocks --model
    K("TESSERACT_MODEL_DIR", Kind::Text),
    K("TESSERACT_JOIN_DUMP", Kind::Text),
    K("TESSERACT_PLASMID_CLUSTERS", Kind::Text),
    K("TESSERACT_DEV_FORK_BATCH", Kind::Text),
    // Organism Model 2.0, C1 SEAM (om2_ledger.h); all unset = release behaviour
    K("TESSERACT_OM2_EVIDENCE", Kind::Binary),
    C("TESSERACT_OM2_SEAM", "0|audit|act"),
    C("TESSERACT_OM2_UNSIZED", "butt|break|sized"),
    C("TESSERACT_OM2_CAP", "keep|size"),
    K("TESSERACT_OM2_PAIRS", Kind::Binary),
    K("TESSERACT_OM2_PLASMID_RULE", Kind::Binary),
    K("TESSERACT_OM2_CLOSEGAPS_POLICY", Kind::Binary),
    K("TESSERACT_OM2_ADMIT", Kind::Text),
    I("TESSERACT_OM2_TANGLE_MAX", 1, kIntMax),
    // Organism Model 2.0, C2 close (om2_close.cpp, om2_alloc.cpp), all default off
    K("TESSERACT_OM2_CLOSE", Kind::Binary),                  // stage [4c/7]
    C("TESSERACT_OM2_FILL", "none|genome|contig"),           // genome is the default when on
    C("TESSERACT_OM2_ALLOC", "off|consensus|phased|prior"),  // phased is the default when on
    K("TESSERACT_OM2_CIRC", Kind::Binary),
    K("TESSERACT_OM2_RRN_PRIOR", Kind::Text),                // <org>.om2rrn, md5-pinned to the model
    K("TESSERACT_OM2_DNAA", Kind::Text),                     // <org>.om2dnaa, md5-pinned to the model
    K("TESSERACT_OM2_FLOW", Kind::Binary),
    // Organism Model 2.0, component C3 (om2_output.cpp, organism_detect.cpp); all default off
    K("TESSERACT_OM2_OUTPUT", Kind::Binary),                 // genome/ owner view, report.json om2 block
    K("TESSERACT_OM2_AGP_EVIDENCE", Kind::Binary),           // per-gap linkage evidence in scaffolds.agp
    C("TESSERACT_OM2_DETECT", "off|warn|gate"),              // check --organism against the reads
    K("TESSERACT_OM2_DETECT_SKETCH", Kind::Text),            // default <model dir>/om2detect.sketch
    // Organism Model 2.0, component C4 CLONAL (om2_clonal.cpp); all default off
    K("TESSERACT_OM2_CLONAL", Kind::Binary),                 // stage [4b2/7] + junction labels
    I("TESSERACT_OM2_CLONAL_K", 1, 32),                      // nearest relatives consulted (5)
    I("TESSERACT_OM2_CLONAL_KMIN", 1, 32),                   // relatives that must agree (3)
    K("TESSERACT_OM2_CLONAL_LAYOUT", Kind::Binary),          // propose joins (1)
    K("TESSERACT_OM2_CLONAL_BREAK", Kind::Binary),           // break contradicted panel-only joins (1)
    C("TESSERACT_OM2_CLONAL_NONPOS", "fill|gap"),            // non-positional elements with bases (fill)
    R("TESSERACT_OM2_CLONAL_DMAX", 0.0, 1.0),                // distance gate for clonal joins (1)
    R("TESSERACT_OM2_CLONAL_CONFIDENT", 0.5, 1.0),           // claim threshold (0.99)
    R("TESSERACT_OM2_CLONAL_CONF_DMAX", 0.0, 1.0),           // distance gate for `confident` (1)
    K("TESSERACT_OM2_CLONAL_CAL", Kind::Text),               // dev-fitted calibration table
    I("TESSERACT_OM2_CLONAL_MINMARK", 1, 32),                // markers placing an end (3)
    K("TESSERACT_OM2_CLONAL_NRP", Kind::Text),               // nearest-relative plasmid sidecar (<org>.om2nrp)
    C("TESSERACT_OM2_CLONAL_WEIGHT", "distance|count"),      // relatives' votes (distance)
    R("TESSERACT_OM2_CLONAL_DSCALE", 1e-6, 0.1),             // vote weight distance scale (0.0005)
    K("TESSERACT_OM2_CLONAL_OVERLAP", Kind::Binary),         // round 2: exact-overlap joins (1)
    K("TESSERACT_OM2_CLONAL_VOUCH", Kind::Binary),           // round 2: relatives vouch over C1 pair DENY at repeat ends (1)
    K("TESSERACT_OM2_CLONAL_BRACKET", Kind::Binary),         // round 2: walk / relatives size bracket (0)
    K("TESSERACT_OM2_CLONAL_RESIZE", Kind::Binary),          // round 2: clonal size check of panel-sized joins (1)
    K("TESSERACT_OM2_CLONAL_WALKSIZE", Kind::Binary),        // round 3: data-driven gap sizing from the isolate's walks (0)
    K("TESSERACT_OM2_CLONAL_KC2", Kind::Binary),             // round 3: KC2, non-positional elements never filled (0)
    K("TESSERACT_OM2_CLONAL_CONFRULE", Kind::Binary),        // round 3: count-based `confident` rule (0)
    I("TESSERACT_OM2_CLONAL_DISCORD", 0, 1000),              // round 3: CONFRULE discord threshold (2)
    K("TESSERACT_OM2_CLONAL_INDEL", Kind::Binary),           // round 3b: end placements tolerate a 250 kb indel (0)
    K("TESSERACT_OM2_CLONAL_UNVERIFIED", Kind::Binary),      // round 3b: break panel-sized layout gaps nobody sizes (0)
    K("TESSERACT_OM2_CLONAL_CAPWALK", Kind::Binary),         // round 3: C1 cap rule sized at the isolated shortest walk (0)
    K("TESSERACT_OM2_CLONAL_STRICT", Kind::Binary),          // round 3c: graph-first emission rule of model junctions (0)
    K("TESSERACT_OM2_CLONAL_STRICT_NR0", Kind::Binary),      // round 3c: STRICT also breaks junctions no relative places (0)
    K("TESSERACT_OM2_CLONAL_ALLOW_CONFIDENT", Kind::Binary), // 1.5: write the `confident` claim (0: KC1, written `supported`)
    // 1.5 layout-only genome view (om2_output.cpp); default off
    K("TESSERACT_OM2_LAYOUT_ONLY", Kind::Binary),            // join bases only at isolate-confirmed junctions
    C("TESSERACT_OM2_LAYOUT_CONFIRM", "isolate|c1pass"),     // what confirms a junction (isolate)
    // Phase 2 (om2/design/EVAL_PLAN_P2.md), W1; all default off, unset = release behaviour
    // emit-A (p2_emit.h)
    K("TESSERACT_P2_LIBGUARD", Kind::Binary),                // R2: library orientation / insert-model guard
    C("TESSERACT_P2_CIRC", "off|close|verify"),              // R3: circle closure (close) + verified claims (verify)
    K("TESSERACT_P2_SPIKEIN", Kind::Binary),                 // F5: PhiX174 spike-in screen, `_spikein` label
    // emit-B (p2_emitb.h)
    K("TESSERACT_P2_ENDS", Kind::Binary),                    // ends.tsv (report-only)
    K("TESSERACT_P2_TIPS_LOWERCASE", Kind::Binary),          // lower-case unsupported tips (case-only)
    K("TESSERACT_P2_SELF_QA", Kind::Binary),                 // k-mer self-QA, report.json p2.self_qa (report-only)
    K("TESSERACT_P2_PROVENANCE", Kind::Binary),              // provenance manifest (report-only)
    K("TESSERACT_P2_DETECT_REPORT", Kind::Binary),           // organism detection as a report (report-only)
    K("TESSERACT_P2_DETECT_SKETCH", Kind::Text),             // its sketch (else OM2_DETECT_SKETCH / MODEL_DIR)
    // read by tesseract-eskape, tesseract-klebsiella and tesseract-get-models
    K("TESSERACT_ASM", Kind::External),
    K("TESSERACT_KP_MODEL", Kind::External),
    K("TESSERACT_MODEL_URL", Kind::External),
    K("TESSERACT_TIME", Kind::External),
};

bool startsWithSpace(const char* v) { return std::isspace(static_cast<unsigned char>(*v)) != 0; }

bool parseInteger(const char* v, long long lo, long long hi, long long& out, std::string& why) {
    if (!*v) { why = "empty value"; return false; }
    if (startsWithSpace(v)) { why = "leading whitespace"; return false; }
    errno = 0;
    char* end = nullptr;
    const long long x = std::strtoll(v, &end, 10);
    if (end == v || *end != '\0') { why = "not a base-10 integer"; return false; }
    if (errno == ERANGE) { why = "out of range"; return false; }
    if (x < lo || x > hi) { why = "out of range"; return false; }
    out = x;
    return true;
}

std::string shortest(double x) {
    char buf[64];
    for (int p = 1; p <= 17; ++p) {
        std::snprintf(buf, sizeof buf, "%.*g", p, x);
        if (std::strtod(buf, nullptr) == x) break;
    }
    return buf;
}

bool isChoice(const char* words, const char* v, long long lo = 0, long long hi = -1) {
    const size_t n = std::strlen(v);
    if (n == 0) return false;
    for (const char* p = words; *p;) {
        const char* bar = std::strchr(p, '|');
        const size_t len = bar ? static_cast<size_t>(bar - p) : std::strlen(p);
        if (len == n && std::strncmp(p, v, n) == 0) return true;
        if (len == 1 && *p == '#' && lo <= hi) {   // an integer in [lo, hi]
            long long x = 0;
            std::string whyNot;
            if (parseInteger(v, lo, hi, x, whyNot)) return true;
        }
        if (!bar) break;
        p = bar + 1;
    }
    return false;
}

const Spec& lookup(const char* name, std::initializer_list<Kind> kinds, const char* reader) {
    const Spec* s = find(name);
    if (!s) {
        std::fprintf(stderr, "internal error: %s is read through env::%s but is not in the "
                             "flag table (src/envflags.cpp)\n", name, reader);
        std::abort();
    }
    for (Kind k : kinds)
        if (s->kind == k) return *s;
    std::fprintf(stderr, "internal error: %s is a %s flag but is read through env::%s\n",
                 name, kindName(s->kind), reader);
    std::abort();
}

std::string errorLine(const Spec& s, const char* v, const std::string& why) {
    return "error: " + std::string(s.name) + "='" + v + "' is not valid: " + why;
}

[[noreturn]] void reject(const Spec& s, const char* v, const std::string& why) {
    std::fflush(stdout);
    std::fprintf(stderr, "%s\n", errorLine(s, v, why).c_str());
    std::exit(2);
}

Parsed parsedOrExit(const Spec& s, const char* v) {
    Parsed p;
    std::string why;
    if (!parse(s, v, p, why)) reject(s, v, why);
    return p;
}

}  // namespace

const char* kindName(Kind k) {
    switch (k) {
        case Kind::Integer:  return "integer";
        case Kind::Real:     return "real";
        case Kind::Switch:   return "switch";
        case Kind::Binary:   return "binary";
        case Kind::Presence: return "presence";
        case Kind::One:      return "one";
        case Kind::Text:     return "text";
        case Kind::Choice:   return "choice";
        case Kind::External: return "external";
    }
    return "?";
}

const Spec* table(std::size_t& count) {
    count = sizeof(kTable) / sizeof(kTable[0]);
    return kTable;
}

const Spec* find(const char* name) {
    if (!name) return nullptr;
    for (const Spec& s : kTable)
        if (std::strcmp(s.name, name) == 0) return &s;
    return nullptr;
}

bool parse(const Spec& s, const char* v, Parsed& out, std::string& why) {
    out = Parsed();
    if (!v) { why = "no value"; return false; }
    switch (s.kind) {
        case Kind::Integer:
            if (!parseInteger(v, s.ilo, s.ihi, out.i, why)) {
                why += " (expected a base-10 integer in [" + std::to_string(s.ilo) + ", " +
                       std::to_string(s.ihi) + "])";
                return false;
            }
            out.canonical = std::to_string(out.i);
            return true;
        case Kind::Real: {
            const std::string expect = " (expected a finite number in [" + shortest(s.rlo) + ", " +
                                       shortest(s.rhi) + "])";
            if (!*v) { why = "empty value" + expect; return false; }
            if (startsWithSpace(v)) { why = "leading whitespace" + expect; return false; }
            char* end = nullptr;
            const double x = std::strtod(v, &end);   // atof(v) is defined as this call
            if (end == v || *end != '\0') { why = "not a number" + expect; return false; }
            if (!std::isfinite(x)) { why = "not finite" + expect; return false; }
            if (x < s.rlo || x > s.rhi) { why = "out of range" + expect; return false; }
            out.r = x;
            out.canonical = shortest(x);
            return true;
        }
        case Kind::Switch: {
            long long x = 0;
            if (!parseInteger(v, kIntMin, kIntMax, x, why)) {
                why += " (expected an integer: 0 = off, any other value = on)";
                return false;
            }
            out.i = x;
            out.on = x != 0;
            out.canonical = std::to_string(x) + (out.on ? " (on)" : " (off)");
            return true;
        }
        case Kind::Binary:
            if (std::strcmp(v, "0") == 0) { out.on = false; out.canonical = "0 (off)"; return true; }
            if (std::strcmp(v, "1") == 0) { out.on = true; out.canonical = "1 (on)"; return true; }
            why = "expected exactly 0 or 1";
            return false;
        case Kind::One:
            if (std::strcmp(v, "1") == 0) { out.on = true; out.canonical = "1 (on)"; return true; }
            why = "expected exactly 1 (any setting turns this on; unset it to turn it off)";
            return false;
        case Kind::Presence:
            out.on = true;
            out.canonical = std::string(v) + " (set)";
            return true;
        case Kind::Text:
        case Kind::External:
            out.canonical = v;
            return true;
        case Kind::Choice:
            if (!isChoice(s.choices, v, s.ilo, s.ihi)) {
                why = std::string("expected one of ") + s.choices;
                if (s.ilo <= s.ihi && std::strchr(s.choices, '#'))
                    why += " (# = a base-10 integer in [" + std::to_string(s.ilo) + ", " +
                           std::to_string(s.ihi) + "])";
                return false;
            }
            out.canonical = v;
            return true;
    }
    why = "unknown kind";
    return false;
}

long long integer(const char* name, long long def) {
    const Spec& s = lookup(name, {Kind::Integer}, "integer");
    const char* v = std::getenv(name);
    return v ? parsedOrExit(s, v).i : def;
}

double real(const char* name, double def) {
    const Spec& s = lookup(name, {Kind::Real}, "real");
    const char* v = std::getenv(name);
    return v ? parsedOrExit(s, v).r : def;
}

bool on(const char* name, bool def) {
    const Spec& s = lookup(name, {Kind::Switch, Kind::Binary}, "on");
    const char* v = std::getenv(name);
    return v ? parsedOrExit(s, v).on : def;
}

bool present(const char* name) {
    const Spec& s = lookup(name, {Kind::Presence, Kind::One}, "present");
    const char* v = std::getenv(name);
    return v ? parsedOrExit(s, v).on : false;
}

const char* text(const char* name) {
    const Spec& s = lookup(name, {Kind::Text, Kind::Choice}, "text");
    const char* v = std::getenv(name);
    if (v) parsedOrExit(s, v);
    return v;
}

bool isSet(const char* name) {
    if (!find(name)) {
        std::fprintf(stderr, "internal error: %s is read through env::isSet but is not in the "
                             "flag table (src/envflags.cpp)\n", name);
        std::abort();
    }
    return std::getenv(name) != nullptr;
}

int validateEnvironment(std::FILE* log) {
    std::vector<std::pair<std::string, std::string>> lines, errors;
    for (char** e = ::environ; e && *e; ++e) {
        if (std::strncmp(*e, "TESSERACT_", 10) != 0) continue;
        const char* eq = std::strchr(*e, '=');
        const std::string name = eq ? std::string(*e, static_cast<size_t>(eq - *e)) : std::string(*e);
        const char* value = eq ? eq + 1 : "";
        const Spec* s = find(name.c_str());
        if (!s) {
            lines.push_back({name, "[config] warning: " + name +
                                   " is set but is not a TesserACT flag; it has no effect"});
            continue;
        }
        if (s->kind == Kind::External) continue;
        Parsed p;
        std::string why;
        if (!parse(*s, value, p, why)) {
            errors.push_back({name, errorLine(*s, value, why)});
            continue;
        }
        lines.push_back({name, "[config] " + name + "=" + p.canonical});
    }
    if (!errors.empty()) {
        std::sort(errors.begin(), errors.end());
        for (const auto& l : errors) std::fprintf(log, "%s\n", l.second.c_str());
        std::fflush(log);
        return 2;
    }
    std::sort(lines.begin(), lines.end());
    for (const auto& l : lines) std::fprintf(log, "%s\n", l.second.c_str());
    std::fflush(log);
    return 0;
}

std::vector<SetVariable> setVariables() {
    std::vector<SetVariable> out;
    for (char** e = ::environ; e && *e; ++e) {
        if (std::strncmp(*e, "TESSERACT_", 10) != 0) continue;
        const char* eq = std::strchr(*e, '=');
        SetVariable v;
        v.name = eq ? std::string(*e, static_cast<size_t>(eq - *e)) : std::string(*e);
        v.value = eq ? std::string(eq + 1) : std::string();
        if (const Spec* s = find(v.name.c_str())) {
            v.registered = true;
            v.kind = kindName(s->kind);
            Parsed p;
            std::string why;
            v.valid = parse(*s, v.value.c_str(), p, why);
            v.canonical = v.valid ? p.canonical : why;
        }
        out.push_back(v);
    }
    std::sort(out.begin(), out.end(), [](const SetVariable& a, const SetVariable& b) { return a.name < b.name; });
    return out;
}

}}  // namespace ts::env
