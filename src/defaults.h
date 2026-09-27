// TesserACT 1.4.0 defaults: the combo3 "F2" configuration (arm cb3s_eb_t3), misassembly-first.
//
// One named constant per default, and one reader per setting. Every read site takes its
// fallback from this header (RELEASE_PLAN 1.4.0 step 3), and main() prints the effective
// values on one "[defaults]" stderr line on every run, so every log records what ran.
//
// What changed from 1.3.0, and the opt-out that restores the 1.3.0 value:
//
//   setting                                   1.3.0   1.4.0        opt-out
//   --tie-ratio                               1.02    3.0          --tie-ratio 1.02
//   TESSERACT_COMMON_PREFIX (walk budget, bp) 3000    250          TESSERACT_COMMON_PREFIX=3000
//   TESSERACT_PREFIX_MIN_BODY (ENDS-A)        off     reach        TESSERACT_PREFIX_MIN_BODY=0
//   TESSERACT_MIN_FALLBACK_DEST               0       1000000000   TESSERACT_MIN_FALLBACK_DEST=0
//   TESSERACT_REQUIRE_SUPPORT_SINGLE          off     on           TESSERACT_REQUIRE_SUPPORT_SINGLE=0
//   TESSERACT_EXCLUDE_SHARED_REPEAT_SUPPORT   off     on           TESSERACT_EXCLUDE_SHARED_REPEAT_SUPPORT=0
//   TESSERACT_FIXES (umbrella of the 21       off     on           TESSERACT_FIXES=0, or
//     output-changing TESSERACT_FIX_* fixes)                       TESSERACT_FIX_<NAME>=0 for one fix
//   TESSERACT_DROPOUT_BRIDGE                  off     on           TESSERACT_DROPOUT_BRIDGE=0
//
// TESSERACT_JOIN_TRACE (logging only) stays off. The run-mode presets are unchanged
// (--mode fast 1.3, careful 1.4, aggressive 1.05 tie ratio), so `careful` is no longer
// stricter than `standard` on the tie ratio.
//
// Each reader parses through the envflags table (strict, per call, never cached): a set
// value always wins, and an unset flag takes the constant below.
#pragma once

#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <string>

#include "envflags.h"

namespace ts {
namespace defaults {

constexpr double kDefaultTieRatio = 3.0;
constexpr double kDefaultCommonPrefix = 250.0;
constexpr const char* kDefaultPrefixMinBody = "reach";
constexpr long long kDefaultMinFallbackDest = 1000000000;
constexpr bool kDefaultRequireSupportSingle = true;
constexpr bool kDefaultExcludeSharedRepeat = true;
constexpr bool kDefaultFixes = true;
constexpr bool kDefaultDropoutBridge = true;

// Common-prefix walk budget in bp (<= 0 disables the walk).
inline double commonPrefix() { return env::real("TESSERACT_COMMON_PREFIX", kDefaultCommonPrefix); }

// ENDS-A minimum body, as written: "reach" or a base count ("0" = off).
inline std::string prefixMinBodyText() {
    const char* e = env::text("TESSERACT_PREFIX_MIN_BODY");
    return e ? std::string(e) : std::string(kDefaultPrefixMinBody);
}

// ENDS-A minimum body as the resolver uses it: 0 = off, -1 = the insert model's reach,
// otherwise that many bases.
inline long prefixMinBody() {
    const std::string v = prefixMinBodyText();
    if (v == "reach") return -1L;
    return std::strtol(v.c_str(), nullptr, 10);   // validated: an integer in [0, INT_MAX]
}

// The three fallback joins are withdrawn below this destination length (0 = fallbacks on).
inline std::size_t minFallbackDest() {
    return static_cast<std::size_t>(env::integer("TESSERACT_MIN_FALLBACK_DEST", kDefaultMinFallbackDest));
}

inline bool requireSupportSingle() {
    return env::on("TESSERACT_REQUIRE_SUPPORT_SINGLE", kDefaultRequireSupportSingle);
}

inline bool excludeSharedRepeatSupport() {
    return env::on("TESSERACT_EXCLUDE_SHARED_REPEAT_SUPPORT", kDefaultExcludeSharedRepeat);
}

// The defect-fix umbrella. A fix's own TESSERACT_FIX_<NAME> (0 or 1) always wins over it.
inline bool fixesUmbrella() { return env::on("TESSERACT_FIXES", kDefaultFixes); }

inline bool dropoutBridge() { return env::on("TESSERACT_DROPOUT_BRIDGE", kDefaultDropoutBridge); }

// The "[defaults]" line: the eight effective values of this run, whether they come from the
// defaults above, the environment or the command line. `tieRatio` is the effective one (after
// --tie-ratio and --mode).
inline std::string line(double tieRatio) {
    char buf[512];
    std::snprintf(buf, sizeof buf,
                  "[defaults] tie_ratio=%.3f common_prefix=%.0f prefix_min_body=%s min_fallback_dest=%zu "
                  "require_support_single=%d exclude_shared_repeat_support=%d fixes=%d dropout_bridge=%d",
                  tieRatio, commonPrefix(), prefixMinBodyText().c_str(), minFallbackDest(),
                  requireSupportSingle() ? 1 : 0, excludeSharedRepeatSupport() ? 1 : 0,
                  fixesUmbrella() ? 1 : 0, dropoutBridge() ? 1 : 0);
    return buf;
}

}  // namespace defaults
}  // namespace ts
