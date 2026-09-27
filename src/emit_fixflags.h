// Switches for the G-emit defect fixes (combo3, 2026-09-25). Default OFF in build_v3; since
// 1.4.0 the umbrella TESSERACT_FIXES is ON by default (src/defaults.h), so every fix below is
// on unless TESSERACT_FIXES=0 or its own TESSERACT_FIX_<NAME>=0.
//
// Each output-changing fix is behind its own TESSERACT_FIX_<NAME> and all of them are
// switched on together by the umbrella TESSERACT_FIXES=1. An explicit per-fix value wins
// over the umbrella, so TESSERACT_FIXES=1 TESSERACT_FIX_TRIM_COPY_GUARD=0 is every fix but
// one. With nothing set the release code path runs.
//
// Values are read on every call -- never cached in a function-local static -- so a batch or
// library caller that changes the environment between runs gets the value it set (T16).
// Values are strict: exactly "0" or "1" (GAPFILL_BACKOFF also takes a base count). build_v3:
// every name is a row of the envflags table (T16), whose readers exit with status 2 on a
// malformed value; main() validates the whole environment before any work, and validate()
// below re-reads every switch at the start of Assembler::run for library callers.
//
// build_v3 package aliases (the campaign package build_v2_merge carried its own
// implementations of two of these fixes; build_v3 keeps ONE implementation of each):
//   TESSERACT_BOUNDARY_SAFE_TRIM=0|1 (Q40)  -> T10, TESSERACT_FIX_BOUNDARY_SAFE_TRIM
//   TESSERACT_GAPFILL_BACKOFF=<B>, 0..200    -> T11's back-off with B bases (package
//                                               meaning: 1 is one base, 0 is off)
// Precedence: the fix's own variable, then the package alias, then the umbrella.
//
// Every name below appears as one complete string literal, because the wrapper guards
// check `strings -a` of the binary for the exact variable name.
#pragma once

#include <cstdio>
#include <cstring>
#include <string>

#include "defaults.h"
#include "envflags.h"

namespace ts {
namespace emitfix {

constexpr const char* kUmbrella          = "TESSERACT_FIXES";
constexpr const char* kPolishSkipN       = "TESSERACT_FIX_POLISH_SKIP_N";         // T02
constexpr const char* kGapfillStrict     = "TESSERACT_FIX_GAPFILL_STRICT_BUDGET"; // T09
constexpr const char* kGapfillBackoff    = "TESSERACT_FIX_GAPFILL_BACKOFF";       // T11
constexpr const char* kGapfillSkipInputN = "TESSERACT_FIX_GAPFILL_SKIP_INPUT_N";  // T27
constexpr const char* kBoundarySafeTrim  = "TESSERACT_FIX_BOUNDARY_SAFE_TRIM";    // T10
constexpr const char* kTrimCopyGuard     = "TESSERACT_FIX_TRIM_COPY_GUARD";       // T22
constexpr const char* kSplitPostprocess  = "TESSERACT_FIX_SPLIT_POSTPROCESS";     // T25
constexpr const char* kAgpGfaV2          = "TESSERACT_FIX_AGP_GFA_V2";            // T26
// Package aliases (build_v2_merge PKG-MERGE Q40 and PKG-GAP).
constexpr const char* kPkgBoundarySafeTrim = "TESSERACT_BOUNDARY_SAFE_TRIM";      // -> T10
constexpr const char* kPkgGapfillBackoff   = "TESSERACT_GAPFILL_BACKOFF";         // -> T11

constexpr const char* kBinarySwitches[] = {kPolishSkipN, kGapfillStrict, kGapfillSkipInputN,
                                           kBoundarySafeTrim, kTrimCopyGuard, kSplitPostprocess,
                                           kAgpGfaV2};

// Back-off budget used when TESSERACT_FIX_GAPFILL_BACKOFF is "1" or the umbrella turns it on.
constexpr int kDefaultBackoff = 60;
constexpr int kMaxBackoff = 200;

inline bool umbrella() { return defaults::fixesUmbrella(); }   // kUmbrella, default ON (1.4.0)

// The package alias of a switch, or nullptr.
inline const char* aliasOf(const char* name) {
    return std::strcmp(name, kBoundarySafeTrim) == 0 ? kPkgBoundarySafeTrim : nullptr;
}

// Effective state of one binary switch: an explicit 0/1 wins, then its package alias, and
// unset follows the umbrella. A malformed value exits with status 2 (envflags).
inline bool enabled(const char* name) {
    if (env::isSet(name)) return env::on(name, false);
    if (const char* alias = aliasOf(name))
        if (env::isSet(alias)) return env::on(alias, false);
    return umbrella();
}

// TESSERACT_FIX_GAPFILL_BACKOFF: "0" -> 0; "1" -> kDefaultBackoff; 2..kMaxBackoff -> that many
// bases. Unset: the package alias TESSERACT_GAPFILL_BACKOFF=<B> (B bases, 0 = off), else
// kDefaultBackoff under the umbrella, else 0. Out-of-range values exit 2 (envflags table).
inline int backoffBases() {
    if (env::isSet(kGapfillBackoff)) {
        const long long b = env::integer(kGapfillBackoff, 0);
        return b == 1 ? kDefaultBackoff : static_cast<int>(b);
    }
    if (env::isSet(kPkgGapfillBackoff)) return static_cast<int>(env::integer(kPkgGapfillBackoff, 0));
    return umbrella() ? kDefaultBackoff : 0;
}

// Reads every switch once, so a malformed value stops the run before any stage reads it
// (env readers exit 2 naming the variable). Returns true; `error` is kept for the API.
inline bool validate(std::string& error) {
    (void)error;
    (void)umbrella();
    for (const char* name : kBinarySwitches) (void)enabled(name);
    (void)backoffBases();
    return true;
}

// The effective state of every G-emit switch, one line, printed once per run.
inline void logSummary() {
    std::fprintf(stderr,
                 "[fixes-emit] umbrella=%d polish_skip_n=%d gapfill_strict_budget=%d gapfill_backoff=%d "
                 "gapfill_skip_input_n=%d boundary_safe_trim=%d trim_copy_guard=%d "
                 "split_postprocess=%d agp_gfa_v2=%d\n",
                 umbrella() ? 1 : 0, enabled(kPolishSkipN) ? 1 : 0, enabled(kGapfillStrict) ? 1 : 0,
                 backoffBases(), enabled(kGapfillSkipInputN) ? 1 : 0,
                 enabled(kBoundarySafeTrim) ? 1 : 0, enabled(kTrimCopyGuard) ? 1 : 0,
                 enabled(kSplitPostprocess) ? 1 : 0, enabled(kAgpGfaV2) ? 1 : 0);
}

}  // namespace emitfix
}  // namespace ts
