// Environment flags for the build_v3 G-graph defect fixes (T01, T08, T28, T29, T30).
//
// Rules these helpers enforce, from combo3 triage T16 and OBJECTIVE amendment A2:
//   * read on every call, never cached in a function-local static, so a library or
//     multi-configuration driver sees the value in force when the stage runs;
//   * parsed strictly: an unset variable means "not set"; anything else, the empty string
//     included (build_v3: the envflags table contract), must be exactly 0 or 1 (or an
//     integer in [0, 2] for the carry gate), or the process exits with status 2 and names
//     the variable -- a typo must never silently mean "off";
//   * every output-changing fix is default OFF behind TESSERACT_FIX_<NAME>, and the
//     umbrella TESSERACT_FIXES=1 turns on every fix whose own variable is unset.
//     TESSERACT_FIX_<NAME>=0 keeps that one fix off even under the umbrella.
//
// Callers pass the full literal variable name so that the name is present in the binary
// and the wrappers' `strings -a BIN | grep -qx NAME` guard can find it.
//
// build_v3 (integration): folded onto src/envflags.h (T16, G-xcut). Every name is a row of
// the flag table; the readers there are strict and per call (a malformed value exits with
// status 2 naming the variable, the release T16 contract).
#pragma once

#include "envflags.h"

namespace ts {
namespace fixflags {

// The umbrella: TESSERACT_FIXES=1 turns on every fix whose own variable is unset.
inline bool umbrella() { return env::on("TESSERACT_FIXES", false); }

// Level of an output-changing defect fix: the fix's own variable when set (a Binary flag
// when maxLevel is 1, an Integer flag in [0, maxLevel] otherwise), else 1 under the
// umbrella, else 0 (release behaviour).
inline int fixLevel(const char* name, int maxLevel = 1) {
    if (env::isSet(name)) {
        if (maxLevel <= 1) return env::on(name, false) ? 1 : 0;
        return static_cast<int>(env::integer(name, 0));
    }
    return umbrella() ? 1 : 0;
}

}  // namespace fixflags
}  // namespace ts
