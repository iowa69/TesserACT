// Test-process environment hygiene.
//
// The library reads its experiment switches from TESSERACT_* environment variables, and
// several are cached in function-local statics on first use. A test binary started from a
// shell that exports an experiment environment (for example the K2 arm: COMMON_PREFIX=0,
// MIN_FALLBACK_DEST=1000000000, ...) therefore ran its fixtures under that configuration:
// test_resolver_prefix, test_resolver_prefix_ownership and test_resolver_read_threads
// aborted, and any test could equally pass because an ambient flag hid the behaviour it
// exists to check.
//
// clearTesseractEnv() removes every TESSERACT_* variable. Call it FIRST in main(), before
// anything reaches the library, then set the flags the test depends on explicitly.
#pragma once

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

// C linkage, matching glibc <unistd.h> (which declares it under _GNU_SOURCE) and platforms
// whose headers do not declare it at all.
extern "C" char** environ;

namespace testenv {

inline size_t clearTesseractEnv() {
    // Collect first: unsetenv() rewrites environ while it is being walked.
    std::vector<std::string> names;
    for (char** e = environ; e != nullptr && *e != nullptr; ++e) {
        const char* eq = std::strchr(*e, '=');
        if (eq == nullptr) continue;
        std::string name(*e, static_cast<size_t>(eq - *e));
        if (name.compare(0, 10, "TESSERACT_") == 0) names.push_back(name);
    }
    for (const std::string& name : names) unsetenv(name.c_str());
    if (!names.empty()) {
        std::fprintf(stderr, "  testenv: cleared %zu ambient TESSERACT_* variable(s)\n",
                     names.size());
    }
    return names.size();
}

// 1.4.0 changed eight defaults (src/defaults.h), six of them environment flags. A test written
// against the release 1.3.0 resolver pins their 1.3.0 values with this, after clearing the
// environment; the common-prefix walk budget and the tie ratio stay each test's own choice.
inline void pinRelease130Flags() {
    setenv("TESSERACT_FIXES", "0", 1);
    setenv("TESSERACT_PREFIX_MIN_BODY", "0", 1);
    setenv("TESSERACT_MIN_FALLBACK_DEST", "0", 1);
    setenv("TESSERACT_REQUIRE_SUPPORT_SINGLE", "0", 1);
    setenv("TESSERACT_EXCLUDE_SHARED_REPEAT_SUPPORT", "0", 1);
    setenv("TESSERACT_DROPOUT_BRIDGE", "0", 1);
}

}  // namespace testenv
