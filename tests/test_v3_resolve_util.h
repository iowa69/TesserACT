// Shared helpers for the build_v3 G-resolve component tests (tests/test_v3_resolve_*.cpp).
// Header only: the Makefile builds tests/test_*.cpp, never this file on its own.
#pragma once

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <random>
#include <string>
#include <utility>
#include <vector>
#include <unistd.h>

#include "graph.h"
#include "resolve.h"
#include "seqio.h"

namespace v3r {

inline int& failures() { static int f = 0; return f; }
inline int& checks() { static int c = 0; return c; }
inline void expect(bool ok, const std::string& what) {
    ++checks();
    std::printf("  %s  %s\n", ok ? "PASS" : "FAIL", what.c_str());
    if (!ok) ++failures();
}
inline int finish(const char* name) {
    std::printf("%s: %d/%d checks failed\n", name, failures(), checks());
    return failures() ? 1 : 0;
}

inline std::string dna(size_t n, unsigned seed) {
    std::mt19937 rng(seed);
    std::string s(n, 'A');
    for (char& b : s) b = "ACGT"[rng() % 4];
    return s;
}
inline std::string rc(const std::string& s) { return ts::reverseComplement(s); }
using Pair = std::pair<std::string, std::string>;

inline ts::SequenceStore store(const std::vector<Pair>& pairs, const std::string& tag) {
    static int serial = 0;
    const auto dir = std::filesystem::temp_directory_path() /
                     ("tesseract-v3r-" + tag + "-" + std::to_string(getpid()) + "-" + std::to_string(serial++));
    std::filesystem::create_directories(dir);
    const auto first = dir / "r1.fa", second = dir / "r2.fa";
    {
        std::ofstream f(first), s(second);
        for (size_t i = 0; i < pairs.size(); ++i) {
            f << ">" << i << "\n" << pairs[i].first << "\n";
            s << ">" << i << "\n" << pairs[i].second << "\n";
        }
    }
    ts::SequenceStore st;
    ts::Library lib;
    lib.r1 = first.string();
    lib.r2 = second.string();
    std::string error;
    if (!st.load({lib}, 1, error)) { std::fprintf(stderr, "load: %s\n", error.c_str()); std::exit(2); }
    std::filesystem::remove_all(dir);
    return st;
}

// Redirects file descriptor 2 into a temporary file between start() and stop().
class StderrCapture {
public:
    void start() {
        std::fflush(stderr);
        file_ = std::tmpfile();
        saved_ = dup(2);
        dup2(fileno(file_), 2);
    }
    std::string stop() {
        std::fflush(stderr);
        dup2(saved_, 2);
        close(saved_);
        std::rewind(file_);
        std::string out;
        char buf[4096];
        size_t n = 0;
        while ((n = std::fread(buf, 1, sizeof(buf), file_)) > 0) out.append(buf, n);
        std::fclose(file_);
        return out;
    }
private:
    FILE* file_ = nullptr;
    int saved_ = -1;
};

// All lines of `log` that start with `tag` (for example "[gapflank]").
inline std::vector<std::string> lines(const std::string& log, const std::string& tag) {
    std::vector<std::string> out;
    size_t p = 0;
    while (p < log.size()) {
        size_t e = log.find('\n', p);
        if (e == std::string::npos) e = log.size();
        const std::string l = log.substr(p, e - p);
        if (l.compare(0, tag.size(), tag) == 0) out.push_back(l);
        p = e + 1;
    }
    return out;
}
// Integer value of `key=` in `line`; `def` when absent.
inline long long field(const std::string& line, const std::string& key, long long def = -1) {
    const std::string pat = " " + key + "=";
    const size_t p = line.find(pat);
    if (p == std::string::npos) return def;
    return std::atoll(line.c_str() + p + pat.size());
}

// Unset every flag that changes resolver decisions, so each test starts from the
// release-default resolver. Function-local statics in resolve.cpp (JOIN_TRACE,
// REQUIRE_SUPPORT_SINGLE, MIN_FALLBACK_DEST, COMMON_PREFIX, ...) cache the value seen by
// the FIRST resolve() call of the process, so tests set them once, here, before any run.
inline void releaseDefaults() {
    setenv("TESSERACT_COMMON_PREFIX", "0", 1);
    // V3R_KEEP_TRACE=1 keeps TESSERACT_JOIN_TRACE for debugging a fixture (stderr only).
    const bool keepTrace = std::getenv("V3R_KEEP_TRACE") != nullptr;
    const char* trace = std::getenv("TESSERACT_JOIN_TRACE");
    const std::string traceValue = trace ? trace : "";
    for (const char* f : {"TESSERACT_ROUTE_DISTANCE", "TESSERACT_EXACT_READ_THREADS", "TESSERACT_PREFIX_SNP_BUBBLES",
                          "TESSERACT_EXCLUDE_SHARED_REPEAT_SUPPORT", "TESSERACT_REQUIRE_SUPPORT_SINGLE",
                          "TESSERACT_SCAF_SUPPORT", "TESSERACT_JOIN_TRACE", "TESSERACT_MIN_FALLBACK_DEST",
                          "TESSERACT_WEIGHTED_RESOLVER_COVERAGE", "TESSERACT_WEIGHTED_ELIGIBLE_COVERAGE",
                          "TESSERACT_NO_UNSPANNED_FALLBACK", "TESSERACT_SHARED_SUPPORT_AUDIT",
                          "TESSERACT_FIXES", "TESSERACT_FIX_GAP_FLANK", "TESSERACT_FIX_REVISIT_GUARD",
                          "TESSERACT_FIX_SCAFFOLD_CYCLE",
                          "TESSERACT_FIX_GAP_ESTIMATE", "TESSERACT_FIX_TRUNC_GUARD", "TESSERACT_FIX_COV_CONTRIB",
                          "TESSERACT_FIX_ROUTE_ORDER", "TESSERACT_FIX_MIRROR_ROUTE"})
        unsetenv(f);
    if (keepTrace && trace) setenv("TESSERACT_JOIN_TRACE", traceValue.c_str(), 1);
}
inline void setFlag(const char* name, const char* value) {
    if (value) setenv(name, value, 1); else unsetenv(name);
}

}  // namespace v3r
