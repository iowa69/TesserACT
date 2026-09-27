// T16 (build_v3, group G-xcut), phase 3: after the call-site retrofit every TESSERACT_* flag
// is read through the strict readers, per call, and never cached.
//
// 1. Strict source scan (tests/env_scan.h): no raw getenv outside envflags.cpp, no static
//    caches a flag value, every registered flag is read.
// 2. End to end through the real PairedResolver on the T16 fallback fixture (chain end
//    A -> repeat R -> {C, D}, no linking pairs, so TESSERACT_MIN_FALLBACK_DEST decides whether
//    A and C are joined):
//      - valid values give the release results;
//      - a value changed between two resolves in ONE process takes effect (release: the first
//        value was kept for the whole process);
//      - a malformed value makes the resolver itself exit 2 (library use, no main()).
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <random>
#include <string>
#include <unistd.h>
#include <vector>

#include "envflags.h"
#include "resolve.h"
#include "env_reject.h"
#include "env_scan.h"

namespace {

int checks = 0;
void check(bool ok, const std::string& label) {
    ++checks;
    if (!ok) { std::fprintf(stderr, "FAIL: %s\n", label.c_str()); std::exit(1); }
}

void clearTesseractEnvironment() {
    std::vector<std::string> names;
    for (char** e = ::environ; e && *e; ++e) {
        if (std::strncmp(*e, "TESSERACT_", 10) != 0) continue;
        const char* eq = std::strchr(*e, '=');
        names.emplace_back(*e, eq ? static_cast<size_t>(eq - *e) : std::strlen(*e));
    }
    for (const std::string& n : names) ::unsetenv(n.c_str());
}

// ---- 2. no caching, end to end through the real resolver (T16 fixture) ----------------------

std::string dna(size_t n, unsigned seed) {
    std::mt19937 r(seed);
    std::string s(n, 'A');
    for (char& c : s) c = "ACGT"[r() % 4];
    return s;
}
void edge(ts::UnitigGraph& g, unsigned a, unsigned b) {
    g.nodes[a].ends[1].push_back({b, 0});
    g.nodes[b].ends[0].push_back({a, 1});
}
struct Result { bool joined = false; long long shortDest = -1; int traceLines = 0; };

Result resolveFixture() {
    ts::UnitigGraph g;
    g.setK(31);
    const size_t ov = 30;
    g.nodes.resize(9);
    const std::string rep = dna(200, 2);
    g.nodes[1].seq = rep;                                                  // R repeat, 50x
    g.nodes[0].seq = dna(300, 1).substr(0, 270) + rep.substr(0, ov);        // A source, 10x
    g.nodes[5].seq = dna(300, 6).substr(0, 270) + rep.substr(0, ov);        // B source, 3x
    g.nodes[2].seq = rep.substr(rep.size() - ov) + dna(270, 3);             // C destination, 10x
    g.nodes[3].seq = rep.substr(rep.size() - ov) + dna(270, 4);             // D destination, 3x
    g.nodes[4].seq = dna(300, 5);                                           // filler
    for (size_t i = 6; i < 9; ++i) g.nodes[i].seq = dna(300, 10 + unsigned(i));
    for (auto& n : g.nodes) n.coverage = 10;
    g.nodes[1].coverage = 50; g.nodes[5].coverage = 3; g.nodes[3].coverage = 3;
    edge(g, 0, 1); edge(g, 5, 1); edge(g, 1, 2); edge(g, 1, 3);
    check(g.validate().empty(), "fixture graph valid");

    // Pairs wholly inside the filler: reads exist but none links A to C or D, so the join
    // falls to the coverage fallback, which MIN_FALLBACK_DEST can withdraw.
    const auto dir = std::filesystem::temp_directory_path() / ("test_envflags_" + std::to_string(::getpid()));
    std::filesystem::create_directories(dir);
    const auto r1 = dir / "r1.fa", r2 = dir / "r2.fa";
    {
        std::ofstream a(r1), b(r2);
        for (int i = 0; i < 20; ++i) {
            a << ">p" << i << "/1\n" << g.nodes[4].seq.substr(10 + i, 70) << '\n';
            b << ">p" << i << "/2\n" << ts::reverseComplement(g.nodes[4].seq.substr(200 + i, 70)) << '\n';
        }
    }
    ts::Library lib;
    lib.r1 = r1.string();
    lib.r2 = r2.string();
    ts::SequenceStore reads;
    std::string err;
    check(reads.load({lib}, 1, err), "fixture reads load");

    std::fflush(stderr);
    std::FILE* cap = std::tmpfile();
    const int saved = ::dup(2);
    ::dup2(::fileno(cap), 2);
    ts::PairedResolver res(g, reads, 1, 2, 1.02, 0.02);
    res.setInsertBounds(0, 2000);
    res.buildSupport();
    std::vector<std::string> contigs;
    std::vector<double> covs;
    res.resolve(contigs, covs);
    std::fflush(stderr);
    ::dup2(saved, 2);
    ::close(saved);
    std::rewind(cap);
    std::string log;
    char buf[4096];
    size_t n = 0;
    while ((n = std::fread(buf, 1, sizeof buf, cap)) > 0) log.append(buf, n);
    std::fclose(cap);
    std::filesystem::remove_all(dir);

    Result out;
    const std::string aCore = g.nodes[0].seq.substr(0, 200), cCore = g.nodes[2].seq.substr(100, 150);
    for (const auto& c : contigs) {
        const std::string rc = ts::reverseComplement(c);
        if ((c.find(aCore) != std::string::npos && c.find(cCore) != std::string::npos) ||
            (rc.find(aCore) != std::string::npos && rc.find(cCore) != std::string::npos))
            out.joined = true;
    }
    const size_t p = log.find("short-dest=");
    if (p != std::string::npos) out.shortDest = std::atoll(log.c_str() + p + 11);
    for (size_t q = 0; (q = log.find("[jointrace]", q)) != std::string::npos; ++q) ++out.traceLines;
    return out;
}

void setOrUnset(const char* name, const char* v) {
    if (v) ::setenv(name, v, 1); else ::unsetenv(name);
}

void resolverNoCaching() {
    clearTesseractEnvironment();
    ::setenv("TESSERACT_COMMON_PREFIX", "0", 1);   // isolate chain arbitration
    ::setenv("TESSERACT_DEBUG_RESOLVE", "1", 1);   // prints the short-dest counter

    // Valid values give the release results (repro_output.txt, verify/T16).
    struct V { const char* value; bool joined; long long shortDest; };
    for (const V& v : std::vector<V>{{nullptr, true, 0}, {"0", true, 0}, {"250", true, 0},
                                     {"301", false, 4}, {"1000000000", false, 4}}) {
        setOrUnset("TESSERACT_MIN_FALLBACK_DEST", v.value);
        const Result r = resolveFixture();
        check(r.joined == v.joined && r.shortDest == v.shortDest,
              std::string("MIN_FALLBACK_DEST=") + (v.value ? v.value : "(unset)") + " gives the release result");
    }

    // Changing the value between two resolves in ONE process takes effect (release: the
    // first value was kept -- joined=0 twice, then joined=1 twice).
    setOrUnset("TESSERACT_MIN_FALLBACK_DEST", "1000000000");
    const Result a1 = resolveFixture();
    setOrUnset("TESSERACT_MIN_FALLBACK_DEST", "0");
    const Result a2 = resolveFixture();
    check(!a1.joined && a1.shortDest == 4 && a2.joined && a2.shortDest == 0, "guard on then off in one process");
    setOrUnset("TESSERACT_MIN_FALLBACK_DEST", "1000000000");
    const Result a3 = resolveFixture();
    check(!a3.joined && a3.shortDest == 4, "guard off then on in one process");
    setOrUnset("TESSERACT_MIN_FALLBACK_DEST", nullptr);

    setOrUnset("TESSERACT_JOIN_TRACE", nullptr);
    const Result t1 = resolveFixture();
    setOrUnset("TESSERACT_JOIN_TRACE", "1");
    const Result t2 = resolveFixture();
    check(t1.traceLines == 0 && t2.traceLines == 4, "JOIN_TRACE enabled mid-process traces the second resolve");

    // The resolver itself refuses a malformed value (library use has no main() to validate).
    // resolveFixture redirects stderr itself, so only the exit status is checked here; the
    // message is checked by readerRejections().
    for (const char* bad : {"1e9", "abc", "-1", "4294967296"}) {
        check(exitsWithEnvError(nullptr, [&] {
                  ::setenv("TESSERACT_MIN_FALLBACK_DEST", bad, 1);
                  (void)resolveFixture();
              }), std::string("resolver exits 2 on MIN_FALLBACK_DEST='") + bad + "'");
    }
    clearTesseractEnvironment();
}

}  // namespace

int main() {
    resolverNoCaching();

    const std::filesystem::path src = envscan::sourceDir(__FILE__);
    check(std::filesystem::exists(src / "envflags.cpp"),
          "source tree found at " + src.string() + " (run from the tree root)");
    const std::vector<std::string> v = envscan::scan(src, true);
    for (const std::string& s : v) std::fprintf(stderr, "  %s\n", s.c_str());
    check(v.empty(), std::to_string(v.size()) + " strict source-scan violation(s) in " + src.string());

    std::printf("env flags retrofit: %d checks passed\n", checks);
    return 0;
}
