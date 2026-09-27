// T16 (build_v3, group G-xcut), phase 0: the strict flag registry in src/envflags.{h,cpp} and
// the startup validation in main(). Holds before and after the phase-3 call-site retrofit
// (which tests/test_envflags_retrofit.cpp covers).
//
//   build/test_envflags                     all checks (run from the tree root, as make does)
//   build/test_envflags --scan DIR [strict] only the source scan, on another tree's src/; run
//                                           it on the release sources to see what it catches
//
// 1. Parse contract: for every registered flag and a corpus of candidate strings, a value
//    the table accepts parses to exactly what the release expression (atoi, atof, atol,
//    strcmp) gave; values that were garbage, wrapped or undefined are rejected.
// 2. Readers: a rejected value exits with status 2 and names the flag (the T16 cases); a
//    changed value is seen by the next read (no caching in the readers).
// 3. validateEnvironment: [config] lines for set flags, warning for unknown names, 2 on error.
// 4. Source scan (tests/env_scan.h), non-strict: every TESSERACT_* literal is registered,
//    readers match kinds, no static caches an env:: read, fork-batch flags are 0/1 flags.
#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <map>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <vector>

#include "envflags.h"
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

// ---- 1. parse contract --------------------------------------------------------------------

const std::vector<const char*> kCorpus = {
    "0", "1", "2", "-1", "7", "+5", "007", "10", "250", "301", "3000", "1000000000",
    "2147483647", "2147483648", "-2147483648", "-2147483649", "3000000000", "4294967295",
    "4294967296", "99999999999999999999", "1e9", "1.5", "0.1", "0.35", "-0.5", "1e15", "1e16",
    "1e18", "1e19", "1e-400", "1e400", "abc", "true", "false", "yes", "on", " 7", "7 ", "7x",
    "", "0x10", "inf", "-inf", "nan", "greedy", "mutual", "Greedy", "/tmp/some path"};

void parseContract() {
    size_t n = 0;
    const ts::env::Spec* table = ts::env::table(n);
    check(n > 60, "flag table populated");
    std::set<std::string> names;
    for (size_t i = 0; i < n; ++i) {
        const ts::env::Spec& s = table[i];
        check(std::strncmp(s.name, "TESSERACT_", 10) == 0, std::string("flag name prefix ") + s.name);
        check(names.insert(s.name).second, std::string("flag listed once ") + s.name);
        check(ts::env::find(s.name) == &s, std::string("find() returns the row for ") + s.name);
        for (const char* v : kCorpus) {
            ts::env::Parsed p;
            std::string why;
            const bool ok = ts::env::parse(s, v, p, why);
            const std::string at = std::string(s.name) + "='" + v + "'";
            if (!ok) { check(!why.empty(), "rejection has a reason " + at); continue; }
            switch (s.kind) {
                case ts::env::Kind::Integer: {
                    const bool isLong = std::strstr(s.name, "COMEMBER_") != nullptr;
                    if (isLong) check(std::atol(v) == p.i, "atol agrees " + at);
                    else {
                        check(p.i >= INT_MIN && p.i <= INT_MAX, "atoi-domain flag stays in int " + at);
                        check(std::atoi(v) == p.i, "atoi agrees " + at);
                    }
                    if (s.ilo >= 0) {
                        check(p.i >= 0, "unsigned target never wraps " + at);
                        check(static_cast<size_t>(std::atoi(v)) == static_cast<size_t>(p.i),
                              "cast to size_t agrees " + at);
                    }
                    break;
                }
                case ts::env::Kind::Real:
                    check(std::isfinite(p.r), "real is finite " + at);
                    check(std::atof(v) == p.r, "atof agrees " + at);
                    break;
                case ts::env::Kind::Switch:
                    check((std::atoi(v) != 0) == p.on, "atoi != 0 agrees " + at);
                    check(!(std::atoi(v) != 0) == !p.on, "default-on switch (!e || atoi != 0) agrees " + at);
                    break;
                case ts::env::Kind::Binary:
                    // Every release spelling of a 0/1 flag agrees on the two accepted values.
                    check((std::strcmp(v, "1") == 0) == p.on, "== \"1\" agrees " + at);
                    check((std::strcmp(v, "0") != 0) == p.on, "!= \"0\" agrees " + at);
                    check((*v != '0') == p.on, "IS_VETO_CHR form agrees " + at);
                    check((*v && *v != '0') == p.on, "PLASMID_COMPLETE form agrees " + at);
                    break;
                case ts::env::Kind::One:
                    check(std::strcmp(v, "1") == 0 && p.on, "one-flag accepts only 1 " + at);
                    break;
                case ts::env::Kind::Presence:
                    check(p.on, "presence flag is on when set " + at);
                    break;
                case ts::env::Kind::Choice: {
                    // A listed word, or (build_v3: a "#" in the list) an integer in [ilo, ihi]
                    // that atoi reads the same way.
                    bool word = false;
                    std::string list = std::string(s.choices) + "|";
                    for (size_t a = 0, b; (b = list.find('|', a)) != std::string::npos; a = b + 1)
                        if (list.compare(a, b - a, v) == 0 && std::string(v) != "#") word = true;
                    const bool numeric = std::strchr(s.choices, '#') != nullptr && !word;
                    if (numeric) {
                        check(p.canonical == v, "integer choice kept verbatim " + at);
                        check(std::atoi(v) >= s.ilo && std::atoi(v) <= s.ihi, "integer choice in range " + at);
                    } else {
                        check(word, "choice word " + at);
                    }
                    break;
                }
                case ts::env::Kind::Text:
                case ts::env::Kind::External:
                    check(p.canonical == v, "text kept verbatim " + at);
                    break;
            }
        }
    }

    // The T16 cases and the silent-default family: all rejected now.
    const std::vector<std::pair<const char*, const char*>> rejected = {
        {"TESSERACT_MIN_FALLBACK_DEST", "1e9"},       // atoi -> 1: guard silently equal to off
        {"TESSERACT_MIN_FALLBACK_DEST", "abc"},       // atoi -> 0
        {"TESSERACT_MIN_FALLBACK_DEST", "4294967296"},// atoi UB (0 on glibc): guard off
        {"TESSERACT_MIN_FALLBACK_DEST", "3000000000"},// atoi UB, wrapped to 1.8e19
        {"TESSERACT_MIN_FALLBACK_DEST", "-1"},        // wrapped to SIZE_MAX
        {"TESSERACT_MIN_FALLBACK_DEST", ""},          // atoi -> 0
        {"TESSERACT_GAPCLOSE_VOTES", "-1"},           // wrapped to 4294967295
        {"TESSERACT_GAPCLOSE", " 7x"},                // atoi -> 7
        {"TESSERACT_JOIN_TRACE", "true"},             // atoi -> 0: trace silently off
        {"TESSERACT_JOIN_TRACE", "yes"},
        {"TESSERACT_REQUIRE_SUPPORT_SINGLE", "on"},
        {"TESSERACT_RC_DOVETAIL", "true"},            // default-ON feature silently turned OFF
        {"TESSERACT_RC_DOVETAIL", ""},
        {"TESSERACT_EXACT_READ_THREADS", "false"},    // stayed ON
        {"TESSERACT_ROUTE_DISTANCE", "2"},            // silently OFF
        {"TESSERACT_OWNED_ANCHOR_PREFIX", "true"},
        {"TESSERACT_GAP_SEQUENCE_COMPETITION", "2"},
        {"TESSERACT_IS_VETO_CHR", "01"},              // first char '0' -> off
        {"TESSERACT_NO_PLASMID_VOUCH", "0"},          // presence: "0" turned it ON
        {"TESSERACT_COMMON_PREFIX", "inf"},
        {"TESSERACT_COMMON_PREFIX", "1e19"},          // cast to size_t undefined
        {"TESSERACT_TIP_ABS_MULT", "1e16"},
        {"TESSERACT_MODEL_MIN_PANEL", "-1"},          // double -> uint32 undefined
        {"TESSERACT_MODEL_MATCH", "Greedy"},          // silently mutual
        {"TESSERACT_COMEMBER_WEAK", "99999999999999999999"},
    };
    for (const auto& r : rejected) {
        const ts::env::Spec* s = ts::env::find(r.first);
        check(s != nullptr, std::string("registered ") + r.first);
        ts::env::Parsed p;
        std::string why;
        check(!ts::env::parse(*s, r.second, p, why), std::string("rejected ") + r.first + "='" + r.second + "'");
    }

    // Every value the campaign harness has ever passed to a release flag (xqueue jobs,
    // wrappers, sweeps; 2026-09-25 grep) and the K2 configuration must still be accepted.
    const std::vector<std::pair<const char*, const char*>> used = {
        {"TESSERACT_COMMON_PREFIX", "0"}, {"TESSERACT_COMMON_PREFIX", "250"},
        {"TESSERACT_COMMON_PREFIX", "1000"}, {"TESSERACT_COMMON_PREFIX", "3000"},
        {"TESSERACT_MIN_FALLBACK_DEST", "1000000000"}, {"TESSERACT_MIN_FALLBACK_DEST", "2000"},
        {"TESSERACT_REQUIRE_SUPPORT_SINGLE", "1"}, {"TESSERACT_EXCLUDE_SHARED_REPEAT_SUPPORT", "1"},
        {"TESSERACT_SHARED_SUPPORT_AUDIT", "1"}, {"TESSERACT_JOIN_TRACE", "1"},
        {"TESSERACT_GAPCLOSE", "10"}, {"TESSERACT_GAPCLOSE_VOTES", "2"}, {"TESSERACT_GAPCLOSE_VOTES", "3"},
        {"TESSERACT_GAP_ORIENTED", "1"}, {"TESSERACT_GAP_SEQUENCE_COMPETITION", "1"},
        {"TESSERACT_EXACT_READ_THREADS", "0"}, {"TESSERACT_OWNED_ANCHOR_PREFIX", "1"},
        {"TESSERACT_WEIGHTED_ELIGIBLE_COVERAGE", "1"}, {"TESSERACT_WEIGHTED_RESOLVER_COVERAGE", "0"},
        {"TESSERACT_DEBUG_RESOLVE", "1"}, {"TESSERACT_NO_UNSPANNED_FALLBACK", "1"},
        {"TESSERACT_MODEL_AUTHOR", "1"}, {"TESSERACT_MODEL_MIN_FRACTION", "0.1"},
        {"TESSERACT_MODEL_MIN_PANEL", "3"}, {"TESSERACT_NO_PLASMID_VOUCH", "1"},
        {"TESSERACT_RC_DOVETAIL", "0"}, {"TESSERACT_GRAPH_PHASES", "1"},
        {"TESSERACT_MODEL_DIR", "/home/x/models100/bundled"}, {"TESSERACT_EC_REQUIRE_UNIQUE_BEST", "1"},
    };
    for (const auto& u : used) {
        const ts::env::Spec* s = ts::env::find(u.first);
        check(s != nullptr, std::string("registered ") + u.first);
        ts::env::Parsed p;
        std::string why;
        check(ts::env::parse(*s, u.second, p, why), std::string("accepted ") + u.first + "='" + u.second + "' " + why);
    }
    check(ts::env::find("TESSERACT_NOT_A_FLAG") == nullptr, "unknown name is not in the table");
}

// ---- 2. readers exit 2 --------------------------------------------------------------------

void readerRejections() {
    struct Case { const char* name; const char* value; int reader; };  // 0 int, 1 real, 2 on, 3 present, 4 text
    const std::vector<Case> cases = {
        {"TESSERACT_MIN_FALLBACK_DEST", "1e9", 0}, {"TESSERACT_MIN_FALLBACK_DEST", "abc", 0},
        {"TESSERACT_MIN_FALLBACK_DEST", "4294967296", 0}, {"TESSERACT_MIN_FALLBACK_DEST", "-1", 0},
        {"TESSERACT_GAPCLOSE_VOTES", "-1", 0}, {"TESSERACT_JOIN_TRACE", "true", 2},
        {"TESSERACT_RC_DOVETAIL", "true", 2}, {"TESSERACT_COMMON_PREFIX", "abc", 1},
        {"TESSERACT_NO_PLASMID_VOUCH", "0", 3}, {"TESSERACT_MODEL_MATCH", "GREEDY", 4},
    };
    for (const Case& c : cases) {
        const bool rejected = exitsWithEnvError(c.name, [&] {
            ::setenv(c.name, c.value, 1);
            switch (c.reader) {
                case 0: (void)ts::env::integer(c.name, 0); break;
                case 1: (void)ts::env::real(c.name, 0.0); break;
                case 2: (void)ts::env::on(c.name, false); break;
                case 3: (void)ts::env::present(c.name); break;
                default: (void)ts::env::text(c.name); break;
            }
        });
        check(rejected, std::string("reader exits 2 on ") + c.name + "='" + c.value + "'");
    }
    // Defaults and valid values come straight back.
    clearTesseractEnvironment();
    check(ts::env::integer("TESSERACT_MIN_FALLBACK_DEST", 0) == 0, "unset -> default");
    check(ts::env::on("TESSERACT_RC_DOVETAIL", true), "unset default-on switch stays on");
    check(ts::env::on("TESSERACT_EXACT_READ_THREADS", true), "unset default-on binary stays on");
    check(ts::env::real("TESSERACT_COMMON_PREFIX", 3000.0) == 3000.0, "unset real -> default");
    check(!ts::env::present("TESSERACT_DEBUG_RESOLVE"), "unset presence is off");
    check(ts::env::text("TESSERACT_JOIN_DUMP") == nullptr, "unset text is null");
    ::setenv("TESSERACT_MIN_FALLBACK_DEST", "1000000000", 1);
    check(ts::env::integer("TESSERACT_MIN_FALLBACK_DEST", 0) == 1000000000, "K2 guard value");
    ::setenv("TESSERACT_JOIN_TRACE", "2", 1);
    check(ts::env::on("TESSERACT_JOIN_TRACE", false), "=2 is on, as atoi != 0 was");
    ::setenv("TESSERACT_RC_DOVETAIL", "0", 1);
    check(!ts::env::on("TESSERACT_RC_DOVETAIL", true), "=0 turns a default-on switch off");
    ::setenv("TESSERACT_DEBUG_RESOLVE", "0", 1);
    check(ts::env::present("TESSERACT_DEBUG_RESOLVE"), "presence flag: any value is on (log only)");

    // 3/4a. No caching at reader level: the second read sees the new value.
    ::setenv("TESSERACT_MIN_FALLBACK_DEST", "0", 1);
    check(ts::env::integer("TESSERACT_MIN_FALLBACK_DEST", 0) == 0, "second read sees the changed value");
    clearTesseractEnvironment();
}

// ---- 3. validateEnvironment ---------------------------------------------------------------

std::string validate(int& rc) {
    std::FILE* log = std::tmpfile();
    rc = ts::env::validateEnvironment(log);
    std::rewind(log);
    std::string out;
    char buf[4096];
    size_t n = 0;
    while ((n = std::fread(buf, 1, sizeof buf, log)) > 0) out.append(buf, n);
    std::fclose(log);
    return out;
}

void validation() {
    clearTesseractEnvironment();
    int rc = -1;
    check(validate(rc).empty() && rc == 0, "empty environment: no lines, status 0");

    // The K2 configuration plus one typo and one wrapper-script variable.
    ::setenv("TESSERACT_COMMON_PREFIX", "0", 1);
    ::setenv("TESSERACT_MIN_FALLBACK_DEST", "1000000000", 1);
    ::setenv("TESSERACT_REQUIRE_SUPPORT_SINGLE", "1", 1);
    ::setenv("TESSERACT_EXCLUDE_SHARED_REPEAT_SUPPORT", "1", 1);
    ::setenv("TESSERACT_MIN_FALLBAK_DEST", "5", 1);
    ::setenv("TESSERACT_ASM", "/usr/bin/tesseract-asm", 1);
    const std::string out = validate(rc);
    check(rc == 0, "K2 configuration validates");
    check(out ==
              "[config] TESSERACT_COMMON_PREFIX=0\n"
              "[config] TESSERACT_EXCLUDE_SHARED_REPEAT_SUPPORT=1 (on)\n"
              "[config] TESSERACT_MIN_FALLBACK_DEST=1000000000\n"
              "[config] warning: TESSERACT_MIN_FALLBAK_DEST is set but is not a TesserACT flag; it has no effect\n"
              "[config] TESSERACT_REQUIRE_SUPPORT_SINGLE=1 (on)\n",
          "[config] lines sorted, parsed values, typo warned, wrapper variable silent:\n" + out);

    ::setenv("TESSERACT_MIN_FALLBACK_DEST", "1e9", 1);
    ::setenv("TESSERACT_JOIN_TRACE", "true", 1);
    const std::string bad = validate(rc);
    check(rc == 2, "malformed values: status 2");
    check(bad ==
              "error: TESSERACT_JOIN_TRACE='true' is not valid: not a base-10 integer (expected an integer: 0 = off, any other value = on)\n"
              "error: TESSERACT_MIN_FALLBACK_DEST='1e9' is not valid: not a base-10 integer (expected a base-10 integer in [0, 2147483647])\n",
          "one error line per malformed flag, and no [config] lines:\n" + bad);
    clearTesseractEnvironment();
}

// ---- 4. source scan -----------------------------------------------------------------------

void scanSelfTest() {
    // The scanner must recognise the release pattern and pass its fixed form.
    const auto dir = std::filesystem::temp_directory_path() / ("test_envflags_scan_" + std::to_string(::getpid()));
    std::filesystem::create_directories(dir);
    std::ofstream(dir / "bad.cpp") <<
        "static const bool joinTrace_ = [] {\n"
        "    const char* e = std::getenv(\"TESSERACT_JOIN_TRACE\");\n"
        "    return e && std::atoi(e) != 0;\n"
        "}();\n"
        "static const size_t minDest{static_cast<size_t>(env::integer(\"TESSERACT_MIN_FALLBACK_DEST\", 0))};\n"
        "static int helper(int x) { return x + env::on(\"TESSERACT_JOIN_TRACE\", false); }\n"
        "static const double z = env::integer(\"TESSERACT_COMMON_PREFIX\", 0);\n"
        "const char* q = \"TESSERACT_NOT_REGISTERED\";\n";
    std::ofstream(dir / "good.cpp") <<
        "// static const bool x = getenv(\"TESSERACT_JOIN_TRACE\");  (a comment)\n"
        "static void computeContigStats(int a);\n"
        "static const char* kName = \"static getenv env::\";\n"
        "void f() { const bool joinTrace = env::on(\"TESSERACT_JOIN_TRACE\", false); (void)joinTrace; }\n"
        "static const int kPlain = 3;\n";
    for (bool strict : {false, true}) {
        const std::vector<std::string> found = envscan::scan(dir, strict);
        auto has = [&](const std::string& needle) {
            return std::any_of(found.begin(), found.end(), [&](const std::string& s) { return s.find(needle) != std::string::npos; });
        };
        const std::string mode = strict ? " (strict)" : " (phase 0)";
        check(has("bad.cpp:1: static caches") == strict, "static lambda over getenv flagged only in strict mode" + mode);
        check(has("bad.cpp:2: raw getenv") == strict, "raw getenv flagged only in strict mode" + mode);
        check(has("bad.cpp:5: static caches"), "brace-initialised static over env:: flagged" + mode);
        check(has("bad.cpp:7: static caches"), "static over env::integer flagged" + mode);
        check(has("bad.cpp:7: TESSERACT_COMMON_PREFIX is a real flag read through env::integer"), "kind mismatch flagged" + mode);
        check(has("bad.cpp:8: TESSERACT_NOT_REGISTERED is not in the flag table"), "unregistered literal flagged" + mode);
        check(!has("bad.cpp:6"), "a static function that reads a flag per call is fine" + mode);
        check(!has("good.cpp"), "comments, string contents and plain statics are not flagged" + mode);
    }
    std::filesystem::remove_all(dir);
}

}  // namespace

int main(int argc, char** argv) {
    if ((argc == 3 || argc == 4) && std::string(argv[1]) == "--scan") {
        const bool strict = argc == 4 && std::string(argv[3]) == "strict";
        const std::vector<std::string> v = envscan::scan(argv[2], strict);
        for (const std::string& s : v) std::printf("%s\n", s.c_str());
        std::printf("%zu violation(s) in %s%s\n", v.size(), argv[2], strict ? " (strict)" : "");
        return v.empty() ? 0 : 1;
    }
    clearTesseractEnvironment();
    parseContract();
    readerRejections();
    validation();
    scanSelfTest();

    const std::filesystem::path src = envscan::sourceDir(__FILE__);
    check(std::filesystem::exists(src / "envflags.cpp"),
          "source tree found at " + src.string() + " (run from the tree root)");
    const std::vector<std::string> v = envscan::scan(src, false);
    for (const std::string& s : v) std::fprintf(stderr, "  %s\n", s.c_str());
    check(v.empty(), std::to_string(v.size()) + " source-scan violation(s) in " + src.string());

    std::printf("env flags: %d checks passed\n", checks);
    return 0;
}
