// Component test for the 1.4.0 defaults (RELEASE_PLAN 1.4.0 step 5.1): src/defaults.h.
//
// 1.4.0 makes the combo3 F2 configuration the default. This pins, with the environment cleared:
//   * each default-resolution helper returns the 1.4.0 value, and the fix switches that follow
//     the umbrella are on;
//   * each opt-out returns the 1.3.0 value, and the whole opt-out list gives the 1.3.0 state;
//   * an explicit per-fix TESSERACT_FIX_<NAME>=0 beats the umbrella default (and =1 beats
//     TESSERACT_FIXES=0);
//   * the run-mode presets are unchanged (fast 1.3, careful 1.4, aggressive 1.05);
//   * the text of the [defaults] line, from the helper and from the built binary.
// The binary is $TESSERACT_ASM_BIN, else ./tesseract-asm; `make componenttest` builds it.
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

#include <fcntl.h>
#include <sys/wait.h>
#include <unistd.h>

#include "assembler.h"
#include "defaults.h"
#include "emit_fixflags.h"
#include "graph_fix_flags.h"
#include "resolve.h"
#include "test_env.h"

namespace {

int checks = 0;
int failures = 0;

void check(bool ok, const std::string& what) {
    ++checks;
    if (!ok) {
        ++failures;
        std::fprintf(stderr, "FAIL: %s\n", what.c_str());
    }
}

using Env = std::vector<std::pair<std::string, std::string>>;

void setAll(const Env& env) {
    for (const auto& kv : env) setenv(kv.first.c_str(), kv.second.c_str(), 1);
}
void unsetAll(const Env& env) {
    for (const auto& kv : env) unsetenv(kv.first.c_str());
}

// The opt-out list that restores 1.3.0 (plus --tie-ratio 1.02 on the command line).
const Env kOptOut130 = {
    {"TESSERACT_FIXES", "0"},
    {"TESSERACT_COMMON_PREFIX", "3000"},
    {"TESSERACT_PREFIX_MIN_BODY", "0"},
    {"TESSERACT_MIN_FALLBACK_DEST", "0"},
    {"TESSERACT_REQUIRE_SUPPORT_SINGLE", "0"},
    {"TESSERACT_EXCLUDE_SHARED_REPEAT_SUPPORT", "0"},
    {"TESSERACT_DROPOUT_BRIDGE", "0"},
};

const char* kLine14 =
    "[defaults] tie_ratio=3.000 common_prefix=250 prefix_min_body=reach min_fallback_dest=1000000000 "
    "require_support_single=1 exclude_shared_repeat_support=1 fixes=1 dropout_bridge=1";
const char* kLine130 =
    "[defaults] tie_ratio=1.020 common_prefix=3000 prefix_min_body=0 min_fallback_dest=0 "
    "require_support_single=0 exclude_shared_repeat_support=0 fixes=0 dropout_bridge=0";

double modeTie(ts::RunMode m, bool userSet = false, double user = 0) {
    ts::AssemblyOptions o;
    o.mode = m;
    if (userSet) { o.tieRatio = user; o.userSetTie = true; }
    o.applyMode();
    return o.tieRatio;
}

struct Run {
    int rc = -1;
    std::string out, err;
};

// The binary with PATH, HOME and exactly `env`; stdout and stderr captured.
Run runBinary(const std::string& bin, const Env& env, const std::vector<std::string>& args, const std::string& dir) {
    const std::string errPath = dir + "/stderr.txt", outPath = dir + "/stdout.txt";
    std::vector<std::string> envStore;
    if (const char* p = std::getenv("PATH")) envStore.push_back(std::string("PATH=") + p);
    if (const char* h = std::getenv("HOME")) envStore.push_back(std::string("HOME=") + h);
    for (const auto& kv : env) envStore.push_back(kv.first + "=" + kv.second);
    std::vector<char*> envp;
    for (std::string& s : envStore) envp.push_back(const_cast<char*>(s.c_str()));
    envp.push_back(nullptr);
    std::vector<char*> argv;
    argv.push_back(const_cast<char*>(bin.c_str()));
    for (const std::string& a : args) argv.push_back(const_cast<char*>(a.c_str()));
    argv.push_back(nullptr);
    Run r;
    const pid_t pid = fork();
    if (pid == 0) {
        const int o = open(outPath.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
        const int e = open(errPath.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
        if (o >= 0) dup2(o, 1);
        if (e >= 0) dup2(e, 2);
        execve(bin.c_str(), argv.data(), envp.data());
        _exit(127);
    }
    int status = 0;
    if (pid < 0 || waitpid(pid, &status, 0) != pid) return r;
    r.rc = WIFEXITED(status) ? WEXITSTATUS(status) : 128 + WTERMSIG(status);
    std::ifstream f(errPath), g(outPath);
    r.err.assign(std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>());
    r.out.assign(std::istreambuf_iterator<char>(g), std::istreambuf_iterator<char>());
    return r;
}

bool hasLine(const std::string& text, const std::string& line) {
    return text.find(line + "\n") != std::string::npos;
}

}  // namespace

int main() {
    // Read before the environment is cleared: TESSERACT_ASM_BIN is itself a TESSERACT_* name.
    const char* e = std::getenv("TESSERACT_ASM_BIN");
    std::string bin = e && *e ? e : "./tesseract-asm";
    testenv::clearTesseractEnv();
    namespace d = ts::defaults;

    // ---- 1. the 1.4.0 defaults, nothing set
    check(d::kDefaultTieRatio == 3.0, "kDefaultTieRatio is 3.0");
    check(ts::AssemblyOptions().tieRatio == d::kDefaultTieRatio, "AssemblyOptions().tieRatio is the 1.4.0 default 3.0");
    check(d::commonPrefix() == 250.0, "common-prefix budget defaults to 250");
    check(d::prefixMinBodyText() == "reach" && d::prefixMinBody() == -1, "prefix min body defaults to reach");
    check(d::minFallbackDest() == 1000000000u, "min fallback dest defaults to 1000000000 (fallbacks off)");
    check(d::requireSupportSingle(), "require-support-single defaults on");
    check(d::excludeSharedRepeatSupport(), "exclude-shared-repeat-support defaults on");
    check(d::fixesUmbrella(), "the TESSERACT_FIXES umbrella defaults on");
    check(d::dropoutBridge(), "the dropout bridge defaults on");
    check(ts::fixflags::umbrella() && ts::emitfix::umbrella(), "the graph and emit umbrella helpers follow the default");
    check(ts::fixflags::fixLevel("TESSERACT_FIX_CARRY_READ_GATE", 2) == 1, "T01 carry-read gate on at level 1 by default");
    check(ts::fixflags::fixLevel("TESSERACT_FIX_HAIRPIN_KEEP") == 1, "T28 on by default");
    check(ts::emitfix::enabled(ts::emitfix::kPolishSkipN), "T02 polish skip-N on by default");
    check(ts::emitfix::enabled(ts::emitfix::kTrimCopyGuard), "T22 trim copy guard on by default");
    check(ts::emitfix::backoffBases() == ts::emitfix::kDefaultBackoff, "T11 back-off at its umbrella default");
    check(ts::gapFlankFixEnabled(), "T03 gap-flank restore on by default");
    check(ts::checkFixCoupling(stderr) == 0, "the default configuration passes the T02/T03 coupling check");
    check(d::line(d::kDefaultTieRatio) == kLine14, "[defaults] line text, 1.4.0 defaults (got " + d::line(3.0) + ")");

    // ---- 2. run-mode presets unchanged (D-e): careful is no longer stricter than standard on tau
    check(modeTie(ts::RunMode::Standard) == 3.0, "--mode standard keeps tie ratio 3.0");
    check(std::fabs(modeTie(ts::RunMode::Fast) - 1.3) < 1e-12, "--mode fast keeps tie ratio 1.3");
    check(std::fabs(modeTie(ts::RunMode::Careful) - 1.4) < 1e-12, "--mode careful keeps tie ratio 1.4");
    check(std::fabs(modeTie(ts::RunMode::Aggressive) - 1.05) < 1e-12, "--mode aggressive keeps tie ratio 1.05");
    check(std::fabs(modeTie(ts::RunMode::Careful, true, 1.02) - 1.02) < 1e-12, "--tie-ratio beats a mode preset");

    // ---- 3. each opt-out gives the 1.3.0 value
    setenv("TESSERACT_COMMON_PREFIX", "3000", 1);
    check(d::commonPrefix() == 3000.0, "TESSERACT_COMMON_PREFIX=3000 restores 3000");
    unsetenv("TESSERACT_COMMON_PREFIX");
    setenv("TESSERACT_PREFIX_MIN_BODY", "0", 1);
    check(d::prefixMinBody() == 0 && d::prefixMinBodyText() == "0", "TESSERACT_PREFIX_MIN_BODY=0 turns ENDS-A off");
    setenv("TESSERACT_PREFIX_MIN_BODY", "400", 1);
    check(d::prefixMinBody() == 400, "TESSERACT_PREFIX_MIN_BODY=400 is 400 bp");
    unsetenv("TESSERACT_PREFIX_MIN_BODY");
    setenv("TESSERACT_MIN_FALLBACK_DEST", "0", 1);
    check(d::minFallbackDest() == 0, "TESSERACT_MIN_FALLBACK_DEST=0 restores the fallbacks");
    unsetenv("TESSERACT_MIN_FALLBACK_DEST");
    setenv("TESSERACT_REQUIRE_SUPPORT_SINGLE", "0", 1);
    check(!d::requireSupportSingle(), "TESSERACT_REQUIRE_SUPPORT_SINGLE=0 turns it off");
    unsetenv("TESSERACT_REQUIRE_SUPPORT_SINGLE");
    setenv("TESSERACT_EXCLUDE_SHARED_REPEAT_SUPPORT", "0", 1);
    check(!d::excludeSharedRepeatSupport(), "TESSERACT_EXCLUDE_SHARED_REPEAT_SUPPORT=0 turns it off");
    unsetenv("TESSERACT_EXCLUDE_SHARED_REPEAT_SUPPORT");
    setenv("TESSERACT_DROPOUT_BRIDGE", "0", 1);
    check(!d::dropoutBridge(), "TESSERACT_DROPOUT_BRIDGE=0 turns the bridge off");
    unsetenv("TESSERACT_DROPOUT_BRIDGE");
    setenv("TESSERACT_FIXES", "0", 1);
    check(!d::fixesUmbrella() && !ts::fixflags::umbrella() && !ts::emitfix::umbrella(), "TESSERACT_FIXES=0 turns the umbrella off");
    check(ts::fixflags::fixLevel("TESSERACT_FIX_CARRY_READ_GATE", 2) == 0, "TESSERACT_FIXES=0: T01 off");
    check(!ts::emitfix::enabled(ts::emitfix::kPolishSkipN) && ts::emitfix::backoffBases() == 0 && !ts::gapFlankFixEnabled(),
          "TESSERACT_FIXES=0: T02, T11 and T03 off");
    unsetenv("TESSERACT_FIXES");

    setAll(kOptOut130);
    check(d::line(1.02) == kLine130, "[defaults] line text, the 1.3.0 opt-out list (got " + d::line(1.02) + ")");
    check(ts::checkFixCoupling(stderr) == 0, "the 1.3.0 opt-out passes the T02/T03 coupling check");
    unsetAll(kOptOut130);
    check(d::line(3.0) == kLine14, "unsetting the opt-outs gives the 1.4.0 line back");

    // ---- 4. an explicit per-fix value beats the umbrella, in both directions
    setenv("TESSERACT_FIX_CARRY_READ_GATE", "0", 1);
    check(ts::fixflags::fixLevel("TESSERACT_FIX_CARRY_READ_GATE", 2) == 0 && d::fixesUmbrella(),
          "TESSERACT_FIX_CARRY_READ_GATE=0 turns T01 off under the default umbrella");
    check(ts::fixflags::fixLevel("TESSERACT_FIX_HAIRPIN_KEEP") == 1, "... and leaves the other graph fixes on");
    unsetenv("TESSERACT_FIX_CARRY_READ_GATE");
    setenv("TESSERACT_FIX_TRIM_COPY_GUARD", "0", 1);
    check(!ts::emitfix::enabled(ts::emitfix::kTrimCopyGuard) && ts::emitfix::enabled(ts::emitfix::kPolishSkipN),
          "TESSERACT_FIX_TRIM_COPY_GUARD=0 turns T22 alone off");
    unsetenv("TESSERACT_FIX_TRIM_COPY_GUARD");
    setenv("TESSERACT_FIX_GAP_FLANK", "0", 1);
    check(!ts::gapFlankFixEnabled(), "TESSERACT_FIX_GAP_FLANK=0 turns T03 off under the default umbrella");
    unsetenv("TESSERACT_FIX_GAP_FLANK");
    setenv("TESSERACT_FIX_GAPFILL_BACKOFF", "0", 1);
    check(ts::emitfix::backoffBases() == 0, "TESSERACT_FIX_GAPFILL_BACKOFF=0 turns T11 off under the default umbrella");
    unsetenv("TESSERACT_FIX_GAPFILL_BACKOFF");
    setenv("TESSERACT_FIXES", "0", 1);
    setenv("TESSERACT_FIX_CARRY_READ_GATE", "1", 1);
    check(ts::fixflags::fixLevel("TESSERACT_FIX_CARRY_READ_GATE", 2) == 1, "TESSERACT_FIX_CARRY_READ_GATE=1 wins over TESSERACT_FIXES=0");
    unsetenv("TESSERACT_FIX_CARRY_READ_GATE");
    unsetenv("TESSERACT_FIXES");

    // ---- 5. the binary prints the line before any work (a missing input then stops it)
    if (access(bin.c_str(), X_OK) != 0) {
        check(false, "the binary " + bin + " is built (make componenttest builds it)");
    } else {
        bin = std::filesystem::absolute(bin).string();
        const char* tmpdir = std::getenv("TMPDIR");
        std::string pattern = std::string(tmpdir && *tmpdir ? tmpdir : "/tmp") + "/tesseract_v14_defaults_XXXXXX";
        std::vector<char> buf(pattern.begin(), pattern.end());
        buf.push_back('\0');
        if (!mkdtemp(buf.data())) { std::perror("mkdtemp"); return 2; }
        const std::string dir = buf.data();
        const std::string missing = dir + "/missing_1.fq";
        const std::vector<std::string> base = {"-1", missing, "-o", dir + "/out", "-q"};

        Run r = runBinary(bin, {}, base, dir);
        check(r.rc == 1 && hasLine(r.err, kLine14), "binary, nothing set: prints the 1.4.0 [defaults] line (rc " +
                                                        std::to_string(r.rc) + ", stderr: " + r.err + ")");
        std::vector<std::string> a130 = base;
        a130.insert(a130.begin(), {"--tie-ratio", "1.02"});
        r = runBinary(bin, kOptOut130, a130, dir);
        check(r.rc == 1 && hasLine(r.err, kLine130), "binary, opt-out list + --tie-ratio 1.02: prints the 1.3.0 line (stderr: " +
                                                         r.err + ")");
        std::vector<std::string> careful = base;
        careful.insert(careful.begin(), {"--mode", "careful"});
        r = runBinary(bin, {}, careful, dir);
        check(r.err.find("[defaults] tie_ratio=1.400 ") != std::string::npos,
              "binary, --mode careful: the line shows the effective tie ratio 1.400 (stderr: " + r.err + ")");
        r = runBinary(bin, {}, {"--help"}, dir);
        check(r.rc == 0, "--help exits 0");
        // --help states the new tie-ratio default and the exact 1.3.0 opt-out list.
        for (const char* want : {"(default: 3.0;", "DEFAULTS (1.4.0)", "REPRODUCING 1.3.0",
                                 "TESSERACT_FIXES=0 TESSERACT_COMMON_PREFIX=3000 TESSERACT_PREFIX_MIN_BODY=0",
                                 "TESSERACT_MIN_FALLBACK_DEST=0 TESSERACT_REQUIRE_SUPPORT_SINGLE=0",
                                 "TESSERACT_EXCLUDE_SHARED_REPEAT_SUPPORT=0 TESSERACT_DROPOUT_BRIDGE=0",
                                 "tesseract-asm --tie-ratio 1.02 ..."})
            check(r.out.find(want) != std::string::npos, std::string("--help says '") + want + "'");
        std::filesystem::remove_all(dir);
    }

    if (failures) {
        std::fprintf(stderr, "test_v14_defaults: %d of %d checks FAILED\n", failures, checks);
        return 1;
    }
    std::printf("test_v14_defaults: %d checks passed\n", checks);
    return 0;
}
