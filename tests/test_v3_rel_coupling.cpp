// Component test for N9 (combo3 DEFECTS.md; RELEASE_PLAN 1.4.0 step 2): the T02/T03 coupling.
//
// T03 (TESSERACT_FIX_GAP_FLANK, alias TESSERACT_GAP_KEEP_FLANK, or the umbrella TESSERACT_FIXES)
// restores the k-1 bases after each open N-gap and writes the estimated gap; T02
// (TESSERACT_FIX_POLISH_SKIP_N) keeps the polisher off those N-runs. Before 1.4.0 an
// environment with T03 on and T02 off was neither refused nor warned. It is now refused at
// startup with exit status 2; every other combination runs.
//
// Two levels, each in forked children so that every case starts from its own environment:
//   1. ts::checkFixCoupling() under each flag combination (the rule itself);
//   2. the built binary (main() owns the call): T03 on + T02 off exits 2 and writes nothing,
//      both on and both off assemble a small synthetic library and exit 0.
// The binary is $TESSERACT_ASM_BIN, else ./tesseract-asm (make runs tests from the tree root);
// `make componenttest` builds it first (N14).
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

#include <fcntl.h>
#include <sys/wait.h>
#include <unistd.h>

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

std::string show(const Env& env) {
    if (env.empty()) return "(no TESSERACT_* set)";
    std::string s;
    for (const auto& kv : env) s += (s.empty() ? "" : " ") + kv.first + "=" + kv.second;
    return s;
}

std::string slurp(const std::string& path) {
    std::ifstream f(path);
    return std::string(std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>());
}

// Level 1: checkFixCoupling() in a child whose environment is exactly `env`.
int couplingIn(const Env& env, const std::string& errPath) {
    const pid_t pid = fork();
    if (pid == 0) {
        testenv::clearTesseractEnv();
        for (const auto& kv : env) setenv(kv.first.c_str(), kv.second.c_str(), 1);
        std::FILE* log = std::fopen(errPath.c_str(), "w");
        const int rc = ts::checkFixCoupling(log ? log : stderr);
        if (log) std::fclose(log);
        _exit(rc);
    }
    int status = 0;
    if (pid < 0 || waitpid(pid, &status, 0) != pid) return -1;
    return WIFEXITED(status) ? WEXITSTATUS(status) : 128 + WTERMSIG(status);
}

// Level 2: the binary, with an environment of PATH, HOME and exactly `env`.
int runBinary(const std::string& bin, const Env& env, const std::vector<std::string>& args,
              const std::string& errPath) {
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
    const pid_t pid = fork();
    if (pid == 0) {
        const int devnull = open("/dev/null", O_WRONLY);
        const int e = open(errPath.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
        if (devnull >= 0) dup2(devnull, 1);
        if (e >= 0) dup2(e, 2);
        execve(bin.c_str(), argv.data(), envp.data());
        _exit(127);
    }
    int status = 0;
    if (pid < 0 || waitpid(pid, &status, 0) != pid) return -1;
    return WIFEXITED(status) ? WEXITSTATUS(status) : 128 + WTERMSIG(status);
}

// A 6 kb random genome and error-free 2x100 pairs (insert 300) at about 40x, deterministic.
void writeLibrary(const std::string& r1, const std::string& r2) {
    unsigned long long x = 0x9E3779B97F4A7C15ULL;
    auto next = [&x]() { x ^= x << 13; x ^= x >> 7; x ^= x << 17; return x; };
    std::string g(6000, 'A');
    for (char& c : g) c = "ACGT"[next() % 4];
    auto rc = [](const std::string& s) {
        std::string r(s.rbegin(), s.rend());
        for (char& c : r) c = c == 'A' ? 'T' : c == 'C' ? 'G' : c == 'G' ? 'C' : 'A';
        return r;
    };
    std::ofstream f1(r1), f2(r2);
    const int readLen = 100, insert = 300, pairs = 1200;
    for (int i = 0; i < pairs; ++i) {
        const size_t s = next() % (g.size() - insert);
        const std::string frag = g.substr(s, insert);
        const std::string a = frag.substr(0, readLen), b = rc(frag.substr(insert - readLen));
        const std::string q(readLen, 'I');
        f1 << "@p" << i << "/1\n" << a << "\n+\n" << q << "\n";
        f2 << "@p" << i << "/2\n" << b << "\n+\n" << q << "\n";
    }
}

}  // namespace

int main() {
    // Read before the environment is cleared: TESSERACT_ASM_BIN is itself a TESSERACT_* name.
    const char* e = std::getenv("TESSERACT_ASM_BIN");
    std::string bin = e && *e ? e : "./tesseract-asm";
    testenv::clearTesseractEnv();
    const char* tmpdir = std::getenv("TMPDIR");
    std::string pattern = std::string(tmpdir && *tmpdir ? tmpdir : "/tmp") + "/tesseract_v3_coupling_XXXXXX";
    std::vector<char> buf(pattern.begin(), pattern.end());
    buf.push_back('\0');
    if (!mkdtemp(buf.data())) { std::perror("mkdtemp"); return 2; }
    const std::string dir = buf.data();
    const std::string err = dir + "/stderr.txt";
    const char* kMessage = "TESSERACT_FIX_GAP_FLANK needs TESSERACT_FIX_POLISH_SKIP_N";

    // ---- level 1: the rule
    struct Case { Env env; int rc; };
    const std::vector<Case> cases = {
        {{}, 0},                                                                  // library defaults
        {{{"TESSERACT_FIXES", "1"}}, 0},                                          // both on (umbrella)
        {{{"TESSERACT_FIXES", "0"}}, 0},                                          // both off
        {{{"TESSERACT_FIX_GAP_FLANK", "1"}, {"TESSERACT_FIX_POLISH_SKIP_N", "0"}}, 2},
        {{{"TESSERACT_GAP_KEEP_FLANK", "1"}, {"TESSERACT_FIX_POLISH_SKIP_N", "0"}}, 2},   // package alias
        {{{"TESSERACT_FIXES", "1"}, {"TESSERACT_FIX_POLISH_SKIP_N", "0"}}, 2},
        {{{"TESSERACT_FIXES", "0"}, {"TESSERACT_FIX_GAP_FLANK", "1"}}, 2},
        {{{"TESSERACT_FIX_GAP_FLANK", "1"}, {"TESSERACT_FIX_POLISH_SKIP_N", "1"}}, 0},   // both on
        {{{"TESSERACT_FIX_GAP_FLANK", "0"}, {"TESSERACT_FIX_POLISH_SKIP_N", "0"}}, 0},   // both off
        {{{"TESSERACT_FIX_GAP_FLANK", "0"}, {"TESSERACT_FIX_POLISH_SKIP_N", "1"}}, 0},   // T02 alone
        {{{"TESSERACT_FIXES", "1"}, {"TESSERACT_FIX_GAP_FLANK", "0"}, {"TESSERACT_FIX_POLISH_SKIP_N", "0"}}, 0},
        {{{"TESSERACT_FIXES", "1"}, {"TESSERACT_GAP_KEEP_FLANK", "0"}, {"TESSERACT_FIX_POLISH_SKIP_N", "0"}}, 0},
    };
    for (const Case& c : cases) {
        const int rc = couplingIn(c.env, err);
        check(rc == c.rc, "checkFixCoupling under " + show(c.env) + " returns " + std::to_string(c.rc) +
                              " (got " + std::to_string(rc) + ")");
        const bool said = slurp(err).find(kMessage) != std::string::npos;
        check(said == (c.rc == 2), "checkFixCoupling under " + show(c.env) +
                                       (c.rc == 2 ? " names both flags" : " prints nothing"));
    }

    // ---- level 2: main() refuses before any work, and runs the accepted combinations
    if (access(bin.c_str(), X_OK) != 0) {
        check(false, "the binary " + bin + " is built (make componenttest builds it; N14)");
    } else {
        bin = std::filesystem::absolute(bin).string();
        const std::string r1 = dir + "/r_1.fq", r2 = dir + "/r_2.fq";
        writeLibrary(r1, r2);
        struct Run { Env env; int rc; const char* label; };
        const std::vector<Run> runs = {
            {{{"TESSERACT_FIX_GAP_FLANK", "1"}, {"TESSERACT_FIX_POLISH_SKIP_N", "0"}}, 2, "t03_only"},
            {{{"TESSERACT_GAP_KEEP_FLANK", "1"}, {"TESSERACT_FIX_POLISH_SKIP_N", "0"}}, 2, "alias_only"},
            {{{"TESSERACT_FIXES", "1"}, {"TESSERACT_FIX_POLISH_SKIP_N", "0"}}, 2, "umbrella_no_t02"},
            {{{"TESSERACT_FIX_GAP_FLANK", "1"}, {"TESSERACT_FIX_POLISH_SKIP_N", "1"}}, 0, "both_on"},
            {{{"TESSERACT_FIX_GAP_FLANK", "0"}, {"TESSERACT_FIX_POLISH_SKIP_N", "0"}}, 0, "both_off"},
            {{}, 0, "defaults"},
        };
        for (const Run& r : runs) {
            const std::string out = dir + "/out_" + r.label;
            const int rc = runBinary(bin, r.env, {"-1", r1, "-2", r2, "-o", out, "-t", "2", "-q"}, err);
            const std::string text = slurp(err);
            check(rc == r.rc, std::string("tesseract-asm under ") + show(r.env) + " exits " + std::to_string(r.rc) +
                                  " (got " + std::to_string(rc) + "; stderr tail: " +
                                  text.substr(text.size() > 300 ? text.size() - 300 : 0) + ")");
            const bool wrote = std::filesystem::exists(out + "/contigs.fasta");
            if (r.rc == 2) {
                check(text.find(kMessage) != std::string::npos,
                      std::string("tesseract-asm under ") + show(r.env) + " names both flags");
                check(!wrote, std::string("tesseract-asm under ") + show(r.env) + " writes no assembly");
            } else {
                check(wrote && std::filesystem::file_size(out + "/contigs.fasta") > 0,
                      std::string("tesseract-asm under ") + show(r.env) + " writes contigs.fasta");
                check(text.find(kMessage) == std::string::npos,
                      std::string("tesseract-asm under ") + show(r.env) + " prints no coupling error");
            }
        }
    }

    std::filesystem::remove_all(dir);
    if (failures) {
        std::fprintf(stderr, "test_v3_rel_coupling (N9): %d of %d checks FAILED\n", failures, checks);
        return 1;
    }
    std::printf("test_v3_rel_coupling (N9): %d checks passed\n", checks);
    return 0;
}
