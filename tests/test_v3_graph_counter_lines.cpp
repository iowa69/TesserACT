// T15 and the A2 counter lines of every G-graph fix (build_v3), end to end.
//
// Runs the real Assembler on a small synthetic paired library, once per configuration,
// each in a forked child (the release caches several flags in function-local statics, so
// one process can only ever see one configuration), with stderr captured. Every rung must
// print, zeros included:
//   [carrygate]  (T01)      [gapclose]  (T15/T08)      [hairpin_keep]  (T28)
//   [simplify_count_all] (T29)                         [joindeadends]  (T36)
// and the run prints [ec_cap_mask] (T30) once. With TESSERACT_GAPCLOSE set, the legacy
// "gap close:" line must appear on every rung too: release 1.3.0 printed it only when a
// rung had voted or closed something, so an idle rung looked like an absent feature.
// Uses only release APIs, so the same source fails on release 1.3.0.
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <random>
#include <sstream>
#include <string>
#include <sys/wait.h>
#include <unistd.h>
#include <vector>

#include "assembler.h"

extern char** environ;

namespace {

int failures = 0;
void check(bool ok, const std::string& what) {
    std::printf("  %s  %s\n", ok ? "PASS" : "FAIL", what.c_str());
    std::fflush(stdout);
    if (!ok) ++failures;
}

void clearTesseractEnv() {
    std::vector<std::string> names;
    for (char** e = environ; e && *e; ++e) {
        if (std::strncmp(*e, "TESSERACT_", 10) == 0) {
            const char* eq = std::strchr(*e, '=');
            names.emplace_back(*e, eq ? static_cast<size_t>(eq - *e) : std::strlen(*e));
        }
    }
    for (const std::string& n : names) unsetenv(n.c_str());
}

std::string rc(const std::string& s) { return ts::reverseComplement(s); }

std::string tmpDir;

void writeLibrary() {
    std::mt19937 rng(15);
    std::string genome(30000, 'A');
    for (char& c : genome) c = "ACGT"[rng() & 3];
    // Two copies of a 600 bp repeat and a thin stretch, so there is something to resolve and
    // something to leave open.
    const std::string rep = genome.substr(5000, 600);
    genome.replace(20000, 600, rep);
    std::ofstream f1(tmpDir + "/r1.fq"), f2(tmpDir + "/r2.fq");
    std::uniform_real_distribution<double> U(0, 1);
    const int L = 150, insert = 400, pairs = 4000;
    int written = 0;
    for (int i = 0; i < pairs; ++i) {
        const size_t p = rng() % (genome.size() - insert);
        if (p > 12000 && p < 12600 && U(rng) < 0.9) continue;   // coverage dropout
        std::string a = genome.substr(p, L);
        std::string b = rc(genome.substr(p + insert - L, L));
        for (char& c : a) if (U(rng) < 0.003) c = "ACGT"[rng() & 3];
        for (char& c : b) if (U(rng) < 0.003) c = "ACGT"[rng() & 3];
        f1 << "@p" << written << "/1\n" << a << "\n+\n" << std::string(a.size(), 'I') << "\n";
        f2 << "@p" << written << "/2\n" << b << "\n+\n" << std::string(b.size(), 'I') << "\n";
        ++written;
    }
}

// A circular replicon with no repeat: once cleaned, every rung is one closed unitig with no
// dead end, so the gap closer has nothing to vote on -- the idle case T15 is about.
void writeCircularLibrary() {
    std::mt19937 rng(151);
    std::string genome(12000, 'A');
    for (char& c : genome) c = "ACGT"[rng() & 3];
    const std::string ring = genome + genome;
    std::ofstream f1(tmpDir + "/c1.fq"), f2(tmpDir + "/c2.fq");
    const int L = 150, insert = 400, pairs = 3000;
    for (int i = 0; i < pairs; ++i) {
        const size_t p = rng() % genome.size();
        const std::string a = ring.substr(p, L);
        const std::string b = rc(ring.substr(p + insert - L, L));
        f1 << "@c" << i << "/1\n" << a << "\n+\n" << std::string(a.size(), 'I') << "\n";
        f2 << "@c" << i << "/2\n" << b << "\n+\n" << std::string(b.size(), 'I') << "\n";
    }
}

struct Run { bool ok = false; std::string log; };

Run assemble(const std::string& tag, const std::vector<std::pair<const char*, const char*>>& env,
             const std::string& lib1 = "r1.fq", const std::string& lib2 = "r2.fq") {
    const std::string logPath = tmpDir + "/" + tag + ".log";
    std::fflush(stdout);
    std::fflush(stderr);
    const pid_t pid = fork();
    if (pid == 0) {
        for (const auto& kv : env) setenv(kv.first, kv.second, 1);
        if (!std::freopen(logPath.c_str(), "w", stderr)) _exit(3);
        ts::AssemblyOptions opt;
        ts::Library lib;
        lib.r1 = tmpDir + "/" + lib1;
        lib.r2 = tmpDir + "/" + lib2;
        opt.libraries.push_back(lib);
        opt.outDir = tmpDir + "/out_" + tag;
        opt.kValues = {21, 33, 55};
        opt.userSetK = true;
        opt.threads = 2;
        opt.emitHtml = false;
        opt.verbose = true;
        ts::Assembler a(opt);
        std::string error;
        const bool ok = a.run(error);
        if (!ok) std::fprintf(stderr, "ERROR: %s\n", error.c_str());
        std::fflush(stderr);
        _exit(ok ? 0 : 1);
    }
    int status = 0;
    waitpid(pid, &status, 0);
    Run r;
    r.ok = WIFEXITED(status) && WEXITSTATUS(status) == 0;
    std::ifstream in(logPath);
    std::stringstream ss;
    ss << in.rdbuf();
    r.log = ss.str();
    return r;
}

size_t count(const std::string& log, const std::string& needle) {
    size_t n = 0;
    std::istringstream in(log);
    std::string line;
    while (std::getline(in, line)) n += line.find(needle) != std::string::npos;
    return n;
}
std::string firstLine(const std::string& log, const std::string& needle) {
    std::istringstream in(log);
    std::string line;
    while (std::getline(in, line))
        if (line.find(needle) != std::string::npos) return line;
    return "";
}

void perRung(const Run& r, const std::string& tag, const std::string& needle, size_t rungs) {
    const size_t n = count(r.log, needle);
    check(n == rungs, tag + ": '" + needle + "' on every rung (" + std::to_string(n) + "/" +
                          std::to_string(rungs) + ")");
}

}  // namespace

int main() {
    clearTesseractEnv();
    char tmpl[] = "/tmp/tess_v3_lines.XXXXXX";
    const char* d = mkdtemp(tmpl);
    if (!d) return 2;
    tmpDir = d;
    writeLibrary();
    writeCircularLibrary();

    struct Config { std::string tag; std::vector<std::pair<const char*, const char*>> env; };
    // 1.4.0: the umbrella TESSERACT_FIXES is on by default (src/defaults.h), so "unset" now
    // reports the fixes on; "release" (TESSERACT_FIXES=0) is the 1.3.0 state the other
    // configurations pin.
    const std::vector<Config> configs = {
        {"unset", {}},
        {"release", {{"TESSERACT_FIXES", "0"}}},
        {"gapclose_legacy", {{"TESSERACT_FIXES", "0"}, {"TESSERACT_GAPCLOSE", "10"}}},
        {"gapclose_oriented", {{"TESSERACT_FIXES", "0"}, {"TESSERACT_GAPCLOSE", "10"}, {"TESSERACT_GAP_ORIENTED", "1"}}},
        {"fixes", {{"TESSERACT_FIXES", "1"}, {"TESSERACT_GAPCLOSE", "10"}}},
    };
    for (const Config& c : configs) {
        std::printf("config %s\n", c.tag.c_str());
        const Run r = assemble(c.tag, c.env);
        check(r.ok, c.tag + ": assembly succeeded");
        const size_t rungs = count(r.log, " count ");
        check(rungs == 3, c.tag + ": three rungs counted");
        perRung(r, c.tag, "[carrygate] enabled=", rungs);
        perRung(r, c.tag, "[hairpin_keep] enabled=", rungs);
        perRung(r, c.tag, "[simplify_count_all] enabled=", rungs);
        perRung(r, c.tag, "[joindeadends] requested=", rungs);
        perRung(r, c.tag, "[gapclose] ", rungs);
        check(count(r.log, "[ec_cap_mask] enabled=") == 1, c.tag + ": [ec_cap_mask] once per run");
        const std::string cg = firstLine(r.log, "k=21  [carrygate]");
        check(cg.find("carriedKmers=0 ") != std::string::npos, c.tag + ": first rung carries nothing");
        const bool fixes = c.tag == "fixes" || c.tag == "unset";
        const std::string want = fixes ? "enabled=1" : "enabled=0";
        check(firstLine(r.log, "[carrygate]").find(want) != std::string::npos &&
                  firstLine(r.log, "[hairpin_keep]").find(want) != std::string::npos &&
                  firstLine(r.log, "[simplify_count_all]").find(want) != std::string::npos &&
                  firstLine(r.log, "[ec_cap_mask]").find(want) != std::string::npos,
              c.tag + ": fixes report " + want);
        if (c.tag == "unset" || c.tag == "release") {
            perRung(r, c.tag, fixes ? "[gapclose] off nominatorFix=1" : "[gapclose] off nominatorFix=0", rungs);
            check(count(r.log, "gap close:") == 0, c.tag + ": no legacy gap-close line when off");
        } else {
            perRung(r, c.tag, "gap close:", rungs);
            perRung(r, c.tag, "[gapclose] minInt=10 votes=2", rungs);
            const std::string gl = firstLine(r.log, "[gapclose] minInt=");
            std::printf("    %s\n", gl.c_str());
            const bool oriented = c.tag == "gapclose_oriented";
            check(gl.find(oriented ? "oriented=1" : "oriented=0") != std::string::npos &&
                      gl.find(fixes ? "nominatorFix=1" : "nominatorFix=0") != std::string::npos,
                  c.tag + ": [gapclose] reports the configuration in force");
        }
    }

    {
        std::printf("config gapclose_idle (circular replicon, nothing to close)\n");
        const Run r = assemble("gapclose_idle", {{"TESSERACT_GAPCLOSE", "10"}}, "c1.fq", "c2.fq");
        check(r.ok, "gapclose_idle: assembly succeeded");
        const size_t rungs = count(r.log, " count ");
        const size_t idle = count(r.log, "voted=0 nominated=0 tested=0 bridged=0");
        std::printf("    %s\n", firstLine(r.log, "[gapclose] minInt=").c_str());
        check(rungs == 3 && idle == rungs, "gapclose_idle: every rung idle (0 voted, 0 bridged)");
        perRung(r, "gapclose_idle", "gap close:", rungs);   // release 1.3.0: 0 of 3
        perRung(r, "gapclose_idle", "[gapclose] minInt=10 votes=2", rungs);
    }

    std::printf("test_v3_graph_counter_lines: %s (%d failure%s)\n", failures ? "FAILED" : "ALL PASSED",
                failures, failures == 1 ? "" : "s");
    const std::string cmd = "rm -rf '" + tmpDir + "'";
    if (std::system(cmd.c_str()) != 0) std::printf("note: could not remove %s\n", tmpDir.c_str());
    return failures ? 1 : 0;
}
