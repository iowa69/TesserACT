// Component test for command-line validation (combo3 G-io: T32) and the end-to-end refusal of a
// truncated gzip input (T18). It drives the built binary, because main() owns the parser.
//
// The binary is $TESSERACT_ASM_BIN, else ./tesseract-asm (make runs tests from the tree root).
// `make componenttest` builds it first (N14); run by hand without it, the test says SKIP and
// exits 0 -- build it first (make tesseract-asm) for the checks to run.
//
// Release 1.3.0 exited 0 on every "refused" case below: it truncated fractional integers, ran a
// -k ladder in the order given (correcting at its first rung), let a second -1/-2/--12 or -k
// replace the first, and loaded the same file twice when it was named twice.
//
// Parse-level cases end in --version, which exits 0 as soon as the parser reaches it: a refused
// option fails before that point, and an accepted one never starts an assembly.
#include <zlib.h>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include <fcntl.h>
#include <sys/wait.h>
#include <unistd.h>
#include "test_env.h"   // build_v3: T40, clear the ambient TESSERACT_* first

namespace {

int checks = 0;
int failures = 0;
std::string bin;
std::string dir;

void check(bool ok, const std::string& what) {
    ++checks;
    if (!ok) {
        ++failures;
        std::fprintf(stderr, "FAIL: %s\n", what.c_str());
    }
}

struct Run {
    int rc = -1;
    std::string err;
};

Run run(const std::vector<std::string>& args) {
    const std::string errPath = dir + "/stderr.txt";
    std::vector<char*> argv;
    argv.push_back(const_cast<char*>(bin.c_str()));
    for (const std::string& a : args) argv.push_back(const_cast<char*>(a.c_str()));
    argv.push_back(nullptr);
    Run r;
    const pid_t pid = fork();
    if (pid == 0) {
        const int devnull = open("/dev/null", O_WRONLY);
        const int e = open(errPath.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
        if (devnull >= 0) dup2(devnull, 1);
        if (e >= 0) dup2(e, 2);
        execv(bin.c_str(), argv.data());
        _exit(127);
    }
    int status = 0;
    if (pid < 0 || waitpid(pid, &status, 0) != pid) return r;
    r.rc = WIFEXITED(status) ? WEXITSTATUS(status) : 128 + WTERMSIG(status);
    std::ifstream f(errPath);
    r.err.assign(std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>());
    return r;
}

std::string show(const std::vector<std::string>& args) {
    std::string s;
    for (const std::string& a : args) {
        if (!s.empty()) s += ' ';
        s += a.rfind(dir, 0) == 0 ? a.substr(dir.size() + 1) : a;
    }
    return s;
}

void refused(const std::vector<std::string>& args, int rc, const char* needle) {
    const Run r = run(args);
    check(r.rc == rc, "'" + show(args) + "' exits " + std::to_string(rc) + " (got " + std::to_string(r.rc) +
                          ", stderr: " + r.err + ")");
    check(r.err.find(needle) != std::string::npos,
          "'" + show(args) + "' says '" + needle + "' (stderr: " + r.err + ")");
}

void accepted(const std::vector<std::string>& args) {
    const Run r = run(args);
    check(r.rc == 0, "'" + show(args) + "' is accepted (got " + std::to_string(r.rc) + ", stderr: " + r.err + ")");
}

std::string fastq(const char* mate, int n) {
    std::string s;
    for (int i = 0; i < n; ++i) {
        s += "@r" + std::to_string(i) + mate + "\nACGTACGTAC\n+\nIIIIIIIIII\n";
    }
    return s;
}

}  // namespace

int main() {
    testenv::clearTesseractEnv();
    const char* e = std::getenv("TESSERACT_ASM_BIN");
    bin = e && *e ? e : "./tesseract-asm";
    if (access(bin.c_str(), X_OK) != 0) {
        std::printf("test_v3_io_cli: SKIP -- %s is not built (run 'make tesseract-asm' first, or set "
                    "TESSERACT_ASM_BIN)\n", bin.c_str());
        return 0;
    }
    bin = std::filesystem::absolute(bin).string();
    const char* tmpdir = std::getenv("TMPDIR");
    std::string pattern = std::string(tmpdir && *tmpdir ? tmpdir : "/tmp") + "/tesseract_v3_cli_XXXXXX";
    std::vector<char> buf(pattern.begin(), pattern.end());
    buf.push_back('\0');
    if (!mkdtemp(buf.data())) { std::perror("mkdtemp"); return 2; }
    dir = buf.data();

    const std::string a = dir + "/a_1.fq", b = dir + "/a_2.fq", link = dir + "/link_1.fq";
    std::ofstream(a) << fastq("/1", 50);
    std::ofstream(b) << fastq("/2", 50);
    std::filesystem::create_symlink(a, link);
    const std::string out = dir + "/out";
    const std::string V = "--version";

    // ---- T32: integers, -k ladder
    refused({"-k", "21.9,55", V}, 2, "expects an integer");
    refused({"-k", "77,21,55", V}, 2, "strictly increasing");
    refused({"-k", "127,21", V}, 2, "strictly increasing");
    refused({"-k", "21,21", V}, 2, "strictly increasing");
    refused({"-k", "21,33", "-k", "55", V}, 2, "more than once");
    refused({"-t", "4.5", V}, 2, "expects an integer");
    refused({"--min-contig", "500.7", V}, 2, "expects an integer");
    refused({"-c", "2.5", V}, 2, "expects an integer");
    accepted({"-k", "21,33,55,77", V});
    accepted({"-k", "21", V});
    accepted({"-t", "4", V});
    accepted({"-t", "4.0", V});
    accepted({"--min-contig", "1e3", V});
    accepted({"--tie-ratio", "2.0", V});       // a real-valued option keeps its fraction
    refused({"-k", "22", V}, 1, "must be odd");  // unchanged release behaviour
    refused({"-k", "129", V}, 2, "must be between");

    // ---- T32: read options
    refused({"-1", a, "-1", b, V}, 2, "-1 was already given");
    refused({"-2", a, "-2", b, V}, 2, "-2 was already given");
    refused({"--12", a, "-2", b, V}, 2, "--12 already gives both mates");
    refused({"-2", b, "--12", a, V}, 2, "-2 was already given");
    refused({"-1", a, "--12", b, V}, 2, "-1 was already given");
    refused({"--12", a, "-1", b, V}, 2, "--12 was already given");
    refused({"--12", a, "--12", b, V}, 2, "--12 was already given");
    refused({"--read1", a, "--read1", b, V}, 2, "-1 was already given");
    accepted({"-1", a, "-2", b, V});
    accepted({"-2", b, "-1", a, V});
    accepted({"--12", a, V});
    accepted({"-s", a, "-s", b, V});
    accepted({"-1", a, "-2", b, "-s", link + "-not-there", V});  // file checks come after parsing

    // ---- T32: the same file twice (checked after parsing, before any read is loaded)
    refused({"-1", a, "-2", a, "-o", out, "-q"}, 2, "name the same file");
    refused({"-1", a, "-2", link, "-o", out, "-q"}, 2, "name the same file");
    refused({"-1", a, "-2", dir + "/./a_1.fq", "-o", out, "-q"}, 2, "name the same file");
    refused({"-s", a, "-s", a, "-o", out, "-q"}, 2, "name the same file");
    refused({"-1", a, "-2", b, "-s", link, "-o", out, "-q"}, 2, "name the same file");
    refused({"--12", a, "-s", a, "-o", out, "-q"}, 2, "name the same file");
    check(!std::filesystem::exists(out + "/contigs.fasta"), "no assembly was written by a refused run");

    // ---- T18 end to end: a gzip input cut at a record boundary is refused
    {
        const std::string text = fastq("", 400);
        const std::string head = fastq("", 250);
        const std::string gz = dir + "/cut.fq.gz";
        gzFile g = gzopen(gz.c_str(), "wb");
        gzwrite(g, text.data(), static_cast<unsigned>(head.size()));
        gzflush(g, Z_SYNC_FLUSH);
        const auto cut = static_cast<uintmax_t>(gzoffset(g));
        gzwrite(g, text.data() + head.size(), static_cast<unsigned>(text.size() - head.size()));
        gzclose(g);
        std::filesystem::resize_file(gz, cut);
        refused({"-s", gz, "-o", dir + "/out_cut", "-q", "-t", "1"}, 1, "unexpected end of file");
    }

    std::filesystem::remove_all(dir);
    if (failures) {
        std::fprintf(stderr, "test_v3_io_cli: %d of %d checks FAILED (binary %s)\n", failures, checks, bin.c_str());
        return 1;
    }
    std::printf("test_v3_io_cli: %d checks passed (binary %s)\n", checks, bin.c_str());
    return 0;
}
