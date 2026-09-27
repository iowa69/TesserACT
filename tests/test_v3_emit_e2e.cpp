// G-emit end to end (Assembler::run on a small synthetic isolate with a repeat at two loci and a
// coverage dropout that leaves a scaffold gap):
//   T34  a re-run into a directory holding another run's optional outputs removes them
//   T17  report.json contig_* / scaffold_gaps describe contigs.fasta / scaffolds.fasta, and
//        report.json carries the gap-filler statistics
//   T26  with TESSERACT_FIX_AGP_GFA_V2=1 every AGP W row names a contigs.fasta record and
//        matches it base for base, and every P-line spells the record it names
//   A2   every G-emit counter line is printed, flags off and on
// Against the release sources the T34 and T17 assertions fail (stale files survive; the report
// counts runs of 10+ N in untrimmed scaffolds and has no gap_fill object).
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <map>
#include <random>
#include <regex>
#include <sstream>
#include <string>
#include <vector>
#include <unistd.h>

#include "assembler.h"
#include "graph.h"
#include "test_env.h"   // build_v3: T40, clear the ambient TESSERACT_* first

namespace fs = std::filesystem;

namespace {
int checks = 0, failures = 0;
void check(bool ok, const std::string& what) {
    ++checks;
    std::printf("%s %s\n", ok ? "ok  " : "FAIL", what.c_str());
    if (!ok) ++failures;
}
std::string rc(const std::string& s) { return ts::reverseComplement(s); }

std::map<std::string, std::string> readFasta(const std::string& path) {
    std::map<std::string, std::string> out;
    std::ifstream f(path);
    std::string line, name;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') { name = line.substr(1); out[name]; }
        else out[name] += line;
    }
    return out;
}
std::string slurp(const std::string& p) {
    std::ifstream f(p);
    std::stringstream ss;
    ss << f.rdbuf();
    return ss.str();
}
long long jsonInt(const std::string& json, const std::string& key) {
    std::smatch m;
    if (std::regex_search(json, m, std::regex("\"" + key + "\": *([0-9]+)"))) return std::stoll(m[1]);
    return -1;
}
size_t nRuns(const std::string& s) {
    size_t n = 0;
    for (size_t i = 0; i < s.size();) {
        if (s[i] != 'N' && s[i] != 'n') { ++i; continue; }
        ++n;
        while (i < s.size() && (s[i] == 'N' || s[i] == 'n')) ++i;
    }
    return n;
}

struct Data { std::string dir, r1, r2; };
Data makeData() {
    std::mt19937 rng(7);
    auto dna = [&](size_t n) { std::string s(n, 'A'); for (char& c : s) c = "ACGT"[rng() % 4]; return s; };
    const std::string rep = dna(700);
    const std::string g = dna(14000) + rep + dna(14000) + rep + dna(14000);
    Data d;
    d.dir = (fs::temp_directory_path() / ("v3emit-e2e-" + std::to_string(getpid()))).string();
    fs::create_directories(d.dir);
    d.r1 = d.dir + "/r1.fa";
    d.r2 = d.dir + "/r2.fa";
    std::ofstream a(d.r1), b(d.r2);
    size_t n = 0;
    // No read touches [21970, 22030), but fragments of 300-500 bp span it: the pairs join the
    // two sides across a gap nothing can fill.
    auto touches = [](size_t b, size_t e) { return b < 22030 && e > 21970; };
    for (int i = 0; i < 9000; ++i) {
        const size_t len = 300 + rng() % 200;
        const size_t p = rng() % (g.size() - len);
        if (touches(p, p + 150) || touches(p + len - 150, p + len)) continue;
        const std::string frag = g.substr(p, len);
        a << ">p" << n << "/1\n" << frag.substr(0, 150) << "\n";
        b << ">p" << n << "/2\n" << rc(frag.substr(len - 150)) << "\n";
        ++n;
    }
    return d;
}

bool assemble(const Data& d, const std::string& out, bool gfa, bool html) {
    ts::AssemblyOptions opt;
    ts::Library lib;
    lib.r1 = d.r1;
    lib.r2 = d.r2;
    opt.libraries.push_back(lib);
    opt.outDir = out;
    opt.threads = 2;
    opt.verbose = false;
    opt.emitGfa = gfa;
    opt.emitHtml = html;
    ts::Assembler asmb(std::move(opt));
    std::string error;
    const bool ok = asmb.run(error);
    if (!ok) std::fprintf(stderr, "assembly failed: %s\n", error.c_str());
    return ok;
}

void clearEnv() {
    for (const char* v : {"TESSERACT_FIXES", "TESSERACT_FIX_POLISH_SKIP_N", "TESSERACT_FIX_GAPFILL_STRICT_BUDGET",
                          "TESSERACT_FIX_GAPFILL_BACKOFF", "TESSERACT_FIX_GAPFILL_SKIP_INPUT_N",
                          "TESSERACT_FIX_BOUNDARY_SAFE_TRIM", "TESSERACT_FIX_TRIM_COPY_GUARD",
                          "TESSERACT_FIX_SPLIT_POSTPROCESS", "TESSERACT_FIX_AGP_GFA_V2"})
        unsetenv(v);
}

// stderr of the process is captured into a file for the counter-line checks.
struct StderrToFile {
    int saved = -1;
    explicit StderrToFile(const std::string& path) {
        std::fflush(stderr);
        saved = dup(2);
        std::FILE* f = std::freopen(path.c_str(), "w", stderr);
        (void)f;
    }
    ~StderrToFile() {
        std::fflush(stderr);
        dup2(saved, 2);
        close(saved);
        clearerr(stderr);
    }
};

void checkReport(const std::string& out, const char* label) {
    const auto contigs = readFasta(out + "/contigs.fasta");
    size_t total = 0;
    for (const auto& kv : contigs) total += kv.second.size();
    size_t gaps = 0;
    if (fs::exists(out + "/scaffolds.fasta"))
        for (const auto& kv : readFasta(out + "/scaffolds.fasta")) gaps += nRuns(kv.second);
    const std::string json = slurp(out + "/report.json");
    check(jsonInt(json, "contig_count") == static_cast<long long>(contigs.size()) &&
              jsonInt(json, "contig_total_length") == static_cast<long long>(total) &&
              jsonInt(json, "scaffold_gaps") == static_cast<long long>(gaps),
          std::string("T17 ") + label + ": report.json contig_count/total/scaffold_gaps match the written files (" +
              std::to_string(contigs.size()) + " records, " + std::to_string(gaps) + " gaps)");
    check(json.find("\"gap_fill\"") != std::string::npos && jsonInt(json, "gaps_seen") >= 0,
          std::string("T17 ") + label + ": report.json carries the gap_fill object");
}
}  // namespace

int main() {
    testenv::clearTesseractEnv();
    clearEnv();
    const Data d = makeData();
    const std::string outA = d.dir + "/outA";
    const std::string markers[] = {"scaffolds.fasta", "scaffolds.agp", "unitigs.fasta", "assembly_graph.gfa",
                                   "report.html"};

    // ---- T34 + T17 + counters, switches off ---------------------------------------------------
    fs::create_directories(outA);
    for (const std::string& m : markers) std::ofstream(outA + "/" + m) << ">STALE\nACGT\n";
    bool ok;
    {
        StderrToFile cap(d.dir + "/a.err");
        ok = assemble(d, outA, /*gfa=*/false, /*html=*/false);
    }
    check(ok, "assembly A ran");
    size_t stale = 0;
    for (const std::string& m : markers)
        if (fs::exists(outA + "/" + m) && slurp(outA + "/" + m).find("STALE") != std::string::npos) ++stale;
    check(stale == 0, "T34: no optional output from the earlier run survives (" + std::to_string(stale) + " stale)");
    check(!fs::exists(outA + "/unitigs.fasta") && !fs::exists(outA + "/assembly_graph.gfa") &&
              !fs::exists(outA + "/report.html"),
          "T34: outputs this run did not write are absent");
    check(fs::exists(outA + "/scaffolds.fasta"), "fixture: the dropout left a scaffold gap");
    checkReport(outA, "off");
    const std::string errA = slurp(d.dir + "/a.err");
    for (const char* line : {"[fixes-emit] umbrella=0", "[polish-n] skip_n=0", "[gapfill-budget] strict=0",
                             "[gapfill-backoff] B=0", "[gapfill-inputn] enabled=0", "[dedup] ",
                             "[outputs] stale_removed=5", "[splitpost] enabled=0", "[boundarytrim] enabled=0",
                             "[trimguard] enabled=0", "[agpgfa] enabled=0"})
        check(errA.find(line) != std::string::npos, std::string("A2 off: counter line '") + line + "' present");

    // ---- every switch on (umbrella) -------------------------------------------------------------
    setenv("TESSERACT_FIXES", "1", 1);
    const std::string outB = d.dir + "/outB";
    {
        StderrToFile cap(d.dir + "/b.err");
        ok = assemble(d, outB, true, false);
    }
    check(ok, "assembly B (TESSERACT_FIXES=1) ran");
    checkReport(outB, "on");
    const std::string errB = slurp(d.dir + "/b.err");
    for (const char* line : {"[fixes-emit] umbrella=1 polish_skip_n=1 gapfill_strict_budget=1 gapfill_backoff=60 "
                             "gapfill_skip_input_n=1 boundary_safe_trim=1 trim_copy_guard=1 split_postprocess=1 "
                             "agp_gfa_v2=1",
                             "[polish-n] skip_n=1", "[gapfill-budget] strict=1", "[gapfill-backoff] B=60",
                             "[gapfill-inputn] enabled=1 provenance=1", "[splitpost] enabled=1",
                             "[boundarytrim] enabled=1", "[trimguard] enabled=1", "[agpgfa] enabled=1"})
        check(errB.find(line) != std::string::npos, std::string("A2 on: counter line '") + line + "' present");
    // T26: the AGP describes contigs.fasta.
    {
        const auto contigs = readFasta(outB + "/contigs.fasta");
        const auto scaf = readFasta(outB + "/scaffolds.fasta");
        std::ifstream agp(outB + "/scaffolds.agp");
        std::string line;
        size_t w = 0, wIn = 0, bad = 0, nRows = 0, alignGenus = 0;
        std::map<std::string, size_t> next;
        while (std::getline(agp, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::istringstream is(line);
            std::string o, t, c;
            size_t b, e, part;
            is >> o >> b >> e >> part >> t;
            const std::string& s = scaf.at(o);
            if (b != next[o] + 1) ++bad;
            next[o] = e;
            const std::string region = s.substr(b - 1, e - b + 1);
            if (t == "N") {
                ++nRows;
                if (line.find("align_genus") != std::string::npos) ++alignGenus;
                if (region.find_first_not_of("Nn") != std::string::npos) ++bad;
                continue;
            }
            ++w;
            size_t cb, ce;
            char orient;
            is >> c >> cb >> ce >> orient;
            auto it = contigs.find(c);
            if (it == contigs.end()) continue;
            ++wIn;
            std::string piece = it->second.substr(cb - 1, ce - cb + 1);
            if (orient == '-') piece = rc(piece);
            if (piece != region) ++bad;
        }
        for (const auto& kv : scaf) if (next[kv.first] != kv.second.size()) ++bad;
        check(w > 0 && wIn == w && bad == 0 && nRows > 0 && alignGenus == 0,
              "T26: every AGP W row names a contigs.fasta record and matches it; tiling complete; no "
              "align_genus without a model");
        // P-lines spell the record they name.
        ts::UnitigGraph dummy;
        std::map<std::string, std::string> seg;
        size_t pLines = 0, pBad = 0;
        int k = 0;
        std::ifstream gfa(outB + "/assembly_graph.gfa");
        std::vector<std::string> lines;
        while (std::getline(gfa, line)) lines.push_back(line);
        for (const std::string& l : lines) {
            if (l.rfind("S\t", 0) == 0) {
                std::istringstream is(l);
                std::string t, id, s;
                is >> t >> id >> s;
                seg[id] = s;
            } else if (l.rfind("L\t", 0) == 0 && !k) {
                k = std::stoi(l.substr(l.rfind('\t') + 1)) + 1;
            }
        }
        for (const std::string& l : lines) {
            if (l.rfind("P\t", 0) != 0) continue;
            ++pLines;
            std::istringstream is(l);
            std::string t, name, walk;
            is >> t >> name >> walk;
            std::string spelled;
            std::stringstream ws(walk);
            std::string step;
            while (std::getline(ws, step, ',')) {
                const std::string id = step.substr(0, step.size() - 1);
                std::string piece = seg[id];
                if (step.back() == '-') piece = rc(piece);
                spelled += spelled.empty() ? piece : piece.substr(static_cast<size_t>(k - 1));
            }
            auto it = contigs.find(name);
            if (it == contigs.end() || it->second != spelled) ++pBad;
        }
        check(pBad == 0, "T26: every P-line (" + std::to_string(pLines) + ") names a contigs.fasta record and spells it");
    }
    clearEnv();
    fs::remove_all(d.dir);
    std::printf("test_v3_emit_e2e: %d checks, %d failures\n", checks, failures);
    return failures ? 1 : 0;
}
