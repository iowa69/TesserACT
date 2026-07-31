#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <thread>

#include "assembler.h"
#include "util.h"

namespace {

constexpr const char* kVersion = "1.0.0";

void usage() {
    std::printf(
        "tessera %s - de novo short-read assembler (de Bruijn, multi-k)\n"
        "\n"
        "USAGE\n"
        "  tessera -1 R1.fq.gz -2 R2.fq.gz -o OUTDIR [options]\n"
        "  tessera --12 interleaved.fq.gz -o OUTDIR\n"
        "  tessera -s single.fq.gz -o OUTDIR\n"
        "\n"
        "INPUT\n"
        "  -1, --read1 FILE        forward reads (FASTQ/FASTA, optionally gzipped)\n"
        "  -2, --read2 FILE        reverse reads\n"
        "      --12 FILE           interleaved paired reads in one file\n"
        "  -s, --single FILE       unpaired reads\n"
        "\n"
        "OUTPUT\n"
        "  -o, --out DIR           output directory (default: tessera_out)\n"
        "      --min-contig N      minimum contig length to report (default: 2*k)\n"
        "\n"
        "ASSEMBLY\n"
        "  -k, --kmers LIST        comma-separated k values, e.g. 21,33,55,77\n"
        "                          (default: chosen from the read length; max 96)\n"
        "  -c, --cutoff N          k-mer abundance cutoff (default: auto-detect)\n"
        "      --min-link N        paired reads needed to trust a join (default: 2)\n"
        "      --tie-ratio F       winning branch must beat the runner-up by F (default: 1.15)\n"
        "      --aggressive        also collapse diverged repeat copies: more\n"
        "                          contiguity, but risks a misassembly per genome\n"
        "      --no-correct        skip read error correction\n"
        "      --no-resolve        skip paired-end repeat resolution\n"
        "      --no-scaffold       skip scaffolding\n"
        "      --no-polish         skip consensus polishing\n"
        "\n"
        "GENERAL\n"
        "  -t, --threads N         worker threads (default: all cores)\n"
        "  -q, --quiet             suppress progress output\n"
        "  -v, --version           print version\n"
        "  -h, --help              print this message\n"
        "\n"
        "EXAMPLES\n"
        "  tessera -1 R1.fq.gz -2 R2.fq.gz -o asm -t 16\n"
        "  tessera -1 R1.fq.gz -2 R2.fq.gz -o asm -k 21,33,55 --min-contig 500\n",
        kVersion);
}


}  // namespace

int main(int argc, char** argv) {
    using namespace ts;

    AssemblyOptions opt;
    Library lib;
    bool haveLib = false;

    if (argc < 2) { usage(); return 1; }

    auto needValue = [&](int& i, const char* flag) -> const char* {
        if (i + 1 >= argc) {
            std::fprintf(stderr, "error: %s requires a value\n", flag);
            std::exit(1);
        }
        return argv[++i];
    };

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "-h" || a == "--help") { usage(); return 0; }
        else if (a == "-v" || a == "--version") { std::printf("tessera %s\n", kVersion); return 0; }
        else if (a == "-1" || a == "--read1") { lib.r1 = needValue(i, "-1"); haveLib = true; }
        else if (a == "-2" || a == "--read2") { lib.r2 = needValue(i, "-2"); haveLib = true; }
        else if (a == "--12") { lib.r1 = needValue(i, "--12"); lib.interleaved = true; haveLib = true; }
        else if (a == "-s" || a == "--single") { lib.r1 = needValue(i, "-s"); haveLib = true; }
        else if (a == "-o" || a == "--out") opt.outDir = needValue(i, "-o");
        else if (a == "--min-contig") opt.minContigLen = static_cast<size_t>(std::atoll(needValue(i, "--min-contig")));
        else if (a == "-k" || a == "--kmers") {
            opt.kValues.clear();
            for (const std::string& tok : util::split(needValue(i, "-k"), ',')) {
                int k = std::atoi(tok.c_str());
                if (k < 5 || k > 96) {
                    std::fprintf(stderr, "error: k must be between 5 and 96 (got %d)\n", k);
                    return 1;
                }
                if (k % 2 == 0) {
                    std::fprintf(stderr, "error: k must be odd to avoid palindromic k-mers (got %d)\n", k);
                    return 1;
                }
                opt.kValues.push_back(k);
            }
        }
        else if (a == "-c" || a == "--cutoff") opt.forcedCutoff = static_cast<uint32_t>(std::atoi(needValue(i, "-c")));
        else if (a == "--min-link") opt.minLinkSupport = std::atoi(needValue(i, "--min-link"));
        else if (a == "--tie-ratio") opt.tieRatio = std::atof(needValue(i, "--tie-ratio"));
        else if (a == "--aggressive") opt.bubbleCoverageLimit = 10.0;
        else if (a == "--no-correct") opt.correctReads = false;
        else if (a == "--no-resolve") opt.resolveRepeats = false;
        else if (a == "--no-scaffold") opt.scaffold = false;
        else if (a == "--no-polish") opt.polish = false;
        else if (a == "-t" || a == "--threads") opt.threads = std::atoi(needValue(i, "-t"));
        else if (a == "-q" || a == "--quiet") opt.verbose = false;
        else {
            std::fprintf(stderr, "error: unknown option '%s'\nRun 'tessera --help' for usage.\n", a.c_str());
            return 1;
        }
    }

    if (!haveLib || lib.r1.empty()) {
        std::fprintf(stderr, "error: no input reads given (use -1/-2, --12, or -s)\n");
        return 1;
    }
    for (const std::string& f : {lib.r1, lib.r2}) {
        if (!f.empty() && !util::fileExists(f)) {
            std::fprintf(stderr, "error: input file not found: %s\n", f.c_str());
            return 1;
        }
    }
    opt.libraries.push_back(lib);

    Assembler asmb(std::move(opt));
    std::string error;
    if (!asmb.run(error)) {
        std::fprintf(stderr, "error: %s\n", error.c_str());
        return 1;
    }

    const AssemblyStats& st = asmb.stats();
    std::fprintf(stderr,
                 "\n"
                 "  contigs      %s\n"
                 "  total length %s bp\n"
                 "  largest      %s bp\n"
                 "  N50          %s bp\n"
                 "  GC           %.2f%%\n"
                 "  mean coverage %.1fx\n"
                 "  elapsed      %.1fs   peak memory %s\n",
                 util::commify(static_cast<long long>(st.contigs)).c_str(),
                 util::commify(static_cast<long long>(st.totalLength)).c_str(),
                 util::commify(static_cast<long long>(st.largest)).c_str(),
                 util::commify(static_cast<long long>(st.n50)).c_str(),
                 st.gcPercent, st.meanCoverage, st.seconds,
                 util::humanBytes(static_cast<double>(util::peakMemoryBytes())).c_str());
    return 0;
}
