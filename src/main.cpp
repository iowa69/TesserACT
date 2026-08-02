#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <thread>

#include "assembler.h"
#include "util.h"

namespace {

constexpr const char* kVersion = "1.1.0";

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
        "ORGANISM MODEL\n"
        "      --organism NAME     organism the reads come from (e.g. klebsiella)\n"
        "      --model FILE        genus model built by tessera-model. Used only at\n"
        "                          junctions no fragment can span: the chromosome is\n"
        "                          reconstructed first from conserved gene order, then\n"
        "                          plasmids are refined against the plasmid table.\n"
        "                          Without it those junctions are left broken rather\n"
        "                          than guessed at.\n"
        "\n"
        "      --is-panel FILE     FASTA of known insertion sequences. Contig ends\n"
        "                          lying inside one are left unjoined: that is where\n"
        "                          this isolate differs from the panel by construction.\n"
        "      --is-sites FILE     table of recurrent insertion sites (flank\n"
        "                          signatures); lets an element be placed on the\n"
        "                          chromosome instead of leaving the end unjoined\n"
        "      --map-polish NAME   polish the finished contigs against a full read\n"
        "                          alignment: bowtie2 | bwa | none (default none)\n"
        "      --mapper-dir DIR    where to find the mapper binaries\n"
        "\n"
        "MODES\n"
        "      --mode NAME         fast | standard (default) | careful | aggressive\n"
        "                          fast       fewer k values, no polishing\n"
        "                          standard   balanced; what the benchmarks use\n"
        "                          careful    denser k ladder, stricter joins, 2 polish passes\n"
        "                          aggressive collapse diverged repeats for maximum contiguity\n"
        "\n"
        "OUTPUT FILES\n"
        "      --no-gfa            skip assembly_graph.gfa\n"
        "      --no-html           skip report.html (report.json is always written)\n"
        "      --unitigs           also write unitigs.fasta (pre-resolution graph)\n"
        "\n"
        "ASSEMBLY\n"
        "  -k, --kmers LIST        comma-separated k values, e.g. 21,33,55,77\n"
        "                          (default: chosen from the read length; max 96)\n"
        "  -c, --cutoff N          k-mer abundance cutoff (default: auto-detect)\n"
        "      --min-link N        paired reads needed to trust a join (default: 2)\n"
        "      --tie-ratio F       winning branch must beat the runner-up by F (default: 1.15)\n"
        "      --link-per-x F      paired support required per unit of median coverage\n"
        "                          before a contested join is taken (default: 0.10)\n"
        "      --aggressive        alias for --mode aggressive\n"
        "      --simplify-rounds N max graph simplification passes per k\n"
        "      --polish-passes N   consensus polishing passes (0 disables)\n"
        "      --no-correct        skip read error correction\n"
        "      --no-qtrim          keep 3' ends that fall below the quality cutoff\n"
        "      --qtrim-quality N   3' trim below mean quality N (default: 20)\n"
        "      --qtrim-window N    window used for the 3' quality scan (default: 4)\n"
        "      --no-resolve        skip paired-end repeat resolution\n"
        "      --no-scaffold       skip scaffolding\n"
        "      --no-gapfill        skip scaffold gap closing\n"
        "      --no-polish         skip consensus polishing\n"
        "\n"
        "GENERAL\n"
        "  -t, --threads N         worker threads (default: all cores)\n"
        "      --max-memory GB     abort cleanly above this much RAM (default: 80%% of total)\n"
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
    std::vector<Library> singles;
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
        // A separate library, not an overwrite of -1. Merging overlapping
        // mates yields a set of single-end fragments plus the pairs that did
        // not merge, and both belong in the same assembly; -s may be repeated.
        else if (a == "-s" || a == "--single") {
            Library sl;
            sl.r1 = needValue(i, "-s");
            singles.push_back(std::move(sl));
            haveLib = true;
        }
        else if (a == "-o" || a == "--out") opt.outDir = needValue(i, "-o");
        else if (a == "--organism") opt.organism = needValue(i, "--organism");
        else if (a == "--model") opt.organismModelPath = needValue(i, "--model");
        else if (a == "--mapper-dir") opt.mapperDir = needValue(i, "--mapper-dir");
        else if (a == "--is-panel") opt.isPanelPath = needValue(i, "--is-panel");
        else if (a == "--is-sites") opt.isSitesPath = needValue(i, "--is-sites");
        else if (a == "--map-polish") {
            const std::string m = needValue(i, "--map-polish");
            if (!parseMapper(m, opt.mapPolisher)) {
                std::fprintf(stderr, "error: unknown mapper '%s' (bowtie2 | bwa | none)\n", m.c_str());
                return 2;
            }
        }
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
            opt.userSetK = true;
        }
        else if (a == "-c" || a == "--cutoff") opt.forcedCutoff = static_cast<uint32_t>(std::atoi(needValue(i, "-c")));
        else if (a == "--min-link") { opt.minLinkSupport = std::atoi(needValue(i, "--min-link")); opt.userSetMinLink = true; }
        else if (a == "--link-per-x") {
            opt.linkSupportPerX = std::atof(needValue(i, "--link-per-x"));
            opt.userSetLinkPerX = true;
        }
        else if (a == "--tie-ratio") { opt.tieRatio = std::atof(needValue(i, "--tie-ratio")); opt.userSetTie = true; }
        else if (a == "--bubble-coverage") { opt.bubbleCoverageLimit = std::atof(needValue(i, "--bubble-coverage")); opt.userSetBubble = true; }
        else if (a == "--simplify-rounds") { opt.simplifyRounds = std::atoi(needValue(i, "--simplify-rounds")); opt.userSetRounds = true; }
        else if (a == "--polish-passes") { opt.polishPasses = std::atoi(needValue(i, "--polish-passes")); opt.userSetPolishPasses = true; }
        else if (a == "--aggressive") opt.mode = RunMode::Aggressive;
        else if (a == "--mode") {
            const std::string m = needValue(i, "--mode");
            if (!parseRunMode(m, opt.mode)) {
                std::fprintf(stderr,
                             "error: unknown mode '%s' (expected fast, standard, careful or aggressive)\n",
                             m.c_str());
                return 1;
            }
        }
        else if (a == "--max-memory") {
            opt.maxMemoryBytes = static_cast<long long>(std::atof(needValue(i, "--max-memory")) * 1073741824.0);
        }
        else if (a == "--no-gfa") opt.emitGfa = false;
        else if (a == "--no-html") opt.emitHtml = false;
        else if (a == "--unitigs") opt.emitUnitigs = true;
        else if (a == "--no-correct") opt.correctReads = false;
        else if (a == "--no-qtrim") opt.qtrim.enabled = false;
        else if (a == "--qtrim-quality") opt.qtrim.meanQuality = std::atoi(needValue(i, "--qtrim-quality"));
        else if (a == "--qtrim-window") opt.qtrim.windowSize = std::atoi(needValue(i, "--qtrim-window"));
        else if (a == "--no-resolve") opt.resolveRepeats = false;
        else if (a == "--no-scaffold") opt.scaffold = false;
        else if (a == "--no-gapfill") opt.gapFill = false;
        else if (a == "--no-polish") opt.polish = false;
        else if (a == "-t" || a == "--threads") opt.threads = std::atoi(needValue(i, "-t"));
        else if (a == "-q" || a == "--quiet") opt.verbose = false;
        else {
            std::fprintf(stderr, "error: unknown option '%s'\nRun 'tessera --help' for usage.\n", a.c_str());
            return 1;
        }
    }

    if (!haveLib || (lib.r1.empty() && singles.empty())) {
        std::fprintf(stderr, "error: no input reads given (use -1/-2, --12, or -s)\n");
        return 1;
    }
    if (lib.r1.empty() && !lib.r2.empty()) {
        std::fprintf(stderr, "error: -2 given without -1\n");
        return 1;
    }
    if (!lib.r1.empty()) opt.libraries.push_back(lib);
    for (Library& sl : singles) opt.libraries.push_back(std::move(sl));
    for (const Library& l : opt.libraries) {
        for (const std::string& f : {l.r1, l.r2}) {
            if (!f.empty() && !util::fileExists(f)) {
                std::fprintf(stderr, "error: input file not found: %s\n", f.c_str());
                return 1;
            }
        }
    }

    Assembler asmb(std::move(opt));
    std::string error;
    if (!asmb.run(error)) {
        std::fprintf(stderr, "error: %s\n", error.c_str());
        return 1;
    }

    const AssemblyStats& st = asmb.stats();
    const AssemblyReport& rep = asmb.report();
    const int finalK = rep.iterations.empty() ? 0 : rep.iterations.back().k;

    std::fprintf(stderr,
                 "\n"
                 "  contigs      %s\n"
                 "  total length %s bp\n"
                 "  largest      %s bp\n"
                 "  N50          %s bp    N90 %s bp    L50 %s\n"
                 "  GC           %.2f%%\n",
                 util::commify(static_cast<long long>(st.contigs)).c_str(),
                 util::commify(static_cast<long long>(st.totalLength)).c_str(),
                 util::commify(static_cast<long long>(st.largest)).c_str(),
                 util::commify(static_cast<long long>(st.n50)).c_str(),
                 util::commify(static_cast<long long>(rep.n90)).c_str(),
                 util::commify(static_cast<long long>(rep.l50)).c_str(),
                 st.gcPercent);

    // k-mer depth and read depth differ by a factor of L/(L-k+1); reporting one
    // as "coverage" invites reading it as the other.
    std::fprintf(stderr, "  k-mer depth  %.1fx (k=%d)\n", st.meanCoverage, finalK);
    if (rep.polishRun && rep.polish.meanDepth > 0) {
        std::fprintf(stderr, "  read depth   %.1fx (measured by mapping reads back)\n",
                     rep.polish.meanDepth);
    }
    if (rep.gapBases) {
        std::fprintf(stderr, "  scaffold gaps %s joins spanning %s N bases\n",
                     util::commify(static_cast<long long>(rep.resolve.scaffoldJoins)).c_str(),
                     util::commify(static_cast<long long>(rep.gapBases)).c_str());
    }
    std::fprintf(stderr, "  elapsed      %.1fs   peak memory %s\n",
                 st.seconds,
                 util::humanBytes(static_cast<double>(util::peakMemoryBytes())).c_str());
    return 0;
}
