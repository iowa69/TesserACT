// Development probe for EVAL_PLAN_P2 W1 R3 (never built by `all`; `make p2probe`).
//
// Runs the assembler's own R3 code (src/p2_emit.cpp: closeWalk, verifyJunctions) on a FINISHED
// assembly directory -- contigs.fasta and assembly_graph.gfa of an F2 (1.5.0) run -- so its graph
// closure and its junction evidence can be compared, circle for circle, with the prototype verdicts
// the plan lets an R3 builder calibrate against (work/p2/panels/jv_verdicts_w1.tsv, EVAL_PLAN_P2
// s7.4). It reads no reference and no truth.
//
// usage: p2_r3_probe <asm_dir> [-1 R1.fq.gz -2 R2.fq.gz] [-t threads]
// prints one TSV row per `_circular` record of contigs.fasta.
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "graph.h"
#include "p2_emit.h"
#include "seqio.h"

namespace {

bool readFasta(const std::string& path, std::vector<std::string>& names, std::vector<std::string>& seqs) {
    std::ifstream in(path);
    if (!in) return false;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            std::string n = line.substr(1);
            const size_t sp = n.find_first_of(" \t");
            if (sp != std::string::npos) n.resize(sp);
            names.push_back(n);
            seqs.emplace_back();
        } else if (!seqs.empty()) {
            for (char c : line) seqs.back() += static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
        }
    }
    return true;
}

std::vector<std::string> split(const std::string& s, char d) {
    std::vector<std::string> out;
    std::string cur;
    for (char c : s) {
        if (c == d) { out.push_back(cur); cur.clear(); }
        else cur += c;
    }
    out.push_back(cur);
    return out;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "usage: p2_r3_probe <asm_dir> [-1 R1 -2 R2] [-t threads]\n");
        return 2;
    }
    const std::string dir = argv[1];
    std::string r1, r2;
    int threads = 4;
    for (int i = 2; i + 1 < argc; i += 2) {
        if (!std::strcmp(argv[i], "-1")) r1 = argv[i + 1];
        else if (!std::strcmp(argv[i], "-2")) r2 = argv[i + 1];
        else if (!std::strcmp(argv[i], "-t")) threads = std::atoi(argv[i + 1]);
    }
    // ---- the graph as the GFA writes it ----------------------------------------------------
    std::ifstream gfa(dir + "/assembly_graph.gfa");
    if (!gfa) { std::fprintf(stderr, "no GFA in %s\n", dir.c_str()); return 1; }
    ts::UnitigGraph g;
    std::map<std::string, std::vector<uint64_t>> paths;
    struct L { uint32_t u; int eu; uint32_t v; int ev; };
    std::vector<L> links;
    int k = 0;
    std::string line;
    while (std::getline(gfa, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> f = split(line, '\t');
        if (f[0] == "S" && f.size() >= 3) {
            const uint32_t u = static_cast<uint32_t>(std::stoul(f[1]));
            if (g.nodes.size() <= u) {
                const size_t old = g.nodes.size();
                g.nodes.resize(u + 1);
                for (size_t x = old; x <= u; ++x) g.nodes[x].deleted = true;
            }
            g.nodes[u].seq = f[2];
            g.nodes[u].deleted = false;
            for (size_t t = 3; t < f.size(); ++t)
                if (f[t].rfind("dp:f:", 0) == 0) g.nodes[u].coverage = std::strtod(f[t].c_str() + 5, nullptr);
        } else if (f[0] == "L" && f.size() >= 6) {
            const uint32_t u = static_cast<uint32_t>(std::stoul(f[1])), v = static_cast<uint32_t>(std::stoul(f[3]));
            links.push_back({u, f[2] == "+" ? 1 : 0, v, f[4] == "+" ? 0 : 1});
            if (!k) k = std::atoi(f[5].c_str()) + 1;
        } else if (f[0] == "P" && f.size() >= 3) {
            std::vector<uint64_t> w;
            for (const std::string& step : split(f[2], ',')) {
                if (step.size() < 2) continue;
                const uint64_t u = std::stoull(step.substr(0, step.size() - 1));
                w.push_back((u << 1) | (step.back() == '-' ? 1ULL : 0ULL));
            }
            paths[f[1]] = w;
        }
    }
    for (const L& l : links) {
        if (l.u >= g.nodes.size() || l.v >= g.nodes.size()) continue;
        g.nodes[l.u].ends[l.eu].push_back({l.v, static_cast<uint8_t>(l.ev)});
        if (!(l.u == l.v && l.eu == l.ev)) g.nodes[l.v].ends[l.ev].push_back({l.u, static_cast<uint8_t>(l.eu)});
    }
    g.setK(k);

    std::vector<std::string> names, seqs;
    if (!readFasta(dir + "/contigs.fasta", names, seqs)) { std::fprintf(stderr, "no contigs.fasta\n"); return 1; }

    // ---- closure per circular record -----------------------------------------------------
    std::vector<size_t> cand;
    std::vector<ts::p2::ClosureResult> cl;
    ts::p2::VerifyBatch b;
    b.records = &seqs;
    std::vector<size_t> selfOv;
    for (size_t i = 0; i < names.size(); ++i) {
        const std::string& n = names[i];
        if (n.size() < 9 || n.compare(n.size() - 9, 9, "_circular") != 0) continue;
        const auto it = paths.find(n);
        ts::p2::ClosureResult r = it == paths.end() ? ts::p2::ClosureResult()
                                                    : ts::p2::closeWalk(seqs[i], it->second, g);
        std::string T;
        size_t ov = 0;
        if (r.closed) T = r.seq;
        else {
            ov = seqs[i].size() >= 40 ? ts::p2::terminalSelfOverlap(seqs[i]) : 0;
            T = seqs[i].substr(0, seqs[i].size() - ov);
        }
        cand.push_back(i);
        cl.push_back(r);
        selfOv.push_back(ov);
        const size_t h = T.size() / 2;
        b.candRecord.push_back(i);
        b.candCircle.push_back(T.substr(h) + T.substr(0, h));
    }
    std::vector<ts::p2::JunctionStats> js(cand.size());
    ts::p2::VerifyRunStats run;
    if (!r1.empty() && !cand.empty()) {
        ts::SequenceStore reads;
        ts::Library lib;
        lib.r1 = r1;
        lib.r2 = r2;
        std::string err;
        if (!reads.load({lib}, threads, err)) { std::fprintf(stderr, "reads: %s\n", err.c_str()); return 1; }
        js = ts::p2::verifyJunctions(b, reads, threads, run);
    }
    const ts::p2::VerifyThresholds th;
    std::printf("# k=%d window=%d-%d window_pairs=%zu reads_placed=%zu index_kmers=%zu\n", k, run.windowLo, run.windowHi,
                run.windowPairs, run.readsPlaced, run.indexKmers);
    std::printf("contig\tg_mode\tg_isolated\tclosed_len\tself_ov\tcircle_len\tspan_reads\tclip_reads\tclip_frac\tctrl_span\t"
                "near_join\texits_near_join\texit_ratio\tJ\tJX02\n");
    for (size_t c = 0; c < cand.size(); ++c) {
        const ts::p2::ClosureResult& r = cl[c];
        const ts::p2::JunctionStats& j = js[c];
        const bool J = ts::p2::passJ(j, th);
        const bool X = ts::p2::passX(j, th);
        std::printf("%s\t%s\t%s\t%s\t%zu\t%zu\t%zu\t%zu\t%.4f\t%zu\t%zu\t%zu\t%.4f\t%s\t%s\n", names[cand[c]].c_str(),
                    ts::p2::closureName(r.mode), r.isolated >= 0 ? std::to_string(r.isolated).c_str() : "",
                    r.closed ? std::to_string(r.seq.size()).c_str() : "", selfOv[c], j.circleLen, j.spanReads,
                    j.clipReads, j.clipFrac, j.ctrlSpan, j.nearJoin, j.exitsNearJoin, j.exitRatio,
                    j.computed ? (J ? "1" : "0") : "", j.computed ? (J && X ? "1" : "0") : "");
    }
    return 0;
}
