// Machine-readable run report. Written with a minimal serialiser rather than a
// JSON library so the project keeps its zero-dependency build.
#include <cstdio>
#include <string>
#include <vector>

#include "report.h"

namespace ts {

namespace {

void esc(const std::string& in, std::string& out) {
    for (char c : in) {
        switch (c) {
            case '"':  out += "\\\""; break;
            case '\\': out += "\\\\"; break;
            case '\n': out += "\\n";  break;
            case '\r': out += "\\r";  break;
            case '\t': out += "\\t";  break;
            default:
                if (static_cast<unsigned char>(c) < 0x20) {
                    char buf[8];
                    std::snprintf(buf, sizeof(buf), "\\u%04x", c);
                    out += buf;
                } else {
                    out += c;
                }
        }
    }
}

struct Writer {
    std::string s;
    int indent = 0;

    void pad() { s.append(static_cast<size_t>(indent) * 4, ' '); }
    void key(const char* k) { pad(); s += '"'; s += k; s += "\": "; }

    void str(const char* k, const std::string& v) {
        key(k); s += '"'; esc(v, s); s += "\",\n";
    }
    void num(const char* k, double v) {
        key(k);
        char buf[64];
        std::snprintf(buf, sizeof(buf), "%.6g", v);
        s += buf; s += ",\n";
    }
    void uint(const char* k, unsigned long long v) {
        key(k); s += std::to_string(v); s += ",\n";
    }
    void boolean(const char* k, bool v) { key(k); s += v ? "true" : "false"; s += ",\n"; }

    void openObj(const char* k) { key(k); s += "{\n"; ++indent; }
    void openArr(const char* k) { key(k); s += "[\n"; ++indent; }
    void close(char c) {
        // Drop the trailing comma the members left behind.
        if (s.size() >= 2 && s[s.size() - 2] == ',') s.erase(s.size() - 2, 1);
        --indent;
        pad();
        s += c;
        s += ",\n";
    }

    template <typename T>
    void arrOfNumbers(const char* k, const std::vector<T>& v, size_t cap = 0) {
        key(k); s += '[';
        const size_t n = (cap && v.size() > cap) ? cap : v.size();
        for (size_t i = 0; i < n; ++i) {
            if (i) s += ',';
            s += std::to_string(v[i]);
        }
        s += "],\n";
    }
};

}  // namespace

bool writeJsonReport(const std::string& path, const AssemblyReport& rep, std::string& error) {
    Writer w;
    w.s += "{\n";
    w.indent = 1;

    w.str("tool", "tessera");
    w.str("version", rep.version);
    w.str("mode", rep.mode);
    w.str("started_at", rep.startedAt);
    w.str("command", rep.command);
    w.uint("threads", static_cast<unsigned long long>(rep.threads));
    w.num("total_seconds", rep.totalSeconds);
    w.num("peak_memory_bytes", rep.peakMemoryBytes);

    w.openObj("input");
    w.uint("reads", rep.reads);
    w.uint("bases", rep.inputBases);
    w.uint("max_read_length", rep.maxReadLength);
    w.uint("quality_trimmed_bases", rep.qualityTrimmedBases);
    w.boolean("paired", rep.paired);
    w.key("files"); w.s += "[";
    for (size_t i = 0; i < rep.inputFiles.size(); ++i) {
        if (i) w.s += ',';
        w.s += '"'; esc(rep.inputFiles[i], w.s); w.s += '"';
    }
    w.s += "],\n";
    w.close('}');

    w.openObj("error_correction");
    w.boolean("run", rep.correctionRun);
    w.uint("k", static_cast<unsigned long long>(rep.correctionK));
    w.uint("reads_examined", rep.correction.readsExamined);
    w.uint("reads_corrected", rep.correction.readsCorrected);
    w.uint("bases_corrected", rep.correction.basesCorrected);
    w.uint("reads_without_anchor", rep.correction.readsUncorrectable);
    w.uint("bases_masked", rep.correction.basesMasked);
    w.num("seconds", rep.correctionSeconds);
    w.close('}');

    w.openArr("iterations");
    for (const KIteration& it : rep.iterations) {
        w.pad(); w.s += "{\n"; ++w.indent;
        w.uint("k", static_cast<unsigned long long>(it.k));
        w.uint("total_kmers", it.totalKmers);
        w.uint("distinct_kmers", it.distinctKmers);
        w.uint("solid_kmers", it.solidKmers);
        w.uint("cutoff", it.cutoff);
        w.num("peak_coverage", it.peakCoverage);
        w.num("median_coverage", it.medianCoverage);
        w.uint("carry_over_contigs", it.carryOverContigs);
        w.uint("unitigs_built", it.unitigsBuilt);
        w.uint("length_built", it.lengthBuilt);
        w.uint("n50_built", it.n50Built);
        w.uint("unitigs_final", it.unitigsFinal);
        w.uint("length_final", it.lengthFinal);
        w.uint("n50_final", it.n50Final);
        w.num("count_seconds", it.countSeconds);
        w.num("graph_seconds", it.graphSeconds);
        w.num("simplify_seconds", it.simplifySeconds);
        w.arrOfNumbers("count_histogram", it.countHistogram);
        w.openArr("simplify_rounds");
        for (const SimplifyRoundStats& r : it.rounds) {
            w.pad(); w.s += "{";
            w.s += "\"round\":" + std::to_string(r.round);
            w.s += ",\"tips\":" + std::to_string(r.tipsRemoved);
            w.s += ",\"bubbles\":" + std::to_string(r.bubblesPopped);
            w.s += ",\"chimeras\":" + std::to_string(r.chimerasRemoved);
            w.s += ",\"isolated\":" + std::to_string(r.isolatedRemoved);
            w.s += ",\"merged\":" + std::to_string(r.merged);
            w.s += ",\"unitigs\":" + std::to_string(r.unitigs);
            w.s += ",\"n50\":" + std::to_string(r.n50);
            w.s += ",\"total_length\":" + std::to_string(r.totalLength);
            w.s += "},\n";
        }
        w.close(']');
        w.close('}');
    }
    w.close(']');

    w.openObj("repeat_resolution");
    w.boolean("run", rep.resolveRun);
    w.uint("reads_anchored", rep.resolve.readsMapped);
    w.uint("linking_pairs", rep.resolve.pairsLinking);
    w.uint("distinct_links", rep.resolve.distinctLinks);
    w.uint("paths_built", rep.resolve.pathsBuilt);
    w.uint("unitigs_joined", rep.resolve.unitigsJoined);
    w.uint("scaffold_joins", rep.resolve.scaffoldJoins);
    w.uint("gap_bases", rep.resolve.gapBases);
    w.num("seconds", rep.resolveSeconds);
    w.openObj("insert_size");
    w.boolean("usable", rep.resolve.insert.usable);
    w.num("mean", rep.resolve.insert.mean);
    w.num("stddev", rep.resolve.insert.stddev);
    w.uint("observations", rep.resolve.insert.observations);
    w.uint("min_plausible", static_cast<unsigned long long>(rep.resolve.insert.minPlausible));
    w.uint("max_plausible", static_cast<unsigned long long>(rep.resolve.insert.maxPlausible));
    w.close('}');
    w.arrOfNumbers("insert_histogram", rep.insertHistogram, 2000);
    w.close('}');

    w.openObj("polishing");
    w.boolean("run", rep.polishRun);
    w.uint("reads_used", rep.polish.readsUsed);
    w.uint("bases_changed", rep.polish.basesChanged);
    w.uint("positions_covered", rep.polish.positionsCovered);
    w.uint("low_coverage_positions", rep.polish.lowCoveragePositions);
    w.num("mean_depth", rep.polish.meanDepth);
    w.num("seconds", rep.polishSeconds);
    w.close('}');

    w.openObj("assembly");
    w.uint("contigs", rep.contigs.size());
    w.uint("total_length", rep.totalLength);
    w.uint("largest", rep.largest);
    w.uint("n50", rep.n50);
    w.uint("n90", rep.n90);
    w.uint("l50", rep.l50);
    w.num("gc_percent", rep.gcPercent);
    w.num("mean_coverage", rep.meanCoverage);
    w.uint("gap_bases", rep.gapBases);
    w.uint("gfa_segments", rep.gfaSegments);
    w.uint("gfa_links", rep.gfaLinks);
    w.openArr("contig_table");
    for (const ContigRecord& c : rep.contigs) {
        w.pad(); w.s += "{";
        w.s += "\"length\":" + std::to_string(c.length);
        char buf[64];
        std::snprintf(buf, sizeof(buf), "%.4f", c.coverage);
        w.s += ",\"coverage\":"; w.s += buf;
        std::snprintf(buf, sizeof(buf), "%.2f", c.gcPercent);
        w.s += ",\"gc\":"; w.s += buf;
        w.s += ",\"gap_bases\":" + std::to_string(c.gapBases);
        w.s += "},\n";
    }
    w.close(']');
    w.close('}');

    if (w.s.size() >= 2 && w.s[w.s.size() - 2] == ',') w.s.erase(w.s.size() - 2, 1);
    w.s += "}\n";

    std::FILE* f = std::fopen(path.c_str(), "wb");
    if (!f) { error = "cannot write " + path; return false; }
    // The final flush happens in fclose, so a disk that fills during it would
    // otherwise be reported as a clean write. writeFasta already gets this right.
    const bool ok = std::fwrite(w.s.data(), 1, w.s.size(), f) == w.s.size();
    const bool closed = std::fclose(f) == 0;
    if (!ok || !closed) { error = "write failed on " + path; return false; }
    return true;
}

}  // namespace ts
