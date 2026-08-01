#include "gfa.h"

#include <cstdio>
#include <string>
#include <vector>

namespace ts {

namespace {

// Leaving a unitig through end 1 means traversing it forward; through end 0
// means traversing its reverse complement. Entering at end 0 is forward.
inline char exitSign(int end) { return end == 1 ? '+' : '-'; }
inline char entrySign(int toEnd) { return toEnd == 0 ? '+' : '-'; }

}  // namespace

bool writeGfa(const std::string& path, const UnitigGraph& g,
              const std::vector<GfaPath>& paths, size_t& segmentCount, size_t& linkCount,
              std::string& error) {
    segmentCount = 0;
    linkCount = 0;

    std::FILE* f = std::fopen(path.c_str(), "wb");
    if (!f) { error = "cannot write " + path; return false; }
    std::setvbuf(f, nullptr, _IOFBF, 1u << 20);

    std::fprintf(f, "H\tVN:Z:1.0\n");

    // Segment names are the node indices, so links and paths can reference them
    // without a second lookup table.
    std::vector<char> live(g.nodes.size(), 0);
    for (uint32_t u = 0; u < g.nodes.size(); ++u) {
        if (g.nodes[u].deleted) continue;
        live[u] = 1;
        const Unitig& n = g.nodes[u];
        std::fprintf(f, "S\t%u\t%s\tLN:i:%zu\tdp:f:%.4f\n",
                     u, n.seq.c_str(), n.seq.size(), n.coverage);
        ++segmentCount;
    }

    const int ov = g.k() - 1;
    for (uint32_t u = 0; u < g.nodes.size(); ++u) {
        if (!live[u]) continue;
        for (int e = 0; e < 2; ++e) {
            for (const Link& l : g.nodes[u].ends[e]) {
                if (!live[l.to]) continue;
                // Every link is stored on both sides; emit each once.
                const uint64_t a = (static_cast<uint64_t>(u) << 1) | static_cast<uint64_t>(e);
                const uint64_t b = (static_cast<uint64_t>(l.to) << 1) | static_cast<uint64_t>(l.toEnd);
                if (a > b) continue;
                std::fprintf(f, "L\t%u\t%c\t%u\t%c\t%dM\n",
                             u, exitSign(e), l.to, entrySign(l.toEnd), ov);
                ++linkCount;
            }
        }
    }

    for (const GfaPath& p : paths) {
        if (p.oriented.empty()) continue;
        std::fprintf(f, "P\t%s\t", p.name.c_str());
        for (size_t i = 0; i < p.oriented.size(); ++i) {
            const uint32_t u = static_cast<uint32_t>(p.oriented[i] >> 1);
            const int orient = static_cast<int>(p.oriented[i] & 1);
            if (i) std::fputc(',', f);
            std::fprintf(f, "%u%c", u, orient == 0 ? '+' : '-');
        }
        std::fputc('\t', f);
        // Overlaps between consecutive steps; a scaffold gap has no overlap, so
        // it is recorded as a jump rather than a (k-1) match.
        for (size_t i = 1; i < p.oriented.size(); ++i) {
            if (i > 1) std::fputc(',', f);
            const int gap = (i < p.gaps.size()) ? p.gaps[i] : 0;
            if (gap > 0) std::fprintf(f, "%dN", gap);
            else std::fprintf(f, "%dM", ov);
        }
        if (p.oriented.size() == 1) std::fputc('*', f);
        std::fputc('\n', f);
    }

    const bool ok = std::ferror(f) == 0;
    std::fclose(f);
    if (!ok) { error = "write failed on " + path; return false; }
    return true;
}

}  // namespace ts
