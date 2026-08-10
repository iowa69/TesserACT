// plasmid_index -- emit the marker set of every plasmid in a panel.
//
// The model stores plasmid evidence as counts and adjacencies: "how many panel plasmids
// carry this marker" and "how often does marker A precede marker B". Neither can answer
// the question a plasmid reconstruction actually needs, which is whether markers A and B
// occur together on ONE real plasmid -- beyond the 30 kb window the adjacency table spans.
//
// That co-membership question is a SET test, not an order test, and that is exactly why it
// is the right shape for plasmids: they are mosaic, so their gene order is not conserved
// and an order-based layout imposes an arrangement that does not exist. Which panel
// plasmid carries this whole collection of markers is a question mosaicism does not
// disturb.
//
// This tool emits the raw material for that test, using the same hash-sampled canonical
// 31-mers the model uses, so the marker space is identical and the two are comparable.
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include <algorithm>

#include "organism.h"

namespace {

void usage() {
    std::fprintf(stderr,
                 "plasmid_index --panel FASTA --out TSV [--min-len N] [--max-len N]\n"
                 "  Writes: accession <TAB> length <TAB> comma-separated marker hashes\n");
}

}  // namespace

int main(int argc, char** argv) {
    std::string panel, out;
    size_t minLen = 1000, maxLen = 300000;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        auto val = [&]() -> std::string {
            if (i + 1 >= argc) { usage(); std::exit(2); }
            return argv[++i];
        };
        if (a == "--panel") panel = val();
        else if (a == "--out") out = val();
        else if (a == "--min-len") minLen = static_cast<size_t>(std::atol(val().c_str()));
        else if (a == "--max-len") maxLen = static_cast<size_t>(std::atol(val().c_str()));
        else { usage(); return 2; }
    }
    if (panel.empty() || out.empty()) { usage(); return 2; }

    std::FILE* in = std::fopen(panel.c_str(), "rb");
    if (!in) { std::fprintf(stderr, "error: cannot read %s\n", panel.c_str()); return 1; }
    std::FILE* o = std::fopen(out.c_str(), "w");
    if (!o) { std::fprintf(stderr, "error: cannot write %s\n", out.c_str()); std::fclose(in); return 1; }

    std::string name, seq;
    size_t kept = 0, total = 0;
    char buf[1 << 16];

    auto emit = [&]() {
        if (name.empty() || seq.size() < minLen || seq.size() > maxLen) return;
        ++total;
        std::vector<uint64_t> markers;
        ts::forEachMarkerKmer(seq, [&](uint64_t km, uint32_t, int) { markers.push_back(km); });
        if (markers.empty()) return;
        // Sorted and deduplicated: a marker occurring twice on one plasmid names a repeat
        // there and says nothing about co-membership beyond its first occurrence.
        std::sort(markers.begin(), markers.end());
        markers.erase(std::unique(markers.begin(), markers.end()), markers.end());
        std::fprintf(o, "%s\t%zu\t", name.c_str(), seq.size());
        for (size_t i = 0; i < markers.size(); ++i) {
            std::fprintf(o, "%s%llu", i ? "," : "",
                         static_cast<unsigned long long>(markers[i]));
        }
        std::fputc('\n', o);
        ++kept;
    };

    while (std::fgets(buf, sizeof(buf), in)) {
        if (buf[0] == '>') {
            emit();
            name.assign(buf + 1);
            const size_t cut = name.find_first_of(" \t\n\r");
            if (cut != std::string::npos) name.resize(cut);
            seq.clear();
            continue;
        }
        for (char* p = buf; *p && *p != '\n' && *p != '\r'; ++p) seq.push_back(*p);
    }
    emit();
    std::fclose(in);
    std::fclose(o);
    std::fprintf(stderr, "indexed %zu of %zu plasmids -> %s\n", kept, total, out.c_str());
    return 0;
}
