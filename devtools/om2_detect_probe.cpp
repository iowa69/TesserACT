// Organism Model 2.0 C3: detection sketch builder and scorer (development tool; `make om2probe`,
// never built by `all`). It calls the same library code the assembler runs for --organism auto
// and TESSERACT_OM2_DETECT (src/organism_detect.cpp).
//
//   om2_detect_probe build OUT.sketch [--hold NAME=PATH]... MODEL.tsm...
//       core markers of each model (single-copy on >= 90% of its panel chromosomes, never on a
//       plasmid set) -> one sketch; records each .tsm md5 and, per --hold, the hold list md5
//   om2_detect_probe show SKETCH
//   om2_detect_probe score SKETCH LIST.tsv          LIST: sp iso r1 r2 (r2 may be -)
//       prints: sp iso reads bases model core core_hit score seconds
#include <cstdio>
#include <cstring>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "organism_detect.h"

using namespace ts::om2;

int main(int argc, char** argv) {
    if (argc < 3) {
        std::fprintf(stderr, "usage: om2_detect_probe build OUT.sketch [--hold NAME=PATH]... MODEL.tsm...\n"
                             "       om2_detect_probe show SKETCH\n"
                             "       om2_detect_probe score SKETCH LIST.tsv\n");
        return 2;
    }
    const std::string cmd = argv[1];
    std::string err;
    if (cmd == "build") {
        const std::string out = argv[2];
        std::map<std::string, std::string> hold;
        DetectSketch sk;
        for (int i = 3; i < argc; ++i) {
            if (std::strcmp(argv[i], "--hold") == 0 && i + 1 < argc) {
                const std::string h = argv[++i];
                const size_t eq = h.find('=');
                if (eq == std::string::npos) { std::fprintf(stderr, "--hold wants NAME=PATH\n"); return 2; }
                hold[h.substr(0, eq)] = h.substr(eq + 1);
                continue;
            }
            DetectSketchModel m;
            if (!sketchFromModel(argv[i], m, err)) { std::fprintf(stderr, "%s\n", err.c_str()); return 2; }
            auto it = hold.find(m.name);
            if (it != hold.end()) {
                m.holdList = it->second;
                m.holdMd5 = md5OfFile(it->second);
                if (m.holdMd5.empty()) { std::fprintf(stderr, "cannot read hold list %s\n", it->second.c_str()); return 2; }
            }
            std::fprintf(stderr, "%s\t%s\tgChr=%u\tdenom=%u\tcore=%zu\ttsm_md5=%s\thold=%s\n", argv[i], m.name.c_str(),
                         m.genomesChr, m.denom, m.core.size(), m.tsmMd5.c_str(), m.holdMd5.c_str());
            sk.models.push_back(std::move(m));
        }
        if (!writeDetectSketch(out, sk, err)) { std::fprintf(stderr, "%s\n", err.c_str()); return 2; }
        return 0;
    }
    DetectSketch sk;
    if (!loadDetectSketch(argv[2], sk, err)) { std::fprintf(stderr, "%s\n", err.c_str()); return 2; }
    if (cmd == "show") {
        std::printf("sketch\t%s\tmd5=%s\tmodels=%zu\n", argv[2], md5OfFile(argv[2]).c_str(), sk.models.size());
        for (const DetectSketchModel& m : sk.models)
            std::printf("%s\tgChr=%u\tdenom=%u\tcore=%zu\ttsm_md5=%s\thold=%s\thold_md5=%s\n", m.name.c_str(), m.genomesChr,
                        m.denom, m.core.size(), m.tsmMd5.c_str(), m.holdList.c_str(), m.holdMd5.c_str());
        return 0;
    }
    if (cmd == "score" && argc >= 4) {
        std::ifstream in(argv[3]);
        std::string ln;
        std::printf("sp\tiso\treads\tbases\tmodel\tcore\tcore_hit\tscore\tseconds\n");
        while (std::getline(in, ln)) {
            std::istringstream ss(ln);
            std::string sp, iso, r1, r2;
            ss >> sp >> iso >> r1 >> r2;
            if (r1.empty()) continue;
            std::vector<std::string> files{r1};
            if (!r2.empty() && r2 != "-") files.push_back(r2);
            DetectResult r;
            if (!scoreReadFiles(sk, files, r, err)) { std::fprintf(stderr, "%s: %s\n", iso.c_str(), err.c_str()); continue; }
            for (const DetectModelScore& s : r.scores)
                std::printf("%s\t%s\t%llu\t%llu\t%s\t%u\t%u\t%.6f\t%.1f\n", sp.c_str(), iso.c_str(),
                            static_cast<unsigned long long>(r.reads), static_cast<unsigned long long>(r.bases),
                            s.name.c_str(), s.core, s.coreHit, s.score, r.seconds);
            std::fflush(stdout);
        }
        return 0;
    }
    std::fprintf(stderr, "unknown command %s\n", cmd.c_str());
    return 2;
}
