// om2_nrp_build -- the nearest-relative plasmid sidecar of the Organism Model 2.0 clonal stage (C4).
//
//   om2_nrp_build MODEL.tsm PLASMID_MAP.tsv PLASMID_DB.fasta OUT.om2nrp [--hold HOLD_LIST] [--scenario NAME]
//
// For every layout-track genome of MODEL (the panel chromosomes the model was built from), the plasmid records of that
// genome (PLASMID_MAP: safe_acc <TAB> accession <TAB> plasmid_record <TAB> length) that the MODEL ITSELF carries as a
// plasmid membership set, with their marker order: canonical 31-mers sampled at the model's density, single-copy on
// the record, as (marker id << 1 | orientation) and start position -- the same as a chromosome layout track.
//
// Leakage by construction: a record is written only if its name is one of the model's own plasmid membership sets, so
// nothing the scenario withheld (build_model.py exclusions, the LTO / LCO2 extra lists) can appear; a genome is written
// only if it is one of the model's tracks. Record names are matched exactly, else with the NZ_ prefix and the version
// stripped on both sides (the curated database keeps GenBank twins of RefSeq records).
//
// Output (text): header lines
//   #om2nrp 1 / #tsm_md5 <md5 of MODEL> / #tracks_sha256 <sha256 of the sorted track names, one per line> /
//   #hold_md5 <md5 of the scenario's hold list, or NA> / #scenario <name> / #denom <model density> /
//   #genomes <n> / #records <n>   (EVAL_PLAN_CLONAL s2.4: .tsm md5, hold-list md5, track-list sha256)
// then one line per record: P <TAB> genome <TAB> record <TAB> length <TAB> n <TAB> oriented,... <TAB> pos,...
// The clonal stage refuses a sidecar whose #tsm_md5 is not the loaded model's.
#include <algorithm>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "om2_alloc.h"
#include "om2_output.h"
#include "organism.h"

using namespace ts;

namespace {
std::string norm(std::string x) {
    if (x.compare(0, 3, "NZ_") == 0) x = x.substr(3);
    const size_t d = x.find('.');
    if (d != std::string::npos) x.resize(d);
    return x;
}
}  // namespace

int main(int argc, char** argv) {
    if (argc < 5) {
        std::fprintf(stderr, "usage: om2_nrp_build MODEL.tsm PLASMID_MAP.tsv PLASMID_DB.fasta OUT.om2nrp "
                             "[--hold HOLD_LIST] [--scenario NAME]\n");
        return 2;
    }
    const std::string modelPath = argv[1], mapPath = argv[2], dbPath = argv[3], outPath = argv[4];
    std::string holdMd5 = "NA", scenario = "NA";
    for (int i = 5; i + 1 < argc; i += 2) {
        const std::string a = argv[i];
        if (a == "--hold") holdMd5 = om2::md5File(argv[i + 1]);
        else if (a == "--scenario") scenario = argv[i + 1];
        else { std::fprintf(stderr, "error: unknown option %s\n", a.c_str()); return 2; }
    }
    if (holdMd5.empty()) { std::fprintf(stderr, "error: cannot read the hold list\n"); return 2; }
    OrganismModel m;
    std::string err;
    if (!m.load(modelPath, err)) { std::fprintf(stderr, "error: %s\n", err.c_str()); return 1; }
    std::set<std::string> tracks;
    for (const LayoutTrack& t : m.tracks()) tracks.insert(t.name);
    std::unordered_map<std::string, std::string> exact, byNorm;   // model membership names
    for (const PlasmidMembership& p : m.plasmids()) { exact[p.name] = p.name; byNorm.emplace(norm(p.name), p.name); }
    // model record name -> genomes carrying it
    std::map<std::string, std::vector<std::string>> owners;
    {
        std::ifstream in(mapPath);
        if (!in) { std::fprintf(stderr, "error: cannot read %s\n", mapPath.c_str()); return 1; }
        std::string l;
        size_t rows = 0, kept = 0;
        while (std::getline(in, l)) {
            std::vector<std::string> c;
            std::stringstream ss(l);
            std::string f;
            while (std::getline(ss, f, '\t')) c.push_back(f);
            if (c.size() < 3 || c[0] == "safe_acc") continue;
            ++rows;
            if (!tracks.count(c[0])) continue;
            std::string name;
            auto e = exact.find(c[2]);
            if (e != exact.end()) name = e->second;
            else {
                auto n = byNorm.find(norm(c[2]));
                if (n != byNorm.end()) name = n->second;
            }
            if (name.empty()) continue;
            auto& v = owners[name];
            if (std::find(v.begin(), v.end(), c[0]) == v.end()) { v.push_back(c[0]); ++kept; }
        }
        std::fprintf(stderr, "plasmid_map: %zu rows, %zu (genome, record) pairs on model tracks and in the model\n", rows, kept);
    }
    // tracks sha256
    std::string trackSha;
    {
        const std::string tmp = outPath + ".tracks.tmp";
        std::ofstream t(tmp);
        for (const std::string& n : tracks) t << n << '\n';
        t.close();
        trackSha = om2::sha256File(tmp);
        std::remove(tmp.c_str());
    }
    const uint64_t denom = m.markerDenom();
    std::vector<std::string> lines;
    std::set<std::string> genomes;
    std::FILE* f = std::fopen(dbPath.c_str(), "rb");
    if (!f) { std::fprintf(stderr, "error: cannot read %s\n", dbPath.c_str()); return 1; }
    std::string cur, name;
    char buf[1 << 16];
    auto flush = [&]() {
        auto it = owners.find(name);
        if (it == owners.end() || cur.empty()) { cur.clear(); return; }
        std::vector<std::pair<uint64_t, std::pair<uint32_t, int>>> hits;
        forEachMarkerKmer(cur, [&](uint64_t km, uint32_t pos, int o) { hits.push_back({km, {pos, o}}); }, denom);
        std::unordered_map<uint64_t, int> occ;
        for (const auto& h : hits) ++occ[h.first];
        std::string ids, pos;
        size_t n = 0;
        for (const auto& h : hits) {
            if (occ[h.first] != 1) continue;
            const uint32_t id = m.markerOf(h.first);
            if (id == UINT32_MAX) continue;
            if (n) { ids += ','; pos += ','; }
            ids += std::to_string((id << 1) | static_cast<uint32_t>(h.second.second));
            pos += std::to_string(h.second.first);
            ++n;
        }
        for (const std::string& g : it->second) {
            genomes.insert(g);
            lines.push_back("P\t" + g + "\t" + name + "\t" + std::to_string(cur.size()) + "\t" + std::to_string(n) + "\t" +
                            (n ? ids : ".") + "\t" + (n ? pos : "."));
        }
        cur.clear();
    };
    while (std::fgets(buf, sizeof buf, f)) {
        if (buf[0] == '>') {
            flush();
            name.assign(buf + 1);
            const size_t cut = name.find_first_of(" \t\n\r");
            if (cut != std::string::npos) name.resize(cut);
            continue;
        }
        if (!owners.count(name)) continue;
        for (char* p = buf; *p && *p != '\n' && *p != '\r'; ++p) cur.push_back(*p);
    }
    flush();
    std::fclose(f);
    std::sort(lines.begin(), lines.end());
    std::ofstream out(outPath);
    out << "#om2nrp\t1\n#tsm_md5\t" << om2::md5File(modelPath) << "\n#tracks_sha256\t" << trackSha << "\n#hold_md5\t"
        << holdMd5 << "\n#scenario\t" << scenario << "\n#denom\t" << denom
        << "\n#genomes\t" << genomes.size() << "\n#records\t" << lines.size() << "\n";
    for (const std::string& l : lines) out << l << '\n';
    out.close();
    std::fprintf(stderr, "wrote %s: %zu records of %zu genomes (of %zu tracks)\n", outPath.c_str(), lines.size(),
                 genomes.size(), tracks.size());
    return 0;
}
