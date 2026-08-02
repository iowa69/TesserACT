// tessera-model: builds a genus prior from closed reference genomes.
//
//   tessera-model --organism klebsiella --out kleb.tsm ref/*.fasta
//   tessera-model --organism klebsiella --out kleb.tsm --exclude ERR123 \
//                 --plasmids plasmid_db.fasta ref/*.fasta
//
// Inputs are classified by file name: `*_chr.fasta` is chromosome, anything
// matching `*_plasmid*` is plasmid, and `--plasmids` takes a multi-record
// plasmid database where every record is its own replicon. The two classes
// are learned into separate tables, because conserved chromosomal gene order
// and mosaic plasmid structure are not the same kind of evidence.
//
// `--exclude` is the whole point of the tool being separate. When the model is
// scored against an isolate, that isolate's own reference must not be in the
// panel, or the result says only that a genome predicts itself. The excluded
// accessions are recorded inside the model file so a run can prove which
// genomes it never saw.
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <unordered_map>
#include <vector>

#include "organism.h"
#include "util.h"

namespace {

struct MarkerHit {
    uint64_t kmer;
    uint32_t pos;
    int orient;
};

// Reads a FASTA file, calling `fn(sequence)` per record so a 385 MB plasmid
// database never has to be held in memory at once.
template <typename F>
bool forEachFastaRecord(const std::string& path, F&& fn) {
    std::FILE* f = std::fopen(path.c_str(), "rb");
    if (!f) return false;
    std::string cur;
    char buf[1 << 16];
    while (std::fgets(buf, sizeof(buf), f)) {
        if (buf[0] == '>') {
            if (!cur.empty()) { fn(cur); cur.clear(); }
            continue;
        }
        for (char* p = buf; *p && *p != '\n' && *p != '\r'; ++p) cur.push_back(*p);
    }
    if (!cur.empty()) fn(cur);
    std::fclose(f);
    return true;
}

// The accession is the leading token of the file name, which is how the panel
// is laid out: ERR11578413_chr.fasta, ERR11578413_plasmid_1.fasta.
std::string accessionOf(const std::string& path) {
    size_t slash = path.find_last_of('/');
    std::string base = slash == std::string::npos ? path : path.substr(slash + 1);
    size_t cut = base.find_first_of("_.");
    return cut == std::string::npos ? base : base.substr(0, cut);
}

bool looksPlasmid(const std::string& path) {
    return path.find("plasmid") != std::string::npos;
}

// Learns one replicon set: markers single-copy within it, and their ordered
// adjacencies. A marker occurring twice names a repeat family rather than a
// locus, and is exactly the kind that cannot settle a junction.
void learnReplicons(ts::OrganismModel& model, const std::vector<std::vector<MarkerHit>>& replicons,
                    ts::Replicon cls) {
    std::unordered_map<uint64_t, uint32_t> occurrences;
    for (const auto& rep : replicons) {
        for (const MarkerHit& h : rep) ++occurrences[h.kmer];
    }

    model.noteGenome(cls);
    std::unordered_map<uint64_t, char> countedHere;

    for (const auto& rep : replicons) {
        std::vector<MarkerHit> single;
        single.reserve(rep.size());
        for (const MarkerHit& h : rep) {
            if (occurrences[h.kmer] == 1) single.push_back(h);
        }
        for (const MarkerHit& h : single) {
            const uint32_t id = model.internMarker(h.kmer);
            if (!countedHere.count(h.kmer)) {
                countedHere[h.kmer] = 1;
                model.noteMarkerGenome(id, cls);
            }
        }
        for (size_t i = 0; i < single.size(); ++i) {
            const uint32_t idA = model.internMarker(single[i].kmer);
            for (int n = 1; n <= ts::kMarkerNeighbours && i + static_cast<size_t>(n) < single.size(); ++n) {
                const MarkerHit& b = single[i + static_cast<size_t>(n)];
                const int32_t dist = static_cast<int32_t>(b.pos) - static_cast<int32_t>(single[i].pos);
                if (dist <= 0 || dist > ts::kMaxMarkerDistance) break;
                const uint32_t idB = model.internMarker(b.kmer);
                // Both the walk and its mirror are recorded, so a query need
                // not know which strand the assembly happens to be on.
                const uint64_t fwdFrom = (static_cast<uint64_t>(idA) << 1) | static_cast<uint64_t>(single[i].orient);
                const uint64_t fwdTo = (static_cast<uint64_t>(idB) << 1) | static_cast<uint64_t>(b.orient);
                model.addObservation(fwdFrom, fwdTo, dist, cls);
                const uint64_t revFrom = (static_cast<uint64_t>(idB) << 1) | static_cast<uint64_t>(1 - b.orient);
                const uint64_t revTo = (static_cast<uint64_t>(idA) << 1) | static_cast<uint64_t>(1 - single[i].orient);
                model.addObservation(revFrom, revTo, dist, cls);
            }
        }
    }
}

void usage() {
    std::fprintf(stderr,
                 "tessera-model -- build a genus prior from closed genomes\n\n"
                 "usage: tessera-model --organism NAME --out FILE [options] FASTA...\n\n"
                 "  --organism NAME     organism the model describes (e.g. klebsiella)\n"
                 "  --out FILE          model file to write\n"
                 "  --exclude ACC       omit this accession from the panel (repeatable)\n"
                 "  --plasmids FILE     multi-record plasmid database; each record is\n"
                 "                      treated as its own plasmid replicon\n"
                 "  --min-support N     chromosomal marker pairs seen in fewer than N\n"
                 "                      genomes are dropped (default 5)\n"
                 "  --min-support-plasmid N   the same for plasmids (default 3, because\n"
                 "                      any one plasmid is carried by far fewer isolates)\n");
}

}  // namespace

int main(int argc, char** argv) {
    std::string organism = "klebsiella", outPath, plasmidDb;
    std::vector<std::string> inputs, excludes;
    uint32_t minSupport = 5, minSupportPls = 3;

    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        auto value = [&](const char* what) -> std::string {
            if (i + 1 >= argc) {
                std::fprintf(stderr, "error: %s needs a value\n", what);
                std::exit(2);
            }
            return argv[++i];
        };
        if (a == "--organism") organism = value("--organism");
        else if (a == "--out") outPath = value("--out");
        else if (a == "--exclude") excludes.push_back(value("--exclude"));
        else if (a == "--plasmids") plasmidDb = value("--plasmids");
        else if (a == "--min-support") minSupport = static_cast<uint32_t>(std::atoi(value("--min-support").c_str()));
        else if (a == "--min-support-plasmid") minSupportPls = static_cast<uint32_t>(std::atoi(value("--min-support-plasmid").c_str()));
        else if (a == "-h" || a == "--help") { usage(); return 0; }
        else if (!a.empty() && a[0] == '-') {
            std::fprintf(stderr, "error: unknown option %s\n", a.c_str());
            return 2;
        } else inputs.push_back(a);
    }

    if (outPath.empty() || (inputs.empty() && plasmidDb.empty())) { usage(); return 2; }

    // Group the panel files by accession and class: a genome's chromosome and
    // its plasmids are learned separately, and marker adjacency must never run
    // across a replicon boundary.
    std::vector<std::string> order;
    std::vector<std::vector<std::string>> chrFiles, plsFiles;
    {
        std::unordered_map<std::string, size_t> at;
        for (const std::string& p : inputs) {
            const std::string acc = accessionOf(p);
            bool skip = false;
            for (const std::string& e : excludes) if (e == acc) { skip = true; break; }
            if (skip) continue;
            auto it = at.find(acc);
            if (it == at.end()) {
                at.emplace(acc, order.size());
                order.push_back(acc);
                chrFiles.emplace_back();
                plsFiles.emplace_back();
                it = at.find(acc);
            }
            (looksPlasmid(p) ? plsFiles : chrFiles)[it->second].push_back(p);
        }
    }

    if (order.empty() && plasmidDb.empty()) {
        std::fprintf(stderr, "error: every input was excluded; nothing to build from\n");
        return 2;
    }

    ts::OrganismModel model;
    model.beginBuild(organism, ts::kMarkerK);
    for (const std::string& e : excludes) model.noteExcluded(e);

    ts::util::Timer timer;
    size_t chrReplicons = 0, plsReplicons = 0;

    auto collect = [&](const std::vector<std::string>& files,
                       std::vector<std::vector<MarkerHit>>& out) {
        for (const std::string& path : files) {
            if (!forEachFastaRecord(path, [&](const std::string& s) {
                    out.emplace_back();
                    auto& rep = out.back();
                    ts::forEachMarkerKmer(s, [&](uint64_t km, uint32_t pos, int orient) {
                        rep.push_back({km, pos, orient});
                    });
                })) {
                std::fprintf(stderr, "warning: cannot read %s\n", path.c_str());
            }
        }
    };

    for (size_t gi = 0; gi < order.size(); ++gi) {
        std::vector<std::vector<MarkerHit>> chrReps, plsReps;
        collect(chrFiles[gi], chrReps);
        collect(plsFiles[gi], plsReps);
        chrReplicons += chrReps.size();
        plsReplicons += plsReps.size();
        if (!chrReps.empty()) learnReplicons(model, chrReps, ts::Replicon::Chromosome);
        // A genome's own plasmids are one plasmid observation, learned against
        // the plasmid table rather than the chromosomal one.
        if (!plsReps.empty()) learnReplicons(model, plsReps, ts::Replicon::Plasmid);

        if ((gi + 1) % 50 == 0 || gi + 1 == order.size()) {
            std::fprintf(stderr, "  %zu/%zu genomes, %zu markers (%.1fs)\n",
                         gi + 1, order.size(), model.markerCount(), timer.elapsed());
        }
    }

    if (!plasmidDb.empty()) {
        // Each record of the database is an independent plasmid: learned on its
        // own so that adjacency support counts distinct plasmids, not contigs
        // of one.
        size_t records = 0;
        const bool ok = forEachFastaRecord(plasmidDb, [&](const std::string& s) {
            std::vector<std::vector<MarkerHit>> one(1);
            ts::forEachMarkerKmer(s, [&](uint64_t km, uint32_t pos, int orient) {
                one[0].push_back({km, pos, orient});
            });
            learnReplicons(model, one, ts::Replicon::Plasmid);
            ++records;
            ++plsReplicons;
            if (records % 500 == 0) {
                std::fprintf(stderr, "  %zu plasmid records (%.1fs)\n", records, timer.elapsed());
            }
        });
        if (!ok) {
            std::fprintf(stderr, "error: cannot read plasmid database %s\n", plasmidDb.c_str());
            return 1;
        }
        std::fprintf(stderr, "  plasmid database: %zu records\n", records);
    }

    model.finalise(minSupport, minSupportPls);

    std::string error;
    if (!model.save(outPath, error)) {
        std::fprintf(stderr, "error: %s\n", error.c_str());
        return 1;
    }

    std::fprintf(stderr,
                 "\n%s model -> %s (%.1fs)\n"
                 "  chromosome: %u genomes, %zu replicons, %zu adjacencies (min support %u)\n"
                 "  plasmid:    %u sets, %zu replicons, %zu adjacencies (min support %u)\n"
                 "  %zu markers, %zu accessions excluded\n",
                 organism.c_str(), outPath.c_str(), timer.elapsed(),
                 model.genomes(ts::Replicon::Chromosome), chrReplicons,
                 model.edgeCount(ts::Replicon::Chromosome), minSupport,
                 model.genomes(ts::Replicon::Plasmid), plsReplicons,
                 model.edgeCount(ts::Replicon::Plasmid), minSupportPls,
                 model.markerCount(), excludes.size());
    return 0;
}
