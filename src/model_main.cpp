// tesseract-model: builds a genus prior from closed reference genomes.
//
//   tesseract-model --organism klebsiella --out kleb.tsm ref/*.fasta
//   tesseract-model --organism klebsiella --out kleb.tsm --exclude ERR123
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
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <unordered_map>
#include <vector>

#include "organism.h"

#include "version.h"
namespace {
using ts::kVersion;
}  // namespace
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

// As above, but hands the record's name to the callback. Needed because a plasmid
// database is one file of many replicons, so the only way to withhold one is by name.
template <typename F>
bool forEachNamedFastaRecord(const std::string& path, F&& fn) {
    std::FILE* f = std::fopen(path.c_str(), "rb");
    if (!f) return false;
    std::string cur, name;
    char buf[1 << 16];
    auto flush = [&]() {
        if (!name.empty() || !cur.empty()) fn(name, cur);
        cur.clear();
    };
    while (std::fgets(buf, sizeof(buf), f)) {
        if (buf[0] == '>') {
            flush();
            name.assign(buf + 1);
            const size_t cut = name.find_first_of(" \t\n\r");
            if (cut != std::string::npos) name.resize(cut);
            continue;
        }
        for (char* p = buf; *p && *p != '\n' && *p != '\r'; ++p) cur.push_back(*p);
    }
    flush();
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
                    ts::Replicon cls, const std::string& trackName = std::string()) {
    std::unordered_map<uint64_t, uint32_t> occurrences;
    for (const auto& rep : replicons) {
        for (const MarkerHit& h : rep) ++occurrences[h.kmer];
    }

    model.noteGenome(cls);
    std::unordered_map<uint64_t, char> countedHere;

    // A layout track is recorded for the longest replicon of the set only. A genome's
    // chromosome is what a contig is laid out against; its plasmids are separate
    // molecules, and a track spanning both would assert an order between sequences that
    // have none.
    size_t longest = 0;
    for (size_t i = 1; i < replicons.size(); ++i) {
        if (replicons[i].size() > replicons[longest].size()) longest = i;
    }

    for (size_t ri = 0; ri < replicons.size(); ++ri) {
        const auto& rep = replicons[ri];
        std::vector<MarkerHit> single;
        single.reserve(rep.size());
        for (const MarkerHit& h : rep) {
            if (occurrences[h.kmer] == 1) single.push_back(h);
        }
        if (!trackName.empty() && ri == longest && !single.empty()) {
            ts::LayoutTrack track;
            track.name = trackName;
            track.oriented.reserve(single.size());
            track.pos.reserve(single.size());
            for (const MarkerHit& h : single) {
                const uint32_t id = model.internMarker(h.kmer);
                track.oriented.push_back((id << 1) | static_cast<uint32_t>(h.orient));
                track.pos.push_back(h.pos);
            }
            model.addTrack(std::move(track));
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

void usage(std::FILE* to = stderr) {
    std::fprintf(to,
                 "tesseract-model -- build a genus prior from closed genomes\n\n"
                 "usage: tesseract-model --organism NAME --out FILE [options] FASTA...\n\n"
                 "  --organism NAME     organism the model describes (e.g. klebsiella)\n"
                 "  --out FILE          model file to write\n"
                 "  --exclude ACC       omit this accession from the panel (repeatable)\n"
                 "  --plasmids FILE     multi-record plasmid database; each record is\n"
                 "                      treated as its own plasmid replicon\n"
                 "  --exclude-plasmids FILE   one accession per line; those records are left\n"
                 "                      out of the plasmid database. Without this, --exclude\n"
                 "                      does NOT protect the plasmid table -- it filters the\n"
                 "                      chromosome inputs only. Measured on this cohort,\n"
                 "                      46%% of test plasmids have a near-identical match in\n"
                 "                      a RefSeq-derived panel, so a model built without an\n"
                 "                      exclusion list scores plasmids against themselves.\n"
                 "  --min-support N     chromosomal marker pairs seen in fewer than N\n"
                 "                      genomes are dropped (default 5)\n"
                 "  --min-support-plasmid N   the same for plasmids (default 3, because\n"
                 "                      any one plasmid is carried by far fewer isolates)\n"
                 "  --marker-density N  keep one k-mer in N as a marker (default 512,\n"
                 "                      about 500 bp spacing). Recorded in the model and\n"
                 "                      applied automatically when it is queried, so build\n"
                 "                      and query cannot disagree. 64 costs roughly 8x the\n"
                 "                      model size and build time and is what plasmid\n"
                 "                      grouping wants: at 512 a 1 kb contig expects two\n"
                 "                      markers and grouping needs three.\n"
                 "  --layout-tracks     also record each panel chromosome's marker order,\n"
                 "                      so an assembly can be laid out against the relative\n"
                 "                      it most resembles instead of joined junction by\n"
                 "                      junction. Costs roughly 80 kB per genome.\n");
}

}  // namespace

int main(int argc, char** argv) {
    // No default organism: a model built from arbitrary genomes and stamped
    // "klebsiella" is worse than one that refuses to build.
    std::string organism, outPath, plasmidDb, plasmidExcludePath;
    std::vector<std::string> inputs, excludes;
    uint32_t minSupport = 5, minSupportPls = 3;
    // Marker sampling denominator, recorded in the model so a query cannot disagree with it.
    // 512 gives roughly 500 bp spacing and is the default the chromosome result was measured
    // at; 64 costs about eight times the model size and build time and is what the plasmid
    // grouping wants -- a 1 kb contig expects two markers at 512 and grouping needs three.
    uint32_t markerDenom = ts::kMarkerSampleDenom;
    bool withTracks = false;
    bool withPlasmidMembership = true;

    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        auto value = [&](const char* what) -> std::string {
            if (i + 1 >= argc) {
                std::fprintf(stderr, "error: %s needs a value\n", what);
                std::exit(2);
            }
            return argv[++i];
        };
        // Range-checked: atoi cannot report failure, so "abc" and an overflowing
        // value both became 0, which means "keep every pair", and a negative one
        // wrapped to about four billion, which drops all of them.
        auto count = [&](const char* what) -> uint32_t {
            const std::string t = value(what);
            char* end = nullptr;
            const long long v = std::strtoll(t.c_str(), &end, 10);
            if (end == t.c_str() || (end && *end != '\0') || v < 1 || v > 1000000) {
                std::fprintf(stderr, "error: %s must be between 1 and 1000000 (got '%s')\n",
                             what, t.c_str());
                std::exit(2);
            }
            return static_cast<uint32_t>(v);
        };
        if (a == "--organism") organism = value("--organism");
        else if (a == "--out") outPath = value("--out");
        else if (a == "--exclude") excludes.push_back(value("--exclude"));
        else if (a == "--plasmids") plasmidDb = value("--plasmids");
        else if (a == "--exclude-plasmids") plasmidExcludePath = value("--exclude-plasmids");
        else if (a == "--layout-tracks") withTracks = true;
        else if (a == "--no-plasmid-membership") withPlasmidMembership = false;
        else if (a == "--min-support") minSupport = count("--min-support");
        else if (a == "--min-support-plasmid") minSupportPls = count("--min-support-plasmid");
        else if (a == "--marker-density") markerDenom = count("--marker-density");
        else if (a == "-h" || a == "--help") { usage(stdout); return 0; }
        else if (a == "-v" || a == "--version") {
            std::printf("tesseract-model %s\n", kVersion);
            return 0;
        }
        else if (!a.empty() && a[0] == '-') {
            std::fprintf(stderr, "error: unknown option %s\n", a.c_str());
            return 2;
        } else inputs.push_back(a);
    }

    if (organism.empty()) {
        std::fprintf(stderr, "error: --organism is required\n");
        return 2;
    }
    if (outPath.empty() || (inputs.empty() && plasmidDb.empty())) { usage(stderr); return 2; }

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
    if (markerDenom == 0) {
        std::fprintf(stderr, "error: --marker-density must be at least 1\n");
        return 2;
    }
    model.beginBuild(organism, ts::kMarkerK);
    model.setMarkerDenom(markerDenom);
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
                    }, markerDenom);
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
        if (!chrReps.empty()) {
            learnReplicons(model, chrReps, ts::Replicon::Chromosome,
                           withTracks ? order[gi] : std::string());
        }
        // A genome's own plasmids are one plasmid observation, learned against
        // the plasmid table rather than the chromosomal one.
        if (!plsReps.empty()) learnReplicons(model, plsReps, ts::Replicon::Plasmid);

        if ((gi + 1) % 50 == 0 || gi + 1 == order.size()) {
            std::fprintf(stderr, "  %zu/%zu genomes, %zu markers (%.1fs)\n",
                         gi + 1, order.size(), model.markerCount(), timer.elapsed());
        }
    }

    if (!plasmidDb.empty()) {
        // Accessions to withhold. The chromosome side has had whole-cluster exclusion
        // since the beginning; the plasmid side had none, and a plasmid panel is far more
        // likely to contain the isolate's own molecule than a chromosome panel is,
        // because plasmids are deposited and redeposited across submissions.
        std::unordered_map<std::string, char> plasmidExclude;
        if (!plasmidExcludePath.empty()) {
            std::FILE* xf = std::fopen(plasmidExcludePath.c_str(), "r");
            if (!xf) {
                std::fprintf(stderr, "error: cannot read --exclude-plasmids %s\n",
                             plasmidExcludePath.c_str());
                return 2;
            }
            char lb[512];
            while (std::fgets(lb, sizeof(lb), xf)) {
                std::string acc(lb);
                while (!acc.empty() && (acc.back() == '\n' || acc.back() == '\r' ||
                                        acc.back() == ' ' || acc.back() == '\t')) acc.pop_back();
                // Accept a whole FASTA header and keep the accession, so an exclusion list
                // can be produced by cutting a column out of a mash or blast report.
                const size_t sp = acc.find_first_of(" \t");
                if (sp != std::string::npos) acc = acc.substr(0, sp);
                if (!acc.empty()) plasmidExclude.emplace(acc, 1);
            }
            std::fclose(xf);
            if (plasmidExclude.empty()) {
                std::fprintf(stderr, "error: --exclude-plasmids %s listed no accessions\n",
                             plasmidExcludePath.c_str());
                return 2;
            }
            for (const auto& kv : plasmidExclude) model.noteExcluded(kv.first);
        }
        size_t skipped = 0;

        // Each record of the database is an independent plasmid: learned on its
        // own so that adjacency support counts distinct plasmids, not contigs
        // of one.
        size_t records = 0;
        const bool ok = forEachNamedFastaRecord(plasmidDb, [&](const std::string& name,
                                                               const std::string& s) {
            if (!plasmidExclude.empty() && plasmidExclude.count(name)) { ++skipped; return; }
            std::vector<std::vector<MarkerHit>> one(1);
            ts::forEachMarkerKmer(s, [&](uint64_t km, uint32_t pos, int orient) {
                one[0].push_back({km, pos, orient});
            }, markerDenom);
            learnReplicons(model, one, ts::Replicon::Plasmid);
            // Record what this plasmid CONTAINS, separately from what follows what on it.
            // Membership is the signal that survives mosaicism: measured on 5,370 contig
            // pairs with known truth, it separates same-plasmid from different-plasmid
            // pairs at AUC 0.857 against 0.416 for a degree-preserving null.
            if (withPlasmidMembership && !one[0].empty()) {
                // Only markers that are SINGLE-COPY on this plasmid. A marker occurring
                // twice here names a repeat rather than a locus, and the query side drops
                // exactly those, so including them would populate the membership table
                // with entries no lookup can ever match -- while interning k-mers the rest
                // of the model never sees and inflating the marker space for nothing.
                std::unordered_map<uint64_t, uint32_t> occ;
                for (const MarkerHit& h : one[0]) ++occ[h.kmer];
                ts::PlasmidMembership pm;
                pm.name = name;
                pm.length = static_cast<uint32_t>(s.size());
                pm.markers.reserve(one[0].size());
                for (const MarkerHit& h : one[0]) {
                    if (occ[h.kmer] != 1) continue;
                    pm.markers.push_back(model.internMarker(h.kmer));
                }
                std::sort(pm.markers.begin(), pm.markers.end());
                pm.markers.erase(std::unique(pm.markers.begin(), pm.markers.end()),
                                 pm.markers.end());
                model.addPlasmid(std::move(pm));
            }
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
        std::fprintf(stderr, "  plasmid database: %zu records (%zu excluded)\n",
                     records, skipped);
    }

    model.finalise(minSupport, minSupportPls);

    // Every input unreadable, or none carrying a usable marker, produces a file
    // that loads without complaint and can never place anything. That is worse
    // than no file, because the run that uses it looks like it worked.
    if (model.markerCount() == 0) {
        std::fprintf(stderr,
                     "error: no markers could be learned -- check the input FASTA paths "
                     "and that they hold closed replicons\n");
        return 1;
    }

    std::string error;
    if (!model.save(outPath, error)) {
        std::fprintf(stderr, "error: %s\n", error.c_str());
        return 1;
    }

    std::fprintf(stderr,
                 "\n%s model -> %s (%.1fs)\n"
                 "  chromosome: %u genomes, %zu replicons, %zu adjacencies (min support %u)\n"
                 "  plasmid:    %u sets, %zu replicons, %zu adjacencies (min support %u)\n"
                 "  layout:     %zu tracks, %zu plasmid membership sets\n"
                 "  %zu markers, %zu accessions excluded\n",
                 organism.c_str(), outPath.c_str(), timer.elapsed(),
                 model.genomes(ts::Replicon::Chromosome), chrReplicons,
                 model.edgeCount(ts::Replicon::Chromosome), minSupport,
                 model.genomes(ts::Replicon::Plasmid), plsReplicons,
                 model.edgeCount(ts::Replicon::Plasmid), minSupportPls,
                 model.trackCount(), model.plasmidCount(), model.markerCount(),
                 model.excluded().size());
    return 0;
}
