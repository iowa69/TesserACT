// join_probe -- run only the model join stage, on contigs that already exist.
//
// Assembling a 160x library takes half an hour, almost all of it spent counting
// k-mers and resolving the graph -- work that is identical no matter what the join
// stage is asked to do. Iterating on join thresholds through the full assembler
// therefore costs about a thousand times more than the change being measured.
//
// This driver loads a finished contig set and a model, runs joinByModel, and reports
// what happened. Seconds instead of half an hour, which is the difference between
// measuring a dozen variants and guessing at one.
//
// It is a development tool, deliberately outside src/ so it never links into the
// assembler, and it makes no attempt to reproduce the rest of the pipeline: the
// contigs it reads have already been through error correction, graph cleaning and
// paired resolution, and the gap closer that would normally run afterwards is absent.
// Numbers from here are comparable to each other, not to a full run.
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include "organism.h"
#include "organism_join.h"
#include "organism_layout.h"
#include "replicon.h"

namespace {

bool readFasta(const std::string& path, std::vector<std::string>& seqs,
               std::vector<std::string>& names) {
    std::FILE* f = std::fopen(path.c_str(), "rb");
    if (!f) return false;
    std::string cur, name;
    char buf[1 << 16];
    auto flush = [&]() {
        if (!name.empty()) {
            seqs.push_back(cur);
            names.push_back(name);
        }
        cur.clear();
    };
    while (std::fgets(buf, sizeof(buf), f)) {
        if (buf[0] == '>') {
            flush();
            name.assign(buf + 1);
            while (!name.empty() && (name.back() == '\n' || name.back() == '\r')) name.pop_back();
            continue;
        }
        for (char* p = buf; *p && *p != '\n' && *p != '\r'; ++p) cur.push_back(*p);
    }
    flush();
    std::fclose(f);
    return true;
}

// Coverage as written by the assembler into the contig name, so the join stage sees
// the same numbers it would in a real run. Absent, every contig gets 1.0, which the
// join stage does not currently use but which must not be left uninitialised.
double covFromName(const std::string& name) {
    const size_t at = name.find("cov=");
    if (at == std::string::npos) return 1.0;
    return std::atof(name.c_str() + at + 4);
}

void usage() {
    std::fprintf(stderr,
                 "join_probe --model FILE --contigs FASTA [--out FASTA] [--k N]\n"
                 "           [--is-panel FASTA] [--is-sites TSV] [--layout] [--quiet]\n");
}

}  // namespace

int main(int argc, char** argv) {
    std::string modelPath, contigPath, outPath, isPanelPath, isSitesPath;
    int k = 127;
    bool verbose = true;
    bool layout = false;

    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        auto val = [&]() -> std::string {
            if (i + 1 >= argc) { usage(); std::exit(2); }
            return argv[++i];
        };
        if (a == "--model") modelPath = val();
        else if (a == "--contigs") contigPath = val();
        else if (a == "--out") outPath = val();
        else if (a == "--is-panel") isPanelPath = val();
        else if (a == "--is-sites") isSitesPath = val();
        else if (a == "--k") k = std::atoi(val().c_str());
        else if (a == "--layout") layout = true;
        else if (a == "--quiet") verbose = false;
        else { usage(); return 2; }
    }
    if (modelPath.empty() || contigPath.empty()) { usage(); return 2; }

    ts::OrganismModel model;
    std::string err;
    if (!model.load(modelPath, err)) {
        std::fprintf(stderr, "error: %s\n", err.c_str());
        return 1;
    }

    std::vector<std::string> contigs, names;
    if (!readFasta(contigPath, contigs, names)) {
        std::fprintf(stderr, "error: cannot read %s\n", contigPath.c_str());
        return 1;
    }
    if (contigs.empty()) {
        std::fprintf(stderr, "error: no contigs in %s\n", contigPath.c_str());
        return 1;
    }
    std::vector<double> covs;
    covs.reserve(names.size());
    for (const std::string& nm : names) covs.push_back(covFromName(nm));

    ts::IsPanel panel;
    ts::IsPanel* panelPtr = nullptr;
    if (!isPanelPath.empty()) {
        if (!panel.load(isPanelPath, err)) {
            std::fprintf(stderr, "error: %s\n", err.c_str());
            return 1;
        }
        panelPtr = &panel;
    }
    if (!isSitesPath.empty()) {
        if (!panel.loadSites(isSitesPath, err)) {
            std::fprintf(stderr, "error: %s\n", err.c_str());
            return 1;
        }
        panelPtr = &panel;
    }

    size_t basesIn = 0;
    for (const std::string& s : contigs) basesIn += s.size();

    const ts::OrganismJoinStats st =
        ts::joinByModel(model, contigs, covs, k, verbose, panelPtr, nullptr);

    ts::LayoutStats ls;
    if (layout) ls = ts::layoutByModel(model, contigs, covs, k, verbose);

    // Classify replicons so the tagged output can be scored for grouping. Read pairs are
    // not available here -- the probe never sees the reads -- so this exercises the panel
    // co-membership path alone, which is exactly the part under test.
    std::vector<char> layoutMembers;
    if (layout) {
        layoutMembers.assign(contigs.size(), 0);
        for (size_t i = 0; i < ls.scaffolds && i < layoutMembers.size(); ++i) layoutMembers[i] = 1;
    }
    const ts::RepliconAssignment rep =
        ts::assignReplicons(contigs, covs, layoutMembers, model, nullptr, verbose);

    size_t basesOut = 0, longest = 0;
    for (const std::string& s : contigs) {
        basesOut += s.size();
        if (s.size() > longest) longest = s.size();
    }

    std::printf("contigs_in\t%zu\ncontigs_out\t%zu\njoins\t%zu\n"
                "chromosome_joins\t%zu\nplasmid_joins\t%zu\ninsertion_site_joins\t%zu\n"
                "markers_placed\t%zu\ncandidates\t%zu\nrejected_weak\t%zu\n"
                "rejected_inconsistent\t%zu\nrejected_not_mutual\t%zu\n"
                "rejected_insertion_sequence\t%zu\noverlap_merges\t%zu\n"
                "overlap_unconfirmed\t%zu\nrounds\t%zu\ngap_bases\t%zu\n"
                "bases_in\t%zu\nbases_out\t%zu\nlongest\t%zu\n",
                st.contigsIn, st.contigsOut, st.joins, st.chromosomeJoins, st.plasmidJoins,
                st.insertionSiteJoins, st.markersFound, st.candidates, st.rejectedWeak,
                st.rejectedInconsistent, st.rejectedNotMutual, st.rejectedInsertionSeq,
                st.overlapMerges, st.overlapUnconfirmed, st.rounds,
                st.gapBases, basesIn, basesOut, longest);
    if (layout) {
        std::printf("layout_run\t%d\nlayout_track\t%s\nlayout_shared_markers\t%zu\n"
                    "layout_placed\t%zu\nlayout_unplaced\t%zu\nlayout_incoherent\t%zu\n"
                    "layout_dropped\t%zu\nlayout_scaffolds\t%zu\nlayout_gap_bases\t%zu\n"
                    "layout_overlap_merges\t%zu\n",
                    ls.run ? 1 : 0, ls.track.empty() ? "-" : ls.track.c_str(),
                    ls.sharedMarkers, ls.placed, ls.unplaced, ls.incoherent, ls.dropped,
                    ls.scaffolds, ls.gapBases, ls.overlapMerges);
        std::printf("replicon_chr\t%zu\nreplicon_plasmid\t%zu\nreplicon_groups\t%zu\n"
                    "replicon_unassigned\t%zu\n",
                    rep.chromosomeContigs, rep.plasmidContigs, rep.plasmidGroups,
                    rep.unassignedContigs);
    }

    if (!outPath.empty()) {
        std::FILE* o = std::fopen(outPath.c_str(), "w");
        if (!o) {
            std::fprintf(stderr, "error: cannot write %s\n", outPath.c_str());
            return 1;
        }
        for (size_t i = 0; i < contigs.size(); ++i) {
            char tag[48];
            tag[0] = '\0';
            if (i < rep.calls.size()) {
                const ts::RepliconCall& rc = rep.calls[i];
                const char* cls = rc.cls == ts::RepliconClass::Chromosome ? "chr"
                                : rc.cls == ts::RepliconClass::Plasmid    ? "plas"
                                                                          : "unk";
                if (rc.cls == ts::RepliconClass::Plasmid && rc.group > 0) {
                    std::snprintf(tag, sizeof(tag), "_%s_%u", cls, rc.group);
                } else {
                    std::snprintf(tag, sizeof(tag), "_%s", cls);
                }
            }
            std::fprintf(o, ">probe_%zu%s len=%zu\n", i, tag, contigs[i].size());
            for (size_t at = 0; at < contigs[i].size(); at += 80) {
                std::fprintf(o, "%.*s\n", static_cast<int>(std::min<size_t>(80, contigs[i].size() - at)),
                             contigs[i].c_str() + at);
            }
        }
        std::fclose(o);
    }
    return 0;
}
