#include "replicon.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <functional>
#include <map>
#include <unordered_map>

namespace ts {

namespace {

// A contig this short carries too few markers and too little depth signal to classify,
// and calling it either way would be guessing. It is reported unassigned.
constexpr size_t kMinClassifiable = 500;

// Depth ratios. The chromosome sits at 1.0 by construction (it defines the mode). Below
// the floor a contig is contamination, a fragment, or a mis-assembly; above the ceiling it
// is multi-copy and cannot be chromosomal.
constexpr double kMultiCopy = 1.75;
constexpr double kLowDepth = 0.35;

// Marker votes must be this lopsided before they decide. The model's own join stage uses
// the same 3:1 rule, and using a different one here would let the two stages disagree
// about what a contig is.
constexpr uint32_t kVoteRatio = 3;

// Pairs linking two contigs before they are called co-resident on one molecule. Low
// absolute counts are dominated by mismapping, and IS elements shared between replicons
// are exactly the sequence that links a plasmid to the chromosome.
constexpr size_t kMinCoResident = 8;

// Markers each side of a co-membership test must carry. Asymmetric on purpose: see the
// comment at the test itself for what the two numbers are doing and why the old symmetric
// requirement of three excluded most of the sequence that needed grouping.
// Measured over 40 fold-0 isolates by plasmids delivered whole in one group: 3/3 gave 1/63,
// 5/1 gave 8/65, 3/1 gave 9/65, 8/1 gave 8/65, 5/2 gave 6/65. The poor side is what matters
// -- raising it to two loses ground, and the rich side is flat from three to eight.
constexpr size_t kMinMarkersWeak = 1;
constexpr size_t kMinMarkersStrong = 3;

}  // namespace

const char* repliconName(RepliconClass c) {
    switch (c) {
        case RepliconClass::Chromosome: return "chromosome";
        case RepliconClass::Plasmid:    return "plasmid";
        default:                        return "unassigned";
    }
}

RepliconAssignment assignReplicons(const std::vector<std::string>& contigs,
                                   const std::vector<double>& covs,
                                   const std::vector<char>& layoutMembers,
                                   const OrganismModel& model,
                                   const ContigEndLinks* links,
                                   bool verbose) {
    RepliconAssignment out;
    const size_t n = contigs.size();
    out.calls.assign(n, RepliconCall{});
    if (n == 0) return out;

    // ---- modal depth --------------------------------------------------------
    // Length-weighted median rather than the mean: the mean is dragged by a handful of
    // very high-copy small plasmids, and the chromosome is most of the sequence, so the
    // length-weighted median lands on it.
    {
        std::vector<std::pair<double, size_t>> byCov;
        byCov.reserve(n);
        size_t total = 0;
        for (size_t i = 0; i < n; ++i) {
            byCov.emplace_back(covs[i], contigs[i].size());
            total += contigs[i].size();
        }
        std::sort(byCov.begin(), byCov.end());
        size_t acc = 0;
        for (const auto& kv : byCov) {
            acc += kv.second;
            if (acc * 2 >= total) { out.modalDepth = kv.first; break; }
        }
    }
    if (out.modalDepth <= 0) out.modalDepth = 1.0;

    // ---- marker votes -------------------------------------------------------
    std::vector<uint32_t> chrVotes(n, 0), plsVotes(n, 0);
    if (model.loaded()) {
        std::unordered_map<uint32_t, uint32_t> seen;
        std::vector<std::pair<size_t, uint32_t>> hits;
        for (size_t c = 0; c < n; ++c) {
            forEachMarkerKmer(contigs[c], [&](uint64_t km, uint32_t, int) {
                const uint32_t id = model.markerOf(km);
                if (id == UINT32_MAX) return;
                hits.emplace_back(c, id);
                ++seen[id];
            }, model.markerDenom());
        }
        for (const auto& h : hits) {
            // A marker occurring more than once in this assembly names a repeat here, not
            // a locus, and cannot say which molecule its contig belongs to.
            if (seen[h.second] != 1) continue;
            const uint32_t gc = model.markerGenomes(h.second, Replicon::Chromosome);
            const uint32_t gp = model.markerGenomes(h.second, Replicon::Plasmid);
            if (gc > 0 && gp == 0) ++chrVotes[h.first];
            else if (gp > 0 && gc == 0) ++plsVotes[h.first];
        }
    }

    // ---- per contig ---------------------------------------------------------
    for (size_t i = 0; i < n; ++i) {
        RepliconCall& call = out.calls[i];
        call.depthRatio = out.modalDepth > 0 ? covs[i] / out.modalDepth : 0;
        if (links && links->usable()) {
            const PairEvidence ev = links->circular(static_cast<uint32_t>(i));
            call.circular = ev.consistent >= links->supportFloor();
        }

        if (contigs[i].size() < kMinClassifiable) {
            call.basis = "too_short";
            continue;
        }
        // Layout placement is the strongest evidence there is: the contig was ordered
        // against a panel chromosome and agreed with it.
        if (i < layoutMembers.size() && layoutMembers[i]) {
            call.cls = RepliconClass::Chromosome;
            call.basis = "layout";
            continue;
        }
        if (chrVotes[i] >= plsVotes[i] * kVoteRatio && chrVotes[i] > 0) {
            // Markers say chromosome, but depth can overrule: a contig at three times the
            // chromosome's depth is not part of a single-copy molecule whatever its
            // markers resemble. Shared mobile elements are exactly how a plasmid contig
            // acquires chromosomal-looking markers.
            if (call.depthRatio >= kMultiCopy) {
                call.cls = RepliconClass::Plasmid;
                call.basis = "depth_over_markers";
            } else {
                call.cls = RepliconClass::Chromosome;
                call.basis = "markers";
            }
            continue;
        }
        if (plsVotes[i] >= chrVotes[i] * kVoteRatio && plsVotes[i] > 0) {
            call.cls = RepliconClass::Plasmid;
            call.basis = "markers";
            continue;
        }
        if (call.depthRatio >= kMultiCopy) {
            call.cls = RepliconClass::Plasmid;
            call.basis = "depth";
            continue;
        }
        if (call.depthRatio <= kLowDepth) {
            call.basis = "low_depth";
            continue;
        }
        call.basis = "no_signal";
    }

    // ---- propagate class along read pairs -----------------------------------
    // Most contigs end up unassigned, and they are the short ones: on KP1000, 69 of 103.
    // That is a classification failure that precedes grouping, because a contig never
    // called plasmid can never be put on a plasmid. The cause is the same sampling arithmetic
    // that broke co-membership -- one marker per 512 31-mers means a 1 kb contig expects two
    // markers and often has none, so it falls through every marker test to "no_signal".
    //
    // Read pairs do not have that problem. A 1 kb contig carries hundreds of fragments whose
    // partners land on its neighbours, and a neighbour that has already been classified can
    // say what this contig is. One round only, and only from contigs classified on their own
    // evidence: iterating would let a single wrong call spread through the assembly, and the
    // point of this pass is to reach contigs with no evidence, not to overrule contigs that
    // have some.
    if (links && links->usable()) {
        const std::vector<RepliconCall> seed = out.calls;
        for (size_t i = 0; i < n; ++i) {
            if (out.calls[i].cls != RepliconClass::Unassigned) continue;
            // Unassigned is two populations, and only one of them is evidence-free. A
            // contig the depth floor rejected carries a verdict -- contamination, or a
            // fragment at a coverage no replicon in this isolate has -- and propagating a
            // class onto it overrules that verdict, which is what the paragraph above says
            // this pass must not do. Without this line a 0.1x contaminant that mismaps
            // enough pairs onto the chromosome is emitted as `_chr`.
            if (std::strcmp(out.calls[i].basis, "low_depth") == 0) continue;
            if (contigs[i].size() < kMinClassifiable) continue;
            size_t chrW = 0, plsW = 0;
            for (size_t j = 0; j < n; ++j) {
                if (i == j || seed[j].cls == RepliconClass::Unassigned) continue;
                const size_t w = links->linkWeight(static_cast<uint32_t>(i),
                                                   static_cast<uint32_t>(j));
                if (w < kMinCoResident) continue;
                (seed[j].cls == RepliconClass::Chromosome ? chrW : plsW) += w;
            }
            // The same 3:1 rule the marker votes use, so the two stages cannot disagree
            // about what counts as decisive.
            if (chrW > 0 && chrW >= plsW * kVoteRatio) {
                out.calls[i].cls = RepliconClass::Chromosome;
                out.calls[i].basis = "pair_propagated";
            } else if (plsW > 0 && plsW >= chrW * kVoteRatio) {
                out.calls[i].cls = RepliconClass::Plasmid;
                out.calls[i].basis = "pair_propagated";
            }
        }
    }

    // ---- a chromosome fragment is not a closed molecule ----------------------
    // The only circularity call that was wrong across 30 isolates was a 2,800 bp contig
    // classified chromosomal, and it was wrong for a reason the classification already
    // knows: a contig assigned to the chromosome is part of a five-megabase molecule, so
    // whatever made its ends look joined was a collapsed repeat rather than a circle. A
    // genuinely circularised chromosome must still be allowed, so the test is relative --
    // the contig has to be most of the chromosomal sequence in the assembly, not merely
    // some of it.
    //
    // This must run AFTER class propagation, not before. It first ran before, and KP7 came
    // back with an 845 bp contig tagged `_chr_circular`: it was unassigned when the guard
    // looked at it and was propagated to chromosome afterwards, keeping a tag the guard
    // would have removed.
    {
        size_t chrLongest = 0;
        for (size_t i = 0; i < n; ++i) {
            if (out.calls[i].cls == RepliconClass::Chromosome) {
                chrLongest = std::max(chrLongest, contigs[i].size());
            }
        }
        for (size_t i = 0; i < n; ++i) {
            if (!out.calls[i].circular) continue;
            if (out.calls[i].cls != RepliconClass::Chromosome) continue;
            if (contigs[i].size() * 2 < chrLongest) out.calls[i].circular = false;
        }
    }

    // ---- group the plasmid contigs into molecules ---------------------------
    // Contigs on one plasmid are linked by read pairs; contigs on different molecules
    // should not be. Union-find over links above the floor. Without pair evidence every
    // plasmid contig becomes its own group, which is honest -- it says the contigs are
    // plasmid but does not claim to know how many plasmids there are.
    std::vector<uint32_t> parent(n);
    for (size_t i = 0; i < n; ++i) parent[i] = static_cast<uint32_t>(i);
    std::function<uint32_t(uint32_t)> find = [&](uint32_t x) {
        while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; }
        return x;
    };
    // Panel co-membership: does one real plasmid in the panel carry markers from both
    // contigs? This is the signal that survives mosaicism, because it asks what occurs
    // together rather than in what order -- and it is the only one of the three candidates
    // that works. Measured against truth on 5,370 contig pairs it separates same-plasmid
    // from different-plasmid at AUC 0.857, and a degree-preserving null of the same graph
    // falls to 0.416, below chance. Depth cannot do this (43% of same-isolate plasmid pairs
    // sit within 1.5x of each other) and read pairs cannot reach across the repeats that
    // separate plasmid contigs.
    std::vector<std::vector<uint32_t>> contigMarkers(n);
    if (model.plasmidCount() > 0) {
        // Occurrence is counted over the WHOLE assembly, not over the plasmid contigs alone.
        // Only markers occurring exactly once here are interned, and that claim is the
        // entire safety argument for co-membership and for the hub: it is what keeps IS
        // copies and other mobile elements out. Counting within the plasmid subset instead
        // made "single copy" mean "single copy among plasmids", so an IS26 marker sitting in
        // twelve chromosomal copies and one plasmid copy passed the filter -- and that is
        // exactly the marker that matches hundreds of panel plasmids and would join two
        // unrelated molecules into one group.
        std::unordered_map<uint32_t, uint32_t> seenM;
        std::vector<std::pair<size_t, uint32_t>> hits;
        for (size_t c = 0; c < n; ++c) {
            const bool keep = out.calls[c].cls == RepliconClass::Plasmid;
            forEachMarkerKmer(contigs[c], [&](uint64_t km, uint32_t, int) {
                const uint32_t id = model.markerOf(km);
                if (id == UINT32_MAX) return;
                if (keep) hits.emplace_back(c, id);
                ++seenM[id];
            }, model.markerDenom());
        }
        for (const auto& h : hits) {
            if (seenM[h.second] == 1) contigMarkers[h.first].push_back(h.second);
        }
        for (auto& v : contigMarkers) std::sort(v.begin(), v.end());
    }
    // marker -> panel plasmids carrying it, built once over the plasmids that matter here.
    std::unordered_map<uint32_t, std::vector<uint32_t>> markerToPlasmid;
    if (model.plasmidCount() > 0) {
        std::unordered_map<uint32_t, char> wanted;
        for (size_t c = 0; c < n; ++c) for (uint32_t m : contigMarkers[c]) wanted[m] = 1;
        for (uint32_t pi = 0; pi < model.plasmids().size(); ++pi) {
            for (uint32_t m : model.plasmids()[pi].markers) {
                if (wanted.count(m)) markerToPlasmid[m].push_back(pi);
            }
        }
    }
    // Both minima are tunable for the same reason the threshold is: the useful values are
    // set by how much sequence a contig needs before its markers mean anything, which is a
    // property of this panel at this sampling rate and is measured, not derived.
    auto envSize = [](const char* name, size_t dflt) -> size_t {
        const char* e = std::getenv(name);
        if (!e) return dflt;
        const long v = std::atol(e);
        return v > 0 ? static_cast<size_t>(v) : dflt;
    };
    const size_t weakMin = envSize("TESSERACT_COMEMBER_WEAK", kMinMarkersWeak);
    const size_t strongMin = envSize("TESSERACT_COMEMBER_STRONG", kMinMarkersStrong);

    auto coMembership = [&](size_t a, size_t b) -> double {
        const std::vector<uint32_t>& ma = contigMarkers[a];
        const std::vector<uint32_t>& mb = contigMarkers[b];
        // Both sides needed three markers, and that requirement was the whole failure.
        // Markers are sampled at one 31-mer in 512, so three of them want about 1.5 kb of
        // panel-matching sequence. Measured on 40 fold-0 isolates: contigs that got a group
        // have median length 9,559 bp, contigs left ungrouped have median 1,001 bp with 70%
        // under 1.5 kb. The short contigs were not being rejected, they were never eligible,
        // and they are most of what a plasmid is missing -- only 1 of 63 multi-contig
        // plasmids ended up wholly inside one group.
        //
        // The asymmetric rule keeps the evidence requirement where it can be met. One side
        // must be marker-rich enough to identify a panel plasmid; the other only has to be
        // consistent with it. A single marker makes that side's fraction 0 or 1, so the
        // threshold then demands the marker be present and the rich side match strongly.
        // What keeps this from joining everything through mobile elements is upstream: only
        // markers occurring exactly once in this assembly are interned, which excludes the
        // IS copies that are the usual way two unrelated replicons look related.
        const size_t lo = std::min(ma.size(), mb.size());
        const size_t hi = std::max(ma.size(), mb.size());
        if (lo < weakMin || hi < strongMin) return 0.0;
        std::unordered_map<uint32_t, uint32_t> hitA;
        for (uint32_t m : ma) {
            auto it = markerToPlasmid.find(m);
            if (it == markerToPlasmid.end()) continue;
            for (uint32_t pi : it->second) ++hitA[pi];
        }
        double best = 0;
        for (const auto& kv : hitA) {
            const std::vector<uint32_t>& pm = model.plasmids()[kv.first].markers;
            size_t nb = 0;
            for (uint32_t m : mb) {
                if (std::binary_search(pm.begin(), pm.end(), m)) ++nb;
            }
            if (nb == 0) continue;
            const double score = std::min(static_cast<double>(kv.second) / ma.size(),
                                          static_cast<double>(nb) / mb.size());
            if (score > best) best = score;
        }
        return best;
    };

    const bool haveMembership = model.plasmidCount() > 0;
    double panelThreshold = 0.5;
    if (const char* e = std::getenv("TESSERACT_COMEMBER_MIN")) {
        const double v = std::atof(e);
        if (v > 0 && v <= 1.0) panelThreshold = v;
    }
    if ((links && links->usable()) || haveMembership) {
        for (size_t a = 0; a < n; ++a) {
            if (out.calls[a].cls != RepliconClass::Plasmid) continue;
            for (size_t b = a + 1; b < n; ++b) {
                if (out.calls[b].cls != RepliconClass::Plasmid) continue;
                // Either signal may carry the pair; both are checked, and co-membership
                // can act where read pairs are silent, which at these distances is most
                // of the time.
                const bool byPairs =
                    links && links->usable() &&
                    links->linkWeight(static_cast<uint32_t>(a), static_cast<uint32_t>(b))
                        >= kMinCoResident;
                // 0.5 sits between the measured medians: same-plasmid pairs score 0.979,
                // different-plasmid pairs 0.333. It is tunable because the threshold is the
                // one knob that moves grouping along its trade -- lower merges more, which
                // buys completeness with homogeneity -- and the right point is a measured
                // question, not a derivable one.
                const bool byPanel = haveMembership && coMembership(a, b) >= panelThreshold;
                if (!byPairs && !byPanel) continue;
                // Depth must also agree: two contigs of one molecule are present in the
                // same number of copies, so a link between contigs at 8x and 45x is a
                // shared repeat rather than co-residence.
                const double ra = out.calls[a].depthRatio, rb = out.calls[b].depthRatio;
                const double hi = std::max(ra, rb), lo = std::min(ra, rb);
                if (lo <= 0 || hi / lo > 2.0) continue;
                parent[find(static_cast<uint32_t>(a))] = find(static_cast<uint32_t>(b));
            }
        }
    }

    // ---- hub grouping: let the panel plasmid be the shared referent ---------
    // The pairwise rule above needs evidence running directly between two contigs, and that
    // is why groups fragment. Adding read pairs cut the ungrouped pile from 260 contigs to
    // 73 without moving the whole-plasmid rate at all: contigs are being grouped, into
    // several groups per plasmid. Two contigs from opposite ends of one molecule may share
    // no markers with each other while both matching the same panel plasmid strongly, and
    // pairwise scoring cannot see that.
    //
    // So each contig names its own best panel plasmid, and contigs naming the same one are
    // joined. The panel plasmid is a hub, not a claim: it does not assert that the isolate
    // carries that plasmid, only that these contigs answer to the same reference.
    //
    // The obvious hazard is a common backbone -- an IncF replicon that half the panel shares
    // would merge unrelated molecules. Three things hold it back: only markers occurring
    // exactly once in this assembly are eligible, so multi-copy mobile elements never enter;
    // the match must clear the same threshold as the pairwise rule; and the depth gate still
    // applies to every union. Measured before being enabled, and off by default until then.
    //
// Measured across 40 fold-0 isolates against a panel clustered at mash d<0.01, the hub is a
// no-op: whole plasmids delivered in one group 9/65 either way, one contig moved. It stays
// here, behind a knob, because the mechanism demonstrably fires -- collisions went 0 -> 5 on
// the first isolate -- and it is the right place to attach a denser panel. It is not on.
    const char* hubEnv = std::getenv("TESSERACT_HUB_GROUP");
    // Value, not presence. Every other boolean knob in this tree parses its value, so a
    // sweep harness that exports TESSERACT_HUB_GROUP=0 for its off arm would otherwise turn
    // the pass ON in both arms and report a flat curve as a measured null result.
    if (haveMembership && hubEnv && std::atoi(hubEnv) != 0) {
        double hubThreshold = panelThreshold;
        if (const char* e = std::getenv("TESSERACT_HUB_MIN")) {
            const double v = std::atof(e);
            if (v > 0 && v <= 1.0) hubThreshold = v;
        }
        // Optional accession -> cluster table, applied at use time rather than baked into
        // the model. The model already stores each panel plasmid's accession, so collapsing
        // the panel needs no rebuild of a 339 MB file -- and the clustering threshold stays
        // a thing that can be swept without regenerating anything.
        //
        // Without it, hubs cannot form: measured on one isolate, all ten plasmid contigs
        // named a panel plasmid and no two named the same one. Argmax over tens of thousands
        // of near-identical accessions is unique per contig. Clustering is what gives two
        // contigs of one molecule something to agree on.
        std::vector<uint32_t> clusterOf;      // panel index -> cluster id, empty if unused
        if (const char* cf = std::getenv("TESSERACT_PLASMID_CLUSTERS")) {
            std::unordered_map<std::string, uint32_t> byName;
            // Every failure below is announced. Silence here is indistinguishable from
            // "clustering applied and merged nothing", which is precisely the result this
            // experiment is trying to tell apart: a mistyped path would report a clean null.
            if (FILE* fh = std::fopen(cf, "r")) {
                char acc[256];
                unsigned long cid = 0;
                while (std::fscanf(fh, "%255s %lu", acc, &cid) == 2) {
                    byName.emplace(acc, static_cast<uint32_t>(cid));
                }
                const bool truncated = std::feof(fh) == 0;
                std::fclose(fh);
                if (truncated) {
                    std::fprintf(stderr, "warning: %s stopped parsing after %zu entries "
                                         "(malformed line); clustering is partial\n",
                                 cf, byName.size());
                }
            } else {
                std::fprintf(stderr, "warning: cannot open TESSERACT_PLASMID_CLUSTERS=%s; "
                                     "hub grouping runs unclustered\n", cf);
            }
            if (!byName.empty()) {
                clusterOf.assign(model.plasmids().size(), UINT32_MAX);
                size_t mapped = 0;
                for (size_t pi = 0; pi < model.plasmids().size(); ++pi) {
                    auto it = byName.find(model.plasmids()[pi].name);
                    if (it != byName.end()) { clusterOf[pi] = it->second; ++mapped; }
                }
                if (verbose) {
                    std::fprintf(stderr, "      hub: %zu of %zu panel plasmids mapped to "
                                         "clusters from %s\n",
                                 mapped, model.plasmids().size(), cf);
                }
                // A table that names almost nothing in this model is a mismatch, not a
                // configuration -- silently hubbing on the handful that did match would be
                // worse than not hubbing at all.
                if (mapped * 2 < model.plasmids().size()) {
                    std::fprintf(stderr, "warning: only %zu of %zu panel plasmids are in "
                                         "%s; ignoring the table\n",
                                 mapped, model.plasmids().size(), cf);
                    clusterOf.clear();
                }
            }
        }
        // With clusters, argmax runs over clusters; without, over accessions, which is the
        // arrangement already measured to produce no collisions.
        //
        // The two key spaces must not overlap. Cluster ids run 1..k and panel indices run
        // 0..size-1 over the same panel, so a plasmid absent from the table -- allowed, the
        // guard above tolerates up to half of them -- would land on the key of an unrelated
        // cluster with the same number. That is not a near miss: `hubOwner` is keyed on this
        // value, so two contigs answering to different referents would be unioned and emitted
        // under one `_plas_<group>`. Unmapped plasmids are offset past every cluster id.
        uint32_t clusterCeiling = 0;
        for (uint32_t cid : clusterOf) {
            if (cid != UINT32_MAX) clusterCeiling = std::max(clusterCeiling, cid);
        }
        auto keyOf = [&, clusterCeiling](uint32_t pi) -> uint32_t {
            if (clusterOf.empty() || pi >= clusterOf.size()) return pi;
            return clusterOf[pi] == UINT32_MAX ? clusterCeiling + 1 + pi : clusterOf[pi];
        };

        size_t hubNamed = 0, hubCollided = 0;
        std::map<uint32_t, std::vector<size_t>> hubOwner;  // hub key -> contigs naming it
        for (size_t a = 0; a < n; ++a) {
            if (out.calls[a].cls != RepliconClass::Plasmid) continue;
            const std::vector<uint32_t>& ma = contigMarkers[a];
            if (ma.size() < strongMin) continue;
            std::unordered_map<uint32_t, uint32_t> hits;
            // A marker shared by several panel plasmids of one cluster must count once, or
            // a large cluster wins every argmax by size alone.
            std::unordered_map<uint32_t, uint32_t> lastMarker;
            for (uint32_t m : ma) {
                auto it = markerToPlasmid.find(m);
                if (it == markerToPlasmid.end()) continue;
                for (uint32_t pi : it->second) {
                    const uint32_t key = keyOf(pi);
                    auto seen = lastMarker.find(key);
                    if (seen != lastMarker.end() && seen->second == m) continue;
                    lastMarker[key] = m;
                    ++hits[key];
                }
            }
            uint32_t bestPi = UINT32_MAX;
            double bestScore = 0;
            for (const auto& kv : hits) {
                const double score = static_cast<double>(kv.second) / ma.size();
                bool better = score > bestScore;
                if (!better && score == bestScore && bestPi != UINT32_MAX) {
                    // Unclustered, the key is a panel index and ties go to the shorter
                    // plasmid: a marker set that fills a small replicon says more than the
                    // same set inside a large one. Clustered, the key is a cluster id and
                    // indexes nothing, so ties go to the lower id purely to stay
                    // deterministic -- unordered_map iteration order is not.
                    better = clusterOf.empty()
                                 ? model.plasmids()[kv.first].length <
                                       model.plasmids()[bestPi].length
                                 : kv.first < bestPi;
                }
                if (better) { bestScore = score; bestPi = kv.first; }
            }
            // The hub needs its own threshold, and not for tuning. At the pairwise threshold
            // it is a no-op by construction: the pairwise score is min(fracA, fracB) on a
            // shared plasmid, so if two contigs each put half their markers on P they have
            // already been joined pairwise. A hub can only add something below that line,
            // where one contig matches P strongly and the other weakly.
            if (bestPi == UINT32_MAX || bestScore < hubThreshold) continue;
            ++hubNamed;
            hubOwner[bestPi].push_back(a);
        }
        // Union within each hub, in depth order. Comparing every contig against whichever
        // one happened to claim the key first threw away the pairs the hub exists to find:
        // a single-copy fragment at 1x taking ownership of a cluster blocked two 3x contigs
        // of one genuine multi-copy plasmid from ever being compared to each other, and the
        // result depended on contig order. Sorting by depth and anchoring each block on its
        // own first member keeps every union inside one 2x window without that dependence.
        for (auto& kv : hubOwner) {
            std::vector<size_t>& members = kv.second;
            if (members.size() < 2) continue;
            std::sort(members.begin(), members.end(), [&](size_t x, size_t y) {
                if (out.calls[x].depthRatio != out.calls[y].depthRatio)
                    return out.calls[x].depthRatio < out.calls[y].depthRatio;
                return x < y;
            });
            size_t anchor = members.front();
            for (size_t idx = 1; idx < members.size(); ++idx) {
                const size_t m = members[idx];
                const double lo = out.calls[anchor].depthRatio;
                const double hi = out.calls[m].depthRatio;
                if (lo <= 0 || hi / lo > 2.0) { anchor = m; continue; }
                ++hubCollided;
                parent[find(static_cast<uint32_t>(m))] = find(static_cast<uint32_t>(anchor));
            }
        }
        if (verbose) {
            std::fprintf(stderr, "      hub: %zu contigs named a panel plasmid, %zu named "
                                 "one another contig had already named\n",
                         hubNamed, hubCollided);
        }
    }
    // A group number is a claim that these contigs belong together, so a contig that was
    // grouped with nothing must not carry one. Numbering every plasmid contig -- which is
    // what this did -- produced output with no bare `_plas` anywhere and a mean of 17.4
    // "plasmid groups" per clinical isolate, a figure that reads as seventeen plasmids and
    // is mostly singletons. Singletons now say `_plas` and nothing more: plasmid, molecule
    // unknown, which is the truth about them.
    std::vector<uint32_t> rootSize(n, 0);
    for (size_t i = 0; i < n; ++i) {
        if (out.calls[i].cls == RepliconClass::Plasmid) ++rootSize[find(static_cast<uint32_t>(i))];
    }
    std::unordered_map<uint32_t, uint32_t> groupOf;
    for (size_t i = 0; i < n; ++i) {
        RepliconCall& call = out.calls[i];
        if (call.cls == RepliconClass::Plasmid) {
            const uint32_t root = find(static_cast<uint32_t>(i));
            if (rootSize[root] >= 2) {
                auto it = groupOf.find(root);
                if (it == groupOf.end()) {
                    const uint32_t g = static_cast<uint32_t>(groupOf.size()) + 1;
                    groupOf.emplace(root, g);
                    call.group = g;
                } else {
                    call.group = it->second;
                }
            }
            ++out.plasmidContigs;
        } else if (call.cls == RepliconClass::Chromosome) {
            ++out.chromosomeContigs;
        } else {
            ++out.unassignedContigs;
        }
        if (call.circular) ++out.circularContigs;
    }
    out.plasmidGroups = groupOf.size();

    if (verbose) {
        std::fprintf(stderr,
                     "      replicons: %zu chromosomal, %zu plasmid contigs in %zu groups, "
                     "%zu unassigned, %zu circular (modal depth %.1fx)\n",
                     out.chromosomeContigs, out.plasmidContigs, out.plasmidGroups,
                     out.unassignedContigs, out.circularContigs, out.modalDepth);
    }
    return out;
}

}  // namespace ts
