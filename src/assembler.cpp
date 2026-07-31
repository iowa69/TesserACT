#include "assembler.h"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <thread>

#include "correct.h"
#include "counter.h"
#include "polish.h"
#include "resolve.h"
#include "util.h"

namespace ts {

Assembler::Assembler(AssemblyOptions opt) : opt_(std::move(opt)) {
    if (opt_.threads <= 0) {
        unsigned hc = std::thread::hardware_concurrency();
        opt_.threads = hc ? static_cast<int>(hc) : 4;
    }
}

std::vector<int> Assembler::resolveKLadder() const {
    if (!opt_.kValues.empty()) return opt_.kValues;
    const int rl = static_cast<int>(reads_.maxReadLength());
    // Each step must stay comfortably below the read length or too few k-mers
    // per read survive to keep the graph connected.
    if (rl >= 200) return {21, 33, 55, 77, 95};
    if (rl >= 140) return {21, 33, 55, 77, 95};
    if (rl >= 100) return {21, 33, 45, 55};
    if (rl >= 70)  return {21, 33, 45};
    if (rl >= 45)  return {17, 25, 31};
    return {15, 21};
}

bool Assembler::iterate(int k, const std::vector<std::string>& carryOver,
                        UnitigGraph& graph, double& meanCoverage, std::string& error) {
    util::Timer t;
    KmerCounter counter(k, opt_.threads);

    // Contigs from the previous, smaller k are trusted sequence: weighting them
    // at the coverage peak keeps them above the abundance cutoff without
    // letting them dominate the histogram used to choose that cutoff.
    const uint32_t carryWeight = carryOver.empty() ? 0 : 4;
    counter.count(reads_, carryOver, carryWeight);

    KmerTable solid;
    counter.extractSolid(opt_.forcedCutoff, solid);
    const CountingStats& cs = counter.stats();

    if (opt_.verbose) {
        std::fprintf(stderr, "  k=%-3d count %.1fs  distinct=%s solid=%s cutoff=%u peak=%.1f\n",
                     k, t.elapsed(), util::commify(static_cast<long long>(cs.distinctKmers)).c_str(),
                     util::commify(static_cast<long long>(cs.solidKmers)).c_str(),
                     cs.cutoff, cs.peakCoverage);
    }
    if (solid.size() == 0) {
        error = "no solid k-mers at k=" + std::to_string(k) +
                "; coverage may be too low or the cutoff too aggressive";
        return false;
    }

    t.reset();
    graph = UnitigGraph::build(solid, k, opt_.threads);
    graph.compact();
    if (opt_.verbose) {
        std::fprintf(stderr, "  k=%-3d graph %.1fs  ", k, t.elapsed());
        graph.stats("built", true);
    }

    meanCoverage = cs.peakCoverage > 0 ? cs.peakCoverage : graph.medianCoverage();
    t.reset();
    graph.simplify(meanCoverage, static_cast<int>(reads_.maxReadLength()), opt_.verbose,
                   opt_.bubbleCoverageLimit);
    if (opt_.verbose) {
        std::fprintf(stderr, "  k=%-3d simplify %.1fs  ", k, t.elapsed());
        graph.stats("simplified", true);
    }
    return true;
}

bool Assembler::run(std::string& error) {
    util::Timer total;

    if (opt_.verbose) std::fprintf(stderr, "[1/5] loading reads\n");
    if (!reads_.load(opt_.libraries, opt_.threads, error)) return false;
    if (reads_.size() == 0) { error = "no reads were loaded"; return false; }
    if (opt_.verbose) {
        std::fprintf(stderr, "      %s reads, %s bases, max length %u%s\n",
                     util::commify(static_cast<long long>(reads_.size())).c_str(),
                     util::commify(static_cast<long long>(reads_.totalBases())).c_str(),
                     reads_.maxReadLength(), reads_.paired() ? ", paired" : "");
    }

    const std::vector<int> ladder = resolveKLadder();
    stats_.kUsed = ladder;

    if (opt_.correctReads && !ladder.empty()) {
        // Correcting against the smallest k in the ladder: short k-mers stay
        // above the abundance cutoff even in thin coverage, so the trusted set
        // covers the genome densely enough to anchor almost every read.
        const int kc = ladder.front();
        if (opt_.verbose) std::fprintf(stderr, "[2/5] read error correction (k=%d)\n", kc);
        util::Timer t;
        KmerCounter cc(kc, opt_.threads);
        cc.count(reads_, {}, 0);
        KmerTable trusted;
        cc.extractSolid(opt_.forcedCutoff, trusted);
        const CorrectionStats rs = correctReads(reads_, trusted, kc, opt_.threads);
        if (opt_.verbose) {
            std::fprintf(stderr,
                         "      %s bases corrected in %s reads (%s had no trusted anchor), %.1fs\n",
                         util::commify(static_cast<long long>(rs.basesCorrected)).c_str(),
                         util::commify(static_cast<long long>(rs.readsCorrected)).c_str(),
                         util::commify(static_cast<long long>(rs.readsUncorrectable)).c_str(),
                         t.elapsed());
        }
    }

    UnitigGraph graph;
    std::vector<std::string> carryOver;
    double meanCoverage = 0;

    if (opt_.verbose) std::fprintf(stderr, "[3/5] de Bruijn iterations\n");
    for (size_t i = 0; i < ladder.size(); ++i) {
        const int k = ladder[i];
        if (k >= static_cast<int>(reads_.maxReadLength())) {
            if (opt_.verbose) {
                std::fprintf(stderr, "  k=%-3d skipped (not shorter than the reads)\n", k);
            }
            continue;
        }
        if (!iterate(k, carryOver, graph, meanCoverage, error)) return false;

        if (i + 1 < ladder.size()) {
            carryOver.clear();
            std::vector<double> covs;
            // Only carry sequence long enough to be meaningful at the next k.
            graph.toContigs(static_cast<size_t>(ladder[i + 1]) + 1, carryOver, covs);
        }
    }

    const int finalK = graph.k();
    size_t minLen = opt_.minContigLen ? opt_.minContigLen : static_cast<size_t>(2 * finalK);

    std::vector<std::string> seqs;
    std::vector<double> covs;

    if (opt_.resolveRepeats && reads_.paired()) {
        if (opt_.verbose) std::fprintf(stderr, "[4/5] paired-end repeat resolution\n");
        util::Timer t;
        PairedResolver resolver(graph, reads_, opt_.threads, opt_.minLinkSupport, opt_.tieRatio);
        resolver.setScaffolding(opt_.scaffold);
        resolver.buildSupport();
        const ResolveStats& rs = resolver.stats();
        if (opt_.verbose) {
            if (rs.insert.usable) {
                std::fprintf(stderr, "      insert size %.0f +/- %.0f from %s pairs\n",
                             rs.insert.mean, rs.insert.stddev,
                             util::commify(static_cast<long long>(rs.insert.observations)).c_str());
            } else {
                std::fprintf(stderr, "      insert size could not be estimated; "
                                     "using graph topology only\n");
            }
            std::fprintf(stderr, "      %s reads anchored, %s linking pairs, %s distinct links\n",
                         util::commify(static_cast<long long>(rs.readsMapped)).c_str(),
                         util::commify(static_cast<long long>(rs.pairsLinking)).c_str(),
                         util::commify(static_cast<long long>(rs.distinctLinks)).c_str());
        }
        resolver.resolve(seqs, covs);
        if (opt_.verbose) {
            std::fprintf(stderr, "      %s paths, %s joins, %.1fs\n",
                         util::commify(static_cast<long long>(resolver.stats().pathsBuilt)).c_str(),
                         util::commify(static_cast<long long>(resolver.stats().unitigsJoined)).c_str(),
                         t.elapsed());
            if (resolver.stats().scaffoldJoins) {
                std::fprintf(stderr, "      %s scaffold joins spanning %s gap bases\n",
                             util::commify(static_cast<long long>(resolver.stats().scaffoldJoins)).c_str(),
                             util::commify(static_cast<long long>(resolver.stats().gapBases)).c_str());
            }
        }
        // Drop anything shorter than the reporting threshold.
        std::vector<std::string> keptSeqs;
        std::vector<double> keptCovs;
        for (size_t i = 0; i < seqs.size(); ++i) {
            if (seqs[i].size() >= minLen) {
                keptSeqs.push_back(std::move(seqs[i]));
                keptCovs.push_back(covs[i]);
            }
        }
        seqs.swap(keptSeqs);
        covs.swap(keptCovs);
    } else {
        if (opt_.verbose) std::fprintf(stderr, "[4/5] collecting contigs\n");
        graph.toContigs(minLen, seqs, covs);
    }

    if (opt_.polish && !seqs.empty()) {
        if (opt_.verbose) std::fprintf(stderr, "[5/5] consensus polishing\n");
        util::Timer t;
        const PolishStats ps = polishContigs(seqs, reads_, opt_.threads, 31, 15, 0.90);
        if (opt_.verbose) {
            std::fprintf(stderr, "      %s reads used, mean depth %.1fx, %s bases corrected, %.1fs\n",
                         util::commify(static_cast<long long>(ps.readsUsed)).c_str(),
                         ps.meanDepth,
                         util::commify(static_cast<long long>(ps.basesChanged)).c_str(),
                         t.elapsed());
        }
    }

    // Longest first, which is both conventional and what most downstream tools
    // assume when they take the top N contigs.
    std::vector<size_t> order(seqs.size());
    for (size_t i = 0; i < order.size(); ++i) order[i] = i;
    std::sort(order.begin(), order.end(),
              [&](size_t a, size_t b) { return seqs[a].size() > seqs[b].size(); });

    std::vector<std::string> outSeqs, outNames;
    outSeqs.reserve(order.size());
    outNames.reserve(order.size());
    size_t totalLen = 0, gc = 0;
    double covWeighted = 0;
    for (size_t rank = 0; rank < order.size(); ++rank) {
        const size_t i = order[rank];
        char name[128];
        std::snprintf(name, sizeof(name), "NODE_%zu_length_%zu_cov_%.4f",
                      rank + 1, seqs[i].size(), covs[i]);
        outNames.emplace_back(name);
        totalLen += seqs[i].size();
        covWeighted += covs[i] * static_cast<double>(seqs[i].size());
        for (char c : seqs[i]) if (c == 'G' || c == 'C') ++gc;
        outSeqs.push_back(std::move(seqs[i]));
    }

    if (!util::makeDirs(opt_.outDir)) {
        error = "cannot create output directory " + opt_.outDir;
        return false;
    }
    const std::string contigPath = opt_.outDir + "/contigs.fasta";
    if (!writeFasta(contigPath, outSeqs, outNames, 80, error)) return false;

    stats_.contigs = outSeqs.size();
    stats_.totalLength = totalLen;
    stats_.largest = outSeqs.empty() ? 0 : outSeqs.front().size();
    stats_.meanCoverage = totalLen ? covWeighted / static_cast<double>(totalLen) : 0;
    stats_.gcPercent = totalLen ? 100.0 * static_cast<double>(gc) / static_cast<double>(totalLen) : 0;
    stats_.seconds = total.elapsed();
    {
        std::vector<size_t> lens;
        lens.reserve(outSeqs.size());
        for (const std::string& s : outSeqs) lens.push_back(s.size());
        std::sort(lens.begin(), lens.end(), std::greater<size_t>());
        size_t acc = 0;
        for (size_t l : lens) {
            acc += l;
            if (acc * 2 >= totalLen) { stats_.n50 = l; break; }
        }
    }

    if (opt_.verbose) std::fprintf(stderr, "      wrote %s\n", contigPath.c_str());
    return true;
}

}  // namespace ts
