#include "assembler.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <thread>

#include "correct.h"
#include "counter.h"
#include "gapfill.h"
#include "gfa.h"
#include "polish.h"
#include "resolve.h"
#include "util.h"

namespace ts {

const char* runModeName(RunMode m) {
    switch (m) {
        case RunMode::Fast:       return "fast";
        case RunMode::Standard:   return "standard";
        case RunMode::Careful:    return "careful";
        case RunMode::Aggressive: return "aggressive";
    }
    return "standard";
}

bool parseRunMode(const std::string& s, RunMode& out) {
    if (s == "fast")       { out = RunMode::Fast; return true; }
    if (s == "standard")   { out = RunMode::Standard; return true; }
    if (s == "careful")    { out = RunMode::Careful; return true; }
    if (s == "aggressive") { out = RunMode::Aggressive; return true; }
    return false;
}

void AssemblyOptions::applyMode() {
    switch (mode) {
        case RunMode::Fast:
            if (!userSetRounds) simplifyRounds = 6;
            if (!userSetMinLink) minLinkSupport = 3;
            if (!userSetTie) tieRatio = 1.3;
            if (!userSetPolishPasses) polishPasses = 0;
            break;
        case RunMode::Standard:
            break;
        case RunMode::Careful:
            // More simplification passes and stricter joins: fewer but better
            // supported decisions, at the cost of some contiguity.
            if (!userSetRounds) simplifyRounds = 24;
            if (!userSetMinLink) minLinkSupport = 3;
            if (!userSetLinkPerX) linkSupportPerX = 0.14;
            if (!userSetTie) tieRatio = 1.4;
            if (!userSetPolishPasses) polishPasses = 2;
            break;
        case RunMode::Aggressive:
            // Collapse diverged repeat copies and accept narrower wins: maximum
            // contiguity, with a real risk of a misassembly.
            if (!userSetBubble) bubbleCoverageLimit = 10.0;
            if (!userSetTie) tieRatio = 1.05;
            if (!userSetMinLink) minLinkSupport = 2;
            if (!userSetLinkPerX) linkSupportPerX = 0.05;
            if (!userSetRounds) simplifyRounds = 16;
            break;
    }
    if (polishPasses <= 0) polish = false;
}

Assembler::Assembler(AssemblyOptions opt) : opt_(std::move(opt)) {
    if (opt_.threads <= 0) {
        unsigned hc = std::thread::hardware_concurrency();
        opt_.threads = hc ? static_cast<int>(hc) : 4;
    }
    opt_.applyMode();
    if (opt_.maxMemoryBytes <= 0) {
        // Leave headroom for the graph and the read store, which live alongside
        // the counting tables.
        const long long total = util::totalMemoryBytes();
        opt_.maxMemoryBytes = total > 0 ? static_cast<long long>(total * 0.8) : 0;
    }
}

std::vector<int> Assembler::resolveKLadder() const {
    if (!opt_.kValues.empty()) return opt_.kValues;

    // Mean, not maximum. Many public runs arrive already trimmed by the
    // submitter, so their reads are variable-length: one 301 bp survivor in a
    // library whose reads average 190 would otherwise pick a ladder most of
    // the reads are too short to contribute a single k-mer to. The mean is
    // what determines how much k-mer each read still carries at the top rung.
    const int rl = reads_.size() ? static_cast<int>(reads_.totalBases() / reads_.size())
                                 : static_cast<int>(reads_.maxReadLength());

    // Each rung must stay comfortably below the read length, or too few k-mers
    // per read survive to keep the graph connected.
    if (opt_.mode == RunMode::Fast) {
        if (rl >= 200) return {21, 55, 99};
        if (rl >= 140) return {21, 55, 77};
        if (rl >= 100) return {21, 45, 55};
        if (rl >= 70)  return {21, 45};
        return {15, 21};
    }
    if (opt_.mode == RunMode::Careful) {
        if (rl >= 200) return {21, 33, 45, 55, 67, 77, 87, 99, 111, 127};
        if (rl >= 165) return {21, 33, 45, 55, 67, 77, 87, 99, 111};
        if (rl >= 140) return {21, 33, 45, 55, 67, 77, 87, 95};
        if (rl >= 100) return {21, 31, 41, 51, 61, 71};
        if (rl >= 70)  return {17, 25, 33, 41, 49};
        return {15, 21, 27};
    }
    // 2x250 libraries carry far more k-mer per read than 2x150, so the top
    // rung climbs with the read length rather than stopping at 95. The rung
    // count stays at five: on the closed-reference panel, raising the top from
    // 95 to 127 took unitig N50 from 142,604 to 207,162 and mismatches from
    // 2.78 to 0.52 per 100 kbp, while a sixth rung only bought runtime.
    if (rl >= 200) return {21, 33, 55, 77, 127};
    if (rl >= 165) return {21, 33, 55, 77, 115};
    if (rl >= 140) return {21, 33, 55, 77, 95};
    if (rl >= 100) return {21, 33, 45, 55};
    if (rl >= 70)  return {21, 33, 45};
    if (rl >= 45)  return {17, 25, 31};
    return {15, 21};
}

bool Assembler::iterate(int k, const std::vector<std::string>& carryOver, UnitigGraph& graph,
                        double& meanCoverage, KIteration& it, std::string& error) {
    util::Timer t;
    it.k = k;
    it.carryOverContigs = carryOver.size();

    KmerCounter counter(k, opt_.threads);
    counter.setMemoryLimit(opt_.maxMemoryBytes);
    // Contigs from the previous, smaller k are trusted sequence: weighting them
    // keeps them above the abundance cutoff without letting them dominate the
    // histogram that cutoff is chosen from.
    const uint32_t carryWeight = carryOver.empty() ? 0 : 4;
    counter.count(reads_, carryOver, carryWeight);

    if (counter.exceededMemory()) {
        char buf[512];
        std::snprintf(buf, sizeof(buf),
                      "ran out of the memory budget (%.1f GB) while counting %d-mers. "
                      "Raise it with --max-memory, raise the abundance cutoff with -c, "
                      "drop the largest k with -k, or use --mode fast",
                      static_cast<double>(opt_.maxMemoryBytes) / 1073741824.0, k);
        error = buf;
        return false;
    }

    KmerTable solid;
    counter.extractSolid(opt_.forcedCutoff, solid);
    const CountingStats& cs = counter.stats();
    it.countSeconds = t.elapsed();
    it.totalKmers = cs.totalKmers;
    it.distinctKmers = cs.distinctKmers;
    it.solidKmers = cs.solidKmers;
    it.cutoff = cs.cutoff;
    it.peakCoverage = cs.peakCoverage;
    // The histogram tail is mostly empty; keep enough to show the error
    // shoulder and the coverage peak.
    const size_t histCap = std::min<size_t>(cs.histogram.size(),
                                            static_cast<size_t>(cs.peakCoverage * 4) + 50);
    it.countHistogram.assign(cs.histogram.begin(),
                             cs.histogram.begin() + static_cast<long>(histCap));

    if (opt_.verbose) {
        std::fprintf(stderr, "  k=%-3d count %.1fs  distinct=%s solid=%s cutoff=%u peak=%.1f\n",
                     k, it.countSeconds,
                     util::commify(static_cast<long long>(cs.distinctKmers)).c_str(),
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
    it.graphSeconds = t.elapsed();
    it.unitigsBuilt = graph.liveCount();
    it.lengthBuilt = graph.totalLength();
    it.n50Built = graph.n50();
    if (opt_.verbose) {
        std::fprintf(stderr, "  k=%-3d graph %.1fs  ", k, it.graphSeconds);
        graph.stats("built", true);
    }

    meanCoverage = cs.peakCoverage > 0 ? cs.peakCoverage : graph.medianCoverage();
    t.reset();
    graph.simplify(meanCoverage, static_cast<int>(reads_.maxReadLength()), opt_.verbose,
                   opt_.bubbleCoverageLimit, opt_.simplifyRounds, &it.rounds);
    it.simplifySeconds = t.elapsed();
    it.unitigsFinal = graph.liveCount();
    it.lengthFinal = graph.totalLength();
    it.n50Final = graph.n50();
    it.medianCoverage = graph.medianCoverage();
    if (opt_.verbose) {
        std::fprintf(stderr, "  k=%-3d simplify %.1fs  ", k, it.simplifySeconds);
        graph.stats("simplified", true);
    }
    return true;
}

bool Assembler::run(std::string& error) {
    util::Timer total;
    report_.version = "1.1.0";
    report_.mode = runModeName(opt_.mode);
    report_.threads = opt_.threads;
    {
        std::time_t tt = std::time(nullptr);
        std::tm tmv{};
        localtime_r(&tt, &tmv);
        char buf[64];
        std::strftime(buf, sizeof(buf), "%Y-%m-%d %H:%M:%S", &tmv);
        report_.startedAt = buf;
    }
    for (const Library& l : opt_.libraries) {
        if (!l.r1.empty()) report_.inputFiles.push_back(l.r1);
        if (!l.r2.empty()) report_.inputFiles.push_back(l.r2);
    }

    if (opt_.verbose) std::fprintf(stderr, "[1/7] loading reads\n");
    reads_.setQualityTrim(opt_.qtrim);
    if (!reads_.load(opt_.libraries, opt_.threads, error)) return false;
    if (reads_.size() == 0) { error = "no reads were loaded"; return false; }
    report_.reads = reads_.size();
    report_.inputBases = reads_.totalBases();
    report_.maxReadLength = reads_.maxReadLength();
    report_.paired = reads_.paired();
    report_.qualityTrimmedBases = reads_.trimmedBases();
    if (opt_.verbose) {
        std::fprintf(stderr, "      %s reads, %s bases, max length %u%s  [mode: %s]\n",
                     util::commify(static_cast<long long>(reads_.size())).c_str(),
                     util::commify(static_cast<long long>(reads_.totalBases())).c_str(),
                     reads_.maxReadLength(), reads_.paired() ? ", paired" : "",
                     runModeName(opt_.mode));
        if (reads_.trimmedBases() > 0) {
            const double pct = 100.0 * static_cast<double>(reads_.trimmedBases()) /
                               static_cast<double>(reads_.trimmedBases() + reads_.totalBases());
            std::fprintf(stderr, "      quality trimmed %s bases from 3' ends (%.1f%%, below Q%d)\n",
                         util::commify(static_cast<long long>(reads_.trimmedBases())).c_str(),
                         pct, opt_.qtrim.meanQuality);
        }
    }

    const std::vector<int> ladder = resolveKLadder();
    stats_.kUsed = ladder;

    if (opt_.correctReads && !ladder.empty()) {
        // Correcting against the smallest k: short k-mers stay above the
        // abundance cutoff even in thin coverage, so the trusted set covers the
        // genome densely enough to anchor almost every read.
        const int kc = ladder.front();
        if (opt_.verbose) std::fprintf(stderr, "[2/7] read error correction (k=%d)\n", kc);
        util::Timer t;
        KmerCounter cc(kc, opt_.threads);
        cc.count(reads_, {}, 0);
        KmerTable trusted;
        cc.extractSolid(opt_.forcedCutoff, trusted);
        report_.correction = correctReads(reads_, trusted, kc, opt_.threads);
        report_.correctionSeconds = t.elapsed();
        report_.correctionK = kc;
        report_.correctionRun = true;
        if (opt_.verbose) {
            std::fprintf(stderr,
                         "      %s bases corrected in %s reads (%s had no trusted anchor), %.1fs\n",
                         util::commify(static_cast<long long>(report_.correction.basesCorrected)).c_str(),
                         util::commify(static_cast<long long>(report_.correction.readsCorrected)).c_str(),
                         util::commify(static_cast<long long>(report_.correction.readsUncorrectable)).c_str(),
                         report_.correctionSeconds);
            // Masking removes bases from every later stage, so it has to be
            // visible: a run that silently drops a large share of the data
            // would otherwise look like a clean assembly at low coverage.
            if (report_.correction.basesMasked > 0) {
                std::fprintf(stderr, "      %s bases masked as unvouchable (%.2f%% of input)\n",
                             util::commify(
                                 static_cast<long long>(report_.correction.basesMasked)).c_str(),
                             100.0 * static_cast<double>(report_.correction.basesMasked) /
                                 static_cast<double>(reads_.totalBases()));
            }
        }
    }

    // If no rung is shorter than the longest read, every iteration below is
    // skipped and the run would finish "successfully" with an empty assembly.
    // Quality trimming makes that reachable from ordinary mistakes -- a
    // Phred+64 file read as Phred+33, or a --qtrim-quality above the data's
    // range -- so say what happened instead of writing an empty FASTA.
    if (ladder.empty() || ladder.front() >= static_cast<int>(reads_.maxReadLength())) {
        char buf[512];
        const int smallest = ladder.empty() ? 0 : ladder.front();
        if (reads_.trimmedBases() > 0) {
            std::snprintf(buf, sizeof(buf),
                          "after quality trimming the longest read is %u bp, shorter than the "
                          "smallest k (%d), so there is nothing to assemble. %s bases (%.1f%%) "
                          "were trimmed -- if that seems too many, the file may be Phred+64 "
                          "(tessera assumes Phred+33) or --qtrim-quality may be set too high; "
                          "--no-qtrim disables trimming entirely",
                          reads_.maxReadLength(), smallest,
                          util::commify(static_cast<long long>(reads_.trimmedBases())).c_str(),
                          100.0 * static_cast<double>(reads_.trimmedBases()) /
                              static_cast<double>(reads_.trimmedBases() + reads_.totalBases()));
        } else {
            std::snprintf(buf, sizeof(buf),
                          "the longest read is %u bp, shorter than the smallest k (%d), so there "
                          "is nothing to assemble; lower it with -k",
                          reads_.maxReadLength(), smallest);
        }
        error = buf;
        return false;
    }

    UnitigGraph graph;
    std::vector<std::string> carryOver;
    double meanCoverage = 0;

    if (opt_.verbose) std::fprintf(stderr, "[3/7] de Bruijn iterations\n");
    for (size_t i = 0; i < ladder.size(); ++i) {
        const int k = ladder[i];
        if (k >= static_cast<int>(reads_.maxReadLength())) {
            if (opt_.verbose) {
                std::fprintf(stderr, "  k=%-3d skipped (not shorter than the reads)\n", k);
            }
            continue;
        }
        KIteration it;
        if (!iterate(k, carryOver, graph, meanCoverage, it, error)) return false;
        report_.iterations.push_back(std::move(it));

        if (i + 1 < ladder.size()) {
            carryOver.clear();
            std::vector<double> covs;
            graph.toContigs(static_cast<size_t>(ladder[i + 1]) + 1, carryOver, covs);
        }
    }

    const int finalK = graph.k();
    const size_t minLen = opt_.minContigLen ? opt_.minContigLen : static_cast<size_t>(2 * finalK);

    std::vector<std::string> seqs;
    std::vector<double> covs;
    std::vector<GfaPath> gfaPaths;

    if (opt_.resolveRepeats && reads_.paired()) {
        if (opt_.verbose) std::fprintf(stderr, "[4/7] paired-end repeat resolution\n");
        util::Timer t;
        PairedResolver resolver(graph, reads_, opt_.threads, opt_.minLinkSupport, opt_.tieRatio,
                                opt_.linkSupportPerX);
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
        report_.resolve = resolver.stats();
        report_.resolveSeconds = t.elapsed();
        report_.resolveRun = true;
        report_.insertHistogram = resolver.insertHistogram();

        // Keep the graph walks alongside the sequences so short contigs drop
        // out of both together.
        const std::vector<ResolvedPath>& rp = resolver.paths();
        std::vector<std::string> keptSeqs;
        std::vector<double> keptCovs;
        for (size_t j = 0; j < seqs.size(); ++j) {
            if (seqs[j].size() < minLen) continue;
            keptSeqs.push_back(std::move(seqs[j]));
            keptCovs.push_back(covs[j]);
            GfaPath gp;
            if (j < rp.size()) {
                gp.oriented = rp[j].oriented;
                gp.gaps = rp[j].gaps;
            }
            gfaPaths.push_back(std::move(gp));
        }
        seqs.swap(keptSeqs);
        covs.swap(keptCovs);

        if (opt_.verbose) {
            std::fprintf(stderr, "      %s paths, %s joins, %.1fs\n",
                         util::commify(static_cast<long long>(report_.resolve.pathsBuilt)).c_str(),
                         util::commify(static_cast<long long>(report_.resolve.unitigsJoined)).c_str(),
                         report_.resolveSeconds);
            if (report_.resolve.scaffoldJoins) {
                std::fprintf(stderr, "      %s scaffold joins spanning %s gap bases\n",
                             util::commify(static_cast<long long>(report_.resolve.scaffoldJoins)).c_str(),
                             util::commify(static_cast<long long>(report_.resolve.gapBases)).c_str());
            }
        }
    } else {
        if (opt_.verbose) std::fprintf(stderr, "[4/7] collecting contigs\n");
        graph.toContigs(minLen, seqs, covs);
    }

    // Scaffolding wrote the joins it trusts as Ns. Those Ns are a statement
    // about the graph, not about the data: the reads that cross the join are
    // still in memory. Closing them here, before polishing, means the bases
    // that come back are polished with everything else.
    if (opt_.gapFill && opt_.scaffold && !seqs.empty()) {
        if (opt_.verbose) std::fprintf(stderr, "[5/7] scaffold gap closing\n");
        const GapFillStats gs = closeGaps(seqs, reads_, opt_.threads, 31, 300);
        report_.gapFill = gs;
        report_.gapFillRun = true;
        if (opt_.verbose && gs.gapsSeen) {
            std::fprintf(stderr,
                         "      %s/%s gaps closed (%s ambiguous, %s unspanned), "
                         "%s Ns -> %s bases, %.1fs (thin=%zu budget=%zu seedgone=%zu tgtgone=%zu depth=%.0f floor=%.1f)\n",
                         util::commify(static_cast<long long>(gs.gapsClosed)).c_str(),
                         util::commify(static_cast<long long>(gs.gapsSeen)).c_str(),
                         util::commify(static_cast<long long>(gs.gapsAmbiguous)).c_str(),
                         util::commify(static_cast<long long>(gs.gapsNoPath)).c_str(),
                         util::commify(static_cast<long long>(gs.nBasesRemoved)).c_str(),
                         util::commify(static_cast<long long>(gs.basesInserted)).c_str(),
                         gs.seconds, gs.gapsThinPool, gs.gapsOutOfBudget,
                         gs.seedBelowFloor, gs.targetBelowFloor, gs.meanLocalDepth, gs.meanFloor);
        }
    }

    // How much of the pileup has to back a base before it replaces what the
    // graph produced.
    //
    // This is deliberately close to unanimity, and the consequence is that the
    // stage corrects nothing at all on a clean bacterial isolate -- it covers
    // 99.9998% of assembly positions at 64x and changes zero bases. That looks
    // like a threshold set too high, and it was tested as one. It is not.
    //
    // Lowering it to 0.70 on a panel isolate made 15 corrections and took
    // mismatches against the closed reference from 0.49 to 0.80 per 100 kbp:
    // every correction it bought was wrong. The positions where the pileup
    // disagrees with the contig are collapsed repeats, where the reads come
    // from several copies at once and the majority belongs to whichever copy is
    // commonest -- not to the locus being polished. Near-unanimity is what
    // keeps the stage from rewriting one repeat copy into another.
    //
    // So polishing is a guard, not a corrector: it earns its 3 seconds by
    // confirming the graph consensus already matches the reads, and it should
    // stay silent. Tunable for experiments; do not lower it on this evidence.
    const double kPolishFraction = std::getenv("TESSERA_POLISH_FRACTION")
                                       ? std::atof(std::getenv("TESSERA_POLISH_FRACTION"))
                                       : 0.90;

    if (opt_.polish && !seqs.empty()) {
        if (opt_.verbose) std::fprintf(stderr, "[6/7] consensus polishing\n");
        util::Timer t;
        for (int pass = 0; pass < opt_.polishPasses; ++pass) {
            const PolishStats ps = polishContigs(seqs, reads_, opt_.threads, 31, 15, kPolishFraction);
            report_.polish.readsUsed = ps.readsUsed;
            report_.polish.basesChanged += ps.basesChanged;
            report_.polish.positionsCovered = ps.positionsCovered;
            report_.polish.lowCoveragePositions = ps.lowCoveragePositions;
            report_.polish.meanDepth = ps.meanDepth;
            if (opt_.verbose) {
                std::fprintf(stderr,
                             "      pass %d: %s reads used, mean depth %.1fx, %s bases corrected\n",
                             pass + 1, util::commify(static_cast<long long>(ps.readsUsed)).c_str(),
                             ps.meanDepth,
                             util::commify(static_cast<long long>(ps.basesChanged)).c_str());
            }
            if (ps.basesChanged == 0) break;   // converged
        }
        report_.polishSeconds = t.elapsed();
        report_.polishRun = true;
    }

    // ---- order and name -------------------------------------------------
    std::vector<size_t> order(seqs.size());
    for (size_t i = 0; i < order.size(); ++i) order[i] = i;
    std::sort(order.begin(), order.end(),
              [&](size_t a, size_t b) { return seqs[a].size() > seqs[b].size(); });

    std::vector<std::string> outSeqs, outNames;
    std::vector<GfaPath> outPaths;
    outSeqs.reserve(order.size());
    outNames.reserve(order.size());
    for (size_t rank = 0; rank < order.size(); ++rank) {
        const size_t i = order[rank];
        char name[160];
        std::snprintf(name, sizeof(name), "NODE_%zu_length_%zu_cov_%.4f",
                      rank + 1, seqs[i].size(), covs[i]);
        outNames.emplace_back(name);

        ContigRecord rec;
        rec.length = seqs[i].size();
        rec.coverage = covs[i];
        size_t gc = 0, ns = 0;
        for (char c : seqs[i]) {
            if (c == 'G' || c == 'C') ++gc;
            else if (c == 'N') ++ns;
        }
        const size_t called = seqs[i].size() - ns;
        rec.gcPercent = called ? 100.0 * static_cast<double>(gc) / static_cast<double>(called) : 0;
        rec.gapBases = ns;
        report_.contigs.push_back(rec);

        if (i < gfaPaths.size() && !gfaPaths[i].oriented.empty()) {
            GfaPath gp = gfaPaths[i];
            gp.name = name;
            outPaths.push_back(std::move(gp));
        }
        outSeqs.push_back(std::move(seqs[i]));
    }
    report_.finalize();

    // ---- write outputs --------------------------------------------------
    if (!util::makeDirs(opt_.outDir)) {
        error = "cannot create output directory " + opt_.outDir;
        return false;
    }
    if (opt_.verbose) std::fprintf(stderr, "[7/7] writing output\n");

    const std::string contigPath = opt_.outDir + "/contigs.fasta";
    if (!writeFasta(contigPath, outSeqs, outNames, 80, error)) return false;
    if (opt_.verbose) std::fprintf(stderr, "      %s\n", contigPath.c_str());

    if (opt_.emitUnitigs) {
        std::vector<std::string> us;
        std::vector<double> ucov;
        graph.toContigs(static_cast<size_t>(finalK), us, ucov);
        std::vector<std::string> un;
        un.reserve(us.size());
        for (size_t i = 0; i < us.size(); ++i) {
            char nm[160];
            std::snprintf(nm, sizeof(nm), "UNITIG_%zu_length_%zu_cov_%.4f",
                          i + 1, us[i].size(), ucov[i]);
            un.emplace_back(nm);
        }
        const std::string up = opt_.outDir + "/unitigs.fasta";
        if (!writeFasta(up, us, un, 80, error)) return false;
        if (opt_.verbose) std::fprintf(stderr, "      %s\n", up.c_str());
    }

    if (opt_.emitGfa) {
        const std::string gp = opt_.outDir + "/assembly_graph.gfa";
        if (!writeGfa(gp, graph, outPaths, report_.gfaSegments, report_.gfaLinks, error)) {
            return false;
        }
        if (opt_.verbose) {
            std::fprintf(stderr, "      %s  (%s segments, %s links)\n", gp.c_str(),
                         util::commify(static_cast<long long>(report_.gfaSegments)).c_str(),
                         util::commify(static_cast<long long>(report_.gfaLinks)).c_str());
        }
    }

    stats_.contigs = outSeqs.size();
    stats_.totalLength = report_.totalLength;
    stats_.largest = report_.largest;
    stats_.n50 = report_.n50;
    stats_.gcPercent = report_.gcPercent;
    stats_.meanCoverage = report_.meanCoverage;
    stats_.seconds = total.elapsed();
    report_.totalSeconds = stats_.seconds;
    report_.peakMemoryBytes = static_cast<double>(util::peakMemoryBytes());

    {
        std::string err;
        const std::string jp = opt_.outDir + "/report.json";
        if (writeJsonReport(jp, report_, err)) {
            if (opt_.verbose) std::fprintf(stderr, "      %s\n", jp.c_str());
        } else if (opt_.verbose) {
            std::fprintf(stderr, "      warning: %s\n", err.c_str());
        }
    }
    if (opt_.emitHtml) {
        std::string err;
        const std::string hp = opt_.outDir + "/report.html";
        if (writeHtmlReport(hp, report_, err)) {
            if (opt_.verbose) std::fprintf(stderr, "      %s\n", hp.c_str());
        } else if (opt_.verbose) {
            std::fprintf(stderr, "      warning: %s\n", err.c_str());
        }
    }
    return true;
}

void AssemblyReport::finalize() {
    totalLength = 0;
    gapBases = 0;
    largest = 0;
    double covWeighted = 0, gcWeighted = 0;
    size_t calledTotal = 0;
    for (const ContigRecord& c : contigs) {
        totalLength += c.length;
        gapBases += c.gapBases;
        largest = std::max(largest, c.length);
        covWeighted += c.coverage * static_cast<double>(c.length);
        const size_t called = c.length - c.gapBases;
        gcWeighted += c.gcPercent * static_cast<double>(called);
        calledTotal += called;
    }
    meanCoverage = totalLength ? covWeighted / static_cast<double>(totalLength) : 0;
    gcPercent = calledTotal ? gcWeighted / static_cast<double>(calledTotal) : 0;

    std::vector<size_t> lens;
    lens.reserve(contigs.size());
    for (const ContigRecord& c : contigs) lens.push_back(c.length);
    std::sort(lens.begin(), lens.end(), std::greater<size_t>());
    size_t acc = 0;
    n50 = n90 = l50 = 0;
    for (size_t i = 0; i < lens.size(); ++i) {
        acc += lens[i];
        if (!n50 && acc * 2 >= totalLength) { n50 = lens[i]; l50 = i + 1; }
        if (!n90 && acc * 10 >= totalLength * 9) { n90 = lens[i]; break; }
    }
}

}  // namespace ts
