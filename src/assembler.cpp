#include "assembler.h"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <ctime>
#include <thread>
#include <unordered_map>
#include <utility>

#include "correct.h"
#include "counter.h"
#include "gapfill.h"
#include "mappolish.h"
#include "organism_join.h"
#include "pairends.h"
#include "replicon.h"
#include "gfa.h"
#include "polish.h"
#include "resolve.h"
#include "util.h"
#include "version.h"

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

// Drop rungs off the TOP of the ladder whose k-mer coverage will be too thin to assemble.
//
// k-mer coverage falls as depth*(L-k+1)/L, so the last rung is always the thinnest, and on
// a 2x151 library at k=95 it retains only 37% of read depth. A 21x library therefore
// reaches the top rung with about 8x k-mer coverage -- against which an abundance cutoff of
// 2 and an error rate of ~0.1% per base leave a graph made largely of noise. The contigs
// this assembler reports come from the final rung, so a shredded final rung is reported
// directly as shredded contigs.
//
// SPAdes stops at k=77 for 151 bp reads. That is not timidity: it is the same arithmetic,
// applied with a fixed read-length rule because SPAdes has no depth estimate available when
// it picks its ladder. Having one lets the decision key off the quantity that actually
// governs the outcome rather than a proxy for it -- a thin 2x250 library should stop early
// and a deep 2x150 library should be allowed to climb, and read length alone cannot tell
// those apart.
//
// The floor is not baked in yet: TESSERA_QC_MIN_KCOV defaults to 0, which disables this
// entirely and reproduces the previous ladder exactly. That default is deliberate -- it
// keeps the pre-QC behaviour available as a baseline while the value is swept.
std::vector<int> Assembler::trimLadderToCoverage(std::vector<int> ladder) const {
    static const double kMinKmerCov = [] {
        const char* e = std::getenv("TESSERA_QC_MIN_KCOV");
        return e ? std::atof(e) : 0.0;
    }();
    static const size_t kMinRungs = [] {
        const char* e = std::getenv("TESSERA_QC_MIN_RUNGS");
        return e ? static_cast<size_t>(std::atoi(e)) : 3u;
    }();
    if (kMinKmerCov <= 0 || !qc_.loaded || ladder.size() <= kMinRungs) return ladder;

    size_t keep = ladder.size();
    while (keep > kMinRungs && qc_.expectedKmerCoverage(ladder[keep - 1]) < kMinKmerCov) --keep;
    if (keep == ladder.size()) return ladder;

    if (opt_.verbose) {
        std::fprintf(stderr,
                     "      QC: dropping %zu top rung(s) -- k=%d would carry only %.1fx "
                     "k-mer coverage at %.1fx read depth (floor %.1f)\n",
                     ladder.size() - keep, ladder.back(),
                     qc_.expectedKmerCoverage(ladder.back()), qc_.readDepth, kMinKmerCov);
    }
    ladder.resize(keep);
    return ladder;
}

std::vector<int> Assembler::resolveKLadder() const {
    // An explicit -k is the user's decision and is never second-guessed, not even by a
    // measurement that disagrees with it.
    if (!opt_.kValues.empty()) return opt_.kValues;
    return trimLadderToCoverage(baseKLadder());
}

std::vector<int> Assembler::baseKLadder() const {

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
    // 127 is the top of what four 64-bit words hold and the most separating
    // rung available: every repeat shorter than it collapses at every smaller
    // k. It needs enough read left over to contribute k-mers, which 165 bp
    // gives -- measured on the closed-reference panel, moving this rung from
    // 115 to 127 left contiguity unchanged and improved both genome fraction
    // (98.945 -> 98.988) and per-base accuracy (0.08 -> 0.06 mismatches per
    // 100 kbp).
    // The 165+ ladder steps 77 -> 127, a jump of 50. Every other step in every other
    // ladder here is between 12 and 22, and the 140-164 ladder steps 77 -> 95. That gap is
    // not a free choice: each rung exists to resolve the repeats the previous one could
    // not, and to hand the next rung a graph built from sequence it has already trusted.
    // Across a 50-k gap the carried contigs have to bridge every repeat between 77 and 127
    // unaided, at a rung where k-mer coverage has dropped by another third.
    //
    // Measured over 491 strains with closed references, split by the top rung their read
    // length selects:
    //
    //     top rung   n     median SPAdes/vanilla NGA50 ratio
    //     95        24     0.78   (TesserACT ahead)
    //     127      458     1.07   (TesserACT behind)
    //
    // SPAdes, on the same 2x250 libraries, uses {21,33,55,77,99,127} -- it does climb to
    // 127, but it gets there through 99. The rung count is what differs, not the ceiling.
    //
    // TESSERA_DENSE_LADDER selects how finely the top of the ladder is sampled:
    //
    //     0   {21,33,55,77,127}                  the original five rungs
    //     1   {21,33,55,77,99,127}               matches what SPAdes uses here (default)
    //     2   {21,33,55,77,99,111,119,127}       four rungs above 77
    //     3   {21,33,55,77,87,99,111,119,127}    also halves the 77->99 step (default)
    //
    // Level 3 is the default, on a full-cohort measurement. A 50-strain comparison had
    // said the extra rungs did nothing (NGA50 5 better / 6 worse, p=0.79) and that reading
    // was wrong -- 40 of those 50 strains tied, so the comparison rested on about ten
    // informative pairs. Repeated over 666 strains, level 3 beats level 1:
    //
    //     NGA50            97 better /  49 worse   p=0.0008
    //     contigs         241 better / 145 worse   p=3e-7
    //     genome fraction 277 better / 178 worse   p=7e-8
    //     misassemblies, mismatches, duplication   no significant difference
    //
    // The cost is real but small: nine count-and-simplify passes instead of six. Measured
    // alone on a quiet machine the denser ladder is actually FASTER (0.92x), because each
    // rung inherits a better-consolidated carry-over and therefore builds a smaller graph;
    // under contention that advantage disappears and it runs somewhat slower. Either way
    // the difference is well inside what the contiguity gain is worth.
    static const int kLadderLevel = [] {
        const char* e = std::getenv("TESSERA_DENSE_LADDER");
        return e ? std::atoi(e) : 3;
    }();
    if (rl >= 165) {
        if (kLadderLevel >= 3) return {21, 33, 55, 77, 87, 99, 111, 119, 127};
        if (kLadderLevel >= 2) return {21, 33, 55, 77, 99, 111, 119, 127};
        if (kLadderLevel >= 1) return {21, 33, 55, 77, 99, 127};
        return {21, 33, 55, 77, 127};
    }
    if (rl >= 140) return {21, 33, 55, 77, 95};
    if (rl >= 100) return {21, 33, 45, 55};
    if (rl >= 70)  return {21, 33, 45};
    if (rl >= 45)  return {17, 25, 31};
    return {15, 21};
}

bool Assembler::iterate(int k, const std::vector<std::string>& carryOver, UnitigGraph& graph,
                        double& meanCoverage, KIteration& it, std::string& error,
                        double prevPeak) {
    util::Timer t;
    it.k = k;
    it.carryOverContigs = carryOver.size();

    KmerCounter counter(k, opt_.threads);
    counter.setMemoryLimit(opt_.maxMemoryBytes);
    // Contigs from the previous, smaller k are trusted sequence: weighting them
    // keeps them above the abundance cutoff without letting them dominate the
    // histogram that cutoff is chosen from.
    //
    // The weight has to scale with coverage. k-mer coverage falls as
    // depth*(L-k+1)/L, so the top of the ladder is always the thinnest rung; on a
    // shallow library it is thin enough that the carried contigs land near the
    // abundance cutoff chosen from the same histogram, get dropped, and the
    // contiguity they existed to preserve is lost with them. A weight of 4 is a
    // nudge at 50x and useless at 13x.
    //
    // Measured over 666 Klebsiella isolates with closed references: below ~25x
    // peak coverage a weight of 32 improved NGA50 on 4/4 strains (+4% to +782%,
    // median +337%), and at 179x it was neutral. But at 39.9x it HALVED NGA50
    // (886,314 -> 429,730), so this must be gated, not applied blanket -- a
    // uniform change would buy a 22% minority's gain with the majority's loss.
    //
    // The gate must key off a number the boost cannot contaminate. Using the
    // PREVIOUS k's observed peak is self-defeating: raising the weight inflates
    // that peak (23.0 at k=55 -> 47.0 at k=77 on ERR6293532), which switches the
    // gate back off at the next rung -- exactly where coverage is thinnest.
    //
    // So project from the FIRST rung, which is counted with no carry-over at all
    // and therefore reflects the library alone:
    //
    //     peak(k) ~ peak(k0) * (L - k + 1) / (L - k0 + 1)
    //
    // Checked on ERR6293532 (L=251, peak(21)=25.0): predicts 13.5 at k=127
    // against 13.35 observed.
    static const uint32_t kCarryWeight = [] {
        const char* e = std::getenv("TESSERA_CARRY_WEIGHT");
        return e ? static_cast<uint32_t>(std::atoi(e)) : 4u;
    }();
    // Threshold on the PROJECTED peak at the final rung. Calibrated against the two
    // strains that bracket the decision:
    //
    //   ERR11578724  peak(k=21)=47  -> projects 25.4 at k=127; boosting gives +337%
    //   ERR5056274   peak(k=21)=61  -> projects ~33   at k=127; boosting COSTS -51%
    //
    // 30 sits between them with margin either side. NOTE: this is calibrated on a
    // handful of strains and wants confirming on a held-out set before it is trusted
    // as a universal default -- the separation is clear but the sample is small.
    // Retuned against 77 changed assemblies from a full-corpus run. Summed NGA50 delta
    // by candidate threshold:  20 -> +60k, 24 -> +194k, 25 -> +265k, 28 -> +39k,
    // 30 -> -120k. The previous value of 30 was net NEGATIVE -- harm concentrated in the
    // 25-30 band (3 better / 8 worse) while the win sits at 20-25 (11/6). The whole
    // 20-28 band is positive, so 25 is a broad optimum, not a knife-edge fit.
    //
    // Magnitude in perspective: +265k summed over 45 boosted strains is a refinement.
    // The substantive fix nearby is the peak-detection repair in counter.cpp, worth
    // +202k on a single strain.
    static const double kCarryBoostBelow = [] {
        const char* e = std::getenv("TESSERA_CARRY_BOOST_BELOW");
        return e ? std::atof(e) : 25.0;
    }();
    // DISABLED (set equal to the ordinary weight). Everything above describes how this
    // boost was calibrated and why it looked like a win; it does not survive being
    // measured against the metric that matters.
    //
    // The first evidence against it was a unit test: "spanned repeat resolved to one
    // contig" began failing at the commit that introduced this, returning 17,231 bp of
    // sequence for a 16,600 bp construct, in two contigs instead of one. That is the boost
    // over-weighting carried contigs until they collapse into each other.
    //
    // Confirmed on 50 closed-reference strains. Turning it off, against the same build
    // with it on: NGA50 10 better / 1 worse, misassemblies 7 better / 1 worse, and no
    // measurable cost anywhere -- mismatches 6/5, duplication 5/4, both indistinguishable
    // from noise. A change that improves contiguity AND misassemblies at once is not a
    // trade being made well; it is a defect being removed.
    //
    // Kept as an environment knob rather than deleted, because the coverage-projection
    // machinery it sits on is sound and worth keeping available: set
    // TESSERA_CARRY_WEIGHT_LOW=32 to restore the previous behaviour exactly.
    static const uint32_t kCarryWeightLow = [] {
        const char* e = std::getenv("TESSERA_CARRY_WEIGHT_LOW");
        return e ? static_cast<uint32_t>(std::atoi(e)) : 4u;
    }();
    // Decide ONCE, from the LAST rung, and apply the decision to every rung.
    //
    // Projecting per-rung is too timid: on a library whose final peak is 24, the
    // early rungs project well above the threshold, so the boost fires only at the
    // last k -- far too late to build the contiguity the later rungs inherit. That
    // is measurable: ERR11578724 scored 8,989 with a per-rung gate against 40,819
    // with the weight applied throughout. What matters is whether the assembly will
    // END thin, because the final rung is where the reported contigs come from.
    //
    // When a QC report is available the projection is unnecessary: it measures read depth
    // directly, from a spectrum built at a small k where the coverage mode is unambiguous,
    // with the error mass already separated out. The chain above starts from THIS
    // assembler's peak estimate at the first rung -- the same estimate that was silently
    // landing on the error shoulder until recently, and that the carry boost itself
    // perturbs. Reading depth from outside the loop removes both problems at once.
    static const bool kQcCarry = [] {
        const char* e = std::getenv("TESSERA_QC_CARRY");
        return e && std::atoi(e) != 0;
    }();
    double projected = prevPeak;
    if (kQcCarry && qc_.loaded && finalK_ > 0 && qc_.expectedKmerCoverage(finalK_) > 0) {
        projected = qc_.expectedKmerCoverage(finalK_);
    } else if (basePeak_ > 0 && baseK_ > 0 && finalK_ > 0) {
        const double L = static_cast<double>(reads_.maxReadLength());
        const double num = L - static_cast<double>(finalK_) + 1.0;
        const double den = L - static_cast<double>(baseK_) + 1.0;
        if (num > 0 && den > 0) projected = basePeak_ * num / den;
    }
    const bool thin = projected > 0 && projected < kCarryBoostBelow;
    const uint32_t chosenWeight = thin ? kCarryWeightLow : kCarryWeight;
    const uint32_t carryWeight = carryOver.empty() ? 0 : chosenWeight;
    if (opt_.verbose && carryWeight && thin) {
        std::fprintf(stderr,
                     "  k=%-3d carry weight raised to %u "
                     "(projected peak at final k=%d is %.1f < %.1f)\n",
                     k, carryWeight, finalK_, projected, kCarryBoostBelow);
    }
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
    // First rung is counted without carry-over, so its peak is the honest library
    // coverage and the only safe basis for projecting the thinner rungs.
    if (basePeak_ <= 0 && carryOver.empty() && cs.peakCoverage > 0) {
        basePeak_ = cs.peakCoverage;
        baseK_ = k;
    }
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
    report_.version = ts::kVersion;
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

    // Check the output directory is usable before spending the run on it. This
    // is checked here rather than at the final stage, so a run that cannot
    // write its output fails before it spends the time to produce one.
    if (!util::makeDirs(opt_.outDir)) {
        error = "cannot create output directory " + opt_.outDir;
        return false;
    }
    // The model, panel and site table are not consulted until stage 4b, by which
    // point a real isolate has cost several minutes of read loading, correction
    // and de Bruijn iteration. A path that cannot be opened is knowable now.
    if (opt_.organismModelPath.empty() &&
        (!opt_.isPanelPath.empty() || !opt_.isSitesPath.empty())) {
        error = "--is-panel and --is-sites are only consulted while joining with a model; "
                "add --model FILE, or drop them";
        return false;
    }
    const std::pair<const char*, const std::string*> namedInputs[] = {
        {"--model", &opt_.organismModelPath},
        {"--is-panel", &opt_.isPanelPath},
        {"--is-sites", &opt_.isSitesPath},
    };
    for (const auto& named : namedInputs) {
        if (named.second->empty()) continue;
        std::FILE* probe = std::fopen(named.second->c_str(), "rb");
        if (!probe) {
            error = std::string("cannot open ") + named.first + " file: " + *named.second;
            return false;
        }
        std::fclose(probe);
    }
    {
        const std::string probe = opt_.outDir + "/.tessera_write_test";
        std::FILE* f = std::fopen(probe.c_str(), "w");
        if (!f) {
            error = "output directory is not writable: " + opt_.outDir;
            return false;
        }
        std::fclose(f);
        std::remove(probe.c_str());
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

    // Load the QC report before the ladder is chosen -- choosing the ladder is one of the
    // decisions it informs, and after resolveKLadder() returns it is already too late.
    //
    // A named-but-unusable report is a hard error, not a warning. The whole point of
    // passing one is that the run behaves differently because of it; silently falling back
    // to the self-derived estimates would produce a run that looks QC-guided, is not, and
    // is indistinguishable in its output from one that is.
    if (!opt_.qcPath.empty()) {
        std::string qcError;
        if (!loadLibraryQC(opt_.qcPath, qc_, qcError)) { error = qcError; return false; }
        if (opt_.verbose) {
            std::fprintf(stderr,
                         "      QC: depth %.1fx, genome %.2f Mb, error rate %.4f%%, "
                         "insert %.0f bp%s\n",
                         qc_.readDepth, qc_.genomeSize / 1e6, qc_.errorRate * 100.0,
                         qc_.insertPeak, qc_.spectrumReliable ? "" : "  [spectrum flagged unreliable]");
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
    // The minimum, not the first: -k takes its values in the order given, so
    // testing front() rejected "-k 127,21" while accepting the same set as
    // "-k 21,127", and called 127 the smallest k while doing it.
    const int smallest = ladder.empty() ? 0 : *std::min_element(ladder.begin(), ladder.end());
    if (ladder.empty() || smallest >= static_cast<int>(reads_.maxReadLength())) {
        char buf[512];
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

    // The last rung shorter than the reads is where the reported contigs come from,
    // so it is the rung whose coverage decides the carry weight for all of them.
    finalK_ = 0;
    for (int kk : ladder) {
        if (kk < static_cast<int>(reads_.maxReadLength())) finalK_ = kk;
    }

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
        const double prevPeak = meanCoverage;   // peak from the previous, smaller k
        if (!iterate(k, carryOver, graph, meanCoverage, it, error, prevPeak)) return false;
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
    // Which contigs the layout stage placed on a chromosome track. Declared out
    // here because the classification that consumes it runs much later, after gap
    // closing and polishing have had their turn at the sequences.
    std::vector<char> layoutMembers;

    if (opt_.resolveRepeats && reads_.paired()) {
        if (opt_.verbose) std::fprintf(stderr, "[4/7] paired-end repeat resolution\n");
        util::Timer t;
        PairedResolver resolver(graph, reads_, opt_.threads, opt_.minLinkSupport, opt_.tieRatio,
                                opt_.linkSupportPerX);
        resolver.setScaffolding(opt_.scaffold);
        // OFF by default: measured, and it makes things worse.
        //
        // The reasoning that motivated it was sound as far as it went -- the fitted window
        // really is mean +/- 4 sd on a heavy-tailed sample, really does work out to about
        // [0, 859] against a mean of 262, and really does admit every candidate. Replacing
        // it with percentiles of the QC histogram narrowed it to [49, 459].
        //
        // On 50 held-out difficult strains that lost, decisively: NGA50 2 better / 20 worse
        // (p=0.001) and contig count 0 better / 45 worse (p=5e-9). Not one strain improved
        // its contiguity.
        //
        // The cause is a property of how the QC histogram is built. It measures fragments
        // by overlapping the two mates, so a fragment longer than the two reads combined
        // contributes nothing to it -- 10.8% of pairs here are recorded as "unknown", and
        // those are exactly the longest ones. Its 99th percentile is therefore the 99th
        // percentile of the SHORT subpopulation, and using it as a hard ceiling discards
        // the long-range links that span repeats, which is the whole reason to consult
        // pairs at all. A wide window admits junk; a truncated one rejects the evidence.
        //
        // Kept behind the flag rather than deleted: the lower bound is measured honestly
        // (short fragments all overlap) and may still be worth something on its own.
        static const bool kQcInsert = [] {
            const char* e = std::getenv("TESSERA_QC_INSERT");
            return e && std::atoi(e) != 0;
        }();
        if (kQcInsert && qc_.loaded && qc_.insertUsable) {
            resolver.setInsertBounds(qc_.insertP1, qc_.insertP99);
            if (opt_.verbose) {
                std::fprintf(stderr,
                             "      insert window from QC: [%d, %d] bp "
                             "(peak %.0f, %zu observations)\n",
                             qc_.insertP1, qc_.insertP99, qc_.insertPeak,
                             qc_.insertObservations);
            }
        }
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

    // What is left unjoined at this point is not weakly supported -- it is
    // unsupported. The repeats that break these contigs are longer than a
    // fragment, so no pair could have spoken to them either way. Where a model
    // of the organism is supplied, those junctions get asked of closed genomes
    // of the same species instead of guessed at from coverage.
    // A model named on the command line that cannot be read is fatal, not a
    // warning. It was only ever printed under --verbose, so a --quiet run
    // produced a vanilla assembly while the user believed the model had guided
    // it -- and the two differ in exactly the places the model exists to settle.
    if (!organismModel_.loaded() && !opt_.organismModelPath.empty()) {
        std::string err;
        if (!organismModel_.load(opt_.organismModelPath, err)) {
            error = err;
            return false;
        }
        // --organism names what the reads are; the model names what it was built
        // from. If they disagree the model is the wrong one for these reads, and
        // its joins would be placed from another species' gene order.
        auto lower = [](std::string v) {
            for (char& c : v) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
            return v;
        };
        if (!opt_.organism.empty() && !organismModel_.organism().empty() &&
            lower(opt_.organism) != lower(organismModel_.organism())) {
            error = "--organism says '" + opt_.organism + "' but " + opt_.organismModelPath +
                    " was built for '" + organismModel_.organism() +
                    "'; use a model for this organism, or drop --organism";
            return false;
        }
    }
    if (!opt_.isPanelPath.empty() && !isPanel_.loaded()) {
        std::string err;
        if (!isPanel_.load(opt_.isPanelPath, err)) {
            error = err;
            return false;
        } else if (opt_.verbose) {
            std::fprintf(stderr, "      insertion-sequence panel: %s copies, %s k-mers\n",
                         util::commify(static_cast<long long>(isPanel_.copies())).c_str(),
                         util::commify(static_cast<long long>(isPanel_.size())).c_str());
        }
    }
    if (!opt_.isSitesPath.empty() && isPanel_.siteCount() == 0) {
        std::string err;
        // Fatal for the same reason --model is: a site table the user asked
        // for and tessera could not read changes where contigs are placed.
        if (!isPanel_.loadSites(opt_.isSitesPath, err)) {
            error = err;
            return false;
        } else if (opt_.verbose) {
            std::fprintf(stderr, "      insertion sites: %s recurrent loci\n",
                         util::commify(static_cast<long long>(isPanel_.siteCount())).c_str());
        }
    }
    if (organismModel_.loaded() && !seqs.empty()) {
        if (opt_.verbose) {
            std::fprintf(stderr, "[4b/7] %s model joining (%u genomes, %u plasmid sets)\n",
                         organismModel_.organism().c_str(),
                         organismModel_.genomes(Replicon::Chromosome),
                         organismModel_.genomes(Replicon::Plasmid));
        }
        std::vector<uint32_t> joinSource;
        report_.organism = joinByModel(organismModel_, seqs, covs, finalK, opt_.verbose,
                                      isPanel_.loaded() ? &isPanel_ : nullptr, &joinSource);
        // joinByModel rebuilds the contig vector: chains become one entry and
        // the survivors come back in a different order. gfaPaths is parallel to
        // it and has to be permuted with it, or a P line names one contig and
        // walks another. A contig built by joining has no single walk, so it
        // gets none.
        if (!joinSource.empty()) {
            std::vector<GfaPath> permuted(joinSource.size());
            for (size_t i = 0; i < joinSource.size(); ++i) {
                const uint32_t from = joinSource[i];
                if (from != UINT32_MAX && from < gfaPaths.size()) permuted[i] = gfaPaths[from];
            }
            gfaPaths.swap(permuted);
        }
        report_.organismRun = true;
        report_.organismName = organismModel_.organism();
        report_.organismGenomes = organismModel_.genomes(Replicon::Chromosome);
        report_.organismPlasmids = organismModel_.genomes(Replicon::Plasmid);
        report_.organismExcluded = organismModel_.excluded();

        // Then lay the result out against the panel chromosome it most resembles, when
        // the model carries layout tracks. Adjacency answers one junction at a time and
        // needs to be right about every one of them; over 666 isolates that reached 23%
        // of chromosomes closed at 90% or better, against 95.6% for the layout, with
        // fewer misassemblies and a longer NGA50.
        //
        // Order runs join-then-layout because joining first makes the pieces longer, so
        // each carries more markers to be placed by. What the layout cannot place is kept
        // exactly as it was: in this cohort that residue is 70% plasmid by base, which is
        // the part a clinical read-out actually wants to look at.
        if (organismModel_.trackCount() > 0 && opt_.layout && seqs.size() > 1) {
            report_.layout = layoutByModel(organismModel_, seqs, covs, finalK, opt_.verbose);
            // gfaPaths cannot survive a layout that ran: a scaffold is several contigs with
            // N between them and has no single walk through the graph. Dropping it is
            // correct; permuting it would name one contig and describe another.
            //
            // Only when it ran, though. layoutByModel returns early -- no located markers,
            // no track sharing 200 of them, nothing placeable -- without touching `seqs`,
            // and on those runs the walks are still in exact correspondence. Clearing
            // regardless cost every P record in assembly_graph.gfa, with both --emit-gfa and
            // --layout on by default and nothing said about it: report.json recorded
            // `layout.run = false`, which reads as "the layout changed nothing".
            //
            // The test is `run && scaffolds`, not a change in contig count: a layout in which
            // every scaffold holds one contig reorders without changing the count. `scaffolds`
            // is incremented only inside the flush that builds the new vector, so it is
            // nonzero exactly when that vector was rebuilt.
            if (report_.layout.run && report_.layout.scaffolds > 0) gfaPaths.clear();
            // layoutByModel emits its scaffolds first and appends everything it could not
            // place, so the leading `scaffolds` entries are the chromosome. Recorded here
            // because the later stages preserve contig count and order but the output sort
            // does not.
            layoutMembers.assign(seqs.size(), 0);
            for (size_t i = 0; i < report_.layout.scaffolds && i < layoutMembers.size(); ++i) {
                layoutMembers[i] = 1;
            }
        }
    }

    // Scaffolding wrote the joins it trusts as Ns. Those Ns are a statement
    // about the graph, not about the data: the reads that cross the join are
    // still in memory. Closing them here, before polishing, means the bases
    // that come back are polished with everything else.
    if (opt_.gapFill && !seqs.empty()) {
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

    if (opt_.mapPolisher != Mapper::None) util::makeDirs(opt_.outDir);
    // The k-mer polisher is blind wherever a base is wrong consistently across
    // every read covering it: those wrong k-mers are themselves solid. A real
    // alignment is independent evidence about every position, so it runs last.
    if (opt_.mapPolisher != Mapper::None && !seqs.empty()) {
        if (opt_.verbose) {
            std::fprintf(stderr, "[6b/7] alignment polishing (%s)\n",
                         mapperName(opt_.mapPolisher));
        }
        MapPolishOptions mp;
        mp.mapper = opt_.mapPolisher;
        mp.mapperDir = opt_.mapperDir;
        mp.workDir = opt_.outDir;
        mp.threads = opt_.threads;
        mp.verbose = opt_.verbose;
        for (const Library& lib : opt_.libraries) {
            if (!lib.r1.empty()) mp.reads1.push_back(lib.r1);
            if (!lib.r2.empty()) mp.reads2.push_back(lib.r2);
        }
        report_.mapPolish = mapPolish(seqs, mp);
        // Printed even under --quiet. Unlike a model, a failed polish leaves
        // the contigs correct rather than differently shaped, so it is not
        // fatal -- but the user did ask for it, and silence is how they end up
        // believing it happened.
        if (!report_.mapPolish.error.empty()) {
            std::fprintf(stderr, "      warning: alignment polishing skipped: %s\n",
                         report_.mapPolish.error.c_str());
        }
    }


    // ---- order and name -------------------------------------------------
    // ---- what molecule is each contig on -----------------------------------
    // Runs last, on the finished sequences, because everything before it can still change
    // them: the gap closer rewrites Ns into bases and the polisher corrects them, and a
    // classification made earlier would describe contigs that no longer exist.
    //
    // Read pairs are re-anchored here rather than carried from the resolver. The
    // resolver's evidence is keyed on unitigs, and the mapping from unitigs to final
    // contigs is destroyed by joining -- a joined contig has no single graph walk. Anchoring
    // again costs one pass over the reads and lands every distance directly in the frame
    // the answer is reported in.
    RepliconAssignment replicons;
    {
        std::unique_ptr<ContigEndLinks> links;
        if (reads_.pairCount() > 0 && report_.resolve.insert.usable && seqs.size() > 1) {
            if (opt_.verbose) std::fprintf(stderr, "      anchoring pairs onto final contigs\n");
            links.reset(new ContigEndLinks(seqs, reads_, report_.resolve.insert,
                                           opt_.threads));
        }
        replicons = assignReplicons(seqs, covs, layoutMembers, organismModel_,
                                    links ? links.get() : nullptr, opt_.verbose);
    }

    // ---- order the output as a genome, not as a length ranking ---------------
    // Sorting by length alone scatters the answer: a chromosomal fragment lands between two
    // plasmids, and the contigs of one plasmid land wherever their sizes put them. The file
    // is the deliverable, so it is ordered the way it is read -- chromosome first, then each
    // plasmid molecule whole and contiguous, then what could not be called.
    //
    // Group numbers are renumbered here so that `_plas_1` is the largest molecule in the
    // isolate and file order agrees with the numbering. They were previously handed out in
    // whatever order the contigs happened to be visited, which is stable but arbitrary.
    // Renumbering happens before `report_.replicons` is stored, so the report, the GFA and
    // the FASTA all quote the same number.
    {
        std::unordered_map<uint32_t, size_t> groupBases;
        for (size_t i = 0; i < seqs.size() && i < replicons.calls.size(); ++i) {
            if (replicons.calls[i].group > 0) groupBases[replicons.calls[i].group] += seqs[i].size();
        }
        std::vector<std::pair<size_t, uint32_t>> byBases;   // (-bases, old group)
        byBases.reserve(groupBases.size());
        for (const auto& kv : groupBases) byBases.emplace_back(kv.second, kv.first);
        std::sort(byBases.begin(), byBases.end(), [](const auto& a, const auto& b) {
            if (a.first != b.first) return a.first > b.first;
            return a.second < b.second;   // total order; unordered_map iteration is not one
        });
        std::unordered_map<uint32_t, uint32_t> renumber;
        for (size_t r = 0; r < byBases.size(); ++r) renumber[byBases[r].second] = static_cast<uint32_t>(r + 1);
        for (RepliconCall& rc : replicons.calls) {
            if (rc.group > 0) rc.group = renumber[rc.group];
        }
    }
    report_.replicons = replicons;

    // Rank of the block a contig belongs to. Chromosome, then plasmids, then unassigned;
    // within the plasmids, every grouped molecule (in the new numbering) before the
    // singletons, whose molecule is unknown and which therefore belong to no block.
    const size_t nGroups = replicons.plasmidGroups;
    auto blockOf = [&](size_t i) -> size_t {
        if (i >= replicons.calls.size()) return nGroups + 2;
        const RepliconCall& rc = replicons.calls[i];
        switch (rc.cls) {
            case RepliconClass::Chromosome: return 0;
            case RepliconClass::Plasmid:    return rc.group > 0 ? rc.group : nGroups + 1;
            default:                        return nGroups + 2;
        }
    };

    // Inside the chromosome block, a laid-out contig keeps its layout order. The layout
    // stage places contigs along a panel chromosome and emits its scaffolds in track
    // order, so for those entries the index IS the position on the reconstructed
    // chromosome -- the one ordering in this pipeline that was actually paid for. Sorting
    // them by length would put a 1.9 Mb piece ahead of the 2.4 Mb piece that precedes it
    // on the genome. It changes nothing when the chromosome closes in a single scaffold,
    // which is the common case (49 of 50 clinical isolates), and it is exactly right in
    // the isolates where it did not, which are the ones a reader is looking at.
    //
    // Gated on the chromosome block rather than on `layoutMembers` alone: a laid-out
    // contig shorter than the classifier's floor is left unassigned, and phase must never
    // disagree with the block it sits in.
    auto phaseOf = [&](size_t i) -> int {
        return (blockOf(i) == 0 && i < layoutMembers.size() && layoutMembers[i]) ? 0 : 1;
    };

    std::vector<size_t> order(seqs.size());
    for (size_t i = 0; i < order.size(); ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
        const size_t ba = blockOf(a), bb = blockOf(b);
        if (ba != bb) return ba < bb;
        const int pa = phaseOf(a), pb = phaseOf(b);
        if (pa != pb) return pa < pb;
        if (pa == 0) return a < b;      // track order
        // Length alone is not a total order and std::sort is not stable, so
        // equal-length contigs would be numbered by whatever order they arrived
        // in. The sequence breaks the tie and makes the naming reproducible.
        if (seqs[a].size() != seqs[b].size()) return seqs[a].size() > seqs[b].size();
        return seqs[a] < seqs[b];
    });

    std::vector<std::string> outSeqs, outNames;
    std::vector<GfaPath> outPaths;
    outSeqs.reserve(order.size());
    outNames.reserve(order.size());
    for (size_t rank = 0; rank < order.size(); ++rank) {
        const size_t i = order[rank];
        // The replicon tag goes on the end of the existing name so nothing that parses
        // the old format breaks. `chr` and `plas` are the call; the trailing number on a
        // plasmid is a GROUP, present only when read pairs tied several contigs to the
        // same molecule, and absent rather than invented when they did not. `circular`
        // means this contig's own two ends are joined by pairs -- a finished molecule,
        // which for a plasmid is the difference between a contig and a result.
        char tag[48];
        tag[0] = '\0';
        if (i < replicons.calls.size()) {
            const RepliconCall& rc = replicons.calls[i];
            const char* cls = rc.cls == RepliconClass::Chromosome ? "chr"
                            : rc.cls == RepliconClass::Plasmid    ? "plas"
                                                                  : "unk";
            if (rc.cls == RepliconClass::Plasmid && rc.group > 0) {
                std::snprintf(tag, sizeof(tag), "_%s_%u%s", cls, rc.group,
                              rc.circular ? "_circular" : "");
            } else {
                std::snprintf(tag, sizeof(tag), "_%s%s", cls,
                              rc.circular ? "_circular" : "");
            }
        }
        char name[224];
        std::snprintf(name, sizeof(name), "NODE_%zu_length_%zu_cov_%.4f%s",
                      rank + 1, seqs[i].size(), covs[i], tag);
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
    report_.command = opt_.commandLine;
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
        } else {
            std::fprintf(stderr, "      warning: %s\n", err.c_str());
        }
    }
    if (opt_.emitHtml) {
        std::string err;
        const std::string hp = opt_.outDir + "/report.html";
        if (writeHtmlReport(hp, report_, err)) {
            if (opt_.verbose) std::fprintf(stderr, "      %s\n", hp.c_str());
        } else {
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
