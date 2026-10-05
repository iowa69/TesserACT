// Component tests for Organism Model 2.0 C3 organism detection (src/organism_detect.*).
//
// Two synthetic "organisms" whose core sets are their sampled 31-mers; reads drawn from one of
// them, from a mixture, and at too low a depth. Checks the score (core markers seen >= 2x),
// the decisions of --organism auto and of TESSERACT_OM2_DETECT=warn|gate, --organism-force,
// the core-marker rule read off a real .tsm, the sketch round trip, FASTQ/FASTA/gz parity with
// in-memory scoring, and the sketch-vs-model md5 check.
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <random>
#include <string>
#include <vector>
#include <unistd.h>
#include <zlib.h>

#include "test_env.h"
#include "organism.h"
#include "organism_detect.h"

using namespace ts;
using namespace ts::om2;

namespace {

int checks = 0, failures = 0;
void check(bool ok, const std::string& what) {
    ++checks;
    if (!ok) { ++failures; std::fprintf(stderr, "FAIL: %s\n", what.c_str()); }
}
std::mt19937_64 rng(77);
std::string randomSeq(size_t n) {
    std::string s(n, 'A');
    for (char& c : s) c = "ACGT"[rng() % 4];
    return s;
}
std::string rc(const std::string& s) {
    std::string r(s.rbegin(), s.rend());
    for (char& c : r) c = c == 'A' ? 'T' : c == 'C' ? 'G' : c == 'G' ? 'C' : 'A';
    return r;
}
std::vector<std::string> reads(const std::string& g, double depth, size_t len = 150) {
    std::vector<std::string> out;
    const size_t n = static_cast<size_t>(depth * g.size() / len);
    for (size_t i = 0; i < n; ++i) {
        const size_t p = rng() % (g.size() - len);
        std::string r = g.substr(p, len);
        out.push_back(rng() % 2 ? rc(r) : r);
    }
    return out;
}
std::vector<uint64_t> sampled(const std::string& g, uint32_t denom) {
    std::vector<uint64_t> v;
    forEachMarkerKmer(g, [&](uint64_t km, uint32_t, int) { v.push_back(km); }, denom);
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
    return v;
}
void writeFastqGz(const std::string& path, const std::vector<std::string>& rs, size_t from, size_t step) {
    gzFile g = gzopen(path.c_str(), "wb");
    for (size_t i = from; i < rs.size(); i += step)
        gzprintf(g, "@r%zu\n%s\n+\n%s\n", i, rs[i].c_str(), std::string(rs[i].size(), 'I').c_str());
    gzclose(g);
}

}  // namespace

int main() {
    testenv::clearTesseractEnv();
    const std::string dir = (std::filesystem::temp_directory_path() / ("tesseract-om2-detect-" + std::to_string(getpid()))).string();
    std::filesystem::create_directories(dir);
    std::FILE* log = std::tmpfile();

    const uint32_t denom = 16;
    const std::string A = randomSeq(150000), B = randomSeq(150000);
    DetectSketch sk;
    sk.models.resize(2);
    sk.models[0].name = "alpha";
    sk.models[0].denom = denom;
    sk.models[0].core = sampled(A, denom);
    sk.models[1].name = "beta";
    sk.models[1].denom = denom;
    sk.models[1].core = sampled(B, denom);
    check(sk.models[0].core.size() > 5000, "enough sampled markers per synthetic organism");

    // in-memory scoring
    const std::vector<std::string> rA = reads(A, 10);
    DetectResult r;
    scoreSequences(sk, rA, r);
    check(r.scores.size() == 2 && r.scores[0].name == "alpha" && r.best() > 0.99 && r.second() < 0.01,
          "reads of alpha: alpha scores > 0.99, beta < 0.01");
    std::vector<std::string> mix = rA;
    for (const std::string& x : reads(B, 2.5)) mix.push_back(x);
    scoreSequences(sk, mix, r);
    check(r.best() > 0.99 && r.second() >= kDetectMixture, "a 25% admixture of beta is seen as a mixture");
    scoreSequences(sk, reads(A, 0.4), r);
    check(r.best() < kDetectAccept, "at 0.4x no model reaches the accept threshold");

    // files: FASTQ.gz pair and FASTA agree exactly with in-memory scoring
    writeFastqGz(dir + "/R1.fq.gz", rA, 0, 2);
    writeFastqGz(dir + "/R2.fq.gz", rA, 1, 2);
    {
        std::ofstream fa(dir + "/reads.fa");
        for (size_t i = 0; i < rA.size(); ++i) fa << ">r" << i << "\n" << rA[i].substr(0, 70) << "\n" << rA[i].substr(70) << "\n";
    }
    DetectResult mem, fq, fa;
    std::string err;
    scoreSequences(sk, rA, mem);
    check(scoreReadFiles(sk, {dir + "/R1.fq.gz", dir + "/R2.fq.gz"}, fq, err), "score FASTQ.gz files: " + err);
    check(scoreReadFiles(sk, {dir + "/reads.fa"}, fa, err), "score a multi-line FASTA file: " + err);
    check(fq.reads == rA.size() && fa.reads == rA.size(), "every read of every file is scanned");
    check(fq.scores.size() == 2 && fq.scores[0].coreHit == mem.scores[0].coreHit &&
              fq.scores[1].coreHit == mem.scores[1].coreHit && fa.scores[0].coreHit == mem.scores[0].coreHit,
          "file scoring equals in-memory scoring exactly");

    // sketch round trip
    check(writeDetectSketch(dir + "/s.sketch", sk, err), "write sketch");
    DetectSketch back;
    check(loadDetectSketch(dir + "/s.sketch", back, err) && back.models.size() == 2 &&
              back.models[0].core == sk.models[0].core && back.models[1].name == "beta" && back.models[1].denom == denom,
          "sketch round trip");
    check(!loadDetectSketch(dir + "/R1.fq.gz", back, err), "a non-sketch file is refused");

    // core rule on a real .tsm: single-copy on >= 90% of chromosomes and on no plasmid set
    {
        OrganismModel m;
        m.beginBuild("gamma", kMarkerK);
        m.setMarkerDenom(denom);
        for (int g = 0; g < 10; ++g) m.noteGenome(Replicon::Chromosome);
        m.noteGenome(Replicon::Plasmid);
        const std::vector<uint64_t> km = sampled(A, denom);
        // marker 0: 10/10 chr, no plasmid -> core; 1: 9/10 -> core; 2: 8/10 -> not; 3: 10/10 + plasmid -> not
        const int chrCount[4] = {10, 9, 8, 10};
        const int plsCount[4] = {0, 0, 0, 1};
        for (int i = 0; i < 4; ++i) {
            const uint32_t id = m.internMarker(km[i]);
            for (int c = 0; c < chrCount[i]; ++c) m.noteMarkerGenome(id, Replicon::Chromosome);
            for (int c = 0; c < plsCount[i]; ++c) m.noteMarkerGenome(id, Replicon::Plasmid);
        }
        m.finalise(1, 1);
        check(m.save(dir + "/gamma.tsm", err), "save a synthetic .tsm: " + err);
        DetectSketchModel sm;
        check(sketchFromModel(dir + "/gamma.tsm", sm, err), "read the core markers of a .tsm: " + err);
        std::vector<uint64_t> want{km[0], km[1]};
        std::sort(want.begin(), want.end());
        check(sm.name == "gamma" && sm.core == want && sm.denom == denom && sm.genomesChr == 10,
              "core = markers on >= 90% of chromosomes and no plasmid");
        check(sm.tsmMd5 == md5OfFile(dir + "/gamma.tsm") && sm.tsmMd5.size() == 32, "the sketch records the .tsm md5");
    }

    // decisions (main()'s entry point)
    setenv("TESSERACT_OM2_DETECT_SKETCH", (dir + "/s.sketch").c_str(), 1);
    const std::vector<std::string> files{dir + "/R1.fq.gz", dir + "/R2.fq.gz"};
    std::string org = "auto";
    check(detectRequested(org), "--organism auto requests detection with the flag unset");
    check(detectOrganism(org, false, files, dir, log) == 0 && org == "alpha" && lastDetection().decision == "auto",
          "--organism auto picks alpha");
    org = "beta";
    check(!detectRequested(org), "an explicit --organism with the flag unset requests nothing");
    setenv("TESSERACT_OM2_DETECT", "gate", 1);
    check(detectRequested(org), "TESSERACT_OM2_DETECT=gate requests detection");
    check(detectOrganism(org, false, files, dir, log) == 2 && lastDetection().decision == "refused",
          "gate refuses --organism beta for alpha reads (exit 2)");
    org = "beta";
    check(detectOrganism(org, true, files, dir, log) == 0 && lastDetection().decision == "forced" && org == "beta",
          "--organism-force runs it anyway, recorded as forced");
    org = "alpha";
    check(detectOrganism(org, false, files, dir, log) == 0 && lastDetection().decision == "accepted" &&
              lastDetection().requestedScore > 0.99,
          "gate accepts --organism alpha");
    // the sketch must come from the model file the run loads
    {
        DetectSketch withMd5 = sk;
        std::ofstream(dir + "/alpha.tsm") << "pretend model";
        withMd5.models[0].tsmMd5 = md5OfFile(dir + "/alpha.tsm");
        writeDetectSketch(dir + "/s.sketch", withMd5, err);
        org = "alpha";
        detectOrganism(org, false, files, dir, log);
        check(checkSketchModel(dir + "/alpha.tsm", log) == 0 && lastDetection().modelCheck == "match", "sketch md5 matches the model");
        std::ofstream(dir + "/other.tsm") << "another model";
        check(checkSketchModel(dir + "/other.tsm", log) == 2 && lastDetection().modelCheck == "mismatch",
              "gate refuses a model the sketch was not built from");
        setenv("TESSERACT_OM2_DETECT", "warn", 1);
        detectOrganism(org, false, files, dir, log);
        check(checkSketchModel(dir + "/other.tsm", log) == 0 && lastDetection().modelCheck == "mismatch",
              "warn only warns about it");
    }
    org = "beta";
    check(detectOrganism(org, false, files, dir, log) == 0 && lastDetection().decision == "warned", "warn: beta is warned about, not refused");
    unsetenv("TESSERACT_OM2_DETECT");
    // low depth: auto finds nothing and assembles without a model
    const std::vector<std::string> thin = reads(A, 0.4);
    writeFastqGz(dir + "/T1.fq.gz", thin, 0, 1);
    org = "auto";
    check(detectOrganism(org, false, {dir + "/T1.fq.gz"}, dir, log) == 0 && org.empty() &&
              lastDetection().decision == "auto_none",
          "--organism auto below the threshold runs without a model");
    // no sketch
    setenv("TESSERACT_OM2_DETECT_SKETCH", (dir + "/missing.sketch").c_str(), 1);
    org = "auto";
    check(detectOrganism(org, false, files, dir, log) == 2, "--organism auto without a sketch is an error");
    setenv("TESSERACT_OM2_DETECT", "warn", 1);
    org = "alpha";
    check(detectOrganism(org, false, files, dir, log) == 0 && lastDetection().decision == "skipped",
          "warn without a sketch warns and continues");
    unsetenv("TESSERACT_OM2_DETECT");
    unsetenv("TESSERACT_OM2_DETECT_SKETCH");
    check(md5OfFile(dir + "/nonexistent").empty(), "md5 of a missing file is empty");

    std::fclose(log);
    std::filesystem::remove_all(dir);
    std::printf("test_om2_detect: %d checks, %d failures\n", checks, failures);
    return failures ? 1 : 0;
}
