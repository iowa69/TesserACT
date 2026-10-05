// Organism Model 2.0, component C3: which organism are these reads? (DESIGN.md C3 "Detection")
//
// Default OFF. Nothing here runs unless `--organism auto` is given or TESSERACT_OM2_DETECT is
// warn or gate. Port of the om2/diverge/tools/detect_org2.cpp prototype: a model's core markers
// are its sampled 31-mers that are single-copy on >= 90% of its panel chromosomes and never on
// a plasmid set; the score of a model is the fraction of its core markers seen >= 2 times in
// the reads. Every read of every input file is scanned (a file prefix misjudges sorted runs;
// om2/diverge/product.md), with the same splitmix64 1/512 sampling the model builder used.
//
// The sketch holds only the core marker hashes of each model (about 18 k k-mers for the seven
// ESKAPEE models, < 200 KB), and records the md5 of every .tsm it was built from and its hold
// list, so a sketch built from other models than the ones a run loads is detected.
//
// Decision (thresholds from om2/diverge/detect/analysis.txt, MEASURED on 329 held-out isolates
// and 142 isolates of 30 other species):
//   score >= 0.85 accepts a model; second-best >= 0.10 prints a mixture warning;
//   --organism auto picks the best model when it is accepted, else runs without a model;
//   an explicit --organism below 0.85 is refused with exit 2 under gate (unless
//   --organism-force), and only warned about under warn.
#pragma once

#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

namespace ts {
namespace om2 {

constexpr double kDetectAccept = 0.85;
constexpr double kDetectMixture = 0.10;

struct DetectSketchModel {
    std::string name;            // organism name recorded in the .tsm
    std::string tsmMd5;          // md5 of the .tsm the core set was taken from
    std::string holdList;        // hold list the .tsm was built with (path as given), or "-"
    std::string holdMd5;         // its md5, or "-"
    uint32_t genomesChr = 0;
    uint32_t denom = 512;
    std::vector<uint64_t> core;  // canonical 31-mers, sorted
};
struct DetectSketch { std::vector<DetectSketchModel> models; };

struct DetectModelScore { std::string name; uint32_t core = 0, coreHit = 0; double score = 0; };

struct DetectResult {
    bool ran = false;
    std::string mode;            // auto | warn | gate
    std::string sketchPath, sketchMd5;
    std::vector<DetectModelScore> scores;   // best first
    std::string requested;       // --organism as given
    std::string chosen;          // organism the run uses after detection ("" = none)
    double requestedScore = -1;  // score of the requested organism (-1: not in the sketch)
    std::string decision;        // accepted | refused | forced | warned | auto | auto_none
    bool mixture = false;
    std::string modelCheck = "not_checked";  // match | mismatch | absent | not_checked
    uint64_t reads = 0, bases = 0;
    double seconds = 0;
    double best() const { return scores.empty() ? -1 : scores[0].score; }
    double second() const { return scores.size() < 2 ? -1 : scores[1].score; }
};

// The result of the detection this process ran (ran == false when none did). Read by the
// output writer for the [om2-out] line and the report.json om2 block.
DetectResult& lastDetection();

// Core markers of one .tsm (the marker section only; the adjacency tables are never read).
bool sketchFromModel(const std::string& tsmPath, DetectSketchModel& out, std::string& error);
bool writeDetectSketch(const std::string& path, const DetectSketch& s, std::string& error);
bool loadDetectSketch(const std::string& path, DetectSketch& s, std::string& error);

// Scores every model of the sketch against the reads of `files` (FASTQ or FASTA, gzipped or
// not; one thread per file). Fills scores (best first), reads, bases and seconds.
bool scoreReadFiles(const DetectSketch& s, const std::vector<std::string>& files, DetectResult& out,
                    std::string& error);
// The same on in-memory sequences (the component tests).
void scoreSequences(const DetectSketch& s, const std::vector<std::string>& seqs, DetectResult& out);

// Whether this run asks for detection: `--organism auto`, or TESSERACT_OM2_DETECT=warn|gate.
bool detectRequested(const std::string& organism);
// Runs detection for main(): may rewrite `organism` (auto -> the detected name, or "" for no
// model). Returns 0 to continue, or the exit status (2) for a refused or impossible request.
// `modelDir` is where --organism looks for NAME.tsm; the sketch defaults to om2detect.sketch there.
int detectOrganism(std::string& organism, bool force, const std::vector<std::string>& files,
                   const std::string& modelDir, std::FILE* log);
// After --organism resolved to a model file: is it the file the sketch was built from?
// Returns 0 to continue, 2 when gate mode must refuse a mismatch.
int checkSketchModel(const std::string& modelPath, std::FILE* log);

std::string md5OfFile(const std::string& path);

}  // namespace om2
}  // namespace ts
