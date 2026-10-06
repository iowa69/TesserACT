// Phase 2, emit-B, feature F10: the provenance manifest (report.json p2.provenance).
//
// What it records (EVAL_PLAN_P2 L-F10a):
//   * every input read file: path as given, real path, size, mtime, md5 of the first MiB (the armx
//     manifest key's md5_head1m) and the SHA-256 of the whole file;
//   * the running binary (/proc/self/exe): real path, size, md5 and SHA-256;
//   * the model and every other file named on the command line (--model/--organism model, --qc,
//     --is-panel, --is-sites) or by a TESSERACT_* text flag (sketches, sidecars): md5 and SHA-256;
//   * every TESSERACT_* variable of the environment (value, parsed value, kind), the full argv, the
//     [defaults] line, the requested and resolved thread counts, mode and k ladder;
//   * versions: TesserACT, the git commit the binary was built from (when the build knew it), the
//     compiler, zlib and the C++ standard.
// Hashing runs on a background thread started when the run starts, so the SHA-256 of a multi-GB read
// set overlaps the assembly instead of adding to it. The thread only reads files.
#pragma once

#include <cstdint>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace ts {
namespace p2 {

struct FileDigest {
    std::string role, path, realpath, md5Head1m, md5, sha256, error;
    uint64_t size = 0;
    long long mtime = 0;
};

struct ProvenanceInput {
    std::vector<std::pair<std::string, std::string>> inputs;      // (role, path): -1, -2, --12, -s
    std::vector<std::pair<std::string, std::string>> namedFiles;  // (role, path)
    std::string defaultsLine;
    std::string mode;
};

class ProvenanceJob {
public:
    ProvenanceJob() = default;
    ProvenanceJob(const ProvenanceJob&) = delete;
    ProvenanceJob& operator=(const ProvenanceJob&) = delete;
    ~ProvenanceJob();
    void start(const ProvenanceInput& in);
    bool started() const { return started_; }
    // Joins the hashing thread and returns the JSON object.
    std::string finishJson(int threadsResolved, const std::vector<int>& kLadder);
    size_t filesHashed() const { return filesHashed_; }
    uint64_t bytesHashed() const { return bytesHashed_; }
    double seconds() const { return seconds_; }

private:
    void work();
    ProvenanceInput in_;
    std::thread t_;
    bool started_ = false;
    std::vector<FileDigest> inputs_, named_, env_;
    FileDigest binary_;
    size_t filesHashed_ = 0;
    uint64_t bytesHashed_ = 0;
    double seconds_ = 0;
};

FileDigest digestFile(const std::string& role, const std::string& path, bool head1m);
std::vector<std::string> processArgv();
std::string gitCommit();

}  // namespace p2
}  // namespace ts
