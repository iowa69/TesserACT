#include "p2_provenance.h"

#include <zlib.h>

#include <chrono>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sys/stat.h>

#include "envflags.h"
#include "p2_json.h"
#include "p2_sha256.h"
#include "version.h"

#ifndef TS_GIT_COMMIT
#define TS_GIT_COMMIT "unknown"
#endif

namespace ts {
namespace p2 {

namespace {

std::string realPath(const std::string& p) {
    char buf[PATH_MAX];
    if (::realpath(p.c_str(), buf)) return buf;
    return std::string();
}

bool isRegularFile(const std::string& p) {
    struct stat st;
    return ::stat(p.c_str(), &st) == 0 && S_ISREG(st.st_mode);
}

std::string digestJson(const FileDigest& d, bool withHead) {
    JsonObj o;
    if (!d.role.empty()) o.str("role", d.role);
    o.str("path", d.path);
    o.str("realpath", d.realpath);
    o.u("size", d.size);
    if (withHead) {
        o.i("mtime", d.mtime);
        o.str("md5_head1m", d.md5Head1m);
    }
    if (!d.md5.empty()) o.str("md5", d.md5);
    o.str("sha256", d.sha256);
    if (!d.error.empty()) o.str("error", d.error);
    return o.done();
}

}  // namespace

std::string gitCommit() { return TS_GIT_COMMIT; }

FileDigest digestFile(const std::string& role, const std::string& path, bool head1m) {
    FileDigest d;
    d.role = role;
    d.path = path;
    d.realpath = realPath(path);
    struct stat st;
    if (::stat(path.c_str(), &st) == 0) d.mtime = static_cast<long long>(st.st_mtime);
    uint64_t bytes = 0;
    if (!sha256File(path, d.sha256, bytes)) {
        d.error = "cannot read";
        return d;
    }
    d.size = bytes;
    if (head1m) d.md5Head1m = md5FileHead(path, 1u << 20);
    else d.md5 = md5File(path);
    return d;
}

std::vector<std::string> processArgv() {
    std::vector<std::string> out;
    std::ifstream f("/proc/self/cmdline", std::ios::binary);
    if (!f) return out;
    std::string all((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
    size_t at = 0;
    while (at < all.size()) {
        const size_t z = all.find('\0', at);
        const size_t e = z == std::string::npos ? all.size() : z;
        out.push_back(all.substr(at, e - at));
        at = e + 1;
    }
    return out;
}

ProvenanceJob::~ProvenanceJob() {
    if (t_.joinable()) t_.join();
}

void ProvenanceJob::start(const ProvenanceInput& in) {
    in_ = in;
    started_ = true;
    t_ = std::thread([this]() { work(); });
}

void ProvenanceJob::work() {
    const auto t0 = std::chrono::steady_clock::now();
    for (const auto& x : in_.inputs) {
        inputs_.push_back(digestFile(x.first, x.second, true));
        ++filesHashed_;
        bytesHashed_ += inputs_.back().size;
    }
    binary_ = digestFile("binary", "/proc/self/exe", false);
    binary_.path = binary_.realpath;
    ++filesHashed_;
    bytesHashed_ += binary_.size;
    for (const auto& x : in_.namedFiles) {
        if (x.second.empty()) continue;
        named_.push_back(digestFile(x.first, x.second, false));
        ++filesHashed_;
        bytesHashed_ += named_.back().size;
    }
    // every TESSERACT_* text flag whose value names a regular file (detection sketches, sidecars)
    for (const env::SetVariable& v : env::setVariables()) {
        if (v.kind != "text" || v.value.empty() || !isRegularFile(v.value)) continue;
        env_.push_back(digestFile(v.name, v.value, false));
        ++filesHashed_;
        bytesHashed_ += env_.back().size;
    }
    seconds_ = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

std::string ProvenanceJob::finishJson(int threadsResolved, const std::vector<int>& kLadder) {
    if (t_.joinable()) t_.join();
    JsonObj o;
    o.str("hash_method", "sha256 = SHA-256 of every byte of the file; md5_head1m = MD5 of its first 1 MiB "
                         "(the armx manifest key); md5 = MD5 of every byte");
    {
        std::vector<std::string> a;
        for (const FileDigest& d : inputs_) a.push_back(digestJson(d, true));
        o.raw("inputs", jsonArray(a));
    }
    o.raw("binary", digestJson(binary_, false));
    {
        std::vector<std::string> a;
        for (const FileDigest& d : named_) a.push_back(digestJson(d, false));
        o.raw("model_and_named_files", jsonArray(a));
    }
    {
        std::vector<std::string> a;
        for (const FileDigest& d : env_) a.push_back(digestJson(d, false));
        o.raw("flag_files", jsonArray(a));
    }
    {
        std::vector<std::string> a, names;
        for (const env::SetVariable& v : env::setVariables()) {
            JsonObj e;
            e.str("name", v.name).str("value", v.value);
            e.b("registered", v.registered);
            if (v.registered) e.str("kind", v.kind).str("parsed", v.canonical).b("valid", v.valid);
            a.push_back(e.done());
            names.push_back(v.name + "=" + v.value);
        }
        o.raw("tesseract_env", jsonArray(a));
        o.strs("non_default_env", names);
    }
    {
        const std::vector<std::string> argv = processArgv();
        o.strs("argv", argv);
        // the options that change behaviour: everything but the inputs, the output and -t/-q
        std::vector<std::string> opts;
        for (size_t i = 1; i < argv.size(); ++i) {
            const std::string& a = argv[i];
            if (a == "-1" || a == "-2" || a == "--read1" || a == "--read2" || a == "--12" || a == "-s" ||
                a == "--single" || a == "-o" || a == "--out" || a == "-t" || a == "--threads") { ++i; continue; }
            if (a == "-q" || a == "--quiet") continue;
            opts.push_back(a);
        }
        o.strs("non_default_options", opts);
    }
    o.str("defaults_line", in_.defaultsLine);
    o.str("mode", in_.mode);
    {
        std::string k = "[";
        for (size_t i = 0; i < kLadder.size(); ++i) k += (i ? ", " : "") + std::to_string(kLadder[i]);
        o.raw("k_ladder", k + "]");
    }
    o.i("threads", threadsResolved);
    {
        JsonObj v;
        v.str("tesseract", kVersion);
        v.str("git_commit", gitCommit());
#ifdef __VERSION__
        v.str("compiler", __VERSION__);
#endif
        v.str("zlib", zlibVersion());
        v.i("cplusplus", static_cast<long long>(__cplusplus));
        o.raw("versions", v.done());
    }
    o.u("files_hashed", filesHashed_);
    o.u("bytes_hashed", bytesHashed_);
    o.num("hash_seconds", seconds_, 3);
    return o.done();
}

}  // namespace p2
}  // namespace ts
