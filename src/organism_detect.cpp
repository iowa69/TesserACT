// Organism Model 2.0, component C3: organism detection. See organism_detect.h.
#include "organism_detect.h"

#include <zlib.h>

#include <algorithm>
#include <chrono>
#include <cstring>
#include <memory>
#include <mutex>
#include <thread>
#include <unordered_map>

#include "envflags.h"
#include "organism.h"
#include "util.h"

namespace ts {
namespace om2 {

namespace {

// ---- MD5 (RFC 1321), for the provenance fields only ----------------------------------------
class Md5 {
    uint32_t a0 = 0x67452301, b0 = 0xefcdab89, c0 = 0x98badcfe, d0 = 0x10325476;
    uint64_t bytes = 0;
    size_t used = 0;
    unsigned char block[64]{};
    static uint32_t rotl(uint32_t x, int c) { return (x << c) | (x >> (32 - c)); }
    void compress() {
        static const uint32_t K[64] = {
            0xd76aa478, 0xe8c7b756, 0x242070db, 0xc1bdceee, 0xf57c0faf, 0x4787c62a, 0xa8304613, 0xfd469501,
            0x698098d8, 0x8b44f7af, 0xffff5bb1, 0x895cd7be, 0x6b901122, 0xfd987193, 0xa679438e, 0x49b40821,
            0xf61e2562, 0xc040b340, 0x265e5a51, 0xe9b6c7aa, 0xd62f105d, 0x02441453, 0xd8a1e681, 0xe7d3fbc8,
            0x21e1cde6, 0xc33707d6, 0xf4d50d87, 0x455a14ed, 0xa9e3e905, 0xfcefa3f8, 0x676f02d9, 0x8d2a4c8a,
            0xfffa3942, 0x8771f681, 0x6d9d6122, 0xfde5380c, 0xa4beea44, 0x4bdecfa9, 0xf6bb4b60, 0xbebfbc70,
            0x289b7ec6, 0xeaa127fa, 0xd4ef3085, 0x04881d05, 0xd9d4d039, 0xe6db99e5, 0x1fa27cf8, 0xc4ac5665,
            0xf4292244, 0x432aff97, 0xab9423a7, 0xfc93a039, 0x655b59c3, 0x8f0ccc92, 0xffeff47d, 0x85845dd1,
            0x6fa87e4f, 0xfe2ce6e0, 0xa3014314, 0x4e0811a1, 0xf7537e82, 0xbd3af235, 0x2ad7d2bb, 0xeb86d391};
        static const int S[64] = {7, 12, 17, 22, 7, 12, 17, 22, 7, 12, 17, 22, 7, 12, 17, 22,
                                  5, 9,  14, 20, 5, 9,  14, 20, 5, 9,  14, 20, 5, 9,  14, 20,
                                  4, 11, 16, 23, 4, 11, 16, 23, 4, 11, 16, 23, 4, 11, 16, 23,
                                  6, 10, 15, 21, 6, 10, 15, 21, 6, 10, 15, 21, 6, 10, 15, 21};
        uint32_t M[16];
        for (int i = 0; i < 16; ++i)
            M[i] = uint32_t(block[4 * i]) | (uint32_t(block[4 * i + 1]) << 8) | (uint32_t(block[4 * i + 2]) << 16) |
                   (uint32_t(block[4 * i + 3]) << 24);
        uint32_t A = a0, B = b0, C = c0, D = d0;
        for (int i = 0; i < 64; ++i) {
            uint32_t F;
            int g;
            if (i < 16) { F = (B & C) | (~B & D); g = i; }
            else if (i < 32) { F = (D & B) | (~D & C); g = (5 * i + 1) % 16; }
            else if (i < 48) { F = B ^ C ^ D; g = (3 * i + 5) % 16; }
            else { F = C ^ (B | ~D); g = (7 * i) % 16; }
            F = F + A + K[i] + M[g];
            A = D; D = C; C = B;
            B = B + rotl(F, S[i]);
        }
        a0 += A; b0 += B; c0 += C; d0 += D;
    }
public:
    void add(const void* p, size_t n) {
        const auto* s = static_cast<const unsigned char*>(p);
        bytes += n;
        while (n) {
            const size_t take = std::min(n, 64 - used);
            std::memcpy(block + used, s, take);
            used += take; s += take; n -= take;
            if (used == 64) { compress(); used = 0; }
        }
    }
    std::string finish() {
        const uint64_t bits = bytes * 8;
        const unsigned char one = 0x80, zero = 0;
        add(&one, 1);
        while (used != 56) add(&zero, 1);
        unsigned char len[8];
        for (int i = 0; i < 8; ++i) len[i] = static_cast<unsigned char>(bits >> (8 * i));
        add(len, 8);
        char out[33];
        const uint32_t h[4] = {a0, b0, c0, d0};
        for (int i = 0; i < 4; ++i)
            for (int b = 0; b < 4; ++b) std::snprintf(out + 8 * i + 2 * b, 3, "%02x", (h[i] >> (8 * b)) & 0xff);
        return std::string(out, 32);
    }
};

template <typename T>
bool rdPod(std::FILE* f, T& v) { return std::fread(&v, sizeof(T), 1, f) == 1; }
template <typename T>
bool wrPod(std::FILE* f, const T& v) { return std::fwrite(&v, sizeof(T), 1, f) == 1; }
bool rdStr(std::FILE* f, std::string& s, uint32_t cap = 1u << 20) {
    uint32_t n = 0;
    if (!rdPod(f, n) || n > cap) return false;
    s.assign(n, '\0');
    return n == 0 || std::fread(&s[0], 1, n, f) == n;
}
bool wrStr(std::FILE* f, const std::string& s) {
    const uint32_t n = static_cast<uint32_t>(s.size());
    return wrPod(f, n) && (n == 0 || std::fwrite(s.data(), 1, n, f) == n);
}

constexpr char kSketchMagic[8] = {'O', 'M', '2', 'D', 'E', 'T', '0', '1'};

std::string lowerCase(std::string v) {
    for (char& c : v) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return v;
}

// One file's counts of the sketch's k-mers. `keys` maps every sketch k-mer to its slot.
struct FileScan {
    std::vector<uint32_t> counts;
    uint64_t reads = 0, bases = 0;
    std::string error;
};

void countSequence(const std::string& seq, const std::vector<uint64_t>& denoms,
                   const std::unordered_map<uint64_t, uint32_t>& slot, std::vector<uint32_t>& counts) {
    // A k-mer is sampled when its hash falls under the denominator's threshold, so a k-mer
    // sampled at 1/512 is sampled at every model of that density. Models of different density
    // are counted in separate passes (every bundled model is 1/512: one pass).
    for (uint64_t d : denoms) {
        forEachMarkerKmer(seq, [&](uint64_t km, uint32_t, int) {
            auto it = slot.find(km);
            if (it != slot.end()) ++counts[it->second];
        }, d);
    }
}

void scanFile(const std::string& path, const std::vector<uint64_t>& denoms,
              const std::unordered_map<uint64_t, uint32_t>& slot, FileScan& out) {
    out.counts.assign(slot.size(), 0);
    gzFile g = gzopen(path.c_str(), "rb");
    if (!g) { out.error = "cannot open " + path; return; }
    gzbuffer(g, 1 << 20);
    std::vector<char> buf(1 << 16);
    std::string line, seq;
    int format = 0;            // 1 FASTQ, 2 FASTA
    uint64_t lineNo = 0;
    bool more = true;
    while (more) {
        // One whole line, whatever its length.
        line.clear();
        bool got = false;
        while (gzgets(g, buf.data(), static_cast<int>(buf.size()))) {
            got = true;
            line += buf.data();
            if (!line.empty() && line.back() == '\n') break;
        }
        if (!got) { more = false; }
        while (!line.empty() && (line.back() == '\n' || line.back() == '\r')) line.pop_back();
        if (!more) break;
        if (format == 0) {
            if (line.empty()) continue;
            format = line[0] == '>' ? 2 : 1;
        }
        if (format == 1) {
            if ((lineNo & 3) == 1) {
                ++out.reads;
                out.bases += line.size();
                countSequence(line, denoms, slot, out.counts);
            }
            ++lineNo;
        } else {
            if (!line.empty() && line[0] == '>') {
                if (!seq.empty()) { ++out.reads; out.bases += seq.size(); countSequence(seq, denoms, slot, out.counts); }
                seq.clear();
            } else {
                seq += line;
            }
        }
    }
    if (format == 2 && !seq.empty()) { ++out.reads; out.bases += seq.size(); countSequence(seq, denoms, slot, out.counts); }
    int zerr = 0;
    const char* msg = gzerror(g, &zerr);
    if (zerr != Z_OK && zerr != Z_STREAM_END) out.error = path + ": " + (msg ? msg : "read error");
    gzclose(g);
}

void buildSlots(const DetectSketch& s, std::unordered_map<uint64_t, uint32_t>& slot, std::vector<uint64_t>& denoms) {
    for (const DetectSketchModel& m : s.models) {
        if (std::find(denoms.begin(), denoms.end(), static_cast<uint64_t>(m.denom)) == denoms.end())
            denoms.push_back(m.denom);
        for (uint64_t km : m.core) slot.emplace(km, static_cast<uint32_t>(slot.size()));
    }
}

void finishScores(const DetectSketch& s, const std::unordered_map<uint64_t, uint32_t>& slot,
                  const std::vector<uint32_t>& counts, DetectResult& out) {
    out.scores.clear();
    for (const DetectSketchModel& m : s.models) {
        DetectModelScore sc;
        sc.name = m.name;
        sc.core = static_cast<uint32_t>(m.core.size());
        for (uint64_t km : m.core) {
            auto it = slot.find(km);
            if (it != slot.end() && counts[it->second] >= 2) ++sc.coreHit;   // prototype: seen >= 2x
        }
        sc.score = sc.core ? static_cast<double>(sc.coreHit) / sc.core : 0.0;
        out.scores.push_back(sc);
    }
    std::stable_sort(out.scores.begin(), out.scores.end(),
                     [](const DetectModelScore& a, const DetectModelScore& b) { return a.score > b.score; });
}

// TESSERACT_OM2_DETECT_SKETCH, else om2detect.sketch next to the models (main() resolves the
// model directory exactly as it does for --organism).
std::string sketchPath(const std::string& modelDir) {
    if (const char* p = env::text("TESSERACT_OM2_DETECT_SKETCH")) return p;
    return modelDir.empty() ? std::string() : modelDir + "/om2detect.sketch";
}

std::string detectMode() {
    const char* m = env::text("TESSERACT_OM2_DETECT");
    return m ? std::string(m) : std::string("off");
}

}  // namespace

DetectResult& lastDetection() {
    static DetectResult r;   // process-wide record of the one detection a run makes; not a cached flag
    return r;
}

std::string md5OfFile(const std::string& path) {
    std::FILE* f = std::fopen(path.c_str(), "rb");
    if (!f) return std::string();
    Md5 h;
    std::vector<char> buf(1 << 20);
    size_t n;
    while ((n = std::fread(buf.data(), 1, buf.size(), f)) > 0) h.add(buf.data(), n);
    const bool bad = std::ferror(f) != 0;
    std::fclose(f);
    return bad ? std::string() : h.finish();
}

bool sketchFromModel(const std::string& tsmPath, DetectSketchModel& out, std::string& error) {
    std::FILE* raw = std::fopen(tsmPath.c_str(), "rb");
    if (!raw) { error = "cannot open model file: " + tsmPath; return false; }
    const std::unique_ptr<std::FILE, int (*)(std::FILE*)> guard(raw, std::fclose);
    std::FILE* f = raw;
    char magic[8];
    if (std::fread(magic, 1, 8, f) != 8 || std::memcmp(magic, "TSMODEL", 7) != 0 || magic[7] < '2' || magic[7] > '5') {
        error = "not a TesserACT model file: " + tsmPath;
        return false;
    }
    const bool hasDensity = magic[7] == '5';
    uint32_t kk = 0, denom = 512, gChr = 0, gPls = 0;
    bool ok = rdPod(f, kk) && kk == static_cast<uint32_t>(kMarkerK);
    if (ok && hasDensity) ok = rdPod(f, denom) && denom > 0;
    ok = ok && rdPod(f, gChr) && rdPod(f, gPls);
    std::string name;
    ok = ok && rdStr(f, name);
    uint32_t nExcluded = 0;
    ok = ok && rdPod(f, nExcluded) && nExcluded < (1u << 24);
    for (uint32_t i = 0; ok && i < nExcluded; ++i) {
        std::string acc;
        ok = rdStr(f, acc);
    }
    uint64_t nMarkers = 0;
    ok = ok && rdPod(f, nMarkers) && nMarkers < (1ULL << 34);
    const double coreT = 0.9 * gChr;
    out = DetectSketchModel();
    for (uint64_t i = 0; ok && i < nMarkers; ++i) {
        uint64_t km = 0;
        uint32_t gc = 0, gp = 0;
        ok = rdPod(f, km) && rdPod(f, gc) && rdPod(f, gp);
        // Core: single-copy on >= 90% of the panel chromosomes and on no plasmid set
        // (detect_org2.cpp cls bit 4, the same expression).
        if (ok && gc >= coreT && gp == 0 && gChr > 0) out.core.push_back(km);
    }
    if (!ok) { error = "truncated or unreadable model file: " + tsmPath; return false; }
    std::sort(out.core.begin(), out.core.end());
    out.core.erase(std::unique(out.core.begin(), out.core.end()), out.core.end());
    out.name = name;
    out.genomesChr = gChr;
    out.denom = denom;
    out.tsmMd5 = md5OfFile(tsmPath);
    out.holdList = "-";
    out.holdMd5 = "-";
    return true;
}

bool writeDetectSketch(const std::string& path, const DetectSketch& s, std::string& error) {
    std::FILE* f = std::fopen(path.c_str(), "wb");
    if (!f) { error = "cannot write " + path; return false; }
    bool ok = std::fwrite(kSketchMagic, 1, 8, f) == 8;
    const uint32_t n = static_cast<uint32_t>(s.models.size());
    ok = ok && wrPod(f, n);
    for (const DetectSketchModel& m : s.models) {
        ok = ok && wrStr(f, m.name) && wrStr(f, m.tsmMd5) && wrStr(f, m.holdList) && wrStr(f, m.holdMd5) &&
             wrPod(f, m.genomesChr) && wrPod(f, m.denom);
        const uint64_t nc = m.core.size();
        ok = ok && wrPod(f, nc);
        if (ok && nc) ok = std::fwrite(m.core.data(), sizeof(uint64_t), nc, f) == nc;
    }
    if (std::fclose(f) != 0) ok = false;
    if (!ok) error = "write failed on " + path;
    return ok;
}

bool loadDetectSketch(const std::string& path, DetectSketch& s, std::string& error) {
    s.models.clear();
    std::FILE* raw = std::fopen(path.c_str(), "rb");
    if (!raw) { error = "cannot open detection sketch: " + path; return false; }
    const std::unique_ptr<std::FILE, int (*)(std::FILE*)> guard(raw, std::fclose);
    char magic[8];
    if (std::fread(magic, 1, 8, raw) != 8 || std::memcmp(magic, kSketchMagic, 8) != 0) {
        error = "not an om2 detection sketch: " + path;
        return false;
    }
    uint32_t n = 0;
    bool ok = rdPod(raw, n) && n <= 4096;
    for (uint32_t i = 0; ok && i < n; ++i) {
        DetectSketchModel m;
        uint64_t nc = 0;
        ok = rdStr(raw, m.name) && rdStr(raw, m.tsmMd5) && rdStr(raw, m.holdList) && rdStr(raw, m.holdMd5) &&
             rdPod(raw, m.genomesChr) && rdPod(raw, m.denom) && m.denom > 0 && rdPod(raw, nc) && nc < (1ULL << 26);
        if (ok && nc) {
            m.core.resize(nc);
            ok = std::fread(m.core.data(), sizeof(uint64_t), nc, raw) == nc;
        }
        if (ok) s.models.push_back(std::move(m));
    }
    if (!ok) { error = "truncated or corrupt detection sketch: " + path; s.models.clear(); return false; }
    return true;
}

bool scoreReadFiles(const DetectSketch& s, const std::vector<std::string>& files, DetectResult& out,
                    std::string& error) {
    const auto t0 = std::chrono::steady_clock::now();
    std::unordered_map<uint64_t, uint32_t> slot;
    std::vector<uint64_t> denoms;
    buildSlots(s, slot, denoms);
    std::vector<FileScan> scans(files.size());
    std::vector<std::thread> pool;
    for (size_t i = 0; i < files.size(); ++i)
        pool.emplace_back([&, i]() { scanFile(files[i], denoms, slot, scans[i]); });
    for (std::thread& t : pool) t.join();
    std::vector<uint32_t> counts(slot.size(), 0);
    out.reads = out.bases = 0;
    for (const FileScan& fs : scans) {
        if (!fs.error.empty()) { error = fs.error; return false; }
        for (size_t k = 0; k < counts.size(); ++k) counts[k] += fs.counts[k];
        out.reads += fs.reads;
        out.bases += fs.bases;
    }
    finishScores(s, slot, counts, out);
    out.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return true;
}

void scoreSequences(const DetectSketch& s, const std::vector<std::string>& seqs, DetectResult& out) {
    std::unordered_map<uint64_t, uint32_t> slot;
    std::vector<uint64_t> denoms;
    buildSlots(s, slot, denoms);
    std::vector<uint32_t> counts(slot.size(), 0);
    out.reads = out.bases = 0;
    for (const std::string& q : seqs) {
        ++out.reads;
        out.bases += q.size();
        countSequence(q, denoms, slot, counts);
    }
    finishScores(s, slot, counts, out);
}

bool detectRequested(const std::string& organism) {
    return organism == "auto" || detectMode() != "off";
}

int detectOrganism(std::string& organism, bool force, const std::vector<std::string>& files,
                   const std::string& modelDir, std::FILE* log) {
    DetectResult& r = lastDetection();
    r = DetectResult();
    const std::string mode = detectMode();
    const bool isAuto = organism == "auto";
    r.mode = isAuto ? "auto" : mode;
    r.requested = organism;
    r.sketchPath = sketchPath(modelDir);
    // Nothing to check: detection only speaks about a model, and none was asked for.
    if (!isAuto && organism.empty()) {
        std::fprintf(log, "[om2-detect] mode=%s requested=none decision=skipped (no --organism)\n", r.mode.c_str());
        return 0;
    }
    DetectSketch sk;
    std::string err;
    if (r.sketchPath.empty() || !loadDetectSketch(r.sketchPath, sk, err)) {
        if (err.empty()) err = "no detection sketch (set TESSERACT_OM2_DETECT_SKETCH or TESSERACT_MODEL_DIR)";
        // auto and gate cannot do what they were asked without a sketch; warn only warns.
        if (isAuto || mode == "gate") {
            std::fprintf(log, "error: organism detection: %s\n", err.c_str());
            return 2;
        }
        std::fprintf(log, "      warning: organism detection skipped: %s\n", err.c_str());
        r.decision = "skipped";
        return 0;
    }
    r.sketchMd5 = md5OfFile(r.sketchPath);
    std::vector<std::string> present;
    for (const std::string& f : files)
        if (!f.empty() && util::fileExists(f)) present.push_back(f);
    if (!scoreReadFiles(sk, present, r, err)) {
        std::fprintf(log, "error: organism detection: %s\n", err.c_str());
        return 2;
    }
    r.ran = true;
    r.mixture = r.second() >= kDetectMixture;
    if (isAuto) {
        if (!r.scores.empty() && r.best() >= kDetectAccept) {
            r.chosen = canonicalOrganism(r.scores[0].name);
            r.decision = "auto";
        } else {
            r.chosen.clear();
            r.decision = "auto_none";
        }
        organism = r.chosen;
    } else {
        const std::string want = lowerCase(canonicalOrganism(organism));
        for (const DetectModelScore& sc : r.scores)
            if (lowerCase(canonicalOrganism(sc.name)) == want) r.requestedScore = sc.score;
        r.chosen = organism;
        if (r.requestedScore >= kDetectAccept) r.decision = "accepted";
        else if (mode == "gate" && !force) r.decision = "refused";
        else r.decision = force && mode == "gate" ? "forced" : "warned";
    }
    std::fprintf(log, "[om2-detect] mode=%s requested=%s best=%s:%.4f second=%s:%.4f requested_score=%.4f "
                      "decision=%s mixture=%d reads=%llu seconds=%.1f sketch_md5=%s\n",
                 r.mode.c_str(), organism.empty() && isAuto ? "auto" : r.requested.c_str(),
                 r.scores.empty() ? "none" : r.scores[0].name.c_str(), r.best(),
                 r.scores.size() < 2 ? "none" : r.scores[1].name.c_str(), r.second(), r.requestedScore,
                 r.decision.c_str(), r.mixture ? 1 : 0, static_cast<unsigned long long>(r.reads), r.seconds,
                 r.sketchMd5.c_str());
    if (r.mixture)
        std::fprintf(log, "      warning: reads also match %s at %.3f (>= %.2f): a mixed or contaminated sample?\n",
                     r.scores[1].name.c_str(), r.second(), kDetectMixture);
    if (r.decision == "auto_none")
        std::fprintf(log, "      warning: --organism auto: no model reaches %.2f (best %s %.3f); assembling "
                          "without a model\n", kDetectAccept, r.scores.empty() ? "none" : r.scores[0].name.c_str(),
                     r.best());
    if (r.decision == "warned" || r.decision == "forced")
        std::fprintf(log, "      warning: --organism %s scores %.3f against these reads (accept >= %.2f; best %s "
                          "%.3f)%s\n", organism.c_str(), r.requestedScore, kDetectAccept,
                     r.scores.empty() ? "none" : r.scores[0].name.c_str(), r.best(),
                     r.decision == "forced" ? "; --organism-force given, continuing" : "");
    if (r.decision == "refused") {
        std::fprintf(log, "error: --organism %s scores %.3f against these reads, below %.2f (best %s %.3f).\n"
                          "  The model would place contigs by another organism's gene order. Use the right\n"
                          "  --organism, drop it, or pass --organism-force to run anyway.\n",
                     organism.c_str(), r.requestedScore, kDetectAccept,
                     r.scores.empty() ? "none" : r.scores[0].name.c_str(), r.best());
        return 2;
    }
    return 0;
}

int checkSketchModel(const std::string& modelPath, std::FILE* log) {
    DetectResult& r = lastDetection();
    if (!r.ran || modelPath.empty()) return 0;
    DetectSketch sk;
    std::string err;
    if (!loadDetectSketch(r.sketchPath, sk, err)) return 0;   // already reported by detectOrganism
    const std::string want = lowerCase(canonicalOrganism(r.chosen));
    const DetectSketchModel* m = nullptr;
    for (const DetectSketchModel& x : sk.models)
        if (lowerCase(canonicalOrganism(x.name)) == want) m = &x;
    if (!m) { r.modelCheck = "absent"; return 0; }
    const std::string md5 = md5OfFile(modelPath);
    r.modelCheck = md5 == m->tsmMd5 ? "match" : "mismatch";
    if (r.modelCheck == "mismatch") {
        const bool gate = r.mode == "gate";
        std::fprintf(log, "%s: the detection sketch was built from %s md5 %s, but this run loads %s md5 %s\n",
                     gate ? "error" : "      warning", m->name.c_str(), m->tsmMd5.c_str(), modelPath.c_str(),
                     md5.c_str());
        if (gate) return 2;
    }
    return 0;
}

}  // namespace om2
}  // namespace ts
