#include "libqc.h"

#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace ts {
namespace {

// A minimal JSON scanner, sufficient for a machine-written report and no more.
//
// The temptation is to grep for `"estimated_depth"` and read the number after it. That
// breaks the moment the same key appears under two parents -- and it already does here:
// `peak` lives under both `insert_size` and other blocks, and `k` is a key in several. A
// scanner that tracks nesting can be asked for a PATH, which is unambiguous. Pulling in a
// JSON library for six numbers would be the worse trade.
class Json {
public:
    explicit Json(const std::string& text) : s_(text) {}

    // Walk to a dotted path (e.g. "scepter.kmer_spectrum.estimated_depth") and return the
    // raw token that follows it. Empty string means "not present", which callers must
    // distinguish from "present and zero".
    std::string value(const std::string& path) const {
        std::vector<std::string> want;
        std::string part;
        std::istringstream parts(path);
        while (std::getline(parts, part, '.')) want.push_back(part);
        if (want.empty()) return "";

        std::vector<std::string> stack;   // current object path
        size_t i = 0;
        std::string pendingKey;
        while (i < s_.size()) {
            const char c = s_[i];
            if (c == '"') {
                const std::string tok = readString(i);      // advances i past the string
                skipSpace(i);
                if (i < s_.size() && s_[i] == ':') {
                    ++i;
                    skipSpace(i);
                    pendingKey = tok;
                    if (i < s_.size() && s_[i] == '{') {
                        stack.push_back(pendingKey);
                        ++i;
                        pendingKey.clear();
                        continue;
                    }
                    if (i < s_.size() && s_[i] == '[') { skipArray(i); pendingKey.clear(); continue; }
                    // A scalar: check whether stack + key is the path we were asked for.
                    if (matches(stack, pendingKey, want)) return readScalar(i);
                    skipScalar(i);
                    pendingKey.clear();
                    continue;
                }
                continue;   // a bare string value; the scalar branch above consumed keys
            }
            if (c == '{') { stack.push_back(""); ++i; continue; }
            if (c == '}') { if (!stack.empty()) stack.pop_back(); ++i; continue; }
            if (c == '[') { skipArray(i); continue; }
            ++i;
        }
        return "";
    }

private:
    const std::string& s_;

    void skipSpace(size_t& i) const {
        while (i < s_.size() && std::isspace(static_cast<unsigned char>(s_[i]))) ++i;
    }
    std::string readString(size_t& i) const {
        ++i;                                     // opening quote
        std::string out;
        while (i < s_.size() && s_[i] != '"') {
            if (s_[i] == '\\' && i + 1 < s_.size()) { out += s_[i + 1]; i += 2; continue; }
            out += s_[i++];
        }
        if (i < s_.size()) ++i;                  // closing quote
        return out;
    }
    void skipArray(size_t& i) const {
        int depth = 0;
        for (; i < s_.size(); ++i) {
            if (s_[i] == '"') { size_t j = i; readString(j); i = j - 1; continue; }
            if (s_[i] == '[') ++depth;
            else if (s_[i] == ']') { --depth; if (depth == 0) { ++i; return; } }
        }
    }
public:
    // Numeric array at a dotted path. Returns empty when the path is absent or holds a
    // scalar; an empty histogram and a missing one are treated identically by callers,
    // which is correct here -- neither can support a percentile.
    std::vector<double> array(const std::string& path) const {
        std::vector<std::string> want;
        std::string part;
        std::istringstream parts(path);
        while (std::getline(parts, part, '.')) want.push_back(part);
        if (want.empty()) return {};

        std::vector<std::string> stack;
        size_t i = 0;
        while (i < s_.size()) {
            const char c = s_[i];
            if (c == '"') {
                const std::string tok = readString(i);
                skipSpace(i);
                if (i < s_.size() && s_[i] == ':') {
                    ++i;
                    skipSpace(i);
                    if (i < s_.size() && s_[i] == '{') { stack.push_back(tok); ++i; continue; }
                    if (i < s_.size() && s_[i] == '[') {
                        if (matches(stack, tok, want)) {
                            const size_t start = i;
                            skipArray(i);
                            return parseNumbers(s_.substr(start + 1, i - start - 2));
                        }
                        skipArray(i);
                        continue;
                    }
                    skipScalar(i);
                    continue;
                }
                continue;
            }
            if (c == '{') { stack.push_back(""); ++i; continue; }
            if (c == '}') { if (!stack.empty()) stack.pop_back(); ++i; continue; }
            if (c == '[') { skipArray(i); continue; }
            ++i;
        }
        return {};
    }

private:
    static std::vector<double> parseNumbers(const std::string& body) {
        std::vector<double> out;
        std::string tok;
        std::istringstream in(body);
        while (std::getline(in, tok, ',')) {
            const size_t a = tok.find_first_not_of(" \t\r\n");
            if (a == std::string::npos) continue;
            out.push_back(std::atof(tok.c_str() + a));
        }
        return out;
    }

    std::string readScalar(size_t& i) const {
        size_t j = i;
        skipScalar(j);
        return s_.substr(i, j - i);
    }
    void skipScalar(size_t& i) const {
        if (i < s_.size() && s_[i] == '"') { readString(i); return; }
        while (i < s_.size() && s_[i] != ',' && s_[i] != '}' && s_[i] != ']') ++i;
    }
    static bool matches(const std::vector<std::string>& stack, const std::string& key,
                        const std::vector<std::string>& want) {
        if (want.empty() || key != want.back()) return false;
        // Match the parent chain from the inside out, ignoring the anonymous frames that
        // an unnamed `{` pushes.
        size_t w = want.size() - 1;
        for (size_t si = stack.size(); si-- > 0 && w > 0;) {
            if (stack[si].empty()) continue;
            if (stack[si] == want[w - 1]) --w;
            else return false;
        }
        return w == 0;
    }
};

double num(const Json& j, const char* path, bool* found = nullptr) {
    const std::string v = j.value(path);
    if (found) *found = !v.empty();
    if (v.empty()) return 0.0;
    return std::atof(v.c_str());
}

bool flag(const Json& j, const char* path) {
    const std::string v = j.value(path);
    return v == "true" || v == "1";
}

}  // namespace

double LibraryQC::expectedKmerCoverage(int k) const {
    if (!loaded || readDepth <= 0 || meanReadLength <= 0) return 0.0;
    const double span = meanReadLength - static_cast<double>(k) + 1.0;
    if (span <= 0) return 0.0;
    return readDepth * span / meanReadLength;
}

bool loadLibraryQC(const std::string& path, LibraryQC& out, std::string& error) {
    std::ifstream in(path, std::ios::binary);
    if (!in) { error = "cannot open QC report: " + path; return false; }
    std::ostringstream buf;
    buf << in.rdbuf();
    const std::string text = buf.str();
    if (text.empty()) { error = "QC report is empty: " + path; return false; }

    const Json j(text);
    LibraryQC qc;
    qc.source = path;

    // The QC tool was renamed fastplus -> scepter, and the JSON block carrying everything
    // below was renamed with it. Reports written under either name are accepted: a rename
    // is no reason to reject a measurement of a library that may no longer be available to
    // re-measure. The block name is resolved ONCE and every field path is built from it,
    // so the two spellings cannot drift apart field by field.
    std::string root = "scepter";
    bool haveDepth = false;
    qc.readDepth = num(j, (root + ".kmer_spectrum.estimated_depth").c_str(), &haveDepth);
    if (!haveDepth) {
        root = "fastplus";
        qc.readDepth = num(j, (root + ".kmer_spectrum.estimated_depth").c_str(), &haveDepth);
    }
    if (!haveDepth) {
        error = "QC report has no scepter.kmer_spectrum block (was it run with "
                "--preset wgs-bacteria?): " + path;
        return false;
    }
    const std::string ks = root + ".kmer_spectrum.";
    const std::string em = root + ".error_model.";
    qc.spectrumReliable = flag(j, (ks + "reliable").c_str());
    qc.spectrumK        = static_cast<int>(num(j, (ks + "k").c_str()));
    qc.peakDepth        = num(j, (ks + "peak_depth").c_str());
    qc.genomeSize       = num(j, (ks + "estimated_genome_size").c_str());
    qc.genomeSizeCore   = num(j, (ks + "estimated_genome_size_core").c_str());
    qc.errorKmerFraction = num(j, (ks + "error_kmer_fraction").c_str());
    qc.repeatFraction   = num(j, (ks + "repeat_fraction").c_str());

    qc.errorModelHasData = flag(j, (em + "has_data").c_str());
    qc.errorRate         = num(j, (em + "error_rate").c_str());
    qc.insertPeak        = num(j, "insert_size.peak");

    // Fit the insert model from the histogram rather than from the reported peak alone.
    // The peak is a single bin; what the resolver actually needs is a plausible RANGE, and
    // the range is what the assembler's own estimator gets wrong -- it fits a standard
    // deviation to a heavy-tailed sample and then opens the window to +/-4 of it, which on
    // a typical library here spans [0, 3.3x the mean] and therefore excludes nothing.
    // Percentiles of the observed distribution cannot blow up that way.
    const std::vector<double> hist = j.array("insert_size.histogram");
    double total = 0;
    for (double v : hist) total += v;
    if (total >= 1000) {
        double running = 0, sum = 0;
        int p1 = -1, p99 = -1;
        for (size_t bin = 0; bin < hist.size(); ++bin) {
            running += hist[bin];
            sum += hist[bin] * static_cast<double>(bin);
            if (p1 < 0 && running >= total * 0.01) p1 = static_cast<int>(bin);
            if (p99 < 0 && running >= total * 0.99) p99 = static_cast<int>(bin);
        }
        const double mean = sum / total;
        double var = 0;
        for (size_t bin = 0; bin < hist.size(); ++bin) {
            const double d = static_cast<double>(bin) - mean;
            var += hist[bin] * d * d;
        }
        qc.insertMean = mean;
        qc.insertSd = std::sqrt(var / total);
        qc.insertP1 = p1 < 0 ? 0 : p1;
        qc.insertP99 = p99 < 0 ? static_cast<int>(hist.size()) : p99;
        qc.insertObservations = static_cast<size_t>(total);
        // The last bin of a fixed-width histogram is an overflow bucket: everything longer
        // than the range lands in it. If that bucket holds the 99th percentile the
        // distribution is not contained by the histogram and no percentile from it can be
        // trusted, so the model is refused rather than reported with a bogus upper bound.
        qc.insertUsable = mean > 0 && qc.insertP99 > qc.insertP1 &&
                          static_cast<size_t>(qc.insertP99) + 1 < hist.size();
    }

    const double reads = num(j, (root + ".run.reads_in").c_str());
    const double bases = num(j, (root + ".run.bases_in").c_str());
    if (reads > 0 && bases > 0) qc.meanReadLength = bases / reads;

    if (qc.readDepth <= 0) {
        error = "QC report carries a non-positive depth estimate: " + path;
        return false;
    }
    qc.loaded = true;
    out = qc;
    return true;
}

}  // namespace ts
