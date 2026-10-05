// Organism Model 2.0, component C3 (SURFACE): output writer. See om2_output.h.
#include "om2_output.h"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstring>
#include <fstream>
#include <map>
#include <sstream>
#include <tuple>
#include <unordered_map>

#include "envflags.h"
#include "om2_clonal.h"
#include "organism_detect.h"
#include "seqio.h"
#include "util.h"

namespace ts {
namespace om2 {

namespace {

// ---- SHA-256, for the model and sidecar identities in report.json -------------------------
class Sha256 {
    uint32_t h[8] = {0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19};
    uint64_t bytes = 0;
    size_t used = 0;
    unsigned char block[64]{};
    static uint32_t rotr(uint32_t x, int n) { return (x >> n) | (x << (32 - n)); }
    void compress() {
        static const uint32_t k[64] = {
            0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
            0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
            0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
            0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
            0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
            0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
            0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
            0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2};
        uint32_t w[64];
        for (int i = 0; i < 16; ++i)
            w[i] = (uint32_t(block[4 * i]) << 24) | (uint32_t(block[4 * i + 1]) << 16) |
                   (uint32_t(block[4 * i + 2]) << 8) | block[4 * i + 3];
        for (int i = 16; i < 64; ++i) {
            const uint32_t a = w[i - 15], b = w[i - 2];
            w[i] = w[i - 16] + (rotr(a, 7) ^ rotr(a, 18) ^ (a >> 3)) + w[i - 7] + (rotr(b, 17) ^ rotr(b, 19) ^ (b >> 10));
        }
        uint32_t a = h[0], b = h[1], c = h[2], d = h[3], e = h[4], f = h[5], g = h[6], z = h[7];
        for (int i = 0; i < 64; ++i) {
            const uint32_t t1 = z + (rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25)) + ((e & f) ^ (~e & g)) + k[i] + w[i];
            const uint32_t t2 = (rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22)) + ((a & b) ^ (a & c) ^ (b & c));
            z = g; g = f; f = e; e = d + t1; d = c; c = b; b = a; a = t1 + t2;
        }
        h[0] += a; h[1] += b; h[2] += c; h[3] += d; h[4] += e; h[5] += f; h[6] += g; h[7] += z;
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
        unsigned char end[8];
        for (int i = 0; i < 8; ++i) end[7 - i] = static_cast<unsigned char>(bits >> (i * 8));
        add(end, 8);
        char out[65];
        for (int i = 0; i < 8; ++i) std::snprintf(out + i * 8, 9, "%08x", h[i]);
        return std::string(out, 64);
    }
};

// ---- JSON (the report.json om2 block, rendered at the report's indentation) ----------------
void jesc(const std::string& in, std::string& out) {
    for (char c : in) {
        switch (c) {
            case '"': out += "\\\""; break;
            case '\\': out += "\\\\"; break;
            case '\n': out += "\\n"; break;
            case '\t': out += "\\t"; break;
            default:
                if (static_cast<unsigned char>(c) < 0x20) {
                    char b[8];
                    std::snprintf(b, sizeof b, "\\u%04x", c);
                    out += b;
                } else {
                    out += c;
                }
        }
    }
}
struct Json {
    std::string s;
    int indent = 1;   // the om2 object sits one level inside the report
    void pad() { s.append(static_cast<size_t>(indent) * 4, ' '); }
    void key(const char* k) { pad(); s += '"'; s += k; s += "\": "; }
    void str(const char* k, const std::string& v) { key(k); s += '"'; jesc(v, s); s += "\",\n"; }
    void num(const char* k, double v) {
        key(k);
        if (!std::isfinite(v)) { s += "null,\n"; return; }
        char b[64];
        std::snprintf(b, sizeof b, "%.6g", v);
        s += b; s += ",\n";
    }
    void numOrNull(const char* k, double v, bool isNull) { if (isNull) { key(k); s += "null,\n"; } else num(k, v); }
    void uint(const char* k, unsigned long long v) { key(k); s += std::to_string(v); s += ",\n"; }
    void boolean(const char* k, bool v) { key(k); s += v ? "true" : "false"; s += ",\n"; }
    void open(const char* k, char c) {
        if (k) key(k); else pad();
        s += c; s += '\n'; ++indent;
    }
    void close(char c) {
        if (s.size() >= 2 && s[s.size() - 2] == ',') s.erase(s.size() - 2, 1);
        --indent;
        pad();
        s += c; s += ",\n";
    }
};

std::string lastN(const std::string& s, size_t n) { return s.size() <= n ? s : s.substr(s.size() - n); }
std::string firstN(const std::string& s, size_t n) { return s.substr(0, std::min(n, s.size())); }

std::string fmtFloat(double v, const char* f = "%.4g") {
    char b[64];
    std::snprintf(b, sizeof b, f, v);
    return b;
}

// A seam is an N-run of a record written by this run.
struct Seam {
    size_t rec = 0, pos = 0, len = 0;
    std::string fl, fr;
    int ji = -1;           // ledger junction, -1 unrecorded
    bool flip = false;     // matched in the reverse orientation
    bool ambiguous = false;
    Admit admit = Admit::Scaffold;
};

struct Seg {
    enum Kind { Comp, Gap, Fill, Merge } kind = Comp;   // Merge: a clonal exact-overlap join (no bases, C4 round 2)
    std::string comp;      // component id (Comp, Fill)
    size_t cb = 0, ce = 0; // 1-based inclusive range on the component (Comp, Fill)
    char ori = '+';
    size_t len = 0;        // bases in the genome record
    int seam = -1;         // Gap/Fill: seam index; Fill of a wrap: -1
    int ji = -1;           // ledger junction behind a Gap or Fill
    std::string text;      // Comp/Fill: the bases as written into genome.fasta
};

struct GRec {
    TagInfo tag;
    size_t rec = 0;                   // source record
    bool whole = true;                // covers its whole source record
    std::vector<Seg> segs;
    std::string name;
    bool circular = false, rotated = false, wrapOpen = false;
    std::string circBasis = "none";
    size_t gaps = 0, fills = 0, allocated = 0, length = 0;
    double expMis = 0;
    bool expNA = false;
};

bool fillValid(const Junction& j) {
    if (j.fillSeq.empty() || j.contigFilled) return false;
    uint32_t at = 0;
    for (const FillSpan& f : j.fillSpans) {
        if (f.off != at || f.len == 0 || f.basis == Basis::Observed) return false;
        at += f.len;
    }
    return at == j.fillSeq.size();
}

// The fill in the orientation it is written, lowercase inside allocation spans.
std::string fillText(const Junction& j) {
    std::string s = j.fillSeq;
    for (char& c : s) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    for (const FillSpan& f : j.fillSpans)
        if (isAllocation(f.basis))
            for (uint32_t i = f.off; i < f.off + f.len && i < s.size(); ++i)
                s[i] = static_cast<char>(std::tolower(static_cast<unsigned char>(s[i])));
    return s;
}

std::string partName(const std::string& name, size_t k, size_t len) {
    std::string out;
    if (name.compare(0, 5, "NODE_") == 0) {
        size_t p = 5;
        while (p < name.size() && std::isdigit(static_cast<unsigned char>(name[p]))) ++p;
        out = name.substr(0, p) + "." + std::to_string(k) + name.substr(p);
        const size_t lp = out.find("_length_");
        if (lp != std::string::npos) {
            size_t q = lp + 8;
            while (q < out.size() && std::isdigit(static_cast<unsigned char>(out[q]))) ++q;
            out = out.substr(0, lp + 8) + std::to_string(len) + out.substr(q);
        }
        return out;
    }
    return name + "_om2part" + std::to_string(k);
}

bool writeText(const std::string& path, const std::string& text, std::string& error) {
    std::FILE* f = std::fopen(path.c_str(), "wb");
    if (!f) { error = "cannot write " + path + ": " + std::strerror(errno); return false; }
    const bool ok = std::fwrite(text.data(), 1, text.size(), f) == text.size();
    const bool closed = std::fclose(f) == 0;
    if (!ok || !closed) { error = "write failed on " + path; return false; }
    return true;
}

std::vector<std::string> splitTab(const std::string& l) {
    std::vector<std::string> v;
    size_t p = 0;
    while (true) {
        const size_t t = l.find('\t', p);
        v.push_back(l.substr(p, t == std::string::npos ? std::string::npos : t - p));
        if (t == std::string::npos) break;
        p = t + 1;
    }
    return v;
}

const char* gapType(const Junction* j) {
    if (j && (j->endA == EndClass::Repeat || j->endB == EndClass::Repeat)) return "repeat";
    return "scaffold";
}
char gapComponent(const Junction* j, size_t len) {
    return (j && j->source == Source::ResolverUnknown100 && len == 100) ? 'U' : 'N';
}

// What would close an open gap, by what the ledger knows about it.
std::string closureClass(const Junction* j) {
    if (!j) return "unrecorded";
    if (j->endA == EndClass::Repeat || j->endB == EndClass::Repeat) return "repeat_bounded";
    if (j->endA == EndClass::DeadEnd || j->endB == EndClass::DeadEnd) return "dead_end";
    if (j->verdict == Verdict::Silent || j->tier == Tier::E) return "panel_only";
    if (j->verdict == Verdict::Unanchored || j->verdict == Verdict::Abstain) return "unanchored";
    return "unique_open";
}
std::string closureAdvice(const std::string& cls, int64_t longest) {
    if (cls == "unrecorded")
        return "not classified: this build has no C1 junction ledger, so the gap's provenance "
               "(pairs, model join, layout) and its end classes are unknown";
    if (cls == "repeat_bounded") {
        const int64_t kb = (std::max<int64_t>(longest, 0) + 1000 + 999) / 1000;
        return "repeat-bounded (IS/rRNA class): a read spanning the repeat walk plus its flanks, >= " +
               std::to_string(kb) + " kb (longest walk class " + std::to_string(longest) + " bp)";
    }
    if (cls == "dead_end") return "no graph path from a flank: more coverage, or long reads across the gap";
    if (cls == "panel_only") return "ordered by the panel only (align_genus): a long read or this isolate's own closed genome";
    if (cls == "unanchored") return "an end could not be anchored in the graph: long reads";
    return "unique flanks without a closing walk (coverage dropout): more coverage";
}

}  // namespace

// ---- names ----------------------------------------------------------------------------------
std::string revcomp(const std::string& s) {
    std::string r(s.rbegin(), s.rend());
    for (char& c : r) {
        switch (c) {
            case 'A': c = 'T'; break; case 'C': c = 'G'; break; case 'G': c = 'C'; break; case 'T': c = 'A'; break;
            case 'a': c = 't'; break; case 'c': c = 'g'; break; case 'g': c = 'c'; break; case 't': c = 'a'; break;
            default: break;
        }
    }
    return r;
}
const char* basisName(Basis b) {
    switch (b) {
        case Basis::Observed: return "OBSERVED";
        case Basis::GraphWalk: return "GRAPH_WALK";
        case Basis::PairPhased: return "PAIR_PHASED";
        case Basis::ThreadPhased: return "THREAD_PHASED";
        case Basis::IsolateRepeatExact: return "ISOLATE_REPEAT_EXACT";
        case Basis::PriorAllocated: return "PRIOR_ALLOCATED";
        case Basis::Multiplicity: return "MULTIPLICITY";
        case Basis::Consensus: return "CONSENSUS";
    }
    return "?";
}
bool isAllocation(Basis b) {
    return b == Basis::PriorAllocated || b == Basis::Multiplicity || b == Basis::Consensus;
}
// sourceName, verdictName, tierName, admitName and endClassName are defined once, in
// om2_evidence.cpp (C1), since the integration (build_om2): C1 and C3 each carried a copy.
// C1's spelling is kept: the calibration class keys (om2_seams.tsv `cls`) are built from it.

TagInfo parseTag(const std::string& name) {
    TagInfo t;
    // The tag follows the coverage field: NODE_<n>_length_<L>_cov_<c><tag>.
    size_t p = name.rfind("_cov_");
    std::string tag;
    if (p != std::string::npos) {
        p += 5;
        while (p < name.size() && (std::isdigit(static_cast<unsigned char>(name[p])) || name[p] == '.')) ++p;
        tag = name.substr(p);
    }
    if (tag.compare(0, 4, "_chr") == 0) t.cls = 'c';
    else if (tag.compare(0, 5, "_plas") == 0) {
        t.cls = 'p';
        size_t q = 5;
        if (q < tag.size() && tag[q] == '_' && q + 1 < tag.size() && std::isdigit(static_cast<unsigned char>(tag[q + 1]))) {
            t.group = static_cast<uint32_t>(std::strtoul(tag.c_str() + q + 1, nullptr, 10));
        }
    } else {
        t.cls = 'u';
    }
    t.circular = tag.size() >= 9 && tag.compare(tag.size() - 9, 9, "_circular") == 0;
    return t;
}

std::string agpEvidence(const Junction* j, bool modelRan) {
    if (!j) return modelRan ? "unspecified" : "paired-ends";
    std::vector<std::string> ev;
    const bool pairs = j->source == Source::Resolver || j->source == Source::ResolverUnknown100 ||
                       j->tier == Tier::B || j->verdict == Verdict::ContraPairsOk;
    if (pairs) ev.push_back("paired-ends");
    if (j->source != Source::Resolver && j->source != Source::ResolverUnknown100) ev.push_back("align_genus");
    // Graph evidence has no AGP 2.1 term of its own; DESIGN C3 maps it to `map`.
    if (j->tier == Tier::A || j->tier == Tier::C) ev.push_back("map");
    if (ev.empty()) return "unspecified";
    std::string out;
    for (size_t i = 0; i < ev.size(); ++i) out += (i ? ";" : "") + ev[i];
    return out;
}

std::string sha256File(const std::string& path) {
    std::FILE* f = std::fopen(path.c_str(), "rb");
    if (!f) return std::string();
    Sha256 h;
    std::vector<char> buf(1 << 20);
    size_t n;
    while ((n = std::fread(buf.data(), 1, buf.size(), f)) > 0) h.add(buf.data(), n);
    const bool bad = std::ferror(f) != 0;
    std::fclose(f);
    return bad ? std::string() : h.finish();
}
// md5File is defined once, in om2_alloc.cpp (C2), since the integration.

bool outputEnabled() { return env::on("TESSERACT_OM2_OUTPUT", false); }
bool agpEvidenceEnabled() { return env::on("TESSERACT_OM2_AGP_EVIDENCE", false); }

std::string validateAgp(const std::string& path) {
    std::ifstream in(path);
    if (!in) return "cannot read " + path;
    std::string l;
    std::map<std::string, std::pair<long long, long long>> last;   // object -> (last end, last part)
    size_t lineNo = 0;
    bool version = false;
    while (std::getline(in, l)) {
        ++lineNo;
        if (l.empty()) continue;
        if (l[0] == '#') { if (l.compare(0, 13, "##agp-version") == 0) version = l.find("2.1") != std::string::npos; continue; }
        const std::vector<std::string> x = splitTab(l);
        const std::string at = path + ":" + std::to_string(lineNo) + ": ";
        if (x.size() != 9) return at + "expected 9 columns, got " + std::to_string(x.size());
        const long long ob = std::atoll(x[1].c_str()), oe = std::atoll(x[2].c_str()), part = std::atoll(x[3].c_str());
        auto& st = last[x[0]];
        if (ob != st.first + 1) return at + "object coordinates are not contiguous";
        if (oe < ob) return at + "object end before its start";
        if (part != st.second + 1) return at + "part numbers are not consecutive";
        st = {oe, part};
        if (x[4] == "N" || x[4] == "U") {
            const long long gl = std::atoll(x[5].c_str());
            if (gl != oe - ob + 1) return at + "gap length differs from its object span";
            if (x[4] == "U" && gl != 100) return at + "a U gap must be 100 bp";
            static const char* kTypes[] = {"scaffold", "contig", "centromere", "short_arm", "heterochromatin",
                                           "telomere", "repeat", "contamination"};
            if (std::find_if(std::begin(kTypes), std::end(kTypes), [&](const char* t) { return x[6] == t; }) ==
                std::end(kTypes))
                return at + "unknown gap type " + x[6];
            if (x[7] != "yes" && x[7] != "no") return at + "linkage must be yes or no";
            if (x[8].empty()) return at + "missing linkage evidence";
            if (x[7] == "yes" && x[8] == "na") return at + "a linked gap needs linkage evidence other than na";
        } else if (x[4] == "W") {
            const long long cb = std::atoll(x[6].c_str()), ce = std::atoll(x[7].c_str());
            if (cb < 1 || ce < cb) return at + "bad component range";
            if (ce - cb != oe - ob) return at + "component span differs from its object span";
            if (x[8] != "+" && x[8] != "-" && x[8] != "?" && x[8] != "0" && x[8] != "na")
                return at + "bad orientation";
        } else {
            return at + "unsupported component type " + x[4];
        }
    }
    if (!version) return path + ": missing ##agp-version 2.1 header";
    return std::string();
}

// ---- the writer -------------------------------------------------------------------------------
SurfaceStats writeSurface(const SurfaceInput& in, std::string& om2Json) {
    SurfaceStats st;
    om2Json.clear();
    st.enabled = outputEnabled();
    st.agpEvidence = agpEvidenceEnabled();
    const DetectResult& det = lastDetection();
    if (det.ran) { st.detectScore = det.best(); st.detectSecond = det.second(); }
    const bool anyLedger = in.ledger && !in.ledger->j.empty();
    if (!st.enabled && !st.agpEvidence && !det.ran && !anyLedger) return st;   // flags off: nothing at all
    if (!in.seqs || !in.names || in.seqs->size() != in.names->size()) {
        st.ok = false;
        st.error = "om2 output: records and names disagree";
        return st;
    }
    const std::vector<std::string>& seqs = *in.seqs;
    const std::vector<std::string>& names = *in.names;
    const std::vector<Junction> noJ;
    const std::vector<Junction>& J = in.ledger ? in.ledger->j : noJ;

    // ---- seams, and which ledger junction each one is -----------------------------------
    std::unordered_map<std::string, std::vector<std::pair<int, bool>>> byFlank;
    std::vector<int> wrapOrdinal(J.size(), -1);
    int wraps = 0;
    for (size_t i = 0; i < J.size(); ++i) {
        const Junction& j = J[i];
        if (j.source == Source::Wrap) { wrapOrdinal[i] = wraps++; continue; }
        if (j.writtenN < 0) continue;
        byFlank[j.flankL32 + "|" + j.flankR32].push_back({static_cast<int>(i), false});
        byFlank[revcomp(j.flankR32) + "|" + revcomp(j.flankL32)].push_back({static_cast<int>(i), true});
    }
    std::vector<Seam> seams;
    std::vector<std::vector<size_t>> seamsOf(seqs.size());
    std::vector<char> jUsed(J.size(), 0);
    // Integration: C1's own location of its rows (anchors and far probes) decides first, so a
    // repeat-bounded junction whose final flanks equal another's is still found; the exact-flank
    // rule below then matches only runs C1 left unowned, and never a row C1 placed elsewhere.
    auto ownerOf = [&](size_t r, size_t k, size_t p0, size_t e0) -> const RunOwner* {
        if (!in.runOwner || r >= in.runOwner->size() || k >= (*in.runOwner)[r].size()) return nullptr;
        const RunOwner& o = (*in.runOwner)[r][k];
        if (o.id < 0 || static_cast<size_t>(o.id) >= J.size() || o.start != static_cast<int64_t>(p0) ||
            o.end != static_cast<int64_t>(e0))
            return nullptr;
        const Junction& j = J[static_cast<size_t>(o.id)];
        if (j.source == Source::Wrap || j.writtenN < 0 || j.contigFilled) return nullptr;
        return &o;
    };
    if (in.runOwner) {
        for (size_t r = 0; r < seqs.size() && r < in.runOwner->size(); ++r)
            for (size_t k = 0; k < (*in.runOwner)[r].size(); ++k) {
                const RunOwner& o = (*in.runOwner)[r][k];
                if (ownerOf(r, k, static_cast<size_t>(o.start), static_cast<size_t>(o.end)))
                    jUsed[static_cast<size_t>(o.id)] = 1;
            }
    }
    for (size_t r = 0; r < seqs.size(); ++r) {
        const std::string& q = seqs[r];
        size_t p = 0;
        size_t runK = 0;
        while (p < q.size()) {
            if (q[p] != 'N' && q[p] != 'n') { ++p; continue; }
            size_t e = p;
            while (e < q.size() && (q[e] == 'N' || q[e] == 'n')) ++e;
            Seam s;
            s.rec = r; s.pos = p; s.len = e - p;
            s.fl = q.substr(p >= 32 ? p - 32 : 0, std::min<size_t>(32, p));
            s.fr = q.substr(e, 32);
            const RunOwner* own = ownerOf(r, runK++, p, e);
            auto it = own ? byFlank.end() : byFlank.find(s.fl + "|" + s.fr);
            if (own) {
                s.ji = static_cast<int>(own->id);
                s.flip = own->orient == '-';
                s.admit = J[static_cast<size_t>(s.ji)].admit;
                ++st.ownedByC1;
            } else if (it != byFlank.end()) {
                // Distinct junctions only: a palindromic flank pair can list one junction twice.
                std::vector<std::pair<int, bool>> c;
                for (const auto& x : it->second)
                    if (std::none_of(c.begin(), c.end(), [&](const std::pair<int, bool>& y) { return y.first == x.first; }))
                        c.push_back(x);
                if (c.size() == 1 && !jUsed[c[0].first]) {
                    s.ji = c[0].first; s.flip = c[0].second; jUsed[s.ji] = 1;
                    s.admit = J[s.ji].admit;
                } else {
                    s.ambiguous = true;
                }
            }
            if (s.ji < 0) { ++st.unrecorded; if (s.ambiguous) ++st.ambiguous; }
            seamsOf[r].push_back(seams.size());
            seams.push_back(std::move(s));
            p = e;
        }
    }
    st.seams = seams.size();

    // ---- scaffolds.fasta: split where the ledger admits a join to the genome view only ----
    // part boundaries per record: [start, end) in record coordinates, and the part's name
    std::vector<std::vector<std::pair<size_t, size_t>>> parts(seqs.size());
    std::vector<std::vector<std::string>> partNames(seqs.size());
    bool anySplit = false;
    for (size_t r = 0; r < seqs.size(); ++r) {
        size_t start = 0;
        for (size_t si : seamsOf[r]) {
            const Seam& s = seams[si];
            if (s.admit == Admit::Scaffold) continue;
            parts[r].push_back({start, s.pos});
            start = s.pos + s.len;
        }
        parts[r].push_back({start, seqs[r].size()});
        if (parts[r].size() == 1) {
            partNames[r].push_back(names[r]);
        } else {
            anySplit = true;
            st.scaffoldSplits += parts[r].size() - 1;
            for (size_t k = 0; k < parts[r].size(); ++k)
                partNames[r].push_back(partName(names[r], k + 1, parts[r][k].second - parts[r][k].first));
        }
    }
    if (anySplit && in.scaffoldsFile) {
        std::vector<std::string> ps, pn;
        for (size_t r = 0; r < seqs.size(); ++r)
            for (size_t k = 0; k < parts[r].size(); ++k) {
                ps.push_back(seqs[r].substr(parts[r][k].first, parts[r][k].second - parts[r][k].first));
                pn.push_back(partNames[r][k]);
            }
        if (!writeFasta(in.outDir + "/scaffolds.fasta", ps, pn, 80, st.error)) { st.ok = false; return st; }
    }

    // ---- scaffolds.agp: object renames at the splits, and per-gap evidence -------------------
    if (in.scaffoldsFile && (anySplit || st.agpEvidence)) {
        const std::string agpPath = in.outDir + "/scaffolds.agp";
        std::ifstream agpIn(agpPath);
        if (agpIn) {
            std::unordered_map<std::string, size_t> recOf;
            for (size_t r = 0; r < names.size(); ++r) recOf[names[r]] = r;
            std::string out, l;
            std::map<std::string, long long> partNo;
            while (std::getline(agpIn, l)) {
                if (l.empty() || l[0] == '#') { out += l + "\n"; continue; }
                std::vector<std::string> x = splitTab(l);
                auto rit = x.size() == 9 ? recOf.find(x[0]) : recOf.end();
                if (rit == recOf.end()) { out += l + "\n"; continue; }
                const size_t r = rit->second;
                const long long ob = std::atoll(x[1].c_str()), oe = std::atoll(x[2].c_str());
                const bool isGap = x[4] == "N" || x[4] == "U";
                const Seam* seam = nullptr;
                if (isGap)
                    for (size_t si : seamsOf[r])
                        if (static_cast<long long>(seams[si].pos) + 1 == ob &&
                            static_cast<long long>(seams[si].len) == oe - ob + 1) seam = &seams[si];
                if (seam && seam->admit != Admit::Scaffold) continue;   // the gap the split removed
                size_t k = 0;
                while (k + 1 < parts[r].size() && static_cast<long long>(parts[r][k + 1].first) < ob) ++k;
                const long long shift = static_cast<long long>(parts[r][k].first);
                x[0] = partNames[r][k];
                x[1] = std::to_string(ob - shift);
                x[2] = std::to_string(oe - shift);
                x[3] = std::to_string(++partNo[x[0]]);
                if (isGap && st.agpEvidence) {
                    const Junction* j = seam && seam->ji >= 0 ? &J[seam->ji] : nullptr;
                    x[6] = gapType(j);
                    x[7] = "yes";
                    x[8] = agpEvidence(j, in.modelRan);
                    ++st.agpEvidenceGaps;
                }
                std::string line;
                for (size_t c = 0; c < x.size(); ++c) line += (c ? "\t" : "") + x[c];
                out += line + "\n";
            }
            agpIn.close();
            if (!writeText(agpPath, out, st.error)) { st.ok = false; return st; }
        }
    }

    std::vector<GRec> recs;
    if (st.enabled) {
        // ---- genome records ------------------------------------------------------------------
        auto addComp = [&](GRec& g, size_t r, size_t a, size_t b) {
            if (b <= a) return;
            size_t k = 0;
            while (k + 1 < parts[r].size() && parts[r][k + 1].first <= a) ++k;
            Seg s;
            s.kind = Seg::Comp;
            s.comp = partNames[r][k];
            s.cb = a - parts[r][k].first + 1;
            s.ce = b - parts[r][k].first;
            s.len = b - a;
            s.text = seqs[r].substr(a, b - a);
            g.segs.push_back(std::move(s));
        };
        auto addFill = [&](GRec& g, int ji, bool flip, int seam) {
            const Junction& j = J[ji];
            Seg s;
            s.kind = Seg::Fill;
            s.comp = "om2fill_" + std::to_string(j.id);
            s.cb = 1; s.ce = j.fillSeq.size();
            s.ori = flip ? '-' : '+';
            s.len = j.fillSeq.size();
            s.seam = seam; s.ji = ji;
            s.text = flip ? revcomp(fillText(j)) : fillText(j);
            g.segs.push_back(std::move(s));
        };
        auto finish = [&](GRec& g) {
            while (!g.segs.empty() && g.segs.back().kind == Seg::Gap) g.segs.pop_back();
            while (!g.segs.empty() && g.segs.front().kind == Seg::Gap) g.segs.erase(g.segs.begin());
            if (!g.segs.empty()) recs.push_back(std::move(g));
        };
        for (size_t r = 0; r < seqs.size(); ++r) {
            GRec g;
            g.tag = parseTag(names[r]);
            g.rec = r;
            size_t pos = 0;
            for (size_t sk = 0; sk < seamsOf[r].size(); ++sk) {
                const size_t si = seamsOf[r][sk];
                const Seam& s = seams[si];
                addComp(g, r, pos, s.pos);
                size_t skip = 0;
                if (s.admit == Admit::Break) {
                    ++st.broken;
                    g.whole = false;
                    GRec next;
                    next.tag = g.tag; next.rec = r; next.whole = false;
                    finish(g);
                    g = std::move(next);
                } else {
                    if (s.admit == Admit::GenomeOnly) ++st.genomeOnlyJoins;
                    // C4 round 2: a clonal exact-overlap join writes the shared bases once. Verified here on the
                    // record itself: the ov bases before the seam equal the ov bases after it, both inside the
                    // pieces next to the seam; otherwise the seam is written as its N-run.
                    size_t ov = (s.ji >= 0 && !fillValid(J[s.ji]) && J[s.ji].mergeOverlap > 0)
                                    ? static_cast<size_t>(J[s.ji].mergeOverlap) : 0;
                    if (ov) {
                        const size_t nextBound = sk + 1 < seamsOf[r].size() ? seams[seamsOf[r][sk + 1]].pos : seqs[r].size();
                        const size_t b0 = s.pos + s.len;
                        bool ok = s.pos >= pos + ov && b0 + ov <= nextBound;
                        for (size_t i = 0; ok && i < ov; ++i) {
                            ok = std::toupper(static_cast<unsigned char>(seqs[r][s.pos - ov + i])) ==
                                 std::toupper(static_cast<unsigned char>(seqs[r][b0 + i]));
                        }
                        if (!ok) ov = 0;
                    }
                    if (s.ji >= 0 && fillValid(J[s.ji])) {
                        addFill(g, s.ji, s.flip, static_cast<int>(si));
                    } else if (ov) {
                        Seg m;
                        m.kind = Seg::Merge;
                        m.len = 0;
                        m.seam = static_cast<int>(si);
                        m.ji = s.ji;
                        g.segs.push_back(std::move(m));
                        skip = ov;
                        if (in.clonal) in.clonal->noteMerge();
                    } else {
                        Seg gap;
                        gap.kind = Seg::Gap;
                        gap.len = s.len;
                        gap.seam = static_cast<int>(si);
                        gap.ji = s.ji;
                        g.segs.push_back(std::move(gap));
                    }
                }
                pos = s.pos + s.len + skip;
            }
            addComp(g, r, pos, seqs[r].size());
            finish(g);
        }

        // ---- wrap junctions: circular only when the wrap is closed with bases ------------
        for (GRec& g : recs) {
            const std::string& q = seqs[g.rec];
            if (g.whole && g.tag.circular) g.circBasis = "pairs";   // ends joined by pairs; wrap not closed
            if (!g.whole || q.size() < 64) continue;
            const std::string l32 = lastN(q, 32), f32 = firstN(q, 32);
            for (size_t i = 0; i < J.size(); ++i) {
                const Junction& j = J[i];
                if (j.source != Source::Wrap || jUsed[i]) continue;
                bool flip;
                if (j.flankL32 == l32 && j.flankR32 == f32) flip = false;
                else if (j.flankL32 == revcomp(f32) && j.flankR32 == revcomp(l32)) flip = true;
                else continue;
                jUsed[i] = 1;
                const bool closed = j.writtenN >= 0 && j.admit != Admit::Break &&
                                    (fillValid(j) || (j.fillSeq.empty() && j.writtenN == 0));
                if (!closed) { g.wrapOpen = true; break; }
                if (fillValid(j)) addFill(g, static_cast<int>(i), flip, -1);
                // Integration: a wrap closed by an overlap of the record's two ends (claimedN = -overlap,
                // no fill) is closed by writing the repeated end once. Refused, and the record left
                // linear, unless the record's last `overlap` bases equal its first.
                if (j.fillSeq.empty() && j.claimedN < 0) {
                    const size_t ov = static_cast<size_t>(-static_cast<int64_t>(j.claimedN));
                    Seg* last = g.segs.empty() ? nullptr : &g.segs.back();
                    size_t L0 = 0;
                    for (const Seg& sg : g.segs) L0 += sg.len;
                    std::string head;
                    for (const Seg& sg : g.segs) { if (head.size() >= ov) break; head += sg.text; }
                    const bool ok = last && last->kind == Seg::Comp && last->ori == '+' && last->len > ov &&
                                    L0 > 2 * ov && head.size() >= ov &&
                                    last->text.compare(last->text.size() - ov, ov, head, 0, ov) == 0;
                    if (!ok) { g.wrapOpen = true; break; }
                    last->ce -= ov;
                    last->len -= ov;
                    last->text.resize(last->text.size() - ov);
                    ++st.overlapTrimmed;
                }
                g.circular = true;
                g.circBasis = j.tier == Tier::B ? "pairs" : "graph";
                // Rotation to dnaA (C2 records the offset per wrap junction, in the junction's
                // orientation of record + wrap fill). Never inside a gap.
                const int w = wrapOrdinal[i];
                int64_t off0 = (w >= 0 && w < 64 && in.ledger) ? in.ledger->rotateOffset[w] : -1;
                size_t L = 0;
                for (const Seg& s : g.segs) L += s.len;
                // Integration: the origin is located on the finished genome record itself. A
                // reverse-strand origin reverse-complements the record first (every segment turned).
                bool located = false;
                if (in.originLocator) {
                    std::string T;
                    T.reserve(L);
                    for (const Seg& sg : g.segs) T += sg.text;
                    for (char& c : T) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
                    const int64_t h = in.originLocator(T);
                    off0 = -1;
                    if (h <= -2 && static_cast<size_t>(-2 - h) < L) {
                        std::reverse(g.segs.begin(), g.segs.end());
                        for (Seg& sg : g.segs) {
                            if (sg.kind != Seg::Gap) sg.ori = sg.ori == '+' ? '-' : '+';
                            sg.text = revcomp(sg.text);
                        }
                        ++st.reversed;
                        off0 = static_cast<int64_t>(L) - 1 - (-2 - h);
                    } else if (h >= 0 && static_cast<size_t>(h) < L) {
                        off0 = h;
                    }
                    located = true;
                }
                if (off0 > 0 && static_cast<size_t>(off0) < L) {
                    const size_t off = flip && !located ? L - static_cast<size_t>(off0) : static_cast<size_t>(off0);
                    size_t at = 0, idx = 0;
                    while (idx < g.segs.size() && at + g.segs[idx].len <= off) at += g.segs[idx++].len;
                    if (idx < g.segs.size() && g.segs[idx].kind == Seg::Gap && off > at) {
                        ++st.rotateRefused;
                    } else if (idx < g.segs.size()) {
                        std::vector<Seg> rot;
                        const size_t d = off - at;
                        if (d > 0) {
                            const Seg& s = g.segs[idx];
                            Seg a = s, b = s;   // a: genome bases [0, d) of s; b: [d, len)
                            if (s.ori == '+') { a.ce = s.cb + d - 1; b.cb = s.cb + d; }
                            else { a.cb = s.ce - d + 1; b.ce = s.ce - d; }
                            a.len = d; b.len = s.len - d;
                            a.text = s.text.substr(0, d); b.text = s.text.substr(d);
                            rot.push_back(b);
                            for (size_t t = idx + 1; t < g.segs.size(); ++t) rot.push_back(g.segs[t]);
                            for (size_t t = 0; t < idx; ++t) rot.push_back(g.segs[t]);
                            rot.push_back(a);
                        } else {
                            for (size_t t = idx; t < g.segs.size(); ++t) rot.push_back(g.segs[t]);
                            for (size_t t = 0; t < idx; ++t) rot.push_back(g.segs[t]);
                        }
                        g.segs.swap(rot);
                        g.rotated = true;
                        ++st.rotated;
                    }
                }
                break;
            }
        }

        // ---- names: chromosome first, then plasmid groups, then the rest (file order) ------
        std::map<uint32_t, size_t> groupSize, groupSeen;
        for (const GRec& g : recs) if (g.tag.cls == 'p' && g.tag.group > 0) ++groupSize[g.tag.group];
        size_t nChr = 0, nUngrouped = 0, nUnplaced = 0;
        for (GRec& g : recs) {
            if (g.tag.cls == 'c') g.name = "chromosome_" + std::to_string(++nChr);
            else if (g.tag.cls == 'p' && g.tag.group > 0) {
                g.name = "plasmid_" + std::to_string(g.tag.group);
                if (groupSize[g.tag.group] > 1) g.name += "_" + std::to_string(++groupSeen[g.tag.group]);
            } else if (g.tag.cls == 'p') g.name = "plasmid_ungrouped_" + std::to_string(++nUngrouped);
            else g.name = "unplaced_" + std::to_string(++nUnplaced);
        }

        // ---- C4: clonal labels of every junction of the view (only with a clonal surface) ----------
        // The contexts are the genome records as written (fills included), N-free next to the junction;
        // a circular record lends its other end.
        std::map<std::pair<size_t, size_t>, ClonalAnnot> clonalAt;   // (record, segment) -> labels
        if (in.clonal) {
            std::vector<std::string> texts;
            texts.reserve(recs.size());
            for (const GRec& g : recs) {
                std::string t;
                for (const Seg& sg : g.segs) {
                    if (sg.kind == Seg::Gap) t.append(sg.len, 'N');
                    else t += sg.text;
                }
                texts.push_back(std::move(t));
            }
            in.clonal->prepare(texts);
            std::vector<ClonalSurface::Req> reqs;
            std::vector<std::pair<size_t, size_t>> where;
            const size_t kCtx = 30000;
            for (size_t gi = 0; gi < recs.size(); ++gi) {
                const GRec& g = recs[gi];
                const std::string& T = texts[gi];
                size_t at = 0;
                for (size_t k = 0; k < g.segs.size(); ++k) {
                    const Seg& sg = g.segs[k];
                    const size_t s0 = at, s1 = at + sg.len;
                    at = s1;
                    if (sg.kind == Seg::Comp) continue;
                    ClonalSurface::Req rq;
                    if (g.circular && T.size() > sg.len) {
                        const size_t L = std::min(kCtx, T.size() - sg.len);
                        const std::string TT = T + T;
                        rq.left = TT.substr(s0 + T.size() - L, L);
                        rq.right = TT.substr(s1, L);
                    } else {
                        rq.left = T.substr(s0 > kCtx ? s0 - kCtx : 0, s0 > kCtx ? kCtx : s0);
                        rq.right = T.substr(s1, std::min(kCtx, T.size() - s1));
                    }
                    rq.j = sg.ji >= 0 ? &J[static_cast<size_t>(sg.ji)] : nullptr;
                    rq.asserted = static_cast<int64_t>(sg.len);
                    rq.sized = sg.kind == Seg::Fill || gapComponent(rq.j, sg.len) != 'U';
                    rq.wrap = rq.j && rq.j->source == Source::Wrap;
                    rq.replicon = g.tag.cls;
                    reqs.push_back(std::move(rq));
                    where.emplace_back(gi, k);
                }
            }
            const std::vector<ClonalAnnot> an = in.clonal->annotate(reqs);
            for (size_t i = 0; i < an.size() && i < where.size(); ++i) clonalAt[where[i]] = an[i];
        }
        auto clonalCols = [&](const ClonalAnnot* a) -> std::string {
            if (!in.clonal) return std::string();
            return "\t" + ClonalSurface::columns(a ? *a : ClonalAnnot());
        };
        auto clonalOf = [&](size_t gi, size_t k) -> const ClonalAnnot* {
            auto it = clonalAt.find({gi, k});
            return it == clonalAt.end() ? nullptr : &it->second;
        };

        // ---- write genome/ --------------------------------------------------------------------
        const std::string gdir = in.outDir + "/genome";
        if (!util::makeDirs(gdir)) { st.ok = false; st.error = "cannot create " + gdir; return st; }
        std::string modelSha8 = "none";
        std::string modelSha, modelMd5;
        if (!in.modelPath.empty()) {
            modelSha = sha256File(in.modelPath);
            modelMd5 = md5File(in.modelPath);
            modelSha8 = in.modelName + "@" + firstN(modelSha, 8);
        }
        std::vector<std::string> gSeq, gHdr, fSeq, fName;
        std::string agp = "##agp-version\t2.1\n"
                          "# TesserACT om2 genome view (genome.fasta). W components are scaffolds.fasta records\n"
                          "# (contigs.fasta records in a gap-free run) or om2fill_<junction> in fills.fasta.\n"
                          "# A gap row is an adjacency asserted without sequence; column 9 says what asserts it.\n";
        std::string bed = "#record\tstart0\tend\tbasis\tjunction\tdetail\n";
        std::string jt =
            "junction\tstatus\tgenome_record\tgenome_start0\tgenome_end\tscaffold_record\tscaffold_pos0\tn_len\t"
            "fill_len\tview\tagp_evidence\tgap_type\tsource\tverdict\ttier\tadmit\tcls\tclaimedN\twrittenN\t"
            "a_piece\ta_tail\tb_piece\tb_tail\tendA\tendB\tcopyA\tcopyB\tgmin\twalks\texhaustive\thairpin\tuA\tuB\t"
            "pairLambda\tpairK\tpairContra\tpanelSupport\tpanelGenomes\ttieMargin\tp_misjoin\tallowCloseGaps\t"
            "contigFilled\tfill_bases\tclosure_class\tflankL32\tflankR32" +
            (in.clonal ? "\t" + ClonalSurface::header() : std::string()) + "\n";
        auto ledgerCols = [&](const Junction* j) {
            if (!j) {
                std::string s;
                for (int c = 0; c < 31; ++c) s += "\t.";
                return s;
            }
            std::string w;
            for (size_t i = 0; i < j->walks.size(); ++i)
                w += (i ? ";" : "") + std::to_string(j->walks[i].len) + ":" + std::to_string(j->walks[i].n);
            if (w.empty()) w = ".";
            std::string fb;
            for (size_t i = 0; i < j->fillSpans.size(); ++i)
                fb += std::string(i ? ";" : "") + basisName(j->fillSpans[i].basis) + ":" +
                      std::to_string(j->fillSpans[i].len);
            if (fb.empty()) fb = ".";
            std::ostringstream o;
            o << '\t' << sourceName(j->source) << '\t' << verdictName(j->verdict) << '\t' << tierName(j->tier)
              << '\t' << admitName(j->admit) << '\t' << (j->cls.empty() ? "." : j->cls) << '\t' << j->claimedN
              << '\t' << j->writtenN << '\t' << j->a.piece << '\t' << (j->a.tail ? 1 : 0) << '\t' << j->b.piece
              << '\t' << (j->b.tail ? 1 : 0) << '\t' << endClassName(j->endA) << '\t' << endClassName(j->endB)
              << '\t' << fmtFloat(j->copyA) << '\t' << fmtFloat(j->copyB) << '\t' << j->gmin << '\t' << w << '\t'
              << (j->exhaustive ? 1 : 0) << '\t' << (j->hairpin ? 1 : 0) << '\t' << j->uA << '\t' << j->uB << '\t'
              << fmtFloat(j->pairLambda) << '\t' << j->pairK << '\t' << j->pairContra << '\t' << j->panelSupport
              << '\t' << j->panelGenomes << '\t' << fmtFloat(j->tieMargin) << '\t'
              << (j->pMisjoin < 0 ? std::string("NA") : fmtFloat(j->pMisjoin)) << '\t' << (j->allowCloseGaps ? 1 : 0)
              << '\t' << (j->contigFilled ? 1 : 0) << '\t' << fb;
            return o.str();
        };
        std::vector<char> seamListed(seams.size(), 0);
        std::string closure;
        size_t unrecordedNo = 0;
        std::map<size_t, std::string> unrecordedId;
        for (size_t i = 0; i < seams.size(); ++i)
            if (seams[i].ji < 0) unrecordedId[i] = "u" + std::to_string(++unrecordedNo);
        std::map<std::string, std::pair<size_t, int64_t>> closeCount;   // per record, reused
        for (GRec& g : recs) {
            std::string seq;
            size_t part = 0;
            closeCount.clear();
            const size_t gIdx = static_cast<size_t>(&g - &recs[0]);
            for (const Seg& s : g.segs) {
                const size_t ob = seq.size() + 1, oe = seq.size() + s.len;
                const Junction* j = s.ji >= 0 ? &J[s.ji] : nullptr;
                const ClonalAnnot* ca = in.clonal ? clonalOf(gIdx, static_cast<size_t>(&s - &g.segs[0])) : nullptr;
                if (ca && s.kind != Seg::Comp) agp += in.clonal->agpComment(j ? std::to_string(j->id) : std::string("unrecorded"), *ca);
                if (s.kind == Seg::Merge) {
                    const Seam& sm = seams[s.seam];
                    seamListed[s.seam] = 1;
                    jt += std::to_string(j->id) + "\tmerged_overlap\t" + g.name + "\t" + std::to_string(ob - 1) + "\t" +
                          std::to_string(ob - 1) + "\t" + names[sm.rec] + "\t" + std::to_string(sm.pos) + "\t" +
                          std::to_string(sm.len) + "\t0\t" +
                          (sm.admit == Admit::GenomeOnly ? "genome_only" : "scaffolds+genome") + "\t.\t." + ledgerCols(j) +
                          "\tmerged:" + std::to_string(j->mergeOverlap) + "\t" + (sm.fl.empty() ? "." : sm.fl) + "\t" +
                          (sm.fr.empty() ? "." : sm.fr) + clonalCols(ca) + "\n";
                    continue;
                }
                if (s.kind == Seg::Gap) {
                    const Seam& sm = seams[s.seam];
                    const char ct = gapComponent(j, s.len);
                    agp += g.name + "\t" + std::to_string(ob) + "\t" + std::to_string(oe) + "\t" +
                           std::to_string(++part) + "\t" + ct + "\t" + std::to_string(s.len) + "\t" + gapType(j) +
                           "\tyes\t" + agpEvidence(j, in.modelRan) + "\n";
                    seq.append(s.len, 'N');
                    ++g.gaps;
                    if (!j || j->pMisjoin < 0) g.expNA = true; else g.expMis += j->pMisjoin;
                    const std::string cc = closureClass(j);
                    auto& c = closeCount[cc];
                    ++c.first;
                    if (j) {
                        int64_t longest = j->gmin;
                        for (const WalkClass& w : j->walks) longest = std::max<int64_t>(longest, w.len);
                        c.second = std::max(c.second, longest);
                    }
                    seamListed[s.seam] = 1;
                    jt += (j ? std::to_string(j->id) : unrecordedId[s.seam]) + "\t" + (j ? "open_gap" : "unrecorded_gap") +
                          "\t" + g.name + "\t" + std::to_string(ob - 1) + "\t" + std::to_string(oe) + "\t" +
                          names[sm.rec] + "\t" + std::to_string(sm.pos) + "\t" + std::to_string(sm.len) + "\t0\t" +
                          (sm.admit == Admit::GenomeOnly ? "genome_only" : "scaffolds+genome") + "\t" +
                          agpEvidence(j, in.modelRan) + "\t" + gapType(j) + ledgerCols(j) + "\t" + cc + "\t" +
                          (sm.fl.empty() ? "." : sm.fl) + "\t" + (sm.fr.empty() ? "." : sm.fr) + clonalCols(ca) + "\n";
                    continue;
                }
                agp += g.name + "\t" + std::to_string(ob) + "\t" + std::to_string(oe) + "\t" + std::to_string(++part) +
                       "\tW\t" + s.comp + "\t" + std::to_string(s.cb) + "\t" + std::to_string(s.ce) + "\t" + s.ori + "\n";
                if (s.kind == Seg::Fill) {
                    ++g.fills;
                    // The fills.fasta record, once per junction, in the junction's own orientation.
                    if (std::find(fName.begin(), fName.end(), s.comp) == fName.end()) {
                        fName.push_back(s.comp);
                        fSeq.push_back(fillText(*j));
                        ++st.fills;
                    }
                    // Basis intervals of the component range this segment carries.
                    for (const FillSpan& f : j->fillSpans) {
                        const size_t a = std::max<size_t>(f.off, s.cb - 1), b = std::min<size_t>(f.off + f.len, s.ce);
                        if (b <= a) continue;
                        size_t g0, g1;
                        if (s.ori == '+') { g0 = ob - 1 + (a - (s.cb - 1)); g1 = ob - 1 + (b - (s.cb - 1)); }
                        else { g0 = ob - 1 + (s.ce - b); g1 = ob - 1 + (s.ce - a); }
                        bed += g.name + "\t" + std::to_string(g0) + "\t" + std::to_string(g1) + "\t" +
                               basisName(f.basis) + "\t" + std::to_string(j->id) + "\t" + s.comp + ":" +
                               std::to_string(a + 1) + "-" + std::to_string(b) + (s.ori == '-' ? ":rc" : "") + "\n";
                        if (isAllocation(f.basis)) g.allocated += b - a;
                    }
                    if (j->pMisjoin < 0) g.expNA = true; else g.expMis += j->pMisjoin;
                    const bool isWrap = j->source == Source::Wrap;
                    const bool firstPiece = s.seam < 0 ? true : !seamListed[s.seam];
                    if (s.seam >= 0) seamListed[s.seam] = 1;
                    if (firstPiece) {
                        const Seam* sm = s.seam >= 0 ? &seams[s.seam] : nullptr;
                        jt += std::to_string(j->id) + "\t" + (isWrap ? "wrap_closed" : "filled_genome") + "\t" + g.name +
                              "\t" + std::to_string(ob - 1) + "\t" + std::to_string(oe) + "\t" + names[g.rec] + "\t" +
                              (sm ? std::to_string(sm->pos) : std::string(".")) + "\t" +
                              (sm ? std::to_string(sm->len) : std::string("0")) + "\t" + std::to_string(j->fillSeq.size()) +
                              "\t" + ((sm && sm->admit == Admit::GenomeOnly) || isWrap ? "genome_only" : "scaffolds+genome") +
                              "\t.\t." + ledgerCols(j) + "\tfilled\t" + (j->flankL32.empty() ? "." : j->flankL32) + "\t" +
                              (j->flankR32.empty() ? "." : j->flankR32) + clonalCols(ca) + "\n";
                    }
                }
                seq += s.text;
            }
            g.length = seq.size();
            st.allocatedBp += g.allocated;
            ++st.records;
            if (g.tag.cls == 'c') ++st.chrRecords;
            if (g.circular) ++st.circular;
            char hdr[512];
            std::snprintf(hdr, sizeof hdr, "%s topology=%s replicon=%s gaps=%zu exp_misjoins=%s allocated_bp=%zu "
                                           "circular_basis=%s model=%s",
                          g.name.c_str(), g.circular ? "circular" : "linear",
                          g.tag.cls == 'c' ? "chr" : g.tag.cls == 'p' ? "plasmid" : "unplaced", g.gaps,
                          g.expNA ? "NA" : fmtFloat(g.expMis, "%.3f").c_str(), g.allocated, g.circBasis.c_str(),
                          modelSha8.c_str());
            gHdr.push_back(hdr);
            gSeq.push_back(std::move(seq));
            closure += g.name + "\tlength=" + std::to_string(g.length) + "\ttopology=" +
                       (g.circular ? "circular" : "linear") + "\tstatus=" +
                       (g.circular && g.gaps == 0 ? "closed" : "open") + "\tgaps=" + std::to_string(g.gaps) +
                       "\tfills=" + std::to_string(g.fills) + "\tallocated_bp=" + std::to_string(g.allocated) +
                       "\tcircular_basis=" + g.circBasis + (g.rotated ? "\trotated=dnaA" : "") + "\n";
            for (const auto& kv : closeCount)
                closure += "    " + std::to_string(kv.second.first) + " gap" + (kv.second.first == 1 ? "" : "s") + " " +
                           kv.first + ": " + closureAdvice(kv.first, kv.second.second) + "\n";
            if (!g.circular) {
                if (g.wrapOpen) closure += "    wrap: a wrap junction is recorded but not closed with bases\n";
                else if (g.circBasis == "pairs")
                    closure += "    wrap: the two ends are joined by read pairs, but no unique graph walk closes "
                               "them with bases; the record is written linear\n";
            }
        }
        // Seams in no genome record (their record was split by a break), and ledger junctions
        // that are not an N-run of the finished records.
        for (size_t i = 0; i < seams.size(); ++i) {
            if (seamListed[i]) continue;
            const Seam& sm = seams[i];
            const Junction* j = sm.ji >= 0 ? &J[sm.ji] : nullptr;
            jt += (j ? std::to_string(j->id) : unrecordedId[i]) + "\tbroken\t.\t.\t.\t" + names[sm.rec] + "\t" +
                  std::to_string(sm.pos) + "\t" + std::to_string(sm.len) + "\t0\tnone\t.\t." + ledgerCols(j) +
                  "\tbroken\t" + (sm.fl.empty() ? "." : sm.fl) + "\t" + (sm.fr.empty() ? "." : sm.fr) +
                  clonalCols(j && in.clonal ? in.clonal->stageRow(j->id) : nullptr) + "\n";
        }
        for (size_t i = 0; i < J.size(); ++i) {
            if (jUsed[i]) continue;
            const Junction& j = J[i];
            const char* status = j.writtenN < 0 ? "broken_by_gate" : j.contigFilled ? "contig_filled"
                                 : j.source == Source::Wrap ? "wrap_unmatched" : "closed_or_absent";
            jt += std::to_string(j.id) + "\t" + status + "\t.\t.\t.\t.\t.\t0\t" + std::to_string(j.fillSeq.size()) +
                  "\tnone\t.\t." + ledgerCols(&j) + "\t.\t" + (j.flankL32.empty() ? "." : j.flankL32) + "\t" +
                  (j.flankR32.empty() ? "." : j.flankR32) + clonalCols(in.clonal ? in.clonal->stageRow(j.id) : nullptr) +
                  "\n";
        }

        // repeat_variants.tsv: inter-copy diversity of multi-copy repeats, as data.
        std::string rv = "family\tcluster\tsite\tjunction\tgenome_record\tgenome_pos0\tsite_len\tstrand\tplaced_allele\tbasis\t"
                         "alleles(allele:carriers_est)\tn_alleles\tminority_carriers_est\tphasing\tcandidate_loci\n";
        if (in.variants) {
            // where each junction's fill landed: junction id -> (record, genome start, orientation, fill length)
            std::map<uint32_t, std::tuple<std::string, size_t, char, size_t>> fillAt;
            for (const GRec& g : recs) {
                size_t at = 0;
                for (const Seg& s : g.segs) {
                    if (s.kind == Seg::Fill && s.cb == 1 && s.ce == J[s.ji].fillSeq.size() && !fillAt.count(J[s.ji].id))
                        fillAt[J[s.ji].id] = std::make_tuple(g.name, at, s.ori, s.len);
                    at += s.len;
                }
            }
            for (const RepeatVariant& v : *in.variants) {
                std::string rec = ".", pos = ".", strand = ".";
                if (v.junction >= 0) {
                    auto it = fillAt.find(static_cast<uint32_t>(v.junction));
                    if (it != fillAt.end() &&
                        static_cast<size_t>(v.fillOffset) + v.siteLen <= std::get<3>(it->second)) {
                        const auto& [rn, at, ori, flen] = it->second;
                        rec = rn;
                        const size_t off = ori == '+' ? v.fillOffset : flen - v.fillOffset - v.siteLen;
                        pos = std::to_string(at + off);
                        strand = std::string(1, ori);
                    }
                }
                std::string al;
                float minority = 0;
                for (size_t i = 0; i < v.alleles.size(); ++i) {
                    al += (i ? ";" : "") + v.alleles[i].first + ":" + fmtFloat(v.alleles[i].second, "%.2f");
                    if (v.alleles[i].first != v.placed) minority += v.alleles[i].second;
                }
                rv += v.family + "\t" + v.cluster + "\t" + v.site + "\t" +
                      (v.junction >= 0 ? std::to_string(v.junction) : std::string(".")) + "\t" + rec + "\t" + pos +
                      "\t" + std::to_string(v.siteLen) + "\t" + strand + "\t" + (v.placed.empty() ? "." : v.placed) + "\t" + (v.placed.empty() ? "." : basisName(v.basis)) +
                      "\t" + (al.empty() ? "." : al) + "\t" + std::to_string(v.alleles.size()) + "\t" +
                      fmtFloat(minority, "%.2f") + "\t" + v.phasing + "\t" + v.loci + "\n";
                ++st.variantRows;
            }
        }

        const bool ledgerKnown = in.ledger != nullptr;
        std::string readme =
            "TesserACT Organism Model 2.0 -- genome view (genome/)\n"
            "=====================================================\n\n"
            "This directory is the OWNER view of the assembly. It is written in addition to, never instead of,\n"
            "the two other views, and each view has its own misassembly budget (om2/design/EVAL_PLAN.md s5.1):\n\n"
            "  contigs.fasta    strict / training view: split at every N-run; only data-proven sequence.\n"
            "                   Genus models are trained on this file.\n"
            "  scaffolds.fasta  verified view: joins whose evidence class is admitted to scaffolds; no\n"
            "                   allocated bases.\n"
            "  genome/genome.fasta  owner view: chromosome first, then plasmid groups, then unplaced records;\n"
            "                   also joins admitted to the genome view only (panel-ordered, evidence tier E,\n"
            "                   `align_genus`), genome-view fills, and labelled allocation.\n\n"
            "Files\n"
            "  genome.fasta          records named chromosome_<n>, plasmid_<group>[_<k>], plasmid_ungrouped_<k>,\n"
            "                        unplaced_<k>. Header fields: topology, replicon, gaps (open N gaps),\n"
            "                        exp_misjoins (sum of calibrated misjoin rates of the record's junctions;\n"
            "                        NA when any junction has no calibrated class), allocated_bp, circular_basis,\n"
            "                        model=<name>@<sha256 prefix>.\n"
            "  genome.agp            AGP 2.1. W components are scaffolds.fasta records (contigs.fasta records\n"
            "                        in a gap-free run) or om2fill_<junction> records of fills.fasta. Gap rows\n"
            "                        carry linkage evidence: paired-ends (read pairs), align_genus (ordered by\n"
            "                        the organism model's panel), map (a graph walk; AGP has no graph term),\n"
            "                        unspecified. gap type `repeat` = an end of the gap is a repeat node.\n"
            "  fills.fasta           bases written into gaps from this isolate's own graph and reads.\n"
            "  genome.mask.bed       the basis label of every filled base (tiles every fill exactly).\n"
            "  junctions.tsv         one row per junction: where it is, what view admitted it, and every\n"
            "                        field of the junction ledger (source, verdict, tier, walks, pairs, ...).\n"
            "  repeat_variants.tsv   inter-copy variant sites of multi-copy repeats (rRNA operons, IS):\n"
            "                        which allele was written at which locus and on what basis, and every\n"
            "                        allele the isolate carries at the site with its estimated carrier copies.\n"
            "                        This table is where the diversity between copies is kept as data,\n"
            "                        including sites whose locus short reads cannot identify.\n"
            "  closure.txt           per record: closed or open, open gaps by class, what would close them.\n\n"
            "Basis labels (genome.mask.bed; lowercase bases in genome.fasta are allocation)\n"
            "  GRAPH_WALK            a unique or forced walk of this isolate's graph, agreed by pairs/threads\n"
            "  PAIR_PHASED           a bubble branch tied to this locus by flank-anchored read pairs (P>=0.99)\n"
            "  THREAD_PHASED         the same, by reads threading the branch and a unique flank\n"
            "  ISOLATE_REPEAT_EXACT  a uniform repeat: every copy identical at read resolution\n"
            "  PRIOR_ALLOCATED       ALLOCATION: ITS class chosen by the leave-clone-out population prior\n"
            "                        (the choice flips without the prior)\n"
            "  MULTIPLICITY          ALLOCATION: the allele forced by the isolate's own copy budget\n"
            "  CONSENSUS             ALLOCATION: the isolate's majority branch; the locus is not identifiable\n"
            "Allocated bases are never written into contigs.fasta or scaffolds.fasta.\n\n"
            "circular_basis\n"
            "  graph   a unique graph walk closes the wrap; the record is written closed (topology=circular)\n"
            "          and, when a dnaA locator was available, starts at dnaA\n"
            "  pairs   read pairs join the two ends, but no walk closes them with bases; topology stays linear\n"
            "  none    no evidence that the record is circular\n\n"
            "What short reads decide, and what they do not (DESIGN.md s11)\n"
            "  decided by this isolate: single-copy adjacency; forward skips; seams a read pair can cross;\n"
            "    unique walks agreed by pairs or threads; variants within flank-anchored pair reach.\n"
            "  decided by the panel, labelled align_genus / tier D-E: order across silent junctions; which exit\n"
            "    of a repeat cluster follows which entry.\n"
            "  allocated, labelled, genome view only: rRNA interior variants beyond pair reach; ITS class by\n"
            "    the population prior; multi-copy branch counts.\n"
            "  never claimed: full operon haplotypes, copy numbers of 3 or more, sequence across dead ends.\n\n";
        readme += ledgerKnown ? "This run: a junction ledger was supplied.\n"
                              : "This run: no junction ledger (C1/C2 are not in this build). Every N-run is kept exactly as\n"
                                "scaffolds.fasta has it and listed as `unrecorded` in junctions.tsv; no fills, no\n"
                                "allocation and no circle closure can occur.\n";

        if (in.clonal) closure = in.clonal->closureHeader() + closure;
        std::string err;
        if (!writeFasta(gdir + "/genome.fasta", gSeq, gHdr, 80, err) ||
            !writeFasta(gdir + "/fills.fasta", fSeq, fName, 80, err) ||
            !writeText(gdir + "/genome.agp", agp, err) || !writeText(gdir + "/genome.mask.bed", bed, err) ||
            !writeText(gdir + "/junctions.tsv", jt, err) || !writeText(gdir + "/repeat_variants.tsv", rv, err) ||
            !writeText(gdir + "/closure.txt", closure, err) || !writeText(gdir + "/README.txt", readme, err)) {
            st.ok = false;
            st.error = err;
            return st;
        }

        // ---- report.json om2 block -----------------------------------------------------------
        Json js;
        js.s = "{\n";
        js.indent = 2;
        js.uint("schema", 1);
        js.boolean("output", true);
        js.boolean("agp_evidence", st.agpEvidence);
        js.boolean("ledger", ledgerKnown);
        js.open("model", '{');
        js.str("name", in.modelName);
        js.str("path", in.modelPath);
        js.str("sha256", modelSha);
        js.str("md5", modelMd5);
        js.close('}');
        js.open("sidecars", '[');
        for (const auto& sc : in.sidecars) {
            js.open(nullptr, '{');
            js.str("label", sc.first);
            js.str("path", sc.second);
            js.str("md5", md5File(sc.second));
            js.close('}');
        }
        js.close(']');
        js.open("records", '[');
        for (const GRec& g : recs) {
            js.open(nullptr, '{');
            js.str("name", g.name);
            js.str("replicon", g.tag.cls == 'c' ? "chr" : g.tag.cls == 'p' ? "plasmid" : "unplaced");
            js.uint("length", g.length);
            js.str("topology", g.circular ? "circular" : "linear");
            js.str("status", g.circular && g.gaps == 0 ? "closed" : "open");
            js.uint("gaps", g.gaps);
            js.uint("fills", g.fills);
            js.uint("allocated_bp", g.allocated);
            js.numOrNull("exp_misjoins", g.expMis, g.expNA);
            js.str("circular_basis", g.circBasis);
            js.boolean("rotated", g.rotated);
            js.close('}');
        }
        js.close(']');
        if (in.clonal) js.s += in.clonal->jsonBlock(js.indent);
        om2Json = js.s;
    }

    // Counters and detection go into the block even when only detection (or only the AGP
    // evidence) ran; the genome part above is present only with TESSERACT_OM2_OUTPUT=1.
    Json js;
    if (om2Json.empty()) {
        js.s = "{\n";
        js.indent = 2;
        js.uint("schema", 1);
        js.boolean("output", false);
        js.boolean("agp_evidence", st.agpEvidence);
        js.boolean("ledger", in.ledger != nullptr);
    } else {
        js.s = om2Json;
        js.indent = 2;
    }
    js.open("counters", '{');
    js.uint("records", st.records);
    js.uint("chr_records", st.chrRecords);
    js.uint("circular", st.circular);
    js.uint("rotated", st.rotated);
    js.uint("genome_only_joins", st.genomeOnlyJoins);
    js.uint("broken", st.broken);
    js.uint("fills", st.fills);
    js.uint("allocated_bp", st.allocatedBp);
    js.uint("agp_evidence_gaps", st.agpEvidenceGaps);
    js.uint("seams", st.seams);
    js.uint("unrecorded", st.unrecorded);
    js.uint("ambiguous", st.ambiguous);
    js.uint("scaffold_splits", st.scaffoldSplits);
    js.uint("variant_rows", st.variantRows);
    js.uint("rotate_refused", st.rotateRefused);
    js.uint("owned_by_c1", st.ownedByC1);
    js.uint("overlap_trimmed", st.overlapTrimmed);
    js.uint("reversed", st.reversed);
    js.close('}');
    if (!in.stageCounters.empty()) {
        js.open("stage_counters", '{');
        for (const auto& sc : in.stageCounters) js.str(sc.first.c_str(), sc.second);
        js.close('}');
    }
    js.open("detection", '{');
    js.boolean("ran", det.ran);
    if (det.ran) {
        js.str("mode", det.mode);
        js.str("decision", det.decision);
        js.str("requested", det.requested);
        js.numOrNull("requested_score", det.requestedScore, det.requestedScore < 0);
        js.str("chosen", det.chosen);
        js.str("best", det.scores.empty() ? "" : det.scores[0].name);
        js.num("best_score", det.best());
        js.str("second", det.scores.size() < 2 ? "" : det.scores[1].name);
        js.numOrNull("second_score", det.second(), det.scores.size() < 2);
        js.boolean("mixture", det.mixture);
        js.str("model_check", det.modelCheck);
        js.str("sketch", det.sketchPath);
        js.str("sketch_md5", det.sketchMd5);
        js.uint("reads", det.reads);
        js.num("seconds", det.seconds);
        js.open("scores", '{');
        for (const DetectModelScore& sc : det.scores) js.num(sc.name.c_str(), sc.score);
        js.close('}');
    }
    js.close('}');
    if (js.s.size() >= 2 && js.s[js.s.size() - 2] == ',') js.s.erase(js.s.size() - 2, 1);
    js.s += "    }";
    om2Json = js.s;
    return st;
}

void logSurfaceCounters(const SurfaceStats& s, std::FILE* log) {
    std::fprintf(log,
                 "[om2-out] enabled=%d records=%zu chr_records=%zu circular=%zu rotated=%zu genome_only_joins=%zu "
                 "fills=%zu allocated_bp=%zu agp_evidence_gaps=%zu detect_score=%s detect_second=%s seams=%zu "
                 "unrecorded=%zu ambiguous=%zu scaffold_splits=%zu broken=%zu variant_rows=%zu agp_evidence=%d "
                 "owned_by_c1=%zu overlap_trimmed=%zu reversed=%zu\n",
                 s.enabled ? 1 : 0, s.records, s.chrRecords, s.circular, s.rotated, s.genomeOnlyJoins, s.fills,
                 s.allocatedBp, s.agpEvidenceGaps, s.detectScore < 0 ? "NA" : fmtFloat(s.detectScore, "%.4f").c_str(),
                 s.detectSecond < 0 ? "NA" : fmtFloat(s.detectSecond, "%.4f").c_str(), s.seams, s.unrecorded,
                 s.ambiguous, s.scaffoldSplits, s.broken, s.variantRows, s.agpEvidence ? 1 : 0, s.ownedByC1,
                 s.overlapTrimmed, s.reversed);
}

}  // namespace om2
}  // namespace ts
