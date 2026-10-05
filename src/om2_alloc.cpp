#include "om2_alloc.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <map>
#include <sstream>

#include "kmer.h"

namespace ts {
namespace om2 {

const char* c2BasisName(Basis b) {
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

// ---- allocation ------------------------------------------------------------------------------

int majorityBranch(const SiteInput& s) {
    int best = -1;
    for (size_t b = 0; b < s.branches.size(); ++b) {
        if (best < 0) { best = static_cast<int>(b); continue; }
        const SiteBranch& x = s.branches[b];
        const SiteBranch& y = s.branches[static_cast<size_t>(best)];
        if (x.carriers > y.carriers || (x.carriers == y.carriers && x.cov > y.cov)) best = static_cast<int>(b);
    }
    return best;
}

double phasePosterior(const std::vector<uint32_t>& n, size_t b, double eps) {
    const size_t m = n.size();
    if (m == 0 || b >= m) return 0.0;
    if (m == 1) return 1.0;
    // P(branch j | counts) is proportional to r^n_j with r = (1-eps)(m-1)/eps.
    const double lr = std::log((1.0 - eps) * static_cast<double>(m - 1) / eps);
    double denom = 0;
    for (size_t j = 0; j < m; ++j) {
        denom += std::exp(lr * (static_cast<double>(n[j]) - static_cast<double>(n[b])));
    }
    return denom > 0 ? 1.0 / denom : 0.0;
}

namespace {

struct Decision { int branch = -1; Basis basis = Basis::Consensus; double marginal = 0; bool exhausted = false; };

// The choice made with no prior: the highest budget marginal when it reaches the bar, else
// the majority branch.
Decision decideNoPrior(const SiteInput& s, const std::vector<int>& rem, int maj, double marginalMin) {
    Decision d;
    long total = 0;
    bool changed = false;
    for (size_t b = 0; b < rem.size(); ++b) {
        total += std::max(0, rem[b]);
        if (rem[b] != s.branches[b].carriers) changed = true;
    }
    if (total <= 0) {
        d.branch = maj;
        d.basis = Basis::Consensus;
        d.exhausted = true;
        return d;
    }
    int best = -1;
    double bestM = -1;
    for (size_t b = 0; b < rem.size(); ++b) {
        const double m = static_cast<double>(std::max(0, rem[b])) / static_cast<double>(total);
        if (m > bestM + 1e-12 || (std::fabs(m - bestM) <= 1e-12 && static_cast<int>(b) == maj)) {
            bestM = m;
            best = static_cast<int>(b);
        }
    }
    if (bestM >= marginalMin) {
        d.branch = best;
        d.marginal = bestM;
        d.basis = (best != maj || changed) ? Basis::Multiplicity : Basis::Consensus;
    } else {
        d.branch = maj;
        d.marginal = static_cast<double>(std::max(0, rem[static_cast<size_t>(maj)])) / static_cast<double>(total);
        d.basis = Basis::Consensus;
    }
    return d;
}

void tally(AllocCounts& c, Basis b) {
    switch (b) {
        case Basis::PairPhased: ++c.pair; break;
        case Basis::ThreadPhased: ++c.thread; break;
        case Basis::PriorAllocated: ++c.prior; break;
        case Basis::Multiplicity: ++c.multiplicity; break;
        default: ++c.consensus; break;
    }
}

}  // namespace

std::vector<SiteChoice> allocateSite(const SiteInput& s, const AllocParams& p, AllocCounts& counts) {
    std::vector<SiteChoice> out(s.loci.size());
    const size_t m = s.branches.size();
    if (m == 0 || p.mode == AllocMode::Off) return out;
    const int maj = majorityBranch(s);
    if (m == 1 || s.complex || p.mode == AllocMode::Consensus) {
        for (SiteChoice& c : out) {
            c.branch = maj;
            c.basis = Basis::Consensus;
            c.noPriorBranch = maj;
            tally(counts, c.basis);
        }
        return out;
    }
    std::vector<int> rem(m);
    for (size_t b = 0; b < m; ++b) rem[b] = s.branches[b].carriers;
    std::vector<char> done(s.loci.size(), 0);

    // A: loci whose own flank-anchored fragments name a branch.
    for (size_t i = 0; i < s.loci.size(); ++i) {
        const SiteLocus& L = s.loci[i];
        if (L.threadN.size() != m || L.pairN.size() != m) continue;
        std::vector<uint32_t> n(m);
        uint32_t tot = 0;
        for (size_t b = 0; b < m; ++b) { n[b] = L.threadN[b] + L.pairN[b]; tot += n[b]; }
        if (tot == 0) continue;
        size_t best = 0;
        for (size_t b = 1; b < m; ++b) if (n[b] > n[best]) best = b;
        const double post = phasePosterior(n, best, p.eps);
        if (post < p.phaseMin) continue;
        SiteChoice& c = out[i];
        c.branch = static_cast<int>(best);
        c.basis = L.threadN[best] > 0 && L.threadN[best] >= L.pairN[best] ? Basis::ThreadPhased : Basis::PairPhased;
        c.posterior = post;
        c.noPriorBranch = c.branch;
        --rem[best];
        done[i] = 1;
        tally(counts, c.basis);
    }

    // B: the leave-clone-out locus prior, ITS-class sites only, strongest first. The prior names an
    // ITS class; the branch within the class is the one with the most budget left (then depth).
    if (p.mode == AllocMode::Prior && s.itsClass) {
        auto classOf = [&](int b) -> int {
            if (b < 0) return -1;
            return s.branchClass.empty() ? b : s.branchClass[static_cast<size_t>(b)];
        };
        std::vector<std::pair<double, size_t>> order;
        for (size_t i = 0; i < s.loci.size(); ++i) {
            const SiteLocus& L = s.loci[i];
            if (done[i] || L.prior.empty() || L.priorGenomes < p.priorGenomesMin) continue;
            const double top = *std::max_element(L.prior.begin(), L.prior.end());
            order.emplace_back(-top, i);
        }
        std::sort(order.begin(), order.end());
        for (const auto& oi : order) {
            const size_t i = oi.second;
            const SiteLocus& L = s.loci[i];
            size_t cp = 0;
            for (size_t c = 1; c < L.prior.size(); ++c) if (L.prior[c] > L.prior[cp]) cp = c;
            if (L.prior[cp] < p.priorMin) continue;
            int bp = -1;
            for (size_t b = 0; b < m; ++b) {
                if (classOf(static_cast<int>(b)) != static_cast<int>(cp) || rem[b] <= 0) continue;
                if (bp < 0 || rem[b] > rem[static_cast<size_t>(bp)] ||
                    (rem[b] == rem[static_cast<size_t>(bp)] && s.branches[b].cov > s.branches[static_cast<size_t>(bp)].cov)) {
                    bp = static_cast<int>(b);
                }
            }
            if (bp < 0) continue;   // no branch of that class left in the budget
            const Decision cf = decideNoPrior(s, rem, maj, p.marginalMin);
            SiteChoice& c = out[i];
            c.noPriorBranch = cf.branch;
            if (classOf(cf.branch) != static_cast<int>(cp)) {
                c.branch = bp;
                c.basis = Basis::PriorAllocated;
                c.posterior = L.prior[cp];
                c.priorFlip = true;
            } else {
                c.branch = cf.branch;
                c.basis = cf.basis;
                c.posterior = cf.marginal;
                ++counts.priorNoFlip;
                if (cf.exhausted) ++counts.budgetExhausted;
            }
            --rem[static_cast<size_t>(c.branch)];
            done[i] = 1;
            tally(counts, c.basis);
        }
    }

    // C: everything else, in the order given.
    for (size_t i = 0; i < s.loci.size(); ++i) {
        if (done[i]) continue;
        const Decision d = decideNoPrior(s, rem, maj, p.marginalMin);
        SiteChoice& c = out[i];
        c.branch = d.branch;
        c.basis = d.basis;
        c.posterior = d.marginal;
        c.noPriorBranch = d.branch;
        if (d.exhausted) ++counts.budgetExhausted;
        --rem[static_cast<size_t>(d.branch)];
        done[i] = 1;
        tally(counts, c.basis);
    }
    return out;
}

// ---- hashing ---------------------------------------------------------------------------------

uint64_t sidecarHash(uint64_t canonical31) { return mix64(canonical31); }

namespace {

struct Md5 {
    uint32_t h[4] = {0x67452301u, 0xefcdab89u, 0x98badcfeu, 0x10325476u};
    unsigned char buf[64];
    size_t used = 0;
    uint64_t bytes = 0;
    uint32_t K[64];
    Md5() {
        for (int i = 0; i < 64; ++i) {
            K[i] = static_cast<uint32_t>(static_cast<uint64_t>(std::fabs(std::sin(static_cast<double>(i + 1))) *
                                                               4294967296.0));
        }
    }
    static uint32_t rotl(uint32_t x, int c) { return (x << c) | (x >> (32 - c)); }
    void block(const unsigned char* q) {
        static const int S[64] = {7, 12, 17, 22, 7, 12, 17, 22, 7, 12, 17, 22, 7, 12, 17, 22,
                                  5, 9,  14, 20, 5, 9,  14, 20, 5, 9,  14, 20, 5, 9,  14, 20,
                                  4, 11, 16, 23, 4, 11, 16, 23, 4, 11, 16, 23, 4, 11, 16, 23,
                                  6, 10, 15, 21, 6, 10, 15, 21, 6, 10, 15, 21, 6, 10, 15, 21};
        uint32_t M[16];
        for (int i = 0; i < 16; ++i) {
            M[i] = static_cast<uint32_t>(q[4 * i]) | (static_cast<uint32_t>(q[4 * i + 1]) << 8) |
                   (static_cast<uint32_t>(q[4 * i + 2]) << 16) | (static_cast<uint32_t>(q[4 * i + 3]) << 24);
        }
        uint32_t A = h[0], B = h[1], C = h[2], D = h[3];
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
        h[0] += A; h[1] += B; h[2] += C; h[3] += D;
    }
    void addRaw(const unsigned char* p, size_t n) {
        while (n) {
            const size_t take = std::min(n, 64 - used);
            std::memcpy(buf + used, p, take);
            used += take; p += take; n -= take;
            if (used == 64) { block(buf); used = 0; }
        }
    }
    void add(const void* p, size_t n) { bytes += n; addRaw(static_cast<const unsigned char*>(p), n); }
    std::string finish() {
        const uint64_t bits = bytes * 8;
        const unsigned char one = 0x80, zero = 0;
        addRaw(&one, 1);
        while (used != 56) addRaw(&zero, 1);
        unsigned char len[8];
        for (int i = 0; i < 8; ++i) len[i] = static_cast<unsigned char>(bits >> (8 * i));
        addRaw(len, 8);
        char out[33];
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) std::snprintf(out + 8 * i + 2 * j, 3, "%02x", (h[i] >> (8 * j)) & 0xffu);
        }
        return std::string(out, 32);
    }
};

// Canonical 31-mers of a sequence: fn(canonical, position, forwardIsCanonical).
template <typename F>
void forEach31(const std::string& s, F&& fn) {
    const int k = 31;
    const uint64_t mask = (1ULL << (2 * k)) - 1;
    uint64_t fwd = 0, rev = 0;
    int valid = 0;
    for (size_t i = 0; i < s.size(); ++i) {
        const int c = baseCode(s[i]);
        if (c < 0) { valid = 0; fwd = rev = 0; continue; }
        fwd = ((fwd << 2) | static_cast<uint64_t>(c)) & mask;
        rev = (rev >> 2) | (static_cast<uint64_t>(3 - c) << (2 * (k - 1)));
        if (++valid < k) continue;
        fn(fwd <= rev ? fwd : rev, i + 1 - static_cast<size_t>(k), fwd <= rev);
    }
}

bool parseHex(const std::string& t, uint64_t& v) {
    if (t.empty() || t.size() > 16) return false;
    char* end = nullptr;
    v = std::strtoull(t.c_str(), &end, 16);
    return end && *end == '\0';
}

// Reads the '#key value' header of a sidecar and checks its format tag and model pin.
bool readHeader(std::istream& in, const char* format, const std::string& modelMd5,
                std::map<std::string, std::string>& hdr, std::string& firstBody, std::string& error) {
    std::string line;
    bool sawFormat = false;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] != '#') { firstBody = line; break; }
        std::istringstream ls(line.substr(1));
        std::string key, value;
        ls >> key;
        std::getline(ls, value);
        const size_t f = value.find_first_not_of(" \t");
        value = f == std::string::npos ? std::string() : value.substr(f);
        if (key == format) sawFormat = true;
        hdr[key] = value;
    }
    if (!sawFormat) { error = std::string("not a ") + format + " file"; return false; }
    const auto it = hdr.find("tsm_md5");
    if (it == hdr.end()) { error = "no tsm_md5 pin in the header"; return false; }
    if (modelMd5.empty()) { error = "no organism model loaded to check the tsm_md5 pin against"; return false; }
    if (it->second != modelMd5) {
        error = "built for model md5 " + it->second + ", the loaded model is " + modelMd5;
        return false;
    }
    return true;
}

std::vector<std::string> splitOn(const std::string& s, char sep) {
    std::vector<std::string> out;
    std::string cur;
    for (char c : s) {
        if (c == sep) { out.push_back(cur); cur.clear(); }
        else cur.push_back(c);
    }
    out.push_back(cur);
    return out;
}

}  // namespace

std::string md5File(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) return std::string();
    Md5 m;
    std::vector<char> buf(1 << 16);
    while (in) {
        in.read(buf.data(), static_cast<std::streamsize>(buf.size()));
        m.add(buf.data(), static_cast<size_t>(in.gcount()));
    }
    if (!in.eof()) return std::string();
    return m.finish();
}

// ---- rRNA locus prior ------------------------------------------------------------------------

bool RrnPrior::load(const std::string& path, const std::string& modelMd5, std::string& error) {
    loaded_ = false;
    std::ifstream in(path);
    if (!in) { error = "cannot open " + path; return false; }
    std::map<std::string, std::string> hdr;
    std::string line;
    if (!readHeader(in, "om2rrn", modelMd5, hdr, line, error)) return false;
    if (hdr.count("flank_denom")) denom_ = std::strtoull(hdr["flank_denom"].c_str(), nullptr, 10);
    if (hdr.count("flank_bp")) flankBp_ = std::atoi(hdr["flank_bp"].c_str());
    if (denom_ == 0) { error = "flank_denom must be positive"; return false; }
    do {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ls(line);
        std::string tag;
        ls >> tag;
        if (tag == "C") {
            int id = -1;
            size_t n = 0;
            std::string hashes;
            ls >> id >> n >> hashes;
            if (id < 0 || id > 4096) { error = "bad class line: " + line; return false; }
            if (classes_.size() <= static_cast<size_t>(id)) classes_.resize(static_cast<size_t>(id) + 1);
            for (const std::string& t : splitOn(hashes, ',')) {
                uint64_t v;
                if (!t.empty() && parseHex(t, v)) classes_[static_cast<size_t>(id)].insert(v);
            }
        } else if (tag == "L") {
            int id = -1;
            uint32_t genomes = 0;
            std::string counts, markers;
            ls >> id >> genomes >> counts >> markers;
            if (id < 0 || id > 100000) { error = "bad locus line: " + line; return false; }
            if (loci_.size() <= static_cast<size_t>(id)) loci_.resize(static_cast<size_t>(id) + 1);
            Locus& L = loci_[static_cast<size_t>(id)];
            L.genomes = genomes;
            for (const std::string& t : splitOn(counts, ';')) {
                const size_t c = t.find(':');
                if (c == std::string::npos) continue;
                L.counts[std::atoi(t.substr(0, c).c_str())] =
                    static_cast<uint32_t>(std::strtoul(t.substr(c + 1).c_str(), nullptr, 10));
            }
            for (const std::string& t : splitOn(markers, ',')) {
                uint64_t v;
                if (!t.empty() && parseHex(t, v)) markerLoci_[v].push_back(static_cast<uint32_t>(id));
            }
        }
    } while (std::getline(in, line));
    loaded_ = !loci_.empty() && !classes_.empty();
    if (!loaded_) error = "no loci or no classes in " + path;
    return loaded_;
}

int RrnPrior::matchLocus(const std::vector<uint64_t>& markers, int* score) const {
    std::unordered_map<uint32_t, int> votes;
    std::unordered_set<uint64_t> seen;
    for (uint64_t h : markers) {
        if (!seen.insert(h).second) continue;
        const auto it = markerLoci_.find(h);
        if (it == markerLoci_.end()) continue;
        for (uint32_t l : it->second) ++votes[l];
    }
    int best = -1, bestV = 0, second = 0;
    for (const auto& kv : votes) {
        if (kv.second > bestV || (kv.second == bestV && static_cast<int>(kv.first) < best)) {
            second = std::max(second, bestV);
            bestV = kv.second;
            best = static_cast<int>(kv.first);
        } else {
            second = std::max(second, kv.second);
        }
    }
    if (score) *score = bestV;
    if (best < 0 || bestV < 3 || bestV < 2 * second) return -1;
    return best;
}

int RrnPrior::classOf(const std::string& seq) const {
    std::vector<int> hits(classes_.size(), 0);
    std::unordered_set<uint64_t> seen;
    forEach31(seq, [&](uint64_t canon, size_t, bool) {
        const uint64_t h = sidecarHash(canon);
        if (!seen.insert(h).second) return;
        for (size_t c = 0; c < classes_.size(); ++c) if (classes_[c].count(h)) ++hits[c];
    });
    int best = -1, bestV = 0, second = 0;
    for (size_t c = 0; c < hits.size(); ++c) {
        if (hits[c] > bestV) { second = bestV; bestV = hits[c]; best = static_cast<int>(c); }
        else second = std::max(second, hits[c]);
    }
    if (best < 0 || bestV < 3 || bestV < 2 * second) return -1;
    return best;
}

uint32_t RrnPrior::locusGenomes(int locus) const {
    if (locus < 0 || static_cast<size_t>(locus) >= loci_.size()) return 0;
    return loci_[static_cast<size_t>(locus)].genomes;
}

uint32_t RrnPrior::classCount(int locus, int cls) const {
    if (locus < 0 || static_cast<size_t>(locus) >= loci_.size()) return 0;
    const auto& c = loci_[static_cast<size_t>(locus)].counts;
    const auto it = c.find(cls);
    return it == c.end() ? 0 : it->second;
}

// ---- dnaA sketch -----------------------------------------------------------------------------

bool DnaaSketch::load(const std::string& path, const std::string& modelMd5, std::string& error) {
    loaded_ = false;
    std::ifstream in(path);
    if (!in) { error = "cannot open " + path; return false; }
    std::map<std::string, std::string> hdr;
    std::string line;
    if (!readHeader(in, "om2dnaa", modelMd5, hdr, line, error)) return false;
    do {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ls(line);
        std::string hex;
        int offset = -1, strand = -1;
        ls >> hex >> offset >> strand;
        uint64_t v;
        if (!parseHex(hex, v) || offset < 0 || (strand != 0 && strand != 1)) {
            error = "bad sketch line: " + line;
            return false;
        }
        Entry e;
        e.offset = offset;
        e.strand = static_cast<uint8_t>(strand);
        entries_[v] = e;
    } while (std::getline(in, line));
    loaded_ = !entries_.empty();
    if (!loaded_) error = "empty sketch " + path;
    return loaded_;
}

DnaaSketch::Hit DnaaSketch::locate(const std::string& rec) const {
    Hit hit;
    const int64_t L = static_cast<int64_t>(rec.size());
    if (!loaded_ || L < 31) return hit;
    std::map<std::pair<int, int64_t>, uint32_t> votes;
    forEach31(rec, [&](uint64_t canon, size_t pos, bool fwdCanon) {
        const auto it = entries_.find(sidecarHash(canon));
        if (it == entries_.end()) return;
        ++hit.hits;
        const int recStrand = fwdCanon ? 0 : 1;
        const int64_t p = static_cast<int64_t>(pos);
        if (recStrand == it->second.strand) {
            votes[{0, ((p - it->second.offset) % L + L) % L}]++;
        } else {
            votes[{1, ((p + 30 + it->second.offset) % L + L) % L}]++;
        }
    });
    // Sum each candidate with its neighbours within 5 bp (indels in the first 1.5 kb).
    std::pair<int, int64_t> best{0, -1};
    uint32_t bestV = 0;
    for (const auto& kv : votes) {
        uint32_t v = 0;
        for (const auto& kw : votes) {
            if (kw.first.first != kv.first.first) continue;
            int64_t d = std::llabs(kw.first.second - kv.first.second);
            d = std::min(d, L - d);
            if (d <= 5) v += kw.second;
        }
        if (v > bestV) { bestV = v; best = kv.first; }
    }
    hit.votes = bestV;
    if (bestV >= 10 && static_cast<double>(bestV) >= 0.6 * static_cast<double>(hit.hits)) {
        hit.found = true;
        hit.reverse = best.first == 1;
        hit.offset = best.second;
    }
    return hit;
}

}  // namespace om2
}  // namespace ts
