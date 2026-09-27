// T36 (build_v3, G-graph): joinDeadEnds in simplify() -- skipped when it cannot join.
//
// joinDeadEnds forces its overlap to exactly k-1, and build() already links every exact
// (k-1) overlap between solid k-mers, so on a graph from build() it never joins anything,
// yet simplify() called it every round and it copied every dead-end unitig. build_v3 runs
// it only when a cheap exact count of (k-1) tail/head pairs finds a candidate, and prints
//   k=<k> [joindeadends] requested= effective= roundsSkipped= pairsFound= calls= joined=
// on every simplify() call.
//
//  1. Invariant (passes before and after): over random graphs from build(), after random
//     deletions and simplify(), no dead-end pair has an exact (k-1) overlap, and neither
//     simplify() nor a direct joinDeadEnds() joins anything. Prints a digest of every
//     final graph so the release and patched builds can be compared (they must match).
//  2. Fallback (passes before and after): two hand-made unlinked unitigs with a (k-1)
//     overlap are still joined and merged by simplify().
//  3. Counter line (fails on release, which prints none), in forked children because the
//     release reads TESSERACT_JOIN_DEADENDS once per process: default, =40 and =0.
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <fstream>
#include <random>
#include <sstream>
#include <string>
#include <sys/wait.h>
#include <unistd.h>
#include <unordered_map>
#include <vector>

#include "counter.h"
#include "graph.h"
#include "kmer.h"

extern char** environ;
using namespace ts;

namespace {

int failures = 0;
void check(bool ok, const std::string& what) {
    std::printf("  %s  %s\n", ok ? "PASS" : "FAIL", what.c_str());
    std::fflush(stdout);
    if (!ok) ++failures;
}

void clearTesseractEnv() {
    std::vector<std::string> names;
    for (char** e = environ; e && *e; ++e) {
        if (std::strncmp(*e, "TESSERACT_", 10) == 0) {
            const char* eq = std::strchr(*e, '=');
            names.emplace_back(*e, eq ? static_cast<size_t>(eq - *e) : std::strlen(*e));
        }
    }
    for (const std::string& n : names) unsetenv(n.c_str());
}

std::string rc(const std::string& s) { return reverseComplement(s); }

size_t bruteOverlapPairs(const UnitigGraph& g) {
    const size_t K1 = static_cast<size_t>(g.k() - 1);
    std::unordered_map<std::string, std::vector<uint32_t>> heads;
    std::vector<std::pair<uint32_t, std::string>> tails;
    for (uint32_t u = 0; u < g.nodes.size(); ++u) {
        if (g.nodes[u].deleted || g.nodes[u].seq.size() < K1) continue;
        const std::string& s = g.nodes[u].seq;
        for (int e = 0; e < 2; ++e) {
            if (!g.nodes[u].ends[e].empty()) continue;
            const std::string head = e == 0 ? s : rc(s);
            const std::string tail = e == 1 ? s : rc(s);
            heads[head.substr(0, K1)].push_back(u);
            tails.emplace_back(u, tail.substr(tail.size() - K1));
        }
    }
    size_t pairs = 0;
    for (const auto& t : tails) {
        auto it = heads.find(t.second);
        if (it == heads.end()) continue;
        for (uint32_t h : it->second) pairs += h != t.first;
    }
    return pairs;
}

uint64_t fnv(uint64_t h, const std::string& s) {
    for (unsigned char c : s) { h ^= c; h *= 1099511628211ULL; }
    return h;
}
uint64_t graphDigest(uint64_t h, const UnitigGraph& g) {
    for (uint32_t u = 0; u < g.nodes.size(); ++u) {
        const Unitig& n = g.nodes[u];
        h = fnv(h, n.deleted ? "D" : n.seq);
        if (n.deleted) continue;
        for (int e = 0; e < 2; ++e) {
            for (const Link& l : n.ends[e])
                h = fnv(h, std::to_string(e) + ":" + std::to_string(l.to) + ":" + std::to_string(l.toEnd) + ";");
        }
    }
    return h;
}

// Random graph from build(): repeats, a homopolymer run, 0.5-2% errors, dropouts.
UnitigGraph randomGraph(std::mt19937_64& rng, int k) {
    const size_t glen = 3000 + rng() % 9000;
    std::string genome;
    for (size_t i = 0; i < glen; ++i) genome += "ACGT"[rng() % 4];
    for (int r = 0; r < 4; ++r) {
        const size_t len = 50 + rng() % 400;
        const size_t from = rng() % (genome.size() - len);
        std::string seg = genome.substr(from, len);
        if (rng() % 2) seg = rc(seg);
        if (rng() % 2) seg[rng() % seg.size()] = "ACGT"[rng() % 4];
        genome.insert(rng() % genome.size(), seg);
    }
    genome.insert(rng() % genome.size(), std::string(20 + rng() % 60, "ACGT"[rng() % 4]));
    const double errRate = 0.005 + 0.015 * static_cast<double>(rng() % 100) / 100.0;
    std::unordered_map<std::string, uint32_t> counts;
    const size_t nReads = genome.size() * (8 + rng() % 30) / 150;
    std::uniform_real_distribution<double> U(0, 1);
    const size_t dropPhase = rng() % 5;
    for (size_t r = 0; r < nReads; ++r) {
        const size_t p = rng() % (genome.size() - 150);
        if ((p / 700) % 5 == dropPhase && U(rng) < 0.8) continue;
        std::string read = genome.substr(p, 150);
        for (char& c : read) if (U(rng) < errRate) c = "ACGT"[rng() % 4];
        if (rng() % 2) read = rc(read);
        for (size_t i = 0; i + static_cast<size_t>(k) <= read.size(); ++i) {
            const std::string km = read.substr(i, static_cast<size_t>(k));
            const std::string r2 = rc(km);
            ++counts[km < r2 ? km : r2];
        }
    }
    const uint32_t cutoff = 1 + static_cast<uint32_t>(rng() % 3);
    KmerTable solid(counts.size());
    for (const auto& kv : counts) {
        if (kv.second < cutoff) continue;
        bool ok = false;
        const Kmer x = stringToKmer(kv.first, k, ok);
        if (ok) solid.put(canonical(x, k), kv.second);
    }
    UnitigGraph g = UnitigGraph::build(solid, k, 1);
    g.compact();
    for (uint32_t u = 0; u < g.nodes.size(); ++u)
        if (!g.nodes[u].deleted && rng() % 7 == 0) g.deleteNode(u);
    g.compact();
    return g;
}

// Runs `fn` in a child with stderr captured; returns the captured text.
template <typename Fn>
std::string captureStderrInChild(const char* envName, const char* envValue, Fn&& fn) {
    char path[] = "/tmp/tess_v3_jde.XXXXXX";
    const int fd = mkstemp(path);
    if (fd < 0) return "";
    std::fflush(stdout);
    std::fflush(stderr);
    const pid_t pid = fork();
    if (pid == 0) {
        if (envName) setenv(envName, envValue, 1);
        dup2(fd, 2);
        fn();
        std::fflush(stderr);
        _exit(0);
    }
    int status = 0;
    waitpid(pid, &status, 0);
    close(fd);
    std::ifstream in(path);
    std::stringstream ss;
    ss << in.rdbuf();
    std::remove(path);
    return ss.str();
}

std::string lineWith(const std::string& text, const std::string& tag) {
    std::istringstream in(text);
    std::string line, last;
    while (std::getline(in, line))
        if (line.find(tag) != std::string::npos) last = line;
    return last;
}

}  // namespace

int main(int argc, char** argv) {
    clearTesseractEnv();
    const int trials = argc > 1 ? std::atoi(argv[1]) : 72;

    // The counter-line checks run first, in forked children, before this process has
    // called simplify(): the release caches TESSERACT_JOIN_DEADENDS in a static on first use.
    auto handMade = []() {
        UnitigGraph h;
        h.setK(33);
        std::mt19937_64 rng(11);
        std::string a;
        for (int i = 0; i < 400; ++i) a += "ACGT"[rng() % 4];
        std::string b = a.substr(a.size() - 32);
        for (int i = 0; i < 400; ++i) b += "ACGT"[rng() % 4];
        Unitig ua; ua.seq = a; ua.coverage = 40; h.nodes.push_back(ua);
        Unitig ub; ub.seq = rc(b); ub.coverage = 40; h.nodes.push_back(ub);   // stored reversed
        return std::make_pair(h, a + b.substr(32));
    };
    // 3. Counter line.
    auto builtGraph = []() {
        std::mt19937_64 rng(0x5eed);
        (void)randomGraph(rng, 21);
        return randomGraph(rng, 33);
    };
    {
        const std::string err = captureStderrInChild(nullptr, nullptr, [&]() {
            UnitigGraph g = builtGraph();
            g.simplify(15.0, 150, true, 0.35, 12, nullptr, 0.0);
        });
        const std::string line = lineWith(err, "[joindeadends]");
        std::printf("    default: %s\n", line.c_str());
        check(line.find("k=33 ") != std::string::npos && line.find("requested=31 effective=32") != std::string::npos &&
                  line.find("pairsFound=0 calls=0 joined=0") != std::string::npos &&
                  line.find("roundsSkipped=0 ") == std::string::npos,
              "default: counter line, every round skipped, 0 pairs, 0 calls");
    }
    {
        const std::string err = captureStderrInChild("TESSERACT_JOIN_DEADENDS", "40", [&]() {
            UnitigGraph g = builtGraph();
            g.simplify(15.0, 150, true, 0.35, 12, nullptr, 0.0);
        });
        const std::string line = lineWith(err, "[joindeadends]");
        std::printf("    =40:     %s\n", line.c_str());
        check(line.find("requested=40 effective=0") != std::string::npos,
              "TESSERACT_JOIN_DEADENDS=40 at k=33: reported as not in effect");
    }
    {
        const std::string err = captureStderrInChild("TESSERACT_JOIN_DEADENDS", "0", [&]() {
            UnitigGraph g = builtGraph();
            g.simplify(15.0, 150, true, 0.35, 12, nullptr, 0.0);
        });
        const std::string line = lineWith(err, "[joindeadends]");
        std::printf("    =0:      %s\n", line.c_str());
        check(line.find("requested=0 effective=0 roundsSkipped=0 pairsFound=0 calls=0 joined=0") !=
                  std::string::npos,
              "TESSERACT_JOIN_DEADENDS=0: off, reported as off");
    }
    {
        const std::string err = captureStderrInChild(nullptr, nullptr, [&]() {
            UnitigGraph h = handMade().first;
            h.simplify(40.0, 150, true, 0.35, 12, nullptr, 0.0);
        });
        const std::string line = lineWith(err, "[joindeadends]");
        std::printf("    fallback:%s\n", line.c_str());
        check(line.find("calls=1 joined=1") != std::string::npos && line.find("pairsFound=0") == std::string::npos,
              "fallback: a found pair runs joinDeadEnds, reported");
    }

    // 1. Invariant, plus a digest to compare release and patched builds.
    {
        std::mt19937_64 rng(0x7336);
        const int ks[] = {21, 33, 55, 77, 99, 127};
        size_t graphs = 0, brute = 0, simplifyJoins = 0, direct = 0, rounds = 0;
        uint64_t digest = 1469598103934665603ULL;
        bool valid = true;
        for (int t = 0; t < trials; ++t) {
            const int k = ks[t % 6];
            UnitigGraph g = randomGraph(rng, k);
            if (g.liveCount() == 0) continue;
            ++graphs;
            brute += bruteOverlapPairs(g);
            std::vector<SimplifyRoundStats> rs;
            g.simplify(15.0, 150, false, 0.35, 12, &rs, 0.0);
            for (const auto& st : rs) { simplifyJoins += st.deadEndsJoined; ++rounds; }
            brute += bruteOverlapPairs(g);
            digest = graphDigest(digest, g);
            for (size_t x : {size_t(1), size_t(8), size_t(31), size_t(k - 1), size_t(k), size_t(1000)})
                direct += g.joinDeadEnds(x);
            valid = valid && g.validate().empty();
        }
        std::printf("  invariant: graphs=%zu rounds=%zu brutePairs=%zu simplifyJoins=%zu directJoins=%zu\n",
                    graphs, rounds, brute, simplifyJoins, direct);
        std::printf("  simplify digest (must equal the release build's): %016llx\n",
                    static_cast<unsigned long long>(digest));
        check(graphs > 0 && brute == 0 && simplifyJoins == 0 && direct == 0 && valid,
              "graphs from build(): no exact (k-1) dead-end pair, nothing ever joined");
    }

    // 2. Fallback: an unlinked (k-1) overlap is still joined, exactly as before.
    {
        auto hm = handMade();
        UnitigGraph h = hm.first;
        std::vector<SimplifyRoundStats> rs;
        h.simplify(40.0, 150, false, 0.35, 12, &rs, 0.0);
        const bool one = h.liveCount() == 1;
        std::string merged;
        for (const auto& n : h.nodes) if (!n.deleted) merged = n.seq;
        check(one && (merged == hm.second || merged == rc(hm.second)) && !rs.empty() &&
                  rs.front().deadEndsJoined == 1 && h.validate().empty(),
              "hand-made unlinked (k-1) overlap: still joined and merged by simplify()");
    }

    std::printf("test_v3_graph_joindeadends: %s (%d failure%s)\n", failures ? "FAILED" : "ALL PASSED",
                failures, failures == 1 ? "" : "s");
    return failures ? 1 : 0;
}
