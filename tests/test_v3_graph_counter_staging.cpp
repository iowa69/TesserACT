// T39 (build_v3, G-graph): KmerCounter staging buffers must not grow with T^2.
//
// Each counting worker used to stage nShards * 512 k-mers (32 bytes each), with
// nShards = nextpow2(max(256, 8T)), allocated and touched before the memory guard could
// run: 16 MiB in total at T=4, 512 MiB at T=64, 2 GiB at T=128. build_v3 caps the shard
// count and holds each worker's staging at 4 MiB. Counts do not depend on sharding.
//
//  1. Counts at T=64 are identical to T=4 (solid-table digest, histogram, cutoff).
//  2. With the memory limit set 384 MiB above the counter's own footprint, a T=64 pass
//     over enough reads to occupy every worker stays under it. Release 1.3.0 stages
//     512 MiB there and trips the guard (the assembler then aborts the run).
// Uses only release APIs, so the same source fails on release 1.3.0.
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <random>
#include <string>
#include <sys/wait.h>
#include <unistd.h>
#include <vector>

#include "counter.h"
#include "seqio.h"
#include "util.h"

extern char** environ;

namespace {

int failures = 0;
void check(bool ok, const std::string& what) {
    std::printf("  %s  %s\n", ok ? "PASS" : "FAIL", what.c_str());
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

struct Digest {
    uint64_t sum = 0, xr = 0;
    size_t size = 0;
    uint32_t cutoff = 0;
    uint64_t total = 0, distinct = 0;
    std::vector<uint64_t> hist;
    bool operator==(const Digest& o) const {
        return sum == o.sum && xr == o.xr && size == o.size && cutoff == o.cutoff && total == o.total &&
               distinct == o.distinct && hist == o.hist;
    }
};

Digest countAt(const ts::SequenceStore& store, int threads, long long limitAbove, bool* exceeded) {
    ts::KmerCounter c(21, threads);
    if (limitAbove > 0) {
        const long long base = ts::util::currentMemoryBytes();
        c.setMemoryLimit(base + limitAbove);
    }
    c.count(store, {}, 0);
    if (exceeded) *exceeded = c.exceededMemory();
    ts::KmerTable solid;
    c.extractSolid(2, solid);
    Digest d;
    solid.forEach([&](ts::Kmer k, uint32_t n) {
        const uint64_t h = ts::kmerHash(k) * 0x9e3779b97f4a7c15ULL + n;
        d.sum += h;
        d.xr ^= h;
    });
    d.size = solid.size();
    d.cutoff = c.stats().cutoff;
    d.total = c.stats().totalKmers;
    d.distinct = c.stats().distinctKmers;
    d.hist = c.stats().histogram;
    return d;
}

}  // namespace

int main() {
    clearTesseractEnv();
    char tmpl[] = "/tmp/tess_v3_staging.XXXXXX";
    const char* dir = mkdtemp(tmpl);
    if (!dir) return 2;
    const std::string path = std::string(dir) + "/reads.fa";
    {
        // 1536 blocks of 1024 reads: each of 64 workers takes ~24 blocks and re-checks the
        // guard on each, so the later checks see every worker's staging allocated. Short
        // reads from a small genome keep the shard tables tiny, so staging is what the
        // memory figure measures.
        std::mt19937 rng(39);
        std::string genome(200000, 'A');
        for (char& c : genome) c = "ACGT"[rng() & 3];
        std::ofstream o(path);
        for (int i = 0; i < 1536 * 1024; ++i) {
            const size_t p = rng() % (genome.size() - 40);
            o << ">r" << i << "\n" << genome.substr(p, 40) << "\n";
        }
    }
    ts::SequenceStore store;
    {
        ts::Library lib;
        lib.r1 = path;
        std::string err;
        ts::QualityTrim qt;
        qt.enabled = false;
        store.setQualityTrim(qt);
        if (!store.load({lib}, 1, err)) { std::printf("load: %s\n", err.c_str()); return 2; }
    }

    // Memory first, each in a fresh child: a freed staging area stays resident in this
    // process's malloc arenas and would hide the next pass's growth from the guard.
    auto limitedPass = [&](int threads) -> int {
        std::fflush(stdout);
        const pid_t pid = fork();
        if (pid == 0) {
            bool exceeded = true;
            (void)countAt(store, threads, 384LL << 20, &exceeded);
            std::printf("    T=%d limited pass: peak RSS %.0f MiB, exceeded=%d\n", threads,
                        static_cast<double>(ts::util::peakMemoryBytes()) / 1048576.0, exceeded ? 1 : 0);
            std::fflush(stdout);
            _exit(exceeded ? 1 : 0);
        }
        int status = 0;
        waitpid(pid, &status, 0);
        return WIFEXITED(status) ? WEXITSTATUS(status) : 2;
    };
    check(limitedPass(64) == 0, "T=64 counting stays within 384 MiB above the counter's footprint");
    check(limitedPass(4) == 0, "control: T=4 within the same limit");

    const Digest d4 = countAt(store, 4, 0, nullptr);
    const Digest d64 = countAt(store, 64, 0, nullptr);
    std::printf("  T=4  solid=%zu cutoff=%u total=%llu | T=64 solid=%zu cutoff=%u total=%llu\n", d4.size,
                d4.cutoff, static_cast<unsigned long long>(d4.total), d64.size, d64.cutoff,
                static_cast<unsigned long long>(d64.total));
    check(d4 == d64 && d4.size > 0, "counts at T=64 identical to T=4");

    std::printf("test_v3_graph_counter_staging: %s (%d failure%s)\n", failures ? "FAILED" : "ALL PASSED",
                failures, failures == 1 ? "" : "s");
    std::remove(path.c_str());
    rmdir(dir);
    return failures ? 1 : 0;
}
