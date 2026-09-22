// Standalone targeted regression test; see work/agent_extension.md for build.
#include "correct.h"
#include "seqio.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <random>
#include <string>
#include <vector>

namespace {
int checks = 0;
void check(bool ok, const char* what) {
    ++checks;
    if (!ok) { std::fprintf(stderr, "FAIL: %s\n", what); std::exit(1); }
}
std::mt19937 rng(2918);
std::string randomSeq(size_t n) {
    std::string s(n, 'A');
    for (char& c : s) c = "ACGT"[rng() % 4];
    return s;
}
std::string rc(std::string s) {
    std::reverse(s.begin(), s.end());
    for (char& c : s) { const int code = ts::baseCode(c); c = code < 0 ? 'N' : ts::codeBase(3 - code); }
    return s;
}
std::string directory;
ts::SequenceStore load(const std::string& a, const std::string& b,
                       const std::string& single = "") {
    std::ofstream(directory + "/a.fa") << ">pair/1\n" << a << '\n';
    std::ofstream(directory + "/b.fa") << ">pair/2\n" << b << '\n';
    ts::Library lib; lib.r1 = directory + "/a.fa"; lib.r2 = directory + "/b.fa";
    std::vector<ts::Library> libs{lib};
    if (!single.empty()) {
        std::ofstream(directory + "/s.fa") << ">single\n" << single << '\n';
        ts::Library sl; sl.r1 = directory + "/s.fa"; libs.push_back(sl);
    }
    ts::SequenceStore reads;
    std::string error;
    check(reads.load(libs, 2, error), "load fixture");
    return reads;
}
size_t kmers(const ts::SequenceStore& reads, size_t r, int k) {
    size_t n = 0;
    ts::forEachKmer(reads, r, k, [&](ts::Kmer, uint32_t) { ++n; });
    return n;
}
}

int main() {
    char templ[] = "/tmp/tesseract-mate-rescue-XXXXXX";
    const char* tmp = mkdtemp(templ);
    check(tmp != nullptr, "temporary directory"); directory = tmp;
    const std::string fragment = randomSeq(400);
    const std::string a = fragment.substr(0, 300), b = rc(fragment.substr(100, 300));
    {
        auto reads = load(a, b, randomSeq(300));
        const size_t bases = reads.totalBases(), size = reads.size();
        reads.maskRange(0, 200, 300);
        const size_t before = kmers(reads, 0, 127);
        const auto stats = ts::rescueMateOverlaps(reads, {0, 0, 2, UINT32_MAX});
        check(stats.basesRescued == 100 && stats.readsRescued == 1, "recover true masked overlap");
        check(reads.decode(0) == a && reads.decode(1) == b, "exact original sequences recovered");
        check(kmers(reads, 0, 127) == before + 100, "recover 100 missing high-k windows");
        check(reads.size() == size && reads.totalBases() == bases && reads.pairCount() == 1 &&
              reads.mateOf(0) == 1 && reads.length(0) == 300, "preserve reads, depth and pair layout");
    }
    {
        auto reads = load(a, b);
        reads.maskRange(1, 200, 300);
        const auto stats = ts::rescueMateOverlaps(reads, {1});
        check(stats.basesRescued == 100 && reads.decode(1) == b, "reverse mate rescue uses correct coordinates");
    }
    {
        auto reads = load(a, b);
        reads.setBase(0, 220, (reads.baseAt(0, 220) + 1) % 4);
        reads.maskRange(0, 200, 300);
        const auto stats = ts::rescueMateOverlaps(reads, {0});
        check(stats.basesRescued == 99 && reads.baseAt(0, 220) < 0, "masked sequencing disagreement stays masked");
    }
    {
        auto reads = load(a, b);
        reads.maskRange(0, 200, 260);
        // R2 positions 140..199 correspond to R1 positions 200..259.
        reads.maskRange(1, 140, 200);
        const auto stats = ts::rescueMateOverlaps(reads, {0, 1});
        check(stats.basesRescued == 0 && reads.baseAt(0, 220) < 0, "two masked mates cannot corroborate each other");
    }
    {
        auto reads = load(a, randomSeq(300));
        reads.maskRange(0, 200, 300);
        check(ts::rescueMateOverlaps(reads, {0}).basesRescued == 0, "unrelated mates do not rescue");
    }
    {
        auto reads = load(std::string(300, 'A'), std::string(300, 'T'));
        reads.maskRange(0, 200, 300);
        check(ts::rescueMateOverlaps(reads, {0}).basesRescued == 0, "repetitive offsets cannot rescue");
    }
    {
        std::string mosaic = randomSeq(300);
        mosaic.replace(0, 50, a.substr(100, 50));
        mosaic.replace(100, 50, a.substr(250, 50));
        auto reads = load(a, rc(mosaic));
        reads.maskRange(0, 0, 100);
        reads.maskRange(0, 150, 250);
        reads.maskRange(1, 0, 150);
        reads.maskRange(1, 200, 250);
        const auto stats = ts::rescueMateOverlaps(reads, {0, 1});
        check(stats.ambiguousOverlaps == 1 && stats.basesRescued == 0,
              "two distinct supported offsets reject the whole overlap");
    }
    {
        auto reads = load(a, b);
        reads.maskRange(0, 120, 300);
        check(ts::rescueMateOverlaps(reads, {0}).basesRescued == 0, "insufficient independent overlap evidence");
    }
    {
        std::string ambiguous = a; ambiguous[220] = 'N';
        auto reads = load(ambiguous, b);
        reads.maskRange(1, 240, 260);
        const auto stats = ts::rescueMateOverlaps(reads, {1});
        check(reads.baseAt(0, 220) < 0, "input N is never considered eligible through raw A");
        check(stats.basesRescued == 20, "input N elsewhere does not erase independent evidence");
    }
    {
        auto reads = load(a, b);
        ts::KmerTable solid(400);
        bool valid = false;
        for (size_t p = 0; p + 21 <= 200; ++p)
            solid.put(ts::canonical(ts::stringToKmer(a.substr(p, 21), 21, valid), 21), 10);
        std::vector<uint32_t> masked;
        const auto correction = ts::correctReads(reads, solid, 21, 2, 8, &masked);
        check(correction.basesMasked > 0 && !masked.empty(), "corrector exports actual mask provenance");
        std::string withN = a; withN[220] = 'N';
        auto ambiguousReads = load(withN, b);
        ts::correctReads(ambiguousReads, solid, 21, 2, 8, &masked);
        check(std::find(masked.begin(), masked.end(), 0) == masked.end(), "corrector excludes input-N read from provenance");
        ts::rescueMateOverlaps(ambiguousReads, masked);
        check(ambiguousReads.baseAt(0, 220) < 0, "input N preserved through complete provenance/rescue flow");
    }
    std::filesystem::remove_all(directory);
    std::printf("mate rescue: %d checks passed\n", checks);
}
