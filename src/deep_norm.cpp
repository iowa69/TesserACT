#include "deep_norm.h"

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <string>
#include <thread>

#include "counter.h"
#include "envflags.h"
#include "kmer.h"
#include "seqio.h"

namespace ts {

DeepNormConfig deepNormConfigFromEnv() {
    DeepNormConfig c;
    c.sel = env::real("TESSERACT_DEEP_NORM", 0.0);
    c.target = env::real("TESSERACT_DEEP_NORM_TARGET", 2.0);
    // A Choice with a number token: "both", or a fraction the table has already checked to
    // be a finite number in [0, 1] (strtod, exactly as env::real parses it).
    if (const char* m = env::text("TESSERACT_DEEP_NORM_MATE")) {
        if (std::strcmp(m, "both") == 0) c.both = true;
        else c.mate = std::strtod(m, nullptr);
    }
    return c;
}

DeepNormVerdict deepNormDecide(const DeepNormConfig& c, double peak, float medA, float medB,
                               bool paired, uint64_t unitIndex) {
    const double sel = c.sel * peak;
    // The medians are compared as floats and widened once, as both source binaries did.
    const double hiM = paired ? std::max(medA, medB) : medA;
    const double loM = paired ? std::min(medA, medB) : medA;
    double m;
    if (c.both) {
        // dn7 (sg_dn30b): both mates at the bar; the keep probability from the lower mate.
        if (loM < sel) return DeepNormVerdict::Skip;
        m = loM;
    } else {
        // dev/salfix b3b6d41 (dn8): one mate at the bar, the other at `mate` x the bar; the
        // keep probability from the higher mate.
        if (hiM < sel || loM < c.mate * sel) return DeepNormVerdict::Skip;
        m = hiM;
    }
    const double keep = c.target * peak / m;
    const uint64_t h = mix64(unitIndex ^ 0x5a17f1c3ULL);
    const double u = static_cast<double>(h >> 11) / 9007199254740992.0;
    return u < keep ? DeepNormVerdict::Keep : DeepNormVerdict::Thin;
}

std::vector<float> deepNormMedians(const SequenceStore& reads, const KmerTable& trusted, int k,
                                   int threads) {
    const size_t n = reads.size();
    std::vector<float> med(n, 0.0f);
    auto work = [&](size_t t0, size_t step) {
        std::vector<uint32_t> cs;
        for (size_t r = t0; r < n; r += step) {
            const int len = static_cast<int>(reads.length(r));
            if (len < k) continue;
            cs.clear();
            Kmer f = 0;
            int run = 0;
            for (int i = 0; i < len; ++i) {
                const int c = reads.baseAt(r, static_cast<uint32_t>(i));
                if (c < 0) { run = 0; continue; }
                f = pushBack(f, c, k);
                if (++run >= k && (i % 4) == 0) cs.push_back(trusted.get(canonical(f, k)));
            }
            if (cs.empty()) continue;
            std::nth_element(cs.begin(), cs.begin() + cs.size() / 2, cs.end());
            med[r] = static_cast<float>(cs[cs.size() / 2]);
        }
    };
    const size_t T = static_cast<size_t>(std::max(1, threads));
    std::vector<std::thread> pool;
    for (size_t t = 0; t < T; ++t) pool.emplace_back(work, t, T);
    for (auto& th : pool) th.join();
    return med;
}

DeepNormStats deepNormApply(SequenceStore& reads, const KmerTable& trusted, int k, int threads,
                            double peak, const DeepNormConfig& c) {
    DeepNormStats s;
    s.peak = peak;
    if (!c.enabled() || !(peak > 0)) return s;
    s.ran = true;
    const std::vector<float> med = deepNormMedians(reads, trusted, k, threads);
    const size_t n = reads.size();
    for (size_t r = 0; r < n; ++r) {
        const bool paired = reads.hasMate(r);
        if (paired && (r & 1)) continue;   // decide once per pair, at the even mate
        const size_t mt = paired ? reads.mateOf(r) : r;
        const DeepNormVerdict v = deepNormDecide(c, peak, med[r], med[mt], paired,
                                                 static_cast<uint64_t>(paired ? r / 2 : r));
        if (v == DeepNormVerdict::Skip) continue;
        ++(paired ? s.selectedPairs : s.selectedSingle);
        if (v == DeepNormVerdict::Keep) continue;
        ++(paired ? s.thinnedPairs : s.thinnedSingle);
        reads.maskRange(r, 0, static_cast<uint32_t>(reads.length(r)));
        if (paired) reads.maskRange(mt, 0, static_cast<uint32_t>(reads.length(mt)));
    }
    return s;
}

void deepNormPrintCounts(std::FILE* out, const DeepNormConfig& c, const DeepNormStats& s) {
    std::string rule = "off";
    if (c.enabled()) {
        if (c.both) {
            rule = "both";
        } else {
            char buf[64];
            std::snprintf(buf, sizeof buf, "mate:%.17g", c.mate);
            rule = buf;
        }
    }
    std::fprintf(out, "  [deep_norm_counts] enabled=%d ran=%d rule=%s selected_pairs=%zu "
                      "thinned_pairs=%zu selected_single=%zu thinned_single=%zu masked_reads=%zu\n",
                 c.enabled() ? 1 : 0, s.ran ? 1 : 0, rule.c_str(), s.selectedPairs, s.thinnedPairs,
                 s.selectedSingle, s.thinnedSingle, s.maskedReads());
}

}  // namespace ts
