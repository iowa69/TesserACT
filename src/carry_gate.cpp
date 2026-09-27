#include "carry_gate.h"

#include <algorithm>

#include "kmer.h"

namespace ts {

namespace {

enum : uint8_t { kNotSolid = 0, kReadAbsent = 1, kReadSupported = 2 };

}  // namespace

CarryGateStats gateCarryOnlyJunctions(const std::vector<std::string>& carry, KmerTable& solid,
                                      int k, uint32_t carryWeight, uint32_t cutoff, int mode) {
    CarryGateStats st;
    if (k < 2 || k > kMaxK || carryWeight == 0 || carry.empty() || solid.size() == 0) return st;
    const size_t K = static_cast<size_t>(k);

    // Every carried occurrence, canonical, sorted: occurrences(x) is an equal_range. A sorted
    // vector costs 32 bytes per carried k-mer, a third of a KmerTable at its load factor.
    std::vector<Kmer> carried;
    {
        size_t total = 0;
        for (const std::string& c : carry) {
            if (c.size() >= K) total += c.size() - K + 1;
        }
        carried.reserve(total);
    }
    for (const std::string& c : carry) {
        Kmer fwd = 0, rc = 0;
        int valid = 0;
        for (size_t p = 0; p < c.size(); ++p) {
            const int b = baseCode(c[p]);
            if (b < 0) { valid = 0; fwd = 0; rc = 0; continue; }
            fwd = pushBack(fwd, b, k);
            rc = pushFrontRc(rc, b, k);
            if (++valid >= k) carried.push_back(fwd < rc ? fwd : rc);
        }
    }
    std::sort(carried.begin(), carried.end());

    auto occurrences = [&](const Kmer& x) -> uint64_t {
        const auto r = std::equal_range(carried.begin(), carried.end(), x);
        return static_cast<uint64_t>(r.second - r.first);
    };
    // Reads' share of a solid k-mer's count; 0 when absent from the table.
    auto readCount = [&](const Kmer& canon) -> uint64_t {
        const uint64_t s = solid.get(canon);
        if (s == 0) return 0;
        const uint64_t w = static_cast<uint64_t>(carryWeight) * occurrences(canon);
        return s > w ? s - w : 0;
    };
    const uint64_t competitorFloor = std::max<uint64_t>(2, cutoff);

    std::vector<uint8_t> state;
    std::vector<Kmer> canon;
    std::vector<Kmer> withheld;   // one entry per withheld occurrence

    for (const std::string& c : carry) {
        if (c.size() < K) continue;
        const size_t n = c.size() - K + 1;
        state.assign(n, kNotSolid);
        canon.assign(n, Kmer());
        {
            Kmer fwd = 0, rc = 0;
            int valid = 0;
            for (size_t p = 0; p < c.size(); ++p) {
                const int b = baseCode(c[p]);
                if (b < 0) { valid = 0; fwd = 0; rc = 0; continue; }
                fwd = pushBack(fwd, b, k);
                rc = pushFrontRc(rc, b, k);
                if (++valid < k) continue;
                const size_t i = p + 1 - K;
                const Kmer x = fwd < rc ? fwd : rc;
                canon[i] = x;
                ++st.carriedKmers;
                if (solid.get(x) == 0) continue;            // below the cutoff: already a break
                if (readCount(x) == 0) {
                    state[i] = kReadAbsent;
                    ++st.readAbsent;
                } else {
                    state[i] = kReadSupported;
                }
            }
        }

        for (size_t i = 0; i < n;) {
            if (state[i] != kReadAbsent) { ++i; continue; }
            size_t j = i;
            while (j + 1 < n && state[j + 1] == kReadAbsent) ++j;
            const size_t next = j + 1;
            // Interior and flanked on both sides by read-supported k-mers.
            if (i == 0 || j + 1 >= n || state[i - 1] != kReadSupported ||
                state[j + 1] != kReadSupported) {
                i = next;
                continue;
            }
            ++st.runs;
            const size_t len = j - i + 1;
            if (len > K - 1) {
                ++st.runsTooLong;
                i = next;
                continue;
            }
            // Left flank: the k-mer at i-1, read along the contig. Does some successor other
            // than the carried one have read support?
            bool left = false, right = false;
            {
                bool ok = false;
                const Kmer x = stringToKmer(c.substr(i - 1, K), k, ok);
                const int carriedNext = baseCode(c[i - 1 + K]);
                for (int b = 0; ok && b < 4 && !left; ++b) {
                    if (b == carriedNext) continue;
                    if (readCount(canonical(pushBack(x, b, k), k)) >= competitorFloor) left = true;
                }
            }
            // Right flank: the k-mer at j+1. Does some predecessor other than the carried one?
            {
                bool ok = false;
                const Kmer z = stringToKmer(c.substr(j + 1, K), k, ok);
                const int carriedPrev = baseCode(c[j]);
                for (int b = 0; ok && b < 4 && !right; ++b) {
                    if (b == carriedPrev) continue;
                    if (readCount(canonical(pushFront(z, b, k), k)) >= competitorFloor) right = true;
                }
            }
            const bool fire = mode >= 2 ? (left && right) : (left || right);
            if (!fire) {
                ++st.runsNoCompetitor;
                i = next;
                continue;
            }
            ++st.candidates;
            if (mode > 0) {
                ++st.runsGated;
                for (size_t p = i; p <= j; ++p) withheld.push_back(canon[p]);
            }
            i = next;
        }
    }

    for (const Kmer& x : withheld) {
        const uint32_t s = solid.get(x);
        if (s == 0) continue;
        const uint32_t reduced = s > carryWeight ? s - carryWeight : 0;
        if (reduced < cutoff || reduced == 0) {
            solid.erase(x);
            ++st.kmersWithheld;
        } else {
            solid.put(x, reduced);
        }
    }
    return st;
}

}  // namespace ts
