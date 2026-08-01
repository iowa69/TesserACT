// K-mer encoding: two bits per base packed into a little-endian multiword
// integer. Four 64-bit words give 256 bits, so k can reach 128.
//
// Three words (k <= 96) was enough for 150 bp reads, but on 2x250 MiSeq
// libraries -- which is what most of the closed-reference panel is -- the
// unitig N50 was still climbing steeply at the top of the ladder, so k, not
// the data, was the binding constraint. The fourth word costs one more
// comparison per k-mer operation and eight more bytes per table entry, and
// buys every repeat between 96 and 128 bases.
//
// Words are ordered least-significant first, and the most recently appended
// base occupies the lowest bits. Because 64 is even, no 2-bit group ever
// straddles a word boundary, which keeps every operation here branch-light.
#pragma once

#include <cstdint>
#include <cstring>
#include <string>

namespace ts {

constexpr int kKmerWords = 4;
constexpr int kKmerBits = 64 * kKmerWords;
constexpr int kMaxK = kKmerBits / 2;   // 96

struct Kmer {
    uint64_t w[kKmerWords];

    Kmer() : w{} {}
    // Implicit so `Kmer x = 0;` and rolling-window resets read naturally; the
    // value is taken as the low word, matching the packed layout.
    Kmer(uint64_t v) : w{v} {}

    bool operator==(const Kmer& o) const {
        for (int i = 0; i < kKmerWords; ++i) {
            if (w[i] != o.w[i]) return false;
        }
        return true;
    }
    bool operator!=(const Kmer& o) const { return !(*this == o); }
    bool operator<(const Kmer& o) const {
        for (int i = kKmerWords - 1; i >= 0; --i) {
            if (w[i] != o.w[i]) return w[i] < o.w[i];
        }
        return false;
    }
    bool operator<=(const Kmer& o) const { return !(o < *this); }
};

// A=0, C=1, G=2, T=3; anything else is invalid and breaks the k-mer run.
inline int baseCode(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return -1;
    }
}

inline char codeBase(int c) { return "ACGT"[c & 3]; }

inline Kmer shiftLeft2(Kmer a) {
    Kmer r;
    for (int i = kKmerWords - 1; i > 0; --i) r.w[i] = (a.w[i] << 2) | (a.w[i - 1] >> 62);
    r.w[0] = a.w[0] << 2;
    return r;
}

inline Kmer shiftRight2(Kmer a) {
    Kmer r;
    for (int i = 0; i < kKmerWords - 1; ++i) r.w[i] = (a.w[i] >> 2) | (a.w[i + 1] << 62);
    r.w[kKmerWords - 1] = a.w[kKmerWords - 1] >> 2;
    return r;
}

inline Kmer shiftRightN(Kmer a, int n) {
    Kmer r;
    if (n >= kKmerBits || n < 0) return r;
    const int words = n / 64;
    const int bits = n % 64;
    for (int i = 0; i < kKmerWords; ++i) {
        const int src = i + words;
        uint64_t v = (src < kKmerWords) ? a.w[src] : 0;
        if (bits) {
            v >>= bits;
            const uint64_t hi = (src + 1 < kKmerWords) ? a.w[src + 1] : 0;
            v |= hi << (64 - bits);
        }
        r.w[i] = v;
    }
    return r;
}

// A k-mer holding `code` at base position `pos`, counting from the low end.
inline Kmer baseAtPosition(int code, int pos) {
    Kmer r;
    const int bit = 2 * pos;
    if (bit >= 0 && bit < kKmerBits) {
        r.w[bit / 64] = static_cast<uint64_t>(code & 3) << (bit % 64);
    }
    return r;
}

inline Kmer bitOr(Kmer a, Kmer b) {
    Kmer r;
    for (int i = 0; i < kKmerWords; ++i) r.w[i] = a.w[i] | b.w[i];
    return r;
}

inline void maskToK(Kmer& a, int k) {
    const int bits = 2 * k;
    for (int i = 0; i < kKmerWords; ++i) {
        const int lo = i * 64;
        if (bits <= lo) a.w[i] = 0;
        else if (bits < lo + 64) a.w[i] &= (uint64_t{1} << (bits - lo)) - 1;
    }
}

inline int lowBase(Kmer a) { return static_cast<int>(a.w[0] & 3); }

// Append a base on the right, dropping the leftmost base.
inline Kmer pushBack(Kmer km, int code, int k) {
    Kmer r = shiftLeft2(km);
    r.w[0] |= static_cast<uint64_t>(code & 3);
    maskToK(r, k);
    return r;
}

// Append a base on the left, dropping the rightmost base.
inline Kmer pushFront(Kmer km, int code, int k) {
    return bitOr(shiftRight2(km), baseAtPosition(code, k - 1));
}

// Append the complement of a base on the left of the reverse-complement strand.
inline Kmer pushFrontRc(Kmer rc, int code, int k) {
    return bitOr(shiftRight2(rc), baseAtPosition(3 - (code & 3), k - 1));
}

inline uint64_t reverse2BitGroups64(uint64_t v) {
    v = ((v >> 2) & 0x3333333333333333ULL) | ((v & 0x3333333333333333ULL) << 2);
    v = ((v >> 4) & 0x0F0F0F0F0F0F0F0FULL) | ((v & 0x0F0F0F0F0F0F0F0FULL) << 4);
    v = ((v >> 8) & 0x00FF00FF00FF00FFULL) | ((v & 0x00FF00FF00FF00FFULL) << 8);
    v = ((v >> 16) & 0x0000FFFF0000FFFFULL) | ((v & 0x0000FFFF0000FFFFULL) << 16);
    return (v >> 32) | (v << 32);
}

inline Kmer reverseComplement(Kmer km, int k) {
    // Complementing a 2-bit code is a bitwise NOT; reversing the base order is
    // a reversal of 2-bit groups across the whole width, after which the result
    // is shifted down so the k-mer occupies the low 2k bits again.
    Kmer c;
    for (int i = 0; i < kKmerWords; ++i) c.w[i] = ~km.w[i];
    Kmer rev;
    for (int i = 0; i < kKmerWords; ++i) {
        rev.w[i] = reverse2BitGroups64(c.w[kKmerWords - 1 - i]);
    }
    return shiftRightN(rev, kKmerBits - 2 * k);
}

inline Kmer canonical(Kmer km, int k) {
    const Kmer rc = reverseComplement(km, k);
    return km < rc ? km : rc;
}

inline bool isCanonical(Kmer km, int k) { return km <= reverseComplement(km, k); }

inline std::string kmerToString(Kmer km, int k) {
    std::string s(static_cast<size_t>(k), 'A');
    for (int i = k - 1; i >= 0; --i) {
        s[static_cast<size_t>(i)] = codeBase(lowBase(km));
        km = shiftRight2(km);
    }
    return s;
}

inline Kmer stringToKmer(const std::string& s, int k, bool& ok) {
    Kmer km;
    ok = true;
    for (int i = 0; i < k; ++i) {
        const int c = baseCode(s[static_cast<size_t>(i)]);
        if (c < 0) { ok = false; return Kmer(); }
        km = shiftLeft2(km);
        km.w[0] |= static_cast<uint64_t>(c);
    }
    return km;
}

inline uint64_t mix64(uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

inline uint64_t kmerHash(Kmer km) {
    uint64_t h = 0;
    for (int i = kKmerWords - 1; i >= 0; --i) h = mix64(km.w[i] ^ h);
    return h;
}

// Standard containers keyed by Kmer need this explicitly.
struct KmerHasher {
    size_t operator()(const Kmer& km) const noexcept {
        return static_cast<size_t>(kmerHash(km));
    }
};

// The (k-1)-mer prefix and suffix, used as the graph's node identities.
inline Kmer prefixOf(Kmer km, int) { return shiftRight2(km); }
inline Kmer suffixOf(Kmer km, int k) {
    Kmer r = km;
    maskToK(r, k - 1);
    return r;
}

}  // namespace ts
