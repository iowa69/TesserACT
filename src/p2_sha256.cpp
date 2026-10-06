#include "p2_sha256.h"

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <vector>

namespace ts {
namespace p2 {

namespace {

constexpr uint32_t kShaK[64] = {
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2};

inline uint32_t rotr(uint32_t x, int n) { return (x >> n) | (x << (32 - n)); }

// MD5 (RFC 1321), only for the armx-compatible md5 of a file's first MiB.
class Md5 {
public:
    void add(const void* p, size_t n) {
        const auto* s = static_cast<const unsigned char*>(p);
        bytes_ += n;
        while (n) {
            const size_t take = std::min(n, static_cast<size_t>(64) - used_);
            std::memcpy(block_ + used_, s, take);
            used_ += take; s += take; n -= take;
            if (used_ == 64) { compress(); used_ = 0; }
        }
    }
    std::string finish() {
        const uint64_t bits = bytes_ * 8;
        const unsigned char one = 0x80, zero = 0;
        add(&one, 1);
        while (used_ != 56) add(&zero, 1);
        unsigned char len[8];
        for (int i = 0; i < 8; ++i) len[i] = static_cast<unsigned char>(bits >> (8 * i));
        add(len, 8);
        char out[33];
        const uint32_t h[4] = {a_, b_, c_, d_};
        for (int i = 0; i < 4; ++i)
            for (int b = 0; b < 4; ++b) std::snprintf(out + 8 * i + 2 * b, 3, "%02x", (h[i] >> (8 * b)) & 0xffu);
        return std::string(out, 32);
    }

private:
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
            M[i] = uint32_t(block_[4 * i]) | (uint32_t(block_[4 * i + 1]) << 8) |
                   (uint32_t(block_[4 * i + 2]) << 16) | (uint32_t(block_[4 * i + 3]) << 24);
        uint32_t A = a_, B = b_, C = c_, D = d_;
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
        a_ += A; b_ += B; c_ += C; d_ += D;
    }
    uint32_t a_ = 0x67452301, b_ = 0xefcdab89, c_ = 0x98badcfe, d_ = 0x10325476;
    uint64_t bytes_ = 0;
    size_t used_ = 0;
    unsigned char block_[64]{};
};

}  // namespace

Sha256::Sha256() {
    static const uint32_t init[8] = {0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
                                     0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19};
    std::memcpy(h_, init, sizeof(h_));
}

void Sha256::block(const uint8_t* p) {
    uint32_t w[64];
    for (int i = 0; i < 16; ++i)
        w[i] = (uint32_t(p[4 * i]) << 24) | (uint32_t(p[4 * i + 1]) << 16) | (uint32_t(p[4 * i + 2]) << 8) |
               uint32_t(p[4 * i + 3]);
    for (int i = 16; i < 64; ++i) {
        const uint32_t s0 = rotr(w[i - 15], 7) ^ rotr(w[i - 15], 18) ^ (w[i - 15] >> 3);
        const uint32_t s1 = rotr(w[i - 2], 17) ^ rotr(w[i - 2], 19) ^ (w[i - 2] >> 10);
        w[i] = w[i - 16] + s0 + w[i - 7] + s1;
    }
    uint32_t a = h_[0], b = h_[1], c = h_[2], d = h_[3], e = h_[4], f = h_[5], g = h_[6], h = h_[7];
    for (int i = 0; i < 64; ++i) {
        const uint32_t S1 = rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25);
        const uint32_t ch = (e & f) ^ (~e & g);
        const uint32_t t1 = h + S1 + ch + kShaK[i] + w[i];
        const uint32_t S0 = rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22);
        const uint32_t mj = (a & b) ^ (a & c) ^ (b & c);
        const uint32_t t2 = S0 + mj;
        h = g; g = f; f = e; e = d + t1; d = c; c = b; b = a; a = t1 + t2;
    }
    h_[0] += a; h_[1] += b; h_[2] += c; h_[3] += d; h_[4] += e; h_[5] += f; h_[6] += g; h_[7] += h;
}

void Sha256::add(const void* data, size_t len) {
    const auto* p = static_cast<const uint8_t*>(data);
    bytes_ += len;
    if (fill_) {
        const size_t take = std::min(len, static_cast<size_t>(64) - fill_);
        std::memcpy(buf_ + fill_, p, take);
        fill_ += take; p += take; len -= take;
        if (fill_ == 64) { block(buf_); fill_ = 0; }
    }
    while (len >= 64) { block(p); p += 64; len -= 64; }
    if (len) { std::memcpy(buf_, p, len); fill_ = len; }
}

std::string Sha256::finish() {
    const uint64_t bits = bytes_ * 8;
    uint8_t pad[72];
    size_t n = 0;
    pad[n++] = 0x80;
    while ((fill_ + n) % 64 != 56) pad[n++] = 0;
    for (int i = 7; i >= 0; --i) pad[n++] = static_cast<uint8_t>(bits >> (8 * i));
    const uint64_t keep = bytes_;
    add(pad, n);
    bytes_ = keep;
    char out[65];
    for (int i = 0; i < 8; ++i) std::snprintf(out + 8 * i, 9, "%08x", h_[i]);
    return std::string(out, 64);
}

std::string sha256Hex(const std::string& data) {
    Sha256 h;
    h.add(data.data(), data.size());
    return h.finish();
}

bool sha256File(const std::string& path, std::string& hex, uint64_t& bytes) {
    hex.clear();
    bytes = 0;
    std::FILE* f = std::fopen(path.c_str(), "rb");
    if (!f) return false;
    Sha256 h;
    std::vector<char> buf(4u << 20);
    size_t n;
    uint64_t total = 0;
    while ((n = std::fread(buf.data(), 1, buf.size(), f)) > 0) {
        h.add(buf.data(), n);
        total += n;
    }
    const bool bad = std::ferror(f) != 0;
    std::fclose(f);
    if (bad) return false;
    hex = h.finish();
    bytes = total;
    return true;
}

std::string md5FileHead(const std::string& path, size_t limit) {
    std::FILE* f = std::fopen(path.c_str(), "rb");
    if (!f) return std::string();
    Md5 h;
    std::vector<char> buf(1u << 20);
    size_t left = limit;
    while (left > 0) {
        const size_t want = std::min(left, buf.size());
        const size_t n = std::fread(buf.data(), 1, want, f);
        if (n == 0) break;
        h.add(buf.data(), n);
        left -= n;
    }
    const bool bad = std::ferror(f) != 0;
    std::fclose(f);
    return bad ? std::string() : h.finish();
}

std::string md5File(const std::string& path) { return md5FileHead(path, SIZE_MAX); }

}  // namespace p2
}  // namespace ts
