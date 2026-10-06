// Phase 2, emit-B (F10 provenance): SHA-256 of files and buffers (FIPS 180-4), self-contained.
//
// The provenance block records the SHA-256 of every input file and of the running binary, so an
// output directory can be tied to the exact bytes it was made from. No library outside the C++
// standard library and zlib is linked by TesserACT, so the digest is implemented here; it is
// checked against the FIPS 180-4 / NIST test vectors in tests/test_p2_emitb.cpp.
#pragma once

#include <cstddef>
#include <cstdint>
#include <string>

namespace ts {
namespace p2 {

class Sha256 {
public:
    Sha256();
    void add(const void* data, size_t len);
    // Lower-case hex digest (64 characters). The object must not be used afterwards.
    std::string finish();

private:
    void block(const uint8_t* p);
    uint32_t h_[8];
    uint8_t buf_[64];
    size_t fill_ = 0;
    uint64_t bytes_ = 0;
};

std::string sha256Hex(const std::string& data);

// SHA-256 and size of a whole file, streamed in 4 MiB blocks. Returns false (and leaves the
// outputs empty / 0) when the file cannot be opened or read.
bool sha256File(const std::string& path, std::string& hex, uint64_t& bytes);

// MD5 of the first `limit` bytes of a file (the armx manifest's md5_head1m uses 1 MiB), and
// of a whole file. Empty string when the file cannot be read.
std::string md5FileHead(const std::string& path, size_t limit);
std::string md5File(const std::string& path);

}  // namespace p2
}  // namespace ts
