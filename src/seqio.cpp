#include "seqio.h"

#include <zlib.h>

#include <algorithm>
#include <atomic>
#include <cerrno>
#include <cstdio>
#include <cstring>
#include <thread>

namespace ts {
namespace {

// Large block reads with an in-buffer newline scan. gzgets() is roughly an
// order of magnitude slower because it copies byte by byte through the zlib
// state machine.
constexpr size_t kIoBuffer = 1u << 22;    // 4 MiB
constexpr size_t kZlibBuffer = 1u << 20;  // 1 MiB internal inflate buffer
constexpr size_t kWriteBuffer = 1u << 20;

// Table form of baseCode(), so the packing loop is a load instead of a switch.
struct BaseLut {
    int8_t v[256];
    BaseLut() {
        for (int i = 0; i < 256; ++i) v[i] = static_cast<int8_t>(baseCode(static_cast<char>(i)));
    }
};
const BaseLut kLut;

uint64_t nameHashOf(const char* s, size_t n) {
    size_t e = 0;
    while (e < n && s[e] != ' ' && s[e] != '\t') ++e;
    // Mate suffixes differ between the two files of a pair by design.
    if (e >= 2 && s[e - 2] == '/' && (s[e - 1] == '1' || s[e - 1] == '2')) e -= 2;
    uint64_t h = 1469598103934665603ULL;
    for (size_t i = 0; i < e; ++i) {
        h ^= static_cast<unsigned char>(s[i]);
        h *= 1099511628211ULL;
    }
    return h;
}

std::string trimName(const std::string& s) {
    size_t e = s.find_first_of(" \t");
    return e == std::string::npos ? s : s.substr(0, e);
}

class GzReader {
public:
    GzReader() = default;
    ~GzReader() { close(); }
    GzReader(const GzReader&) = delete;
    GzReader& operator=(const GzReader&) = delete;

    bool open(const std::string& path, std::string& error) {
        gz_ = gzopen(path.c_str(), "rb");
        if (!gz_) {
            error = "cannot open '" + path + "': " + std::strerror(errno);
            return false;
        }
        gzbuffer(gz_, static_cast<unsigned>(kZlibBuffer));
        buf_.resize(kIoBuffer);
        pos_ = end_ = 0;
        eof_ = failed_ = false;
        return true;
    }

    void close() {
        if (gz_) {
            gzclose(gz_);
            gz_ = nullptr;
        }
    }

    bool failed() const { return failed_; }

    std::string ioMessage() const {
        if (!gz_) return "read error";
        int code = 0;
        const char* m = gzerror(gz_, &code);
        return (m && *m) ? std::string(m) : std::string("read error");
    }

    // `out` points into the internal buffer and is invalidated by the next call.
    bool nextLine(const char*& out, size_t& len) {
        for (;;) {
            if (pos_ < end_) {
                char* base = buf_.data();
                const void* nl = std::memchr(base + pos_, '\n', end_ - pos_);
                if (nl) {
                    size_t idx = static_cast<size_t>(static_cast<const char*>(nl) - base);
                    out = base + pos_;
                    len = idx - pos_;
                    if (len && out[len - 1] == '\r') --len;
                    pos_ = idx + 1;
                    return true;
                }
            }
            if (eof_) {
                if (pos_ < end_) {
                    out = buf_.data() + pos_;
                    len = end_ - pos_;
                    if (len && out[len - 1] == '\r') --len;
                    pos_ = end_;
                    return true;
                }
                return false;
            }
            if (pos_ > 0) {
                std::memmove(buf_.data(), buf_.data() + pos_, end_ - pos_);
                end_ -= pos_;
                pos_ = 0;
            }
            if (end_ == buf_.size()) buf_.resize(buf_.size() * 2);  // single very long FASTA line
            int got = gzread(gz_, buf_.data() + end_, static_cast<unsigned>(buf_.size() - end_));
            if (got < 0) {
                failed_ = true;
                eof_ = true;
            } else if (got == 0) {
                eof_ = true;
            } else {
                end_ += static_cast<size_t>(got);
            }
        }
    }

private:
    gzFile gz_ = nullptr;
    std::vector<char> buf_;
    size_t pos_ = 0;
    size_t end_ = 0;
    bool eof_ = false;
    bool failed_ = false;
};

// Walks a FASTA or FASTQ stream (format detected from the first record) and
// calls fn(name, seq, seqLen, qual) once per record. `seq` and `qual` are only
// valid inside the callback, and `qual` is null for FASTA. fn returns false to
// abort, having already filled `error`.
template <typename Fn>
bool forEachRecord(GzReader& r, const std::string& path, std::string& error, Fn&& fn) {
    const char* line = nullptr;
    size_t len = 0;
    std::string name;
    std::string seqBuf;
    name.reserve(128);
    seqBuf.reserve(1024);

    do {
        if (!r.nextLine(line, len)) {
            if (r.failed()) {
                error = "read error in '" + path + "': " + r.ioMessage();
                return false;
            }
            return true;  // empty file
        }
    } while (len == 0);

    if (line[0] == '@') {
        for (;;) {
            if (len == 0 || line[0] != '@') {
                error = "malformed FASTQ in '" + path + "': expected a '@' header line";
                return false;
            }
            name.assign(line + 1, len - 1);
            if (!r.nextLine(line, len)) {
                error = "truncated FASTQ record '" + trimName(name) + "' in '" + path + "'";
                return false;
            }
            // Copied because reading the quality line invalidates `line`.
            seqBuf.assign(line, len);
            const size_t seqLen = seqBuf.size();
            if (!r.nextLine(line, len)) {
                error = "truncated FASTQ record '" + trimName(name) + "' in '" + path + "'";
                return false;
            }
            if (len == 0 || line[0] != '+') {
                error = "malformed FASTQ record '" + trimName(name) + "' in '" + path +
                        "': expected a '+' separator";
                return false;
            }
            if (!r.nextLine(line, len)) {
                error = "truncated FASTQ record '" + trimName(name) + "' in '" + path + "'";
                return false;
            }
            if (len != seqLen) {
                error = "malformed FASTQ record '" + trimName(name) + "' in '" + path +
                        "': quality length differs from sequence length";
                return false;
            }
            if (!fn(name, seqBuf.data(), seqLen, line)) return false;
            if (!r.nextLine(line, len)) break;
        }
    } else if (line[0] == '>') {
        bool more = true;
        while (more) {
            name.assign(line + 1, len - 1);
            seqBuf.clear();
            for (;;) {
                more = r.nextLine(line, len);
                if (!more || (len > 0 && line[0] == '>')) break;
                seqBuf.append(line, len);
            }
            if (!fn(name, seqBuf.data(), seqBuf.size(), nullptr)) return false;
        }
    } else {
        error = "unrecognised format in '" + path + "': expected FASTA ('>') or FASTQ ('@')";
        return false;
    }

    if (r.failed()) {
        error = "read error in '" + path + "': " + r.ioMessage();
        return false;
    }
    return true;
}

}  // namespace

// Length of the read once its 3' end has been cut back to where quality holds
// up. Illumina reads decay towards the 3' end, and on 2x250/2x300 MiSeq runs
// the last stretch can fall to Q3 -- a ~50% error rate, i.e. noise. Those bases
// are worse than useless: the corrector will happily rewrite them into
// plausible genomic sequence, and the fabrication enters the graph as trusted
// evidence. Cutting them before anything is counted is the only place the fix
// is cheap, because the trusted k-mer set is built from this same pass.
//
// The window is walked from the 3' end, the same rule fastp's --cut_tail uses,
// so the two tools agree on where a read stops being trustworthy.
uint32_t qualityTrimmedLength(const char* qual, size_t len, const QualityTrim& qt) {
    if (!qt.enabled || qual == nullptr || len == 0) return static_cast<uint32_t>(len);
    const int w = qt.windowSize;
    if (w <= 0 || static_cast<size_t>(w) > len) return static_cast<uint32_t>(len);

    const int threshold = qt.meanQuality;
    size_t end = len;
    int sum = 0;
    for (size_t i = end - static_cast<size_t>(w); i < end; ++i) {
        sum += qual[i] - qt.phredOffset;
    }
    // Shrink one base at a time; the window slides with the end, so the sum is
    // maintained rather than recomputed.
    while (end > static_cast<size_t>(w)) {
        if (sum >= threshold * w) break;
        sum -= qual[end - 1] - qt.phredOffset;
        sum += qual[end - 1 - static_cast<size_t>(w)] - qt.phredOffset;
        --end;
    }
    if (end == static_cast<size_t>(w) && sum < threshold * w) end = 0;
    return static_cast<uint32_t>(end);
}

namespace {

// Error path only: recovers a record's name so the message can point at it.
std::string recordNameAt(const std::string& path, size_t index) {
    GzReader r;
    std::string err;
    if (!r.open(path, err)) return "<unreadable>";
    std::string found = "<missing>";
    size_t i = 0;
    forEachRecord(r, path, err, [&](const std::string& name, const char*, size_t, const char*) {
        if (i++ == index) {
            found = trimName(name);
            return false;
        }
        return true;
    });
    return found;
}

// Writes one read into the shared 2-bit buffer. Reads are packed back to back,
// so the first and last word of a read can be shared with a neighbouring read
// that another thread owns; those two get an atomic OR while the fully owned
// interior words get plain stores. Everything starts zeroed, so OR is enough.
void packSequence(uint64_t* data, uint64_t* amb, uint64_t startBit, const char* seq, uint32_t len) {
    uint32_t p = 0;
    while (p < len) {
        const uint64_t bit = startBit + p;
        const size_t w = static_cast<size_t>(bit >> 5);
        const unsigned off = static_cast<unsigned>(bit & 31);
        unsigned take = 32 - off;
        if (take > len - p) take = len - p;
        uint64_t acc = 0;
        for (unsigned t = 0; t < take; ++t) {
            int c = kLut.v[static_cast<unsigned char>(seq[p + t])];
            if (c < 0) {
                const uint64_t ab = bit + t;
                __atomic_fetch_or(&amb[ab >> 6], 1ULL << (ab & 63), __ATOMIC_RELAXED);
                c = 0;
            }
            acc |= static_cast<uint64_t>(c) << ((off + t) * 2);
        }
        if (off == 0 && take == 32) {
            data[w] = acc;
        } else {
            __atomic_fetch_or(&data[w], acc, __ATOMIC_RELAXED);
        }
        p += take;
    }
}

template <typename Fn>
void runParallel(size_t nTasks, int threads, Fn&& fn) {
    if (nTasks == 0) return;
    size_t nt = threads > 0 ? static_cast<size_t>(threads) : 1;
    if (nt > nTasks) nt = nTasks;
    if (nt <= 1) {
        for (size_t i = 0; i < nTasks; ++i) fn(i);
        return;
    }
    std::atomic<size_t> next{0};
    std::vector<std::thread> pool;
    pool.reserve(nt);
    for (size_t t = 0; t < nt; ++t) {
        pool.emplace_back([&] {
            for (;;) {
                size_t i = next.fetch_add(1, std::memory_order_relaxed);
                if (i >= nTasks) break;
                fn(i);
            }
        });
    }
    for (auto& th : pool) th.join();
}

// One input file. Records of a file always land on an arithmetic progression of
// global read indices, which is what lets the second pass run file-parallel:
// stride 2 for the two members of a pair, stride 1 otherwise.
struct FileSlot {
    std::string path;
    bool needNames = false;
    size_t base = 0;
    size_t stride = 1;
    size_t count = 0;
    uint64_t bases = 0;
    uint64_t trimmedBases = 0;
    uint32_t maxLen = 0;
    std::vector<uint32_t> lengths;
    std::vector<uint64_t> nameHashes;
};

enum LibMode { kSingle, kPaired, kInterleaved };

struct LibPlan {
    size_t slotA = 0;
    size_t slotB = 0;
    LibMode mode = kSingle;
};

}  // namespace

bool SequenceStore::load(const std::vector<Library>& libs, int threads, std::string& error) {
    error.clear();
    data_.clear();
    ambiguous_.clear();
    offsets_.clear();
    totalBases_ = 0;
    maxLen_ = 0;
    paired_ = false;
    pairedReads_ = 0;
    trimmedBases_ = 0;

    std::vector<FileSlot> slots;
    std::vector<LibPlan> plans;
    slots.reserve(libs.size() * 2);
    plans.reserve(libs.size());

    for (const Library& lib : libs) {
        if (lib.r1.empty()) {
            error = "library has no r1 file";
            return false;
        }
        LibPlan plan;
        plan.slotA = slots.size();
        if (!lib.r2.empty()) {
            plan.mode = kPaired;
            slots.emplace_back();
            slots.back().path = lib.r1;
            slots.back().needNames = true;
            plan.slotB = slots.size();
            slots.emplace_back();
            slots.back().path = lib.r2;
            slots.back().needNames = true;
        } else if (lib.interleaved) {
            plan.mode = kInterleaved;
            plan.slotB = plan.slotA;
            slots.emplace_back();
            slots.back().path = lib.r1;
            slots.back().needNames = true;
        } else {
            plan.mode = kSingle;
            plan.slotB = plan.slotA;
            slots.emplace_back();
            slots.back().path = lib.r1;
        }
        plans.push_back(plan);
    }

    // Pass 1: count records and bases so every buffer is allocated exactly once.
    std::vector<std::string> slotErr(slots.size());
    runParallel(slots.size(), threads, [&](size_t s) {
        FileSlot& slot = slots[s];
        GzReader r;
        if (!r.open(slot.path, slotErr[s])) return;
        forEachRecord(r, slot.path, slotErr[s],
                      [&](const std::string& name, const char*, size_t rawLen, const char* qual) {
                          if (rawLen > 0xFFFFFFFFull) {
                              slotErr[s] = "record '" + trimName(name) + "' in '" + slot.path +
                                           "' is longer than 4 Gbp";
                              return false;
                          }
                          // Must match pass 2 exactly or the buffers mis-size.
                          const uint32_t len = qualityTrimmedLength(qual, rawLen, qtrim_);
                          slot.trimmedBases += rawLen - len;
                          slot.lengths.push_back(len);
                          slot.bases += len;
                          if (len > slot.maxLen) slot.maxLen = len;
                          if (slot.needNames) slot.nameHashes.push_back(nameHashOf(name.data(), name.size()));
                          return true;
                      });
    });
    for (const std::string& e : slotErr) {
        if (!e.empty()) {
            error = e;
            return false;
        }
    }

    for (FileSlot& slot : slots) {
        slot.count = slot.lengths.size();
        trimmedBases_ += slot.trimmedBases;
    }

    // Pairing sanity checks, then the global read index layout.
    //
    // Paired libraries are laid out first, so every read below `pairedReads_`
    // has its mate at i^1 and everything above is unpaired. That lets one run
    // mix the two -- which is exactly what merging overlapping mates produces:
    // a pile of single-end fragments plus the pairs that did not merge.
    size_t nextRead = 0;
    for (int phase = 0; phase < 2; ++phase) {
    for (size_t li = 0; li < plans.size(); ++li) {
        LibPlan& plan = plans[li];
        const bool isSingle = plan.mode == kSingle;
        if ((phase == 0) == isSingle) continue;
        if (phase == 1 && nextRead > 0 && pairedReads_ == 0) pairedReads_ = nextRead;
        FileSlot& a = slots[plan.slotA];
        if (plan.mode == kPaired) {
            FileSlot& b = slots[plan.slotB];
            if (a.count != b.count) {
                error = "paired files are out of sync: '" + a.path + "' has " +
                        std::to_string(a.count) + " records but '" + b.path + "' has " +
                        std::to_string(b.count);
                return false;
            }
            for (size_t j = 0; j < a.count; ++j) {
                if (a.nameHashes[j] != b.nameHashes[j]) {
                    error = "paired files are out of sync at record " + std::to_string(j + 1) +
                            ": '" + a.path + "' has '" + recordNameAt(a.path, j) + "' but '" +
                            b.path + "' has '" + recordNameAt(b.path, j) + "'";
                    return false;
                }
            }
            a.base = nextRead;
            a.stride = 2;
            b.base = nextRead + 1;
            b.stride = 2;
            nextRead += a.count * 2;
        } else if (plan.mode == kInterleaved) {
            if (a.count % 2 != 0) {
                error = "interleaved file '" + a.path + "' has an odd number of records (" +
                        std::to_string(a.count) + "); its last read has no mate";
                return false;
            }
            for (size_t j = 0; j + 1 < a.count; j += 2) {
                if (a.nameHashes[j] != a.nameHashes[j + 1]) {
                    error = "interleaved file '" + a.path + "' is out of sync at record " +
                            std::to_string(j + 1) + ": '" + recordNameAt(a.path, j) +
                            "' is not mated with '" + recordNameAt(a.path, j + 1) + "'";
                    return false;
                }
            }
            a.base = nextRead;
            a.stride = 1;
            nextRead += a.count;
        } else {
            a.base = nextRead;
            a.stride = 1;
            nextRead += a.count;
        }
    }

    }
    if (pairedReads_ == 0) pairedReads_ = nextRead;   // no single-end library at all
    const size_t nReads = nextRead;
    uint64_t total = 0;
    for (const FileSlot& slot : slots) {
        total += slot.bases;
        if (slot.maxLen > maxLen_) maxLen_ = slot.maxLen;
    }
    totalBases_ = total;
    paired_ = pairedReads_ > 0;

    offsets_.assign(nReads + 1, 0);
    for (const FileSlot& slot : slots) {
        for (size_t j = 0; j < slot.count; ++j) {
            offsets_[slot.base + j * slot.stride + 1] = slot.lengths[j];
        }
    }
    for (size_t i = 0; i < nReads; ++i) offsets_[i + 1] += offsets_[i];

    for (FileSlot& slot : slots) {
        slot.lengths.clear();
        slot.lengths.shrink_to_fit();
        slot.nameHashes.clear();
        slot.nameHashes.shrink_to_fit();
    }

    data_.assign(static_cast<size_t>((totalBases_ + 31) / 32), 0);
    ambiguous_.assign(static_cast<size_t>((totalBases_ + 63) / 64), 0);

    // Pass 2: file-parallel packing into offsets that are already known.
    uint64_t* data = data_.data();
    uint64_t* amb = ambiguous_.data();
    runParallel(slots.size(), threads, [&](size_t s) {
        FileSlot& slot = slots[s];
        GzReader r;
        if (!r.open(slot.path, slotErr[s])) return;
        size_t j = 0;
        forEachRecord(r, slot.path, slotErr[s],
                      [&](const std::string& name, const char* seq, size_t rawLen, const char* qual) {
                          if (j >= slot.count) {
                              slotErr[s] = "'" + slot.path + "' changed while it was being read";
                              return false;
                          }
                          const size_t len = qualityTrimmedLength(qual, rawLen, qtrim_);
                          const size_t read = slot.base + j * slot.stride;
                          const uint64_t start = offsets_[read];
                          if (offsets_[read + 1] - start != len) {
                              slotErr[s] = "'" + slot.path + "' changed while it was being read (" +
                                           "record '" + trimName(name) + "')";
                              return false;
                          }
                          packSequence(data, amb, start, seq, static_cast<uint32_t>(len));
                          ++j;
                          return true;
                      });
        if (slotErr[s].empty() && j != slot.count) {
            slotErr[s] = "'" + slot.path + "' changed while it was being read";
        }
    });
    for (const std::string& e : slotErr) {
        if (!e.empty()) {
            error = e;
            data_.clear();
            ambiguous_.clear();
            offsets_.clear();
            totalBases_ = 0;
            maxLen_ = 0;
            paired_ = false;
            pairedReads_ = 0;
            return false;
        }
    }
    return true;
}

void SequenceStore::decode(size_t read, std::string& out) const {
    const uint32_t len = length(read);
    out.resize(len);
    const uint64_t start = offsets_[read];
    const bool hasAmb = !ambiguous_.empty();
    for (uint32_t p = 0; p < len; ++p) {
        const uint64_t bit = start + p;
        if (hasAmb && ((ambiguous_[bit >> 6] >> (bit & 63)) & 1ULL)) {
            out[p] = 'N';
            continue;
        }
        out[p] = codeBase(static_cast<int>((data_[bit >> 5] >> ((bit & 31) * 2)) & 3));
    }
}

void SequenceStore::setBase(size_t read, uint32_t pos, int code) {
    const uint64_t bit = offsets_[read] + pos;
    const unsigned shift = static_cast<unsigned>(bit & 31) * 2;
    uint64_t& word = data_[bit >> 5];
    word = (word & ~(static_cast<uint64_t>(3) << shift)) |
           (static_cast<uint64_t>(code & 3) << shift);
    if (!ambiguous_.empty()) ambiguous_[bit >> 6] &= ~(1ULL << (bit & 63));
}

void SequenceStore::maskRange(size_t read, uint32_t from, uint32_t to) {
    const uint32_t len = length(read);
    if (from >= to || from >= len) return;
    if (to > len) to = len;
    // load() always sizes the bitmap, so this only fires for an empty store,
    // where the range guard above has already returned. Kept as a cheap
    // guarantee that the indexing below is in bounds.
    if (ambiguous_.empty()) {
        ambiguous_.assign(static_cast<size_t>((totalBases_ + 63) / 64), 0);
    }
    const uint64_t base = offsets_[read];
    for (uint32_t p = from; p < to; ++p) {
        const uint64_t bit = base + p;
        ambiguous_[bit >> 6] |= 1ULL << (bit & 63);
    }
}

bool writeFasta(const std::string& path, const std::vector<std::string>& seqs,
                const std::vector<std::string>& names, int lineWidth, std::string& error) {
    error.clear();
    if (!names.empty() && names.size() != seqs.size()) {
        error = "writeFasta: " + std::to_string(seqs.size()) + " sequences but " +
                std::to_string(names.size()) + " names";
        return false;
    }
    std::FILE* f = std::fopen(path.c_str(), "wb");
    if (!f) {
        error = "cannot open '" + path + "' for writing: " + std::strerror(errno);
        return false;
    }

    std::string buf;
    buf.reserve(kWriteBuffer + (1u << 16));
    auto flush = [&]() {
        if (buf.empty()) return true;
        bool ok = std::fwrite(buf.data(), 1, buf.size(), f) == buf.size();
        buf.clear();
        return ok;
    };

    const size_t width = lineWidth > 0 ? static_cast<size_t>(lineWidth) : 0;
    for (size_t i = 0; i < seqs.size(); ++i) {
        buf += '>';
        if (names.empty()) {
            buf += "contig_";
            buf += std::to_string(i + 1);
        } else {
            buf += names[i];
        }
        buf += '\n';
        const std::string& s = seqs[i];
        if (width == 0 || s.size() <= width) {
            buf += s;
            buf += '\n';
        } else {
            for (size_t p = 0; p < s.size(); p += width) {
                buf.append(s, p, std::min(width, s.size() - p));
                buf += '\n';
            }
        }
        if (buf.size() >= kWriteBuffer && !flush()) {
            error = "write error on '" + path + "': " + std::strerror(errno);
            std::fclose(f);
            return false;
        }
    }
    if (!flush()) {
        error = "write error on '" + path + "': " + std::strerror(errno);
        std::fclose(f);
        return false;
    }
    if (std::fclose(f) != 0) {
        error = "write error on '" + path + "': " + std::strerror(errno);
        return false;
    }
    return true;
}

}  // namespace ts
