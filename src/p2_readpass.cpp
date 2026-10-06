#include "p2_readpass.h"

#include <zlib.h>

#include <algorithm>
#include <chrono>
#include <condition_variable>
#include <cstring>
#include <deque>
#include <mutex>
#include <thread>

namespace ts {
namespace p2 {

namespace {

// A batch: record sequences back to back, `ends` holding each record's end offset.
struct Batch {
    std::string data;
    std::vector<uint32_t> ends;
};

class BatchQueue {
public:
    explicit BatchQueue(size_t cap) : cap_(cap) {}
    void push(std::unique_ptr<Batch> b) {
        std::unique_lock<std::mutex> l(m_);
        notFull_.wait(l, [&] { return q_.size() < cap_; });
        q_.push_back(std::move(b));
        notEmpty_.notify_one();
    }
    // nullptr once every producer is done and the queue is empty.
    std::unique_ptr<Batch> pop() {
        std::unique_lock<std::mutex> l(m_);
        notEmpty_.wait(l, [&] { return !q_.empty() || producers_ == 0; });
        if (q_.empty()) return nullptr;
        std::unique_ptr<Batch> b = std::move(q_.front());
        q_.pop_front();
        notFull_.notify_one();
        return b;
    }
    void addProducer() { std::lock_guard<std::mutex> l(m_); ++producers_; }
    void producerDone() {
        std::lock_guard<std::mutex> l(m_);
        --producers_;
        notEmpty_.notify_all();
    }

private:
    std::mutex m_;
    std::condition_variable notFull_, notEmpty_;
    std::deque<std::unique_ptr<Batch>> q_;
    size_t cap_;
    int producers_ = 0;
};

constexpr size_t kBatchBytes = 4u << 20;

struct FileResult {
    uint64_t reads = 0, bases = 0, wanted = 0;
    std::string error;
};

inline bool wantedRead(const ReadPassSinks& k, uint64_t idx) {
    const std::vector<uint64_t>& b = *k.wantBits;
    const uint64_t w = idx >> 6;
    return w < b.size() && ((b[w] >> (idx & 63)) & 1ULL);
}

// Decompresses one file and pushes its record sequences in batches; hands the wanted records to
// sinks.onWanted (after the 3' quality trim for FASTQ).
void readFile(const std::string& path, size_t fileIdx, BatchQueue& q, FileResult& out, ReadPassSinks& sinks) {
    gzFile g = gzopen(path.c_str(), "rb");
    if (!g) { out.error = "cannot open " + path; return; }
    gzbuffer(g, 1 << 20);
    const bool collect = sinks.wantBits && sinks.onWanted && fileIdx < sinks.layout.size() &&
                         sinks.layout[fileIdx].second > 0;
    const uint64_t base = collect ? sinks.layout[fileIdx].first : 0;
    const uint64_t stride = collect ? sinks.layout[fileIdx].second : 0;
    std::vector<char> buf(1u << 22);
    std::string carry;           // an incomplete line from the previous block
    std::string fastaSeq;        // the FASTA record being assembled
    std::string lastSeq;         // the FASTQ sequence line awaiting its quality line
    bool lastWanted = false;
    int format = 0;              // 1 FASTQ, 2 FASTA
    uint64_t lineNo = 0;
    std::unique_ptr<Batch> b(new Batch);
    auto emit = [&](const char* s, size_t n) {
        ++out.reads;
        out.bases += n;
        b->data.append(s, n);
        b->ends.push_back(static_cast<uint32_t>(b->data.size()));
        if (b->data.size() >= kBatchBytes) {
            q.push(std::move(b));
            b.reset(new Batch);
        }
    };
    auto line = [&](const char* s, size_t n) {
        while (n && (s[n - 1] == '\r')) --n;
        if (format == 0) {
            if (n == 0) return;
            format = s[0] == '>' ? 2 : 1;
        }
        if (format == 1) {
            const uint64_t phase = lineNo & 3;
            if (phase == 1) {
                const uint64_t rec = out.reads;       // record index in this file
                emit(s, n);
                lastWanted = collect && wantedRead(sinks, base + rec * stride);
                if (lastWanted) lastSeq.assign(s, n);
            } else if (phase == 3 && lastWanted) {
                size_t e = std::min(n, lastSeq.size());
                while (e > 0 && static_cast<int>(static_cast<unsigned char>(s[e - 1])) - 33 < sinks.trimQ) --e;
                sinks.onWanted(base + (out.reads - 1) * stride, lastSeq.data(), e);
                ++out.wanted;
                lastWanted = false;
            }
            ++lineNo;
        } else {
            if (n && s[0] == '>') {
                if (!fastaSeq.empty()) {
                    const uint64_t rec = out.reads;
                    emit(fastaSeq.data(), fastaSeq.size());
                    if (collect && wantedRead(sinks, base + rec * stride)) {
                        sinks.onWanted(base + rec * stride, fastaSeq.data(), fastaSeq.size());
                        ++out.wanted;
                    }
                }
                fastaSeq.clear();
            } else {
                fastaSeq.append(s, n);
            }
        }
    };
    for (;;) {
        const int got = gzread(g, buf.data(), static_cast<unsigned>(buf.size()));
        if (got < 0) {
            int zerr = 0;
            const char* msg = gzerror(g, &zerr);
            out.error = path + ": " + (msg ? msg : "read error");
            break;
        }
        if (got == 0) break;
        const char* p = buf.data();
        const char* end = p + got;
        while (p < end) {
            const char* nl = static_cast<const char*>(std::memchr(p, '\n', static_cast<size_t>(end - p)));
            if (!nl) { carry.append(p, static_cast<size_t>(end - p)); break; }
            if (!carry.empty()) {
                carry.append(p, static_cast<size_t>(nl - p));
                line(carry.data(), carry.size());
                carry.clear();
            } else {
                line(p, static_cast<size_t>(nl - p));
            }
            p = nl + 1;
        }
    }
    if (out.error.empty()) {
        if (!carry.empty()) line(carry.data(), carry.size());
        if (format == 2 && !fastaSeq.empty()) {
            const uint64_t rec = out.reads;
            emit(fastaSeq.data(), fastaSeq.size());
            if (collect && wantedRead(sinks, base + rec * stride)) {
                sinks.onWanted(base + rec * stride, fastaSeq.data(), fastaSeq.size());
                ++out.wanted;
            }
        }
        if (!b->ends.empty()) q.push(std::move(b));
    }
    gzclose(g);
}

// Per-worker state: shard buffers for the sampled table.
struct Worker {
    std::vector<std::vector<uint64_t>> buf;
    explicit Worker(bool sampled) {
        if (sampled) buf.assign(ShardedCountTable::kShards, std::vector<uint64_t>());
    }
};

constexpr size_t kFlush = 1024;

inline void countOne(const char* s, size_t n, ReadPassSinks& k, Worker& w) {
    forEachCanon31(s, n, [&](size_t, uint64_t c) {
        if (k.sampled && sqaSampled(c, k.sampleS)) {
            const size_t sh = ShardedCountTable::shardOf(c);
            std::vector<uint64_t>& v = w.buf[sh];
            v.push_back(c);
            if (v.size() >= kFlush) { k.sampled->addReadBatch(sh, v.data(), v.size()); v.clear(); }
        }
        if (k.exact || k.detectSlot) {
            const uint64_t h = mix64(c);
            if (k.exact) k.exact->hitHashed(c, h);
            if (k.detectSlot) {
                uint32_t times = 0;
                for (uint64_t thr : k.detectThresholds) if (h <= thr) ++times;
                if (times) {
                    auto it = k.detectSlot->find(c);
                    if (it != k.detectSlot->end())
                        k.detectCounts[it->second].fetch_add(times, std::memory_order_relaxed);
                }
            }
        }
    });
}

void flushWorker(ReadPassSinks& k, Worker& w) {
    if (!k.sampled) return;
    for (size_t sh = 0; sh < w.buf.size(); ++sh) {
        if (!w.buf[sh].empty()) k.sampled->addReadBatch(sh, w.buf[sh].data(), w.buf[sh].size());
        w.buf[sh].clear();
    }
}

}  // namespace

bool runReadPass(const std::vector<std::string>& files, int threads, ReadPassSinks& sinks, ReadPassStats& st) {
    const auto t0 = std::chrono::steady_clock::now();
    st = ReadPassStats();
    st.files = files.size();
    if (threads < 1) threads = 1;
    BatchQueue q(static_cast<size_t>(threads) * 2 + 2);
    std::vector<FileResult> results(files.size());
    std::vector<std::thread> readers, workers;
    for (size_t i = 0; i < files.size(); ++i) q.addProducer();
    for (size_t i = 0; i < files.size(); ++i) {
        readers.emplace_back([&, i]() {
            readFile(files[i], i, q, results[i], sinks);
            q.producerDone();
        });
    }
    for (int t = 0; t < threads; ++t) {
        workers.emplace_back([&]() {
            Worker w(sinks.sampled != nullptr);
            while (std::unique_ptr<Batch> b = q.pop()) {
                uint32_t from = 0;
                for (uint32_t e : b->ends) {
                    countOne(b->data.data() + from, e - from, sinks, w);
                    from = e;
                }
            }
            flushWorker(sinks, w);
        });
    }
    for (std::thread& t : readers) t.join();
    for (std::thread& t : workers) t.join();
    for (const FileResult& r : results) {
        st.reads += r.reads;
        st.bases += r.bases;
        st.wanted += r.wanted;
        st.fileRecords.push_back(r.reads);
        if (st.error.empty() && !r.error.empty()) st.error = r.error;
    }
    st.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return st.error.empty();
}

void countSequences(const std::vector<std::string>& seqs, ReadPassSinks& sinks, ReadPassStats& st) {
    Worker w(sinks.sampled != nullptr);
    for (const std::string& s : seqs) {
        ++st.reads;
        st.bases += s.size();
        countOne(s.data(), s.size(), sinks, w);
    }
    flushWorker(sinks, w);
}

}  // namespace p2
}  // namespace ts
