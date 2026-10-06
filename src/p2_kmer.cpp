#include "p2_kmer.h"

#include <algorithm>

namespace ts {
namespace p2 {

uint64_t canon31At(const std::string& s, size_t pos) {
    if (pos + kK > s.size()) return UINT64_MAX;
    uint64_t fw = 0, rv = 0;
    constexpr int shift = 2 * (kK - 1);
    for (int i = 0; i < kK; ++i) {
        const int c = baseCode2(s[pos + static_cast<size_t>(i)]);
        if (c < 0) return UINT64_MAX;
        fw = ((fw << 2) | static_cast<uint64_t>(c)) & kKMask;
        rv = (rv >> 2) | (static_cast<uint64_t>(3 - c) << shift);
    }
    return fw < rv ? fw : rv;
}

std::string reverseComplement(const std::string& s) {
    std::string r(s.rbegin(), s.rend());
    for (char& c : r) {
        switch (c) {
            case 'A': c = 'T'; break;
            case 'C': c = 'G'; break;
            case 'G': c = 'C'; break;
            case 'T': c = 'A'; break;
            case 'a': c = 't'; break;
            case 'c': c = 'g'; break;
            case 'g': c = 'c'; break;
            case 't': c = 'a'; break;
            default: break;
        }
    }
    return r;
}

// ---- FixedCountTable ---------------------------------------------------------------------------

void FixedCountTable::addKey(uint64_t canon) { pending_.push_back(canon); }

size_t FixedCountTable::slotOfH(uint64_t canon, uint64_t h) const {
    size_t i = static_cast<size_t>(h) & mask_;
    while (keys_[i] && keys_[i] != canon + 1) i = (i + 1) & mask_;
    return i;
}

void FixedCountTable::freeze(bool filter) {
    std::sort(pending_.begin(), pending_.end());
    pending_.erase(std::unique(pending_.begin(), pending_.end()), pending_.end());
    size_t cap = 16;
    while (cap < pending_.size() * 2 + 16) cap <<= 1;
    keys_.assign(cap, 0);
    counts_.reset(new std::atomic<uint32_t>[cap]);
    for (size_t i = 0; i < cap; ++i) counts_[i].store(0, std::memory_order_relaxed);
    mask_ = cap - 1;
    n_ = 0;
    if (filter) filter_.assign((size_t(1) << 24) / 64, 0);
    else filter_.clear();
    for (uint64_t k : pending_) {
        const uint64_t h = mix64(k);
        const size_t i = slotOfH(k, h);
        if (!keys_[i]) { keys_[i] = k + 1; ++n_; }
        if (filter) {
            const uint64_t top = h >> 40;
            filter_[top >> 6] |= 1ULL << (top & 63);
        }
    }
    pending_.clear();
    pending_.shrink_to_fit();
    frozen_ = true;
}

void FixedCountTable::resetCounts() {
    for (size_t i = 0; i < keys_.size(); ++i) counts_[i].store(0, std::memory_order_relaxed);
}

bool FixedCountTable::contains(uint64_t canon) const {
    if (!frozen_ || keys_.empty()) return false;
    return keys_[slotOf(canon)] != 0;
}

void FixedCountTable::hitHashed(uint64_t canon, uint64_t h) {
    if (!frozen_ || n_ == 0) return;
    if (!filter_.empty()) {
        const uint64_t top = h >> 40;
        if (!((filter_[top >> 6] >> (top & 63)) & 1ULL)) return;
    }
    const size_t i = slotOfH(canon, h);
    if (!keys_[i]) return;
    uint32_t cur = counts_[i].load(std::memory_order_relaxed);
    while (cur != UINT32_MAX &&
           !counts_[i].compare_exchange_weak(cur, cur + 1, std::memory_order_relaxed)) {
    }
}

void FixedCountTable::hitSerial(uint64_t canon) {
    if (!frozen_ || n_ == 0) return;
    const size_t i = slotOf(canon);
    if (!keys_[i]) return;
    const uint32_t cur = counts_[i].load(std::memory_order_relaxed);
    if (cur != UINT32_MAX) counts_[i].store(cur + 1, std::memory_order_relaxed);
}

uint32_t FixedCountTable::count(uint64_t canon) const {
    if (!frozen_ || n_ == 0) return 0;
    const size_t i = slotOf(canon);
    return keys_[i] ? counts_[i].load(std::memory_order_relaxed) : 0;
}

// ---- ShardedCountTable -------------------------------------------------------------------------

ShardedCountTable::ShardedCountTable() : shards_(kShards) {
    for (Shard& s : shards_) {
        s.slots.assign(size_t(1) << 14, Slot{0, 0, 0, 0});
        s.mask = s.slots.size() - 1;
    }
}

void ShardedCountTable::grow(Shard& s) {
    std::vector<Slot> old;
    old.swap(s.slots);
    s.slots.assign(old.size() * 2, Slot{0, 0, 0, 0});
    s.mask = s.slots.size() - 1;
    for (const Slot& x : old) {
        if (!x.key) continue;
        size_t i = slotHash(x.key - 1) & s.mask;
        while (s.slots[i].key) i = (i + 1) & s.mask;
        s.slots[i] = x;
    }
}

ShardedCountTable::Slot& ShardedCountTable::insert(Shard& s, uint64_t canon) {
    if ((s.n + 1) * 10 > s.slots.size() * 7) grow(s);
    size_t i = slotHash(canon) & s.mask;
    while (s.slots[i].key && s.slots[i].key != canon + 1) i = (i + 1) & s.mask;
    if (!s.slots[i].key) { s.slots[i].key = canon + 1; ++s.n; }
    return s.slots[i];
}

void ShardedCountTable::addAssembly(uint64_t canon) {
    Slot& x = insert(shards_[shardOf(canon)], canon);
    if (x.ac < 65535) ++x.ac;
}

void ShardedCountTable::addReadBatch(size_t shard, const uint64_t* keys, size_t n) {
    Shard& s = shards_[shard];
    std::lock_guard<std::mutex> lock(s.m);
    for (size_t j = 0; j < n; ++j) {
        Slot& x = insert(s, keys[j]);
        if (x.rc != UINT32_MAX) ++x.rc;
    }
}

size_t ShardedCountTable::distinct() const {
    size_t n = 0;
    for (const Shard& s : shards_) n += s.n;
    return n;
}

size_t ShardedCountTable::memoryBytes() const {
    size_t n = 0;
    for (const Shard& s : shards_) n += s.slots.capacity() * sizeof(Slot);
    return n;
}

uint32_t ShardedCountTable::readCount(uint64_t canon) const {
    const Shard& s = shards_[shardOf(canon)];
    size_t i = slotHash(canon) & s.mask;
    while (s.slots[i].key && s.slots[i].key != canon + 1) i = (i + 1) & s.mask;
    return s.slots[i].key ? s.slots[i].rc : 0;
}

}  // namespace p2
}  // namespace ts
