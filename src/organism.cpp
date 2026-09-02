#include "organism.h"

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <memory>

namespace ts {

namespace {

// Version 5 records the marker sampling denominator the model was built at; version 4
// appended per-plasmid marker membership; version 3 appended chromosome layout tracks;
// version 2 has neither. All still load, and a model missing a section simply does not offer
// the stage that consumes it -- adjacency joining works without tracks, and layout works
// without membership. Rebuilding a fold model costs minutes and there are ten of them on
// disk, so making older files unreadable would be a real cost for no benefit.
//
// The density field is the one addition that cannot be defaulted wrongly and left alone: a
// query run at a different density than the model shares almost none of its markers and
// reports nothing, without erroring. Files older than version 5 predate the option and were
// all built at 512, which is what they are read as.
constexpr char kMagic[8] = {'T', 'S', 'M', 'O', 'D', 'E', 'L', '5'};
constexpr char kMagicV4[8] = {'T', 'S', 'M', 'O', 'D', 'E', 'L', '4'};
constexpr char kMagicV3[8] = {'T', 'S', 'M', 'O', 'D', 'E', 'L', '3'};
constexpr char kMagicV2[8] = {'T', 'S', 'M', 'O', 'D', 'E', 'L', '2'};

template <typename T>
bool writePod(std::FILE* f, const T& v) {
    return std::fwrite(&v, sizeof(T), 1, f) == 1;
}
template <typename T>
bool readPod(std::FILE* f, T& v) {
    return std::fread(&v, sizeof(T), 1, f) == 1;
}

}  // namespace

void OrganismModel::beginBuild(const std::string& organism, int k) {
    organism_ = organism;
    k_ = k;
    genomesChr_ = genomesPls_ = 0;
    index_.clear();
    for (int i = 0; i < 2; ++i) {
        genomesPerMarker_[i].clear();
        edges_[i].clear();
        pending_[i].clear();
    }
}

uint32_t OrganismModel::internMarker(uint64_t canonicalKmer) {
    auto it = index_.find(canonicalKmer);
    if (it != index_.end()) return it->second;
    const uint32_t id = static_cast<uint32_t>(genomesPerMarker_[0].size());
    index_.emplace(canonicalKmer, id);
    genomesPerMarker_[0].push_back(0);
    genomesPerMarker_[1].push_back(0);
    return id;
}

void OrganismModel::addObservation(uint64_t fromOriented, uint64_t toOriented, int32_t dist,
                                   Replicon r) {
    Pending& p = pending_[static_cast<int>(r)][(fromOriented << 32) | toOriented];
    if (p.count < static_cast<uint32_t>(kDistSamples)) p.dist[p.count] = dist;
    ++p.count;
}

void OrganismModel::finalise(uint32_t minSupportChr, uint32_t minSupportPls) {
    for (int c = 0; c < 2; ++c) {
        const uint32_t minSupport = c == 0 ? minSupportChr : minSupportPls;
        edges_[c].clear();
        edges_[c].reserve(pending_[c].size());
        for (auto& kv : pending_[c]) {
            Pending& p = kv.second;
            if (p.count < minSupport) continue;
            const int n = static_cast<int>(std::min<uint32_t>(p.count, kDistSamples));
            std::sort(p.dist, p.dist + n);
            MarkerEdge e;
            e.support = p.count;
            e.medianDist = p.dist[n / 2];
            edges_[c].emplace(kv.first, e);
        }
        pending_[c].clear();
    }
    loaded_ = true;
}

bool OrganismModel::save(const std::string& path, std::string& error) const {
    std::FILE* f = std::fopen(path.c_str(), "wb");
    if (!f) { error = "cannot write model file: " + path; return false; }
    bool ok = std::fwrite(kMagic, 1, 8, f) == 8;

    const uint32_t kk = static_cast<uint32_t>(k_);
    ok = ok && writePod(f, kk) && writePod(f, markerDenom_) &&
         writePod(f, genomesChr_) && writePod(f, genomesPls_);

    const uint32_t nameLen = static_cast<uint32_t>(organism_.size());
    ok = ok && writePod(f, nameLen) &&
         (nameLen == 0 || std::fwrite(organism_.data(), 1, nameLen, f) == nameLen);

    const uint32_t nExcluded = static_cast<uint32_t>(excluded_.size());
    ok = ok && writePod(f, nExcluded);
    for (const std::string& acc : excluded_) {
        const uint32_t n = static_cast<uint32_t>(acc.size());
        ok = ok && writePod(f, n) && (n == 0 || std::fwrite(acc.data(), 1, n, f) == n);
    }

    const uint64_t nMarkers = genomesPerMarker_[0].size();
    ok = ok && writePod(f, nMarkers);
    // Markers are written in id order so ids survive the round trip.
    std::vector<uint64_t> byId(nMarkers, 0);
    for (const auto& kv : index_) byId[kv.second] = kv.first;
    for (uint64_t i = 0; i < nMarkers && ok; ++i) {
        ok = writePod(f, byId[i]) && writePod(f, genomesPerMarker_[0][i]) &&
             writePod(f, genomesPerMarker_[1][i]);
    }

    for (int c = 0; c < 2 && ok; ++c) {
        const uint64_t nEdges = edges_[c].size();
        ok = writePod(f, nEdges);
        for (const auto& kv : edges_[c]) {
            if (!ok) break;
            ok = writePod(f, kv.first) && writePod(f, kv.second.support) &&
                 writePod(f, kv.second.medianDist);
        }
    }

    const uint64_t nTracks = tracks_.size();
    ok = ok && writePod(f, nTracks);
    for (const LayoutTrack& t : tracks_) {
        if (!ok) break;
        const uint32_t nameLen2 = static_cast<uint32_t>(t.name.size());
        ok = writePod(f, nameLen2) &&
             (nameLen2 == 0 || std::fwrite(t.name.data(), 1, nameLen2, f) == nameLen2);
        const uint64_t n = t.oriented.size();
        ok = ok && writePod(f, n);
        for (uint64_t i = 0; i < n && ok; ++i) {
            ok = writePod(f, t.oriented[i]) && writePod(f, t.pos[i]);
        }
    }

    const uint64_t nPls = plasmids_.size();
    ok = ok && writePod(f, nPls);
    for (const PlasmidMembership& m : plasmids_) {
        if (!ok) break;
        const uint32_t nl = static_cast<uint32_t>(m.name.size());
        ok = writePod(f, nl) && (nl == 0 || std::fwrite(m.name.data(), 1, nl, f) == nl) &&
             writePod(f, m.length);
        const uint64_t nm = m.markers.size();
        ok = ok && writePod(f, nm);
        for (uint64_t i = 0; i < nm && ok; ++i) ok = writePod(f, m.markers[i]);
    }

    std::fclose(f);
    if (!ok) error = "short write on model file: " + path;
    return ok;
}

bool OrganismModel::checkHeader(const std::string& path, std::string& error) {
    std::FILE* raw = std::fopen(path.c_str(), "rb");
    if (!raw) { error = "cannot open model file: " + path; return false; }
    const std::unique_ptr<std::FILE, int (*)(std::FILE*)> guard(raw, std::fclose);

    char magic[8];
    if (std::fread(magic, 1, 8, raw) != 8) {
        error = "not a TesserACT model file (too short): " + path;
        return false;
    }
    const bool known = std::memcmp(magic, kMagic, 8) == 0 ||
                       std::memcmp(magic, kMagicV4, 8) == 0 ||
                       std::memcmp(magic, kMagicV3, 8) == 0 ||
                       std::memcmp(magic, kMagicV2, 8) == 0;
    if (!known) {
        error = "not a TesserACT model file (or a model from an older version): " + path;
        return false;
    }
    // The declared marker count against the bytes that could possibly hold it. A download
    // cut short keeps a valid magic and a plausible header, so the magic alone would pass
    // it; the size is what exposes the truncation without reading the whole file.
    uint32_t kk = 0, denom = 0, gchr = 0, gpls = 0, nameLen = 0, nExcluded = 0;
    bool ok = readPod(raw, kk);
    if (ok && std::memcmp(magic, kMagic, 8) == 0) ok = readPod(raw, denom);
    ok = ok && readPod(raw, gchr) && readPod(raw, gpls) && readPod(raw, nameLen);
    if (ok && nameLen > 0) ok = std::fseek(raw, nameLen, SEEK_CUR) == 0;
    ok = ok && readPod(raw, nExcluded);
    for (uint32_t i = 0; i < nExcluded && ok; ++i) {
        uint32_t n = 0;
        ok = readPod(raw, n) && (n == 0 || std::fseek(raw, n, SEEK_CUR) == 0);
    }
    uint64_t nMarkers = 0;
    ok = ok && readPod(raw, nMarkers);
    if (!ok) { error = "truncated or corrupt model file: " + path; return false; }

    long here = std::ftell(raw);
    if (std::fseek(raw, 0, SEEK_END) == 0 && here >= 0) {
        const long end = std::ftell(raw);
        // 16 bytes per marker record, and the tables that follow are larger still, so this
        // is a floor rather than an estimate: anything smaller cannot be a whole model.
        if (end >= 0 && static_cast<uint64_t>(end - here) < nMarkers * 16ULL) {
            error = "truncated or corrupt model file: " + path;
            return false;
        }
    }
    return true;
}

bool OrganismModel::load(const std::string& path, std::string& error) {
    std::FILE* raw = std::fopen(path.c_str(), "rb");
    if (!raw) { error = "cannot open model file: " + path; return false; }
    // Held by the guard so the handle is closed even if an allocation below
    // throws: every count in this file comes off disk, and a corrupt one used
    // to reach resize() directly.
    const std::unique_ptr<std::FILE, int (*)(std::FILE*)> guard(raw, std::fclose);
    std::FILE* f = raw;

    char magic[8];
    if (std::fread(magic, 1, 8, f) != 8) {
        error = "not a TesserACT model file: " + path;
        return false;
    }
    const bool hasDensity = std::memcmp(magic, kMagic, 8) == 0;
    const bool hasPlasmids = hasDensity || std::memcmp(magic, kMagicV4, 8) == 0;
    const bool hasTracks = hasPlasmids || std::memcmp(magic, kMagicV3, 8) == 0;
    if (!hasTracks && std::memcmp(magic, kMagicV2, 8) != 0) {
        error = "not a TesserACT model file (or a model from an older version): " + path;
        return false;
    }

    // The file's own length, used to sanity-check the counts it declares. A
    // count larger than the bytes that could hold it is corruption, and must be
    // rejected before it reaches an allocation.
    long fileSize = 0;
    if (std::fseek(f, 0, SEEK_END) == 0) {
        fileSize = std::ftell(f);
        if (std::fseek(f, 8, SEEK_SET) != 0) fileSize = 0;
    }
    const auto plausible = [&](uint64_t count, uint64_t bytesEach) {
        if (fileSize <= 0) return true;              // unseekable: fall back to short-read checks
        const uint64_t left = static_cast<uint64_t>(fileSize);
        return bytesEach == 0 || count <= left / bytesEach;
    };

    uint32_t kk = 0;
    bool ok = readPod(f, kk);
    k_ = static_cast<int>(kk);
    markerDenom_ = kMarkerSampleDenom;
    if (ok && hasDensity) {
        uint32_t denom = 0;
        ok = readPod(f, denom);
        if (ok && denom == 0) {
            error = "model declares a marker sampling denominator of zero: " + path;
            return false;
        }
        if (ok) markerDenom_ = denom;
    }
    ok = ok && readPod(f, genomesChr_) && readPod(f, genomesPls_);

    uint32_t nameLen = 0;
    ok = ok && readPod(f, nameLen) && plausible(nameLen, 1);
    if (ok) {
        organism_.assign(nameLen, '\0');
        ok = nameLen == 0 || std::fread(&organism_[0], 1, nameLen, f) == nameLen;
    }

    uint32_t nExcluded = 0;
    ok = ok && readPod(f, nExcluded) && plausible(nExcluded, 4);
    excluded_.clear();
    for (uint32_t i = 0; i < nExcluded && ok; ++i) {
        uint32_t n = 0;
        ok = readPod(f, n) && plausible(n, 1);
        if (!ok) break;
        std::string acc(n, '\0');
        ok = ok && (n == 0 || std::fread(&acc[0], 1, n, f) == n);
        if (ok) excluded_.push_back(acc);
    }

    uint64_t nMarkers = 0;
    // 8 bytes of k-mer plus two 4-byte counts per marker.
    ok = ok && readPod(f, nMarkers) && plausible(nMarkers, 16);
    index_.clear();
    if (ok) {
        index_.reserve(nMarkers * 2);
        genomesPerMarker_[0].resize(nMarkers);
        genomesPerMarker_[1].resize(nMarkers);
    }
    for (uint64_t i = 0; i < nMarkers && ok; ++i) {
        uint64_t km = 0;
        uint32_t gc = 0, gp = 0;
        ok = readPod(f, km) && readPod(f, gc) && readPod(f, gp);
        if (!ok) break;
        index_.emplace(km, static_cast<uint32_t>(i));
        genomesPerMarker_[0][i] = gc;
        genomesPerMarker_[1][i] = gp;
    }

    for (int c = 0; c < 2 && ok; ++c) {
        uint64_t nEdges = 0;
        // 8 bytes of key plus support and median distance.
        ok = readPod(f, nEdges) && plausible(nEdges, 16);
        edges_[c].clear();
        if (ok) edges_[c].reserve(nEdges * 2);
        for (uint64_t i = 0; i < nEdges && ok; ++i) {
            uint64_t key = 0;
            MarkerEdge e;
            ok = readPod(f, key) && readPod(f, e.support) && readPod(f, e.medianDist);
            if (ok) edges_[c].emplace(key, e);
        }
    }

    tracks_.clear();
    if (ok && hasTracks) {
        uint64_t nTracks = 0;
        // Each track costs at least a length field per entry, so 8 bytes is the floor.
        ok = readPod(f, nTracks) && plausible(nTracks, 8);
        for (uint64_t t = 0; t < nTracks && ok; ++t) {
            LayoutTrack track;
            uint32_t nameLen2 = 0;
            ok = readPod(f, nameLen2) && plausible(nameLen2, 1);
            if (!ok) break;
            track.name.assign(nameLen2, '\0');
            ok = nameLen2 == 0 || std::fread(&track.name[0], 1, nameLen2, f) == nameLen2;
            uint64_t n = 0;
            ok = ok && readPod(f, n) && plausible(n, 8);
            if (!ok) break;
            track.oriented.resize(n);
            track.pos.resize(n);
            for (uint64_t i = 0; i < n && ok; ++i) {
                ok = readPod(f, track.oriented[i]) && readPod(f, track.pos[i]);
            }
            if (ok) tracks_.push_back(std::move(track));
        }
    }

    plasmids_.clear();
    if (ok && hasPlasmids) {
        uint64_t nPls = 0;
        ok = readPod(f, nPls) && plausible(nPls, 12);
        for (uint64_t t = 0; t < nPls && ok; ++t) {
            PlasmidMembership m;
            uint32_t nl = 0;
            ok = readPod(f, nl) && plausible(nl, 1);
            if (!ok) break;
            m.name.assign(nl, '\0');
            ok = (nl == 0 || std::fread(&m.name[0], 1, nl, f) == nl) && readPod(f, m.length);
            uint64_t nm = 0;
            ok = ok && readPod(f, nm) && plausible(nm, 4);
            if (!ok) break;
            m.markers.resize(nm);
            for (uint64_t i = 0; i < nm && ok; ++i) ok = readPod(f, m.markers[i]);
            if (ok) plasmids_.push_back(std::move(m));
        }
    }

    if (!ok) {
        error = "truncated or corrupt model file: " + path;
        return false;
    }
    loaded_ = true;
    return true;
}

}  // namespace ts
