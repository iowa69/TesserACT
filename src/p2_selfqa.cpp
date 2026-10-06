#include "p2_selfqa.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>

#include "p2_json.h"

namespace ts {
namespace p2 {

const char* const kSqaBinName[kSqaBins] = {"lt0.1x", "0.1-0.5x", "0.5-1.5x", "1.5-3x", "3-10x", "gt10x"};

namespace {
constexpr uint32_t kHistMax = 100000;
const double kEdges[kSqaBins + 1] = {0, 0.1, 0.5, 1.5, 3, 10, 1e18};
}  // namespace

std::string kmerString(uint64_t canon) {
    std::string s(kK, 'A');
    for (int i = kK - 1; i >= 0; --i) {
        s[static_cast<size_t>(i)] = "ACGT"[canon & 3];
        canon >>= 2;
    }
    return s;
}

double sqaQv(uint64_t bad, uint64_t positions) {
    const double p = positions ? static_cast<double>(bad) / static_cast<double>(positions) : 0;
    if (p <= 0) return 99.0;
    const double e = 1 - std::pow(1 - p, 1.0 / kK);
    return std::min(99.0, -10 * std::log10(e));
}

void SelfQa::prepare(const std::vector<std::string>& records, const std::vector<std::string>& names,
                     const std::vector<std::string>& spike) {
    recs_.clear();
    names_.clear();
    asmPos_ = asmBases_ = 0;
    for (size_t i = 0; i < records.size(); ++i) {
        if (records[i].size() < kSqaMinRecord) continue;
        recs_.push_back(&records[i]);
        names_.push_back(i < names.size() ? names[i] : std::string("record_") + std::to_string(i + 1));
        asmBases_ += records[i].size();
        forEachCanon31(records[i].data(), records[i].size(), [&](size_t, uint64_t c) {
            if (!sqaSampled(c, kSqaSample)) return;
            table_.addAssembly(c);
            ++asmPos_;
        });
    }
    spike_.clear();
    for (const std::string& s : spike)
        forEachCanon31(s.data(), s.size(), [&](size_t, uint64_t c) {
            if (sqaSampled(c, kSqaSample)) spike_.push_back(c);
        });
    std::sort(spike_.begin(), spike_.end());
    spike_.erase(std::unique(spike_.begin(), spike_.end()), spike_.end());
}

SelfQaResult SelfQa::finish(uint64_t reads, uint64_t readBases) const {
    const auto t0 = std::chrono::steady_clock::now();
    SelfQaResult r;
    r.ran = true;
    r.asmRecords = recs_.size();
    r.asmBases = asmBases_;
    r.reads = reads;
    r.readBases = readBases;
    r.asmPos = asmPos_;
    r.tableDistinct = table_.distinct();
    r.tableBytes = table_.memoryBytes();

    // histogram of read counts (sampled distinct k-mers seen in the reads), first local minimum from 2,
    // and the peak above it (kqa.cpp, unchanged)
    std::vector<uint64_t> h(kHistMax + 1, 0);
    table_.forEach([&](uint64_t, uint32_t rc, uint16_t) {
        if (rc) ++h[std::min(rc, kHistMax)];
    });
    uint32_t valley = 2;
    while (valley + 1 < 1000 && h[valley] >= h[valley + 1]) ++valley;
    uint32_t peak = valley;
    for (uint32_t c = valley; c < 2000; ++c) if (h[c] > h[peak]) peak = c;
    r.valley = valley;
    r.peak = peak;

    auto isSpike = [&](uint64_t k) { return std::binary_search(spike_.begin(), spike_.end(), k); };
    std::vector<uint32_t> spikeCounts;
    table_.forEach([&](uint64_t key, uint32_t c, uint16_t a) {
        if (a) {
            if (c == 0) r.asmRead0 += a;
            if (c < valley) r.asmBelowValley += a;
        }
        if (c >= valley && !spike_.empty() && isSpike(key)) {
            ++r.spikeSolid;
            spikeCounts.push_back(c);
            if (a) ++r.spikeInAsm;
            return;
        }
        if (c >= valley) {
            const double x = static_cast<double>(c) / peak;
            int b = 0;
            while (!(x >= kEdges[b] && x < kEdges[b + 1])) ++b;
            ++r.solid;
            ++r.solidBin[b];
            if (!a && x >= 3) r.missingHi.push_back({key, c});
            if (a) ++r.solidFound; else ++r.missBin[b];
        }
    });
    r.completeness = r.solid ? static_cast<double>(r.solidFound) / static_cast<double>(r.solid) : 0;
    r.qvRead0 = sqaQv(r.asmRead0, r.asmPos);
    r.qvSolid = sqaQv(r.asmBelowValley, r.asmPos);
    r.spikeRef = spike_.size();
    std::sort(spikeCounts.begin(), spikeCounts.end());
    r.spikeMedianCn = spikeCounts.empty() ? 0 : static_cast<double>(spikeCounts[spikeCounts.size() / 2]) / peak;
    r.alarm = r.missingBpEst(kSqaBins - 1) >= kSqaAlarmBp;
    std::sort(r.missingHi.begin(), r.missingHi.end(), [](const auto& a, const auto& b) {
        if (a.second != b.second) return a.second > b.second;
        return a.first < b.first;
    });

    // per record (kqa .contigs.tsv) and read-absent windows (kqa .absent.bed)
    for (size_t ci = 0; ci < recs_.size(); ++ci) {
        const std::string& s = *recs_[ci];
        SqaContig q;
        q.name = names_[ci];
        q.len = s.size();
        std::vector<uint32_t> cs;
        long ws = -1, we = -1, cnt = 0;
        forEachCanon31(s.data(), s.size(), [&](size_t p, uint64_t k) {
            if (!sqaSampled(k, kSqaSample)) return;
            const uint32_t rc = table_.readCount(k);
            ++q.sampled;
            if (rc == 0) ++q.read0;
            if (rc < valley) ++q.belowValley;
            cs.push_back(rc);
            if (rc == 0) {
                if (ws >= 0 && static_cast<long>(p) <= we + 200) { we = static_cast<long>(p) + kK; ++cnt; }
                else {
                    if (ws >= 0) { r.absentWindows.push_back({q.name, {ws, we}}); r.absentCounts.push_back(cnt); }
                    ws = static_cast<long>(p); we = ws + kK; cnt = 1;
                }
            }
        });
        if (ws >= 0) { r.absentWindows.push_back({q.name, {ws, we}}); r.absentCounts.push_back(cnt); }
        if (!cs.empty()) {
            std::nth_element(cs.begin(), cs.begin() + static_cast<long>(cs.size() / 2), cs.end());
            q.medianCn = static_cast<double>(cs[cs.size() / 2]) / peak;
        }
        r.contigs.push_back(q);
    }
    r.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return r;
}

std::string selfQaSummaryTsv(const SelfQaResult& r) {
    std::string o;
    char b[256];
    auto kv = [&](const char* k, const std::string& v) { o += k; o += '\t'; o += v; o += '\n'; };
    auto u = [&](uint64_t x) { return std::to_string(x); };
    auto f = [&](double x, int p) { std::snprintf(b, sizeof b, "%.*f", p, x); return std::string(b); };
    kv("k", u(kK));
    kv("sample", u(kSqaSample));
    kv("asm_records", u(r.asmRecords));
    kv("asm_bases", u(r.asmBases));
    kv("read_bases", u(r.readBases));
    kv("reads", u(r.reads));
    kv("valley", u(r.valley));
    kv("peak", u(r.peak));
    kv("solid_sampled", u(r.solid));
    kv("completeness", f(r.completeness, 5));
    kv("asm_pos_sampled", u(r.asmPos));
    kv("asm_pos_read0", u(r.asmRead0));
    kv("asm_pos_below_valley", u(r.asmBelowValley));
    kv("qv_read0", f(r.qvRead0, 2));
    kv("qv_solid", f(r.qvSolid, 2));
    kv("spike_kmers_ref", u(r.spikeRef));
    kv("spike_solid", u(r.spikeSolid));
    kv("spike_in_asm", u(r.spikeInAsm));
    kv("spike_median_cn", f(r.spikeMedianCn, 2));
    for (int i = 0; i < kSqaBins; ++i) {
        kv((std::string("solid_") + kSqaBinName[i]).c_str(), u(r.solidBin[i]));
        kv((std::string("missing_") + kSqaBinName[i]).c_str(), u(r.missBin[i]));
        kv((std::string("missing_bp_est_") + kSqaBinName[i]).c_str(), u(r.missingBpEst(i)));
    }
    kv("asm_bases_ge_minlen", u(r.asmBases));
    kv("lost_replicon_alarm", r.alarm ? "1" : "0");
    return o;
}

std::string selfQaVerdict(const SelfQaResult& r) {
    if (!r.ran) return "self-QA did not run" + (r.error.empty() ? std::string() : ": " + r.error);
    char b[512];
    const double hiKb = static_cast<double>(r.missingBpEst(kSqaBins - 1)) / 1000.0;
    const double midKb = static_cast<double>(r.missingBpEst(kSqaBins - 2)) / 1000.0;
    if (r.alarm) {
        std::snprintf(b, sizeof b,
                      "~%.1f kb of sequence present in the reads at >10x the genome depth is absent from the "
                      "assembly: likely a lost small plasmid, rRNA operon or other high-copy element "
                      "(3-10x: ~%.1f kb absent)", hiKb, midKb);
    } else {
        std::snprintf(b, sizeof b,
                      "no lost-replicon alarm: ~%.1f kb of >10x and ~%.1f kb of 3-10x read sequence absent "
                      "(alarm at >= %.1f kb of >10x)", hiKb, midKb, static_cast<double>(kSqaAlarmBp) / 1000.0);
    }
    return b;
}

std::string selfQaJson(const SelfQaResult& r) {
    JsonObj o;
    o.b("ran", r.ran);
    if (!r.ran) {
        o.str("error", r.error);
        return o.done();
    }
    o.str("method", "kqa port: canonical 31-mers hash-sampled 1 in 8; contigs.fasta records >= 500 bp; "
                    "every base of the raw input reads; spike-in reference PhiX174 NC_001422.1");
    o.u("k", kK);
    o.u("sample", kSqaSample);
    o.u("min_record_length", kSqaMinRecord);
    o.u("asm_records", r.asmRecords);
    o.u("asm_bases", r.asmBases);
    o.u("reads", r.reads);
    o.u("read_bases", r.readBases);
    o.u("valley", r.valley);
    o.u("peak", r.peak);
    o.u("solid_sampled", r.solid);
    o.num("completeness", r.completeness, 5);
    o.u("asm_pos_sampled", r.asmPos);
    o.u("asm_pos_read0", r.asmRead0);
    o.u("asm_pos_below_valley", r.asmBelowValley);
    o.num("qv_read0", r.qvRead0, 2);
    o.num("qv_solid", r.qvSolid, 2);
    {
        JsonObj m;
        for (int i = 0; i < kSqaBins; ++i) {
            JsonObj bin;
            bin.u("solid", r.solidBin[i]).u("missing", r.missBin[i]).u("missing_bp_est", r.missingBpEst(i));
            m.raw(kSqaBinName[i], bin.done());
        }
        o.raw("missing_by_depth", m.done());
    }
    {
        JsonObj s;
        s.str("reference", "PhiX174 NC_001422.1");
        s.u("kmers_ref", r.spikeRef).u("solid", r.spikeSolid).u("in_asm", r.spikeInAsm);
        s.num("median_cn", r.spikeMedianCn, 2);
        o.raw("spike_in", s.done());
    }
    {
        JsonObj a;
        a.str("rule", "missing_bp_est_gt10x >= 1000 (EVAL_PLAN_P2 L-SQA-b, kqa flag A)");
        a.u("missing_bp_est_gt10x", r.missingBpEst(kSqaBins - 1));
        a.u("missing_bp_est_3_10x", r.missingBpEst(kSqaBins - 2));
        a.b("fired", r.alarm);
        o.raw("lost_replicon_alarm", a.done());
    }
    o.str("verdict", selfQaVerdict(r));
    o.u("table_distinct", r.tableDistinct);
    o.u("table_bytes", r.tableBytes);
    o.num("finish_seconds", r.seconds, 3);
    return o.done();
}

bool writeSelfQaFiles(const std::string& dir, const SelfQaResult& r, std::string& error) {
    auto put = [&](const std::string& name, const std::string& body) {
        const std::string p = dir + "/" + name;
        std::FILE* f = std::fopen(p.c_str(), "wb");
        if (!f) { error = "cannot write " + p; return false; }
        const bool ok = std::fwrite(body.data(), 1, body.size(), f) == body.size();
        if (std::fclose(f) != 0 || !ok) { error = "write failed on " + p; return false; }
        return true;
    };
    if (!put("p2_self_qa.tsv", selfQaSummaryTsv(r))) return false;
    std::string c = "contig\tlen\tsampled\tread0\tbelow_valley\tmedian_cn\n";
    char b[128];
    for (const SqaContig& q : r.contigs) {
        std::snprintf(b, sizeof b, "\t%zu\t%llu\t%llu\t%llu\t%.3f\n", q.len, static_cast<unsigned long long>(q.sampled),
                      static_cast<unsigned long long>(q.read0), static_cast<unsigned long long>(q.belowValley),
                      q.medianCn);
        c += q.name + b;
    }
    if (!put("p2_self_qa.contigs.tsv", c)) return false;
    std::string m = "kmer\tread_count\tdepth_ratio\n";
    for (const auto& x : r.missingHi) {
        std::snprintf(b, sizeof b, "\t%u\t%.2f\n", x.second, r.peak ? static_cast<double>(x.second) / r.peak : 0.0);
        m += kmerString(x.first) + b;
    }
    if (!put("p2_self_qa.missing_hi.tsv", m)) return false;
    std::string a;
    for (size_t i = 0; i < r.absentWindows.size(); ++i) {
        std::snprintf(b, sizeof b, "\t%ld\t%ld\t%ld\n", r.absentWindows[i].second.first,
                      r.absentWindows[i].second.second, r.absentCounts[i]);
        a += r.absentWindows[i].first + b;
    }
    return put("p2_self_qa.absent.bed", a);
}

}  // namespace p2
}  // namespace ts
