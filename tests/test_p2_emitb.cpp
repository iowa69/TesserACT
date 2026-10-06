// Component tests for phase 2 emit-B (src/p2_*): SHA-256/MD5 known answers, the embedded PhiX
// reference, the k-mer tables, the raw read pass (thread-count invariance, FASTQ/FASTA/gz parity),
// the kqa self-QA port (alarm, completeness, QV, spike-ins), the contig-end audit (tip rule,
// lower case, non-unique tails, pair partners), the detection report (parity with
// om2::scoreReadFiles), the provenance digests, and the flags-off "p2" block.
#include <algorithm>
#include <atomic>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <map>
#include <random>
#include <set>
#include <string>
#include <thread>
#include <vector>
#include <unistd.h>
#include <zlib.h>

#include "test_env.h"
#include "envflags.h"
#include "organism.h"
#include "organism_detect.h"
#include "p2_emitb.h"
#include "p2_ends.h"
#include "p2_kmer.h"
#include "p2_phix.h"
#include "p2_provenance.h"
#include "p2_readpass.h"
#include "p2_selfqa.h"
#include "p2_sha256.h"
#include "seqio.h"

using namespace ts;
using namespace ts::p2;

namespace {

int checks = 0, failures = 0;
void check(bool ok, const std::string& what) {
    ++checks;
    if (!ok) { ++failures; std::fprintf(stderr, "FAIL: %s\n", what.c_str()); }
}
std::mt19937_64 rng(2026);
std::string randomSeq(size_t n) {
    std::string s(n, 'A');
    for (char& c : s) c = "ACGT"[rng() % 4];
    return s;
}
std::string rc(const std::string& s) { return ts::p2::reverseComplement(s); }

void writeFile(const std::string& p, const std::string& body) {
    std::ofstream f(p, std::ios::binary);
    f << body;
}
void writeFastqGz(const std::string& p, const std::vector<std::string>& seqs, const std::string& tag) {
    gzFile g = gzopen(p.c_str(), "wb");
    for (size_t i = 0; i < seqs.size(); ++i) {
        const std::string rec = "@r" + std::to_string(i) + tag + "\n" + seqs[i] + "\n+\n" + std::string(seqs[i].size(), 'I') + "\n";
        gzwrite(g, rec.data(), static_cast<unsigned>(rec.size()));
    }
    gzclose(g);
}
// pairs (r1, r2) of fragments of g; insert uniform in [lo, hi]; circular g wraps
void pairsOf(const std::string& g, size_t n, size_t len, size_t lo, size_t hi, bool circular,
             std::vector<std::string>& r1, std::vector<std::string>& r2) {
    const std::string gg = circular ? g + g.substr(0, hi + 1) : g;
    for (size_t i = 0; i < n; ++i) {
        const size_t ins = lo + rng() % (hi - lo + 1);
        const size_t p = rng() % ((circular ? g.size() : g.size() - ins) + 0);
        const std::string frag = gg.substr(p, ins);
        if (frag.size() < ins) continue;
        if (rng() % 2) { r1.push_back(frag.substr(0, len)); r2.push_back(rc(frag.substr(ins - len))); }
        else { r2.push_back(frag.substr(0, len)); r1.push_back(rc(frag.substr(ins - len))); }
    }
}

void testDigests(const std::string& dir) {
    check(sha256Hex("") == "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855", "sha256('')");
    check(sha256Hex("abc") == "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad", "sha256(abc)");
    check(sha256Hex("abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq") ==
              "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1", "sha256(448-bit vector)");
    {
        Sha256 h;
        const std::string a(1000, 'a');
        for (int i = 0; i < 1000; ++i) h.add(a.data(), a.size());
        check(h.finish() == "cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e046d39ccc7112cd0", "sha256(1e6 x 'a')");
    }
    {   // odd chunking equals one shot
        const std::string m = randomSeq(10007);
        Sha256 h;
        size_t at = 0, step = 1;
        while (at < m.size()) { const size_t n = std::min(step, m.size() - at); h.add(m.data() + at, n); at += n; step = step * 3 + 1; }
        check(h.finish() == sha256Hex(m), "sha256 streaming == one shot");
    }
    writeFile(dir + "/abc.txt", "abc");
    std::string hex;
    uint64_t bytes = 0;
    check(sha256File(dir + "/abc.txt", hex, bytes) && bytes == 3 &&
              hex == "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad", "sha256File(abc)");
    check(md5FileHead(dir + "/abc.txt", 1u << 20) == "900150983cd24fb0d6963f7d28e17f72", "md5 head of abc");
    check(md5FileHead(dir + "/abc.txt", 2) == "187ef4436122d1cc2f40dc2b92f0eba0", "md5 of the first 2 bytes (ab)");
    check(md5File(dir + "/abc.txt") == "900150983cd24fb0d6963f7d28e17f72", "md5File(abc)");
    check(!sha256File(dir + "/missing.txt", hex, bytes) && hex.empty(), "missing file refused");
    // the embedded PhiX174
    const std::string phix(kPhiX174);
    writeFile(dir + "/phix.txt", phix);
    check(phix.size() == 5386, "PhiX174 length 5386");
    check(md5File(dir + "/phix.txt") == kPhiX174Md5, "PhiX174 md5 matches NC_001422.1 (phix.fa)");
    // provenance digests
    const FileDigest d = digestFile("-1", dir + "/abc.txt", true);
    check(d.size == 3 && d.sha256 == "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad" &&
              d.md5Head1m == "900150983cd24fb0d6963f7d28e17f72" && !d.realpath.empty() && d.error.empty(),
          "digestFile: size, sha256, md5_head1m, realpath");
    const FileDigest e = digestFile("x", dir + "/missing.txt", false);
    check(!e.error.empty() && e.sha256.empty(), "digestFile on a missing file records an error");
}

void testKmers() {
    const std::string s = randomSeq(500);
    size_t n = 0;
    bool same = true;
    forEachCanon31(s.data(), s.size(), [&](size_t p, uint64_t c) {
        ++n;
        if (canon31At(s, p) != c) same = false;
    });
    check(n == 470 && same, "forEachCanon31 == canon31At at every position");
    const std::string r = rc(s);
    check(canon31At(s, 0) == canon31At(r, r.size() - 31), "canonical k-mer is strand independent");
    std::string low = s;
    for (char& c : low) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    check(canon31At(low, 17) == canon31At(s, 17), "lower case gives the same k-mer");
    std::string withN = s;
    withN[40] = 'N';
    size_t m = 0;
    forEachCanon31(withN.data(), withN.size(), [&](size_t, uint64_t) { ++m; });
    check(m == 470 - 31, "an N breaks every k-mer that covers it");
    check(canon31At(withN, 20) == UINT64_MAX, "canon31At over N is invalid");
    check(kmerString(canon31At(std::string(31, 'A'), 0)) == std::string(31, 'A'), "kmerString(AAA..)");

    FixedCountTable t;
    t.addKey(5); t.addKey(9); t.addKey(5);
    t.freeze();
    check(t.size() == 2 && t.contains(5) && !t.contains(7), "fixed table keys");
    std::vector<std::thread> pool;
    for (int i = 0; i < 4; ++i) pool.emplace_back([&]() { for (int j = 0; j < 1000; ++j) { t.hit(5); t.hit(7); } });
    for (auto& x : pool) x.join();
    check(t.count(5) == 4000 && t.count(7) == 0 && t.count(9) == 0, "fixed table: concurrent hits counted exactly");

    ShardedCountTable st;
    st.addAssembly(11);
    st.addAssembly(11);
    std::vector<std::thread> pool2;
    for (int i = 0; i < 4; ++i)
        pool2.emplace_back([&, i]() {
            std::vector<uint64_t> keys;
            for (uint64_t k = 0; k < 20000; ++k) keys.push_back(k * 7919 + static_cast<uint64_t>(i % 2));
            std::vector<std::vector<uint64_t>> by(ShardedCountTable::kShards);
            for (uint64_t k : keys) by[ShardedCountTable::shardOf(k)].push_back(k);
            for (size_t sh = 0; sh < by.size(); ++sh) st.addReadBatch(sh, by[sh].data(), by[sh].size());
        });
    for (auto& x : pool2) x.join();
    uint64_t total = 0;
    uint16_t ac11 = 0;
    st.forEach([&](uint64_t k, uint32_t rcount, uint16_t ac) { total += rcount; if (k == 11) ac11 = ac; });
    check(total == 80000 && st.readCount(7919) == 2 && st.readCount(7920) == 2 && ac11 == 2 && st.distinct() == 40001,
          "sharded table: exact counts under concurrent batches");
}

void testReadPass(const std::string& dir) {
    const std::string g = randomSeq(20000);
    std::vector<std::string> r1, r2;
    pairsOf(g, 4000, 100, 250, 350, false, r1, r2);
    writeFastqGz(dir + "/rp_1.fq.gz", r1, "/1");
    // the second file as multi-line FASTA, uncompressed
    std::string fa;
    for (size_t i = 0; i < r2.size(); ++i) fa += ">r" + std::to_string(i) + "/2\n" + r2[i].substr(0, 60) + "\n" + r2[i].substr(60) + "\n";
    writeFile(dir + "/rp_2.fa", fa);
    std::vector<std::string> all = r1;
    all.insert(all.end(), r2.begin(), r2.end());

    std::vector<uint64_t> keys;
    for (size_t p = 0; p + 31 <= 2000; p += 3) keys.push_back(canon31At(g, p));
    auto run = [&](int threads, bool inMemory, std::vector<uint32_t>& counts, std::vector<uint64_t>& hist, ReadPassStats& st) {
        FixedCountTable exact;
        for (uint64_t k : keys) exact.addKey(k);
        exact.freeze();
        ShardedCountTable sampled;
        ReadPassSinks s;
        s.exact = &exact;
        s.sampled = &sampled;
        if (inMemory) countSequences(all, s, st);
        else runReadPass({dir + "/rp_1.fq.gz", dir + "/rp_2.fa"}, threads, s, st);
        counts.clear();
        for (uint64_t k : keys) counts.push_back(exact.count(k));
        hist.assign(64, 0);
        sampled.forEach([&](uint64_t, uint32_t c, uint16_t) { ++hist[std::min<uint32_t>(c, 63)]; });
    };
    std::vector<uint32_t> c1, c4, cm;
    std::vector<uint64_t> h1, h4, hm;
    ReadPassStats s1, s4, sm;
    run(1, false, c1, h1, s1);
    run(4, false, c4, h4, s4);
    run(1, true, cm, hm, sm);
    check(s1.error.empty() && s1.reads == all.size() && s1.bases == all.size() * 100, "read pass: every record of FASTQ.gz + FASTA read");
    check(c1 == c4 && h1 == h4, "read pass: identical counts at 1 and 4 threads");
    check(c1 == cm && h1 == hm && sm.reads == s1.reads, "read pass on files == in-memory counting");
    ReadPassStats bad;
    ReadPassSinks none;
    FixedCountTable ex;
    ex.addKey(1);
    ex.freeze();
    none.exact = &ex;
    check(!runReadPass({dir + "/nonexistent.fq.gz"}, 2, none, bad) && !bad.error.empty(), "read pass: a missing file is an error");
}

// the kqa port against a brute-force recomputation of its own definitions
void testSelfQa() {
    const std::string chrom = randomSeq(60000);
    const std::string plasmid = randomSeq(4000);     // high copy, left out of the assembly
    std::vector<std::string> reads;
    auto sample = [&](const std::string& g, double depth) {
        const size_t n = static_cast<size_t>(depth * g.size() / 100);
        for (size_t i = 0; i < n; ++i) {
            std::string r = g.substr(rng() % (g.size() - 100), 100);
            if (rng() % 50 == 0) r[rng() % 100] = "ACGT"[rng() % 4];   // occasional error
            reads.push_back(rng() % 2 ? rc(r) : r);
        }
    };
    sample(chrom, 40);
    sample(plasmid, 600);
    std::string phix(kPhiX174);
    sample(phix, 30);
    auto runSqa = [&](const std::vector<std::string>& recs) {
        SelfQa q;
        q.prepare(recs, {}, {phix});
        ReadPassSinks s;
        s.sampled = &q.table();
        ReadPassStats st;
        countSequences(reads, s, st);
        return q.finish(st.reads, st.bases);
    };
    const SelfQaResult a = runSqa({chrom});
    check(a.ran && a.asmRecords == 1 && a.asmBases == chrom.size(), "self-QA: one record");
    check(a.peak >= 25 && a.peak <= 45, "self-QA: spectrum peak near the 31-mer depth (" + std::to_string(a.peak) + ")");
    check(a.alarm && a.missingBpEst(kSqaBins - 1) >= 3000 && a.missingBpEst(kSqaBins - 1) <= 5000,
          "self-QA: a 4-kb 600x plasmid absent from the assembly fires the alarm (" +
              std::to_string(a.missingBpEst(kSqaBins - 1)) + " bp)");
    // (the last bases of the chromosome are never at the end of a sampled read: a k-mer or two read-absent)
    check(a.asmRead0 <= 3 && a.qvRead0 > 50 && a.completeness > 0.85 && a.completeness < 0.99,
          "self-QA: QV high, completeness < 1 without the plasmid (qv " + std::to_string(a.qvRead0) + ", completeness " +
              std::to_string(a.completeness) + ", valley " + std::to_string(a.valley) + ", read0 " + std::to_string(a.asmRead0) + ")");
    check(a.spikeRef > 500 && a.spikeSolid > 0.9 * a.spikeRef && a.spikeInAsm == 0, "self-QA: PhiX spike-in k-mers solid and absent");
    // brute force: missing >10x sampled k-mers
    {
        std::set<uint64_t> asmK;
        forEachCanon31(chrom.data(), chrom.size(), [&](size_t, uint64_t c) { if (sqaSampled(c, 8)) asmK.insert(c); });
        std::set<uint64_t> spikeK;
        forEachCanon31(phix.data(), phix.size(), [&](size_t, uint64_t c) { if (sqaSampled(c, 8)) spikeK.insert(c); });
        std::map<uint64_t, uint32_t> cnt;
        for (const std::string& r : reads)
            forEachCanon31(r.data(), r.size(), [&](size_t, uint64_t c) { if (sqaSampled(c, 8)) ++cnt[c]; });
        uint64_t miss10 = 0, solid = 0, found = 0;
        for (const auto& kv : cnt) {
            if (kv.second < a.valley || spikeK.count(kv.first)) continue;
            ++solid;
            if (asmK.count(kv.first)) ++found;
            else if (static_cast<double>(kv.second) / a.peak >= 10) ++miss10;
        }
        check(miss10 == a.missBin[kSqaBins - 1] && solid == a.solid && found == a.solidFound,
              "self-QA: missing >10x, solid and found equal a brute-force count");
    }
    const SelfQaResult b = runSqa({chrom, plasmid});
    check(!b.alarm && b.missingBpEst(kSqaBins - 1) < 1000 && b.completeness > a.completeness, "self-QA: no alarm with the plasmid assembled");
    std::string bad = chrom;
    for (size_t p = 1000; p < 60000; p += 997) bad[p] = bad[p] == 'A' ? 'C' : 'A';
    const SelfQaResult c = runSqa({bad});
    check(c.qvRead0 < 60 && c.asmRead0 > 0, "self-QA: substitutions lower the QV (" + std::to_string(c.qvRead0) + ")");
    check(sqaQv(0, 100) == 99.0 && sqaQv(1, 0) == 99.0, "QV cap");
    const std::string tsv = selfQaSummaryTsv(a);
    check(tsv.find("missing_bp_est_gt10x\t") != std::string::npos && tsv.find("lost_replicon_alarm\t1") != std::string::npos,
          "self-QA summary keys (kqa names)");
    const std::string js = selfQaJson(a);
    check(js.find("\"lost_replicon_alarm\": {") != std::string::npos && js.find("\"fired\": true") != std::string::npos,
          "self-QA JSON carries the alarm");
    const SelfQaResult tiny = runSqa({randomSeq(499)});
    check(tiny.asmRecords == 0, "self-QA: records under 500 bp are not part of the assembly (kqa -m 500)");
}

void testTipRule() {
    check(medianOf({}) == 0 && medianOf({3, 1, 2}) == 2 && medianOf({4, 1, 3, 2}) == 2.5, "median (statistics.median)");
    check(tipFromCounts({0, 0, 1, 9, 9}, 20) == 3, "tip: three unsupported, then a vouched k-mer (thr 5)");
    check(tipFromCounts({4, 9}, 20) == 1, "tip: count below 0.25 x depth is unsupported");
    check(tipFromCounts({2, 0}, 4) == 0, "tip: floor of 2 reads vouches at low depth");
    check(tipFromCounts({1, 1, 1}, 0) == 3, "tip: nothing vouched -> every scanned position");
    check(tipFromCounts({}, 30) == 0, "tip: no k-mers");
}

void testEnds(const std::string& dir) {
    // A genome cut into records A | gap | B; C = a repeat copy shared by two records' ends.
    const std::string g = randomSeq(9000);
    const std::string A = g.substr(0, 4000);
    const std::string B = g.substr(4150, 4850);
    std::vector<std::string> r1, r2;
    pairsOf(g, 9000, 100, 280, 320, false, r1, r2);
    writeFastqGz(dir + "/e_1.fq.gz", r1, "/1");
    writeFastqGz(dir + "/e_2.fq.gz", r2, "/2");
    SequenceStore store;
    std::string err;
    Library lib;
    lib.r1 = dir + "/e_1.fq.gz";
    lib.r2 = dir + "/e_2.fq.gz";
    QualityTrim qt;
    qt.enabled = false;
    store.setQualityTrim(qt);
    check(store.load({lib}, 2, err) && store.pairCount() > 0, "ends: load the pairs (" + err + ")");

    // record A with an error 5 bases from its right end; B with a 120-bp tail shared with D
    std::string Aerr = A;
    Aerr[A.size() - 5] = Aerr[A.size() - 5] == 'A' ? 'C' : 'A';
    const std::string shared = B.substr(B.size() - 120);
    const std::string D = randomSeq(900) + shared;
    std::vector<std::string> seqs = {Aerr, B, D, randomSeq(300)};
    std::vector<std::string> names = {"NODE_1_length_4000_cov_30.0_chr", "NODE_2_length_4850_cov_30.0_chr",
                                      "NODE_3_length_1020_cov_30.0_unk", "NODE_4_length_300_cov_30.0_unk"};
    EndsInput in;
    in.seqs = &seqs;
    in.names = &names;
    std::vector<double> covs = {30, 30, 30, 30};
    in.covs = &covs;
    in.reads = &store;
    in.insert.mean = 300;
    in.insert.stddev = 12;
    in.insert.maxPlausible = 400;
    in.insert.usable = true;
    in.threads = 3;
    EndsAudit audit(in);
    audit.computeStructure();
    FixedCountTable raw;
    audit.addRawKeys(raw);
    raw.freeze();
    ReadPassSinks s;
    s.exact = &raw;
    s.wantBits = &audit.poolBits();
    s.onWanted = [&](uint64_t idx, const char* q, size_t n) { audit.countPoolRead(idx, q, n); };
    std::vector<std::string> files;
    readFileLayout({lib}, files, s.layout);
    check(files.size() == 2 && s.layout[0] == std::make_pair(uint64_t(0), uint64_t(2)) &&
              s.layout[1] == std::make_pair(uint64_t(1), uint64_t(2)), "ends: store layout of a paired library");
    ReadPassStats st;
    check(runReadPass(files, 2, s, st) && st.wanted > 0 && st.fileRecords.size() == 2 &&
              st.fileRecords[0] * 2 == store.size(), "ends: raw pass collects the locus-pool reads");
    audit.finishTips(raw, "raw");
    std::vector<EndRow> rows = audit.rows();
    check(rows.size() == 6, "ends: two rows per record >= 500 bp (300-bp record has none)");
    auto row = [&](size_t rec, char e) -> const EndRow& {
        for (const EndRow& r : rows) if (r.record == rec && r.end == e) return r;
        return rows[0];
    };
    check(row(0, 'R').tip == 5, "ends: an error 5 bp from the end gives a 5-bp unsupported tip (" + std::to_string(row(0, 'R').tip) + ")");
    check(row(0, 'L').tip > 0, "ends: the genome start is a dead end the reads thin out at");
    check(row(1, 'L').tip == 0, "ends: a supported end has no tip");
    check(row(1, 'R').nonuniqueTail == 120 && row(2, 'R').nonuniqueTail == 120 && row(1, 'L').nonuniqueTail == 0,
          "ends: a 120-bp tail shared with another record is the non-unique tail");
    check(row(0, 'R').partnerClass == "unique" && row(0, 'R').bestPartner.find("NODE_2_length_4850_cov_30.0_chr:L:") == 0,
          "ends: A's right end has B's left end as its unique partner (" + row(0, 'R').partnerClass + " " + row(0, 'R').bestPartner + ")");
    check(row(1, 'L').partnerClass == "unique" && row(1, 'L').outward > 0, "ends: and the reverse");
    check(row(0, 'R').cn > 0.99 && row(0, 'R').cn < 1.01 && row(0, 'R').anchor, "ends: copy number 1 anchor");
    // lower case
    EndsStats es = audit.stats();
    std::vector<std::string> edits;
    const std::vector<std::string> lc = lowercaseTips(seqs, names, rows, es, edits);
    size_t lower = 0;
    for (const std::string& q : lc) for (char c : q) if (std::islower(static_cast<unsigned char>(c))) ++lower;
    size_t tips = 0;
    for (const EndRow& r : rows) tips += r.tip;
    check(lower == tips && es.lowercasedBases == tips && edits.size() == es.tipEnds, "lower case: exactly the tip bases, one edit row per tipped end");
    bool runsOk = true;
    for (const EndRow& r : rows) {
        const std::string& q = lc[r.record];
        if (r.end == 'R') {
            for (size_t i = 0; i < q.size(); ++i) {
                const bool isLow = std::islower(static_cast<unsigned char>(q[i])) != 0;
                if (i >= q.size() - r.tip && !isLow) runsOk = false;
            }
        }
    }
    std::string upper = lc[0];
    for (char& c : upper) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    check(runsOk && upper == seqs[0] && lc[3] == seqs[3], "lower case: suffix run of tip bases; upper-cased record unchanged");
    const std::string tsv = endsTsv(rows, names, seqs);
    check(tsv.find("record\tend\tlength\ttopology\tcn\tanchor") == 0 &&
              std::count(tsv.begin(), tsv.end(), '\n') == 7, "ends.tsv: header + one line per end");
    // the same audit from the store's reads (the fallback) agrees on this error-free library
    {
        EndsAudit b(in);
        b.computeStructure();
        b.countPoolsFromStore();
        b.finishTips(raw, "store_corrected");
        bool sameTips = true;
        for (size_t i = 0; i < rows.size(); ++i) sameTips = sameTips && b.rows()[i].tip == rows[i].tip && b.rows()[i].pool == rows[i].pool;
        check(sameTips && b.stats().poolSource == "store_corrected", "ends: store-read pools give the same tips");
    }
    check(row(0, 'R').tipStatus == "audited" && row(0, 'R').pool >= kEndsMinPool, "ends: anchor ends with a pool are audited");
    // no usable insert model: partner class na
    EndsInput in2 = in;
    in2.insert.usable = false;
    EndsAudit a2(in2);
    a2.computeStructure();
    check(a2.rows()[0].partnerClass == "na" && a2.stats().na == 6, "ends: partners na without an insert model");
}

void testDetect(const std::string& dir) {
    // two synthetic organisms; the sketch's cores are their 1/512-sampled 31-mers (denser for B)
    const std::string ga = randomSeq(400000), gb = randomSeq(300000);
    om2::DetectSketch sk;
    om2::DetectSketchModel ma, mb;
    ma.name = "kpneumoniae"; ma.denom = 512; ma.tsmMd5 = "a"; ma.holdList = "-"; ma.holdMd5 = "-";
    mb.name = "saureus"; mb.denom = 64; mb.tsmMd5 = "b"; mb.holdList = "-"; mb.holdMd5 = "-";
    forEachMarkerKmer(ga, [&](uint64_t km, uint32_t, int) { ma.core.push_back(km); }, 512);
    forEachMarkerKmer(gb, [&](uint64_t km, uint32_t, int) { mb.core.push_back(km); }, 64);
    for (auto* m : {&ma, &mb}) { std::sort(m->core.begin(), m->core.end()); m->core.erase(std::unique(m->core.begin(), m->core.end()), m->core.end()); }
    sk.models = {ma, mb};
    std::string err;
    check(om2::writeDetectSketch(dir + "/d.sketch", sk, err), "detect: write sketch");
    std::vector<std::string> r1, r2;
    pairsOf(ga, 12000, 150, 300, 400, false, r1, r2);
    pairsOf(gb, 300, 150, 300, 400, false, r1, r2);   // a little of B
    writeFastqGz(dir + "/d_1.fq.gz", r1, "/1");
    writeFastqGz(dir + "/d_2.fq.gz", r2, "/2");
    om2::DetectResult ref;
    check(om2::scoreReadFiles(sk, {dir + "/d_1.fq.gz", dir + "/d_2.fq.gz"}, ref, err), "detect: reference scoring");

    setenv("TESSERACT_P2_DETECT_REPORT", "1", 1);
    setenv("TESSERACT_P2_DETECT_SKETCH", (dir + "/d.sketch").c_str(), 1);
    EmitBState st;
    st.opt = EmitBOptions::fromEnv();
    unsetenv("TESSERACT_P2_DETECT_REPORT");
    unsetenv("TESSERACT_P2_DETECT_SKETCH");
    check(st.opt.detect && !st.opt.endsAudit() && !st.opt.selfQa && st.opt.detectSketch == dir + "/d.sketch", "detect: options");
    const std::vector<std::string> recs = {randomSeq(1000)}, names = {"NODE_1"};
    EmitBRecords in;
    in.outDir = dir;
    in.seqs = &recs;
    in.names = &names;
    in.covs = {30};
    in.threads = 3;
    Library dl;
    dl.r1 = dir + "/d_1.fq.gz";
    dl.r2 = dir + "/d_2.fq.gz";
    in.libraries = {dl};
    std::vector<std::string> written;
    const bool replaced = runEmitBRecords(in, st, written);
    check(!replaced && st.det.ran && st.det.scores.size() == 2, "detect: ran, nothing rewritten");
    bool same = st.det.scores.size() == ref.scores.size();
    for (size_t i = 0; same && i < ref.scores.size(); ++i)
        same = st.det.scores[i].name == ref.scores[i].name && st.det.scores[i].coreHit == ref.scores[i].coreHit &&
               st.det.scores[i].core == ref.scores[i].core;
    check(same, "detect: scores identical to om2::scoreReadFiles (both densities)");
    check(st.det.accepted() && st.det.scores[0].name == "kpneumoniae", "detect: the right call, accepted");
    const std::string js = p2BlockJson(st);
    check(js.find("\"organism_detect\": {\"mode\": \"report_only\"") != std::string::npos &&
              js.find("\"model_applied\": false") != std::string::npos && js.find("\"call\": \"kpneumoniae\"") != std::string::npos,
          "detect: report-only block");
    // no sketch: recorded, no crash
    EmitBState s2;
    s2.opt.detect = true;
    runEmitBRecords(in, s2, written);
    check(!s2.det.ran && !s2.det.error.empty(), "detect: a missing sketch is reported, not fatal");
}

void testBlockAndFlags(const std::string& dir) {
    EmitBState off;
    off.opt = EmitBOptions::fromEnv();
    check(!off.opt.any(), "flags off by default");
    const std::string js = p2BlockJson(off);
    check(js == "{\"enabled\": 0, \"counters\": {\"ends_written\": 0, \"end_rows\": 0, \"tip_ends\": 0, \"tip_bases\": 0, "
                "\"lowercased_bases\": 0, \"self_qa_run\": 0, \"self_qa_alarm\": 0, \"provenance_files\": 0, "
                "\"detect_run\": 0, \"read_pass_reads\": 0}}",
          "flags off: the p2 block holds only enabled=0 and zero counters: " + js);
    std::FILE* log = std::tmpfile();
    logEmitBCounters(off, log);
    std::rewind(log);
    std::string lines;
    char buf[1024];
    while (std::fgets(buf, sizeof buf, log)) lines += buf;
    std::fclose(log);
    check(lines.find("[p2-ends] enabled=0 lowercase=0 records=0") != std::string::npos &&
              lines.find("[p2-sqa] enabled=0 ran=0") != std::string::npos && lines.find("[p2-prov] enabled=0 files=0 bytes=0") != std::string::npos &&
              lines.find("[p2-detect] enabled=0 ran=0") != std::string::npos && lines.find("[p2-readpass] enabled=0 files=0 reads=0 bases=0") != std::string::npos,
          "flags off: every counter line with enabled=0 and zeros");
    std::vector<std::string> written;
    const std::vector<std::string> recs = {randomSeq(800)}, names = {"NODE_1"};
    EmitBRecords in;
    in.outDir = dir + "/off_out";
    std::filesystem::create_directories(in.outDir);
    in.seqs = &recs;
    in.names = &names;
    check(!runEmitBRecords(in, off, written) && std::filesystem::is_empty(in.outDir), "flags off: no file written, nothing rewritten");
    // a malformed value is rejected by the registry
    const env::Spec* sp = env::find("TESSERACT_P2_TIPS_LOWERCASE");
    env::Parsed p;
    std::string why;
    check(sp && !env::parse(*sp, "yes", p, why) && env::parse(*sp, "1", p, why) && p.on, "flag registry: binary 0/1");
    // setVariables sees a set flag
    setenv("TESSERACT_P2_PROVENANCE", "1", 1);
    bool seen = false;
    for (const env::SetVariable& v : env::setVariables())
        if (v.name == "TESSERACT_P2_PROVENANCE" && v.registered && v.valid && v.kind == "binary") seen = true;
    unsetenv("TESSERACT_P2_PROVENANCE");
    check(seen, "setVariables lists a set TESSERACT_* flag with its kind");
}

void testProvenance(const std::string& dir) {
    writeFile(dir + "/in1.fq", "@a\nACGT\n+\nIIII\n");
    ProvenanceInput in;
    in.inputs = {{"-1", dir + "/in1.fq"}};
    in.namedFiles = {{"model", ""}, {"--qc", dir + "/abc.txt"}};
    in.defaultsLine = "[defaults] x";
    in.mode = "standard";
    std::string js;
    {
        ProvenanceJob job;
        job.start(in);
        js = job.finishJson(4, {21, 33, 55});
        check(job.filesHashed() == 3, "provenance: input, binary and one named file hashed");
    }
    std::string sha;
    uint64_t n = 0;
    sha256File(dir + "/in1.fq", sha, n);
    check(js.find("\"sha256\": \"" + sha + "\"") != std::string::npos && js.find("\"size\": 15") != std::string::npos,
          "provenance: input sha256 and size");
    check(js.find("\"binary\": {\"role\": \"binary\"") != std::string::npos && js.find("\"k_ladder\": [21, 33, 55]") != std::string::npos &&
              js.find("\"threads\": 4") != std::string::npos && js.find("\"git_commit\"") != std::string::npos &&
              js.find("\"argv\": [") != std::string::npos,
          "provenance: binary, k ladder, threads, versions, argv");
    { ProvenanceJob never; }   // destroying an unstarted job is fine
    {
        ProvenanceJob abandoned;
        abandoned.start(in);   // destructor joins
    }
    check(true, "provenance: a started job is joined on destruction");
}

}  // namespace

int main() {
    testenv::clearTesseractEnv();
    const std::string dir = (std::filesystem::temp_directory_path() / ("tesseract-p2-emitb-" + std::to_string(getpid()))).string();
    std::filesystem::create_directories(dir);
    testDigests(dir);
    testKmers();
    testReadPass(dir);
    testSelfQa();
    testTipRule();
    testEnds(dir);
    testDetect(dir);
    testBlockAndFlags(dir);
    testProvenance(dir);
    std::filesystem::remove_all(dir);
    std::printf("test_p2_emitb: %d checks, %d failures\n", checks, failures);
    return failures ? 1 : 0;
}
