#include "p2_emitb.h"

#include <algorithm>
#include <atomic>
#include <memory>
#include <numeric>
#include <unordered_map>

#include "envflags.h"
#include "p2_json.h"
#include "p2_phix.h"
#include "util.h"

namespace ts {
namespace p2 {

namespace {

bool putFile(const std::string& path, const std::string& body, std::string& error) {
    std::FILE* f = std::fopen(path.c_str(), "wb");
    if (!f) { error = "cannot write " + path; return false; }
    const bool ok = std::fwrite(body.data(), 1, body.size(), f) == body.size();
    if (std::fclose(f) != 0 || !ok) { error = "write failed on " + path; return false; }
    return true;
}

}  // namespace

EmitBOptions EmitBOptions::fromEnv() {
    EmitBOptions o;
    o.ends = env::on("TESSERACT_P2_ENDS", false);
    o.lowercase = env::on("TESSERACT_P2_TIPS_LOWERCASE", false);
    o.selfQa = env::on("TESSERACT_P2_SELF_QA", false);
    o.provenance = env::on("TESSERACT_P2_PROVENANCE", false);
    o.detect = env::on("TESSERACT_P2_DETECT_REPORT", false);
    if (o.detect) {
        if (const char* p = env::text("TESSERACT_P2_DETECT_SKETCH")) o.detectSketch = p;
        else if (const char* q = env::text("TESSERACT_OM2_DETECT_SKETCH")) o.detectSketch = q;
        else if (const char* d = env::text("TESSERACT_MODEL_DIR")) o.detectSketch = std::string(d) + "/om2detect.sketch";
    }
    return o;
}

void readFileLayout(const std::vector<Library>& libs, std::vector<std::string>& files,
                    std::vector<std::pair<uint64_t, uint64_t>>& layout) {
    files.clear();
    layout.clear();
    for (const Library& l : libs) {
        if (!l.r1.empty()) files.push_back(l.r1);
        if (!l.r2.empty()) files.push_back(l.r2);
    }
    layout.assign(files.size(), {0, 0});
    // SequenceStore::load: a paired library's R1 record j is read 2j and its R2 record read 2j+1; an
    // interleaved or single-end file's record j is read j. With more than one library the bases depend
    // on record counts not known before reading, so nothing is matched.
    if (libs.size() != 1 || files.empty()) return;
    if (!libs[0].r2.empty() && files.size() == 2) layout = {{0, 2}, {1, 2}};
    else layout[0] = {0, 1};
}

bool runEmitBRecords(const EmitBRecords& in, EmitBState& st, std::vector<std::string>& written) {
    const EmitBOptions& o = st.opt;
    if (!o.needsReadPass() || !in.seqs) return false;

    // ---- consumers of the raw read pass --------------------------------------------------------
    std::unique_ptr<EndsAudit> audit;
    FixedCountTable endKmers;
    ReadPassSinks sinks;
    if (o.endsAudit()) {
        EndsInput ei;
        ei.seqs = in.seqs;
        ei.names = in.names;
        ei.covs = &in.covs;
        ei.reads = in.reads;
        ei.insert = in.insert;
        ei.threads = in.threads;
        audit.reset(new EndsAudit(ei));
        audit->computeStructure();
        audit->addRawKeys(endKmers);
        endKmers.freeze();
        sinks.exact = &endKmers;
        // the locus pools' raw reads, fetched by their store index while the files stream past
        sinks.wantBits = &audit->poolBits();
        EndsAudit* a = audit.get();
        sinks.onWanted = [a](uint64_t idx, const char* seq, size_t len) { a->countPoolRead(idx, seq, len); };
    }
    std::vector<std::string> readFiles;
    readFileLayout(in.libraries, readFiles, sinks.layout);
    std::unique_ptr<SelfQa> sqa;
    if (o.selfQa) {
        sqa.reset(new SelfQa);
        sqa->prepare(*in.seqs, *in.names, {std::string(kPhiX174)});
        sinks.sampled = &sqa->table();
        sinks.sampleS = kSqaSample;
    }
    om2::DetectSketch sketch;
    std::unordered_map<uint64_t, uint32_t> slot;
    std::unique_ptr<std::atomic<uint32_t>[]> slotCounts;
    if (o.detect) {
        st.det.sketchPath = o.detectSketch;
        std::string err;
        if (o.detectSketch.empty()) {
            st.det.error = "no detection sketch (set TESSERACT_P2_DETECT_SKETCH, TESSERACT_OM2_DETECT_SKETCH or "
                           "TESSERACT_MODEL_DIR)";
        } else if (!om2::loadDetectSketch(o.detectSketch, sketch, err)) {
            st.det.error = err;
        } else {
            st.det.sketchMd5 = om2::md5OfFile(o.detectSketch);
            // organism_detect.cpp buildSlots(): one slot per distinct core k-mer, in model order
            std::vector<uint64_t> denoms;
            for (const om2::DetectSketchModel& m : sketch.models) {
                if (std::find(denoms.begin(), denoms.end(), static_cast<uint64_t>(m.denom)) == denoms.end())
                    denoms.push_back(m.denom);
                for (uint64_t km : m.core) slot.emplace(km, static_cast<uint32_t>(slot.size()));
            }
            slotCounts.reset(new std::atomic<uint32_t>[slot.size() + 1]);
            for (size_t i = 0; i <= slot.size(); ++i) slotCounts[i].store(0, std::memory_order_relaxed);
            sinks.detectSlot = &slot;
            for (uint64_t d : denoms) sinks.detectThresholds.push_back(~0ULL / (d ? d : 512));
            sinks.detectCounts = slotCounts.get();
        }
    }

    // ---- the pass ----------------------------------------------------------------------------
    bool passOk = false;
    if (sinks.any()) {
        passOk = runReadPass(readFiles, in.threads, sinks, st.rp);
        st.readPassRan = true;
    }
    if (!passOk) {
        const std::string why = st.rp.error.empty() ? std::string("read pass did not run") : st.rp.error;
        if (o.endsAudit()) st.endsError = why;
        if (o.selfQa) { st.sqa.ran = false; st.sqa.error = why; }
        if (o.detect && st.det.error.empty()) st.det.error = why;
        return false;
    }

    // ---- results -----------------------------------------------------------------------------
    bool replaced = false;
    if (audit) {
        // the raw records match the store only if every collected file held the records the store did
        bool rawPools = in.reads != nullptr && !sinks.layout.empty();
        for (size_t f = 0; rawPools && f < sinks.layout.size(); ++f) {
            const uint64_t stride = sinks.layout[f].second;
            if (stride == 0) { rawPools = false; break; }
            const uint64_t expect = in.reads->size() / stride;
            if (f >= st.rp.fileRecords.size() || st.rp.fileRecords[f] != expect) rawPools = false;
        }
        if (!rawPools) audit->countPoolsFromStore();
        audit->finishTips(endKmers, rawPools ? "raw" : "store_corrected");
        std::vector<EndRow> rows = audit->rows();
        st.ends = audit->stats();
        if (o.lowercase) {
            written = lowercaseTips(*in.seqs, *in.names, rows, st.ends, st.edits);
            replaced = true;
            std::string body = "record\tfeature\toperation\tstart\tend\tbases\n";
            for (const std::string& e : st.edits) body += e + "\n";
            if (!putFile(in.outDir + "/p2_edits.tsv", body, st.fileError)) {
                // the edit log is part of the feature: without it no lower case is written
                written.clear();
                replaced = false;
                for (EndRow& r : rows) r.lowercased = 0;
                st.ends.lowercasedBases = 0;
                st.edits.clear();
            }
        }
        std::string err;
        if (putFile(in.outDir + "/ends.tsv", endsTsv(rows, *in.names, *in.seqs), err)) {
            st.endsRan = true;
        } else {
            st.endsError = err;
            if (st.fileError.empty()) st.fileError = err;
        }
    }
    if (sqa) {
        st.sqa = sqa->finish(st.rp.reads, st.rp.bases);
        std::string err;
        if (!writeSelfQaFiles(in.outDir, st.sqa, err) && st.fileError.empty()) st.fileError = err;
    }
    if (o.detect && st.det.error.empty()) {
        std::vector<om2::DetectModelScore> sc;
        for (const om2::DetectSketchModel& m : sketch.models) {
            om2::DetectModelScore s;
            s.name = m.name;
            s.core = static_cast<uint32_t>(m.core.size());
            for (uint64_t km : m.core) {
                auto it = slot.find(km);
                if (it != slot.end() && slotCounts[it->second].load(std::memory_order_relaxed) >= 2) ++s.coreHit;
            }
            s.score = s.core ? static_cast<double>(s.coreHit) / s.core : 0.0;
            sc.push_back(s);
        }
        std::vector<size_t> idx(sc.size());
        std::iota(idx.begin(), idx.end(), size_t(0));
        std::stable_sort(idx.begin(), idx.end(), [&](size_t a, size_t b) { return sc[a].score > sc[b].score; });
        for (size_t i : idx) {
            st.det.scores.push_back(sc[i]);
            st.det.tsmMd5.push_back(sketch.models[i].tsmMd5);
            st.det.holdList.push_back(sketch.models[i].holdList + " (md5 " + sketch.models[i].holdMd5 + ")");
        }
        st.det.ran = true;
    }
    return replaced;
}

std::string p2BlockJson(const EmitBState& st) {
    const EmitBOptions& o = st.opt;
    JsonObj p;
    p.u("enabled", o.any() ? 1 : 0);
    {
        JsonObj c;
        c.u("ends_written", st.endsRan ? 1 : 0);
        c.u("end_rows", st.endsRan ? st.ends.ends : 0);
        c.u("tip_ends", st.endsRan ? st.ends.tipEnds : 0);
        c.u("tip_bases", st.endsRan ? st.ends.tipBases : 0);
        c.u("lowercased_bases", st.ends.lowercasedBases);
        c.u("self_qa_run", st.sqa.ran ? 1 : 0);
        c.u("self_qa_alarm", st.sqa.ran && st.sqa.alarm ? 1 : 0);
        c.u("provenance_files", st.provFiles);
        c.u("detect_run", st.det.ran ? 1 : 0);
        c.u("read_pass_reads", st.rp.reads);
        p.raw("counters", c.done());
    }
    if (!o.any()) return p.done();
    {
        JsonObj f;
        f.b("ends", o.ends).b("tips_lowercase", o.lowercase).b("self_qa", o.selfQa).b("provenance", o.provenance)
            .b("detect_report", o.detect);
        p.raw("flags", f.done());
    }
    if (st.readPassRan) {
        JsonObj r;
        r.str("source", "raw input read files, every base (no trimming, no correction)");
        r.u("files", st.rp.files).u("reads", st.rp.reads).u("bases", st.rp.bases);
        if (!st.rp.error.empty()) r.str("error", st.rp.error);
        r.num("seconds", st.rp.seconds, 3);
        p.raw("read_pass", r.done());
    }
    if (o.endsAudit()) {
        if (st.endsRan) {
            p.raw("ends", endsJson(st.ends, o.lowercase));
        } else {
            JsonObj e;
            e.b("ran", false).str("error", st.endsError);
            p.raw("ends", e.done());
        }
    }
    if (o.selfQa) p.raw("self_qa", selfQaJson(st.sqa));
    if (o.provenance) p.raw("provenance", st.provenanceJson.empty() ? std::string("{}") : st.provenanceJson);
    if (o.detect) {
        JsonObj d;
        d.str("mode", "report_only");
        d.b("ran", st.det.ran);
        d.str("sketch_path", st.det.sketchPath);
        d.str("sketch_md5", st.det.sketchMd5);
        if (!st.det.error.empty()) d.str("error", st.det.error);
        if (st.det.ran) {
            std::vector<std::string> ms;
            for (size_t i = 0; i < st.det.scores.size(); ++i) {
                JsonObj m;
                m.str("name", st.det.scores[i].name).u("core", st.det.scores[i].core)
                    .u("core_hit", st.det.scores[i].coreHit).num("score", st.det.scores[i].score, 4)
                    .str("tsm_md5", st.det.tsmMd5[i]).str("hold_list", st.det.holdList[i]);
                ms.push_back(m.done());
            }
            d.raw("models", jsonArray(ms));
            d.str("best", st.det.scores.empty() ? "none" : st.det.scores[0].name);
            d.num("best_score", st.det.best(), 4);
            d.str("second", st.det.scores.size() < 2 ? "none" : st.det.scores[1].name);
            d.num("second_score", st.det.second(), 4);
            d.num("accept_threshold", om2::kDetectAccept, 2);
            d.num("mixture_threshold", om2::kDetectMixture, 2);
            d.b("accepted", st.det.accepted());
            d.str("call", st.det.accepted() ? st.det.scores[0].name : "none");
            d.b("mixture", st.det.mixture());
        }
        d.b("model_applied", false);
        p.raw("organism_detect", d.done());
    }
    if (!st.fileError.empty()) p.str("file_error", st.fileError);
    return p.done();
}

void logEmitBCounters(const EmitBState& st, std::FILE* log) {
    const EmitBOptions& o = st.opt;
    const EndsStats& e = st.ends;
    const bool er = st.endsRan;
    std::fprintf(log,
                 "[p2-ends] enabled=%d lowercase=%d records=%zu ends=%zu tip_ends=%zu tip_bases=%zu "
                 "lowercased_bases=%zu unique=%zu branching=%zu repeat_only=%zu none=%zu na=%zu\n",
                 o.endsAudit() ? 1 : 0, o.lowercase ? 1 : 0, er ? e.records : 0, er ? e.ends : 0,
                 er ? e.tipEnds : 0, er ? e.tipBases : 0, e.lowercasedBases, er ? e.unique : 0,
                 er ? e.branching : 0, er ? e.repeatOnly : 0, er ? e.none : 0, er ? e.na : 0);
    const SelfQaResult& q = st.sqa;
    std::fprintf(log,
                 "[p2-sqa] enabled=%d ran=%d completeness=%.5f qv_read0=%.2f missing_gt10x_bp=%llu "
                 "missing_3_10x_bp=%llu spike_solid=%llu alarm=%d\n",
                 o.selfQa ? 1 : 0, q.ran ? 1 : 0, q.ran ? q.completeness : 0.0, q.ran ? q.qvRead0 : 0.0,
                 static_cast<unsigned long long>(q.ran ? q.missingBpEst(kSqaBins - 1) : 0),
                 static_cast<unsigned long long>(q.ran ? q.missingBpEst(kSqaBins - 2) : 0),
                 static_cast<unsigned long long>(q.ran ? q.spikeSolid : 0), q.ran && q.alarm ? 1 : 0);
    std::fprintf(log, "[p2-prov] enabled=%d files=%zu bytes=%llu\n", o.provenance ? 1 : 0, st.provFiles,
                 static_cast<unsigned long long>(st.provBytes));
    const DetectReport& d = st.det;
    std::fprintf(log, "[p2-detect] enabled=%d ran=%d best=%s best_score=%.4f second=%s second_score=%.4f call=%s "
                      "mixture=%d model_applied=0\n",
                 o.detect ? 1 : 0, d.ran ? 1 : 0, d.ran && !d.scores.empty() ? d.scores[0].name.c_str() : "-",
                 d.ran ? std::max(0.0, d.best()) : 0.0,
                 d.ran && d.scores.size() > 1 ? d.scores[1].name.c_str() : "-", d.ran ? std::max(0.0, d.second()) : 0.0,
                 d.ran && d.accepted() ? d.scores[0].name.c_str() : "-", d.ran && d.mixture() ? 1 : 0);
    std::fprintf(log, "[p2-readpass] enabled=%d files=%zu reads=%llu bases=%llu\n", st.readPassRan ? 1 : 0,
                 st.readPassRan ? st.rp.files : 0, static_cast<unsigned long long>(st.rp.reads),
                 static_cast<unsigned long long>(st.rp.bases));
}

}  // namespace p2
}  // namespace ts
