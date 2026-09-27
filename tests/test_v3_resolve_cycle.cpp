// T13 (build_v3 G-resolve): accepted scaffold joins form simple paths AND simple cycles; the
// emit walk starts only at a free port, so a cycle (a circular element split at >= 2 gaps)
// was rendered chain by chain and every one of its accepted joins silently dropped.
//
//  L  : linear control, two unitigs, one 200 bp hole spanned by pairs -> 1 scaffold, 1 join.
//  C2 : circular 12 kb element as two unitigs with two 200 bp holes -> a 2-cycle.
//  C3 : circular 15 kb element as three unitigs -> a 3-cycle.
//  M  : C2 plus an unrelated linear pair in the same graph.
// release: C2 -> 2 contigs / 0 joins, C3 -> 3 / 0 ([scafcycle] still counts the cycle);
// TESSERACT_FIX_SCAFFOLD_CYCLE=1 (or TESSERACT_FIXES=1): each cycle broken once at its weakest
// join -> C2 1 scaffold / 1 join, C3 1 / 2; every emitted non-N run is genome sequence.
#include "test_v3_resolve_util.h"
#include "test_env.h"   // build_v3: T40, clear the ambient TESSERACT_* first

using namespace v3r;

namespace {
void tile(std::vector<Pair>& out, const std::string& g, bool circular, size_t step = 4, size_t len = 100,
          size_t frag = 320) {
    const std::string gg = circular ? g + g.substr(0, frag) : g;
    const size_t limit = circular ? g.size() : g.size() - frag + 1;
    for (size_t p = 0; p < limit; p += step) out.push_back({gg.substr(p, len), rc(gg.substr(p + frag - len, len))});
}
struct Out {
    std::vector<std::string> seqs;
    ts::ResolveStats st;
    size_t long1k = 0, nRuns = 0;
    std::string log;
};
Out run(const ts::UnitigGraph& g, const ts::SequenceStore& reads, const char* fix, const char* umbrella) {
    releaseDefaults();
    setFlag("TESSERACT_FIX_SCAFFOLD_CYCLE", fix);
    setFlag("TESSERACT_FIXES", umbrella ? umbrella : "0");   // 1.4.0: unset would follow the default umbrella (on)
    StderrCapture cap;
    cap.start();
    ts::PairedResolver r(g, reads, 1, 2, 2.0, 0.10);
    r.setScaffolding(true);
    r.buildSupport();
    Out o;
    std::vector<double> covs;
    r.resolve(o.seqs, covs);
    o.log = cap.stop();
    o.st = r.stats();
    for (const auto& s : o.seqs) {
        if (s.size() >= 1000) ++o.long1k;
        for (size_t i = 0; i < s.size(); ++i)
            if (s[i] == 'N' && (i == 0 || s[i - 1] != 'N')) ++o.nRuns;
    }
    const auto sc = lines(o.log, "[scafcycle]");
    std::printf("    seqs=%zu (>=1kb %zu) nRuns=%zu joins=%zu | %s\n", o.seqs.size(), o.long1k, o.nRuns,
                o.st.scaffoldJoins, sc.empty() ? "(no [scafcycle])" : sc[0].c_str());
    return o;
}
void addNodes(ts::UnitigGraph& g, const std::string& genome, const std::vector<size_t>& offsets, size_t len) {
    for (size_t o : offsets) {
        ts::Unitig nd;
        nd.seq = genome.substr(o, len);
        nd.coverage = 40;
        g.nodes.push_back(nd);
    }
}
bool piecesInGenome(const Out& o, const std::string& genome, bool circular) {
    const std::string gg = circular ? genome + genome : genome;
    const std::string rg = rc(gg);
    for (const auto& s : o.seqs) {
        size_t i = 0;
        while (i < s.size()) {
            if (s[i] == 'N') { ++i; continue; }
            size_t j = i;
            while (j < s.size() && s[j] != 'N') ++j;
            const std::string piece = s.substr(i, j - i);
            if (gg.find(piece) == std::string::npos && rg.find(piece) == std::string::npos) return false;
            i = j;
        }
    }
    return true;
}
long long counter(const Out& o, const char* key) {
    const auto sc = lines(o.log, "[scafcycle]");
    return sc.size() == 1 ? field(sc[0], key) : -99;
}
}  // namespace

int main() {
    testenv::clearTesseractEnv();
    const int k = 31;
    {
        std::printf("[L] linear control\n");
        const std::string G = dna(12000, 101);
        ts::UnitigGraph g; g.setK(k);
        addNodes(g, G, {0, 6000}, 5800);
        std::vector<Pair> p; tile(p, G, false);
        const auto reads = store(p, "L");
        for (const char* fx : {static_cast<const char*>(nullptr), "1"}) {
            const Out o = run(g, reads, fx, nullptr);
            expect(o.long1k == 1 && o.st.scaffoldJoins == 1 && o.nRuns == 1, "L: one scaffold, 1 join, 1 N-run");
            expect(counter(o, "cycles") == 0, "L: [scafcycle] cycles=0");
        }
    }
    {
        std::printf("[C2] circular 12 kb, two unitigs, two holes\n");
        const std::string C = dna(12000, 202);
        ts::UnitigGraph g; g.setK(k);
        addNodes(g, C, {0, 6000}, 5800);
        std::vector<Pair> p; tile(p, C, true);
        const auto reads = store(p, "C2");
        const Out rel = run(g, reads, nullptr, nullptr);
        expect(rel.long1k == 2 && rel.st.scaffoldJoins == 0 && rel.nRuns == 0,
               "C2 release: two unjoined contigs, 0 joins (both accepted joins dropped: the defect)");
        expect(counter(rel, "enabled") == 0 && counter(rel, "cycles") == 1 && counter(rel, "joinsInCycles") == 2 &&
               counter(rel, "joinsDropped") == 2, "C2 release: [scafcycle] enabled=0 cycles=1 joinsInCycles=2 joinsDropped=2");
        for (int umbrella = 0; umbrella < 2; ++umbrella) {
            const Out fix = umbrella ? run(g, reads, nullptr, "1") : run(g, reads, "1", nullptr);
            expect(fix.long1k == 1 && fix.st.scaffoldJoins == 1 && fix.nRuns == 1,
                   std::string("C2 ") + (umbrella ? "TESSERACT_FIXES=1" : "FIX=1") + ": one scaffold, 1 join, 1 N-run");
            expect(counter(fix, "enabled") == 1 && counter(fix, "broken") == 1 && counter(fix, "joinsDropped") == 1,
                   "C2 fixed: [scafcycle] enabled=1 broken=1 joinsDropped=1");
            expect(piecesInGenome(fix, C, true), "C2 fixed: every emitted non-N run is genome sequence");
        }
        const Out optOut = run(g, reads, "0", "1");
        expect(optOut.long1k == 2 && optOut.st.scaffoldJoins == 0, "C2 TESSERACT_FIXES=1 with FIX=0: release behaviour");
    }
    {
        std::printf("[C3] circular 15 kb, three unitigs, three holes\n");
        const std::string C = dna(15000, 303);
        ts::UnitigGraph g; g.setK(k);
        addNodes(g, C, {0, 5000, 10000}, 4800);
        std::vector<Pair> p; tile(p, C, true);
        const auto reads = store(p, "C3");
        const Out rel = run(g, reads, nullptr, nullptr);
        expect(rel.long1k == 3 && rel.st.scaffoldJoins == 0, "C3 release: three unjoined contigs, 0 joins (the defect)");
        const Out fix = run(g, reads, "1", nullptr);
        expect(fix.long1k == 1 && fix.st.scaffoldJoins == 2 && fix.nRuns == 2, "C3 fixed: one scaffold, 2 joins, 2 N-runs");
        expect(piecesInGenome(fix, C, true), "C3 fixed: every emitted non-N run is genome sequence");
    }
    {
        std::printf("[M] 2-cycle plus an unrelated linear pair\n");
        const std::string C = dna(12000, 404), G = dna(12000, 505);
        ts::UnitigGraph g; g.setK(k);
        addNodes(g, C, {0, 6000}, 5800);
        addNodes(g, G, {0, 6000}, 5800);
        std::vector<Pair> p; tile(p, C, true); tile(p, G, false);
        const auto reads = store(p, "M");
        const Out rel = run(g, reads, nullptr, nullptr);
        expect(rel.long1k == 3 && rel.st.scaffoldJoins == 1, "M release: linear pair joined, cycle chains alone");
        const Out fix = run(g, reads, "1", nullptr);
        expect(fix.long1k == 2 && fix.st.scaffoldJoins == 2, "M fixed: linear pair joined and the cycle broken once");
    }
    releaseDefaults();
    return finish("test_v3_resolve_cycle (T13)");
}
