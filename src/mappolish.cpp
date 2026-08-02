#include "mappolish.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <unordered_map>

#include "util.h"

namespace ts {

const char* mapperName(Mapper m) {
    switch (m) {
        case Mapper::Bowtie2: return "bowtie2";
        case Mapper::BwaMem:  return "bwa";
        case Mapper::None:    return "none";
    }
    return "none";
}

bool parseMapper(const std::string& s, Mapper& out) {
    if (s == "none" || s.empty()) { out = Mapper::None; return true; }
    if (s == "bowtie2") { out = Mapper::Bowtie2; return true; }
    if (s == "bwa" || s == "bwa-mem" || s == "bwamem") { out = Mapper::BwaMem; return true; }
    return false;
}

namespace {

// Per-position evidence. Counts are capped well below overflow at any depth a
// bacterial library reaches.
struct Column {
    uint32_t base[4] = {0, 0, 0, 0};
    uint32_t deletion = 0;
    uint32_t depth = 0;
};

inline int codeOf(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return -1;
    }
}

std::string shellQuote(const std::string& s) {
    std::string out = "'";
    for (char c : s) {
        if (c == '\'') out += "'\\''";
        else out.push_back(c);
    }
    out.push_back('\'');
    return out;
}

std::string joinList(const std::vector<std::string>& v) {
    std::string out;
    for (size_t i = 0; i < v.size(); ++i) {
        if (i) out.push_back(',');
        out += shellQuote(v[i]);
    }
    return out;
}

std::string binary(const std::string& dir, const char* name) {
    if (dir.empty()) return name;
    std::string p = dir;
    if (!p.empty() && p.back() != '/') p.push_back('/');
    return shellQuote(p + name);
}

bool writeFasta(const std::string& path, const std::vector<std::string>& contigs) {
    std::FILE* f = std::fopen(path.c_str(), "w");
    if (!f) return false;
    for (size_t i = 0; i < contigs.size(); ++i) {
        std::fprintf(f, ">c%zu\n", i);
        const std::string& s = contigs[i];
        for (size_t p = 0; p < s.size(); p += 80) {
            std::fwrite(s.data() + p, 1, std::min<size_t>(80, s.size() - p), f);
            std::fputc('\n', f);
        }
    }
    std::fclose(f);
    return true;
}

}  // namespace

MapPolishStats mapPolish(std::vector<std::string>& contigs, const MapPolishOptions& opt) {
    MapPolishStats st;
    st.mapper = mapperName(opt.mapper);
    if (opt.mapper == Mapper::None || contigs.empty() || opt.reads1.empty()) return st;

    util::Timer timer;
    // The same lesson the k-mer polisher learned, for the same reason: inside a
    // repeat the pileup is fed by several copies at once, and a simple majority
    // belongs to whichever copy is commonest rather than to the locus being
    // polished. Only near-unanimity distinguishes a real error from a copy that
    // outvoted its neighbour. Tunable for experiments.
    const double minFraction = std::getenv("TESSERA_MAPPOLISH_FRACTION")
                                   ? std::atof(std::getenv("TESSERA_MAPPOLISH_FRACTION"))
                                   : opt.minFraction;
    const int minDepth = std::getenv("TESSERA_MAPPOLISH_DEPTH")
                             ? std::atoi(std::getenv("TESSERA_MAPPOLISH_DEPTH"))
                             : opt.minDepth;
    const std::string ref = opt.workDir + "/mappolish_ref.fasta";
    if (!writeFasta(ref, contigs)) {
        st.error = "cannot write " + ref;
        return st;
    }

    // Build the index and stream SAM back through a pipe: the alignment is
    // consumed a record at a time and never lands on disk.
    std::string cmd;
    const std::string qref = shellQuote(ref);
    const std::string t = std::to_string(opt.threads > 0 ? opt.threads : 1);
    if (opt.mapper == Mapper::Bowtie2) {
        cmd = binary(opt.mapperDir, "bowtie2-build") + " --threads " + t + " -q " + qref + " " +
              qref + " >/dev/null 2>&1 && " + binary(opt.mapperDir, "bowtie2") +
              " -p " + t + " -x " + qref + " -1 " + joinList(opt.reads1);
        if (!opt.reads2.empty()) cmd += " -2 " + joinList(opt.reads2);
        // Local alignment keeps a read whose tail runs off a contig end, which
        // is exactly where the assembly is least certain.
        cmd += " --local --no-unal 2>/dev/null";
    } else {
        cmd = binary(opt.mapperDir, "bwa") + " index " + qref + " >/dev/null 2>&1 && " +
              binary(opt.mapperDir, "bwa") + " mem -t " + t + " " + qref + " " +
              joinList(opt.reads1);
        if (!opt.reads2.empty()) cmd += " " + joinList(opt.reads2);
        cmd += " 2>/dev/null";
    }
    // bwa takes its mates as two file arguments rather than a comma list.
    if (opt.mapper == Mapper::BwaMem && opt.reads1.size() > 1) {
        st.error = "bwa polishing takes a single read pair";
        return st;
    }

    std::FILE* pipe = popen(cmd.c_str(), "r");
    if (!pipe) {
        st.error = std::string("cannot run ") + mapperName(opt.mapper);
        return st;
    }

    std::vector<std::vector<Column>> cols(contigs.size());
    for (size_t i = 0; i < contigs.size(); ++i) cols[i].resize(contigs[i].size());
    // Insertions are sparse, so they are kept per (contig, position, sequence)
    // rather than costing a slot at every base.
    std::unordered_map<uint64_t, std::unordered_map<std::string, uint32_t>> insertions;

    std::string line;
    line.reserve(1 << 16);
    char buf[1 << 16];
    while (std::fgets(buf, sizeof(buf), pipe)) {
        if (buf[0] == '@') continue;
        ++st.reads;

        // qname flag rname pos mapq cigar rnext pnext tlen seq qual
        char* p = buf;
        auto field = [&]() -> char* {
            char* start = p;
            while (*p && *p != '\t') ++p;
            if (*p == '\t') *p++ = '\0';
            return start;
        };
        field();                              // qname
        const int flag = std::atoi(field());
        const char* rname = field();
        const long pos = std::atol(field());  // 1-based
        const int mapq = std::atoi(field());
        const char* cigar = field();
        field(); field(); field();            // rnext pnext tlen
        const char* seq = field();

        if (flag & 0x4) continue;             // unmapped
        if (flag & 0x100 || flag & 0x800) continue;   // secondary or supplementary
        // A read that could sit in several places says nothing about which
        // base belongs here.
        if (mapq < 10) continue;
        if (rname[0] != 'c') continue;
        const size_t ci = static_cast<size_t>(std::atoll(rname + 1));
        if (ci >= contigs.size() || pos < 1) continue;
        ++st.alignedReads;

        size_t rpos = static_cast<size_t>(pos - 1);   // reference cursor
        size_t qpos = 0;                              // query cursor
        std::vector<Column>& col = cols[ci];

        for (const char* c = cigar; *c;) {
            long len = 0;
            while (*c >= '0' && *c <= '9') len = len * 10 + (*c++ - '0');
            const char op = *c++;
            switch (op) {
                case 'M': case '=': case 'X':
                    for (long i = 0; i < len && rpos < col.size(); ++i, ++rpos, ++qpos) {
                        const int b = codeOf(seq[qpos]);
                        if (b >= 0) ++col[rpos].base[b];
                        ++col[rpos].depth;
                    }
                    break;
                case 'D': case 'N':
                    for (long i = 0; i < len && rpos < col.size(); ++i, ++rpos) {
                        ++col[rpos].deletion;
                        ++col[rpos].depth;
                    }
                    break;
                case 'I': {
                    if (rpos > 0 && rpos <= col.size()) {
                        const uint64_t key = (static_cast<uint64_t>(ci) << 32) |
                                             static_cast<uint64_t>(rpos);
                        insertions[key][std::string(seq + qpos, static_cast<size_t>(len))] += 1;
                    }
                    qpos += static_cast<size_t>(len);
                    break;
                }
                case 'S': qpos += static_cast<size_t>(len); break;
                case 'H': case 'P': break;
                default: break;
            }
        }
    }
    const int rc = pclose(pipe);
    if (st.alignedReads == 0) {
        st.error = rc == 0 ? "no reads aligned" : "mapper failed";
        st.seconds = timer.elapsed();
        return st;
    }

    // ---- apply the pileup --------------------------------------------------
    // Rebuilt rather than edited in place, because insertions and deletions
    // change the coordinates of everything after them.
    unsigned long long depthSum = 0;
    for (size_t ci = 0; ci < contigs.size(); ++ci) {
        const std::vector<Column>& col = cols[ci];
        std::string out;
        out.reserve(contigs[ci].size() + 64);

        for (size_t i = 0; i < contigs[ci].size(); ++i) {
            const Column& c = col[i];
            if (c.depth) { ++st.positions; depthSum += c.depth; }

            char keep = contigs[ci][i];
            if (static_cast<int>(c.depth) >= minDepth) {
                // Ns are what the gap filler could not close; the aligner may
                // still have placed reads across them, and any call there is
                // an improvement on nothing.
                int bestB = -1;
                uint32_t bestN = 0;
                for (int b = 0; b < 4; ++b) {
                    if (c.base[b] > bestN) { bestN = c.base[b]; bestB = b; }
                }
                const double frac = static_cast<double>(bestN) / static_cast<double>(c.depth);
                const double delFrac = static_cast<double>(c.deletion) / static_cast<double>(c.depth);

                if (delFrac >= minFraction) {
                    ++st.deletions;
                    continue;                      // the column is dropped
                }
                if (bestB >= 0 && frac >= minFraction) {
                    const char call = "ACGT"[bestB];
                    if (call != keep) {
                        ++st.substitutions;
                        keep = call;
                    }
                }
            }
            out.push_back(keep);

            const uint64_t key = (static_cast<uint64_t>(ci) << 32) | static_cast<uint64_t>(i + 1);
            auto it = insertions.find(key);
            if (it != insertions.end() && static_cast<int>(c.depth) >= minDepth) {
                const std::string* bestSeq = nullptr;
                uint32_t bestCount = 0;
                for (const auto& kv : it->second) {
                    if (kv.second > bestCount) { bestCount = kv.second; bestSeq = &kv.first; }
                }
                if (bestSeq &&
                    static_cast<double>(bestCount) >= minFraction * static_cast<double>(c.depth)) {
                    out += *bestSeq;
                    ++st.insertions;
                }
            }
        }
        contigs[ci].swap(out);
    }

    st.meanDepth = st.positions ? static_cast<double>(depthSum) / static_cast<double>(st.positions) : 0;
    st.seconds = timer.elapsed();
    st.ran = true;

    if (opt.verbose) {
        std::fprintf(stderr,
                     "      %s: %s reads aligned, mean depth %.1fx, "
                     "%zu substitutions %zu insertions %zu deletions, %.1fs\n",
                     mapperName(opt.mapper),
                     util::commify(static_cast<long long>(st.alignedReads)).c_str(),
                     st.meanDepth, st.substitutions, st.insertions, st.deletions, st.seconds);
    }
    return st;
}

}  // namespace ts
