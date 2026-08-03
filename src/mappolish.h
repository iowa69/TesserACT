// Polishing contigs against a full read alignment.
//
// The assembler's own polisher works from the k-mer spectrum, which is blind
// in exactly one place: where a base is wrong in a way that is *consistent*
// across the reads covering it, the wrong k-mers are themselves solid. A real
// aligner places each read against the finished contig, including reads the
// graph never resolved cleanly, and the resulting pileup is independent
// evidence about every position.
//
// The mapper is chosen by the caller because bowtie2 and bwa-mem do not fail
// the same way. bowtie2 is end-to-end by default and conservative about
// gapped alignments; bwa-mem soft-clips and finds indels more readily but
// scatters multi-mapping reads differently. Which is better for polishing is
// an empirical question about this data, not a matter of preference, so both
// are supported and measured.
#pragma once

#include <string>
#include <vector>

namespace ts {

enum class Mapper { None, Bowtie2, BwaMem };

const char* mapperName(Mapper m);
bool parseMapper(const std::string& s, Mapper& out);

struct MapPolishStats {
    bool ran = false;
    std::string mapper;
    size_t reads = 0;             // alignment records consumed
    size_t alignedReads = 0;
    size_t malformed = 0;         // records rejected as unparseable
    size_t positions = 0;         // contig bases covered at all
    size_t substitutions = 0;
    size_t deletions = 0;
    size_t insertions = 0;
    double meanDepth = 0;
    double seconds = 0;
    std::string error;            // non-empty if the pass could not run
};

struct MapPolishOptions {
    Mapper mapper = Mapper::None;
    std::string mapperDir;        // where to find the binaries, if not on PATH
    std::string workDir;          // scratch for the index and temporary FASTA
    std::vector<std::string> reads1, reads2;
    int threads = 1;
    // A base is rewritten only when the pileup is both deep enough to mean
    // something and lopsided enough to be unambiguous. These defaults are
    // measured, not chosen: on ERR11578413 a 0.70 threshold took mismatches
    // per 100 kbp from 0.08 to 0.36, and 0.90 to 0.32, because inside a repeat
    // the pileup is fed by several copies at once and a majority belongs to
    // whichever copy is commonest. Only at 0.99 does the pass stop doing harm.
    int minDepth = 10;
    double minFraction = 0.99;
    bool verbose = true;
};

// Rewrites `contigs` in place. Never touches a position the alignment does not
// cover, so a failure to map leaves the assembly exactly as it was.
MapPolishStats mapPolish(std::vector<std::string>& contigs, const MapPolishOptions& opt);

}  // namespace ts
