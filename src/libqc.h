#pragma once

#include <string>

namespace ts {

// What a QC pass over the raw reads can tell the assembler before it counts a single
// k-mer. Every field here is a MEASUREMENT of the library, not a parameter of the run:
// nothing in this struct should ever be set from a command line except the path it was
// read from.
//
// The assembler currently reconstructs versions of several of these from its own k-mer
// histograms, which is a harder problem than it looks -- the histogram it reads is
// contaminated by the very error k-mers it is trying to measure around, and at the top of
// the k ladder the real coverage mode can be smaller than the error shoulder. A QC pass
// that has already separated the error mass from the genomic mass at a small, safe k can
// hand over the answer directly.
struct LibraryQC {
    bool loaded = false;
    std::string source;

    // --- k-mer spectrum -----------------------------------------------------------
    bool   spectrumReliable = false;   // the QC tool's own verdict; false means ignore
    int    spectrumK = 0;              // k the spectrum was built at (21 in practice)
    double peakDepth = 0;              // k-mer coverage mode at spectrumK
    double genomeSize = 0;             // bp, whole library including plasmids
    double genomeSizeCore = 0;         // bp, excluding the repeat/multi-copy mass
    double readDepth = 0;              // x, READ depth -- not k-mer depth
    double errorKmerFraction = 0;      // share of DISTINCT k-mers that are errors
    double repeatFraction = 0;

    // --- error model --------------------------------------------------------------
    // Measured by comparing overlapping mates against each other, so it is an observed
    // substitution rate rather than the instrument's own opinion of its Q scores.
    bool   errorModelHasData = false;
    double errorRate = 0;

    // --- library geometry ---------------------------------------------------------
    // Derived from the QC tool's insert-size histogram, which is built by overlapping
    // mates on the RAW reads. The assembler builds its own from reads anchored onto the
    // unitig graph, which is circular in a way that matters: the pairs it can place are
    // the ones the current graph already explains, so the model it fits is biased toward
    // the graph it is about to be used to improve.
    bool   insertUsable = false;
    double insertPeak = 0;             // bp, modal fragment length
    double insertMean = 0;             // bp, mean over the trimmed histogram
    double insertSd = 0;
    int    insertP1 = 0;               // 1st percentile, bp
    int    insertP99 = 0;              // 99th percentile, bp
    size_t insertObservations = 0;
    double meanReadLength = 0;         // bases / reads, over the whole input

    // Expected k-mer coverage at a given k, from read depth and read length:
    //
    //     cov(k) = depth * (L - k + 1) / L
    //
    // This is the quantity that actually governs whether a rung of the ladder is
    // assemblable, and it is knowable in advance. The assembler's own estimate of it is
    // only available AFTER counting that rung -- too late to decide whether to count it.
    // Returns 0 when the QC data cannot support the calculation.
    double expectedKmerCoverage(int k) const;
};

// Parse a fastplus/scepter JSON report. Returns false and sets `error` when the file is
// missing, malformed, or carries no k-mer spectrum -- callers must not proceed with a
// half-populated struct, since a zeroed field is indistinguishable from a real measurement
// of zero and would silently disable exactly the decisions this exists to inform.
bool loadLibraryQC(const std::string& path, LibraryQC& out, std::string& error);

}  // namespace ts
