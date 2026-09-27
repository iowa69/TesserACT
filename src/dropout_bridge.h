// Bridging facing dead ends across read-coverage dropouts (TESSERACT_DROPOUT_BRIDGE=1,
// default OFF; combo2 PKG-BRIDGE).
//
// A de Bruijn graph breaks wherever the reads are too thin for the k-mers to pass the
// solidity cutoff, even though individual reads still run straight across the thin
// stretch. Such a dropout leaves two dead ends facing each other with nothing in the
// graph to join them. This stage looks for reads that leave one dead end's terminal
// window and, further along the same read, enter another dead end's terminal window,
// and bridges the pair when enough independent fragments agree on the bases between.
//
// It runs once, on the final-rung simplified graph, before the paired-end resolver.
//
// Reads are examined THROUGH the error corrector's mask (raw bases), unlike the design's
// first draft: correctReads masks from wherever correction stalls to the read's end, so a
// read spanning a dropout almost always has the far flank masked. On the 2x301 sentinel a
// masked scan found 3 nominations from 25,804 anchor-hitting reads; the raw scan finds 188.
// Every inserted base is either a flank base already in the graph or a consensus base
// carried by at least `minReads` distinct fragments that hold unique anchors of both
// flanks. The bridge unitig follows closeGapsByOverlap's insertBridge convention.
#pragma once

#include <cstddef>
#include <cstdint>
#include <string>

#include "graph.h"
#include "seqio.h"

namespace ts {

struct DropoutBridgeStats {
    int k = 0;
    uint32_t minReads = 0;
    double medianCoverage = 0;
    size_t deadEnds = 0;          // live unitig ends with no link
    size_t eligible = 0;          // dead ends passing length/coverage AND keeping >= 1 anchor
    size_t endsNoAnchor = 0;      // passed length/coverage but lost every anchor
    size_t anchors = 0;           // distinct window A-mers kept (unique in windows and graph)
    size_t anchorsDropped = 0;    // distinct window A-mers dropped (repeated)
    size_t readsHit = 0;          // reads with at least one anchor hit
    size_t maskedHits = 0;        // anchor hits whose 31 read bases include an EC-masked base
    size_t nominations = 0;       // read-level nominations (before any dedup)
    size_t pairsNominated = 0;    // distinct unordered end pairs with >= 1 nomination
    size_t refMinReads = 0;       // fewer than minReads distinct pair ids
    size_t refDup = 0;            // fewer than minReads distinct fragments after start dedup
    size_t refModal = 0;          // modal gap below 2/3 of fragments or below minReads
    size_t refConsensus = 0;      // some gap position below 2/3 agreement
    size_t refCompeting = 0;      // an end has another partner with >= 2 distinct pairs
    size_t refOverlap = 0;        // negative gap whose overlap disagrees between the flanks
    size_t refLowComplexity = 0;  // windows + consensus low-complexity
    size_t refSelf = 0;           // both ends on one unitig (circularising; not done in round 1)
    size_t bridged = 0;
    size_t overlapBridged = 0;    // bridged with gap < 0
    size_t bpAdded = 0;           // consensus bases inserted (sum of max(gap, 0))
    size_t rawGapVotes = 0;       // consensus votes of accepted bridges read through the EC mask
    size_t merged = 0;            // merges done by the post-bridge compact()
    double rssDeltaMb = 0;
    double seconds = 0;
    std::string error;            // non-empty only if the graph failed validate() afterwards
    bool enabled = false;         // build_v3: false for the idle line a run without the flag prints
};

// Bridges dropouts in `g` in place. `threads` workers scan the reads. `minReads` is the
// evidence floor (distinct pairs and distinct fragment starts). With `trace`, prints one
// "[dropoutbridge-join]" line per bridge and one "[dropoutbridge-refused]" line per end
// pair that reached minReads distinct pairs but was refused.
DropoutBridgeStats bridgeDropouts(UnitigGraph& g, const SequenceStore& reads, int threads,
                                  uint32_t minReads, bool trace);

// The single "[dropoutbridge] ..." counter line (no trailing newline). build_v3: printed on
// every run (OBJECTIVE A2); a run without TESSERACT_DROPOUT_BRIDGE=1 prints the zero line with
// the trailing field enabled=0.
std::string formatDropoutBridgeStats(const DropoutBridgeStats& s);

}  // namespace ts
