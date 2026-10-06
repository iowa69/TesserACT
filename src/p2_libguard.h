// Phase 2, W1 item R2 (EVAL_PLAN_P2): the library orientation and insert-model guard.
//
// TESSERACT_P2_LIBGUARD=1 (default off). The paired resolver learns its insert model from
// pairs whose two mates anchor on one unitig, and assumes an inward-facing (FR) library
// everywhere else: in the scaffolder, in the repeat-resolution support table and in the
// circularity call (pairends.cpp). A mate-pair (RF, outward) library or a library whose mates
// do not belong together therefore turns into confident-looking joins and circles. The
// backward audit measured 54 of the 59 false F2 scaffold joins and 21 false circles from one
// mate-pair library (99.3 % outward pairs, FR model fitted on 0.27 % of pairs).
//
// The guard measures, on the same-unitig pairs the resolver already anchors:
//   inward   the forward-strand mate starts at or before the reverse-strand mate (FR);
//   outward  the forward-strand mate starts at or after the END of the reverse-strand mate
//            (the two reads face away from each other and do not overlap: RF);
//   dovetail the rest of the opposite-strand pairs (overlapping and offset: adapter
//            read-through in an FR library looks exactly like this, so it counts as neither);
//   same     both mates on one strand.
// and fires when
//   outward / (outward + inward) > 0.50 (evaluated with >= 1000 such pairs), or
//   frFit = insert-model observations / read pairs < 0.25 (the 508-input census definition:
//   report.json insert_size.observations over input pairs; median library 0.95).
// When it fires the run makes no pair-based joins and no pair-based circle calls: the
// resolver's pair support table is emptied (repeats are crossed only where the graph or a
// single read thread says so), scaffolding is off, and the final-contig pair links
// (ContigEndLinks: circularity and pair-based plasmid grouping) are not built. On a library
// where it does not fire nothing changes, byte for byte.
#pragma once

#include <cstddef>

namespace ts {

struct LibGuardStats {
    bool enabled = false;       // TESSERACT_P2_LIBGUARD
    bool evaluated = false;     // the resolver ran on a paired library
    bool fired = false;
    size_t pairs = 0;           // read pairs in the library
    size_t sameUnitig = 0;      // pairs with both mates anchored on one unitig
    size_t inward = 0;
    size_t outward = 0;
    size_t dovetail = 0;
    size_t sameStrand = 0;
    size_t frObservations = 0;  // pairs the FR insert model was fitted on
    double outwardFrac = 0;     // outward / (outward + inward); 0 when not evaluable
    double frFit = 0;           // frObservations / pairs
    bool outwardEvaluable = false;
    int reason = 0;             // bit 0: outward > 0.50; bit 1: FR fit < 0.25
};

namespace libguard {

constexpr double kMaxOutwardFrac = 0.50;   // EVAL_PLAN_P2 / FINAL R2: "> 50 % of pairs point outward"
constexpr double kMinFrFit = 0.25;         // "the FR model fits < 25 % of pairs"
constexpr size_t kMinOrientedPairs = 1000; // the insert fit's own minimum sample

// Classifies one same-unitig pair from its two anchors. `fwdPos`/`fwdLen` describe the mate
// that runs along the unitig's forward strand, `revPos`/`revLen` the other one. Returns
// 0 inward, 1 outward, 2 dovetail.
inline int classifyOpposite(long fwdPos, long fwdLen, long revPos, long revLen) {
    (void)fwdLen;
    if (fwdPos <= revPos) return 0;
    if (fwdPos >= revPos + revLen) return 1;
    return 2;
}

// Fills outwardFrac, frFit, outwardEvaluable, reason and fired from the counts.
inline void decide(LibGuardStats& s) {
    s.evaluated = true;
    const size_t oriented = s.inward + s.outward;
    s.outwardEvaluable = oriented >= kMinOrientedPairs;
    s.outwardFrac = oriented ? static_cast<double>(s.outward) / static_cast<double>(oriented) : 0.0;
    s.frFit = s.pairs ? static_cast<double>(s.frObservations) / static_cast<double>(s.pairs) : 0.0;
    s.reason = 0;
    if (s.outwardEvaluable && s.outwardFrac > kMaxOutwardFrac) s.reason |= 1;
    if (s.pairs > 0 && s.frFit < kMinFrFit) s.reason |= 2;
    s.fired = s.enabled && s.reason != 0;
}

inline const char* verdict(const LibGuardStats& s) {
    if (!s.enabled) return "off";
    if (!s.evaluated) return "not_evaluated";
    if (!s.fired) return "pass";
    if (s.reason == 3) return "fired_outward_and_fr_fit";
    return (s.reason & 1) ? "fired_outward" : "fired_fr_fit";
}

}  // namespace libguard
}  // namespace ts
