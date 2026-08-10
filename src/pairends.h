// Paired-end evidence at contig ends.
//
// The paired resolver uses read pairs to walk the unitig graph, and then its evidence is
// thrown away: `PairedResolver` goes out of scope before the model join and the layout
// stage ever run. Those two stages decide what joins to which, at contig ends, with no
// paired evidence at all -- on a prior learned from other genomes and on sequence overlap.
//
// This rebuilds the evidence in the frame where those decisions are actually made. Reads
// are anchored onto the FINAL contigs rather than the unitigs, which avoids the lossy
// unitig->contig mapping (a contig built by joining has no single graph walk, and the
// mapping is cleared before layout anyway) and lands every distance directly in contig
// coordinates.
//
// What it can and cannot say, stated plainly because it bounds every use below: a pair
// only speaks about a junction it can span. On this cohort fragments run about 217 +/- 110
// bp while the repeats that break the assembly are 500 bp and up, so for MOST junctions the
// model stage exists precisely because no pair reaches across. Paired evidence is decisive
// on the short-gap and overlap tail, silent elsewhere, and silence must never be read as a
// veto.
#pragma once

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

#include "resolve.h"     // InsertModel
#include "seqio.h"

namespace ts {

// A contig end. Port 2c+1 is contig c's right end, 2c+0 its left -- the same convention
// the model join uses, so the two stages speak the same language.
inline uint32_t endPort(uint32_t contig, int end) { return (contig << 1) | (end & 1); }

struct PairEvidence {
    size_t total = 0;        // pairs linking these two ends at all
    size_t consistent = 0;   // ... whose implied gap agrees with the one proposed
    int32_t medianGap = 0;   // gap the pairs themselves imply
    bool spannable = false;  // could a fragment have reached across? if not, silence
                             // means nothing and the caller must not treat it as a veto
};

class ContigEndLinks {
public:
    // Anchors every read onto `contigs` and accumulates links between contig ends.
    // `insert` supplies the fragment distribution; without a usable one nothing is built
    // and every query returns `spannable = false`.
    ContigEndLinks(const std::vector<std::string>& contigs, const SequenceStore& reads,
                   const InsertModel& insert, int threads, int anchorK = 31);

    bool usable() const { return usable_; }

    // Evidence that `exitPort` continues into `entryPort` across a gap of `predictedGap`.
    PairEvidence query(uint32_t exitPort, uint32_t entryPort, int32_t predictedGap) const;

    // Evidence that a contig's two ends meet -- that it is a circle. For a plasmid this is
    // the difference between a contig and a finished replicon, so it is worth reporting on
    // its own. Note what it does NOT mean: a circular plasmid carrying an unresolved repeat
    // also shows head-to-tail pairs, so this says "the ends join", not "the molecule is
    // complete".
    PairEvidence circular(uint32_t contig) const;

    // Pairs linking two contigs anywhere, not just at their ends. Contigs on one plasmid
    // link to each other; contigs on different molecules should not.
    size_t linkWeight(uint32_t contigA, uint32_t contigB) const;

    // Support floor scaled to depth, so a 20x library and a 200x library are not held to
    // the same absolute count. Same shape as the scaffolder's floor.
    size_t supportFloor() const;

private:
    struct Link { std::vector<int32_t> spans; };

    bool usable_ = false;
    InsertModel insert_;
    double meanDepth_ = 0;
    std::vector<size_t> contigLen_;
    std::unordered_map<uint64_t, Link> ends_;        // (portA<<32)|portB -> spans
    std::unordered_map<uint64_t, size_t> contigLinks_;  // (a<<32)|b -> pair count
};

}  // namespace ts
