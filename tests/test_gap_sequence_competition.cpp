// Default-off, two-phase sequence-compatible gap arbitration.
#define main originalGapEvidenceTests
#include "test_gap_evidence.cpp"
#undef main
#include <algorithm>

namespace {
const char* competitionFlag = "TESSERACT_GAP_SEQUENCE_COMPETITION";
std::vector<std::string> canonicalSequences(ts::UnitigGraph g) {
    g.compact();
    std::vector<std::string> seqs;
    for (const auto& node : g.nodes)
        if (!node.deleted) seqs.push_back(std::min(node.seq, rc(node.seq)));
    std::sort(seqs.begin(), seqs.end());
    return seqs;
}
}

int main() {
    const char* prior = std::getenv(competitionFlag);
    const bool hadPrior = prior != nullptr;
    const std::string priorValue = prior ? prior : "";
    unsetenv(competitionFlag);
    if (originalGapEvidenceTests()) return 1;
    const std::string a = sequence(180);
    const std::string overlap = a.substr(160);
    const std::string b = overlap + sequence(160);
    const std::string c = sequence(180);
    auto base = graph({a, b, c});
    base.setK(31);
    const std::pair<std::string, std::string> ab{a.substr(100, 60), rc(b.substr(20, 60))};
    const std::pair<std::string, std::string> ac{ab.first, rc(c.substr(20, 60))};
    const auto input = reads({ab, ab, ac, ac});
    const auto oldNominations = nominate(base, input);
    auto oldGraph = base;
    check(oldNominations.empty() && oldGraph.closeGapsByOverlap(10, &oldNominations) == 0,
          "flag absent preserves original early competition veto");
    setenv(competitionFlag, "0", 1);
    check(nominate(base, input) == oldNominations, "flag value zero preserves v1");
    setenv(competitionFlag, "1", 1);
    const auto nominations = nominate(base, input);
    check(nominations.size() == 2, "v2 retains both mate-supported partners for sequence evaluation");
    auto good = base;
    check(good.closeGapsByOverlap(10, &nominations, nullptr, true) == 1,
          "unclosable competitor does not block unique compatible pair");
    check(good.validate().empty(), "bridge preserves bidirected invariant");
    auto expected = graph({a + b.substr(20), c});
    check(canonicalSequences(good) == canonicalSequences(expected), "bridge spells the expected biological sequence");

    const std::string other = overlap + sequence(160);
    auto competing = graph({a, b, other}); competing.setK(31);
    const std::pair<std::string, std::string> ad{ab.first, rc(other.substr(20, 60))};
    const auto ambiguous = nominate(competing, reads({ab, ab, ad, ad}));
    check(ambiguous.size() == 2, "both actual alternatives reach sequence arbitration");
    const auto before = canonicalSequences(competing);
    check(competing.closeGapsByOverlap(10, &ambiguous, nullptr, true) == 0,
          "two actual compatible partners remain unresolved");
    check(competing.validate().empty() && canonicalSequences(competing) == before,
          "ambiguous graph remains unchanged");
    auto legacy = graph({a, b, other}); legacy.setK(31);
    check(legacy.closeGapsByOverlap(10, &ambiguous, nullptr, false) == 1,
          "explicit default graph option preserves legacy greedy behavior");

    auto reordered = graph({c, b, a}); reordered.setK(31);
    const auto reorderedNominations = nominate(reordered, input);
    check(reordered.closeGapsByOverlap(10, &reorderedNominations, nullptr, true) == 1 &&
          reordered.validate().empty() && canonicalSequences(reordered) == canonicalSequences(good),
          "node order preserves selected biological closure");
    auto flipped = graph({rc(a), b, rc(c)}); flipped.setK(31);
    const auto flippedNominations = nominate(flipped, input);
    check(flipped.closeGapsByOverlap(10, &flippedNominations, nullptr, true) == 1 &&
          flipped.validate().empty() && canonicalSequences(flipped) == canonicalSequences(good),
          "unitig reverse-complement representation preserves closure");

    const std::string motif = "ACGTTGCAAC";
    const std::string periodicA = sequence(160) + motif + motif;
    const std::string periodicB = motif + motif + sequence(160);
    auto multipleLengths = graph({periodicA, periodicB}); multipleLengths.setK(31);
    std::set<std::pair<uint64_t, uint64_t>> onePartner{{1, 2}};
    check(multipleLengths.closeGapsByOverlap(10, &onePartner, nullptr, true) == 1,
          "two overlap lengths and reverse traversal count as one partner");
    check(multipleLengths.nodes.back().seq.size() == 40 && multipleLengths.validate().empty(),
          "same partner keeps largest valid overlap");

    const std::string d = sequence(180), e = d.substr(160) + sequence(160);
    auto disjoint = graph({a, b, d, e}); disjoint.setK(31);
    std::set<std::pair<uint64_t, uint64_t>> twoPairs{{1, 2}, {5, 6}};
    check(disjoint.closeGapsByOverlap(10, &twoPairs, nullptr, true) == 2 && disjoint.validate().empty(),
          "all disjoint mutually unique pairs can be inserted");
    if (hadPrior) setenv(competitionFlag, priorValue.c_str(), 1); else unsetenv(competitionFlag);
    std::printf("gap sequence competition: %d checks, %d failures\n", checks, failures);
    return failures ? 1 : 0;
}
