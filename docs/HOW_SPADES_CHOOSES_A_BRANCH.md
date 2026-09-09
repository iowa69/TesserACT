# What SPAdes actually does at a branch, read from its source

TesserACT loses one metric to SPAdes reproducibly, on two independent panels:

| panel | arm | misassemblies, TesserACT vs SPAdes | p |
|---|---|---|---|
| 666 *K. pneumoniae*, no model | vanilla | 151 / 218 (win/loss) | 0.002 |
| 151 *S. aureus*, non-clonal | `--organism` | 26 / 42 (win/loss) | 0.048 |

Same direction, different organism, different model state. That reproducibility is the
reason to look for one algorithmic cause rather than tuning a threshold.

## The shape of the error is not distinctive

Classifying every extensive misassembly in the *S. aureus* cohort by the size of the
relocation QUAST reports:

| relocation size | TesserACT | SPAdes |
|---|---|---|
| 1–5 kb | 80 (48.2%) | 52 (46.0%) |
| 5–20 kb | 13 (7.8%) | 13 (11.5%) |
| 20–100 kb | 15 (9.0%) | 11 (9.7%) |
| >=100 kb | 58 (34.9%) | 37 (32.7%) |
| **total** | **166** | **113** |

The distributions are the same shape. We are not making a *different kind* of mistake
that a special case could catch — we make the same kinds about 1.47x as often. A uniform
multiplier across every size class points at one decision rule applied too readily, not
at a structural blind spot.

## Where SPAdes refuses and we do not

`PairedResolver::resolve` in `src/resolve.cpp` reaches a chain end, enumerates the
continuations, and scores each by read pairs. When the best score does not clear
`linkBar` it does not stop. It tries three fallbacks in order:

1. `matchingAgrees` — a 2-in/2-out repeat short enough that no pair can span it, where
   the two possible perfect matchings are scored against each other and the intended one
   must dominate by 3x. Structurally sound: a wrong pairing has to beat the right one
   twice.
2. `allSameDestination` -> `pickByInterior` — every candidate ends on the same oriented
   unitig, so the destination was never in doubt and only the route through the repeat
   is being chosen. Also sound: the join itself is certain.
3. **`pickByCoverage`** — take the candidate whose terminal unitig sits within 25% of the
   chain's own depth, provided every rival is at least 60% off it.

The third is different in kind from the first two. It decides *where the chain goes* on
depth similarity alone, with no paired evidence and no structural constraint. The code's
own comment says so: it "fires only where the paired reads had nothing to say, so it has
no second opinion to check itself against."

SPAdes has written this heuristic twice, and ships it disabled for isolates.

* `SimpleCoverageExtensionChooser` (`extension_chooser.hpp:284`). Enabled by
  `pe_params.info:87` -> `simple_coverage_resolver { enabled false }`, and by
  `extenders_logic.cpp:508`, which only constructs it when `pset.multi_path_extend` is
  also set. `multi_path_extend` is `false` in `pe_params.info:18` and `true` only in
  `rna_mode.info:132`.
* `CoordinatedCoverageExtensionChooser` (`extension_chooser.hpp:1346`). Gated on
  `use_coordinated_coverage`, which is `false` in `pe_params.info:78` and `true` only in
  `meta_mode.info:202`.

`isolate_mode.info` sets `mode isolate` and includes `careful_mode.info`, which touches
only graph simplification — bulge removal, tip conditions — and overrides neither flag.
So for `--isolate` and for the default bacterial pipeline, **SPAdes performs no
coverage-based branch selection at all.**

What it uses instead is `SimpleExtensionChooser`, and its refusal is explicit
(`ExcludingExtensionChooser::FindFilteredEdges`, `extension_chooser.hpp:427`):

```cpp
AlternativeContainer weights = FindWeights(path, edges, to_exclude);
auto max_weight = (--weights.end())->first;
EdgeContainer top = FindPossibleEdges(weights, max_weight);
EdgeContainer result;
if (CheckThreshold(max_weight)) {     // below the weight threshold -> empty
    result = top;                     // ties within prior_coeff_ -> more than one
}
return result;
```

`FindPossibleEdges` keeps every candidate scoring above `max_weight / prior_coeff_`, so a
near-tie returns two edges, and `SimpleExtender` extends only on exactly one. Two ways to
stop, both of them a refusal: the evidence is too weak, or it does not separate the
candidates. There is no third branch.

That is the answer to "how does SPAdes make the smart decision": at an ambiguous branch
it does not make one. It leaves the sequence as a separate contig and accepts the shorter
N50. Its coverage heuristics exist, were implemented, and are switched off for isolate
data — which is a considered judgement by their authors that on a genome at roughly
uniform depth, "this candidate is at my coverage" separates almost nothing.

Note what SPAdes' disabled chooser does when it *is* enabled, because it is not what
`pickByCoverage` does. It requires exactly two candidates; it walks back along the path to
find a vertex of in-degree 2 and refuses if there is none; it fetches the *sibling entry*
into that vertex; it refuses if the two entries have similar coverage or the two exits do;
and then it matches by **rank**, not by value — if my entry is the deeper of the two
entries, take the deeper of the two exits. It is a perfect matching over a 2-in/2-out
structure, scored on coverage. `matchingAgrees` in `resolve.cpp` is that same matching
scored on read pairs. `pickByCoverage` is neither: it looks at no sibling, requires no
structure, and compares absolute depths.

## Status

The mechanism above is read from the source and the shipped configs, and is not in
question. Whether removing `pickByCoverage` is a net improvement for TesserACT is a
separate empirical question — it buys contiguity on libraries whose mates overlap, which
is most of this panel, and that is why it was added. That measurement is in progress on a
27-isolate subset (20 drawn at random, plus the 10 isolates where TesserACT is both more
misassembled and shorter than SPAdes) and this document will record the result, whichever
way it goes.
