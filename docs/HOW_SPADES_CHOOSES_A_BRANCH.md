# What SPAdes actually does at a branch, read from its source

TesserACT loses one metric to SPAdes reproducibly, on two independent panels:

| panel | arm | misassemblies, TesserACT vs SPAdes | p |
|---|---|---|---|
| 666 *K. pneumoniae*, no model | vanilla | 151 / 218 (win/loss) | 0.002 |
| 151 *S. aureus*, non-clonal | `--organism` | 26 / 42 (win/loss) | 0.048 |

Same direction, different organism, different model state. That reproducibility is the
reason to look for an algorithmic cause rather than tuning a threshold — but the *S.
aureus* row is weaker than it looks, and the next section takes it apart before drawing
anything from it.

## How big the gap really is, once contiguity is controlled

**Updated at the full cohort (n=180). The earlier version of this section, written at
n=151, concluded that at matched contiguity there was no significant difference. Twenty-nine
more isolates reversed that, and the conclusion below is the one that stands.**

Raw counts over the complete cohort: 223 misassemblies to SPAdes' 150, an excess of 73.
Splitting by how much more contiguous our assembly is than SPAdes' on the same isolate:

| NGA50, model vs SPAdes | n | model | SPAdes | excess | win/loss | p |
|---|---|---|---|---|---|---|
| SPAdes ahead (< −10%) | 35 | 33 | 27 | +6 | 9/10 | 0.702 |
| **matched (−10%..+10%)** | **48** | **86** | **70** | **+16** | **5/17** | **0.022** |
| we lead 10–50% | 50 | 44 | 25 | +19 | 7/15 | 0.072 |
| we lead >50% | 47 | 60 | 28 | +32 | 8/11 | 0.445 |

Two things are true at once, and the n=151 reading got the second one wrong.

**Part of the gap is bought contiguity.** The excess grows monotonically with how far ahead
we are — +6, +16, +19, +32 — and the largest single block sits in the band where we are more
than 50% more contiguous. A longer contig crosses more junctions and has more chances to
cross one wrongly.

**But a residual gap survives matching.** In the contiguity-matched band the excess is +16
over 48 isolates at 5 wins to 17, **p=0.022**. That is a real difference in decision quality
at equal contiguity, and it agrees with the 666-isolate *K. pneumoniae* benchmark
(151/218, p=0.002) rather than contradicting it. At n=151 the same band read 5/14, p=0.067,
and was reported here as "not significant"; that was a sample-size statement, not a result.

The heavy tail is still worth naming, because it inflates the headline: GCF046742145v1 alone
contributes 29 of the 73 (30 against SPAdes' 1) while leading NGA50 by 269% — there SPAdes
reached 67 kb against our 248 kb and collected one misassembly by producing fragments too
short to be wrong. GCF003184985v1 adds 5 more (13 v 8). But removing them no longer removes
the effect, because the matched band does not depend on them.

Relocation sizes show we have no distinctive blind spot — the two distributions are the same
shape, so this is one decision made too readily rather than a specific structure we mishandle:

| relocation size | TesserACT | SPAdes |
|---|---|---|
| 1–5 kb | 80 (48.2%) | 52 (46.0%) |
| 5–20 kb | 13 (7.8%) | 13 (11.5%) |
| 20–100 kb | 15 (9.0%) | 11 (9.7%) |
| >=100 kb | 58 (34.9%) | 37 (32.7%) |

(measured at n=151; the shape does not change at 180.)

## What the model costs here

The model is not innocent in this. Against the no-model arm over the same 180 isolates it
buys a great deal of contiguity and pays for it in exactly this coin:

| model vs base | median | win/loss | p |
|---|---|---|---|
| NGA50 | 247,188 vs 216,080 (**+14.4%**) | **101/0** | 2.7e-18 |
| NG50 | 253,984 vs 223,999 (+13.4%) | 103/0 | 1.3e-18 |
| genome fraction | 98.83 vs 98.80 | 97/28 | 6.1e-07 |
| # misassemblies | 0 vs 0 | **2/18** | **0.0011** |

101 wins and no losses on NGA50 is as clean as a result gets. It comes with 18 isolates
made worse on misassemblies against 2 made better. The model's adjacency table is telling
the resolver to make joins the paired reads would not have made on their own, and most of
those joins are right.

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

## Which of our fallbacks actually fires

`resolve.cpp` tries three fallbacks when paired support does not clear `linkBar`, plus a
fourth escape inside the tie branch. Instrumented counts, per isolate:

| isolate | matched | sameDest | depthPick | tieEscape | enum hit cap |
|---|---|---|---|---|---|
| GCF026547035v1 | 0 | 47 | 0 | (not counted) | 0 |
| GCF022693245v1 | 0 | 5 | 5 | 17 | 0 |
| GCF005931015v2 | 8 | 18 | 5 | 8 | 0 |
| GCF054392195v1 | 10 | 35 | **38** | 98 | 0 |

`pickByCoverage` reads 0 on the first isolate and 38 on the fourth, where it is the
largest fallback. One isolate is not enough to call a code path dead, and an earlier
draft of this document did exactly that.

`enumHitCap` is 0 throughout, so no decision here was taken over a candidate list
truncated at `kMaxCandidates` — `allSameDestination` is seeing the real candidate set,
not an artefact of the cap.

## Ablations

Screen: the 14 isolates carrying >=2 model-arm misassemblies, excluding GCF022832835v1
(61 against SPAdes' 58 — both assemblers fail there, so it is the isolate). Control is
the campaign's existing `model` arm; the probe binary reproduces it byte-for-byte
(md5 `77374e51a82ca5a8bb6035572156eea4` on GCF026547035v1), so no control was re-run and
no campaign output was touched.

**`--no-gapfill`** (n=10 scored). Gap closing is the one join source that converts a
scaffold adjacency into contig-level sequence before the N-split, and its 3,000 bp
ceiling matches the 1–5 kb relocation bucket.

| metric | model | nogapfill | win/loss | p |
|---|---|---|---|---|
| NGA50 | 154,834 | 150,699 | 6/0 | 0.036 |
| # contigs | 86 | 94 | 10/0 | 0.006 |
| genome fraction | 99 | 99 | 8/2 | 0.067 |
| misassemblies, total over the ten | **65** | **56** | 1/4 | n/a |

Gap filling accounts for 14% of the misassemblies on the worst isolates and pays for them
with contiguity on every isolate tested. A real but minor contributor, and the trade is
not clearly worth taking. **Not the cause.**

**`TESSERACT_NO_DEPTH_PICK`** — removing `pickByCoverage` entirely.

| isolate set | n | NGA50 identical | genome fraction identical | misassemblies |
|---|---|---|---|---|
| the 14 worst-misassembly isolates | 14 | 14/14 | — | 73 -> 75 |
| unselected isolates | 9 | 9/9 | 9/9 | identical 9/9 |

On 23 isolates the fallback changes **nothing measurable**, to the base pair, and where it
does move anything it moves it the wrong way. It is not idle: it fires 5, 5 and 38 times
on the three isolates instrumented above, and the contigs it produces genuinely differ
(md5 differs on GCF054392195v1 and GCF046742145v1). It simply arrives at the same
assembly.

That makes it neither a cause of misassemblies nor a source of contiguity: roughly forty
lines of heuristic that the code's own comment justified on reasoning ("a chain running at
45x continues into sequence at 45x") that measurement does not support. Removing it is a
simplification, not an optimisation, and it lands TesserACT where SPAdes already is --
they wrote the same heuristic and disabled it.

## Status

The source reading above — that SPAdes ships both coverage choosers disabled for isolates
and refuses at an ambiguous branch — is settled and not in question.

What is *not* established is that this costs TesserACT anything measurable on *S. aureus*.
At matched contiguity the difference is not significant, and the headline p=0.048 does not
survive dropping one isolate. Two candidate causes have been tested and neither explains
the gap: gap filling accounts for 14%, and `pickByCoverage` was wrongly written off from a
single isolate and is still under test.

The defensible conclusion today is narrower than the one this document originally reached:
TesserACT and SPAdes sit at different points on the same contiguity/correctness trade,
TesserACT further along it, and most of the misassembly difference on this panel is that
displacement rather than a worse decision rule. Whether a residual difference exists needs
the 666-isolate *K. pneumoniae* panel or the remaining six ESKAPEE cohorts, not this one.
