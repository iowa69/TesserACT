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

The raw counts are 177 misassemblies to SPAdes' 118 over 151 isolates, a 1.50x ratio.
Two things cut that down, and both were found after the first version of this document
claimed the multiplier was uniform. It is not.

**Most of it is bought contiguity.** Splitting the cohort by how much more contiguous our
assembly is than SPAdes' on the same isolate:

| NGA50, model vs SPAdes | n | model | SPAdes | excess | win/loss | p |
|---|---|---|---|---|---|---|
| SPAdes ahead (< −10%) | 30 | 27 | 22 | +5 | 7/7 | 0.754 |
| **matched (−10%..+10%)** | **43** | **82** | **70** | **+12** | **5/14** | **0.067** |
| we lead 10–50% | 40 | 25 | 13 | +12 | 7/12 | 0.184 |
| we lead >50% | 38 | 43 | 13 | +30 | 7/9 | 0.587 |

Half the excess (+30 of +59) sits in the band where we are more than 50% more contiguous
than SPAdes. In the contiguity-matched band the difference is **not significant**. A
longer contig crosses more junctions and therefore has more chances to cross one wrongly;
that part of the gap is the price of the contiguity, not a defect.

**And the significance is fragile.** One isolate carries half the excess:

| isolate | model | SPAdes | NGA50, model vs SPAdes |
|---|---|---|---|
| GCF046742145v1 | 30 | 1 | 248,358 vs 67,253 (**+269%**) |
| GCF022693245v1 | 5 | 0 | 58,027 vs 112,053 (−48%) |
| GCF013836745v1 | 4 | 0 | 76,446 vs 65,471 (+17%) |

| cohort | excess | p |
|---|---|---|
| all 151 | +59 | **0.048** |
| minus GCF046742145v1 | +30 | 0.071 (ns) |
| minus that and GCF022693245v1 | +25 | 0.102 (ns) |

Drop one isolate and the result stops being significant. On GCF046742145v1 SPAdes did not
really assemble the genome — NGA50 67 kb against our 248 kb — and collected one
misassembly by producing fragments too short to be wrong. That is not SPAdes making a
better decision; it is SPAdes making fewer decisions.

So the honest statement for *S. aureus* is: **at matched contiguity there is no
significant difference**, and the headline p=0.048 rests on a single isolate. The
666-isolate *K. pneumoniae* benchmark (151/218, p=0.002) is the more robust evidence that
a real residual effect exists; this 151-isolate panel is too small to resolve it, and the
first version of this document overstated what it showed.

Relocation sizes are still worth recording, because they show we do not have a distinctive
blind spot — the two distributions are the same shape:

| relocation size | TesserACT | SPAdes |
|---|---|---|
| 1–5 kb | 80 (48.2%) | 52 (46.0%) |
| 5–20 kb | 13 (7.8%) | 13 (11.5%) |
| 20–100 kb | 15 (9.0%) | 11 (9.7%) |
| >=100 kb | 58 (34.9%) | 37 (32.7%) |
| **total** | **166** | **113** |

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

`TESSERACT_NO_DEPTH_PICK` is running on the same 14 and will be recorded here.

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
