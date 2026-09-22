# Experiment ledger

Every change tried against TesserACT during the seven-organism ESKAPEE campaign, what it
was measured at, and what was decided. Negative results are kept deliberately: the cost of
re-testing a dead idea is the same as testing a live one.

**Protocol for everything below.** Each arm is a contrast against a control that reproduces
the then-shipped binary **byte-for-byte** (verified by md5: `77374e51…` pre-AD,
`4dee0032…` post-AD). Cohorts are dereplicated with `mash` at `d <= 0.0005` so no two
isolates are clonal. Scoring is QUAST `-s` with **only the `_broken` column read**, each
isolate against **its own closed RefSeq reference**, so scaffold gaps are never counted as
contiguity. Quality gates discard an isolate when the reference and the reads disagree
(>5% unaligned in every arm, or genome fraction <80%) — a database error must not be
recorded as an assembly error.

---

## 1. Shipped, with the measurement that justified it

| change | effect | n | evidence |
|---|---|---|---|
| **A — k ladder from the read-length distribution** | NGA50 **+26%** on severely truncated ladders | 16 | 13 better / 0 worse, p=0.002 |
| **D — tip-clipper absolute branch** | genome fraction **18/6** (p=0.009), contigs **−5.8%** (16/7, p=0.015) | 38 | — |
| **AD together** | genome fraction 22/7 (**p=0.001**), largest alignment **+14.5%** (p=0.038), contigs **−11.4%** (p=0.031), **no structural cost** (misassemblies 3/2, duplication 5/4) | 37 | commit `5b84c95` |
| **Plasmid copy-number exemption** | whole-plasmid recovery **33% → 44%**, low-copy control band unmoved exactly as predicted | 39 | +1 misassembly total |
| **1,500 bp replicon floor + corroboration + propagation guard** | pooled precision **0.225 → 0.712**, recall 1.000 → 0.794 | 185 | — |

**A, in detail.** `baseKLadder` binned on the **mean** read length with a cliff (>=165 → nine
rungs to k=127; >=140 → five to 95; >=100 → four to 55). The rule was *anti-correlated with
what it protected*: the library retaining the **most** k-mer mass at k=127 (35%) was capped
at **k=55**, while a fixed-length 151 bp library retaining the **least** (17%) was given 95.
Submitter-trimmed libraries are bimodal — a low mean with a tail at the instrument's full
cycle count. Now the top rung is the largest retaining 35% of the k=21 k-mer mass.
Independent confirmation from a natural experiment: of 68 isolates at max read length 151,
the 12 given the k=55 ladder score **0.815** against SPAdes where the 56 given k=95 score
**1.330** (Cliff d −0.562, p=0.0024).

**D, in detail.** SPAdes runs two tip rules; TesserACT shipped only the relative one at
3.5×. SPAdes' second branch reaches 10× with **no competitor requirement**
(`tip_clipper.hpp`, `"{ tc_lb 10., cb auto }"`). A tip 250–850 bp at 2–7× beside a
competitor at 3× survived our rule forever, and every survivor is a junction that blocks
compaction — i.e. a contig break.

---

## 2. Rejected on measurement

| change | result | why it fails |
|---|---|---|
| **B — trusted carry-over** (SPAdes' `TrustedContigs`) | NGA50 **4/11**, genome fraction 10/22 | SPAdes fills coverage from a **separate stream** (`construction.cpp:110`, *"Has to be separate stream for not counting it in coverage"*). Here the **count IS the coverage**, so a k-mer admitted at the cutoff reads as error-level to every rule downstream and is pruned. The two mechanisms are not interchangeable. |
| **E — fitted per-rung EC threshold** | NGA50 **6/23** (p=6.9e-04), mismatches **+103%** | The idea is right — SPAdes fits this per rung and its values swing 7/10/4/12/2 across one ladder where ours is fixed. My approximation of the fit is wrong, and a bad fitted threshold is worse than a crude fixed one. Wired in before checking it reproduced SPAdes' values; that was the wrong order. |
| **ABC** | NGA50 10/11, p=0.835 | A's gain exactly cancelled by B's harm — also a useful check that the factorial composes as the individual arms predict. |
| **Unicycler best-rung selection** | NGA50 **0 better / 5 worse** | Unicycler builds each k **independently** and cleans before scoring, so its "best k" is a real choice among peers. Our rungs are **iterative with carry-over** — the last rung *is* the accumulation. Choosing an earlier one discards work. |
| **Naive ladder union (concat)** | duplication **2.99**, total length 7.9 Mb for a 2.8 Mb genome | Re-emits sequence three times. Matches the literature: MetaHipMer's 12-way coassembly went duplication 1.2 → 7.3 with genome fraction halved; `dedupe` only recovers 7.3 → 4.3. |
| **Adapter read-through trimming** | condition is **3.3%** of pairs | Implemented and correct; the target barely exists. 69.2% of pairs have a *perfect* overlap but only 3.3% have a fragment shorter than R1. |
| **Mate merging (via fastplus/scepter)** | genome fraction **94.94 → 81.07** on one isolate | +2–5% on two isolates, catastrophic on the third. |
| **Shorter k ladder (≤77)** | +2% on one isolate, **−4% and −10%** on the others | Our tall ladder earns its keep. |
| **`pickByCoverage` (depth-continuity fallback)** | every metric **identical to the base pair** on 23 isolates | Fires up to 38×/isolate and changes the contigs, but arrives at the same assembly. ~40 lines justified on reasoning measurement does not support. SPAdes ships the same heuristic **disabled**. |
| **Scaffolder symmetric-target fix** | joins **5 → 5** | `reach` bounds the window to ~600 bp, so it widened the target from 2 unitigs to 4. The terminus restriction it attacked is **load-bearing**: a pair linking to the *middle* of another chain is evidence of a repeat, not of a join at that chain's end. Reverted (96 lines). |
| **Gap closing on the raw graph, sequence only** | **37,685** bridges against a ceiling of 804 | Paired nomination is not optional — it is the other half of the conjunction. |

---

### The ladder union, settled at n=38

The hypothesis was: *"we emit contigs from the final graph only, so the highest rung decides
what survives — maybe the solution is to emit the union over its whole ladder, which is why
a high top-k doesn't cost it anything."*

Implemented as `rescueLadderContigs()`: after the final rung, walk back down the ladder and
append any contig a smaller k resolved that the final graph no longer contains, deduplicated
on 31-mer containment, longest-first, with each rescued contig registered as emitted so two
near-copies from adjacent rungs cannot both pass. Run as arm **C** of the change factorial
over 38 *S. aureus* isolates against each isolate's own closed genome.

**Result: 37 of the 38 output files are byte-identical to base0.** Across the whole panel
the union rescued **2 contigs, 5,697 bp, on one isolate**:

| metric | base0 | arm C |
|---|---|---|
| # contigs | 76 | 78 |
| total length | 2,893,898 | 2,899,595 |
| genome fraction | 98.540 | **98.707** (+0.167 pp) |
| NGA50 / NG50 | 125,253 | 125,253 (unchanged) |
| duplication | 1.005 | 1.005 |
| # misassemblies | 0 | 0 |

So the mechanism is *correct* — the two contigs it recovered are real genome, placed without
a misassembly and without inflating duplication — and **almost always inert.**

**Why, measured rather than argued:** the premise is false for TesserACT. The highest rung
does not decide alone, because our ladder is **iterative with carry-over** — every rung's
contigs are fed into the next rung's k-mer count, so what a small k resolved is already
carried into the final graph. There is nothing left for a union to rescue. This is the same
reason Unicycler's best-rung selection lost 0/5 here: for a carry-over ladder the last rung
*is* the accumulation.

The half of the hypothesis that was right is the other half — *"a high top-k doesn't cost
anything"*. It doesn't, and the way to cash that in is not to union the ladder but to
**choose the ladder better**, which is arm A, and with the tip rule (D) it is the change
that shipped:

| metric (n=38) | base0 | **AD** | win/loss | p |
|---|---|---|---|---|
| NGA50 | 244,916 | **259,636** (+6.0%) | 13/4 | 0.037 |
| NG50 | 247,582 | **269,365** (+8.8%) | 12/4 | 0.052 |
| genome fraction | 98.828 | **98.856** | 22/8 | **0.0021** |
| # contigs | 34.5 | **31.0** | 17/9 | 0.022 |
| LGA50 | 4.5 | **4.0** | **9/0** | **0.0092** |
| # misassemblies | 0 | 0 | 3/3 | 0.75 |
| duplication | 1.002 | 1.002 | 6/4 | 0.54 |

Nine isolates improve LGA50 and none worsen. Misassemblies and duplication do not move —
this contiguity is not bought with correctness.

Arm C is kept behind `--ladder-union` (off by default): it is free, it has never harmed a
metric, and on the one isolate in 38 where the ladder does drop something it recovers it.
It is not a lever.

---

### Duplication: where our only hard metric loss came from, and the fix

On the complete *S. aureus* cohort, `model` loses QUAST duplication ratio to SPAdes
**1.003 vs 1.000, 1 win / 135 losses, p=1.3e-23** — one of only two metrics we lose.
Duplication counts a reference base covered by two contigs twice in the numerator and once
in the denominator, so the question is literally *which contigs sit on top of each other*.
Decomposed from QUAST's own alignment tables, 185 isolates, both arms:

| | TesserACT `model` | SPAdes |
|---|---|---|
| total duplicated reference bases | 2,255,627 | **322,080** |
| median per isolate | 6,271 | **1,014** |
| **contigs wholly contained in another** | **1,057,040 bp (46.9%), 856 contigs** | **7,183 bp (2.2%), 11 contigs** |
| partial, end-abutting overlaps | 1,198,587 bp | 314,897 bp |

Two separate causes, in almost equal parts.

**Cause 1 — we emit contigs that are exact substrings of other contigs; SPAdes does not.**

| | contigs >=500 bp | exact substring of a longer contig |
|---|---|---|
| TesserACT | 1,411 | **87 (6.2%)**, 3.5 per isolate |
| SPAdes | 1,554 | **0 (0.0%)** |

Concretely: `NODE_44_length_1615_cov_74.5_unk` in GCF000010465v1 is an exact substring of
`NODE_25_length_568084_cov_39.9_chr_circular`, at 1.87x the chromosomal depth — a two-copy
repeat. The resolver resolves the repeat into its flanking chains **and then also emits the
standalone repeat unitig**. Present in `base` too (43.6%), so it is core emission, not the
model.

**Cause 2 — our contig-end overlaps are individually far longer than SPAdes'.** Almost the
same number of overlapping pairs (2,598 vs 2,642) carrying 3.8x the bases:

| overlap length | model pairs / bases | SPAdes pairs / bases |
|---|---|---|
| <= 126 bp | 1,659 / 75,810 | 2,358 / 153,806 |
| 127-1000 bp | 764 / 319,608 | 989 / 142,852 |
| **> 1000 bp** | **540 / 800,088 (67%)** | **1 / 1,026** |

**70% of SPAdes' duplicated bases are a single overlap of exactly its final K** — the
unavoidable de Bruijn overhang. **90.6% of ours sit in overlaps longer than 2k**, 90% of
them terminal on *both* contigs, at >=99% identity, with lengths that recur to the base pair
across unrelated isolates at canonical *S. aureus* IS sizes (1324 = IS256, ~1512 = IS1181,
790 = IS431). At the 203 loci where both assemblers duplicate, our overlap is a median
1,324 bp against SPAdes' 77 bp — the same junction, 17x the overhang.

Adopting SPAdes' k ladder changes nothing (5,128 -> 4,974 duplicated bases per isolate,
matched 18 isolates), so **the k ladder is not the mechanism.**

**The fix, and the fix we deliberately did not take.** The obvious reading is "stop
extending into the repeat", which is what produces SPAdes' clean exactly-K signature. That
is rejected: extending *through* repeats is precisely where our NGA50 advantage (249 kb vs
177 kb) comes from, and copying SPAdes here would trade away the win to fix the loss.
Instead, keep the extension and drop the redundant second copy:

* `--dedup-contained` — drop a contig that is an exact substring, forward or
  reverse-complement, of a longer contig. Cannot lose sequence by definition. Near-copies
  that share a core but differ are real variant sequence, so this is a substring test and
  not a k-mer containment test. Plasmid calls and circular contigs are never dropped: a call
  is a result, not merely bases.
* `--trim-overlap N` — trim an exact contig-end dovetail longer than the final k, from the
  **shorter** partner only, so the longer contig keeps it and the assembly still contains
  every base once.

**Measured, 44 isolates, against each isolate's own closed genome:**

| metric | control | dedup | trim | **dedup+trim** | SPAdes |
|---|---|---|---|---|---|
| duplication | 1.0020 | 1.0010 | 1.0010 | **1.0000** | 1.0000 |
| # contigs | 39.5 | 34.0 | 38.5 | **34.0** | 39.5 |
| genome fraction | 98.783 | 98.783 | 98.737 | **98.703** | 98.272 |
| NGA50 | 249,267 | 249,267 | 249,267 | **249,267** | 176,653 |
| LGA50 | 4 | 4 | 4 | **4** | 5 |

Against SPAdes, the loss becomes a dead heat and **no win is given up**:

| metric | control vs SPAdes | dedup+trim vs SPAdes |
|---|---|---|
| **duplication** | **0/35, p=2.6e-07** | **9/6, p=0.887** |
| genome fraction | 41/3, p=2.9e-08 | 39/5, p=2.3e-07 |
| NGA50 | 34/9, p=8.5e-05 | 34/9, p=8.9e-05 |
| # contigs | 20/21, p=0.86 | 28/13, p=0.13 |
| LGA50 | 29/9, p=0.028 | 29/9, p=0.028 |
| mismatches /100 kbp | 28/13 | 28/13 |

**The price, stated honestly.** Some of what QUAST calls our duplication is legitimately
covering *both* copies of a two-copy repeat, so removing the second copy removes that
coverage. Measured exchange rate:

| variant | duplicated bases removed | unique reference coverage lost | ratio |
|---|---|---|---|
| dedup | 128,641 | 19,675 | 6.5 : 1 |
| trim | 480,157 | 64,995 | 7.4 : 1 |
| dedup+trim | 582,999 | 106,149 | 5.5 : 1 |

A threshold sweep (trim only overlaps > 300 / 600 / 1200 / 2400 bp) found **no sweet spot** —
the rate is flat at 2.4-2.5 : 1 across the range and only reaches 3.2 : 1 at >2400, where
duplication no longer reaches parity. There was nothing to tune, only a decision. It is
taken because the genome-fraction lead survives at p=2.3e-07 and every other win is
untouched to the decimal.

### The `model` arm's `_cov_` field is not per-contig coverage

Found while testing the above. In the `model` arm all pieces of one scaffold inherit the
scaffold's single depth: in GCF002025145v1 all 15 chromosomal contigs report
`cov_59.8445`, a 4,570 bp piece and a 663,293 bp piece alike. Distinct cov values per
contig: `model` 25/39, `base` **41/41**, SPAdes **63/63**. So it is introduced by the
scaffold gap split, which reuses the scaffold's `cov` for every piece.

It affects no QUAST metric, but it is the shipped `--organism` path and anyone parsing the
header — for binning, contamination screening, or plasmid copy number — gets a constant.
Queued behind the two metric fixes.

---

## 3. Diagnoses that redirected the work

**The corrector is not the gap.** Assembling SPAdes' own BayesHammer-corrected reads with
TesserACT closes **13% of the gap on one isolate and 0.5% on another**, and our corrector
adds nothing on top of theirs. This **falsifies** the standing conclusion in
`WHY_SPADES_SOMETIMES_WINS.md`, which now carries a superseded banner.

**We degrade faster than SPAdes as graphs get harder.** Strongest single predictor of
losing is scaffold-gap count (Spearman **−0.402**, p≈1e-8); the NGA50 ratio falls
**1.60 → 1.47 → 1.12 → 0.90** across difficulty quartiles. *E. faecium* is that prediction
landing on a whole organism: median NGA50 45 kb against *S. aureus*'s 247 kb, and the
misassembly deficit goes from 29/53 (p=0.011) to **14/75 (p=2.2e-09)**.

**The resolver is not the limit.** On the worst isolates the continuation counters read
`no-candidate=984` of 1,018 chain ends (**97%**) — the link-support bar, tie rule and
coverage fallback are never reached. There is nothing to decide.

**The dead ends are intrinsic.** They stay flat across the whole ladder (1,338 → 632) while
unitigs collapse tenfold (6,096 → 623). Not created by climbing to k=127, and N50 improves
to the top rung, so stopping earlier only costs contiguity.

**`joinDeadEnds` never fires.** 0 joins across **36 graphs and ~10,000 dead ends**. Its
`minOverlap = k-1` requirement (`graph.cpp:917`) is unsatisfiable in practice. The
misassembly risk the source review flagged is therefore theoretical — and so is its
benefit. `filterByReadDepth` is likewise never called.

**The model is not reckless.** On *E. faecium* it makes **9,872 joins and gets ~75 wrong —
99.24% correct** — while winning NGA50 on 90 isolates and losing on **none**. The
per-isolate ratio (1 better / 40 worse on misassemblies) is volume arithmetic at ~80 joins
per isolate, not carelessness. Tightening the bar would suppress thousands of correct joins
to remove 75 errors.

---

### The misassembly gap, decomposed — and a previous reading withdrawn

`HOW_SPADES_CHOOSES_A_BRANCH.md` reported a residual at matched NGA50: 48 isolates,
excess +16, 5/17, **p=0.022**, measured at n=180 on the pre-AD binary. **On the current
binary it does not reproduce** (p=0.074 raw, 0.089 excluding the four isolates where both
assemblers fail). More importantly the control itself was wrong: a sub-2 kb wrong tail
costs essentially no contiguity, so a defect of that shape hides *inside* the matched band
rather than showing up outside it. NGA50-matching cannot see it.

The control that does work is the length of the **shorter alignment block flanking the
breakpoint**. Re-derived independently from QUAST's own stdout, n=177, arms mtime-gated to
one binary:

| shorter flank | model | base | SPAdes | excess | p (model v SPAdes) |
|---|---|---|---|---|---|
| **< 2 kb** | 192 | 185 | 100 | **+92** | **1.9e-04** |
| 2 - 20 kb | 137 | 120 | 140 | −3 | 0.74 |
| **>= 20 kb** | 71 | 59 | 40 | **+31** | 3.3e-04 |

Two disjoint classes with a dead zone between them, and they behave oppositely.

**The >=20 kb class is bought contiguity plus the model's joins.** It concentrates where we
lead, and model (71) > base (59) > SPAdes (40). A legitimate trade, to be argued rather
than fixed.

**The <2 kb class is a genuine core defect.** It is *identical* in `base` and `model`
(192 v 185, 5/7, **p=0.50**), so it is not the organism model. And in the 29 isolates where
**SPAdes is the more contiguous assembler**, it is still +11 at 2/10, **p=0.031** — so it is
not bought contiguity either. 81% of those short blocks sit at a contig terminus (SPAdes:
30%), their lengths cluster in two humps at 0-250 bp and 1,000-1,500 bp — IS-element scale —
and the sequence they belong to is a median **209 kb** away (SPAdes: 6.7 kb). We append a
short dispersed repeat from the wrong copy and stop there.

**Also ruled out, each with its measurement.** Contig-level coverage as a misassembly
signature: median coverage ratio is **1.000 in both** misassembled and correct contigs, and
identical once length-matched — do not build a filter on it. Plasmid/chromosome chimerism:
6 translocations in the whole cohort. Gap closing or short-range error: local misassemblies
are flat and 7 *in our favour*. A relocation-size blind spot: the distributions are the same
shape. Every cohort covariate — depth, read length, contiguity advantage, contig count —
all p > 0.06.

**The change under test.** `TESSERACT_MIN_FALLBACK_DEST=2000`: withdraw the three
`PairedResolver` fallbacks (`matchingAgrees`, `allSameDestination`, `pickByCoverage`) when
the continuation they would take lands on a unitig shorter than 2 kb. Paired support that
clears `linkBar` outright is untouched. 24 isolates — the 14 carrying 155 of the 192 short-
flank events, plus 10 drawn at random as a contiguity-neutrality check.

Control validation for this panel is unusually strong: **24 of 24 `z0` arms (the same
experimental binary, every new feature off) are byte-identical to the campaign's own `model`
arm**, so any difference is attributable to the guard alone. The guard fires on 10 of the 24,
changing contig counts by 0 to +6 — the expected shape, since a refused join leaves two
pieces.

**Result: a clean negative. The guard does nothing.**

| metric | control (z0) | short-dest guard | win/loss | p |
|---|---|---|---|---|
| # misassemblies (total) | **306** | **304** | 1/0 | — |
| **<2 kb flank class** | **166** | **164** | 1/0 | — |
| 2–20 kb class | 107 | 107 | 0/0 | — |
| >=20 kb class | 33 | 33 | 0/0 | — |
| # contigs | 40.0 | 41.0 | **1/9** | **0.0125** |
| NGA50 | 173,695 | 164,810 | 0/2 | — |

Two misassemblies removed from the class it was designed for, out of 166, while costing
contigs on nine isolates and −10.4% NGA50 on one. **The inferred mechanism is wrong.** The
fallback path fires — contig counts move on 10 of 24 isolates — but the joins it withdraws
are not the joins that produce these misassemblies.

What that rules out: the `<2 kb`-flank misassemblies are **not** created by
`matchingAgrees` / `allSameDestination` / `pickByCoverage` authorising an extension onto a
short unitig. They survive the withdrawal of all three. So either the join clears `linkBar`
on genuine paired support that cannot distinguish two copies of a ~1.3 kb repeat, or the
chimera is already in the graph before resolution, or it is introduced by scaffolding or gap
filling — all of which are in `base`, which carries this class identically.

Recorded because the sequence matters: the agent that found the decomposition proposed
instrumenting the last extension step *first* and ablating second, and going straight to the
ablation cost a panel to learn the instrumentation was the necessary step.

### Stage ablation: the resolver owns 82% of the excess

Eight carriers holding 152 of the 166 `<2 kb` events, each stage disabled in turn:

| arm | `<2 kb` | total misassemblies | NGA50 | contigs |
|---|---|---|---|---|
| control (full pipeline) | **154** | 285 | 68,606 | 39 |
| `--no-scaffold` | 153 (−1) | 283 | 68,606 | 40 |
| `--no-gapfill` | 147 (−7) | 276 | 64,257 | 42 |
| **`--no-resolve` (graph only)** | **100 (−54)** | 219 | 61,661 | 66 |
| **SPAdes** | **88** | 213 | 67,493 | 47 |

**Our de Bruijn graph is very nearly as clean as SPAdes': 100 against 88.** The excess of 66
over SPAdes decomposes as **resolver 54, gap filling 7, graph 12**. Scaffolding contributes
nothing. It is concentrated: GCF031190615v1 24 → 50, GCF046742145v1 9 → 26 (SPAdes: 0),
GCF003184985v1 2 → 11 — three isolates carrying +52 of the +54.

The resolver is also what buys the contiguity on these isolates (NGA50 61,661 → 68,606,
contigs 66 → 39), so this is a bar-setting problem, not a code path to delete.

### Threshold tuning is ruled out as a family

Three independent knobs, none of which touches the class:

| knob | range tested | `<2 kb` events | cost |
|---|---|---|---|
| withdraw all three fallbacks | on/off | 166 → **164** | contigs 1/9, p=0.0125; −10.4% NGA50 on one isolate |
| `--min-link` (linkBar floor) | 6 → 8 → 12 → 16 | 147 → **147** | contigs 39 → 49 |
| `tieRatio` | 1.02 → 1.40 | 218 → **212** (13/8, p=0.322) | — |

The `tieRatio` row is free evidence: our own `careful` mode already runs at **1.40**, close to
SPAdes' `priority_coeff` of **1.50**, and is scored on all 185 isolates:

| arm | tieRatio | `<2 kb` | 2–20 kb | >=20 kb | NGA50 | contigs | genome fraction |
|---|---|---|---|---|---|---|---|
| base | 1.02 | 211 | 149 | 70 | 222,906 | 42 | 98.844 |
| model | 1.02 | 218 | 165 | 83 | 250,176 | 36 | 98.862 |
| careful | **1.40** | **212** | 148 | 67 | 246,903 | 38 | 98.859 |
| aggressive | 1.05 | 220 | 187 | 98 | 249,571 | 36 | 99.059 |
| SPAdes | **1.50** | **123** | 163 | 49 | 184,157 | 40 | 98.333 |

**What the three negatives say together:** these joins are *strongly supported* (a bar of 16
does not touch them) and *not near-ties* (a 1.4 discrimination ratio does not touch them).
The wrong candidate wins **decisively**. So the resolver's excess is not bad arbitration
between two visible options. Either the correct continuation is absent from the candidate
set, or the support for the wrong one is genuine but uninformative — pairs anchored on
sequence that is itself repetitive, which would support either copy equally.

One version of that was checked and ruled out: mates landing *inside* the traversed repeat do
not inflate the score. `scoreCandidate(tailFirst[m], terminal, ...)` scores chain-member to
**destination** pairs only; the interior route contributes only its length.

### What the breakpoint blocks actually are

Measured directly (31-mer probes from each `<2 kb` block, counted against the isolate's own
closed reference and against our own GFA segment depths):

| isolate | block len | copies in reference | block depth | assembly median | ratio |
|---|---|---|---|---|---|
| GCF003184985v1 | 1,063 | **7.0** | 697 | 70.7 | **9.9x** |
| GCF003184985v1 | 1,134 | **7.0** | 637 | 70.7 | **9.0x** |
| GCF003184985v1 | 1,700 | 0.5 | 76 | 70.7 | 1.1x |
| GCF022693245v1 | 1,975 | **5.0** | 182 | 29.6 | **6.2x** |
| GCF022693245v1 | 1,061 | 0.5 | 74 | 29.6 | 2.5x |

They **are** dispersed repeats — 5 to 7 copies in the reference — but they sit at **6-10x
median depth, far ABOVE the 1.6x repeat guard**, not below it. The hypothesis that they
escape by coverage dilution (which the code's own comment predicts) is **false**. Whatever
places them, it is not a guard they slip past on depth.

### Fifth negative: refusing fallbacks whose route no pair can span

Reasoning from the table above: a ~1,063 bp interior against an insert window of 204+/-140
(`maxPlausible` 764) cannot be spanned by any pair, so every candidate scores zero, `best <
linkBar` holds, and the decision falls to the coverage fallback. Implemented as
`TESSERACT_NO_UNSPANNED_FALLBACK`: refuse when the chosen route's interior exceeds
`insert_.maxPlausible`. Structural, not depth-based, so the dilution that defeats `isRepeat`
cannot defeat it.

**It fires hard and changes nothing.** Debug counters show 12, 10, 18, 21, 36 refusals per
isolate — comparable to or exceeding the coverage-fallback counts it replaces. Over 20 panel
isolates:

| arm | `<2 kb` | total misassemblies | NGA50 | contigs |
|---|---|---|---|---|
| control | 166 | 304 | 179,662 | 37.5 |
| no-unspanned-fallback | **166** | **304** | **179,662** | **37.5** |

Identical on 18 of 20 isolates. The refused joins are reconstructed downstream — by
scaffolding and gap filling, or by reaching the same junction from the opposite end.

### Five negatives, and what they jointly establish

| hypothesis | result |
|---|---|
| fallbacks authorise extension onto a short unitig | 166 -> 164 |
| paired-support bar too low (`--min-link` 6 -> 16) | 154 -> 153 -> 153 -> 155 |
| tie discrimination too loose (`tieRatio` 1.02 -> 1.40) | 218 -> 212, p=0.32 |
| blocks are repeats diluted below the 1.6x guard | **false** — they sit at 6-10x |
| fallbacks on routes no pair can span | 166 -> 166, output identical |

Two of these withdrew the fallback path by different tests and **neither moved the metric**,
while `--no-resolve` removes 54 of the 66 excess. That is only consistent with the
misassemblies arising on the resolver's **main, supported path** — joins that clear
`linkBar` on genuine paired evidence which cannot distinguish one copy of a 5-7 copy repeat
from another. No local rule in this code path separates them, because locally there is
nothing to separate: the evidence is real, it is just not specific.

**This is where guessing stops.** Every further attempt must start from instrumentation that
ties an individual join to an individual QUAST breakpoint — log per extension the candidate
count, best/second scores, destination and interior identity, and the evidence class, then
intersect with the breakpoint list. Five panels have now been spent testing plausible rules
blind; the sixth should not be.

---

### Why six resolver guards in a row did nothing: mutual choice

Instrumenting every accepted join (`TESSERACT_JOIN_TRACE`) on carrier GCF003184985v1 ended a
long run of failed experiments by showing the premise of all of them was wrong.

**381 joins traced. `by-coverage=0` — not one was fallback-decided.** So the three guards
that withdrew fallbacks (#12, #15, and the short-destination test) had nothing to withdraw on
this isolate. **Zero joins had an interior exceeding `maxPlausible`** (max 415 bp against a
762 bp window), so the unspanned-route theory was false as well.

What the trace did show: **51 of 381 joins (13%) have a single candidate**, and 25 of those
are accepted on **zero paired support**, every one routed through a unitig above the 1.6x
repeat threshold. `linkBar` and `tieRatio` are both gated on `cands.size() > 1`, so a lone
candidate had never been checked at all. That looked like the defect.

It is not. Guarding it (`TESSERACT_REQUIRE_SUPPORT_SINGLE`) fires 22-37 times per isolate and
changes **nothing**: over 14 isolates misassemblies 266 -> 266, `<2 kb` class 137 -> 137,
NGA50 identical on 14/14, contigs identical on 13/14.

Differencing the two traces explains it. The 25 refused joins are exactly the 25 that vanish,
**0 of them reappear in the reverse direction**, and all 25 share one destination — unitig 41,
nominated by five different sources (29, 33, 99, 173, 181), each with zero support. That is a
collapsed-repeat hub. And `resolve.cpp:1015` already requires:

```cpp
// Only join when both ends independently chose each other.
if (!back.ok || back.chainB != c || back.endB != e) continue;
```

Only one of the five can be chosen back, so four were discarded by arbitration before any
guard ran. 198 of 202 accepted joins are symmetric, confirming the rule binds. **The guards
were deleting candidates that had already lost.**

**What this establishes.** The misassembling joins are *mutually chosen* and *genuinely
supported* — both ends independently pick each other, on real read pairs, and are still wrong.
No local evidence rule separates them, because locally the evidence is correct: the pairs do
link those two sequences, in some copy of the repeat. That is why six threshold and guard
experiments in a row moved the metric by 0-2 events.

**Stop condition, restated.** Do not add a seventh guard to `bestContinuation`. The next step
must map individual accepted joins to individual QUAST breakpoints — the trace format now
supports it — so a wrong join is identified from the reference rather than from a rule
somebody hoped would separate them.

---

## 4. In flight

**Paired-read gap closer.** SPAdes' `late_gapcloser`, reimplemented as its own stage after
simplification. Two independent tests must both pass: an exact overlap of [10, k−2] bases
(the missing bases are already carried in the two flanks' terminal overhangs — nothing is
invented), **and** ≥2 read pairs vouching for that specific pair of ends, gathered through
a tip-neighbourhood map that walks back 5 kb along non-branching paths.

The false blocker is gone: `graph.cpp:917` claimed shorter overlaps require *"giving Link
its own overlap field first"*. They do not — SPAdes inserts an **edge carrying the missing
bases**, so every link still splices exactly k−1 and `mergeInto` is untouched.

**Pilot, two isolates, three arms each, QUAST against each isolate's own closed genome.**
Control is the probe binary with the feature off; it reproduces the shipped AD binary
byte-for-byte (md5 `4dee003211c1f1238c9912daed942891` on GCF026547035v1), so the three arms
differ only in the gap closer.

*GCF018093065v1* — the isolate with real gap-close activity (56 bridges under nomination,
337 under sequence alone):

| metric | control | sequence-only | **seq + paired vote** |
|---|---|---|---|
| # contigs | 166 | 139 | **146** |
| NGA50 | 35,731 | 40,474 (+13.3%) | **40,362 (+13.0%)** |
| NG50 | 35,731 | 40,489 | **40,444** |
| genome fraction | 98.030 | 98.161 | **98.124** |
| duplication | 1.001 | 1.000 | **1.000** |
| LGA50 | 21 | 19 | **19** |
| **# misassemblies** | **0** | **2** | **0** |

The nomination keeps **97.7% of the NGA50 gain and 74% of the contig reduction, and pays
none of the misassembly cost.** Mismatches and indels per 100 kbp are unchanged to two
decimals across all three arms, so nothing is being invented — these are joins, not edits.

*GCF022693245v1* — a quiet isolate (3 bridges under nomination, 5 under sequence alone):

| metric | control | sequence-only | seq + paired vote |
|---|---|---|---|
| # contigs | 136 | 133 | 135 |
| NGA50 | 58,027 | 58,625 | 58,027 |
| NG50 | 58,625 | 63,344 | 58,625 |
| largest alignment | 153,466 | 207,955 | 153,466 |
| # misassemblies | 5 | 5 | 5 |

Here the nomination is *too* conservative: sequence alone found five joins, gained NG50 and
a 54 kb longer largest alignment, and cost nothing, while the paired test admitted three and
gained nothing. That is the honest counterweight to the first table — requiring two vouching
pairs does discard real joins on libraries where the graph carries thin paired evidence at
dead ends (this one is 55x, 250 bp, against 196x, 301 bp for the first).

### The 16-isolate panel

14 isolates drawn at random from those whose `model` arm is already on AD, plus the two
pilots as a reproduction check. The pilots reproduced **byte-for-byte** — same contig
counts, same split points — so the panel config is provably the pilot config.

Control is each isolate's own post-AD `model` arm (15 of 16) or the probe with the feature
off (1, whose campaign arm is still pre-AD). All three arms scored by QUAST against that
isolate's own closed genome.

| metric | control | seq-only | seq+vote | sq w/l | p | nm w/l | p |
|---|---|---|---|---|---|---|---|
| # contigs (total) | 1,439 | **1,395** | **1,401** | 5/2 | 0.13 | 5/1 | 0.075 |
| NGA50 (median) | 256,222 | 275,040 | 275,040 | **4/0** | — | **3/0** | — |
| NG50 (median) | 275,040 | 306,160 | 306,160 | 4/0 | — | 3/0 | — |
| genome fraction | 98.986 | 98.989 | 98.989 | **6/0** | 0.036 | 4/1 | — |
| duplication | 1.003 | 1.003 | 1.003 | 2/0 | — | 2/0 | — |
| # misassemblies (total) | **95** | **101** | **98** | 0/3 | — | 0/2 | — |
| mismatches /100 kbp | 0.355 | 0.355 | 0.355 | 2/0 | — | 1/1 | — |

**Nine of the sixteen isolates are completely unchanged by either arm.** Everything the gap
closer does, it does on seven, and the whole NGA50 prize sits on three:

| isolate | misasm c/sq/nm | contigs c/sq/nm | NGA50 sq | NGA50 nm |
|---|---|---|---|---|
| GCF046742135v1 | 1/**3**/**3** | 73/69/69 | **+62.85%** | **+62.85%** |
| GCF900607305v1 | 0/0/0 | 91/79/80 | **+51.24%** | **+51.24%** |
| GCF018093065v1 | 0/**2**/0 | 166/139/146 | +13.27% | +12.96% |
| GCF022693245v1 | 5/5/5 | 136/133/135 | +1.03% | +0.00% |
| GCF046742145v1 | 27/**29**/**28** | 630/628/627 | 0 | 0 |
| GCF009912455v1 | 0/0/0 | 26/**29**/**27** | 0 | 0 |
| GCF003945405v1 | 1/1/1 | 86/**87**/86 | 0 | 0 |

**Verdict, and it is narrower than the pilot suggested.** The paired vote keeps all three
large NGA50 wins and costs **+3 misassemblies against sequence-only's +6** — it halves the
price, it does not abolish it. On GCF018093065v1 it removes both added misassemblies, as the
pilot showed. On GCF046742135v1 it removes neither, and that isolate gains **+62.85% NGA50**:
that is the contiguity/correctness trade showing up inside the gap closer, not a defect in
the vote. Note too that closing a graph gap can *raise* the contig count (26 → 27, 86 → 87),
because the join changes what the scaffolder can lay out downstream of it.

**Against SPAdes on the same 16**, with the gap closer on: NGA50 **12/4** (median 275,039 vs
181,475), misassemblies **98 vs 64**. The panel reproduces the whole-cohort picture in
miniature — we win contiguity decisively and lose misassemblies, and the excess is
concentrated where we are far more contiguous (GCF046742145v1: 28 against SPAdes' 1, at
248 kb NGA50 against SPAdes' 67 kb).

**Decision: the vote ships as a gate, not a tie-breaker.** Confirmed at n=40 below, where
the case is far stronger than at n=16.

### The gap closer at n=40 — the vote is load-bearing, not a refinement

| metric | control | seq-only | **seq+vote** | sq w/l | p | nm w/l | p |
|---|---|---|---|---|---|---|---|
| NGA50 | 236,904 | 247,915 | **247,915** | 10/1 | 0.0145 | **9/0** | **0.0092** |
| NG50 | 248,228 | 263,780 | **263,780** | 13/1 | 0.0039 | **11/0** | **0.0039** |
| # contigs | 37.5 | 36.5 | **36.5** | 12/5 | 0.17 | **14/1** | **0.0015** |
| genome fraction | 98.952 | 98.945 | 98.945 | 11/3 | 0.13 | 10/3 | 0.21 |
| duplication | 1.0030 | 1.0030 | 1.0030 | 3/3 | 1.0 | 2/2 | — |
| **# misassemblies (total)** | **198** | **250 (+52)** | **207 (+9)** | 1/6 | 0.035 | 0/5 | — |

**Sixteen of the forty isolates move at all**; NGA50 gains reach +75%, +63%, +51%, +47%,
+31%, +28%.

The vote's value is not the marginal misassembly count — it is that **sequence alone blows
up catastrophically on isolates the vote holds.** GCF010364725v2 goes 4 → **47**
misassemblies without the vote and 4 → 7 with it; that one isolate is most of the +52.
GCF038025155v1 goes 84 → **153** contigs and **−9.87%** NGA50 without the vote, and 84
contigs at **+10.25%** with it. Requiring two read pairs to vouch for the specific pair of
ends is what makes the difference between a gap closer and a chimera generator.

Residual cost of the vote: **+9 misassemblies over 40 isolates, on 5 isolates, none
improved.** Against NGA50 9/0 and contigs 14/1 that is the same contiguity/correctness
trade the rest of the assembler sits on — GCF046742135v1 alone pays +2 for +62.85% NGA50.

**Ceiling, measured before building** (the check skipped before arm E): across 36 graphs,
**median isolate has 1 candidate**, 7 of 36 have ≥10, 3 have ≥50, and one isolate holds
**73% of the entire prize**. This will transform a few isolates and leave most untouched.

---

## 7. Round 2: three divergent agents, four hypotheses, four negatives

Agents were sent to read Shovill/Unicycler/bactopia, SKESA/Velvet, and the
alignment-based gap closers (minimap2, Sealer, GapFiller, Pilon/Racon). Each produced a
candidate mechanism for the contig-count gap (166 vs SPAdes' 110 on `GCF018093065v1`).
**All four were tested and none survived.** Recorded so they are not re-proposed.

| hypothesis | source | verdict |
|---|---|---|
| **Our own cleaning creates the dead ends.** SKESA never deletes anything — all its conservatism is query-time filtering. Our `deleteNode` destroys links, and `removeLocallyWeak`/`removeErroneousConnections` delete *interior* unitigs with no dead-end guard; deleting an interior unitig creates two dead ends. | SKESA/Velvet agent, ranked its own #1 | **FALSIFIED.** Instrumented `totalDeadEnds()` on the raw graph before any cleaning: **51,788 at k=21, falling to 1,338 after twelve rounds.** Simplification removes 97% of dead ends. The survivors are residue, not damage. |
| **Dead ends are single-base bubbles at repeat-copy variant sites**, branches of length `2(k-1)+1` = 253 bp differing at one centre base, declined by the bubble popper because loser depth (26–387) exceeds `meanCoverage × 0.35` = 20.9. | alignment agent, from tracing 7 junctions | **DOES NOT REPRODUCE.** 9 segments of length 253 in the whole graph; **2** two-branch converging bubbles total. The proposed targeted 1-base bubble pop would act on a population of 2. |
| **Short repeat-copy fragments inflate our contig count**; SPAdes deletes every isolated edge up to `max(RL,150)+k` unconditionally while ours is coverage-gated. | SPAdes source review, follow-up | **NO.** SPAdes has fewer contigs at **every** length threshold *and* more total sequence: ≥500 bp 110/2.81 Mb vs our 166/2.75 Mb; ≥5 kb **62/2.74 Mb vs our 98/2.63 Mb**. The 500–2000 bp populations are comparable (37 vs 44). No size filter touches this. |
| **`gfa_connector`-style reconnection** — rebuild at lower k (41) with `min_count 2`, walk 2 kb out of each contig end, keep only paths landing on another contig. | SKESA/Velvet **and** wrapper agents, independently | **Not built.** Ceiling measured first: median isolate **1** closable pair, 7 of 36 with ≥10. Two agents converging made it the most attractive candidate; the ceiling says it would transform ~3 isolates and leave 29 untouched. |

**What did hold**, and is the most useful thing from the round: **45% of dead-end segments
are sequence that already appears verbatim elsewhere in the assembly** (18 of 40 sampled),
median length 170 bp, 71% shorter than 2k. They are not missing sequence — they are
unplaceable repeat copies. That explains why the scaffolder finds nothing (a repeat copy
links everywhere) and why the gap closer found only 32 joins (these are not gaps).

**The residual gap is therefore genuine contiguity in long contigs**, not emission, not
filtering, not dead-end rescue. SPAdes assembles this genome into materially longer pieces
from identical reads.

### Other findings worth keeping

* **Our abundance cutoff implements SKESA's `max(2, depth/50)` but feeds it the k-mer peak
  where SKESA feeds read depth** (`counter.cpp:365`). Ours therefore scales *inversely with
  k*: on a 196x library the first rung showed `cutoff=4 peak=189`, while at k=127 the peak
  falls to ~31 and the scaling never engages. Untested — and note our own measured table
  says a *stricter* top-rung cutoff is worse (cutoff 5 → NGA50 39,855; cutoff 2 → 132,059),
  so the two pieces of evidence disagree and only measurement settles it.
* **Subsampling to ~150x is folklore.** Shovill's default moved 50 → 100 → 150 with no
  stated rationale, and the literature contradicts it for de Bruijn assemblers: GAGE-B chose
  250x *because* 100x was worse; GABenchToB found MiSeq "not susceptible to oversampling";
  Radai et al. (39,000 SPAdes assemblies) found depth positively associated with contiguity.
  SPAdes' own author: *"deep sequencing does not present any problem for assembly algorithms
  relying on the de Bruijn graph per se."*
* **Short-read polishing of short-read assemblies is contraindicated** — Wick 2023 found
  stock Pilon never reduced errors and often added them. Shovill only survives it by
  crippling Pilon (`--fix bases --minmq 60 --mindepth 0.25` plus `samclip`).
* **SKESA clips k bases off both ends of every seed** so contig ends are verified from both
  directions. Our chain ends sit exactly where assembly failed — the least-verified base.
* **Velvet's second-order link derivation** (if A↔B and B↔C, place C relative to A with
  summed variance) reaches ends with **zero** direct links — 85 of our 279. Untested.

### A binary swap mid-campaign silently mixes a cohort

Recorded because it nearly produced a meaningless table and nothing in the output would
have shown it.

AD was made the campaign default at 2026-09-10 18:32, while *A. baumannii* was already
assembling. The monitor then reported `abaumannii COHORT READY: 143/150 scored` — correctly,
by its own definition. But those 143 were **78 isolates assembled by the pre-AD binary and
66 by AD**: two different assemblers averaged into one significance table.

| organism | pre-AD arms | AD arms |
|---|---|---|
| *S. aureus* | 76 | 109 |
| *E. faecium* | 124 | 0 |
| *A. baumannii* | **78** | **66** |

The fix is mtime-gated: rebuild only arms older than the swap. A blanket rebuild wastes half
the work; a blanket skip leaves the cohort mixed. Both are wrong, in opposite directions.

**The general rule this implies:** any binary change during the campaign invalidates every
cohort *currently in flight*, not merely the ones not yet started, and "cohort ready" must
mean "assembled by one binary" rather than "n isolates scored". This will apply again if the
gap-closer work ships.

---

## 5. Campaign state  (2026-09-12 12:46)

| organism | cohort | TesserACT arms on AD | SPAdes benchmark |
|---|---|---|---|
| *S. aureus* | 185 | **185/185 all four arms ✅** | 185/185 ✅ |
| *E. faecium* | 124 | 27/124, running | 124/124 ✅ |
| *A. baumannii* | 150 | 72/150, queued | 150/150 ✅ |
| *P. aeruginosa* | 188 | queued | 68/188 |
| *E. cloacae* | 195 | queued | 0/195 |
| *E. coli* | 200 | queued | 0/200 |
| *K. pneumoniae* | 200 | queued | 0/200 |

*S. aureus* is the first cohort where every arm came from one binary, so its table is clean
by construction rather than by filtering.

All **1,242 libraries stored and `gzip -t` verified**, so everything remaining is
compute-bound. SPAdes arms are never re-run.

**A 22-hour silent stall, worth recording.** The chain that was to run *E. faecium* and
*A. baumannii* after *S. aureus* was `while pgrep -f "rerun_ad.sh saureus"; do sleep 120;
done; ./rerun_all.sh`. That waiter's own command line contains the pattern, so `pgrep -f`
matched **itself** and the loop never exited. *S. aureus* finished at 01:28 and the next two
organisms did not start until the waiter was killed by PID at 08:47. Nothing in any log
showed it and the monitor correctly reported every driver alive — because the waiter *was*
alive. The general rule: **never poll for a previous stage; use `wait` on the PID, or put
both stages in one script.** `rerun_all.sh` already serialises its organisms under one
`flock`, so the outer wait loop was never needed in the first place.

---|---|---|---|
| *S. aureus* | 185 | 124/185 on AD, re-running | 185/185 ✅ |
| *E. faecium* | 124 | 0/124 on AD, queued behind *S. aureus* | 124/124 ✅ |
| *A. baumannii* | 150 | 72/150 on AD, queued | 148/150 |
| *P. aeruginosa* | 188 | queued | 0/188 |
| *E. cloacae* | 195 | queued | 0/195 |
| *E. coli* | 200 | queued | 0/200 |
| *K. pneumoniae* | 200 | queued | 0/200 |

(state at 2026-09-11 13:10)

All **1,242 libraries stored and `gzip -t` verified** on the SSD, so everything remaining
is compute-bound. SPAdes arms are never re-run — they are unaffected by our changes and
cost about as much as all four TesserACT arms combined.

---

## 5b. Power loss, 2026-09-13 18:35 — what survived and what did not

Three reboots (18:38, 18:41, then 06:52 on 09-14) and the machine came back with
**`/media/iowa/u` not mounted**, which is where every FASTQ lives. The disk itself was
healthy — `/dev/sda1`, label `u`, ext4, present in `lsblk` — it simply had no fstab entry
and udisks only auto-mounts on desktop login. `udisksctl mount -b /dev/sda1` restored it at
373 GB used / 518 GB free, identical to before.

Full audit after remount:

| check | result |
|---|---|
| FASTQ libraries | **1,242 / 1,242 present** |
| contigs.fasta truncated (last byte not newline) | **0** |
| QUAST dirs missing `report.tsv` | **0** |
| SPAdes arms with empty contigs | **0** |
| TesserACT arms without `report.json` | **4** (all *A. baumannii*, all in flight at 18:35) |
| leftover `spades_work` | 5 (all *P. aeruginosa*, none had finished) |
| `tesseract-asm` md5 | `ce7e757827ce7198b5820ab878236362`, unchanged |
| `repo/src` fixes + docs | all present |

Damage was confined **entirely to work in flight at the moment of the cut**. The four
part-written arms were deleted so they rebuild, and the five `spades_work` directories
removed. Nothing else needed repair.

Two design choices paid for themselves here. Writing `report.json` **last** makes
"incomplete" trivially detectable — an arm with contigs but no report is a partial write, and
that single test found all four. And never overwriting an existing arm means a restart costs
only what was actually lost.

Scratch binaries under `/tmp` (`tesseract-z`, `tesseract-z2`) were wiped, as `/tmp` is
cleared on boot. They are rebuildable from `repo/src`, which survived.

**Restart order: the whole SPAdes benchmark first, then the TesserACT arms.** SPAdes is
unaffected by our binary, so finishing it first leaves every TesserACT arm to be built by
**one** binary in a single pass — which is the only way a cohort is a valid paired
comparison, and it means a decision to ship the pending fixes no longer risks splitting a
cohort. 641 SPAdes runs remain (*P. aeruginosa* 46, *E. cloacae* 195, *E. coli* 200,
*K. pneumoniae* 200), chained inside one script under one `flock` — never by polling for the
previous stage.

---

## 6. Open questions

1. ~~Does the paired nomination remove the misassembly cost?~~ **Settled at n=40: the vote
   is a gate, not a tie-breaker.** NGA50 9/0 (p=0.0092), NG50 11/0, contigs 14/1 (p=0.0015),
   at +9 misassemblies over 40 isolates. Sequence alone costs +52 and blows up one isolate
   from 4 to 47. Ready to ship; held only until the campaign's cohorts close.
2. ~~Does the ladder union recover anything?~~ **Settled at n=38: 37 of 38 outputs
   byte-identical to control.** The carry-over already does its job.
3. ~~Does the duplication loss have a mechanism?~~ **Settled: two, in equal parts, both
   fixed.** `dedup+trim` takes duplication vs SPAdes from 0/35 (p=2.6e-07) to 9/6 (p=0.887)
   and gives up no win. Ready to ship.
4. ~~Is the misassembly gap fully explained by bought contiguity?~~ **No.** It splits by
   breakpoint flank length into a `>=20 kb` class that is bought contiguity and a `<2 kb`
   class that is not — identical in `base` and `model`, and significant even where SPAdes is
   the more contiguous assembler. **Localised to the resolver's main path (54 of the 66
   excess), and threshold tuning is ruled out as a family.** Next: instrument the joins.
5. **Do the models still work after AD?** **Yes, on *S. aureus*, decisively** — `model` vs
   `base` at n=160 on one binary: NGA50 90/4 (p=2.0e-13), contigs 141/2 (p=6.2e-25), LGA50
   69/3, genome fraction 82/33. No re-refinement needed there. Unanswered for the other six.
6. Does AD hold across all seven organisms, or only where ladders were truncated?
7. Does the *E. faecium* plasmid recovery deficit (86.7% vs SPAdes' 95.0%, whole 8.0% vs
   14.6%) close with the copy-number exemption? Never tested on the plasmid-densest cohort —
   121 of 124 isolates, 625 plasmids. *E. faecium* is re-running on AD now, so this becomes
   answerable shortly.
8. The `model` arm's `_cov_` header field is the scaffold's depth, not the contig's. Affects
   no QUAST metric; affects every user who parses it.
9. `models-v1` has never been released, so `tesseract-get-models` fails for every user today.
