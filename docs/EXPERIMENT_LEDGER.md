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

## 4. In flight

**Paired-read gap closer.** SPAdes' `late_gapcloser`, reimplemented as its own stage after
simplification. Two independent tests must both pass: an exact overlap of [10, k−2] bases
(the missing bases are already carried in the two flanks' terminal overhangs — nothing is
invented), **and** ≥2 read pairs vouching for that specific pair of ends, gathered through
a tip-neighbourhood map that walks back 5 kb along non-branching paths.

The false blocker is gone: `graph.cpp:917` claimed shorter overlaps require *"giving Link
its own overlap field first"*. They do not — SPAdes inserts an **edge carrying the missing
bases**, so every link still splices exactly k−1 and `mergeInto` is untouched.

Sequence-only result (nomination pending): **32 bridged, contigs 166 → 139, NGA50 +13.3%,
genome fraction +0.13 pp, duplication 1.001 → 1.000 — at a cost of 0 → 2 misassemblies.**
The nomination is expected to remove that cost.

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

## 5. Campaign state

| organism | cohort | TesserACT arms | SPAdes benchmark |
|---|---|---|---|
| *S. aureus* | 185 | re-running on AD | 185/185 ✅ |
| *E. faecium* | 124 | pre-AD, needs re-run | 124/124 ✅ |
| *A. baumannii* | 150 | in progress | ~127/150 |
| *P. aeruginosa* | 188 | queued | queued |
| *E. cloacae* | 195 | queued | queued |
| *E. coli* | 200 | queued | queued |
| *K. pneumoniae* | 200 | queued | queued |

All **1,242 libraries stored and `gzip -t` verified** on the SSD, so everything remaining
is compute-bound. SPAdes arms are never re-run — they are unaffected by our changes and
cost about as much as all four TesserACT arms combined.

---

## 6. Open questions

1. Does the paired nomination remove the 2 misassemblies while keeping the 27 contigs?
2. Does AD hold across all seven organisms, or only where ladders were truncated?
3. Does the *E. faecium* plasmid recovery deficit (86.7% vs SPAdes' 95.0%, whole 8.0% vs
   14.6%) close with the copy-number exemption? It has never been tested on the
   plasmid-densest cohort — 121 of 124 isolates, 625 plasmids.
4. The residual misassembly gap at matched contiguity (p=0.022 at n=180) still has no
   identified mechanism.
