# Every variant tested, with its parameters and result

One row per change that has ever been built and measured. `n` is the number of isolates it
was measured on, always paired and always scored by QUAST against **each isolate's own closed
reference genome**. Win/loss counts are isolates, not metrics. Incomplete rows are marked as
such rather than omitted.

Status key — **SHIPPED**: on by default in the current binary. **DEFAULT-ON (v2)**: on by
default in the candidate binary, validated but not yet in a full campaign. **OFF**: built and
gated, not default. **REJECTED**: measured and abandoned. **INERT**: works, changes nothing.

---

## 1. k-ladder and graph construction

| # | variant | switch / parameter | value | n | result | status |
|---|---|---|---|---|---|---|
| A | k ladder from read-length **distribution** instead of mean-binned cliff | `TESSERACT_MIN_KMER_MASS` | 0.35 | 16 | NGA50 **+26%** on truncated ladders, 13 better / 0 worse, p=0.002 | **SHIPPED** |
| D | tip-clipper **absolute** branch (SPAdes `tc_lb 10`) | `TESSERACT_TIP_ABS_MULT` | 10.0 | 38 | genome fraction 18/6 (p=0.009), contigs −5.8% (16/7, p=0.015) | **SHIPPED** |
| AD | A + D together | both | — | 38 | NGA50 **+6.0%** (13/4, p=0.037), NG50 +8.8%, genome fraction 22/8 (**p=0.0021**), contigs 17/9, **LGA50 9/0 (p=0.0092)**, misassemblies 3/3, duplication 6/4 | **SHIPPED** |
| B | trusted carry-over (SPAdes `TrustedContigs`) | `TESSERACT_TRUSTED_CARRY` | 1 | 38 | NGA50 **4/11**, genome fraction 10/22, largest aln 2/7 (p=0.044) | **REJECTED** |
| E | per-rung EC threshold fitted from that rung's histogram | `TESSERACT_FITTED_EC` | 1 | 38 | NGA50 **6/23** (p=6.9e-04), contigs 3/30 (p=1.5e-04), mismatches **+103%** | **REJECTED** |
| C | ladder union — emit what a smaller k resolved and the final graph lost | `--ladder-union` | off | 38 | **37/38 outputs byte-identical.** One isolate: +2 contigs, +0.167 pp genome fraction, NGA50 unchanged | **INERT**, kept behind flag |
| — | Unicycler-style best-rung selection | — | — | 5 | **0 better / 5 worse** | **REJECTED** |
| — | naive ladder union (concatenate rungs) | — | — | 3 | duplication **2.99**, 7.9 Mb for a 2.8 Mb genome | **REJECTED** |
| — | shorter k ladder (cap 77) | — | — | 3 | +2% on one, **−4% and −10%** on the others | **REJECTED** |
| — | ABC combined | — | — | 38 | NGA50 10/11, p=0.835 — A's gain exactly cancelled by B's harm | **REJECTED** |

## 2. Output emission — the duplication fix

| # | variant | switch | value | n | result | status |
|---|---|---|---|---|---|---|
| 7 | drop contigs that are exact substrings of a longer contig | `--dedup-contained` | **on** | 44 + 124 | 6.2% of our contigs are exact substrings; SPAdes emits **0 of 1554** | **DEFAULT-ON (v2)** |
| 8 | trim exact terminal dovetails from the shorter partner | `--trim-overlap` | **final k** | 44 + 124 | 90.6% of our duplicated bases were in overlaps **>2k**; SPAdes puts 70% of its in a single overlap of exactly K | **DEFAULT-ON (v2)** |
| 7+8 | both, *model* arm, S. aureus | — | — | 44 | duplication vs SPAdes **0/35 (p=2.6e-07) → 9/6 (p=0.887)**; every other metric unchanged to the decimal | **DEFAULT-ON** |
| 7+8 | both, *model* arm, E. faecium | — | — | 124 | duplication Cliff **+0.97 → +0.05**; contigs win strengthens 98/24 → **113/11** | **DEFAULT-ON** |
| 7+8 | both, **vanilla** arm | — | — | 185 + 124 | **Contigs tie → WIN** (77/104→112/70; 63/60→105/17), **Indel tie → WIN** (efm), duplication 4/155→28/42 and 0/124→32/46, **medians now identical to SPAdes** | **DEFAULT-ON** |
| 8b | trim threshold sweep | `--trim-overlap` | 300/600/1200/2400 | 20 | **no sweet spot** — exchange rate flat at 2.4–2.5 duplicated bases removed per base of unique coverage lost | use final k |

## 3. Gap closing — the strongest unshipped change

| # | variant | switch | value | n | result | status |
|---|---|---|---|---|---|---|
| 9 | gap close on the **raw** graph, sequence only | — | — | 1 | **37,685** bridges against a ceiling of 804 | **REJECTED** |
| 10 | gap close after simplification, **sequence only** | `TESSERACT_GAPCLOSE=10`, `VOTES=0` | — | 40 | NGA50 10/1 (p=0.015), NG50 13/1, contigs 12/5, **misassemblies 198 → 250 (+52)**; one isolate 4 → **47** | **REJECTED** |
| 11 | **gap close + paired-read nomination** | `TESSERACT_GAPCLOSE=10`, `VOTES=2` | 2 votes, 5 kb tip neighbourhood, 31-mer index | 40 | **NGA50 9/0 (p=0.0092)**, **NG50 11/0 (p=0.0039)**, **contigs 14/1 (p=0.0015)**, genome fraction 10/3, duplication 2/2, **misassemblies +9 over 40 isolates** | **OFF — best unshipped candidate** |

Gains on individual isolates reach **+75%, +63%, +51%, +47%, +31%, +28% NGA50**; 16 of 40
isolates move at all. The paired vote is load-bearing, not a refinement: without it one
isolate goes 4 → 47 misassemblies and another loses 9.87% NGA50 while gaining 69 contigs.

## 4. Repeat resolver — the misassembly blocker (five consecutive negatives)

| # | variant | switch | value | n | result | status |
|---|---|---|---|---|---|---|
| 12 | withdraw all three fallbacks at short destinations | `TESSERACT_MIN_FALLBACK_DEST` | 2000 | 24 | `<2 kb` class **166 → 164**; contigs 1/9 worse (p=0.0125), one isolate −10.4% NGA50 | **REJECTED** |
| 13 | raise paired-support bar | `--min-link` | 8 / 12 / 16 | 14 | `<2 kb` class **154 → 153 → 153 → 155**; NGA50 −7.3%, contigs 32 → 36 | **REJECTED** |
| 14 | raise tie-discrimination ratio | `--tie-ratio` | 1.40 (vs 1.02) | 185 | `<2 kb` class **218 → 212**, p=0.322 | **REJECTED** |
| 15 | refuse fallbacks on routes no read pair can span | `TESSERACT_NO_UNSPANNED_FALLBACK` | 1 | 20 | fires 12–36×/isolate, **output identical on 18 of 20** | **REJECTED** |
| 16 | remove `pickByCoverage` entirely | `TESSERACT_NO_DEPTH_PICK` | 1 | 23 | every metric **identical to the base pair** | **INERT** (simplification, not optimisation) |
| — | scaffolder symmetric-target fix | — | — | 1 | joins **5 → 5**; terminus restriction is load-bearing | **REJECTED**, reverted 96 lines |
| — | disable gap filling | `--no-gapfill` | — | 10 | misassemblies −14% on the worst isolates, contigs **10/0 worse**, NGA50 6/0 worse | **REJECTED** |

**What the five negatives jointly establish:** the wrong joins are *strongly supported* (a bar
of 16 does not touch them) and *not near-ties* (1.40 does not touch them), so the wrong
candidate wins **decisively**. The breakpoint blocks are dispersed repeats with **5–7 copies**
in the reference sitting at **6–10× median depth** — far above the 1.6× repeat guard, so they
do not escape it by dilution. Stage ablation: our raw graph gives 100 such events vs SPAdes'
88, and the **resolver adds 54 of the 66 excess**.

## 5. Read handling

| # | variant | n | result | status |
|---|---|---|---|---|
| — | adapter read-through trimming | — | condition is **3.3%** of pairs; 69.2% overlap perfectly but only 3.3% have a fragment shorter than R1 | **REJECTED** |
| — | mate merging (fastplus/scepter) | 3 | +2–5% on two isolates, genome fraction **94.94 → 81.07** on the third | **REJECTED** |
| — | BayesHammer-corrected reads as input | — | closes 13% / 0.5% of the gap | **the corrector is not the gap** |

## 6. Run modes — parameter sets, not code changes

| mode | simplifyRounds | minLinkSupport | linkSupportPerX | tieRatio | bubbleCoverageLimit | polishPasses |
|---|---|---|---|---|---|---|
| fast | 6 | 3 | 0.02 | 1.30 | default | 0 |
| **standard (vanilla)** | **12** | **2** | **0.02** | **1.02** | default | **1** |
| careful | 24 | 3 | 0.14 | 1.40 | default | 2 |
| aggressive | 16 | 2 | 0.05 | 1.05 | 10.0 | 1 |

Losses against SPAdes out of 30 cells (3 organisms × 10 metrics), before the duplication fix:

| mode | losses | which |
|---|---|---|
| careful | **6** | duplication ×3, misassemblies ×3 |
| model | 6 | duplication ×3, misassemblies ×3 |
| base | 8 | + LGA50 and contigs on *A. baumannii* |
| aggressive | 10 | + mismatches ×3, indels ×1 |

Misassembly losses worsen monotonically with contiguity — base 26/51, careful 26/52,
model 25/62, aggressive 26/84 — so no point on this curve reaches zero losses. The mode is
not the lever.

## 7. Replicon / plasmid

| variant | n | result | status |
|---|---|---|---|
| plasmid copy-number exemption | 39 | whole-plasmid recovery **33% → 44%**, low-copy control band unmoved | **SHIPPED** (+1 misassembly total) |
| 1,500 bp replicon floor + corroboration + propagation guard | 185 | pooled precision **0.225 → 0.712**, recall 1.000 → 0.794 | **SHIPPED** |

---

## Ranking of what is not yet default

1. **Gap closer + paired nomination (#11)** — the only unshipped change with a significant
   multi-metric win: NGA50 9/0, NG50 11/0, contigs 14/1. Cost: +9 misassemblies over 40
   isolates. Measured on *S. aureus* only; needs a second organism before it ships.
2. **dedup + trim (#7+8)** — already default in v2, proven on two organisms and on both the
   model and vanilla arms.
3. Everything else is shipped, rejected, or inert.
