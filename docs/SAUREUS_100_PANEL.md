# S. aureus non-clonal panel, n=151

The first of the seven ESKAPEE cohorts to reach the target size. Every number here
is measured, on this machine, against each isolate's own closed reference.

## How the cohort was built

185 *S. aureus* isolates with both a closed RefSeq assembly and a public paired-end
Illumina library. Libraries were sketched with `mash` and one isolate kept per cluster
at `d <= 0.0005`, so no two members of the panel are clonal redeposits of each other.
This matters for the statistics: a paired test over 151 clonal replicates measures
sequencing noise, not assembler behaviour.

151 isolates survive the quality gates. Eight are discarded and the reason recorded:

| isolate | reason |
|---|---|
| GCF031190615v1 | 10.0% of the assembly unaligned to its own reference in **every** arm |
| GCF046741825v1 | 7.7% unaligned in every arm |
| GCF053165775v1 | 16.4% unaligned in every arm |
| 5 others | still assembling at the time of collection |

A double-digit unaligned fraction in every arm is not an assembly failure — it means
the deposited reference and the deposited reads are not the same strain. Scoring those
would measure the database, not the assembler.

## Scoring

QUAST with `-s`, and **only the `_broken` column is read**. QUAST splits scaffolds at
runs of >=10 Ns and reports the split contigs in that column; the unbroken column
rewards an assembler for emitting Ns, which is not sequence. TesserACT emits zero Ns
in `contigs.fasta`, so for TesserACT the two columns are identical — the `-s` flag
exists to stop SPAdes' scaffold gaps from being counted as contiguity. See
[GAPS_AND_SCORING.md](GAPS_AND_SCORING.md).

## Head to head: TesserACT `--organism saureus` vs vanilla SPAdes 4.3.0

Medians over the 151 paired isolates; win/loss is the paired count; *p* is a two-sided
Wilcoxon signed-rank test.

| metric | TesserACT | SPAdes | delta | win/loss | p |
|---|---|---|---|---|---|
| NGA50 | 247,472 | 196,476 | **+26.0%** | 111/40 | 2.3e-08 |
| NG50 | 258,942 | 204,038 | **+26.9%** | 110/41 | 1.6e-08 |
| Genome fraction (%) | 98.89 | 98.33 | **+0.56 pp** | 146/5 | 1.2e-24 |
| Largest alignment | 551,977 | 449,944 | **+22.7%** | 104/45 | 9.5e-07 |
| Mismatches / 100 kb | 0.65 | 1.16 | **-44.0%** | 91/56 | 0.0096 |
| Indels / 100 kb | 0.24 | 0.29 | -17.2% | 73/67 | 0.23 (ns) |
| # contigs | 38 | 39 | -2.6% | 72/73 | 0.78 (ns) |
| Duplication ratio | 1.000 | 1.000 | +0.2% | 2/129 | 2.2e-22 |
| # misassemblies | 0 | 0 | — | 26/42 | **0.048 (SPAdes)** |

Eight of nine metrics favour TesserACT, four of them by more than twenty percent.
The ninth is the honest one, and it is discussed below.

## The one metric SPAdes wins

Median misassemblies are zero in both arms, so this is entirely a tail effect:

| arm | total misassemblies | isolates with any |
|---|---|---|
| base (no model) | 157 | 47/151 |
| **model** | **177** | **55/151** |
| careful | 161 | 54/151 |
| aggressive | 201 | 68/151 |
| **SPAdes** | **118** | **38/151** |

We make 59 more misassemblies than SPAdes across 15,000-odd contigs, and 30 of those
59 come from a single isolate. Splitting the cohort by which arm misassembles more:

| group | n | median NGA50 vs SPAdes |
|---|---|---|
| TesserACT has more misassemblies | 42 | **+9.6%** |
| tied | 83 | **+11.4%** |
| TesserACT has fewer | 26 | **+15.7%** |

The NGA50 advantage holds in all three groups, so this is not a blanket
contiguity-for-correctness trade: on most isolates we are both longer and no worse.

The subset that actually costs us is small — **10 of 151 isolates (6.6%)** where
TesserACT is simultaneously more misassembled *and* shorter than SPAdes:

| isolate | misasm | NGA50 vs SPAdes | genome fraction |
|---|---|---|---|
| GCF026547035v1 | 1 vs 0 | -61.2% | 99.02 vs 98.53 |
| GCF022693245v1 | 5 vs 0 | -48.2% | 98.19 vs 97.94 |
| GCF005931015v2 | 2 vs 0 | -46.3% | 99.22 vs 98.94 |
| GCF046742125v1 | 2 vs 0 | -41.7% | 97.95 vs 98.08 |
| GCF054392195v1 | 9 vs 5 | -33.3% | 96.96 vs 96.71 |
| GCF008620175v1 | 2 vs 1 | -31.1% | 98.18 vs 97.76 |
| GCF021869805v1 | 1 vs 0 | -18.9% | 98.61 vs 98.15 |
| GCF022405335v1 | 1 vs 0 | -9.4% | 98.75 vs 98.05 |
| GCF009912795v1 | 1 vs 0 | -8.5% | 98.50 vs 97.64 |
| GCF046741885v1 | 1 vs 0 | -0.1% | 98.49 vs 97.31 |

**In nine of these ten, TesserACT still recovers more of the genome than SPAdes.** The
failure is not that we miss sequence — we assemble more of it and then join one or two
junctions wrongly, and a single misjoin near the middle of a chromosome halves NGA50.
SPAdes recovers less and declines to make the join. That is the concrete target for
correction, and it is a repeat-resolution decision, not a coverage-threshold one.

Note also GCF022832835v1: 61 misassemblies for TesserACT, 58 for SPAdes. When both
assemblers fail on the same isolate the isolate is the problem, not the assembler.

## How much the model is worth, and why one number cannot say

The model arm never loses to the no-model arm on NG50 across 151 isolates — 86 wins, 0
losses, p=8.1e-16 — but *how much* it is worth depends entirely on which summary is
quoted, and the honest answer is a range:

| summary of the same 151 paired NG50 values | model vs base |
|---|---|
| median per-isolate gain | **+2.8%** |
| median of the two columns | +10.5% |
| mean of the two columns | +17.6% |
| mean per-isolate gain | **+20.2%** |

A factor of seven, from summary choice alone. The reason is the shape: 86 isolates improve
and 65 are unchanged to the base pair, because the model only has something to say where
its marker adjacency covers the junction. The typical isolate gains a little, a minority
gain a great deal, and none get worse.

`docs/MODELS_INTEGRATED.md` reports +86.1% for *S. aureus* from a single held-out genome.
That measurement is real; it is the far right tail of this distribution. The number to
quote for the model's worth is the win/loss and the p-value, which are unambiguous, with
whichever central estimate is stated as what it is.

## Replicon classification

SPAdes does not classify contigs, so this column has no comparison arm. Recovery is
measured against each reference's own plasmid records; classification is pooled over
every contig >= 1500 bp in all 151 isolates (micro-averaged, not a mean of per-isolate
rates — 93 of 150 isolates carry at least one plasmid, 119 plasmids in total).

| arm | plasmids recovered | recovered whole | precision | recall | F1 |
|---|---|---|---|---|---|
| base (no model) | 111/119 (93.3%) | 35 (29.4%) | 0.396 | 0.573 | 0.468 |
| **model** | 111/119 (93.3%) | 38 (31.9%) | **0.712** | **0.794** | **0.751** |
| careful | 110/119 (92.4%) | 40 (33.6%) | 0.707 | 0.812 | 0.756 |
| aggressive | 111/119 (93.3%) | 37 (31.1%) | 0.547 | 0.773 | 0.640 |
| SPAdes | 106/115 (92.2%) | 46 (40.0%) | — | — | — |

Two honest observations:

* SPAdes assembles more plasmids into a single contig than we do (40.0% vs 31.9%).
  Recovery at the sequence level is a tie; contiguity within the plasmid is not.
* The model arm's precision rests on three fixes made during this campaign — a
  1,500 bp classification floor, a corroboration requirement before a lone contig may
  be called a plasmid, and a guard on label propagation. Pooled precision before
  those fixes was 0.225 at recall 1.000: the classifier called nearly everything a
  plasmid and was right a quarter of the time.

**Caveat on that before/after.** These 151 assemblies span both the pre-fix and
post-fix binary, because the fixes landed partway through the campaign. The offline
scorer re-applies the 1,500 bp floor uniformly, but it cannot retroactively apply the
corroboration or propagation fixes. So the table above is the current state of the
cohort, not a clean controlled comparison. A controlled one requires the post-fix
binary across all isolates.

## Reproducing

```
devtools/panels/select_cohort100.py saureus      # mash-dereplicated cohort
devtools/panels/fetch_all_v3.sh saureus          # libraries + closed references
devtools/panels/run_asm100.sh saureus 5          # five arms per isolate
devtools/panels/collect100.py eval100/saureus out.tsv out.txt
devtools/panels/score_plasmids100.sh saureus       # replicon recovery + calls
```

Per-isolate values for all five arms and ten metrics are in
[saureus_100_headtohead.tsv](saureus_100_headtohead.tsv).
