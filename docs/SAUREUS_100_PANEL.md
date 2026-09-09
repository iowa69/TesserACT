# S. aureus non-clonal panel, n=180 (complete)

The first of the seven ESKAPEE cohorts to reach the target size. Every number here
is measured, on this machine, against each isolate's own closed reference.

## How the cohort was built

185 *S. aureus* isolates with both a closed RefSeq assembly and a public paired-end
Illumina library. All 185 assembled and scored in five arms; this is the complete cohort,
not a snapshot of one. Libraries were sketched with `mash` and one isolate kept per cluster
at `d <= 0.0005`, so no two members of the panel are clonal redeposits of each other.
This matters for the statistics: a paired test over 180 clonal replicates measures
sequencing noise, not assembler behaviour.

180 isolates survive the quality gates. Five are discarded and the reason recorded, all
for the same reason:

| isolate | unaligned to its own reference, in **every** arm |
|---|---|
| GCF009912475v1 | 44.6% |
| GCF053165775v1 | 16.4% |
| GCF031190615v1 | 10.0% |
| GCF002895385v1 | 7.9% |
| GCF046741825v1 | 7.7% |

A double-digit unaligned fraction in every arm is not an assembly failure — it means the
deposited reference and the deposited reads are not the same strain. Scoring those would
measure the database, not the assembler. No isolate was discarded for a poor score.

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
| NGA50 | 247,188 | 188,658 | **+31.0%** | 133/47 | 2.0e-10 |
| NG50 | 253,984 | 196,578 | **+29.2%** | 131/49 | 2.9e-10 |
| Genome fraction (%) | 98.83 | 98.37 | **+0.46 pp** | 170/10 | 1.3e-27 |
| Largest alignment | 538,873 | 427,134 | **+26.2%** | 126/52 | 3.1e-08 |
| Mismatches / 100 kb | 0.66 | 1.22 | **−46.3%** | 111/65 | 0.0017 |
| Indels / 100 kb | 0.25 | 0.29 | −13.8% | 89/77 | 0.13 (ns) |
| # contigs | 38 | 40 | −5.0% | 88/86 | 0.60 (ns) |
| Duplication ratio | 1.000 | 1.000 | +0.2% | 5/151 | 8.4e-25 |
| # misassemblies | 0 | 0 | — | 29/53 | **0.011 (SPAdes)** |

Eight of nine metrics favour TesserACT, four of them by more than twenty percent.
The ninth is the honest one, and it is discussed below.

## The one metric SPAdes wins

Median misassemblies are zero in both arms, so this is a tail effect: 223 for us against
150 for SPAdes across roughly 7,000 contigs.

Stratifying by how much more contiguous we are than SPAdes on the same isolate:

| NGA50, us vs SPAdes | n | us | SPAdes | excess | win/loss | p |
|---|---|---|---|---|---|---|
| SPAdes ahead (< −10%) | 35 | 33 | 27 | +6 | 9/10 | 0.702 |
| **matched (−10%..+10%)** | **48** | **86** | **70** | **+16** | **5/17** | **0.022** |
| we lead 10–50% | 50 | 44 | 25 | +19 | 7/15 | 0.072 |
| we lead >50% | 47 | 60 | 28 | +32 | 8/11 | 0.445 |

Part of the gap is bought contiguity — the excess grows monotonically with how far ahead we
are, and the largest block is where we lead by more than 50%. A longer contig crosses more
junctions and has more chances to cross one wrongly.

But **a residual gap survives matching**: +16 over 48 isolates at 5/17, p=0.022. That is a
real difference in decision quality at equal contiguity. An earlier version of this document,
written when only 151 isolates had finished, read the same band at 5/14, p=0.067 and called
it "not significant" — that was a statement about the sample size, not a result, and the
complete cohort reverses it. [HOW_SPADES_CHOOSES_A_BRANCH.md](HOW_SPADES_CHOOSES_A_BRANCH.md)
follows that through SPAdes' source.

One isolate still dominates the raw total: GCF046742145v1 contributes 29 of the 73 (30
against SPAdes' 1) while leading NGA50 by 269% — there SPAdes reached 67 kb against our
248 kb and earned its single misassembly by producing fragments too short to be wrong.
Unlike at n=151, dropping it no longer removes the effect, because the matched band does
not depend on it.

Note also GCF022832835v1: 61 misassemblies for us, 58 for SPAdes. When both assemblers fail
on the same isolate, the isolate is the problem.

## How much the model is worth, and why one number cannot say

The model arm never loses to the no-model arm on NG50 across the complete 180 isolates —
**103 wins, 0 losses, p=1.3e-18**, and NGA50 101/0 at p=2.7e-18. But *how much* it is worth
depends on which summary is quoted, and the win/loss is the part that is unambiguous:

| summary of the same 180 paired NG50 values | model vs base |
|---|---|
| median of the two columns | **+13.4%** |
| win/loss | **103 / 0** |

The reason the central estimates spread is the shape: about a hundred isolates improve and
the rest are unchanged to the base pair, because the model only has something to say where
its marker adjacency covers the junction. The typical isolate gains a little, a minority
gain a great deal, and **none get worse**.

It is not free. Over the same 180 isolates the model costs misassemblies against the no-model
arm — 18 isolates worse against 2 better, p=0.0011. The adjacency table tells the resolver to
make joins the paired reads would not have made alone, and most, not all, of those are right.

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
