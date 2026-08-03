# Closed-reference benchmark — seven mixed isolates

> **Scope.** Vanilla tessera, **no organism model**, on seven isolates of mixed
> chemistry (2x150, 2x250 and 2x300). Run under an earlier k ceiling of 96.
> This page measures the model-free assembler on heterogeneous libraries; it is
> not the head-to-head that the README quotes, which uses the genus model on a
> *Klebsiella* panel. For that, see [KLEBSIELLA_PANEL.md](KLEBSIELLA_PANEL.md)
> and the table in the README.

Real Illumina reads assembled against **complete (closed) reference genomes from
the same isolate**, so accuracy is measured rather than inferred. Each isolate
has its reads and its finished genome deposited under one BioSample; the
references are Unicycler hybrid (Illumina + ONT/PacBio) assemblies. Two of them
(INF164, KSB1_7J) come from the Wick/Holt *Klebsiella* set.

## How to read this table

SPAdes' `contigs.fasta` is unscaffolded; tessera's is scaffolded. Comparing the
two directly flatters tessera, because N-gap joins inflate NGA50 without
resolving anything in the graph. **Every number below is contig level**:
tessera's scaffolds are split at N-runs >= 10 (QUAST `-s`, the `_broken`
assembly) before scoring. QUAST 5.3.0, default `--min-contig 500`.

## Results

| dataset | tessera NGA50 | SPAdes NGA50 | delta | tessera misasm | SPAdes misasm | tessera mism/100kb | SPAdes mism/100kb |
|---|---|---|---|---|---|---|---|
| kpn_RHBSTW00472 | 81,684 | 107,055 | -23.7% | 1 | 2 | 0.74 | 1.83 |
| kpn_RHBSTW00433 | 58,411 | 87,634 | -33.3% | 3 | 0 | 1.27 | 2.10 |
| kpn_RHBSTW00128 | 127,786 | 173,782 | -26.5% | 0 | 0 | 0.00 | 1.87 |
| kpn_1GR13 | 67,879 | 295,825 | -77.1% | 0 | 4 | 0.52 | 0.69 |
| eco_RHB07C12 | 96,312 | 162,825 | -40.8% | 1 | 1 | 0.55 | 5.31 |
| kpn_INF164 | 63,852 | 117,408 | -45.6% | 3 | 3 | 0.20 | 1.42 |
| kpn_KSB17J | 172,683 | 172,174 | +0.3% | 0 | 0 | 0.43 | 0.02 |

**Without a model, SPAdes leads contiguity on six of these seven, with one tie.**
That is the position on this panel and it should not be quoted any other way.
Contiguity without a model is bounded by what a fragment can span: a 217 bp
insert cannot settle a repeat longer than itself, and tessera stops at those
junctions rather than guessing. Supplying a genus model is what changes this
axis, and is measured separately.

Where tessera leads on this panel:

- **Base accuracy** — mean 0.53 mismatches/100 kbp against SPAdes' 1.89, better
  on 6 of 7, 3.6x better overall.
- **Misassemblies** — 8 against 10 in total, with zero local misassemblies
  throughout.
- **Speed and memory** — 2.4-12.6x faster on this panel, on roughly half to a
  third of the RAM.

| dataset | tessera | SPAdes |
|---|---|---|
| kpn_RHBSTW00472 | 76s / 2.75 GB | 217s / 4.30 GB |
| kpn_RHBSTW00433 | 71s / 1.77 GB | 198s / 4.25 GB |
| kpn_RHBSTW00128 | 63s / 2.77 GB | 192s / 4.30 GB |
| kpn_1GR13 | 86s / 2.76 GB | 1086s / 5.57 GB |
| eco_RHB07C12 | 70s / 2.75 GB | 197s / 4.28 GB |
| kpn_INF164 | 79s / 1.75 GB | 187s / 4.33 GB |
| kpn_KSB17J | 58s / 2.84 GB | 166s / 4.30 GB |

## Two defaults this panel set

**Error correction requires corroboration.** `kpn_1GR13` is 2x300 MiSeq whose R2
quality collapses to Q3.6 by position 281 — a ~44% per-base error rate, which is
noise. Raw noise is harmless, because unique k-mers die at the abundance cutoff.
*Corrected* noise is not: it enters the graph indistinguishable from evidence.
A corrector free to spend a budget of `len/10` substitutions, walking whichever
path keeps k-mers solid, will rewrite such a tail into plausible genomic
sequence. tessera therefore requires the k-mers following a correction to stay
solid for four more steps, and masks stretches it cannot vouch for rather than
rewriting them. On this dataset the difference is NGA50 38,459 against 67,886.

**The abundance cutoff excludes the saturation bin.** The count histogram folds
everything at or above 100,000 into its last bin, so that bin's index is a floor
rather than a count. Weighting it by that index lets it outvote the real
coverage mode at small k, where collapsed repeats and low-complexity sequence
pile up there: on this dataset that gives `cutoff=787, peak=100000, solid=1,538`
out of 27M distinct 21-mers. Since error correction builds its trusted set at
the smallest k, a cutoff chosen that way leaves correction with nothing to work
from. Excluding the saturation bin gives `cutoff=7, peak=56, solid=5.6M`.

Both defaults together, on the dataset that is most sensitive to them, from raw
reads:

| | histogram valley, uncorroborated correction | shipped defaults |
|---|---|---|
| NGA50 | 38,459 | 67,886 (+76%) |
| misassemblies | 1 | 0 |
| mismatches/100 kbp | 0.99 | 0.52 |
| largest contig | 166,175 | 292,926 |
| peak RSS | 9.16 GB | 2.76 GB |
| wall clock | 93.5s | 86s |

tessera also quality-trims 3' ends as reads are loaded (`--no-qtrim`,
`--qtrim-quality`, `--qtrim-window`), before the count that builds the trusted
k-mer set. On binned-quality HiSeq/NovaSeq data this is a no-op; on 2x250 and
2x300 MiSeq it is what makes the rest of the pipeline behave.

## What the k ladder does and does not reach

On this panel the ladder plateaus below the ceiling it was run under: on
`kpn_1GR13` the ladders `21,33,55`, `21,33,55,77` and `21,33,55,77,95` give N50
76,933 / 72,168 / 77,131. `kpn_KSB17J`, the one dataset where SPAdes itself only
used k<=55, is also the one tessera ties.

That plateau is a property of 2x150 and 2x300 libraries and does not generalise
to 2x250, where a read carries far more k-mer. On the *Klebsiella* panel, which
is mostly 2x250, unitig N50 was still climbing steeply at k=95 and raising the
top rung to 127 took it from 142,604 to 207,162 on one isolate. This is why the
ladder's top rung tracks read length rather than stopping at 95 — see
[KLEBSIELLA_PANEL.md](KLEBSIELLA_PANEL.md).
