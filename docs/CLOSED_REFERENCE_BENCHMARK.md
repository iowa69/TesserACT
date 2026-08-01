# Closed-reference benchmark

Real Illumina reads assembled against **complete (closed) reference genomes from
the same isolate**, so accuracy can be measured rather than inferred. Seven
isolates, each with the reads and the finished genome deposited under one
BioSample; references are Unicycler hybrid (Illumina + ONT/PacBio) assemblies.
Two of them (INF164, KSB1_7J) come from the Wick/Holt *Klebsiella* set.

Reproduce with `closed/run_closed_benchmark.sh` and `closed/rerun_v3.sh`.

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

**Contiguity: SPAdes wins 6 of 7, one tie.** This is the honest position and it
should not be quoted any other way. An earlier version of this document claimed
a tessera contiguity win; that comparison was scaffolds against contigs and did
not survive being redone properly.

Where tessera is ahead:

- **Base accuracy** — mean 0.53 mismatches/100 kbp against SPAdes' 1.89, better
  on 6 of 7, and 3.6x better overall.
- **Misassemblies** — 8 against 10 in total, with zero local misassemblies
  throughout.
- **Speed and memory** — 2.4-12.6x faster, on roughly half to a third the RAM.

| dataset | tessera | SPAdes |
|---|---|---|
| kpn_RHBSTW00472 | 76s / 2.75 GB | 217s / 4.30 GB |
| kpn_RHBSTW00433 | 71s / 1.77 GB | 198s / 4.25 GB |
| kpn_RHBSTW00128 | 63s / 2.77 GB | 192s / 4.30 GB |
| kpn_1GR13 | 86s / 2.76 GB | 1086s / 5.57 GB |
| eco_RHB07C12 | 70s / 2.75 GB | 197s / 4.28 GB |
| kpn_INF164 | 79s / 1.75 GB | 187s / 4.33 GB |
| kpn_KSB17J | 58s / 2.84 GB | 166s / 4.30 GB |

## What this benchmark found

Two real defects, both fixed, both regression-tested:

**The corrector fabricated sequence.** `kpn_1GR13` is 2x300 MiSeq whose R2
quality collapses to Q3.6 by position 281 -- a ~44% per-base error rate, i.e.
noise. The corrector had no corroboration requirement and a budget of `len/10`
substitutions, so it rewrote those tails into plausible genomic sequence by
walking whichever path kept k-mers solid. Raw noise is harmless (unique k-mers
die at the abundance cutoff); *corrected* noise is not, because it enters the
graph indistinguishable from evidence. Disabling correction entirely scored
better than running it (NGA50 58,161 vs 38,459), which is what exposed it.
Corrections now require the following k-mers to stay solid for four more steps,
and stretches that cannot be vouched for are masked rather than rewritten.

**The abundance cutoff was chosen from a saturation artefact.** The count
histogram folds everything at or above 100,000 into its last bin, so that bin's
index is a floor, not a count. The peak search weighted it by that index anyway.
At small k, where collapsed repeats and low-complexity sequence pile up there,
it outweighed the real coverage mode: `cutoff=787, peak=100000, solid=1,538` out
of 27M distinct 21-mers. Since error correction builds its trusted set at the
smallest k, correction had effectively been running against an empty set on this
dataset. Excluding the saturation bin restores `cutoff=7, peak=56, solid=5.6M`.

Combined effect on the dataset that exposed them, from raw reads:

| | before | after |
|---|---|---|
| NGA50 | 38,459 | 67,886 (+76%) |
| misassemblies | 1 | 0 |
| mismatches/100 kbp | 0.99 | 0.52 |
| largest contig | 166,175 | 292,926 |
| peak RSS | 9.16 GB | 2.76 GB |
| wall clock | 93.5s | 86s |

tessera also now quality-trims 3' ends as reads are loaded (`--no-qtrim`,
`--qtrim-quality`, `--qtrim-window`), before the count that builds the trusted
k-mer set. On binned-quality HiSeq/NovaSeq data this is a no-op; on 2x250/2x300
MiSeq it is what makes the rest of the pipeline behave.

## Known limitation

The remaining contiguity gap is repeat resolution in the graph, not k. Extending
the k ladder does not close it -- on `kpn_1GR13` the ladders `21,33,55`,
`21,33,55,77` and `21,33,55,77,95` give N50 76,933 / 72,168 / 77,131, i.e. it
has plateaued well below the k<=96 ceiling this panel was run under.
`kpn_KSB17J`, the one dataset where SPAdes itself only used k<=55, is also the
one tessera ties.

> **Superseded in part.** An earlier version of this paragraph went further and
> said widening the k-mer representation to reach k=127 "would not help". That
> was a claim about these seven libraries generalised past its evidence. On the
> 212-isolate *Klebsiella* panel (`KLEBSIELLA_PANEL.md`), which is mostly 2x250
> MiSeq, the unitig N50 was still climbing steeply at k=95 and widening to
> k<=128 took it from 142,604 to 207,162 on one isolate. The plateau above is
> real for 2x150 and 2x300 data and says nothing about 2x250. The larger finding
> from that panel -- that the abundance cutoff, not k or repeat resolution, was
> the dominant limiter -- applies here too and postdates every number on this
> page.
