# Where SPAdes still beats TesserACT, and what it is not

> **Scope.** Model-free and model-guided TesserACT against vanilla SPAdes 4.3.0 on the
> same reads, 47 non-clonal *S. aureus* isolates each with its own closed reference,
> QUAST 5.3.0, contig level, zero-N output on both sides. SPAdes wins contig NGA50 on
> **14 of 47 (30%)**. This page is about those 14.

Most of what follows is elimination. That is the useful part: it says where not to look.

## What it is not

Each of these was tested directly, on the isolates where the gap is largest, and each
is refuted.

| hypothesis | test | result |
|---|---|---|
| Repeat content / IS elements | contig ends mapped back to the reference | breaks sit on **AT-rich coverage dropouts**, not repeats: 56-69% of break windows are >9 GC points below the genome mean, against 5% of random windows |
| The k ladder stops too low | per-rung `n50_final` from `report.json` | the ladder reaches k=127 and still loses; on the isolate we win it reaches the same ceiling and gets 129,919 |
| Paired-end resolution is too strict | break sites vs SPAdes contigs | **76% of break sites fall inside continuous SPAdes sequence** built from the same reads. The reads carry it; the graph never got it. Nothing for the resolver to decline |
| Read correction masks too much | `--no-correct` on the three worst | NGA50 moves **<0.5%** (17,440→17,532; 35,068→35,068) |
| Library contamination | solid k-mers / genome length | SPAdes wins 38% of libraries above 1.3x excess and 26% below. Weak, and the highest-excess library (2.31x) is one we win |
| Merging overlapping mates | `fastp --merge`, 89% and 97% merged | **worse on both**: 17,440→14,825 and 35,068→20,965. Confirms the rejection recorded in KLEBSIELLA_PANEL.md, now on the library type that most favours merging |

## What it is

The only covariate that separates the two groups is **insert size relative to read
length**: median 1.364 where SPAdes wins against 1.513 where it does not. At 1.14 the
mates overlap by roughly 85%, so a pair carries barely more information than one read.

And the sequence SPAdes recovers is **sub-threshold k-mers**. Running with the
abundance cutoff forced to 1:

| isolate | `-c 1` | default | SPAdes | genome fraction, default |
|---|---|---|---|---|
| GCF010364725v2 | 2,874 | **5,781** | 14,366 | 94.9% |
| GCF003111725v1 | 7,044 | **17,440** | 28,223 | 95.6% |
| GCF016591995v1 | **51,105** | 35,068 | **51,105** | 98.6% |

On the third isolate `-c 1` reaches NGA50 **51,105 -- exactly SPAdes' figure**, +46% over
the default, with better genome fraction (98.83 against 98.56) and fewer contigs (165
against 214). The missing sequence is there, in k-mers seen once, and it is recoverable.

On the other two the same setting is a disaster: genome fraction falls to 89.3%. The
split is by library quality -- `-c 1` helps the clean library and wrecks the noisy ones.
On a clean library a count-1 k-mer is real low-coverage sequence; on a noisy one it is
an error.

## What would fix it, and why it is not done here

Not a lower cutoff. A rescue that admits a sub-threshold k-mer only when it **bridges**
two solid ones, so the decision is made on graph topology rather than on a count.

That is harder than it sounds and is the reason it is not in this release. A single
substitution error produces a chain of k count-1 k-mers, and so does a genuine
low-coverage stretch; in both cases the interior k-mers' immediate neighbours are also
sub-threshold. A rule that requires both neighbours solid therefore rescues almost
nothing, and a rule that accepts short sub-threshold chains between solid regions is a
graph operation, not a counting one -- and is exactly where error k-mers return.

The evidence here is three isolates. `counter.cpp` already carries a measured table for
this threshold (cutoff 5 / 3 / 2 giving NGA50 39,855 / 86,779 / 132,059), and changing a
default that was set that carefully needs more than three isolates against it. The
larger cohort now running is what should settle it.

## The honest summary

TesserACT wins contig NGA50 on 70% of these isolates and loses on 30%. The losses are
concentrated on libraries whose fragments barely exceed their reads, and the mechanism
is low-coverage sequence discarded at the abundance threshold rather than anything about
repeats, k, pairing, correction or contamination. Every one of those was tested and
none of them is the cause.
