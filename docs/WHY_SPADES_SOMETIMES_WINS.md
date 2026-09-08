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

And the sequence SPAdes recovers is removed by TesserACT **twice**, which is why it took
three experiments to see it.

TesserACT has two gates that discard the same low-coverage bases:

1. the read corrector masks stretches it cannot vouch for -- 923,741 bases (0.57%) on a
   clean library, **16,631,849 bases (10.4%)** on a noisy one
2. the abundance cutoff then discards k-mers seen fewer than `cutoff` times

Either gate alone removes the sequence, so opening either alone changes nothing. Measured
on GCF016591995v1, contig NGA50:

| configuration | NGA50 | contigs | genome fraction | mismatches/100kb |
|---|---|---|---|---|
| default (correction on, cutoff auto) | 35,068 | 214 | 98.561 | 0.46 |
| `--no-correct` (correction off, cutoff auto) | 35,068 | - | - | - |
| cutoff 1, correction **on** | 35,068 | 218 | 98.602 | 0.74 |
| cutoff 1, correction **off** | **51,105** | 165 | 98.832 | - |
| SPAdes | 51,105 | 139 | 98.286 | 0.60 |

Only the fourth row moves, and it needs both gates open at once. That row reaches SPAdes'
figure exactly.

**A defect made this hard to see, and it is fixed.** `-c N` used to set the abundance
cutoff *and* the trusted set the read corrector anchors on, because both read
`opt_.forcedCutoff`. At `-c 1` every k-mer became trusted, so no k-mer run ever broke and
the corrector silently did nothing -- every `-c 1` run on disk reports "0 bases corrected
in 0 reads". So `-c 1` was never a cutoff experiment; it was a cutoff-and-correction
experiment. The trusted set now has its own `--trust-cutoff` and never inherits `-c`.

The earlier version of this page concluded that the missing sequence was count-1 k-mers.
That was measured with the confounded flag and is wrong: with the corrector actually
running, cutoff 1 changes NGA50 by nothing at all.

## The complete 2x2, and where the gap actually is

With the corrector's trusted set decoupled from `-c`, and the mask threshold exposed as
`--mask-min-run`, the two gates can finally be varied independently. On GCF016591995v1:

| configuration | NGA50 | contigs | genome fraction | misassemblies | mismatches/100kb |
|---|---|---|---|---|---|
| default (mask on, cutoff auto) | 35,068 | 214 | 98.561 | 0 | **0.46** |
| mask off, cutoff auto | 35,068 | 212 | 98.584 | 0 | **0.46** |
| mask off + cutoff 1 | **48,447** | 168 | 98.806 | 2 | **2.48** |
| SPAdes | **51,105** | 139 | 98.286 | 1 | **0.60** |

Neither gate alone changes anything. Opening both buys +38% contiguity and costs 5.4x the
mismatch rate.

**SPAdes reaches the same contiguity at 0.60 mismatches, and that is the whole gap.** It
is not more permissive than TesserACT-with-both-gates-open; it is cleaner upstream.
BayesHammer will not promote a read unless every base of it is covered by a solid k-mer
(`hammer/expander.cpp`, `covered_by_solid`), so the reads reaching its graph are
trustworthy and a permissive cutoff is safe on them. TesserACT masks the unvouchable tail
and leaves the rest of the read in play, so opening the cutoff admits the errors still in
those reads along with the real sequence.

The rule that follows is read-level acceptance rather than tail masking: judge the whole
read by whether it is fully supported, instead of trimming the part that is not. That is
a change to correct.cpp, it is not made here, and the flags needed to measure it
(`--mask-min-run`, `--trust-cutoff`, `-c`) are now independent so it can be.

## It is not the thresholds. It is the corrector.

Masking harder does nothing: `--mask-min-run 1` masks the same 0.57% of input as the
default 8, because every unvouchable stretch is already longer than 8 bases. With the
cutoff open as well it still gives 35,068.

So the whole +38% in "mask off + cutoff 1" comes from the masked tails themselves. Those
tails are not noise -- they carry the real low-coverage sequence SPAdes recovers. They
also carry the errors, which is why keeping them raw costs 5.4x the mismatch rate.

That leaves one difference, and it is not a threshold:

* TesserACT's corrector will not change a base unless the k-mers following the correction
  stay solid for four more steps (`kMinCorroboration`, correct.cpp), and MASKS whatever it
  cannot vouch for. The reasoning is recorded in the file and is sound in itself: a
  corrector free to rewrite a degraded tail turns noise into plausible fiction, and that
  fiction enters the graph as evidence.
* BayesHammer CORRECTS those tails instead, by clustering each k-mer against solid centres
  under a Bayesian sub-cluster model (`hammer/kmer_cluster.cpp`), and only then decides
  what is solid.

So SPAdes keeps the sequence and keeps it clean: 51,105 at 0.60 mismatches, where
TesserACT must choose between 35,068 at 0.46 and 48,447 at 2.48. The missing capability is
a corrector strong enough to rescue a low-coverage tail rather than having to discard it.

That is a real piece of work, not a parameter, and it is the honest end of this line of
enquiry. No threshold in TesserACT closes this gap, and every threshold that was tried is
recorded above so the same ground is not covered twice.

## What implementing it would actually cost

The cheap form of BayesHammer's idea is: admit a count-1 k-mer only when it has no
Hamming-1 neighbour at high count -- an error is a near-miss of something abundant, a
genuine low-coverage k-mer is not. That is the discrimination a flat threshold cannot make
and it is why SPAdes can keep count-1 k-mers safely.

It does not fit TesserACT's counting as written, and the obstacle is memory rather than
logic. `extractSolid` (counter.cpp) compacts the sharded tables into the solid table and
frees each shard as it goes; counter.h states plainly that this "is what keeps peak memory
near the size of the genome rather than the size of the error cloud". A Hamming-1 lookup
needs the whole k-mer space available while deciding, because sharding is by hash and a
k-mer's 3k neighbours land in unrelated shards. So the rescue pass needs both structures
resident at once.

On these libraries the error cloud is 50-88% of distinct k-mers (from SPAdes' own
BayesHammer accounting on the same reads), so peak memory roughly doubles: 2-3 GB today,
5-7 GB with the rescue. That is a deliberate design decision being reversed, not an
oversight, and it should be taken as its own piece of work with its own measurement --
not slipped in alongside a benchmark run.

## What would fix it, and why it is not done here

Opening both gates recovers the contiguity and costs base accuracy: on the four isolates
where it was measured with correction off, mismatches per 100 kbp rose 2-12x
(0.72 -> 8.62, 0.69 -> 5.37, 0.00 -> 3.64, 2.22 -> 4.16) and misassemblies rose. That is
the trade, and it is not obviously worth taking -- TesserACT's base accuracy is one of the
things it currently wins on.

What SPAdes does instead is worth stating precisely, because it is not "keep everything".
BayesHammer filters at the READ level before any k-mer table exists, and it will not
promote a read unless every base of it is covered by a solid k-mer
(`hammer/expander.cpp`, `covered_by_solid`). On our own libraries that removes 50-88% of
distinct k-mers. What survives at count 1 in SPAdes' graph is therefore disproportionately
real low-coverage sequence, and SPAdes keeps it, applying its remaining error removal to
graph structures with a fitted per-dataset bound rather than to raw counts.

TesserACT's corrector is the analogous stage and it currently MASKS what it cannot vouch
for rather than dropping the read. The open question -- and the smallest experiment that
would settle it -- is whether masking less aggressively while keeping the cutoff at 2
recovers the contiguity without the mismatch cost. That is one flag away now that the
trusted set is decoupled, and it is what the cohort should be used for.

## The honest summary

TesserACT wins contig NGA50 on 70% of these isolates and loses on 30%. The losses
concentrate on libraries whose fragments barely exceed their reads. The mechanism is
low-coverage sequence removed twice over -- once by the corrector's masking and again by
the abundance cutoff -- and recovering it costs base accuracy, which is an axis TesserACT
currently wins. Repeats, the k ladder, pairing rules and contamination were each tested
directly and none of them is the cause.
