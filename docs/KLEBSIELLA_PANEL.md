# The 212-isolate *Klebsiella pneumoniae* panel

Every isolate here has **paired Illumina reads and a complete closed reference
genome for the same DNA** — chromosome plus every plasmid, from NCBI CP
accessions. That is what makes it possible to measure completeness and accuracy
rather than infer them from contiguity alone.

Both assemblers get the same `fastplus`-trimmed reads, and **every number is
contig level**: QUAST runs with `-s`, so tessera's scaffolds are split at their
N runs before scoring and neither tool is credited for a gap it did not close.
Comparing tessera's scaffolds against SPAdes' contigs flatters tessera and is
not done here.

Run with `kle_bench/run_batch.sh batchN.txt TAG`; score with
`kle_bench/headtohead.py TAG`.

## Library characteristics, and why they matter

These are mostly 2×250 MiSeq runs whose **fragments are shorter than the read
pair** — around 230 bp of insert against ~190 bp reads after trimming. Two
consequences run through everything below:

* **Paired-end resolution has almost nothing to work with.** A pair only links
  two unitigs if a junction falls between the two anchor positions, which is a
  window of about `insert − readLength` ≈ 36 bp. Measured: 2,595 linking pairs
  out of 702,392 — 0.37%, and in line with that arithmetic. This is the
  library's information content, not a defect in the anchoring.
* **k is doing the work instead**, which is why the k ceiling and the abundance
  cutoff turned out to matter far more here than the pairing logic.

Roughly 86% of pairs will merge into a single ~238 bp read. Merging was tried
and is *not* used by default: it halves k-mer depth by removing the overlap's
double-counting, which starves the top rungs of the ladder. On one isolate it
cost 73% of contig NGA50. tessera does accept the merged output
(`-1 unmerged_1 -2 unmerged_2 -s merged`) for anyone who wants it.

## What the panel found

Three defects, all fixed, in rough order of how much they cost.

### The abundance cutoff was set at the histogram valley

Much the largest. The valley is where the error shoulder meets the coverage
mode, so it looks like the principled place to threshold — but on a real genome
the counts in between are not empty. They hold the AT-rich prophages and
genomic islands, the low-copy plasmids, the k-mers flanking repeats.

The diagnosis came from diffing against SPAdes: of 83,765 bp tessera failed to
recover on one isolate, SPAdes recovered 55,663 bp, and **every large stretch
was GC 29–42% against a 57% genome**. Composition bias like that points at a
coverage threshold, not at an algorithm.

| cutoff | NGA50 | genome fraction | misassemblies | contigs |
|---|---|---|---|---|
| 5 (valley) | 39,855 | 97.43% | 1 | 328 |
| 3 | 86,779 | 98.15% | 0 | 195 |
| **2** | **132,059** | **98.41%** | 1 | **157** |
| *SPAdes* | *81,707* | *98.58%* | *3* | *188* |

The error k-mers a low cutoff admits do not survive anyway — they enter the
graph as tips, bubbles and erroneous connections, and simplification removes
them *by topology*, which is evidence a count threshold does not have. Be
permissive at the k-mer level and strict at the topology level.

### The k ceiling was set for 150 bp reads

The 192-bit k-mer capped k at 96. On 2×250 libraries the unitig N50 was still
climbing steeply there. Widening to 256 bits (k ≤ 128) and topping the ladder at
127 took unitig N50 from 142,604 to 207,162 on one isolate and mismatches from
2.78 to 0.52 per 100 kbp.

This directly contradicts what the older seven-isolate benchmark concluded, and
both are correct for their data: on those 2×150 and 2×300 libraries the ladder
had genuinely plateaued below k=96, so widening would have bought nothing. The
conclusion was right about that panel and wrong as a general claim.

### The ladder keyed off maximum read length

Many public runs arrive already trimmed by the submitter, so their reads are
variable-length. One 301 bp survivor in a library averaging 190 bp was
selecting a ladder most of the reads were too short to contribute a k-mer to.
It now uses the mean.

## Batch 1: the session's effect, contig level

| isolate | before | after | change |
|---|---|---|---|
| ERR11578909 | 35,768 | 132,059 | +269% |
| ERR11578757 | 72,904 | 265,033 | +264% |
| ERR11579004 | 78,313 | 240,856 | +208% |
| ERR11578413 | 101,486 | 265,335 | +161% |
| ERR11578712 | 73,879 | 173,713 | +135% |
| ERR11578485 | 49,939 | 102,944 | +106% |
| ERR11578347 | 45,342 | 88,184 | +94% |
| ERR11578610 | 99,487 | 168,365 | +69% |
| ERR2631558 | 162,079 | 265,037 | +64% |
| ERR11578154 | 200,724 | 304,809 | +52% |
| ERR11578086 | 224,616 | 282,003 | +26% |
| ERR12791282 | 137,038 | 143,286 | +5% |
| ERR10447223 | 373,371 | 387,438 | +4% |
| ERR11579075 | 367,615 | 377,674 | +3% |
| ERR11578240 | 304,377 | 304,165 | −0% |
| ERR11578571 | 348,740 | 283,971 | **−19%** |
| ERR11578837 | 219,607 | 168,470 | **−23%** |

Median NGA50 137,038 → 265,033. Genome fraction up on all 17, median 98.64% →
98.78%. Misassemblies 9 → 6.

The two regressions are the honest cost of the permissive cutoff: on those
isolates the extra low-coverage k-mers add graph branching faster than
simplification removes it, halving unitig N50 (160,804 → 80,136 on ERR11578837).
Raising `--simplify-rounds` does not recover them — tested, no change.

## Against vanilla SPAdes

Fifteen isolates, same trimmed reads, contig level throughout.

| | tessera | SPAdes |
|---|---|---|
| NGA50 wins | 4/15 | 11/15 |
| median NGA50 ratio | **0.88×** | — |
| genome fraction | 98.71% | 98.85% |
| **misassemblies (total)** | **6** | **11** |
| wall clock | ~100 s | ~180–260 s |

**SPAdes still leads contiguity by about 12% at the median.** tessera makes
roughly half the structural errors and runs 2–3× faster on a third to a half the
memory. That is the position; it should not be quoted any other way.

## Where the remaining gap is

Repeat resolution, and specifically what to do when the read pairs say nothing.
With the graph now well connected, the resolver's dominant failure changed from
"no candidate continuation exists" (2,620 chain ends before, 636 after) to
"candidates exist and no paired evidence separates them" (616, against 239
successful joins). On these overlapping-mate libraries that is the normal case,
so anything that decides those without inventing joins is worth more than any
further work on the pairing itself.
