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
  two unitigs if a junction falls between the two anchor positions, a window of
  about `insert − readLength` ≈ 36 bp. Measured: 2,595 linking pairs out of
  702,392 — 0.37%, in line with that arithmetic. This is the library's
  information content, not a defect in the anchoring.
* **k and the abundance cutoff do the work instead**, which is why they turned
  out to matter far more here than the pairing logic.

Roughly 86% of pairs will merge into a single ~238 bp read. Merging was tried
and is **not** used: it halves k-mer depth by removing the overlap's
double-counting, which starves the top of the ladder — on one isolate it cost
73% of contig NGA50. tessera does accept merged output
(`-1 unmerged_1 -2 unmerged_2 -s merged`) for anyone who wants it.

## What the panel found

In descending order of how much it cost. None of it was where the search
started.

### 1. The abundance cutoff was set at the histogram valley

The valley is where the error shoulder meets the coverage mode, so it looks
like the principled place to threshold. But on a real genome the counts in
between are not empty: they hold the AT-rich prophages and genomic islands, the
low-copy plasmids, the k-mers flanking repeats.

The diagnosis came from diffing against SPAdes. Of 83,765 bp tessera failed to
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
them *by topology*, which is evidence a count threshold does not have.

### 2. The trimming preset was throwing the assembly away

`fastplus --preset assembly` ran `cut_right` at Q20 and enabled overlap
correction. Both cost contiguity, and together they explained the entire gap to
plain fastp — to the base.

| trimming | NGA50 |
|---|---|
| `cut_right` Q20 + correction *(what it did)* | 132,059 |
| `cut_tail` Q10 + correction | 149,192 |
| `cut_tail` Q10, no correction | **160,991** |
| fastp's own defaults | **160,991** |

`cut_right` truncates a read at the first sliding window below the threshold,
discarding 13–24% of all bases. Overlap correction rewrites disagreeing bases
toward the higher-quality mate — wrong whenever the mates disagree because they
came off different repeat copies. On one isolate the fix was +95% NGA50.

Same lesson as the cutoff, in a different component: an assembler removes errors
by graph topology, pooling evidence across every read covering a locus; a
trimmer decides per read, in advance, with far less information.

### 3. The k ceiling was set for 150 bp reads

The 192-bit k-mer capped k at 96. On 2×250 libraries the unitig N50 was still
climbing steeply there. Widening to 256 bits (k ≤ 128) and topping the ladder at
127 took unitig N50 from 142,604 to 207,162 on one isolate and mismatches from
2.78 to 0.52 per 100 kbp.

This contradicts what the older seven-isolate benchmark concluded, and both are
correct for their data: on those 2×150 and 2×300 libraries the ladder had
genuinely plateaued below k=96. The conclusion was right about that panel and
wrong as a general claim.

### 4. Smaller findings

* **The ladder keyed off maximum read length.** Many public runs arrive already
  trimmed by the submitter, so one 301 bp survivor in a library averaging 190 bp
  selected a ladder most reads were too short to contribute to. Now uses the mean.
* **A chain end with candidates but no paired support was rejected outright**,
  even when every candidate ended on the same unitig and only the route through
  the repeat differed. The destination is not in doubt there. Taking those joins
  is +26% and +23% on two isolates, with misassemblies unchanged.
* **Scaffold gaps had no closing stage.** Local reassembly from the reads
  crossing each gap closes them; the corrector's masking had been hiding the
  very reads needed, so gap closing reads through the mask.
* **`minLinkSupport` was a flat count.** Pairs crossing a junction scale with
  coverage, so 2 is a real bar at 40x and almost none at 100x. It became visible
  the moment better trimming raised the depth: contiguity rose and so did
  long-range chimeric joins -- 4.4 Mb and 1.6 Mb inconsistencies inside single
  contigs. It now scales with the median unitig coverage, and the factor was
  measured across all twenty isolates rather than guessed:

  | factor | median NGA50 | misassemblies |
  |---|---|---|
  | flat 2 | 303,624 | 12 |
  | 0.07 | 291,265 | 7 |
  | **0.10** | **277,799** | **3** |

  The trade is close to linear, so this is a judgement rather than an optimum.
  0.10 is the default because a chimeric contig is the worse error for what
  these assemblies get used for -- it misplaces resistance genes relative to
  their genomic context and corrupts cgMLST, while a shorter contig only costs
  contiguity. The knob follows the run modes: careful 0.14, standard 0.10,
  aggressive 0.05.

## Batch 1: the session's effect, contig level

| | before | after |
|---|---|---|
| median NGA50 | 181,402 | **277,799** (+53%) |
| median genome fraction | 98.67% | **98.79%** |
| median mismatches per 100 kbp | — | **0.11** |
| misassemblies across 20 isolates | 10 | **3** |

Seven isolates more than doubled: ERR11578909 +317%, ERR11578757 +275%,
ERR11578712 +255%, ERR11579004 +225%, ERR2631558 +133%, ERR11578347 +119%,
ERR11578485 +106%. Two regressed: ERR11578571 −25% and ERR11578240 −4%, both
under the depth-scaled support bar that took the panel from twelve
misassemblies to three.

## Does `--mode aggressive` buy the contiguity back?

The two isolates furthest behind SPAdes are also the two that gave up most to
the depth-scaled support bar, so the obvious question is whether the mode built
for that trade recovers it. Sometimes, and sometimes it is simply worse:

| isolate | standard | aggressive |
|---|---|---|
| ERR11578413 | 184,767 · 0 misasm · 0.08 mm | **262,699** · 0 misasm · 0.21 mm |
| ERR11578610 | **170,619** · **0** misasm · **0.22** mm | 162,645 · **3** misasm · **4.72** mm |

On the first it recovers the whole gap for nothing -- higher genome fraction
too. On the second it loses on every axis at once: shorter contigs, three
misassemblies, and twenty times the mismatch rate. Which of the two an isolate
will behave like is not knowable in advance, which is the argument for the
cautious default and for keeping the aggressive mode available rather than
tuning towards it.

## Across the whole panel

Batch 1 is scored against SPAdes; the rest measure whether the improvements hold
outside the twenty isolates they were developed on. They do.

| | batch 1 (20) | batches 3+ (127) |
|---|---|---|
| median contig NGA50 | 277,799 | 231,529 |
| median genome fraction | 98.79% | 98.76% |
| misassemblies per genome | 0.15 | 0.49 |
| median mismatches / 100 kbp | 0.11 | 0.17 |

Restricting to the 109 isolates whose reads and reference agree (under 1 kb
unaligned) barely moves it: NGA50 234,832, genome fraction 98.78%, 0.15
mismatches. So the defaults were not overfitted to batch 1 -- the numbers are
the same shape on 127 isolates that had no part in choosing them.

## Reference-free cross-check

QUAST asks how close an assembly is to the truth. Mapping the reads back asks
how much of the evidence it accounts for, which is the question that survives
when there is no closed reference. Both assemblers explain essentially all of
it, so neither is discarding data:

| isolate | tessera | SPAdes |
|---|---|---|
| ERR11578909 | 99.89% | 99.95% |
| ERR11578837 | 99.87% | 99.94% |
| ERR11578347 | 99.87% | 99.88% |
| ERR11578086 | 99.94% | 99.95% |

## Runtime

Graph construction was profiled with `TESSERA_GRAPH_PHASES=1`. The unitig walk
was 78–84% of it, against 13–16% for the start classification. The walk
parallelises exactly — the walks are disjoint by construction — taking graph
build from 8.5 s to 3.4 s per k-rung with byte-identical output.

## Five things that looked right and were not

Each was measured, rejected, and the measurement recorded in the code so the
same plausible argument does not get made twice.

* **Merging overlapping mates.** +32% on the hardest isolate, −73% on another.
  It halves k-mer depth.
* **More simplification rounds.** A 4× mismatch improvement on the first isolate
  tried; reproduced on neither of the next two, identical NGA50 *and* identical
  mismatches on both.
* **Widening the repeat-enumeration budget.** Changed no output at all — those
  chain ends were graph dead ends, not budget exhaustion. Kept anyway, as the
  correct separation of "what continuations exist" from "which do the reads
  prefer", but with no measured gain claimed.
* **Polishing being blind to repeats.** It covers 99.9998% of assembly positions
  at 64x. It sees essentially everything.
* **Polishing's threshold being too strict.** Lowering it from 0.90 to 0.70 made
  15 corrections and took mismatches from 0.49 to 0.80 per 100 kbp; every
  correction was wrong. The positions where the pileup disagrees with the contig
  are collapsed repeats, where the read majority belongs to whichever copy is
  commonest rather than to the locus being polished. Near-unanimity is what stops
  the stage rewriting one repeat copy into another, so polishing is a guard, not
  a corrector, and should stay silent.

## A caution about the panel itself

Eight of batch 2's twenty isolates have ≥1 kb of assembly that does not align to
their own closed reference, up to 53 kb — contamination or imperfect
strain/reference pairing in the public dataset, not an assembler defect. Those
eight carry 11 of that batch's 19 misassemblies and a mean 2.19 mismatches per
100 kbp against 0.27 for the twelve cleanly paired ones. Scoring should report
the clean subset separately.

## The remaining contiguity gap, and what it is not

SPAdes leads contig NGA50 by about 14% at the median. Eleven explanations were
measured against the closed reference on the isolate furthest behind
(ERR5056466: ours 472,260, SPAdes 739,807) and none of them is the cause:

| tested | result |
|---|---|
| SPAdes' read correction | tessera on their BayesHammer reads: 461,007 |
| k ceiling | k<=95: 268,882 against k=127's 472,260 |
| bubble popping | disabled entirely: NGA50 identical |
| chimera removal | 0.12 / 0.30 / 0.50: all 472,260 |
| repeat-gated support bar | four extra misassemblies, no better trade |
| support-bar level | flat 2 / 0.07 / 0.10 trade linearly |
| matching dominance | 3.0 and 1.5 both fire exactly ten times |
| risk appetite | `--mode aggressive`: 460,761, below the default |
| multi-k carry weight | 4 / 16 / 48: 472,260 / 451,613 / 414,314 |
| merged fragments as spanning evidence | 361,437, worse than pairs alone |
| unfair comparison | both have zero Ns, duplication ratio 1.000 |
| circularity / layout constraint | built and measured; see below |

The junctions SPAdes crosses and we do not are short branch nodes, 219-398 bp,
in-degree two and out-degree two, whose sequence occurs exactly twice in the
genome. Every 40-mer of it is in our graph. No pair can span them -- read + node
+ read is about 708 bp against a 355 +/- 119 fragment distribution -- and the
matching rule confirms there is nothing to weigh: the intended assignment scores
zero.

### The circularity constraint: built, measured, and not the answer either

A bacterial replicon is a circle, so picking the wrong pairing at a two-copy
repeat splits it into two disjoint cycles while the right one leaves it whole.
That is a property of the whole layout, which is exactly why local paired
scoring cannot reach it -- so it looked like the missing ingredient.

It was implemented twice and neither version helps.

The first was a standalone pass over unplaced 2-in/2-out repeat nodes, choosing
the pairing that leaves fewer connected components. It never fired: of 58 nodes
with the right shape, **none** had all four neighbours as free chain ends. The
junctions are not single repeat nodes between two chains, they are small tangles
of several short nodes, which is what path enumeration already walks through.

The second put the same test inside candidate evaluation, as a last resort when
no local evidence separates the options: prefer a candidate whose destination is
not already in this chain's component, since taking one that is would close the
replicon early. It fires eight times on the worst isolate and lifts joins from
538 to 546 -- and changes NGA50 on **none** of four isolates measured, nor any
misassembly count. The joins it makes are between small contigs that do not
reach N50.

A third version did the joint optimisation the second one implied: collect every
ambiguous junction, enumerate all 2^m combinations of pairings, count the
components each leaves, and act only when the minimum is achieved by exactly one
combination. It collected **zero** junctions.

All three fail the same way, and that is the real finding. Each is built on the
resolver's abstraction -- a chain end, a repeat, another chain end -- and at the
junctions that matter the terminals are not clean chain ends at all. They are
tangles of several short unitigs, some unanchorable, some mid-chain. The first
version found 58 repeat nodes with the right shape and none with four free chain
ends; the third found no ambiguous junction with two clean chain-end
destinations.

So the gap is not reachable from the chain abstraction, by any decision rule
layered on top of it. Closing it means path extension over the unitig graph
itself, carrying global constraints as it goes -- which is exSPAnder's design,
and is why that is a rewrite of the resolver rather than an addition to it.
Three measured attempts say so, rather than one argument.
