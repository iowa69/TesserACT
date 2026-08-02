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

## Head to head with vanilla SPAdes, two independent batches

Forty isolates, same fastplus-trimmed reads for both assemblers, contig level
throughout. Batch 3 played no part in choosing any parameter.

| | batch 1 (20) | batch 3 (20) |
|---|---|---|
| median NGA50 ratio | 0.86x | 0.87x |
| NGA50 wins | 4/20 | 1/20 |
| genome fraction | 98.78% vs 98.94% | 98.74% vs 98.89% |
| **misassemblies** | **3 vs 6** | **9 vs 20** |
| wall clock | ~70 s vs ~340 s | same |

**SPAdes leads contiguity by about 15% at the median, on both batches.** That is
the position and it is stable: the ratio measured on isolates used for tuning
and on isolates that were not is the same to within a percentage point.

Against that, tessera makes **less than half the structural errors** -- 12
against 26 across the forty -- and wins per-base accuracy on most isolates,
while running four to five times faster.

The mechanism behind both halves of that trade is the same one traced in the
elimination study above: SPAdes joins at ambiguous repeats where the paired
evidence is zero, and is usually but not always right. Its misassembly rate is
what "usually" costs.

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

That looked like proof the chain abstraction was the limiter, and that closing
the gap needed path extension over the unitig graph itself -- a rewrite rather
than an addition.

**That was wrong, and building it is what showed so.** An alternative resolver
that walks the graph edge by edge, with a per-unitig traversal budget from
coverage so a repeat can be used as often as its depth says it occurs, is worse
on five isolates:

| isolate | chain NGA50 | extend NGA50 | misassemblies |
|---|---|---|---|
| ERR5056466 | 451,301 | 362,493 | 0 -> 0 |
| ERR11578413 | 184,767 | 135,007 | 0 -> 3 |
| ERR11578909 | 149,192 | 139,096 | 0 -> 5 |
| ERR11578347 | 99,464 | 92,200 | 0 -> 9 |
| ERR11578086 | 282,003 | 188,136 | 0 -> 14 |

Lower contiguity on every one, and thirty-one misassemblies where the chain
resolver makes none. It does recover more of the genome -- 0.4 to 0.8 points of
genome fraction, on one isolate beating SPAdes -- but pays for it with exactly
those errors, which is the same guessing trade seen everywhere else in this
study rather than a free win.

Adding look-ahead so the walk pools evidence the way chain enumeration does
changed almost nothing (362,493 against 361,930). Spending the traversal budget
per path rather than globally produced an assembly 3.4x the size of the genome.

So the chain abstraction is not what limits contiguity here. Where the limit
actually lies is, after all of this, still unidentified -- which is a more
honest place to stop than a confident architectural claim that measurement
contradicts.

## Where the limit actually lies

The claim above -- that the limit was unidentified -- lasted until the
resolver was asked what it was rejecting rather than what it was choosing.

On ERR11578413, of 498 chain ends:

    no-partner=461  low-support=35  not-mutual=0  taken=2

and of 1,635 continuation decisions:

    ok=665  no-candidate=405  low-support=387  tie=0  by-coverage=125

Ninety-three percent of chain ends have **no paired partner at all**, and 405
continuations have **no graph candidate at all**. Neither number is a decision
the resolver got wrong; both are the absence of anything to decide. Every knob
tried over this study -- the link bar, the tie ratio, matching dominance, risk
appetite -- adjudicates the 387 low-support cases, which is a quarter of the
problem at most.

Aligning what the assembly fails to cover against the reads shows what sits at
those dead ends. 44,662 bp of ERR11578413 is covered by reads and missing from
the assembly, and the largest pieces are eight stretches of 3.5-4.0 kb at ~52%
GC and 1.3-2.4x typical depth: the rRNA operons. They are not absent from the
assembly. They are shattered:

| | tessera | SPAdes |
|---|---|---|
| operon at chr 4,151,374-4,155,420 | five contigs, 365-428 bp, cov 212-280x | one contig, 2,801 bp, cov 211x |

Two mechanisms were tested and both are innocent. Disabling read correction
entirely leaves the gap boundaries **byte-identical**, so the corrector is not
collapsing the copies. Bubble popping never fires there at all -- zero bubbles
in the late simplification rounds -- because an operon tangle is not a
two-branch bubble, and neither raising the length cap fourfold and tenfold nor
adding a rule that collapses bubbles whose both sides are multi-copy changes a
single contig. Both knobs were reverted as inert.

What is left is arithmetic. The fragments in this library are 230-355 bp; the
repeats are ~4,000 bp. `scoreCandidate` accepts a continuation only when the
fragment length it implies is plausible, so a pair can never vouch for a walk
across an operon -- not because the bar is too high, but because no such
fragment exists. **The remaining contiguity is not reachable from paired-end
evidence on this library by any decision rule.**

SPAdes crosses them anyway, on coverage multiplicity rather than evidence, and
the price is visible in the same reports: on ERR11578413, contig level,

| | tessera | SPAdes |
|---|---|---|
| NGA50 | 234,964 | 320,313 |
| genome fraction | 98.95 | 99.14 |
| contigs | 65 | 53 |
| misassemblies | 0 | 0 |
| mismatches per 100 kbp | **0.08** | 0.72 |

Nine times the per-base error rate for 36% more NGA50. Across the 40-isolate
head-to-head the same trade shows up as 26 misassemblies against our 12. And
when that heuristic was implemented here -- the extension resolver above -- it
reproduced the trade exactly: more genome fraction, thirty-one new
misassemblies.

So the gap is real, it is measured, and it is a choice rather than a defect.
Closing it means buying contiguity with error. That is worth stating plainly
rather than continuing to search for a decision rule that cannot exist.

## The scaffold splice, and why it was reverted

One defect was real. The scaffold stage refused any partner reachable
through the graph, on the grounds that chain extension had handled it; chain
extension had in fact rejected those junctions for want of support on a single
from->terminal pair, while this stage pools every pair landing within a
fragment of the chain end. The evidence was being dropped twice and the
junction joined by nobody. Letting it act there, splicing the interior unitigs
rather than a run of Ns and only when both ends independently see the same
unique path, gained contiguity on five isolates with no misassembly.

Widening the test to every batch-1 isolate with a baseline reversed the
verdict:

| isolate | dNGA50 | misassemblies | mismatches/100 kbp |
|---|---|---|---|
| ERR10447223 | +9,754 | 0 -> 1 | 0.11 -> **0.46** |
| ERR11578413 | +2,910 | 0 -> 0 | 0.08 -> 0.08 |
| ERR11578571 | 0 | 1 -> 1 | 0.78 -> 0.78 |
| ERR11578610 | 0 | 0 -> **1** | 0.22 -> 0.22 |
| ERR11578909 | 0 | 0 -> 0 | 0.28 -> 0.28 |
| ERR5056466 | +509 | 0 -> 0 | 0.12 -> 0.12 |

+13,173 bp of NGA50 across six isolates, bought with two misassemblies and a
fourfold per-base accuracy regression on one of them -- and on ERR11578610,
bought with no contiguity gain at all.

The obvious repair does not work. A pair can only vouch for a junction it could
physically span, so requiring the spliced interior to fit inside the largest
plausible fragment should have removed the coincidental support. It changes
nothing: both regressing isolates come out byte-identical, because the
offending interiors were already within reach. The joins are spannable and
still wrong.

So the change was reverted. It is the same trade this document argues against
everywhere else, and the fact that it was our own change does not make the
trade a better one. Five isolates were not enough to see it; the honest sample
here is every isolate with a baseline.
