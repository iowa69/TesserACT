# Where TesserACT loses, by how much, and why

All five arms on the complete 180-isolate non-clonal *S. aureus* panel, contig level
(QUAST `-s`, `_broken` column), each isolate against its own closed reference.

The arms are: **vanilla** (no model), **`--organism`** (the shipped default with a genus
model), **`careful`** and **`aggressive`** (presets on top of the model), and **SPAdes 4.3.0**
with its own error correction and automatic k.

## The whole table

| metric (median) | vanilla | `--organism` | careful | aggressive | SPAdes |
|---|---|---|---|---|---|
| NGA50 | 216,080 | **247,188** | 246,817 | 240,484 | 188,658 |
| NG50 | 224,000 | **253,984** | 247,921 | 245,751 | 196,578 |
| Genome fraction (%) | 98.80 | 98.83 | 98.83 | **99.01** | 98.37 |
| Largest alignment | 482,713 | **538,873** | 509,955 | 510,502 | 427,134 |
| # contigs | 42 | **38** | 40 | 39.5 | 40 |
| # misassemblies | 0 | 0 | 0 | 0.5 | 0 |
| Mismatches /100 kb | 0.59 | 0.66 | **0.56** | 1.33 | 1.22 |
| Indels /100 kb | 0.22 | 0.25 | **0.21** | 0.32 | 0.29 |
| Duplication ratio | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |
| **total misassemblies** | 201 | 223 | 216 | 272 | **150** |
| **total local misassemblies** | 119 | 144 | 140 | 151 | 129 |
| **isolates with any misassembly** | 61 | 70 | 70 | 90 | **47** |

### Each arm against SPAdes (win = the TesserACT arm is better on that isolate)

| metric | vanilla | `--organism` | careful | aggressive |
|---|---|---|---|---|
| NGA50 | 113/67 (0.008) | **133/47** (2e-10) | 122/58 (4e-06) | 123/57 (8e-06) |
| Genome fraction | 170/10 (4e-27) | 170/10 (1e-27) | 170/10 (1e-27) | **172/8** (5e-29) |
| Largest alignment | 110/67 (0.035) | **126/52** (3e-08) | 119/60 (1e-04) | 122/57 (6e-06) |
| # contigs | 63/111 (0.014) | 88/86 (ns) | 80/91 (ns) | 79/90 (ns) |
| # misassemblies | 30/45 (ns) | 29/53 (0.011) | 25/51 (0.017) | 24/75 (6e-06) |
| # local misassemblies | 49/26 (0.034) | 45/37 (ns) | 43/33 (ns) | 43/41 (ns) |
| Mismatches /100 kb | 118/57 (6e-05) | 111/65 (0.002) | 119/58 (5e-06) | 69/108 (0.006) |
| Duplication ratio | 4/151 (4e-25) | 5/151 (8e-25) | 3/152 (3e-26) | 3/171 (4e-29) |

### Each option against vanilla (what each one actually changes)

| metric | `--organism` | careful | aggressive | SPAdes |
|---|---|---|---|---|
| NGA50 | **101/0** (3e-18) | 102/42 (1e-05) | 119/37 (1e-06) | 67/113 (0.008) |
| Genome fraction | 97/28 (6e-07) | 97/68 (0.003) | **170/9** (4e-30) | 10/170 (4e-27) |
| # contigs | **157/1** (1e-27) | 107/52 (4e-08) | 122/46 (2e-11) | 111/63 (0.014) |
| # misassemblies | 2/18 (0.001) | 19/30 (ns) | 9/59 (4e-08) | 45/30 (ns) |
| Mismatches /100 kb | 15/42 (4e-05) | 68/55 (ns) | 20/153 (4e-22) | 57/118 (6e-05) |
| Duplication ratio | 52/13 (6e-06) | 42/34 (ns) | 20/130 (2e-19) | 151/4 (4e-25) |

Read the presets honestly. **`aggressive` buys genome fraction and pays for it everywhere
else**: +170/9 on genome fraction against vanilla, but 9/59 on misassemblies and 20/153 on
mismatches. **`careful` is the accuracy arm** — best mismatches and indels of any arm, at
some contiguity. The default `--organism` is the contiguity arm and never loses NGA50 to
vanilla on a single isolate out of 180.

## The two losses, sized

| loss vs SPAdes | isolates losing | median loss | worst | total across cohort |
|---|---|---|---|---|
| Duplication ratio | 151 / 180 | **0.0020** | 0.452 | 0.964 |
| # misassemblies | 53 / 180 | 1 | 29 | 114 |

**The duplication loss is statistically overwhelming and practically small.** p=8e-25 on
151 losses, and the median loss is 0.0020 — about **5,600 extra bases on a 2.8 Mb genome**,
0.2%. Only **5 of 180** isolates exceed 1% extra sequence, and the worst of them
(GCF046742145v1 at 1.458) is the same isolate that contributes 29 of the 114 excess
misassemblies. A p-value that small on an effect that small is a statement about
consistency, not importance: we are almost always very slightly more redundant.

**The misassembly loss is smaller than the raw count suggests but is real.** See
[HOW_SPADES_CHOOSES_A_BRANCH.md](HOW_SPADES_CHOOSES_A_BRANCH.md): part is the price of
contiguity, and a residual survives contiguity-matching at p=0.022.

## Where the advantage actually goes: paired reach

**This section was wrong in its first version and is corrected here. The correction is
recorded rather than quietly edited, because the mistake is instructive.**

Define **paired reach** = fitted insert size − read length. It is how far past its own read
a pair can vouch for; when the fragment is no longer than the read, the mate says nothing the
read did not already say.

The first version took "read length" from the assembly log's `max length`. These libraries
are not fixed-length — median mean read length across the panel is **202 bp against a median
max of 251**, and 36 of 180 isolates have a max more than 20% above their mean. That inflated
the read length and invented a 23-isolate band of "overlapping mates" which, measured against
the mean, is **one isolate**. (`baseKLadder()` in `assembler.cpp` already uses the mean for
exactly this reason, and its comment describes the trap precisely. The assembler was right;
the analysis was not.)

Recomputed against mean read length:

| library | n | med insert | med read | our NGA50 | SPAdes | ratio | win/loss | p |
|---|---|---|---|---|---|---|---|---|
| mates overlap (reach <= 0) | 1 | 126 | 145 | — | — | 4.20x | 1/0 | — |
| **reach 0–60 bp** | **38** | 248 | 201 | — | — | **1.00x** | 25/13 | **0.131 (ns)** |
| reach 60–120 bp | 58 | 328 | 237 | — | — | 1.11x | 43/15 | 6.2e-04 |
| reach > 120 bp | 83 | 368 | 149 | — | — | 1.29x | 64/19 | 3.0e-07 |

Spearman across all 180 isolates, reach against our NGA50 ratio: **+0.193, p=0.0085** —
identical whether reach is computed from mean or max read length, so the *trend* was never in
doubt; only the band boundaries were.

**On 38 of 180 isolates our entire margin over SPAdes disappears.** Not a loss — a parity we
should not be at, given we are 1.29x ahead when the pairs have reach. Losing and winning
isolates have the same coverage (94.6x against 96.8x median) and the same read length.

### It is not adapter read-through

The first version attributed this to adapter read-through and proposed trimming it. Measured
on the raw FASTQ, stride-sampled across the whole file rather than its head:

| | GCF010364725v2 |
|---|---|
| pairs with a **perfect, zero-mismatch** overlap | **69.2%** |
| median overlap length | 127 bp of a 131 bp read |
| pairs whose fragment is shorter than R1 (true read-through) | **3.3%** |

The overlap is real and near-universal; read-through is rare. The fragments are barely longer
than the reads, so **a pair is two near-duplicate copies of the same ~120 bp of genome**.
Paired reach is zero by construction. That is a property of the library, not a defect in the
resolver, and no amount of adapter trimming addresses it.

`scepter` (fastplus) agrees independently and fails all three of the worst libraries outright:
insert peaks of **78, 35 and 35 bp** against 301- and 351-cycle reads, and on the worst one
reported quality is 8.5 Q units optimistic against its empirical error rate.

## Four fixes tested on the three worst isolates, four negatives

The three: GCF010364725v2 (56x, mean read 136), GCF046268025v1 (140x, 163), GCF045347525v1
(53x, 187). NGA50 for every arm:

| arm | GCF010364725v2 | GCF046268025v1 | GCF045347525v1 |
|---|---|---|---|
| vanilla | 5,653 | 48,982 | 77,166 |
| `--organism` | 5,781 | 53,876 | 77,166 |
| careful | 5,834 | 51,574 | 77,166 |
| aggressive | 5,781 | 51,574 | 77,166 |
| mates merged (scepter) | **1,726** | 56,613 | 78,519 |
| ladder capped at k=77 | 5,910 | 51,632 | **69,186** |
| **SPAdes** | **14,366** | **158,959** | **106,210** |

1. **Adapter trimming.** Implemented (overlap detection, `--no-overlap-trim`, activates on a
   stride sample). Correct, and inert: the condition is 3.3% of pairs, below any sane
   threshold. Not committed — a fix for a problem this rare does not earn a place in the
   loader.
2. **Mate merging.** Two isolates gain 2–5%; the third loses genome fraction **94.94% ->
   81.07%** and doubles its contig count. Not shippable as a default.
3. **Shorter k ladder.** Helps the worst isolate by 2%, costs the other two 4% and 10%. Our
   tall ladder is earning its keep; effective coverage falling from 28.8x to 20.9x across the
   top rungs is a price worth paying.
4. **Library depth** — does not separate winners from losers at all.

Note the third column: **all four TesserACT arms return exactly 77,166**, to the base pair.
Model, presets, resolution settings — none of them move it. Whatever caps us on that isolate
is upstream of repeat resolution, in the graph itself.

## The corrector hypothesis, tested and falsified

`WHY_SPADES_SOMETIMES_WINS.md` concluded, from the threshold side, that the gap is the read
corrector rather than any cutoff. That was the last standing explanation, and it is wrong.

The test: run SPAdes' BayesHammer alone (`--only-error-correction`) and assemble **its own
corrected reads** with TesserACT, once with our corrector off (`--no-correct`, a clean
head-to-head of the two correctors over identical downstream code) and once with it on.

| | ours | BH reads, ours OFF | BH reads, ours ON | SPAdes |
|---|---|---|---|---|
| GCF010364725v2 NGA50 | 5,781 | **6,868** | 6,868 | **14,366** |
| — genome fraction | 94.94 | 95.55 | 95.41 | 97.18 |
| — contigs | 704 | 677 | 674 | **329** |
| GCF045347525v1 NGA50 | 77,166 | **77,302** | 77,302 | **106,210** |
| — contigs | 71 | 68 | 69 | 66 |

BayesHammer's reads close **13% of the gap on one isolate and 0.5% on the other**. And the
two BH columns are the same to within noise: our corrector adds nothing on top of theirs, and
removing ours costs nothing. The reads are not the difference.

## Where the difference actually is

Given identical corrected reads, the contig length profiles diverge:

| | contigs | total | N50 | largest | >10 kb | >1 kb |
|---|---|---|---|---|---|---|
| TesserACT | 865 | 2,709,511 | 7,015 | 41,743 | 57 | 536 |
| SPAdes | 353 | 2,725,898 | **14,494** | **79,543** | **79** | **295** |

Same reads, same genome recovered, and we emit 241 more contigs over 1 kb while producing 22
fewer over 10 kb. Our own unitig graph saturates: N50 by rung is 3,197 -> 5,719 -> 5,929 ->
5,956 at k=21/33/45/55 and stops improving after k=45.

Both assemblers do multi-k with graph carry-over (`carryOver` in `Assembler::iterate`), so
that architectural difference is not it either. What is left is graph **simplification** --
tip clipping, bulge removal, erroneous-connection removal, relative-coverage removal -- where
SPAdes has a long-tuned pipeline and we have tips, bubbles and chimeras. On a 40x library with
a 155 +/- 93 bp insert, that is where the remaining 2x lives.

That is not a lever to flip. It is a body of work, and naming it honestly is more useful than
proposing a fix that has not been measured.

## What to fix, in order of expected value

1. **Plasmid copy number in the repeat test.** The one improvement here that is measured and
   works: whole-plasmid recovery 33% -> 44% on a 39-isolate subset, with the low-copy control
   band unmoved and one extra misassembly across 39 isolates. See
   [PLASMID_COPY_NUMBER.md](PLASMID_COPY_NUMBER.md).
2. **Graph simplification on low-coverage, short-insert libraries.** Localised above but not
   attempted. The measurable target is the 38 isolates at reach 0-60 bp sitting at 1.00x
   against SPAdes where the rest of the cohort is 1.29x ahead.
3. **The residual misassembly gap at matched contiguity** (p=0.022), mechanism still unknown.

### Keeping this in proportion

The three isolates above were chosen as the worst in the panel. Over all 180, TesserACT wins
NGA50 by 31.0% (133/47, p=2e-10) and eight of nine metrics. What this section describes is the
tail, not the centre.

## What was tested and did *not* explain a loss

Recorded so the same ground is not covered twice:

* **`pickByCoverage`**, the resolver's depth-continuity fallback — ablated on 23 isolates,
  every metric identical to the base pair. It fires up to 38 times per isolate and changes
  the contigs, but arrives at the same assembly.
* **Gap filling** (`--no-gapfill`) — removes 9 of 65 misassemblies on the worst 10 isolates
  and costs NGA50 6/0 (p=0.036) and contigs 10/0 (p=0.006). A real contributor at 14%, and
  a bad trade.
* **Library depth** — does not separate winning from losing isolates (94.6x vs 96.8x).
* **The read corrector** — BayesHammer's own corrected reads close 13% and 0.5% of the gap on
  the two isolates tested. Falsifies the standing conclusion in
  [WHY_SPADES_SOMETIMES_WINS.md](WHY_SPADES_SOMETIMES_WINS.md).
* **Preset switching by library type** — `aggressive` looked better than `--organism` in the
  overlapping-mate band by medians (178,905 vs 162,059) but the paired test is 7/8, p=0.589.
  There is no support for an automatic preset rule, and the median difference was a
  small-sample artefact.
