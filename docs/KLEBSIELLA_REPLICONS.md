# Replicon calls on *Klebsiella pneumoniae*

> **Scope.** TesserACT **with** a genus model carrying chromosome layout tracks and an
> Enterobacterales plasmid panel. This page measures what the model adds beyond contiguity:
> which molecule each contig belongs to, which contigs belong together, and which contigs are
> whole replicons. For the model-free head-to-head against SPAdes see
> [KLEBSIELLA_PANEL.md](KLEBSIELLA_PANEL.md).

A pile of contigs is not what a clinical read-out needs. It needs "this is the chromosome,
these are plasmid 1 and plasmid 2, and this is what did not fit anywhere", with the evidence
attached so a reader can tell an observation from an inference.

The chromosome half of that is finished. The plasmid half is not, and this page says so with
numbers rather than adjectives.

---

## 1. Data and leakage control

**Panel.** 2,933 closed *K. pneumoniae* genomes; 2,221 contribute chromosome layout tracks,
51,789 records contribute the plasmid panel (Enterobacterales plasmids under 300 kb).

**Test.** 666 isolates with paired reads *and* a closed reference for the same DNA — the only
collection where the right answer is known for both chromosome and plasmids.

**Leakage.** Panel and test were clustered together by mash distance under single linkage at
d < 0.001 and **whole clusters** were excluded, not individual genomes: dropping the matched
genome leaves its cluster-mates, which carry the same information. Five leave-cluster-out
folds.

Plasmids needed a separate audit and it changed the design. **46.2% of test plasmids had a
near-identical match in the panel**, affecting 448 of 680 isolates. Without intervention the
plasmid model would have been scored largely on its ability to recognise sequence it had
memorised. 4,942 panel plasmids are withheld per fold.

---

## 2. Chromosome

| | |
|---|---|
| chromosome in one contig at ≥ 90% | **95.2 – 96.4%** of the 666, across folds |
| median chromosome-in-one-contig | 98.8 – 99.0% |
| NGA50 | 2.82 Mb |
| genome fraction | 99.23% |
| misassemblies | 4 |

On 50 clinical isolates that contributed nothing to the model and have no closed reference,
48 of 50 (96.0%) return a chromosome in a single contig of ≥ 5.0 Mb, mean longest/total
0.929. The two measurements are not identical — one is scored against a true chromosome, the
other is a 5 Mb contig inside a 5.7 Mb assembly — but they agree.

---

## 3. Plasmids

### 3.1 Classification

Four independent signals, none sufficient alone: layout placement on a chromosome track,
panel markers, depth against modal coverage, and read pairs between contigs. Contigs no
signal reaches are reported `_unk` rather than guessed at.

### 3.2 Grouping

Grouping asks a *set* question — does one panel plasmid carry markers from both contigs —
because that survives the mosaicism that defeats order-based evidence. Against truth on 5,370
contig pairs it separates same-plasmid from different-plasmid at **AUC 0.857**; a
degree-preserving null of the same graph falls to 0.416, below chance.

The delivered result is weaker than that number suggests:

| 30 fold-0 isolates | |
|---|---|
| multi-contig plasmids delivered in one group | **7 / 44 (15.9%)** |
| per-isolate median homogeneity | 1.000 |
| per-isolate median completeness | 0.215 |

The groups are pure and fragmented. That is one failure, not two, and its cause is a single
number: markers are canonical 31-mers sampled at **1 in 512**, so a 1 kb contig expects two
markers and grouping needs three. The contigs that fail to group have a median length of
about 1.2 kb, with 62.5% under 1.5 kb — they return zero before any threshold is consulted.

Read pairs reach some of them. Adding pair-link propagation cut the ungrouped pile from 260
contigs to 73 and moved the whole-plasmid rate not at all: more contigs get a group, and the
groups stay split.

### 3.3 Closure

| 30 fold-0 isolates | |
|---|---|
| contigs tagged `_circular`, non-chromosomal | 11 |
| of those, complete reference replicons | 10 (**90.9%** precision) |
| plasmids assembled into a single contig | 14 of 60 |
| of those, carrying the tag | 10 (71.4% recall) |

Every correct call overhangs its reference plasmid by 126 bp — k−1 at k=127, the assembler's
junction overlap, sitting exactly where a circularised contig should carry it.

The eleventh call is worth naming rather than rounding away. `ERR11578929 NODE_10`, 5,563 bp
at 45×, does not align to its isolate's closed reference anywhere; it aligns end to end to
panel plasmid `NZ_JABXPX010000008.1` in two blocks meeting at the origin, with the junction
overhang. It is most likely a small plasmid the reference assembly does not carry. It is
scored as a miss anyway — crediting a call because our own panel agrees with it would be
scoring the model against itself.

The recall ceiling is not detection: only 14 of 60 plasmids assembled into a single contig at
all, and the other 46 were never candidates for the tag.

---

## 4. Output

The replicon call is appended to the contig **name**, not placed after the space — an earlier
version put it in the description, where aligners and scorers silently ignored it and two
arms with visibly different tags scored identically.

| suffix | meaning |
|---|---|
| `_chr` | chromosomal |
| `_plas` | plasmid, molecule unknown |
| `_plas_<n>` | plasmid, grouped with the other contigs carrying the same `<n>` |
| `_unk` | no signal reached it |
| `_circular` | its two ends are joined by read pairs: the contig is the whole molecule |

The file is ordered as a genome: chromosome first — in layout order where the layout placed
it, so the pieces read along the reconstructed chromosome rather than by size — then each
plasmid molecule whole and contiguous, then the plasmid contigs whose molecule is unknown,
then the unassigned. Longest first within a block, sequence as the final tie-break, so the
order is total and the numbering reproducible. Group numbers are per-isolate and assigned by
total group length, so `_plas_1` is the largest molecule in that isolate.

---

## 4b. What layout costs, and how to decline it

Every figure above measures how much chromosome lands in one record. That is what this project
set out to deliver, but on its own it says nothing about whether the order inside that record
is right, and layout takes its order from a different genome. Scored by QUAST against the full
reference on the same 105 held-out isolates, with the assembly also re-scored after splitting
it back into contigs at the N runs layout writes:

| | misassemblies | median NGA50 |
|---|---|---|
| SPAdes input | 121 | 242,291 |
| after join, gaps split | 358 | 315,291 |
| after join + layout | 584 | 2,959,825 |

The join stage buys 1.3× contiguity for roughly three times the misassemblies. Layout buys a
further 9.4× for another 226, and **every one of those 226 sits at a gap layout inserted** —
splitting the output at its N runs removes all of them and returns the assembly to 358. The
input has none at all, because it has no gaps.

Two things follow. The first is that these are not silent errors: they are localised at
declared N runs, which are visibly not sequence, and a gap's *length* is copied from the panel
genome, so it is wrong by construction rather than by mistake — QUAST charges 640 gap-size
extensive misassemblies on this cohort for that alone (700 for the dense build). The second is that the choice is the
user's and the tool already offers it. `--no-layout` stops before this stage, and splitting a
finished assembly at `N{10,}` recovers the conservative version from an already-laid-out run.

Take layout when you want a chromosome to look at, and decline it when you need order you can
defend per junction — for a structural-variant call, a rearrangement claim, or anything where
a 5 Mb record implying observed contiguity would mislead. The default is on, because for the
question this tool exists to answer it is the right trade; it is documented here so that it is
a choice rather than a surprise.

## 5. What this does not do

- It does not close most plasmids. 14 of 60 assembled into one contig; grouping holds 15.9%
  of the multi-contig ones in a single group.
- It leaves a long tail unassigned — on one clinical isolate, 69 of 103 contigs. They are
  reported as unknown rather than guessed at, which is the right behaviour for a question
  about mobilizable content, but it is not a small residue.
- It does not get both halves at once. Every plasmid limitation above traces to the marker
  sampling rate: at one 31-mer in 512 a 1 kb contig expects two markers and grouping needs
  three. Sampling eight times denser is the fix, and it works — whole plasmids 13.8% → 21.3%,
  grouping completeness 0.215 → 0.306, per-isolate homogeneity unchanged at 1.000, so the
  extra grouping is not over-merging. It also costs the chromosome: over the same 104
  isolates, chromosome-in-one-contig at ≥90% falls from 98.1% to 94.2%, with five isolates
  dropping below the line and one rising above it. The falls are steep (99.2 → 86.9) and
  share a signature — the layout places fewer contigs and displaces some.

  That was traced to a parameter that did not scale. `kMarkerNeighbours` is 16, chosen so
  that at ~500 bp spacing the adjacency reaches ~8 kb and covers the rRNA operons; at ~62 bp
  spacing the same 16 reach ~1 kb. Holding the reach constant instead — 1 in 128 with 64
  neighbours, again ~8.2 kb — recovers four of the five broken isolates exactly, and beats
  the 1-in-64 build on both halves:

  | | 1/512, 16 nb | 1/64, 16 nb | 1/128, 64 nb |
  |---|---|---|---|
  | plasmids delivered whole | 9/65 (13.8%) | 13/61 (21.3%) | **15/63 (23.8%)** |
  | grouping completeness, per isolate | 0.219 | 0.306 | **0.325** |
  | grouping homogeneity, per isolate | 1.000 | 1.000 | 1.000 |
  | chromosome ≥90% | **98.1%** | 94.2% | 95.2% |
  | chromosome ≥98% | 64.8% | 59.6% | **69.5%** |
  | median chr_best | 98.6 | 98.5 | **98.7** |

  The ≥90% line still favours the default: the reach-corrected build breaks a different set
  of four isolates while rescuing one badly broken one, so the tail moved rather than shrank.
  512 with 16 neighbours stays the default; 1/128 with 64 is the build for anyone who wants
  the plasmid half.

  Two candidate mechanisms for the residual were then tested, and both are closed. Neither is
  the overlap-rejection rule, which the displaced counts made the obvious suspect: displaced
  placements do rise with density, but the isolates that lost closure have a median of 3
  against 2 for those that did not, and the worst loser of all has 2. The other was the
  scaffold split threshold `kMaxTrackGap`, and testing it turned up something that matters
  more than the parameter. Raising it recovers the split isolates on `chr_best` — 80.3 → 98.1
  on one — while leaving NGA50 and largest alignment identical to the base pair and adding
  scaffold relocations. The recovery was the same aligned bases relabelled from three records
  into one. The constant stays at 400 kb.

  That result applies to the table above, because `chr_best` and the ≥90% line are the same
  kind of measure and the two arms do not scaffold equally: the default places 12.44 contigs
  per isolate, the reach-corrected build 10.04. Scored by QUAST over all 105 held-out
  isolates, paired, the arms are close to parity — the dense build is more complete on 92 of
  them and more contiguous on 60, the default has fewer misassemblies on 48 against 30 (584
  against 674 cohort-wide), and median NGA50 differs by 1,547 bp on a 5.3 Mb chromosome. The
  ≥90% row overstates a gap that alignment shows as a trade. Quote the two together.
- Hub grouping over a clustered panel does not rescue it either. Clustering the 56,731-plasmid panel
  by mash single linkage at d < 0.01 — the last threshold before percolation; the largest
  component holds 2.1% there against 29.6% at d < 0.02 — gives 16,495 clusters, every panel
  plasmid in the model maps to one, and contigs start naming the same cluster where before
  every argmax was unique. Across 40 isolates it moves one contig. Clustering changes *which*
  panel entry a contig names, not *whether* it can name one. The pass ships behind
  `TESSERACT_HUB_GROUP`, off.
