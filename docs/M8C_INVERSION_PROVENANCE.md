# M8C — Is the layout score detecting assembly error, or biology?

**Answer: neither. It is detecting genuine ambiguity. 24 of 26 inversion events (92.3%, Wilson
95% CI 75.9–97.9%) are flanked by inverted repeat copies that exceed the library's anchor reach
by one to two orders of magnitude, so both orientations are exactly consistent with the data.
Under the pre-registered rule (A > 70%), the flag ships as an ambiguity / structural-variant
warning worded "large-scale orientation not determined by the data" — never "assembly likely
incorrect".**

## Control gate — passed, so silence means something

An undetectable discriminator and a genuine absence of evidence produce identical silence. The
experiment was gated on telling them apart, and the gate passed:

- injected inversions recovered **3/3**, breakpoints within +3..+8 and −11..−12 bp of the
  injection sites;
- the read test was decisive at **6/6** injected breakpoints, median 26 qualifying fragments;
- **power symmetry**: reads simulated *from the injected assembly* flip all 6 to `asm`, so the
  test is not a "reference" stamp;
- false positives **0/20** (Wilson 0.0–16.1%) — and notably 12 of 20 alternative junctions did
  receive crossing alignments, with maximum two-sided reach **3 bp**. The 50 bp anchor rule
  rejected every one. That is the leakage failure mode, caught.

## The result

| class | meaning | n | % | Wilson 95% CI |
|---|---|---|---|---|
| **A** | repeat-flanked, no discriminating read possible → **undecidable from this library** | **24** | **92.3%** | **75.9–97.9%** |
| B | unique sequence, reads support the reference → assembly error | 2 | 7.7% | 2.1–24.1% |
| C | unique sequence, reads support the assembly | 0 | 0.0% | 0.0–12.9% |
| D | repeat-flanked but reads decide | 0 | 0.0% | 0.0–12.9% |
| E | blocks span >1 reference record → replicon misassignment | 0 | 0.0% | 0.0–12.9% |

At breakpoint level (n = 52):

| quantity | value | Wilson 95% CI |
|---|---|---|
| repeat-flanked breakpoints with **any** discriminating fragment | **0 / 39** | 0.0–9.0% |
| breakpoints supporting the assembly orientation, anywhere | **0 / 52** | 0.0–6.9% |
| qualifying fragments, reference vs assembly | **378 vs 0** | — |
| **median two-sided anchor reach** | **0 bp** | — |

The arithmetic predicted this before the experiment ran: ~48 bp of pair reach against rRNA
operons of ~5,000 bp and IS elements of ~800–2,500 bp. Measured shortfalls run from 28 bp to
4,687 bp. The closest miss in the whole study was **58 bp short**.

## The plasmid hypothesis — tested, and not supported at this scale

The experiment was extended to align against `ref_full` (chromosome **and** plasmids) rather than
the chromosome-only `ref_chr` the original labels used, on the theory that plasmid contigs with
nowhere correct to land would migrate to homologous IS copies and be miscounted as inversions.

**That did not happen.** Class E = **0/26**. All 26 minority-strand blocks are chromosome-internal;
of 298 primary ≥100 kb blocks in the cohort, 38 are plasmid and every one stays inside a single
plasmid record. **At ≥100 kb, aligning to `ref_full` changes no inversion call.**

It is *not* refuted below that scale. The sensitivity ladder gives 2/108 strains with
multi-replicon contigs at ≥50 kb, 15 at ≥5 kb, and **32/108 (29.6%, 21.8–38.8%) with no block
floor**. Chromosome–plasmid chimerism is real in these assemblies and sits entirely beneath the
block size the inherited method can see.

## What the preflight changed

**P1 — the scorer is order-sensitive, and two earlier claims were over-scoped.** All 20
single-contig flips changed V (median |ΔV| 0.0129, max 0.154); the earlier "provable no-op,
residual 6e−15" was correctly scoped to *whole-molecule* reverse complement only. Two corrections
follow. That 6e−15 **does not reproduce** on delivered multi-contig layouts — whole-assembly RC
gives ΔV up to 1.6e−2, because the concatenation length is not a multiple of the 5 kb window, so
anyone comparing layouts at the 1e−2 level must pad or circularly window first. And magnitude
tracks contig length: a plasmid-sized contig moves V by roughly the noise floor, so **V cannot
adjudicate plasmid-contig orientation at all**.

**V is a weak detector and must not be default-on.** Sensitivity **8/20 = 40%** (21.9–61.3%),
precision 8/8 (67.6–100%), with a 12/100 inversion rate among unflagged controls. The largest
missed inversions carry 1.0–1.9 Mb of minority-strand sequence at V = 0.89–0.93. V may flag; it
may equally not.

**P2 — the original labels have three defects, all confirmed from source.** They were computed
against `ref_chr` (`grep ref_full` over every `.py` returns zero hits); they aligned **only the
single longest contig**, which produces false negatives by construction — ERR5056647's inversion
lives on `probe_2`, never submitted to the aligner; and `inv_frac` is a strand-mass ratio that
never looks at coordinates and silently scores whole-contig reverse complement as 0.0 (14 of 129
strains). No "inverted" threshold exists anywhere in the code, and `inv.tsv` is read by nothing.
The 129-strain file is a sample, not a flag set.

## Blind re-derivation: 9/10 concordant, and the disagreement was the blind agent being right

An independent agent, shown no prior calls and instructed to default to "cannot reproduce",
re-derived ten events from raw data with its own alignments, annotations and read mapping.
Concordance **9/10** (59.6–98.2%), `cannot_reproduce = 0`.

The single disagreement (ERR11578845, main = A, blind = B) is adopted as **B**. The main pipeline
annotated the breakpoint as a 2,455 bp IS element and mechanically applied the "anchor ≥50 bp
outside the repeat" rule, demanding a 2,554 bp span from a 677 bp library — a **forced zero**. But
its own annotation recorded zero partner copies: the element is **single-copy**, therefore not a
repeat, and cannot mediate recombination. The blind path mapped reads to the *actual delivered
contig junction* and found **0 of 1,938 reads** spanning it, against 17 reads and 9 pairs
supporting reference continuity. Where the two differ, the blind test examined the delivered
sequence and the main pipeline reconstructed the junction from alignment coordinates. **The blind
result is the stronger test.**

**A third finding, from the blind path only, matters beyond this experiment.** The delivered files
are **scaffolds**: 23–84 N-runs per strain totalling 40–78 kb, and the ~8 largest runs sit exactly
on the 8 rRNA operons. **Six of twenty junctions examined are pure N-gaps** — the assembly asserts
no sequence there at all, only an ordering. No read can test a junction that contains no sequence,
and calling such a thing "an inversion" overstates what was claimed.

## Decision

**Ship as an ambiguity flag.** Wording: *"large-scale orientation not determined by the data"*.
Never "assembly likely incorrect" — on this evidence that phrasing would be wrong 92% of the time.

**Not promoted to default-on.** The `B > 70%` branch is excluded under every reading, and 40%
sensitivity independently disqualifies it as a gate.

**Per-isolate class annotation ships with the flag**, per the "mixed" clause: the two class-B
isolates (ERR11578845, ERR6293666) are genuine assembly errors and are worth chasing on their own.

**What would resolve class A: one long read spanning the operon.** This is undecidable from *this
library*, not undecidable. A single nanopore or PacBio read crossing a 5 kb rRNA operon settles
every one of the 24 events immediately.

## What this eliminates

- *"V flags assembly errors"* — refuted. 92.3% of what it flags is not an error.
- *"V flags biology"* — not established either; the data cannot distinguish the two states.
- *"Plasmid mis-mapping explains the inversion labels"* — tested and not supported at ≥100 kb.
- *"The orientation flip is a provable no-op"* — over-scoped; true only for whole-molecule RC,
  and even then only when the length divides the window evenly.
- *"chr_best-style coverage is a usable outcome for layout"* — reaffirmed as false; every one of
  these 26 events is invisible to it.
