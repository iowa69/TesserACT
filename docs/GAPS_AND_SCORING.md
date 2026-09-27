# Gaps, N, and how the misassembly count depends on where you split

> **Scope.** Why `contigs.fasta` no longer contains N, what that changed and what it
> did not, and a scoring artifact that made the organism model look four times worse
> than it is. Measured on 30 held-out *S. aureus* isolates with real Illumina reads
> and their own closed references, QUAST 5.3.0.

## The output changed

`contigs.fasta` is now split at every gap and contains **no N at all**. The ordered,
oriented form is kept in `scaffolds.fasta`, and `scaffolds.agp` (AGP 2.1) records
order, orientation and gap length with per-gap provenance. Nothing is lost; the
adjacency is stated in the AGP rather than asserted as sequence.

Only model-guided runs are affected. Without a model the assembler leaves unspanned
junctions broken rather than guessing, so its output was already 0.008% N.

## What that did to the numbers

Nothing, at contig level — which is the point:

| | scaffolded `contigs.fasta`, scored with QUAST `-s` | zero-N `contigs.fasta` |
|---|---|---|
| NGA50 | 440,260 | 440,260 |
| # contigs | 25 | 25 |
| # misassemblies | 3 | 3 |
| genome fraction | 99.084 | 99.084 |
| N per 100 kbp | 0.00 *(only after QUAST split it)* | **0.00 in the delivered file** |

The contig-level numbers were always contig-level. What changed is that the file we
ship now matches the numbers we quote about it, and QUAST reports
`nothing was broken` rather than silently rescuing the comparison.

## The scoring artifact, which is the part worth reading

QUAST `-s` splits an assembly at runs of **10 or more** N. Gaps shorter than that are
left joined, and a wrong adjacency inside one is counted as a misassembly *within a
contig* — which is what it looks like, because at that point it is one.

The organism model emits many such gaps. Splitting at every N run instead of at runs
of ten removes them from the sequence, and with them the misassemblies they carried.
Same assemblies, same reads, same models, same scorer — only the split threshold
differs:

| cohort sum, n=30 | split at N-runs ≥ 10 | split at every N run |
|---|---|---|
| base (no model) | 7 | 7 |
| **model** | **51** | **11** |

**40 of the 51 were sub-threshold gaps.** The `base` column is identical in both,
which is what confirms the comparison rather than the split moved.

This is not a way of hiding misassemblies. A gap is an assertion of adjacency with no
sequence behind it; scoring it as though the flanking sequence were contiguous
measures a claim the assembly never made. The honest reading is that the model's
residual structural error against a no-model baseline is **11 against 7**, not
51 against 7 — and that any comparison of this assembler which uses QUAST's default
`-s` threshold on scaffolded output overstates its misassembly count roughly
fourfold.

Anyone reproducing earlier numbers from this repository should know that
[`REAL_WORLD_RESULTS.md`](REAL_WORLD_RESULTS.md) and
[`MODELS_INTEGRATED.md`](MODELS_INTEGRATED.md) were measured with `-s` on scaffolded
output, so their misassembly columns are on the other side of this artifact.

## What is genuinely unclosable

Zero N was reached by declining to assert gaps, not by filling them. Of the residual
N in a model-guided *S. aureus* assembly:

* **69.4% lies in runs longer than 3,000 bp**, which local gap-filling never attempts
  (`kMaxGapLen` in `src/gapfill.cpp`). Those are rRNA-operon-scale repeats, and no
  150 bp read pair resolves them.
* Of the 30.6% that is attempted, closure recovers roughly 800 bp of real sequence
  per assembly.

"Fill the gaps instead" is therefore not available for most of them. Splitting is not
a workaround for a weak gap-filler; it is the correct representation of an adjacency
the reads do not support.

## What 1.4.0 changed in gaps and in the reports

These fixes are on by default since 1.4.0. `TESSERACT_FIXES=0` restores the 1.3.0 behaviour
described above; the defect ids are those of the combo3 register.

**The statistics describe the files as written (T17).** In 1.3.0 `report.json` and the
summary were computed from the scaffolds before they were split at N-runs and before the
terminal-overlap trim. On one *E. faecium* isolate that overstated `contig_total_length` by
109,288 bp. Now the `contig_*` fields are computed from the `contigs.fasta` records, and
`scaffold_gaps` counts every N-run in `scaffolds.fasta`, a 1-N run included. `report.json`
also gains `trimmed_overlaps` and `trimmed_overlap_bases`, and a `gap_fill` object: gaps seen
and closed, why the rest stayed open (ambiguous, no path, thin read pool, out of budget),
searches truncated or capped, and back-off closures. A `--quiet` run used to keep no record
of any of it. This part is unconditional.

**An N-run now carries the estimated gap (T03, T20).** In 1.3.0 the resolver wrote each
scaffold gap as an N-run standing in for the first k−1 bases of the unitig after it, and
those real bases were lost. Now the gap filler still sees that layout (so its target k-mer
lies in full-depth sequence), and every gap it leaves open is then rewritten:

* the k−1 bases come back;
* the N-run becomes the estimated true gap length;
* a verified exact overlap between the two flanks (at least 4 bases, not contradicted by the
  estimate) becomes 1 N, with the overlap trimmed once;
* an estimate of 0 or less with no verified overlap has no measurable length and is written
  as **100 N**, meaning "length unknown".

The estimate still runs short: T20 is only partly fixed. The 100-N gaps are written as
ordinary N rows in the AGP, not as type U, because gap provenance is not tracked yet.

**The polisher leaves N-runs alone (T02).** In 1.3.0 the polisher could overwrite a short
N-run with read bases, turning a join the reads never spanned into contiguous sequence. It
now skips N positions. T02 and T03 go together: `tesseract-asm` refuses to start with T03
on and T02 off.

**AGP and GFA describe the final records (T26).** In 1.3.0 no W row of `scaffolds.agp` named a
`contigs.fasta` record, because both files were built from records before post-processing.
Now:

* every W row names a `contigs.fasta` record, with its final coordinates;
* a piece below the length floor still sits inside its scaffold, and its W row names
  `<scaffold>_unlisted_<n>`, a component that `contigs.fasta` does not hold;
* the gap evidence is `paired-ends` for resolver gaps and `unspecified` once a model or the
  layout has joined contigs, because gap origin is not tracked;
* GFA P-lines cover the contigs that are an exact walk of graph segments. A GFA1 P-line
  cannot express "an N-gap, then a segment with L bases trimmed", so those P-lines are left
  out and counted on the log as `pLineOverspell`.

**What this does to scoring.** The restored flank bases lengthen the pieces after gaps. Some
of them then cross QUAST's 500 bp floor, which adds to the contig count (INFERRED; it is one
of the contiguity costs of 1.4.0, see the README).
