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
