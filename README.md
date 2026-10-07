# TesserACT

A de novo short-read assembler for bacterial isolates: de Bruijn graphs over a multi-k
ladder, paired-end repeat resolution, verified circular replicons, and optional organism
models for the junctions no read pair can span. It reads paired FASTQ (gzipped or not) and
needs nothing beyond zlib and a C++17 compiler.

Its priority is **correctness first**: on held-out isolates with closed reference genomes it
makes fewer misassemblies and far fewer base errors than SPAdes 4.3.0, at equal or better
contiguity (see [Accuracy](#accuracy)).

---

## Quick start

```sh
tesseract -1 sample_R1.fq.gz -2 sample_R2.fq.gz -t 32 -o sample
```

That is the whole command. It runs in three steps and stops at the first one that fails:

1. **Read check.** Both files are complete (gzip integrity), well-formed FASTQ in Phred+33, and
   pair record for record. It measures read length, depth, adapter read-through and
   two-colour poly-G tails, and checks free disk and memory. Problems are caught *before* an
   hour of assembly, with a sentence saying how to fix them.
2. **Assembly**, vanilla: no organism model, no prediction from other genomes. Everything in
   the output was assembled from these reads.
3. **Output check.** Every file exists and parses, the FASTA is valid, the size is plausible,
   and three self-checks are summarised: what the reads contain that the assembly does not
   (k-mer self-check), which replicons are verified circles, and whether the library behaved
   like a normal paired-end library.

The result is in `sample/`:

| file | what it is |
|---|---|
| `contigs.fasta` | **the assembly** -- use this |
| `scaffolds.fasta` | contigs joined across gaps of estimated size (written only when such gaps remain) |
| `SUMMARY.txt` | what was checked, every warning, and what to do about it |
| `preflight.txt` | the read check alone |
| `report.html` / `report.json` | every number and decision of the run |
| `assembly_graph.gfa` | the assembly graph (Bandage-compatible) |
| `replicons.tsv` | each circular candidate and the evidence for or against closing it |
| `ends.tsv` | why each contig ends where it does |
| `p2_self_qa.tsv` | k-mer completeness, QV and sequence missing from the assembly |
| `assembly.log` | the assembler's log |

Exit status 0 means every check passed; warnings are listed in `SUMMARY.txt`. Exit status 1
means a check failed (nothing is assembled if the *read* check fails; add `--force` to try
anyway).

To check reads without assembling:

```sh
tesseract check -1 sample_R1.fq.gz -2 sample_R2.fq.gz
```

### Using an organism model (manual, optional)

```sh
tesseract-get-models                                   # once: downloads the models (~1.7 GB)
tesseract -1 R1.fq.gz -2 R2.fq.gz -t 32 -o sample --organism kpneumoniae
```

Models exist for `kpneumoniae` (`klebsiella`), `ecoli`, `ecloacae`, `abaumannii`,
`paeruginosa`, `saureus`, `efaecium` and `senterica` (`salmonella`). A model is used **only when
you ask for it**. It lays out the chromosome from conserved gene order in thousands of closed
genomes of the species, which raises contiguity (median NGA50 134 -> 156 kb on 329 held-out
isolates) but roughly doubles misassemblies, because the order it uses was observed in other
genomes, not in your reads. Use it when contiguity matters more than certainty; keep the
vanilla assembly as your reference. If the reads look like one of these organisms, the vanilla
`SUMMARY.txt` says so, as a hint.

### Many samples

```sh
for r1 in reads/*_R1.fastq.gz; do
    s=$(basename "$r1" _R1.fastq.gz)
    tesseract -1 "$r1" -2 "reads/${s}_R2.fastq.gz" -t 32 -o "assemblies/$s" || echo "$s: see assemblies/$s/SUMMARY.txt"
done
```

File names may contain spaces. Each sample's `SUMMARY.txt` tells you whether it needs a look.

---

## Install

```sh
git clone https://github.com/iowa69/TesserACT.git
cd TesserACT
./install.sh
```

Run in a terminal, `./install.sh` checks the four tools it needs and names the one command
that installs any that are missing, asks where to put things, builds, puts the commands on
your PATH and offers the optional organism models (the default answer is no). Nothing goes
system-wide and no password is needed. From a script it takes the defaults silently;
`--prefix DIR` or `--conda-env NAME` choose the destination.

You get `tesseract` (the command above), `tesseract-asm` (the assembler itself),
`tesseract-get-models`, and the organism front ends `tesseract-eskape` and
`tesseract-klebsiella`. Optional: `minimap2` on the PATH adds a read map-back to the output
check.

To build without installing: `make -j` (needs zlib headers), then `make check`.

---

## Getting the best assembly from your data

`tesseract` applies the settings that were measured to be best on 273 development isolates
of 30+ species and every common library type; you do not need to tune anything. This section
says what it does for each kind of data, and what the warnings in `SUMMARY.txt` mean.

**Illumina 2x150 (NextSeq, NovaSeq, HiSeq), 30-150x.** The common case, and the one the
defaults are tuned on. Nothing to set.

**MiSeq 2x250 and 2x300.** The k ladder grows with the read length automatically (up to
k = 127). Longer reads resolve short repeats that 2x150 cannot, but the inserts of these
libraries are usually *shorter* (278-327 bp measured), so long repeats (rRNA operons, IS
elements) are no more resolvable than with 2x150. No special setting improved them on the
panel: a stronger connector cleaning that helped one 2x301 project made no difference on 41
long-read panel isolates and cost misassemblies on 2x150, so it is not used.

**Low depth (below 30x).** Reported as a warning; below 15x as a strong one. Expect more
breaks where coverage drops (AT- or GC-rich stretches). If a very thin library (around 10x)
assembles almost nothing with the default k ladder, `tesseract` retries automatically with
k = 21,33,55, keeps the first attempt in `attempt1_default_k/`, and says so in the summary. At
normal depth the short ladder is worse, so it is never used there. More sequencing is the only
real fix.

**Very deep runs (above ~150x).** Handled: depth beyond ~150x does not improve the result, and
runtime and memory grow with it; above 400x the read check warns. Small multicopy plasmids --
ColE-type plasmids present at 100-500 copies per cell, i.e. thousands of x -- are a special
case: their reads are thinned to about twice the genome depth before assembly, so their
sequencing errors no longer swamp the graph. This recovers small replicons that earlier
versions lost (see [What 1.6 changes](#what-16-changes)). The k-mer self-check reports any
high-copy sequence that is still missing.

**Two-colour instruments (NextSeq, NovaSeq).** Poly-G tails from dark cycles are measured; above
2% of read pairs the check suggests `fastp --trim_poly_g`.

**Untrimmed adapters.** Above 5% of read pairs with adapter read-through the check suggests
trimming (fastp, cutadapt). TesserACT trims low-quality 3' ends itself but not adapters.

**Mate-pair or unusual libraries.** The library-orientation check measures pair orientation
on the assembly graph. If the library is not a standard forward-reverse paired-end library
(more than half the pairs point outward, or fewer than a quarter fit an FR insert), pair-based
joins and circle calls are switched off for that run and the summary says so. On the
development panel one mate-pair library was the source of 54 of 59 false scaffold joins and
17 false circles; this check removes them and changes nothing on normal libraries.

**PhiX spike-in.** Contigs that are PhiX174 are labelled `_spikein` and never called plasmid;
the summary tells you to remove them before submission. A full PhiX contig appeared in 30 of
508 assemblies of public data.

**Mixed or contaminated samples.** An assembly much larger than the expected genome size is
flagged (give `--organism` or `--genome-size` for an exact expectation).

**When the chromosome must be in as few pieces as possible.** Use an organism model
(`--organism`, above), or `--mode aggressive`, which collapses diverged repeat copies (about one
extra misassembly per genome). Both trade certainty for contiguity; the vanilla assembly does
not.

**Closing the chromosome.** With short reads alone, about half of the remaining breaks in a
typical *Klebsiella* chromosome sit in long identical repeats (IS elements, rRNA operons) that
no read pair spans: no short-read assembler can place them from the reads. Relatives conserve
the *order* of these loci but not their content, which is what the organism models use. A
modest long-read set (about 20x) would decide 95% of these junctions.

**Reproducibility.** Output is identical across thread counts. `report.json` records the
binary, the input file checksums and every non-default setting.

---

## Reading the output

### Contig names

```
NODE_1_length_5456968_cov_34.9404_chr
NODE_10_length_86253_cov_75.7539_plas_1
NODE_6_length_24538_cov_96.8381_plas_circular
NODE_31_length_5386_cov_412.1_unk_spikein
```

| suffix | meaning |
|---|---|
| `_chr`, `_plas`, `_plas_<n>`, `_unk` | replicon call: chromosome, plasmid (grouped by `<n>` when known), unknown |
| `_circular` | a **verified** circle: the graph closes it, it is isolated, and reads span the join. The sequence is written once, without the duplicated (k-1)-bp end older versions carried |
| `_spikein` | PhiX174 control DNA |

Lowercase bases at a contig end are bases the reads at that end do not support well (see
`ends.tsv`); they are kept, not trimmed, so you can decide.

### SUMMARY.txt

Each line is `OK`, `WARN` or `FAIL` with a sentence. The ones worth knowing:

| line | meaning | what to do |
|---|---|---|
| `self-check ... not in the assembly` | sequence present in the reads at >10x genome depth is missing (typically a small high-copy plasmid, rRNA or IS copies) | `p2_self_qa.missing_hi.tsv` lists it |
| `circles ... written linear` | a circle candidate failed verification | it is written linear; the reason is in `replicons.tsv` |
| `library ...` | the pairs do not behave like a paired-end library | check the library type; pair-based joins were switched off |
| `size ...` | total length far from the expected genome size | contamination, a mixed sample, or the wrong expectation |
| `map-back ...` (with minimap2) | fewer than 95% of reads map back | some of the genome is not assembled |

---

## Options

### tesseract

| option | |
|---|---|
| `-1`, `-2 FILE` | paired reads |
| `-o DIR` | output directory (must not already hold an assembly) |
| `-t N` | threads (default: all cores) |
| `--organism NAME` | use this organism model (manual; see above) |
| `--genome-size MB` | expected genome size for the depth and size checks |
| `--force` | assemble even if the read check fails |
| anything else | passed to `tesseract-asm` unchanged, e.g. `--mode careful`, `--min-contig 500` |

`tesseract check -1 R1 -2 R2` runs the read check only. `tesseract --version` prints both versions.

### tesseract-asm

`tesseract-asm` is the assembler itself, for pipelines that want direct control. On its own it
keeps the 1.5.0 defaults exactly (byte-identical output); the 1.6 additions are switched on by
the environment settings below, which `tesseract` sets for you.

| option | |
|---|---|
| `-1 / -2 / --12 / -s FILE` | paired, interleaved or single-end reads |
| `-o DIR`, `-t N` | output directory, threads |
| `--organism NAME` | organism model (`TESSERACT_MODEL_DIR` says where the `.tsm` files are) |
| `--mode NAME` | `fast`, `standard` (default), `careful`, `aggressive` |
| `--min-contig N` | shortest contig reported (default 2k) |
| `-k LIST`, `-c N` | k ladder and abundance cutoff (default: automatic) |
| `--qc FILE` | a scepter QC report of the same reads |
| `--map-polish bwa\|bowtie2` | polish against a full read alignment (`--mapper-dir` where to find the mapper) |
| `--max-memory GB` | abort cleanly above this (default 80% of RAM) |
| `--no-gfa`, `--no-html`, `--unitigs` | output files |

`tesseract-asm --help` lists every option, including the repeat-resolution, simplification,
gap-filling and polishing settings; each of them is exercised by the release test.

### The 1.6 profile

`tesseract` sets these for every run; a value you set in the environment wins, so any one can
be switched off (`=0`, or `TESSERACT_P2_CIRC=off`):

| setting | does |
|---|---|
| `TESSERACT_P2_CIRC=verify` | close circles in the graph, verify the join with reads, write verified circles without the duplicated end and failed ones as linear |
| `TESSERACT_P2_LIBGUARD=1` | library-orientation check (no pair joins for non-FR libraries) |
| `TESSERACT_P2_SPIKEIN=1` | PhiX screen and `_spikein` label |
| `TESSERACT_DEEP_NORM=30`, `TESSERACT_DEEP_NORM_MATE=both` | thin read pairs from replicons at >30x the genome depth (small multicopy plasmids) |
| `TESSERACT_P2_ENDS=1`, `TESSERACT_P2_TIPS_LOWERCASE=1` | end audit (`ends.tsv`) and lowercase unsupported tip bases |
| `TESSERACT_P2_SELF_QA=1` | k-mer self-check (`p2_self_qa.*`) |
| `TESSERACT_P2_PROVENANCE=1` | input checksums, binary and settings in `report.json` |
| `TESSERACT_P2_DETECT_REPORT=1` | organism detection, reported only (never applies a model) |

---

## What 1.6 changes

**One command, with checks.** `tesseract -1 R1 -2 R2 -t N -o OUT` checks the reads before
assembling (truncated gzip, broken FASTQ, unpaired or misordered mates, Phred+64 qualities,
adapters, poly-G, depth, disk, memory), assembles, and checks every output, with a plain-language
`SUMMARY.txt`. The organism model is applied only when you ask for it (`--organism`). For very
thin data it retries automatically with a short k ladder.

**Verified circular replicons.** A circle is claimed only when the graph closes it, the component
is isolated, and reads span the join; the duplicated (k-1)-bp end that 1.5.0 (and SPAdes) wrote
is removed. On 273 development isolates with closed references: circle claims exact as written
2 -> 81, false circle claims 65 -> 3, claim precision 0.63 -> 0.97.

**Small multicopy plasmids.** Reads of replicons at more than 30x the genome depth (ColE-type
plasmids at hundreds of copies) are thinned before assembly, so their sequencing errors no
longer turn them into a tangle. On 133 development isolates: small replicons (< 10 kb)
assembled whole 28 -> 38, present 63 -> 77. The same step stops a failure on thin libraries
where an ultra-abundant contaminant set the abundance cutoff: a ~10x *S. aureus* library went
from 0 contigs to a 2.71 Mb assembly.

**Labels and self-checks.** PhiX spike-ins are labelled `_spikein`; non-paired-end libraries
(mate-pair, wrong orientation) are detected and their pair joins disabled; contig ends are
audited (`ends.tsv`, unsupported tip bases in lowercase); a k-mer self-check reports sequence
present in the reads but missing from the assembly; `report.json` records input checksums,
binary and settings; organism detection is reported as a hint.

**Correctness unchanged.** Misassemblies are identical on every development panel, isolate for
isolate (108 diverse isolates: 350 = 350; 189 ESKAPEE: 935 = 935; 27 *Salmonella*: 4 = 4;
133 isolates with plasmid recovery: 354 = 354), and NGA50 is never lower. Through the full
`tesseract` command, 13 isolates covering 2x100 to 2x300 reads and 8x to 258x depth give the
same misassemblies (4 = 4) and NGA50 as 1.5.0, with more genome where small plasmids were lost.
`tesseract-asm` on its own still produces the 1.5.0 output byte for byte.

What did not change: the organism models (models-v2) and their cost. A redesigned model that
would order chromosomes from relatives was tested against true junctions and rejected for
this release (it added 0.4 misjoins per isolate).

---

## Accuracy

Held-out isolates with closed reference genomes, QUAST 5.3.0, contigs >= 500 bp. "Correct
reference" excludes isolates whose public reference is a different strain from the reads.

| | TesserACT 1.6 (same as 1.5.0) | SPAdes 4.3.0 |
|---|---|---|
| ESKAPEE, 140 fresh isolates: misassemblies (correct reference) | 48 | 114 |
| ESKAPEE, 140 fresh isolates: mismatches per 100 kb | 0.5 | 2.1 |
| ESKAPEE, 140 fresh isolates: median NGA50 | 146 kb | 144 kb |
| 82 isolates of 30 other species: misassemblies (correct reference) | 31 | 85 |
| *Salmonella*, 30 isolates: misassemblies | 4 | 10 |

Where SPAdes still does better: it keeps rRNA operons and IS copies as collapsed consensus
contigs (TesserACT leaves them in the graph: full-length 16S in 38 vs 270 of 329 assemblies),
and it closes a few more single-read breaks in very low-coverage stretches. Full benchmark
history, earlier release notes and design notes: [docs/HISTORY.md](docs/HISTORY.md).

---

## Limitations

* **Short reads only.** No long-read or hybrid assembly.
* **Long identical repeats stay open** unless an organism model is used: rRNA operons
  (~5 kb) and IS elements (1-2.5 kb) longer than the library's insert size. No short-read
  assembler can place them from the reads alone.
* **Bacterial isolates.** Metagenomes and eukaryotes are untested.
* **Organism detection** covers the seven ESKAPEE species; for *Salmonella* and others give
  `--organism` or `--genome-size` yourself.
* **Gap fills from a capped search (open defect N8).** Rare (5 fills in 248 development
  isolates); counted on the `[gapfill-budget]` line of `assembly.log` as `cappedAccepted=`.

## Citation and licence

See `CITATION.cff` and `LICENSE`.
