# TesserACT

A de novo short-read genome assembler: de Bruijn graphs over a multi-k ladder, with
paired-end repeat resolution and an optional genus model for the junctions no read pair
can span.

Built for bacterial isolates. Reads FASTQ or FASTA, gzipped or plain, paired, interleaved
or single-end. No dependencies beyond zlib and a C++17 compiler.

---

## Benchmark

666 *Klebsiella pneumoniae* isolates, every one with a closed reference genome. Both
assemblers ran on raw untrimmed reads with default settings — SPAdes 4.3.0 with its own
error correction on and automatic k selection, TesserACT with no model. Scored by
QUAST 5.2.0 in a single invocation per strain, so no comparison crosses tool versions.
Paired Wilcoxon signed-rank; W/L counts are per strain.

| Metric | TesserACT | SPAdes | W / L | |
|---|---|---|---|---|
| NGA50 | **262,365** | 252,009 | 389 / 275 | **win** (p=7e-4) |
| NG50 | **265,589** | 253,289 | 384 / 280 | **win** (p=8e-4) |
| Largest contig | 614,847 | 627,713 | 386 / 279 | **win** (p=0.03) |
| Genome fraction | **99.2 %** | 98.9 % | 626 / 39 | **win** (p=7e-96) |
| Mismatches /100 kb | **0.5** | 0.8 | 484 / 177 | **win** (p=5e-20) |
| LGA50 | 8 | 8 | 261 / 225 | tie (p=0.67) |
| Contigs | 79 | 79 | 284 / 372 | loss (p=8e-4) |
| Misassemblies | 1 | 0 | 151 / 218 | loss (p=0.002) |
| Duplication ratio | 1.003 | 1.000 | 9 / 594 | loss (p=8e-94) |

Read that honestly: TesserACT reconstructs **more** of the genome and gets **more of the
bases right**, in **longer** contigs. It pays for that with slightly more redundant
sequence and slightly more misassemblies. It is not better in every way, and the table
above is the whole table, not the flattering half of it.

Two rows need care, because the columns disagree:

* **Largest contig.** The medians favour SPAdes (627,713 against 614,847) while the paired
  comparison favours TesserACT (386 strains better, median gain +489). Both are true. The
  win/loss column asks "on this isolate, which assembler did better?" and the median column
  asks "what does a typical assembly look like?" They part company when the two tools lose
  on different strains, which is exactly what happens here.
* **Contigs.** The medians are equal at 79, yet SPAdes wins 372 strains to 284. A minority
  of isolates carry a large share of that loss.

Per-strain numbers, not just summaries, are what these were computed from — ask if you
want them for a comparison of your own.

---

## Install

```sh
git clone https://github.com/iowa69/TesserACT.git
cd TesserACT
make -j                     # needs zlib headers; on conda, CPATH=$CONDA_PREFIX/include
make test                   # 21 end-to-end checks on synthetic genomes
```

Produces `./tessera`.

## Quick start

```sh
# Paired-end isolate, all cores
tessera -1 reads_R1.fq.gz -2 reads_R2.fq.gz -o out/

# Interleaved, 8 threads
tessera --12 reads.fq.gz -o out/ -t 8

# Maximum contiguity, accepting more misassembly risk
tessera -1 R1.fq.gz -2 R2.fq.gz -o out/ --mode aggressive

# Maximum caution: more simplification passes, stricter joins
tessera -1 R1.fq.gz -2 R2.fq.gz -o out/ --mode careful

# With a genus model, for junctions no fragment spans
tessera -1 R1.fq.gz -2 R2.fq.gz -o out/ --organism klebsiella --model kleb.tsm

# Hand it a QC report from scepter (see below)
tessera -1 R1.fq.gz -2 R2.fq.gz -o out/ --qc sample.json
```

Outputs `contigs.fasta`, `assembly_graph.gfa`, and `report.html` / `report.json` —
per-k-rung statistics, the coverage histogram, and every decision the run made.

---

## How it works

**Multi-k ladder.** Small k keeps the graph connected where coverage is thin; large k
separates repeats. TesserACT does both, carrying the contigs of each rung forward into the
next as trusted sequence. For 2×250 reads the ladder is
`21, 33, 55, 77, 87, 99, 111, 119, 127`; shorter reads get proportionally shorter ladders.

The rung *spacing* matters as much as the ceiling. An earlier ladder stepped 77 → 127 in
one jump of 50 while every other step was 12–22, and closing that gap is worth
97 better / 49 worse on NGA50 across the cohort.

**Abundance cutoff.** Chosen per rung from the k-mer count histogram, deliberately
permissive: erroneous k-mers that slip past a count threshold are removed later by graph
topology, which can see structure a bare count cannot. Cutting hard here costs real
sequence and saves little.

**Graph simplification.** Tip removal, bubble popping, erroneous-connection removal and
weak-link pruning, over repeated rounds whose thresholds ramp up as the graph settles.

**Paired-end resolution.** Reads are anchored to unitigs at k=31 — short enough that a
single sequencing error does not invalidate every k-mer in the read, which is what
happens at the graph's own k. Branches are then resolved by the fragment lengths the
pairs imply.

**Genus model (optional).** At junctions no fragment can span, a model built from closed
genomes of the same genus supplies the ordering evidence the reads cannot. Without one,
those junctions are left broken rather than guessed at.

---

## Working with scepter

[scepter](https://github.com/iowa69/scepter) is the companion QC and preprocessing tool.
Run it first and hand the report to the assembler:

```sh
scepter -i R1.fq.gz -I R2.fq.gz --qc-only --preset wgs-bacteria -j sample.json
tessera -1 R1.fq.gz -2 R2.fq.gz -o out/ --qc sample.json
```

The report carries read depth, an estimated genome size, a substitution rate measured by
comparing overlapping mates against each other, and the insert-size distribution — all
available *before* the first k-mer is counted. Measured against the 666 closed references,
the genome-size estimate lands within 10 % on 620 of them and the depth estimate tracks
truth at Spearman 0.95.

What that currently buys you: the run reports what the library actually is, and the
assembler will refuse a QC file that does not match its reads. The decisions those numbers
could drive are implemented but **off by default**, because none of them has yet beaten
what the assembler already infers from its own graph. They are exposed as `TESSERA_QC_*`
environment variables for anyone who wants to experiment. Accurate measurement in search
of a use is worth shipping as exactly that, and not as a feature.

---

## Genus models

`--model` takes a model file and `--organism` names the genus it must match. `make` builds
the model builder alongside the assembler:

```sh
tessera-model --organism klebsiella --out kleb.tsm --layout-tracks \
              --plasmids plasmid_panel.fna --exclude-plasmids withheld.txt \
              --marker-density 64 \
              closed/*.fasta
```

A model records, from closed genomes of one genus, which canonical 31-mers occur near each
other and on what kind of replicon. Markers are sampled by hash at a fixed rate, so the same
loci are picked in every genome and in the assembly under test. `--layout-tracks` also stores
each panel chromosome's marker order whole, which lets an assembly be ordered against the
relative it most resembles instead of joined junction by junction. `--plasmids` adds a
plasmid panel, from which the replicon calls below are derived.

**Density.** `--marker-density N` keeps one k-mer in N; the default is 512, about 500 bp
apart. The value is written into the model and applied automatically when it is queried, so
build and query cannot disagree — which matters because they fail silently rather than
loudly: sampling is by hash threshold, so a query at the wrong density shares almost none of
the model's markers and simply reports nothing.

Denser is not strictly better, and the trade is measured. On 104 *Klebsiella* isolates with
closed references and a leave-cluster-out model:

| | 1 in 512 | 1 in 64 |
|---|---|---|
| multi-contig plasmids delivered in one group | 13.8% | **21.3%** |
| plasmid grouping completeness (per-isolate median) | 0.215 | **0.306** |
| chromosome in one contig at ≥ 90% | **98.1%** | 94.2% |
| model size / build time | 339 MB / 4 min | 2.4 GB / 20 min |

512 is the default because the chromosome result is the one that is finished. 64 buys plasmid
grouping, which is marker-starved at 512 — a 1 kb contig expects two markers there and
grouping needs three — and costs five isolates their chromosome out of 104.

The chromosome cost looks like a parameter rather than a law: `kMarkerNeighbours` is 16,
chosen so that at ~500 bp spacing the adjacency table reaches ~8 kb and covers the rRNA
operons and IS elements that defeat paired evidence. At 1 in 64 the spacing is ~62 bp and
those sixteen neighbours reach ~1 kb. Scaling the neighbour count with the density is
untested.

**Leakage.** A model must never contain the genome being assembled, or anything close enough
to stand in for it. `--exclude` drops chromosomes by accession and `--exclude-plasmids` takes
a list for the plasmid panel — and the second is not optional in practice. On one *Klebsiella*
cohort, 46% of test plasmids had a near-identical match in a RefSeq-derived panel, so a model
built without the list scores plasmids against themselves. Exclude whole mash-distance
clusters rather than individual accessions: dropping only the matched genome leaves its
cluster-mates behind, and they carry the same information.

A model is consulted **only** at junctions no read pair spans. It cannot override
evidence the reads provide, so it changes contiguity without changing what the data say.

---

## Replicon calls

With a model carrying a plasmid panel, every contig is called chromosomal, plasmid or
unknown, and the call is appended to the contig **name** — not placed after the space, where
aligners and scorers silently ignore it:

```
NODE_1_length_5456968_cov_34.9404_chr
NODE_10_length_86253_cov_75.7539_plas_1
NODE_22_length_35126_cov_91.8331_plas
NODE_44_length_1203_cov_9.1100_unk
NODE_6_length_24538_cov_96.8381_plas_8_circular
```

| suffix | meaning |
|---|---|
| `_chr` | chromosomal |
| `_plas` | plasmid, molecule unknown |
| `_plas_<n>` | plasmid, grouped with the other contigs carrying the same `<n>` |
| `_unk` | no signal reached it — reported as unknown rather than guessed at |
| `_circular` | its two ends are joined by read pairs: the contig is the whole molecule |

Four independent signals decide it, none sufficient alone: layout placement on a chromosome
track, panel markers, depth against the modal coverage, and read pairs running between
contigs. Group numbers are per-isolate — `_plas_1` in two genomes is not the same plasmid —
and are assigned by total group length, so `_plas_1` is the largest molecule in that isolate.

**The file is ordered as a genome**, not as a length ranking: the chromosome first, then each
plasmid molecule whole and contiguous, then the plasmid contigs whose molecule is unknown,
then the unassigned. Longest first within a block, with the sequence itself as the final
tie-break so the ordering is total and the numbering reproducible.

What this is worth, on 666 *Klebsiella* isolates with closed references and leave-cluster-out
models: the chromosome lands in a single contig at ≥90% for 95–96% of them, and on 50 unseen
clinical isolates 48 of 50 return a chromosome in one contig of ≥5.0 Mb.

Plasmids are a weaker result and should be read as one. 14 of 60 assembled into a single
contig at all; grouping holds 9 of 65 multi-contig plasmids in one group; a `_circular` tag
is correct 10 times in 11. The groups that are reported are pure (per-isolate homogeneity
1.000) and fragmented (completeness 0.215). The limit is marker density rather than the
panel: at one 31-mer in 512, a 1 kb contig expects two markers and grouping needs three.

---

## Options

| | |
|---|---|
| `-1 / -2 / --12 / -s` | input reads |
| `-o DIR` | output directory |
| `-t N` | threads (default: all cores) |
| `-k LIST` | override the k ladder, e.g. `21,33,55,77` |
| `-c N` | force the abundance cutoff (default: automatic) |
| `--mode NAME` | `fast`, `standard` (default), `careful`, `aggressive` |
| `--organism` / `--model` | genus model, see above |
| `--qc FILE` | scepter QC report |
| `--min-contig N` | shortest contig to report (default 2k) |
| `--max-memory GB` | counting-table budget (default 80 % of RAM) |

`tessera --help` lists the rest, including the repeat-resolution and polishing knobs.

---

## Limitations

* **Short reads only.** No long-read or hybrid support.
* **Repeats longer than the top k stay unresolved** unless paired reads span them or a
  model covers them. In *Klebsiella* that means rRNA operons (~5 kb) and IS elements
  (1–2.5 kb). Raising k further does not help — 99.7 % of the genome is already unique at
  k=99, and the rest is far longer than any k a 250 bp read can support.
* **Duplication ratio runs ~0.3 % above SPAdes** (1.003 against 1.000), consistently and
  across every configuration tried. The cause is not yet identified.
* **Bacterial isolates.** Metagenomes and eukaryotes are untested.

## Licence

See `LICENSE`.
