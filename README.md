<div align="center">

<img src="docs/logo.svg" alt="TesserACT" width="132" height="132">

# TesserACT

**A de novo genome assembler for bacterial isolates.**
Short reads in, chromosomes out — in ~8,600 lines of C++17 with no dependencies
beyond a compiler, zlib and pthreads.

[Install](#install) · [Quick start](#quick-start) · [How it works](#how-it-works) · [Benchmark](#benchmark) · [Options](#options)

</div>

---

## What it does

TesserACT assembles haploid bacterial isolates from paired-end Illumina reads. It
builds a de Bruijn graph over a ladder of k-mer sizes, cleans it, resolves
repeats with paired-end evidence, closes gaps and polishes — then writes the
contigs, the assembly graph, and a report of what every stage decided.

**It declines to guess.** Where a repeat is longer than a sequencing fragment
there is no paired evidence that can settle it, and tessera stops there rather
than joining on a hunch. That is why it makes about half the misassemblies and a
third the per-base errors of the usual alternative.

**And it can do better than stop.** Bacterial genomes of one species are not
arbitrary strings — the sequence flanking a repeat is largely conserved. Give it a model built from closed genomes of the same organism and it settles
those junctions from independent evidence. On 40 *Klebsiella pneumoniae*
isolates with complete reference genomes it then makes fewer structural errors
than SPAdes, is about twice as accurate per base, recovers more of the genome,
and draws level on contiguity.

| On 40 isolates with closed references | TesserACT + model | SPAdes |
| --- | --- | --- |
| Median contig NGA50 | 304,672 | **319,598** |
| Median NGA50 ratio to SPAdes | **1.02x** | 1.00x |
| Isolates with the longer NGA50 | **25 / 40** | 15 / 40 |
| Total misassemblies | **21** | 26 |
| Mismatches / 100 kbp | **0.23** | 0.44 |
| Genome fraction | **99.23%** | 98.95% |

Contiguity is a tie, and the two rows that look contradictory are both real:
TesserACT has the longer contigs on most isolates, while SPAdes wins the median
because the isolates it wins, it wins by more. Every figure is contig level and
comes from [`docs/final_headtohead.tsv`](docs/final_headtohead.tsv).

## Install

**With conda** — the easiest route, and it brings the optional tools too:

```sh
git clone https://github.com/iowa69/TesserACT.git && cd TesserACT
./install.sh --conda-env tesseract
conda activate tesseract
```

**Without conda** — nothing is needed but a C++17 compiler and zlib:

```sh
git clone https://github.com/iowa69/TesserACT.git && cd TesserACT
./install.sh                      # into the active conda env, else ~/.local/bin
./install.sh --prefix /usr/local  # or wherever you like
```

**Or just build it:**

```sh
make -j
./tessera --help
```

The optional extras — `mash` for choosing a model panel, `quast` for scoring
against a reference, `bandage` for looking at the graph, `bwa`/`bowtie2` for
alignment polishing — are listed in `conda/environment.yml`.

## Quick start

Assemble a pair of FASTQ files:

```sh
tessera -1 reads_1.fq.gz -2 reads_2.fq.gz -o assembly
```

That writes `assembly/contigs.fasta`, `assembly/assembly_graph.gfa` and
`assembly/report.html`. On a 5 Mb genome at 60x it takes about a minute on a
laptop and around 3 GB of memory.

Assembling *Klebsiella*, and you have closed genomes to learn from? Build the
model once and reuse it:

```sh
tessera-model --organism klebsiella --out kleb.tsm references/*.fasta
tessera --organism klebsiella --model kleb.tsm \
        -1 reads_1.fq.gz -2 reads_2.fq.gz -o assembly
```

Trimming first is worth it, and `fastplus` has a preset for exactly this:

```sh
fastplus --preset assembly -i R1.fq.gz -I R2.fq.gz --out-dir trimmed/
```

## Design

The pipeline, in the order it runs:

1. **Load reads.** Every read is 2-bit packed and kept in memory for the whole
   run; ambiguous bases get a side bitmap. Mates are stored adjacently, so the
   mate of read `i` is `i^1`.
2. **Read error correction.** A k-mer spectrum pass at the smallest k in the
   ladder: each base that breaks a run of trusted k-mers is replaced with the
   one that restores it (`--no-correct` skips this).
3. **The multi-k loop.** For each k in the ladder:
   - parallel sharded k-mer counting,
   - automatic abundance-cutoff selection from the count histogram,
   - compacted (unitig) graph construction,
   - iterative simplification: tip clipping, bubble popping, erroneous-connection
     removal and isolated-node removal,
   - the resulting contigs are carried into the next, larger k as extra evidence.
4. **Paired-end repeat resolution.** Reads are anchored to unitigs with k=31
   probes, a fragment-length model is learned from pairs landing on one unitig,
   and chains of unitigs are extended by mutual-best scoring of whole paths
   through repeats, with a fragment-length consistency test (`--no-resolve`).
5. **Scaffolding**, which runs as the tail of the same stage. Chain ends with
   strong paired support but no path through the graph are joined across a gap
   of `N`s sized from the fragment model (`--no-scaffold`).
6. **Consensus polishing.** Reads are re-anchored to the finished sequence and a
   per-position majority fixes residual substitutions (`--no-polish`).
7. **Output.** Contigs are written longest-first to `OUTDIR/contigs.fasta`.

k-mers are packed into four 64-bit words (256 bits), so **k can go up to 127**
and must be odd.

## Building from source

Requirements:

* a C++17 compiler (CI builds and runs `make check` under both g++ and
  clang++)
* zlib (`libz-dev` / `zlib-devel`)
* pthreads

Nothing else — no CMake, no Boost, no external k-mer counter.

### Build targets

```sh
make -j                 # optimised build -> ./tessera
make native             # same, plus -march=native (rebuilds from scratch)
make debug              # -O0 -g3
make asan               # AddressSanitizer + UBSan
make test               # end-to-end test suite (needs python3)
make unittest           # C++ unit tests
make check              # unit tests, then the end-to-end suite
make install PREFIX=~/.local
make clean
```

`make install` honours `PREFIX` (default `/usr/local`) and `DESTDIR`;
`make uninstall` removes the installed binary. The build is warning-free with
`-Wall -Wextra`.

## Options

```
INPUT
  -1, --read1 FILE        forward reads (FASTQ/FASTA, optionally gzipped)
  -2, --read2 FILE        reverse reads
      --12 FILE           interleaved paired reads in one file
  -s, --single FILE       unpaired reads

OUTPUT
  -o, --out DIR           output directory (default: tessera_out)
      --min-contig N      minimum contig length to report (default: 2*k)

ASSEMBLY
  -k, --kmers LIST        comma-separated k values, e.g. 21,33,55,77
                          (default: chosen from the read length; max 127)
  -c, --cutoff N          k-mer abundance cutoff (default: auto-detect)
      --min-link N        paired reads needed to trust a join (default: 2)
      --tie-ratio F       winning branch must beat the runner-up by F (default: 1.15)
      --aggressive        also collapse diverged repeat copies: more
                          contiguity, but risks a misassembly per genome
      --no-correct        skip read error correction
      --no-resolve        skip paired-end repeat resolution
      --no-scaffold       skip scaffolding
      --no-gapfill        skip gap closing inside scaffolds
      --no-polish         skip consensus polishing
      --polish-passes N   consensus polishing passes (default: 1)
      --simplify-rounds N graph simplification rounds (default: 12)
      --link-per-x F      extra paired support required per unit of median
                          coverage before a contested join is taken
                          (default: 0.10)
      --bubble-coverage F a bubble branch below this fraction of its twin's
                          coverage is dropped (default: 0.35)

READ TRIMMING
      --no-qtrim          do not quality-trim read 3' ends
      --qtrim-quality N   3' trim threshold (default: 20)
      --qtrim-window N    3' trim window (default: 4)

ORGANISM MODEL
      --organism NAME     organism the reads come from (e.g. klebsiella)
      --model FILE        genus model built by tessera-model, consulted only at
                          junctions no fragment can span
      --is-panel FILE     FASTA of known insertion sequences; contig ends lying
                          inside one are left unjoined
      --is-sites FILE     table of recurrent insertion sites, which lets an
                          element be placed rather than only avoided

A --model or --is-panel that cannot be read is a fatal error, not a warning:
the assembly it would have produced differs from the one without it in exactly
the junctions the model exists to settle.

POLISHING
      --map-polish NAME   polish contigs against a full read alignment:
                          bowtie2 | bwa | none (default none)
      --mapper-dir DIR    where to find the mapper binaries

GENERAL
  -t, --threads N         worker threads (default: all cores)
      --max-memory G      memory ceiling in GB (default: 80% of physical RAM)
      --no-gfa            do not write assembly_graph.gfa
      --no-html           do not write report.html
      --unitigs           also write unitigs.fasta
  -q, --quiet             suppress progress output
  -v, --version           print version
  -h, --help              print this message
```

Every k given to `-k` must be odd (even k admits palindromic k-mers, which have
no well-defined canonical form) and between 5 and 127. Invalid values are
rejected with a message on stderr and a non-zero exit.

## Output files

| Path | Contents |
| --- | --- |
| `OUTDIR/contigs.fasta` | the assembly, longest contig first, wrapped at 80 columns |
| `OUTDIR/assembly_graph.gfa` | the simplified unitig graph in GFA1, with a `P` line per contig |
| `OUTDIR/report.html` | self-contained run report (disable with `--no-html`) |
| `OUTDIR/report.json` | the same data, machine-readable; always written |
| `OUTDIR/unitigs.fasta` | pre-resolution unitigs, only with `--unitigs` |

The GFA opens directly in [Bandage](https://rrwick.github.io/Bandage/). Segments
are unitigs carrying a `dp:f:` depth tag, links carry the `(k-1)M` overlap the de
Bruijn graph implies, and each contig appears as a `P` path recording exactly
which oriented unitigs it walks through — so a contig that stops at a repeat is
visible as a path ending at a branch. Bandage's own `info` on a real assembly
reports every edge overlap as exactly `k-1` and an overlap-free total matching
`contigs.fasta`, which is a useful independent check that the file is sound.

Headers follow the SPAdes convention:

```
>NODE_1_length_5186_cov_26.3083
```

`NODE_<rank>` counts from 1 in output order, `length` is the sequence length in
bases (including any `N` gap bases), and `cov` is the mean k-mer coverage of the
unitigs the contig was built from. If scaffolding joined two contigs, the gap
appears as a run of `N`s inside the sequence. A run summary — contig count,
total length, largest, N50, GC, mean coverage, elapsed time and peak RSS — is
printed to stderr at the end.

## Run modes

`--mode` picks a preset; any flag you set explicitly still wins over it.

| Mode | k ladder (150 bp reads; longer reads top out at 127) | Simplify rounds | Joins | Polish | Use when |
| --- | --- | --- | --- | --- | --- |
| `fast` | 21,55,77 | 6 | strict | off | triage, or a first look at many isolates |
| `standard` | 21,33,55,77,95 | 12 | balanced | 1 pass | the default; what the benchmarks use |
| `careful` | 21,33,45,55,67,77,87,95 | 24 | strictest | 2 passes | when a wrong join costs more than a break |
| `aggressive` | same as `standard` | 16 | loosest, collapses diverged repeats | 1 pass | maximum contiguity, accepting misassembly risk |

Measured on the simulated *S. aureus* set (2.82 Mb, 100x, 5 threads), the modes
trade off as intended:

| Mode | Contigs | N50 | Wall clock |
| --- | --- | --- | --- |
| `fast` | 57 | 199,297 | 42 s |
| `standard` | 43 | 294,925 | 78 s |
| `careful` | 44 | 248,041 | 102 s |
| `aggressive` | 36 | 304,753 | 67 s |

`careful` deliberately lands below `standard` on N50: it refuses joins that
`standard` accepts. `aggressive` buys its extra contiguity by collapsing
diverged repeat copies, which on the benchmark panel costs about one misassembly
per genome — see [Behaviour worth knowing](#behaviour-worth-knowing).

## How it works

### The multi-k ladder

A single k is always a compromise: small k keeps the graph connected through
low-coverage and error-ridden regions but collapses every repeat shorter than k,
while large k separates repeats but shatters wherever coverage dips. tessera
runs the whole assembly at each k in an increasing ladder — by default a ladder
chosen from the observed read length (for 150 bp reads: 21, 33, 55, 77, 95; k
values not shorter than the reads are skipped) — and the contigs from iteration
*i* are fed into the counter at iteration *i+1* as extra sequence, weighted as
several observations each.

That carry-forward is the point of the ladder. A region assembled at k=21 on
thin coverage yields a contig whose 95-mers would otherwise be too rare to pass
the abundance cutoff; injecting the contig puts them comfortably above it, so
the large-k graph keeps the connectivity the small-k graph found while gaining
the repeat resolution that only large k gives.

### Choosing the abundance cutoff

Sequencing errors create k-mers that occur once or twice; genuine k-mers occur
at roughly the sequencing depth. The histogram of "how many distinct k-mers have
count *c*" therefore has a huge spike at low *c*, a valley, and a peak at the
coverage mode. tessera finds the peak by maximising `c * histogram[c]` — which
weights by bases rather than by distinct k-mers, so the error spike cannot win —
then takes the minimum of the histogram between 1 and the peak as the cutoff,
clamped to at least 2 and at most a quarter of the peak so a noisy histogram can
never throw away real sequence. `-c` overrides the whole procedure.

### Graph simplification

The graph is bidirected: each unitig is a double-stranded sequence with two
ends, and each link records which end of the neighbour it attaches to. Four
operations run in a loop until nothing changes, with thresholds that ramp up
over rounds so the obviously spurious goes first and genuinely low-coverage
sequence gets a chance to be joined into something defensible:

* **Tip clipping** — a short dead-end branch is removed when the junction it
  hangs off offers a better-supported alternative.
* **Bubble popping** — two near-identical parallel paths between the same
  branch points are collapsed onto the better-supported one.
* **Erroneous-connection removal** — a short, low-coverage unitig that connects
  two otherwise well-supported regions is cut.
* **Isolated-node removal** — short, low-coverage, entirely disconnected
  fragments are dropped.

After every operation the graph is re-compacted (maximal non-branching chains
merged into single unitigs). `UnitigGraph::validate()` checks the bidirected
link invariant and is asserted clean by the unit tests after each operation.

### Why read anchoring uses k=31

The repeat resolver has to know where each read sits in the graph. The obvious
choice is to look up the graph's own k-mers, but at k=95 that fails badly: a
single sequencing error invalidates every k-mer overlapping it, so one error in
a 150 bp read destroys 95 of its 56 possible 95-mers — nearly all of them. The
reads that matter most for resolution are exactly the ones spanning a junction,
and losing them defeats the stage.

So anchoring uses its own index at k=31 (or the graph's k, whichever is
smaller). A 31-mer index is specific enough to place a read uniquely in a
bacterial genome, and short enough that a read with one or two errors still has
plenty of clean 31-mers. Each read is probed at up to 12 positions spread across
its length, every probe votes for a (unitig, position, orientation) triple, and
a read is anchored only if at least two probes agree — which also rejects reads
that hit a repeat, since a 31-mer occurring in more than one place is marked
ambiguous at index time and ignored.

### Mutual-best chain extension

Anchored pairs are accumulated into a table of oriented unitig → oriented
unitig links, each carrying the fragment span it implies *excluding* whatever
sequence lies between the two unitigs. The fragment-length distribution itself
is learned from pairs whose two reads land on the same unitig, with the top and
bottom 1% trimmed before fitting so chimeric pairs cannot inflate the standard
deviation.

Resolution then starts with every uniquely-anchorable unitig as a chain of one
and repeatedly joins chains. For each chain end, every path out of it that
terminates on an anchorable unitig within fragment reach is enumerated —
intermediate nodes being unanchorable collapsed repeats — and each candidate is
scored by the paired evidence supporting it. Crucially, a pair only counts for a
candidate if the fragment length it implies *once that candidate's intermediate
sequence is inserted* is plausible under the learned model. This is what
distinguishes the correct way through a repeat from a coincidental link: the
same read pair supports one exit and not another purely by arithmetic on the
fragment length. Scoring accumulates over every member of the chain still within
fragment reach of the boundary, not just the last unitig, so support builds up
over a whole contig tail.

Two chains are only joined when both ends independently pick each other
(mutual best). Chains grow round by round until no join is found.

### Why the tie-break is deliberately conservative

A join is refused when the best candidate does not beat the runner-up by
`--tie-ratio` (default 1.15), or has fewer than `--min-link` supporting pairs
(default 2). That threshold is intentionally cautious. A near tie is not a
missing feature — it is the data saying the repeat is genuinely unresolved at
this fragment length. Forcing a choice would join the wrong flanks, and unlike a
break, a misassembly is not detectable downstream and silently corrupts every
analysis built on the assembly. tessera reports a fragmented but correct
assembly instead: an unresolved repeat costs contiguity — it shows up as extra
contigs — rather than correctness.

## Assembling with a genus model

Paired reads settle a junction only when a fragment can span it. On a typical
*K. pneumoniae* library the fragments run about 217 bp while the genome carries
hundreds of repeat copies between 500 bp and 5 kb — rRNA operons, IS elements,
duplicated loci. At those junctions there is no paired evidence to weigh, so an
assembler either stops there or guesses from coverage.

tessera takes a third route: ask closed genomes of the same organism. Core gene
order is strongly conserved within a species, so a panel of finished genomes can
say which flank normally follows which, and at what distance.

### Building a model

```sh
tessera-model --organism klebsiella --out kleb.tsm \
              --plasmids plasmid_db.fasta \
              references/*.fasta
```

The builder samples canonical 31-mers at roughly 500 bp spacing, keeps those
occurring exactly once in a replicon set — a marker seen twice names a repeat
family, which is exactly what cannot settle a junction — and records how often
each ordered pair of markers is observed across the panel, and at what median
distance.

Chromosome and plasmid are learned into separate tables. Core chromosomal gene
order is close to a rule; plasmids are mosaic, recombine, and vary in copy
number, so pooling the two would let plasmid rearrangements vote on chromosomal
junctions. Inputs are classified by file name (`*_chr.fasta` versus
`*_plasmid*`), and `--plasmids` takes a multi-record plasmid database in which
each record is its own replicon.

`--exclude ACC` omits an accession from the panel and records the omission
inside the model file. Every assembly that uses the model repeats the list in
`report.json` under `organism_model.excluded_accessions`, so a leave-one-out
evaluation can be checked from its own output rather than taken on trust:

```console
$ jq '.organism_model.excluded_accessions' out/report.json
[
  "ERR10447223"
]
```

That block also records how many markers were placed, how many joins were made
on the chromosome and on plasmids, and how many candidate joins each rejection
rule threw out.

### Assembling with it

```sh
tessera --organism klebsiella --model kleb.tsm --is-panel is_elements.fna \
        -1 R1.fq.gz -2 R2.fq.gz -o out
```

The model stage runs after paired-end resolution, on exactly the junctions the
reads could not reach, and in two passes: the chromosome first, then plasmids.
A join is taken only when several independent marker pairs agree, a substantial
share of the panel genomes carrying both markers support it, and the implied gap
is consistent across them. The gap is sized from the panel's median distance,
and the existing gap-closing stage then fills it with real read-derived
sequence. Everything failing those tests is left broken.

### Choosing a strain-relevant panel

A whole-species panel describes what the organism does on average. Where the
isolate carries a genuine rearrangement — an IS insertion, an inversion — the
average is confidently wrong about it. Selecting the panel per isolate fixes
that: sketch the closed-genome corpus once, sketch the isolate's own draft
assembly, and learn only from its nearest neighbours, which in practice are its
own sequence type and close relatives.

```sh
mash sketch -s 5000 -o corpus -l chromosome_list.txt
mash dist corpus.msh draft.fasta | sort -k3,3g | head -400 | cut -f1 > near.txt
tessera-model --organism klebsiella --out strain.tsm $(cat near.txt)
```

The number of neighbours is a dial rather than a constant: fewer makes the
evidence more relevant to this strain, more makes it better supported.

### Tuning for contiguity without a model

Two knobs move the conservative default toward contiguity, and they are
independent of each other:

```sh
tessera -1 R1.fq.gz -2 R2.fq.gz -o out \
        --link-per-x 0.04 --tie-ratio 1.05
```

`--link-per-x` sets how much paired support a contested join must have per unit
of median coverage, and `--tie-ratio` how far the winning branch must beat the
runner-up. Lowering both takes more joins. Measured on five closed-reference
isolates, this moves median contig NGA50 from 261,117 to 273,595 and genome
fraction from 98.94% to 98.97%, for one additional misassembly across the set.

Separately, a contig may absorb up to 3 kb beyond its end when the graph offers
only one way forward into a multi-entrance repeat. Every continuation begins
with that same sequence, so appending it decides nothing; it recovers the repeat
copy at that locus, which would otherwise be represented once in the assembly
and counted as missing at every other copy. Measured across 37 closed-reference
isolates this raises genome fraction from 98.83% to 99.02% with the duplication
ratio unchanged at 1.000, so it is **on by default**; `TESSERA_COMMON_PREFIX=0`
disables it and any other value sets the budget in bases.

The join-bar flags are off by default because the default trades contiguity for
certainty.
Whether that trade is right depends on what the assembly is for: a contiguity
number, or a sequence you intend to trust base by base.

### The insertion-sequence panel

`--is-panel` takes a FASTA of known insertion sequences. A contig ending inside
a mobile element ends there *because* this isolate carries an element the panel
need not have, so the panel's adjacency across that point describes a genome
without the insertion. Those ends are left unjoined. A panel can be produced
from any IS annotation — extracting the predicted element coordinates from a
set of closed genomes is enough.

## Benchmark

Per-isolate figures behind everything below:
[`docs/final_headtohead.tsv`](docs/final_headtohead.tsv). Two further pages
cover the model-free assembler in more detail — what each default is worth on
the *Klebsiella* panel ([`docs/KLEBSIELLA_PANEL.md`](docs/KLEBSIELLA_PANEL.md))
and how it behaves on mixed chemistry
([`docs/CLOSED_REFERENCE_BENCHMARK.md`](docs/CLOSED_REFERENCE_BENCHMARK.md)).
Both measure tessera **without** a model, so their contiguity figures are the
0.89x column above, not the model one.

### Real isolates with closed reference genomes

The benchmark is *Klebsiella pneumoniae* clinical isolates for which both paired
Illumina reads and a complete, closed reference genome are available. Every
figure is **contig level** — QUAST is run with `-s` and the `_broken` column is
reported, which splits assemblies at runs of `N`, so scaffolded output is never
compared against another tool's contigs. SPAdes 4.2.0 was run on the same
trimmed reads.

The model column uses a panel chosen per isolate: the closed-genome corpus is
mash-sketched once, the isolate's own draft assembly retrieves its nearest
neighbours, and the model is learned from those. That corpus is a separate
collection from the isolates being tested, so no assembly is scored against a
model that has seen its own genome.

Measured on 40 isolates:

| Metric | tessera | tessera + model | SPAdes |
| --- | --- | --- | --- |
| Median contig NGA50 | 268,292 | 304,672 | **319,598** |
| Median NGA50 ratio to SPAdes | 0.89x | **1.02x** | 1.00x |
| Isolates leading SPAdes on NGA50 | 11 | **25 / 40** | — |
| Median contigs | 78 | 70 | **65** |
| Total misassemblies | **19** | 21 | 26 |
| Assemblies with zero misassemblies | **27** | 26 | 24 |
| Median mismatches / 100 kbp | **0.22** | 0.23 | 0.44 |
| Median genome fraction (%) | 99.22 | **99.23** | 98.95 |

The two contiguity rows disagree on purpose. The model has the longer NGA50 on
25 of 40 isolates and its median *ratio* to SPAdes is 1.02x, but SPAdes holds
the higher *median* NGA50, because the isolates SPAdes wins it wins by a wider
margin. Neither statement is the whole picture, so both are given.

Read across rather than down. Without a model tessera is the more conservative
assembler: fewer misassemblies, half the mismatch rate, higher genome fraction,
and 0.89x the contiguity, because it declines joins it cannot support. The model
buys back most of that contiguity -- 0.89x to 1.02x at the median ratio, and the
longer NGA50 on 25 of 40 isolates -- while keeping the accuracy: still fewer
misassemblies than SPAdes, still about half the mismatch rate, still the higher
genome fraction. It reports more contigs than SPAdes, which is what a higher
genome fraction looks like when the extra sequence is short: the pieces are
counted, and they are also present.

Both columns use the same trimmed reads, the same QUAST invocation and the same
closed references, and the model is always learned from a corpus that excludes
the genome being assembled.

### Plasmids

Genome fraction is dominated by the chromosome, so a plasmid can be missing
entirely and barely move it. Scored per replicon across 78 plasmids in the
panel:

| | Recovered >=99% | >=90% | Below 50% | Median |
| --- | --- | --- | --- | --- |
| tessera | 59 / 78 | **78** | **0** | 99.6% |
| SPAdes | **68 / 78** | **78** | **0** | 99.7% |

Small plasmids are the hard case, and for a reason worth knowing: a 2.5 kb
multicopy plasmid can sit at twenty times the genome's median depth, where its
own sequencing errors clear an abundance cutoff chosen for the chromosome. Those
errors become real unitigs and shatter the replicon -- one came out of the graph
as fifty-seven pieces. tessera judges each unitig against its own neighbours as
well as against the genome, which recovers most of them; `TESSERA_LOCAL_WEAK=0`
turns that off and costs 18 of the 78.

### Reconstructing the chromosome

Contiguity says how long the pieces are, not how much of the chromosome is
present. Aligning every block back to the closed reference answers the second
question, for one isolate:

| | chromosome covered | breaks >=100 bp | in a repeat | in unique sequence |
| --- | --- | --- | --- | --- |
| tessera | 99.72% | 46 | 40 | 6 |
| tessera + model | **99.96%** | **5** | **0** | 5 |
| SPAdes | 99.93% | 9 | 4 | 5 |

With the model the chromosome is essentially complete and no remaining break
lies in a repeat. The five that survive total 1,650 bp in one stretch at 23-30%
GC against a 57% genome, where read depth is 3.2x with 7% of positions
uncovered — SPAdes breaks at the same coordinates. `bench/breakpoints.py`
reproduces this analysis for any isolate.

### Why the model changes the numbers

Fragments in these libraries average about 217 bp. The genomes carry a mean of
355 repeat copies of 500 bp or longer — 8 rRNA operons of ~5.1 kb each, a median
of 39 insertion sequences of 0.8–2 kb, and hundreds of shorter duplications. Not
one of them can be spanned by a read pair, so at those junctions there is no
paired evidence to weigh. The model supplies independent evidence there and
nowhere else.

### Speed

tessera runs roughly 4–5x faster than SPAdes on the same isolate and threads.
Building a model from a few hundred genomes takes seconds; from a corpus of
several thousand, about a minute, and a model is reusable across runs.

## Limitations

* **k is capped at 127** (odd values only). k-mers are packed into four 64-bit words. For 150 bp
  reads the ladder tops out at 95; 2x250 libraries carry enough k-mer per read
  to use 127, which is where the top rung goes automatically.
* **Repeats longer than the fragment length cannot be resolved from the reads
  alone.** No paired-end library can span them, and tessera leaves them unjoined
  rather than guessing. Where a model of the organism is supplied it can settle
  some of those junctions from independent evidence; without one, longer inserts
  are the only fix.
* **Scaffolding is paired-end only and conservative.** Gaps are sized from the
  fragment model and require mutually-best support scaled with depth --
  `min(10, max(3, median coverage x 0.06))` pairs, so a floor of three;
  there is no iterative gap closing, and no scaffolding across anything wider
  than one fragment.
* **No mate-pair or long-read support**, and only one library per run.
* **Haploid bacterial isolates only.** There is no diploid/heterozygosity
  handling, no plasmid copy-number modelling, and no metagenome mode: the
  coverage model assumes a single genome at roughly uniform depth. Eukaryotic
  genomes are out of scope, both for the coverage model and for memory.
* **Memory is roughly proportional to the number of distinct k-mers**, plus the
  reads themselves (2 bits per base, held for the entire run). The error cloud
  before the abundance cutoff dominates peak usage.
* **Single-machine, in-memory only.** No checkpointing, no resume, no
  distribution.

## Behaviour worth knowing

These are deliberate choices rather than accidents, described here so the
output is not surprising:

* **Erroneous-connection removal is length-limited.** Cutting long low-coverage
  connectors risks severing real sequence, so the operation is confined to short
  ones; a long chimeric connector is left in place and splits the components it
  touches.
* **Bubble popping is restrained by default.** A bubble side is discarded only
  when its coverage falls well below the mean, which confines the operation to
  sequencing-error bubbles. Bubbles between genuinely diverged repeat copies sit
  at comparable depth and are kept, because collapsing them fuses distinct loci.
  `--aggressive` lifts the limit and collapses them, trading that risk for
  contiguity.
* **The consensus polisher is a guard, not a corrector.** It requires
  near-unanimity before rewriting a base, and on clean isolate data it is
  expected to change nothing. Inside a repeat the pileup is fed by several
  copies at once, so a simple majority would rewrite one copy into another.
* **Alignment polishing is off by default** for the same reason, and uses a
  0.99 agreement threshold when enabled.
* **One library per run.** Passing `-1`/`-2` more than once replaces the earlier
  value rather than adding a second library.

## Testing

```sh
make check      # unit, end-to-end and flag suites
```

`make unittest` builds `tests/test_units.cpp` against the project objects and
runs 277,712 assertions covering k-mer encoding/decoding, reverse complement,
canonical form and `pushFront` at k = 15, 31, 63, 77, 95 and 127 — spanning all
four 64-bit words of the packed representation — the rolling
forward/reverse-complement window, the open-addressed k-mer table (100,000
random operations checked against `std::unordered_map`, which is what validates
its backward-shift deletion), the banded sequence identity, and the graph's
bidirected link invariant after every simplification operation.

`make flagcheck` runs `tests/check_flags.sh`, which asserts that every flag the
help text advertises parses and runs, that invalid values are rejected, that
`--min-contig` demonstrably changes the output, that runs stay deterministic and
thread-invariant, and that the model and insertion-sequence paths actually
engage rather than being silently ignored.

`make test` runs `tests/run_tests.sh`, which needs only bash and python3. It
generates synthetic genomes and reads itself — no external read simulator — in a
temporary directory it removes on exit, and prints one `PASS`/`FAIL` line per
test. It covers exact reconstruction of an error-free genome, multiple
chromosomes, tolerance of 1% substitution errors, resolution of a repeat that
the fragments span, non-misassembly of a repeat they cannot span, all four input
shapes (paired, interleaved, single, FASTA) and gzipped vs plain, `-k` and
missing-file validation, `--min-contig`, determinism, thread invariance, output
format, and graceful handling of empty and tiny inputs. The whole suite runs in
a few seconds; it also passes cleanly under `make asan`.

The misassembly test deserves a note: it builds a genome `A R B R C` with a
5 kb identical repeat and 350 bp fragments, so the repeat cannot be spanned, and
then asserts that every assembled contig is an exact substring of the source.
Checking only that each contig's k-mers exist in the source is not enough — a
contig that jumps from one copy of the repeat to the wrong flank is built
entirely from real k-mers.

## License

MIT — see [LICENSE](LICENSE).
