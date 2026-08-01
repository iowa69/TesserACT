# tessera

A de novo short-read genome assembler written from scratch in C++17. tessera
builds a de Bruijn graph over a ladder of k-mer sizes, simplifies it, resolves
repeats with paired-end information, and polishes the result — the same overall
shape as SPAdes, in under 4,000 lines with no dependencies beyond a C++17
compiler, zlib and pthreads.

It is aimed at haploid bacterial isolates sequenced with paired-end Illumina
reads.

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

k-mers are packed into three 64-bit words (192 bits), so **k can go up to 96**
and must be odd.

## Requirements

* a C++17 compiler (developed and tested with g++; `make CXX=clang++` should
  work but is not part of the tested configuration)
* zlib (`libz-dev` / `zlib-devel`)
* pthreads

Nothing else — no CMake, no Boost, no external k-mer counter.

## Build

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

## Quick start

```sh
tessera -1 reads_1.fq.gz -2 reads_2.fq.gz -o asm -t 16
```

Assembled contigs land in `asm/contigs.fasta`. Other input shapes:

```sh
tessera --12 interleaved.fq.gz -o asm         # one interleaved file
tessera -s single.fq.gz -o asm                # unpaired reads
tessera -1 R1.fa -2 R2.fa -o asm              # FASTA is accepted too
```

Input may be FASTQ or FASTA, plain or gzipped; the format is detected from the
first record and gzip from the file contents, not the extension.

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
                          (default: chosen from the read length; max 96)
  -c, --cutoff N          k-mer abundance cutoff (default: auto-detect)
      --min-link N        paired reads needed to trust a join (default: 2)
      --tie-ratio F       winning branch must beat the runner-up by F (default: 1.15)
      --aggressive        also collapse diverged repeat copies: more
                          contiguity, but risks a misassembly per genome
      --no-correct        skip read error correction
      --no-resolve        skip paired-end repeat resolution
      --no-scaffold       skip scaffolding
      --no-polish         skip consensus polishing

GENERAL
  -t, --threads N         worker threads (default: all cores)
  -q, --quiet             suppress progress output
  -v, --version           print version
  -h, --help              print this message
```

Every k given to `-k` must be odd (even k admits palindromic k-mers, which have
no well-defined canonical form) and between 5 and 96. Invalid values are
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

| Mode | k ladder (150 bp reads) | Simplify rounds | Joins | Polish | Use when |
| --- | --- | --- | --- | --- | --- |
| `fast` | 21,55,77 | 6 | strict | off | triage, or a first look at many isolates |
| `standard` | 21,33,55,77,95 | 12 | balanced | 1 pass | the default; what the benchmarks use |
| `careful` | 21,33,45,55,67,77,87,95 | 24 | strictest | 2 passes | when a wrong join costs more than a break |
| `aggressive` | 21,33,55,77,95 | 16 | loosest, collapses diverged repeats | 1 pass | maximum contiguity, accepting misassembly risk |

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
per genome — see [Known issues](#known-issues).

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

## Benchmark

The reference benchmark is *E. coli* K-12 MG1655 with simulated 2×150 bp
Illumina reads at 100× and an insert size of 350 ± 50, assembled on 16 threads
and evaluated with QUAST 5.x against the reference, alongside SPAdes 4.2.0 run
with `--isolate`. Reported metrics are NGA50, contigs ≥ 500 bp, largest contig,
genome fraction, misassemblies and local misassemblies, mismatches and indels
per 100 kbp, duplication ratio, runtime and peak RAM.

Results are being re-measured against the current code and will be inserted
here:

### Setup

Four bacterial references spanning 33-57% GC and 2.8-5.7 Mb, including a
multi-replicon isolate. Reads simulated with ART (`HS25`, 2x150 bp, 100x,
fragment 350 +/- 50 bp). Both assemblers ran with 16 threads on the same
machine; assemblies were scored with QUAST 5.x against the reference. SPAdes
4.2.0 ran in `--isolate` mode, the setting its documentation recommends for
high-coverage bacterial isolates. `SPAdes-scaf` is SPAdes' scaffold output,
shown for completeness — the contig columns are the like-for-like comparison.

### Contiguity — NGA50 (bp, higher is better)

| Genome | Size | GC | tessera | SPAdes | SPAdes-scaf |
| --- | --- | --- | --- | --- | --- |
| *E. coli* K-12 MG1655 | 4.64 Mb | 50.8% | **133,088** | 133,048 | 176,676 |
| *S. aureus* NCTC 8325 | 2.82 Mb | 32.9% | **294,925** | 247,992 | 448,842 |
| *K. pneumoniae* HS11286 | 5.68 Mb | 57.1% | **214,410** | 181,239 | 192,569 |
| *L. monocytogenes* EGD-e | 2.94 Mb | 38.0% | **496,800** | 496,766 | 496,766 |

tessera's contigs beat SPAdes' contigs on all four genomes, and beat SPAdes'
*scaffolds* on two of the four.

### Correctness

| Genome | Misassemblies | Mismatches /100 kbp | Indels /100 kbp | Duplication |
| --- | --- | --- | --- | --- |
| *E. coli* | 0 vs 0 | **0.04** vs 0.31 | 0.20 vs 0.18 | 1.000 vs 1.000 |
| *S. aureus* | 0 vs 0 | **0.00** vs 0.32 | 0.22 vs 0.18 | 1.000 vs 1.000 |
| *K. pneumoniae* | 0 vs 0 | **0.00** vs 0.09 | **0.05** vs 0.11 | 1.000 vs 1.000 |
| *L. monocytogenes* | 0 vs 0 | **0.00** vs 0.03 | 0.00 vs 0.00 | 1.000 vs 1.000 |

Values are tessera vs SPAdes contigs. Neither assembler misassembled anything on
this panel. tessera's substitution rate is at or below SPAdes' everywhere and
zero on three of the four genomes.

### Completeness

| Genome | Genome fraction % | Contigs >=500 bp | Largest contig |
| --- | --- | --- | --- |
| *E. coli* | 98.313 vs 98.344 | **80** vs 83 | **327,141** vs 327,108 |
| *S. aureus* | **98.913** vs 98.884 | 26 vs **24** | **995,613** vs 881,901 |
| *K. pneumoniae* | **98.399** vs 98.360 | **72** vs 77 | **502,601** vs 326,874 |
| *L. monocytogenes* | 98.959 vs **99.003** | **13** vs 14 | 887,084 vs **887,423** |

### Resources

| Genome | Wall clock | Peak RAM |
| --- | --- | --- |
| *E. coli* | **68.7 s** vs 185.8 s | **2.59 GB** vs 6.04 GB |
| *S. aureus* | **37.8 s** vs 109.5 s | **1.45 GB** vs 3.68 GB |
| *K. pneumoniae* | **79.3 s** vs 248.7 s | **2.83 GB** vs 7.38 GB |
| *L. monocytogenes* | **40.5 s** vs 109.8 s | **1.48 GB** vs 3.87 GB |

tessera is 2.7-3.1x faster and uses 2.3-2.6x less memory.

### Reproducing

```sh
bash bench/simulate.sh        # fetch references, simulate reads
bash bench/run_benchmark.sh   # assemble with both tools, score with QUAST
```

### Caveat

These are simulated reads with uniform coverage and no adapter or contaminant
content. Real data has coverage bias, chimeras and library artefacts that the
two assemblers handle differently, so read this as a controlled comparison of
assembly algorithms rather than a prediction of real-world performance.


Note that `make test` does not reproduce these figures: it uses small synthetic
genomes and checks correctness, not assembly performance.

## Limitations

* **k is capped at 96.** k-mers are packed into three 64-bit words. SPAdes goes
  to 127. In practice k=95 is past the point of diminishing returns for 150 bp
  reads, but longer reads would benefit from more.
* **Repeats longer than the fragment length cannot be resolved.** No paired-end
  library can span them, and tessera deliberately leaves them unjoined rather
  than guessing (see the tie-break discussion above). Longer inserts are the
  only fix.
* **Scaffolding is paired-end only and conservative.** Gaps are sized from the
  fragment model and require mutually-best support from at least five pairs;
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

## Known issues

These are real, observed behaviours of the current code, kept here rather than
in the prose above so they are not mistaken for design:

* **Erroneous-connection removal is length-limited to `2k`.** A low-coverage
  chimeric connector longer than that survives the whole simplification schedule
  and splits both components it touches. The threshold is deliberately tight —
  cutting long connectors risks severing real sequence — but it is tighter than
  the operation's description suggests.
* **Bubble popping is deliberately restrained.** By default a bubble side is
  only discarded when its coverage falls below 35% of the mean, which confines
  the operation to sequencing-error bubbles. On the benchmark panel the bubbles
  that survive the abundance cutoff are almost all genuinely diverged repeat
  copies at comparable depth, so the default pops nothing there. `--aggressive`
  lifts the limit and collapses them too. Measured on the panel, that raises
  NGA50 substantially on two genomes (*E. coli* 133,088 to 174,375; *S. aureus*
  294,925 to 304,753) but lowers it on a third (*K. pneumoniae* 214,410 to
  203,855) and introduces one misassembly each on *E. coli* and *S. aureus*,
  where the default has none anywhere. Contiguity is not worth a fused locus by
  default, so it is opt-in.
* **Only one library per run.** Passing `-1`/`-2` more than once overwrites the
  earlier value rather than adding a second library.

## Testing

```sh
make check      # both suites
```

`make unittest` builds `tests/test_units.cpp` against the project objects and
runs ~140,000 assertions covering k-mer encoding/decoding, reverse complement,
canonical form and `pushFront` at k = 15, 31, 63, 77 and 95 — spanning all three
64-bit words of the packed representation — the rolling
forward/reverse-complement window, the open-addressed k-mer table (100,000
random operations checked against `std::unordered_map`, which is what validates
its backward-shift deletion), the banded sequence identity, and the graph's
bidirected link invariant after every simplification operation.

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
