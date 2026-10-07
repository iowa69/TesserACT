# TesserACT: release history and benchmarks (through 1.5.0)

This page keeps the README of 1.5.0 unchanged, for its release notes, benchmark tables and design notes. The current manual is [README.md](../README.md).


A de novo short-read genome assembler: de Bruijn graphs over a multi-k ladder, with
paired-end repeat resolution and an optional genus model for the junctions no read pair
can span.

Built for bacterial isolates. Reads FASTQ or FASTA, gzipped or plain, paired, interleaved
or single-end. No dependencies beyond zlib and a C++17 compiler.

---

## What 1.5.0 changes

**The default assembly does not change.** 1.5.0 writes the same `contigs.fasta`,
`scaffolds.fasta`, `scaffolds.agp` and `assembly_graph.gfa`, byte for byte, as 1.4.0, and
`--organism` runs the same 1.4.0 model path. Measured byte for byte against 1.4.0 on 9 isolates
(one of each ESKAPEE species, *M. tuberculosis* and *S. enterica*), and with `--organism` on 4.

* **A Salmonella model.** `senterica.tsm` (1,530 closed *S. enterica* chromosomes and the
  Enterobacterales plasmid database) joins the seven ESKAPEE models in the `models-v2` release.
  Select it with `--organism senterica` (`salmonella` is accepted) or
  `tesseract-eskape --preset salmonella`. Like every model it is optional and has a cost: on 15
  development isolates (held-out not yet measured) it raised the median NGA50 from 150 to 267 kb
  and the misassemblies from 3 to 10. See [`models/README.md`](models/README.md).
* The Organism Model 2.0 code (junction evidence, gap filling from the isolate's own graph, a
  genome view with a layout-only mode, the nearest-relative layout) is included, every part of it
  off by default and not reachable from the command line. It is still in development: on the
  development isolates its layout-only mode came within 6 misassemblies of the default
  (123 against 117 on 175 ESKAPEE isolates), which missed its own pre-set bar of at most 5 % more,
  so it is not offered in this release.

## What 1.4.0 changed

**Misassemblies first.** TesserACT assemblies are the raw material for genus models, and a
misassembly built into a model is repeated in every assembly the model later completes. So
1.4.0 puts correctness ahead of contiguity. It fixes 42 defects found in a systematic audit
of 1.3.0, and it makes the configuration with the fewest misassemblies of every
configuration measured (called F2 in the campaign records) the default. That costs some
contiguity, and the cost is stated below with the gain.

The measurements come from the combo3 campaign: 4 panels of real Illumina libraries, each
isolate scored by QUAST against its own closed reference, `--min-contig 500`. The full record
is kept with the campaign files (`combo3/FINAL.md`, `combo3/DEFECTS.md`), not in this
repository.

### New defaults

| Setting | 1.3.0 | 1.4.0 | What it does |
|---|---|---|---|
| `--tie-ratio` | 1.02 | **3.0** | The winning branch needs 3× the runner-up's pair support, which removes paired-join errors |
| `TESSERACT_COMMON_PREFIX` | 3000 | **250** | A short evidence-free walk at chain ends, in bp |
| `TESSERACT_PREFIX_MIN_BODY` | off | **reach** | That walk starts only from chain ends whose body reaches the insert size |
| `TESSERACT_MIN_FALLBACK_DEST` | 0 | **1000000000** | Withdraws the three fallback joins |
| `TESSERACT_REQUIRE_SUPPORT_SINGLE` | off | **on** | A lone join candidate still needs pair support |
| `TESSERACT_EXCLUDE_SHARED_REPEAT_SUPPORT` | off | **on** | A repeat that supports every destination loses its votes |
| `TESSERACT_FIXES` | off | **on** | The 21 output-changing defect fixes (below) |
| `TESSERACT_DROPOUT_BRIDGE` | off | **on** | Bridges a coverage dropout between two dead ends when at least 2 read pairs span it |

Every run prints the values in force on one line of its log:

```
[defaults] tie_ratio=3.000 common_prefix=250 prefix_min_body=reach min_fallback_dest=1000000000 require_support_single=1 exclude_shared_repeat_support=1 fixes=1 dropout_bridge=1
```

The run modes keep their 1.3.0 tie ratios (`fast` 1.3, `careful` 1.4, `aggressive` 1.05). So
`--mode careful` still adds simplification rounds and polishing passes, but it is no longer
stricter than `standard` on the tie ratio. It will be re-measured for 1.4.1.

### Reproducing 1.3.0

The 1.3.0 behaviour is one command line away. Measured: the four output files
(`contigs.fasta`, `scaffolds.fasta`, `assembly_graph.gfa`, `scaffolds.agp`) are byte-identical
to release 1.3.0 on 3 real isolates and on the test fixtures (`tests/golden_130.md5`).

```sh
TESSERACT_FIXES=0 TESSERACT_COMMON_PREFIX=3000 TESSERACT_PREFIX_MIN_BODY=0 \
TESSERACT_MIN_FALLBACK_DEST=0 TESSERACT_REQUIRE_SUPPORT_SINGLE=0 \
TESSERACT_EXCLUDE_SHARED_REPEAT_SUPPORT=0 TESSERACT_DROPOUT_BRIDGE=0 \
  tesseract-asm --tie-ratio 1.02 -1 R1.fq.gz -2 R2.fq.gz -o out/
```

A single fix is turned off with `TESSERACT_FIX_<NAME>=0`, for example
`TESSERACT_FIX_REVISIT_GUARD=0` (T07) or `TESSERACT_FIX_CARRY_READ_GATE=0` (T01). An explicit
per-fix value always wins over the umbrella. One combination is refused: the gap-flank restore
(T03) with the polisher's N-skip (T02) off exits with status 2, because the polisher could
then overwrite gaps that T03 has restored.

### Results

On two fresh panels that no development decision had seen, against SPAdes 4.3.0:

* **Fresh ESKAPEE isolates (140):** fewer misassemblies (636 against 681), but not
  significantly at the pre-registered level (one-sided p = 0.0495 against α = 0.025). Because
  that test failed, everything else on this panel is descriptive only: genome fraction is
  higher and assembly size error is lower.
* **Fresh 28-species diversity panel (82 isolates):** fewer misassemblies (133 against 187;
  −0.64 per isolate, one-sided p = 0.0005 in the pre-registered analysis, which drops the
  one isolate SPAdes could not assemble). This result is **descriptive, not confirmatory**:
  that analysis needed two protocol deviations, and they were not ratified.
* **Contiguity is not better than SPAdes** on either panel. The NGA50 and contig-count
  comparisons are inconclusive.

These results do not show that 1.4.0 is better than SPAdes on ESKAPEE isolates. They also
do not show a confirmed improvement over K2, the best configuration of the 1.3 series
(`--tie-ratio 2.0 TESSERACT_COMMON_PREFIX=0` with the fallbacks withdrawn and both support
rules on; it was never a shipped default). On the fresh ESKAPEE panel K2 and F2 tie on
misassemblies: 637 against 636.

Misassembly totals on all four panels (lower is better). The development panels (dev_c and
esk140b) are the ones F2 was selected and checked on:

| | dev_c (108) | esk140b (140) | esk140c (140, fresh) | divfresh (82, fresh) |
|---|---|---|---|---|
| **1.4.0 defaults (F2)** | **350** | **628** | **636** | **133** |
| K2 | 367 | 643 | 637 | 155 |
| SPAdes 4.3.0 | 411 | 679 | 681 | 187 |

Of the 22 fewer divfresh misassemblies against K2, 15 come from one low-coverage isolate (see
the cost below). On fresh data, the gains over K2 that did replicate are genome fraction
(+0.13 points on esk140c, +0.07 on divfresh), size error, and local misassemblies (314
against 337, and 107 against 119).

The 1.3.0 defaults had significantly more misassemblies than SPAdes on the campaign's two
earlier development panels (ESK-140 and DIV-146; `combo/FINAL.md` §3.2), as the 1.3.0 table
below shows for the 30-species panel. F2 itself was not scored head-to-head against 1.3.0.

### What it costs

Measured against K2 on the two fresh panels:

| | esk140c (140 ESKAPEE) | divfresh (82, 28 species) |
|---|---|---|
| Contigs | +3.3 per isolate (mean) | +9.7 mean, +3 median per isolate |
| NGA50, geometric mean ratio | ×0.940 | ×0.931 |
| Duplication ratio | +0.0002 | +0.0002 |

The campaign's legacy scorer prints the NGA50 cells as "WIN" and "tie". Both are
significant losses, as the rank sums and the ratios show; the scorer takes the direction
from win/loss counts that are nearly equal here.

Two of the defaults account for most of the cost, and both are kept on deliberately:

* **The revisit guard (T07).** Alone it costs NGA50 ×0.96, 2.5 contigs per isolate and
  0.01 points of genome fraction on the development panel, without removing an extensive
  misassembly. It is kept on because it removes local misassemblies (232 → 223) and
  mismatches (26 isolates better, 0 worse).
  `TESSERACT_FIX_REVISIT_GUARD=0` turns it off.
* **The carry-read gate (T01) on low-coverage libraries.** T01 accounts for 11 of the 12
  fewer misassemblies that the fixes give on the development panel. On one low-coverage
  *S. maltophilia* isolate (k-mer depth 25 at k = 99) it takes misassemblies from 18 to 3,
  but costs 2.6 points of genome fraction (90.58 % → 87.99 %) and about 260 more contigs.
  Whether the 15 removed events were real chimeras has not been checked on raw reads.
  `TESSERACT_FIX_CARRY_READ_GATE=0` turns it off. A low-coverage guard is planned for 1.4.1.

Run time is about ×1.07 that of K2 and peak memory about ×1.01 (development panel).

### The 42 defect fixes

45 defects were triaged. 42 are fixed, T20 only in part; 5 of the 42 (T04, T05, T19, T44,
T45) are in the evaluation harness, outside this repository. 3 were refuted as filed (T06,
T12, T23). The ids are those of the campaign's defect register.

**Output-changing, on by default through `TESSERACT_FIXES`** (each also has its own
`TESSERACT_FIX_<NAME>` switch):

| Id | Defect | Switch |
|---|---|---|
| T01 | Carried k-mers made read-free junctions solid ("inside one unitig" chimeras) | `CARRY_READ_GATE` |
| T02 | The polisher overwrote scaffold N-runs | `POLISH_SKIP_N` |
| T03 | k−1 real bases were dropped after every open N-gap | `GAP_FLANK` |
| T07 | A tandem loop was collapsed by the lone-candidate path | `REVISIT_GUARD` |
| T08 | The legacy gap-close nominator continued from the wrong end | `GAP_NOMINATOR` |
| T09 | A gap-fill search cut off by its budget was accepted as unique | `GAPFILL_STRICT_BUDGET` |
| T10 | The terminal-overlap trim could cut both copies of an equal-length pair | `BOUNDARY_SAFE_TRIM` |
| T11 | The gap filler anchored on the exact terminal k-mers, with no back-off | `GAPFILL_BACKOFF` |
| T13 | A scaffold cycle lost every join in it | `SCAFFOLD_CYCLE` |
| T20 | The gap-length estimate was one-sided (partly fixed: its short bias remains) | `GAP_ESTIMATE` |
| T21 | A rival discarded by a truncated enumeration counted as absent | `TRUNC_GUARD` |
| T22 | The terminal-overlap trim ignored copy number | `TRIM_COPY_GUARD` |
| T24 | Contig coverage counted the k−1 overlap of every unitig | `COV_CONTRIB` |
| T25 | Split pieces skipped the per-record filters and deduplication | `SPLIT_POSTPROCESS` |
| T26 | AGP and GFA path records were built from records before post-processing | `AGP_GFA_V2` |
| T27 | The gap filler read input N bases as A | `GAPFILL_SKIP_INPUT_N` |
| T28 | A hairpin self-link was lost in a merge | `HAIRPIN_KEEP` |
| T29 | Graph simplification stopped without counting every removal | `SIMPLIFY_COUNT_ALL` |
| T30 | Read correction left the rest of a capped read unmasked | `EC_CAP_MASK` |
| T31 | Route-distance allocation summed in thread order | `ROUTE_ORDER` |
| R6 | Mutual joins did not check that both ends chose the same route | `MIRROR_ROUTE` |

**Unconditional, and output-neutral** (logs, reports and refusal of invalid input):

* T14, T15, T35, T36: every run prints the resolver, gap-close, repeat-threshold and
  dead-end counters, zeros included.
* T16: every `TESSERACT_*` variable is validated at startup, and a malformed value exits 2
  with the variable named. 1.3.0 read `1e9` as 1 and `true` as 0 without a word.
* T17: `report.json` and the summary describe the records actually written. `report.json`
  gains a `gap_fill` object.
* T18: a truncated gzip input is refused.
* T32: conflicting read options and bad `-k` lists are refused.
* T33: FASTQ blank lines and `.1`/`.2` mate suffixes are accepted.
* T34: re-running into an output directory leaves no stale files.
* T37, T38, T39: faster deduplication, no repeated gap-filler walk, and thread buffers that
  no longer grow with the square of `-t`.

**Tests and packaging:** T40 to T43. Tests no longer read the caller's `TESSERACT_*`
environment; the python tests run in `make check`; the conda recipe is checked by
`make recipecheck`.

**New in 1.4.0 beyond the register:** the T02/T03 coupling is enforced at startup (N9), and
`make componenttest` builds the binary its end-to-end tests need (N14).

**Known open items.** Gap fills from a search that hit its solution cap are still accepted
(N8; see Limitations). The coverage labels of split pieces use a different convention from
the resolver's (N10; names only). The T01 low-coverage guard (N20) is 1.4.1 work.

---

## Benchmark

### 30 species, 146 isolates (1.3.0)

Measured with the 1.3.0 defaults, before the 1.4.0 changes above.

A generalisation panel outside the ESKAPE pathogens: 30 bacterial species from 1.6 Mb
(*H. pylori*) to 7.7 Mb (*B. cenocepacia*) and 27–67 % GC, up to five isolates each, every
isolate paired by BioSample with its own closed reference genome. Raw untrimmed Illumina
reads (100–301 bp), both assemblers at default settings — SPAdes 4.3.0 with its own error
correction and automatic k, TesserACT with no model. QUAST 5.2.0, `--min-contig 500`.
Paired Wilcoxon signed-rank on per-isolate differences, Holm-corrected across the six rows;
win / tie / loss counts are per isolate. (NGA50 is undefined for one isolate, hence n=145.)

| Metric | TesserACT | SPAdes | win / tie / loss | |
|---|---|---|---|---|
| NGA50 | **149,610** | 139,749 | 95 / 0 / 50 | **win** (p=0.005) |
| Genome fraction | **98.46 %** | 97.93 % | 131 / 0 / 15 | **win** (p=1e-17) |
| Contigs | **68** | 84 | 109 / 0 / 37 | **win** (p=3e-7) |
| Assembly size error | **56,475 bp** | 71,170 bp | 114 / 0 / 32 | **win** (p=2e-10) |
| Duplication ratio | 1.0000 | 1.0000 | 40 / 79 / 27 | tie (p=0.72) |
| Misassemblies | 1 | 0 | 29 / 56 / 61 | loss (p=0.005) |

Misassemblies are the one row still lost: 328 events against SPAdes' 182 across the panel.
Nearly all of the excess is *relocations* (305 against 161; inversions 17 against 13,
translocations 6 against 8) — a contig that joins two correctly assembled blocks across a
repeat it has collapsed by one copy. At most of those junctions the repeat is as long as the
sequencing fragment, so the reads are equally consistent with the right and the wrong join;
SPAdes avoids them by stopping there, which is part of why it emits more contigs. Per-base
accuracy is higher than SPAdes' on the same panel (17,970 against 20,491 mismatches in total).
It is being worked on, not hidden.

What 1.3.0 changed to get here, each measured on this panel:

* **Head-to-head overlaps.** The terminal-overlap trimmer only ever probed each contig's
  forward 3′ end, so a 5′/5′ overlap — a quarter of all exact terminal overlaps — was never
  seen. Duplication went from a loss (p=0.001) to a tie.
* **The repeat threshold's baseline.** The single-copy depth the resolver compares against
  was an unweighted median over graph nodes. On some libraries a cloud of short, shallow
  fragments outvotes the chromosome — on one *S. enterica* isolate 33 % of the nodes held
  2 % of the sequence and set the baseline at 5.7× for a genome sequenced at 35×, so almost
  everything counted as repeat and the resolver joined nearly nothing. The baseline is now a
  length-weighted median. NGA50 went from a tie to a win.
* **Single-read threading.** Where one read spans a whole branch point, its exact path through
  the graph now counts as evidence for the join. Against the release without it: NGA50
  better on 27 isolates and worse on 1, fewer contigs on 65 against 4, genome fraction better
  on 39 against 10. It is not free — misassemblies rose on 4 isolates and fell on 1, five
  events net across 146 — and that is the row this release still loses.

Each of the three can be switched off for comparison with the environment variables
`TESSERACT_RC_DOVETAIL=0`, `TESSERACT_WEIGHTED_RESOLVER_COVERAGE=0` and
`TESSERACT_EXACT_READ_THREADS=0`.

### 666 *Klebsiella pneumoniae* (1.2.x)

Measured on the 1.2 series and not yet re-run on 1.3.0.

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

### With a genus model: 185 non-clonal *S. aureus*

The panel above is vanilla against vanilla. This one is `--organism saureus`, on a cohort
dereplicated with `mash` at `d <= 0.0005` so that no two isolates are clonal redeposits —
which is what makes a paired test meaningful. Every arm of all 185 isolates was built by a
single binary, so the comparison is clean by construction rather than by filtering; the
cohort is described in [docs/SAUREUS_100_PANEL.md](docs/SAUREUS_100_PANEL.md).

| Metric | TesserACT | SPAdes | W / L | |
|---|---|---|---|---|
| NGA50 | **250,176** | 184,157 | 143 / 41 | **win** (+35.8%, p=1e-13) |
| NG50 | **270,913** | 193,052 | 146 / 39 | **win** (+40.3%, p=2e-14) |
| Genome fraction | **98.86 %** | 98.33 % | 174 / 11 | **win** (p=2e-28) |
| Largest alignment | **546,622** | 423,362 | 141 / 41 | **win** (+29.1%, p=2e-12) |
| LGA50 | **4** | 5 | 107 / 30 | **win** (p=5e-06) |
| Mismatches /100 kb | **0.67** | 1.32 | 113 / 66 | **win** (−49.2%, p=4e-04) |
| Contigs | **36** | 40 | 102 / 77 | **win** (p=8e-03) |
| Indels /100 kb | 0.25 | 0.31 | 91 / 79 | tie (p=0.06) |
| Duplication ratio | 1.003 | 1.000 | 4 / 152 | loss (p=3e-25) |
| Misassemblies | 0 | 0 | 25 / 62 | loss (p=2e-04) |

Seven wins, one tie, two losses — more of the genome, more of the bases right, in longer and
fewer contigs, paid for in redundancy and misassemblies.

**Both losses have been taken apart, and one of them is solved.**

*Duplication* had two causes in almost equal parts. We emitted contigs that are exact
substrings of longer contigs — 6.2% of ours against 0 of SPAdes' 1,554 — because a resolved
repeat was written both inside its flanking chain and again standalone. And our contig-end
overlaps ran far longer than SPAdes': 70% of SPAdes' duplicated bases sit in a single overlap
of exactly its final K, while 90.6% of ours sat in overlaps longer than 2k, at IS-element
lengths. Dropping the contained copy and trimming the redundant dovetail takes duplication
against SPAdes from 0 wins / 35 losses (p=3e-07) to **9 / 6, p=0.887** — a dead heat — while
every other win above is unchanged to the decimal. It costs ~0.08 pp of genome fraction,
because some of that redundancy was legitimately covering both copies of a two-copy repeat.

*Misassemblies* split by the length of the shorter alignment block flanking each breakpoint
into two disjoint classes with a dead zone between them: a `>= 20 kb` class that is the price
of the contiguity, and a `< 2 kb` class that is not — it is identical in the no-model arm and
is significant even on the isolates where SPAdes is the more contiguous assembler. Stage
ablation puts our raw graph at 100 such events against SPAdes' 88, with our repeat resolver
adding 54 more. Five candidate mechanisms have been tested and all five fail, including every
threshold in the resolver. The analysis, the negatives and the stop condition are in
[docs/EXPERIMENT_LEDGER.md](docs/EXPERIMENT_LEDGER.md).

An earlier version of this section reported a misassembly residual surviving contiguity
matching at p=0.022. That measurement was made on an older binary, **does not reproduce**
(p=0.074), and the control itself was wrong: a sub-2 kb wrong tail costs no contiguity, so a
defect of that shape hides inside the matched band instead of showing up outside it.

**What the model is worth.** Against the same binary with no model, over the same 185
isolates: NGA50 **104 wins / 4 losses** (+12.2%, p=1e-15), NG50 109/4, contigs **164/2**
(p=1e-28), LGA50 84/3, genome fraction 99/36. It costs misassemblies (5/26, p=3e-04) and a
little accuracy (mismatches 14/51). Without a model TesserACT still beats SPAdes on NGA50
(124/60, +21.0%), genome fraction (174/11) and mismatches (120/59), and ties on contigs.

**Replicons.** 117 of the 185 isolates carry plasmids, 146 in total. Every arm and SPAdes
recover ~93% of them at >=50% coverage. SPAdes assembles 39.7% of them into a single contig
against our 28.8%; that gap is a copy-number effect with an identified cause and a measured
fix ([docs/PLASMID_COPY_NUMBER.md](docs/PLASMID_COPY_NUMBER.md)). Contig classification has
no SPAdes arm, since SPAdes does not classify: pooled over every contig >=1500 bp, the model
arm reaches precision 0.700 / recall 0.803 (F1 0.748) against 0.406 / 0.549 with no model.

### With a genus model on 1.4.0: all seven ESKAPEE organisms

The same comparison has now been made for all seven organisms, with the 1.4.0 defaults and a
leave-clone-out model for each, on 47 held-out isolates per organism (329 in all). The model
lowers the median contig count from 99 to 86 and raises the median NGA50 from 134 to 156 kb.
It also raises misassemblies in every organism: from 154 to 290 on the 306 isolates whose
closed reference matches the reads. Of those isolates, 89 get more misassemblies, 2 get fewer
and 215 are unchanged. With or without a model, no isolate had 90% of its chromosome in one
correct block. The downloadable `models-v2` files are whole-panel builds of the same models
with nothing withheld. The per-organism table, and what was and was not measured on the
downloads, are in [models/README.md](models/README.md).

So the models are **optional and opt-in**. TesserACT is complete without them:
`tesseract-asm` loads one only when `--organism` is given, `./install.sh` asks before
downloading them and defaults to no, a `tesseract-eskape` preset uses one only if it is
installed, and `tesseract-klebsiella` uses one only with `--with-model`. An improved organism
model is in development.

---

## Install

```sh
git clone https://github.com/iowa69/TesserACT.git
cd TesserACT
./install.sh
```

Run in a terminal, `./install.sh` walks through it: it checks the four tools it needs and
names the one command that installs any that are missing, asks where to put things, builds,
puts the commands on your PATH, and offers the optional organism models (the default answer
is no). Nothing goes system-wide and no password is needed.

It is also the non-interactive installer -- run from a script, a CI job or a pipe there is
no terminal to ask, so it takes the defaults and says nothing. `--no-prompt` forces that in
a terminal too, `--guided` forces the questions anywhere, and `--prefix DIR` or
`--conda-env NAME` choose the destination outright.

You end up with four commands: `tesseract-asm`, `tesseract-klebsiella`, `tesseract-eskape` and
`tesseract-get-models`. The model builder, `tesseract-model`, is built by `make model` and is not
installed.

To build without installing:

```sh
make -j                     # needs zlib headers; on conda, CPATH=$CONDA_PREFIX/include
make test                   # end-to-end checks on synthetic genomes
```

Produces `./tesseract-asm`.

## Quick start

For *Klebsiella*, one command does everything -- it builds the assembler if needed and
assembles every read pair it is given. It uses no model unless you add `--with-model`, which
fetches and checks the *K. pneumoniae* model once and prints what it costs
([above](#with-a-genus-model-on-140-all-seven-eskapee-organisms)):

```sh
./tesseract-klebsiella reads/                  # a directory of pairs
./tesseract-klebsiella sample_R1.fq.gz         # the mate is found automatically
./tesseract-klebsiella --with-model reads/     # with the organism model (optional)
```

Interrupt it and run it again; it picks up where it stopped. Everything below is the general
interface.

```sh
## Paired-end isolate, all cores
tesseract-asm -1 reads_R1.fq.gz -2 reads_R2.fq.gz -o out/

## Interleaved, 8 threads
tesseract-asm --12 reads.fq.gz -o out/ -t 8

## Maximum contiguity, accepting more misassembly risk
tesseract-asm -1 R1.fq.gz -2 R2.fq.gz -o out/ --mode aggressive

## More simplification passes and 2 polishing passes. Its tie ratio (1.4) is below the
## 1.4.0 default (3.0), so it is no longer the stricter mode on joins
tesseract-asm -1 R1.fq.gz -2 R2.fq.gz -o out/ --mode careful

## With an organism model, for junctions no fragment spans. Reads kpneumoniae.tsm from
## ~/.tesseract/models (tesseract-get-models puts the seven ESKAPEE models and senterica
## there), or from the directory TESSERACT_MODEL_DIR names; `klebsiella` is accepted for
## kpneumoniae and `salmonella` for senterica
tesseract-asm -1 R1.fq.gz -2 R2.fq.gz -o out/ --organism kpneumoniae


## Hand it a QC report from scepter (see below)
tesseract-asm -1 R1.fq.gz -2 R2.fq.gz -o out/ --qc sample.json
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
tesseract-asm -1 R1.fq.gz -2 R2.fq.gz -o out/ --qc sample.json
```

The report carries read depth, an estimated genome size, a substitution rate measured by
comparing overlapping mates against each other, and the insert-size distribution — all
available *before* the first k-mer is counted. Measured against the 666 closed references,
the genome-size estimate lands within 10 % on 620 of them and the depth estimate tracks
truth at Spearman 0.95.

What that currently buys you: the run reports what the library actually is, and the
assembler will refuse a QC file that does not match its reads. The decisions those numbers
could drive are implemented but **off by default**, because none of them has yet beaten
what the assembler already infers from its own graph. They are exposed as `TESSERACT_QC_*`
environment variables for anyone who wants to experiment. Accurate measurement in search
of a use is worth shipping as exactly that, and not as a feature.

---

## Genus models

`--organism NAME` selects a model by organism: the assembler reads `NAME.tsm` from
`~/.tesseract/models`, or from the directory `TESSERACT_MODEL_DIR` names, and the model must
have been built for that organism. The names are `saureus`, `efaecium`, `abaumannii`,
`paeruginosa`, `ecloacae`, `ecoli`, `kpneumoniae` and, since 1.5.0, `senterica`; `klebsiella` is
accepted as another name for `kpneumoniae`, and `salmonella` for `senterica`.
`tesseract-get-models` downloads and checks the eight models from the `models-v2` release, and `tesseract-eskape` presets and `tesseract-klebsiella --with-model`
select theirs the same way. The models are optional: each buys contiguity and costs
misassemblies (see "With a genus model on 1.4.0" above). There is no option that takes a model
file: `--model FILE` is kept for building and validating the bundled models and is refused
unless `TESSERACT_MODEL_AUTHOR` is set.

The 1.2 *Klebsiella* models are attached to the
[releases](https://github.com/iowa69/TesserACT/releases) — see
[`models/README.md`](models/README.md) for which one to take. To build your own, `make model`
builds the model builder, and `TESSERACT_MODEL_DIR` points the assembler at the result:

```sh
tesseract-model --organism kpneumoniae --out mymodels/kpneumoniae.tsm --layout-tracks \
              --plasmids plasmid_panel.fna --exclude-plasmids withheld.txt \
              --marker-density 64 \
              closed/*.fasta
TESSERACT_MODEL_DIR=mymodels tesseract-asm -1 R1.fq.gz -2 R2.fq.gz -o out/ --organism kpneumoniae
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

Most of that chromosome cost is a parameter rather than a law, and the parameter is
`kMarkerNeighbours` — 16, chosen so that at ~500 bp spacing the adjacency table reaches ~8 kb
and covers the rRNA operons and IS elements that defeat paired evidence. At 1 in 64 the
spacing is ~62 bp and those sixteen neighbours reach ~1 kb. Holding the *reach* constant while
raising density — 1 in 128 with 64 neighbours, again ~8.2 kb — recovers four of the five
isolates the 1-in-64 build broke, exactly (99.2 → 86.9 → 99.2, and three like it), and lands
ahead of the default at the tighter thresholds:

| | 1/512, 16 nb | 1/64, 16 nb | 1/128, 64 nb |
|---|---|---|---|
| plasmids delivered whole | 13.8% | 21.3% | **23.8%** |
| chromosome ≥ 90% | **98.1%** | 94.2% | 95.2% |
| chromosome ≥ 98% | 64.8% | 59.6% | **69.5%** |
| median chr_best | 98.6 | 98.5 | **98.7** |

`kMarkerNeighbours` is build-side only — the query looks edges up rather than regenerating
them — so a model built with a different count is read by the stock binary. It is a
compile-time constant (`TS_MARKER_NEIGHBOURS`) rather than a flag, because the value that is
right for a given density is not something a user should have to guess at.

The ≥90% line still favours the default, and that is not smoothed over: the denser build
breaks a different set of four isolates while rescuing one badly broken one, so the tail moved
rather than shrank. A second mechanism is visible in the logs and is not reach — layout
placements that collide and get discarded track density, not neighbour count.

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
| `--organism NAME` | optional organism model `NAME.tsm`, see above: more contiguity, more misassemblies (`TESSERACT_MODEL_DIR` says where) |
| `--qc FILE` | scepter QC report |
| `--min-contig N` | shortest contig to report (default 2k) |
| `--max-memory GB` | counting-table budget (default 80 % of RAM) |

`tesseract-asm --help` lists the rest, including the repeat-resolution and polishing knobs.

---

## Limitations

* **Short reads only.** No long-read or hybrid support.
* **Repeats longer than the top k stay unresolved** unless paired reads span them or a
  model covers them. In *Klebsiella* that means rRNA operons (~5 kb) and IS elements
  (1–2.5 kb). Raising k further does not help — 99.7 % of the genome is already unique at
  k=99, and the rest is far longer than any k a 250 bp read can support.
* **Bacterial isolates.** Metagenomes and eukaryotes are untested.
* **Low-coverage libraries lose some genome fraction to the carry-read gate (T01; N20,
  open).** On one *S. maltophilia* isolate at k-mer depth 25 it cost 2.6 points while removing
  15 of 18 misassemblies. `TESSERACT_FIX_CARRY_READ_GATE=0` trades them back. A low-coverage
  guard is planned for 1.4.1.
* **Gap fills from a search that hit its solution cap (defect N8, open).** When the gap
  filler's search stops at its cap of 24 candidate fills, a fill can still be accepted by
  dominance, although the search did not see every alternative. This is the sibling of the
  expansion-budget case that `TESSERACT_FIX_GAPFILL_STRICT_BUDGET` (T09) refuses. It is rare
  but can emit unverified sequence. Measured on the development panels: 5 such fills in 248
  isolates with the configuration that 1.4.0 makes the default (6 isolates were affected
  across the configurations measured). Every run reports the count as `cappedAccepted=` on
  the `[gapfill-budget]` line of its log. A fix will first get its own default-off flag.

## Licence

See `LICENSE`.
