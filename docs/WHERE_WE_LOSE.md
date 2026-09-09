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

Define **paired reach** = fitted insert size − read length. It is how far past its own read a
pair can vouch for. When the fragment is shorter than the read, the mates overlap completely
and reach goes to zero or below — the pair says nothing the read did not already say.

TesserACT resolves repeats from paired evidence. SPAdes' `SimpleExtensionChooser` does too,
but refuses at an ambiguous branch instead of continuing. So the two should converge exactly
where paired evidence vanishes, and they do:

| library | n | our NGA50 | SPAdes | ratio | win/loss | p |
|---|---|---|---|---|---|---|
| **mates overlap (reach <= 0)** | **23** | **162,059** | **158,959** | **1.02x** | **15/8** | **0.354 (ns)** |
| reach 0–60 bp | 28 | 228,496 | 176,652 | 1.29x | 20/8 | 0.023 |
| reach 60–120 bp | 53 | 272,584 | 231,262 | 1.18x | 39/14 | 0.001 |
| reach > 120 bp | 76 | 288,882 | 189,584 | **1.52x** | 59/17 | 7.7e-07 |

Spearman across all 180 isolates, reach against our NGA50 ratio: **+0.193, p=0.0085.**
The trend is continuous, not an artefact of the binning.

**On 23 of 180 isolates — 13% — our entire margin over SPAdes disappears.** Not a loss; a
parity we should not be at, given we are 1.5x ahead when the library has reach. Losing
isolates and winning isolates have the same coverage (94.6x vs 96.8x median) and the same
read length. Reach separates them; depth does not.

### These libraries are also physically damaged, and we do not repair them

In that band the median fitted insert is **225 bp against a 301 bp read**. Every read runs
about 76 bp past the end of its own fragment, into adapter. That shows up in our own logs:

| band | n | median insert | median read length | 3' quality-trimmed | masked unvouchable |
|---|---|---|---|---|---|
| mates overlap | 23 | 225 | 301 | 0.90% | **3.15%** |
| reach 0–60 | 28 | 266 | 251 | 0.20% | 1.74% |
| reach 60–120 | 53 | 333 | 251 | 0.00% | **0.23%** |
| reach > 120 | 76 | 368 | 151 | 0.10% | 1.10% |

Three to fourteen times more sequence masked as unvouchable than the clean bands. The worst
cases:

| isolate | insert | read length | overhang | masked | NGA50 vs SPAdes |
|---|---|---|---|---|---|
| GCF010364725v2 | 154 | 301 | 147 bp | 4.68% | **0.40x** |
| GCF046268025v1 | 175 | 301 | 126 bp | 1.85% | **0.34x** |
| GCF038024725v1 | 227 | 351 | 124 bp | **22.31%** | 1.06x |
| GCF045347525v1 | 214 | 351 | 137 bp | 4.98% | 0.73x |

The corrector is doing its job — it correctly finds that those trailing bases have no k-mer
support and masks them. But masking a 130 bp tail leaves a fragment, not a read, and the
information that the tail was *adapter* is never used.

**TesserACT does no adapter detection at all.** There is no occurrence of `AGATCGGAAGAGC`
or any adapter logic anywhere in `src/`. What exists is only 3' quality trimming, which
removes 0.9% where 25–40% of the read is adapter.

And the detection is free: `LibraryQC` already carries `insertPeak` and `meanReadLength`,
and the QC stage **already overlaps the mates** to build the insert histogram and to measure
the substitution rate. Everything needed to notice `insertPeak < meanReadLength` is computed
and then not acted on.

## What to fix, in order of expected value

1. **Adapter read-through / mate merging.** Detect `insertPeak < meanReadLength` and either
   hard-trim each read to the fragment or merge the overlapping mates into one consensus
   read. Merging is strictly better: it removes the adapter *and* corrects the overlap by
   consensus, which is the same comparison `libqc` already performs. Target: 23 of 180
   isolates (13%) currently at parity with SPAdes where the rest of the cohort is 1.2–1.5x
   ahead.
2. **Plasmid copy number in the repeat test.** Separately measured and implemented; whole-
   plasmid recovery 33% -> 44% on a 39-isolate subset. See
   [PLASMID_COPY_NUMBER.md](PLASMID_COPY_NUMBER.md).
3. **The residual misassembly gap at matched contiguity** (p=0.022). No mechanism identified
   yet; `pickByCoverage` and gap filling were both tested and neither explains it.

## What was tested and did *not* explain a loss

Recorded so the same ground is not covered twice:

* **`pickByCoverage`**, the resolver's depth-continuity fallback — ablated on 23 isolates,
  every metric identical to the base pair. It fires up to 38 times per isolate and changes
  the contigs, but arrives at the same assembly.
* **Gap filling** (`--no-gapfill`) — removes 9 of 65 misassemblies on the worst 10 isolates
  and costs NGA50 6/0 (p=0.036) and contigs 10/0 (p=0.006). A real contributor at 14%, and
  a bad trade.
* **Library depth** — does not separate winning from losing isolates (94.6x vs 96.8x).
* **Preset switching by library type** — `aggressive` looked better than `--organism` in the
  overlapping-mate band by medians (178,905 vs 162,059) but the paired test is 7/8, p=0.589.
  There is no support for an automatic preset rule, and the median difference was a
  small-sample artefact.
