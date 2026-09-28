# Released models

Models are attached to the [releases](https://github.com/iowa69/TesserACT/releases) rather
than committed. There are two sets, and neither belongs in a git history:

- **`models-v2`**, the seven ESKAPEE organism models, one `<organism>.tsm` each, about 1.4 GB
  together. `tesseract-get-models` downloads them into `~/.tesseract/models` and checks each
  against [`models.sha256`](../models.sha256). The assembler selects one with
  `--organism <organism>` and looks for it there, or in the directory `TESSERACT_MODEL_DIR`
  names. `tesseract-eskape --preset <organism>` and `tesseract-klebsiella --with-model` use the
  same file.
- **The 1.2 *Klebsiella* models**, 339 MB and 2.9 GB, attached to the `v1.2.0` release.
  `tesseract-klebsiella --plasmid` downloads the plasmid one for itself.

**The models are optional and opt-in.** TesserACT is complete without them, and nothing
downloads or loads one unless you ask: `./install.sh` asks before downloading them and defaults
to no.

**A model buys contiguity and costs misassemblies.** With the 1.4.0 defaults, on 329 held-out
isolates, adding the organism model cut the median contig count from 99 to 86 and raised the
median NGA50 from 134 to 156 kb. It also raised misassemblies from 154 to 290 on the 306
isolates whose closed reference matches the reads. The assembler loads a model only when
`--organism` is given. The measurement is in [What a model costs](#what-a-model-costs).

## models-v2: the seven ESKAPEE models

The release is named `models-v2`, and there is no `models-v1`. The `tesseract-get-models` of
1.3.0 fetches from a `models-v1` release and checks against a different checksum list, so
publishing these files as `models-v1` would make every 1.3.0 install download 1.4 GB only to
reject each file.

| asset | organism | bytes | size | chromosomes | plasmid records |
|---|---|---:|---:|---:|---:|
| `abaumannii.tsm` | *Acinetobacter baumannii* | 96,477,839 | 96.5 MB | 831 | 3,487 |
| `ecloacae.tsm` | *Enterobacter cloacae* | 296,827,515 | 296.8 MB | 934 | 41,443 |
| `ecoli.tsm` | *Escherichia coli* | 300,556,888 | 300.6 MB | 1,507 | 41,336 |
| `efaecium.tsm` | *Enterococcus faecium* | 63,509,962 | 63.5 MB | 811 | 5,676 |
| `kpneumoniae.tsm` | *Klebsiella pneumoniae* | 349,390,015 | 349.4 MB | 2,238 | 41,585 |
| `paeruginosa.tsm` | *Pseudomonas aeruginosa* | 214,016,711 | 214.0 MB | 1,358 | 1,146 |
| `saureus.tsm` | *Staphylococcus aureus* | 92,588,552 | 92.6 MB | 1,473 | 4,135 |
| **all seven** | | **1,413,367,482** | **1.41 GB** | | |

These are the sha256 checksums. They are the same as [`models.sha256`](../models.sha256) and
the release's `SHA256SUMS` asset:

```
8d4de675a3c89209d1349796b46c607c09aacea6bc3858396a98ef549f796ec5  abaumannii.tsm
bd175e9ab2877f95383f1329a314933c607ecb71344252e9db6bc26e3727dc30  ecloacae.tsm
d8702ec9ab45aac72f4a33d03475ddff2708901fd6d3ef0c1bedf5399bb5ae9d  ecoli.tsm
825f883c37e2566d0376b1abb321de29ab477b070adb466fa7e08f93411aae1a  efaecium.tsm
c668812cc7ef81f3d34e4be9b9ac38c40f5b56d980bd741ea2a89bd8c43fa8a0  kpneumoniae.tsm
ed47f4035b0a473c4d2091e9707d31bdd9b6bbff941613360e6e7470884773a3  paeruginosa.tsm
86e6ec7f1afbf1a391ce41bfbcb5267f5eef68a90fbc2a0033e4f96637ede216  saureus.tsm
```

`tesseract-get-models` fetches all seven, or only the ones named (`tesseract-get-models saureus`).
It checks each file against that list and skips any it has already verified. To fetch by hand,
download from the [`models-v2` release](https://github.com/iowa69/TesserACT/releases/tag/models-v2)
and run `sha256sum --check --ignore-missing SHA256SUMS`.

Each model is `TSMODEL5`, k = 31, with one layout track per chromosome:

- **Chromosomes.** Every selected genome of the organism's clonal panel. For *E. coli* this is
  the capped panel.
- **Plasmids.** PLSDB 2025 (`2024_05_31_v2`) and the MOB-suite database 3.1.8, curated by
  [`devtools/panels/curate_plasmids.py`](../devtools/panels/curate_plasmids.py) and split by
  genus.
- **Parameters.** `tesseract-model --marker-density 512 --min-support 5 --min-support-plasmid 3
  --layout-tracks`.

### Whole-panel builds, not the builds that were measured

The measurement below was made with **leave-clone-out** builds. For each organism, the build
withheld the whole clonal cluster of every test isolate, together with every plasmid record of
those clusters. So no measured gain can come from a model that had already seen the answer.

The `models-v2` files are **whole-panel builds**. They use the same chromosome panel, the same
plasmid panel and the same parameters, with **nothing withheld**. A released model has no test
isolate to protect, so it holds everything the measured build held, and more. That includes
the closed genomes of the isolates the measurement used, so these files cannot be scored on
those isolates.

**The numbers below describe the leave-clone-out builds, not these files.** Two things were
checked on the released files:

- Rebuilding each leave-clone-out model with the same builder and panels reproduces the
  measured model byte for byte.
- Each released model loads through `--organism` and runs its model stage on a real isolate of
  its organism.

That shows the files work, not how accurate they are.

### What a model costs

This compares the 1.4.0 defaults with and without each organism's leave-clone-out model. There
are 47 held-out isolates per organism, 329 in all, and each has a closed genome and its own
Illumina reads. QUAST 5.3.0 scored every `contigs.fasta` against the closed genome
(`-s --min-contig 500`).

The misassembly counts leave out 23 isolates, whose closed genome is not the strain that was
sequenced: 9 *K. pneumoniae*, 11 *P. aeruginosa*, 2 *E. coli* and 1 *A. baumannii*. TesserACT
and SPAdes both disagree with the reference on each of them, at more than 25 times the median
mismatch rate, and misassemblies counted against a wrong reference are not the assembler's.
The other rows use all 329.

| | no model | with model |
|---|---:|---:|
| misassemblies, the 306 isolates whose reference matches the reads | **154** | 290 |
| of those 306, isolates with more / fewer / the same misassemblies | | 89 / 2 / 215 |
| misassemblies, all 329 isolates | **1,571** | 1,743 |
| contigs, median | 99 | **86** |
| NGA50, median | 134 kb | **156 kb** |
| genome fraction, mean | 96.85% | **96.91%** |
| mismatches per 100 kb, median | **0.53** | 0.57 |

| organism | isolates (reference matches) | misassemblies | median NGA50, kb | median contigs |
|---|---|---:|---:|---:|
| *K. pneumoniae* | 47 (38) | 70 → 93 | 155 → 204 | 94 → 76 |
| *A. baumannii* | 47 (46) | 14 → 44 | 97 → 115 | 116 → 99 |
| *E. cloacae* | 47 (47) | 7 → 24 | 208 → 209 | 77 → 64 |
| *E. coli* | 47 (45) | 23 → 45 | 117 → 133 | 162 → 138 |
| *E. faecium* | 47 (47) | 23 → 47 | 49 → 50 | 161 → 149 |
| *P. aeruginosa* | 47 (36) | 5 → 16 | 181 → 210 | 97 → 86 |
| *S. aureus* | 47 (47) | 12 → 21 | 186 → 214 | 39 → 32 |

Each arrow runs from no model to model. All seven organisms move the same way: fewer, longer
contigs and more misassemblies.

With or without a model, no isolate had 90% of its chromosome in one correct block.

TesserACT puts misassemblies first, so models are opt-in. `tesseract-asm` never loads a model
on its own, `./install.sh` asks before downloading them and defaults to no,
`tesseract-eskape --preset <organism>` uses the model only if it is installed, and
`tesseract-klebsiella` uses one only with `--with-model`, `--model FILE` or `--plasmid`. To
assemble one of these organisms without a model while the models are installed, run
`tesseract-asm` without `--organism`: the presets add no other option. An improved organism
model is in development.

Both arms were run with the pre-release build of the configuration that 1.4.0 made the
default. Without a model, 1.4.0 writes byte-identical assemblies to that build on every isolate
where this was checked.

These figures replace the model gains of earlier releases, "+25% to +70% longer contigs, and
fewer misassemblies". Those were measured with older defaults and older models, in
[`docs/REAL_WORLD_RESULTS.md`](../docs/REAL_WORLD_RESULTS.md) and
[`docs/ESKAPEE_MODELS.md`](../docs/ESKAPEE_MODELS.md), and they do not hold for 1.4.0.

## Using one

For the seven ESKAPEE organisms:

```sh
tesseract-get-models                                            # optional; once, about 1.4 GB
tesseract-asm -1 R1.fq.gz -2 R2.fq.gz -o out/ --organism saureus
tesseract-eskape --preset saureus -1 R1.fq.gz -2 R2.fq.gz -o out/   # the same model, via the preset
```

`klebsiella` is accepted as another name for `kpneumoniae`, and reads the same
`kpneumoniae.tsm`. The model must have been built for the organism named, which is a guard
against pointing a *Klebsiella* model at something else by accident; the 1.2 models were built
as `klebsiella`, which counts as `kpneumoniae`. The assembler does not take a model file on
the command line: `--model FILE` is refused unless `TESSERACT_MODEL_AUTHOR` is set, which is
for building and validating the bundled models.

`tesseract-klebsiella` runs without a model unless asked. `--with-model` uses the same
`kpneumoniae.tsm`, and fetches it with `tesseract-get-models` if it is not installed;
`--model FILE` takes another model file, and `--plasmid` the 1.2 plasmid model below. The
helper hands any of them to the assembler itself and prints what a model costs. To run the
assembler directly with a 1.2 model, put it where `--organism` looks for it, under the name
`kpneumoniae.tsm`:

```sh
mkdir -p kleb-model
ln -s "$PWD/tesseract-klebsiella-default-v1.2.0.tsm" kleb-model/kpneumoniae.tsm
TESSERACT_MODEL_DIR=kleb-model tesseract-asm -1 R1.fq.gz -2 R2.fq.gz -o out/ --organism kpneumoniae
```

## The 1.2 *Klebsiella* models

Kept for compatibility. Since 1.4.0 `tesseract-klebsiella` no longer uses the default model on
its own; `--plasmid` still fetches the plasmid model. Neither was measured with the 1.4.0
defaults: the figures below are from the 1.2 series.

| asset | sampling | download | unpacked |
|---|---|---|---|
| `tesseract-klebsiella-default-v1.2.0.tsm` | 1 marker in 512, 16 neighbours | 339 MB | — |
| `tesseract-klebsiella-plasmid-v1.2.0.tsm.zst` | 1 marker in 128, 64 neighbours | 1.3 GB | 2.9 GB |

The default model is uploaded as-is and works straight away. The plasmid model is
zstd-compressed only because GitHub caps a release asset at 2 GB and the raw file is 2.9 GB:

```sh
zstd -d tesseract-klebsiella-plasmid-v1.2.0.tsm.zst -o kleb-plasmid.tsm
```

Check your download first — a truncated 1.3 GB file will otherwise fail later with a confusing
error about the file not being a model:

```sh
sha256sum --check --ignore-missing SHA256SUMS
```

### Which one

**Take the default unless you specifically want plasmid grouping.** Every figure quoted in
[`docs/KLEBSIELLA_REPLICONS.md`](../docs/KLEBSIELLA_REPLICONS.md) was measured with it.

| | default | plasmid |
|---|---|---|
| multi-contig plasmids delivered whole | 13.8% | **23.8%** |
| grouping completeness, per isolate | 0.219 | **0.325** |
| grouping homogeneity, per isolate | 1.000 | 1.000 |
| chromosome in one contig ≥90% | **98.1%** | 95.2% |
| chromosome in one contig ≥98% | 64.8% | **69.5%** |

The plasmid model nearly doubles whole-plasmid delivery. The chromosome rows above are the
share of isolates whose chromosome lands in one *record*, and that measure rewards a model for
scaffolding aggressively whether or not the resulting order is right — the default model places
12.44 contigs per isolate against the plasmid model's 10.04. So the same 105 held-out isolates
were also scored by alignment, with QUAST against the full reference:

| | default | plasmid |
|---|---|---|
| misassemblies, cohort total | **584** | 674 |
| median NGA50 | 2,959,825 | 2,958,278 |
| median genome fraction | 99.186% | **99.215%** |
| isolates won, misassemblies | **48** | 30 (27 tied) |
| isolates won, NGA50 | 35 | **60** (10 tied) |
| isolates won, genome fraction | 10 | **92** (3 tied) |

By alignment the two are close to parity: the plasmid model reconstructs more of the genome on
92 of 105 isolates and is more contiguous on 60, and pays for it with 90 extra misassemblies
across the cohort. Median NGA50 differs by 1,547 bp on a 5.3 Mb chromosome. That is a genuine
trade, and it is a good deal narrower than the ≥90% row on its own suggests — quote the two
together rather than either alone.

The README's *Genus models* section explains why the two cannot currently be a single model —
the adjacency reach has to scale with the sampling density, and doing that recovers most but
not all of the cost.

### Provenance and leakage

Both are **fold-0 leave-cluster-out** builds: 2,221 closed *K. pneumoniae* chromosomes and
51,789 Enterobacterales plasmid records, with the fold-0 mash clusters withheld and 4,942 panel
plasmids excluded from the plasmid panel. That exclusion is what makes the published numbers
honest — 46% of test plasmids had a near-identical panel match, so without it the model would
have been scored on recall of memorised sequence.

It also means these models have never seen roughly a fifth of public *K. pneumoniae*
diversity. For production use on arbitrary isolates a model built over the whole panel would be
marginally stronger; these are released because they are the ones every published figure was
measured with.

## What is in the file

The format is plain and documented by `src/organism.cpp`. A model stores canonical 31-mers
with per-replicon-class counts, the marker adjacency tables, one marker-order track per panel
chromosome, and per-plasmid marker membership sets. Nothing is obfuscated: the k-mers decode
directly back to sequence, which is public RefSeq in any case.

Build your own with `tesseract-model` (`make model`) — see the README. About four minutes at
the default density, twenty at the dense one. Name the file `<organism>.tsm` and point
`TESSERACT_MODEL_DIR` at its directory to use it.
