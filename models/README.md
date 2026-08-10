# Released models

Models are attached to the [releases](https://github.com/iowa69/TesserACT/releases) rather
than committed: 339 MB and 2.9 GB is not something that belongs in a git history.

| asset | sampling | download | unpacked |
|---|---|---|---|
| `tessera-klebsiella-default-v1.2.0.tsm` | 1 marker in 512, 16 neighbours | 339 MB | — |
| `tessera-klebsiella-plasmid-v1.2.0.tsm.zst` | 1 marker in 128, 64 neighbours | 1.3 GB | 2.9 GB |

The default model is uploaded as-is and works straight away. The plasmid model is
zstd-compressed only because GitHub caps a release asset at 2 GB and the raw file is 2.9 GB:

```sh
zstd -d tessera-klebsiella-plasmid-v1.2.0.tsm.zst -o kleb-plasmid.tsm
```

Check your download first — a truncated 1.3 GB file will otherwise fail later with a confusing
error about the file not being a model:

```sh
sha256sum --check --ignore-missing SHA256SUMS
```

## Using one

```sh
tessera -1 R1.fq.gz -2 R2.fq.gz -o out/ --organism klebsiella \
        --model tessera-klebsiella-default-v1.2.0.tsm
```

`--organism` must match the name the model was built with, which is a guard against pointing a
*Klebsiella* model at something else by accident.

## Which one

**Take the default unless you specifically want plasmid grouping.** Every figure quoted in
[`docs/KLEBSIELLA_REPLICONS.md`](../docs/KLEBSIELLA_REPLICONS.md) was measured with it, and it
is the better chromosome model.

| | default | plasmid |
|---|---|---|
| multi-contig plasmids delivered whole | 13.8% | **23.8%** |
| grouping completeness, per isolate | 0.219 | **0.325** |
| grouping homogeneity, per isolate | 1.000 | 1.000 |
| chromosome in one contig ≥90% | **98.1%** | 95.2% |
| chromosome in one contig ≥98% | 64.8% | **69.5%** |

The plasmid model nearly doubles whole-plasmid delivery and costs a tail of the chromosome
result: it breaks four isolates in a held-out hundred while rescuing one badly broken one. The
README's *Genus models* section explains why the two cannot currently be a single model — the
adjacency reach has to scale with the sampling density, and doing that recovers most but not
all of the cost.

## Provenance and leakage

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

Build your own with `tessera-model` — see the README. About four minutes at the default
density, twenty at the dense one.
