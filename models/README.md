# Released models

Models are attached to the [releases](https://github.com/iowa69/TesserACT/releases) rather
than committed: the smaller one is 339 MB and the larger 2.9 GB, well past what belongs in a
git history.

They are **encrypted**. The asset is public — a release asset on a public repository can be
downloaded by anyone with the URL — so encryption is the only thing that restricts the models
themselves. Without the key the download is an opaque blob. The key is issued separately;
open an issue or get in touch.

Be clear about what that does and does not protect. It stops the models being picked up and
reused by anyone who happens across the repository. It does not protect them from anyone who
has been given the key, because this assembler is open source and reads the format: a
key-holder can extract everything in the file. And it protects the *model* — the panel
curation, the leakage exclusions, the density and support thresholds — not the underlying
sequence, which is public RefSeq.

## Unlocking

```sh
export TESSERA_MODEL_KEY=/path/to/keyfile          # or omit and gpg will prompt
./unlock.sh tessera-klebsiella-default-v1.2.0.tsm.zst.gpg kleb.tsm
tessera -1 R1.fq.gz -2 R2.fq.gz -o out/ --organism klebsiella --model kleb.tsm
```

`unlock.sh` verifies `SHA256SUMS` when it is next to the asset, decrypts, decompresses, and
moves the result into place only on success — a failed decryption leaves nothing behind
rather than a partial file the assembler would later reject with a confusing error.

Needs `gpg` and `zstd`. Both are in every distribution's base repositories and in conda.

## What is in the release

| asset | sampling | download | unpacked |
|---|---|---|---|
| `tessera-klebsiella-default` | 1 marker in 512, 16 neighbours | 132 MB | 339 MB |
| `tessera-klebsiella-plasmid` | 1 marker in 128, 64 neighbours | 1.3 GB | 2.9 GB |

**Take the default unless you specifically want plasmid grouping.** Every figure quoted in
[`docs/KLEBSIELLA_REPLICONS.md`](../docs/KLEBSIELLA_REPLICONS.md) was measured with it, and it
is the better model for the chromosome: 98.1% of a 105-isolate held-out set get their
chromosome in one contig at ≥90%, against 95.2% for the plasmid build.

The plasmid build nearly doubles whole-plasmid delivery (23.8% of multi-contig plasmids in one
group against 13.8%) and raises grouping completeness from 0.219 to 0.325, at the cost of that
chromosome tail. The README's *Genus models* section has the full comparison and the reason
the two cannot currently be one model.

## Provenance and leakage

Both are **fold-0 leave-cluster-out** models: built from 2,221 closed *K. pneumoniae*
chromosomes and 51,789 Enterobacterales plasmid records, with the fold-0 mash clusters
withheld — 4,942 panel plasmids among them. That is what makes the benchmark numbers honest,
and it also means these models have never seen roughly a fifth of the public *K. pneumoniae*
diversity. For production use on arbitrary isolates a model built over the whole panel would
be marginally stronger; these are the ones every published figure was measured with, which is
why they are the ones released.

Build your own with `tessera-model` — see the README. It takes about four minutes for the
default density and twenty for the dense one.
