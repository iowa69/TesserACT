# ESKAPEE models — six new organism models, measured on held-out genomes

**All six improve contig NG50. Two are free; four cost misassemblies. *S. aureus* is the best of
them by a wide margin, on both axes.**

| species | panel | genome | base contig NG50 | model contig NG50 | gain | misassemblies | scaffold inflation |
|---|---|---|---|---|---|---|---|
| ***S. aureus*** | 137 | 2.94 Mb | 210,470 | **391,724** | **+86.1%** | **2 → 2** | 8.9× |
| *A. baumannii* | 131 | 4.12 Mb | 181,909 | 281,219 | +54.6% | 1 → 4 | 17.0× |
| *E. faecium* | 125 | 2.88 Mb | 49,722 | 67,484 | +35.7% | 0 → 3 | 43.9× |
| *P. aeruginosa* | 157 | 7.11 Mb | 269,010 | 341,819 | +27.1% | 0 → 5 | 20.2× |
| *E. coli* | 124 | 5.63 Mb | 101,287 | 125,434 | +23.8% | 1 → 6 | 43.5× |
| *E. cloacae* | 126 | 4.89 Mb | 236,092 | 275,188 | +16.6% | **0 → 0** | 19.1× |

Four held-out genomes per species, each **removed from the panel before training**, reads
simulated at 60× (wgsim, seed 11). NG50 against the true genome length.

## Panels

Complete chromosomes only, mash-dereplicated at d ≥ 0.0005 so that outbreak clones do not
dominate — public collections are heavily skewed toward them, and 100 near-identical ST239
*S. aureus* genomes would teach one arrangement rather than a species. *S. aureus* was given the
widest net (600 listed, 300 downloaded, 137 kept) as the priority organism.

## Scaffold inflation is why both numbers are reported

The rightmost column is scaffold NG50 divided by contig NG50 — the factor by which reporting the
record as written, N included, would overstate the result. It runs from 8.9× to **43.9×**.

On *E. faecium* the model emits a 2.8 Mb "chromosome" whose contig NG50 is 67 kb: the record is
overwhelmingly N. Reported as scaffold NG50, *E. faecium* would look like the **best** model here.
It is one of the weaker ones. This is the same error that inflated the Klebsiella headline 7.9×,
and reporting both columns is the fix.

## The cost is species-specific and must be measured, not predicted

An earlier reading of the first four species suggested that larger, more repeat-dense genomes
gain least and pay most. **The full six do not support it.** *E. cloacae* is 4.89 Mb and costs
nothing; *E. faecium* is 2.88 Mb and costs three. Genome size tracks the gain loosely and does
not predict the misassembly cost at all. There is no shortcut here: each organism has to be
measured on held-out data, which is what this design is for.

## Recommendation per organism

- ***S. aureus*: ship, enabled by the preset.** +86.1% for zero misassemblies and the lowest
  inflation of the six. This is a good model.
- ***E. cloacae*: ship.** Modest +16.6%, but free.
- ***A. baumannii*, *E. faecium*, *P. aeruginosa*, *E. coli*: ship, but as an opt-in.** Each buys
  real contiguity and pays real misassemblies (+3 to +5 across four genomes). That is a trade a
  user should make deliberately, and the preset documentation says so rather than hiding it.

## Caveats

Simulated reads — uniform coverage, no adapters, no GC bias, flat 0.5% error, no indels — are
easier than real ones, so every gain here is an upper bound. Four held-out genomes per species is
a small sample; the per-genome spread is wide (*S. aureus* ranged +68.6% to +103.7%). The
misassembly proxy is a collinearity check on minimap2 primary alignments, not QUAST.

**A build defect found and fixed along the way.** `accessionOf()` in `model_main.cpp` cuts a
filename at the first `_` or `.`, so panel files named `GCF_001518735.1.fasta` all resolved to
the accession `GCF` and collapsed into a single genome. The build reported success — 1 genome,
0 adjacencies, 1 layout track, 69 kB — with no error. Three of the six panels were affected
(125, 100 and 33 files) and would have produced silently useless models.
