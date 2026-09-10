# E. faecium non-clonal panel, n=113 (complete)

The second of the seven ESKAPEE cohorts. Same protocol as
[SAUREUS_100_PANEL.md](SAUREUS_100_PANEL.md): 124 isolates with both a closed RefSeq
assembly and a public paired-end Illumina library, dereplicated with `mash` at
`d <= 0.0005` so no two are clonal, all 124 assembled in five arms, QUAST `-s` with only
the `_broken` column read, each isolate against its own closed reference.

113 survive the quality gates. Nine are discarded, including two whose genome fraction is
73.1% and 78.5% **in every arm** and one with 25.0% of its assembly unaligned to its own
reference — those are reference/library mismatches, not assembly failures.

## Head to head: `--organism efaecium` vs vanilla SPAdes 4.3.0

| metric | TesserACT | SPAdes | delta | win/loss | p |
|---|---|---|---|---|---|
| NGA50 | 45,128 | 37,972 | **+18.8%** | 100/13 | 4.7e-14 |
| NG50 | 45,835 | 37,972 | **+20.7%** | 101/12 | 2.8e-14 |
| **Genome fraction (%)** | **95.86** | **93.51** | **+2.35 pp** | **110/3** | **4.2e-20** |
| Largest alignment | 151,804 | 145,687 | +4.2% | 83/30 | 9.5e-05 |
| Mismatches / 100 kb | 1.40 | 2.10 | −33.3% | 76/36 | 0.00089 |
| Indels / 100 kb | 0.36 | 0.51 | −29.4% | 68/42 | 0.077 (ns) |
| # local misassemblies | 0 | 1 | −100% | 64/23 | 0.00016 |
| # contigs | 181 | 187 | −3.2% | 69/40 | 0.012 |
| **# misassemblies** | **2** | **0** | — | **14/75** | **2.2e-09** |
| **Duplication ratio** | 1.01 | 1.00 | +0.9% | **0/111** | 6.1e-20 |

Seven wins, two losses. The genome-fraction margin is **five times** the S. aureus one
(+2.35 pp against +0.46 pp) — on a genome this repetitive, recovering the sequence at all
is the harder problem, and that is where the model earns most.

## The misassembly deficit scales with genome difficulty

| panel | median NGA50 | misassemblies, us vs SPAdes | p |
|---|---|---|---|
| *S. aureus*, n=180 | 247,188 | 29/53 | 0.011 |
| *E. faecium*, n=113 | 45,128 | **14/75** | **2.2e-09** |

*E. faecium* assembles to a fifth of the contiguity of *S. aureus* on the same protocol —
it is a hard-graph organism, dense in IS elements. The deficit that reads as a marginal
tail effect on *S. aureus* is decisive here.

This is not a surprise, it is a prediction landing. Data-mining the 180-isolate *S.
aureus* panel found the strongest single predictor of losing to SPAdes was scaffold-gap
count (Spearman −0.402, p≈1e-8), with the NGA50 ratio falling **1.60 → 1.47 → 1.12 →
0.90** across difficulty quartiles: *TesserACT degrades faster than SPAdes as the graph
gets harder.* *E. faecium* is that prediction applied to a whole organism.

## What the model is worth here, and what it costs

Against the same binary with no model, over the same 113 isolates:

| model vs no model | median | win/loss | p |
|---|---|---|---|
| NGA50 | 45,128 vs 41,224 (+9.5%) | **90 / 0** | 1.8e-16 |
| NG50 | 45,835 vs 41,224 (+11.2%) | **92 / 0** | 8.3e-17 |
| Genome fraction | 95.86 vs 95.84 | 81 / 26 | 7e-11 |
| **# misassemblies** | 2 vs 1 | **1 / 40** | **7e-08** |

**The model never loses contiguity on a single isolate out of 113.** It also makes 40 of
them worse on misassemblies and 1 better.

Compare *S. aureus*, where the same contrast was 2 better / 18 worse. The cost is an
order of magnitude larger here, and it scales with the same axis as everything else: the
marker-adjacency table licenses joins the paired reads would not make alone, and on a
repeat-dense genome more of those joins are wrong.

The conclusion is that the models need **re-tuning for hard genomes, not retirement**. The
shape of the fix is known, because the same problem was already solved once for the
replicon classifier: requiring a second, independent signal before acting took its pooled
precision from 0.225 to 0.712. An adjacency-only join wants the same treatment.

## Reproducing

```
devtools/panels/select_cohort100.py efaecium
devtools/panels/fetch_all_v3.sh efaecium
devtools/panels/run_asm100.sh efaecium 5
devtools/panels/collect100.py eval100/efaecium out.tsv out.txt
```

Per-isolate values for all five arms and ten metrics:
[efaecium_100_headtohead.tsv](efaecium_100_headtohead.tsv).
