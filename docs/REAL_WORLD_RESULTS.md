# Real-world ESKAPEE evaluation — recovered results, 279 isolates

**Stopped mid-run and consolidated. 279 of 360 isolates scored, every one a real clinical
isolate with a closed genome as truth and its own Illumina reads as input, against models
trained under leave-cluster-out.**

| species | n | median cov | base contig NG50 | model contig NG50 | gain | median | misassemblies | improved | genome fraction | scaffold inflation |
|---|---|---|---|---|---|---|---|---|---|---|
| *A. baumannii* | 42 | 30× | 138,287 | 234,593 | **+69.6%** | +62.1% | 70 → 48 | **42/42** | 97.3 → 99.6% | 22× |
| *P. aeruginosa* | 49 | 25× | 176,333 | 295,855 | **+67.8%** | +58.0% | 101 → 54 | **49/49** | 98.3 → 99.7% | 23× |
| *E. cloacae* | 53 | 44× | 215,275 | 321,318 | **+49.3%** | +35.7% | 157 → 99 | 46/53 | 98.5 → 99.7% | 16× |
| *S. aureus* | 28 | 55× | 215,947 | 307,779 | **+42.5%** | +50.9% | 44 → 23 | 26/28 | 98.1 → 99.7% | 12× |
| *E. faecium* | 47 | 51× | 35,245 | 45,699 | **+29.7%** | +26.7% | 175 → 136 | **47/47** | 94.5 → 99.1% | 62× |
| *E. coli* | 60 | 27× | 102,892 | 128,781 | **+25.2%** | +23.0% | 215 → 154 | 59/60 | 93.9 → 98.4% | 43× |
| **all** | **279** | — | — | — | — | — | **762 → 514 (−33%)** | **269/279 (96%)** | — | — |

## Three results, all measured

**Contiguity improves everywhere.** Mean contig NG50 gains of +25% to +70%; 269 of 279 isolates
(96%) improved individually.

**Misassemblies fall everywhere: 762 → 514, −33%.** This was the opposite of my prediction, made
twice. Simulated reads had suggested the model *adds* 0–6 misassemblies per species.

**Genome fraction rises everywhere**, most where it was worst: *E. coli* 93.9 → 98.4%,
*E. faecium* 94.5 → 99.1%. The model is not trading completeness for contiguity — it recovers
more of the genome as well as ordering it better.

## Misassembly composition

| species | inversion | relocation | translocation |
|---|---|---|---|
| *S. aureus* | 14 | 7 | 2 |
| *E. faecium* | 28 | 77 | 31 |
| *A. baumannii* | 25 | 21 | 2 |
| *P. aeruginosa* | 27 | 25 | 2 |
| *E. cloacae* | 13 | 79 | 7 |
| *E. coli* | 71 | 67 | 16 |

Inversions preserve every base; relocations displace sequence; translocations mix replicons.
*E. faecium* stands out with 31 translocations and 77 relocations — consistent with its very
large plasmid complement (1,104 plasmid records in its panel, by far the most) and its 62×
scaffold inflation. It is the organism where this approach works least well and needs most care.

## What is provisional

Reads for these were fetched with a coverage cap of ~100× — an optimisation that should not have
existed, since real coverage spans 24×–267× and that spread is the condition being measured. The
uncapped re-run was 565/720 files in when work stopped, so the **median coverage column above
(25×–55×) reflects capped depth for some isolates and full depth for others**. Direction and
ranking are robust; exact magnitudes will move on a clean uncapped run.

*S. aureus* n = 28 rather than ~60 because it had the most high-coverage runs and therefore lost
the most isolates when the cap was removed for refetching.

## Traps that these numbers avoid

- **contig NG50, not scaffold NG50.** Scaffold inflation runs 12×–62× here. Reported as scaffold
  NG50, *E. faecium* — the weakest model — would look like the best, showing ~2.8 Mb instead of
  45,699 bp.
- **Leave-cluster-out.** Every panel genome within mash d < 0.0005 of a test isolate was excluded
  from training. For *S. aureus* that removed 22 of 137 panel genomes; without it the model has
  effectively seen the answer.
- **Real reads, not simulated.** wgsim inverted the species ranking entirely and reversed the
  misassembly finding.
