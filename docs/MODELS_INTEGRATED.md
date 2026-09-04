# The integrated ESKAPEE models

Seven organisms, each a single model carrying **both** the chromosome panel (layout tracks) and
the plasmid database, with every held-out isolate's own plasmids withheld from training.

| organism | chromosome panel | plasmid records | held-out plasmids withheld | model |
|---|---|---|---|---|
| *S. aureus* | 133 genomes | 318 | 2 | 14.6 MB |
| *E. faecium* | 121 | 1,080 | 24 | 12.6 MB |
| *A. baumannii* | 127 | 345 | 14 | 22.9 MB |
| *P. aeruginosa* | 153 | 89 | 2 | 38.7 MB |
| *E. cloacae* | 122 | 464 | 9 | 38.6 MB |
| *E. coli* | 120 | 206 | 7 | 25.2 MB |
| *K. pneumoniae* | 2,221 | (existing release model) | — | 338.8 MB |

`--exclude-plasmids` matches the **exact record name**, not an accession prefix. A list of bare
accessions excludes nothing and the build still reports success, so each test plasmid would be
scored against itself — worth 46% of test plasmids on the Klebsiella panel. The lists here are
built from full record names and each build log reports the count actually excluded.

## Chromosome performance (held out of the panel before training)

| organism | base contig NG50 | model contig NG50 | gain | misassemblies |
|---|---|---|---|---|
| *S. aureus* | 210,470 | **391,724** | **+86.1%** | 2 → 2 |
| *A. baumannii* | 181,909 | 281,219 | +54.6% | 1 → 4 |
| *E. faecium* | 49,722 | 67,484 | +35.7% | 0 → 3 |
| *P. aeruginosa* | 269,010 | 341,819 | +27.1% | 0 → 5 |
| *E. coli* | 101,287 | 125,434 | +23.8% | 1 → 6 |
| *E. cloacae* | 236,092 | 275,188 | +16.6% | 0 → 0 |

NG50 is contig NG50 — measured after splitting at runs of ≥10 N, against the true genome length.
Scaffold NG50 would overstate these by 8.9× to 43.9×.

## Plasmid performance, and the test that had to be redone

Plasmid-contig tagging on a held-out *S. aureus* genome carrying a 54 kb plasmid:

| | precision | recall |
|---|---|---|
| **integrated model, realistic copy number** | **0.913** | **0.808** |
| same model, uniform-coverage simulation | 0.200 | 0.042 |

**The first measurement was invalid and nearly became a conclusion.** `wgsim` samples uniformly
across whatever it is given, so simulating from a concatenated chromosome-plus-plasmid file puts
the plasmid at 1× chromosome depth. Copy number is the strongest signal the replicon classifier
has, and the simulation had removed it. Re-simulated with the chromosome at 60× and the plasmid
at 300×, precision and recall go from 0.20/0.04 to 0.91/0.81.

A denser build (`--marker-density 128`) gives slightly better precision (0.947) and worse recall
(0.720). The default density wins on F1 (0.858 against 0.818) and is what ships.

## Clinical impact of the residual misassemblies

58 breakpoints across seven species were scanned against NCBI AMR, ResFinder, VFDB and MEGARes.
**No acquired resistance gene and no virulence gene is disrupted.** Serious genes are present near
breakpoints in Klebsiella — `blaOXA-1`, `blaCTX-M-15`, `aac(6')-Ib-cr`, `dfrA14`, `tet(D)` — and
all are intact, sitting inside the inverted block rather than across its boundary.

The only clinically legible damage: 2 of 42 Klebsiella breakpoints split the **intrinsic**
chromosomal class A beta-lactamase (blaSHV/LEN/OKP family). Anyone calling that gene from a
scaffolded assembly should confirm it against the unscaffolded contigs.

See `../MISASSEMBLY_ANATOMY.md` for why the misassemblies happen — they are prophage and IS
boundaries inherited from the panel track, not invented by the assembler — and for the two fixes
that would address them.

## Installing

```sh
tesseract-get-models            # ~491 MB, checksum-verified, safe to re-run
tesseract-eskape --preset saureus -1 R1.fastq.gz -2 R2.fastq.gz -o result
```
