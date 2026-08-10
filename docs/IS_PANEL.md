# Klebsiella insertion-sequence panel and recurrent sites

Two files, both derived from mapping the 5,970 ISFinder reference elements onto 2,933
closed *Klebsiella pneumoniae* chromosomes: **57,119 IS copies located across 2,910 of
them**.

## `klebsiella_is_panel.fasta` — for `--is-panel`

4,964 distinct sequences, dereplicated from every located copy rather than one reference
per element. That choice is measured, not stylistic: held out on 1,462 genomes, a panel
built from all copies detects a copy at 90-92% identity to its reference with mean k-mer
density **0.997**, while a one-reference-per-element panel drops to **0.481** for copies at
92-94% identity and misses a third of them outright.

A contig end whose 600 bp window is >=30% panel k-mers has its joins silenced. That
threshold separates three cleanly distinct regimes: host DNA fires at 0.0001 (false-fire
rate 0.010%, 29 windows in 290,389), a window straddling an element edge at 0.49, and a
window inside one at 0.997.

## `klebsiella_is_sites.tsv` — for `--is-sites`

183 recurrent insertion sites: two unique flanking 31-mers, the median element length
between them, and how many panel genomes share the site.

Only near-fixed sites are included, and the filter matters. Across the panel, IS copies hit
the same host loci repeatedly -- 3,705 loci are shared by two or more genomes and account
for 69% of the census -- but the median locus is **occupied in only 0.008** of the genomes
carrying its flanks. A site that is bimodal across the species cannot size a join, because
the panel's own genomes disagree about whether anything is there. Conditioning on lineage
instead of species keeps the sites where a clone is consistent.

The yield is deliberately uneven and lands where it is needed: **ST258 gets about 5
placeable loci per genome, ST37 about 0** -- and ST37 already closes at 98.8% while ST258
manages 61.8%.

## Why the veto now applies to the chromosome

It used to run on plasmids only, on the reasoning that chromosomal sites recur and so an
end on one can be placed rather than refused. They recur as *sites* and not as
*insertions*: per family, the chance that a panel chromosome lacks the element at a given
locus runs from 83.5% to 99.9%. The panel's adjacency across such a point therefore
describes a genome without the insertion, and acting on it deletes this isolate's element.

Copy-number stability does not rescue any family. IS1182 is the steadiest in the census --
copy-number CV 0.05 among carriers, never a private insertion -- and the panel still lacks
it at its locus 85% of the time, because it is a lineage marker rather than species
furniture.

## What it buys

Measured over the 666 closed-reference cohort, against the same pipeline without it:

| | chromosome NGA50 | chr misassemblies | plasmid chimeric | plasmid collinear |
|---|---|---|---|---|
| without | 2.84 Mb | 4 | 43 | 97.9% |
| with | 2.82 Mb | 4 | **27** | **99.3%** |

Free on the chromosome, and it halves the plasmid contigs carrying foreign sequence.
