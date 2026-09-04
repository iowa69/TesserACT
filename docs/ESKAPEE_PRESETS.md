> ## UPDATE — six trained models now exist
>
> This page was written when only *Klebsiella* had a model and the measured answer was that
> parameter tuning does nothing. **That finding stands**: no tuning flag beats the defaults on any
> ESKAPEE organism. What has changed is that models were then built for the other six, from
> panels of 124–157 mash-dereplicated closed chromosomes, and tested on genomes held out of the
> panel before training:
>
> | organism | contig NG50 gain | misassembly cost |
> |---|---|---|
> | *S. aureus* | **+86.1%** | none |
> | *A. baumannii* | +54.6% | +3 |
> | *E. faecium* | +35.7% | +3 |
> | *P. aeruginosa* | +27.1% | +5 |
> | *E. coli* | +23.8% | +5 |
> | *E. cloacae* | +16.6% | none |
>
> Every preset now selects its model when one is installed. See `ESKAPEE_MODELS.md` for the full
> table, the held-out design, and the scaffold-inflation column — reported as contig NG50
> throughout, because scaffold NG50 would overstate these by 8.9× to 43.9×.

# `tesseract-eskape` presets — what was measured, and what we ship

TesserACT 1.2.5 · presets frozen 2026-09 · seven ESKAPEE organisms

---

## Read this first

**Six of the seven presets change nothing.** They pass no extra flags to the
assembler. That is not an oversight and it is not a placeholder — it is the
measured result. We tried every tuning knob the assembler exposes on three
real closed genomes per species, and on six of the seven species not one knob
made the assembly meaningfully better than leaving the settings alone. Several
made it worse. The honest preset for those six is *use the defaults*.

The preset name still works for all seven, so you can write the same command
whatever you are assembling. Only *Klebsiella pneumoniae* actually adds
something, and what it adds is a **trained model file** — a data asset built
from thousands of finished *Klebsiella* genomes — not a tuning flag.

**Everything below was measured on SIMULATED reads.** See
[Why the real-world gain will be smaller](#why-the-real-world-gain-will-be-smaller).

---

## The table

Genome size and GC are medians of the three genomes tested per species.
NG50 is a contiguity measure: **half of the true genome sits in pieces this
long or longer.** Bigger is better. It is computed against the *true* genome
length, not the assembly length, so an assembler cannot inflate it by emitting
junk. All NG50 values below are **contig** NG50 — see
[the metric footnote](#the-metric-that-nearly-fooled-us-contigs-vs-scaffolds).

| Organism | Preset name | Genome | GC | Flags the preset adds | Default NG50 | Preset NG50 | Beats default? |
|---|---|---:|---:|---|---:|---:|:--|
| *Klebsiella pneumoniae* | `kpneumoniae` | 5.28 Mb | 57.5 % | `--organism klebsiella --model <model>.tsm` | 325,336 bp | **409,059 bp** | **YES — +25.7 %** |
| *Staphylococcus aureus* | `saureus` | 3.05 Mb | 32.8 % | *(none — defaults)* | 249,909 bp | 249,909 bp | no |
| *Enterococcus faecium* | `efaecium` | 2.98 Mb | 37.9 % | *(none — defaults)* | 47,639 bp | 47,639 bp | no |
| *Acinetobacter baumannii* | `abaumannii` | 4.16 Mb | 39.0 % | *(none — defaults)* | 160,376 bp | 160,376 bp | no |
| *Pseudomonas aeruginosa* | `paeruginosa` | 7.06 Mb | 65.8 % | *(none — defaults)* | 284,592 bp | 284,592 bp | no |
| *Enterobacter cloacae* | `ecloacae` | 5.32 Mb | 54.8 % | *(none — defaults)* | 211,541 bp | 211,541 bp | no |
| *Escherichia coli* | `ecoli` | 5.15 Mb | 50.8 % | *(none — defaults)* | 123,506 bp | 123,506 bp | no |

The *Klebsiella* row is measured on **three held-out genomes** that were never
used to choose the setting. The other six rows are the tuning-set means over
three genomes each; since nothing beat the default there was nothing to hold
out and validate.

### What we tried and rejected, per organism

The acceptance bar was **+10 % mean NG50 with no increase in the misassembly
proxy.** Best non-default configuration found, per species:

| Organism | Best alternative tried | Effect on NG50 | Verdict |
|---|---|---:|---|
| *S. aureus* | `--mode aggressive` | **−8.4 %** | rejected; also raised misassembly proxy 10.3 → 13.7 |
| *E. faecium* | `--mode aggressive` | −0.07 % | rejected; zero misassemblies, but inflated assembly to 101.4 % of true length |
| *A. baumannii* | `--mode careful` | ±0.0 % (exact tie on all 3 genomes) | rejected; `--mode aggressive` was −7.5 % |
| *P. aeruginosa* | `--mode aggressive` | −0.18 % | rejected; misassembly proxy 18.3 → 19.3 |
| *E. cloacae* | `--mode aggressive` | +2.87 % on tuning set | **rejected — did not reproduce**, see below |
| *E. coli* | `--mode aggressive` | +0.08 % | rejected; misassembly proxy 0.33 → 2.33 (7×) |

Flags tested across all seven species: `--mode careful`, `--mode aggressive`,
`--simplify-rounds 6`, `--bubble-coverage 0.6`, `--min-link 4`,
`--min-link 1 --link-per-x 0.01`, `--tie-ratio 1.0`. Every one of them was
accepted by the assembler (none were silently dropped — we checked by feeding
it a bogus flag, which it rejects with an error). Several were **exact
no-ops**: `--simplify-rounds 6` produced a byte-identical `contigs.fasta` to
the default on every *S. aureus* genome and on 2 of 3 *E. coli* genomes,
because graph simplification already converges before round 6.

### The one we shipped and then un-shipped

*E. cloacae* with `--mode aggressive` looked like a small win on the three
tuning genomes (+2.87 %). It was then re-tested on three genomes it had never
seen, and the sign **reversed**: +10.0 %, −25.2 %, 0.0 % on the individual
genomes, **−9.2 % on the mean**, with the misassembly proxy rising from 0.00
to 0.67 and assembly length inflating past the true genome length. A +2.9 %
effect measured on three genomes was noise. `ecloacae` ships defaults.

This is the reason the bar was +10 % and not +2 %.

---

## The metric that nearly fooled us: contigs vs scaffolds

The first version of this table claimed the *Klebsiella* preset raised NG50
from 560,263 bp to 5,277,836 bp — **a 9.4× gain**. That number was wrong, and
it is worth explaining exactly how, because the same trap catches published
assembly papers.

The model produces one chromosome-length sequence per genome. But that
sequence is a **scaffold**, not a contig: it contains runs of the letter `N`,
which mean *"we believe some DNA goes here and we know roughly how much, but
we did not read it."* The model's scaffolds carried 46,000–57,000 N's each.
The default assemblies carried 0, 1 and 0 N's. Comparing an N-padded scaffold
against N-free contigs is comparing two different things.

Split the sequences at every run of ≥10 N's — that is, ask *how long are the
pieces we actually read* — and the honest numbers are:

| Metric | Default | Klebsiella preset | Ratio |
|---|---:|---:|---:|
| **Contig** NG50, held-out genomes | 325,336 bp | 409,059 bp | **1.26×** ← *what we ship* |
| **Contig** NG50, tuning genomes | 560,263 bp | 669,408 bp | 1.19× |
| *Scaffold* NG50, held-out genomes | 325,336 bp | 5,247,330 bp | *16.1× — do not quote this* |

So the gain is real, it holds on genomes the setting was never tuned on, and
it clears the +10 % bar on **every individual held-out genome** (+56.0 %,
+15.0 %, +13.7 %) — but it is about a quarter, not an order of magnitude.

Two further honest notes on the *Klebsiella* preset:

* **Where the ~46 kb of N comes from.** The model orders the contigs by
  comparing them to a panel of 2,221 finished *Klebsiella* chromosomes. The
  assembler's own log says it plainly: *"21 joins (21 chromosomal) spanning
  45,862 gap bases."* That order is **inferred from other genomes, not
  observed in your reads.** It is documented, intended behaviour, and you can
  switch it off with `--no-layout` if that distinction matters for your work
  (outbreak analysis, structural-variant calling). We checked the order is
  correct on held-out data — the biggest scaffold aligns to its reference as a
  single collinear block — but "correct on three genomes" is not "always
  correct".
* **It is not a 60× trick.** Tested at 30×, 60× and 100× read depth, the
  preset clears the bar at all three, and helps *most* at low coverage
  (+93.7 % at 30×, +56.0 % at 60×, +49.7 % at 100×). The model output barely
  changes with depth; the gain shrinks only because the default gets better
  with more reads.

---

## Why the real-world gain will be smaller

Every number on this page came from **simulated reads** — generated with
`wgsim`, not produced by a sequencer. Simulated reads are much easier to
assemble than real ones:

* coverage is perfectly uniform; real coverage dips in GC-rich and GC-poor regions
* there are no adapter sequences, no PCR duplicates, no optical duplicates
* there is no contamination, no host DNA, no cross-barcode bleed
* the error model is a flat 0.5 % substitution rate with **no indels at all**
* there is no true within-sample variation

So **every gain quoted here is an upper bound.** The *Klebsiella* +25.7 %
should be read as "at most about a quarter, probably less on your data." And
for the six default-shipping species, an upper bound of zero is a fairly
strong negative result: if a knob cannot help on easy data, it is unlikely to
rescue hard data.

---

## Everything else you should know before trusting this

* **n = 3 genomes per species**, one random seed. Small. The *E. cloacae*
  reversal above is exactly what n = 3 does to you. Treat differences under
  ~10 % as noise; that is why the bar was set there.
* **The "misassembly proxy" is not QUAST.** It aligns contigs back to the
  source genome with `minimap2 -x asm10`, keeps only primary alignments
  (`tp:A:P`) and blocks ≥1 kb, and flags a contig whose pieces land on
  different references, on both strands, out of order, or >10 kb apart, with
  an exemption for contigs crossing the circular origin. The **absolute**
  counts are not real misassembly counts; only the comparison *between*
  configurations is meaningful, and all configurations were scored with the
  same code. It cannot see small indels or local rearrangements that minimap2
  absorbs into one alignment block, and it cannot detect a wrong-but-plausible
  join that the *Klebsiella* panel also contains.
* **Only *Klebsiella* has a trained model.** For the other six, the
  model-guided layout step that produces the *Klebsiella* gain simply does not
  run. The single largest available improvement for *A. baumannii*,
  *P. aeruginosa* and the rest is **a trained model for that genus**, not a
  flag preset. That is a data-collection project, not a tuning project.
* **The contiguity ceiling here is physics, not settings.** With 150 bp reads
  and a ~350 bp insert, any repeat longer than the insert — rRNA operons, IS
  elements, transposases — breaks the assembly graph and no threshold knob can
  span it. That is why loosening the join thresholds (`--min-link 1
  --link-per-x 0.01`) produced *byte-identical output* to the default on
  *A. baumannii* and *E. coli*. The contigs are not join-limited; they are
  repeat-limited.
* **Reference sets were chromosome-only**, so plasmid contigs in an assembly
  count toward total length but cannot be validated by the alignment. This
  affects the length-vs-truth percentages, not the NG50 comparison.
* **Reproducing this:** reads were simulated as
  `wgsim -N (bp*60/300) -1 150 -2 150 -d 350 -s 40 -e 0.005 -r 0 -R 0 -X 0 -S 11`,
  assembled with `tesseract-asm -t 3` run strictly serially. Across the sweep
  and the held-out verification, every assembly exited 0; no run failed and no
  flag was rejected.

---

## What this means for you, in one paragraph

Run `tesseract-eskape` with the preset for your organism and do not think
about it any further. For six of the seven, the preset is the assembler's own
default, which we checked is already the best available setting — you are not
missing a magic flag, because there isn't one. If you work on *Klebsiella
pneumoniae*, get the model file: it is worth roughly a quarter more
contiguity, and more than that if your coverage is low. If you need to know
which bases were actually sequenced rather than inferred from other genomes,
add `--no-layout`.
