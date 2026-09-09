# Why SPAdes assembles more plasmids whole, and what it costs us

The *S. aureus* panel gives TesserACT eight of nine assembly metrics. It also shows a
plasmid deficit that is **not** explained by the contiguity trade, because it runs the
other way: SPAdes is less contiguous overall and yet delivers more plasmids in one piece.

| | TesserACT `--organism` | SPAdes |
|---|---|---|
| reference plasmids recovered (>=50% covered) | 111/119 (93.3%) | 106/115 (92.2%) |
| **recovered whole (>=95% in one contig)** | **38 (32%)** | **49 (41%)** |

We find the sequence. We do not put it in one contig.

## The deficit is entirely a copy-number effect

Every reference plasmid in the panel — 119 of them across 93 plasmid-bearing isolates —
scored by the depth of our own contigs covering it, relative to the median depth of our
contigs on that isolate's chromosome:

| plasmid depth / chromosome depth | n | TesserACT whole | SPAdes whole | median contigs, TesserACT | SPAdes |
|---|---|---|---|---|---|
| **< 1.6x** | 49 | **18 (37%)** | 16 (33%) | 1 | 1 |
| 1.6 – 3x | 35 | 12 (34%) | 15 (43%) | 3 | 2 |
| 3 – 8x | 20 | 7 (35%) | 9 (45%) | 2 | 2 |
| **>= 8x** | 15 | **1 (7%)** | **9 (60%)** | **6** | **1** |
| total | 119 | 38 (32%) | 49 (41%) | | |

Below 1.6x **we beat SPAdes**. At or above 8x we recover one plasmid whole out of fifteen
where SPAdes recovers nine, and our median contig count per plasmid goes to six while
theirs stays at one. The entire 11-plasmid deficit is in the bands at or above 1.6x, and
eight of those eleven are in the >=8x band alone.

The threshold is not a coincidence. `resolve.cpp`:

```cpp
const double repeatThreshold = medianCoverage_ * 1.6;
auto isRepeat = [&](uint32_t u) {
    return medianCoverage_ > 0 && g_.nodes[u].coverage > repeatThreshold;
};
```

and a unitig that answers true to that "cannot seed a chain, terminate a path, or be
traversed -- it is emitted verbatim and walls off both neighbours." A plasmid at four
copies per cell is indistinguishable, by depth alone, from a four-copy repeat. So it is
excluded from chain building and leaves the assembler as loose unitigs.

**59% of the panel's reference plasmids sit at or above that 1.6x threshold**, so this is
the common case, not an edge case.

The existing comment beside `repeatThreshold` already says this is happening, and puts the
cost at 4.8% of one assembly including "a 68 kb unitig at 1.78x median that was the whole
of that assembly's NGA50". It records two attempted fixes, both rejected on measurement:
raising the multiplier to 2.4 (frees plasmids, stops genuine two-copy repeats resolving)
and requiring a branching end (no NGA50 change, more misassemblies). Both stay inside the
same evidence — depth and graph shape — and the comment's own conclusion is that
"identifying the risky joins needs a structural signal rather than a depth one".

## The signal the model already has, and was not being asked for

`OrganismModel` stores, per sampled canonical k-mer, how many panel chromosomes and how
many panel plasmids carry that marker exactly once. A marker seen on panel plasmids and
never on a panel chromosome is direct evidence of replicon class, and it is evidence of a
different kind from depth: it does not care how many copies of the plasmid this cell
carries.

`resolve.cpp` had no reference to `OrganismModel` at all. The model is loaded at
`assembler.cpp:778`, *after* the resolver has already run at line 674 — deliberately, so
that a 67–277 MB model is not added to the memory peak of the k-mer ladder.

So the change is not "use the model in the resolver", which would undo that memory
decision. It is to read only the part that answers this one question:

* `OrganismModel::loadExclusiveMarkers` reads the file's marker section and stops there,
  keeping only markers exclusive to one class (1 = plasmid-only, 2 = chromosome-only).
  The adjacency tables and layout tracks behind it are never allocated.
* `PairedResolver::setExclusiveMarkers` takes that map; `isRepeat` exempts a unitig that
  is at least 1,000 bp, carries at least two plasmid-exclusive markers, and carries **no**
  chromosome-exclusive marker.

The rule is deliberately one-sided. A wrong exemption lets a real repeat be traversed, so
a repeat family shared between chromosome and plasmid — which contributes no exclusive
markers either way — is never exempted, and a single chromosome-exclusive marker vetoes.

`TESSERACT_NO_PLASMID_VOUCH=1` disables it for A/B measurement.

## Two attempts, and what the first one taught

**First version: per unitig. It vouched nothing.** The rule asked each high-depth unitig
for at least two plasmid-exclusive markers and no chromosome-exclusive one. On both
isolates tested it fired zero times, and the output was byte-identical to the control --
which at least confirms the change is inert when it does not fire.

Reading the model header directly explains why:

```
magic TSMODEL5   k=31   denom=512
panel chromosomes=984   plasmid sets=3657   excluded accessions=1304
markers=61,891   chr-only 42,776 (69.1%)   plasmid-only 15,563 (25.1%)   both 3,552 (5.7%)
```

**The sampling is one canonical k-mer in 512**, not one in 64. A 1,000 bp unitig
contributes about two sampled k-mers in total, so asking it for two plasmid-exclusive
markers asks for very nearly all of them. The threshold was calibrated against a
denominator the model does not use.

**Second version: per connected group of high-depth unitigs.** This matches the structure
as well as the sampling. A multi-copy plasmid is not one deep unitig; it is a cluster of
them, joined to each other and cut off from the chromosome by the very depth test at issue.
A 20 kb plasmid contributes about 40 sampled k-mers however many pieces it is in.

The veto stays one-sided: a single chromosome-exclusive marker anywhere in the group
disqualifies the whole group, and a repeat family shared between chromosome and plasmid
contributes no exclusive markers either way and is never exempted.

## First measured isolate: it works, and it has a cost

GCF038428225v1 carries three plasmids at 10.6x, 17.1x and 23.4x chromosome depth — the
extreme end of the table above. TesserACT recovered one of the three whole; SPAdes
recovered all three.

| plasmid | length | depth | model | vouch | SPAdes |
|---|---|---|---|---|---|
| NZ_CP131656.1 | 4,397 | 10.6x | 1 contig | 1 contig | 1 contig |
| NZ_CP131657.1 | 2,920 | 17.1x | **2 contigs** (73% in the largest) | **1 contig** | 1 contig |
| NZ_CP131658.1 | 1,491 | 23.4x | **4 contigs** (36% in the largest, 83% covered) | **1 contig** (99% covered) | 1 contig |

Whole-plasmid recovery on this isolate goes **1/3 to 3/3**, matching SPAdes. The rest of
the assembly:

| metric | model | vouch | SPAdes |
|---|---|---|---|
| NGA50 | 227,814 | 227,814 | 153,346 |
| genome fraction | 98.897 | **98.930** | 97.858 |
| duplication ratio | 1.005 | 1.006 | 1.000 |
| **# misassemblies** | **0** | **1** | 0 |
| # local misassemblies | 1 | 2 | 0 |
| total length | 2,847,493 | 2,853,278 | 2,805,223 |

**The cost is real and has an identifiable cause.** One recovered plasmid came out at
2.03x its true length: a circular molecule traversed twice. That is the exact risk of
exempting a unitig from `isRepeat` — the exemption that lets a plasmid be assembled also
lets it be assembled twice — and it is what the extra misassembly and the +5,785 bases of
total length are.

So this is not free. Two plasmids delivered whole against one extra misassembly and one
extra local misassembly, on one isolate.

## Status

Under measurement on 39 isolates: 27 that carry a plasmid at >=3x chromosome depth (where
the change should act) and 12 whose plasmids are all below 1.6x (where the prediction,
recorded before the result, is that **nothing moves at all**). Whichever way it comes out
goes here, including the duplication and misassembly columns, which are the ones that
decide whether this ships.

The doubling is separately fixable — a circular contig whose ends overlap by one period is
detectable — but that is a second change and is not attempted until this one is measured.
