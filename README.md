# TesserACT

## What it is

TesserACT is a de novo genome assembler for **bacterial isolates** sequenced with
**paired-end Illumina reads**. One command takes your two FASTQ files and returns the
assembled genome:

- it **checks the reads first**: complete files, valid FASTQ, mates in the same order,
  depth, adapters, poly-G tails. Problems are reported before any time is spent assembling;
- it **assembles** them: de Bruijn graphs over a ladder of k-mer sizes, repeat resolution
  with read pairs, gap closing and polishing;
- it **writes circular replicons** (plasmids, sometimes the chromosome) only when the reads
  prove the circle closes;
- it **checks the result** and tells you in plain language what to look at.

Its special touch is the **organism models**. For *Klebsiella pneumoniae*, *E. coli*,
*Enterobacter cloacae*, *Acinetobacter baumannii*, *Pseudomonas aeruginosa*,
*Staphylococcus aureus*, *Enterococcus faecium* and *Salmonella enterica*, a model built
from thousands of closed genomes of the species can lay out your chromosome along the
conserved gene order of its relatives, joining pieces that short reads alone cannot place.
Models are optional and used only when you ask.

## How to install it

With conda (Linux and macOS):

```sh
conda install -c conda-forge -c bioconda tesseract-assembler
```

The package is called `tesseract-assembler`; the command it installs is `tesseract`.

From source (needs a C++17 compiler, make and zlib):

```sh
git clone https://github.com/iowa69/TesserACT.git
cd TesserACT
./install.sh
```

`install.sh` builds TesserACT, puts `tesseract` on your PATH (it asks where; nothing needs
administrator rights), and offers to download the organism models.

Check it works:

```sh
tesseract --version
```

## How to run it

```sh
tesseract -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz -t 16 -o sample
```

That is a complete, vanilla assembly: built only from your reads, nothing from other
genomes. `-t` is the number of CPU threads (default: all of them); `-o` is a new output
folder. The genome is in **`sample/contigs.fasta`**, and **`sample/SUMMARY.txt`** lists
what was checked and anything worth a look.

To check the reads without assembling:

```sh
tesseract check -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz
```

Many samples:

```sh
for r1 in reads/*_R1.fastq.gz; do
    s=$(basename "$r1" _R1.fastq.gz)
    tesseract -1 "$r1" -2 "reads/${s}_R2.fastq.gz" -t 16 -o "assemblies/$s"
done
```

## How to use the options

### Organism models

Download once (all of them, about 1.7 GB, or only the ones you need):

```sh
tesseract models download                 # all eight
tesseract models download kpneumoniae     # just one
tesseract models list                     # what is available and installed
```

Then name the organism:

```sh
tesseract -1 R1.fastq.gz -2 R2.fastq.gz -t 16 -o sample --organism kpneumoniae
```

| `--organism` | species |
|---|---|
| `kpneumoniae` (or `klebsiella`) | *Klebsiella pneumoniae* |
| `ecoli` | *Escherichia coli* |
| `ecloacae` | *Enterobacter cloacae* complex |
| `abaumannii` | *Acinetobacter baumannii* |
| `paeruginosa` | *Pseudomonas aeruginosa* |
| `saureus` | *Staphylococcus aureus* |
| `efaecium` | *Enterococcus faecium* |
| `senterica` (or `salmonella`) | *Salmonella enterica* |

With a model the chromosome comes out in fewer, longer pieces. The order of some joins
then comes from the species' relatives rather than from your reads, so a few can be wrong
where your isolate differs from them. Use the vanilla assembly when every join must be
backed by your own reads, and the model when contiguity matters most. If the reads look
like one of these species, the vanilla `SUMMARY.txt` says so.

Models go to `~/.tesseract/models`; set `TESSERACT_MODEL_DIR` to keep them elsewhere.

### Everyday options

| option | use it to |
|---|---|
| `-t N` | set the number of threads |
| `--genome-size MB` | give the expected genome size, for better depth and size checks (set automatically by `--organism`) |
| `--min-contig N` | report only contigs of at least N bp (default: about 2 x the largest k) |
| `--max-memory GB` | stop cleanly above this much memory (default: 80% of RAM) |

### Assembly modes

| option | what it does |
|---|---|
| `--mode standard` | the default |
| `--mode fast` | fewer k-mer sizes and no polishing: quicker, a little more fragmented |
| `--mode careful` | more k-mer sizes, more graph simplification, two polishing passes |
| `--mode aggressive` | merges diverged repeat copies for longer contigs, at a risk of about one wrong join per genome |
| `--map-polish bwa` (or `bowtie2`) | polish the final contigs against a full read alignment (needs the aligner installed) |

`tesseract --help` lists the main options; `tesseract --help-all` lists every assembler
setting (k-mer sizes, cutoffs, repeat resolution, scaffolding, gap filling, polishing).

### What the run does for you

You do not need to tune anything for the usual cases; `tesseract` adapts to the data:

- **any read length** from 2x100 to 2x300: the k-mer sizes follow the read length;
- **very deep data and small multicopy plasmids**: reads from replicons at many times the
  genome depth are thinned before assembly, so small plasmids are recovered whole;
- **very thin data** (about 10x): if the standard k-mer sizes assemble almost nothing, it
  retries with smaller ones and says so;
- **unusual libraries** (mate-pair, wrong orientation): detected, and read-pair joins are
  switched off for that run;
- **PhiX spike-in**: PhiX contigs are labelled `_spikein` so you can remove them.

### The output

| file | content |
|---|---|
| `contigs.fasta` | the assembly |
| `scaffolds.fasta` | contigs joined across gaps of estimated size (only when such gaps remain) |
| `SUMMARY.txt` | every check, as `OK` / `WARN` / `FAIL`, with what to do |
| `report.html` | the whole run, for a browser; `report.json` for scripts |
| `assembly_graph.gfa` | the assembly graph, for Bandage |
| `replicons.tsv` | each circular candidate and why it was or was not closed |
| `ends.tsv` | why each contig ends where it does |
| `p2_self_qa.tsv` | how much of the read content made it into the assembly |

Contig names end with what TesserACT knows about them: `_circular` (a closed replicon,
verified by reads), `_spikein` (PhiX), and with a model `_chr`, `_plas` or `_unk`
(chromosome, plasmid, unknown). Lowercase bases at a contig end are bases the reads there
support only weakly.

Exit status 0 means every check passed (warnings are in `SUMMARY.txt`); 1 means a check
failed. If the read check fails, nothing is assembled; `--force` assembles anyway.

## Citation and licence

Please cite TesserACT as described in [`CITATION.cff`](CITATION.cff). MIT licence, see
[`LICENSE`](LICENSE).
