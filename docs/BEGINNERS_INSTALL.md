# Installing `tesseract-eskape`

This guide assumes you have never installed a piece of bioinformatics software
before. Every unfamiliar word is explained the first time it appears.

---

## 0. What you are installing, and why there are two pieces

**`tesseract-asm`** is the *assembler*. Sequencing machines do not read a
genome from end to end; they read millions of short fragments a few hundred
letters long. An assembler is the program that puts those fragments back
together into long stretches of DNA. `tesseract-asm` is the real work.

**`tesseract-eskape`** is the small script in this folder. It is a *wrapper*:
it does not assemble anything itself, it just checks your files, picks the
right settings for your organism, and then calls `tesseract-asm` for you. It
exists so that you can type one short command instead of remembering flags.

You need both. Install the assembler first, then drop the wrapper next to it.

---

## 1. The one-liner (if you already have conda)

**conda** is a program that installs scientific software and all the other
software it depends on, without you having to compile anything. If you already
have it, this is the whole installation:

```bash
conda create -n tesseract -c conda-forge -c bioconda tesseract-assembler
conda activate tesseract
```

Line by line:

* `conda create -n tesseract` — make a new, isolated *environment* (a private
  folder of software) called `tesseract`. Environments keep projects from
  breaking each other.
* `-c conda-forge -c bioconda` — the two *channels* (software repositories)
  that the package and its dependencies live in. Order matters; keep it as
  written.
* `tesseract-assembler` — the package name. Note it is **not** the same as the
  command name, which is `tesseract-asm`.
* `conda activate tesseract` — step into that environment. **You must run this
  in every new terminal window** before the command will be found.

### Honest caveat

The conda recipe for this package lives in the source tree at
`conda/meta.yaml` and declares version 1.2.5. I could **not** verify from this
machine that the package is actually published on bioconda yet — the network
was not available to check. If the command above fails with

```
PackagesNotFoundError: The following packages are not available from current channels: tesseract-assembler
```

then it has not been published, and that is not your fault. Use section 3
(build from source) instead. It is only two commands longer.

---

## 2. If you do not have conda

Check first — you may already have it:

```bash
conda --version
```

If that prints something like `conda 26.3.2`, skip to section 1.
If it prints `command not found`, install **Miniforge**, a small, free,
no-strings-attached conda distribution:

**Linux (most servers and HPC clusters):**

```bash
curl -L -o Miniforge3.sh "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
bash Miniforge3.sh -b -p "$HOME/miniforge3"
"$HOME/miniforge3/bin/conda" init bash
```

**macOS (Apple Silicon — M1/M2/M3/M4):**

```bash
curl -L -o Miniforge3.sh "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-arm64.sh"
bash Miniforge3.sh -b -p "$HOME/miniforge3"
"$HOME/miniforge3/bin/conda" init zsh
```

Then **close your terminal and open a new one.** This is not optional — the
`init` step edits your shell's start-up file, and that file is only read when a
terminal starts. Now `conda --version` should work, and you can go back to
section 1.

**Windows:** the assembler is a Unix program. Install *WSL* (Windows Subsystem
for Linux) from the Microsoft Store — it gives you a real Ubuntu terminal
inside Windows — then follow the Linux instructions inside it.

**No admin rights?** Everything above installs into your own home folder and
needs no administrator or `sudo`. That is why we use Miniforge rather than a
system package manager.

---

## 3. If conda will not work: build from source

The assembler needs almost nothing: a C++ compiler and zlib. If you have a
copy of the source tree, this is the whole build:

```bash
cd /path/to/TesserACT
./install.sh
```

By default that installs into your active conda environment, or into
`~/.local/bin` if you are not in one. To choose the location yourself:

```bash
./install.sh --prefix "$HOME/tesseract"
export PATH="$HOME/tesseract/bin:$PATH"
```

The second line tells your shell where to look for commands. To make it stick,
add that same line to the end of `~/.bashrc` (Linux) or `~/.zshrc` (macOS).

If `./install.sh` complains **`no C++ compiler found`**, install one:
`conda install -c conda-forge cxx-compiler`, or on Ubuntu
`sudo apt install build-essential`.

---

## 4. Install the wrapper

The wrapper is the file `tesseract-eskape` sitting next to this document.
Copy it somewhere on your `PATH` (the list of folders your shell searches for
commands) and make it executable:

```bash
mkdir -p "$HOME/.local/bin"
cp tesseract-eskape "$HOME/.local/bin/"
chmod +x "$HOME/.local/bin/tesseract-eskape"
export PATH="$HOME/.local/bin:$PATH"
```

`chmod +x` means "mark this file as a program that may be run". Without it you
get `Permission denied`.

You do not have to do any of this — you can always run it in place as
`./tesseract-eskape` from the folder it is in.

### Telling the wrapper where the assembler is

The wrapper looks for `tesseract-asm` at a built-in path. If yours is
elsewhere (very likely, if you installed it yourself), tell it:

```bash
export TESSERACT_ASM="$(command -v tesseract-asm)"
```

Add that line to `~/.bashrc` so you do not have to repeat it. If the path is
wrong, the wrapper says so clearly before it starts anything slow — it will
not fail halfway through an hour-long run.

---

## 5. The Klebsiella model file (only if you work on *Klebsiella pneumoniae*)

*Klebsiella pneumoniae* is the one organism with a **trained model** — a data
file, built from thousands of finished *Klebsiella* genomes, that helps the
assembler put the pieces in the right order. It is the only preset in this
package that changes anything (see `PRESETS.md`).

The model is a large file (a few hundred megabytes) and is **not** part of the
software install. Point the wrapper at your copy:

```bash
export TESSERACT_KP_MODEL=/path/to/fold0_mem.tsm
```

or pass it per-run with `--model /path/to/fold0_mem.tsm`.

If you do not have the model, do not worry: run without `--preset` (or with
any other preset) and you get a perfectly normal assembly. You just do not get
the Klebsiella contiguity bonus. The wrapper tells you this rather than
failing silently.

The other six organisms need no model file at all. Nothing to download.

---

## 6. Check that it worked

Run these three commands. All three should finish instantly and print
something sensible.

```bash
tesseract-asm --version
```
Expected: `TesserACT 1.2.5` and an author line. **Verified on this machine.**

```bash
tesseract-eskape --list
```
Expected: a table of the seven organism presets. **Verified on this machine.**

```bash
tesseract-eskape --help
```
Expected: a help page with a copyable example. **Verified on this machine.**

### The real test: assemble something

If you have a pair of read files, do a dry run first. `--dry-run` checks all
your files and prints the exact command it *would* run, then stops — so you
find out about a typo in one second instead of ten minutes.

```bash
tesseract-eskape --preset saureus \
    -1 reads_R1.fastq.gz \
    -2 reads_R2.fastq.gz \
    -o test_result --dry-run
```

Happy with what it printed? Remove `--dry-run` and run it for real:

```bash
tesseract-eskape --preset saureus \
    -1 reads_R1.fastq.gz \
    -2 reads_R2.fastq.gz \
    -o test_result -t 4
```

A typical bacterial genome takes a few minutes on 4 cores and about 2 GB of
memory. When it is done, your assembled genome is the file
`test_result/contigs.fasta`, and the wrapper prints that path for you.

Check it is not empty:

```bash
grep -c '^>' test_result/contigs.fasta
```

That counts the sequences. A good short-read bacterial assembly is typically
tens to a few hundred sequences. `0` means something went wrong — read
`test_result/tesseract-eskape.log`, which holds everything the assembler said.

---

## 7. Words you just met

| Word | What it means |
|---|---|
| **read** | one short DNA fragment (~150 letters) off the sequencer |
| **FASTQ** | the file format holding reads. Usually ends `.fastq.gz` or `.fq.gz` |
| **R1 / R2** | the two halves of a *paired* read. Two files, always used together |
| **gz** | the file is compressed. Leave it compressed; the tools read it directly |
| **contig** | one continuous stretch of assembled DNA |
| **FASTA** | the file format holding contigs. The output, `contigs.fasta` |
| **assembly** | the whole set of contigs — your reconstructed genome |
| **environment** | a private folder of installed software, managed by conda |
| **PATH** | the list of folders your shell searches when you type a command |
| **preset** | a named bundle of settings, here one per organism |

---

## 8. Common problems

| What you see | What it means | What to do |
|---|---|---|
| `command not found: tesseract-eskape` | not on your `PATH` | run it as `./tesseract-eskape`, or redo section 4 |
| `command not found: conda` | conda not installed, or you did not open a new terminal | section 2 |
| `PackagesNotFoundError` | the package is not on that channel | build from source, section 3 |
| `I cannot find the assembler program at: ...` | wrapper does not know where `tesseract-asm` is | `export TESSERACT_ASM=$(command -v tesseract-asm)` |
| `File not found: ... check the spelling` | typo in a filename | `ls` the folder; use Tab to autocomplete names |
| `Permission denied` | file not marked runnable | `chmod +x tesseract-eskape` |
| `You gave the same file for -1 and -2` | you passed R1 twice | find the matching `_R2` file |
| assembler killed, no message | ran out of memory | use a machine with more RAM, or `--max-memory 8` |
