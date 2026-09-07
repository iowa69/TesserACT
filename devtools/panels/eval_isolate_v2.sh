#!/usr/bin/env bash
# eval_isolate_v2.sh ORG SAFE_ACC RUN OUTDIR [MODEL:NAME ...]
#
# As eval_isolate.sh, plus a vanilla SPAdes run off the SAME reads in the SAME
# download pass. Fetching an isolate's FASTQ once and running every assembler on
# it is the only way to keep peak read storage at one isolate while still
# comparing tools, and it also guarantees both assemblers saw byte-identical input.
#
# "Vanilla" is taken from the repo's own convention (README head-to-head):
#   "Both assemblers ran on raw untrimmed reads with default settings -- SPAdes
#    4.3.0 with its own error correction on and automatic k selection, TesserACT
#    with no model."
# So: no --isolate, no --careful, no pre-trimming, no k list. Adding --isolate or
# --careful would be a tuned SPAdes, and comparing tuned-ours against tuned-theirs
# is a different claim than the one the repo makes.
#
# Scoring pairs like with like: SPAdes' contigs.fasta is unscaffolded, so it is
# compared against TesserACT's N-split (_broken) row. QUAST -s produces both.
set -euo pipefail
org="$1"; safe="$2"; run="$3"; outdir="$4"; shift 4
root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HOME/miniconda3/etc/profile.d/conda.sh"; conda activate tesseract
asm="$root/../TesserACT-main/tesseract-asm"
threads="${THREADS:-8}"
spades_mem="${SPADES_MEM:-12}"
mkdir -p "$outdir"

ref=$(python3 - "$root/genomes/$org" "$safe" <<'PY'
import glob, os, sys
raw, safe = sys.argv[1], sys.argv[2]
for d in glob.glob(os.path.join(raw, "batch.*.d", "ncbi_dataset", "data", "GC*")):
    if os.path.basename(d).replace("_", "").replace(".", "v") == safe:
        f = glob.glob(os.path.join(d, "*.fna"))
        if f:
            print(f[0]); break
PY
)
[ -n "$ref" ] && [ -s "$ref" ] || { echo "  no reference for $safe" >&2; exit 1; }

names=(base)
for spec in "$@"; do names+=("${spec##*:}"); done
names+=(spades)

need=0
for n in "${names[@]}"; do [ -s "$outdir/$n/contigs.fasta" ] || need=1; done

if [ "$need" = 1 ]; then
    rd="$outdir/reads"; mkdir -p "$rd"
    if [ ! -s "$rd/${run}_1.fastq" ]; then
        ok=0
        for try in 1 2 3; do
            fasterq-dump --split-3 -e "$threads" -O "$rd" -t "$rd" "$run" >/dev/null 2>&1 && { ok=1; break; }
            sleep 15
        done
        [ "$ok" = 1 ] || { echo "  fasterq-dump failed $run" >&2; rm -rf "$rd"; exit 1; }
    fi
    [ -s "$rd/${run}_1.fastq" ] && [ -s "$rd/${run}_2.fastq" ] \
        || { echo "  $run not paired" >&2; rm -rf "$rd"; exit 1; }

    for spec in base "$@"; do
        if [ "$spec" = base ]; then n=base; extra=()
        else n="${spec##*:}"; extra=(--organism "$org" --model "${spec%%:*}"); fi
        [ -s "$outdir/$n/contigs.fasta" ] && continue
        "$asm" -1 "$rd/${run}_1.fastq" -2 "$rd/${run}_2.fastq" -o "$outdir/$n" \
               -t "$threads" "${extra[@]}" > "$outdir/$n.log" 2>&1 \
            || echo "  assembly failed: $n" >&2
    done

    if [ ! -s "$outdir/spades/contigs.fasta" ]; then
        # SPAdes writes a large working tree; only contigs/scaffolds are kept below.
        conda run -n spades spades.py -1 "$rd/${run}_1.fastq" -2 "$rd/${run}_2.fastq" \
            -o "$outdir/spades_work" -t "$threads" -m "$spades_mem" \
            > "$outdir/spades.log" 2>&1 || echo "  spades failed" >&2
        if [ -s "$outdir/spades_work/contigs.fasta" ]; then
            mkdir -p "$outdir/spades"
            cp "$outdir/spades_work/contigs.fasta" "$outdir/spades/contigs.fasta"
            [ -s "$outdir/spades_work/scaffolds.fasta" ] && \
                cp "$outdir/spades_work/scaffolds.fasta" "$outdir/spades/scaffolds.fasta"
        fi
        rm -rf "$outdir/spades_work"
    fi
    rm -rf "$rd"
fi

for n in "${names[@]}"; do
    [ -s "$outdir/$n/contigs.fasta" ] || continue
    [ -s "$outdir/quast_$n/report.tsv" ] && continue
    conda run -n quast quast.py -s -t "$threads" -R "$ref" --min-contig 500 \
        -o "$outdir/quast_$n" "$outdir/$n/contigs.fasta" >/dev/null 2>&1 || true
done
echo "  done $safe"
