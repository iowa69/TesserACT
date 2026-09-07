#!/usr/bin/env bash
# eval_isolate.sh ORG SAFE_ACC RUN OUTDIR MODEL[:NAME]...
#
# One held-out isolate, one read download, N assemblies: a model-free baseline
# plus one run per model handed to it. Reads are fetched once, all assemblies run
# off them, and the FASTQ is deleted before any scoring happens -- so peak disk is
# one isolate's reads, not the cohort's.
#
# Scoring is QUAST -s and only the *_broken* row is used: TesserACT scaffolds, and
# the unsplit row would credit it for N-gaps it never closed (8.9x-43.9x
# inflation, docs/ESKAPEE_MODELS.md).
set -euo pipefail
org="$1"; safe="$2"; run="$3"; outdir="$4"; shift 4
root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
asm="$root/../TesserACT-main/tesseract-asm"
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate tesseract
threads="${THREADS:-32}"
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

need=0
for spec in base "$@"; do
    name="${spec##*:}"; [ "$spec" = base ] && name=base
    [ -s "$outdir/$name/contigs.fasta" ] || need=1
done

if [ "$need" = 1 ]; then
    rd="$outdir/reads"; mkdir -p "$rd"
    if [ ! -s "$rd/${run}_1.fastq" ]; then
        ok=0
        for try in 1 2 3; do
            if fasterq-dump --split-3 -e "$threads" -O "$rd" -t "$rd" "$run" >/dev/null 2>&1; then ok=1; break; fi
            sleep 15
        done
        [ "$ok" = 1 ] || { echo "  fasterq-dump failed for $run" >&2; rm -rf "$rd"; exit 1; }
    fi
    if [ ! -s "$rd/${run}_1.fastq" ] || [ ! -s "$rd/${run}_2.fastq" ]; then
        echo "  $run is not paired after dump" >&2; rm -rf "$rd"; exit 1
    fi

    for spec in base "$@"; do
        if [ "$spec" = base ]; then
            name=base; extra=()
        else
            model="${spec%%:*}"; name="${spec##*:}"
            extra=(--organism "$org" --model "$model")
        fi
        [ -s "$outdir/$name/contigs.fasta" ] && continue
        "$asm" -1 "$rd/${run}_1.fastq" -2 "$rd/${run}_2.fastq" -o "$outdir/$name" \
               -t "$threads" "${extra[@]}" > "$outdir/$name.log" 2>&1 \
            || echo "  assembly failed: $name" >&2
    done
    rm -rf "$rd"
fi

for spec in base "$@"; do
    name="${spec##*:}"; [ "$spec" = base ] && name=base
    [ -s "$outdir/$name/contigs.fasta" ] || continue
    [ -s "$outdir/quast_$name/report.tsv" ] && continue
    conda run -n quast quast.py -s -t "$threads" -R "$ref" --min-contig 500 \
        -o "$outdir/quast_$name" "$outdir/$name/contigs.fasta" >/dev/null 2>&1 || true
done
echo "  done $safe"
