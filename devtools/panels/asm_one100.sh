#!/usr/bin/env bash
# asm_one100.sh ORG SAFE RUN — one isolate, five arms, scored against its own closed genome.
set -euo pipefail
org="$1"; safe="$2"; run="$3"
root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "$root"
source "$HOME/miniconda3/etc/profile.d/conda.sh"; conda activate tesseract
# The model arms go through --organism, which after the lockdown resolves only from
# TESSERACT_MODEL_DIR (or ~/.tesseract/models). Set here so the script is
# self-contained and exercises exactly the interface a user has.
export TESSERACT_MODEL_DIR="${TESSERACT_MODEL_DIR:-$root/models100/bundled}"
asm="$root/../repo/tesseract-asm"
store="/media/iowa/u/tesseract_fastq/$org"
out="eval100/$org/$safe"; mkdir -p "$out"
R1="$store/${run}_1.fastq.gz"; R2="$store/${run}_2.fastq.gz"
[ -s "$R1" ] && [ -s "$R2" ] || { echo "  no stored reads for $safe ($run)" >&2; exit 1; }
th="${THREADS:-8}"

ref=$(python3 - "$root/genomes/$org" "$safe" <<'PY'
import glob, os, sys
raw, safe = sys.argv[1], sys.argv[2]
for d in glob.glob(os.path.join(raw, "batch.*.d", "ncbi_dataset", "data", "GC*")):
    if os.path.basename(d).replace("_","").replace(".","v") == safe:
        f = glob.glob(os.path.join(d, "*.fna"))
        if f: print(f[0]); break
PY
)
[ -n "$ref" ] && [ -s "$ref" ] || { echo "  no closed reference for $safe" >&2; exit 1; }

run_arm(){ # name, extra args...
  local n="$1"; shift
  [ -s "$out/$n/contigs.fasta" ] && return 0
  "$asm" -1 "$R1" -2 "$R2" -o "$out/$n" -t "$th" "$@" > "$out/$n.log" 2>&1 \
    || echo "  arm failed: $safe/$n" >&2
}
run_arm base
run_arm model      --organism "$org"
run_arm careful    --organism "$org" --mode careful
run_arm aggressive --organism "$org" --mode aggressive

if [ ! -s "$out/spades/contigs.fasta" ]; then
  for mem in "${SPADES_MEM:-14}" $(( ${SPADES_MEM:-14} * 2 )); do
    rm -rf "$out/spades_work"
    conda run -n spades spades.py -1 "$R1" -2 "$R2" -o "$out/spades_work" \
      -t "$th" -m "$mem" > "$out/spades.log" 2>&1 && break
  done
  if [ -s "$out/spades_work/contigs.fasta" ]; then
    mkdir -p "$out/spades"; cp "$out/spades_work/contigs.fasta" "$out/spades/contigs.fasta"
    [ -s "$out/spades_work/scaffolds.fasta" ] && cp "$out/spades_work/scaffolds.fasta" "$out/spades/scaffolds.fasta"
  fi
  rm -rf "$out/spades_work"
fi

for n in base model careful aggressive spades; do
  [ -s "$out/$n/contigs.fasta" ] || continue
  [ -s "$out/quast_$n/report.tsv" ] && continue
  conda run -n quast quast.py -s -t "$th" -R "$ref" --min-contig 500 \
    -o "$out/quast_$n" "$out/$n/contigs.fasta" >/dev/null 2>&1 || true
done
echo "  done $safe"
