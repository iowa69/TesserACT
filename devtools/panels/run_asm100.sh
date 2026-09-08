#!/usr/bin/env bash
# run_asm100.sh ORG [JOBS] — five arms per isolate off the stored FASTQ.
#
# Reads come from the SSD store and are NEVER deleted or modified here; this script
# only reads them. Every arm writes to its own directory and an arm already present
# is skipped, so a re-run costs only what is missing and nothing already computed is
# overwritten.
#
# The model arms go through --organism with TESSERACT_MODEL_DIR pointed at the
# bundled directory: that is the only interface users have after the lockdown, so the
# benchmark exercises the shipped path rather than a developer one.
set -euo pipefail
org="$1"; JOBS="${2:-3}"
root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "$root"
source "$HOME/miniconda3/etc/profile.d/conda.sh"; conda activate tesseract
export TESSERACT_MODEL_DIR="$root/models100/bundled"
PER=$(( ${TOTAL_THREADS:-24} / JOBS )); [ "$PER" -lt 1 ] && PER=1
export THREADS="$PER" SPADES_MEM="${SPADES_MEM:-14}"
mkdir -p "eval100/$org"

tail -n +2 "eval100/$org.cohort.tsv" | awk -F'\t' '{print $1"\t"$3}' \
| xargs -P "$JOBS" -n 1 -d '\n' -I{} bash -c '
    set -- $(printf "%s" "{}")
    safe="$1"; run="$2"
    ./asm_one100.sh '"$org"' "$safe" "$run" || echo "  SKIPPED $safe"
  '
echo "$org: $(ls -d eval100/$org/*/quast_* 2>/dev/null | wc -l) scored"
