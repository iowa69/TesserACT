#!/usr/bin/env bash
# run_definitive.sh ORG [JOBS]
#
# The shipping run. Four arms per held-out isolate, one read download:
#
#   base        TesserACT, no model            -- what the assembler does alone
#   small       TesserACT + ~137-genome panel  -- the panel size the shipped models use
#   clonal      TesserACT + full clonal panel  -- the panel this work builds
#   spades      vanilla SPAdes 4.3.0           -- the external reference point
#
# "small" is here because it is the open question that decides what ships: the shipped
# ESKAPEE models were built on ~137 dereplicated genomes and docs/REAL_WORLD_RESULTS.md
# reports misassemblies FALLING with them, while the 879/1268-genome panels measured here
# show them rising 4 -> 30. Either panel breadth costs structural accuracy, or that
# documented result does not reproduce. Running both panels against the same isolates and
# the same reads is the only way to tell, and it has to be settled before publishing
# models built on the larger panels.
#
# SPAdes is vanilla per README.md: raw untrimmed reads, default settings, its own error
# correction, automatic k. Its contigs.fasta is unscaffolded; TesserACT's contigs.fasta is
# now zero-N, so both are contig-level and directly comparable with no column to choose.
set -euo pipefail
org="$1"; JOBS="${2:-4}"
root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "$root"
source "$HOME/miniconda3/etc/profile.d/conda.sh"; conda activate tesseract
PER=$(( ${TOTAL_THREADS:-32} / JOBS )); [ "$PER" -lt 1 ] && PER=1
export THREADS="$PER" SPADES_MEM="${SPADES_MEM:-12}"
mkdir -p "eval/${org}_final"

tail -n +2 "eval/$org.holdout.tsv" | awk -F'\t' '{print $1"\t"$5}' \
| xargs -P "$JOBS" -n 1 -d '\n' -I{} bash -c '
    set -- $(printf "%s" "{}")
    ./eval_isolate_v2.sh '"$org"' "$1" "$2" "eval/'"$org"'_final/$1" \
        "models/'"$org"'_v2/small_lco.tsm:small" \
        "models/'"$org"'_v2/clonal_lco.tsm:clonal" || echo "  SKIPPED $1"
  '
python3 collect_multi.py "eval/${org}_final" "eval/$org.final.tsv" "eval/$org.holdout.tsv" \
  | tee "eval/$org.final.summary.txt"
