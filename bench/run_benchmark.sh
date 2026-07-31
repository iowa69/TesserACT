#!/usr/bin/env bash
# Assemble the benchmark panel with tessera and SPAdes, then score both with
# QUAST against the reference.
#
# Requires: spades.py and quast.py on PATH (or set SPADES=/QUAST= to point at
# them), and reads produced by simulate.sh.
#
# One genome at a time, deleting SPAdes' per-k working directories as soon as
# each run finishes, so peak disk stays a few GB rather than tens.
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TESSERA="${TESSERA:-$HERE/../tessera}"
SPADES="${SPADES:-spades.py}"
QUAST="${QUAST:-quast.py}"
THREADS="${THREADS:-16}"

READS="$HERE/reads"
REFS="$HERE/refs"
RESULTS="$HERE/results"
WORK="$HERE/work"
mkdir -p "$RESULTS" "$WORK"

[[ -x "$TESSERA" ]] || { echo "error: tessera binary not found at $TESSERA (run make first)" >&2; exit 1; }

for g in "${@:-ecoli saureus kpneu lmono}"; do
  R1="$READS/${g}_1.fq.gz"
  R2="$READS/${g}_2.fq.gz"
  [[ -s "$R1" && -s "$R2" ]] || { echo "error: reads for $g missing; run simulate.sh" >&2; continue; }

  REF="$RESULTS/${g}_ref.fna"
  zcat "$REFS/${g}.fna.gz" > "$REF"
  echo "=================== $g ==================="

  if [[ ! -s "$RESULTS/${g}_tessera.fasta" ]]; then
    echo "[tessera] $g"
    /usr/bin/time -f "TESSERA_TIME %e s  MEM %M KB" \
      "$TESSERA" -1 "$R1" -2 "$R2" -o "$WORK/ts_$g" -t "$THREADS" -q \
      > "$RESULTS/${g}_tessera.log" 2>&1
    tail -1 "$RESULTS/${g}_tessera.log"
    cp "$WORK/ts_$g/contigs.fasta" "$RESULTS/${g}_tessera.fasta" 2>/dev/null
    rm -rf "$WORK/ts_$g"
  fi

  if [[ ! -s "$RESULTS/${g}_spades.fasta" ]] && command -v "$SPADES" >/dev/null; then
    echo "[spades] $g"
    /usr/bin/time -f "SPADES_TIME %e s  MEM %M KB" \
      "$SPADES" --isolate -1 "$R1" -2 "$R2" -o "$WORK/sp_$g" -t "$THREADS" -m 24 \
      > "$RESULTS/${g}_spades.log" 2>&1
    grep -E "SPADES_TIME" "$RESULTS/${g}_spades.log" | tail -1
    cp "$WORK/sp_$g/contigs.fasta"   "$RESULTS/${g}_spades.fasta"      2>/dev/null
    cp "$WORK/sp_$g/scaffolds.fasta" "$RESULTS/${g}_spades_scaf.fasta" 2>/dev/null
    rm -rf "$WORK/sp_$g"
  fi

  echo "[quast] $g"
  ASSEMBLIES=("$RESULTS/${g}_tessera.fasta")
  LABELS="tessera"
  for extra in spades spades_scaf; do
    if [[ -s "$RESULTS/${g}_${extra}.fasta" ]]; then
      ASSEMBLIES+=("$RESULTS/${g}_${extra}.fasta")
      LABELS="$LABELS,$extra"
    fi
  done
  "$QUAST" -r "$REF" -o "$RESULTS/quast_$g" -t 8 --silent -l "$LABELS" \
    "${ASSEMBLIES[@]}" > /dev/null 2>&1
  grep -E "^(# contigs +|Largest contig|NGA50|Genome fraction|# misassemblies|# mismatches|# indels|Duplication)" \
    "$RESULTS/quast_$g/report.txt"
  rm -f "$REF"
done

rmdir "$WORK" 2>/dev/null || true
echo "=================== done ==================="
