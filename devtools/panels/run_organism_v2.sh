#!/usr/bin/env bash
# run_organism.sh ORG TAXID SOURCE MINLEN MAXLEN [NHOLD]
#
# The whole per-organism pipeline. Genomes are kept; FASTQ is not -- each held-out
# isolate's reads are fetched, assembled against every model, and deleted before
# the next isolate starts, so peak read storage is one isolate rather than a cohort.
set -euo pipefail
org="$1"; taxid="$2"; src="$3"; minlen="$4"; maxlen="$5"; nhold="${6:-30}"
root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$root"
# Background shells do not read .bashrc, so conda is activated explicitly.
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate tesseract
BIN=../TesserACT-main/tesseract-model
export THREADS="${THREADS:-32}"

step() { echo; echo "######## $org: $* ########"; }

step "metadata"
[ -s "meta/$org.tsv" ] || python3 fetch_meta.py "$taxid" "$src" "meta/$org.tsv"
tail -n +2 "meta/$org.tsv" | cut -f1 > "meta/$org.acc"

step "genomes"
[ -d "panel/$org/chr" ] && [ "$(ls "panel/$org/chr" | wc -l)" -gt 0 ] \
    || ./prep_genomes.sh "$org" "meta/$org.acc"
[ -s "panel/$org/plasmid_map.tsv" ] || python3 plasmid_map.py "genomes/$org" "panel/$org/plasmid_map.tsv"

step "reads available"
# Two passes. The first pairs reads to a genome through a shared BioSample; the
# second joins on strain name through ENA, which indexes runs the BioSample key
# misses because curated collections often deposit assembly and reads separately.
[ -s "meta/$org.reads.tsv" ] || python3 find_reads.py "meta/$org.tsv" "meta/$org.reads.tsv"
[ -s "meta/$org.reads2.tsv" ] || python3 find_reads_ena.py "meta/$org.tsv" "$taxid" \
    "meta/$org.reads.tsv" "meta/$org.reads2.tsv"
READS="meta/$org.reads2.tsv"

step "panels"
for spec in "1:panel_derep" "16:panel_clonal"; do
    cap=${spec%%:*}; tag=${spec##*:}
    [ -s "panel/$org/$tag.tsv" ] || python3 select_panel.py --org "$org" \
        --meta "meta/$org.tsv" --chrdir "panel/$org/chr" --out "panel/$org" \
        --cap "$cap" --tag "$tag" --reads "$READS" \
        --min-len "$minlen" --max-len "$maxlen"
done

step "hold-out"
mkdir -p eval
[ -s "eval/$org.holdout.tsv" ] || python3 choose_holdout.py --panel "panel/$org" \
    --tag panel_clonal --reads "$READS" --n "$nhold" \
    --out "eval/$org.holdout.tsv"

step "models"
mkdir -p "models/$org"
CL=$(cat "eval/$org.holdout.tsv.clusters"); AC=$(cat "eval/$org.holdout.tsv.accs")
for spec in "derep_lco:panel_derep:c" "clonal_lco:panel_clonal:c" \
            "derep_lio:panel_derep:a" "clonal_lio:panel_clonal:a"; do
    name=${spec%%:*}; rest=${spec#*:}; tag=${rest%%:*}; mode=${rest##*:}
    [ -s "models/$org/$name.tsm" ] && continue
    if [ "$mode" = c ]; then H=(--hold-clusters "$CL"); else H=(--hold-accs "$AC"); fi
    python3 build_model.py --org "$org" --panel "panel/$org" --tag "$tag" \
        --out "models/$org/$name.tsm" --binary "$BIN" "${H[@]}"
done

step "evaluation"
mkdir -p "eval/$org"
# Four isolates at a time on eight threads each, rather than one on thirty-two.
# Assembly scales sublinearly with threads, so the wider layout finishes the cohort
# in roughly a third of the wall clock at the same total CPU. Each worker still
# deletes its own reads before the next isolate, so peak read storage is four
# isolates, not the cohort.
JOBS="${JOBS:-4}"
PER=$(( ${TOTAL_THREADS:-32} / JOBS )); [ "$PER" -lt 1 ] && PER=1
export THREADS="$PER"
tail -n +2 "eval/$org.holdout.tsv" \
  | awk -F'\t' '{print $1"\t"$5}' \
  | xargs -P "$JOBS" -n 1 -d '\n' -I{} bash -c '
        set -- $(printf "%s" "{}")
        safe="$1"; run="$2"
        ./eval_isolate_v2.sh '"$org"' "$safe" "$run" "eval/'"$org"'/$safe" \
            "models/'"$org"'/derep_lco.tsm:derep_lco" \
            "models/'"$org"'/clonal_lco.tsm:clonal_lco" \
            "models/'"$org"'/derep_lio.tsm:derep_lio" \
            "models/'"$org"'/clonal_lio.tsm:clonal_lio" \
            || echo "  SKIPPED $safe"
    '

step "results"
python3 collect_multi.py "eval/$org" "eval/$org.results.tsv" "eval/$org.holdout.tsv" | tee "eval/$org.summary.txt"

step "shippable model (whole panel, nothing withheld)"
[ -s "models/$org/ship.tsm" ] || python3 build_model.py --org "$org" --panel "panel/$org" \
    --tag panel_clonal --out "models/$org/ship.tsm" --binary "$BIN"
ls -lh "models/$org/"
