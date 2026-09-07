#!/usr/bin/env bash
# sweep_fraction.sh ORG OUTROOT ISOLATE_LIST
#
# Tests the clonality claim under the gate that actually governs it.
#
# organism_join.cpp admits a candidate adjacency only if
#     e->support >= minFraction * genomes_carrying_both_markers   (kMinFractionChr 0.5)
# and separately if
#     support >= minPanel                                          (kMinPanelChr 20)
#
# A junction specific to one lineage is carried by that lineage's genomes and no
# others, so its share among genomes carrying both flanking markers is small by
# construction -- around 0.02 for a 16-genome lineage in an 800-genome panel. It
# is discarded before the vote. No amount of clonal depth changes that, which is
# why a clonal panel added 205,814 adjacencies to the model file and produced
# FEWER joins than the dereplicated one.
#
# So the panel and the threshold are crossed, 2x2, rather than changing both at
# once: "a lower threshold helps" and "clonal depth helps" are different claims.
# The absolute bar (minPanel 20) is deliberately left alone -- it is the real
# evidence requirement, and a lineage must still field 20 genomes to be heard.
set -euo pipefail
org="$1"; outroot="$2"; list="$3"
root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HOME/miniconda3/etc/profile.d/conda.sh"; conda activate tesseract
asm="$root/../TesserACT-main/tesseract-asm"
threads="${THREADS:-8}"

while IFS=$'\t' read -r safe run; do
    [ -n "${safe:-}" ] || continue
    out="$outroot/$safe"; mkdir -p "$out"
    ref=$(python3 - "$root/genomes/$org" "$safe" <<'PY'
import glob, os, sys
raw, safe = sys.argv[1], sys.argv[2]
for d in glob.glob(os.path.join(raw, "batch.*.d", "ncbi_dataset", "data", "GC*")):
    if os.path.basename(d).replace("_","").replace(".","v") == safe:
        f = glob.glob(os.path.join(d, "*.fna"))
        if f: print(f[0]); break
PY
)
    [ -n "$ref" ] || { echo "  no reference $safe" >&2; continue; }

    need=0
    for cfg in derep_f50 derep_f10 deep_f50 deep_f10 deep_f02 deep_f02p10 derep_f02p10; do
        [ -s "$out/$cfg/contigs.fasta" ] || need=1
    done
    if [ "$need" = 1 ]; then
        rd="$out/reads"; mkdir -p "$rd"
        if [ ! -s "$rd/${run}_1.fastq" ]; then
            fasterq-dump --split-3 -e "$threads" -O "$rd" -t "$rd" "$run" >/dev/null 2>&1 \
                || { echo "  dump failed $run" >&2; rm -rf "$rd"; continue; }
        fi
        # minFraction alone changes nothing: it admits candidates that then die on the
        # ABSOLUTE bar (kMinPanelChr 20), so joins stayed at 17 across 0.5 and 0.1 while
        # rejected_weak went 0 -> 6/8. And 0.1 is still far above what a lineage can reach:
        # 20 supporting genomes in an ~800-genome panel is a share of 0.025. So the low
        # arms drop the share below that AND relax the absolute bar, which is the only
        # configuration in which lineage-specific structure can be used at all.
        for spec in "derep_f50:derep_lio:0.5:20" "derep_f10:derep_lio:0.1:20" \
                    "deep_f50:deep_lio:0.5:20"  "deep_f10:deep_lio:0.1:20" \
                    "deep_f02:deep_lio:0.02:20" "deep_f02p10:deep_lio:0.02:10" \
                    "derep_f02p10:derep_lio:0.02:10"; do
            cfg=${spec%%:*}; rest=${spec#*:}; mdl=${rest%%:*}
            rest2=${rest#*:}; frac=${rest2%%:*}; panel=${rest2##*:}
            [ -s "$out/$cfg/contigs.fasta" ] && continue
            TESSERACT_MODEL_MIN_PANEL="$panel" \
            TESSERACT_MODEL_MIN_FRACTION="$frac" \
            "$asm" -1 "$rd/${run}_1.fastq" -2 "$rd/${run}_2.fastq" -o "$out/$cfg" \
                   -t "$threads" --organism "$org" --model "$root/models/$org/$mdl.tsm" \
                   > "$out/$cfg.log" 2>&1 || echo "  failed $safe/$cfg" >&2
        done
        rm -rf "$rd"
    fi
    for cfg in derep_f50 derep_f10 deep_f50 deep_f10 deep_f02 deep_f02p10 derep_f02p10; do
        [ -s "$out/$cfg/contigs.fasta" ] || continue
        [ -s "$out/quast_$cfg/report.tsv" ] && continue
        conda run -n quast quast.py -s -t "$threads" -R "$ref" --min-contig 500 \
            -o "$out/quast_$cfg" "$out/$cfg/contigs.fasta" >/dev/null 2>&1 || true
    done
    echo "  swept $safe"
done < "$list"
