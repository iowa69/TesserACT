#!/usr/bin/env bash
# score_plasmids100.sh ORG — plasmid recovery + replicon-call accuracy for every
# completed isolate of a cohort. Read-only over the assemblies; safe to run while the
# campaign is still assembling, and re-running only adds rows for new isolates.
set -euo pipefail
org="$1"
root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "$root"
source "$HOME/miniconda3/etc/profile.d/conda.sh"; conda activate tesseract
out="eval100/$org.plasmids.tsv"
tmp="$out.tmp.$$"
printf "isolate\tarm\tn_plasmid\trecovered\twhole\tprecision\trecall\ttp_fp_fn\n" > "$tmp"
for d in eval100/"$org"/*/; do
    [ -d "$d" ] || continue
    iso=$(basename "$d")
    ref=$(python3 - "genomes/$org" "$iso" <<'PY'
import glob, os, sys
for d in glob.glob(os.path.join(sys.argv[1],"batch.*.d","ncbi_dataset","data","GC*")):
    if os.path.basename(d).replace("_","").replace(".","v")==sys.argv[2]:
        f=glob.glob(os.path.join(d,"*.fna"))
        if f: print(f[0]); break
PY
)
    [ -n "$ref" ] || continue
    for arm in base model careful aggressive spades; do
        f="$d/$arm/contigs.fasta"
        [ -s "$f" ] || continue
        python3 plasmid_recovery.py --ref "$ref" --asm "$f" --arm "$arm" \
            --isolate "$iso" --threads 2 >> "$tmp" 2>/dev/null || true
    done
done
mv -f "$tmp" "$out"
echo "$org: $(( $(wc -l < "$out") - 1 )) rows -> $out"
