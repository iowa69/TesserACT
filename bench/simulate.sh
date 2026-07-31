#!/usr/bin/env bash
# Fetch the benchmark reference genomes and simulate Illumina reads from them.
#
# Requires: curl, gzip, and ART (art_illumina). pigz is used when available.
#   conda create -n benchtools -c bioconda -c conda-forge art quast
#
# Genomes are simulated one at a time and compressed immediately, so peak disk
# stays near 1.2 GB rather than the ~4 GB the uncompressed set would need.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REFS="$HERE/refs"
READS="$HERE/reads"
ART="${ART:-art_illumina}"

COV=${COV:-100}
LEN=${LEN:-150}
FRAG=${FRAG:-350}
FRAGSD=${FRAGSD:-50}

command -v "$ART" >/dev/null || { echo "error: art_illumina not found; set ART=/path/to/art_illumina" >&2; exit 1; }
mkdir -p "$REFS" "$READS"

# RefSeq assembly paths for the four panel genomes.
declare -A GENOMES=(
  [ecoli]="GCF/000/005/845/GCF_000005845.2_ASM584v2"
  [saureus]="GCF/000/013/425/GCF_000013425.1_ASM1342v1"
  [kpneu]="GCF/000/240/185/GCF_000240185.1_ASM24018v2"
  [lmono]="GCF/000/196/035/GCF_000196035.1_ASM19603v1"
)

for g in "${!GENOMES[@]}"; do
  if [[ -s "$REFS/${g}.fna.gz" ]]; then continue; fi
  path="${GENOMES[$g]}"
  base=$(basename "$path")
  echo "[fetch] $g"
  curl -sfL --retry 3 -o "$REFS/${g}.fna.gz" \
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/${path}/${base}_genomic.fna.gz"
  gzip -t "$REFS/${g}.fna.gz"
done

for g in ecoli saureus kpneu lmono; do
  if [[ -s "$READS/${g}_1.fq.gz" && -s "$READS/${g}_2.fq.gz" ]]; then
    echo "[skip] $g already simulated"
    continue
  fi
  echo "[sim] $g  ${COV}x  ${LEN}bp PE  frag=${FRAG}+/-${FRAGSD}"
  zcat "$REFS/${g}.fna.gz" > "$READS/${g}.fna"
  "$ART" -ss HS25 -i "$READS/${g}.fna" -p -l "$LEN" -f "$COV" \
         -m "$FRAG" -s "$FRAGSD" -rs 42 -na -o "$READS/${g}_" > "$READS/${g}.artlog" 2>&1
  rm -f "$READS/${g}.fna"
  if command -v pigz >/dev/null; then
    pigz -f "$READS/${g}_1.fq" "$READS/${g}_2.fq"
  else
    gzip -f "$READS/${g}_1.fq" "$READS/${g}_2.fq"
  fi
  echo "[done] $g"
done

ls -lh "$READS"/*.fq.gz
