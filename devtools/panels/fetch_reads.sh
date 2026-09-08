#!/usr/bin/env bash
# fetch_reads.sh ORG [JOBS] — fetch a cohort's reads to the SSD, gzipped, once.
#
# Reads are KEPT this time (the earlier design deleted them per isolate to bound disk).
# Storing them costs ~250 MB/isolate gzipped, ~175 GB for 700, against 891 GB free --
# and it means an assembler can be re-run, or a new arm added, without refetching 700
# libraries from SRA.
#
# Resumable and non-destructive: an isolate whose .gz pair already exists is skipped,
# and a partial fetch is written to a .tmp directory and only moved into place once
# both mates are complete. Nothing already stored is overwritten.
set -euo pipefail
org="$1"; JOBS="${2:-4}"
root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "$root"
source "$HOME/miniconda3/etc/profile.d/conda.sh"; conda activate tesseract
STORE="/media/iowa/u/tesseract_fastq/$org"
mkdir -p "$STORE"

tail -n +2 "eval100/$org.cohort.tsv" | awk -F'\t' '{print $1"\t"$3}' \
| xargs -P "$JOBS" -n 1 -d '\n' -I{} bash -c '
    set -- $(printf "%s" "{}")
    safe="$1"; run="$2"; store="'"$STORE"'"
    if [ -s "$store/${run}_1.fastq.gz" ] && [ -s "$store/${run}_2.fastq.gz" ]; then exit 0; fi
    tmp="$store/.tmp.$run"; rm -rf "$tmp"; mkdir -p "$tmp"
    ok=0
    for try in 1 2 3; do
        if fasterq-dump --split-3 -e 4 -O "$tmp" -t "$tmp" "$run" >/dev/null 2>&1; then ok=1; break; fi
        sleep 20
    done
    if [ "$ok" != 1 ] || [ ! -s "$tmp/${run}_1.fastq" ] || [ ! -s "$tmp/${run}_2.fastq" ]; then
        echo "FETCHFAIL $safe $run" >&2; rm -rf "$tmp"; exit 0
    fi
    pigz -p 4 -3 "$tmp/${run}_1.fastq" "$tmp/${run}_2.fastq"
    mv -n "$tmp/${run}_1.fastq.gz" "$store/${run}_1.fastq.gz"
    mv -n "$tmp/${run}_2.fastq.gz" "$store/${run}_2.fastq.gz"
    rm -rf "$tmp"
    echo "ok $safe $run"
  '
n=$(ls "$STORE"/*_1.fastq.gz 2>/dev/null | wc -l)
echo "$org: $n/$(( $(wc -l < "eval100/$org.cohort.tsv") - 1 )) libraries stored, $(du -sh "$STORE" | cut -f1)"
