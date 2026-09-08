#!/usr/bin/env bash
# campaign100.sh — the 1242-isolate, five-arm campaign.
#
# Runs alongside the fetch: an isolate whose reads are not yet stored is skipped with a
# message and picked up on a later pass, so assembly does not wait for the whole fetch.
# Two passes over every organism, then a repair pass, because an isolate skipped for
# missing reads on pass 1 is usually available by pass 3.
#
# Nothing is overwritten: every arm is skipped if its contigs.fasta already exists.
set -euo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "$root"
exec 9>"$root/.campaign100v2.lock"; flock -n 9 || { echo "campaign100 already running"; exit 1; }
ORGS="saureus efaecium abaumannii paeruginosa ecloacae ecoli kpneumoniae"
for pass in 1 2 3 4; do
    echo "================ PASS $pass ================"
    for org in $ORGS; do
        [ -s "models100/bundled/$org.tsm" ] || { echo "== $org: no model yet, skipping"; continue; }
        avail=$(ls /media/iowa/u/tesseract_fastq/$org/*_1.fastq.gz 2>/dev/null | wc -l)
        [ "$avail" -gt 0 ] || { echo "== $org: no reads stored yet"; continue; }
        done_n=0
        for d in "eval100/$org"/*/; do
            [ -d "$d" ] || continue
            [ "$(ls -d "$d"/quast_* 2>/dev/null | wc -l)" -ge 5 ] && done_n=$((done_n+1))
        done
        want=$(( $(wc -l < "eval100/$org.cohort.tsv") - 1 ))
        [ "$done_n" -ge "$want" ] && { echo "== $org complete ($done_n/$want)"; continue; }
        echo "== $org pass $pass: $done_n/$want done, $avail libraries stored"
        TOTAL_THREADS=30 SPADES_MEM=8 ./run_asm100.sh "$org" 5 2>&1 | tail -3 || true
    done
done
echo "CAMPAIGN100 PASSES COMPLETE"
