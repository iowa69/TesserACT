#!/usr/bin/env bash
# finish_campaign.sh — the fixed protocol, applied identically to every organism.
#
# PROTOCOL (per organism, no deviation):
#   1. 30 held-out isolates, stratified clonal/singleton, depth band 30-200x,
#      real paired Illumina under the same BioSample or an exact strain match.
#   2. Leave-clone-out: the isolate's whole mash cluster withheld, chromosomes AND
#      plasmids, expanded through redundancy-group membership.
#   3. Five arms off ONE read download, reads deleted before the next isolate:
#        base        no model
#        small       ~137-genome panel  (the size the shipped models use)
#        clonal      full clonal panel  (this work)
#        clonal_f10  clonal + TESSERACT_MODEL_MIN_FRACTION=0.1
#        spades      vanilla SPAdes 4.3.0, raw reads, defaults
#   4. QUAST 5.3.0, -s, --min-contig 500, against the isolate's own closed genome.
#      Every TesserACT arm emits zero-N contigs.fasta, so all arms are contig-level.
#   5. Report medians AND paired win/loss, every metric, regressions included.
#
# Concurrency is fixed at 3 workers / 6 GB SPAdes: 4 workers at 10 GB exhausted a
# 62 GB box earlier and the harness killed the run.
set -euo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "$root"
source "$HOME/miniconda3/etc/profile.d/conda.sh"; conda activate tesseract

# Wait out any driver already running rather than competing with it.
while pgrep -f 'run_definitive_v[2].sh' >/dev/null; do sleep 60; done

# efaecium and abaumannii are revisited to repair SPAdes runs that died at -m 6.
# eval_isolate_v3.sh skips every arm already present, so only the missing SPAdes
# runs are redone -- at the cost of refetching those isolates' reads, which is the
# price of not keeping FASTQ around.
for org in efaecium abaumannii paeruginosa ecloacae ecoli kpneumoniae; do
    done_n=0
    for d in "eval/${org}_final"/*/; do
        [ -d "$d" ] || continue
        [ "$(ls -d "$d"/quast_* 2>/dev/null | wc -l)" -ge 5 ] && done_n=$((done_n+1))
    done
    if [ "$done_n" -ge 29 ]; then echo "==> $org already complete ($done_n/30)"; else
        echo "################ $org ################"
        # Two workers at 14 GB (28 GB of SPAdes) rather than three at 6 GB. The lower
# concurrency is the point: a SPAdes crash does not fail loudly, it quietly
# shrinks that column's cohort.
        TOTAL_THREADS=24 SPADES_MEM=14 ./run_definitive_v2.sh "$org" 2 2>&1 | tail -8 || true
    fi
    python3 collect_multi.py "eval/${org}_final" "eval/$org.final.tsv" \
        "eval/$org.holdout.tsv" > "eval/$org.final.summary.txt" 2>&1 || true
    cp -f "eval/$org.final.tsv" "../repo/docs/${org}_headtohead.tsv" 2>/dev/null || true
    echo "==> $org collected"
done
# Repair pass. An isolate can be dropped mid-run -- a killed worker, a refused
# download, or (once) a reads directory deleted out from under an in-flight fetch.
# eval_isolate_v3.sh skips every arm already on disk, so a second pass costs only the
# isolates that are actually missing, and without it the cohort silently ends up at 29
# or 28 with no record of which isolate went or why.
echo "==> repair pass"
for org in efaecium abaumannii paeruginosa ecloacae ecoli kpneumoniae; do
    n=0
    for d in "eval/${org}_final"/*/; do
        [ -d "$d" ] || continue
        [ "$(ls -d "$d"/quast_* 2>/dev/null | wc -l)" -ge 5 ] && n=$((n+1))
    done
    [ "$n" -ge 30 ] && continue
    echo "==> repairing $org ($n/30 complete)"
    TOTAL_THREADS=24 SPADES_MEM=14 ./run_definitive_v2.sh "$org" 2 2>&1 | tail -5 || true
    python3 collect_multi.py "eval/${org}_final" "eval/$org.final.tsv" \
        "eval/$org.holdout.tsv" > "eval/$org.final.summary.txt" 2>&1 || true
    cp -f "eval/$org.final.tsv" "../repo/docs/${org}_headtohead.tsv" 2>/dev/null || true
done
echo "CAMPAIGN COMPLETE: all seven organisms measured"
