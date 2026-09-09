#!/usr/bin/env bash
set -euo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "$root"
exec 9>"$root/.fetch3.lock"; flock -n 9 || { echo "fetch already running"; exit 1; }
for org in saureus efaecium abaumannii paeruginosa ecloacae ecoli kpneumoniae; do
    echo "######## fetching $org ########"
    ./fetch_reads_v2.sh "$org" 2 2>&1 | grep -E 'FETCHFAIL|libraries stored' || true
done
echo "FETCH COMPLETE"
df -h /media/iowa/u | tail -1
