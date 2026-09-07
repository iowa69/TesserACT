#!/usr/bin/env bash
# build_ship_models.sh PANEL_TAG OUTDIR
#
# Builds the models users actually download: the whole panel, nothing withheld.
#
# The evaluation models withhold clusters so a measured gain cannot come from having
# seen the answer. A shipped model has no test isolate to protect and should see
# everything -- models/README.md notes exactly this about the released Klebsiella
# models: "these models have never seen roughly a fifth of public K. pneumoniae
# diversity ... for production use a model built over the whole panel would be
# marginally stronger; these are released because they are the ones every published
# figure was measured with."
#
# That is a defensible choice for a paper and the wrong one for a release. Here the
# two are separated: measured numbers come from the withheld builds, and the download
# is the whole-panel build. The distinction is stated in the release notes rather
# than left for a reader to infer.
set -euo pipefail
tag="${1:?usage: build_ship_models.sh PANEL_TAG OUTDIR}"; out="${2:?}"
root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "$root"
source "$HOME/miniconda3/etc/profile.d/conda.sh"; conda activate tesseract
mkdir -p "$out"
for org in saureus efaecium abaumannii paeruginosa ecloacae ecoli kpneumoniae; do
    t="$tag"
    # E. coli's clonal panel is capped; every other organism uses the plain tag.
    [ "$org" = ecoli ] && [ "$tag" = panel_clonal ] && t=panel_clonal_capped
    [ -s "panel/$org/$t.tsv" ] || { echo "skip $org: no $t.tsv" >&2; continue; }
    [ -s "$out/$org.tsm" ] && { echo "have $org.tsm"; continue; }
    echo "==> $org from $t"
    python3 build_model.py --org "$org" --panel "panel/$org" --tag "$t" \
        --out "$out/$org.tsm" --binary ../TesserACT-main/tesseract-model 2>&1 \
        | grep -E 'panel |chromosome:|plasmid: |layout:|markers,'
done
echo
echo "==> manifest"
( cd "$out" && sha256sum *.tsm > models.sha256 && cat models.sha256 )
du -ch "$out"/*.tsm | tail -1
