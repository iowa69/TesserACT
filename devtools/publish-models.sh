#!/usr/bin/env bash
# Publish the organism models as the `models-v2` release assets.
#
#   devtools/publish-models.sh [MODEL_DIR]
#
# MODEL_DIR defaults to the checked copy kept outside git at
# ../../release/models-v2 -- the models are about 1.4 GB and belong in a release, not a history.
#
# The release is models-v2, not models-v1: the tesseract-get-models of 1.3.0 fetches models-v1
# and checks against a different checksum list, so it would download 1.4 GB and reject it all.
#
# Needs gh already authenticated (`gh auth status`). It never reads a token from the
# environment or an argument: the credential stays wherever gh keeps it.
#
# Safe to re-run: it creates the release only if missing, and --clobber replaces an asset
# rather than failing on a name that already exists.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="iowa69/TesserACT"
tag="models-v2"
dir="${1:-$here/../../release/models-v2}"
manifest="$here/../models.sha256"

command -v gh >/dev/null 2>&1 || { echo "gh is not installed" >&2; exit 1; }
gh auth status >/dev/null 2>&1 || { echo "gh is not authenticated -- run: gh auth login" >&2; exit 1; }
[ -d "$dir" ] || { echo "no such directory: $dir" >&2; exit 1; }
[ -r "$manifest" ] || { echo "cannot read $manifest" >&2; exit 1; }

# Verified before anything is uploaded, not after. tesseract-get-models checks every download
# against this same list and deletes what does not match, so publishing a file that fails here
# produces a release nobody can install from -- and the failure surfaces on the user's machine,
# not ours.
echo "==> verifying $dir against models.sha256"
( cd "$dir" && sha256sum -c "$manifest" ) || { echo "checksum mismatch -- not publishing" >&2; exit 1; }

if gh release view "$tag" -R "$repo" >/dev/null 2>&1; then
    echo "==> release $tag exists; uploading into it"
else
    echo "==> creating release $tag"
    gh release create "$tag" -R "$repo" --target main \
        --title "Organism models (models-v2)" \
        --notes "The seven ESKAPEE organism models that tesseract-get-models fetches and tesseract-eskape loads by preset. Each carries both the chromosome panel and the plasmid database for its organism. Install with: tesseract-get-models -- they land in ~/.tesseract/models and are verified against models.sha256. They are optional: TesserACT assembles without them, and without one tesseract-eskape says so and uses the assembler defaults. A model buys contiguity and costs misassemblies; see models/README.md for the measured cost. Giovanni Lorenzin, IOWA-BioTech"
fi

while read -r _sha name; do
    [ -n "${_sha:-}" ] || continue
    case "$_sha" in \#*) continue ;; esac
    echo "==> uploading $name"
    gh release upload "$tag" -R "$repo" "$dir/$name" --clobber
done < "$manifest"

echo
echo "==> done. Check it end to end with:"
echo "    rm -rf /tmp/mtest && tesseract-get-models --dir /tmp/mtest"
