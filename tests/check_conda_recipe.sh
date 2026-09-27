#!/usr/bin/env bash
# The conda recipe (conda/), checked against this source tree without conda or network.
#
# The recipe once pinned the previous release's tarball, installed files that tarball did not
# have, and tested a tool the package no longer ships and a help flag the help no longer
# lists; none of it failed anything here, because nothing ran the recipe. This does what
# conda-build would, with the compile step being the binary `make` already built:
#   1. meta.yaml's version is the version tesseract-asm reports;
#   2. every file build.sh installs exists in the tree (build.sh runs under set -e);
#   3. build.sh completes into a scratch PREFIX;
#   4. every meta.yaml test command exits 0 with only that PREFIX/bin (and the system) on PATH;
#   5. conda/run_test.sh -- the package test, including its 6 kb assembly -- passes;
#   6. in a git checkout that has the tag v<version>: the pinned sha256 equals
#      `git archive --prefix=TesserACT-<version>/ v<version> | gzip -n | sha256sum`, the method
#      that reproduces the GitHub tarball hash pinned for v1.2.5 exactly. Elsewhere it says so.
set -uo pipefail
while IFS= read -r v; do unset "$v"; done < <(compgen -e | grep '^TESSERACT_' || true)

ROOT=${SRCROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}   # the tree built and packaged
RECIPE=${RECIPE:-$ROOT/conda}
GITREPO=${GITREPO:-$ROOT}   # where to look for the release tag (read only)
TMP=$(mktemp -d "${TMPDIR:-/tmp}/tesseract-recipe.XXXXXXXX")
trap 'rm -rf "$TMP"' EXIT
pass=0; fail=0
ok()  { printf '  ok    %s\n' "$1"; pass=$((pass + 1)); }
bad() { printf '  FAIL  %s -- %s\n' "$1" "$2"; fail=$((fail + 1)); }

meta=$RECIPE/meta.yaml
version=$(sed -n 's/^{% set version = "\([^"]*\)" %}.*/\1/p' "$meta")
sha=$(sed -n 's/^{% set sha256 = "\([0-9a-f]*\)" %}.*/\1/p' "$meta")
binver=$("$ROOT/tesseract-asm" --version 2>/dev/null | sed -n '1s/^TesserACT //p')
echo "conda recipe $RECIPE: version=$version sha256=${sha:0:12}... binary=$binver"

# 1. version
[ -n "$version" ] && [ "$version" = "$binver" ] && ok "meta.yaml version $version is the binary's" \
    || bad "meta.yaml version" "recipe says '$version', tesseract-asm --version says '$binver'"

# 2. installed files exist
missing=""
while read -r src; do
    [ -e "$ROOT/$src" ] || missing="$missing $src"
done < <(sed -n 's/^install -m [0-7]* \([^ ]*\) .*/\1/p' "$RECIPE/build.sh")
[ -z "$missing" ] && ok "every file build.sh installs exists" || bad "build.sh installs" "missing:$missing"

# 3. build.sh into a scratch prefix (its `make` finds everything built)
P=$TMP/prefix; mkdir -p "$P"
if (cd "$ROOT" && PREFIX="$P" CPU_COUNT=1 bash "$RECIPE/build.sh") > "$TMP/build.log" 2>&1; then
    ok "build.sh completes ($(ls "$P/bin" | tr '\n' ' '))"
else
    bad "build.sh" "$(tail -1 "$TMP/build.log")"
fi

# The package test calls `python`; a system with only python3 still has to run it.
SHIM=$TMP/shim; mkdir -p "$SHIM"
# build_v3: decide on the PATH the package test really gets (a `python` elsewhere on the
# caller's PATH, e.g. a conda base, is not on TPATH and used to leave run_test.sh without one).
PATH="$P/bin:/usr/bin:/bin" command -v python >/dev/null 2>&1 || ln -s "$(command -v python3)" "$SHIM/python"
TPATH="$P/bin:$SHIM:/usr/bin:/bin"

# 4. test commands
n=0
while IFS= read -r cmd; do
    n=$((n + 1))
    if (cd "$TMP" && PATH="$TPATH" bash -c "$cmd") > "$TMP/cmd.log" 2>&1; then
        ok "test command: $cmd"
    else
        bad "test command: $cmd" "exit $? ($(tail -1 "$TMP/cmd.log"))"
    fi
done < <(sed -n '/^test:/,/^  requires:/p' "$meta" | sed -n 's/^ *- \(tesseract.*\)$/\1/p')
[ "$n" -gt 0 ] || bad "test commands" "none found in meta.yaml"

# 5. the package test
mkdir -p "$TMP/pkgtest"
if (cd "$TMP/pkgtest" && PATH="$TPATH" bash "$RECIPE/run_test.sh") > "$TMP/pkgtest.log" 2>&1; then
    ok "run_test.sh: $(grep -m1 '^assembled' "$TMP/pkgtest.log")"
else
    bad "run_test.sh" "$(tail -1 "$TMP/pkgtest.log")"
fi

# 6. tarball hash, where the tag is available
if git -C "$GITREPO" rev-parse -q --verify "refs/tags/v$version" >/dev/null 2>&1; then
    got=$(git -C "$GITREPO" archive --format=tar --prefix="TesserACT-$version/" "v$version" | gzip -n | sha256sum | cut -d' ' -f1)
    [ "$got" = "$sha" ] && ok "sha256 is the v$version tarball's" \
        || bad "sha256" "recipe pins ${sha:0:12}..., git archive v$version gives ${got:0:12}..."
else
    echo "  --    sha256 not checked: no git tag v$version here (compare against the GitHub archive)"
fi

printf '\n%d passed, %d failed\n' "$pass" "$fail"
[ "$fail" -eq 0 ]
