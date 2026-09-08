#!/usr/bin/env bash
# conda-build entry point. Kept separate from the inline script in meta.yaml so
# the same steps work when building from a local checkout.
set -euo pipefail

make -j"${CPU_COUNT:-4}" CXX="${CXX:-g++}"

mkdir -p "${PREFIX}/bin"
install -m 0755 tesseract-asm "${PREFIX}/bin/tesseract-asm"
# tesseract-model is not packaged; see the Makefile comment on `all`.

# The Klebsiella runner ships too, as `tesseract-klebsiella`. Without it an installed package
# is the assembler and nothing else: the user still has to find the release page, download a
# 339 MB model by hand and work out the flags. The script finds the model, checks it and
# assembles, which is the difference between "installed" and "usable".
#
# It resolves the assembler as $(dirname $0)/tesseract-asm, which in an installed tree is the
# TesserACT sitting beside it -- so the installed copy uses the installed binary, not whatever
# happens to be on PATH.
install -m 0755 tesseract-klebsiella "${PREFIX}/bin/tesseract-klebsiella"

# The same two commands install.sh installs, for the same reason: tesseract-eskape is the
# entry point every organism preset is documented under, and tesseract-get-models is how its
# models arrive. A package without them installs cleanly and then has no command by the name
# the docs use.
#
# Not in the 1.2.5 package -- that tarball is already cut and bioconda PR #68728 is queued
# against its sha256, so this ships in 1.3 rather than invalidating an open PR.
install -m 0755 tesseract-eskape "${PREFIX}/bin/tesseract-eskape"
install -m 0755 tesseract-get-models "${PREFIX}/bin/tesseract-get-models"
# tesseract-get-models resolves the checksum list as $(dirname $0)/models.sha256, so it has to
# sit beside the script: without it every downloaded model fails verification and is deleted.
install -m 0644 models.sha256 "${PREFIX}/bin/models.sha256"
