#!/usr/bin/env bash
# conda-build entry point. Kept separate from the inline script in meta.yaml so
# the same steps work when building from a local checkout.
set -euo pipefail

make -j"${CPU_COUNT:-4}" CXX="${CXX:-g++}"

# One command on the PATH, `tesseract`; the engine, the model fetcher (with its checksum list)
# and the organism-detection sketch go to libexec/tesseract and share/tesseract. The layout is
# the Makefile's, so a source install and the package are the same.
make install PREFIX="${PREFIX}"
