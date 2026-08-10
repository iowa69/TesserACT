#!/usr/bin/env bash
# unlock.sh <asset.tsm.zst.gpg> [output.tsm]
#
# Decrypts and decompresses a released TesserACT model.
#
# Release models are distributed encrypted. The asset itself is public; the key is not, and
# is issued separately. Without the key the file is a 1.5 GB opaque blob -- that is
# deliberate, and it is the only thing restricting the models, since a GitHub release asset
# on a public repository can be downloaded by anyone with the URL.
#
#   ./unlock.sh tessera-klebsiella-default-v1.2.0.tsm.zst.gpg kleb.tsm
#   tessera -1 R1.fq.gz -2 R2.fq.gz -o out/ --organism klebsiella --model kleb.tsm
#
# The passphrase is read from $TESSERA_MODEL_KEY (a file path) if set, otherwise gpg prompts.
# Passing it on the command line would put it in your shell history, so there is no flag for
# it.
set -euo pipefail

SRC=${1:?usage: unlock.sh <asset.tsm.zst.gpg> [output.tsm]}
DST=${2:-$(basename "$SRC" .zst.gpg)}

command -v gpg  >/dev/null || { echo "gpg is required" >&2; exit 1; }
command -v zstd >/dev/null || { echo "zstd is required" >&2; exit 1; }
[ -s "$SRC" ] || { echo "no such asset: $SRC" >&2; exit 1; }

# Verify against SHA256SUMS when it is sitting next to the asset. A truncated download is the
# most common failure and it is worth catching before spending the decryption.
SUMS=$(dirname "$SRC")/SHA256SUMS
if [ -s "$SUMS" ]; then
  if (cd "$(dirname "$SRC")" && sha256sum --check --ignore-missing --status SHA256SUMS); then
    echo "checksum ok"
  else
    echo "CHECKSUM FAILED for $SRC -- re-download before going further" >&2
    exit 1
  fi
fi

gpg_args=(--batch --yes --decrypt)
if [ -n "${TESSERA_MODEL_KEY:-}" ]; then
  [ -s "$TESSERA_MODEL_KEY" ] || { echo "TESSERA_MODEL_KEY is set but empty: $TESSERA_MODEL_KEY" >&2; exit 1; }
  gpg_args+=(--passphrase-file "$TESSERA_MODEL_KEY")
else
  gpg_args=(--decrypt)   # let gpg prompt on the tty
fi

# Straight to a temporary file and moved into place on success: a failed decryption otherwise
# leaves a partial .tsm that the assembler would try to load and reject with a confusing error
# about the file not being a model.
tmp=$DST.part
trap 'rm -f "$tmp"' EXIT
gpg "${gpg_args[@]}" "$SRC" | zstd -d -q -o "$tmp" -f
mv "$tmp" "$DST"
trap - EXIT

echo "wrote $DST ($(du -h "$DST" | cut -f1))"
echo "check it loads:  tessera --model $DST --organism klebsiella -1 <reads> -2 <reads> -o out/"
