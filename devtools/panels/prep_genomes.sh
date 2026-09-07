#!/usr/bin/env bash
# prep_genomes.sh ORG ACC_LIST
#
# Downloads the listed assemblies and splits every one into a chromosome FASTA
# and its plasmid records.
#
# Panel file names carry the accession as their leading token, and
# `accessionOf()` in src/model_main.cpp cuts that token at the FIRST '_' or '.'.
# So GCF_001518735.1_chr.fasta resolves to the accession "GCF" -- every panel
# file collapsing onto one genome, which the build reports as success. This is
# the defect ESKAPEE_MODELS.md records; it was fixed by renaming, not in code,
# so the renaming is mandatory: GCF_001518735.1 -> GCF001518735v1.
set -euo pipefail
org="$1"; acclist="$2"
root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
raw="$root/genomes/$org"; out="$root/panel/$org"
mkdir -p "$raw" "$out/chr"
: > "$out/plasmids.fna"
: > "$out/split.log"

n=$(wc -l < "$acclist")
echo "==> $org: downloading $n assemblies"
# Batched: one datasets call per 500 accessions keeps a failure from costing the
# whole set and keeps the zip a sane size.
split -l 500 -d "$acclist" "$raw/batch."
# Globbed on exactly two digits. Plain batch.* also matches the batch.NN.zip
# archives and the batch.NN.d extraction directories this loop creates, so on a
# re-run a directory was handed to --inputfile and datasets answered with its
# usage text -- the script worked once and failed on resume, which is when it
# actually matters, because resume is the whole point of an interrupted fetch.
for b in "$raw"/batch.[0-9][0-9]; do
    z="$b.zip"
    [ -s "$z" ] || datasets download genome accession --inputfile "$b" \
        --include genome --no-progressbar --filename "$z"
    rm -rf "$b.d" && mkdir -p "$b.d" && unzip -q -o "$z" -d "$b.d"
done

echo "==> $org: splitting chromosome / plasmid"
python3 - "$raw" "$out" <<'PY'
import glob, os, re, sys
raw, out = sys.argv[1], sys.argv[2]
plas = open(os.path.join(out, "plasmids.fna"), "w")
log  = open(os.path.join(out, "split.log"), "w")
nchr = npls = ngen = 0
for d in sorted(glob.glob(os.path.join(raw, "batch.*.d", "ncbi_dataset", "data", "GC*"))):
    acc = os.path.basename(d)
    safe = acc.replace("_", "").replace(".", "v")   # GCF_001518735.1 -> GCF001518735v1
    fnas = glob.glob(os.path.join(d, "*.fna"))
    if not fnas:
        log.write(f"{acc}\tNO_FASTA\n"); continue
    recs, name, seq = [], None, []
    for fn in fnas:
        with open(fn) as fh:
            for line in fh:
                if line.startswith(">"):
                    if name is not None: recs.append((name, "".join(seq)))
                    name, seq = line[1:].rstrip("\n"), []
                else:
                    seq.append(line.strip())
    if name is not None: recs.append((name, "".join(seq)))
    # A record is a plasmid if its definition line says so. Everything else on a
    # complete assembly is chromosomal; unplaced material does not exist at this
    # assembly level.
    chrs = [(n, s) for n, s in recs if "plasmid" not in n.lower()]
    pls  = [(n, s) for n, s in recs if "plasmid" in n.lower()]
    if not chrs:
        log.write(f"{acc}\tNO_CHROMOSOME\n"); continue
    with open(os.path.join(out, "chr", f"{safe}_chr.fasta"), "w") as fh:
        for n, s in chrs:
            fh.write(f">{n}\n")
            for i in range(0, len(s), 80): fh.write(s[i:i+80] + "\n")
    nchr += len(chrs); ngen += 1
    for n, s in pls:
        # The record name must be the bare accession: --exclude-plasmids matches
        # the first whitespace-delimited token of the header, exactly.
        tok = n.split()[0]
        plas.write(f">{tok}\n")
        for i in range(0, len(s), 80): plas.write(s[i:i+80] + "\n")
        npls += 1
    log.write(f"{acc}\t{safe}\tchr={len(chrs)}\tplasmid={len(pls)}\n")
plas.close(); log.close()
print(f"  {ngen} genomes, {nchr} chromosome records, {npls} plasmid records")
PY
echo "==> $org: done"
