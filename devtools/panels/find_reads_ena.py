#!/usr/bin/env python3
"""Second-pass read linkage via ENA, joining on strain name as well as BioSample.

The first pass (find_reads.py) pairs a genome with reads only when both sit under
the same BioSample, which covers 13-26% of complete genomes. Curated collections
frequently deposit the assembly and the raw reads under different samples, so the
pairing exists but the key does not.

ENA's portal API returns run_accession, sample_accession, strain and base_count
for an entire taxon in a single request, which makes a strain-name join cheap.

A wrong pairing is worse than a missing one -- it scores an assembly against
another isolate's genome -- so the join is deliberately conservative:

  * exact match on normalised strain (case and punctuation folded),
  * the normalised strain must identify exactly ONE assembly in our set, so
    lineage labels like "USA300" that name hundreds of isolates are rejected,
  * strains shorter than 4 characters, or on a stop-list of non-identifying
    labels, are rejected outright,
  * BioSample matches always win over strain matches.

Pairings are still checked downstream: an assembly scored against the wrong
reference shows up as a collapsed genome fraction, and collect_multi.py gates on
it. docs/KLEBSIELLA_PANEL.md records the same failure in the existing panel --
eight isolates with up to 53 kb not aligning to their own reference.
"""
import csv, re, sys, urllib.parse, urllib.request

ENA = "https://www.ebi.ac.uk/ena/portal/api/search"
STOP = {"wildtype", "wild type", "typestrain", "type strain", "reference", "na", "none",
        "unknown", "atcc", "nctc", "dsm", "clinical", "isolate", "control", "n/a", "missing"}


def norm(s):
    return re.sub(r"[^a-z0-9]", "", (s or "").lower())


def ena_runs(taxid, min_bases=150_000_000):
    q = (f'tax_eq({taxid}) AND library_layout="PAIRED" AND '
         f'instrument_platform="ILLUMINA" AND library_strategy="WGS"')
    data = urllib.parse.urlencode({
        "result": "read_run", "query": q,
        "fields": "run_accession,sample_accession,strain,isolate,base_count,"
                  "instrument_model,read_count",
        "format": "tsv", "limit": "0"}).encode()
    with urllib.request.urlopen(urllib.request.Request(ENA, data=data), timeout=900) as r:
        text = r.read().decode("utf-8", "replace")
    out = []
    for row in csv.DictReader(text.splitlines(), delimiter="\t"):
        try:
            b = int(row.get("base_count") or 0)
        except ValueError:
            continue
        if b < min_bases:
            continue
        out.append(row)
    return out


def main(metafile, taxid, existing, out):
    meta = list(csv.DictReader(open(metafile), delimiter="\t"))
    have = {}
    try:
        for r in csv.DictReader(open(existing), delimiter="\t"):
            have[r["genome"]] = r
    except FileNotFoundError:
        pass

    bs2acc = {r["biosample"]: r["accession"] for r in meta if r["biosample"]}

    # A normalised strain is usable only if it names exactly one ISOLATE. Ambiguity
    # is judged on distinct BioSamples rather than distinct accessions: GenBank and
    # RefSeq carry the same assembly under GCA_ and GCF_ accessions sharing one
    # BioSample, so an accession-level test marks every strain in a GenBank-inclusive
    # set as ambiguous and discards the lot. That is what reduced E. faecium -- the
    # organism most short of reads -- to zero strain matches while the RefSeq-only
    # organisms were unaffected.
    bystrain, seen, dup = {}, {}, set()
    for r in meta:
        s = norm(r.get("strain", ""))
        if len(s) < 4 or s in STOP:
            continue
        bs = r["biosample"]
        if s in seen and seen[s] != bs:
            dup.add(s)
        seen[s] = bs
        # Prefer the RefSeq copy of a GCA/GCF pair.
        if s not in bystrain or r["accession"].startswith("GCF_"):
            bystrain[s] = r["accession"]
    for s in dup:
        bystrain.pop(s, None)
    sys.stderr.write(f"{len(bystrain)} assemblies identified by a unique strain name "
                     f"({len(dup)} names rejected as ambiguous)\n")

    sys.stderr.write(f"querying ENA for taxid {taxid} ...\n")
    runs = ena_runs(taxid)
    sys.stderr.write(f"{len(runs)} Illumina paired WGS runs of sufficient depth\n")

    found, via = {}, {"biosample": 0, "strain": 0}
    for r in runs:
        acc = bs2acc.get(r["sample_accession"])
        how = "biosample"
        if not acc:
            s = norm(r.get("strain") or r.get("isolate") or "")
            acc = bystrain.get(s)
            how = "strain"
        if not acc:
            continue
        b = int(r["base_count"])
        cur = found.get(acc)
        # BioSample evidence outranks a strain match; then prefer the deeper run.
        if cur and (cur["how"] == "biosample" and how == "strain"):
            continue
        if cur and cur["how"] == how and cur["bases"] >= b:
            continue
        found[acc] = {"genome": acc, "biosample": r["sample_accession"],
                      "run": r["run_accession"], "bases": b,
                      "spots": int(r.get("read_count") or 0),
                      "readlen": (b // int(r["read_count"]) if r.get("read_count") and
                                  int(r["read_count"]) else 0),
                      "model": r.get("instrument_model", ""), "strategy": "WGS",
                      "how": how}
    for f in found.values():
        via[f["how"]] += 1

    merged = dict(have)
    added = 0
    for acc, f in found.items():
        if acc not in merged:
            merged[acc] = f
            added += 1
    cols = ["genome", "biosample", "run", "spots", "bases", "readlen", "model", "strategy"]
    with open(out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for g in sorted(merged):
            fh.write("\t".join(str(merged[g].get(c, "")) for c in cols) + "\n")
    sys.stderr.write(f"ENA matched {len(found)} genomes "
                     f"({via['biosample']} by BioSample, {via['strain']} by strain); "
                     f"{len(have)} known, {added} new -> {len(merged)} total in {out}\n")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
