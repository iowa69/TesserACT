#!/usr/bin/env python3
"""Find Illumina runs deposited under the same BioSample as a closed genome.

The point of the pairing is that reads and reference come from the same DNA, so
completeness and accuracy are measured rather than inferred -- the design of
docs/CLOSED_REFERENCE_BENCHMARK.md. It replaces the simulated-read evaluation
the ESKAPEE models were scored with, which docs/ESKAPEE_MODELS.md already flags
as an upper bound: wgsim gives uniform coverage, no adapters, no GC bias, flat
0.5% error and no indels.
"""
import csv, io, re, sys, time, urllib.parse, urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def fetch(url, tries=5):
    for a in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return r.read().decode("utf-8", "replace")
        except Exception:
            if a == tries - 1: raise
            time.sleep(2 * (a + 1))


def post(url, data, tries=5):
    for a in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, data=data), timeout=180) as r:
                return r.read().decode("utf-8", "replace")
        except Exception:
            if a == tries - 1: raise
            time.sleep(2 * (a + 1))


def runinfo(biosamples):
    """One esearch per batch of BioSamples, then runinfo CSV for the hits."""
    # BSPL is the SRA BioSample field. "[BioSample]" is not a valid SRA field and
    # esearch answers PhraseNotFound with Count 0 -- a silent empty result, not an
    # error, which is why the first pass reported that no genome had reads.
    term = " OR ".join(f"{b}[BSPL]" for b in biosamples)
    data = urllib.parse.urlencode({"db": "sra", "term": term, "retmax": "5000",
                                   "usehistory": "y"}).encode()
    xml = post(f"{EUTILS}/esearch.fcgi", data)
    if "<WebEnv>" not in xml:
        return []
    env = xml.split("<WebEnv>")[1].split("</WebEnv>")[0]
    key = xml.split("<QueryKey>")[1].split("</QueryKey>")[0]
    q2 = urllib.parse.urlencode({"db": "sra", "WebEnv": env, "query_key": key,
                                 "rettype": "runinfo", "retmode": "csv"})
    csvtext = fetch(f"{EUTILS}/efetch.fcgi?{q2}")
    return parse_runinfo(csvtext)


PLATFORMS = {"ILLUMINA", "OXFORD_NANOPORE", "PACBIO_SMRT", "ION_TORRENT",
             "BGISEQ", "CAPILLARY", "ELEMENT", "ULTIMA", "DNBSEQ", "COMPLETE_GENOMICS"}


def parse_runinfo(text):
    """Parse runinfo CSV positionally, locating fields by their value shape.

    efetch returns runinfo for a WebEnv query with NO header line, and the column
    offsets are not stable across rows. Handing it to csv.DictReader silently
    consumes the first run as the header and every lookup then returns None --
    which is why this reported that no S. aureus genome had reads. Fields are
    found by pattern instead, so a shifted column cannot produce a silent empty
    result again.
    """
    out = []
    for row in csv.reader(io.StringIO(text)):
        if len(row) < 20:
            continue
        if not re.match(r"^[DES]RR\d+$", row[0] or ""):
            continue
        bs = next((v for v in row if re.match(r"^SAM[NED]", v or "")), "")
        plat = next((v for v in row if v in PLATFORMS), "")
        layout = "PAIRED" if "PAIRED" in row else ("SINGLE" if "SINGLE" in row else "")
        strat = next((v for v in row if v in
                      ("WGS", "OTHER", "WGA", "CLONE", "AMPLICON", "RNA-Seq")), "")
        try:
            spots, bases = int(row[3] or 0), int(row[4] or 0)
        except ValueError:
            continue
        model = ""
        try:
            model = row[row.index(plat) + 1] if plat in row else ""
        except (ValueError, IndexError):
            pass
        out.append({"Run": row[0], "spots": spots, "bases": bases, "BioSample": bs,
                    "Platform": plat, "LibraryLayout": layout,
                    "LibraryStrategy": strat, "Model": model})
    return out


def main(metafile, out, want_accs=None):
    rows = list(csv.DictReader(open(metafile), delimiter="\t"))
    if want_accs:
        want = set(open(want_accs).read().split())
        rows = [r for r in rows
                if r["accession"].replace("_", "").replace(".", "v") in want]
    bs2acc = {r["biosample"]: r["accession"] for r in rows if r["biosample"]}
    bss = sorted(bs2acc)
    sys.stderr.write(f"{len(bss)} biosamples to query\n")

    found = []
    for i in range(0, len(bss), 100):
        batch = bss[i:i + 100]
        try:
            recs = runinfo(batch)
        except Exception as e:
            sys.stderr.write(f"  batch {i}: {e}\n"); continue
        for r in recs:
            if r.get("Platform") != "ILLUMINA": continue
            if r.get("LibraryLayout") != "PAIRED": continue
            if r.get("LibraryStrategy") not in ("WGS", "OTHER", "Other"): continue
            bs = r.get("BioSample", "")
            if bs not in bs2acc: continue
            bases, spots = r["bases"], r["spots"]
            if bases < 150_000_000:      # under ~30x on a 5 Mb genome
                continue
            found.append({
                "genome": bs2acc[bs], "biosample": bs, "run": r.get("Run", ""),
                "spots": spots, "bases": bases,
                "readlen": (bases // spots // 2) if spots else 0,
                "model": r.get("Model", ""), "strategy": r.get("LibraryStrategy", ""),
            })
        sys.stderr.write(f"  {i+len(batch)}/{len(bss)} biosamples, {len(found)} runs\n")
        time.sleep(0.4)

    # One run per genome: the deepest, which is the one most likely to assemble.
    best = {}
    for f in found:
        if f["genome"] not in best or f["bases"] > best[f["genome"]]["bases"]:
            best[f["genome"]] = f
    cols = ["genome", "biosample", "run", "spots", "bases", "readlen", "model", "strategy"]
    with open(out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for g in sorted(best):
            fh.write("\t".join(str(best[g][c]) for c in cols) + "\n")
    sys.stderr.write(f"{len(best)} genomes with paired Illumina reads -> {out}\n")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else None)
