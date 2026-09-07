#!/usr/bin/env python3
"""Compose a tesseract-model build from a selected panel.

Two modes matter:

  ship  -- every selected genome. This is the model users get. models/README.md
           notes the released Klebsiella models were fold-0 builds and that "for
           production use on arbitrary isolates a model built over the whole
           panel would be marginally stronger"; the shipped model should be that
           stronger one, with the measured numbers coming from the eval builds.

  eval  -- whole clonal clusters withheld, chromosomes AND their plasmids. Any
           gain measured against a withheld clone cannot come from having
           memorised it. Dropping only the matched genome leaves its clonemates
           behind carrying the same information, which is the leak
           docs/KLEBSIELLA_PANEL.md warns about.
"""
import argparse, csv, os, subprocess, sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--org", required=True)
    ap.add_argument("--panel", required=True, help="panel/<org> directory")
    ap.add_argument("--out", required=True)
    ap.add_argument("--tag", default="panel")
    ap.add_argument("--binary", default="./tesseract-model")
    ap.add_argument("--hold-clusters", default="",
                    help="LCO: comma-separated clonal cluster ids to withhold whole")
    ap.add_argument("--hold-accs", default="",
                    help="LIO: comma-separated safe accessions to withhold, clonemates kept")
    ap.add_argument("--marker-density", type=int, default=512)
    ap.add_argument("--min-support", type=int, default=5)
    ap.add_argument("--min-support-plasmid", type=int, default=3)
    ap.add_argument("--no-tracks", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    allrows = list(csv.DictReader(open(os.path.join(a.panel, a.tag + ".tsv")), delimiter="\t"))
    rows = [r for r in allrows if r["selected"] == "1"]

    # Redundancy-group expansion. An accession within d_dup of a held-out isolate is that
    # isolate redeposited under another accession, and its plasmid records can be
    # byte-identical to the held-out genome's own. Building the exclusion set from the
    # SELECTED rows alone leaves those deposits in the plasmid database: measured,
    # GCF022405375v1 sits at mash d=0 from hold-out GCF022405275v1, is not selected, and
    # contributed a record with the same md5 as the hold-out's. 6 of 30 S. aureus
    # hold-outs retained 8 of their own plasmid records that way.
    members = {}
    dpath = os.path.join(a.panel, "dup_members.tsv")
    if os.path.exists(dpath):
        for r in csv.DictReader(open(dpath), delimiter="\t"):
            members.setdefault(r["representative"], set()).add(r["member"])
    hold = {c.strip() for c in a.hold_clusters.split(",") if c.strip()}
    holdacc = {c.strip() for c in a.hold_accs.split(",") if c.strip()}
    if hold and holdacc:
        sys.exit("error: --hold-clusters and --hold-accs are different protocols; pick one")

    def withheld(r):
        return r["cluster"] in hold if hold else r["safe_acc"] in holdacc

    keep = [r for r in rows if not withheld(r)]
    drop = [r for r in rows if withheld(r)]

    # Every accession that must not appear anywhere in training: the withheld panel rows,
    # any unselected row that falls in a withheld cluster, and every redundancy-group
    # member of all of those.
    dropped_safe = set()
    for r in allrows:
        if withheld(r):
            dropped_safe.add(r["safe_acc"])
            dropped_safe |= members.get(r["safe_acc"], set())
    for s_ in list(dropped_safe):
        dropped_safe |= members.get(s_, set())
    if not keep:
        sys.exit("error: every panel genome was withheld")

    chrdir = os.path.join(a.panel, "chr")
    files = [os.path.join(chrdir, r["safe_acc"] + "_chr.fasta") for r in keep]
    files = [f for f in files if os.path.exists(f)]

    # Plasmids of every withheld genome, by exact record name. The whole cluster's
    # plasmids go, not just the test isolate's.
    xpls = os.path.join(a.out + ".exclude_plasmids.txt")
    npl = 0
    pmap = os.path.join(a.panel, "plasmid_map.tsv")
    with open(xpls, "w") as fh:
        if os.path.exists(pmap) and dropped_safe:
            for r in csv.DictReader(open(pmap), delimiter="\t"):
                if r["safe_acc"] in dropped_safe:
                    fh.write(r["plasmid_record"] + "\n"); npl += 1

    cmd = [a.binary, "--organism", a.org, "--out", a.out,
           "--marker-density", str(a.marker_density),
           "--min-support", str(a.min_support),
           "--min-support-plasmid", str(a.min_support_plasmid)]
    if not a.no_tracks:
        cmd.append("--layout-tracks")
    plasdb = os.path.join(a.panel, "plasmids.fna")
    if os.path.exists(plasdb) and os.path.getsize(plasdb) > 0:
        cmd += ["--plasmids", plasdb]
        if npl:
            cmd += ["--exclude-plasmids", xpls]
    for r in drop:
        cmd += ["--exclude", r["safe_acc"]]
    cmd += files

    sys.stderr.write(f"{a.org}: panel {len(rows)} selected, {len(drop)} withheld "
                     f"({len(hold)} clusters), {len(files)} chromosome files, "
                     f"{npl} plasmid records withheld\n")
    if a.dry_run:
        sys.stderr.write(" ".join(cmd[:40]) + f" ... [{len(files)} files]\n")
        return
    subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()
