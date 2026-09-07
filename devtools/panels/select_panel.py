#!/usr/bin/env python3
"""Choose a panel that keeps clonal structure instead of dereplicating it away.

The shipped panels were mash-dereplicated at d >= 0.0005, one genome per
cluster, on the argument that "100 near-identical ST239 S. aureus genomes would
teach one arrangement rather than a species" (docs/ESKAPEE_MODELS.md). That is
right if the model must describe a species. It is the wrong objective here,
because the model is not asked to describe a species -- it is asked to place the
contigs of one clinical isolate, and clinical isolates are drawn from a handful
of high-prevalence lineages rather than uniformly from the species.

Dereplication costs two things that are mechanically identifiable in the code:

  1. Layout tracks (--layout-tracks) lay an assembly out against "the relative it
     most resembles". Dereplication deletes precisely the close relatives.
  2. --min-support 5 drops chromosomal adjacencies seen in fewer than five panel
     genomes. A junction carried only by one lineage -- a prophage or IS
     insertion specific to ST258 -- is present in exactly one genome of a
     dereplicated panel and is therefore DROPPED at finalise(). With twenty
     clonemates it clears the bar and survives.

So the panel keeps clonal structure, with a per-cluster cap so that one outbreak
clone cannot swamp the panel -- the original concern, addressed by bounding the
dominance rather than by discarding the lineage.

"Non-redundant" is enforced separately and strictly: genomes within d_dup of each
other are the same isolate redeposited, and exactly one survives.

Leakage is handled at evaluation time, not here: folds withhold whole clonal
clusters, so a gain measured against a held-out clone cannot come from having
memorised it. See build_folds.py.
"""
import argparse, csv, os, subprocess, sys
from collections import defaultdict


class DSU:
    def __init__(self, n): self.p = list(range(n))
    def find(self, a):
        while self.p[a] != a:
            self.p[a] = self.p[self.p[a]]; a = self.p[a]
        return a
    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb: self.p[ra] = rb


def clusters_from(names, edges, thresh):
    idx = {n: i for i, n in enumerate(names)}
    d = DSU(len(names))
    for a, b, dist in edges:
        if dist <= thresh:
            d.union(idx[a], idx[b])
    out = defaultdict(list)
    for n in names:
        out[d.find(idx[n])].append(n)
    return list(out.values())


CLINICAL = ("blood", "urine", "sputum", "wound", "clinical", "patient", "hospital",
            "csf", "abscess", "catheter", "respiratory", "stool", "rectal", "pus",
            "swab", "bronch", "tracheal", "peritoneal", "bile", "human")


def is_clinical(row):
    h = row.get("host", "").lower()
    s = row.get("isolation_source", "").lower()
    if "homo sapiens" in h or "human" in h:
        return True
    return any(t in s for t in CLINICAL)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--org", required=True)
    ap.add_argument("--meta", required=True)
    ap.add_argument("--chrdir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--target", type=int, default=500)
    ap.add_argument("--cap", type=int, default=16,
                    help="max genomes per clonal cluster; 1 reproduces the "
                         "shipped dereplication, 0 keeps every non-redundant genome")
    ap.add_argument("--tag", default="panel")
    ap.add_argument("--d-dup", type=float, default=0.00002,
                    help="at or below this two assemblies are the same isolate")
    ap.add_argument("--d-clone", type=float, default=0.0005,
                    help="clonal cluster; the threshold the shipped panels dereplicated at")
    ap.add_argument("--threads", type=int, default=32)
    ap.add_argument("--reads", default="",
                    help="reads TSV; used ONLY to break ties between redeposits "
                         "of the same isolate, never to weight the panel")
    ap.add_argument("--min-len", type=int, default=0)
    ap.add_argument("--max-len", type=int, default=0)
    a = ap.parse_args()

    meta = {}
    for r in csv.DictReader(open(a.meta), delimiter="\t"):
        meta[r["accession"].replace("_", "").replace(".", "v")] = r

    files = sorted(f for f in os.listdir(a.chrdir) if f.endswith("_chr.fasta"))
    names = [f[:-len("_chr.fasta")] for f in files]

    # Size QC. A "complete" assembly whose length is far off the species is a
    # mislabel or a contaminated submission, and it would contribute adjacencies
    # that are not this organism's.
    keep = []
    for n in names:
        m = meta.get(n)
        if not m: continue
        try: L = int(m["total_len"])
        except ValueError: continue
        if a.min_len and L < a.min_len: continue
        if a.max_len and L > a.max_len: continue
        keep.append(n)
    sys.stderr.write(f"{a.org}: {len(names)} chromosomes, {len(keep)} pass size QC\n")
    names = keep

    os.makedirs(a.out, exist_ok=True)
    lst = os.path.join(a.out, "sketch_input.txt")
    with open(lst, "w") as fh:
        for n in names: fh.write(os.path.join(a.chrdir, n + "_chr.fasta") + "\n")

    msh = os.path.join(a.out, "panel.msh")
    if not os.path.exists(msh):
        sys.stderr.write(f"{a.org}: mash sketch\n")
        subprocess.run(["mash", "sketch", "-p", str(a.threads), "-s", "10000", "-k", "21",
                        "-l", lst, "-o", msh[:-4]], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    tri = os.path.join(a.out, "triangle.tsv")
    if not os.path.exists(tri):
        sys.stderr.write(f"{a.org}: mash triangle ({len(names)} genomes)\n")
        with open(tri, "w") as fh:
            subprocess.run(["mash", "triangle", "-p", str(a.threads), msh],
                           check=True, stdout=fh, stderr=subprocess.DEVNULL)

    # mash triangle: line i holds the distances from genome i to genomes 0..i-1.
    edges, order = [], []
    with open(tri) as fh:
        next(fh)
        for line in fh:
            p = line.rstrip("\n").split("\t")
            nm = os.path.basename(p[0])[:-len("_chr.fasta")]
            for j, v in enumerate(p[1:]):
                d = float(v)
                if d <= a.d_clone:            # only near edges matter for either cut
                    edges.append((nm, order[j], d))
            order.append(nm)

    dup = clusters_from(names, edges, a.d_dup)
    sys.stderr.write(f"{a.org}: {len(names)} -> {len(dup)} non-redundant (d<={a.d_dup})\n")

    havereads = set()
    if a.reads and os.path.exists(a.reads):
        for r in csv.DictReader(open(a.reads), delimiter="\t"):
            havereads.add(r["genome"].replace("_", "").replace(".", "v"))

    def rank(n):
        m = meta.get(n, {})
        # Clinical first: the panel is meant to mirror what a diagnostic lab sees.
        # Then RefSeq's own curation, then the most recent deposit.
        cat = m.get("refseq_category", "")
        return (0 if is_clinical(m) else 1,
                0 if cat in ("reference genome", "representative genome") else 1,
                m.get("release_date", ""))

    # Choosing which accession represents a redeposit group: prefer the copy that
    # has public Illumina reads. Members of a group are the SAME isolate at
    # d<=d_dup, so this changes which accession stands for it and nothing about
    # what the panel contains -- but it decides whether that isolate can ever be
    # evaluated on real reads. Without it only 19 of 811 E. faecium panel genomes
    # were evaluable: its complete genomes come largely from clonal outbreak
    # studies whose members collapse into one another, and the representative was
    # being chosen without reference to whether reads existed.
    def rank_dup(n):
        return (0 if n in havereads else 1,) + rank(n)

    reps = [sorted(g, key=rank_dup)[0] for g in dup]
    dupsize = {sorted(g, key=rank_dup)[0]: len(g) for g in dup}

    clon = clusters_from(reps, [e for e in edges if e[0] in set(reps) and e[1] in set(reps)],
                         a.d_clone)
    clon.sort(key=len, reverse=True)
    sys.stderr.write(f"{a.org}: {len(reps)} non-redundant -> {len(clon)} clonal clusters "
                     f"(d<={a.d_clone}); largest {len(clon[0])}\n")

    # The cap is set from the mechanism, not from a target count. --min-support
    # drops a chromosomal adjacency seen in fewer than N panel genomes (N=5 by
    # default), so a lineage-specific junction needs at least N clonemates to
    # survive finalise(). A cap of 1 -- one genome per cluster -- is precisely the
    # shipped dereplication, and under it EVERY lineage-specific junction is seen
    # once and dropped. The default 16 is 3x min-support: enough margin that the
    # junction survives, while bounding any one outbreak clone's share of the
    # panel, which was the original and correct objection to keeping them all.
    cap = a.cap if a.cap > 0 else 10**9
    chosen = []
    for c in clon:
        chosen.extend(sorted(c, key=rank)[:cap])
    if len(chosen) < a.target:
        sys.stderr.write(f"{a.org}: WARNING only {len(chosen)} genomes, below target {a.target}\n")
    sys.stderr.write(f"{a.org}: cap={a.cap} -> panel of {len(chosen)}\n")

    cl_of = {}
    for i, c in enumerate(clon):
        for n in c: cl_of[n] = i
    sel = set(chosen)
    with open(os.path.join(a.out, a.tag + ".tsv"), "w") as fh:
        fh.write("safe_acc\taccession\tcluster\tcluster_size\tdup_size\tselected\t"
                 "clinical\thost\tisolation_source\ttotal_len\n")
        for n in reps:
            m = meta.get(n, {})
            fh.write(f"{n}\t{m.get('accession','')}\t{cl_of.get(n,-1)}\t"
                     f"{len(clon[cl_of[n]]) if n in cl_of else 0}\t{dupsize.get(n,1)}\t"
                     f"{1 if n in sel else 0}\t{1 if is_clinical(m) else 0}\t"
                     f"{m.get('host','')}\t{m.get('isolation_source','')}\t{m.get('total_len','')}\n")
    with open(os.path.join(a.out, a.tag + "_selected.txt"), "w") as fh:
        for n in sorted(sel): fh.write(n + "\n")

    # Redundancy-group membership, written once per panel directory. Withholding an
    # isolate has to withhold every accession within d_dup of it -- those are the SAME
    # isolate redeposited, and one of them contributes plasmid records byte-identical to
    # the held-out genome's own. The panel tsv lists only group REPRESENTATIVES, so an
    # exclusion built from it alone leaves the non-representative deposits in training.
    with open(os.path.join(a.out, "dup_members.tsv"), "w") as fh:
        fh.write("representative\tmember\n")
        for g in dup:
            rep = sorted(g, key=rank_dup)[0]
            for m in g:
                fh.write(f"{rep}\t{m}\n")

    nclin = sum(1 for n in sel if is_clinical(meta.get(n, {})))
    ncl = len({cl_of[n] for n in sel if n in cl_of})
    sys.stderr.write(f"{a.org}: SELECTED {len(sel)} genomes over {ncl} clonal clusters, "
                     f"{nclin} clinical ({100*nclin/max(1,len(sel)):.0f}%)\n")


if __name__ == "__main__":
    main()
