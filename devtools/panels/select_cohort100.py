#!/usr/bin/env python3
"""Select 100 NON-CLONAL isolates per organism for the definitive campaign.

Non-clonal is enforced structurally: at most one isolate per mash cluster at
d <= 0.0005, so no two isolates in a cohort are clonemates of each other. That is
what makes 100 replicates 100 independent observations rather than a reweighted
outbreak, and it is the condition a significance test needs.

Selection within a cluster prefers, in order:
  1. an isolate already assembled in the 30-isolate run -- those arms can be reused
     instead of recomputed, and reusing them changes nothing about the cohort's
     independence
  2. clinical provenance
  3. depth nearest 60x
Read length must be >= 100 and depth in [30, 200]: outside that band the library is
not what a diagnostic laboratory produces and the comparison is not about the
assembler.
"""
import argparse, csv, math, os, sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--org", required=True)
    ap.add_argument("--panel", required=True)
    ap.add_argument("--reads", required=True)
    ap.add_argument("--n", type=int, default=100)
    ap.add_argument("--reuse-dir", default="", help="existing eval dir whose isolates to prefer")
    ap.add_argument("--min-depth", type=float, default=30.0)
    ap.add_argument("--max-depth", type=float, default=200.0)
    ap.add_argument("--min-readlen", type=int, default=100)
    ap.add_argument("--target-depth", type=float, default=60.0)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    panel = list(csv.DictReader(open(os.path.join(a.panel, "panel_clonal.tsv")), delimiter="\t"))
    reads = {r["genome"].replace("_", "").replace(".", "v"): r
             for r in csv.DictReader(open(a.reads), delimiter="\t")}
    reuse = set()
    if a.reuse_dir and os.path.isdir(a.reuse_dir):
        reuse = {d for d in os.listdir(a.reuse_dir)
                 if os.path.isdir(os.path.join(a.reuse_dir, d))}

    by_cluster = {}
    for r in panel:
        s = r["safe_acc"]
        rd = reads.get(s)
        if not rd:
            continue
        try:
            glen = float(r["total_len"]); bases = int(rd["bases"]); rl = int(rd["readlen"])
        except (ValueError, KeyError):
            continue
        if glen <= 0 or rl < a.min_readlen:
            continue
        depth = bases / glen
        if not (a.min_depth <= depth <= a.max_depth):
            continue
        key = (0 if s in reuse else 1,
               0 if r["clinical"] == "1" else 1,
               abs(math.log(depth / a.target_depth)))
        c = r["cluster"]
        if c not in by_cluster or key < by_cluster[c][0]:
            by_cluster[c] = (key, s, rd, r, depth)

    # Rank clusters so the cohort keeps reusable and clinical isolates first, but the
    # cohort is still one-per-cluster whichever way it is truncated.
    cands = sorted(by_cluster.values(), key=lambda t: t[0])
    if len(cands) < a.n:
        sys.stderr.write(f"{a.org}: WARNING only {len(cands)} non-clonal candidates, "
                         f"wanted {a.n}\n")
    chosen = cands[:a.n]

    with open(a.out, "w") as fh:
        fh.write("safe_acc\tcluster\trun\tbases\treadlen\tdepth\tclinical\tmodel\treused\n")
        for _, s, rd, r, depth in chosen:
            fh.write(f"{s}\t{r['cluster']}\t{rd['run']}\t{rd['bases']}\t{rd['readlen']}\t"
                     f"{depth:.1f}\t{r['clinical']}\t{rd.get('model','')}\t"
                     f"{1 if s in reuse else 0}\n")
    open(a.out + ".clusters", "w").write(",".join(sorted({r["cluster"] for _, _, _, r, _ in chosen})))
    open(a.out + ".accs", "w").write(",".join(s for _, s, _, _, _ in chosen))
    nre = sum(1 for _, s, _, _, _ in chosen if s in reuse)
    ncl = sum(1 for _, _, _, r, _ in chosen if r["clinical"] == "1")
    sys.stderr.write(f"{a.org}: {len(chosen)} non-clonal isolates in {len(chosen)} distinct "
                     f"clusters, {ncl} clinical, {nre} reusable from the 30-isolate run\n")


if __name__ == "__main__":
    main()
