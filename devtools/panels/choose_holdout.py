#!/usr/bin/env python3
"""Pick held-out isolates, stratified by how clonal their lineage is.

Two protocols are run over the same isolates and the same assemblies, because
they answer different questions and only one of them is usually reported:

  LCO -- leave-CLONE-out. The isolate's whole clonal cluster is withheld,
         chromosomes and plasmids. This is the honest generalisation test and
         it is what docs/KLEBSIELLA_PANEL.md demands ("exclude whole
         mash-distance clusters rather than individual accessions"). It also
         DEFINES AWAY the effect being measured: with every clonemate removed,
         a clonal panel cannot use clonemates. LCO is therefore a lower bound
         on what clonality buys, not a measurement of it.

  LIO -- leave-ISOLATE-out. The test isolate's own chromosome and its own
         plasmids are withheld; its clonemates -- which are different isolates,
         independently deposited -- stay. This is the deployed situation: a
         clinical isolate of a prevalent lineage arrives and the panel already
         holds that lineage. It is only honest if the isolate's own sequence is
         genuinely gone, which --exclude plus --exclude-plasmids enforces and
         the model file records.

Reporting one without the other would be a choice about which number flatters.
Both are produced.

Isolates are stratified: half from clusters with enough clonemates to clear
--min-support, half from singletons, so the two regimes are separable in the
results rather than averaged together.
"""
import argparse, csv, os, sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", required=True)
    ap.add_argument("--tag", default="panel_clonal")
    ap.add_argument("--reads", required=True)
    ap.add_argument("--n", type=int, default=30)
    ap.add_argument("--min-support", type=int, default=5)
    ap.add_argument("--max-withheld", type=float, default=0.15)
    ap.add_argument("--min-readlen", type=int, default=100)
    ap.add_argument("--min-depth", type=float, default=30.0)
    ap.add_argument("--max-depth", type=float, default=200.0)
    ap.add_argument("--target-depth", type=float, default=60.0)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    panel = list(csv.DictReader(open(os.path.join(a.panel, a.tag + ".tsv")), delimiter="\t"))
    sel = [r for r in panel if r["selected"] == "1"]
    by_safe = {r["safe_acc"]: r for r in sel}
    size = {}
    for r in sel:
        size[r["cluster"]] = size.get(r["cluster"], 0) + 1

    reads = {}
    for r in csv.DictReader(open(a.reads), delimiter="\t"):
        if int(r["readlen"]) < a.min_readlen:
            continue
        reads[r["genome"].replace("_", "").replace(".", "v")] = r

    # Depth band. Ranking candidates by "deepest library first" selects the tail of
    # the distribution -- here 221x to 2889x -- which is neither what a diagnostic
    # laboratory produces nor what the assembler should be characterised on, and a
    # 2889x library on a 2.8 Mb genome is 8 Gbp of reads per assembly.
    # docs/REAL_WORLD_RESULTS.md records the converse error: a ~100x fetch cap that
    # "should not have existed, since real coverage spans 24x-267x and that spread is
    # the condition being measured". So: keep the natural spread inside a realistic
    # band, and rank by closeness to a typical depth rather than by extremity.
    import math
    cand = []
    for s in by_safe:
        if s not in reads:
            continue
        try:
            glen = float(by_safe[s]["total_len"])
        except (ValueError, KeyError):
            continue
        if glen <= 0:
            continue
        depth = int(reads[s]["bases"]) / glen
        if not (a.min_depth <= depth <= a.max_depth):
            continue
        cand.append((s, reads[s], depth))
    # "Clonal" means the lineage has enough panel members that withholding this
    # one isolate still leaves min-support satisfiable.
    clonal = [c for c in cand if size.get(by_safe[c[0]]["cluster"], 1) > a.min_support]
    single = [c for c in cand if size.get(by_safe[c[0]]["cluster"], 1) == 1]
    sys.stderr.write(f"{len(sel)} panel genomes; {len(cand)} with reads in "
                     f"{a.min_depth:.0f}-{a.max_depth:.0f}x "
                     f"({len(clonal)} clonal lineages, {len(single)} singletons)\n")

    def key(t):
        return (0 if by_safe[t[0]]["clinical"] == "1" else 1,
                abs(math.log(t[2] / a.target_depth)))
    clonal.sort(key=key)
    single.sort(key=key)

    budget = int(len(sel) * a.max_withheld)
    half = a.n // 2
    chosen, used, spent = [], set(), 0
    for pool, want, stratum in ((clonal, half, "clonal"), (single, a.n - half, "singleton")):
        got = 0
        for s, rd, depth in pool:
            if got >= want:
                break
            c = by_safe[s]["cluster"]
            if c in used:
                continue
            if spent + size.get(c, 1) > budget:
                continue
            chosen.append((s, c, rd, stratum, depth))
            used.add(c); spent += size.get(c, 1); got += 1
        if got < want:
            sys.stderr.write(f"  only {got}/{want} available in stratum '{stratum}'\n")

    with open(a.out, "w") as fh:
        fh.write("safe_acc\tcluster\tcluster_size\tstratum\trun\tbases\treadlen\t"
                 "model\tclinical\tdepth\n")
        for s, c, rd, st, dp in chosen:
            fh.write(f"{s}\t{c}\t{size.get(c,1)}\t{st}\t{rd['run']}\t{rd['bases']}\t"
                     f"{rd['readlen']}\t{rd['model']}\t{by_safe[s]['clinical']}\t{dp:.1f}\n")
    open(a.out + ".clusters", "w").write(",".join(sorted(used)))
    open(a.out + ".accs", "w").write(",".join(s for s, _, _, _, _ in chosen))
    nc = sum(1 for _, _, _, st, _ in chosen if st == "clonal")
    sys.stderr.write(f"held out {len(chosen)} isolates ({nc} clonal, {len(chosen)-nc} singleton) "
                     f"in {len(used)} clusters, costing {spent}/{len(sel)} panel genomes "
                     f"({100*spent/len(sel):.1f}%) under LCO\n")


if __name__ == "__main__":
    main()
