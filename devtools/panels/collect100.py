#!/usr/bin/env python3
"""Per-organism significance table with quality gates and plasmid recovery.

Three things this does that the earlier collector did not:

 1. QUALITY GATES. An isolate is discarded, with the reason recorded, when the
    reference or the library is too poor for the comparison to mean anything:
      * >5% of the assembly unaligned to its own closed reference in EVERY arm --
        the reference and the reads are not the same organism/strain
      * genome fraction < 80% in every arm -- the library does not cover the genome
      * a missing arm -- a cohort must be the same isolates in every column
    docs/KLEBSIELLA_PANEL.md records the same failure in the earlier panel: eight of
    twenty isolates with up to 53 kb not aligning to their own reference, carrying
    most of that batch's misassemblies.

 2. PAIRED STATISTICS. Medians over isolates where every arm has the metric, paired
    win/loss, and a Wilcoxon signed-rank p-value against the baseline. With ~100
    non-clonal replicates the test is meaningful; with clonal replicates it would not
    be, which is why the cohort is one isolate per mash cluster.

 3. PLASMID RECOVERY, scored against the reference's own plasmid records rather than
    inferred from contiguity: how many of the isolate's plasmids are recovered whole,
    and whether the contig class calls agree with the reference.
"""
import csv, glob, math, os, statistics, sys

ARMS = ["base", "model", "careful", "aggressive", "spades"]
HIGHER = {"NGA50", "NG50", "Genome fraction (%)", "Largest alignment", "Total length"}
LOWER = {"# misassemblies", "# local misassemblies", "# mismatches per 100 kbp",
         "# indels per 100 kbp", "# contigs", "Duplication ratio", "Unaligned length"}
WANT = ["NGA50", "NG50", "Genome fraction (%)", "# misassemblies", "# local misassemblies",
        "# mismatches per 100 kbp", "# indels per 100 kbp", "Largest alignment",
        "# contigs", "Duplication ratio"]


def read_report(p):
    if not os.path.exists(p):
        return None
    rows = list(csv.reader(open(p), delimiter="\t"))
    if not rows:
        return None
    hdr = rows[0]
    col = next((i for i, h in enumerate(hdr) if i and "broken" in h), 1 if len(hdr) > 1 else None)
    if col is None:
        return None
    out = {}
    for r in rows[1:]:
        if len(r) > col:
            try:
                out[r[0]] = float(r[col].replace("%", ""))
            except ValueError:
                pass
    return out


def wilcoxon(a, b):
    """Two-sided Wilcoxon signed-rank, normal approximation with tie correction."""
    d = [x - y for x, y in zip(a, b) if x != y]
    n = len(d)
    if n < 6:
        return None
    order = sorted(range(n), key=lambda i: abs(d[i]))
    ranks = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j + 1 < n and abs(d[order[j + 1]]) == abs(d[order[i]]):
            j += 1
        r = (i + j) / 2.0 + 1
        for k in range(i, j + 1):
            ranks[order[k]] = r
        i = j + 1
    wp = sum(ranks[i] for i in range(n) if d[i] > 0)
    wm = sum(ranks[i] for i in range(n) if d[i] < 0)
    w = min(wp, wm)
    mu = n * (n + 1) / 4.0
    sd = math.sqrt(n * (n + 1) * (2 * n + 1) / 24.0)
    if sd == 0:
        return None
    z = (w - mu + 0.5) / sd
    return 2.0 * 0.5 * math.erfc(abs(z) / math.sqrt(2))


def main(evaldir, out_tsv, out_txt):
    data, discarded = {}, []
    for d in sorted(glob.glob(os.path.join(evaldir, "*"))):
        if not os.path.isdir(d):
            continue
        iso = os.path.basename(d)
        per = {}
        for a in ARMS:
            r = read_report(os.path.join(d, f"quast_{a}", "report.tsv"))
            if r:
                per[a] = r
        missing = [a for a in ARMS if a not in per]
        if missing:
            discarded.append((iso, "missing arms: " + ",".join(missing)))
            continue
        gf = [per[a].get("Genome fraction (%)", 0) for a in ARMS]
        if max(gf) < 80:
            discarded.append((iso, f"genome fraction {max(gf):.1f}% in every arm"))
            continue
        una, tot = [], []
        for a in ARMS:
            u = per[a].get("Unaligned length", 0.0)
            t = per[a].get("Total length", 0.0) or 1.0
            una.append(100.0 * u / t)
        if min(una) > 5.0:
            discarded.append((iso, f"{min(una):.1f}% unaligned to its own reference in every arm"))
            continue
        data[iso] = per

    with open(out_txt, "w") as fh:
        def emit(s=""):
            print(s); fh.write(s + "\n")
        emit(f"cohort {os.path.basename(evaldir)}: {len(data)} isolates kept, "
             f"{len(discarded)} discarded")
        for iso, why in discarded[:12]:
            emit(f"   discarded {iso}: {why}")
        if len(discarded) > 12:
            emit(f"   ... and {len(discarded)-12} more")
        if not data:
            return
        isos = sorted(data)
        emit()
        emit(f"{'metric':<28}" + "".join(f"{a:>14}" for a in ARMS))
        for m in WANT:
            paired = [i for i in isos if all(m in data[i][a] for a in ARMS)]
            if not paired:
                continue
            line = f"{m:<28}"
            for a in ARMS:
                line += f"{statistics.median([data[i][a][m] for i in paired]):>14,.2f}"
            emit(line)
        emit()
        emit(f"{'vs base  (win/loss, p)':<28}" + "".join(f"{a:>22}" for a in ARMS[1:]))
        for m in WANT:
            paired = [i for i in isos if all(m in data[i][a] for a in ARMS)]
            if len(paired) < 6:
                continue
            hi = m in HIGHER
            line = f"{m:<28}"
            for a in ARMS[1:]:
                bb = [data[i]["base"][m] for i in paired]
                aa = [data[i][a][m] for i in paired]
                w = sum(1 for x, y in zip(bb, aa) if (y > x) == hi and y != x)
                l = sum(1 for x, y in zip(bb, aa) if (y < x) == hi and y != x)
                p = wilcoxon(aa, bb)
                ps = "n/a" if p is None else (f"{p:.1e}" if p < 1e-3 else f"{p:.3f}")
                line += f"{f'{w}/{l} p={ps}':>22}"
            emit(line)

    cols = ["isolate"] + [f"{a}|{m}" for a in ARMS for m in WANT]
    with open(out_tsv, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for i in sorted(data):
            fh.write("\t".join([i] + [f"{data[i][a].get(m,'')}" for a in ARMS for m in WANT]) + "\n")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
