#!/usr/bin/env python3
"""Aggregate one organism's evaluation: every model config against the baseline,
split by whether the isolate's lineage is represented in the panel.

Only the *_broken* QUAST row is read (scaffolds split at N-runs >= 10). Paired
per-isolate win/loss counts are reported next to medians because the two can
disagree honestly -- docs/MODELS_INTEGRATED.md has a case where medians favour
one model and the paired comparison the other, and quoting either alone is a
choice about which flatters.
"""
import csv, glob, os, statistics, sys

HIGHER_BETTER = {"NGA50", "NG50", "NA50", "N50", "Genome fraction (%)", "Largest alignment"}
LOWER_BETTER = {"# misassemblies", "# local misassemblies", "# mismatches per 100 kbp",
                "# indels per 100 kbp", "# contigs", "Duplication ratio"}
# "# local misassemblies" is here because leaving it out hid a real regression: it rises
# from a cohort sum of 6 to 18-23 with the model, and was invisible in every table until
# it was asked for by name.
WANT = ["NGA50", "NG50", "Genome fraction (%)", "# misassemblies", "# local misassemblies",
        "# mismatches per 100 kbp", "# indels per 100 kbp", "Largest alignment",
        "# contigs", "Duplication ratio"]


def read_report(path):
    if not os.path.exists(path):
        return None
    rows = list(csv.reader(open(path), delimiter="\t"))
    if not rows:
        return None
    hdr, col = rows[0], None
    for i, name in enumerate(hdr):
        if i and "broken" in name:
            col = i
    if col is None:
        col = 1 if len(hdr) > 1 else None
    if col is None:
        return None
    out = {}
    for r in rows[1:]:
        if r and r[0] in WANT and len(r) > col:
            try:
                out[r[0]] = float(r[col].replace("%", ""))
            except ValueError:
                pass
    return out


def main(evaldir, out, holdout):
    strat = {}
    if os.path.exists(holdout):
        for r in csv.DictReader(open(holdout), delimiter="\t"):
            strat[r["safe_acc"]] = r["stratum"]

    data, configs = {}, set()
    for d in sorted(glob.glob(os.path.join(evaldir, "*"))):
        if not os.path.isdir(d):
            continue
        iso = os.path.basename(d)
        per = {}
        for q in glob.glob(os.path.join(d, "quast_*")):
            cfg = os.path.basename(q)[len("quast_"):]
            rep = read_report(os.path.join(q, "report.tsv"))
            if rep:
                per[cfg] = rep; configs.add(cfg)
        if "base" in per and len(per) > 1:
            data[iso] = per
    if not data:
        sys.stderr.write("no scored isolates yet\n"); return
    configs = ["base"] + sorted(c for c in configs if c != "base")

    with open(out, "w") as fh:
        cols = ["isolate", "stratum"] + [f"{c}|{m}" for c in configs for m in WANT]
        fh.write("\t".join(cols) + "\n")
        for iso in sorted(data):
            row = [iso, strat.get(iso, "?")]
            for c in configs:
                for m in WANT:
                    row.append(f"{data[iso].get(c, {}).get(m, '')}")
            fh.write("\t".join(row) + "\n")

    def block(title, isos):
        if not isos:
            return
        # Medians are computed over the isolates where EVERY config has the metric, and
        # the n is printed. Taking each column's median over whatever isolates happen to
        # have it compares different cohorts: SPAdes crashed on 13 of 30 A. baumannii
        # isolates, so its column was a median over 17 and the baseline's over 30, and
        # the 13 missing were exactly the chemistry where TesserACT does worst. That
        # inflated a reported gap by roughly 60%. Per-column coverage is printed too, so
        # a crashed arm is visible rather than silently changing the denominator.
        print(f"\n=== {title} (n={len(isos)}) ===")
        cov = {c: sum(1 for i in isos if c in data[i]) for c in configs}
        print(f"{'metric':<28}" + "".join(f"{c:>17}" for c in configs))
        print(f"{'(isolates with output)':<28}" + "".join(f"{cov[c]:>17}" for c in configs))
        for m in WANT:
            paired = [i for i in isos
                      if all(c in data[i] and m in data[i][c] for c in configs)]
            line = f"{m:<28}"
            for c in configs:
                v = [data[i][c][m] for i in paired]
                line += f"{statistics.median(v):>17,.2f}" if v else f"{'-':>17}"
            if len(paired) != len(isos):
                line += f"   [paired n={len(paired)}]"
            print(line)
        print(f"\n{'win/loss vs base':<28}" + "".join(f"{c:>17}" for c in configs))
        for m in WANT:
            if m not in HIGHER_BETTER and m not in LOWER_BETTER:
                continue
            hi = m in HIGHER_BETTER
            line = f"{m:<28}"
            for c in configs:
                if c == "base":
                    line += f"{'-':>17}"; continue
                w = l = 0
                for i in isos:
                    if c not in data[i] or m not in data[i][c] or m not in data[i]["base"]:
                        continue
                    a, b = data[i]["base"][m], data[i][c][m]
                    if a == b: continue
                    better = b > a if hi else b < a
                    w, l = (w + 1, l) if better else (w, l + 1)
                line += f"{f'{w}/{l}':>17}"
            print(line)

    allisos = sorted(data)
    block("ALL held-out isolates", allisos)
    block("lineage REPRESENTED in panel (clonal)", [i for i in allisos if strat.get(i) == "clonal"])
    block("lineage ABSENT from panel (singleton)", [i for i in allisos if strat.get(i) == "singleton"])


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else "")
