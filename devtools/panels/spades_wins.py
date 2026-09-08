#!/usr/bin/env python3
"""Find the isolates where SPAdes beats every TesserACT arm, and describe them.

The point is not the count but the pattern: if the losses share a property -- a
library geometry, a coverage regime, a genome feature -- that property is the lead for
a code fix. A scattered set of losses means there is no single cause to chase.

Reports, per organism:
  * how often SPAdes wins each metric against the BEST TesserACT arm
  * the covariates of those isolates against the ones we win: insert size relative to
    read length, depth, read length, reference plasmid count, genome fraction
  * the isolates themselves, so the follow-up analysis has a work list
"""
import csv, glob, json, os, statistics, sys

ARMS = ["base", "model", "careful", "aggressive"]
HIGHER = {"NGA50", "NG50", "Genome fraction (%)", "Largest alignment"}


def load_quast(d, arm):
    p = os.path.join(d, f"quast_{arm}", "report.tsv")
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


def covariates(d, arm="base"):
    p = os.path.join(d, arm, "report.json")
    if not os.path.exists(p):
        return {}
    try:
        j = json.load(open(p))
    except Exception:
        return {}
    inp = j.get("input", {}) or {}
    rr = j.get("repeat_resolution", {}) or {}
    ins = (rr.get("insert_size") or {})
    reads = inp.get("reads") or 0
    bases = inp.get("bases") or 0
    meanlen = (bases / reads) if reads else 0
    return {
        "mean_read_len": meanlen,
        "max_read_len": inp.get("max_read_length", 0),
        "insert_mean": ins.get("mean", 0),
        "ins_over_len": (ins.get("mean", 0) / meanlen) if meanlen else 0,
        "linking_pairs": rr.get("linking_pairs", 0),
        "reads_anchored": rr.get("reads_anchored", 0),
        "masked_frac": (inp.get("quality_trimmed_bases", 0) / bases) if bases else 0,
    }


def main(org, evaldir, out_tsv):
    rows = []
    for d in sorted(glob.glob(os.path.join(evaldir, "*"))):
        if not os.path.isdir(d):
            continue
        sp = load_quast(d, "spades")
        arms = {a: load_quast(d, a) for a in ARMS}
        if not sp or any(v is None for v in arms.values()):
            continue
        best = {}
        for m in ("NGA50", "Genome fraction (%)", "# misassemblies", "# contigs"):
            vals = [(a, arms[a].get(m)) for a in ARMS if arms[a].get(m) is not None]
            if not vals:
                continue
            best[m] = max(vals, key=lambda t: t[1])[1] if m in HIGHER else min(vals, key=lambda t: t[1])[1]
        if "NGA50" not in best or "NGA50" not in sp:
            continue
        cov = covariates(d)
        rows.append({
            "isolate": os.path.basename(d),
            "spades_nga50": sp["NGA50"], "best_nga50": best["NGA50"],
            "ratio": sp["NGA50"] / best["NGA50"] if best["NGA50"] else 0,
            "spades_wins_nga50": int(sp["NGA50"] > best["NGA50"]),
            "spades_gf": sp.get("Genome fraction (%)", 0), "best_gf": best.get("Genome fraction (%)", 0),
            "spades_mis": sp.get("# misassemblies", 0), "best_mis": best.get("# misassemblies", 0),
            **cov,
        })
    if not rows:
        print(f"{org}: no isolate has all arms scored yet")
        return
    with open(out_tsv, "w") as fh:
        cols = list(rows[0].keys())
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    wins = [r for r in rows if r["spades_wins_nga50"]]
    loss = [r for r in rows if not r["spades_wins_nga50"]]
    print(f"{org}: n={len(rows)}   SPAdes wins NGA50 on {len(wins)}  ({100*len(wins)/len(rows):.0f}%)")
    if not wins or not loss:
        return
    print(f"\n{'covariate':<18}{'SPAdes wins':>14}{'we win':>12}{'ratio':>9}")
    for k in ("ins_over_len", "insert_mean", "mean_read_len", "max_read_len",
              "linking_pairs", "masked_frac"):
        a = statistics.median([r[k] for r in wins if r.get(k)])
        b = statistics.median([r[k] for r in loss if r.get(k)])
        print(f"{k:<18}{a:>14.3f}{b:>12.3f}{(a/b if b else 0):>9.2f}")
    print(f"\nworst 10 losses (SPAdes/best NGA50):")
    for r in sorted(wins, key=lambda r: -r["ratio"])[:10]:
        print(f"  {r['isolate']:<18} ratio {r['ratio']:>5.2f}  spades {r['spades_nga50']:>9,.0f} "
              f"best {r['best_nga50']:>9,.0f}  ins/len {r.get('ins_over_len',0):.2f}  "
              f"maxlen {r.get('max_read_len',0):.0f}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
