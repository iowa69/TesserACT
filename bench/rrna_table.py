#!/usr/bin/env python3
"""Where the unspannable repeats are, and what they look like.

The assembler's contiguity is decided at repeats longer than a fragment. This
builds the table of what those actually are in K. pneumoniae: rRNA operons
located by barrnap, their positions, spacing, GC content and the GC of their
flanks, plus every other repeat found by self-alignment.

Output:
  rrna_operons.tsv   one row per operon: isolate, replicon, start, end, genes, GC
  rrna_summary.tsv   per isolate: operon count, median spacing, GC contrast
  repeats.tsv        every repeat copy >=500 bp from self-alignment, with length/identity
"""
import os
import re
import subprocess
import sys
from collections import defaultdict

R = "/media/iowa/WD_BLACK/kle_bench"
BARRNAP = "/home/iowa/miniconda3/envs/bsi/bin/barrnap"
# barrnap is a Perl script; without its own lib on PERL5LIB it aborts at load.
import glob as _glob
_ENV = "/home/iowa/miniconda3/envs/bsi"
_ENVIRON = dict(os.environ)
_ENVIRON["PERL5LIB"] = ":".join(_glob.glob(f"{_ENV}/lib/perl5/*/"))
_ENVIRON["PATH"] = f"{_ENV}/bin:" + _ENVIRON.get("PATH", "")
MM = "/home/iowa/miniconda3/envs/benchtools/bin/minimap2"
THREADS = 4


def load_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def gc(s):
    if not s:
        return 0.0
    g = sum(1 for c in s if c in "GCgc")
    n = sum(1 for c in s if c in "ACGTacgt")
    return 100.0 * g / n if n else 0.0


def operons_from_gff(gff, contig_len):
    """Group barrnap rRNA hits into operons: 16S/23S/5S within 6 kb of each other."""
    hits = []
    for line in gff.splitlines():
        if line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 9:
            continue
        m = re.search(r"Name=([^;]+)", f[8])
        gene = m.group(1) if m else "rRNA"
        hits.append((f[0], int(f[3]), int(f[4]), gene, f[6]))
    hits.sort(key=lambda h: (h[0], h[1]))

    operons = []
    cur = []
    for h in hits:
        if cur and h[0] == cur[-1][0] and h[1] - cur[-1][2] <= 6000:
            cur.append(h)
        else:
            if cur:
                operons.append(cur)
            cur = [h]
    if cur:
        operons.append(cur)
    return operons


def main():
    ids = sorted(f[:-4] for f in os.listdir(f"{R}/refs") if f.endswith(".fna"))
    limit = int(sys.argv[1]) if len(sys.argv) > 1 else len(ids)
    ids = ids[:limit]

    op_rows, sum_rows, rep_rows = [], [], []
    for n, s in enumerate(ids, 1):
        ref = f"{R}/refs/{s}.fna"
        seqs = load_fasta(ref)
        if not seqs:
            continue
        genome_gc = gc("".join(seqs.values()))

        try:
            gff = subprocess.run([BARRNAP, "--threads", str(THREADS), "--quiet", ref],
                                 capture_output=True, text=True, timeout=900, env=_ENVIRON).stdout
        except Exception:
            continue

        ops = operons_from_gff(gff, seqs)
        # Only full operons: a lone 5S is not what breaks an assembly.
        full = [o for o in ops if len(o) >= 2]
        starts = []
        for o in full:
            contig = o[0][0]
            beg, end = o[0][1], o[-1][2]
            seq = seqs.get(contig, "")
            body = seq[beg - 1:end]
            left = seq[max(0, beg - 1 - 2000):beg - 1]
            right = seq[end:end + 2000]
            genes = ",".join(sorted({h[3].split("_")[0] for h in o}))
            op_rows.append((s, contig, beg, end, end - beg + 1, len(o), genes,
                            f"{gc(body):.1f}", f"{gc(left):.1f}", f"{gc(right):.1f}"))
            if contig == max(seqs, key=lambda k: len(seqs[k])):
                starts.append(beg)

        starts.sort()
        gaps = [starts[i + 1] - starts[i] for i in range(len(starts) - 1)]
        med_gap = sorted(gaps)[len(gaps) // 2] if gaps else 0
        op_gc = [float(r[7]) for r in op_rows if r[0] == s]
        sum_rows.append((s, len(seqs), sum(len(v) for v in seqs.values()), len(full),
                         med_gap, f"{genome_gc:.2f}",
                         f"{sum(op_gc)/len(op_gc):.1f}" if op_gc else "-"))

        # Repeat inventory by self-alignment.
        try:
            out = subprocess.run([MM, "-x", "asm10", "-t", str(THREADS), "-DP", ref, ref],
                                 capture_output=True, text=True, timeout=900).stdout
            buckets = defaultdict(int)
            longest = 0
            for line in out.splitlines():
                f = line.split("\t")
                if len(f) < 12:
                    continue
                qs, qe, ts = int(f[2]), int(f[3]), int(f[7])
                if f[0] == f[5] and abs(qs - ts) < 50:
                    continue
                L = qe - qs
                if L < 500:
                    continue
                longest = max(longest, L)
                key = ("0.5-1kb" if L < 1000 else "1-2kb" if L < 2000
                       else "2-5kb" if L < 5000 else ">=5kb")
                buckets[key] += 1
            rep_rows.append((s, buckets["0.5-1kb"], buckets["1-2kb"], buckets["2-5kb"],
                             buckets[">=5kb"], longest))
        except Exception:
            pass

        if n % 10 == 0 or n == len(ids):
            print(f"  {n}/{len(ids)} isolates", flush=True)

    with open(f"{R}/rrna_operons.tsv", "w") as fh:
        fh.write("isolate\treplicon\tstart\tend\tlength\tgenes_n\tgenes\tGC_operon\tGC_left2kb\tGC_right2kb\n")
        for r in op_rows:
            fh.write("\t".join(str(x) for x in r) + "\n")

    with open(f"{R}/rrna_summary.tsv", "w") as fh:
        fh.write("isolate\treplicons\ttotal_bp\toperons\tmedian_spacing_bp\tGC_genome\tGC_operons_mean\n")
        for r in sum_rows:
            fh.write("\t".join(str(x) for x in r) + "\n")

    with open(f"{R}/repeats.tsv", "w") as fh:
        fh.write("isolate\tn_0.5_1kb\tn_1_2kb\tn_2_5kb\tn_ge5kb\tlongest_bp\n")
        for r in rep_rows:
            fh.write("\t".join(str(x) for x in r) + "\n")

    # ---- summary to stdout -------------------------------------------------
    if sum_rows:
        ops = sorted(int(r[3]) for r in sum_rows)
        gcs = sorted(float(r[5]) for r in sum_rows)
        print(f"\n{len(sum_rows)} isolates")
        print(f"  rRNA operons per genome: median {ops[len(ops)//2]}, range {ops[0]}-{ops[-1]}")
        print(f"  genome GC: median {gcs[len(gcs)//2]:.2f}%, range {gcs[0]:.2f}-{gcs[-1]:.2f}%")
        if op_rows:
            L = sorted(r[4] for r in op_rows)
            og = sorted(float(r[7]) for r in op_rows)
            print(f"  operon length: median {L[len(L)//2]:,} bp, range {L[0]:,}-{L[-1]:,}")
            print(f"  operon GC: median {og[len(og)//2]:.1f}% vs genome {gcs[len(gcs)//2]:.2f}%")
    if rep_rows:
        tot = [sum(r[i] for r in rep_rows) / len(rep_rows) for i in (1, 2, 3, 4)]
        print(f"  repeat copies per genome (mean): "
              f"0.5-1kb {tot[0]:.0f}, 1-2kb {tot[1]:.0f}, 2-5kb {tot[2]:.0f}, >=5kb {tot[3]:.0f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
