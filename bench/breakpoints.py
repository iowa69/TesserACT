#!/usr/bin/env python3
"""Why the chromosome is still in pieces.

Aligns an assembly to its closed chromosome, orders the alignments along the
reference, and reports every remaining break: how wide it is, what sequence
sits in it, and whether that sequence is a known repeat, an rRNA operon or an
insertion sequence. The reference is used only to *diagnose*; nothing here
feeds back into an assembly.

  breakpoints.py ISOLATE TAG [TAG2 ...]
"""
import os
import subprocess
import sys
from collections import Counter

R = "/media/iowa/WD_BLACK/kle_bench"
MM = "/home/iowa/miniconda3/envs/benchtools/bin/minimap2"


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
    g = sum(1 for c in s if c in "GCgc")
    n = sum(1 for c in s if c in "ACGTacgt")
    return 100.0 * g / n if n else 0.0


def repeat_intervals(ref_path):
    """Self-alignment intervals >=500 bp: the sequence that is duplicated."""
    out = subprocess.run([MM, "-x", "asm10", "-t", "6", "-DP", ref_path, ref_path],
                         capture_output=True, text=True).stdout
    iv = []
    for line in out.splitlines():
        f = line.split("\t")
        if len(f) < 12:
            continue
        qs, qe, ts = int(f[2]), int(f[3]), int(f[7])
        if f[0] == f[5] and abs(qs - ts) < 50:
            continue
        if qe - qs >= 500:
            iv.append((f[0], qs, qe))
    return iv


def overlaps(iv, contig, a, b):
    for c, s, e in iv:
        if c == contig and not (e < a or s > b):
            return True
    return False


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        return 2
    iso, tags = sys.argv[1], sys.argv[2:]
    ref = f"{R}/refs/{iso}.fna"
    seqs = load_fasta(ref)
    chrom = max(seqs, key=lambda k: len(seqs[k]))
    clen = len(seqs[chrom])
    print(f"{iso}: chromosome {chrom} {clen:,} bp\n")

    iv = repeat_intervals(ref)

    for tag in tags:
        asm = f"{R}/asm/{iso}_{tag}.fasta"
        if not os.path.exists(asm):
            print(f"  {tag}: no assembly")
            continue
        out = subprocess.run([MM, "-x", "asm10", "-t", "6", ref, asm],
                             capture_output=True, text=True).stdout
        # Every alignment block, not one per contig: a long contig legitimately
        # aligns to several stretches of the chromosome, and keeping only its
        # best would report the rest as missing.
        blocks = []
        names = set()
        for line in out.splitlines():
            f = line.split("\t")
            if len(f) < 12 or f[5] != chrom:
                continue
            if int(f[9]) < 200:            # ignore incidental hits
                continue
            blocks.append((int(f[7]), int(f[8]), f[0]))
            names.add(f[0])
        blocks.sort()
        if not blocks:
            print(f"  {tag}: nothing aligned")
            continue

        covered = 0
        gaps = []
        prev_end = 0
        for s, e, name in blocks:
            if s > prev_end:
                gaps.append((prev_end, s))
            covered += max(0, e - max(s, prev_end))
            prev_end = max(prev_end, e)
        if prev_end < clen:
            gaps.append((prev_end, clen))

        real = [g for g in gaps if g[1] - g[0] >= 100]
        kinds = Counter()
        widths = []
        for a, b in real:
            widths.append(b - a)
            kinds["repeat" if overlaps(iv, chrom, a, b) else "unique"] += 1

        print(f"  {tag}: {len(names)} contigs / {len(blocks)} blocks on the chromosome, "
              f"{100.0*covered/clen:.2f}% covered, {len(real)} breaks >=100 bp")
        if widths:
            widths.sort()
            print(f"      break width median {widths[len(widths)//2]:,} bp, "
                  f"total {sum(widths):,} bp")
            print(f"      in a repeat: {kinds['repeat']}   in unique sequence: {kinds['unique']}")
            print(f"      widest: " + ", ".join(f"{w:,}" for w in widths[-5:]))
            for a, b in sorted(real, key=lambda g: -(g[1]-g[0]))[:6]:
                seg = seqs[chrom][a:b]
                print(f"        {a:>9,}-{b:<9,} {b-a:>7,} bp  GC {gc(seg):4.1f}%  "
                      f"{'repeat' if overlaps(iv, chrom, a, b) else 'unique'}")
        print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
