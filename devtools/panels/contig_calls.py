#!/usr/bin/env python3
"""Dump per-contig replicon calls with truth, so accuracy can be cut by contig length.

One row per contig: length, coverage, the call TesserACT made, and the replicon it
actually aligns to in the isolate's own closed reference. That makes the question
"above what length is the call trustworthy?" answerable from data instead of assumed.
"""
import argparse, os, re, subprocess, sys, tempfile


def read_fasta(path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:].rstrip("\n"), []
            else:
                seq.append(line.strip())
    if name is not None:
        yield name, "".join(seq)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", required=True)
    ap.add_argument("--asm", required=True)
    ap.add_argument("--arm", required=True)
    ap.add_argument("--isolate", required=True)
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--min-id", type=float, default=0.90)
    # Contigs shorter than this are scored as unclassified regardless of the label
    # the assembly carries. The in-binary floor was raised from 500 to 1500 bp partway
    # through this campaign, so assemblies on disk carry two conventions; re-gating here
    # makes the measurement independent of which binary produced the file.
    ap.add_argument("--min-classify", type=int, default=1500)
    a = ap.parse_args()

    plasmids, chroms = set(), set()
    for n, s in read_fasta(a.ref):
        (plasmids if "plasmid" in n.lower() else chroms).add(n.split()[0])

    with tempfile.TemporaryDirectory() as td:
        paf = os.path.join(td, "a.paf")
        with open(paf, "w") as fh:
            subprocess.run(["minimap2", "-x", "asm5", "-t", str(a.threads), "--secondary=no",
                            a.ref, a.asm], stdout=fh, stderr=subprocess.DEVNULL, check=False)
        best = {}
        for line in open(paf):
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            q, qlen, tname = f[0], int(f[1]), f[5]
            ts, te, nm, blk = int(f[7]), int(f[8]), int(f[9]), int(f[10])
            if blk == 0 or nm / blk < a.min_id:
                continue
            span = te - ts
            if q not in best or span > best[q][1]:
                best[q] = (tname, span, qlen)

    for name, seq in read_fasta(a.asm):
        q = name.split()[0]
        L = len(seq)
        mcov = re.search(r"_cov_([0-9.]+)", q)
        cov = mcov.group(1) if mcov else ""
        m = re.search(r"_(chr|plas|unk)(_\d+)?(_circular)?$", q)
        called = m.group(1) if m else "none"
        if L < a.min_classify and called in ("chr", "plas"):
            called = "unk"
        t = best.get(q)
        if not t:
            truth = "unaligned"
        elif t[0] in plasmids:
            truth = "plas"
        elif t[0] in chroms:
            truth = "chr"
        else:
            truth = "other"
        print(f"{a.isolate}\t{a.arm}\t{q}\t{L}\t{cov}\t{called}\t{truth}")


if __name__ == "__main__":
    main()
