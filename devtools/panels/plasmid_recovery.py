#!/usr/bin/env python3
"""Score plasmid recovery against the isolate's OWN closed reference.

Contiguity statistics say nothing about whether a plasmid was recovered: NGA50 is a
chromosome-scale number and a genome can score well while losing every plasmid. This
measures the replicons directly.

Per plasmid record in the reference:
  * covered      -- fraction of the plasmid covered by aligned assembly contigs
  * whole        -- recovered in ONE contig covering >= 95% of it
  * n_contigs    -- how many contigs it took

And, because TesserACT labels each contig with a replicon call in its NAME
(_chr / _plas / _plas_<n> / _unk / _circular), the call is scored against truth:
  * a contig aligning mostly to a reference plasmid should be called plasmid
  * a contig aligning mostly to the chromosome should be called chromosomal
Precision and recall are reported for the plasmid class. SPAdes makes no such call,
so that half is blank for it -- which is itself the point: the classification is a
TesserACT feature and has to be scored, not assumed.
"""
import argparse, os, re, subprocess, sys, tempfile


def read_fasta(path):
    name, seq = None, []
    op = open
    with op(path) as fh:
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

    refs = list(read_fasta(a.ref))
    plasmids = {n.split()[0]: len(s) for n, s in refs if "plasmid" in n.lower()}
    chroms = {n.split()[0]: len(s) for n, s in refs if "plasmid" not in n.lower()}
    if not plasmids:
        print(f"{a.isolate}\t{a.arm}\t0\t0\t0\t\t\t")   # no plasmids in this reference
        return

    with tempfile.TemporaryDirectory() as td:
        paf = os.path.join(td, "aln.paf")
        with open(paf, "w") as fh:
            subprocess.run(["minimap2", "-x", "asm5", "-t", str(a.threads), "--secondary=no",
                            a.ref, a.asm], stdout=fh, stderr=subprocess.DEVNULL, check=False)
        # per reference record: covered intervals; per query contig: best target
        cov = {r: [] for r in list(plasmids) + list(chroms)}
        best = {}
        for line in open(paf):
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            q, qlen, ts, tname, tlen, tstart, tend = f[0], int(f[1]), int(f[2]), f[5], int(f[6]), int(f[7]), int(f[8])
            nmatch, blk = int(f[9]), int(f[10])
            if blk == 0 or nmatch / blk < a.min_id:
                continue
            if tname in cov:
                cov[tname].append((tstart, tend))
            span = tend - tstart
            if q not in best or span > best[q][1]:
                best[q] = (tname, span, qlen)

    def merged(iv):
        if not iv:
            return 0
        iv = sorted(iv); tot = 0; cs, ce = iv[0]
        for s, e in iv[1:]:
            if s > ce:
                tot += ce - cs; cs, ce = s, e
            else:
                ce = max(ce, e)
        return tot + ce - cs

    n_pl = len(plasmids)
    rec = whole = 0
    for p, plen in plasmids.items():
        c = merged(cov[p]) / plen if plen else 0
        if c >= 0.5:
            rec += 1
        ctgs = [q for q, (t, _, _) in best.items() if t == p]
        if c >= 0.95 and len(ctgs) == 1:
            whole += 1

    # replicon-call accuracy, only meaningful for arms that make the call
    tp = fp = fn = 0
    called_any = False
    for q, (t, _, _) in best.items():
        is_plasmid_truth = t in plasmids
        qlen = best[q][2]
        m = re.search(r"_(chr|plas|unk)(_\d+)?(_circular)?$", q)
        if qlen < a.min_classify:
            m = None
        if m:
            called_any = True
            called_plasmid = m.group(1) == "plas"
        else:
            continue
        if is_plasmid_truth and called_plasmid: tp += 1
        elif not is_plasmid_truth and called_plasmid: fp += 1
        elif is_plasmid_truth and not called_plasmid: fn += 1
    prec = f"{tp/(tp+fp):.3f}" if called_any and (tp + fp) else ""
    reca = f"{tp/(tp+fn):.3f}" if called_any and (tp + fn) else ""
    print(f"{a.isolate}\t{a.arm}\t{n_pl}\t{rec}\t{whole}\t{prec}\t{reca}\t{tp}/{fp}/{fn}")


if __name__ == "__main__":
    main()
