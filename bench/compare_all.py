#!/usr/bin/env python3
"""Every configuration against SPAdes, on the isolates where all are scored."""
import os, statistics, sys
R="/media/iowa/WD_BLACK/kle_bench"
TAGS=sys.argv[1:] or ["vanilla","vcp","combo","kmodel","mcp","spades"]
def read(p):
    if not os.path.exists(p): return None
    return {f[0]: f[-1] for f in (l.rstrip("\n").split("\t") for l in open(p)) if len(f)>=2}
def g(d,k,c=float):
    try: return c(d[k])
    except Exception: return None
ids=sorted(f[:-4] for f in os.listdir(f"{R}/refs") if f.endswith(".fna"))
rows={t:[] for t in TAGS}
common=[]
for s in ids:
    r={t:read(f"{R}/quast/{s}_{t}/report.tsv") for t in TAGS}
    if any(r[t] is None for t in TAGS): continue
    common.append(s)
    for t in TAGS: rows[t].append(r[t])
if not common:
    print("no isolate has all tags scored yet"); raise SystemExit(1)
print(f"{len(common)} isolates with every configuration scored\n")
hdr=f"{'metric':<32}" + "".join(t.rjust(11) for t in TAGS)
print(hdr); print("-"*len(hdr))
def line(label,key,cast,agg,fmt):
    out=f"{label:<32}"
    for t in TAGS:
        vals=[g(d,key,cast) for d in rows[t] if g(d,key,cast) is not None]
        v=(statistics.median(vals) if agg=="med" else sum(vals)) if vals else float("nan")
        out+=fmt.format(v).rjust(11)
    print(out)
line("median contig NGA50","NGA50",int,"med","{:,.0f}")
line("median contigs","# contigs",int,"med","{:,.0f}")
line("total misassemblies","# misassemblies",int,"sum","{:,.0f}")
line("median mismatches /100kbp","# mismatches per 100 kbp",float,"med","{:.2f}")
line("median genome fraction (%)","Genome fraction (%)",float,"med","{:.2f}")
sp=rows["spades"] if "spades" in TAGS else None
if sp:
    out=f"{'median ratio to SPAdes':<32}"
    for t in TAGS:
        rr=[g(a,"NGA50",int)/g(b,"NGA50",int) for a,b in zip(rows[t],sp) if g(b,"NGA50",int)]
        out+=(f"{statistics.median(rr):.2f}x").rjust(11)
    print(out)
    out=f"{'isolates leading SPAdes':<32}"
    for t in TAGS:
        out+=str(sum(1 for a,b in zip(rows[t],sp) if g(a,"NGA50",int)>g(b,"NGA50",int))).rjust(11)
    print(out)
