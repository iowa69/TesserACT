#!/usr/bin/env python3
"""Reconstruct vanilla resolver merge histories from true accepted-join traces.

Requires unchanged graph IDs and the default (unweighted, no-model) depth rule.
Sequence searches map connector intervals to final contigs without relying on
post-splitting NODE numbering. Missing exact matches remain explicitly unmapped.
"""
import argparse
import json
import re
from pathlib import Path
from terminal_repeat_breakpoints import fasta, reverse_complement

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--assembly', type=Path, required=True)
p.add_argument('--trace', type=Path, required=True)
p.add_argument('--quast', type=Path)
p.add_argument('--out', type=Path, required=True)
a = p.parse_args()
seqs, coverage, overlaps = {}, {}, []
for line in (a.assembly / 'assembly_graph.gfa').open():
    f = line.rstrip().split('\t')
    if f[0] == 'S':
        seqs[int(f[1])] = f[2]
        coverage[int(f[1])] = float(next(x for x in f[3:] if x.startswith('dp:f:')).split(':')[-1])
    elif f[0] == 'L' and f[5].endswith('M'):
        overlaps.append(int(f[5][:-1]))
ov = sorted(overlaps)[len(overlaps) // 2]
k = ov + 1
covs = sorted(coverage[u] for u, s in seqs.items() if len(s) >= 2 * k)
med = covs[len(covs) // 2]
anchors = [u for u in sorted(seqs) if len(seqs[u]) >= 2 * k and coverage[u] <= 1.6 * med]
chains = {i: [u * 2] for i, u in enumerate(anchors)}
trace = a.trace.read_text()
bounds = [int(x) for x in re.findall(r'maxPl=(\d+)', trace)]
reach = bounds[0] if bounds else 1000
contigs = fasta(a.assembly / 'contigs.fasta')
errors = []
if a.quast:
    gff_files = list((a.quast / 'contigs_reports').glob('*.misassemblies.gff'))
    if not gff_files:
        p.error('no detailed QUAST *.misassemblies.gff files found')
    for gff in gff_files:
        for line in gff.read_text().splitlines():
            if line.startswith('#') or not line:
                continue
            f = line.split('\t')
            errors.append((f[0], int(f[3]) - 1, int(f[4]), f[8]))

flip = lambda chain: [x ^ 1 for x in chain[::-1]]
def sequence(oid):
    s = seqs[oid // 2]
    return reverse_complement(s) if oid & 1 else s

results = []
for line in trace.splitlines():
    if not line.startswith('[acceptedjoin]'):
        continue
    v = dict(re.findall(r'(\w+)=([^\s]+)', line))
    ca, cb, ea, eb = (int(v[x]) for x in ('chainA', 'chainB', 'endA', 'endB'))
    route = list(map(int, v['route'].split(',')))
    left = chains[ca] if ea == 1 else flip(chains[ca])
    right = chains[cb] if eb == 0 else flip(chains[cb])
    if left[-1] != route[0] or right[0] != route[-1]:
        raise ValueError('trace/graph anchor mismatch; verify coverage mode and graph identity')
    repeats = []
    for side, history in [('A', left), ('B', flip(right))]:
        distance = 0
        for oid in reversed(history):
            if distance > reach:
                break
            if coverage[oid // 2] > 1.6 * med:
                repeats.append(dict(side=side, source=oid, distance=distance,
                                    length=len(seqs[oid // 2]), coverage=coverage[oid // 2]))
            distance += len(seqs[oid // 2]) - ov
    first = sequence(route[0])[-200:]
    middle = ''.join(sequence(oid)[ov:] for oid in route[1:-1])
    last = sequence(route[-1])[ov:ov + 200]
    window = first + middle + last
    lo = max(0, len(first) - ov)
    hi = len(first) + len(middle) + ov
    matches = []
    for name, contig in contigs.items():
        for reverse in (False, True):
            needle = reverse_complement(window) if reverse else window
            start = contig.find(needle)
            while start >= 0:
                jlo, jhi = (len(window) - hi, len(window) - lo) if reverse else (lo, hi)
                b, e = start + jlo, start + jhi
                overlap = [dict(start=x + 1, end=y, annotation=note)
                           for c, x, y, note in errors if c == name and x <= e and y >= b]
                matches.append(dict(contig=name, start=b + 1, end=e,
                                    reverse=reverse, breakpoints=overlap))
                start = contig.find(needle, start + 1)
    results.append(dict(v, repeat_history=repeats, matches=matches))
    chains[ca] = left + route[1:-1] + right
    chains[cb] = []
a.out.parent.mkdir(parents=True, exist_ok=True)
a.out.write_text(json.dumps(dict(k=k, median_depth=med, reach=reach, events=results), indent=2) + '\n')
print(json.dumps(dict(accepted=len(results), repeat_history=sum(bool(x['repeat_history']) for x in results),
                     mapped=sum(bool(x['matches']) for x in results),
                     breakpoint_events=sum(any(m['breakpoints'] for m in x['matches']) for x in results))))
