#!/usr/bin/env python3
"""Map QUAST breakpoints to exact terminal repeat sequence in emitted GFA paths.

References are used by QUAST only; this is a diagnostic, never an assembler input.
A terminal repeat is a graph segment above 1.6 times graph median depth. The
reported interval is validated against the actual scaffold FASTA before use.
Example:
  terminal_repeat_breakpoints.py --assembly DIR --quast QUAST_DIR --out FILE.tsv
"""
import argparse
import csv
import re
import statistics
from pathlib import Path


def fasta(path):
    result, name, seq = {}, None, []
    with path.open() as handle:
        for line in handle:
            if line.startswith('>'):
                if name is not None:
                    result[name] = ''.join(seq).upper()
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line.strip())
    if name is not None:
        result[name] = ''.join(seq).upper()
    return result


def reverse_complement(seq):
    return seq.translate(str.maketrans('ACGTN', 'TGCAN'))[::-1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--assembly', type=Path, required=True)
    parser.add_argument('--quast', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    segments, depths, paths, overlaps = {}, {}, {}, []
    with (args.assembly / 'assembly_graph.gfa').open() as handle:
        for line in handle:
            row = line.rstrip().split('\t')
            if row[0] == 'S':
                segments[row[1]] = row[2].upper()
                depths[row[1]] = next((float(x.split(':')[-1]) for x in row[3:]
                                      if x.startswith(('dp:f:', 'DP:f:'))), 0.0)
            elif row[0] == 'P':
                paths[row[1]] = row[2].split(',')
            elif row[0] == 'L' and re.fullmatch(r'\d+M', row[5]):
                overlaps.append(int(row[5][:-1]))
    if not overlaps:
        parser.error('GFA does not provide graph overlap')
    overlap = int(statistics.median(overlaps))
    k = overlap + 1
    covs = sorted(depths[n] for n, seq in segments.items() if len(seq) >= 2 * k)
    median = covs[len(covs) // 2] if covs else 0.0
    contigs = fasta(args.assembly / 'contigs.fasta')
    scaffold_file = args.assembly / 'scaffolds.fasta'
    scaffolds = fasta(scaffold_file if scaffold_file.exists() else args.assembly / 'contigs.fasta')

    def sequence(node):
        seq = segments[node[:-1]]
        return reverse_complement(seq) if node.endswith('-') else seq

    terminal_intervals = {}
    for name, walk in paths.items():
        if name not in scaffolds:
            continue
        scaffold = scaffolds[name]
        intervals = []
        for head in (True, False):
            ordered = walk if head else list(reversed(walk))
            run = []
            for node in ordered:
                if depths[node[:-1]] <= 1.6 * median:
                    break
                run.append(node)
            if not run:
                continue
            if not head:
                run.reverse()
            spelled = sequence(run[0]) + ''.join(sequence(x)[overlap:] for x in run[1:])
            if head and scaffold.startswith(spelled):
                intervals.append(('head', 0, len(spelled), run))
            elif not head and scaffold.endswith(spelled):
                intervals.append(('tail', len(scaffold) - len(spelled), len(scaffold), run))
        terminal_intervals[name] = intervals

    # Splitting scaffolds renames/reorders contigs. Match exact sequence, not NODE
    # numbers. Repeated contained contigs may have several placements; retain all.
    placements = {}
    for name, seq in contigs.items():
        hits = []
        for scaffold_name, scaffold_seq in scaffolds.items():
            offset = scaffold_seq.find(seq)
            if offset >= 0:
                hits.append((scaffold_name, offset))
        placements[name] = hits

    rows = []
    gff_files = sorted((args.quast / 'contigs_reports').glob('*.misassemblies.gff'))
    if not gff_files:
        parser.error('no QUAST *.misassemblies.gff file found')
    for gff in gff_files:
        for line in gff.read_text().splitlines():
            if not line or line.startswith('#'):
                continue
            fields = line.split('\t')
            name, start, end = fields[0], int(fields[3]), int(fields[4])
            if name not in contigs:
                continue
            matched = []
            for scaffold, offset in placements[name]:
                for side, lo, hi, run in terminal_intervals.get(scaffold, []):
                    # QUAST places a breakpoint on an alignment-block boundary,
                    # which can fall anywhere in the graph's shared k-1 bases.
                    if start - 1 + offset <= hi + overlap and end + offset >= lo - overlap:
                        matched.append((scaffold, side, lo - offset + 1, hi - offset, run))
            if not matched:
                matched = [('', '', '', '', [])]
            for scaffold, side, lo, hi, run in matched:
                rows.append(dict(contig=name, contig_length=len(contigs[name]),
                                 breakpoint_start=start, breakpoint_end=end,
                                 terminal_repeat=bool(run), scaffold=scaffold, side=side,
                                 repeat_start=lo, repeat_end=hi, graph_nodes=','.join(run),
                                 graph_depths=','.join(f'{depths[x[:-1]]:.3f}' for x in run),
                                 graph_median_depth=f'{median:.3f}', k=k,
                                 annotation=fields[8]))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    columns = ['contig', 'contig_length', 'breakpoint_start', 'breakpoint_end',
               'terminal_repeat', 'scaffold', 'side', 'repeat_start', 'repeat_end',
               'graph_nodes', 'graph_depths', 'graph_median_depth', 'k', 'annotation']
    with args.out.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    events = {(r['contig'], r['breakpoint_start'], r['breakpoint_end']) for r in rows}
    terminal = {(r['contig'], r['breakpoint_start'], r['breakpoint_end'])
                for r in rows if r['terminal_repeat']}
    print(f'{len(terminal)}/{len(events)} QUAST breakpoints overlap validated terminal repeat context')
    print(args.out)


if __name__ == '__main__':
    main()
