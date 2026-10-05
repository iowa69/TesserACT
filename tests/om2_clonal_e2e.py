#!/usr/bin/env python3
"""om2_clonal_e2e.py -- end-to-end test of the Organism Model 2.0 clonal stage (C4) on synthetic clonal panels.

A synthetic species: an ancestor chromosome with three copies of a 5 kb operon-like repeat at fixed loci (every
genome: positional) and a 1.3 kb IS element. The isolate's lineage carries the IS at two sites; far lineages carry it
elsewhere. A panel of 8 clonal relatives and 12 far genomes is built into a model with layout tracks by
tesseract-model; four isolates are simulated and assembled with the om2 stages and TESSERACT_OM2_CLONAL=1:

  clone        a clonemate of the lineage (every IS site shared with its relatives)
  novel_is     the clone with an extra IS no panel genome has   -> no bases may be written into that junction
  missing_is   the clone without an IS its relatives carry        -> no element may be asserted there
  rearranged   the clone with an inversion whose one breakpoint is an IS copy and the other unique sequence
               -> the relatives' order must not be imposed across the rearranged record
Every isolate also carries the lineage's 50 kb plasmid (two copies of a plasmid IS, 2x copy number); the clonal
relatives carry it in the plasmid database, far lineages carry another. With --nrp-builder the nearest-relative plasmid
sidecar is built and TESSERACT_OM2_CLONAL_NRP set: plasmid records are chained on the relatives' plasmids.
Each isolate is also assembled with --no-layout, so the clonal stage does all of the chromosome ordering.

Checks (on genome/ against each isolate's truth):
  - every junction row carries the clonal columns (I8) and the [om2-clonal] counter line is printed (I9)
  - closure.txt lists the nearest relatives; the nearest are the clonal lineage (never a far genome)
  - every CLONAL join is TRUE: both 300-bp flanks map uniquely to the truth in order and orientation, the asserted
    gap (N-run or fill) within max(1000, 10 %) of the true gap
  - no non-positional fill in novel_is at the novel site; missing_is carries no N-run or fill at the missing site
  - clone: the chromosome is one genome record covering >= 98 % of the truth
Exit 0 = all pass. Prints one line per check.
"""
import argparse
import gzip
import os
import random
import re
import subprocess
import sys

COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def rc(s):
    return s.translate(COMP)[::-1]


def rnd(n, r):
    return ''.join(r.choice('ACGT') for _ in range(n))


def mutate(s, rate, r):
    if rate <= 0:
        return s
    b = list(s)
    n = int(len(b) * rate)
    for p in r.sample(range(len(b)), n):
        b[p] = r.choice([c for c in 'ACGT' if c != b[p]])
    return ''.join(b)


def write_fa(path, name, seq):
    with open(path, 'w') as f:
        f.write('>' + name + '\n')
        for i in range(0, len(seq), 80):
            f.write(seq[i:i + 80] + '\n')


def reads(genome, out, depth, r, insert=350, sd=35, rl=150, err=0.001, plasmid=None, pdepth=0):
    mols = [(genome, depth)] + ([(plasmid, pdepth)] if plasmid else [])
    with gzip.open(out + '_1.fq.gz', 'wt', compresslevel=1) as f1, gzip.open(out + '_2.fq.gz', 'wt', compresslevel=1) as f2:
      i = 0
      for genome, depth in mols:
        L = len(genome)
        g2 = genome + genome[:2000]   # circular
        n = int(L * depth / (2 * rl))
        for _ in range(n):
            i += 1
            fl = max(rl + 20, int(r.gauss(insert, sd)))
            st = r.randrange(L)
            frag = g2[st:st + fl]
            if r.random() < 0.5:
                frag = rc(frag)
            a, b = frag[:rl], rc(frag)[:rl]

            def noisy(s):
                s = list(s)
                for j in range(len(s)):
                    if r.random() < err:
                        s[j] = r.choice([c for c in 'ACGT' if c != s[j]])
                return ''.join(s)
            q = 'I' * rl
            f1.write(f'@r{i}/1\n{noisy(a)}\n+\n{q}\n')
            f2.write(f'@r{i}/2\n{noisy(b)}\n+\n{q}\n')


def read_fasta(p):
    d, name, buf = {}, None, []
    for line in open(p):
        if line.startswith('>'):
            if name is not None:
                d[name] = ''.join(buf)
            name = line[1:].split()[0]
            buf = []
        else:
            buf.append(line.strip())
    if name is not None:
        d[name] = ''.join(buf)
    return d


class Truth:
    """Exact locator on circular truth replicons (both strands). find -> [(replicon, pos, strand)]."""

    def __init__(self, seqs, k=300):
        self.recs = [(s.upper(), (s + s[:k]).upper()) for s in (seqs if isinstance(seqs, list) else [seqs])]
        self.L = len(self.recs[0][0])
        self.k = k

    def length(self, i):
        return len(self.recs[i][0])

    def find(self, s):
        s = s.upper()
        hits = []
        for ri, (seq, g2) in enumerate(self.recs):
            for strand, q in (('+', s), ('-', rc(s))):
                i = g2.find(q)
                while i >= 0 and i < len(seq):
                    hits.append((ri, i, strand))
                    i = g2.find(q, i + 1)
        return hits


def junction_truth(truth, left, right, asserted):
    """left: genome bases ending at the junction, right: starting at it. TRUE/FALSE/UNSCORABLE."""
    k = 300
    if len(left) < k or len(right) < k:
        return 'UNSCORABLE', None
    a = truth.find(left[-k:])
    b = truth.find(right[:k])
    if len(a) != 1 or len(b) != 1:
        return 'UNSCORABLE', None
    (ra, pa, sa), (rb, pb, sb) = a[0], b[0]
    if sa != sb or ra != rb:
        return 'FALSE', None
    L = truth.length(ra)
    if sa == '+':
        gap = (pb - (pa + k)) % L
    else:
        # on the minus strand the left flank's rc lies to the RIGHT of the right flank's rc
        gap = (pa - (pb + k)) % L
    if gap > L // 2:
        gap -= L
    tol = max(1000, abs(gap) // 10)
    return ('TRUE' if abs(gap - asserted) <= tol else 'FALSE'), gap


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--asm', required=True)
    ap.add_argument('--model', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--threads', default='2')
    ap.add_argument('--invariants', default='')
    ap.add_argument('--keep', action='store_true')
    ap.add_argument('--nrp-builder', default='')
    ap.add_argument('--extra-env', default='', help='round 3b: KEY=VALUE[,KEY=VALUE...] added to the assembler environment')
    a = ap.parse_args()
    r = random.Random(20260928)
    os.makedirs(a.out, exist_ok=True)
    fails = 0

    def check(ok, what, detail=''):
        nonlocal fails
        print(('PASS  ' if ok else 'FAIL  ') + what + ('  ' + detail if detail else ''))
        if not ok:
            fails += 1

    # ---- species --------------------------------------------------------------------------------------------
    anc = rnd(260000, r)
    R = rnd(5000, r)
    IS = rnd(1300, r)
    OPS = [40000, 125000, 205000]
    S1, S2, NOVEL, M = 70000, 165000, 100000, 92000
    FAR_IS = [20000, 150000, 235000]

    def build(is_sites, ops=OPS, extra=None, drop=None, invert=None, mut=0.0):
        """Insert operons and IS copies (ancestor coordinates), right to left; then mutate."""
        ins = [(p, 'R') for p in ops] + [(p, 'IS') for p in is_sites if p != drop]
        if extra is not None:
            ins.append((extra, 'IS'))
        s = anc
        if invert is not None:
            x, y = invert
            s = s[:x] + rc(s[x:y]) + s[y:]
        for p, what in sorted(ins, key=lambda t: -t[0]):
            s = s[:p] + (mutate(R, 0.001, random.Random(p)) if what == 'R' else IS) + s[p:]
        return mutate(s, mut, r)

    lineage = build([S1, S2])
    panel = os.path.join(a.out, 'panel')
    os.makedirs(panel, exist_ok=True)
    # the lineage plasmid: 50 kb with two copies of a plasmid IS; far lineages carry another plasmid
    PIS = rnd(1200, r)
    pb = rnd(50000, r)
    plasmid = pb[:15000] + PIS + pb[15000:35000] + PIS + pb[35000:]
    other = rnd(45000, r)
    names = []
    pdb = open(os.path.join(a.out, 'plasmids.fasta'), 'w')
    pmap = open(os.path.join(a.out, 'plasmid_map.tsv'), 'w')
    pmap.write('safe_acc\taccession\tplasmid_record\tlength\n')
    for i in range(8):
        nm = f'CLONE{i + 1}'
        write_fa(os.path.join(panel, nm + '_chr.fasta'), nm, mutate(lineage, 0.0002 * (i + 1), r))
        pl = mutate(plasmid, 0.0002 * (i + 1), r)
        pdb.write(f'>PL_{nm}.1\n{pl}\n')
        pmap.write(f'{nm}\t{nm}\tNZ_PL_{nm}.1\t{len(pl)}\n')   # RefSeq twin name: matched after normalisation
        names.append(nm)
    for i in range(12):
        nm = f'FAR{i + 1}'
        write_fa(os.path.join(panel, nm + '_chr.fasta'), nm, build(FAR_IS, mut=0.012))
        pl = mutate(other, 0.01, r)
        pdb.write(f'>PL_{nm}.1\n{pl}\n')
        pmap.write(f'{nm}\t{nm}\tPL_{nm}.1\t{len(pl)}\n')
        names.append(nm)
    pdb.close()
    pmap.close()
    tsm = os.path.join(a.out, 'synth.tsm')
    rc_ = subprocess.run([a.model, '--organism', 'synth', '--out', tsm, '--layout-tracks', '--min-support', '2',
                          '--marker-density', '64', '--plasmids', os.path.join(a.out, 'plasmids.fasta')] +
                         [os.path.join(panel, n + '_chr.fasta') for n in names],
                         capture_output=True, text=True)
    check(rc_.returncode == 0 and os.path.exists(tsm), 'model built with layout tracks', rc_.stderr.strip().splitlines()[-1] if rc_.stderr else '')
    if fails:
        return 1
    nrp = ''
    if a.nrp_builder:
        nrp = os.path.join(a.out, 'synth.om2nrp')
        p0 = subprocess.run([a.nrp_builder, tsm, os.path.join(a.out, 'plasmid_map.tsv'),
                             os.path.join(a.out, 'plasmids.fasta'), nrp, '--scenario', 'synthetic'], capture_output=True, text=True)
        txt = open(nrp).read() if os.path.exists(nrp) else ''
        check(p0.returncode == 0 and '#records\t20' in txt, 'plasmid sidecar: 20 relatives\' plasmids', p0.stderr.strip().splitlines()[-1] if p0.stderr else '')

    # ---- isolates ---------------------------------------------------------------------------------------------
    # The isolates descend from CLONE1's lineage genome with a clonal amount of change (0.0001).
    isolates = {
        'clone': build([S1, S2], mut=0.0001),
        'novel_is': build([S1, S2], extra=NOVEL, mut=0.0001),
        'missing_is': build([S1, S2], drop=S2, mut=0.0001),
        'rearranged': build([S1, S2], invert=(S1, M), mut=0.0001),
    }
    env = dict(os.environ)
    for k in list(env):
        if k.startswith('TESSERACT_'):
            del env[k]
    env.update({
        'TESSERACT_MODEL_AUTHOR': '1',
        'TESSERACT_OM2_SEAM': 'act', 'TESSERACT_OM2_CLOSE': '1', 'TESSERACT_OM2_FILL': 'contig',
        'TESSERACT_OM2_ALLOC': 'phased', 'TESSERACT_OM2_CIRC': '1', 'TESSERACT_OM2_OUTPUT': '1',
        'TESSERACT_OM2_AGP_EVIDENCE': '1', 'TESSERACT_NO_PLASMID_VOUCH': '1',
        'TESSERACT_OM2_CLONAL': '1',
    })
    if nrp:
        env['TESSERACT_OM2_CLONAL_NRP'] = nrp
    # round 3b: extra settings (e.g. the candidate's round-3 flags), KEY=VALUE[,KEY=VALUE...]
    for kv in filter(None, (a.extra_env or '').split(',')):
        k, v = kv.split('=', 1)
        env[k] = v
    iso_plasmid = mutate(plasmid, 0.0001, r)
    runs = []
    for iso, g in isolates.items():
        runs.append((iso, g, []))
        runs.append((iso + '_nolayout', g, ['--no-layout']))   # the clonal stage does all of the ordering
    for iso, g, extra in runs:
        d = os.path.join(a.out, iso)
        os.makedirs(d, exist_ok=True)
        write_fa(os.path.join(d, 'truth.fasta'), iso, g)
        base_iso = iso.replace('_nolayout', '')
        reads(g, os.path.join(d, 'r'), 60, random.Random(sum(map(ord, base_iso))), plasmid=iso_plasmid, pdepth=120)
        cmd = [a.asm, '--organism', 'synth', '--model', tsm, '-1', os.path.join(d, 'r_1.fq.gz'),
               '-2', os.path.join(d, 'r_2.fq.gz'), '-o', os.path.join(d, 'asm'), '-t', a.threads] + extra
        p = subprocess.run(cmd, env=env, capture_output=True, text=True)
        open(os.path.join(d, 'asm.log'), 'w').write(p.stderr)
        check(p.returncode == 0, f'{iso}: assembler exit 0', f'rc={p.returncode}')
        if p.returncode != 0:
            continue
        line = [l for l in p.stderr.splitlines() if l.startswith('[om2-clonal]')]
        check(len(line) == 1 and 'enabled=1' in line[0], f'{iso}: [om2-clonal] counter line (I9)',
              line[0][:220] if line else 'missing')
        gdir = os.path.join(d, 'asm', 'genome')
        clos = open(os.path.join(gdir, 'closure.txt')).read()
        near = re.findall(r'^#nearest\t(\d+)\t(\S+)\t', clos, re.M)
        check(len(near) >= 5 and all(n.startswith('CLONE') for _, n in near[:5]),
              f'{iso}: the 5 nearest relatives are the clonal lineage', ','.join(n for _, n in near[:6]))
        G = read_fasta(os.path.join(gdir, 'genome.fasta'))
        rows = [l.rstrip('\n').split('\t') for l in open(os.path.join(gdir, 'junctions.tsv'))]
        hdr, rows = rows[0], rows[1:]
        col = {c: i for i, c in enumerate(hdr)}
        check('clonal_class' in col and 'claim' in col and 'nr_agree' in col and 'graph_edge_a' in col
              and 'pair_contra' in col, f'{iso}: junctions.tsv clonal columns (I8)')
        truth = Truth([g, iso_plasmid])
        n_cl = n_true = n_false = n_uns = 0
        realized = [x for x in rows if x[col['genome_record']] not in ('.', '')]
        unsized_conf = 0
        for x in realized:
            if x[col['claim']] == 'confident' and x[col['clonal_class']] != '.' and x[col['gap_type']] == 'U':
                unsized_conf += 1
            if x[col['source']] != 'clonal':
                continue
            n_cl += 1
            rec = G[x[col['genome_record']]]
            s0, e0 = int(x[col['genome_start0']]), int(x[col['genome_end']])
            left, right = rec[:s0].replace('N', ''), rec[e0:]
            # the flanks must be N-free next to the junction
            li = rec.rfind('N', 0, s0)
            left = rec[li + 1:s0] if li >= 0 else rec[:s0]
            ri = rec.find('N', e0)
            right = rec[e0:ri] if ri >= 0 else rec[e0:]
            asserted = e0 - s0
            v, gap = junction_truth(truth, left, right, asserted)
            if v == 'TRUE':
                n_true += 1
            elif v == 'FALSE':
                n_false += 1
                print(f'      FALSE clonal join {x[col["junction"]]}: asserted {asserted} true_gap {gap} '
                      f'class {x[col["clonal_class"]]} claim {x[col["claim"]]} nr {x[col["nr_agree"]]}/{x[col["nr_placed"]]}')
            else:
                n_uns += 1
        check(n_false == 0, f'{iso}: every clonal join is TRUE', f'clonal={n_cl} true={n_true} false={n_false} unscorable={n_uns}')
        check(unsized_conf == 0, f'{iso}: no unsized junction claimed confident')
        chr_recs = {k: v for k, v in G.items() if k.startswith('chromosome_')}
        longest = max((len(v) for v in chr_recs.values()), default=0)
        print(f'      {iso}: chromosome records {len(chr_recs)} longest {longest} / truth {len(g)}; '
              f'clonal joins {n_cl}; {line[0][13:160] if line else ""}')
        # the plasmid: records whose first N-free 300-mer lies on the plasmid truth
        pl_len, pl_joins = 0, 0
        for nm, v in G.items():
            q = v.replace('N', '')[:300]
            h = truth.find(q) if len(q) == 300 else []
            if h and h[0][0] == 1:
                pl_len = max(pl_len, len(v))
        for x in realized:
            if x[col['source']] == 'clonal' and x[col['clonal_class']] == 'plasmid':
                pl_joins += 1
        print(f'      {iso}: plasmid longest record {pl_len} / {len(iso_plasmid)}; clonal plasmid joins {pl_joins}')
        if base_iso == 'clone':
            check(longest >= 0.98 * len(g), f'{iso}: chromosome in one genome record >= 98 %', f'{longest}/{len(g)}')
        if base_iso == 'novel_is':
            # the novel IS sits at NOVEL (ancestor coordinates) -> find its flanks in the truth and make sure no fill
            # of a non-positional clonal junction was written there
            bad = [x for x in realized if x[col['source']] == 'clonal' and x[col['clonal_class']] == 'non_positional'
                   and x[col['status']] == 'filled_genome' and int(x[col['nr_agree']]) < int(x[col['nr_placed']] or 0) // 2 + 1]
            check(not bad, 'novel_is: no non-positional fill without the relatives carrying it', f'{len(bad)} rows')
        if base_iso == 'missing_is':
            # flanks of the missing S2 site are adjacent in the truth; no genome record may assert >= 500 bases between them
            ok = True
            for x in realized:
                if x[col['source']] not in ('clonal',):
                    continue
                if x[col['clonal_class']] == 'non_positional' and x[col['nr_gapdiff']] not in ('0', '.') and int(x[col['n_len']] or 0) > 500:
                    rec = G[x[col['genome_record']]]
                    s0, e0 = int(x[col['genome_start0']]), int(x[col['genome_end']])
                    v, gap = junction_truth(truth, rec[max(0, s0 - 300):s0], rec[e0:e0 + 300], e0 - s0)
                    if v == 'FALSE':
                        ok = False
            check(ok, 'missing_is: no element asserted where the isolate lacks it')
        if a.invariants:
            p2 = subprocess.run([sys.executable, a.invariants, os.path.join(d, 'asm')], capture_output=True, text=True)
            check(p2.returncode == 0, f'{iso}: invariants.py I2/I3/I4', p2.stdout.strip().splitlines()[0] if p2.stdout else p2.stderr[-200:])
    print(f'om2_clonal_e2e: {"FAILED " + str(fails) if fails else "all passed"}')
    return 1 if fails else 0


if __name__ == '__main__':
    sys.exit(main())
