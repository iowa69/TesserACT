#!/usr/bin/env bash
#
# TesserACT end-to-end test suite.
#
# Everything is generated from scratch in a temporary directory: synthetic
# genomes and reads come from the embedded python helper, so the suite has no
# dependency on a read simulator or a reference data set. Run with
#   make test          (or)      bash tests/run_tests.sh
#
# Exits non-zero if any test fails.

set -u

# Start from release defaults: drop every ambient TESSERACT_* so an exported experiment
# environment (an arm's COMMON_PREFIX=0, MIN_FALLBACK_DEST=..., ...) can neither fail nor
# mask the checks below. A test that needs a flag sets it on that one command line.
while IFS= read -r v; do unset "$v"; done < <(compgen -e | grep '^TESSERACT_' || true)

# The model tests exercise the author path: building a model and assembling against
# it. --model only accepts a path under TESSERACT_MODEL_AUTHOR, because a model is
# a curated artifact and --organism selects a bundled one for everybody else.
export TESSERACT_MODEL_AUTHOR=1

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
TESSERACT=$ROOT/tesseract-asm

if [ ! -x "$TESSERACT" ]; then
    echo "error: $TESSERACT not built -- run 'make' first" >&2
    exit 1
fi
if ! command -v python3 >/dev/null 2>&1; then
    echo "error: python3 is required to generate the test fixtures" >&2
    exit 1
fi

TMP=$(mktemp -d "${TMPDIR:-/tmp}/tesseract-test.XXXXXXXX")
trap 'rm -rf "$TMP"' EXIT INT TERM

GEN="$TMP/gen.py"
LOG="$TMP/tesseract.log"

PASSED=0
FAILED=0
SKIPPED=0
XFAILED=0

pass() { printf 'PASS  %-46s %s\n' "$1" "${2:-}"; PASSED=$((PASSED + 1)); }
fail() { printf 'FAIL  %-46s %s\n' "$1" "${2:-}"; FAILED=$((FAILED + 1)); }
# A test that could not run is reported as such, never as a pass.
skip() { printf 'SKIP  %-46s %s\n' "$1" "${2:-}"; SKIPPED=$((SKIPPED + 1)); }

# known NAME DEFECT CONDITION_EXIT_STATUS DETAIL -- a check of the corrected behaviour for a
# tracked, still-open defect. It passes once the fix is in the binary (fixes are default-off
# and enabled together by TESSERACT_FIXES=1, which the caller sets on that command line);
# until then it prints XFAIL and does not fail the suite. RUN_TESTS_STRICT_XFAIL=1 makes an
# XFAIL a failure, for the build in which every fix is expected to have landed.
known() {
    if [ "$3" -eq 0 ]; then
        pass "$1" "[$2 fixed] ${4:-}"
    elif [ "${RUN_TESTS_STRICT_XFAIL:-0}" = 1 ]; then
        fail "$1" "[$2 open] ${4:-}"
    else
        printf 'XFAIL %-46s %s\n' "$1" "[$2 open] ${4:-}"; XFAILED=$((XFAILED + 1))
    fi
}

# same_outputs DIR_A DIR_B -- prints the written outputs that differ between two runs (by
# presence or by a single byte), nothing when they agree. Every file a run writes that
# should be reproducible, not only contigs.fasta: a thread-dependent scaffold, AGP or graph
# tie-break leaves contigs.fasta alone.
same_outputs() {
    local f diff=""
    for f in contigs.fasta scaffolds.fasta scaffolds.agp assembly_graph.gfa; do
        if [ -e "$1/$f" ] || [ -e "$2/$f" ]; then
            cmp -s "$1/$f" "$2/$f" || diff="$diff $f"
        fi
    done
    printf '%s' "$diff"
}
# present_outputs DIR -- which of those files a run wrote (so a pass shows what it covered).
present_outputs() {
    local f out=""
    for f in contigs.fasta scaffolds.fasta scaffolds.agp assembly_graph.gfa; do
        [ -e "$1/$f" ] && out="$out${out:+,}${f}"
    done
    printf '%s' "$out"
}

# check NAME CONDITION_EXIT_STATUS DETAIL
check() {
    if [ "$2" -eq 0 ]; then pass "$1" "${3:-}"; else fail "$1" "${3:-}"; fi
}

TIMEOUT=""
if command -v timeout >/dev/null 2>&1; then TIMEOUT="timeout 180"; fi

# asm OUTDIR [TesserACT args...] -- always quiet, stdout+stderr captured in $LOG
asm() {
    local out=$1
    shift
    rm -rf "$out"
    $TIMEOUT "$TESSERACT" "$@" -o "$out" -q >"$LOG" 2>&1
}

# absdiff A B
absdiff() { if [ "$1" -ge "$2" ]; then echo $(( $1 - $2 )); else echo $(( $2 - $1 )); fi; }

cat > "$GEN" <<'PYEOF'
#!/usr/bin/env python3
"""Fixture generator and assembly checker for the TesserACT test suite.

Subcommands print a one-line summary and exit non-zero when a check fails.
"""
import gzip
import json
import random
import re
import sys

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def rc(s):
    return s.translate(COMP)[::-1]


def opener(path):
    return gzip.open if path.endswith(".gz") else open


def write_fasta(path, records, width=70):
    with opener(path)(path, "wt") as fh:
        for name, seq in records:
            fh.write(">%s\n" % name)
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + "\n")


def read_fasta(path):
    recs, name, buf = [], None, []
    with opener(path)(path, "rt") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line[0] == ">":
                if name is not None:
                    recs.append((name, "".join(buf)))
                name, buf = line[1:], []
            else:
                buf.append(line)
    if name is not None:
        recs.append((name, "".join(buf)))
    return recs


def rand_seq(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def cmd_genome(a):
    """genome OUT SEED LEN [LEN...] -- unrelated random sequences, one per record."""
    out, seed, lens = a[0], int(a[1]), [int(x) for x in a[2:]]
    rng = random.Random(seed)
    write_fasta(out, [("chr%d" % (i + 1), rand_seq(rng, L)) for i, L in enumerate(lens)])
    print("lengths=%s" % ",".join(str(L) for L in lens))


def cmd_repeat_genome(a):
    """repeat_genome OUT SEED UNIQUE REPEAT -- layout A R B R C with R identical."""
    out, seed, uniq, rep = a[0], int(a[1]), int(a[2]), int(a[3])
    rng = random.Random(seed)
    A, B, C, R = (rand_seq(rng, uniq), rand_seq(rng, uniq),
                  rand_seq(rng, uniq), rand_seq(rng, rep))
    write_fasta(out, [("chr1", A + R + B + R + C)])
    print("length=%d" % (3 * uniq + 2 * rep))


def mutate(rng, seq, err):
    if err <= 0:
        return seq
    out = list(seq)
    for i, b in enumerate(out):
        if rng.random() < err:
            out[i] = rng.choice([x for x in "ACGT" if x != b])
    return "".join(out)


def cmd_reads(a):
    """reads GENOME PREFIX COV READLEN INSERT SD ERR SEED MODE FORMAT gz|plain"""
    (gpath, prefix, cov, rlen, ins, sd, err, seed, mode, fmt, gz) = (
        a[0], a[1], float(a[2]), int(a[3]), float(a[4]), float(a[5]),
        float(a[6]), int(a[7]), a[8], a[9], a[10])
    rng = random.Random(seed)
    pairs = []
    minfrag = rlen + 20
    for _, seq in read_fasta(gpath):
        L = len(seq)
        want = int(cov * L / (2.0 * rlen))
        made, guard = 0, 0
        while made < want and guard < want * 40:
            guard += 1
            frag = int(rng.gauss(ins, sd))
            if frag < minfrag:
                continue
            # Fragment starts may hang off either end and are then clipped, so
            # the terminal bases get interior-like depth. Without this the first
            # and last few hundred bases are covered by almost no fragments and
            # contig ends are ragged for reasons that have nothing to do with
            # the assembler.
            s = rng.randint(-(frag - 1), L - 1)
            lo, hi = max(0, s), min(L, s + frag)
            if hi - lo < minfrag:
                continue
            f = seq[lo:hi]
            pairs.append((mutate(rng, f[:rlen], err), mutate(rng, rc(f[-rlen:]), err)))
            made += 1
    rng.shuffle(pairs)

    ext = ("fa" if fmt == "fasta" else "fq") + (".gz" if gz == "gz" else "")
    op = gzip.open if gz == "gz" else open

    def emit(fh, name, seq):
        if fmt == "fasta":
            fh.write(">%s\n%s\n" % (name, seq))
        else:
            fh.write("@%s\n%s\n+\n%s\n" % (name, seq, "I" * len(seq)))

    if mode == "paired":
        with op("%s_1.%s" % (prefix, ext), "wt") as f1, \
             op("%s_2.%s" % (prefix, ext), "wt") as f2:
            for i, (r1, r2) in enumerate(pairs):
                emit(f1, "r%d/1" % i, r1)
                emit(f2, "r%d/2" % i, r2)
    elif mode == "interleaved":
        with op("%s_12.%s" % (prefix, ext), "wt") as fh:
            for i, (r1, r2) in enumerate(pairs):
                emit(fh, "r%d/1" % i, r1)
                emit(fh, "r%d/2" % i, r2)
    else:
        with op("%s_s.%s" % (prefix, ext), "wt") as fh:
            for i, (r1, r2) in enumerate(pairs):
                emit(fh, "r%da" % i, r1)
                emit(fh, "r%db" % i, r2)
    print("pairs=%d" % len(pairs))


def cmd_tiny(a):
    """tiny OUT -- a handful of short reads, far too little to assemble."""
    rng = random.Random(99)
    with open(a[0], "wt") as fh:
        for i in range(3):
            s = rand_seq(rng, 60)
            fh.write("@t%d\n%s\n+\n%s\n" % (i, s, "I" * len(s)))
    print("reads=3")


def contigs_of(path):
    return [s for _, s in read_fasta(path)]


def cmd_stats(a):
    """stats CONTIGS -- shell-evalable summary."""
    lens = sorted((len(s) for s in contigs_of(a[0])), reverse=True)
    total = sum(lens)
    acc, n50 = 0, 0
    for L in lens:
        acc += L
        if acc * 2 >= total:
            n50 = L
            break
    print("n=%d total=%d largest=%d n50=%d" % (len(lens), total, lens[0] if lens else 0, n50))


def cmd_exact(a):
    """exact CONTIGS GENOME -- one contig per source record, identical up to rc."""
    seqs, src = contigs_of(a[0]), [s for _, s in read_fasta(a[1])]
    used = set()
    for g in src:
        for j, c in enumerate(seqs):
            if j not in used and (c == g or rc(c) == g):
                used.add(j)
                break
    print("contigs=%d sources=%d exact=%d" % (len(seqs), len(src), len(used)))
    if len(used) != len(src) or len(seqs) != len(src):
        sys.exit(1)


def pieces_of(seq):
    """Scaffold gaps are runs of N; each side of a gap is checked separately."""
    return [p for p in re.split("N+", seq) if p]


def cmd_substr(a):
    """substr CONTIGS GENOME -- every contig is a contiguous piece of the source.

    This is the misassembly check: a contig that joins the wrong flanks of a
    repeat is built entirely from real source k-mers but is not a substring of
    the source anywhere.
    """
    seqs = contigs_of(a[0])
    src = "".join(s for _, s in read_fasta(a[1]))
    both = src + "\x00" + rc(src)
    bad = sum(1 for c in seqs for p in pieces_of(c) if p not in both)
    gaps = sum(c.count("N") for c in seqs)
    print("contigs=%d off_reference=%d gap_bases=%d" % (len(seqs), bad, gaps))
    if bad:
        sys.exit(1)


def kmerset(seq, k):
    return set(seq[i:i + k] for i in range(len(seq) - k + 1))


def cmd_kmercheck(a):
    """kmercheck CONTIGS GENOME K -- no contig k-mer may be absent from the source."""
    k = int(a[2])
    src = "".join(s for _, s in read_fasta(a[1]))
    ref = kmerset(src, k) | kmerset(rc(src), k)
    tot = bad = 0
    for c in contigs_of(a[0]):
        for p in pieces_of(c):
            for i in range(len(p) - k + 1):
                tot += 1
                if p[i:i + k] not in ref:
                    bad += 1
    frac = 1.0 if tot == 0 else 1.0 - bad / float(tot)
    print("kmers=%d foreign=%d frac=%.6f" % (tot, bad, frac))
    if bad:
        sys.exit(1)


def cmd_identity(a):
    """identity CONTIGS GENOME MIN -- base identity of the best contig vs the source.

    The contig is anchored to the source with an exact 31-mer, then compared
    base by base along that diagonal; the denominator is the longer of the two,
    so missing or extra sequence counts against the score.
    """
    k = 31
    src = [s for _, s in read_fasta(a[1])][0]
    idx = {}
    for i in range(len(src) - k + 1):
        idx.setdefault(src[i:i + k], i)
    best = 0.0
    for c in contigs_of(a[0]):
        for cand in (c, rc(c)):
            off = None
            for i in range(len(cand) - k + 1):
                p = idx.get(cand[i:i + k])
                if p is not None:
                    off = p - i
                    break
            if off is None:
                continue
            m = sum(1 for i, b in enumerate(cand)
                    if 0 <= i + off < len(src) and src[i + off] == b)
            best = max(best, m / float(max(len(cand), len(src))))
    print("identity=%.5f" % best)
    if best < float(a[2]):
        sys.exit(1)


def canon_set(path):
    return sorted(min(s, rc(s)) for s in contigs_of(path))


def cmd_sameset(a):
    """sameset A B -- identical contig multisets, orientation-insensitive."""
    x, y = canon_set(a[0]), canon_set(a[1])
    print("a=%d b=%d equal=%s" % (len(x), len(y), x == y))
    if x != y:
        sys.exit(1)


def cmd_format(a):
    """format CONTIGS -- header syntax, stated length, longest-first ordering."""
    names, lens = [], []
    with open(a[0]) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                names.append(line[1:])
                lens.append(0)
            elif line:
                lens[-1] += len(line)
    # The replicon tag is appended, so the rank/length/coverage contract is unchanged and
    # still checked; only a suffix is newly permitted. The permitted tags are enumerated
    # rather than left open, so a malformed or unexpected one still fails:
    #   _chr | _unk | _plas | _plas_<n>, each optionally followed by _circular
    pat = re.compile(r"^NODE_(\d+)_length_(\d+)_cov_(\d+\.\d+)"
                     r"(?:_(?:chr|unk|plas(?:_\d+)?)(?:_circular)?)?$")
    bad_name = bad_len = bad_rank = 0
    for i, (nm, L) in enumerate(zip(names, lens)):
        m = pat.match(nm)
        if not m:
            bad_name += 1
            continue
        if int(m.group(2)) != L:
            bad_len += 1
        if int(m.group(1)) != i + 1:
            bad_rank += 1
    # Length descending WITHIN a replicon block, not across the whole file. The output is
    # ordered as a genome -- chromosome, then each plasmid molecule contiguous, then what
    # could not be called -- so a chromosomal fragment legitimately precedes a longer plasmid.
    # With no model every contig lands in one block and this reduces to the original global
    # check, which is the case this suite exercises; the per-block form is the real contract.
    tags = [re.search(r"_(?:chr|unk|plas(?:_\d+)?)(?:_circular)?$", nm) for nm in names]
    blocks = [m.group(0).replace("_circular", "") if m else "" for m in tags]
    unsorted = sum(1 for i in range(1, len(lens))
                   if blocks[i] == blocks[i - 1] and lens[i] > lens[i - 1])
    print("contigs=%d bad_header=%d bad_length=%d bad_rank=%d unsorted=%d"
          % (len(names), bad_name, bad_len, bad_rank, unsorted))
    if not names or bad_name or bad_len or bad_rank or unsorted:
        sys.exit(1)


def n50_of(lens):
    total = sum(lens)
    acc = 0
    for l in sorted(lens, reverse=True):
        acc += l
        if acc * 2 >= total:
            return l
    return 0


def contigstats_verdict(out, report_path=None):
    """(detail line, list of mismatches) for one output directory; see cmd_contigstats."""
    import os
    rep = json.load(open(report_path or os.path.join(out, "report.json")))["assembly"]
    recs = [s for _, s in read_fasta(os.path.join(out, "contigs.fasta"))]
    lens = [len(s) for s in recs]
    n_in_contigs = sum(s.upper().count("N") for s in recs)
    sp = os.path.join(out, "scaffolds.fasta")
    runs = [len(m.group(0)) for _, s in (read_fasta(sp) if os.path.exists(sp) else [])
            for m in re.finditer("[Nn]+", s)]
    agp = os.path.join(out, "scaffolds.agp")
    agp_n = 0
    if os.path.exists(agp):
        with open(agp) as fh:
            agp_n = sum(1 for l in fh if not l.startswith("#") and l.split("\t")[4:5] == ["N"])
    n50 = n50_of(lens)
    want = {"scaffold_gaps": len(runs), "contig_count": len(lens),
            "contig_total_length": sum(lens), "contig_n50": n50,
            "contig_largest": max(lens, default=0)}
    bad = [k for k, v in want.items() if rep.get(k) != v]
    notes = ["  %s: report says %s, the files say %s" % (k, rep.get(k), want[k]) for k in bad]
    if n_in_contigs:
        bad.append("N_in_contigs.fasta")
        notes.append("  contigs.fasta carries %d N bases" % n_in_contigs)
    if agp_n != len(runs):
        bad.append("agp_N_rows")
        notes.append("  scaffolds.agp has %d N rows for %d N runs" % (agp_n, len(runs)))
    detail = ("gaps=%d short_gaps=%d contigs=%d contig_n50=%d scaffold_n50=%d mismatched=%s"
              % (len(runs), sum(1 for r in runs if r < 10), len(lens), n50, rep.get("n50", 0),
                 ",".join(bad) if bad else "none"))
    return detail, bad, notes


def cmd_contigstats(a):
    """contigstats OUTDIR [REPORT_JSON] -- report.json describes the files actually written.

    The scaffold figures (n50, largest) are raised by asserting an order across a gap, which
    assembles no additional base; the contig figures are what QUAST and every published
    comparison report. Quoting the first as if it were the second has already misled a
    comparison in this project, so the invariant is pinned here against the files:
      contig_count / contig_total_length / contig_n50 / contig_largest == the records of
        contigs.fasta, taken as written. Since 1.3.0 that file is split at every N run and
        then trimmed, so it is never re-split here: the earlier form of this check split it
        at runs of 10+ N, found none, and so compared nothing (and failed every correct
        report of a gapped assembly, which is why no fixture had a gap).
      scaffold_gaps == the number of N runs of any length in scaffolds.fasta (0 without it);
      contigs.fasta carries no N, and scaffolds.agp has one N row per N run.
    The detail line prints gaps=, so a run that produced no gap is visible as gaps=0.
    """
    detail, bad, notes = contigstats_verdict(a[0], a[1] if len(a) > 1 else None)
    print(detail)
    if bad:
        for n in notes:
            print(n)
        sys.exit(1)


def cmd_contigstats_selftest(a):
    """contigstats_selftest DIR -- the contigstats check itself, on hand-written output dirs.

    Every case is a directory a run could write, with a report that either describes it
    (the check must pass) or is stale in the way report.json has been stale (it must fail):
    a 1-N gap counted by a >=10-N rule, a gap the polisher overwrote, a terminal trim that
    the totals still include, an N left in contigs.fasta, an AGP that disagrees.
    """
    import os
    rng = random.Random(4242)
    A, B, C = rand_seq(rng, 1200), rand_seq(rng, 900), rand_seq(rng, 700)

    def stats(lens, gaps):
        return {"assembly": {"contig_count": len(lens), "contig_total_length": sum(lens),
                             "contig_n50": n50_of(lens), "contig_largest": max(lens, default=0),
                             "scaffold_gaps": gaps, "n50": n50_of(lens)}}

    def agp_rows(scaffolds):
        rows = []
        for name, s in scaffolds:
            pos, part = 1, 1
            for m in re.finditer("[ACGT]+|N+", s):
                L = len(m.group(0))
                if m.group(0)[0] == "N":
                    rows.append("%s\t%d\t%d\t%d\tN\t%d\tscaffold\tyes\tpaired-ends"
                                % (name, pos, pos + L - 1, part, L))
                else:
                    rows.append("%s\t%d\t%d\t%d\tW\t%s_%d\t1\t%d\t+"
                                % (name, pos, pos + L - 1, part, name, part, L))
                pos += L
                part += 1
        return rows

    def case(name, contigs, scaffolds, report, agp=None):
        d = os.path.join(a[0], name)
        os.makedirs(d, exist_ok=True)
        write_fasta(os.path.join(d, "contigs.fasta"),
                    [("NODE_%d" % (i + 1), s) for i, s in enumerate(contigs)])
        if scaffolds is not None:
            recs = [("SCAF_%d" % (i + 1), s) for i, s in enumerate(scaffolds)]
            write_fasta(os.path.join(d, "scaffolds.fasta"), recs)
            with open(os.path.join(d, "scaffolds.agp"), "w") as fh:
                fh.write("##agp-version\t2.1\n")
                for r in (agp if agp is not None else agp_rows(recs)):
                    fh.write(r + "\n")
        json.dump(report, open(os.path.join(d, "report.json"), "w"))
        return d

    L = lambda *xs: [len(x) for x in xs]
    trimmed = C[:-30]
    cases = [
        # name, contigs.fasta, scaffolds.fasta, report, expected verdict (0 = consistent)
        ("no_gap", [A, B, C], None, stats(L(A, B, C), 0), 0),
        ("gap_12N", [A, B, C], [A + "N" * 12 + B, C], stats(L(A, B, C), 1), 0),
        ("gap_1N", [A, B, C], [A + "N" + B, C], stats(L(A, B, C), 1), 0),
        # The release report: pre-split scaffolds under a >=10-N rule, so the 1-N join is
        # one contig of len(A)+1+len(B) and no gap.
        ("gap_1N_stale_10N_rule", [A, B, C], [A + "N" + B, C],
         stats([len(A) + 1 + len(B), len(C)], 0), 1),
        # Contig figures right, but the short run not counted as a scaffold gap.
        ("gap_2N_stale_gap_count", [A, B, C], [A + "NN" + B, C], stats(L(A, B, C), 0), 1),
        # Polisher overwrote the gap with read bases: one record, no N run anywhere.
        ("polished_over_N", [A + "G" * 5 + B, C], [A + "G" * 5 + B, C],
         stats(L(A + "G" * 5 + B, C), 0), 0),
        ("polished_over_N_stale", [A + "G" * 5 + B, C], [A + "G" * 5 + B, C],
         stats(L(A, B, C), 1), 1),
        # Terminal-overlap trim removed 30 bp after the statistics were taken.
        ("trimmed", [A, B, trimmed], None, stats(L(A, B, trimmed), 0), 0),
        ("trimmed_stale", [A, B, trimmed], None, stats(L(A, B, C), 0), 1),
        ("N_in_contigs", [A + "N" * 12 + B, C], None, stats(L(A + "N" * 12 + B, C), 0), 1),
        ("agp_missing_N_row", [A, B, C], [A + "N" * 12 + B, C], stats(L(A, B, C), 1), 1),
    ]
    wrong = []
    for name, contigs, scaffolds, report, expect in cases:
        agp = None
        if name == "agp_missing_N_row":
            agp = [r for r in agp_rows([("SCAF_1", scaffolds[0]), ("SCAF_2", scaffolds[1])])
                   if "\tN\t" not in r]
        d = case(name, contigs, scaffolds, report, agp)
        _, bad, _ = contigstats_verdict(d)
        got = 1 if bad else 0
        if got != expect:
            wrong.append("%s(expected %s, got %s)" % (name, "fail" if expect else "pass",
                                                       "fail" if got else "pass"))
    print("cases=%d misjudged=%d%s" % (len(cases), len(wrong),
                                        (" " + ",".join(wrong)) if wrong else ""))
    if wrong:
        sys.exit(1)


def cmd_holereads(a):
    """holereads GENOME PREFIX COV READLEN INSERT SD SEED HOLE_LO HOLE_HI

    Error-free pairs as `reads ... paired fastq gz`, except that no READ overlaps
    [HOLE_LO, HOLE_HI): fragments may still span it, with both mates outside. A coverage
    dropout that pairs bridge is the one situation that makes the assembler write a scaffold
    gap, which none of the other fixtures produce.
    """
    gpath, prefix = a[0], a[1]
    cov, rlen, ins, sd = float(a[2]), int(a[3]), float(a[4]), float(a[5])
    seed, hlo, hhi = int(a[6]), int(a[7]), int(a[8])
    rng = random.Random(seed)
    pairs = []
    for _, seq in read_fasta(gpath):
        L = len(seq)
        want = int(cov * L / (2.0 * rlen))
        made = guard = 0
        while made < want and guard < want * 40:
            guard += 1
            frag = int(rng.gauss(ins, sd))
            if frag < rlen + 20:
                continue
            s = rng.randint(-(frag - 1), L - 1)
            lo, hi = max(0, s), min(L, s + frag)
            if hi - lo < rlen + 20:
                continue
            if any(x < hhi and y > hlo for x, y in ((lo, lo + rlen), (hi - rlen, hi))):
                continue
            f = seq[lo:hi]
            pairs.append((f[:rlen], rc(f[-rlen:])))
            made += 1
    rng.shuffle(pairs)
    with gzip.open(prefix + "_1.fq.gz", "wt") as f1, gzip.open(prefix + "_2.fq.gz", "wt") as f2:
        for i, (r1, r2) in enumerate(pairs):
            f1.write("@r%d/1\n%s\n+\n%s\n" % (i, r1, "I" * len(r1)))
            f2.write("@r%d/2\n%s\n+\n%s\n" % (i, r2, "I" * len(r2)))
    print("pairs=%d" % len(pairs))


def cmd_truncgz(a):
    """truncgz IN_FQ_GZ OUT KEEP_RECORDS -- a gzip stream cut after KEEP_RECORDS records.

    The deflate stream is flushed at a record boundary and then simply stops: no final
    block, no CRC/ISIZE trailer -- a download that died mid-file. `gzip -t` rejects it, and
    everything it does contain decompresses to whole FASTQ records.
    """
    import zlib
    with gzip.open(a[0], "rt") as fh:
        lines = fh.read().split("\n")
    keep = int(a[2]) * 4
    body = ("\n".join(lines[:keep]) + "\n").encode()
    z = zlib.compressobj(6, zlib.DEFLATED, 31)
    data = z.compress(body) + z.flush(zlib.Z_SYNC_FLUSH)
    with open(a[1], "wb") as out:
        out.write(data)
    print("records=%d bytes=%d" % (keep // 4, len(data)))


def cmd_dotnames(a):
    """dotnames IN_FQ OUT_FQ MATE -- rename '@rN/M' to '@SRR1.N.M' (fastq-dump --readids)."""
    with open(a[0]) as src, open(a[1], "w") as dst:
        for i, line in enumerate(src):
            if i % 4 == 0:
                line = re.sub(r"^@r(\d+)/\d", r"@SRR1.\1." + a[2], line)
            dst.write(line)


if __name__ == "__main__":
    table = {"genome": cmd_genome, "repeat_genome": cmd_repeat_genome, "reads": cmd_reads,
             "tiny": cmd_tiny, "stats": cmd_stats, "exact": cmd_exact, "substr": cmd_substr,
             "kmercheck": cmd_kmercheck, "identity": cmd_identity, "sameset": cmd_sameset,
             "format": cmd_format, "contigstats": cmd_contigstats,
             "contigstats_selftest": cmd_contigstats_selftest, "holereads": cmd_holereads,
             "truncgz": cmd_truncgz, "dotnames": cmd_dotnames}
    table[sys.argv[1]](sys.argv[2:])
PYEOF

gen() { python3 "$GEN" "$@"; }

echo "TesserACT test suite"
echo "  binary   $TESSERACT"
echo "  workdir  $TMP"
echo

# ---------------------------------------------------------------------------
# 1. Trivial reconstruction: error-free pairs from a 20 kb sequence
# ---------------------------------------------------------------------------
D=$TMP/t1; mkdir -p "$D"
gen genome "$D/g.fa" 1 20000 >/dev/null
gen reads "$D/g.fa" "$D/r" 50 150 350 30 0 1001 paired fastq gz >/dev/null
if asm "$D/out" -1 "$D/r_1.fq.gz" -2 "$D/r_2.fq.gz" -t 4; then
    detail=$(gen exact "$D/out/contigs.fasta" "$D/g.fa"); rc=$?
    check "trivial 20kb reconstruction" $rc "$detail"
else
    fail "trivial 20kb reconstruction" "TesserACT exited $? ($(tail -1 "$LOG"))"
fi

# ---------------------------------------------------------------------------
# 2. Two unrelated chromosomes in one run
# ---------------------------------------------------------------------------
D=$TMP/t2; mkdir -p "$D"
gen genome "$D/g.fa" 2 10000 10000 >/dev/null
gen reads "$D/g.fa" "$D/r" 50 150 350 30 0 1002 paired fastq gz >/dev/null
if asm "$D/out" -1 "$D/r_1.fq.gz" -2 "$D/r_2.fq.gz" -t 4; then
    detail=$(gen exact "$D/out/contigs.fasta" "$D/g.fa"); rc=$?
    check "two chromosomes -> two contigs" $rc "$detail"
else
    fail "two chromosomes -> two contigs" "TesserACT exited $? ($(tail -1 "$LOG"))"
fi

# ---------------------------------------------------------------------------
# 3. Error tolerance: 1% substitutions
# ---------------------------------------------------------------------------
D=$TMP/t3; mkdir -p "$D"
gen genome "$D/g.fa" 3 20000 >/dev/null
gen reads "$D/g.fa" "$D/r" 50 150 350 30 0.01 1003 paired fastq gz >/dev/null
if asm "$D/out" -1 "$D/r_1.fq.gz" -2 "$D/r_2.fq.gz" -t 4; then
    eval "$(gen stats "$D/out/contigs.fasta")"
    detail=$(gen identity "$D/out/contigs.fasta" "$D/g.fa" 0.999); rc=$?
    [ "$n" -eq 1 ] || rc=1
    check "1% substitution errors" $rc "contigs=$n $detail"
else
    fail "1% substitution errors" "TesserACT exited $? ($(tail -1 "$LOG"))"
fi

# ---------------------------------------------------------------------------
# 4. Repeat resolution: 800 bp repeat, 1200 bp fragments span it
# ---------------------------------------------------------------------------
D=$TMP/t4; mkdir -p "$D"
gen repeat_genome "$D/g.fa" 4 5000 800 >/dev/null
gen reads "$D/g.fa" "$D/r" 60 150 1200 60 0 1004 paired fastq gz >/dev/null
EXPECT=16600
if asm "$D/out" -1 "$D/r_1.fq.gz" -2 "$D/r_2.fq.gz" -t 4; then
    eval "$(gen stats "$D/out/contigs.fasta")"
    delta=$(absdiff "$total" "$EXPECT")
    rc=0
    [ "$n" -eq 1 ] || rc=1
    [ $((delta * 100)) -le "$EXPECT" ] || rc=1
    gen substr "$D/out/contigs.fasta" "$D/g.fa" >/dev/null || rc=1
    check "spanned repeat resolved to one contig" $rc \
          "contigs=$n length=$total expected=$EXPECT delta=$delta"
else
    fail "spanned repeat resolved to one contig" "TesserACT exited $? ($(tail -1 "$LOG"))"
fi

# ---------------------------------------------------------------------------
# 5. Unresolvable repeat must not be misassembled
#    5 kb repeat, 350 bp fragments -- no pair can span it.
# ---------------------------------------------------------------------------
D=$TMP/t5; mkdir -p "$D"
gen repeat_genome "$D/g.fa" 5 5000 5000 >/dev/null
gen reads "$D/g.fa" "$D/r" 60 150 350 30 0 1005 paired fastq gz >/dev/null
GLEN=25000
if asm "$D/out" -1 "$D/r_1.fq.gz" -2 "$D/r_2.fq.gz" -t 4; then
    eval "$(gen stats "$D/out/contigs.fasta")"
    rc=0
    sub=$(gen substr "$D/out/contigs.fasta" "$D/g.fa") || rc=1
    km=$(gen kmercheck "$D/out/contigs.fasta" "$D/g.fa" 31) || rc=1
    # No duplication blowup, and nothing long enough to have jumped the repeat.
    [ "$total" -le $((GLEN * 12 / 10)) ] || rc=1
    [ "$largest" -lt $((GLEN * 8 / 10)) ] || rc=1
    check "unspanned repeat left unjoined" $rc "$sub $km total=$total largest=$largest"
else
    fail "unspanned repeat left unjoined" "TesserACT exited $? ($(tail -1 "$LOG"))"
fi

# ---------------------------------------------------------------------------
# 6-9. Input handling. One 15 kb genome, five encodings of the same reads.
# ---------------------------------------------------------------------------
D=$TMP/t6; mkdir -p "$D"
gen genome "$D/g.fa" 6 15000 >/dev/null
gen reads "$D/g.fa" "$D/r" 50 150 350 30 0 1006 paired fastq gz >/dev/null
gen reads "$D/g.fa" "$D/i" 50 150 350 30 0 1006 interleaved fastq gz >/dev/null
gen reads "$D/g.fa" "$D/s" 50 150 350 30 0 1006 single fastq gz >/dev/null
gen reads "$D/g.fa" "$D/a" 50 150 350 30 0 1006 paired fasta plain >/dev/null
gen reads "$D/g.fa" "$D/p" 50 150 350 30 0 1006 paired fastq plain >/dev/null

asm "$D/paired" -1 "$D/r_1.fq.gz" -2 "$D/r_2.fq.gz" -t 4
paired_rc=$?

# 6. single-end
if asm "$D/single" -s "$D/s_s.fq.gz" -t 4; then
    eval "$(gen stats "$D/single/contigs.fasta")"
    rc=0
    [ "$n" -ge 1 ] || rc=1
    [ "$total" -ge 14000 ] || rc=1
    gen substr "$D/single/contigs.fasta" "$D/g.fa" >/dev/null || rc=1
    check "single-end input (-s)" $rc "contigs=$n total=$total"
else
    fail "single-end input (-s)" "TesserACT exited $? ($(tail -1 "$LOG"))"
fi

# 7. interleaved
if [ $paired_rc -eq 0 ] && asm "$D/inter" --12 "$D/i_12.fq.gz" -t 4; then
    detail=$(gen sameset "$D/paired/contigs.fasta" "$D/inter/contigs.fasta"); rc=$?
    check "interleaved (--12) matches two-file run" $rc "$detail"
else
    fail "interleaved (--12) matches two-file run" "TesserACT exited non-zero ($(tail -1 "$LOG"))"
fi

# 8. FASTA input
if [ $paired_rc -eq 0 ] && asm "$D/fasta" -1 "$D/a_1.fa" -2 "$D/a_2.fa" -t 4; then
    detail=$(gen sameset "$D/paired/contigs.fasta" "$D/fasta/contigs.fasta"); rc=$?
    check "FASTA input accepted" $rc "$detail"
else
    fail "FASTA input accepted" "TesserACT exited non-zero ($(tail -1 "$LOG"))"
fi

# 9. gzipped vs plain
if [ $paired_rc -eq 0 ] && asm "$D/plain" -1 "$D/p_1.fq" -2 "$D/p_2.fq" -t 4; then
    diffs=$(same_outputs "$D/paired" "$D/plain")
    if [ -z "$diffs" ]; then
        pass "gzipped and plain input identical" "byte-identical $(present_outputs "$D/plain")"
    else
        fail "gzipped and plain input identical" "outputs differ:$diffs"
    fi
else
    fail "gzipped and plain input identical" "TesserACT exited non-zero ($(tail -1 "$LOG"))"
fi

# ---------------------------------------------------------------------------
# 10. -k validation
# ---------------------------------------------------------------------------
kbad() { # kbad LABEL KSPEC
    local err out
    err=$($TIMEOUT "$TESSERACT" -1 "$TMP/t6/r_1.fq.gz" -2 "$TMP/t6/r_2.fq.gz" \
          -o "$TMP/t10" -k "$2" 2>&1 >/dev/null)
    local status=$?
    if [ $status -ne 0 ] && [ -n "$err" ]; then
        pass "reject $1" "exit=$status \"$err\""
    else
        fail "reject $1" "exit=$status stderr=\"$err\""
    fi
}
kbad "even k (-k 22)" 22
kbad "k above the 127 cap (-k 129)" 129
kbad "k < 5 (-k 3)" 3

# The boundary value the cap allows must still assemble. k-mers are packed into
# four 64-bit words, so 127 is the largest odd k the representation holds.
if asm "$TMP/t10ok" -1 "$TMP/t6/r_1.fq.gz" -2 "$TMP/t6/r_2.fq.gz" -k 127 -t 4; then
    eval "$(gen stats "$TMP/t10ok/contigs.fasta")"
    check "largest legal k (-k 127) accepted" $([ "$n" -ge 1 ] && echo 0 || echo 1) \
          "contigs=$n total=$total"
else
    fail "largest legal k (-k 127) accepted" "TesserACT exited non-zero ($(tail -1 "$LOG"))"
fi

# ---------------------------------------------------------------------------
# 11. Missing input file
# ---------------------------------------------------------------------------
err=$($TIMEOUT "$TESSERACT" -1 "$TMP/does-not-exist.fq.gz" -o "$TMP/t11" 2>&1 >/dev/null)
status=$?
if [ $status -ne 0 ] && printf '%s' "$err" | grep -qi "not found"; then
    pass "missing input file reported" "exit=$status \"$err\""
else
    fail "missing input file reported" "exit=$status stderr=\"$err\""
fi

# ---------------------------------------------------------------------------
# 12-15 reuse the 4-contig assembly from the unresolvable-repeat genome.
# ---------------------------------------------------------------------------
D=$TMP/t5
eval "$(gen stats "$D/out/contigs.fasta")"
base_n=$n; base_total=$total; base_n50=$n50; base_largest=$largest

# 12. --min-contig
rc=0
if asm "$D/min_hi" -1 "$D/r_1.fq.gz" -2 "$D/r_2.fq.gz" -t 4 --min-contig "$base_largest"; then
    eval "$(gen stats "$D/min_hi/contigs.fasta")"
    hi_n=$n
    # Only contigs at least as long as the longest one survive.
    [ "$hi_n" -ge 1 ] || rc=1
    [ "$hi_n" -lt "$base_n" ] || rc=1
else
    rc=1; hi_n="?"
fi
if asm "$D/min_hu" -1 "$D/r_1.fq.gz" -2 "$D/r_2.fq.gz" -t 4 --min-contig 1000000; then
    eval "$(gen stats "$D/min_hu/contigs.fasta")"
    huge_n=$n
    [ "$huge_n" -eq 0 ] || rc=1
else
    rc=1; huge_n="?"
fi
check "--min-contig filters output" $rc \
      "unfiltered=$base_n min=$base_largest -> $hi_n, min=1000000 -> $huge_n"

# 13. Determinism
rc=0
asm "$D/det1" -1 "$D/r_1.fq.gz" -2 "$D/r_2.fq.gz" -t 4 || rc=1
asm "$D/det2" -1 "$D/r_1.fq.gz" -2 "$D/r_2.fq.gz" -t 4 || rc=1
if [ $rc -eq 0 ] && [ -z "$(same_outputs "$D/det1" "$D/det2")" ]; then
    pass "determinism (same threads, two runs)" "byte-identical $(present_outputs "$D/det1")"
else
    fail "determinism (same threads, two runs)" \
         "run failed or outputs differ:$(same_outputs "$D/det1" "$D/det2")"
fi

# 14. Thread invariance
rc=0
asm "$D/th1" -1 "$D/r_1.fq.gz" -2 "$D/r_2.fq.gz" -t 1 || rc=1
asm "$D/th8" -1 "$D/r_1.fq.gz" -2 "$D/r_2.fq.gz" -t 8 || rc=1
if [ $rc -eq 0 ]; then
    eval "$(gen stats "$D/th1/contigs.fasta")"; t1_total=$total; t1_n50=$n50
    eval "$(gen stats "$D/th8/contigs.fasta")"; t8_total=$total; t8_n50=$n50
    # Byte equality, not total and N50. Those two agree between assemblies that differ in
    # contig order or in which bases sit where, so the weaker check passes on output the
    # thread count demonstrably changed: injecting a tie-break that depends on opt_.threads
    # left this test green while the "(byte-identical)" note it prints quietly disappeared.
    # A run that reorders its own output by thread count is not thread-invariant, and this is
    # the test that is supposed to say so.
    # Every written output, not only contigs.fasta, and -t 4 (det1 above) as well as -t 8.
    th_diff="$(same_outputs "$D/th1" "$D/th8")"
    [ "$rc" -eq 0 ] && [ -d "$D/det1" ] && th_diff="$th_diff$(same_outputs "$D/th1" "$D/det1" | sed 's/ / t4:/g')"
    if [ -z "$th_diff" ]; then
        pass "thread invariance (-t 1 vs -t 4 vs -t 8)" \
             "total=$t1_total n50=$t1_n50 (byte-identical $(present_outputs "$D/th1"))"
    elif [ "$t1_total" -eq "$t8_total" ] && [ "$t1_n50" -eq "$t8_n50" ]; then
        fail "thread invariance (-t 1 vs -t 4 vs -t 8)" \
             "same total ($t1_total) and n50 ($t1_n50) but differ byte for byte:$th_diff"
    else
        fail "thread invariance (-t 1 vs -t 4 vs -t 8)" \
             "t1: total=$t1_total n50=$t1_n50 / t8: total=$t8_total n50=$t8_n50"
    fi
else
    fail "thread invariance (-t 1 vs -t 4 vs -t 8)" "a run exited non-zero"
fi

# 15. Output format
detail=$(gen format "$D/out/contigs.fasta"); rc=$?
check "contigs.fasta header format and order" $rc "$detail"

# 15b. The contig statistics in report.json describe the contigs actually written.
#
# n50 and largest are scaffold figures whenever the assembly carries gaps: layout asserts an
# order across a gap and fills it with N, which raises both without assembling another base.
# The summary printed "N50 3,034,032" on a real S. aureus isolate whose contig N50 was 86,836,
# a 35x gap, and next to it "0 joins spanning 97,073 N bases" -- a join count taken from one
# stage while the N came from another. Both numbers are now derived from the sequences, and
# this pins them to the file so they cannot drift apart again.
detail=$(gen contigstats "$D/out"); rc=$?
check "report.json contig stats match contigs.fasta" $rc "$detail"

# ---------------------------------------------------------------------------
# 15c. assembly_graph.gfa is well formed, and keeps its P records
#
# The suite had no GFA coverage at all, which is how a real defect survived: gfaPaths was
# cleared unconditionally before layout, but layout returns early on several paths without
# touching the contigs -- so on those runs the walks were still valid and were discarded
# anyway, and assembly_graph.gfa came out with S and L records and not one P. Nothing failed
# and nothing warned; report.json said layout.run=false, which reads as "changed nothing".
#
# This checks both halves: the file parses and no link or path names a segment that is not
# there, and at least one P record survives a run where layout did not rearrange anything.
# ---------------------------------------------------------------------------
D=$TMP/t15c; mkdir -p "$D"
# An unspanned repeat, so the assembly comes out in several pieces. The layout stage is
# guarded on seqs.size() > 1, so a genome that assembles into one contig never reaches the
# code this test exists to cover -- which is the second way this test was written wrong.
gen repeat_genome "$D/g.fa" 42 5000 5000 >/dev/null
gen reads "$D/g.fa" "$D/r" 60 150 350 30 0 1042 paired fastq gz >/dev/null
# A model carrying layout tracks, built from an UNRELATED genome. That combination is the
# one the defect needed: layout is attempted, finds far too few shared markers, and returns
# without touching the contigs -- so the graph walks are still valid and must survive. A run
# with no model never enters that code at all, so testing without one looks like a pass
# whether the bug is present or not. This test was written that way first, and reintroducing
# the bug did not fail it.
gen repeat_genome "$D/other.fa" 99 5000 5000 >/dev/null
# No silent fallback. An earlier version of this test swallowed a failed model build and
# carried on without --model, which drops it straight back into the no-model case where the
# layout code is never entered -- so it passed whether the defect was present or not. If the
# model cannot be built the test has nothing to say and must report that, not pass.
$TIMEOUT "$ROOT/tesseract-model" --organism testus --out "$D/tracks.tsm" \
    --layout-tracks --min-support 1 "$D/other.fa" > "$D/model.log" 2>&1
if [ ! -s "$D/tracks.tsm" ]; then
    fail "assembly_graph.gfa well formed, P records kept" \
         "could not build the track model: $(tail -1 "$D/model.log")"
elif asm "$D/out" -1 "$D/r_1.fq.gz" -2 "$D/r_2.fq.gz" -t 2 \
         --organism testus --model "$D/tracks.tsm"; then
    detail=$(python3 - "$D/out/assembly_graph.gfa" <<'PYGFA'
import sys
segs, dangling, pathbad, paths, links, malformed = set(), 0, 0, 0, 0, 0
for line in open(sys.argv[1]):
    f = line.rstrip("\n").split("\t")
    if f[0] == "S":
        if len(f) < 3: malformed += 1
        else: segs.add(f[1])
    elif f[0] == "L":
        links += 1
        if len(f) < 6: malformed += 1
    elif f[0] == "P":
        paths += 1
        if len(f) < 3: malformed += 1
# Second pass: references can only be checked once every segment is known.
for line in open(sys.argv[1]):
    f = line.rstrip("\n").split("\t")
    if f[0] == "L" and len(f) >= 6:
        if f[1] not in segs or f[3] not in segs: dangling += 1
    elif f[0] == "P" and len(f) >= 3:
        if any(s[:-1] not in segs for s in f[2].split(",") if s): pathbad += 1
print("segments=%d links=%d paths=%d malformed=%d dangling=%d pathbad=%d"
      % (len(segs), links, paths, malformed, dangling, pathbad))
sys.exit(1 if (malformed or dangling or pathbad or not segs or not paths) else 0)
PYGFA
); rc=$?
    check "assembly_graph.gfa well formed, P records kept" $rc "$detail"
    # Same invariant as 15b, but on the run that goes through layout -- the only stage that
    # scaffolds, so the only one that can produce the gaps the contig figures exist to expose.
    # The detail line prints the gap count, so a run where layout happened to join nothing is
    # visible as gaps=0 rather than passing as if the gap arithmetic had been exercised. (It
    # does print gaps=0 here; test 15d below is the fixture that actually writes a gap.)
    detail=$(gen contigstats "$D/out"); rc=$?
    check "contig stats match after layout" $rc "$detail"
else
    fail "assembly_graph.gfa well formed, P records kept" "tesseract-asm exited $?"
fi

# ---------------------------------------------------------------------------
# 15d. A scaffold gap: report.json, scaffolds.fasta/.agp and thread invariance
#
# No other fixture writes a gap (every contigstats line above says gaps=0), so the gap
# arithmetic of report.json, the N rows of the AGP and the thread invariance of
# scaffolds.fasta/.agp were never exercised. Here no read covers a 200 bp window of a 40 kb
# genome while 600 bp fragments span it: the pairs bridge the dropout and the assembler
# writes one scaffold gap (293 N from release 1.3.0).
# ---------------------------------------------------------------------------
D=$TMP/t15d; mkdir -p "$D"
gen genome "$D/g.fa" 7 40000 >/dev/null
gen holereads "$D/g.fa" "$D/r" 60 150 600 40 1007 20000 20200 >/dev/null
if asm "$D/out" -1 "$D/r_1.fq.gz" -2 "$D/r_2.fq.gz" -t 4; then
    detail=$(gen contigstats "$D/out"); rc=$?
    # The fixture is only worth something while it produces a gap.
    case "$detail" in gaps=0\ *) rc=1; detail="$detail (fixture produced no gap)";; esac
    sub=$(gen substr "$D/out/contigs.fasta" "$D/g.fa") || rc=1
    check "gapped scaffold: stats, AGP and pieces" $rc "$detail $sub"
    if asm "$D/th1" -1 "$D/r_1.fq.gz" -2 "$D/r_2.fq.gz" -t 1; then
        diffs=$(same_outputs "$D/out" "$D/th1")
        if [ -z "$diffs" ]; then
            pass "gapped scaffold thread invariance (-t 1/-t 4)" \
                 "byte-identical $(present_outputs "$D/th1")"
        else
            fail "gapped scaffold thread invariance (-t 1/-t 4)" "differ:$diffs"
        fi
    else
        fail "gapped scaffold thread invariance (-t 1/-t 4)" "TesserACT exited non-zero ($(tail -1 "$LOG"))"
    fi
else
    fail "gapped scaffold: stats, AGP and pieces" "TesserACT exited non-zero ($(tail -1 "$LOG"))"
fi

# 15e. The contigstats check itself: hand-written output directories covering what reads
# cannot easily force -- a 1-N and a 12-N gap, a gap the polisher overwrote, a terminal
# trim, an N left in contigs.fasta, an AGP without its N row -- each with a report that
# describes it (must pass) or a stale one (must fail).
mkdir -p "$TMP/t15e"
detail=$(gen contigstats_selftest "$TMP/t15e"); rc=$?
check "contigstats check judges stale reports" $rc "$detail"

# ---------------------------------------------------------------------------
# 15b. Model marker density round-trips through the file
#
# Build and query must agree on the sampling denominator exactly. Sampling is by hash
# threshold, so a query at 512 against a model built at 64 shares an eighth of the markers
# and reports almost nothing -- without erroring, which is what makes it worth a test. The
# model records its own density from version 5 on; this checks that the flag reaches the
# build, that the denser model really is denser, and that both files load.
# ---------------------------------------------------------------------------
MODELBIN=$ROOT/tesseract-model
if [ -x "$MODELBIN" ]; then
    gen genome "$TMP/m.fa" 15 40000 >/dev/null
    m512=$($TIMEOUT "$MODELBIN" --organism testus --out "$TMP/m512.tsm" "$TMP/m.fa" 2>&1 |
           grep -oE '[0-9]+ markers' | tail -1 | cut -d' ' -f1)
    m64=$($TIMEOUT "$MODELBIN" --organism testus --out "$TMP/m64.tsm" --marker-density 64 \
          "$TMP/m.fa" 2>&1 | grep -oE '[0-9]+ markers' | tail -1 | cut -d' ' -f1)
    if [ -s "$TMP/m512.tsm" ] && [ -s "$TMP/m64.tsm" ] &&
       [ -n "$m512" ] && [ -n "$m64" ] && [ "$m64" -gt "$m512" ]; then
        # Both must be usable by this one binary, which is the whole point of recording the
        # density in the file rather than compiling it in.
        gen reads "$TMP/m.fa" "$TMP/mr" 50 150 350 30 0 1015 paired fastq gz >/dev/null
        ok512=0; ok64=0
        # -t as everywhere else: without it these two ran on every core of the machine.
        $TIMEOUT "$TESSERACT" -1 "$TMP/mr_1.fq.gz" -2 "$TMP/mr_2.fq.gz" -o "$TMP/t15b_512" -t 4 \
            --organism testus --model "$TMP/m512.tsm" >/dev/null 2>&1 && ok512=1
        $TIMEOUT "$TESSERACT" -1 "$TMP/mr_1.fq.gz" -2 "$TMP/mr_2.fq.gz" -o "$TMP/t15b_64" -t 4 \
            --organism testus --model "$TMP/m64.tsm" >/dev/null 2>&1 && ok64=1
        if [ "$ok512" -eq 1 ] && [ "$ok64" -eq 1 ]; then
            pass "model marker density round-trips" \
                 "markers 1/512=$m512 1/64=$m64 (both models loaded and assembled)"
        else
            fail "model marker density round-trips" \
                 "assembly failed: 512=$ok512 64=$ok64"
        fi
    else
        fail "model marker density round-trips" \
             "markers 1/512=${m512:-none} 1/64=${m64:-none} (denser must yield more)"
    fi
else
    skip "model marker density round-trips" "no tesseract-model built (make test builds it)"
fi

# ---------------------------------------------------------------------------
# 16. Empty and tiny input
# ---------------------------------------------------------------------------
: > "$TMP/empty.fq"
err=$($TIMEOUT "$TESSERACT" -1 "$TMP/empty.fq" -o "$TMP/t16a" 2>&1 >/dev/null)
status=$?
if [ $status -ne 0 ] && [ $status -ne 124 ] && [ -n "$err" ]; then
    pass "empty input rejected cleanly" "exit=$status \"$(printf '%s' "$err" | tail -1)\""
else
    fail "empty input rejected cleanly" "exit=$status stderr=\"$err\""
fi

gen tiny "$TMP/tiny.fq" >/dev/null
err=$($TIMEOUT "$TESSERACT" -s "$TMP/tiny.fq" -o "$TMP/t16b" -t 4 2>&1 >/dev/null)
status=$?
if [ $status -ne 0 ] && [ $status -ne 124 ] && [ -n "$err" ]; then
    pass "tiny input rejected cleanly" "exit=$status \"$(printf '%s' "$err" | tail -1)\""
else
    fail "tiny input rejected cleanly" "exit=$status stderr=\"$err\""
fi

# ---------------------------------------------------------------------------
# 17. Input handling of tracked open defects (XFAIL until fixed)
#
# Each runs with TESSERACT_FIXES=1, which enables every default-off fix, and checks the
# corrected behaviour; the owning group turns its XFAIL into a PASS when the fix lands.
# T18 truncated gzip accepted; T32 conflicting read options / bad -k lists accepted;
# T33 benign FASTQ variants rejected; T34 a reused -o keeps stale outputs.
# ---------------------------------------------------------------------------
D=$TMP/t17; mkdir -p "$D"
S6=$TMP/t6
# fixed ASM_ARGS... -- run with the fixes enabled; status in $?, stderr in $LOG
fixed() { TESSERACT_FIXES=1 $TIMEOUT "$TESSERACT" "$@" -q >"$LOG" 2>&1; }

gen truncgz "$S6/s_s.fq.gz" "$D/trunc.fq.gz" 4000 >/dev/null
rm -rf "$D/o1"; fixed -s "$D/trunc.fq.gz" -o "$D/o1" -t 4; st=$?
known "truncated gzip rejected" T18 $([ $st -ne 0 ] && [ $st -ne 124 ] && echo 0 || echo 1) \
      "exit=$st $(tail -1 "$LOG" | cut -c1-80)"

rejects() { # rejects LABEL ASM_ARGS... -- the option combination must be refused
    local label=$1; shift
    rm -rf "$D/o2"; fixed "$@" -o "$D/o2" -t 4; local st=$?
    known "$label" T32 $([ $st -ne 0 ] && [ $st -ne 124 ] && [ -s "$LOG" ] && echo 0 || echo 1) \
          "exit=$st $(tail -1 "$LOG" | cut -c1-80)"
}
rejects "reject repeated -1" -1 "$S6/r_1.fq.gz" -1 "$S6/r_2.fq.gz"
rejects "reject the same file as -1 and -2" -1 "$S6/r_1.fq.gz" -2 "$S6/r_1.fq.gz"
# (mates of one file each, so the release run is not stopped by an unrelated name mismatch)
rejects "reject --12 together with -2" --12 "$S6/r_1.fq.gz" -2 "$S6/r_2.fq.gz"
rejects "reject unsorted -k (77,21,55)" -1 "$S6/r_1.fq.gz" -2 "$S6/r_2.fq.gz" -k 77,21,55
rejects "reject duplicate -k (21,21)" -1 "$S6/r_1.fq.gz" -2 "$S6/r_2.fq.gz" -k 21,21
rejects "reject fractional -k (21.9,55)" -1 "$S6/r_1.fq.gz" -2 "$S6/r_2.fq.gz" -k 21.9,55

cp "$S6/p_1.fq" "$D/b_1.fq"; cp "$S6/p_2.fq" "$D/b_2.fq"; echo >> "$D/b_1.fq"; echo >> "$D/b_2.fq"
rm -rf "$D/o3"; fixed -1 "$D/b_1.fq" -2 "$D/b_2.fq" -o "$D/o3" -t 4; st=$?
rc=1; [ $st -eq 0 ] && gen sameset "$S6/paired/contigs.fasta" "$D/o3/contigs.fasta" >/dev/null && rc=0
known "FASTQ with a trailing blank line accepted" T33 $rc "exit=$st $(tail -1 "$LOG" | cut -c1-80)"

gen dotnames "$S6/p_1.fq" "$D/d_1.fq" 1 >/dev/null; gen dotnames "$S6/p_2.fq" "$D/d_2.fq" 2 >/dev/null
rm -rf "$D/o4"; fixed -1 "$D/d_1.fq" -2 "$D/d_2.fq" -o "$D/o4" -t 4; st=$?
rc=1; [ $st -eq 0 ] && gen sameset "$S6/paired/contigs.fasta" "$D/o4/contigs.fasta" >/dev/null && rc=0
known "'.1'/'.2' mate suffixes accepted" T33 $rc "exit=$st $(tail -1 "$LOG" | cut -c1-80)"

# A gap-free run into a directory holding an earlier run's optional outputs: none of them
# may survive into the new result (deleted or rewritten).
rm -rf "$D/o5"; mkdir -p "$D/o5"
for f in scaffolds.fasta scaffolds.agp unitigs.fasta; do printf '>stale\nACGT\n' > "$D/o5/$f"; done
fixed -1 "$S6/r_1.fq.gz" -2 "$S6/r_2.fq.gz" -o "$D/o5" -t 4; st=$?
stale=$(grep -l '^>stale' "$D/o5/scaffolds.fasta" "$D/o5/scaffolds.agp" "$D/o5/unitigs.fasta" 2>/dev/null | xargs -r -n1 basename | tr '\n' ' ')
known "reused -o leaves no stale outputs" T34 $([ $st -eq 0 ] && [ -z "$stale" ] && echo 0 || echo 1) \
      "exit=$st stale=[${stale% }]"

# ---------------------------------------------------------------------------
# 18. 1.4.0 defaults (RELEASE_PLAN 1.4.0 step 5.2)
#
# A run with nothing set uses the combo3 F2 configuration and says so on its counter lines.
# The documented opt-out list restores release 1.3.0 byte for byte: tests/golden_130.md5
# holds what release 1.3.0 (work/tesseract-1.3.0, md5 4af92d1f) wrote on the t5 and t15d
# fixtures above, generated once. On t15d (the scaffold gap) the 1.4.0 defaults write
# different files, so the comparison is not vacuous.
# ---------------------------------------------------------------------------
D=$TMP/t18; mkdir -p "$D"
if asm "$D/def" -1 "$TMP/t15d/r_1.fq.gz" -2 "$TMP/t15d/r_2.fq.gz" -t 4; then
    miss=""
    rf=$(grep -m1 '^\[resolveflags\] ' "$LOG")
    for tok in run=1 tieRatio=3.000 minFallbackDest=1000000000 requireSupportSingle=1 excludeSharedRepeat=1; do
        case " $rf " in *" $tok "*) ;; *) miss="$miss resolveflags:$tok";; esac
    done
    grep -q '^\[fixes-emit\] umbrella=1 ' "$LOG" || miss="$miss fixes-emit:umbrella=1"
    grep -q '^\[dropoutbridge\] .* enabled=1$' "$LOG" || miss="$miss dropoutbridge:enabled=1"
    grep -q '^\[endext\] budget=250 ' "$LOG" || miss="$miss endext:budget=250"
    grep -qx '\[defaults\] tie_ratio=3.000 common_prefix=250 prefix_min_body=reach min_fallback_dest=1000000000 require_support_single=1 exclude_shared_repeat_support=1 fixes=1 dropout_bridge=1' "$LOG" \
        || miss="$miss [defaults]"
    if [ -z "$miss" ]; then
        pass "default run reports the 1.4.0 defaults" "resolveflags, fixes-emit, dropoutbridge, endext, [defaults]"
    else
        fail "default run reports the 1.4.0 defaults" "missing:$miss"
    fi
else
    fail "default run reports the 1.4.0 defaults" "TesserACT exited non-zero ($(tail -1 "$LOG"))"
fi

# One version everywhere (N19): --version, report.json, report.html and src/version.h.
if [ -s "$D/def/report.json" ]; then
    vbin=$("$TESSERACT" --version 2>/dev/null | sed -n '1s/^TesserACT //p')
    vsrc=$(sed -n 's/.*kVersion = "\([^"]*\)".*/\1/p' "$ROOT/src/version.h")
    vjson=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["version"])' "$D/def/report.json" 2>/dev/null)
    vhtml=$(grep -c "Generated by tesseract-asm $vsrc" "$D/def/report.html" 2>/dev/null)
    if [ -n "$vsrc" ] && [ "$vbin" = "$vsrc" ] && [ "$vjson" = "$vsrc" ] && [ "${vhtml:-0}" -ge 1 ]; then
        pass "one version in --version, report.json, report.html" "$vsrc"
    else
        fail "one version in --version, report.json, report.html" \
             "src=$vsrc --version=$vbin report.json=$vjson report.html_matches=${vhtml:-0}"
    fi
else
    fail "one version in --version, report.json, report.html" "no report.json from the default run"
fi

# optout130 OUTDIR [TesserACT args...] -- the 1.3.0 opt-out list from --help, on this command only
optout130() {
    local out=$1; shift
    rm -rf "$out"
    TESSERACT_FIXES=0 TESSERACT_COMMON_PREFIX=3000 TESSERACT_PREFIX_MIN_BODY=0 \
    TESSERACT_MIN_FALLBACK_DEST=0 TESSERACT_REQUIRE_SUPPORT_SINGLE=0 \
    TESSERACT_EXCLUDE_SHARED_REPEAT_SUPPORT=0 TESSERACT_DROPOUT_BRIDGE=0 \
        $TIMEOUT "$TESSERACT" --tie-ratio 1.02 "$@" -o "$out" -q >"$LOG" 2>&1
}
GOLD=$ROOT/tests/golden_130.md5
for fx in t5 t15d; do
    if [ ! -r "$GOLD" ]; then
        fail "1.3.0 opt-out reproduces release 1.3.0 ($fx)" "missing $GOLD"
    elif optout130 "$D/o_$fx" -1 "$TMP/$fx/r_1.fq.gz" -2 "$TMP/$fx/r_2.fq.gz" -t 4; then
        bad=""; n=0
        while read -r name f want; do
            [ "$name" = "$fx" ] || continue
            n=$((n + 1))
            if [ "$want" = absent ]; then
                [ -e "$D/o_$fx/$f" ] && bad="$bad $f(written)"
            else
                got=$(md5sum < "$D/o_$fx/$f" 2>/dev/null | cut -c1-32)
                [ "$got" = "$want" ] || bad="$bad $f"
            fi
        done < <(grep -v '^#' "$GOLD")
        [ "$n" -eq 4 ] || bad="$bad (golden has $n rows for $fx)"
        if [ -z "$bad" ]; then
            pass "1.3.0 opt-out reproduces release 1.3.0 ($fx)" "4/4 files as release 1.3.0 wrote them"
        else
            fail "1.3.0 opt-out reproduces release 1.3.0 ($fx)" "differ:$bad"
        fi
    else
        fail "1.3.0 opt-out reproduces release 1.3.0 ($fx)" "TesserACT exited non-zero ($(tail -1 "$LOG"))"
    fi
done

# ---------------------------------------------------------------------------
# 19. A model selected by organism: directly, and through each helper (N22)
#
# Since 1.3.0 --model is refused unless TESSERACT_MODEL_AUTHOR is set, and tesseract-eskape and
# tesseract-klebsiella still passed it, so every model run through them failed at startup.
# They also passed `--organism klebsiella`, while the installer, tesseract-get-models and
# models.sha256 name the model kpneumoniae.tsm. Nothing here caught it: every model test above
# runs under the TESSERACT_MODEL_AUTHOR this suite exports. Here the model is built with it and
# every run is made without it -- the environment a user has -- and passes only when the run
# succeeds and the model stage ran: report.json says organism_model.run for the model named,
# and the log carries the stage's "model joining" line.
#
# Models are opt-in (1.4.0): tesseract-klebsiella runs without one unless given --with-model,
# --model or --plasmid, even when a model is installed, and a tesseract-eskape preset whose
# model is not installed carries on with the defaults instead of stopping.
# ---------------------------------------------------------------------------
D=$TMP/t19; mkdir -p "$D/models" "$D/cache" "$D/elsewhere"
gen repeat_genome "$D/g.fa" 43 5000 5000 >/dev/null
gen reads "$D/g.fa" "$D/r" 60 150 350 30 0 1043 paired fastq gz >/dev/null
# n22model ORGANISM FILE -- a small track model of the fixture genome, built under
# TESSERACT_MODEL_AUTHOR (the build is the only step that has it)
n22model() {
    TESSERACT_MODEL_AUTHOR=1 $TIMEOUT "$ROOT/tesseract-model" --organism "$1" --out "$2" \
        --layout-tracks --min-support 1 "$D/g.fa" >"$2.log" 2>&1 && [ -s "$2" ]
}
# n22user CMD... -- run as a user would: TESSERACT_MODEL_AUTHOR removed from the environment
n22user() { env -u TESSERACT_MODEL_AUTHOR $TIMEOUT "$@"; }
# modelran NAME STATUS OUTDIR MODEL_ORGANISM LOGFILE
modelran() {
    local detail rc
    if [ "$2" -ne 0 ]; then
        # tesseract-klebsiella leaves a failed isolate in OUTDIR.part
        fail "$1" "exit=$2 ($(cat "$5" "$3.part/run.log" 2>/dev/null | grep -m1 -i 'error'))"
        return
    fi
    detail=$(python3 - "$3/report.json" "$4" <<'PYMODEL'
import json, sys
try:
    m = json.load(open(sys.argv[1]))["organism_model"]
except Exception as e:
    print("no report.json organism_model (%s)" % e); sys.exit(1)
print("organism_model.run=%s organism=%s genomes=%s" % (m.get("run"), m.get("organism"), m.get("genomes")))
sys.exit(0 if m.get("run") is True and m.get("organism") == sys.argv[2] else 1)
PYMODEL
); rc=$?
    grep -q "model joining" "$5" 2>/dev/null || { rc=1; detail="$detail, no 'model joining' line in the log"; }
    check "$1" $rc "$detail"
}
# nomodelran NAME STATUS OUTDIR LOGFILE -- the run succeeded and no model was loaded: report.json
# says organism_model.run is not true, its command has no --organism, and there is no model stage
nomodelran() {
    local detail rc
    if [ "$2" -ne 0 ]; then
        fail "$1" "exit=$2 ($(cat "$4" "$3.part/run.log" 2>/dev/null | grep -m1 -i 'error'))"
        return
    fi
    detail=$(python3 - "$3/report.json" <<'PYNOMODEL'
import json, sys
try:
    r = json.load(open(sys.argv[1]))
except Exception as e:
    print("no report.json (%s)" % e); sys.exit(1)
m = r.get("organism_model") or {}
cmd = r.get("command", "")
print("organism_model.run=%s, --organism in command: %s" % (m.get("run"), "--organism" in cmd))
sys.exit(1 if m.get("run") is True or "--organism" in cmd else 0)
PYNOMODEL
); rc=$?
    grep -q "model joining" "$4" 2>/dev/null && { rc=1; detail="$detail, a 'model joining' line in the log"; }
    check "$1" $rc "$detail"
}
if ! n22model kpneumoniae "$D/models/kpneumoniae.tsm" || ! n22model ecoli "$D/models/ecoli.tsm" ||
   ! n22model klebsiella "$D/cache/tesseract-klebsiella-default-v1.2.0.tsm"; then
    fail "--organism model selection (N22)" "could not build the fixture models (logs: $D/models, $D/cache)"
else
    cp "$D/models/kpneumoniae.tsm" "$D/elsewhere/my-kleb-model.tsm"
    R1=$D/r_1.fq.gz; R2=$D/r_2.fq.gz

    # tesseract-asm itself, by the canonical name and by the alias
    rm -rf "$D/asm_kp"
    TESSERACT_MODEL_DIR="$D/models" n22user "$TESSERACT" -1 "$R1" -2 "$R2" -o "$D/asm_kp" -t 2 \
        --organism kpneumoniae >"$D/asm_kp.log" 2>&1
    modelran "--organism kpneumoniae via TESSERACT_MODEL_DIR" $? "$D/asm_kp" kpneumoniae "$D/asm_kp.log"
    rm -rf "$D/asm_kleb"
    TESSERACT_MODEL_DIR="$D/models" n22user "$TESSERACT" -1 "$R1" -2 "$R2" -o "$D/asm_kleb" -t 2 \
        --organism klebsiella >"$D/asm_kleb.log" 2>&1
    modelran "--organism klebsiella reads kpneumoniae.tsm" $? "$D/asm_kleb" kpneumoniae "$D/asm_kleb.log"

    # tesseract-eskape: the Klebsiella preset from the model folder and from --model FILE under
    # another name, and one of the other six presets
    TESSERACT_MODEL_DIR="$D/models" n22user "$ROOT/tesseract-eskape" --preset kpneumoniae \
        -1 "$R1" -2 "$R2" -o "$D/esk_kp" -t 2 >"$D/esk_kp.out" 2>&1
    modelran "tesseract-eskape --preset kpneumoniae" $? "$D/esk_kp" kpneumoniae "$D/esk_kp/tesseract-eskape.log"
    n22user "$ROOT/tesseract-eskape" --preset klebsiella --model "$D/elsewhere/my-kleb-model.tsm" \
        -1 "$R1" -2 "$R2" -o "$D/esk_file" -t 2 >"$D/esk_file.out" 2>&1
    modelran "tesseract-eskape --model FILE" $? "$D/esk_file" kpneumoniae "$D/esk_file/tesseract-eskape.log"
    TESSERACT_MODEL_DIR="$D/models" n22user "$ROOT/tesseract-eskape" --preset ecoli \
        -1 "$R1" -2 "$R2" -o "$D/esk_ec" -t 2 >"$D/esk_ec.out" 2>&1
    modelran "tesseract-eskape --preset ecoli" $? "$D/esk_ec" ecoli "$D/esk_ec/tesseract-eskape.log"

    # tesseract-klebsiella. By default no model, although kpneumoniae.tsm is installed where
    # --with-model looks; --with-model uses that file and prints its measured cost; --model FILE
    # implies --with-model (here a 1.2-named file, built as `klebsiella` as the 1.2 models were)
    TESSERACT_MODEL_DIR="$D/models" n22user "$ROOT/tesseract-klebsiella" "$R1" "$R2" -o "$D/kl_def" \
        -t 2 >"$D/kl_def.out" 2>&1
    nomodelran "tesseract-klebsiella: no model by default" $? "$D/kl_def/r" "$D/kl_def/r/run.log"
    TESSERACT_MODEL_DIR="$D/models" n22user "$ROOT/tesseract-klebsiella" "$R1" "$R2" -o "$D/kl_with" \
        -t 2 --with-model >"$D/kl_with.out" 2>&1; st=$?
    if [ "$st" -eq 0 ] && ! grep -q 'misassemblies 70 -> 93' "$D/kl_with.out"; then
        fail "tesseract-klebsiella --with-model" "the run did not print the model's measured cost"
    else
        modelran "tesseract-klebsiella --with-model" $st "$D/kl_with/r" kpneumoniae "$D/kl_with/r/run.log"
    fi
    n22user "$ROOT/tesseract-klebsiella" "$R1" "$R2" -o "$D/kl_file" -t 2 \
        --model "$D/cache/tesseract-klebsiella-default-v1.2.0.tsm" >"$D/kl_file.out" 2>&1
    modelran "tesseract-klebsiella --model FILE" $? "$D/kl_file/r" klebsiella "$D/kl_file/r/run.log"

    # tesseract-eskape: the Klebsiella preset without its model installed carries on with the
    # defaults (it used to stop), and a model named explicitly but missing still stops it
    mkdir -p "$D/nomodels"
    TESSERACT_MODEL_DIR="$D/nomodels" n22user "$ROOT/tesseract-eskape" \
        --preset kpneumoniae -1 "$R1" -2 "$R2" -o "$D/esk_none" --dry-run >"$D/esk_none.out" 2>&1; st=$?
    run=$(grep -m1 '^Running: ' "$D/esk_none.out")
    rc=0
    [ "$st" -eq 0 ] || rc=1
    case "$run" in *--organism*|"") rc=1 ;; esac
    grep -q '^Preset: .*no model installed' "$D/esk_none.out" || rc=1
    n22user "$ROOT/tesseract-eskape" --preset kpneumoniae --model "$D/nomodels/absent.tsm" \
        -1 "$R1" -2 "$R2" -o "$D/esk_absent" --dry-run >"$D/esk_absent.out" 2>&1 && rc=1
    check "tesseract-eskape preset without its model" $rc \
          "dry-run exit=$st, no --organism; an explicit missing --model still stops it"

    # Neither helper hands the assembler --model any more.
    TESSERACT_MODEL_DIR="$D/models" n22user "$ROOT/tesseract-eskape" --preset kpneumoniae \
        -1 "$R1" -2 "$R2" -o "$D/esk_dry" --dry-run >"$D/esk_dry.out" 2>&1; st=$?
    run=$(grep -m1 '^Running: ' "$D/esk_dry.out")
    case "$run" in
        *" --model "*) rc=1 ;;
        *"--organism kpneumoniae"*) rc=0 ;;
        *) rc=1 ;;
    esac
    [ "$st" -eq 0 ] || rc=1
    grep -q -- '--model' "$D/kl_file/r/report.json" 2>/dev/null && rc=1
    check "helpers pass --organism, never --model" $rc \
          "dry-run exit=$st, $(printf '%s' "$run" | grep -o -- '--organism [a-z]*')"
fi

# ---------------------------------------------------------------------------
# 19b. 1.5.0: the Salmonella model (senterica, alias salmonella), and the gated CLI surface
#
# senterica.tsm is selected like the seven ESKAPEE models: `--organism salmonella` reads
# senterica.tsm, `tesseract-eskape --preset salmonella` uses it when installed and prints its
# cost as measured on development isolates. Without an om2 variable a model run writes no genome/ directory, as in
# 1.4.0; the experimental TESSERACT_OM2_LAYOUT_ONLY=1 writes the layout-only view and leaves
# contigs.fasta as it was.
# ---------------------------------------------------------------------------
if [ -s "$D/g.fa" ] && n22model senterica "$D/models/senterica.tsm"; then
    R1=$D/r_1.fq.gz; R2=$D/r_2.fq.gz
    rm -rf "$D/asm_salm"
    TESSERACT_MODEL_DIR="$D/models" n22user "$TESSERACT" -1 "$R1" -2 "$R2" -o "$D/asm_salm" -t 2 \
        --organism salmonella >"$D/asm_salm.log" 2>&1
    modelran "--organism salmonella reads senterica.tsm" $? "$D/asm_salm" senterica "$D/asm_salm.log"
    [ ! -e "$D/asm_salm/genome" ]; check "no genome/ view without --layout-view" $? "1.4.0 output set"
    rm -rf "$D/esk_salm"
    TESSERACT_MODEL_DIR="$D/models" n22user "$ROOT/tesseract-eskape" --preset Salmonella \
        -1 "$R1" -2 "$R2" -o "$D/esk_salm" -t 2 >"$D/esk_salm.out" 2>&1
    st=$?
    modelran "tesseract-eskape --preset Salmonella" $st "$D/esk_salm" senterica "$D/esk_salm/tesseract-eskape.log"
    grep -q 'misassemblies 3 -> 10 on 15 development isolates' "$D/esk_salm.out"; check "the senterica preset prints its dev-measured cost" $? ""
    # 1.5.0 ships no --layout-view (its dev gate G1 failed by one misassembly, RELEASE_PLAN_1.5.md);
    # the organism-detection tokens stay behind TESSERACT_OM2_DETECT
    rc=0
    n22user "$TESSERACT" -1 "$R1" -2 "$R2" -o "$D/lv_opt" --organism senterica --layout-view >"$D/lv_opt.log" 2>&1 && rc=1
    grep -q "unknown option '--layout-view'" "$D/lv_opt.log" || rc=1
    n22user "$TESSERACT" -1 "$R1" -2 "$R2" -o "$D/lv_force" --organism-force >"$D/lv_force.log" 2>&1 && rc=1
    grep -q "unknown option '--organism-force'" "$D/lv_force.log" || rc=1
    check "no --layout-view option; --organism-force needs TESSERACT_OM2_DETECT" $rc ""
    # the experimental environment path of the layout-only writer still runs end to end
    rm -rf "$D/lv"
    TESSERACT_OM2_LAYOUT_ONLY=1 TESSERACT_MODEL_DIR="$D/models" n22user "$TESSERACT" -1 "$R1" -2 "$R2" -o "$D/lv" -t 2 \
        --organism salmonella >"$D/lv.log" 2>&1; st=$?
    rc=0; detail="exit=$st"
    [ "$st" -eq 0 ] || rc=1
    for f in genome/genome.fasta genome/layout.agp genome/junctions.tsv contigs.fasta; do
        [ -s "$D/lv/$f" ] || { rc=1; detail="$detail missing:$f"; }
    done
    grep -q '^\[om2-layout\] layout_only=1' "$D/lv.log" || { rc=1; detail="$detail no [om2-layout] line"; }
    cmp -s "$D/lv/contigs.fasta" "$D/asm_salm/contigs.fasta" || { rc=1; detail="$detail contigs.fasta changed"; }
    python3 - "$D/lv/report.json" <<'PYLV' || { rc=1; detail="$detail report.json"; }
import json, sys
r = json.load(open(sys.argv[1]))
assert (r.get("om2") or {}).get("layout_only", {}).get("confirm") == "isolate"
PYLV
    check "TESSERACT_OM2_LAYOUT_ONLY=1 (experimental): layout.agp written, contigs.fasta unchanged" $rc "$detail"
else
    fail "senterica fixture model" "could not build it (log: $D/models/senterica.tsm.log)"
fi

# ---------------------------------------------------------------------------
# 20. Where models come from, and that nothing fetches one unasked
#
# tesseract-get-models fetches from the models-v2 release (models-v1 is what 1.3.0 installs look
# for, with a different checksum list), keeps only what matches models.sha256 and skips what it
# already has; tesseract-klebsiella --with-model installs kpneumoniae.tsm through it; and the
# non-interactive install downloads no model and says models are optional. No network: a curl
# stand-in on PATH records each URL and serves files from a local folder, or answers 404.
# ---------------------------------------------------------------------------
D=$TMP/t20; mkdir -p "$D/shim" "$D/serve" "$D/dest" "$D/mk" "$D/home"
cat > "$D/shim/curl" <<'SHIM'
#!/usr/bin/env bash
out=""; url=""
while [ $# -gt 0 ]; do
    case "$1" in
        -o|-w|--retry|--retry-delay) [ "$1" = -o ] && out=$2; shift 2 ;;
        -*) shift ;;
        *) url=$1; shift ;;
    esac
done
printf '%s\n' "$url" >> "${CURL_LOG:?}"
f="${CURL_SERVE:-/nonexistent}/${url##*/}"
if [ -f "$f" ]; then cp "$f" "$out"; printf 200; exit 0; fi
printf 404; exit 22
SHIM
chmod +x "$D/shim/curl"
printf '#!/bin/sh\nexit 0\n' > "$D/mk/make"; chmod +x "$D/mk/make"

# (a) the published URL, for every model in the manifest
rm -f "$D/urls_a"
env -u TESSERACT_MODEL_URL PATH="$D/shim:$PATH" CURL_LOG="$D/urls_a" \
    "$ROOT/tesseract-get-models" --dir "$D/dest" >"$D/get_a.out" 2>&1; st=$?
want=$(awk '!/^#/ && NF==2 {print "https://github.com/iowa69/TesserACT/releases/download/models-v2/" $2}' "$ROOT/models.sha256")
rc=0
[ "$st" -ne 0 ] || rc=1                      # every asset 404s here, which must fail
[ -n "$want" ] && [ "$(cat "$D/urls_a" 2>/dev/null)" = "$want" ] || rc=1
grep -q 'not on the release page' "$D/get_a.out" || rc=1
[ -z "$(ls -A "$D/dest")" ] || rc=1          # nothing left behind, .part files included
check "tesseract-get-models fetches from models-v2" $rc \
      "$(grep -c . "$D/urls_a" 2>/dev/null) URLs, $(head -1 "$D/urls_a" 2>/dev/null | sed 's|/[^/]*$||')"

# (b) keeps what verifies, discards what does not, then skips what it has
mkdir -p "$D/gm"
cp "$ROOT/tesseract-get-models" "$D/gm/"
printf 'model bytes\n' > "$D/serve/good.tsm"; printf 'truncated\n' > "$D/serve/bad.tsm"
{ echo "$(sha256sum < "$D/serve/good.tsm" | cut -c1-64)  good.tsm"
  echo "$(printf 'the real thing\n' | sha256sum | cut -c1-64)  bad.tsm"; } > "$D/gm/models.sha256"
rm -f "$D/urls_b"
PATH="$D/shim:$PATH" CURL_LOG="$D/urls_b" CURL_SERVE="$D/serve" \
    "$D/gm/tesseract-get-models" --dir "$D/dest_b" >"$D/get_b.out" 2>&1; st1=$?
PATH="$D/shim:$PATH" CURL_LOG="$D/urls_b" CURL_SERVE="$D/serve" \
    "$D/gm/tesseract-get-models" --dir "$D/dest_b" good >"$D/get_b2.out" 2>&1; st2=$?
rc=0
[ "$st1" -ne 0 ] && [ "$st2" -eq 0 ] || rc=1
cmp -s "$D/serve/good.tsm" "$D/dest_b/good.tsm" || rc=1
[ ! -e "$D/dest_b/bad.tsm" ] && [ ! -e "$D/dest_b/bad.tsm.part" ] || rc=1
grep -q 'have  good.tsm' "$D/get_b2.out" || rc=1
[ "$(grep -c . "$D/urls_b")" -eq 2 ] || rc=1  # the second run fetched nothing
check "tesseract-get-models verifies, discards, skips" $rc \
      "exit $st1 then $st2, $(grep -c . "$D/urls_b") fetches"

# (c) tesseract-klebsiella --with-model installs kpneumoniae.tsm through tesseract-get-models and
# assembles with it. A copy of the helpers with a manifest for the fixture model stands in for
# an install; the fixture model is the one section 19 built.
KPFIX=$TMP/t19/models/kpneumoniae.tsm
if [ -s "$KPFIX" ]; then
    mkdir -p "$D/kit" "$D/serve_kp"
    cp "$ROOT/tesseract-klebsiella" "$ROOT/tesseract-get-models" "$D/kit/"
    ln -s "$TESSERACT" "$D/kit/tesseract-asm"
    cp "$KPFIX" "$D/serve_kp/kpneumoniae.tsm"
    echo "$(sha256sum < "$KPFIX" | cut -c1-64)  kpneumoniae.tsm" > "$D/kit/models.sha256"
    rm -f "$D/urls_c"
    TESSERACT_MODEL_DIR="$D/kfetch" PATH="$D/shim:$PATH" CURL_LOG="$D/urls_c" CURL_SERVE="$D/serve_kp" \
        n22user "$D/kit/tesseract-klebsiella" "$TMP/t19/r_1.fq.gz" "$TMP/t19/r_2.fq.gz" \
        -o "$D/kl_fetch" -t 2 --with-model >"$D/kl_fetch.out" 2>&1; st=$?
    if [ "$(cat "$D/urls_c" 2>/dev/null)" != "https://github.com/iowa69/TesserACT/releases/download/models-v2/kpneumoniae.tsm" ] ||
       ! cmp -s "$KPFIX" "$D/kfetch/kpneumoniae.tsm"; then
        fail "tesseract-klebsiella --with-model fetches it" "fetched: $(tr '\n' ' ' < "$D/urls_c" 2>/dev/null)"
    else
        modelran "tesseract-klebsiella --with-model fetches it" $st "$D/kl_fetch/r" kpneumoniae \
                 "$D/kl_fetch/r/run.log"
    fi
else
    fail "tesseract-klebsiella --with-model fetches it" "no fixture model from section 19"
fi

# (d) the non-interactive install: it installs the four commands and the manifest, downloads no
# model, and says the models are optional. `make` is a stand-in, so nothing is rebuilt here.
rm -f "$D/urls_d"; rm -rf "$D/inst"
env -u TESSERACT_MODEL_DIR HOME="$D/home" PATH="$D/mk:$D/shim:$PATH" CURL_LOG="$D/urls_d" \
    bash "$ROOT/install.sh" --prefix "$D/inst" </dev/null >"$D/install.out" 2>&1; st=$?
rc=0; miss=""
[ "$st" -eq 0 ] || rc=1
for f in tesseract-asm tesseract-klebsiella tesseract-eskape tesseract-get-models models.sha256; do
    [ -e "$D/inst/bin/$f" ] || miss="$miss $f"
done
[ -z "$miss" ] || rc=1
[ ! -e "$D/urls_d" ] && [ ! -e "$D/home/.tesseract" ] || rc=1
grep -q 'Organism models are optional' "$D/install.out" || rc=1
grep -q 'fetches the model' "$D/install.out" && rc=1
check "install.sh --prefix: no model, says optional" $rc "exit=$st${miss:+ missing:$miss}"

# (e) the help texts say a model is optional and what it costs
rc=0
"$TESSERACT" --help 2>&1 | grep -q 'roughly doubled misassemblies' || rc=1
"$ROOT/tesseract-eskape" --list | grep -q 'Models are optional and opt-in' || rc=1
"$ROOT/tesseract-klebsiella" --help | grep -q -- '--with-model' || rc=1
"$ROOT/tesseract-get-models" --help | grep -q 'The models are optional' || rc=1
check "help texts: models optional, with their cost" $rc "tesseract-asm, -eskape, -klebsiella, -get-models"

echo
echo "-----------------------------------------------------------------------"
printf '%d passed, %d failed, %d known open defects (XFAIL), %d skipped\n' \
       "$PASSED" "$FAILED" "$XFAILED" "$SKIPPED"
[ "$FAILED" -eq 0 ] || exit 1
exit 0
