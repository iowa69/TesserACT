#!/usr/bin/env bash
# Every flag the help text and README advertise must exist, parse, and take
# effect. A flag that is documented but silently ignored is worse than one that
# is missing, because nothing fails.
set -uo pipefail
BIN=${1:-./tessera}
MODELBIN=${2:-./tessera-model}
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
pass=0; fail=0

ok()   { printf "  ok    %s\n" "$1"; pass=$((pass+1)); }
bad()  { printf "  FAIL  %s -- %s\n" "$1" "$2"; fail=$((fail+1)); }

# --- a tiny but real dataset ------------------------------------------------
python3 - "$TMP" <<'PY'
import random, sys, gzip, os
random.seed(7)
d = sys.argv[1]
g = "".join(random.choice("ACGT") for _ in range(60000))
rep = "".join(random.choice("ACGT") for _ in range(800))
g = g[:20000] + rep + g[20000:40000] + rep + g[40000:]
open(os.path.join(d, "ref.fa"), "w").write(">ref\n" + g + "\n")
def rc(s): return s[::-1].translate(str.maketrans("ACGT", "TGCA"))
r1, r2 = [], []
for i in range(12000):
    p = random.randint(0, len(g) - 400)
    frag = g[p:p + random.randint(280, 360)]
    r1.append(frag[:150]); r2.append(rc(frag[-150:]))
for name, rs in (("R1", r1), ("R2", r2)):
    with gzip.open(os.path.join(d, f"{name}.fq.gz"), "wt") as fh:
        for i, s in enumerate(rs):
            fh.write(f"@r{i}\n{s}\n+\n{'I'*len(s)}\n")
PY

R1="$TMP/R1.fq.gz"; R2="$TMP/R2.fq.gz"

run () {   # run OUTDIR extra-args...
  local out="$1"; shift
  "$BIN" -1 "$R1" -2 "$R2" -o "$out" -t 4 "$@" > "$out.log" 2>&1
}

# --- 1. flags parse and the run succeeds ------------------------------------
for spec in \
  "--mode fast" "--mode standard" "--mode careful" "--mode aggressive" \
  "--aggressive" "--no-correct" "--no-resolve" "--no-scaffold" "--no-polish" \
  "--min-contig 500" "-k 21,33" "-k 21,55,127" "-c 2" "--min-link 3" "--tie-ratio 1.3" \
  "--quiet" "--map-polish none"
do
  d="$TMP/f$(echo "$spec" | tr -d ' /,.-')"
  mkdir -p "$d"
  if run "$d" $spec && [ -s "$d/contigs.fasta" ]; then ok "$spec"; else bad "$spec" "run failed"; fi
done

# --- 2. flags that must be rejected ----------------------------------------
for spec in "-k 4" "-k 200" "--mode nonsense" "--map-polish nonsense" "--nosuchflag" \
            "--tie-ratio 0" "--min-link -5" "-t 0" "--simplify-rounds 0" \
            "--polish-passes -1" "--qtrim-quality 99" "--tie-ratio abc"; do
  d="$TMP/r$(echo "$spec" | tr -d ' /,.-')"; mkdir -p "$d"
  if run "$d" $spec; then bad "reject $spec" "accepted an invalid value"; else ok "reject $spec"; fi
done

# --- 3. flags must actually take effect ------------------------------------
mkdir -p "$TMP/base" "$TMP/minc"
run "$TMP/base"; run "$TMP/minc" --min-contig 100000
a=$(grep -c '^>' "$TMP/base/contigs.fasta" 2>/dev/null); a=${a:-0}
b=$(grep -c '^>' "$TMP/minc/contigs.fasta" 2>/dev/null); b=${b:-0}
[ "$b" -lt "$a" ] && ok "--min-contig changes output ($a -> $b)" || bad "--min-contig" "no effect ($a -> $b)"

[ -s "$TMP/base/assembly_graph.gfa" ] && ok "GFA written by default" || bad "GFA" "missing"
[ -s "$TMP/base/report.html" ] && ok "HTML report written" || bad "report.html" "missing"
[ -s "$TMP/base/report.json" ] && ok "JSON report written" || bad "report.json" "missing"

# determinism
mkdir -p "$TMP/d1" "$TMP/d2"
run "$TMP/d1"; run "$TMP/d2"
cmp -s "$TMP/d1/contigs.fasta" "$TMP/d2/contigs.fasta" \
  && ok "deterministic across runs" || bad "determinism" "outputs differ"

# thread invariance
mkdir -p "$TMP/t1" "$TMP/t8"
"$BIN" -1 "$R1" -2 "$R2" -o "$TMP/t1" -t 1 >/dev/null 2>&1
"$BIN" -1 "$R1" -2 "$R2" -o "$TMP/t8" -t 8 >/dev/null 2>&1
cmp -s "$TMP/t1/contigs.fasta" "$TMP/t8/contigs.fasta" \
  && ok "thread-invariant" || bad "thread invariance" "outputs differ"

# --- 4. the model path ------------------------------------------------------
if [ -x "$MODELBIN" ]; then
  mkdir -p "$TMP/mdl"
  cp "$TMP/ref.fa" "$TMP/mdl/A_chr.fasta"
  python3 - "$TMP" <<'PY'
import random, sys, os
random.seed(11)
d = sys.argv[1]
g = open(os.path.join(d, "ref.fa")).read().split("\n", 1)[1].replace("\n", "")
for n in "BCDEFGHIJ":
    s = list(g)
    for _ in range(len(g) // 200):        # a few percent divergence
        i = random.randrange(len(s)); s[i] = random.choice("ACGT")
    open(os.path.join(d, "mdl", f"{n}_chr.fasta"), "w").write(f">{n}\n" + "".join(s) + "\n")
PY
  if "$MODELBIN" --organism test --out "$TMP/test.tsm" --min-support 2 "$TMP"/mdl/*.fasta \
       > "$TMP/model.log" 2>&1 && [ -s "$TMP/test.tsm" ]; then
    ok "tessera-model builds a model"
  else
    bad "tessera-model" "did not produce a model"
  fi
  if grep -q "excluded" "$TMP/model.log"; then ok "model reports exclusions"; else bad "model exclusions" "not reported"; fi

  "$MODELBIN" --organism test --out "$TMP/excl.tsm" --exclude A --min-support 2 \
      "$TMP"/mdl/*.fasta > "$TMP/excl.log" 2>&1
  grep -q "1 accessions excluded" "$TMP/excl.log" \
    && ok "--exclude is honoured and recorded" || bad "--exclude" "not recorded"

  mkdir -p "$TMP/withmodel"
  if run "$TMP/withmodel" --organism test --model "$TMP/test.tsm" \
     && grep -q "model joining" "$TMP/withmodel.log"; then
    ok "--model runs the model stage"
  else
    bad "--model" "stage did not run"
  fi

  # A model the user named and tessera cannot read is fatal. It used to warn
  # and carry on, and the warning only appeared under --verbose, so a --quiet
  # run returned a vanilla assembly while the user believed the model had
  # guided it -- and the two differ in exactly the junctions the model exists
  # to settle.
  mkdir -p "$TMP/badmodel"
  if run "$TMP/badmodel" --model "$TMP/does-not-exist.tsm"; then
    bad "missing model" "should have failed, not assembled without the model"
  elif grep -qi "cannot open model file" "$TMP/badmodel.log"; then
    ok "an unreadable --model is a fatal error"
  else
    bad "missing model" "failed without saying the model could not be read"
  fi

  # IS panel
  head -2 "$TMP/ref.fa" > "$TMP/is.fna"
  mkdir -p "$TMP/withis"
  if run "$TMP/withis" --organism test --model "$TMP/test.tsm" --is-panel "$TMP/is.fna" \
     && grep -q "insertion-sequence panel" "$TMP/withis.log"; then
    ok "--is-panel loads"
  else
    bad "--is-panel" "not loaded"
  fi
fi

printf "\n%d passed, %d failed\n" "$pass" "$fail"
[ "$fail" -eq 0 ]
