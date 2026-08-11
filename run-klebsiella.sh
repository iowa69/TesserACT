#!/usr/bin/env bash
# run-klebsiella.sh -- assemble Klebsiella reads and say what each contig is, in one command.
#
#   ./run-klebsiella.sh reads/                     every read pair in the directory
#   ./run-klebsiella.sh sample_R1.fq.gz            the mate is found automatically
#   ./run-klebsiella.sh sample_R1.fq.gz sample_R2.fq.gz
#
# Builds the assembler if it is not built, downloads the model if it is not cached, checks
# what it downloaded, assembles, and leaves one FASTA per isolate with every contig labelled
# chromosome, plasmid or unknown. Nothing else to run and nothing to configure.
#
# Everything it does is idempotent: interrupt it and run it again and it picks up where it
# stopped, because each isolate is assembled into a scratch directory and moved into place
# only when it finishes. That is also why a killed run never leaves a half-written FASTA that
# a later run would mistake for a finished one.
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Installed as `tessera-klebsiella`, run from a checkout as `./run-klebsiella.sh`. The help
# should name whichever one the user actually typed.
SELF="$(basename "${BASH_SOURCE[0]}")"
# From a checkout you type ./run-klebsiella.sh; installed on PATH you just type the command.
# The help should show whichever is actually true.
case "$SELF" in run-klebsiella.sh) INVOKE="./$SELF" ;; *) INVOKE="$SELF" ;; esac
REPO=iowa69/TesserACT
# The release the MODELS are attached to. Deliberately separate from the assembler's own
# version: a patch release of the binary does not re-upload a 1.3 GB model, so this pin
# moves only when the models actually change.
MODEL_RELEASE=v1.2.0

# The models live outside the repository -- they are 339 MB and 2.9 GB. Cached per user so a
# second project does not download them again.
CACHE="${TESSERA_MODEL_DIR:-${XDG_CACHE_HOME:-$HOME/.cache}/tessera}"

# Pinned: a download is checked against these before it is used. A truncated model otherwise
# fails much later with a confusing complaint about the file not being a model.
SHA_DEFAULT=cf40b986899a50ddf7b888e75762d1764b646e6d264dde13b041b6aadc9a0606
SHA_PLASMID=308a097f2c39bd5961a00ed1d603419109ad462f747578a503016907215b5880

OUT=tessera-out
THREADS=$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)
JOBS=1
MODEL=""
WANT=default
KEEP_GOING=0

c_info=$'\033[36m'; c_ok=$'\033[32m'; c_warn=$'\033[33m'; c_err=$'\033[31m'; c_off=$'\033[0m'
[ -t 1 ] || { c_info=""; c_ok=""; c_warn=""; c_err=""; c_off=""; }
say()  { printf '%s==>%s %s\n' "$c_info" "$c_off" "$*"; }
ok()   { printf '%s  ok%s %s\n' "$c_ok" "$c_off" "$*"; }
warn() { printf '%swarn:%s %s\n' "$c_warn" "$c_off" "$*" >&2; }
die()  { printf '%serror:%s %s\n' "$c_err" "$c_off" "$*" >&2; exit 1; }

usage() {
    sed -n '2,12p' "${BASH_SOURCE[0]}" | sed -e 's/^# \{0,1\}//' -e "s|\./run-klebsiella\.sh|$INVOKE|g" -e "s|run-klebsiella\.sh|$SELF|g"
    cat <<EOF

Options:
  -o DIR        where to put the results        (default: $OUT)
  -t N          threads per isolate             (default: every core, here $THREADS)
  -j N          isolates at once                (default: 1)
  --plasmid     use the plasmid-oriented model instead of the default
  --model FILE  use this model file; skips the download entirely
  --keep-going  carry on if one isolate fails
  -h, --help    this

Which model: take the default unless you specifically want plasmid grouping. It is the
better chromosome model. See models/README.md for the numbers behind that sentence.
EOF
}

INPUTS=()
while [ $# -gt 0 ]; do
    case "$1" in
        -o) OUT="${2:?-o needs a directory}"; shift 2 ;;
        -t) THREADS="${2:?-t needs a number}"; shift 2 ;;
        -j) JOBS="${2:?-j needs a number}"; shift 2 ;;
        --plasmid) WANT=plasmid; shift ;;
        --model) MODEL="${2:?--model needs a file}"; shift 2 ;;
        --keep-going) KEEP_GOING=1; shift ;;
        -h|--help) usage; exit 0 ;;
        -*) die "unknown option: $1   (try --help)" ;;
        *) INPUTS+=("$1"); shift ;;
    esac
done

[ ${#INPUTS[@]} -gt 0 ] || { usage; exit 2; }
case "$THREADS" in ''|*[!0-9]*) die "-t wants a number, got '$THREADS'" ;; esac
case "$JOBS"    in ''|*[!0-9]*) die "-j wants a number, got '$JOBS'" ;; esac
[ "$THREADS" -ge 1 ] || die "-t must be at least 1"
[ "$JOBS" -ge 1 ]    || die "-j must be at least 1"

# ---------------------------------------------------------------------------
# 1. the assembler
# ---------------------------------------------------------------------------
TESSERA=""
for cand in "$HERE/tessera" "$(command -v tessera 2>/dev/null || true)"; do
    [ -n "$cand" ] && [ -x "$cand" ] && { TESSERA=$cand; break; }
done
if [ -z "$TESSERA" ]; then
    say "no tessera binary yet, building it"
    command -v make >/dev/null || die "make is not installed, so the assembler cannot be built"
    cxx="${CXX:-g++}"
    command -v "$cxx" >/dev/null || die "no C++ compiler found (install g++, or set CXX)"
    echo 'int main(){return 0;}' | "$cxx" -std=c++17 -x c++ - -o /dev/null 2>/dev/null \
        || die "$cxx cannot compile C++17. Install a newer g++ or clang."
    if ! echo '#include <zlib.h>
int main(){return 0;}' | "$cxx" -x c++ - -o /dev/null -lz 2>/dev/null; then
        die "zlib headers are missing.
      Debian/Ubuntu:  sudo apt install zlib1g-dev
      Fedora/RHEL:    sudo dnf install zlib-devel
      conda:          conda install -c conda-forge zlib   (then CPATH=\$CONDA_PREFIX/include)"
    fi
    make -C "$HERE" -j"$THREADS" CXX="$cxx" >/dev/null || die "the build failed; run 'make' in $HERE to see why"
    TESSERA="$HERE/tessera"
    [ -x "$TESSERA" ] || die "the build reported success but produced no binary"
fi
ok "assembler: $("$TESSERA" --version 2>/dev/null || echo "$TESSERA")"

# ---------------------------------------------------------------------------
# 2. work out what to assemble
# ---------------------------------------------------------------------------
# Mate naming in the wild: _R1/_R2, _1/_2, .R1./.R2., any of .fq .fastq, gzipped or not.
mate_of() {
    local r1=$1 r2
    for pat in 'R1:R2' '_1:_2' '.1.:.2.'; do
        local a=${pat%%:*} b=${pat##*:}
        case "$r1" in
            *"$a"*) r2=$(printf '%s' "$r1" | sed "s|\(.*\)$a|\1$b|"); [ -s "$r2" ] && { printf '%s' "$r2"; return 0; } ;;
        esac
    done
    return 1
}
sample_of() {
    basename "$1" \
      | sed -E 's/\.(fastq|fq)(\.gz)?$//; s/[._]R?1$//; s/[._](R1|read1|1)([._].*)?$//'
}

declare -a S_NAME=() S_R1=() S_R2=()
add_pair() {
    local r1=$1 r2=$2 name; name=$(sample_of "$r1")
    for existing in ${S_NAME[@]+"${S_NAME[@]}"}; do
        [ "$existing" = "$name" ] && return 0     # same sample reached two ways
    done
    S_NAME+=("$name"); S_R1+=("$r1"); S_R2+=("$r2")
}

if [ ${#INPUTS[@]} -eq 2 ] && [ -f "${INPUTS[0]}" ] && [ -f "${INPUTS[1]}" ]; then
    add_pair "${INPUTS[0]}" "${INPUTS[1]}"
else
    for inp in "${INPUTS[@]}"; do
        if [ -d "$inp" ]; then
            found=0
            while IFS= read -r r1; do
                if r2=$(mate_of "$r1"); then add_pair "$r1" "$r2"; found=$((found+1)); fi
            done < <(find "$inp" -maxdepth 1 -type f \
                        \( -name '*_R1*.f*q*' -o -name '*_1.f*q*' -o -name '*.R1.f*q*' \) | sort)
            if [ "$found" -eq 0 ]; then
                warn "no read pairs found in $inp"
                echo "      it contains:" >&2
                ls -1 "$inp" 2>/dev/null | head -8 | sed 's/^/        /' >&2
                echo "      expected names like  SAMPLE_R1.fastq.gz  and  SAMPLE_R2.fastq.gz" >&2
            fi
        elif [ -f "$inp" ]; then
            if r2=$(mate_of "$inp"); then
                add_pair "$inp" "$r2"
            else
                die "found $inp but not its mate.
      Looked for the same name with R2 (or _2) in place of R1 (or _1).
      Pass both files explicitly if they are named differently."
            fi
        else
            die "no such file or directory: $inp"
        fi
    done
fi

[ ${#S_NAME[@]} -gt 0 ] || die "nothing to assemble"
say "${#S_NAME[@]} isolate(s) to assemble, ${THREADS} thread(s) each, ${JOBS} at a time"

# ---------------------------------------------------------------------------
# 3. the model, once the inputs are known to be good
# ---------------------------------------------------------------------------
fetch() {   # fetch <url> <destination>
    # A progress bar is worth having on a terminal and is 280 kB of noise in a log file, so
    # it depends on where the output is going.
    local q=(--progress-bar); [ -t 1 ] || q=(--silent --show-error)
    if command -v curl >/dev/null; then
        curl -fL --retry 3 --retry-delay 2 "${q[@]}" -o "$2" "$1"
    elif command -v wget >/dev/null; then
        local w=(--show-progress); [ -t 1 ] || w=()
        wget -q "${w[@]}" -O "$2" "$1"
    else
        die "neither curl nor wget is installed, so the model cannot be downloaded.
      Download it by hand from https://github.com/$REPO/releases/tag/$MODEL_RELEASE
      and pass it with --model"
    fi
}

sha_of() {
    if command -v sha256sum >/dev/null; then sha256sum "$1" | cut -d' ' -f1
    elif command -v shasum   >/dev/null; then shasum -a 256 "$1" | cut -d' ' -f1
    else echo ""; fi
}

if [ -z "$MODEL" ]; then
    if [ "$WANT" = plasmid ]; then
        asset="tessera-klebsiella-plasmid-$MODEL_RELEASE.tsm.zst"; want_sha=$SHA_PLASMID
        MODEL="$CACHE/tessera-klebsiella-plasmid-$MODEL_RELEASE.tsm"
    else
        asset="tessera-klebsiella-default-$MODEL_RELEASE.tsm"; want_sha=$SHA_DEFAULT
        MODEL="$CACHE/$asset"
    fi

    if [ ! -s "$MODEL" ]; then
        mkdir -p "$CACHE" || die "cannot create the model cache at $CACHE"
        need_mb=$([ "$WANT" = plasmid ] && echo 4300 || echo 400)
        free_mb=$(df -Pm "$CACHE" 2>/dev/null | awk 'NR==2{print $4}')
        if [ -n "$free_mb" ] && [ "$free_mb" -lt "$need_mb" ]; then
            die "not enough space in $CACHE: need about ${need_mb} MB, have ${free_mb} MB.
      Point TESSERA_MODEL_DIR somewhere roomier and run again."
        fi
        url="https://github.com/$REPO/releases/download/$MODEL_RELEASE/$asset"
        say "downloading the $WANT model once into $CACHE"
        tmp="$CACHE/.$asset.part"
        rm -f "$tmp"
        fetch "$url" "$tmp" || { rm -f "$tmp"; die "the download failed. Check the network, or fetch it by hand:
      $url
      then pass it with --model"; }

        got=$(sha_of "$tmp")
        if [ -z "$got" ]; then
            warn "no sha256 tool available, so the download could not be checked"
        elif [ "$got" != "$want_sha" ]; then
            rm -f "$tmp"
            die "the download is corrupt (checksum does not match). Run again; if it keeps
      happening, fetch it by hand from $url"
        else
            ok "checksum verified"
        fi

        if [ "$WANT" = plasmid ]; then
            command -v zstd >/dev/null || { rm -f "$tmp"; die "zstd is needed to expand the plasmid model.
      Debian/Ubuntu: sudo apt install zstd     conda: conda install -c conda-forge zstd
      (the default model needs no such step -- drop --plasmid to use it)"; }
            say "expanding it (2.9 GB once expanded)"
            zstd -d -q -f "$tmp" -o "$MODEL.part" || { rm -f "$tmp" "$MODEL.part"; die "could not expand the model"; }
            mv "$MODEL.part" "$MODEL"; rm -f "$tmp"
        else
            mv "$tmp" "$MODEL"
        fi
    fi
fi
[ -s "$MODEL" ] || die "model file is missing or empty: $MODEL"
ok "model: $MODEL"


# ---------------------------------------------------------------------------
# 4. assemble
# ---------------------------------------------------------------------------
mkdir -p "$OUT" || die "cannot create the output directory $OUT"
SUMMARY="$OUT/summary.tsv"

run_one() {
    local name=$1 r1=$2 r2=$3 d="$OUT/$name"
    if [ ! -s "$d/contigs.fasta" ]; then
        rm -rf "$d.part"; mkdir -p "$d.part"
        if ! "$TESSERA" -1 "$r1" -2 "$r2" -o "$d.part" -t "$THREADS" \
                --organism klebsiella --model "$MODEL" > "$d.part/run.log" 2>&1 \
           || [ ! -s "$d.part/contigs.fasta" ]; then
            printf '%sfailed:%s %s -- see %s\n' "$c_err" "$c_off" "$name" "$d.part/run.log" >&2
            return 1
        fi
        mv "$d.part" "$d"
    fi
    awk -v s="$name" '
      /^>/ { n++
             if ($0 ~ /_chr/)        chr++
             else if ($0 ~ /_plas/) { plas++; if (match($0, /_plas_[0-9]+/)) g[substr($0,RSTART,RLENGTH)]=1 }
             else if ($0 ~ /_unk/)   unk++
             if ($0 ~ /_circular/)   circ++
             if (len > max) max = len; tot += len; len = 0; next }
      { len += length($0) }
      END { if (len > max) max = len; tot += len
            ng = 0; for (k in g) ng++
            printf "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", s, n, chr, plas, ng, unk, circ, max, tot }
    ' "$d/contigs.fasta"
}

printf 'sample\tcontigs\tchr\tplas\tplasmid_groups\tunk\tcircular\tlongest\ttotal\n' > "$SUMMARY"
failed=0; done_n=0
i=0
while [ $i -lt ${#S_NAME[@]} ]; do
    batch=0
    while [ $batch -lt "$JOBS" ] && [ $i -lt ${#S_NAME[@]} ]; do
        name=${S_NAME[$i]}
        printf '%s  ..%s %s\n' "$c_info" "$c_off" "$name"
        run_one "$name" "${S_R1[$i]}" "${S_R2[$i]}" >> "$SUMMARY" &
        i=$((i+1)); batch=$((batch+1))
    done
    for _ in $(seq 1 $batch); do
        if wait -n; then done_n=$((done_n+1)); else
            failed=$((failed+1))
            [ "$KEEP_GOING" -eq 1 ] || { warn "stopping after a failure (pass --keep-going to carry on)"; break 2; }
        fi
    done
done

# Sort the summary so the order does not depend on which isolate finished first.
if [ -s "$SUMMARY" ]; then
    { head -1 "$SUMMARY"; tail -n +2 "$SUMMARY" | sort; } > "$SUMMARY.sorted" && mv "$SUMMARY.sorted" "$SUMMARY"
fi

# ---------------------------------------------------------------------------
# 5. say what happened
# ---------------------------------------------------------------------------
echo
if [ "$failed" -gt 0 ]; then
    warn "$failed isolate(s) failed; $done_n finished"
else
    ok "$done_n isolate(s) assembled"
fi
cat <<EOF

Results are in $OUT/
  <sample>/contigs.fasta      the assembly, every contig labelled
  <sample>/report.html        what the run decided, rung by rung
  summary.tsv                 one line per isolate

Contig names end in what the contig is:
  _chr        chromosomal
  _plas       plasmid, which molecule is unknown
  _plas_<n>   plasmid, grouped with the others carrying the same <n>
  _unk        no signal reached it -- reported rather than guessed at
  _circular   the contig's two ends join: it is the whole molecule

The file is ordered as a genome: chromosome first, then each plasmid molecule
whole and contiguous, then the unassigned.
EOF
[ "$failed" -eq 0 ] || exit 1
exit 0
