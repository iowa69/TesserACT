#!/usr/bin/env bash
# Install tessera.
#
#   ./install.sh                 build and install into the active conda env,
#                                or into ~/.local/bin if conda is not active
#   ./install.sh --prefix DIR    install into DIR/bin
#   ./install.sh --conda-env NAME
#                                create (or reuse) that conda environment from
#                                conda/environment.yml, then install into it
#
# tessera needs only a C++17 compiler, zlib and pthreads. The conda environment
# is for the tools around it -- mash to choose a model panel, QUAST to score
# against a reference, Bandage to look at the graph.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
prefix=""
env_name=""

while [ $# -gt 0 ]; do
    case "$1" in
        --prefix)     prefix="${2:?--prefix needs a directory}"; shift 2 ;;
        --conda-env)  env_name="${2:?--conda-env needs a name}"; shift 2 ;;
        -h|--help)    sed -n '2,16p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *)            echo "unknown option: $1" >&2; exit 2 ;;
    esac
done

say() { printf '\033[36m==>\033[0m %s\n' "$*"; }
die() { printf '\033[31merror:\033[0m %s\n' "$*" >&2; exit 1; }

# ---- create the environment if asked -------------------------------------
if [ -n "$env_name" ]; then
    command -v conda >/dev/null 2>&1 || die "conda not found on PATH"
    if conda env list | awk '{print $1}' | grep -qx "$env_name"; then
        say "reusing conda environment '$env_name'"
    else
        say "creating conda environment '$env_name'"
        conda env create -n "$env_name" -f "$here/conda/environment.yml"
    fi
    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "$env_name"
fi

# ---- where to install ------------------------------------------------------
if [ -z "$prefix" ]; then
    if [ -n "${CONDA_PREFIX:-}" ]; then
        prefix="$CONDA_PREFIX"
        say "installing into the active conda environment: $prefix"
    else
        prefix="$HOME/.local"
        say "no conda environment active; installing into $prefix"
    fi
fi

# ---- check the compiler ----------------------------------------------------
cxx="${CXX:-g++}"
command -v "$cxx" >/dev/null 2>&1 || die "no C++ compiler found (set CXX, or install g++)"
if ! echo 'int main(){return 0;}' | "$cxx" -std=c++17 -x c++ - -o /dev/null 2>/dev/null; then
    die "$cxx does not accept -std=c++17"
fi

# ---- build -----------------------------------------------------------------
say "building with $cxx"
make -C "$here" -j"$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)" CXX="$cxx"

# ---- install ---------------------------------------------------------------
mkdir -p "$prefix/bin"
install -m 0755 "$here/tessera" "$prefix/bin/tessera"
install -m 0755 "$here/tessera-model" "$prefix/bin/tessera-model"
say "installed tessera and tessera-model into $prefix/bin"

# ---- verify ----------------------------------------------------------------
if ! "$prefix/bin/tessera" --version >/dev/null 2>&1; then
    die "the installed binary does not run"
fi
say "$("$prefix/bin/tessera" --version)"

case ":$PATH:" in
    *":$prefix/bin:"*) ;;
    *) printf '\n\033[33mnote:\033[0m %s is not on your PATH. Add it with:\n      export PATH="%s/bin:$PATH"\n' "$prefix/bin" "$prefix" ;;
esac

cat <<EOF

Next:
  tessera -1 reads_1.fq.gz -2 reads_2.fq.gz -o assembly

For Klebsiella, build a model once and reuse it:
  tessera-model --organism klebsiella --out kleb.tsm references/*.fasta
  tessera --organism klebsiella --model kleb.tsm -1 R1.fq.gz -2 R2.fq.gz -o out
EOF
