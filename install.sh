#!/usr/bin/env bash
# Install TesserACT.
#
#   ./install.sh                 in a terminal, walks you through it step by step;
#                                anywhere else (a script, CI, a pipe) installs
#                                straight into the active conda env, or into
#                                ~/.local if conda is not active
#   ./install.sh --guided        force the step-by-step version
#   ./install.sh --no-prompt     force the silent version, never asks anything
#   ./install.sh --prefix DIR    install into DIR/bin
#   ./install.sh --conda-env NAME
#                                create (or reuse) that conda environment from
#                                conda/environment.yml, then install into it
#
# TesserACT needs only a C++17 compiler, zlib and pthreads. The conda environment
# is for the tools around it -- mash to choose a model panel, QUAST to score
# against a reference, Bandage to look at the graph.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
prefix=""
env_name=""
guided=""
nargs=$#

while [ $# -gt 0 ]; do
    case "$1" in
        --prefix)     prefix="${2:?--prefix needs a directory}"; shift 2 ;;
        --conda-env)  env_name="${2:?--conda-env needs a name}"; shift 2 ;;
        --guided)     guided=1; shift ;;
        --no-prompt|--yes) guided=0; shift ;;
        -h|--help)    sed -n '2,17p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *)            echo "unknown option: $1" >&2; exit 2 ;;
    esac
done

# Guided by default, but only when there is a person there to answer. A script, a CI job or
# `curl ... | bash` has no terminal on stdin, so it gets exactly the old silent behaviour and
# nothing ever blocks waiting for an answer that cannot come. Passing any flag also means the
# caller already knows what they want, so it stays silent then too.
if [ -z "$guided" ]; then
    if [ "$nargs" -eq 0 ] && [ -t 0 ] && [ -t 1 ]; then guided=1; else guided=0; fi
fi

say() { printf '\033[36m==>\033[0m %s\n' "$*"; }
die() { printf '\033[31merror:\033[0m %s\n' "$*" >&2; exit 1; }
ok()  { printf '  \033[32m[ok]\033[0m   %s\n' "$*"; }
bad() { printf '  \033[31m[--]\033[0m   %s\n' "$*"; }
step() { printf '\n\033[1m%s\033[0m\n' "$*"; }

# Reads one line from the terminal, not from stdin: stdin may be the pipe this script was
# fetched through, and reading the rest of the script as an answer is a memorable way to
# break an install. Returns the default on EOF or an empty line -- and `read` returning
# non-zero at EOF must not trip `set -e`, hence the `|| true`.
ask() {
    local promptline="$1" default="$2" reply=""
    printf '%s [%s]: ' "$promptline" "$default" > /dev/tty
    IFS= read -r reply < /dev/tty || true
    printf '%s' "${reply:-$default}"
}

yesno() {
    local promptline="$1" default="$2" reply
    reply="$(ask "$promptline" "$default")"
    # tr rather than ${reply,,}: that expansion is bash 4 and this script runs under whatever
    # /usr/bin/env bash finds, which on macOS is still 3.2 -- where it is a parse error, so it
    # would kill the script outright, guided or not.
    reply="$(printf '%s' "$reply" | tr '[:upper:]' '[:lower:]')"
    case "$reply" in y|yes) return 0 ;; *) return 1 ;; esac
}

# Plain ASCII on purpose. This is the first thing a new user sees, and it is seen over ssh,
# in a Windows terminal and in whatever locale the machine happens to have; box-drawing
# characters turn into mojibake in all three.
banner() {
cat <<'BANNER'

        #####################################
        ##                                 ##
        ##\  ###########################  /##
        ##   ##                       ##   ##
        ##   ##   T E S S E R A C T   ##   ##
        ##   ##                       ##   ##
        ##/  ###########################  \##
        ##                                 ##
        #####################################

  A short-read assembler for bacterial genomes.
  Giovanni Lorenzin, IOWA-BioTech

BANNER
}

[ "$guided" = 1 ] && banner

# ---- what this machine has -------------------------------------------------
# Checked up front and all at once. The alternative -- failing at the first missing tool --
# makes someone install one thing, re-run, and discover the next one; five minutes each time.
if [ "$guided" = 1 ]; then
    step "Step 1 of 5: checking what this computer already has"
    missing=()
    for tool in g++ make curl sha256sum; do
        if command -v "$tool" >/dev/null 2>&1; then ok "$tool"; else bad "$tool -- missing"; missing+=("$tool"); fi
    done
    if [ "${#missing[@]}" -gt 0 ]; then
        printf '\n'
        say "Some tools are missing. On Ubuntu or Debian, copy this line and run it:"
        printf '\n      sudo apt install build-essential curl coreutils\n\n'
        say "On a Mac: xcode-select --install"
        die "install those first, then run ./install.sh again"
    fi
fi

# ---- create the environment if asked -------------------------------------
if [ -n "$env_name" ]; then
    command -v conda >/dev/null 2>&1 || die "conda not found on PATH"
    # Captured and matched in the shell rather than piped into `grep -q`. Under
    # `set -o pipefail` -- which this script sets -- grep -q exits at the first match and
    # closes the pipe, the producer takes SIGPIPE and exits 141, and the pipeline as a whole
    # reports failure even though the match succeeded. The `if` would then take the wrong
    # branch and try to create an environment that already exists, which conda refuses,
    # which `set -e` turns into a dead install. It needs a long list to fire -- around ten
    # thousand lines, so not at any plausible number of conda environments -- but the
    # construct is wrong regardless and costs nothing to remove.
    existing=$(conda env list | awk '{print $1}')
    if [[ $'\n'"$existing"$'\n' == *$'\n'"$env_name"$'\n'* ]]; then
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
        default_prefix="$CONDA_PREFIX"
    else
        default_prefix="$HOME/.local"
    fi
    if [ "$guided" = 1 ]; then
        step "Step 2 of 5: where to put it"
        if [ -n "${CONDA_PREFIX:-}" ]; then
            printf '  You have a conda environment active, so that is the natural home.\n'
        else
            printf '  No conda environment is active, so this goes in your own home directory.\n'
            printf '  Nothing is installed system-wide and no password is needed.\n'
        fi
        printf '\n'
        prefix="$(ask "  Install into" "$default_prefix")"
    else
        prefix="$default_prefix"
        if [ -n "${CONDA_PREFIX:-}" ]; then
            say "installing into the active conda environment: $prefix"
        else
            say "no conda environment active; installing into $prefix"
        fi
    fi
fi

# ---- check the compiler ----------------------------------------------------
cxx="${CXX:-g++}"
command -v "$cxx" >/dev/null 2>&1 || die "no C++ compiler found (set CXX, or install g++)"
if ! echo 'int main(){return 0;}' | "$cxx" -std=c++17 -x c++ - -o /dev/null 2>/dev/null; then
    die "$cxx does not accept -std=c++17"
fi

# ---- build -----------------------------------------------------------------
[ "$guided" = 1 ] && step "Step 3 of 5: building the assembler (a minute or two)"
say "building with $cxx"
make -C "$here" -j"$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)" CXX="$cxx"

# ---- install ---------------------------------------------------------------
[ "$guided" = 1 ] && step "Step 4 of 5: installing the commands"
mkdir -p "$prefix/bin"
install -m 0755 "$here/tesseract-asm"        "$prefix/bin/tesseract-asm"
install -m 0755 "$here/tesseract-model"      "$prefix/bin/tesseract-model"
install -m 0755 "$here/tesseract-klebsiella" "$prefix/bin/tesseract-klebsiella"
# tesseract-eskape and tesseract-get-models are the whole ESKAPEE path: one command per
# organism, and the fetcher that puts the models where it looks for them. Leaving them out
# of the install -- as this script did until now -- meant a finished install where the
# documented commands did not exist, and no error saying so.
install -m 0755 "$here/tesseract-eskape"     "$prefix/bin/tesseract-eskape"
install -m 0755 "$here/tesseract-get-models" "$prefix/bin/tesseract-get-models"
# tesseract-get-models reads the checksum list from beside itself, so the list has to travel
# with it or every model it downloads fails verification and is deleted.
install -m 0644 "$here/models.sha256"        "$prefix/bin/models.sha256"
say "installed tesseract-asm, tesseract-model, tesseract-klebsiella, tesseract-eskape and tesseract-get-models into $prefix/bin"

# ---- verify ----------------------------------------------------------------
if ! "$prefix/bin/tesseract-asm" --version >/dev/null 2>&1; then
    die "the installed binary does not run"
fi
say "$("$prefix/bin/tesseract-asm" --version)"

on_path=1
case ":$PATH:" in
    *":$prefix/bin:"*) ;;
    *) on_path=0 ;;
esac

if [ "$on_path" = 0 ]; then
    if [ "$guided" = 1 ]; then
        rc="$HOME/.bashrc"
        [ -n "${ZSH_VERSION:-}" ] && rc="$HOME/.zshrc"
        case "${SHELL:-}" in */zsh) rc="$HOME/.zshrc" ;; esac
        printf '\n'
        printf '  The commands are installed, but your shell does not know where to find them yet.\n'
        printf '  This is one line added to %s.\n\n' "$rc"
        if yesno "  Add it for me?" "yes"; then
            printf '\n# added by TesserACT install.sh\nexport PATH="%s/bin:$PATH"\n' "$prefix" >> "$rc"
            ok "added to $rc -- it applies to new terminals, or run: source $rc"
        else
            printf '\n  Then add this yourself, to %s:\n      export PATH="%s/bin:$PATH"\n' "$rc" "$prefix"
        fi
    else
        printf '\n\033[33mnote:\033[0m %s is not on your PATH. Add it with:\n      export PATH="%s/bin:$PATH"\n' "$prefix/bin" "$prefix"
    fi
fi

# ---- models ----------------------------------------------------------------
if [ "$guided" = 1 ]; then
    step "Step 5 of 5: the organism models"
    printf '  TesserACT assembles without these. They are what lets it lay the pieces out\n'
    printf '  against a closed genome of the same species, for the seven clinical bugs:\n'
    printf '  Klebsiella, E. coli, Enterobacter, Acinetobacter, Pseudomonas, S. aureus,\n'
    printf '  Enterococcus. About 150 MB in total, downloaded once.\n\n'
    if yesno "  Download them now?" "yes"; then
        printf '\n'
        # Not fatal. A failed download is a network problem, not an install problem: the
        # assembler is already installed and working, and the fetch is one command to retry.
        if "$prefix/bin/tesseract-get-models"; then
            ok "models are in ${TESSERACT_MODEL_DIR:-$HOME/.tesseract/models}"
        else
            printf '\n'
            say "the download did not finish. Nothing is broken -- retry any time with:"
            printf '      tesseract-get-models\n'
        fi
    else
        say "skipped. When you want them: tesseract-get-models"
    fi
fi

if [ "$guided" = 1 ]; then
cat <<EOF

$(printf '\033[1mDone.\033[0m') Assemble a genome like this:

  tesseract-eskape --preset kpneumoniae -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz -o my_result

Your assembled genome is then the file  my_result/contigs.fasta

  tesseract-eskape --list     the preset name for each organism
  tesseract-eskape --help     every option, explained
EOF
else
cat <<EOF

Next:
  tesseract-asm -1 reads_1.fq.gz -2 reads_2.fq.gz -o assembly

For Klebsiella there is nothing to set up -- this fetches the model and runs everything:
  ./tesseract-klebsiella reads/          every read pair in the directory
  ./tesseract-klebsiella R1.fq.gz        the mate is found automatically

For the other six ESKAPEE organisms:
  tesseract-get-models                   fetch the models (about 150 MB, once)
  tesseract-eskape --list                the preset name for each organism

To build a model of your own instead:
  tesseract-model --organism klebsiella --out kleb.tsm references/*.fasta
  tesseract-asm --organism klebsiella --model kleb.tsm -1 R1.fq.gz -2 R2.fq.gz -o out
EOF
fi
