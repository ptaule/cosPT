#!/usr/bin/env bash
#$ -S /bin/bash                            # shell to be used by SGE
#$ -cwd                                    # run in the Current Working Directory
#$ -l h_vmem=3G,h_cpu=0:59:00,h_fsize=2G   # hardware requirements - reserve enough memory!
#$ -M petter.taule@tum.de
#$ -m ae # a="MAIL_AT_ABORT", e="MAIL_AT_EXIT", combination `-m ae` is allowed

set -e
set -x

sourcedir=$(pwd)
tempdir="$(mktemp -d -t cosPT_XXXXXXXXXX)"

cleanup() {
    rm -rf "$tempdir"
    exit 1
}

# Run cleanup on EXIT
trap cleanup EXIT

exe=build/cosPT

# Set some default values:
k_a_idx=-1
k_b_idx=-1
k_c_idx=-1
n_evals=0
n_cores=4
ini_file=
log_file=

usage()
{
  echo "Usage: build_and_run [ --k_a_idx=INDEX ] [ --k_b_idx=INDEX ]
                     [ --k_c_idx=INDEX ] [ --n_evals=NUM ] [ --n_cores=NUM]
                     ini_file log_file"
  exit 2
}

parsed_arguments=$(getopt -n build_and_run -o a:b:c:n:p: \
    --long k_a_idx:,k_b_idx:,k_c_idx:,n_evals:,n_cores: -- "$@")
valid_arguments=$?
if [ "$valid_arguments" != "0" ]; then
  usage
fi

eval set -- "$parsed_arguments"
while :
do
  case "$1" in
    -a | --k_a_idx) k_a_idx="$2"   ; shift 2 ;;
    -b | --k_b_idx) k_b_idx="$2"   ; shift 2 ;;
    -c | --k_c_idx) k_c_idx="$2"   ; shift 2 ;;
    -n | --n_evals) n_evals=$2 ; shift 2 ;;
    -p | --n_cores) n_cores=$2 ; shift 2 ;;
    # -- means the end of the arguments; drop this, and break out of the while loop
    --) shift; break ;;
    # If invalid options were passed, then getopt should have reported an error,
    # which we checked as VALID_ARGUMENTS when getopt was called...
    *) echo "Unexpected option: $1 - this should not happen."
       usage ;;
  esac
done

if [ $# -ne 2 ]; then
  usage
fi

ini_file="$1"
log_file="$2"

cp -r -t "$tempdir" "$sourcedir"/src "$sourcedir"/include \
    "$sourcedir"/main.cpp "$sourcedir"/Makefile

cd "$sourcedir"
git_sha=$(git rev-parse HEAD)

cd "$tempdir"
# Create src/version_.cpp with git sha from source directory, and add this file
# to the source files listed in the makefile
printf "#include \"../include/version.hpp\"\nstd::string build_git_sha =  \"%s\";" \
    $git_sha > src/version_.cpp
sed -i '/SRC_FILES/ s/\\$/version_.cpp \\/' Makefile
make clean
make -j
export LD_LIBRARY_PATH=/space/ge52sir/local/lib/

"$tempdir"/$exe --k_a_idx $k_a_idx --k_b_idx $k_b_idx --k_c_idx $k_c_idx \
    --n_evals $n_evals --n_cores $n_cores "$ini_file" > "$log_file"
