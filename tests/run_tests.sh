#!/usr/bin/env bash

set -e
set -x

cleanup() {
    rm -rf "$tempdir"
}

# Run cleanup on EXIT
trap cleanup EXIT

# Set some default values:
n_cores=1
log_file=

usage()
{
  echo "Usage: run_tests.sh [ --n_cores NUM] LOG_FILE"
  exit 2
}

PARSED_ARGUMENTS=$(getopt -n run_tests -o hp: --long help,n_cores: -- "$@")
VALID_ARGUMENTS=$?
if [ "$VALID_ARGUMENTS" != "0" ]; then
  usage
fi

eval set -- "$PARSED_ARGUMENTS"
while :
do
  case "$1" in
    -h | --help) usage ; shift 2 ;;
    -p | --n_cores) n_cores=$2 ; shift 2 ;;
    # -- means the end of the arguments; drop this, and break out of the while loop
    --) shift; break ;;
    # If invalid options were passed, then getopt should have reported an error,
    # which we checked as VALID_ARGUMENTS when getopt was called...
    *) echo "Unexpected option: $1 - this should not happen."
       usage ;;
  esac
done

if [ $# -ne 1 ]; then
  usage
fi

log_file=$1


exe=build/cosPT
isapprox=tests/isapprox.jl

sourcedir=$(pwd)
tempdir="$(mktemp -d -t cosPTintTest_XXXXXX)"

mkdir -p "$tempdir"/local
mkdir -p "$tempdir"/input
mkdir -p "$tempdir"/output
mkdir -p "$tempdir"/tests

cp -r -t "$tempdir" "$sourcedir"/src "$sourcedir"/include \
    "$sourcedir"/main.cpp "$sourcedir"/Makefile
cp -r -t "$tempdir/tests/" "$sourcedir"/tests/data \
    "$sourcedir"/tests/ini "$sourcedir"/tests/isapprox.jl
cp -t "$tempdir"/input "$sourcedir"/input/k_*.dat

# Remember git sha from source directory
cd "$sourcedir"
git_sha=$(git rev-parse HEAD)

cd "$tempdir"
# Create src/version_.cpp with git sha from source directory, and add this file
# to the source files listed in the makefile
printf "#include \"../include/version.hpp\"\nstd::string build_git_sha =  \"%s\";" \
    $git_sha > src/version_.cpp
sed -i '/SRC_FILES/ s/\\$/version_.cpp \\/' Makefile
make clean
make -j $n_cores
#export LD_LIBRARY_PATH=path to libs


for f in {L1,L2}; do
    mkdir -p "$tempdir"/output/eds_spt_bs/"$f";
done

{
    echo "Integration test"
    date
} > "$log_file"

for k_a in {000,005,010,015,020,025,030,035,040,045,050,055}; do
    {
        for f in {L1,L2,L2_sh,rsd_L1,rsd_L2,rsd_L2_sh,rsd_ir_L1,rsd_ir_L2,rsd_ir_L2_sh}; do
            {
            mkdir -p "$tempdir"/output/eds_spt_ps/"$f";
            "$exe" --k_a_idx $k_a --n_cores $n_cores "$tempdir"/tests/ini/eds_spt_ps_"$f".cfg
            } >> "$log_file"
        done
    }
done

for k_a in {000,010,020,030,040,050,060,070,080,090,100}; do
    {
        for f in {L1,L2}; do
            {
            mkdir -p "$tempdir"/output/eds_spt_bs/"$f";
            "$exe" --k_a_idx $k_a --k_b_idx $k_a --k_c_idx $k_a --n_cores $n_cores \
                "$tempdir"/tests/ini/eds_spt_bs_"$f".cfg
            } >> "$log_file"
        done
    }
done

for f in output/*/*/; do
    cat "$f"/* > "$f"/total.dat
done

for m in {L1,L2,L2_sh,rsd_L1,rsd_L2,rsd_L2_sh,rsd_ir_L1,rsd_ir_L2,rsd_ir_L2_sh}; do
    {
        julia "$isapprox" --col_A 3 --col_err_A 4 --col_B 3 --col_err_B 4 \
            "$tempdir"/tests/data/eds_spt_ps/"$m".dat "$tempdir"/output/eds_spt_ps/"$m"/total.dat
    } >> "$log_file"
done

for m in {L1,L2}; do
    {
        julia "$isapprox" --col_A 3 --col_err_A 4 --col_B 3 --col_err_B 4 \
            "$tempdir"/tests/data/eds_spt_bs/"$m".dat "$tempdir"/output/eds_spt_bs/"$m"/total.dat
    } >> "$log_file"
done
