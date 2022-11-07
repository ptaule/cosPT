#!/usr/bin/env bash
set -e
set -x

cleanup() {
    rm -rf "$tempdir"
}

# Run cleanup on EXIT
trap cleanup EXIT

# Set some default values:
n_cores=4
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
tempdir="$(mktemp -d -t cosPT_XXXXXXXXXX)"

mkdir -p "$tempdir"/local
mkdir -p "$tempdir"/output
mkdir -p "$tempdir"/ini
mkdir -p "$tempdir"/tests

cp -r -t "$tempdir" "$sourcedir"/src "$sourcedir"/include \
    "$sourcedir"/main.cpp "$sourcedir"/Makefile
cp -r -t "$tempdir"/ini "$sourcedir"/ini/tests/*
cp -r -t "$tempdir"/tests "$sourcedir"/tests/data

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
make -j
#export LD_LIBRARY_PATH=path to libs

mkdir -p "$tempdir"/output/eds_spt_ps/L1
mkdir -p "$tempdir"/output/eds_spt_ps/L2
mkdir -p "$tempdir"/output/eds_spt_ps/L2_sh
mkdir -p "$tempdir"/output/eds_spt_bs/L1
mkdir -p "$tempdir"/output/eds_spt_bs/L2

{
    echo "Integration test"
    date
} > "$log_file"

for k_a in {000,005,010,015,020,025,030,035,040,045,050,055}; do
    {
        "$tempdir"/"$exe" --k_a_idx $k_a --n_cores $n_cores "$tempdir"/ini/eds_spt_ps_L1.cfg
        "$tempdir"/"$exe" --k_a_idx $k_a --n_cores $n_cores "$tempdir"/ini/eds_spt_ps_L2.cfg
        "$tempdir"/"$exe" --k_a_idx $k_a --n_cores $n_cores "$tempdir"/ini/eds_spt_ps_L2_sh.cfg
    } >> "$log_file"
done

for k_a in {000,010,020,030,040,050,060,070,080,090,100}; do
    {
        "$tempdir"/"$exe" --k_a_idx $k_a --k_b_idx $k_a --k_c_idx $k_a --n_cores $n_cores \
            "$tempdir"/ini/eds_spt_bs_L1.cfg
        "$tempdir"/"$exe" --k_a_idx $k_a --k_b_idx $k_a --k_c_idx $k_a --n_cores $n_cores \
            "$tempdir"/ini/eds_spt_bs_L2.cfg
    } >> "$log_file"
done

for f in output/*/*/; do
    cat "$f"/* > "$f"/total.dat
done

{
    julia "$sourcedir/$isapprox" --col_A 3 --col_err_A 4 --col_B 3 --col_err_B 4 \
        "$tempdir"/tests/data/eds_spt_ps/L1/total.dat "$tempdir"/output/eds_spt_ps/L1/total.dat
    julia "$sourcedir/$isapprox" --col_A 3 --col_err_A 4 --col_B 3 --col_err_B 4 \
        "$tempdir"/tests/data/eds_spt_ps/L2/total.dat "$tempdir"/output/eds_spt_ps/L2/total.dat
    julia "$sourcedir/$isapprox" --col_A 3 --col_err_A 4 --col_B 3 --col_err_B 4 \
        "$tempdir"/tests/data/eds_spt_ps/L2_sh/total.dat "$tempdir"/output/eds_spt_ps/L2_sh/total.dat
    julia "$sourcedir/$isapprox" --col_A 5 --col_err_A 6 --col_B 5 --col_err_B 6 \
        "$tempdir"/tests/data/eds_spt_bs/L1/total.dat "$tempdir"/output/eds_spt_bs/L1/total.dat
    julia "$sourcedir/$isapprox" --col_A 5 --col_err_A 6 --col_B 5 --col_err_B 6 \
        "$tempdir"/tests/data/eds_spt_bs/L2/total.dat "$tempdir"/output/eds_spt_bs/L2/total.dat
} >> "$log_file"
