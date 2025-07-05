#!/usr/bin/env bash

set -e

cleanup() {
    git worktree remove -f "$tempdir"
}

# Run cleanup on EXIT
trap cleanup EXIT

# Set some default values:
n_cores=1
log_file=
use_python=false

usage() {
    echo "Usage: run_tests.sh [ --n_cores NUM] LOG_FILE"
    exit 2
}

PARSED_ARGUMENTS=$(getopt -n run_tests -o hn:p --long help,n_cores:,python -- "$@")
VALID_ARGUMENTS=$?
if [ "$VALID_ARGUMENTS" != "0" ]; then
    usage
fi

eval set -- "$PARSED_ARGUMENTS"
while :; do
    case "$1" in
    -h | --help)
        usage
        shift 1
        ;;
    -n | --n_cores)
        n_cores=$2
        shift 2
        ;;
    -p | --python)
        use_python=true
        shift 1
        ;;
    # -- means the end of the arguments; drop this, and break out of the while loop
    --)
        shift
        break
        ;;
    # If invalid options were passed, then getopt should have reported an error,
    # which we checked as VALID_ARGUMENTS when getopt was called...
    *)
        echo "Unexpected option: $1 - this should not happen."
        usage
        ;;
    esac
done

if [ $# -ne 1 ]; then
    usage
fi

log_file=$1

exe=build/cosPT
compare_prog=julia
isapprox=tests/isapprox.jl
if [[ "$use_python" == "true" ]]; then
    compare_prog=python
    isapprox=tests/isapprox.py
fi

sourcedir=$(pwd)
tempdir="$(mktemp -d -t cosPTintTest_XXXXXX)"

git worktree add "$tempdir"

# Make binary
cd "$tempdir"
#export LD_LIBRARY_PATH=path to libs
make clean
make -j "$n_cores"

{
    echo "Integration test"
    date
} >"$log_file"

for k_a in {000,005,010,015,020,025,030,035,040,045,050,055}; do
    {
        for f in {L1,L2,L2_sh,rsd_L1,rsd_L2,rsd_L2_sh,rsd_ir_L1,rsd_ir_L2,rsd_ir_L2_sh,bias_ir_L1}; do
            {
                "$exe" --k_a_idx "$k_a" --n_cores "$n_cores" "$tempdir"/tests/ini/eds_spt_ps_"$f".cfg
            } >>"$log_file"
        done
    }
done

for k_a in {000,010,020,030,040,050,060,070,080,090,100}; do
    {
        for f in {L1,L2}; do
            {
                "$exe" --k_a_idx "$k_a" --k_b_idx "$k_a" --k_c_idx "$k_a" --n_cores $n_cores \
                    "$tempdir"/tests/ini/eds_spt_bs_"$f".cfg
            } >>"$log_file"
        done
    }
done

for k_a in {000,010,020,030,040,050,060,070,080}; do
    {
        "$exe" --k_a_idx "$k_a" --n_cores "$n_cores" "$tempdir"/tests/ini/quijote_Mnu_0p1eV_L1.cfg
    } >>"$log_file"
done

for k_a in {030,045}; do
    {
        "$exe" --k_a_idx "$k_a" --n_cores "$n_cores" "$tempdir"/tests/ini/quijote_Mnu_0p1eV_L2.cfg
    } >>"$log_file"
done

for f in output/*/*/; do
    cat "$f"/* >"$f"/total.dat
done

# Comparison
for m in {L1,L2,L2_sh,rsd_L1,rsd_L2,rsd_L2_sh,rsd_ir_L1,rsd_ir_L2,rsd_ir_L2_sh,bias_ir_L1}; do
    {
        "$compare_prog" "$isapprox" --col_A 3 --col_err_A 4 --col_B 3 --col_err_B 4 \
            "$tempdir"/tests/data/eds_spt_ps/"$m".dat "$tempdir"/output/eds_spt_ps/"$m"/total.dat
    } >>"$log_file"
done

for m in {L1,L2}; do
    {
        "$compare_prog" "$isapprox" --col_A 3 --col_err_A 4 --col_B 3 --col_err_B 4 \
            "$tempdir"/tests/data/eds_spt_bs/"$m".dat "$tempdir"/output/eds_spt_bs/"$m"/total.dat
    } >>"$log_file"
done

for m in {L1,L2}; do
    {
        "$compare_prog" "$isapprox" --col_A 3 --col_err_A 4 --col_B 3 --col_err_B 4 \
            "$tempdir"/tests/data/quijote_Mnu_0p1eV/"$m".dat "$tempdir"/output/quijote_Mnu_0p1eV/"$m"/total.dat
        "$compare_prog" "$isapprox" --col_A 5 --col_err_A 6 --col_B 5 --col_err_B 6 \
            "$tempdir"/tests/data/quijote_Mnu_0p1eV/"$m".dat "$tempdir"/output/quijote_Mnu_0p1eV/"$m"/total.dat
        "$compare_prog" "$isapprox" --col_A 7 --col_err_A 8 --col_B 7 --col_err_B 8 \
            "$tempdir"/tests/data/quijote_Mnu_0p1eV/"$m".dat "$tempdir"/output/quijote_Mnu_0p1eV/"$m"/total.dat
    } >>"$log_file"
done
