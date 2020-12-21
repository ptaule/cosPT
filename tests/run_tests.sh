#!/usr/bin/env bash
set -e
set -x

cleanup() {
    rm -rf $tempdir
}

# Run cleanup on EXIT
trap cleanup EXIT

# Set some default values:
N_CORES=4
LOG_FILE=unset

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
    -p | --n_cores) N_CORES=$2 ; shift 2 ;;
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

LOG_FILE=$1

julia=/space/ge52sir/bin/julia
isapprox=/home/t30/all/ge52sir/non_linear_PS/tests/isapprox.jl

exe=build/cosPT
GSL="gsl-2.6"
CUBA="Cuba-4.2.1"
LIBCONFIG="libconfig-1.7.2"
LIBPATH="/space/ge52sir/"

sourcedir=$(pwd)
tempdir="$(mktemp -d)"

mkdir -p $tempdir/local
mkdir -p $tempdir/output
mkdir -p $tempdir/ini
mkdir -p $tempdir/tests

cp -r $LIBPATH/$GSL $tempdir
cp -r $LIBPATH/$CUBA $tempdir
cp -r $LIBPATH/$LIBCONFIG $tempdir
cp -r -t $tempdir $sourcedir/src $sourcedir/include \
    $sourcedir/main.cpp $sourcedir/Makefile
cp -r -t $tempdir/ini $sourcedir/ini/tests/*
cp -r -t $tempdir/tests $sourcedir/tests/data

cd $tempdir/$LIBCONFIG
autoreconf -f -i
./configure --prefix=$tempdir/local/
make clean
make -j && make -j install

cd $tempdir/$CUBA
./configure --prefix=$tempdir/local/
make clean
make -j install

cd $tempdir/$GSL
./configure --prefix=$tempdir/local/
make clean
make -j && make -j install

cd $sourcedir
git_sha=$(git rev-parse HEAD)
cd $tempdir
printf "#include \"../include/version.hpp\"\\nstd::string build_git_sha =  \"${git_sha}\";" \
    > src/version_.cpp

cd $tempdir
mkdir obj
mkdir build
make clean
make -j cluster
export LD_LIBRARY_PATH=$tempdir/local/lib/

mkdir -p $tempdir/output/eds_spt_ps/L1
mkdir -p $tempdir/output/eds_spt_ps/L2
mkdir -p $tempdir/output/eds_spt_bs/L1
mkdir -p $tempdir/output/eds_spt_bs/L2

echo "Integration test" > $LOG_FILE
echo "$(date)" >> $LOG_FILE

for K_A in {000,005,010,015,020,025,030,035,040,045,050,055}; do
    $tempdir/$exe --k_a $K_A --n_cores $N_CORES $tempdir/ini/eds_spt_ps_L1.cfg \
        >> $LOG_FILE
    $tempdir/$exe --k_a $K_A --n_cores $N_CORES $tempdir/ini/eds_spt_ps_L2.cfg \
        >> $LOG_FILE
done

for K_A in {000,010,020,030,040,050,060,070,080,090,100}; do
    $tempdir/$exe --k_a $K_A --k_b $K_A --k_c $K_A --n_cores $N_CORES \
        $tempdir/ini/eds_spt_bs_L1.cfg >> $LOG_FILE
    $tempdir/$exe --k_a $K_A --k_b $K_A --k_c $K_A --n_cores $N_CORES \
        $tempdir/ini/eds_spt_bs_L2.cfg >> $LOG_FILE
done

for f in output/*/*/; do
    cat $f/* > $f/total.dat
done

$julia $isapprox --col_A 3 --col_err_A 4 --col_B 3 --col_err_B 4 \
    $tempdir/tests/data/eds_spt_ps/L1/total.dat $tempdir/output/eds_spt_ps/L1/total.dat \
    >> $LOG_FILE
$julia $isapprox --col_A 3 --col_err_A 4 --col_B 3 --col_err_B 4 \
    $tempdir/tests/data/eds_spt_ps/L2/total.dat $tempdir/output/eds_spt_ps/L2/total.dat \
    >> $LOG_FILE
$julia $isapprox --col_A 5 --col_err_A 6 --col_B 5 --col_err_B 6 \
    $tempdir/tests/data/eds_spt_bs/L1/total.dat $tempdir/output/eds_spt_bs/L1/total.dat \
    >> $LOG_FILE
$julia $isapprox --col_A 5 --col_err_A 6 --col_B 5 --col_err_B 6 \
    $tempdir/tests/data/eds_spt_bs/L2/total.dat $tempdir/output/eds_spt_bs/L2/total.dat \
    >> $LOG_FILE
