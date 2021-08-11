#!/usr/bin/env bash
set -e
set -x

cleanup() {
    rm -rf $tempdir
    exit $1
}

# Run cleanup on EXIT
trap cleanup EXIT

# Set some default values:
K_A_IDX=-1
K_B_IDX=-1
K_C_IDX=-1
N_EVALS=0
N_CORES=4
INI_FILE=unset
LOG_FILE=unset

usage()
{
  echo "Usage: build_and_run [ --k_a_idx=INDEX ] [ --k_b_idx=INDEX ]
                     [ --k_c_idx=INDEX ] [ --n_evals=NUM ] [ --n_cores=NUM]
                     ini_file log_file"
  exit 2
}

PARSED_ARGUMENTS=$(getopt -n build_and_run -o a:b:c:n:p: --long k_a_idx:,k_b_idx:,k_c_idx:,n_evals:,n_cores: -- "$@")
VALID_ARGUMENTS=$?
if [ "$VALID_ARGUMENTS" != "0" ]; then
  usage
fi

eval set -- "$PARSED_ARGUMENTS"
while :
do
  case "$1" in
    -a | --k_a_idx) K_A_IDX="$2"   ; shift 2 ;;
    -b | --k_b_idx) K_B_IDX="$2"   ; shift 2 ;;
    -c | --k_c_idx) K_C_IDX="$2"   ; shift 2 ;;
    -n | --n_evals) N_EVALS=$2 ; shift 2 ;;
    -p | --n_cores) N_CORES=$2 ; shift 2 ;;
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


INI_FILE="$1"
LOG_FILE="$2"

exe=build/cosPT
GSL="gsl-2.6"
CUBA="Cuba-4.2.1"
LIBCONFIG="libconfig-1.7.2"
LIBPATH="/space/ge52sir/"

sourcedir=$(pwd)
tempdir="$(mktemp -d)"

cp -r $LIBPATH/$GSL $tempdir
cp -r $LIBPATH/$CUBA $tempdir
cp -r $LIBPATH/$LIBCONFIG $tempdir
cp -r -t $tempdir $sourcedir/src $sourcedir/include \
    $sourcedir/main.cpp $sourcedir/Makefile

mkdir -p $tempdir/local/

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

$tempdir/$exe --k_a_idx $K_A_IDX --k_b_idx $K_B_IDX --k_c_idx $K_C_IDX --n_evals $N_EVALS --n_cores $N_CORES $INI_FILE > $LOG_FILE
