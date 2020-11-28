#!/usr/bin/env bash

set -e
set -x

# Set some default values:
K_A=-1
K_B=-1
K_C=-1
N_EVALS=0
N_CORES=4
INI_FILE=unset
LOG_FILE=unset

usage()
{
  echo "Usage: build_and_run [ --k_a INDEX ] [ --k_b INDEX ]
                     [ --k_c INDEX ] [ --n_evals NUM ] [ --n_cores NUM]
                     ini_file log_file"
  exit 2
}

clean() {
    rm -rf $tempdir
    exit $1
}

PARSED_ARGUMENTS=$(getopt -n build_and_run -o a:b:c:n:p: --long k_a:,k_b:,k_c:,n_evals:,n_cores: -- "$@")
VALID_ARGUMENTS=$?
if [ "$VALID_ARGUMENTS" != "0" ]; then
  usage
fi

eval set -- "$PARSED_ARGUMENTS"
while :
do
  case "$1" in
    -a | --k_a)     K_A="$2"   ; shift 2 ;;
    -b | --k_b)     K_B="$2"   ; shift 2 ;;
    -c | --k_c)     K_C="$2"   ; shift 2 ;;
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
CUBA="Cuba-4.2"
LIBCONFIG="libconfig-1.7.2"
LIBPATH="/space/ge52sir/"

sourcedir=$(pwd)
tempdir="$(mktemp -d)"

cp -r $LIBPATH/$GSL $tempdir || clean 1
cp -r $LIBPATH/$CUBA $tempdir || clean 1
cp -r $LIBPATH/$LIBCONFIG $tempdir || clean 1
cp -r -t $tempdir $sourcedir/src $sourcedir/include \
    $sourcedir/main.cpp $sourcedir/Makefile || clean 1

mkdir -p $tempdir/local/

cd $tempdir/$LIBCONFIG || clean 1
autoreconf -f -i || clean 1
./configure --prefix=$tempdir/local/ || clean 1
make clean || clean 1
make -j || clean 1
make -j install || clean 1

cd $tempdir/$CUBA || clean 1
./configure --prefix=$tempdir/local/ || clean 1
make clean || clean 1
make -j || clean 1
make -j install || clean 1

cd $tempdir/$GSL || clean 1
./configure --prefix=$tempdir/local/ || clean 1
make clean || clean 1
make -j || clean 1
make -j install || clean 1

cd $sourcedir || clean 1
git_sha=$(git rev-parse HEAD) || clean 1
cd $tempdir
printf "#include \"../include/version.hpp\"\\nstd::string build_git_sha =  \"${git_sha}\";" \
    > src/version_.cpp || clean 1

cd $tempdir || clean 1
mkdir obj || clean 1
mkdir build || clean 1
make clean || clean 1
make -j cluster || clean 1
export LD_LIBRARY_PATH=$tempdir/local/lib/

$tempdir/$exe --k_a $K_A --k_b $K_B --k_c $K_C --n_evals $N_EVALS --n_cores $N_CORES $INI_FILE > $LOG_FILE || clean 1

clean 0