#!/usr/bin/env bash

julia_depot_path=/space/ge52sir/.julia
log_dir="/space/ge52sir/sge_output"

qsub -N $1 -pe smp 4 -cwd -e "${log_dir}"/error/ -o "${log_dir}"/output/ \
    -v JULIA_DEPOT_PATH=$julia_depot_path tests/run_tests.sh $2
