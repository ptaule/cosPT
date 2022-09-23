#!/usr/bin/env julia

#   cont_sub.jl
#
#   Created by Petter Taule on 17.11.2020
#   Copyright (c) 2020 Petter Taule. All rights reserved.


using ArgParse, Printf

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--bispectrum"
            help = "Bispectrum: set k_b_idx = k_c_idx = k_a_idx"
            action = :store_true
        "--longrun"
            help = "Use queue = longrun"
            action = :store_true
        "--job_prefix"
            help = "SGE job prefix"
            arg_type = String
            default = "ps"
        "--n_evals_step"
            help = "n_evals step"
            arg_type = Float64
            default = 1e4
        "--n_evals_max"
            help = "max. n_evals"
            arg_type = Float64
            default = 1e5
        "--n_cores"
            help = "number of cores to ask for"
            arg_type = Int
            default = 4
        "--log_dir"
            help = "log directory"
            arg_type = String
            default = "/space/ge52sir/sge_output/"
        "--script"
            help = "script to run"
            arg_type = String
            default = "run.sh"
        "ini_file"
            help = "ini file"
            arg_type = String
            required = true
        "k_idx_lower"
            help = "k index lower"
            arg_type = Int
            required = true
        "k_idx_upper"
            help = "k index upper (set to k_idx_lower if not given)"
            arg_type = Int
            required = false
    end

    return parse_args(s)
end

function continuous_submit(k_a_idx::Int,
                           k_b_idx::Int,
                           k_c_idx::Int,
                           n_evals_arr::Vector,
                           n_cores::Int,
                           job_prefix::String,
                           script::String,
                           log_dir::String,
                           ini_file::String,
                           longrun::Bool
                          )
    for n_evals in n_evals_arr
        n_evals_str = @sprintf "%.1e" n_evals
        job_name = job_prefix * "_" * string(k_a_idx) * "_" * n_evals_str
        log_file = log_dir * job_name * ".log"

        # Wait for the last job to finish
        while (occursin(job_prefix * "_$(k_a_idx)_", read(`qstat_formatted.sh`, String)))
            sleep(30)
        end

        qsub_cmd = ``
        if longrun
            qsub_cmd = `qsub -q longrun -N $job_name -pe smp $n_cores
            -e $log_dir/error/ -o $log_dir/output/ $script --k_a_idx=$k_a_idx
            --k_b_idx=$k_b_idx --k_c_idx=$k_c_idx --n_evals $n_evals $ini_file $log_file`
        else
            qsub_cmd = `qsub -N $job_name -pe smp $n_cores -e $log_dir/error/
            -o $log_dir/output/ $script --k_a_idx=$k_a_idx --k_b_idx=$k_b_idx
            --k_c_idx=$k_c_idx --n_evals $n_evals $ini_file $log_file`
        end
        run(pipeline(qsub_cmd, devnull))
    end
end



function main()
    parse_args = parse_commandline()

    job_prefix = parse_args["job_prefix"]
    step       = parse_args["n_evals_step"]
    max        = parse_args["n_evals_max"]

    n_cores = parse_args["n_cores"]
    log_dir = parse_args["log_dir"]
    script  = parse_args["script"]

    k_idx_lower = parse_args["k_idx_lower"]
    k_idx_upper = parse_args["k_idx_upper"]

    if (k_idx_upper == nothing)
        k_idx_upper = k_idx_lower
    end

    k_b_idx = -1
    k_c_idx = -1

    ini_file = pwd() * "/" * parse_args["ini_file"]

    n_evals_arr = [i*step for i=1 : (max/step)]

    # Wait for for-loop to finish before continuing (including forked processes)
    @sync for k_a_idx = k_idx_lower : k_idx_upper
        if (parse_args["bispectrum"])
            # Spawn a new process
            @async continuous_submit(k_a_idx, k_a_idx, k_a_idx, n_evals_arr,
                                     n_cores, job_prefix, script, log_dir,
                                     ini_file, parse_args["longrun"])
        else
            # Spawn a new process
            @async continuous_submit(k_a_idx, k_b_idx, k_c_idx, n_evals_arr,
                                     n_cores, job_prefix, script, log_dir,
                                     ini_file, parse_args["longrun"])
        end
    end
end


main()
