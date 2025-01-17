using DelimitedFiles, ArgParse, Measurements


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--col_A"
            help = "column in file A to compare"
            arg_type = Int
            default = 1
        "--col_err_A"
            help = "corresponding error in file A"
            arg_type = Int
        "--col_B"
            help = "column in file B to compare"
            arg_type = Int
            default = 1
        "--col_err_B"
            help = "corresponding error in file B"
            arg_type = Int
        "--verbose", "-v"
            help = "be verbose"
            action = :store_true
        "fileA"
            help = "File A"
            arg_type = String
            required = true
        "fileB"
            help = "File B"
            arg_type = String
            required = true
    end

    return parse_args(s)
end

function main()
    parse_args = parse_commandline()

    file_a = parse_args["fileA"]
    file_b = parse_args["fileB"]

    a_data = readdlm(file_a, comments=true)
    b_data = readdlm(file_b, comments=true)

    col_a = parse_args["col_A"]
    col_err_a = parse_args["col_err_A"]
    col_b = parse_args["col_B"]
    col_err_b = parse_args["col_err_B"]

    a = measurement.(a_data[:,col_a], (isnothing(col_err_a) ? 0 : a_data[:,col_err_a]))
    b = measurement.(b_data[:,col_b], (isnothing(col_err_b) ? 0 : b_data[:,col_err_b]))

    N = size(a)[1]

    if N != size(b)[1]
        @error "Column dimensions in file A and B differ"
    end

    println("Comparing column $col_a from $file_a and column $col_b from $file_b:")
    if a == b
        println("Colums are equal")
        return
    end

    zscore = stdscore.(a,b)

    for i = 1:N
        z = zscore[i]
        if (parse_args["verbose"] || abs(z) > 1)
            println("Row $i: a = $(a[i]), b = $(b[i]), zscore = $z")
        end
    end
    println("done.")
end

main()
