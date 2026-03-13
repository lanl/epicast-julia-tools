#!/usr/bin/env julia

import Pkg

Pkg.activate("Epicast")
using ArgParse
using Epicast
using DataFrames
using CSV
using Base
using DelimitedFiles
using Glob
using TOML
using JSON
# ============================================================================ #
function parse_args(args)
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--in-dir", "-i"
            default="/home/jkitson/work/epicast/results/scenarios/local"
            help="Path to directory containing the input bin files"
        "--out-dir", "-o"
            default=""
            help="Path to directory to save the resulting CSV to"
    end

    return ArgParse.parse_args(args, s)
end
# ============================================================================ #
function read_bins(dir, id)
    if isdir("$dir/$id")
        tmp = nothing
        foreach(readdir("$dir/$id")) do f
            if endswith(f, ".bin")
                tmp = "$dir/$id/$f"
            end
        end
        return tmp
    end
end
# ============================================================================ #
fear_cols = Dict(
    "total" => ("total", "case_counts"),
    "mentalstate_fear-disease" => ("total", "pro_mitigation"),
    "mentalstate_fear-mitigation" => ("total", "anti_mitigation"),
    "status_prodomal" => ("total", "symptomatic"),
    "status_pox" => ("total", "asymptomatic"),
    "status_latent" => ("total", "presymptomatic"),
    "status_immune" => ("total", "immune"),
    "mentalstate_fear-mitigation" => ("total", "anti_mitigation"),
    "behavior_withdrawn-spont" => ("total", "withdrawn_spont"),
    "behavior_withdrawn-sick" => ("total", "withdrawn_sick"),
    "behavior_withdrawn-hosp" => ("total", "withdrawn_hosp"),
    "broadcaster_fear-spreading" => ("media_broadcaster", "broadcaster_spreading"),
    "broadcaster_fear-countering" => ("media_broadcaster", "broadcaster_countering"),
)
total_and_new = Dict(
    "" => Epicast.total_cases,
    "new_" => Epicast.new_cases
)
# ---------------------------------------------------------------------------- #
function to_df(data::Epicast.RunData;
        col_to_demo=fear_cols, reducers=total_and_new)
    table_like = Dict()
    for (col, (demo, short_col)) in col_to_demo
        for (n, freduce) in reducers
            table_like["$(n)$short_col"] = Epicast.aggregate(data, col; freduce=freduce, demo=demo)
        end
    end
    return DataFrames.DataFrame(table_like)
end
# ============================================================================ #
function main(args)
    in_dir = args["in-dir"]
    bin_files = glob("$in_dir/*/*.bin")
    for f in bin_files
        if f != nothing
            println("Reading $f")
            cur_dir = dirname(f)
            (run, ext) = splitext(basename(f))
            rd = Epicast.read_runfile(f)
            df = to_df(rd)

            f_toml = "$(cur_dir)/$(run)_used_params.toml"
            toml = TOML.parsefile(f_toml)
            log = JSON.Parser.parsefile("$(cur_dir)/$(run)_log.json")
            m = match(r"social_data_dir=\\\"([\w\/]+)\\\" ",
                      log["cmd-string"])
            if m != nothing
                toml["social_data_dir"] = m[1]
            end

            out_dir = args["out-dir"]
            if out_dir == ""
                out_dir = dirname(f)
            end

            f_csv = joinpath(out_dir,
                    basename(replace(f, r".bin" => ".csv")))

            println("  Saving CSV version to $f_csv")
            CSV.write(f_csv, df)

            #m = match(r"run_([0-9]+)", f_toml)
            #if m != nothing
            #    f_toml = "run_$(m[1])_used_params.toml"
            #end
            open(joinpath(out_dir, basename(f_toml)), "w") do io
                TOML.print(io, toml)
            end
        end
    end
end

main(parse_args(ARGS))
