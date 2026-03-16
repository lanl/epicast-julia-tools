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
    "attitude_anti-mitigation" => ("total", "anti_mitigation"),
    "attitude_pro-mitigation" => ("total", "pro_mitigation"),
    "status_susceptible" => ("total", "susceptible"),
    "status_presymptomatic" => ("total", "presymptomatic"),
    "status_asymptomatic" => ("total", "asymptomatic"),
    "status_symptomatic" => ("total", "symptomatic"),
    "status_immune" => ("total", "immune"),
    "behavior_withdrawn-fear" => ("total", "withdrawn_spont"),
    "behavior_withdrawn-sick" => ("total", "withdrawn_sick"),
    "behavior_withdrawn-hosp" => ("total", "withdrawn_hosp"),
    #"broadcaster_fear-spreading" => ("media_broadcaster", "broadcaster_spreading"),
    #"broadcaster_fear-countering" => ("media_broadcaster", "broadcaster_countering"),
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
# using OrderedCollections   # preserves order, useful for debugging

"""Parse a TOML file while ignoring any later duplicate keys."""
function clean_toml_keep_first(path::String)
    # We’ll read the file line‑by‑line, build a new string without repeats.
    seen = Set{String}()                 # keys we have already kept
    keep = String[]                      # lines that survive
    for line in eachline(path)
        # Very naive key detection – works for simple “key = value” lines.
        # If you need full TOML syntax (tables, arrays‑of‑tables, etc.) you’ll have
        # to use a proper tokenizer or a more sophisticated regex.
        m = match(r"^\s*([A-Za-z0-9_-]+)\s*=", line)
        if m === nothing
            push!(keep, line)            # comment, blank line, table header, etc.
        else
            key = m.captures[1]
            if key in seen
                @debug "Dropping duplicate key `$key`"
                continue                # skip this line
            else
                push!(keep, line)
                push!(seen, key)
            end
        end
    end
    cleaned = join(keep, "\n")
    return cleaned
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

            rd = nothing
            try
                rd = Epicast.read_runfile(f)
            catch e
                println("  Error reading $f")
                continue
            end
                
            df = to_df(rd)

            f_toml = "$(cur_dir)/$(run)_used_params.toml"
            # EpiCast parameter files currently have some duplicates, we take the first value
            # though they should all be the same
            cleaned = clean_toml_keep_first(f_toml)
            toml = TOML.parse(cleaned)

            #log = JSON.parsefile("$(cur_dir)/$(run)_log.json")
            #m = match(r"social_data_dir=\\\"([\w\/]+)\\\" ",
            #          log["cmd-string"])
            #if m != nothing
            #    toml["social_data_dir"] = m[1]
            #end

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
