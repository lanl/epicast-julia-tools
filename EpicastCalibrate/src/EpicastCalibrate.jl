module EpicastCalibrate

# Copyright (C) 2025. Triad National Security, LLC. All rights reserved.
# This program was produced under U.S. Government contract 89233218CNA000001
# for Los Alamos National Laboratory (LANL), which is operated by Triad National
# Security, LLC for the U.S. Department of Energy/National Nuclear Security
# Administration. All rights in the program are reserved by Triad National
# Security, LLC, and the U.S. Department of Energy/National Nuclear Security
# Administration. The Government is granted for itself and others acting on its
# behalf a nonexclusive, paid-up, irrevocable worldwide license in this material
# to reproduce, prepare. derivative works, distribute copies to the public,
# perform publicly and display publicly, and to permit others to do so.

using Epicast, TOML, NPZ
using Epicast.EpicastTables
using DelimitedFiles, PyPlot, Statistics, Sobol

const UNDER_REPORTING_FACTOR = 3
# ============================================================================ #
function write_policy(io::IO, scale::Vector{<:AbstractFloat}, counties::Vector{<:Integer})
    N = length(counties)
    for ctx in ["work","social"]
        println(io, ctx, " = [")
        for k in eachindex(counties)
            print(io, "{scale = $(scale[k]), fips = [$(counties[k])]}")
            if k < N
                println(io, ",")
            else
                println(io, "")
            end
        end
        println(io, "]")
    end
end
# ============================================================================ #
@inline function rand_bt(n::Integer, min_v::AbstractFloat, max_v::AbstractFloat)
    return rand(n) .* (max_v - min_v) .+ min_v
end
# ============================================================================ #
function fixed_index_cases(x::Integer=1)
    return [111,5,84,5,5,24,5,81,63,8,5,5,5,56,5,5,23,28,5,6,5,8,11,13,22,5,36,
        5,5,5,5,5,10] .* x
end
# ============================================================================ #
function generate_case_file(param_file::AbstractString,
    counties::Vector{<:Integer}, pop::Vector{<:Integer}, im::Vector{Float64},
    run_n::Integer, odir::AbstractString=dirname(ifile))

    # in ascending sorted FIPS order for NM counties
    n_idx = fixed_index_cases(UNDER_REPORTING_FACTOR)

    ofile = joinpath(odir, 
        replace(
            basename(param_file),
            ".toml" => "_run_" * lpad(run_n, 3, '0') * ".cases"
        )
    )

    open(ofile, "w") do io
        for (k, county) in enumerate(counties)
            println(io, county, " ", n_idx[k], " ",
                round(Int, im[k] * pop[k])
            )
        end
    end
end
# ============================================================================ #
function generate_toml_file(ifile::AbstractString, counties::Vector{<:Integer},
    par::Vector{Float64}, run_n::Integer, odir::AbstractString=dirname(ifile))
    
    param_str = read(ifile, String)

    tmp = replace(param_str, r"run_number = \d+" => "run_number = $(run_n)")

    tmp = replace(tmp, r"p_trans = \d?\.?\d+" => "p_trans = $(par[1])")
    tmp = replace(tmp, r"p_asymptomatic = \d?\.?\d+" => "p_asymptomatic = $(par[2])")
    tmp = replace(tmp, r"rel_trans_asymptomatic = \d?\.?\d+" => "rel_trans_asymptomatic = $(par[3])")
    tmp = replace(tmp, r"withdrawal_scalar = \d?\.?\d+" => "withdrawal_scalar = $(par[4])")

    ofile = joinpath(odir, 
        replace(basename(ifile), ".toml" => "_run_" * lpad(run_n, 3, '0') * ".toml")
    )

    idx_case_file = joinpath(
        "/vast/home/palexander/sandbox/nm_multi-param_sweep-all_params",
        # "/Users/palexander/Documents/emerge+radium/testing_results/nm_multi-param_sweep-all/test",
        replace(basename(ofile), ".toml" => ".cases")
    )

    tmp = replace(tmp, r"index_case_file = [^\n]+" => "index_case_file = \"$(idx_case_file)\"")

    open(ofile, "w") do io
        print(io, tmp)
        scale = par[5:(5 + 33 - 1)]
        write_policy(io, scale, counties)
    end
end
# ============================================================================ #
function generate_param_files(param_file::AbstractString,
    county_file::AbstractString, n::Integer, odir::AbstractString)

    cnty_data = readdlm(county_file, ' ', Int)

    seq = SobolSeq(
        vcat([0.075, 0.05, 0.05, 0.0], fill(0.15, 33), fill(0.00, 33)),
        vcat([0.350, 0.95, 0.95, 2.0], fill(1.00, 33), fill(0.50, 33))
    )

    skip(seq, n)

    cache = zeros(70)
    # out = zeros(70,n)
    for k in 1:n

        next!(seq, cache)
        cache .+= (0.0001 .* randn(70))

        cache .= max.(cache, 0.0)

        generate_toml_file(param_file, cnty_data[:,1], cache[1:(33 + 4)],
            (k-1) * 5, odir)

        generate_case_file(param_file, cnty_data[:,1], cnty_data[:,2],
            cache[(33 + 5):end], (k-1) * 5, odir)

    end

    # return out
end
# ============================================================================ #
function read_case_file_immunity(ifile::AbstractString)
    d = readdlm(ifile, ' ', Int)
    ks = sortperm(d[:,1])
    return d[ks,3]
end
# ============================================================================ #
get_run_number(ifile::AbstractString) = parse(Int, match(r".*run_(\d+).*", ifile)[1])
# ============================================================================ #
function write_npz2(data_dirs::Vector{<:AbstractString},
    param_dirs::Vector{<:AbstractString}, counties::Vector{<:Integer},
    pop::Vector{<:Integer}, ofile::AbstractString)

    n_geo = 33
    n_global = 4 # [p_trans, p_asymp, rel_tans_asymp, withdrawl]
    n_local = 2  # [work/social_scale, initial_immune]
    n_rep = 5

    cols = vcat(
        ["p_trans","p_asymptomatic","rel_trans_asymptomatic","withdrawal_scalar"],
        map(x -> "iisf_" * string(x), counties),
        map(x -> "immune_" * string(x), counties),
    )

    param_files = map(x -> Epicast.find_files(x, r".*\.toml"), param_dirs)
    run_n = collect(Iterators.flatten(map(x -> map(get_run_number, x), param_files)))
        # map(x -> get_run_number(x) + 1000, param_files[2]))

    ks = sortperm(run_n)
    run_n = run_n[ks]
    param_files = reduce(vcat, param_files)
    param_files .= param_files[ks]

    params = zeros(Float64, length(param_files) * n_rep, n_global + n_local * n_geo)
    out = zeros(Float64, length(param_files) * n_rep, 250, n_geo)

    inc = 1

    for (k,file) in enumerate(param_files)

        p = TOML.parsefile(file)

        p_trans = p["p_trans"]
        p_asymptomatic = p["p_asymptomatic"]
        rel_trans_asymptomatic = p["rel_trans_asymptomatic"]
        withdrawal_scalar = p["withdrawal_scalar"]

        policies = sort(p["policies"]["0"]["work"], lt=(a,b)->a["fips"][1] < b["fips"][1])
        work_scale = get.(policies, "scale", NaN)

        case_file = replace(param_files[k], ".toml" => ".cases")
        
        @assert(isfile(case_file), "$(case_file)")
        im = read_case_file_immunity(case_file)

        run = run_n[k]

        for j in run:(run+4)
            data_dir = data_dirs[1]
            data_file = joinpath(data_dir, "run_" * lpad(j, 3, '0') * ".bin")

            data = Epicast.preprocess!(
                Epicast.aggregate(Epicast.County, Epicast.read_runfile(data_file)),
                "total",
                smooth=true,
                diff=true,
                get_denom=Epicast.noop_denom
            )

            out[inc,:,:] .= data.data["total"]

            params[inc,1] = p_trans
            params[inc,2] = p_asymptomatic
            params[inc,3] = rel_trans_asymptomatic
            params[inc,4] = withdrawal_scalar

            idx = n_global + n_geo
            params[inc,(n_global + 1):idx] .= work_scale
            params[inc,idx+1:end] .= im ./ pop

            inc += 1
        end
    end

    npzwrite(ofile, out)
    npzwrite(replace(ofile, ".npz" => "_params.npz"), params)

    open(replace(ofile, ".npz" => "_params.names"), "w") do io
        for col in cols
            println(io, col)
        end
    end
end
# ============================================================================ #
function write_npz(::Type{G}, idir::AbstractString, ofile::AbstractString,
    avg::Bool=true) where {G<:AbstractGeo}

    files = Epicast.find_files(idir, r".*\.bin")

    data = Epicast.preprocess!(
        Epicast.aggregate(G, Epicast.read_runfile(files[1])),
        "total",
        smooth=true,
        diff=true,
        get_denom=noop_denom
    )

    len, n_loc = size(data.data["total"])

    p_trans = zeros(Float64, length(files))

    for k in 1:length(files)
        toml_file = replace(files[k], ".bin" => "_used_params.toml")

        p = TOML.parsefile(toml_file)

        p_trans[k] = p["p_trans"]
    end

    if avg
        up = sort!(unique(p_trans))
    else
        up = p_trans
    end
    
    out = zeros(Float64, length(up), len, n_loc)

    @show(length(p_trans), length(up))

    for k in 1:length(up)
        idx = avg ? findall(isequal(up[k]), p_trans) : [k]
        for j in idx
            tmp = Epicast.preprocess!(
                Epicast.aggregate(G, Epicast.read_runfile(files[j])),
                "total",
                smooth=true,
                diff=true,
                get_denom=noop_denom
            )
            out[k,:,:] .+= tmp.data["total"]
        end

        avg && (out[k,:,:] ./= length(idx))
    end

    npzwrite(ofile, out)

    npzwrite(replace(ofile, ".npz" => "_params.npz"), up)

    return out, up
end
# ============================================================================ #
function generate_files_from_params(param_file::AbstractString,
    counties::Vector{<:Integer}, pop::Vector{<:Integer}, par::Vector{<:Real},
    ofile::AbstractString)

    param_str = read(param_file, String)

    tmp = replace(param_str, r"p_trans = \d?\.?\d+" => "p_trans = $(par[1])")
    tmp = replace(tmp, r"p_asymptomatic = \d?\.?\d+" => "p_asymptomatic = $(par[2])")
    tmp = replace(tmp, r"rel_trans_asymptomatic = \d?\.?\d+" => "rel_trans_asymptomatic = $(par[3])")
    tmp = replace(tmp, r"withdrawal_scalar = \d?\.?\d+" => "withdrawal_scalar = $(par[4])")

    idx_case_file = replace(ofile, ".toml" => ".cases")
    tmp = replace(tmp, r"index_case_file = [^\n]+" => "index_case_file = \"$(idx_case_file)\"")

    last = (5+length(counties)-1)

    open(ofile, "w") do io
        print(io, tmp)
        write_policy(io, par[5:last], counties)
    end

    n_idx = fixed_index_cases(3)

    open(idx_case_file, "w") do io
        for k in eachindex(counties)
            println(io, counties[k], ' ', n_idx[k], ' ', round(Int, par[last+k] * pop[k]))
        end
    end

end
# ============================================================================ #
function load_all_data(idir::AbstractString)
    files = Epicast.find_files(idir, r".*\.bin$")
    out = zeros(Float64, 250, 33, length(files))
    fips = Int[]
    for (k,file) in enumerate(files)
        data = Epicast.preprocess!(
            Epicast.aggregate(County, Epicast.read_rundata(file)),
            "total",
            diff=true,
            smooth=true,
            get_denom=Epicast.noop_denom
        )
        out[:,:,k] .= data.data["total"]
        if isempty(fips)
            fips = sort!(collect(keys(data.data.fips_index)))
        end
    end

    # n = div(length(files), 5)
    # tmp = zeros(250, 33, n)
    # inc = 1
    # for k in 1:5:length(files)
    #     tmp[:,:,inc] .= mean(out[:,:,k:k+4], dims=3)
    #     inc += 1
    # end

    return out, fips
end
# ============================================================================ #
function comparison_plot(idir::AbstractString, obs_file::AbstractString,
    single_plot::Bool=true)

    obs = npzread(obs_file)

    data, fips = load_all_data(idir)

    # @show(sum(data.data["total"]), sum(obs))

    if single_plot
        h, ax = subplots(1,1)

        cols = plt.rcParams["axes.prop_cycle"].by_key()["color"]
        ax.set_prop_cycle("color", cols)

        # ax.plot(data.data["total"])# ./ 10)

        for k in 1:33
            mn = dropdims(mean(data[:,k,:], dims=2), dims=2)
            sd = dropdims(std(data[:,k,:], dims=2), dims=2)
            hp = ax.fill_between(0:249, mn .- sd, mn .+ sd, alpha=0.2)
            ax.plot(0:249, mn, color=hp.get_facecolor()[1:3])
        end
        ax.set_prop_cycle("color", cols)

        # ax.plot(obs, "--")

        ax.spines["right"].set_visible(false)
        ax.spines["top"].set_visible(false)

        ax.set_xlabel("Simulation day (0 = 2020-09-13)", fontsize=14)
        ax.set_ylabel("Newly exposed per day", fontsize=14)

        h.tight_layout()

    else

        h, ax = subplots(7,5)
        h.set_size_inches((12,12))

        for k in 1:33
            # mn = dropdims(mean(data[:,k,:], dims=2), dims=2)
            # sd = dropdims(std(data[:,k,:], dims=2), dims=2)
            # hp = ax[k].fill_between(0:249, mn .- sd, mn .+ sd, alpha=0.2)
            # ax[k].plot(0:249, mn, color=hp.get_facecolor()[1:3])

            ax[k].plot(0:249, data[:,k,:], color="C0")

            ax[k].plot(obs[:,k], color="C1")

            ax[k].set_title("County $(fips[k])", fontsize=14)
            ax[k].spines["right"].set_visible(false)
            ax[k].spines["top"].set_visible(false)
        end

        ax[end-1].set_visible(false)
        ax[end].set_visible(false)
        h.tight_layout()

        n = length(ax[27].lines)
        hl = map(k -> ax[27].lines[k], [1,n])
        ax[27].legend(hl, ["Epicast", "Observed"], fontsize=14, frameon=false,
            loc="upper left", bbox_to_anchor=(1.0,1.0))

    end

    return h, ax
end
# ============================================================================ #
params_from_string(str::AbstractString) = parse.(Float64, split(str, r"\s+"))
# ============================================================================ #
end # module EpicastCalibrate
