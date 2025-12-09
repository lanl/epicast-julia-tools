module EpicastCalibrate

using Epicast, TOML, NPZ
using Epicast.EpicastTables
using DelimitedFiles
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
function rs_rand(n::Integer, min_v::AbstractFloat, max_v::AbstractFloat, epsilon::AbstractFloat=1e-4)

    out = Vector{Float64}(undef, n)
    out[1] = (rand() * (max_v - min_v)) + min_v

    for k = 2:n

        mn = -Inf
        r = 0.0
        attempt = 0
        best = mn
        while mn < epsilon && attempt < 100
            r = (rand() * (max_v - min_v)) + min_v
            mn = minimum(x -> abs.(r - x), view(out, 1:(k-1)))
            attempt += 1
            best = max(mn, best)
        end

        attempt >= 100 && @warn("Failed to generate sample w/in 500 attepmts (best = $(best))")
    
        out[k] = r
    end

    return out
end
# ============================================================================ #
function generate_case_files(param_file::AbstractString,
    counties::Vector{<:Integer}, pop::Vector{<:Integer}, n::Integer,
    odir::AbstractString=dirname(ifile), n_idx_case::Integer=20,
    n_rep::Integer=5)

    # fixed # of index cases across all runs, but uniform? or random small %?
    if true
        # in ascending sorted FIPS order for NM counties
        n_idx = [111,5,84,5,5,24,5,81,63,8,5,5,5,56,5,5,23,28,5,6,5,8,11,13,22,5,36,5,5,5,5,5,10]
    else
        n_idx = fill(n_idx_case, length(counties))
    end

    for k = 1:n

        iim = rs_rand(length(counties), 0.0, 0.5, 1e-4)

        ofile = joinpath(odir, 
            replace(
                basename(param_file),
                ".toml" => "_run_" * lpad((k-1) * n_rep, 3, '0') * ".cases"
            )
        )

        open(ofile, "w") do io
            for (k, county) in enumerate(counties)
                println(io, county, " ", n_idx[k], " ",
                    round(Int, iim[k] * pop[k])
                )
            end
        end

    end
end
# ============================================================================ #
function generate_toml_files(ifile::AbstractString, counties::Vector{<:Integer},
    n::Integer, odir::AbstractString=dirname(ifile), n_rep::Integer=5)
    
    param_str = read(ifile, String)

    N = round.(Int, log10(n) + 1)

    p_trans = round.(rs_rand(n, 0.075, 0.35, 10.0^-N), digits=5)
    p_asymptomatic = round.(rs_rand(n, 0.05, 0.95, 10.0^-N), digits=5)
    rel_trans_asymptomatic = round.(rs_rand(n, 0.05, 0.95, 10.0^-N), digits=5)
    withdrawal_scalar = round.(rs_rand(n, 0.0, 2.0, (10.0^-N) * 2), digits=5)

    for k = 1:n
    
        tmp = replace(param_str, r"run_number = \d+" => "run_number = $(k-1)")

        tmp = replace(tmp, r"p_trans = \d?\.?\d+" => "p_trans = $(p_trans[k])")
        tmp = replace(tmp, r"p_asymptomatic = \d?\.?\d+" => "p_asymptomatic = $(p_asymptomatic[k])")
        tmp = replace(tmp, r"rel_trans_asymptomatic = \d?\.?\d+" => "rel_trans_asymptomatic = $(rel_trans_asymptomatic[k])")
        tmp = replace(tmp, r"withdrawal_scalar = \d?\.?\d+" => "withdrawal_scalar = $(withdrawal_scalar[k])")

        ofile = joinpath(odir, 
            replace(basename(ifile), ".toml" => "_run_" * lpad((k-1) * n_rep, 3, '0') * ".toml")
        )

        idx_case_file = joinpath(
            "/vast/home/palexander/sandbox/nm_multi-param_sweep-all_params",
            # "/Users/palexander/Documents/emerge+radium/testing_results/nm_multi-param_sweep-all/test",
            replace(basename(ofile), ".toml" => ".cases")
        )

        tmp = replace(tmp, r"index_case_file = [^\n]+" => "index_case_file = \"$(idx_case_file)\"")

        open(ofile, "w") do io
            print(io, tmp)
            # scale = round.((rand(length(counties)) .* 0.85) .+ 0.15, digits=4)
            scale = round.((rs_rand(length(counties), 0.0, 1.0, 1e-4) .* 0.85) .+ 0.15, digits=4)
            write_policy(io, scale, counties)
        end

    end
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
    n_global = 5 # 4 + state-level im
    n_local = 2
    n_rep = 5

    param_files = map(x -> Epicast.find_files(x, r".*\.toml"), param_dirs)
    run_n = collect(Iterators.flatten(map(x -> map(get_run_number, x), param_files)))
        # map(x -> get_run_number(x) + 1000, param_files[2]))

    param_files = reduce(vcat, param_files)
    param_files .= param_files[sortperm(run_n)]

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

        # run = get_run_number(file) * 5
        run = run_n[k]# * 5

        for j in run:(run+4)
            data_dir = j < 5000 ? data_dirs[1] : data_dirs[2]
            data_file = joinpath(data_dir, "run_" * lpad(j, 3, '0') * ".bin")
            data = Epicast.preprocess!(
                Epicast.aggregate(Epicast.County, Epicast.read_runfile(data_file)),
                "total",
                smooth=true,
                diff=true,
                get_denom=Epicast.noop_denom
            )

            out[inc,:,:] .+= data.data["total"]

            params[inc,1] = p_trans
            params[inc,2] = p_asymptomatic
            params[inc,3] = rel_trans_asymptomatic
            params[inc,4] = withdrawal_scalar
            params[inc,5] = sum(im) / sum(pop)
            idx = n_global + n_geo
            params[inc,(n_global + 1):idx] .= work_scale
            params[inc,idx+1:end] .= im ./ pop

            inc += 1
        end
    end

    npzwrite(ofile, out)
    npzwrite(replace(ofile, ".npz" => "_params.npz"), params)

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
end # module EpicastCalibrate
