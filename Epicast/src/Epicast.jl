module Epicast

using DelimitedFiles, PyPlot, Colors, Statistics

# ============================================================================ #
function epi_plot(idir::AbstractString, run::Integer=2)

    # # of new symptomatic cases per day total (across tracts)
    n_symp, n_symp_age, n_symp_hh = parse_logfile(idir, run)


    # pop: population of each tract
    # n_infected: # of agents in infected state (inc | pro | pox) per day per tract
    # n_symp: cumulative # of symptoimatic cases per day per tract
    pop, n_infected, n_symp_tract = parse_attackfiles(idir, run)

    h, ax = subplots(1, 4)
    h.set_size_inches((15,5))

    foreach(ax) do cax
        cax.spines["top"].set_visible(false)
        cax.spines["right"].set_visible(false)
        cax.set_xlabel("Simulation day", fontsize=14)
    end

    t = 1:length(n_symp)

    ax[1].plot(t, n_symp, linewidth=3)

    ax[1].set_ylabel("New symptomatic cases", fontsize=14)
    ax[1].set_title("Total new cases", fontsize=16)

    for k in 1:size(n_symp_age, 2)
        ax[2].plot(t, n_symp_age[:,k], linewidth=2, label="age $(k-1)")
    end

    ax[2].set_title("New cases by age", fontsize=16)
    ax[2].legend(frameon=false, fontsize=12)

    for k in 1:size(n_symp_hh, 2)
        lab = k == 1 ? "household size $(k)" : string(k)
        ax[3].plot(t, n_symp_hh[:,k], linewidth=2, label=lab)
    end

    ax[3].set_title("... by household size", fontsize=16)
    ax[3].legend(frameon=false, fontsize=12)


    ax[4].plot(t, n_symp_tract ./ pop', linewidth=1.5)

    ax[4].set_ylabel("Proportion of agents", fontsize=14)
    ax[4].set_title("# symptomatic per tract", fontsize=16)


    h.tight_layout()

    # return h, ax
    return n_symp, n_symp_age, n_symp_hh

end
# ============================================================================ #
function parse_logfile(idir::AbstractString, run::Integer=2)

    ifile = joinpath(idir, "Log" * string(run))

    d1 = Vector{Vector{Float64}}(undef, 0)
    d2 = Vector{Vector{Float64}}(undef, 0)

    for line in eachline(ifile)
        if startswith(line, r"\s*\d+\.\d\s+")
            raw = map(split(strip(line), r"\s+")) do x
                y = tryparse(Float64, x)
                return y == nothing ? NaN : y
            end
            push!(d1, raw)
        elseif startswith(line, r"S\s+\d+\.\d\s+")
            raw = map(split(line, r"\s+")[2:end]) do x
                y = tryparse(Float64, x)
                return y == nothing ? -1 : y
            end
            push!(d2, raw)
        end
    end

    t1 = hcat(d1...)'
    t2 = hcat(d2...)'

    # log files can have repeated time points due to restarts, when multiple
    # time points appear in the file use only the latest
    len = Int(t1[end,1])
    n_symptomatic = zeros(len)
    n_symp_age = zeros(len, 5)
    n_symp_hh = zeros(len, size(t2,2) - 6)

    for k in 1:size(t1, 1)
        idx = floor(Int, t1[k,1])
        if idx > 0 && idx == t1[k,1]
            idx2 = findfirst(isequal(idx), t2[:,1])

            n_symptomatic[idx] = t1[k,2]
            n_symp_age[idx,:] = t2[idx2,2:6]
            n_symp_hh[idx,:] = t2[idx2,7:end]
        end
    end
    
    return n_symptomatic, n_symp_age, n_symp_hh
end
# ============================================================================ #
function attack_cmp_timepoint(a::AbstractString, b::AbstractString)
    at = tryparse(Int, splitext(a)[2][2:end])
    bt = tryparse(Int, splitext(b)[2][2:end])
    return at < bt
end
# ---------------------------------------------------------------------------- #
function parse_attackfiles(idir::AbstractString, run::Integer=2)
    pat = Regex("Attack" * string(run) * "\\.\\d{3}\$")
    files = find_files(idir, pat)
    sort!(files, lt=attack_cmp_timepoint)
    
    data = map(parse_attackfile, files)
    
    n_tract = Int(maximum(vcat(getindex.(data, :, 1)...))) + 1

    pop = zeros(n_tract)
    # # agents in infected state: time x tract
    n_infected = fill(0.0, (length(data), n_tract))

    # cumulative # of symptomatic agents: time x tract
    n_symp = fill(0.0, (length(data), n_tract))

    for k in eachindex(data)
        idx = Int.(data[k][:,1]) .+ 1
        pop[idx] .= data[k][:,3]
        n_infected[k,idx] .= sum(data[k][:,5:7], dims=2)
        n_symp[k,idx] .= data[k][:,4]
    end

    return pop, n_infected, n_symp
end
# ---------------------------------------------------------------------------- #
function parse_attackfile(ifile::AbstractString)
    return readdlm(ifile, skipstart=1)
end
# ============================================================================ #
function run_index(names::AbstractVector{<:AbstractString})
    return Dict{String,Int}(y => x for (x,y) in enumerate(names))
end
# ============================================================================ #
struct EpicastTable{N}
    index::Dict{String,Int}
    data::Array{UInt16,N}
end
Base.getindex(x::EpicastTable{2}, s::AbstractString) = x.data[:, x.index[s]]
Base.getindex(x::EpicastTable{3}, s::AbstractString) = x.data[:, x.index[s], :]
function match_columns(f::Function, x::EpicastTable)
    idx = Vector{Int}(undef, 0)
    for key in keys(x.index)
        f(key) && push!(idx, x.index[key])
    end
    return sort!(idx)
end
filter_columns(f::Function, x::EpicastTable) = sort!(filter(f, collect(keys(x.index))))
# ============================================================================ #
struct RunData
    demog::EpicastTable{2}
    data::EpicastTable{3}
end
rundata(x::RunData, s::AbstractString) = x.data[s]
demographics(x::RunData, s::AbstractString) = x.demog[s]
n_timepoint(x::RunData) = size(x.data.data, 3)
is_demographic(x::RunData, s::AbstractString) = haskey(x.demog.index, s)
function data_groups(x::RunData)
    names = Set{String}()
    for k in keys(x.data.index)
        push!(names, split(k, "_")[1])
    end
    return names
end
# ============================================================================ #
function read_runfile(ifile::AbstractString)
    return open(ifile, "r") do io
        nrow = read(io, UInt64)
        ncol = read(io, UInt64)
        n_pt = read(io, UInt64)
        ncol_demog = read(io, UInt64)
        hdr_len = read(io, UInt64)

        # @show(Int(nrow), Int(ncol), Int(n_pt), Int(ncol_demog), Int(hdr_len))

        names_buf = Vector{UInt8}(undef, hdr_len)
        read!(io, names_buf)
        names = split(String(names_buf), '\0', keepempty=false)

        # demographics for of each tract
        demo = Matrix{UInt16}(undef, nrow, ncol_demog)
        read!(io, demo)

        pos = position(io)
        seekend(io)    
        nbytes = position(io) - pos

        @assert((nbytes % (nrow * ncol * sizeof(UInt16))) == 0,
            "invalid data block size!")
        
        seek(io, pos)

        data = Array{UInt16,3}(undef, nrow, ncol, n_pt)

        for k in 1:n_pt
            read!(io, view(data, :, :, k))
        end

        return RunData(
            EpicastTable{2}(run_index(names[5:(5+ncol_demog-1)]), demo),
            EpicastTable{3}(run_index(names), data),
        )
    end
end
# ============================================================================ #
total_cases(x::Matrix{<:Real}) = dropdims(sum(x, dims=1),dims=1)
# ============================================================================ #
function new_cases(x::Matrix{<:Real}, f::Function=sum)
    out = dropdims(f(x, dims=1),dims=1)
    out[2:end] .= diff(out)
    return out
end
mean_new_cases(x) = new_cases(x, mean)
# ============================================================================ #
function plot_all(data::RunData, ofile::AbstractString; plot_status::Bool=false)
    PyPlot.ioff()
    h, ax = subplots(3, 3)
    h.set_size_inches((16,18))

    plot_run(data, "age", freduce=mean_new_cases, normalize=true,
        title="Age", ylab="Daily cases per tract", h=h, ax=ax[1,1])

    plot_run(data, "household", freduce=mean_new_cases, normalize=true, dropname=true,
        title="Household size", ylab="Daily cases per tract", h=h, ax=ax[1,2])

    plot_run(data, "race", freduce=mean_new_cases, normalize=true, dropname=true,
        title="Race", ylab="Daily cases per tract", h=h, ax=ax[1,3])

    plot_run(data, "ethnicity", freduce=mean_new_cases, normalize=true, dropname=true,
        title="Ethnicity", ylab="Daily cases per tract", h=h, ax=ax[2,1])

    plot_run(data, "hospitalized", freduce=mean_new_cases, normalize=false, dropname=true,
        title="Hospitalization", ylab="Daily incidence per tract", h=h, ax=ax[2,2])

    plot_run(data, "icu", freduce=mean_new_cases, normalize=false, dropname=true,
        title="ICU admitance", ylab="Daily incidence per tract", h=h, ax=ax[2,3])

    plot_run(data, "ventilated", freduce=mean_new_cases, normalize=false, dropname=true,
        title="Ventilation", ylab="Daily incidence per tract", h=h, ax=ax[3,1])
    
    plot_run(data, "dead", freduce=mean_new_cases, normalize=false, dropname=true,
        title="Deaths", ylab="Daily incidence per tract", h=h, ax=ax[3,2])

    plot_run(data, "infection-src", freduce=mean_new_cases, normalize=false, dropname=true,
        title="Infection context", ylab="Daily incidence per tract", h=h, ax=ax[3,3])

    h.tight_layout()

    h.savefig(ofile, dpi=200)

    PyPlot.close(h)

    if plot_status
        h, ax = plot_run(data, "status", freduce=total_cases, normalize=false, dropname=true,
            title="Agent status", ylab="# of agents")

        h.tight_layout()

        ext = splitext(ofile)[2]
        h.savefig(replace(ofile, ext => "_agent-status" * ext), dpi=200)
        PyPlot.close(h)
    end

    PyPlot.ion()
end
# ============================================================================ #
function plot_run(data::RunData, name::AbstractString; freduce=x->total_cases,
    normalize::Bool=false, dropname::Bool=false, ylab::String="",
    title::String="", h=nothing, ax=nothing)

    cols = filter_columns(x -> startswith(x, name), data.data)
    if h == nothing || ax == nothing
        h, ax = subplots(1, 1)
        h.set_size_inches((8,6))
    end

    t = 1:n_timepoint(data)

    rm = dropname ? name * "_" => "" : "_" => " "

    colors = map(col -> (red(col), green(col), blue(col)), 
        distinguishable_colors(length(cols), [RGB(1,1,1), RGB(0,0,0)],
            dropseed=true)
    )
    
    for (k, col) in enumerate(cols)
        dat = Float64.(rundata(data, col))
        if normalize && is_demographic(data, col)
            tmp2 = demographics(data, col)
            dat ./= tmp2
            replace!(x -> isnan(x) || isinf(x) ? 0.0 : x, dat)
        end
        tmp = freduce(dat)
        ax.plot(t, tmp, label=replace(col, rm), color=colors[k])
    end

    ax.spines["top"].set_visible(false)
    ax.spines["right"].set_visible(false)
    ax.set_xlabel("Simulation day", fontsize=14)

    !isempty(ylab) && ax.set_ylabel(ylab, fontsize=14)
    !isempty(title) && ax.set_title(title, fontsize=18)

    ax.legend(fontsize=14, frameon=false)

    h.tight_layout()

    return h, ax
end
# ============================================================================ #
find_files(dir::AbstractString, re=r".*") = return do_match(dir, re, isfile)
# ---------------------------------------------------------------------------- #
find_directories(dir::AbstractString, re=r".*") = return do_match(dir, re, isdir)
# ============================================================================ #
function do_match(dir::AbstractString, re::Regex, f::Function)
    if !isdir(dir)
        error("Input is not a vaild directory path")
    end
    files = [joinpath(dir, x) for x in readdir(dir)]
    return filter(x->occursin(re, x) && f(x), files)
end
# ============================================================================ #
end
