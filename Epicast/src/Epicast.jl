module Epicast

using DelimitedFiles, PyPlot, Colors, Statistics

const TRACT2STATE = 10^9
const TRACT2COUNTY = 10^6
const COUNTY2STATE = 10^3

export TRACT2STATE, TRACT2COUNTY, COUNTY2STATE
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
struct EpicastTable{T,N}
    index::Dict{String,Int}
    data::Array{T,N}
end
Base.getindex(x::EpicastTable{<:Any,2}, s::AbstractString) = view(x.data, :, x.index[s])#x.data[:, x.index[s]]
Base.getindex(x::EpicastTable{<:Any,3}, s::AbstractString) = view(x.data, :, x.index[s], :)#x.data[:, x.index[s], :]
function match_columns(f::Function, x::EpicastTable)
    idx = Vector{Int}(undef, 0)
    for key in keys(x.index)
        f(key) && push!(idx, x.index[key])
    end
    return sort!(idx)
end
filter_columns(f::Function, x::EpicastTable) = sort!(filter(f, collect(keys(x.index))))
# ============================================================================ #
abstract type AbstractRunData end
# ============================================================================ #
struct RunData{T} <: AbstractRunData
    demog::EpicastTable{T,2}
    data::EpicastTable{T,3}
    fips::Vector{UInt64}
    runno::String
end
rundata(x::RunData, s::AbstractString) = x.data[s]
has_data(x::RunData, s::AbstractString) = haskey(x.data.index, s)
demographics(x::RunData, s::AbstractString) = x.demog[s]
n_timepoint(x::RunData) = size(x.data.data, 3)
has_demographic(x::RunData, s::AbstractString) = haskey(x.demog.index, s)
function data_groups(x::RunData)
    names = Set{String}()
    for k in keys(x.data.index)
        push!(names, split(k, "_")[1])
    end
    return names
end
column_names(x::RunData) = collect(keys(x.data.index))
demographic_names(x::RunData) = collect(keys(x.demog.index))
# ============================================================================ #
function case_count!(out::AbstractVector, x::RunData, var::AbstractString,
    idx::AbstractVector{<:Integer})

    if has_demographic(x, var)
        # number of new cases relative to total size of that demographic
        out .= new_cases(view(rundata(x, var), idx, :))
        out ./= sum(view(demographics(x, var), idx))
        out .*= 1e5
    else
        out .= total_cases(view(rundata(x, var), idx, :))
    end
    return out
end
# ============================================================================ #
function group_by(x::RunData, name::AbstractString, conv::Integer,
    fcases!::Function=case_count!)

    if has_data(x, name)
        cols = [name]
    else
        cols = filter_columns(x -> startswith(x, name), x.data)
    end

    dat, grp = group_by(x, cols, conv, fcases!)

    return dropdmins(dat, dims=3), grp
end
# ============================================================================ #
function group_by(x::RunData, conv::Integer, fcases!::Function=case_count!)

    cols = sort!(collect(keys(x.data.index)))
    return group_by(x, cols, conv, fcases!)
end
# ============================================================================ #
function group_by(x::RunData, cols::Vector{<:AbstractString}, conv::Integer,
    fcases!::Function=case_count!)

    fips_conv = div.(x.fips, conv)
    unique_grps = sort!(unique(fips_conv))

    out = Array{Float64,3}(undef, n_timepoint(x), length(unique_grps), length(cols))    

    for k in eachindex(unique_grps)
        idx = findall(isequal(unique_grps[k]), fips_conv)
        for j in eachindex(cols)
            fcases!(view(out,:,k,j), x, cols[j], idx)
        end
    end

    return out, unique_grps
end
# ============================================================================ #
function read_runfile_header(io::IO, ::Type{T}=UInt32) where T <: Integer
    nrow = read(io, UInt64)
    ncol = read(io, UInt64)
    n_pt = read(io, UInt64)
    ncol_demog = read(io, UInt64)
    demo_len = read(io, UInt64)
    col_len = read(io, UInt64)

    # @show(Int(nrow), Int(ncol), Int(n_pt), Int(ncol_demog), Int(hdr_len))

    demo_buf = Vector{UInt8}(undef, demo_len)
    read!(io, demo_buf)
    demo_names = split(String(demo_buf), '\0', keepempty=false)

    col_buf = Vector{UInt8}(undef, col_len)
    read!(io, col_buf)
    col_names = split(String(col_buf), '\0', keepempty=false)

    # FIPS code for each tract
    fips = Vector{UInt64}(undef, nrow)
    read!(io, fips)

    # demographics for of each tract
    demo = Matrix{T}(undef, nrow, ncol_demog)
    read!(io, demo)

    return nrow, ncol, n_pt, col_names, fips,
        EpicastTable{T,2}(run_index(demo_names), demo)
end
# ============================================================================ #
function read_runfile(ifile::AbstractString, ::Type{T}=UInt32) where T <: Integer   
    m = match(r".*run_(\d+)\.bin$", ifile)
    m == nothing && error("failed to parse run number, is this a transitions file?")
    run = m[1]
    return open(ifile, "r") do io

        nrow, ncol, n_pt, col_names, fips, demo = read_runfile_header(io, T)

        pos = position(io)
        seekend(io)    
        nbytes = position(io) - pos

        @assert((nbytes % (nrow * ncol * sizeof(UInt32))) == 0,
            "invalid data block size!")
        
        seek(io, pos)

        data = Array{T,3}(undef, nrow, ncol, n_pt)
        read!(io, data)

        return RunData{T}(
            demo,
            EpicastTable{T,3}(run_index(col_names), data),
            fips,
            m[1]
        )
    end
end
# ============================================================================ #
struct AgentTransition
    agent_id::UInt64
    location_id::UInt64
    timestep::UInt16
    context::UInt8
    state::UInt8
    variant::UInt8
end
# ---------------------------------------------------------------------------- #
home_state(a::AgentTransition) = a.agent_id >> 58
agent_id(a::AgentTransition) = a.agent_id & ~(UInt64(0b0111111) << 58)
tract_fips(a::AgentTransition) = a.location_id >> 8
tract_community(a::AgentTransition) = a.location_id & 0xff

# (tract, community) in which transition to infected occured
infection_location(a::AgentTransition) = tract_fips(a), tract_community(a)
# ============================================================================ #
function read_agent_transitions(io::IO)
    pos = position(io)
    seekend(io)
    nbyte = position(io) - pos
    seek(io, pos)

    if (nbyte % sizeof(AgentTransition)) != 0
        @warn("File contains inomplete transition packets")
    end

    n_packet = div(nbyte, sizeof(AgentTransition))
    
    data = Vector{AgentTransition}(undef, n_packet)
    read!(io, data)
    
    return data
end
# ============================================================================ #
function RunData(demog::EpicastTable{2}, data::Vector{AgentTransition},
    fips::Vector{UInt64}, n_pt::Integer, run::AbstractString)

    mp = Dict{UInt64,Int}(id => k for (k,id) in enumerate(fips))

    tmp = zeros(UInt32, length(fips), 1, n_pt)

    for d in data
        r = mp[tract_fips(d)]
        s = div(d.timestep, 2) + 1

        # count files store data as cumulative sum, so compute that as we go
        view(tmp, r, 1, s:n_pt) .+= 1
    end

    return RunData(demog, EpicastTable{3}(run_index(["total"]), tmp), fips, run)
end
# ============================================================================ #
struct EventData <: AbstractRunData
    demog::EpicastTable{2}
    events::Vector{AgentTransition}
    fips::Vector{UInt64}
    n_pt::UInt64
    run::String
end
# ============================================================================ #
function read_eventfile(::Type{T}, ifile::AbstractString) where T<:AbstractRunData   
    m = match(r".*run_(\d+)\.events\.bin$", ifile)
    m == nothing && error("failed to parse run number, is this a counts file?")
    run = m[1]

    return open(ifile, "r") do io

        nrow, ncol, n_pt, col_names, fips, demog = read_runfile_header(io)

        events = read_agent_transitions(io)

        return T(demog, events, fips, n_pt, m[1])
    end
end

read_eventfile(ifile::AbstractString) = read_eventfile(RunData, ifile)
# ============================================================================ #
total_cases(x::AbstractVector{<:Real}) = x
total_cases(x::AbstractMatrix{<:Real}) = dropdims(sum(x, dims=1),dims=1)
# ============================================================================ #
function new_cases(x::AbstractMatrix{<:Real}, f::Function=sum)
    out = dropdims(f(x, dims=1),dims=1)
    out[2:end] .= diff(out)
    return out
end
function new_cases(x::AbstractVector{<:Real}, f::Function=sum)
    out = copy(x)
    out[2:end] .= diff(out)
    return out
end
mean_new_cases(x) = new_cases(x, mean)
# ============================================================================ #
function plot_all_multi(data::RunData, odir::AbstractString)
    PyPlot.ioff()

    h, ax = plot_run(data, "total", freduce=mean_new_cases, normalize=true,
        ylab="Normalized new cases per tract")

    ax.get_legend().remove()
    h.tight_layout()

    h.savefig(joinpath(odir, "total-cases_$(data.runno).png"), dpi=200)

    PyPlot.close(h)

    if has_data(data, "age_0")
    
        h, ax = subplots(1, 4)
        h.set_size_inches((16, 5))

        plot_run(data, "age", freduce=mean_new_cases, normalize=true,
            title="Age", ylab="Normalized new cases per tract", h=h, ax=ax[1])

        plot_run(data, "household", freduce=mean_new_cases, normalize=true, dropname=true,
            title="Household size", ylab="", h=h, ax=ax[2])

        plot_run(data, "race", freduce=mean_new_cases, normalize=true, dropname=true,
            title="Race", ylab="", h=h, ax=ax[3])

        plot_run(data, "ethnicity", freduce=mean_new_cases, normalize=true, dropname=true,
            title="Ethnicity", ylab="", h=h, ax=ax[4])

        h.tight_layout()

        h.savefig(joinpath(odir, "demographic-cases_$(data.runno).png"), dpi=200)

        PyPlot.close(h)
    end

    if has_data(data, "hospitalized_age0")
    
        h, ax = subplots(1,4)
        h.set_size_inches((16, 5))

        plot_run(data, "hospitalized", freduce=mean_new_cases, normalize=false, dropname=true,
            title="Hospitalization", ylab="Daily incidence per tract", h=h, ax=ax[1])

        plot_run(data, "icu", freduce=mean_new_cases, normalize=false, dropname=true,
            title="ICU admitance", ylab="", h=h, ax=ax[2])

        plot_run(data, "ventilated", freduce=mean_new_cases, normalize=false, dropname=true,
            title="Ventilation", ylab="", h=h, ax=ax[3])
        
        plot_run(data, "dead", freduce=mean_new_cases, normalize=false, dropname=true,
            title="Deaths", ylab="", h=h, ax=ax[4])

        h.tight_layout()

        h.savefig(joinpath(odir, "treatment-status_$(data.runno).png"), dpi=200)

        PyPlot.close(h)
    end

    if has_data(data, "infection-src_family")

        h, ax = plot_run(data, "infection-src", freduce=mean_new_cases, normalize=false, dropname=true,
            title="Infection context", ylab="Daily cases per tract")

        h.tight_layout()

        h.savefig(joinpath(odir, "infection-context_$(data.runno).png"), dpi=200)

        PyPlot.close(h)
    end

    if has_data(data, "status_incubation")
        h, ax = plot_run(data, "status", freduce=total_cases, normalize=false, dropname=true,
            title="Infection status", ylab="# of agents")

        h.tight_layout()

        h.savefig(joinpath(odir, "infection-status_$(data.runno).png"), dpi=200)

        PyPlot.close(h)

    end

    PyPlot.ion()

    return nothing
end
# ============================================================================ #
function plot_all(data::RunData, ofile::AbstractString; plot_status::Bool=false)
    PyPlot.ioff()
    h, ax = subplots(3, 3)
    h.set_size_inches((20,18))

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
        title="Infection context", ylab="Daily cases per tract", h=h, ax=ax[3,3])

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

    return nothing
end
# ============================================================================ #
function plot_states(data::RunData, ofile::AbstractString; title::AbstractString="")
    
    # 10^9 -> FIPS tract to state conversion
    dat, lab = group_by(data, "total", TRACT2STATE, case_count!)

    PyPlot.ioff()
    h, ax = subplots(1,1)
    h.set_size_inches((14,9))

    colors = map(col -> (red(col), green(col), blue(col)), 
        distinguishable_colors(length(lab), [RGB(1,1,1), RGB(0,0,0)],
            dropseed=true)
    )

    foreach(1:length(lab)) do k
        ax.plot(dat[:,1,k] .* 1e5, color=colors[k], label=STATE_FIPS[lab[k]])
    end
    
    ax.legend(frameon=false)

    ax.set_xlabel("Simulation day", fontsize=14)
    ax.set_ylabel("Daily infections per 100k residents", fontsize=14)

    ax.spines["top"].set_visible(false)
    ax.spines["right"].set_visible(false)
    
    # "Daily new cases per state (22 states, ~136M agents)"
    !isempty(title) && ax.set_title(title, fontsize=18)

    h.tight_layout()

    h.savefig(ofile, dpi=200)

    close(h)

    PyPlot.ion()

    return nothing
end
# ============================================================================ #
function plot_run(data::RunData, name::AbstractString; freduce=total_cases,
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
        if normalize && has_demographic(data, col)
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
const STATE_FIPS = Dict(
    1 => "AL", 2 => "AK", 4 => "AZ", 5 => "AR", 6 => "CA", 8 => "CO", 9 => "CT",
    10 => "DE", 11 => "DC", 12 => "FL", 13 => "GA", 15 => "HI", 16 => "ID",
    17 => "IL", 18 => "IN", 19 => "IA", 20 => "KS", 21 => "KY", 22 => "LA",
    23 => "ME", 24 => "MD", 25 => "MA", 26 => "MI", 27 => "MN", 28 => "MS",
    29 => "MO", 30 => "MT", 31 => "NE", 32 => "NV", 33 => "NH", 34 => "NJ",
    35 => "NM", 36 => "NY", 37 => "NC", 38 => "ND", 39 => "OH", 40 => "OK",
    41 => "OR", 42 => "PA", 44 => "RI", 45 => "SC", 46 => "SD", 47 => "TN",
    48 => "TX", 49 => "UT", 50 => "VT", 51 => "VA", 53 => "WA", 54 => "WV",
    55 => "WI", 56 => "WY"
)
end
