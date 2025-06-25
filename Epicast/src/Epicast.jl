module Epicast

using DelimitedFiles, PyPlot, Colors, Statistics

using EpicastTables

import Base

const SignedType = Union{AbstractFloat,Signed}

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
function run_index(items::AbstractVector{T}) where T
    return Dict{T,Int}(y => x for (x,y) in enumerate(items))
end
# ============================================================================ #
abstract type AbstractRunData end
# ============================================================================ #
struct RunData{G<:AbstractGeo,T1<:Real,T2<:Real,I<:Integer} <: AbstractRunData
    demog::FIPSTable{G,T1,2}
    data::FIPSTable{G,T2,3}
    fips::Vector{I}
    runno::String
end
# ============================================================================ #
rundata(x::RunData, s::AbstractString) = x.data[s]
has_data(x::RunData, s::AbstractString) = EpicastTables.has_var(x.data, s)
demographics(x::RunData, s::AbstractString) = x.demog[s]
n_timepoint(x::RunData) = size(x.data.data, 1)
has_demographic(x::RunData, s::AbstractString) = EpicastTables.has_var(x.demog, s)
function data_groups(x::RunData)
    names = Set{String}()
    for k in keys(x.data.index)
        push!(names, split(k, "_")[1])
    end
    return names
end
column_names(x::RunData) = EpicastTables.all_vars(x.data)
demographic_names(x::RunData) = EpicastTables.all_vars(x.demog)
# ============================================================================ #
# noop if geographies are the same and datatype is float
aggregate(::Type{G}, d::RunData{G,<:Real,Float64}, f::Function=sum) where {G<:AbstractGeo} = d
# ---------------------------------------------------------------------------- #
# convert to Float64 if not already
function aggregate(::Type{G}, d::RunData{G,<:Real,<:Real}, f::Function=sum) where {G<:AbstractGeo}
    return RunData(d.demog, convert_datatype(Float64, d.data), d.fips, d.runno)
end
# ---------------------------------------------------------------------------- #
function aggregate(::Type{G}, d::RunData, f::Function=sum) where {G<:AbstractGeo}
    # I don't think it would ever make sense to "mean" the demographic data when
    # aggregating... we'll see
    demog = EpicastTables.aggregate(G, d.demog, sum)
    data = EpicastTables.aggregate(G, d.data, f)
    return RunData(demog, data, EpicastTables.all_fips(data), d.runno)
end
# ============================================================================ #
function smooth!(d::RunData{G,L,T}, var::AbstractString, n::Integer=7,
    f::Function=mean) where {G,L,T<:AbstractFloat}

    EpicastTables.smooth!(d.data, var, n, f)
    return d
end
# ---------------------------------------------------------------------------- #
function smooth!(d::RunData{G,L,T}, vars::AbstractVector{<:AbstractString},
    n::Integer=7, f::Function=mean) where {G,L,T<:AbstractFloat}

    foreach(var -> EpicastTables.smooth!(d.data, var, n, f), vars)
    return d
end
# ---------------------------------------------------------------------------- #
function smooth(d::RunData, vars::AbstractVector{<:AbstractString},
    n::Integer=7, f::Function=mean)

    return smooth!(convert_datatype(Float64, d), var, n, f)
end
# ============================================================================ #
function diff!(d::RunData{G,L,T}, var::AbstractString) where {G,L,T<:SignedType}
    tmp = d.data[var]
    N = size(tmp,1)
    tmp[2:N,:] .= Base.diff(tmp, dims=1)
    return d
end
# ---------------------------------------------------------------------------- #
function diff!(d::RunData{G,L,T}, vars::AbstractVector{<:AbstractString}) where {T<:SignedType,G,L}
    foreach(var -> diff!(d, var), vars)
    return d
end
# ---------------------------------------------------------------------------- #
function diff(d::RunData, vars::AbstractVector{<:AbstractString},
    ::Type{T}=Float64) where T<:SignedType

    return diff!(convert_datatype(T, d), vars)
end
# ============================================================================ #
function default_denom(d::RunData, var::AbstractString)
    if has_demographic(d, denom)
        denom_data = d.demog[denom]
    elseif has_data(d, denom)
        denom_data = d.data[denom]
    else
        error("variable \"$(denom)\" does not exist in given RunData")
    end

    return denom_data
end
# ============================================================================ #
function normalize!(d::RunData{G,L,T}, var::AbstractString,
    get_denom::Function=default_denom) where {G,L,T<:AbstractFloat}

    d.data[var] ./= reshape(get_denom(d, var), 1, :)
    return d
end
# ============================================================================ #
function preprocess!(d::RunData{G,L,T}, var::AbstractString; smooth::Bool=false,
    diff::Bool=false, get_denom::Function=x->default_denom(d, var)) where {G,L,T<:AbstractFloat}

    smooth && smooth!(d, var)
    diff && diff!(d, var)

    return normalize!(d, var, get_denom)
end
# ---------------------------------------------------------------------------- #
function preprocess!(d::RunData{G,L,T}, vars::AbstractVector{<:AbstractString};
    smooth::Bool=false, diff::Bool=false,
    get_denom::Function=x->default_denom(d, var)) where {G,L,T<:AbstractFloat}

    for var in vars
        preprocess!(d, var, smooth=smooth, diff=diff, get_denom=get_denom)
    end

    return d
end
# ---------------------------------------------------------------------------- #
function preprocess!(d::RunData{G,L,T}, fmatch::Function; smooth::Bool=false,
    diff::Bool=false, get_denom::Function=x->default_denom(d, var)) where {G,L,T<:AbstractFloat}

    for name in column_names(d)
        fmatch(name) && preprocess!(d, name, smooth=smooth, diff=diff,
            get_denom=get_denom)
    end

    return d
end
# ---------------------------------------------------------------------------- #
function preprocess!(d::RunData{G,L,T}, pat::Regex; smooth::Bool=false,
    diff::Bool=false, get_denom::Function=x->default_denom(d, var)) where {G,L,T<:AbstractFloat}

    preprocess!(d, x -> occursin(pat, x), smooth=smooth, diff=diff,
        get_denom=get_denom)

    return d
end
# ---------------------------------------------------------------------------- #
function preprocess(::Type{G}, d::RunData, args...) where {G<:AbstractGeo}
    return preprocess!(aggregate(G, d), args...)
end
# ============================================================================ #
function case_count!(out::AbstractVector, x::RunData, var::AbstractString,
    idx::AbstractVector{<:Integer}, norm::AbstractString="")

    if norm == ""
        if has_demographic(x, var)
            norm = var
        end
    end


    if has_demographic(x, norm)
        # number of new cases relative to total size of that demographic
        out .= new_cases(view(rundata(x, var), idx, :))
        out ./= sum(view(demographics(x, norm), idx))
        #out .*= 1e5
    else
        out .= total_cases(view(rundata(x, var), idx, :))
    end
    return out
end
# ============================================================================ #
function group_by(x::RunData, name::AbstractString, conv::Integer,
    fcases!::Function=case_count!, norms::Dict{String,String}=Dict{String,String}())

    if has_data(x, name)
        cols = [name]
    else
        cols = filter_vars(x -> startswith(x, name), x.data)
    end

    dat, grp = group_by(x, cols, conv, fcases!, norms)

    return dropdmins(dat, dims=3), grp
end
# ============================================================================ #
function group_by(x::RunData, conv::Integer, fcases!::Function=case_count!,
    norms::Dict{String,String}=Dict{String,String}())

    cols = sort!(collect(keys(x.data.index)))
    return group_by(x, cols, conv, fcases!, norms)
end
# ============================================================================ #
function group_by(x::RunData, cols::Vector{<:AbstractString}, conv::Integer,
    fcases!::Function=case_count!, norms::Dict{String,String}=Dict{String,String}())

    fips_conv = div.(x.fips, conv)
    unique_grps = sort!(unique(fips_conv))

    out = Array{Float64,3}(undef, n_timepoint(x), length(unique_grps), length(cols))

    for k in eachindex(unique_grps)
        idx = findall(isequal(unique_grps[k]), fips_conv)
        for j in eachindex(cols)
          fcases!(view(out,:,k,j), x, cols[j], idx, get(norms, cols[j], ""))
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
    demo_names = map(string, split(String(demo_buf), '\0', keepempty=false))

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
        FIPSTable{Tract,T,2}(run_index(fips), run_index(demo_names), demo)
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

        return RunData{Tract,T,T,UInt64}(
            demo,
            FIPSTable{Tract,T,3}(
                run_index(fips),
                run_index(map(string, col_names)),
                permutedims(data, (3,1,2))
            ),
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

# (tract, community) in which transition occured
transition_location(a::AgentTransition) = tract_fips(a), tract_community(a)
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
function RunData(demog::FIPSTable{G,T,2}, data::Vector{AgentTransition},
    fips::Vector{UInt64}, n_pt::Integer, run::AbstractString) where {G,T}

    mp = Dict{UInt64,Int}(id => k for (k,id) in enumerate(fips))

    tmp = zeros(T, n_pt, length(fips), 1)

    # count files do not include index cases
    for d in filter(x -> x.state == 0x01 && x.context != 0xff, data)
        r = mp[tract_fips(d)]
        s = div(d.timestep, 2) + 1

        # count files store data as cumulative sum, so compute that as we go
        view(tmp, s:n_pt, r, 1) .+= 1
    end

    return RunData{G,T,T,UInt64}(
        demog,
        FIPSTable{G,T,3}(run_index(fips), run_index(["total"]), tmp),
        fips,
        run
    )
end
# ============================================================================ #
struct EventData <: AbstractRunData
    demog::FIPSTable{Tract,UInt32,2}
    events::Vector{AgentTransition}
    fips::Vector{UInt64}
    n_pt::UInt64
    run::String
end
# ---------------------------------------------------------------------------- #
function Base.:(==)(a::EventData, b::EventData)
    return a.demog == b.demog && a.events == b.events && a.fips == b.fips &&
        a.n_pt == b.n_pt
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
function read_rundata(ifile::AbstractString, ::Type{T}=UInt32) where T<:Integer
    return endswith(ifile, ".events.bin") ? read_eventfile(ifile) :
        read_runfile(ifile, T)
end
# ============================================================================ #
total_cases(x::AbstractVector{<:Real}) = x
total_cases(x::AbstractMatrix{<:Real}) = dropdims(sum(x, dims=2),dims=2)
# ============================================================================ #
function new_cases(x::AbstractMatrix{<:Real}, f::Function=sum)
    out = dropdims(f(x, dims=2),dims=2)
    out[2:end] .= Base.diff(out)
    return out
end
function new_cases(x::AbstractVector{<:Real}, f::Function=sum)
    out = copy(x)
    out[2:end] .= Base.diff(out)
    return out
end
mean_new_cases(x) = new_cases(x, mean)
# ============================================================================ #
function plot_all_multi(data::RunData, odir::AbstractString; ext::AbstractString="png")
    PyPlot.ioff()

    h, ax = plot_run(data, "total", freduce=mean_new_cases, normalize=true,
        ylab="Normalized new cases per tract")

    ax.get_legend().remove()
    h.tight_layout()

    ext = startswith(ext, '.') ? ext[2:end] : ext

    h.savefig(joinpath(odir, "total-cases_$(data.runno)." * ext), dpi=200)

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

        h.savefig(joinpath(odir, "demographic-cases_$(data.runno)." * ext), dpi=200)

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

        h.savefig(joinpath(odir, "treatment-status_$(data.runno)." * ext), dpi=200)

        PyPlot.close(h)
    end

    if has_data(data, "infection-src_family")

        h, ax = plot_run(data, "infection-src", freduce=mean_new_cases, normalize=false, dropname=true,
            title="Infection context", ylab="Daily cases per tract")

        h.tight_layout()

        h.savefig(joinpath(odir, "infection-context_$(data.runno)." * ext), dpi=200)

        PyPlot.close(h)
    end

    if has_data(data, "status_incubation")
        h, ax = plot_run(data, "status", freduce=total_cases, normalize=false, dropname=true,
            title="Infection status", ylab="# of agents")

        h.tight_layout()

        h.savefig(joinpath(odir, "infection-status_$(data.runno)." * ext), dpi=200)

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
    dat, lab = group_by(data, "total", EpicastTables.TRACT2STATE, case_count!)

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
function plot_run(data::Vector{RunData}, name::AbstractString; args...)

    h, ax = (nothing,nothing)

    for d in data
        h, ax = plot_run(d, name; h=h, ax=ax, args...)
    end

    return h, ax
end
# ============================================================================ #
function plot_runs(data, name::AbstractString; freduce=total_cases, reduce_cols=nothing,
    normalize::Bool=false, demo::AbstractString="", dropname::Bool=false, ylab::String="",
    title::String="", h=nothing, ax=nothing, style::Function=(h, ax) -> nothing,
    dataset_names::AbstractVector{<:AbstractString}=nothing)

    cols = map(d -> filter_vars(x -> startswith(x, name), d.data), data)
    if h == nothing || ax == nothing
        h, ax = subplots(1, 1)
        h.set_size_inches((8,6))
    end

    t = 1:maximum(map(n_timepoint, data))

    rm = dropname ? name * "_" => "" : "_" => " "

    n_colors = sum(map(length, cols))
    if reduce_cols != nothing
        n_colors = sum(map(x -> 1, cols))
    end
    colors = map(col -> (red(col), green(col), blue(col)),
                 distinguishable_colors(n_colors, [RGB(1,1,1), RGB(0,0,0)], dropseed=true))

    plot_col(k, dataset, dat, name, col) = begin
        d = demo
        if 0 == cmp(demo, "")
            d = col
        end
        if normalize && has_demographic(dataset, d)
            tmp2 = sum(demographics(dataset, d))
            dat ./= tmp2
            #dat .*= 1e5

            replace!(x -> isnan(x) || isinf(x) ? 0.0 : x, dat)
        end
        tmp = freduce(dat)
        label = replace(col, rm)
        if dataset_names != nothing
            if length(colors) != length(data)
                label = "$(name)_$label"
            else
                label = name
            end
        end
        ax.plot(t, tmp, label=label, color=colors[k])
    end

    k = 1
    for (j, dataset_cols) in enumerate(cols)
        dataset = data[j]
        n = dataset_names[j]
        if reduce_cols == nothing
            for (c, col) in enumerate(dataset_cols)
                d = Float64.(rundata(dataset, col))
                plot_col(k, dataset, d, n, col)
                k += 1
            end
        else
            d = reduce_col(map(col -> Float64.(rundata(dataset, col)), dataset_cols))
            plot_col(k, dataset, d, n, name)
            k += 1
        end
    end

    ax.spines["top"].set_visible(false)
    ax.spines["right"].set_visible(false)
    ax.set_xlabel("Simulation day", fontsize=14)

    !isempty(ylab) && ax.set_ylabel(ylab, fontsize=14)
    !isempty(title) && ax.set_title(title, fontsize=18)

    ax.legend(fontsize=14, frameon=false)
    style(h, ax)

    h.tight_layout()

    return h, ax
end
# ============================================================================ #
function plot_run(data::RunData, name::AbstractString; freduce=total_cases,
    normalize::Bool=false, dropname::Bool=false, ylab::String="",
    title::String="", h=nothing, ax=nothing, style::Function=(h, ax) -> nothing,
    agg_level::Type{<:AbstractGeo}=State, smooth::Bool=false)

    cols = filter_vars(x -> startswith(x, name), data.data)
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
    tmp = EpicastTables.aggregate(agg_level, data.data)
    for (k, col) in enumerate(cols)
        if smooth
            EpicastTables.smooth!(tmp, col, 7, mean)
        end
        cases = freduce(tmp[col])
        if normalize && has_demographic(data, name)
            pop = EpicastTables.aggregate(agg_level, data.demog)[name]
            cases ./= (pop ./ 1e5)
            replace!(x -> isnan(x) || isinf(x) ? 0.0 : x, cases)
        end
        ax.plot(t, cases, label=replace(col, rm), color=colors[k])
    end

    ax.spines["top"].set_visible(false)
    ax.spines["right"].set_visible(false)
    ax.set_xlabel("Simulation day", fontsize=14)

    !isempty(ylab) && ax.set_ylabel(ylab, fontsize=14)
    !isempty(title) && ax.set_title(title, fontsize=18)

    ax.legend(fontsize=14, frameon=false)
    style(h, ax)

    h.tight_layout()

    return h, ax
end
# ============================================================================ #
function plot_infection_src(data::RunData)
    fields = ["family","work","school","neighborhood-cluster","neighborhood"]
    total = Epicast.new_cases(Epicast.rundata(data, "total"))
    d = zeros(Float64, Epicast.n_timepoint(data), length(fields))
    for k in eachindex(fields)
        d[:,k] .= Epicast.new_cases(Epicast.rundata(data, "infection-src_" * fields[k]))
    end
    d ./= total

    h, ax = subplots(1,1)

    h.set_size_inches((10,6))

    ax.plot(d, linewidth=2.5)

    ax[:spines]["right"].set_visible(false)
    ax[:spines]["top"].set_visible(false)

    ax[:spines]["left"].set_linewidth(2.5)
    ax[:spines]["bottom"].set_linewidth(2.5)

    ax.xaxis.set_tick_params(width=2.5)
    ax.yaxis.set_tick_params(width=2.5)

    ax.legend(fields, fontsize=12, frameon=false, loc="upper left",
        bbox_to_anchor=(0.9,1.0))

    ax.set_ylabel("Proportion of total cases", fontsize=14)
    ax.set_xlabel("Simulation time (days)", fontsize=14)

    h.tight_layout()

    return h, ax
end
# ============================================================================ #
function aggregate(data::RunData, col::AbstractString;
    freduce=total_cases, demo::AbstractString="")

    dat = Float64.(rundata(data, col))
    if demo != "" && has_demographic(data, demo)
        tmp2 = demographics(data, demo)
        #dat .*= 1e5

        replace!(x -> isnan(x) || isinf(x) ? 0.0 : x, dat)
        replace!(x -> isnan(x) || isinf(x) ? 0.0 : x, tmp2)
        return freduce(dat) / sum(tmp2)
    else
        return freduce(dat)
    end
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
