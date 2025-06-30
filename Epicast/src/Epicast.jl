module Epicast

using Statistics

using EpicastTables

import Base

const SignedType = Union{AbstractFloat,Signed}

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
function run_index(items::AbstractVector{T}) where T
    return Dict{T,Int}(y => x for (x,y) in enumerate(items))
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
# ============================================================================ #
end
