module EpicastTables

export FIPSTable, aggregate_state, aggregate_county, aggregate_tract,
    all_states, all_counties, all_tracts, filter_vars

export AbstractGeo, BlockGroup, Tract, County, State
# ============================================================================ #
const BG2TRACT = 10
const BG2COUNTY = 10^7
const BG2STATE = 10^10
const TRACT2STATE = 10^9
const TRACT2COUNTY = 10^6
const COUNTY2STATE = 10^3
# ============================================================================ #
abstract type AbstractGeo end
struct BlockGroup <: AbstractGeo end
struct Tract <: AbstractGeo end
struct County <: AbstractGeo end
struct State <: AbstractGeo end
# ============================================================================ #
struct FIPSTable{G<:AbstractGeo,T,N}
    fips_index::Dict{UInt64,Int}
    var_index::Dict{String,Int}
    data::Array{T,N}# time x fips x var | fips x var | fips x 1
end
# ---------------------------------------------------------------------------- #
FIPSTable(::Val{N}=Val(1)) where N = FIPSTable{Tract,Float64,N}(Dict{UInt64,Int}(), Dict{String,Int}(), Array{Float64,N}(undef, tuple(zeros(Int,N)...)))
# ---------------------------------------------------------------------------- #
function FIPSTable(::Type{G}, data::Array{T,1}, fips_index::Dict{UInt64,Int},
    var::AbstractString) where {T,G<:AbstractGeo}

    return FIPSTable{G,T,1}(fips_index, Dict{String,Int}(var => 1), data)
end
# ---------------------------------------------------------------------------- #
function FIPSTable(::Type{G}, data::Array{T,N}, fips_index::Dict{UInt64,Int},
    var_index::Dict{String,Int}) where {T,N,G<:AbstractGeo}

    return FIPSTable{G,T,N}(fips_index, var_index, data)
end
# ---------------------------------------------------------------------------- #
function FIPSTable(::Type{G}, data::Array{T,N}, fips::AbstractVector{<:Integer},
    vars::AbstractVector{<:AbstractString}) where {T,N,G<:AbstractGeo}

    fips_index = Dict{UInt64,Int}(UInt64(x) => k for (k,x) in enumerate(fips))
    var_index = Dict{String,Int}(x => k for (k,x) in enumerate(vars))

    return FIPSTable{G,T,N}(fips_index, var_index, data)
end
# ---------------------------------------------------------------------------- #
all_fips(tbl::FIPSTable) = collect(keys(tbl.fips_index))
has_fips(tbl::FIPSTable, fips::Integer) = haskey(tbl.fips_index, fips)
all_vars(tbl::FIPSTable) = collect(keys(tbl.var_index))
has_var(tbl::FIPSTable, var::AbstractString) = haskey(tbl.var_index, var)
filter_vars(f::Function, tbl::FIPSTable) = sort!(filter(f, all_vars(tbl)))

data_array(tbl::FIPSTable) = tbl.data
# --------------------------------------------------------------------------- #
@inline function check_var(tbl::FIPSTable{G,T,1}, var::AbstractString) where {G,T}
    @assert(haskey(tbl.var_index, var), "variable $(var) not found")
    return true
end
Base.getindex(tbl::FIPSTable{G,T,1}, fips::Integer) where {T,G} = tbl.data[tbl.fips_index[fips]]
Base.getindex(tbl::FIPSTable{G,T,1}, fips::Integer, var::AbstractString) where {T,G} = check_var(tbl, var) && return tbl.data[tbl.fips_index[fips]]
Base.getindex(tbl::FIPSTable{G,T,1}, var::AbstractString) where {T,G} = check_var(tbl, var) && return tbl.data
Base.getindex(tbl::FIPSTable{G,T,1}, k::Integer, fips::Integer, var::AbstractString) where {T,G} = check_var(tbl, var) && return tbl.data[tbl.fips_index[fips]]

Base.getindex(tbl::FIPSTable{G,T,2}, fips::Integer) where {T,G} = view(tbl.data, tbl.fips_index[fips], :)
Base.getindex(tbl::FIPSTable{G,T,2}, var::AbstractString) where {T,G} = view(tbl.data, :, tbl.var_index[var])
Base.getindex(tbl::FIPSTable{G,T,2}, fips::Integer, var::AbstractString) where {T,G} = tbl.data[tbl.fips_index[fips], tbl.var_index[var]]
function Base.getindex(tbl::FIPSTable{G,Float64,2}, k::Integer, fips::Integer, var::AbstractString) where G
    return haskey(tbl.fips_index, fips) ?
        tbl.data[tbl.fips_index[fips], tbl.var_index[var]] :
        NaN
end

Base.getindex(tbl::FIPSTable{G,T,3}, fips::Integer) where {T,G} = view(tbl.data, :, tbl.fips_index[fips], :)
Base.getindex(tbl::FIPSTable{G,T,3}, fips::Integer, var::AbstractString) where {T,G} = view(tbl.data, :, tbl.fips_index[fips], tbl.var_index[var])
Base.getindex(tbl::FIPSTable{G,T,3}, var::AbstractString) where {T,G} = view(tbl.data, :, :, tbl.var_index[var])
Base.getindex(tbl::FIPSTable{G,T,3}, k::Integer, fips::Integer, var::AbstractString) where {T,G} = tbl.data[k, tbl.fips_index[fips], tbl.var_index[var]]
# ============================================================================ #
geo_conversion(::Type{T}, ::Type{T}) where T = 1
geo_conversion(::Type{County}, ::Type{State}) = COUNTY2STATE
geo_conversion(::Type{Tract}, ::Type{State}) = TRACT2STATE
geo_conversion(::Type{BlockGroup}, ::Type{State}) = BG2STATE
# ---------------------------------------------------------------------------- #
geo_conversion(::Type{State}, ::Type{County}) = error("Cannot convert states to counties")
geo_conversion(::Type{Tract}, ::Type{County}) = TRACT2COUNTY
geo_conversion(::Type{BlockGroup}, ::Type{County}) = BG2COUNTY
# ---------------------------------------------------------------------------- #
geo_conversion(::Type{State}, ::Type{Tract}) = error("Cannot convert states to tracts")
geo_conversion(::Type{County}, ::Type{Tract}) = error("Cannot convert counties to tracts")
geo_conversion(::Type{BlockGroup}, ::Type{Tract}) = BG2TRACT
# ============================================================================ #
function aggregate!(out::Array{T,3}, k::Integer, tbl::FIPSTable{G,L,3},
    fips::AbstractVector{<:Integer}, do_avg::Bool) where {T<:Number, G<:AbstractGeo,L<:Number}

    for x in fips
        out[:,k,:] .+= tbl[x]
    end
    if do_avg
        out[:,k,:] ./= length(fips)
    end
    return out
end
# ---------------------------------------------------------------------------- #
function aggregate!(out::Array{T,2}, k::Integer, tbl::FIPSTable{G,L,2},
    fips::AbstractVector{<:Integer}, do_avg::Bool) where {T<:Number, G<:AbstractGeo,L<:Number}

    for x in fips
        out[k,:] .+= tbl[x]
    end
    if do_avg
        out[k,:] ./= length(fips)
    end
    return out
end
# ---------------------------------------------------------------------------- #
function aggregate!(out::Array{T,1}, k::Integer, tbl::FIPSTable{G,L,1},
    fips::AbstractVector{<:Integer}, do_avg::Bool) where {T<:Number, G<:AbstractGeo,L<:Number}

    for x in fips
        out[k] += tbl[x]
    end
    if do_avg
        out[k] /= length(fips)
    end
    return out
end
# ---------------------------------------------------------------------------- #
function aggregate(::Type{Go}, tbl::FIPSTable{Gi,T,N}, do_avg::Bool) where {Go,Gi,T,N}
    conv = geo_conversion(Gi, Go)
    conv == 1 && return tbl
    
    fips = all_fips(tbl)
    fips_conv = div.(fips, conv)
    unique_grps = sort!(unique(fips_conv))

    siz = collect(size(tbl.data))
    idx = N == 3 ? 2 : 1
    siz[idx] = length(unique_grps)

    data = zeros(Float64, siz...)

    for k in eachindex(unique_grps)
        idx = findall(isequal(unique_grps[k]), fips_conv)
        aggregate!(data, k, tbl, fips[idx], do_avg)
    end

    fips_idx = Dict{UInt64,Int}(x => k for (k,x) in enumerate(unique_grps))

    return FIPSTable{Go,Float64,N}(fips_idx, tbl.var_index, data)
end
# ============================================================================ #
aggregate_state(tbl::FIPSTable, avg::Bool) = aggregate(State, tbl, avg)
aggregate_county(tbl::FIPSTable, avg::Bool) = aggregate(County, tbl, avg)
aggregate_tract(tbl::FIPSTable, avg::Bool) = aggregate(Tract, tbl, avg)
# ============================================================================ #
function all_geo(::Type{T}, tbl::FIPSTable{G}) where {G,T<:AbstractGeo}
    return sort!(unique(div.(all_fips(tbl), geo_conversion(G, T))))
end
all_states(tbl::FIPSTable) = all_geo(State, tbl)
all_counties(tbl::FIPSTable) = all_geo(County, tbl)
all_tracts(tbl::FIPSTable) = all_geo(Tract, tbl)
# ============================================================================ #
end # module EpicastTables
