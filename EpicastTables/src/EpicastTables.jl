module EpicastTables

export FIPSTable, aggregate_state, aggregate_county, all_states, all_counties,
filter_columns
# ============================================================================ #
const TRACT2STATE = 10^9
const TRACT2COUNTY = 10^6
const COUNTY2STATE = 10^3
# ============================================================================ #
struct FIPSTable{T,N}
    fips_index::Dict{UInt64,Int}
    var_index::Dict{String,Int}
    data::Array{T,N}# time x fips x var | fips x var | fips x 1
end
# ---------------------------------------------------------------------------- #
FIPSTable(::Val{N}=Val(1)) where N = FIPSTable{N}(Dict{UInt64,Int}(), Dict{String,Int}(), Array{Float64,N}(undef, 0, 0))
# ---------------------------------------------------------------------------- #
function FIPSTable(data::Array{T,1}, fips_index::Dict{UInt64,Int},
    var::AbstractString) where T

    return FIPSTable{T,1}(fips_index, Dict{String,Int}(var => 1), data)
end
# ---------------------------------------------------------------------------- #
function FIPSTable(data::Array{T,N}, fips_index::Dict{UInt64,Int},
    var_index::Dict{String,Int}) where {T,N}

    return FIPSTable{T,N}(fips_index, var_index, data)
end
# ---------------------------------------------------------------------------- #
function FIPSTable(data::Array{T,N}, fips::AbstractVector{<:Integer},
    vars::AbstractVector{<:AbstractString}) where {T,N}

    fips_index = Dict{UInt64,Int}(UInt64(x) => k for (k,x) in enumerate(fips))
    var_index = Dict{String,Int}(x => k for (k,x) in enumerate(vars))

    return FIPSTable{T,N}(fips_index, var_index, data)
end
# ---------------------------------------------------------------------------- #
all_fips(tbl::FIPSTable) = collect(keys(tbl.fips_index))
has_fips(tbl::FIPSTable, fips::Integer) = haskey(tbl.fips_index, fips)
all_vars(tbl::FIPSTable) = collect(keys(tbl.var_index))
has_var(tbl::FIPSTable, var::AbstractString) = haskey(tbl.var_index, var)
data_array(tbl::FIPSTable) = tbl.data
# ---------------------------------------------------------------------------- #
Base.getindex(tbl::FIPSTable{T,1}, fips::Integer) where T = tbl.data[tbl.fips_index[fips]]
Base.getindex(tbl::FIPSTable{T,1}, fips::Integer, var::AbstractString) where T = tbl.data[tbl.fips_index[fips]]
function Base.getindex(tbl::FIPSTable{T,1}, var::AbstractString) where T
    @assert(haskey(tbl.var_index, var), "variable $(var) not found")
    return tbl.data
end

Base.getindex(tbl::FIPSTable{T,2}, fips::Integer) where T = view(tbl.data, tbl.fips_index[fips], :)
Base.getindex(tbl::FIPSTable{T,2}, var::AbstractString) where T = view(tbl.data, :, tbl.var_index[var])
Base.getindex(tbl::FIPSTable{T,2}, fips::Integer, var::AbstractString) where T = tbl.data[tbl.fips_index[fips], tbl.var_index[var]]

Base.getindex(tbl::FIPSTable{T,3}, fips::Integer) where T = view(tbl.data, :, tbl.fips_index[fips], :)
Base.getindex(tbl::FIPSTable{T,3}, fips::Integer, var::AbstractString) where T = view(tbl.data, :, tbl.fips_index[fips], tbl.var_index[var])
Base.getindex(tbl::FIPSTable{T,3}, var::AbstractString) where T = view(tbl.data, :, :, tbl.var_index[var])

filter_columns(f::Function, tbl::FIPSTable) = sort!(filter(f, all_vars(tbl)))
# ============================================================================ #
is_state(fips::Integer) = 0 < fips < 57
is_county(fips::Integer) = 3 < ceil(Int, log10(fips)) < 6

# NOTE: there is overlap b/t block groups and tract in terms of the number of
# digits, so care should be taken with is_tract()
is_tract(fips::Integer) = 9 < ceil(Int, log10(fips)) < 12
# ============================================================================ #
function state_conversion(tbl::FIPSTable)
    tmp = keys(tbl.fips_index)
    isempty(tmp) && return 0

    fips = first(tmp)
    is_state(fips) && return 1
    is_county(fips) && return COUNTY2STATE
    return TRACT2STATE
end
# ---------------------------------------------------------------------------- #
function county_conversion(tbl::FIPSTable)
    tmp = keys(tbl.fips_index)
    isempty(tmp) && return 0

    fips = first(tmp)

    is_state(fips) && return 0
    is_county(fips) && return 1
    return TRACT2COUNTY
end
# ============================================================================ #
function aggregate!(out::Array{T,3}, k::Integer, tbl::FIPSTable{3},
    fips::AbstractVector{<:Integer}) where T <: Number

    for x in fips
        out[:,k,:] .+= tbl[x]
    end
    return out
end
# ---------------------------------------------------------------------------- #
function aggregate!(out::Array{T,2}, k::Integer, tbl::FIPSTable{2},
    fips::AbstractVector{<:Integer}) where T <: Number

    for x in fips
        out[k,:] .+= tbl[x]
    end
    return out
end
# ---------------------------------------------------------------------------- #
function aggregate!(out::Array{T,1}, k::Integer, tbl::FIPSTable{1},
    fips::AbstractVector{<:Integer}) where T <: Number

    for x in fips
        out[k] .+= tbl[x]
    end
    return out
end
# ---------------------------------------------------------------------------- #
function aggregate(tbl::FIPSTable{T,N}, conv::Integer) where {T,N}
    conv == 1 && return tbl
    
    fips = all_fips(tbl)
    fips_conv = div.(fips, conv)
    unique_grps = sort!(unique(fips_conv))

    siz = collect(size(tbl.data))
    idx = N == 3 ? 2 : 1
    siz[idx] = length(unique_grps)

    data = zeros(T, siz...)

    for k in eachindex(unique_grps)
        idx = findall(isequal(unique_grps[k]), fips_conv)
        aggregate!(data, k, tbl, fips[idx])
    end

    fips_idx = Dict{UInt64,Int}(x => k for (k,x) in enumerate(unique_grps))

    return FIPSTable{N}(fips_idx, tbl.var_index, data)
end
# ============================================================================ #
function aggregate_state(tbl::FIPSTable)
    tmp = keys(tbl.fips_index)
    isempty(tmp) && return tbl

    return aggregate(tbl, state_conversion(tbl))
end
# ---------------------------------------------------------------------------- #
function aggregate_county(tbl::FIPSTable)
    tmp = keys(tbl.fips_index)
    isempty(tmp) && return tbl

    return aggregate(tbl, county_conversion(tbl))
end
# ============================================================================ #
function all_states(tbl::FIPSTable)
    return sort!(unique(div.(all_fips(tbl), state_conversion(tbl))))
end
function all_counties(tbl::FIPSTable)
    return sort!(unique(div.(all_fips(tbl), county_conversion(tbl))))
end
# ============================================================================ #
end # module EpicastTables
