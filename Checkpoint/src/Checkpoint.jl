module Checkpoint

using SHA, Mmap

const IntTuple{N} = NTuple{N,Int32}
const SPaSMVector = NTuple{3,Float64}
const SPaSMIVector = NTuple{3,Int32}
# ============================================================================ #
function init_zero(x::T) where T
    Base.memset(pointer_from_objref(x), 0, sizeof(T))
    return x
end
# ---------------------------------------------------------------------------- #
struct CellData
    tract::UInt32
    community::UInt32
    n_workgroup::UInt32
    n_family::UInt32

    nc_0::IntTuple{6}
    nc_0_race::IntTuple{8}
    nc_0_ethnicity::IntTuple{3}
    nc_0_HH::IntTuple{20}
    nc_tot::IntTuple{6}
    nc_tot_race::IntTuple{8}
    nc_tot_ethnicity::IntTuple{3}
    nc_tot_HH::IntTuple{20}
    nc_ill::IntTuple{6}
    nc_ill_race::IntTuple{8}
    nc_ill_ethnicity::IntTuple{3}

    nc_tot2::IntTuple{6}
    nc_tot2_race::IntTuple{8}
    nc_tot2_ethnicity::IntTuple{3}
    nc_ill2::IntTuple{6}
    nc_ill2_race::IntTuple{8}
    nc_ill2_ethnicity::IntTuple{3}
    nc_hosp::IntTuple{6}
    nc_icu::IntTuple{6}
    nc_vent::IntTuple{6}
    nc_dead::IntTuple{6}

    inf_source::NTuple{12,UInt16}

    work_scale::Float64
    social_scale::Float64
    travel_scale::Float64
    bus_open::Float64

    schools_open::UInt8

end
# ---------------------------------------------------------------------------- #
struct Particle
    type::Int32
    tag::Int32

    r::SPaSMVector

    home::SPaSMIVector
    work::SPaSMIVector
    travel::SPaSMIVector

    trip_timer::UInt8
    nbor_all::UInt8
    nborhood::UInt8
    workgroup::UInt8

    small_workgroup::Int8

    naics_code::UInt16
    family::UInt16
    status::UInt16

    hh_size::UInt8
    vacc_tier::UInt8
    vacc_timer::UInt8

    treatment_timer::Int8
    school::Int8
    race::Int8
    ethnicity::Int8
    employment::Int8

    strain2::UInt8

    prob::NTuple{2,Float64}

    p_family::Float64
    p_school::Float64
    p_work::Float64
    p_nc::Float64
    p_hood::Float64
    p_bar::Float64
    p_school_mix::Float64
    p_bc::Float64

end
# ---------------------------------------------------------------------------- #
Base.zero(::Type{NTuple{N,T}}) where {N,T} = tuple(zeros(T, N)...)

init_struct(::Type{T}) where T = T(zero.(fieldtypes(T))...)
# ---------------------------------------------------------------------------- #
mutable struct Community
    cell_data::CellData
    particles::Vector{Particle}
end
# ---------------------------------------------------------------------------- #
function Community(n_agent::Integer)
    return Community(
        init_struct(CellData),
        map(x -> init_struct(Particle), 1:n_agent)
    )
end
# ---------------------------------------------------------------------------- #
function Base.show(io::IO, x::CellData)
    println(io, "tract: ", Int(x.tract), ", comm: ", Int(x.community),
        ", nwg: ", Int(x.n_workgroup), ", nfam: ", Int(x.n_family))
    println(io, "  age: ", x.nc_0, ", race: ", x.nc_0_race,
        ", eth: ", x.nc_0_ethnicity)
    println(io, "  hh: ", x.nc_0_HH)
end
# ---------------------------------------------------------------------------- #
function read_checkpoint(ifile::AbstractString)
    return open(ifile, "r") do io
        n_community = read(io, UInt64)
        cd_size = read(io, UInt64)
        pt_size = read(io, UInt64)

        @assert(cd_size == sizeof(CellData))
        @assert(pt_size == sizeof(Particle))

        data = Vector{Community}(undef, n_community)
        for k = 1:n_community
            n_agent = read(io, UInt64)
            data[k] = Community(n_agent)
            ref = Ref(data[k].cell_data)
            read!(io, ref)
            data[k].cell_data = ref[]
            read!(io, data[k].particles)
        end

        return data
    end
end
# ============================================================================ #
struct_hash(::Type{T}) where T = sha2_256(join(map(string, fieldnames(T))))
# ---------------------------------------------------------------------------- #
struct Tract
    index::UInt64
    n_agent::UInt64
    fips_code::UInt64
end
# ---------------------------------------------------------------------------- #
function read_tracts(ifile::AbstractString)
    hash = struct_hash(Tract)
    return open(ifile, "r") do io
        current_hash = read(io, 32)
        @assert(hash == current_hash, "tract struct layout does not match!")
        n_tract = read(io, UInt64)
        tract_size = read(io, UInt64)
        @assert(tract_size == sizeof(Tract), "tract struct size mismatch")
        
        return mmap(io, Vector{Tract}, (n_tract,), grow=false, shared=false)
    end
end
# ---------------------------------------------------------------------------- #
function read_workerflow(ifile::AbstractString)
    from, to, n = open(ifile, "r") do io
        n_tract = read(io, UInt64)
        n_entry = read(io, UInt64)
        from = Vector{UInt64}(undef, n_entry)
        to = Vector{UInt64}(undef, n_entry)
        n = Vector{UInt32}(undef, n_entry)
        for k = 1:n_entry
            from[k] = read(io, UInt64)
            to[k] = read(io, UInt64)
            n[k] = read(io, UInt32)
        end
        return from, to, n
    end

    return Dict{Tuple{Int,Int}, Int}((Int(f),Int(t)) => n for (f,t,n) in
        zip(from,to,n))
end
# ---------------------------------------------------------------------------- #
function get_column(d::Dict{Tuple{Int,Int}, Int}, idx::Integer)
    col = Dict{Int,Int}()
    for (k,v) in d
        if k[2] == idx
            col[k[1]] = d[k]
        end
    end
    return col
end
# ---------------------------------------------------------------------------- #
function get_row(d::Dict{Tuple{Int,Int}, Int}, idx::Integer)
    row = Dict{Int,Int}()
    for (k,v) in d
        if k[1] == idx
            row[k[2]] = d[k]
        end
    end
    return row
end
# ---------------------------------------------------------------------------- #
vec_sum(d::Dict{Int,Int}, rm::Integer) = sum(values(d)) - d[rm]
# ============================================================================ #
function community_counts(tr_file::AbstractString, wf_file::AbstractString)
    tr = read_tracts(tr_file)
    wf = read_workerflow(wf_file)

    n_res = zeros(UInt8, length(tr))
    n_wrk = zeros(UInt8, length(tr))

    for k in eachindex(tr)
        fips = Int(tr[k].fips_code)
        col = get_column(wf, fips)
        n_res[k] = max(1, round(UInt8, tr[k].n_agent / 2000))
        n_wrk[k] = isempty(col) ? 0x00 : round(UInt8, vec_sum(col, fips) / 2000)
    end

    return tr, n_res, n_wrk
end
# ---------------------------------------------------------------------------- #
function stats(c::CellData)
    if c.n_family == 0
        return (0,1,0)
    else
        return (1, 0, c.nc_0[end])
    end
end
# ---------------------------------------------------------------------------- #
function checkpoint_community_counts(ck_file::AbstractString)
    ck = read_checkpoint(ck_file)
    out = Dict{Int,Tuple{Int,Int,Int}}()
    for k in eachindex(ck)
        if ck[k].cell_data.tract != 0xffffffff
            tract = Int(ck[k].cell_data.tract)
            cs = stats(ck[k].cell_data)
            if !haskey(out, tract)
                out[tract] = cs
            else
                out[tract] = out[tract] .+ cs
            end
        end
    end
    return out
end
# ---------------------------------------------------------------------------- #
function validate_checkpoint(ck_file::AbstractString, tr_file::AbstractString,
    wf_file::AbstractString)

    tr, n_res, n_wrk = community_counts(tr_file, wf_file)

    ck = checkpoint_community_counts(ck_file)

    @assert(length(ck) == length(tr), "Mismatch in # of tracts")

    n = 0

    for k in keys(ck)
        idx = k + 1
        tmp = (n_res[idx], n_wrk[idx], Int(tr[idx].n_agent))
        if ck[k] != tmp
            @warn("FAIL - ", idx, " ", Int(tr[idx].fips_code), " - ck = ", ck[k], " db = ", tmp)
            n += 1
        end
    end

    if n == 0
        @info("SUCCESS - it appears that the test passed...")
    end

    return n
end
# ============================================================================ #
end
#=
TODO:
    1) compare age and household distributions from DB vs CK (at the tract level?)
    2) could also compate race, ethnicity, etc. distributions for v2
=#