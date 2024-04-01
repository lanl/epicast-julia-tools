module Checkpoint

using Workerflow, UrbanPop

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
function read_checkpoint_file(ifile::AbstractString)
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
const N_AGENT_PER_COMMUNITY = 2000
const CELL_UNUSED = 0xffffffff # marks a community/cell as unused
# ============================================================================ #
function community_counts(tr_file::AbstractString, wf_file::AbstractString)
    tr = UrbanPop.read_tract_file(tr_file)
    wf = Workerflow.read_workerflow_file(wf_file)

    n_res = zeros(UInt8, length(tr))
    n_wrk = zeros(UInt8, length(tr))

    for k in eachindex(tr)
        fips = Int(tr[k].fips_code)
        col = Workerflow.get_column(wf, fips)
        n_res[k] = max(1, round(UInt8, tr[k].n_agent / N_AGENT_PER_COMMUNITY))
        
        n_wrk[k] = isempty(col) ? 0x00 :
            round(UInt8, Workerflow.vec_sum(col, fips) / N_AGENT_PER_COMMUNITY)
    end

    return tr, n_res, n_wrk
end
# ---------------------------------------------------------------------------- #
function stats(c::CellData, t::UrbanPop.Tract)
    if c.n_family == 0
        # worker only community
        return 0x00, 0x01, UrbanPop.TractMarginals(t.n_agent, t.fips_code)
    else
        return 0x01, 0x00, UrbanPop.TractMarginals(Int[c.nc_0...], Int[c.nc_0_HH...],
            t.n_agent, t.fips_code)
    end
end
# ---------------------------------------------------------------------------- #
function Base.:+(a::UrbanPop.TractMarginals, b::UrbanPop.TractMarginals)
    @assert(a.fips_code == b.fips_code, "FIPS code mismatch")
    a.age .+= b.age
    a.household_size .+= b.household_size
    return a
end
# ---------------------------------------------------------------------------- #
function Base.:(==)(a::UrbanPop.TractMarginals, b::UrbanPop.TractMarginals)
    return (a.fips_code == b.fips_code) && (a.age == b.age) &&
        (a.household_size == b.household_size)
end
Base.:(!=)(a::UrbanPop.TractMarginals, b::UrbanPop.TractMarginals) = !(a == b)
# ---------------------------------------------------------------------------- #
function checkpoint_community_counts(ck_file::AbstractString, tr::Vector{UrbanPop.Tract})
    ck = read_checkpoint_file(ck_file)

    out = Dict{Int,Tuple{Int,Int,UrbanPop.TractMarginals}}()

    for k in eachindex(ck)
        if ck[k].cell_data.tract != CELL_UNUSED
            idx = Int(ck[k].cell_data.tract)
            @assert(0 <= idx < length(tr), "Invalid tract idx: $(idx)")
            n_res, n_wrk, mrgn = stats(ck[k].cell_data, tr[idx + 1])
            tract = tr[idx + 1].fips_code
            if !haskey(out, tract)
                out[tract] = (n_res, n_wrk, mrgn)
            else
                out[tract] = (out[tract][1] + n_res, out[tract][2] + n_wrk,
                    out[tract][3] + mrgn)
            end
        end
    end
    return out
end
# ---------------------------------------------------------------------------- #
function validate_checkpoint(ck_file::AbstractString, tr_file::AbstractString,
    wf_file::AbstractString, ag_file::AbstractString)

    tr, n_res, n_wrk = community_counts(tr_file, wf_file)

    mrgn = UrbanPop.tract_marginals(ag_file)

    ck = checkpoint_community_counts(ck_file, tr)

    @assert(length(ck) == length(mrgn), "Mismatch in # of tracts")
    @assert(keys(ck) == keys(mrgn), "Mismatch in tract FIPS codes")

    fips = getfield.(tr, :fips_code)

    n = 0

    for k in keys(ck)
        idx = findfirst(isequal(k), fips)

        tmp = (n_res[idx], n_wrk[idx], mrgn[k])
        
        if ck[k] != tmp
            println("FAIL @ ", k)
            println("    ck = ", ck[k][3], "\n    db = ", mrgn[k])
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