module EpicastTypes

export AbstractCellData, AbstractParticle, CellData, CellDataOrig, Particle, 
    ParticleOrig, Community, init_struct

const IntTuple{N} = NTuple{N,Int32}
const SPaSMVector = NTuple{3,Float64}
const SPaSMIVector = NTuple{3,Int32}
# ============================================================================ #
abstract type AbstractCellData end
abstract type AbstractParticle end
# ============================================================================ #
struct CellData <: AbstractCellData
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
struct CellDataOrig <: AbstractCellData
    tract::UInt32
    community::UInt32
    n_workgroup::UInt32
    n_family::UInt32

    nc_0::IntTuple{6}
    nc_tot::IntTuple{6}
    nc_ill::IntTuple{6}

    nc_tot2::IntTuple{6}
    nc_ill2::IntTuple{6}

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
# ============================================================================ #
struct Particle <: AbstractParticle
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
struct ParticleOrig <: AbstractParticle
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

    naics_code::UInt16
    family::UInt16
    status::UInt16

    hh_size::UInt8
    vacc_tier::UInt8
    vacc_timer::UInt8

    treatment_timer::Int8
    school::Int8

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
# ============================================================================ #
Base.zero(::Type{NTuple{N,T}}) where {N,T} = tuple(zeros(T, N)...)
init_struct(::Type{T}) where T = T(zero.(fieldtypes(T))...)
# ============================================================================ #
mutable struct Community{C<:AbstractCellData,P<:AbstractParticle}
    cell_data::C
    particles::Vector{P}
end
# ---------------------------------------------------------------------------- #
function Community(::Type{C}, ::Type{P}, n_agent::Integer) where {C<:AbstractCellData, P<:AbstractParticle}
    return Community(
        init_struct(C),
        map(x -> init_struct(P), 1:n_agent)
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
function Base.show(io::IO, x::CellDataOrig)
    println(io, "tract: ", Int(x.tract), ", comm: ", Int(x.community),
        ", nwg: ", Int(x.n_workgroup), ", nfam: ", Int(x.n_family))
    println(io, "  age: ", x.nc_0)
end
# ---------------------------------------------------------------------------- #
function Base.:(==)(a::CellData, b::CellData)
    return a.tract == b.tract && a.community == b.community &&
        a.nc_0 == b.nc_0 && a.nc_0_HH == b.nc_0_HH
end
function Base.:(==)(a::Community, b::Community)
    return a.cell_data == b.cell_data && a.particles == b.particles
end
# ============================================================================ #
end
