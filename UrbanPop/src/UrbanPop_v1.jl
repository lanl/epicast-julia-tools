module UrbanPop

# Copyright (C) 2025. Triad National Security, LLC. All rights reserved.
# This program was produced under U.S. Government contract 89233218CNA000001
# for Los Alamos National Laboratory (LANL), which is operated by Triad National
# Security, LLC for the U.S. Department of Energy/National Nuclear Security
# Administration. All rights in the program are reserved by Triad National
# Security, LLC, and the U.S. Department of Energy/National Nuclear Security
# Administration. The Government is granted for itself and others acting on its
# behalf a nonexclusive, paid-up, irrevocable worldwide license in this material
# to reproduce, prepare. derivative works, distribute copies to the public,
# perform publicly and display publicly, and to permit others to do so.

using Mmap, SHA

using Arrow, Tables, Printf, Query
# ============================================================================ #
struct Agent
    pums_id::UInt64
    fips_code::UInt64

    household_id::UInt32
    person_id::UInt32    
    
    household_income::UInt32
    
    person_commute_time::Int16

    household_size::UInt8
    household_type::UInt8
    householder_age::UInt8
    household_kids::UInt8
    household_workers::UInt8
    household_nonworkers::UInt8
    household_adult_workers::UInt8
    household_adult_nonworkers::UInt8

    person_age::UInt8
    person_sex::UInt8
    person_race::UInt8
    person_ethnicity::UInt8
    person_commute_mode::UInt8
end
# ---------------------------------------------------------------------------- #
function Agent()
    return Agent(
        rand(UInt64), rand(UInt64), rand(UInt32), rand(UInt32), rand(UInt32),
        rand(Int16), rand(UInt8), rand(UInt8), rand(UInt8), rand(UInt8),
        rand(UInt8), rand(UInt8), rand(UInt8), rand(UInt8), rand(UInt8),
        rand(UInt8), rand(UInt8), rand(UInt8), rand(UInt8)
    )

end
# ---------------------------------------------------------------------------- #
Base.show(io::IO, ::MIME"text/plain", agents::AbstractVector{Agent}) = Base.show(stdout, agents)
function Base.show(io::IO, agents::AbstractVector{Agent})    
    for agent in agents
        Base.show(IOContext(io, :compact => true), agent)
        print(io, "\n")
    end
    print(io, "\n")
end
# ---------------------------------------------------------------------------- #
function Base.show(io::IO, a::Agent)
    prefix, default_suffix = !get(io, :compact, false) ? ("  ", "\n") : ("",", ")
    n = 1
    fields = fieldnames(Agent)
    N = length(fields)
    for f in fields
        suffix = n < N ? default_suffix : ""
        print(io, prefix, f, ": ", getfield(a, f), suffix)
        n += 1
    end
end
# ============================================================================ #
function copy_agent(a::Agent; args...)
    return Agent(
        get(args, :pums_id, a.pums_id),
        get(args, :fips_code, a.fips_code),

        get(args, :household_id, a.household_id),
        get(args, :person_id, a.person_id),
        
        get(args, :household_income, a.household_income),
        
        get(args, :person_commute_time, a.person_commute_time),

        get(args, :household_size, a.household_size),
        get(args, :household_type, a.household_type),
        get(args, :householder_age, a.householder_age),
        get(args, :household_kids, a.household_kids),
        get(args, :household_workers, a.household_workers),
        get(args, :household_nonworkers, a.household_nonworkers),
        get(args, :household_adult_workers, a.household_adult_workers),
        get(args, :household_adult_nonworkers, a.household_adult_nonworkers),

        get(args, :person_age, a.person_age),
        get(args, :person_sex, a.person_sex),
        get(args, :person_race, a.person_race),
        get(args, :person_ethnicity, a.person_ethnicity),
        get(args, :person_commute_mode, a.person_commute_mode)
    )
end
# ============================================================================ #
# for 12 digit FIPS block group codes
@inline fips_tract(x::Integer) = div(x, 10) - (div(x, 10^7) * 10^6)
function count_tracts(data::AbstractVector{Agent})
    s = Set{UInt64}()
    @inbounds for k in eachindex(data)
        push!(s, fips_tract(data[k].fips_code))
    end
    return length(s)
end
function count_households(data::AbstractVector{Agent})
    s = Set{UInt64}()
    @inbounds for k in eachindex(data)
        push!(s, data[k].household_id)
    end
    return length(s)
end
# ============================================================================ #
agent_hash() = sha2_256(join(map(string, fieldnames(Agent))))
# ============================================================================ #
function write(ofile::AbstractString, data::AbstractVector{Agent})

    open(ofile, "w") do io
        n = length(data)
        nbytes = sizeof(Agent)

        # header:
        # 1) sha2_256 hash of field names joined as strings
        #   e3d7a281efe8c09fd88feba73914c4aea22f2b530d7b6d4cfee12c74fd5ba8a3
        # 2) number of entires as UInt64
        # 3) size of each entry/struct in bytes as a UInt64 (to maintain alignment)
        Base.write(io, agent_hash())
        Base.write(io, UInt64(n))
        Base.write(io, UInt16(nbytes))
        Base.write(io, UInt16(count_tracts(data)))
        Base.write(io, UInt32(count_households(data)))
        for k in 1:n
            unsafe_write(io, Ref(data[k]), nbytes)
        end
    end

end
# ============================================================================ #
function memmap(ifile::AbstractString)

    version = agent_hash()

    open(ifile, "r") do io
        tmp = read(io, 32)
        @assert(tmp == version, "Version mismatch between current struct layout and version on disk!")
        n = Base.read(io, UInt64)
        # skip the next 8 bytes (strust size, n-tract, n-household)
        seek(io, position(io) + 8)
        return Mmap.mmap(io, Vector{Agent}, (n,), grow=false)
    end
end
# ============================================================================ #
function is_contiguous(x::AbstractVector)
    # group indices of identical elements of <x>
    idxs = 1:length(x) |> @groupby(x[_]) |> collect

    # if an element appears in a contiguous region of <x>, then the indices
    # for that element are just a unit range min_index:max_index
    return all(idxs) do idx
        mn, mx = extrema(idx)
        return idx == mn:mx
    end
end
# ============================================================================ #
function household_type(x::AbstractString)
    if x == "" # group quarters
        type = 0x00
    elseif x == "dwg_single_fam_detach"
        type = 0x01
    elseif x == "dwg_single_fam_attach"
        type = 0x02
    elseif x == "dwg_2_unit"
        type = 0x03
    elseif x == "dwg_3_4_unit"
        type = 0x04
    elseif x == "dwg_5_9_unit"
        type = 0x05
    elseif x == "dwg_10_19_unit"
        type = 0x06
    elseif x == "dwg_20_49_unit"
        type = 0x07
    elseif x == "dwg_GE50_unit"
        type = 0x08
    elseif x == "dwg_mob_home"
        type = 0x09
    elseif x == "dwg_other"
        type = 0x0a
    else
        error("Unknown household_type: \"$(x)\"")
    end
    return type
end
household_kids(x::AbstractString) = x == "yes" ? 0x01 : 0x00
household_income(x::Real) = x < 0 ? UInt32(0) : UInt32(x)
person_sex(x::AbstractString) = x == "male" ? 0x00 : 0x01
function person_race(x::AbstractString)
    if x == "race_white"
        race = 0x00
    elseif x == "race_blk_af_amer"
        race = 0x01
    elseif x == "race_asian"
        race = 0x02
    elseif x == "race_native_amer"
        race = 0x03
    elseif x == "race_pac_island"
        race = 0x04
    elseif x == "race_other"
        race = 0x05
    elseif x == "race_mult"
        race = 0x06
    end
    return race
end
person_ethnicity(x::AbstractString) = x == "no" ? 0x00 : 0x01
function person_commute_mode(x::AbstractString)
    if x == "" # non-worker
        mode = 0x00
    elseif x == "car_truck_van"
        mode = 0x01
    elseif x == "public_transportation"
        mode = 0x02
    elseif x == "bicycle"
        mode = 0x03
    elseif x == "walked"
        mode = 0x04
    elseif x == "motorcycle"
        mode = 0x05
    elseif x == "taxicab"
        mode = 0x06
    elseif x == "other"
        mode = 0x07
    elseif x == "wfh"
        mode = 0x08
    end
    return mode
end
function parse_pumsid(x::AbstractString)
    pums, hh = split(x, '-')
    # pums id format: yyyyXXddddd...
    # where yyyy is 4-digit year, XX is 00 (<2018) or HU/GQ (>=2018), and ddd... 
    # is the "person serial number"
    return parse(UInt64, pums[vcat(1:4, 7:length(pums))]), parse(UInt64, hh)
end
# ============================================================================ #
function convert_row(row::Tables.ColumnsRow, hh_count::UInt32, 
    last_hh_id::UInt64)

    # we are parsing h_id (not p_id) as h_id ends with '-' and a household id,
    # whereas p_id ends with '-' and household id + person id
    pums_id, hh_id = parse_pumsid(row[:h_id])
    hh_count += hh_id != last_hh_id ? UInt32(1) : UInt32(0)

    return Agent(
        pums_id,
        parse(UInt64, row[:geoid]),

        UInt32(hh_count),
        UInt32(0), # agent id get set at the very end of convert_feather_dir()

        household_income(row[:hh_income]),

        Int16(row[:pr_commute]),

        UInt8(row[:hh_size]),
        household_type(row[:hh_dwg]),
        UInt8(row[:hh_age]),        
        household_kids(row[:hh_has_kids]),        
        UInt8(row[:hh_nb_wrks]),
        UInt8(row[:hh_nb_non_wrks]),
        UInt8(row[:hh_nb_adult_wrks]),
        UInt8(row[:hh_nb_adult_non_wrks]),
        UInt8(row[:pr_age]),
        person_sex(row[:pr_sex]),
        person_race(row[:pr_race]),
        person_ethnicity(row[:pr_hsplat]),
        person_commute_mode(row[:pr_travel])
    ), hh_count, hh_id
end
# ============================================================================ #
function convert_feather(ifile::AbstractString, hh_count::UInt32, 
    last_hh_id::UInt64)

    tbl = Arrow.Table(ifile)

    # make sure all agents from a given household appear together (i.e. in a 
    # contiguous slice)
    @assert(is_contiguous(tbl[:h_id]), "Households are not grouped in table \"$(ifile)\"")

    nagent = Tables.rowcount(tbl)

    data = Vector{Agent}(undef, nagent)
    for (k, row) in enumerate(Tables.rows(tbl))
        data[k], hh_count, last_hh_id = convert_row(row, hh_count, last_hh_id)
    end

    return data, hh_count, last_hh_id
end
# ============================================================================ #
function convert_feather_dir(idir::AbstractString, odir::AbstractString)
    files = find_files(idir, r".*\.feather$")
    fips = map(files) do x
        m = match(r".*syp_(\d+)\.feather$", x)
        @assert(m != nothing, "Failed to parse filename \"$(x)\"")

        # FIPS code is SSCCC (s = state, c = county)
        return parse(Int, m[1])
    end

    state = floor.(Int, fips ./ 1e3)
    @assert(all(isequal(state[1]), state), "not all files have the same FIPS state code!")

    # sort files by fips ascending
    ks = sortperm(fips)
    files .= files[ks]

    ofile = joinpath(odir, @sprintf("%02d.agents.bin", state[1]))

    # we keep a count of the total number households for assigning
    # state-unique household ids, the actualy value get overwritten in the loop
    # at the end of this function, but setting up unique ids in the 
    # convert_feather() loops allows them to be easily changed at the end
    hh_count = UInt32(0)
    
    # place holder so we know when a "new" household should be counted
    last_hh_id = UInt64(0)

    agents = Vector{Agent}(undef, 0)
    N = length(files)
    for (k, file) in enumerate(files)
        data, hh_count, last_hh_id = convert_feather(file, hh_count, last_hh_id)
        append!(agents, data)
        println("[DONE]: $(k)/$(N) \"$(basename(file))\"")
    end

    @info("Total agents in state $(state[1]): $(length(agents))")

    # sort agents fips code then by household id
    sort!(agents, by = x -> (x.fips_code, x.household_id));

    # yes, this is kind of henious, but our lives will be *MUCH* easier if
    # the database is sorted by fips_code, household_id, & person_id *and* it's
    # better to leave the Agent struct as immutable (thus the copy...)
    hhid = -1
    hhid_last = -1
    for k in eachindex(agents)
        person_id = k - 1
        hhid += agents[k].household_id != hhid_last ? 1 : 0
        hhid_last = agents[k].household_id

        agents[k] = copy_agent(agents[k],
            person_id = person_id, 
            household_id = hhid
        )
    end

    write(ofile, agents)

    return ofile
end
# ============================================================================ #
function find_files(idir::AbstractString, re::Regex=r".*")
    out = String[]
    for name in readdir(idir)
        if match(re, name) != nothing
            push!(out, joinpath(idir, name))
        end
    end
    return out
end
# ============================================================================ #
end # module UrbanPop