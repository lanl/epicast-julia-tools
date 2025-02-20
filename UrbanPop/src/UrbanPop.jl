module UrbanPop

using Mmap, SHA

using Arrow, Tables, Printf, Query

using EpicastTables

const NAICS_DIGITS = 3
# ============================================================================ #
struct Agent
    pums_id::UInt64
    fips_code::UInt64

    household_id::UInt32
    person_id::UInt32
    
    household_income::UInt32
    
    person_commute_time::Int16

    person_naics::UInt16

    household_size::UInt8
    household_type::UInt8
    householder_age::UInt8
    household_kids::UInt8
    household_workers::UInt8
    household_nonworkers::UInt8
    household_adult_workers::UInt8
    household_adult_nonworkers::UInt8

    household_arrangement::UInt8
    household_tenure::UInt8
    household_vehicles::UInt8

    person_age::UInt8
    person_sex::UInt8
    person_race::UInt8
    person_ethnicity::UInt8
    person_commute_mode::UInt8

    person_commute_occupancy::UInt8
    person_ipr::UInt8
    person_employment::UInt8
    person_school_grade::UInt8
end
# ---------------------------------------------------------------------------- #
function Agent()
    return Agent(
        rand(UInt64), rand(UInt64), rand(UInt32), rand(UInt32), rand(UInt32),
        rand(Int16), rand(UInt16), rand(UInt8), rand(UInt8), rand(UInt8),
        rand(UInt8), rand(UInt8), rand(UInt8), rand(UInt8), rand(UInt8),
        rand(UInt8), rand(UInt8), rand(UInt8), rand(UInt8), rand(UInt8),
        rand(UInt8), rand(UInt8), rand(UInt8), rand(UInt8), rand(UInt8),
        rand(UInt8), rand(UInt8)
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

        get(args, :person_naics, a.person_naics),

        get(args, :household_size, a.household_size),
        get(args, :household_type, a.household_type),
        get(args, :householder_age, a.householder_age),
        get(args, :household_kids, a.household_kids),
        get(args, :household_workers, a.household_workers),
        get(args, :household_nonworkers, a.household_nonworkers),
        get(args, :household_adult_workers, a.household_adult_workers),
        get(args, :household_adult_nonworkers, a.household_adult_nonworkers),

        get(args, :household_arrangement, a.household_arrangement),
        get(args, :household_tenure, a.household_tenure),
        get(args, :household_vehicles, a.household_vehicles),

        get(args, :person_age, a.person_age),
        get(args, :person_sex, a.person_sex),
        get(args, :person_race, a.person_race),
        get(args, :person_ethnicity, a.person_ethnicity),
        get(args, :person_commute_mode, a.person_commute_mode),

        get(args, :person_commute_occupancy, a.person_commute_occupancy),
        get(args, :person_ipr, a.person_ipr),
        get(args, :person_employment, a.person_employment),
        get(args, :person_school_grade, a.person_school_grade)
    )
end
# ============================================================================ #
# for 12 digit FIPS block group codes
# @inline fips_tract2(x::Integer) = div(x, 10) - (div(x, 10^7) * 10^6)
@inline fips_tract(x::Integer) = div(x, 10)
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
struct_hash(::Type{T}) where T = sha2_256(join(map(string, fieldnames(T))))
agent_hash() = struct_hash(Agent)
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
function read_header(ifile::AbstractString)
    version = agent_hash()

    return open(ifile, "r") do io
        tmp = read(io, 32)
        @assert(tmp == version, "Version mismatch between current struct layout and version on disk!")
        nagent = Base.read(io, UInt64)
        nbytes = Base.read(io, UInt16)
        ntract = Base.read(io, UInt16)
        nhousehold = Base.read(io, UInt32) 
        return nagent, nbytes, ntract, nhousehold
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
    elseif x == "single_fam_detach"
        type = 0x01
    elseif x == "single_fam_attach"
        type = 0x02
    elseif x == "2_unit"
        type = 0x03
    elseif x == "3_4_unit"
        type = 0x04
    elseif x == "5_9_unit"
        type = 0x05
    elseif x == "10_19_unit"
        type = 0x06
    elseif x == "20_49_unit"
        type = 0x07
    elseif x == "GE50_unit"
        type = 0x08
    elseif x == "mob_home"
        type = 0x09
    elseif x == "other"
        type = 0x0a
    else
        error("Unknown household_type: \"$(x)\"")
    end
    return type
end
household_kids(x::AbstractString) = x == "yes" ? 0x01 : 0x00
household_income(x::Real) = x < 0 ? UInt32(0) : round(UInt32, x)
person_sex(x::AbstractString) = x == "male" ? 0x00 : 0x01
function person_race(x::AbstractString)
    if x == "white"
        race = 0x00
    elseif x == "blk_af_amer"
        race = 0x01
    elseif x == "asian"
        race = 0x02
    elseif x == "native_amer"
        race = 0x03
    elseif x == "pac_island"
        race = 0x04
    elseif x == "other"
        race = 0x05
    elseif x == "mult"
        race = 0x06
    else
        error("Unknown race: \"$(x)\"")
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
function household_arrangement(x::AbstractString)
    if x == ""
        out = 0x00
    elseif x == "married"
        out = 0x01
    elseif x == "male_no_spouse"
        out = 0x02
    elseif x == "female_no_spouse"
        out = 0x03
    elseif x == "alone"
        out = 0x04
    elseif x == "not_alone"
        out = 0x05
    else
        error("Invalid arrangement: \"$(x)\"")
    end
    return out
end
function household_tenure(x::AbstractString)
    if x == ""
        out = 0x00
    elseif x == "own"
        out = 0x01
    elseif x == "rent"
        out = 0x02
    elseif x == "other"
        out = 0x03
    else
        error("Invalid tenure type: \"$(x)\"")
    end
    return out
end
function household_vehicles(x::Union{AbstractString, Missing})
    if ismissing(x)
        out = 0x07
    elseif isempty(x)
        out = 0x00
    elseif startswith(x, "GE")
        out = 0x06
    else
        out = tryparse(UInt8, x)
        out == nothing && error("Invalid vehicle #: \"$(x)\"")
    end
    return out
end
function person_commute_occupancy(x::AbstractString)
    if x == ""
        out = 0x00
    elseif x == "drove_alone"
        out = 0x01
    elseif x == "carpooled"
        out = 0x02
    else
        error("Invalid vehicle occupancy: \"$(x)\"")
    end
    return out
end
function person_ipr(x::AbstractString)
    if x == ""
        out = 0x00
    elseif x == "L050"
        out = 0x01
    elseif x == "050_099"
        out = 0x02
    elseif x == "100_124"
        out = 0x03
    elseif x == "125_149"
        out = 0x04
    elseif x == "150_184"
        out = 0x05
    elseif x == "185_199"
        out = 0x06
    elseif x == "GE200"
        out = 0x07
    else
        error("Invalid ipr: \"$(x)\"")
    end
    return out
end
function person_employment(x::AbstractString)
    if x == ""
        out = 0x00
    elseif x == "not.in.force"
        out = 0x01
    elseif x == "unemp"
        out = 0x02
    elseif x == "employed"
        out = 0x03
    elseif x == "mil"
        out = 0x04
    else
        error("Invalid employment: \"$(x)\"")
    end
    return out
end
function person_school_grade(x::AbstractString)
    if x == ""
        out = 0x00
    elseif x == "preschl"
        out = 0x01
    elseif x == "kind"
        out = 0x02
    elseif x == "1st"
        out = 0x03
    elseif x == "2nd"
        out = 0x04
    elseif x == "3rd"
        out = 0x05
    elseif x == "4th"
        out = 0x06
    elseif x == "5th"
        out = 0x07
    elseif x == "6th"
        out = 0x08
    elseif x == "7th"
        out = 0x09
    elseif x == "8th"
        out = 0x0a # 10
    elseif x == "9th"
        out = 0x0b
    elseif x == "10th"
        out = 0x0c
    elseif x == "11th"
        out = 0x0d
    elseif x == "12th"
        out = 0x0e
    elseif x == "undergrad"
        out = 0x0f
    elseif x == "grad"
        out = 0x10 # 16
    else
        error("Invalid grate: \"$(x)\"")
    end
    return out
end
function person_naics(x::AbstractString, employment::AbstractString)

    # all agents w/ employment status of military have a naics of 0,
    # so use the closest one (national security)
    employment == "mil" && return UInt16(928)

    isempty(x) && return UInt16(0)

    # only keep maximum 3-digit specificity
    if length(x) > NAICS_DIGITS
        x = x[1:NAICS_DIGITS]
    end

    #= LETTERS W/IN NAICS CODES
    M = Multiple NAICS codes
    P = Part of a NAICS code - NAICS code split between two or more Census codes
    S = Not specified Industry in NAICS sector - Specific to Census codes only
    Z = Exception to NAICS code - Part of NAICS industry has own Census code

    Source: https://www.census.gov/topics/employment/industry-occupation/guidance/code-lists.html
    =#
    naics = tryparse(Int, replace(x, r"\D+"=>"0"))

    out = 0
    if naics != nothing
        if 0 < naics < 10^(NAICS_DIGITS-1)
            # make sure we have all NAICS_DIGITS digits, trailing 0's do NOT 
            # change the meaning of NAICS codes
            out = naics * 10^((NAICS_DIGITS-1) - floor(Int,log10(naics)))
        else
            out = naics
        end
        if !(10^(NAICS_DIGITS-1) <= out < 10^NAICS_DIGITS)
            @show(out, naics, 10^NAICS_DIGITS)
            error("Failed to properly detirmine naics code: \"$(x)\"")
        end
    else
        @warn("Failed to parse NAICS code \"$(x)\"")
    end
    return UInt16(out)
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

        person_naics(row[:pr_naics], row[:pr_emp_stat]),

        UInt8(row[:hh_size]),
        household_type(row[:hh_dwg]),
        UInt8(row[:hh_age]),        
        household_kids(row[:hh_has_kids]),        
        UInt8(row[:hh_nb_wrks]),
        UInt8(row[:hh_nb_non_wrks]),
        UInt8(row[:hh_nb_adult_wrks]),
        UInt8(row[:hh_nb_adult_non_wrks]),

        household_arrangement(row[:hh_living_arrangement]),
        household_tenure(row[:hh_tenure]),
        household_vehicles(row[:hh_vehicles]),

        UInt8(row[:pr_age]),
        person_sex(row[:pr_sex]),
        person_race(row[:pr_race]),
        person_ethnicity(row[:pr_hsplat]),
        person_commute_mode(row[:pr_travel]),

        person_commute_occupancy(row[:pr_veh_occ]),
        person_ipr(row[:pr_ipr]),
        person_employment(row[:pr_emp_stat]),
        person_school_grade(row[:pr_grade])
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
    # state-unique household ids, the actual values get overwritten in the loop
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
function convert_subdirs(idir::AbstractString, odir::AbstractString)
    dir_paths = filter(readdir(idir, join=true)) do x
        return match(r"\d\d_[A-Z]{2}", basename(x)) != nothing
    end

    ofiles = Vector{String}(undef, length(dir_paths))
    for k in eachindex(dir_paths)
        ofiles[k] = convert_feather_dir(dir_paths[k], odir)
        write_tract_file(ofiles[k])
    end
    cmd = "cd $(odir); tar cvf - *.bin | zstd -9 - > us-demographics.bin.zst; cd $(pwd())"
    return ofiles, cmd
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
struct Tract
    index::UInt64
    n_agent::UInt64
    fips_code::UInt64
end
# ============================================================================ #
Base.isless(a::Tract, b::Tract) = a.fips_code < b.fips_code
# ============================================================================ #
function write_tract_file(ifile::AbstractString)
    m = match(r".*/(\d{2})\.agents\.bin$", ifile)

    m == nothing && error("Invalid agents.bin file name \"$(ifile)\"")

    # state FIPS code
    state = m[1]

    nagent, _, ntract, nhousehold = read_header(ifile)

    data = Vector{Tract}(undef, ntract)

    raw = memmap(ifile)

    idx = 1
    last_tract = fips_tract(raw[1].fips_code)
    last_idx = 0

    for k in eachindex(raw)
        this_tract = fips_tract(raw[k].fips_code)
        if  this_tract != last_tract
            data[idx] = Tract(last_idx, k - last_idx - 1, last_tract)
            last_idx = k - 1
            idx += 1
            last_tract = this_tract
        end
    end

    # don't forget the last one... (don't subtract 1 as length(raw) is not an 
    # index, unlike <k> in the above loop)
    data[idx] = Tract(last_idx, length(raw) - last_idx, last_tract)

    ofile = joinpath(dirname(ifile), state * ".tracts.bin")

    open(ofile, "w") do io
        Base.write(io, struct_hash(Tract))
        Base.write(io, UInt64(length(data)))
        Base.write(io, UInt64(sizeof(Tract)))
        Base.write(io, data)
    end

    return ofile
end
# agent_db % tar cvf - *.bin | zstd -9 - > us-demographics.bin.zst
# ============================================================================ #
function read_tract_file(ifile::AbstractString)
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
# ============================================================================ #
function all_tract_data(idir::AbstractString)
    files = find_files(idir, r".*\.tracts\.bin$")
    all_tracts = Vector{Int64}(undef, 0)
    all_pop = Vector{Int64}(undef, 0)
    for file in files
        raw = read_tract_file(file)
        append!(all_tracts, Int.(getfield.(raw, :fips_code)))
        append!(all_pop, Int.(getfield.(raw, :n_agent)))
    end
    return all_tracts, all_pop
end
# ============================================================================ #
function naics_lookup(ifo::Dict{String,<:Any}, naics::Integer, ndigit=NAICS_DIGITS)
    naics_str = string(naics)
    if !haskey(ifo, naics_str)

        if naics_str[end] == '0'
            # check to see if there is matching code w/ trailing letters
            # (or a letter followed by digits or letters)
            pat = Regex("^" * naics_str[1:(ndigit-1)] * "[A-Z]\\w*\$")
            matches = filter!(collect(keys(ifo))) do key
                match(pat, key) != nothing
            end
            if isempty(matches)
                # punt, remove the trailing 0 and try again for a coarser match
                return naics_lookup(ifo, div(naics, 10), ndigit)
            else
                return map(matches) do key
                    key => ifo[key]
                end
            end
        else
            matches = filter!(collect(keys(ifo))) do key
                startswith(key, naics_str)
            end
            return map(matches) do key
                key => ifo[key]
            end
        end
    else
        return [naics_str => ifo[naics_str]]
    end
end
# ============================================================================ #
# WARNING: don't use this function unless you are SURE you know what you're
# doing... (i.e. you renamed a struct field but otherwise changed *NOTHING*)
function swap_hash(ifile::AbstractString)
    open(ifile, "a") do io
        current_hash = agent_hash()
        version = read(io, 32)
        if version != current_hash
            seekstart(io)
            Base.write(io, current_hash)
        end
    end
    return nothing
end
# ============================================================================ #
function tract_list(idir::AbstractString, ofile::AbstractString)
    total = 0
    open(ofile, "w") do io
        for file in find_files(idir, r".*\.tracts\.bin$")
            data = open(file, "r") do io
                seek(io, 32)
                n_tract = read(io, UInt64)
                seek(io, position(io) + sizeof(UInt64))
                data = Vector{Tract}(undef, n_tract)
                read!(io, data)
                return data
            end
            total += length(data)
            for tract in data
                println(io, tract.fips_code)
            end
        end
    end
    return total
end
# ============================================================================ #
function dist_str(x::AbstractVector{T}) where {T<:Real}
    denom = x[end] != 0 ? x[end] : T(1)
    return "[" * join(map(n -> @sprintf("%0.2f", n), x[1:end-1] ./ denom), ", ") * "]"
end
# ============================================================================ #
mutable struct TractMarginals
    age::Vector{Int}
    household_size::Vector{Int}
    n_agent::Int
    fips_code::Int
end
# ---------------------------------------------------------------------------- #
TractMarginals(n, fips) = TractMarginals(zeros(Int, 6), zeros(Int, 20), n, fips)
# ---------------------------------------------------------------------------- #
function Base.show(io::IO, t::TractMarginals)
    print(io, t.fips_code, ", ", t.n_agent, ", ", dist_str(t.age), ", ", dist_str(t.household_size))
end
# ---------------------------------------------------------------------------- #
function Base.:+(a::TractMarginals, b::TractMarginals)
    @assert(a.fips_code == b.fips_code, "FIPS code mismatch")
    a.n_agent += b.n_agent
    a.age .+= b.age
    a.household_size .+= b.household_size
    return a
end
# ---------------------------------------------------------------------------- #
function Base.:(==)(a::TractMarginals, b::TractMarginals)
    return (a.fips_code == b.fips_code) && (a.age == b.age) &&
        (a.household_size == b.household_size)
end
Base.:(!=)(a::TractMarginals, b::TractMarginals) = !(a == b)
# ---------------------------------------------------------------------------- #
function join_fill(x::Vector{<:Integer}, c::Char=' ')
    tmp = map(n->@sprintf("%5d", n), x)
    return join(tmp, c)
end
# ---------------------------------------------------------------------------- #
function Base.println(io::IO, t::TractMarginals, n::Integer=0)
    TRACT2COUNTY = 10^6
    t.n_agent < 1 && @warn("tract $(t.fips_code) has 0 agents?")
    county_fips = div(t.fips_code, TRACT2COUNTY)
    print(io, @sprintf("%5d", n), " ", @sprintf("%5d", t.n_agent), " ", 0, " ",
        county_fips, " ",
        @sprintf("%6d", t.fips_code - (county_fips * TRACT2COUNTY)), " ")
    
    hh_size = zeros(Int, 7)
    hh_size[1:6] .= t.household_size[1:6]
    hh_size[7] = sum(t.household_size[7:end-1])

    println(io, join_fill(t.age[1:end-1], ' '), " ", join_fill(hh_size, ' '))
end
# ============================================================================ #
function age_group(x::Agent)
    if x.person_age <= 5
        return 0
    elseif x.person_age <= 17
        return 1
    elseif x.person_age <= 29
        return 2
    elseif x.person_age <= 64
        return 3
    else # 65+
        return 4
    end
end
# ============================================================================ #
function household_size(x::Agent)
    if x.household_size > 19
        return 19
    else
        return Int(x.household_size)
    end
end
# ============================================================================ #
function tract_marginals(ifile::AbstractString, agents_by_size::Bool=false)

    raw = memmap(ifile)

    out = Dict{Int64,TractMarginals}()

    last_tract = 0
    n_agent = 0
    last_hh_id = -1

    for k in 1:length(raw)
        # tract fips code
        fips = Int(div(raw[k].fips_code, 10))
        if fips != last_tract
            out[fips] = TractMarginals(0, fips)
            if last_tract > 0
                out[last_tract].n_agent = n_agent
                n_agent = 0
            end
            last_tract = fips
        end
        tmp = out[fips]

        age_index = age_group(raw[k]) + 1
        tmp.age[age_index] += 1

        hh_id = Int(raw[k].household_id)
        
        hh_index = household_size(raw[k])
        if agents_by_size
            # add hosuehold size (count of agents per hh size, not hh count)        
            tmp.household_size[hh_index] += 1
        elseif last_hh_id != hh_id
            tmp.household_size[hh_index] += 1
        end
        last_hh_id = hh_id
        n_agent += 1
    end

    out[last_tract].n_agent = n_agent

    for (k,v) in out
        out[k].age[end] = sum(out[k].age)
        out[k].household_size[end] = sum(out[k].household_size)
    end

    return out
end
# ============================================================================ #
function write_tract_marginals(idir::AbstractString, odir::AbstractString)
    files = find_files(idir, r".*\.agents\.bin")
    n = 0
    for file in files
        m = match(r"(\d{1,2})\.agents\.bin", basename(file))
        m == nothing && error("failed to parse filename $(file)")
        ofile = joinpath(odir, m[1] * ".dat")
        tmp = sort!(collect(values(tract_marginals(file, false))), lt=(x,y)->x.fips_code<y.fips_code)
        
        open(ofile, "w") do io
            println(io, length(tmp))
            foreach((x,c)->println(io, x, c), tmp, n:(n+length(tmp)-1))
        end
        n += length(tmp)
    end

    return n
end
# ============================================================================ #
function locale_counts(ifile::AbstractString)
    tracts = div.(getfield.(memmap(ifile), :fips_code), 10)
    tmp = Set(tracts)
    n_tract = length(tmp)
    n_county = length(Set(div.(tracts, 10^6)))
    return n_county, n_tract
end
# ============================================================================ #
@inline agent_state(id::Integer) = id >> 58

@inline agent_id(id::Integer) = id & ~(UInt64(0b111111) << 58)

function agent_id_lt(a::Integer, b::Integer)
    return agent_state(a) == agent_state(b) ? agent_id(a) < agent_id(b) :
        agent_state(a) < agent_state(b)
end
# ============================================================================ #
function fetch_agent_demographics(agent_ids::AbstractVector{T},
    agent_dir::AbstractString) where T<:Integer

    ks = sortperm(agent_ids, lt=agent_id_lt)
    ids = agent_ids[ks]
    all_states = agent_state.(ids)
    states = unique(all_states)

    # out = Dict{T, Agent}()
    out = Vector{Agent}(undef, length(agent_ids))
    kf = 1
    for state in states
        raw = memmap(joinpath(agent_dir, lpad(state, 2, '0') * ".agents.bin"))
        kl = searchsortedlast(all_states, state)
        for k in kf:kl
            out[k] = raw[agent_id(ids[k]) + 1]
        end
        kf = kl + 1
    end

    return out[invperm(ks)]
end
# ============================================================================ #
@inline bgstr2geoid(x::AbstractString, conv::Integer) = div(parse(Int, x), conv)
# ============================================================================ #
function dt_nt_get_populations(idir::AbstractString; role::AbstractString="",
    state::Integer=35)
    return dt_nt_get_populations(BlockGroup, idir, role=role, state=state)
end
# ---------------------------------------------------------------------------- #
function dt_nt_get_populations(::Type{G}, idir::AbstractString; role::AbstractString="",
    state::Integer=35) where {G<:AbstractGeo}

    files = find_files(idir, r".*\.feather$")

    var_index = Dict{String,Int}("delta_pop"=>1)

    g_conv = EpicastTables.geo_conversion(BlockGroup, G)
    state_conv = EpicastTables.geo_conversion(G, State)

    # [nighttime (dest), daytime (orig)]
    out = [Dict{UInt64,Int}(), Dict{UInt64,Int}()]
    for file in files
        tbl = Arrow.Table(file)
        idx = isempty(role) ? collect(1:length(tbl[:role])) : findall(isequal(role), tbl[:role])

        # [nighttime (dest), daytime (orig)]
        for (k,field) in enumerate([:orig_geoid, :dest_geoid])
            all_geo = map(x -> bgstr2geoid(x, g_conv), tbl[field][idx])
            geo = filter!(unique(all_geo)) do g
                return div(g, state_conv) == state
            end

            for g in geo
                out[k][g] = get(out, g, 0) + count(isequal(g), all_geo)
            end
        end
    end

    # @show(out[1][350010018002], out[2][350010018002])

    all_fips = sort!(union(collect(keys(out[1])), collect(keys(out[2]))))

    fips_index = Dict{UInt64,Int}(v => k for (k,v) in enumerate(all_fips))

    data = zeros(Float64, length(fips_index), 1)

    for fips in all_fips
        # daytime - nighttime
        data[fips_index[fips],1] = (get(out[2], fips, 0) - get(out[1], fips, 0))# /
            #get(out[1], fips, 0)
    end

    return FIPSTable{G,Float64,2}(fips_index, var_index, data)
end
# ============================================================================ #
end # module UrbanPop