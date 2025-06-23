module LDTM

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

using CSV, Arrow
using UrbanPop
# ============================================================================ #
function get_item(line::AbstractString, idx::Integer, ::Type{T}=Int) where T <: Real
    n = 0
    ks = 1
    ke = 1
    while n < idx && ke < length(line)
        if line[ke] == ','
            n += 1
            ks = n < idx ? ke + 1 : ks
        end
        ke += n < idx ? 1 : -1
    end    
    return parse(T, strip(line[ks:ke]))
end
get_float_item(line::AbstractString, idx::Integer) = Int(get_item(line, idx, Float64))
# ============================================================================ #
function get_tract_state(line::AbstractString, idx::Integer)
    tract = get_item(line, idx, Int)
    return div(tract, 10^9)
end
# ============================================================================ #
function split_csv(ifile::AbstractString, oname::AbstractString,
    field::AbstractString="trOState", f_get_item::Function=get_item)

    odir = joinpath(dirname(ifile), oname)
    
    !isdir(odir) && mkpath(odir)

    prefix, _ = splitext(basename(ifile))

    prefix = replace(prefix, r"_part\d"=>"")

    open(ifile, "r") do io
        hdr_str = readline(io)
        names = split(hdr_str, ',')

        state_idx = findfirst(isequal(field), names)
        state_idx == nothing && error("Failed to find field \"$(field)\"")

        last_state = 0
        os = Base.DevNull()

        for line in eachline(io)
            state = f_get_item(line, state_idx)

            if state != last_state
                close(os)
                if last_state > 0
                    println("DONE: state $(last_state)")
                end

                ofile = joinpath(odir, prefix * "_" * lpad(state, 2, '0') * ".csv")
                if isfile(ofile)
                    os = open(ofile, append=true)
                else
                    os = open(ofile, "w")
                    println(os, hdr_str)
                end
                last_state = state
            end

            println(os, line)
        end

        close(os)
    end

    return nothing
end
# ============================================================================ #
struct Trip
    hhid::UInt32 # epicast-urbanpop household id 0 - N-1 where N is the number
                 # of households in the state in which this household resides

    destination::UInt16 # state + county FIPS

    # ENUMS
    trip_id::UInt8 # trip number for this household in a Jan-Dec year
    trip_month::UInt8 # month in which this trip occured (1 - 12)
    trip_nagent::UInt8 # \# of agents from the household that go on this trip
        # 0 => 1 agent
        # 1 => 2 agents
        # 2 => 3 agents
        # 3 => 4+ agents
    trip_purpose::UInt8 # purpose of the trip
        # 0 => personal buisness
        # 1 => visiting friends / relatives
        # 2 => leisure
        # 3 => commute
        # 4 => employment buisness
    trip_duration::UInt8 # duration of the trip
        # 0 => day trip (0 nights)
        # 1 => 2-3 days (1-2 nights)
        # 2 => 4-7 days (3-6 nights)
        # 3 => 8+ days (7+ nights)
end
# ---------------------------------------------------------------------------- #
function Trip(row::CSV.Row)
    return Trip(

    )
end
# ============================================================================ #
function csv2dict(f::Function, csvfile::AbstractString, ::Type{K}, ::Type{V}) where {K,V}
    return Dict{K,V}( f(row) for row in CSV.Rows(csvfile) )
end
# ============================================================================ #
# UrbanPop h_id => BG FIPS code for a single county
function h_id_map(upfile::AbstractString)
    tbl = Arrow.Table(upfile)
    return Dict{String,Int}(
        k => parse(Int, v) for (k,v) in zip(tbl[:h_id], tbl[:geoid])
    )
end
# ============================================================================ #
# UrbanPop h_id => BG FIPS code for an entire state
function h_id_map_all(idir::AbstractString)
    files = UrbanPop.find_files(idir, r"\.feather$")
    out = Dict{String,Int}()
    for file in files
        merge!(out, h_id_map(file))
    end
    return out
end
# ============================================================================ #
# LDTM zone => county FIPS
function zone_map(ifile::AbstractString)
    return csv2dict(ifile, Int, Int) do row
        return parse(Int, row.id) => (parse(Int, row.STATEFP) * 1000 + parse(Int, row.COUNTYFP))
    end
end
# ============================================================================ #
# LDTM hhid => UrbanPop h_id
function hhid_map(popfile::AbstractString)
    return csv2dict(popfile, Int, String) do row
        return parse(Int, row.hhid) => String(row.h_id)
    end
end
# ============================================================================ #
function create_state_db(ldtm_dir::AbstractString, urbanpop_dir::AbstractString,
    ldtm_zone_map::AbstractString)



end
# ============================================================================ #
end # module LDTM
