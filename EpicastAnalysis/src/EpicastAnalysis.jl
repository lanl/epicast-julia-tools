module EpicastAnalysis

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

using Epicast, UrbanPop
export DictMatrix, row_labels, column_labels
# ============================================================================ #
struct DictMatrix{I1,I2,T}
    index1::Dict{I1,Int}
    index2::Dict{I2,Int}
    data::Matrix{T}
end
# DictMatrix(i1::Dict{I1,Int}, i2::Dict{I2,Int}, d::Matrix{T}) where {I1,I2,T} = DictMatrix{I1,I2,T}(i1,i2,d)
data(d::DictMatrix) = d.data
Base.size(d::DictMatrix) = Base.size(d.data)
Base.size(d::DictMatrix, k::Integer) = Base.size(d.data, k)
Base.getindex(d::DictMatrix, r, c) = d.data[d.index1[r], d.index2[c]]
Base.getindex(d::DictMatrix, r, ::Colon) = d.data[d.index1[r], 1:size(d.data,2)]
Base.getindex(d::DictMatrix, ::Colon, c) = d.data[1:size(d.data,1), d.index2[c]]
function dict_labels(idx::Dict{T,Int}) where T
    return collect(keys(idx))[sortperm(collect(values(idx)))]
end
row_labels(d::DictMatrix) = dict_labels(d.index1)
column_labels(d::DictMatrix) = dict_labels(d.index2)
# ============================================================================ #
function count_unique(x::AbstractVector{T}) where T
    return [u => sum(isequal(u), x) for u in sort!(unique(x))]
end
# ============================================================================ #
function state_demographics(agent_dir::AbstractString, state::Integer)
    return UrbanPop.memmap(joinpath(agent_dir, lpad(state, 2, '0') * ".agents.bin"))
end
# ============================================================================ #
function get_demographic(agent_dir::AbstractString, states::AbstractVector{<:Integer}, demog::Symbol)
    raw = state_demographics(agent_dir, states[1])
    out = Dict(count_unique(getfield.(raw, demog)))
    for k in 2:length(states)
        raw = state_demographics(agent_dir, states[k])
        tmp = Dict(count_unique(getfield.(raw, demog)))
        for s in union(keys(out), keys(tmp))
            out[s] = get(out, s, 0) + get(tmp, s, 0)
        end
    end
    return sort!([k => v for (k,v) in out], lt=(x,y)->x.first<y.first)
end
# ============================================================================ #
function aggregate_demographic_by(raw::Vector{UrbanPop.Agent}, demog::Symbol, conv::Integer=10)

    u_tract = sort!(unique(div.(getfield.(raw, :fips_code), conv)))
    u_demog = sort!(unique(getfield.(raw, demog)))
    tract_idx = Dict{UInt64,Int}(v => k for (k,v) in enumerate(u_tract))
    demog_idx = Dict{UInt8,Int}(v => k for (k,v) in enumerate(u_demog))

    data = zeros(Int, length(u_tract), length(u_demog))

    for agent in raw
        tract = div.(agent.fips_code, 10)
        d = getfield(agent, demog)
        data[tract_idx[tract], demog_idx[d]] += 1
    end

    return DictMatrix(tract_idx, demog_idx, data)
end
# ============================================================================ #
function main(ifile::AbstractString)

    agent_dir = "/Users/palexander/Documents/emerge+radium/input-data/agent_db"

    data = Epicast.read_eventfile(Epicast.EventData, ifile)

    states = unique(UrbanPop.agent_state.(getfield.(data.events, :agent_id)))
    ids = getfield.(filter(x -> x.state == 0x01 && x.context != 0xff, data.events), :agent_id)

    demog = UrbanPop.fetch_agent_demographics(ids, agent_dir)

    tmp = count_unique(getfield.(demog, :person_race))
    if length(tmp) < 7
        tmp = vcat(tmp[1:4], 0x04 => 0, tmp[5:end])
    end

    tmp2 = get_demographic(agent_dir, states, :person_race)

    return tmp, tmp2
end
# ============================================================================ #
end # module EpicastAnalysis
