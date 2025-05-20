#!/usr/bin/env julia

import Pkg

Pkg.activate("UrbanPop")
using UrbanPop
using TOML
using Printf

@inline fips_state(x) = floor.(Int, x ./ 1e9)
function get_state_offsets(tract_fips::AbstractVector{<:UInt64},
        tract_pops::AbstractVector{T}) where T<:Integer
    states = fips_state(tract_fips)
    n_states = length(unique(states))
    offsets = zeros(T, n_states, 3)

    offset = 0
    last_state = 0
    state_idx = 0
    state_pop = 0
    for (state, pop) in zip(states, tract_pops)
        if last_state != state
            state_idx += 1
            state_pop = 0

            offsets[state_idx, :] = [state, offset, state_pop]
            last_state = state
        end

        state_pop += pop
        offset += pop
        offsets[state_idx, 3] = state_pop
    end

    return offsets
end

function region_list_to_pop(regions, region_pops, region_states)
    tmp = map(x -> @sprintf("hhs_%.02d.toml", x), regions)
    pop = sum(map(r -> region_pops[r][1], tmp))
    prop = sum(map(r -> region_pops[r][2], tmp))
    states = reduce(vcat, map(r -> region_states[r], tmp))
    return Dict("population" => pop, "proportion of US population" => prop, "geography" => sort(states))
end

function write_toml(data, filepath)
    open(filepath, "w") do f
        TOML.print(f, data)
    end
end

function main()
    regions = map(n -> @sprintf("hhs_%.02d.toml", n), range(1,10))

    region_states = Dict()
    for r in regions
        tmp = TOML.parsefile("../params/$r")
        region_states[r] = tmp["geography"]
    end

    all_tracts, all_pop = UrbanPop.all_tract_data("/vast/home/jkitson/shared/input_data/urbanpop_v2/")
    offsets = get_state_offsets(all_tracts .% UInt , all_pop .% UInt)
    state_pops = Dict()
    for i in range(1, 51)
        (state, offset, pop) = offsets[i, :]
        state_pops[state .% Int] = pop .% Int
    end

    total_pop = sum(all_pop)
    region_pops = Dict()
    for (k, v) in region_states
        tmp = sum(map(s -> state_pops[s], v))
        region_pops[k] = (tmp, tmp / total_pop)
    end

    half = region_list_to_pop([10,9,8,7], region_pops, region_states)
    write_toml(half, "half_us.toml")

    quarter = region_list_to_pop([2,3,4,5], region_pops, region_states)
    write_toml(quarter, "quarter_us.toml")
end

main()
