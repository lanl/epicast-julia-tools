#!/usr/bin/env julia

import Pkg

Pkg.activate("UrbanPop")
using ArgParse
using UrbanPop
using Graphs
using Printf

Id = UInt64

function parse_args(args)
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--in-dir", "-i"
            default="../data"
            help="Path to file containing UrbanPop data"
        "--out-dir", "-o"
            default=""
            help="Path to save resulting data to (defaults to in_dir)"
        "--ave-degree", "-k"
            default=100
            arg_type=Int
            help="Average degree to use in generated network"
        "--beta", "-b"
            default=0.1
            arg_type=Real
            help="Rewiring probability to use when generating network"
        "--states", "-s"
            default=[]
            arg_type=Int
            nargs='+'
            help="A list of state FIPS codes to include in graph (defaults to using all states)"
    end

    args = ArgParse.parse_args(args, s)

    if args["out-dir"] == ""
        args["out-dir"] = args["in-dir"]
    end

    return args
end

@inline fips_state(x) = floor.(Int, x ./ 1e9)
function get_state_offsets(tract_fips::AbstractVector{<:UInt64},
        tract_pops::AbstractVector{T}, states::AbstractVector{<:Int}
        ) where T<:Integer
    tract_states = fips_state(tract_fips)
    if 0 == length(states)
        states = unique(tract_states)
    else
        states = intersect(states, tract_states)
    end
    println(states)

    n_states = length(states)
    offsets = zeros(T, n_states, 3)

    offset = 0
    last_state = 0
    state_idx = 0
    state_pop = 0
    for (state, pop) in zip(tract_states, tract_pops)
        if !in(state, states)
            continue
        end
        if last_state != state
            state_idx += 1
            state_pop = 0

            offsets[state_idx, :] = [state, offset, state_pop]
            last_state = state
            println("state $state has $state_pop people")
        end

        state_pop += pop
        offset += pop
        offsets[state_idx, 3] = state_pop
    end
    println("state $last_state has $state_pop people")

    return offsets, offset
end

n_state_bits = 6
state_shift = sizeof(Id)*8 - n_state_bits
@inline function state_node_range(offset::T, pop::T) where T<:Integer
    return range(offset + 1, offset+pop) .% T
end
function get_agent_ids(total_pop::T,
        state_offsets::Matrix{T}) where T<:Integer
    ids = Vector{T}(range(1, total_pop))

    n_states = size(state_offsets)[1]
    for r in range(1, n_states)
        (state, offset, pop) = state_offsets[r, :]
        state_start = (state << state_shift)
        ids[state_node_range(offset, pop)] =
            range(state_start, state_start+pop-1)

        first = ids[offset + 1] .% Int128
        last = ids[offset+pop] .% Int128
        println("state $state has $pop people with ids: [$first, $last]")
    end

    return ids
end

function get_edgelist(graph::AbstractGraph{T}) where T<:Integer
    return reduce(hcat, [[Graphs.src(e) .% T,
                          Graphs.dst(e) .% T]
                         for e in Graphs.edges(graph)])
end

function get_edge_offsets(graph::AbstractGraph{T},
        nodes::AbstractVector{T}) where T<:Integer
    degrees = Graphs.degree(graph, nodes)
    person_edge_offsets = cumsum(degrees) - degrees
    append!(person_edge_offsets, sum(degrees) .% T)
    return person_edge_offsets
end

function get_dsts(graph::AbstractGraph{T},
        nodes::AbstractVector{T}) where T<:Integer
    return reduce(vcat, Graphs.SimpleGraphs.adj(graph)[nodes])
end

function print_summary(graph::AbstractGraph{T}) where T<:Integer
    g_type = summary(graph)
    ne = Graphs.ne(graph)
    nv = Graphs.nv(graph)
    println("Generated $g_type with $ne edges, $nv nodes")
end

function write_header!(out_stream::IOStream,
        graph::AbstractGraph{T},
        all_pops::AbstractVector{<:UInt64}) where T<:Integer
    ne = Graphs.ne(graph)
    nv = Graphs.nv(graph)
    write(out_stream, ne .% T)
    write(out_stream, nv .% T)

    offsets = get_edge_offsets(graph)
    write(out_stream, offsets .% T)
end

function write_state(out_file::AbstractString,
        graph::AbstractGraph{T}, state::T,
        offset::T, n_nodes::T,
        agent_ids::AbstractVector{T}) where T<:Integer
    open(out_file, "w") do stream
        state_nodes = state_node_range(offset, n_nodes)
        edge_offsets = get_edge_offsets(graph, state_nodes)

        n_edges = edge_offsets[end]
        println("State $state: Saving $n_edges edges, $n_nodes nodes to $out_file")
        write(stream, n_edges .% T)
        write(stream, n_nodes .% T)
        write(stream, edge_offsets .% T)

        dsts = agent_ids[get_dsts(graph, state_nodes)]
        write(stream, dsts)
    end
end

function main(args)
    all_tracts, all_pop = UrbanPop.all_tract_data(args["in-dir"])
    all_tracts = all_tracts .% UInt64
    all_pop = all_pop .% Id
    total_pop = sum(all_pop) .% Id
    state_offsets, used_pop = get_state_offsets(all_tracts, all_pop, args["states"])
    println("Total pop: $total_pop, used pop: $used_pop")

    agent_ids = get_agent_ids(used_pop, state_offsets)

    g = Graphs.newman_watts_strogatz(Id(used_pop), args["ave-degree"], args["beta"])
    #print_summary(g)

    n_states = size(state_offsets)[1]
    out_dir = args["out-dir"]
    for r in range(1, n_states)
        (state, offset, pop) = state_offsets[r, :]
        out_file = @sprintf("%02d.social.bin", state)
        write_state("$out_dir/$out_file", g, state, offset,
                    pop, agent_ids)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(parse_args(ARGS))
end
