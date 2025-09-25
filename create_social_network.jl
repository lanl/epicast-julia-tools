#!/usr/bin/env julia

import Pkg

Pkg.activate("UrbanPop")
using ArgParse
using UrbanPop
using Graphs
using Printf

AgentId = UInt32
EdgeId = UInt64

# ============================================================================ #
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
# ============================================================================ #
@inline fips_state(x) = floor.(Int, x ./ 1e9)
# ---------------------------------------------------------------------------- #
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
            println("state $state has $state_pop people")
            state_idx += 1
            state_pop = 0

            offsets[state_idx, :] = [state, offset, state_pop]
            last_state = state
        end

        state_pop += pop
        offset += pop
        offsets[state_idx, 3] = state_pop
    end
    println("state $last_state has $state_pop people")

    return offsets, offset
end
# ---------------------------------------------------------------------------- #
function get_state_offsets(in_dir::AbstractString, states::AbstractVector{<:Integer})
    all_tracts, all_pop = UrbanPop.all_tract_data(in_dir)
    all_tracts = all_tracts .% UInt64
    all_pop = all_pop .% AgentId
    total_pop = sum(all_pop) % AgentId
    state_offsets, used_pop = get_state_offsets(all_tracts, all_pop, states)

    used_pop = used_pop % AgentId
    println("States: $states, total pop: $total_pop, used pop: $used_pop")

    return state_offsets, used_pop
end
# ============================================================================ #
n_state_bits = 6
state_shift = sizeof(AgentId)*8 - n_state_bits
@inline function state_node_range(offset::T, pop::T) where T<:Integer
    return range(offset + 1, offset+pop) .% T
end
# ---------------------------------------------------------------------------- #
@inline id_to_state(id) = id .>> state_shift
# ---------------------------------------------------------------------------- #
@inline id_to_local_idx(id) = id .- (id_to_state(id) .<< state_shift)
# ---------------------------------------------------------------------------- #
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
# ============================================================================ #
function get_edgelist(graph::AbstractGraph{T}) where T<:Integer
    return reduce(hcat, [[Graphs.src(e) .% T,
                          Graphs.dst(e) .% T]
                         for e in Graphs.edges(graph)])
end
# ---------------------------------------------------------------------------- #
function get_edge_offsets(graph::AbstractGraph{T},
        nodes::AbstractVector{T}) where T<:Integer
    degrees = Graphs.degree(graph, nodes)
    person_edge_offsets = cumsum(degrees) - degrees
    append!(person_edge_offsets, sum(degrees) .% EdgeId)
    return person_edge_offsets
end
# ============================================================================ #
function get_dsts(graph::AbstractGraph{T},
        nodes::AbstractVector{T}) where T<:Integer
    return reduce(vcat, Graphs.SimpleGraphs.adj(graph)[nodes])
end
# ============================================================================ #
function print_summary(graph::AbstractGraph{T}) where T<:Integer
    g_type = summary(graph)
    ne = Graphs.ne(graph)
    nv = Graphs.nv(graph)
    println("Generated $g_type with $ne edges, $nv nodes")
end
# ============================================================================ #
function write_state(out_file::AbstractString,
        graph::AbstractGraph{T}, state::T,
        offset::T, n_nodes::T, agent_ids::AbstractVector{T},
        states::AbstractSet{T}) where T<:Integer
    open(out_file, "w") do stream
        state_nodes = state_node_range(offset, n_nodes)
        edge_offsets = get_edge_offsets(graph, state_nodes)

        n_edges = edge_offsets[end]
        println("State $state: Saving $n_edges edges, $n_nodes nodes to $out_file")
        write(stream, n_edges .% EdgeId)
        write(stream, n_nodes .% T)
        write(stream, edge_offsets .% EdgeId)

        dsts = agent_ids[get_dsts(graph, state_nodes)]
        dst_states = Set(id_to_state(dsts))
        if states != dst_states
            extra_states = setdiff(dst_states, states)
            println("Error: states $extra_states present in output while not selected")
            exit()
        end

        write(stream, dsts)
    end
end
# ---------------------------------------------------------------------------- #
function write_social_network(out_dir::AbstractString,
    graph::AbstractGraph{T}, state_offsets::AbstractMatrix{T}
    ) where T<:Integer

    used_pop = Graphs.nv(graph) .% T
    agent_ids = get_agent_ids(used_pop, state_offsets)
    n_states = size(state_offsets)[1]
    states = Set(state_offsets[:,1])
    for r in range(1, n_states)
        (state, offset, pop) = state_offsets[r, :]
        out_file = @sprintf("%02d.social.bin", state)
        write_state("$out_dir/$out_file", graph, state, offset,
                    pop, agent_ids, states)
    end
end
# ============================================================================ #
HEADER_LENGTH = sizeof(EdgeId) + sizeof(AgentId)
# ----------------------------------------------------------------------------  #
function read_header(in_path::AbstractString)
    n_edges, n_nodes = (0, 0)
    open(in_path, "r") do stream
        n_edges = read(stream, EdgeId)
        n_nodes = read(stream, AgentId)
    end
    return n_edges, n_nodes
end
# ---------------------------------------------------------------------------- #
function read_person_friends(in_path::AbstractString, local_idx::Integer)
    n_edges, n_nodes = read_header(in_path)
    friends = nothing
    open(in_path, "r") do stream
        seek(stream, HEADER_LENGTH + local_idx * sizeof(EdgeId))
        start_offset = read(stream, EdgeId)
        end_offset = read(stream, EdgeId)

        friends = zeros(AgentId, end_offset - start_offset)
        seek(stream, HEADER_LENGTH + (n_nodes + 1) * sizeof(EdgeId)
            + start_offset * sizeof(AgentId))
        read!(stream, friends)
    end

    return friends
end
# ============================================================================ #
function path_to_state(path::AbstractString)
    return parse(Int32, basename(path)[1:2])
end
# ============================================================================ #
#struct Person
#    id::AgentId
#    edge_offset::EdgeId
#    n_edges::AgentId
#end
mutable struct StateNetwork
    people_offsets::Array{EdgeId,1}
    edges::Array{AgentId,1}
    n_edges::EdgeId
    n_nodes::AgentId
    node_offset::AgentId
    fips_code::Int8
end
# ---------------------------------------------------------------------------- #
function global_to_local(s::StateNetwork, global_idx::Integer)
    return global_idx % AgentId - s.node_offset + 1
end
# ---------------------------------------------------------------------------- #
function Base.getindex(s::StateNetwork, global_idx::Integer)
    local_idx = global_to_local(s, global_idx)

    start = s.people_offsets[local_idx] + 1
    stop = s.people_offsets[local_idx + 1] + 1

    return s.edges[start:stop]
end
# ---------------------------------------------------------------------------- #
function read_state_network_header(in_path::AbstractString)
    n_edges, n_nodes = read_header(in_path)

    return StateNetwork(
        zeros(EdgeId, n_nodes + 1),
        zeros(AgentId, n_edges),
        n_edges,
        n_nodes,
        0,
        path_to_state(in_path)
    )
end
# ---------------------------------------------------------------------------- #
function read_offsets!(this::StateNetwork, in_path::AbstractString)
    open(in_path, "r") do stream
        seek(stream, HEADER_LENGTH)
        read!(stream, this.people_offsets)
    end
end
# ---------------------------------------------------------------------------- #
function read_edges!(this::StateNetwork, in_path::AbstractString)
    open(in_path, "r") do stream
        seek(stream, HEADER_LENGTH + (this.n_nodes + 1) * sizeof(EdgeId))
        read!(stream, this.edges)
    end
end
# ---------------------------------------------------------------------------- #
function read_state_network!(this::StateNetwork, in_path::AbstractString)
    read_offsets!(this, in_path)
    read_edges!(this, in_path)
end
# ---------------------------------------------------------------------------- #
function read_state_network(in_path::AbstractString; header_only::Bool=false)
    this = read_state_network_header(in_path)

    if !header_only
        read_state_network!(this, in_path)
    end

    return this
end
# ============================================================================ #
mutable struct SocialNetwork
    states::Vector{StateNetwork}
    state_indices::Dict{Int8,UInt8}
    n_edges::EdgeId
    n_nodes::AgentId
end
# ---------------------------------------------------------------------------- #
function Base.getindex(this::SocialNetwork, state::Int8)
    return this.states[this.state_indices[state]]
end
# ---------------------------------------------------------------------------- #
function Base.getindex(this::SocialNetwork, states::AbstractVector{Int8})
    return this.states[this.state_indices[states]]
end
# ---------------------------------------------------------------------------- #
function convert_to_global!(this::SocialNetwork, state::StateNetwork)
    state_offsets = Dict(k => this.states[i].node_offset for (k, i) in this.state_indices)
    state_ids = id_to_state(state.edges) .% Int8
    edge_offsets = getindex.(Ref(state_offsets), state_ids)
    state.edges = id_to_local_idx(state.edges) + edge_offsets
end
# ---------------------------------------------------------------------------- #
function read_social_network(in_dir::AbstractString;
    states::AbstractVector{<:Integer}=Vector{Int8}(), header_only::Bool=false,
    to_global::Bool=false)
    files = readdir(in_dir, join=true)
    files = filter(f -> occursin(r"\d\d\.social\.bin", f), files)
    if length(states) > 0
        states = Set(states)
        files = filter(f -> path_to_state(f) in states, files)
    end

    states = map(f -> read_state_network(f; header_only=header_only), files)
    state_indices = Dict(states[i].fips_code => i for i in 1:length(states))

    node_offset = 0
    for s in states
        s.node_offset = node_offset
        node_offset += s.n_nodes
    end

    n_edges = sum(map(s -> s.n_edges, values(states)))
    n_nodes = sum(map(s -> s.n_nodes, values(states)))
    this = SocialNetwork(states, state_indices, n_edges, n_nodes)

    if !header_only && to_global
        for s in states
            convert_to_global!(this, s)
        end
    end

    return this
end
# ============================================================================ #
function main(args)
    state_offsets, used_pop = get_state_offsets(args["in-dir"], args["states"])
    g = Graphs.newman_watts_strogatz(AgentId(used_pop), args["ave-degree"], args["beta"])
    #print_summary(g)

    write_social_network(args["out-dir"], g, state_offsets)
end
# ---------------------------------------------------------------------------- #
if abspath(PROGRAM_FILE) == @__FILE__
    main(parse_args(ARGS))
end
