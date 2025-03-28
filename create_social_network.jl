#!/usr/bin/env julia

import Pkg

Pkg.activate("UrbanPop")
using UrbanPop
using Graphs

function get_edgelist(graph::AbstractGraph{T}) where T<:Integer
    return reduce(hcat, [[(Graphs.src(e) - 1) .% T,
                          (Graphs.dst(e) - 1) .% T]
                         for e in Graphs.edges(graph)])
end

function get_edge_offsets(graph::AbstractGraph{T}
    ) where T<:Integer
    degrees = Graphs.degree(graph)
    person_edge_offsets = cumsum(degrees) - degrees
    append!(person_edge_offsets, 2*Graphs.ne(graph) .% T)
    return person_edge_offsets
end

function get_dsts(graph::AbstractGraph{T}) where T<:Integer
    return reduce(vcat, Graphs.SimpleGraphs.adj(graph))
end

function write_header!(out_stream::IOStream,
        graph::AbstractGraph{T},
        all_pops::AbstractVector{<:Int64}) where T<:Integer
    g_type = summary(graph)
    ne = Graphs.ne(graph)
    nv = Graphs.nv(graph)
    println("Generated $g_type with $ne edges, $nv nodes")
    write(out_stream, ne .% T)
    write(out_stream, nv .% T)

    offsets = get_edge_offsets(graph)
    write(out_stream, offsets .% T)
end

function main()
    in_dir = "../data"
    tracts_file = "$in_dir/49.tracts.bin"
    agents_file = "$in_dir/49.agents.bin"
    out_file = "$in_dir/49.social.bin"

    all_tracts, all_pop = UrbanPop.all_tract_data(in_dir)
    total_pop = sum(all_pop)

    g = Graphs.newman_watts_strogatz(UInt64(total_pop), 100, 0.1)

    # Todo: write seperate file for each state, similar to
    # how tracts and agents work
    open(out_file, "w") do stream
        write_header!(stream, g, all_pop)

        dsts = get_dsts(g)
        write(stream, dsts)
    end
end

main()
