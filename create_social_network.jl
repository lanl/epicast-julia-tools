#!/usr/bin/env julia

import Pkg

Pkg.activate("UrbanPop")
using UrbanPop
using Graphs

function write_header!(out_stream::IOStream,
        graph::AbstractGraph,
        all_pops::AbstractVector{<:Int64})
    n_edges = Graphs.degree(graph)
    person_edge_offsets = cumsum(n_edges) - n_edges
    tract_person_offsets = cumsum(all_pop) - all_pop
    tract_edge_offsets = [person_edge_offsets[i+1] for i in tract_person_offsets]

    write(out_stream, tract_edge_offsets)
    write(out_stream, "\n");
end

function get_edgelist(graph::AbstractGraph)
    return reduce(vcat, [[Graphs.src(e), Graphs.dst(e)] for e in Graphs.edges(graph)])
end

function main()
    in_dir = "../data"
    tracts_file = "$in_dir/49.tracts.bin"
    agents_file = "$in_dir/49.agents.bin"
    out_file = "$in_dir/49.social.bin"

    all_tracts, all_pop = UrbanPop.all_tract_data(in_dir)
    total_pop = sum(all_pop)

    g = Graphs.newman_watts_strogatz(total_pop, 100, 0.1)

    # Todo: write seperate file for each state, similar to
    # how tracts and agents work
    open(out_file, "w") do out_stream
        write_header!(out_stream, g, all_pop)

        edges = get_edgelist(g)
        write(out_stream, edges)
    end
end

main()
