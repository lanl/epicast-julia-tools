module Workerflow

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

using DelimitedFiles, CSV
using EpicastTables
# ---------------------------------------------------------------------------- #
function convert_workerflow(ifile::AbstractString, ofile::AbstractString)

    d = readdlm(ifile, '\t', UInt32)

    # remove trailing 0 row (used to mark the EOF?)
    if all(iszero, d[end,:])
        d = d[1:end-1,:]
    end

    # total number of tracts in this file
    n_tract = length(unique(vcat(d[:,1], d[:,2])))
    
    # total number of "entries" (i.e. non-zero matrix elements)
    n_entry = size(d, 1)

    open(ofile, "w") do io
        write(io, UInt64(n_tract))
        write(io, UInt64(n_entry))

        # each "entry" is: from to number..., so reshape to the correct order
        write(io, reshape(d', :, 1))
    end

    return ofile
end
# ---------------------------------------------------------------------------- #
function wf_sort(x, y, mp)
    from1 = mp[x[1]]
    from2 = mp[y[1]]
    return from1 == from2 ? isless(mp[x[2]], mp[y[2]]) : isless(from1, from2)
end
# ---------------------------------------------------------------------------- #
# now works for both ascii and binary workerflow files, however, wf_file and
# tract_file ***MUST*** use the same "ids" or correct results
function convert_workerflow2(wf_file::AbstractString, tract_file::AbstractString,
    ofile::AbstractString)

    if endswith(wf_file, ".dat")
        d = readdlm(wf_file, '\t', UInt32, '\n')
    elseif endswith(wf_file, ".bin")
        nrow = Int(stat(wf_file).size / (sizeof(UInt32) * 3))
        d = zeros(UInt32, 3, nrow)
        open(wf_file, "r") do io
            read!(io, d)
        end
        d = permutedims(d, (2,1))
    else
        error("Unrecognized wf file extension: \"$(wf_file)\"")
    end

    # remove trailing 0 row (used to mark the EOF?)
    if all(iszero, d[end,:])
        d = d[1:end-1,:]
    end

    # total number of "entries" (i.e. non-zero matrix elements)
    n_entry = size(d, 1)

    ifo = readdlm(tract_file, Int, skipstart=1)

    mp = Dict{UInt32,UInt64}()

    for k in 1:size(ifo, 1)
        id = UInt32(ifo[k,1])

        # full tract fips = county ips * 10^6 + tact fips
        fips = ifo[k,4] * 10^6 + ifo[k,5]
        mp[id] = UInt64(fips)
    end

    d = sortslices(d, dims=1, lt=(x,y)->wf_sort(x, y, mp))

    # total number of tracts in this file
    n_tract = length(unique(vcat(d[:,1], d[:,2])))

    open(ofile, "w") do io
        write(io, UInt64(n_tract))
        write(io, UInt64(n_entry))
        # convert "myID"s to FIPS codes
        for k in 1:size(d, 1)
            write(io, mp[d[k,1]]) # UInt64
            write(io, mp[d[k,2]]) # UInt64
            write(io, d[k,3])     # UInt32
        end
    end

end

# ---------------------------------------------------------------------------- #
function find_files(idir::AbstractString, re::Regex=r".*")
    out = String[]
    for name in readdir(idir)
        path = joinpath(idir, name)
        if isfile(path) && match(re, name) != nothing
            push!(out, path)
        end
    end
    return out
end
# ---------------------------------------------------------------------------- #
tract_sort(x, y) = x[1] == y[1] ? isless(x[2], y[2]) : isless(x[1], y[1])
# ---------------------------------------------------------------------------- #
struct LodesMatrixXfm
    adj::Matrix{Bool}
    tract_list::Set{Int}
    xfm::Dict{Int,Vector{Int}}
end
# ---------------------------------------------------------------------------- #
function xfm(x::LodesMatrixXfm, bg_fips::Integer)
    tract = div(bg_fips, 10)
    out = Int[]
    
    if in(tract, x.tract_list)
        # block groups belongs to a UP tract based on FIPS hierarchy
        push!(out, tract)
    end

    if haskey(x.xfm, bg_fips)
        # block group also contributes to another UP tract based on cross-walk
        append!(out, x.xfm[bg_fips])
    end
    return isempty(out) ? [-1] : out
end
# ---------------------------------------------------------------------------- #
# from and to are block-group fips codes
function xfm_entry(x::LodesMatrixXfm, from_bg::Integer, to_bg::Integer)
    bg2state = 10^10
    from = [-1]
    to = [-1]
    if x.adj[div(from_bg, bg2state), div(to_bg, bg2state)]
        # from and to states are considered ~ "neighbors"
        from = xfm(x, from_bg)
        to = xfm(x, to_bg)
    end

    return from, to
end
# ---------------------------------------------------------------------------- #
function get_lodes_tracts(idir::AbstractString)
    files = find_files(idir, r".*_JT00_\d{4}\.csv\.gz")
    from = Set{Int}()
    to = Set{Int}()
    for file in files
        add_tracts!(file, from, to)
    end

    return from, to
end
# ---------------------------------------------------------------------------- #
function get_lodes_matrix(idir::AbstractString, x::LodesMatrixXfm, thr::Integer)
    files = find_files(idir, r".*_JT00_\d{4}\.csv\.gz")
    data = Dict{Tuple{Int,Int},Int}()
    for file in files
        add_file!(data, file, x)
        println("DONE: ", basename(file))
    end

    if thr > 0
        filter!(data) do k
            k.second > thr
        end
    end

    return data
end
# ---------------------------------------------------------------------------- #
function get_lodes_mapping(idir::AbstractString, x::LodesMatrixXfm, thr::Integer)
    files = find_files(idir, r".*_JT00_\d{4}\.csv\.gz")
    data = Dict{Int,Vector{Int}}()
    for file in files
        add_mapping!(data, file, x)
        println("DONE: ", basename(file))
    end

    for (k,v) in data
        data[k] = unique(v)
    end

    return data
end
# ---------------------------------------------------------------------------- #
function load_nhgis_crosswalk(ifile::AbstractString, tracts::Vector{<:Integer})
    csv = CSV.File(ifile)

    mp = Dict{Int, Vector{Int}}()
    xw_src = csv.bg2010ge
    xw_dst = csv.tr2020ge
    xw_tr = Set(xw_dst)

    # n1, n2, n3, n4 = (0,0,0,0)

    for tr in tracts
        # UP tract directy exists in crosswalk as a combination of 2010
        # block-groups 
        if in(tr, xw_tr)
            idx = findall(isequal(tr), xw_dst)
            # n1 += !isempty(idx)
        else
            # was the UP tract split into multiple tracts for 2020 census?
            # e.g. 36053030600 -> [36053030601, 36053030602]
            y = div(tr, 100)
            idx = findall(xw_dst) do x
                div(x, 100) == y
            end
            # n2 += !isempty(idx)
        end
        
        if isempty(idx)
            # no good candidates still, see if any 2010 block-groups appear to
            # be members of this urbanpop tract (could be a 2010 tract that was
            # removed in 2020, but was there in 2019 ACS), if so use those
            idx = findall(xw_src) do x
                div(x, 10) == tr
            end
            # n3 += !isempty(idx)
            if isempty(idx)
                # if not... find the 2020 tract w/ the closest FIPS code and
                # issue a warning (in practice this only happens for Ogala
                # Lakota County, SD, in which case data from a tract w/in the
                # same county will be used instead)
                cnty = div(tr, 10^6)
                mn, kmn = findmin(xw_dst) do x
                    div(x, 10^6) == cnty ? abs(x - tr) : +Inf
                end
                isinf(mn) && error("Completely failed to map tract $(tr)")
                idx = [kmn]
                @warn("Failed for tract $(tr), using $(xw_dst[kmn])")
                # n4 += !isempty(idx)
            end
            
        end

        # map from 2010 block-group to one or more UrbanPop tract
        for k in idx
            if !haskey(mp, csv.bg2010ge[k])
                mp[csv.bg2010ge[k]] = Int[]
            end
            push!(mp[csv.bg2010ge[k]], tr)
        end
    end

    # println("step1: $(n1), step2: $(n2), step3: $(n3), step4: $(n4), sum = $(n1+n2+n3+n4)")

    return mp
end
# ---------------------------------------------------------------------------- #
function convert_lodes_workerflow(lodes_dir::AbstractString, xwlk_file::AbstractString,
    up_tracts::Vector{<:Integer}, adj::Matrix{Bool}, ofile::AbstractString; thr::Integer=0)

    # tract FIPS for given LODES data
    from, to = get_lodes_tracts(lodes_dir)

    # tracts that we need a crosswalk mapping for
    xw_tracts = filter(up_tracts) do x
        return !in(x, from) || !in(x, to)
    end

    # step1: 67, step2: 5, step3: 5, step4: 8, sum = 85
    # construct xwlk mapping 2010 block-group -> UP tract
    mp = load_nhgis_crosswalk(xwlk_file, xw_tracts)


    # general transform object, that given 2010 block-groups FIPS will return
    # the UrbanPop tract FIPS to which that block-groups "belongs"
    xfm = LodesMatrixXfm(adj, Set(up_tracts), mp)

    data = get_lodes_matrix(lodes_dir, xfm, thr)

    # @time data = get_lodes_mapping(lodes_dir, xfm, thr)
    # return data

    tracts = sort!(collect(keys(data)), lt=tract_sort)

    n_tract = length(Set(vcat(getindex.(tracts, 1), getindex.(tracts, 2))))

    open(ofile, "w") do io
        write(io, UInt64(n_tract))
        write(io, UInt64(length(data)))
        for k in tracts
            write(io, UInt64(k[1]))
            write(io, UInt64(k[2]))
            write(io, UInt32(data[k]))
        end
    end

    return nothing
end
# ---------------------------------------------------------------------------- #
function tract_match(fip::Integer, fips::Vector{<:Integer})
    if fip in fips
        return [fip]
    else
        y = div(fip, 100)
        return filter(fips) do x
            div(x, 100) == y
        end
    end
end
# ---------------------------------------------------------------------------- #
function add_tracts!(ifile::AbstractString, from::Set{Int}, to::Set{Int})
    block2tract = 10^4
    csv = CSV.File(ifile)
    union!(from, Set(div.(csv.h_geocode, block2tract)))
    union!(to, Set(div.(csv.w_geocode, block2tract)))
    return nothing
end
# ---------------------------------------------------------------------------- #
function add_mapping!(d::Dict{Int,Vector{Int}}, ifile::AbstractString, x::LodesMatrixXfm)
    block2bg = 10^3

    csv = CSV.File(ifile)

    for k in 1:length(csv)
        # h_ = home, w_ = work
        from, to = xfm_entry(x, div(csv.h_geocode[k], block2bg),
                div(csv.w_geocode[k], block2bg))

        if any(>(0), from)
            if haskey(d, csv.h_geocode[k])
                append!(d[csv.h_geocode[k]], filter(>(0), from))
            else
                d[csv.h_geocode[k]] = filter(>(0), from)
            end
        end

        if any(>(0), to)
            if haskey(d, csv.w_geocode[k])
                append!(d[csv.w_geocode[k]], to)
            else
                d[csv.w_geocode[k]] = filter(>(0), to)
            end
        end

    end

    return d
end
# ---------------------------------------------------------------------------- #
function add_file!(d::Dict{Tuple{Int,Int},Int}, ifile::AbstractString, x::LodesMatrixXfm)
    block2bg = 10^3

    csv = CSV.File(ifile)

    for k in 1:length(csv)
        # h_ = home, w_ = work
        from, to = xfm_entry(x, div(csv.h_geocode[k], block2bg),
                div(csv.w_geocode[k], block2bg))

        for i in eachindex(from)
            for j in eachindex(to)
                if from[i] > 0 && to[j] > 0
                    key = (from[i], to[j])
                    n = get(d, key, 0)
                    d[key] = n + csv.S000[k]
                end
            end
        end
    end

    return d
end
# ---------------------------------------------------------------------------- #
function read_workerflow_file(ifile::AbstractString)
    from, to, n = open(ifile, "r") do io
        n_tract = read(io, UInt64)
        n_entry = read(io, UInt64)
        from = Vector{UInt64}(undef, n_entry)
        to = Vector{UInt64}(undef, n_entry)
        n = Vector{UInt32}(undef, n_entry)
        for k = 1:n_entry
            from[k] = read(io, UInt64)
            to[k] = read(io, UInt64)
            n[k] = read(io, UInt32)
        end
        return from, to, n
    end

    return Dict{Tuple{Int,Int}, Int}((Int(f),Int(t)) => n for (f,t,n) in
        zip(from,to,n))
end
# ---------------------------------------------------------------------------- #
struct Slice{N1,N2} end
# ---------------------------------------------------------------------------- #
function get_slice(::Type{Slice{N1,N2}}, d::Dict{Tuple{Int,Int}, Int}, idx::Integer) where {N1,N2}
    slice = Dict{Int,Int}()
    for (k,v) in d
        if k[N1] == idx
            slice[k[N2]] = d[k]
        end
    end
    return slice
end
# ---------------------------------------------------------------------------- #
get_column(d::Dict{Tuple{Int,Int}, Int}, idx::Integer) = get_slice(Slice{2,1}, d, idx)
get_row(d::Dict{Tuple{Int,Int}, Int}, idx::Integer) = get_slice(Slice{1,2}, d, idx)
# ---------------------------------------------------------------------------- #
vec_sum(d::Dict{Int,Int}, rm::Integer) = sum(values(d)) - d[rm]
# ---------------------------------------------------------------------------- #
function all_tracts(d::Dict{Tuple{Int,Int},Int})
    tracts = Set{Int}()
    for (k,v) in d
        push!(tracts, k[1])
        push!(tracts, k[2])
    end
    return tracts
end
# ---------------------------------------------------------------------------- #
function net_flow(d::Dict{Tuple{Int,Int},Int}, tracts::Set{<:Integer})

    flow = zeros(Float64, length(tracts),1)
    fips_index = Dict{UInt64,Int}()
    k = 1
    for tr in tracts
        inflow = sum(values(get_column(d, tr)))
        outflow = sum(values(get_row(d, tr)))
        flow[k,1] = inflow - outflow
        fips_index[tr] = k
        k += 1
    end

    var_index = Dict{String,Int}("delta_pop" => 1)

    return FIPSTable{Tract,Float64,2}(fips_index, var_index, flow)
end
# ---------------------------------------------------------------------------- #
# list of neighboring and very nearby* states for each state
# *NOTE: nearby means commuting is reasonable b/t these states
function read_neighbor_file(ifile::AbstractString)
    d = Dict{String,Vector{String}}()
    for line in eachline(ifile)
        list = split(line, " ", keepempty=false)
        d[list[1]] = length(list) > 1 ? list[2:end] : String[]
    end
    return d
end
# ---------------------------------------------------------------------------- #
# POSTAL Abbrv -> FIPS code
function fips_map(ifile::AbstractString)
    d = Dict{String,Int}()
    for (k, st) in enumerate(eachline(ifile))
        if !isempty(st)
            d[st] = k
        end
    end
    return d
end
# ---------------------------------------------------------------------------- #
# FIPS indexed state adjacency matrix
function state_adjacency_matrix(neighbor_file::AbstractString, fips_file::AbstractString)

    neighbors = read_neighbor_file(neighbor_file)
    fips = fips_map(fips_file)

    adjm = zeros(Bool, 56, 56)

    for k in 1:56
        adjm[k,k] = true
    end

    for (k, v) in neighbors
        if length(v) > 1
            r = fips[k]
            try
                c = map(x->fips[x], v)
                adjm[r, c] .= true
            catch err
                @show(k, v)
                rethrow(err)
            end
        end
    end

    return adjm
end
# ---------------------------------------------------------------------------- #
end
