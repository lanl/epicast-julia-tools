module Workerflow

using DelimitedFiles, CSV
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
function convert_workerflow3(idir::AbstractString, ofile::AbstractString; thr::Integer=0)
    files = find_files(idir, r".*_JT00_\d{4}\.csv\.gz")
    data = Dict{Tuple{Int,Int},Int}()
    for file in files
        try
            add_file!(data, file)
        catch err
            println("[FAIL]: \"$(file)\"")
            rethrow(err)
        end
        println("DONE: ", basename(file))
    end

    if thr > 0
        filter!(data) do k
            k.second > thr
        end
    end

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
function add_file!(d::Dict{Tuple{Int,Int},Int}, ifile::AbstractString)
    block2tract = 10000

    csv = CSV.File(ifile)

    for k in 1:length(csv)
        # h_ = home, w_ = work
        key = (div(csv.h_geocode[k], block2tract),
            div(csv.w_geocode[k], block2tract))
        n = get(d, key, 0)
        d[key] = n + csv.S000[k]
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
end
