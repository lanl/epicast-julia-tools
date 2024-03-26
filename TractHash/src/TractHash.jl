module TractHash

using DelimitedFiles, Random

const LENGTH = 2^17

function hash_u64(n::UInt64)
    a::UInt64 = n
    a = ~a + a << 21
    a =  a ⊻ a >> 24
    a =  a + a << 3 + a << 8
    a =  a ⊻ a >> 14
    a =  a + a << 2 + a << 4
    a =  a ⊻ a >> 28
    a =  a + a << 31
    return a
end

function hash_u64_2(x::UInt64)
    a::UInt64 = x
    a = (a ⊻ (a >> 30)) * UInt64(0xbf58476d1ce4e5b9)
    a = (a ⊻ (a >> 27)) * UInt64(0x94d049bb133111eb)
    a = a ⊻ (a >> 31)
    return a
end

function tract2index(tract::Int64)
    hsh = hash_u64_2(UInt64(tract))
    idx = ((Base.bitcast(Int64, hsh) & (LENGTH - 1)) + 1)::Int64
    return idx
end

f1(d::Dict, key::Int) = d[key]
f2(d::Vector{Int}, key::Int) = searchsortedfirst(d, key)

function test()
    ifile = "/Users/palexander/Documents/emerge+radium/input-data/tracts.txt"

    d = vec(readdlm(ifile, ',', Int))
    filter!(d) do x
        return 35000000000 < x < 36000000000
    end
    @assert(issorted(d), "Nope...")
    # sort!(d)
    # tmp = map(d) do x
    #     # idx1 =  Base.hashindex(x, 2^17)[1]
    #     idx2 = tract2index(x)
    #     # @assert(idx1 == idx2, "$(x)")
    #     return idx2
    # end

    # @time (tmp = tract2index.(d))

    # @show(length(tmp))
    # @show(length(unique(tmp)))

    # @assert(length(tmp) == length(unique(tmp)), "Fail!!")

    mp = Dict( v => k for (k,v) in enumerate(d))
    return mp
    mem = Vector{Int}(undef, length(d))
    for k = 1:3
        test = shuffle(d)
        print("dict: ")
        @time for k in eachindex(test)
            mem[k] = f1(mp, test[k])
        end
        print("search: ")
        @time for k in eachindex(test)
            mem[k] = f2(d, test[k])
        end
    end
    
end

end # module TractHash
