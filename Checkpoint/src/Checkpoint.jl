module Checkpoint

using DelimitedFiles, Dates
using Workerflow, UrbanPop, EpicastTypes

# ============================================================================ #
const N_AGENT_PER_COMMUNITY = 2000
const CELL_UNUSED = 0xffffffff # marks a community/cell as unused
# ============================================================================ #
struct CheckpointHeader
    n_cell::UInt64
    n_comm::UInt64

    n_agent_per_comm::UInt32
    n_tract::UInt32
    n_hub::UInt32

    n_processor::UInt16
    celldata_size::UInt16
    particle_size::UInt16
    parameters_size::UInt16
    tract_size::UInt16
    urbanpop_version::UInt16
end
CheckpointHeader() = CheckpointHeader(0,0,0,0,0,0,0,0,0,0,0)
function Base.show(io::IO, x::CheckpointHeader)
    println(io, Int(x.n_cell), ", ", Int(x.n_comm), ", ",
        Int(x.n_agent_per_comm), ", ", Int(x.n_tract), ", ", Int(x.n_hub),
        ", ", Int(x.n_processor))
end
# ============================================================================ #
struct CountyCases
    cumulative_cases::UInt32
    new_cases::UInt32
    cumulative_deaths::UInt32
    new_deaths::UInt32
end
# ============================================================================ #
function remove_unused!(c::Vector{Community{A,B}}) where {A,B}
    return filter!(c) do x
        return x.cell_data.tract != CELL_UNUSED && x.cell_data.community != CELL_UNUSED
    end
end
# ============================================================================ #
function check_communities(c::Vector{Community{A,B}}) where {A,B}
    for k in eachindex(c)
        tr = Int(c[k].cell_data.tract)
        cm = Int(c[k].cell_data.community)
        if (c[k].cell_data.n_family > length(c[k].particles)) ||
                (k > 1 && (tr == 0 || cm == 0)) ||
                (c[k].cell_data.n_family == 0 && length(c[k].particles) > 0)
            
            print(c[k].cell_data)
            println("  n_agent: $(length(c[k].particles)), k = $(k)")
        end
    end
end
# ============================================================================ #
function is_approx_equal(d1::Vector{Community{A,B}}, d2::Vector{Community{A,B}})  where {A,B}
    tmp = map(d1, d2) do a,b
        return a.cell_data == b.cell_data && length(a.particles) == length(b.particles)
    end
    
    return all(tmp)
end
# ============================================================================ #
read_checkpoint_file(ifile::AbstractString) = read_checkpoint_file(CellData, Particle, ifile)
# ---------------------------------------------------------------------------- #
function read_checkpoint_file(::Type{C}, ::Type{P}, ifile::AbstractString) where {C<:AbstractCellData, P<:AbstractParticle}

    return open(ifile, "r") do io
        seek(io, sizeof(Int64))
        ref = Ref(CheckpointHeader())
        read!(io, ref)
        hdr = ref[]

        @assert(hdr.celldata_size == sizeof(C), "$(hdr.celldata_size) != $(sizeof(C))")
        @assert(hdr.particle_size == sizeof(P))
        @assert(hdr.tract_size == sizeof(UrbanPop.Tract))
        @assert(hdr.urbanpop_version == 2)

        println(hdr)

        seek(io, position(io) + hdr.parameters_size)

        # tract data block
        tracts = Vector{UrbanPop.Tract}(undef, hdr.n_tract)
        read!(io, tracts)
        # seek(io, position(io) + hdr.n_tract * sizeof(UrbanPop.Tract))

        # index cases and counties block
        # index_counties = Vector{UInt32}(undef, hdr.n_hub)
        # index_cases = Vector{UInt32}(undef, hdr.n_hub)
        # read!(io, index_counties)
        # read!(io, index_cases)
        seek(io, position(io) + hdr.n_hub * (sizeof(UInt16) + sizeof(CountyCases)))

        # flight_matrix block
        seek(io, position(io) + 57 * 57 * sizeof(UInt32))

        # variant prevalence block
        seek(io, position(io) + 57 * sizeof(Float64))

        wrk_sch = Vector{UInt8}(undef, 14)
        sch_sch = Vector{UInt8}(undef, 14)

        read!(io, wrk_sch)
        read!(io, sch_sch)

        println("work_schedule = ", String(wrk_sch), ", school_schedule = ", String(sch_sch))

        policy_str_len = read(io, UInt64)
        policy_str = Vector{UInt8}(undef, policy_str_len)

        read!(io, policy_str)

        println("poilcies: \n", String(policy_str))

        # RNG state block
        rng_state = Vector{UInt64}(undef, hdr.n_processor)
        read!(io, rng_state)

        foreach(rng_state) do x
            println(Int.(reinterpret(NTuple{4,UInt16}, x)))
        end

        # cell_offsets = Vector{UInt64}(undef, hdr.n_comm);
        # read!(io, cell_offsets)
        # # return Int.(cell_offsets)

        # first_block = cell_offsets[1]
        first_block = read(io, UInt64)
        
        @show(Int(first_block))

        seek(io, first_block)

        data = Vector{Community{C,P}}(undef, hdr.n_comm)
        for k = 1:hdr.n_comm
            n_agent = read(io, UInt64)
            if n_agent > 10000
                @show(Int(n_agent))
                continue
            end
            data[k] = Community(C, P, n_agent)
            ref = Ref(data[k].cell_data)
            read!(io, ref)

            data[k].cell_data = ref[]
            read!(io, data[k].particles)
        end

        return data, tracts
    end
end
# ============================================================================ #
function read_epicast_tract_file(tr_file::AbstractString, ::Val{true})
    tr, mrgn = read_epicast_tract_file(tr_file, Val(false))
    sort!(tr)
    return tr, mrgn
end
# ============================================================================ #
function read_epicast_tract_file(tr_file::AbstractString, ::Val{false})
    data = open(tr_file, "r") do io
        n_tract = tryparse(Int, readline(io))
        @assert(n_tract != nothing, "failed to parse # of tracts")
        data = readdlm(io, Int)
        @assert(size(data, 1) == n_tract, "error reading data, wrong # of tracts")
        return data
    end

    tr = Vector{UrbanPop.Tract}(undef, size(data, 1))
    mrgn = Dict{Int, UrbanPop.TractMarginals}()
    
    for k in 1:size(data, 1)
        fips = data[k,4] * 10^6 + data[k,5]
        tr[k] = UrbanPop.Tract(
            0,
            data[k,2],
            fips
        )

        tmp = UrbanPop.TractMarginals(data[k,2], 
            fips)
        tmp.age[1:5] .= data[k, 6:10]
        tmp.age[6] = sum(tmp.age[1:5])
        tmp.household_size[1:7] .= (data[k, 11:17] .* (1:7))
        tmp.household_size[20] = sum(tmp.household_size[1:19])

        mrgn[fips] = tmp
    end

    return tr, mrgn
end
# ============================================================================ #
function community_counts(wf_file::AbstractString, tr::Vector{UrbanPop.Tract})
    
    wf = Workerflow.read_workerflow_file(wf_file)

    n_res = zeros(UInt8, length(tr))
    n_wrk = zeros(UInt8, length(tr))

    for k in eachindex(tr)
        fips = Int(tr[k].fips_code)
        col = Workerflow.get_column(wf, fips)
        n_res[k] = max(1, round(UInt8, tr[k].n_agent / N_AGENT_PER_COMMUNITY))

        n_wrk[k] = (isempty(col) || !haskey(col, fips)) ? 0x00 :
            round(UInt8, Workerflow.vec_sum(col, fips) / N_AGENT_PER_COMMUNITY)
    end

    return n_res, n_wrk
end
# ============================================================================ #
function stats(c::CellData, t::UrbanPop.Tract)
    n_agent = Int(c.nc_0[end])
    if c.n_family == 0
        # worker only community
        return 0x00, 0x01, UrbanPop.TractMarginals(n_agent, t.fips_code)
    else
        return 0x01, 0x00, UrbanPop.TractMarginals(Int[c.nc_0...], Int[c.nc_0_HH...],
            n_agent, t.fips_code)
    end
end
# ---------------------------------------------------------------------------- #
function stats(c::CellDataOrig, t::UrbanPop.Tract)
    n_agent = Int(c.nc_0[end])
    if c.n_family == 0
        # worker only community
        return 0x00, 0x01, UrbanPop.TractMarginals(n_agent, t.fips_code)
    else
        return 0x01, 0x00, UrbanPop.TractMarginals(Int[c.nc_0...], zeros(Int, 20),
            n_agent, t.fips_code)
    end
end
# ============================================================================ #
function checkpoint_community_counts(ck_file::AbstractString,
    tr::Vector{UrbanPop.Tract}, ::Val{B}) where {B}

    C, P = B ? (CellData, Particle) : (CellDataOrig, ParticleOrig)
    ck = read_checkpoint_file(C, P, ck_file)

    out = Dict{Int,Tuple{Int,Int,UrbanPop.TractMarginals}}()

    for k in eachindex(ck)
        if ck[k].cell_data.tract != CELL_UNUSED
            idx = Int(ck[k].cell_data.tract)
            @assert(0 <= idx < length(tr), "Invalid tract idx: $(idx)")
            n_res, n_wrk, mrgn = stats(ck[k].cell_data, tr[idx + 1])
            tract = tr[idx + 1].fips_code
            if !haskey(out, tract)
                out[tract] = (n_res, n_wrk, mrgn)
            else
                out[tract] = (out[tract][1] + n_res, out[tract][2] + n_wrk,
                    out[tract][3] + mrgn)
            end
        end
    end
    return out
end
# ============================================================================ #
function print_difference(io::IO, s::Set)
    if !isempty(s)
        println(io, collect(s))
    end
end
# ============================================================================ #
function compare_checkpoints(ck1::AbstractString, ck2::AbstractString,
    tr1::AbstractString, tr2::AbstractString, ::Val{B1}, ::Val{B2}) where {B1,B2}

    t1, _ = read_epicast_tract_file(tr1, Val(B1))
    t2, _ = read_epicast_tract_file(tr2, Val(B2))

    c1 = checkpoint_community_counts(ck1, t1, Val(B1))
    c2 = checkpoint_community_counts(ck2, t2, Val(B2))

    keys_use = intersect(keys(c1), keys(c2))

    @info("# of tracts in common = $(length(keys_use))")

    n = 0

    ofile = joinpath(@__DIR__, "..", "..",
        Dates.format(Dates.now(), "YYYYmmdd_HHMMSS") * "-compare_checkpoint.log")

    io = open(ofile, "w")

    println(io, "# checkpoint 1: \"$(ck1)\"")
    println(io, "# checkpoint 2: \"$(ck2)\"")

    for k in keys_use
        if c1[k] != c2[k]
            println(io, k)
            if c1[k][1] != c2[k][1] || c1[k][2] != c2[k][2]
                println(io, "    c1: ", c1[k][1:2], ", c2: ", c2[k][1:2])
            end
            println(io, "    c1: ", c1[k][3], "\n    c2: ", c2[k][3])
            n += 1
        end
    end

    println(io, "# in checkpoint 1 but not 2:")
    print_difference(io, setdiff(keys(c1), keys(c2)))
    println(io, "# in checkpoint 2 but not 1:")
    print_difference(io, setdiff(keys(c2), keys(c1)))

    close(io)

    return n
end
# ============================================================================ #
function validate_checkpoint(ck_file::AbstractString, tr_file::AbstractString,
    wf_file::AbstractString, ag_file::AbstractString="", ::Val{B}=Val{true}) where {B}

    if !isempty(ag_file)
        tr = UrbanPop.read_tract_file(tr_file)
        mrgn = UrbanPop.tract_marginals(ag_file)
    else
        tr, mrgn = read_epicast_tract_file(tr_file, Val(B))
    end

    n_res, n_wrk = community_counts(wf_file, tr)

    if B
        ck = checkpoint_community_counts(CellData, Particle, ck_file, tr)
    else
        ck = checkpoint_community_counts(CellDataOrig, ParticleOrig, ck_file, tr)
    end

    keys_use = intersect(keys(ck), keys(mrgn))
    if length(ck) != length(mrgn)
        @warn("Mismatch in # of tracts, $(length(ck)) vs. $(length(mrgn))")
    elseif keys(ck) != keys(mrgn)
        @warn("Mismatch in tract FIPS codes")
    end

    fips = getfield.(tr, :fips_code)

    n = 0

    ofile = joinpath(@__DIR__, "..", "..",
        Dates.format(Dates.now(), "YYYYmmdd_HHMMSS") * "-validate_checkpoint.log")

    io = open(ofile, "w")

    for k in keys_use
        idx = findfirst(isequal(k), fips)

        tmp = (n_res[idx], n_wrk[idx], mrgn[k])
        
        if ck[k] != tmp
            println(io, k)
            if ck[k][1] != n_res[idx] || ck[k][2] != n_wrk[idx]
                println(io, "    (n_res, n_wrk) - ck: ", ck[k][1:2], ", db: ",
                    Int[n_res[idx], n_wrk[idx]])
            end
            println(io, "    ck: ", ck[k][3], "\n    db: ", mrgn[k])
            n += 1
        end
    end

    println(io, "# in checkpoint but not database:")
    print_difference(io, setdiff(keys(ck), keys(mrgn)))
    println(io, "# in database but not checkpoint:")
    print_difference(io, setdiff(keys(mrgn), keys(ck)))

    close(io)

    if n == 0
        @info("SUCCESS - it appears that the test passed...")
    end

    return n
end
# ============================================================================ #
function count_daygroups(c::Community{CellData,Particle})

    sg_count = zeros(Int, 6)
    sg_ngrps = zeros(Int, 6)
    sg_grps = Set{Int}()
    wg_count = zeros(Int, 1000)
    wg_ngrps = zeros(Int, 1000)
    wg_grps = Set{Int}()

    for pt in c.particles
        school = pt.school
        naics = pt.naics_code
        if school > 0
            sg_count[school] += 1
            if !in(pt.daygroup, sg_grps)
                sg_ngrps[school] += 1
                push!(sg_grps, pt.daygroup)
            end
        elseif naics > 0
            pt.employment < 1 && @warn("Agent $(pt.agent_id) naics > 0 by unemployed?")
            wg_count[naics] += 1
            if !in(pt.daygroup, wg_grps)
                wg_ngrps[naics] += 1
                push!(wg_grps, pt.daygroup)
            end
        end
    end

    return sg_count ./ sg_ngrps, wg_count ./ wg_ngrps
end
# ============================================================================ #
function validate_daygroups(cp_file::AbstractString, tr_file::AbstractString)

    tracts = UrbanPop.read_tract_file(tr_file)

    data = read_checkpoint_file(cp_file)



end
# ============================================================================ #
is_student(x) = x.naics_code != 611 && x.school != 0
is_teacher(pt) = pt.naics_code == 611
function school_stats(cp_file::AbstractString, tr_file::AbstractString)

    data = read_checkpoint_file(cp_file)
    tr = UrbanPop.read_tract_file(tr_file)
    res = Dict{Int,Dict{Int,Vector{Int}}}()

    for j in eachindex(data)
        fips = div(tr[data[j].cell_data.tract+1].fips_code,10^6)
        if !haskey(res, fips)
            res[fips] = Dict{Int,Vector{Int}}()
        end
        for k in eachindex(data[j].particles)
            if is_student(data[j].particles[k])
                id = data[j].particles[k].school
                field = min(abs(id),6)
                tmp = get(res[fips], field, zeros(Int,2))
                if id > 0
                    tmp[1] += 1
                else
                    tmp[2] += 1
                end
                res[fips][field] = tmp
            end
        end
    end

    out =  Dict{Int,Dict{Int,Float64}}()
    for (cnty, d) in res
        out[cnty] = Dict( id => v[1] / (v[1]+v[2]) for (id,v) in d)
    end

    return out
end
# ============================================================================ #
is_working(x) = x.employment > 0 && (x.work == floor.(Int, x.r))
function work_stats(cp_file::AbstractString)
    data = read_checkpoint_file(cp_file)
    out = Dict{Int,NamedTuple{(:w,:n),Tuple{Int,Int}}}()

    for j in eachindex(data)
        for k in eachindex(data[j].particles)
            nc = div(data[j].particles[k].naics_code, 10)
            if nc > 0 && is_working(data[j].particles[k])
                tmp = get(out, nc, (w=0,n=0))
                out[nc] = (w=tmp.w + 1, n=tmp.n)
            elseif nc > 0
                tmp = get(out, nc, (w=0,n=0))
                out[nc] = (w=tmp.w, n=tmp.n+1)
            end
        end
    end

    return Dict(k => v.w / (v.w + v.n) for (k,v) in out)
end
# ============================================================================ #
end
#=
TODO:
    * [x] compare age and household distributions from DB vs CK (at the tract level?)
    * [ ] could also compare race, ethnicity, etc. distributions for v2
=#