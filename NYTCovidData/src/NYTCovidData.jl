module NYTCovidData

using CSV, Dates, Random

using UrbanPop

# ============================================================================ #
# up_dir = "/Users/palexander/Documents/emerge+radium/input-data/agent_db"
function create_db(idir::AbstractString, up_dir::AbstractString,
    ofile::AbstractString)

    tract_fips, tract_pop = UrbanPop.all_tract_data(up_dir)
    
    county_data = Dict{Int,Int}()
    for k in eachindex(tract_fips)
        cnty = div(tract_fips[k], 10^6)
        county_data[cnty] = tract_pop[k] + get(county_data, cnty, 0)
    end

    cnty = sort!(unique(div.(tract_fips, 10^6)))
    upc = Dict{Int,Int}(v => k for (k,v) in enumerate(cnty))

    # return upc, county_data

    years = ["2020", "2021", "2022", "2023"]
    data = Vector{Matrix{Datum}}(undef, length(years))

    for (k, yr) in enumerate(years)
        data[k], ref = read_file_data(joinpath(idir, "us-counties-$(yr).csv"), upc, county_data)
    end

    all_data = permutedims(vcat(data...), (2,1))

    open(ofile, "w") do io

        write(io, UInt64(size(all_data, 1))) # # of counties
        write(io, UInt64(size(all_data, 2))) # # of time points        
        write(io, Vector{UInt8}("2020-01-21\0")) # reference date
        write(io, UInt16.(cnty)) # county FIPS as UInt16 in order
        
        write(io, all_data)
    end

    return ofile
    
    # return vcat(data...), cnty, Dates.Date("2020-01-21")
end
# ============================================================================ #
struct Datum
    cases::UInt32
    deaths::UInt32
end
Datum(x::Integer,::Missing) = Datum(x, 0)
Datum(::Missing,::Missing) = Datum(0, 0)
Datum(;args...) = Datum(get(args, :cases, 0), get(args, :deaths, 0))
# ---------------------------------------------------------------------------- #
Base.length(::Datum) = 1
Base.broadcastable(d::Datum) = Ref(d)
# ---------------------------------------------------------------------------- #
# sat_add(a::UInt32, b::UInt32) = a == 0xffffffff || b == 0xffffffff ? 0xffffffff : a + b
# ---------------------------------------------------------------------------- #
function Base.:+(a::Datum, b::Datum)
    # return Datum(sat_add(a.cases, b.cases), sat_add(a.deaths, b.deaths))
    return Datum(a.cases + b.cases, a.deaths + b.deaths)
end
# ---------------------------------------------------------------------------- #
Base.zero(::Type{Datum}) = Datum(0,0)
# ---------------------------------------------------------------------------- #
Base.iszero(x::Datum) = x.cases == 0 && x.deaths == 0
# ============================================================================ #
function read_file_data(ifile::AbstractString, cnty::Dict{Int,Int}, pop::Dict{Int,Int})

    csv = CSV.File(ifile)

    fips_map = Dict{String,Int}()
    k = 1
    while length(fips_map) < 51
        if !ismissing(csv.fips[k]) && !haskey(fips_map, csv.state[k]) &&
            div(csv.fips[k], 10^3) < 57
            
            fips_map[csv.state[k]] = div(csv.fips[k], 10^3)
        end
        k += 1
    end

    ref = csv.date[1]
    n_day = (csv.date[end] - ref).value + 1

    data = zeros(Datum, n_day, length(cnty))
    unknown = zeros(Datum, n_day, 57)

    t0 = time()

    for k in 1:length(csv)
        j = (csv.date[k] - ref).value + 1
        dat = Datum(csv.cases[k], csv.deaths[k])
        
        if haskey(cnty, csv.fips[k])
            data[j,cnty[csv.fips[k]]] = dat

        elseif csv.county[k] == "Unknown" && haskey(fips_map, csv.state[k]) &&
            !iszero(dat)

            state = fips_map[csv.state[k]]
            unknown[j,state] = unknown[j,state] + dat
        end
    end
    
    t1 = time()
    println("Step 1: $(t1 - t0)")

    cnty_fips = collect(keys(cnty))

    # Main.@infiltrate

    for k = 1:57
        if !iszero(sum(view(unknown, :, k)))
            # this state has unallocated cases/deaths
            cntys = filter(x->div(x,10^3)==k, cnty_fips)
            cols = Int[cnty[x] for x in cntys]
            c_pop = Float64[pop[x] for x in cntys]
            for j = 1:n_day
                tmp = view(data, j:size(data, 1), cols)
                if unknown[j,k].cases > 0
                    allocate_unknown!(tmp, unknown[j,k], c_pop, :cases, Datum(1,0))
                end
                if unknown[j,k].deaths > 0
                    allocate_unknown!(tmp, unknown[j,k], c_pop, :deaths, Datum(0,1))
                end
            end

        end
    end

    t2 = time()
    println("Step 2: $(t2 - t1)")

    # ===== New York City
    # allocate cases assigned to "New York City" to it's 5 constituent counties
    # by population (assuming identical case prevalence in each county)
    nyc_fips = [36005, 36047, 36061, 36081, 36085]
    allocate_many!(data, csv, cnty, pop, "New York City", "New York", nyc_fips, ref)

    # ===== Missouri
    # allocate all cases assigned to "Kansas City" MO to Jackson county, as most
    # of Kansas City falls w/in Jackson
    allocate_single!(data, csv, "Kansas City", "Missouri", cnty[29095], ref)

    # allocate all cases assigned to "Joplin" MO to Jasper county, as most
    # of Joplin falls w/in Jasper
    allocate_single!(data, csv, "Joplin", "Missouri", cnty[29097], ref)

    # ===== Alaska
    # "Bristol Bay plus Lake and Peninsula" -> Bristol Bay, Lake and Peninsula
    ak_fips1 = [2060, 2164]
    allocate_many!(data, csv, cnty, pop, "Bristol Bay plus Lake and Peninsula",
        "Alaska", ak_fips1, ref)

    # "Yakutat plus Hoonah-Angoon" -> Yakutat, Hoonah-Angoon
    ak_fips2 = [2282, 2105]
    allocate_many!(data, csv, cnty, pop, "Yakutat plus Hoonah-Angoon", "Alaska",
        ak_fips2, ref)

    # NOTE:
    # prior to 2019-01-02 Chugach and Copper River Census Areas were part of the
    # "Valdez-Cordova Census Area" - FIPS code 2261, which the NY Times data 
    # uses in lieu of the newer, split census areas. 2261 is present and
    # correctly handled by UrbanPop, so no need for any reassignment

    t3 = time()
    println("Step 3: $(t3 - t2)")

    return data, ref
end
# ============================================================================ #
function allocate_many!(data::Matrix{Datum}, csv::CSV.File, cnty::Dict{Int,Int},
    pop::Dict{Int,Int}, name::AbstractString, state::AbstractString, 
    fips::Vector{<:Integer}, ref::Dates.Date)
    
    data_idx = Int[cnty[x] for x in fips]
    pops = Float64[pop[x] for x in fips]
    pops ./= sum(pops)

    return allocate_non_county!(data, csv, name, state, data_idx, pops, ref)
end
# ============================================================================ #
function allocate_single!(data::Matrix{Datum}, csv::CSV.File, name::AbstractString,
     state::AbstractString, data_idx::Int, ref::Dates.Date)

    return allocate_non_county!(data, csv, name, state, [data_idx], [1.0], ref)
end
# ============================================================================ #
function allocate_non_county!(data::Matrix{Datum}, csv::CSV.File, name::AbstractString,
    state::AbstractString, data_idx::Vector{Int}, pop::Vector{Float64}, ref::Dates.Date)

    idx = findall(1:length(csv.county)) do k
        csv.county[k] == name && csv.state[k] == state
    end

    for j in idx
        t = (csv.date[j] - ref).value + 1
        for k in eachindex(data_idx)
            d = Datum(round(UInt64, csv.cases[j] * pop[k]),
                round(UInt64, csv.deaths[j] * pop[k]))
            data[t, data_idx[k]] += d
        end
    end

    return data
end
# ============================================================================ #
function cumsum_i!(x::AbstractVector{<:Number})
    for k in 2:length(x)
        @inbounds x[k] += x[k-1]
    end
    return x
end
# ============================================================================ #
function allocate_unknown!(data::AbstractMatrix{Datum}, un::Datum,
    pop::Vector{Float64}, field::Symbol, tmp::Datum)

    idx = findall(isequal(0), getfield.(view(data, 1, :), field))
    idx = isempty(idx) ? collect(1:size(data, 2)) : idx

    pops = pop[idx]
    pops ./= sum(pops)
    cumsum_i!(pops)

    return allocate_data!(view(data, :, idx), getfield.(un, field), pops, tmp)
end
# ============================================================================ #
function allocate_data!(data::AbstractMatrix{Datum}, N::Integer,
    cdf::Vector{Float64}, tmp::Datum)

    for k in 1:N
        j = inv_sample(cdf)
        @inbounds data[:, j] .+= tmp
    end

    return data
end
# ============================================================================ #
function inv_sample(x::AbstractVector{<:AbstractFloat}, d::Float64=rand())

    for k in eachindex(x)
        @inbounds d < x[k] && return k
    end
    return -1
end
# ============================================================================ #
function read_covid_data(ifile::AbstractString)

    open(ifile, "r") do io
        n_cnty = read(io, UInt64) # # of counties
        n_tp = read(io, UInt64) # # of time points
        date = String(read(io, 11)[1:10]) # reference date

        # fips code for each county
        county_fips = Vector{UInt16}(undef, n_cnty)
        read!(io, county_fips)
        
        data = Matrix{Datum}(undef, n_cnty, n_tp)
        read!(io, data)

        return data, county_fips, date
    end

end
# ============================================================================ #
end
