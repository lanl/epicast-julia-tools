module NYTCovidData

using CSV, Dates, Random

using UrbanPop

# ============================================================================ #
# up_dir = "/Users/palexander/Documents/emerge+radium/input-data/agent_db"
function create_db(idir::AbstractString, up_dir::AbstractString,
    ofile::AbstractString="")

    tract_fips, tract_pop = UrbanPop.all_tract_data(up_dir)
    
    county_data = Dict{Int,Int}()
    for k in eachindex(tract_fips)
        cnty = div(tract_fips[k], 10^6)
        county_data[cnty] = tract_pop[k] + get(county_data, cnty, 0)
    end

    # return county_data
    cnty = sort!(unique(div.(tract_fips, 10^6)))
    upc = Dict{Int,Int}(v => k for (k,v) in enumerate(cnty))

    data = Vector{Matrix{Datum}}(undef, 4)
    years = ["2020", "2021", "2022", "2023"]
    for (k, yr) in enumerate(years)
        data[k], ref = read_file_data(joinpath(idir, "us-counties-$(yr).csv"), upc, county_data)
    end
    
    return vcat(data...), cnty, Dates.Date("2020-01-21")
end
# ============================================================================ #
struct Datum
    cases::UInt32
    deaths::UInt32
end
Datum(x::Integer,::Missing) = Datum(x, 0xffffffff)
Datum(::Missing,::Missing) = Datum(0xffffffff, 0xffffffff)
Datum(;args...) = Datum(get(args, :cases, 0), get(args, :deaths, 0))
# ---------------------------------------------------------------------------- #
Base.length(::Datum) = 1
Base.broadcastable(d::Datum) = Ref(d)
# ---------------------------------------------------------------------------- #
sat_add(a::UInt32, b::UInt32) = a == 0xffffffff || b == 0xffffffff ? 0xffffffff : a + b
# ---------------------------------------------------------------------------- #
function Base.:+(a::Datum, b::Datum)
    return Datum(sat_add(a.cases, b.cases), sat_add(a.deaths, b.deaths))
end
# ---------------------------------------------------------------------------- #
Base.zero(::Type{Datum}) = Datum(0,0)
# ---------------------------------------------------------------------------- #
is_null(x::Datum) = x.cases == 0 && x.deaths == 0
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

    for k in 1:length(csv)
        j = (csv.date[k] - ref).value + 1
        dat = Datum(csv.cases[k], csv.deaths[k])
        
        if haskey(cnty, csv.fips[k])
            data[j,cnty[csv.fips[k]]] = dat

        elseif csv.county[k] == "Unknown" && haskey(fips_map, csv.state[k]) &&
            !is_null(dat)

            state = fips_map[csv.state[k]]
            unknown[j,state] = unknown[j,state] + dat
        end
    end

    cnty_fips = collect(keys(cnty))

    for k = 1:57
        if !is_null(sum(view(unknown, :, k)))
            # this state has unallocated cases/deaths
            cntys = filter(x->div(x,10^3)==k, cnty_fips)
            cols = Int[cnty[x] for x in cntys]
            c_pop = Int[pop[x] for x in cntys]
            for j = 1:n_day
                tmp = view(data, j, cols)
                if unknown[j,k].cases > 0
                    allocate_unknown!(tmp, unknown[j,k], c_pop, :cases)
                end
                if unknown[j,k].deaths > 0
                    allocate_unknown!(tmp, unknown[j,k], c_pop, :deaths)
                end
            end

        end
    end

    # NOTE TODO FIXME WARNING ERROR
    # handel New York City, Kansas City, MO, etc.

    return data, ref
end
# ============================================================================ #
function allocate_unknown!(data::AbstractVector{Datum}, un::Datum,
    pop::Vector{<:Integer}, field::Symbol)

    idx = findall(isequal(0), getfield.(data, field))

    x = getfield(un, field)

    if x > length(idx)
        pops = Vector{Float64}(pop[idx])
        pops ./= sum(pops)

        for i in eachindex(idx)
            d = round(UInt32, x * pops[i])
            data[idx[i]] += Datum(;field => d)
        end
    else
        tidx = shuffle!(idx)[1:x]
        data[tidx] .+= Datum(;field => 1)
    end

    return data
end
# ============================================================================ #
end
