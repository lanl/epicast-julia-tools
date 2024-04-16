module Epicast

using DelimitedFiles, PyPlot

# ============================================================================ #
function epi_plot(idir::AbstractString, run::Integer=2)

    # # of new symptomatic cases per day total (across tracts)
    n_symp, n_symp_age, n_symp_hh = parse_logfile(idir, run)


    # pop: population of each tract
    # n_infected: # of agents in infected state (inc | pro | pox) per day per tract
    # n_symp: cumulative # of symptoimatic cases per day per tract
    pop, n_infected, n_symp_tract = parse_attackfiles(idir, run)

    h, ax = subplots(1, 4)
    h.set_size_inches((15,5))

    foreach(ax) do cax
        cax.spines["top"].set_visible(false)
        cax.spines["right"].set_visible(false)
        cax.set_xlabel("Simulation day", fontsize=14)
    end

    t = 1:length(n_symp)

    ax[1].plot(t, n_symp, linewidth=3)

    ax[1].set_ylabel("New symptomatic cases", fontsize=14)
    ax[1].set_title("Total new cases", fontsize=16)

    for k in 1:size(n_symp_age, 2)
        ax[2].plot(t, n_symp_age[:,k], linewidth=2, label="age $(k-1)")
    end

    ax[2].set_title("New cases by age", fontsize=16)
    ax[2].legend(frameon=false, fontsize=12)

    for k in 1:size(n_symp_hh, 2)
        lab = k == 1 ? "household size $(k)" : string(k)
        ax[3].plot(t, n_symp_hh[:,k], linewidth=2, label=lab)
    end

    ax[3].set_title("... by household size", fontsize=16)
    ax[3].legend(frameon=false, fontsize=12)


    ax[4].plot(t, n_symp_tract ./ pop', linewidth=1.5)

    ax[4].set_ylabel("Proportion of agents", fontsize=14)
    ax[4].set_title("# symptomatic per tract", fontsize=16)


    h.tight_layout()

    return h, ax

end
# ============================================================================ #
function parse_logfile(idir::AbstractString, run::Integer=2)

    ifile = joinpath(idir, "Log" * string(run))

    d1 = Vector{Vector{Float64}}(undef, 0)
    d2 = Vector{Vector{Float64}}(undef, 0)

    for line in eachline(ifile)
        if startswith(line, r"\s*\d+\.\d\s+")
            raw = map(split(strip(line), r"\s+")) do x
                y = tryparse(Float64, x)
                return y == nothing ? NaN : y
            end
            push!(d1, raw)
        elseif startswith(line, r"S\s+\d+\.\d\s+")
            raw = map(split(line, r"\s+")[2:end]) do x
                y = tryparse(Float64, x)
                return y == nothing ? -1 : y
            end
            push!(d2, raw)
        end
    end

    t1 = hcat(d1...)'
    t2 = hcat(d2...)'

    idx = findall(t1[:,1]) do x
        x > 0 && floor(Int, x) == x
    end

    idx2 = findall(t2[:,1]) do x
        x > 0 && floor(Int, x) == x
    end

    # # new symptomaitc cases per day
    n_symptomatic = t1[idx,2]

    # # new symptomaitc cases per day for each age group
    n_symp_age = t2[idx2,2:6]
    
    # # new symptomaitc cases per day for each household size
    n_symp_hh = t2[idx2,6:end]
    
    return n_symptomatic, n_symp_age, n_symp_hh
end
# ============================================================================ #
function attack_cmp_timepoint(a::AbstractString, b::AbstractString)
    at = tryparse(Int, splitext(a)[2][2:end])
    bt = tryparse(Int, splitext(b)[2][2:end])
    return at < bt
end
# ---------------------------------------------------------------------------- #
function parse_attackfiles(idir::AbstractString, run::Integer=2)
    pat = Regex("Attack" * string(run) * "\\.\\d{3}")
    files = find_files(idir, pat)
    sort!(files, lt=attack_cmp_timepoint)
    
    data = map(parse_attackfile, files)
    
    n_tract = Int(maximum(vcat(getindex.(data, :, 1)...))) + 1

    pop = zeros(n_tract)
    # # agents in infected state: time x tract
    n_infected = fill(0.0, (length(data), n_tract))

    # cumulative # of symptomatic agents: time x tract
    n_symp = fill(0.0, (length(data), n_tract))

    for k in eachindex(data)
        idx = Int.(data[k][:,1]) .+ 1
        pop[idx] .= data[k][:,3]
        n_infected[k,idx] .= sum(data[k][:,5:7], dims=2)
        n_symp[k,idx] .= data[k][:,4]
    end

    return pop, n_infected, n_symp
end
# ---------------------------------------------------------------------------- #
function parse_attackfile(ifile::AbstractString)
    return readdlm(ifile, skipstart=1)
end
# ============================================================================ #
find_files(dir::AbstractString, re=r".*") = return do_match(dir, re, isfile)
# ---------------------------------------------------------------------------- #
find_directories(dir::AbstractString, re=r".*") = return do_match(dir, re, isdir)
# ============================================================================ #
function do_match(dir::AbstractString, re::Regex, f::Function)
    if !isdir(dir)
        error("Input is not a vaild directory path")
    end
    files = [joinpath(dir, x) for x in readdir(dir)]
    return filter(x->occursin(re, x) && f(x), files)
end
# ============================================================================ #
end
