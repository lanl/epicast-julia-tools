module EpicastPlot

using Colors, PyPlot
using Epicast, EpicastTables
# ============================================================================ #
function plot_all_multi(data::Epicast.RunData, odir::AbstractString;
    ext::AbstractString="png")

    PyPlot.ioff()

    h, ax = plot_run(data, "total", freduce=mean_new_cases, normalize=true,
        ylab="Normalized new cases per tract")

    ax.get_legend().remove()
    h.tight_layout()

    ext = startswith(ext, '.') ? ext[2:end] : ext

    h.savefig(joinpath(odir, "total-cases_$(data.runno)." * ext), dpi=200)

    PyPlot.close(h)

    if Epicast.has_data(data, "age_0")

        h, ax = subplots(1, 4)
        h.set_size_inches((16, 5))

        plot_run(data, "age", freduce=mean_new_cases, normalize=true,
            title="Age", ylab="Normalized new cases per tract", h=h, ax=ax[1])

        plot_run(data, "household", freduce=mean_new_cases, normalize=true, dropname=true,
            title="Household size", ylab="", h=h, ax=ax[2])

        plot_run(data, "race", freduce=mean_new_cases, normalize=true, dropname=true,
            title="Race", ylab="", h=h, ax=ax[3])

        plot_run(data, "ethnicity", freduce=mean_new_cases, normalize=true, dropname=true,
            title="Ethnicity", ylab="", h=h, ax=ax[4])

        h.tight_layout()

        h.savefig(joinpath(odir, "demographic-cases_$(data.runno)." * ext), dpi=200)

        PyPlot.close(h)
    end

    if Epicast.has_data(data, "hospitalized_age0")

        h, ax = subplots(1,4)
        h.set_size_inches((16, 5))

        plot_run(data, "hospitalized", freduce=mean_new_cases, normalize=false, dropname=true,
            title="Hospitalization", ylab="Daily incidence per tract", h=h, ax=ax[1])

        plot_run(data, "icu", freduce=mean_new_cases, normalize=false, dropname=true,
            title="ICU admitance", ylab="", h=h, ax=ax[2])

        plot_run(data, "ventilated", freduce=mean_new_cases, normalize=false, dropname=true,
            title="Ventilation", ylab="", h=h, ax=ax[3])

        plot_run(data, "dead", freduce=mean_new_cases, normalize=false, dropname=true,
            title="Deaths", ylab="", h=h, ax=ax[4])

        h.tight_layout()

        h.savefig(joinpath(odir, "treatment-status_$(data.runno)." * ext), dpi=200)

        PyPlot.close(h)
    end

    if Epicast.has_data(data, "infection-src_family")

        h, ax = plot_run(data, "infection-src", freduce=mean_new_cases, normalize=false, dropname=true,
            title="Infection context", ylab="Daily cases per tract")

        h.tight_layout()

        h.savefig(joinpath(odir, "infection-context_$(data.runno)." * ext), dpi=200)

        PyPlot.close(h)
    end

    if Epicast.has_data(data, "status_incubation")
        h, ax = plot_run(data, "status", freduce=total_cases, normalize=false, dropname=true,
            title="Infection status", ylab="# of agents")

        h.tight_layout()

        h.savefig(joinpath(odir, "infection-status_$(data.runno)." * ext), dpi=200)

        PyPlot.close(h)

    end

    PyPlot.ion()

    return nothing
end
# ============================================================================ #
function plot_all(data::Epicast.RunData, ofile::AbstractString; plot_status::Bool=false)
    PyPlot.ioff()
    h, ax = subplots(3, 3)
    h.set_size_inches((20,18))

    plot_run(data, "age", freduce=mean_new_cases, normalize=true,
        title="Age", ylab="Daily cases per tract", h=h, ax=ax[1,1])

    plot_run(data, "household", freduce=mean_new_cases, normalize=true, dropname=true,
        title="Household size", ylab="Daily cases per tract", h=h, ax=ax[1,2])

    plot_run(data, "race", freduce=mean_new_cases, normalize=true, dropname=true,
        title="Race", ylab="Daily cases per tract", h=h, ax=ax[1,3])

    plot_run(data, "ethnicity", freduce=mean_new_cases, normalize=true, dropname=true,
        title="Ethnicity", ylab="Daily cases per tract", h=h, ax=ax[2,1])

    plot_run(data, "hospitalized", freduce=mean_new_cases, normalize=false, dropname=true,
        title="Hospitalization", ylab="Daily incidence per tract", h=h, ax=ax[2,2])

    plot_run(data, "icu", freduce=mean_new_cases, normalize=false, dropname=true,
        title="ICU admitance", ylab="Daily incidence per tract", h=h, ax=ax[2,3])

    plot_run(data, "ventilated", freduce=mean_new_cases, normalize=false, dropname=true,
        title="Ventilation", ylab="Daily incidence per tract", h=h, ax=ax[3,1])

    plot_run(data, "dead", freduce=mean_new_cases, normalize=false, dropname=true,
        title="Deaths", ylab="Daily incidence per tract", h=h, ax=ax[3,2])

    plot_run(data, "infection-src", freduce=mean_new_cases, normalize=false, dropname=true,
        title="Infection context", ylab="Daily cases per tract", h=h, ax=ax[3,3])

    h.tight_layout()

    h.savefig(ofile, dpi=200)

    PyPlot.close(h)

    if plot_status
        h, ax = plot_run(data, "status", freduce=total_cases, normalize=false, dropname=true,
            title="Agent status", ylab="# of agents")

        h.tight_layout()

        ext = splitext(ofile)[2]
        h.savefig(replace(ofile, ext => "_agent-status" * ext), dpi=200)
        PyPlot.close(h)
    end

    PyPlot.ion()

    return nothing
end
# ============================================================================ #
function plot_states(data::Epicast.RunData, ofile::AbstractString; title::AbstractString="")

    # 10^9 -> FIPS tract to state conversion
    dat, lab = Epicast.group_by(data, "total", EpicastTables.TRACT2STATE, case_count!)

    PyPlot.ioff()
    h, ax = subplots(1,1)
    h.set_size_inches((14,9))

    colors = map(col -> (red(col), green(col), blue(col)),
        distinguishable_colors(length(lab), [RGB(1,1,1), RGB(0,0,0)],
            dropseed=true)
    )

    foreach(1:length(lab)) do k
        ax.plot(dat[:,1,k] .* 1e5, color=colors[k], label=STATE_FIPS[lab[k]])
    end

    ax.legend(frameon=false)

    ax.set_xlabel("Simulation day", fontsize=14)
    ax.set_ylabel("Daily infections per 100k residents", fontsize=14)

    ax.spines["top"].set_visible(false)
    ax.spines["right"].set_visible(false)

    # "Daily new cases per state (22 states, ~136M agents)"
    !isempty(title) && ax.set_title(title, fontsize=18)

    h.tight_layout()

    h.savefig(ofile, dpi=200)

    close(h)

    PyPlot.ion()

    return nothing
end
# ============================================================================ #
function plot_run(data::Vector{Epicast.RunData}, name::AbstractString; args...)

    h, ax = (nothing,nothing)

    for d in data
        h, ax = plot_run(d, name; h=h, ax=ax, args...)
    end

    return h, ax
end
# ============================================================================ #
function plot_runs(data, name::AbstractString; freduce=total_cases, reduce_cols=nothing,
    normalize::Bool=false, demo::AbstractString="", dropname::Bool=false, ylab::String="",
    title::String="", h=nothing, ax=nothing, style::Function=(h, ax) -> nothing,
    dataset_names::AbstractVector{<:AbstractString}=nothing)

    cols = map(d -> filter_vars(x -> startswith(x, name), d.data), data)
    if h == nothing || ax == nothing
        h, ax = subplots(1, 1)
        h.set_size_inches((8,6))
    end

    t = 1:maximum(map(Epicast.n_timepoint, data))

    rm = dropname ? name * "_" => "" : "_" => " "

    n_colors = sum(map(length, cols))
    if reduce_cols != nothing
        n_colors = sum(map(x -> 1, cols))
    end
    colors = map(col -> (red(col), green(col), blue(col)),
                 distinguishable_colors(n_colors, [RGB(1,1,1), RGB(0,0,0)], dropseed=true))

    plot_col(k, dataset, dat, name, col) = begin
        d = demo
        if 0 == cmp(demo, "")
            d = col
        end
        if normalize && Epicast.has_demographic(dataset, d)
            tmp2 = sum(Epicast.demographics(dataset, d))
            dat ./= tmp2
            #dat .*= 1e5

            replace!(x -> isnan(x) || isinf(x) ? 0.0 : x, dat)
        end
        tmp = freduce(dat)
        label = replace(col, rm)
        if dataset_names != nothing
            if length(colors) != length(data)
                label = "$(name)_$label"
            else
                label = name
            end
        end
        ax.plot(t, tmp, label=label, color=colors[k])
    end

    k = 1
    for (j, dataset_cols) in enumerate(cols)
        dataset = data[j]
        n = dataset_names[j]
        if reduce_cols == nothing
            for (c, col) in enumerate(dataset_cols)
                d = Float64.(rundata(dataset, col))
                plot_col(k, dataset, d, n, col)
                k += 1
            end
        else
            d = reduce_cols(map(col -> Float64.(rundata(dataset, col)), dataset_cols))
            plot_col(k, dataset, d, n, name)
            k += 1
        end
    end

    ax.spines["top"].set_visible(false)
    ax.spines["right"].set_visible(false)
    ax.set_xlabel("Simulation day", fontsize=14)

    !isempty(ylab) && ax.set_ylabel(ylab, fontsize=14)
    !isempty(title) && ax.set_title(title, fontsize=18)

    ax.legend(fontsize=14, frameon=false)
    style(h, ax)

    h.tight_layout()

    return h, ax
end
# ============================================================================ #
function plot_run(data::Epicast.RunData, name::AbstractString; freduce=total_cases,
    normalize::Bool=false, dropname::Bool=false, ylab::String="",
    title::String="", h=nothing, ax=nothing, style::Function=(h, ax) -> nothing,
    agg_level::Type{<:AbstractGeo}=State, smooth::Bool=false)

    cols = filter_vars(x -> startswith(x, name), data.data)
    if h == nothing || ax == nothing
        h, ax = subplots(1, 1)
        h.set_size_inches((8,6))
    end

    t = 1:Epicast.n_timepoint(data)

    rm = dropname ? name * "_" => "" : "_" => " "

    colors = map(col -> (red(col), green(col), blue(col)),
        distinguishable_colors(length(cols), [RGB(1,1,1), RGB(0,0,0)],
            dropseed=true)
    )
    tmp = EpicastTables.aggregate(agg_level, data.data)
    for (k, col) in enumerate(cols)
        if smooth
            EpicastTables.smooth!(tmp, col, 7, mean)
        end
        cases = freduce(tmp[col])
        if normalize && has_demographic(data, name)
            pop = EpicastTables.aggregate(agg_level, data.demog)[name]
            cases ./= (pop ./ 1e5)
            replace!(x -> isnan(x) || isinf(x) ? 0.0 : x, cases)
        end
        ax.plot(t, cases, label=replace(col, rm), color=colors[k])
    end

    ax.spines["top"].set_visible(false)
    ax.spines["right"].set_visible(false)
    ax.set_xlabel("Simulation day", fontsize=14)

    !isempty(ylab) && ax.set_ylabel(ylab, fontsize=14)
    !isempty(title) && ax.set_title(title, fontsize=18)

    ax.legend(fontsize=14, frameon=false)
    style(h, ax)

    h.tight_layout()

    return h, ax
end
# ============================================================================ #
function plot_infection_src(data::Epicast.RunData)
    
    fields = ["family","work","school","neighborhood-cluster","neighborhood"]
    total = Epicast.new_cases(Epicast.rundata(data, "total"))
    
    d = zeros(Float64, Epicast.n_timepoint(data), length(fields))
    for k in eachindex(fields)
        d[:,k] .= Epicast.new_cases(Epicast.rundata(data, "infection-src_" * fields[k]))
    end
   
    d ./= total

    h, ax = subplots(1,1)

    h.set_size_inches((10,6))

    ax.plot(d, linewidth=2.5)

    ax[:spines]["right"].set_visible(false)
    ax[:spines]["top"].set_visible(false)

    ax[:spines]["left"].set_linewidth(2.5)
    ax[:spines]["bottom"].set_linewidth(2.5)

    ax.xaxis.set_tick_params(width=2.5)
    ax.yaxis.set_tick_params(width=2.5)

    ax.legend(fields, fontsize=12, frameon=false, loc="upper left",
        bbox_to_anchor=(0.9,1.0))

    ax.set_ylabel("Proportion of total cases", fontsize=14)
    ax.set_xlabel("Simulation time (days)", fontsize=14)

    h.tight_layout()

    return h, ax
end
# ============================================================================ #
end # module EpicastPlot
