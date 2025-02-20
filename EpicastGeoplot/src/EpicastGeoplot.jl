module EpicastGeoplot

using Shapefile, PyPlot, Colors, Printf, Statistics, Dates

import EpicastTables

using Epicast, EpicastTables

using PyCall

export Polygon, Point, County, Tract, State, BlockGroup

include("./normalize.jl")

const animation = PyNULL()

function __init__()
    copy!(animation, pyimport("matplotlib.animation"))
end

const DATADIR = joinpath(@__DIR__, "..", "assets", "geo-data", "data")

const STATE_OUTLINE_WIDTH = 0.5
# ============================================================================ #
function case_count!(out::AbstractVector, x::Epicast.RunData,
    var::AbstractString, idx::AbstractVector{<:Integer}, norm::AbstractString="")

    in = view(Epicast.rundata(x, var), :, idx)

    # Defaults; can be overriden by passing norm manually
    if norm == ""
        if Epicast.has_demographic(x, var)
            norm = ("new", var)

        elseif endswith(var, r"_age\d")
            # number of new [hospitialized, icu, ventiated, dead] for given age
            # group relative to that age group's total population (for each
            # location)
            norm = "age_" * var[end]

        elseif startswith(var, "infection-src")
            # % of infections attributable to given context
            norm = "new_cases_prop"

        elseif startswith(var, "status_")
            # % of total agents w/in geography that are in given state
            norm = "prop"

        elseif startswith(var, "broadcaster_")
            # % of broadcasters w/in geography that are in given state
            norm = ("total", "media_broadcaster")
        else
            # default: total counts per location
            norm = "default"
        end
    end

    if length(norm) == 2 && norm[1] == "new" && Epicast.has_demographic(x, norm[2])
        # number of new cases relative to total size of that demographic
        out .= Epicast.new_cases(in)
        out ./= sum(view(Epicast.demographics(x, norm[2]), idx))
        out .*= 1e5

    elseif length(norm) == 2 && norm[1] == "total" && Epicast.has_demographic(x, norm[2])
        # number of total cases relative to total size of that demographic
        out .= Epicast.total_cases(in)
        out ./= sum(view(Epicast.demographics(x, norm[2]), idx))
        out .*= 1e5

    elseif endswith(norm, r"age_\d")
        # number of new [hospitialized, icu, ventiated, dead] for given age
        # group relative to that age group's total population (for each
        # location)
        age = "age_" * norm[end]

        out .= Epicast.new_cases(in)
        out ./= sum(view(Epicast.demographics(x, age), idx))
        out .*= 1e5

    elseif norm == "new_cases_prop"
        # % of infections attributable to given context
        out .= Epicast.new_cases(in)
        out ./= Epicast.new_cases(view(Epicast.rundata(x, "total"), :, idx))

    elseif norm == "change_prop"
        # change in state as a % of total agents w/in geography
        out .= Epicast.new_cases(in)
        out ./= sum(view(Epicast.demographics(x, "total"), idx, :))
        out .*= 1e5

    elseif norm == "prop"
        # % of total agents w/in geography that are in given state
        out .= Epicast.total_cases(in)
        out ./= sum(view(Epicast.demographics(x, "total"), idx))
        out .*= 1e5

    elseif norm == "change"
        # change in state as a % of total agents w/in geography
        out .= Epicast.new_cases(in)

    else
        # default: total counts per location
        out .= Epicast.total_cases(in)
    end

    return out
end
# ============================================================================ #
function data_dict(::Type{T}, data::Epicast.RunData, names::Vector{<:AbstractString},
    norms::Dict{String,String}=Dict{String,String}()) where T<:AbstractGeo

    conv = to_geo(T)

    if conv > 1
        # dat is: time x fips x var
        dat, grp = Epicast.group_by(data, names, conv, case_count!, norms)
    else
        # time x fips x var
        dat = Array{Float64,3}(undef, Epicast.n_timepoint(data),
            length(data.fips), length(names))

        # this is pretty slow, and kinda inefficient, but it's nice to reuse
        # case_count!()
        for k in 1:size(dat, 2)
            for j in 1:size(dat, 3)
                case_count!(view(dat, :, k, j), data, names[j], [k], get(norms, names[j], ""))
            end
        end

        grp = data.fips
    end

    replace!(dat, NaN=>0)

    idx = Dict{UInt64,Int}(grp[k] => k for k in eachindex(grp))
    var = Dict{String,Int}(names[k] => k for k in eachindex(names))
    return FIPSTable(T, dat, idx, var)
end
# ============================================================================ #
function area(pts::AbstractVector{Shapefile.Point})
    n = length(pts)
    area = 0.0
    @inbounds for k in eachindex(pts)
        pt1 = pts[k]
        pt2 = pts[(k % n) + 1]
        area += pt1.x * pt2.y - pt1.y * pt2.x
    end
    return abs(area * 0.5)
end
# ============================================================================ #
function largest_part(poly::Shapefile.Polygon)
    ks = poly.parts[end]+1
    ke = length(poly.points)
    mx = area(view(poly.points, ks:ke))

    @inbounds for k = 1:(length(poly.parts)-1)
        ks_tmp = poly.parts[k]+1
        ke_tmp = poly.parts[k+1]
        tmp = area(view(poly.points, ks_tmp:ke_tmp))
        if tmp > mx
            mx = tmp
            ks = ks_tmp
            ke = ke_tmp
        end
    end
    return ks, ke
end
# ============================================================================ #
@inline long2merc(x::Real) = x * (pi/180.0)
@inline lat2merc(y::Real) = log(tan((pi/4.0) + (y * (pi/180) / 2.0)))
@inline point2vec(p::Shapefile.Point) = [long2merc(p.x), lat2merc(p.y)]
# ============================================================================ #
function centroid(poly::Shapefile.Polygon)

    a = 0.0
    c = [0.0, 0.0]
    ks, ke = largest_part(poly)

    @inbounds for k in ks:(ke-1)
        p = point2vec(poly.points[k])
        n = point2vec(poly.points[k+1])
        tmp = p[1] * n[2] - n[1] * p[2]
        a += tmp
        c .+= (p .+ n) .* tmp
    end
    a /= 2
    return c ./ (6 * a)
end
# ============================================================================ #
abstract type AbstractShape end
struct Polygon <: AbstractShape end
struct Point <: AbstractShape end

state_color(::Type{Polygon}) = "white"
state_color(::Type{Point}) = "black"
# ============================================================================ #
shape_file(::Type{BlockGroup}) = joinpath(DATADIR, "cb_2019_us_bg_500k.shp")
shape_file(::Type{Tract}) = joinpath(DATADIR, "cb_2019_us_tract_500k.shp")
shape_file(::Type{County}) = joinpath(DATADIR, "cb_2019_us_county_5m.shp")
shape_file(::Type{State}) = joinpath(DATADIR, "cb_2019_us_state_500k.shp")
to_state(::Type{BlockGroup}) = EpicastTables.BG2STATE
to_state(::Type{Tract}) = EpicastTables.TRACT2STATE
to_state(::Type{County}) = EpicastTables.COUNTY2STATE
to_geo(::Type{County}) = EpicastTables.TRACT2COUNTY
to_geo(::Type{Tract}) = 1
to_geo(::Type{BlockGroup}) = 1
geo_name(::Type{County}) = "county"
geo_name(::Type{Tract}) = "tract"
geo_name(::Type{BlockGroup}) = "block-group"
default_shape(::Type{County}) = Polygon
default_shape(::Type{Tract}) = Point
default_shape(::Type{BlockGroup}) = Point
default_outline(::Type{BlockGroup}) = Tract
default_outline(::Type{Tract}) = County
default_outline(::Type{County}) = State
# ============================================================================ #
struct GeoplotData{T<:AbstractGeo,N}
    shps::Vector{Shapefile.Polygon}
    fips::Vector{UInt64} # FIPS codes for shapes in <shps>
    data::FIPSTable{T,Float64,N}
end

data_matrix(g::GeoplotData, var::AbstractString) = g.data[var]
data_matrix(g::GeoplotData) = data_matrix(g, first(keys(g.data.var_index)))
EpicastTables.all_states(g::GeoplotData) = Set(all_states(g.data))
n_geo(g::GeoplotData) = size(data_matrix(g), 2)
n_shape(g::GeoplotData) = length(g.shps)
n_timepoint(g::GeoplotData) = size(data_matrix(g), 1)
n_state(g::GeoplotData) = length(all_states(g)) #length(g.states)
column_names(g::GeoplotData) = collect(keys(g.data.var_index))
# ============================================================================ #
outline_shapes(::Type{Tract}, d::GeoplotData{BlockGroup}) = shape_file(Tract), Set(all_tracts(d.data))
outline_shapes(::Type{County}, d::GeoplotData{BlockGroup}) = shape_file(County), Set(all_counties(d.data))
outline_shapes(::Type{State}, d::GeoplotData{BlockGroup}) = shape_file(State), Set(all_states(d.data))
outline_shapes(::Type{County}, d::GeoplotData{Tract}) = shape_file(County), Set(all_counties(d.data))
outline_shapes(::Type{State}, d::GeoplotData{Tract}) = shape_file(State), Set(all_states(d.data))
outline_shapes(::Type{County}, d::GeoplotData{County}) = shape_file(County), Set(all_counties(d.data))
outline_shapes(::Type{State}, d::GeoplotData{County}) = shape_file(State), Set(all_states(d.data))
# ============================================================================ #
function geoplot_data(::Type{T}, data_file::AbstractString,
    prefix::AbstractString="", ::Type{S}=UInt32;
    norms::Dict{String,String}=Dict{String,String}()
    ) where {T<:AbstractGeo, S <:Integer}

    if endswith(data_file, ".events.bin")
        data = Epicast.read_eventfile(Epicast.RunData, data_file)
    else
        data = Epicast.read_runfile(data_file, S)
    end

    if isempty(prefix)
        names = sort!(Epicast.column_names(data))
    elseif EpicastTables.has_var(data.data, prefix)
        names = [prefix]
    else
        names = Epicast.filter_columns(x -> startswith(x, prefix), data.data)
    end

    return geoplot_data(T, data, names; norms=norms)
end
# ---------------------------------------------------------------------------- #
function geoplot_data(::Type{T}, data_file::AbstractString, pat::Regex,
    ::Type{S}=UInt32; norms::Dict{String,String}=Dict{String,String}()
    ) where{T <:AbstractGeo, S <: Integer}

    if endswith(data_file, ".events.bin")
        data = Epicast.read_eventfile(Epicast.RunData, data_file)
    else
        data = Epicast.read_runfile(data_file, S)
    end

    names = Epicast.filter_columns(x -> match(pat, x) != nothing, data.data)

    return geoplot_data(T, data, names; norms=norms)
end
# ---------------------------------------------------------------------------- #
function filter_shapes!(shps::Vector{Shapefile.Polygon}, fips::Vector{<:Integer},
    data::FIPSTable)

    idx = findall(!in(EpicastTables.all_fips(data)), fips)
    deleteat!(shps, idx)
    deleteat!(fips, idx)

    return shps, fips
end
# ---------------------------------------------------------------------------- #
function geoplot_data(::Type{T}, all_data::Epicast.RunData,
    names::AbstractVector{<:AbstractString}, ::Type{S}=UInt32;
    norms::Dict{String,String}=Dict{String,String}(),
    ) where{T <:AbstractGeo, S <: Integer}
    data = data_dict(T, all_data, names, norms)

    states = Set(all_states(all_data.data))

    shps, fips = load_polygons(shape_file(T), states, to_state(T))

    filter_shapes!(shps, fips, data)

    return GeoplotData{T,3}(shps, fips, data)
end
# ---------------------------------------------------------------------------- #
aggregate_to(::Type{County}, data::FIPSTable) = aggregate_county(data, true)
aggregate_to(::Type{Tract}, data::FIPSTable) = aggregate_tract(data, true)
aggregate_to(::Type{BlockGroup}, data::FIPSTable) = data
# ---------------------------------------------------------------------------- #
function geoplot_data(data::FIPSTable{G,Float64,N}) where {G<:AbstractGeo,N}
    states = Set(all_states(data))
    shps, fips = load_polygons(shape_file(G), states, to_state(G))
    tmp = aggregate_to(G, data)
    filter_shapes!(shps, fips, tmp)
    return GeoplotData{G,N}(shps, fips, tmp)
end
# ---------------------------------------------------------------------------- #
function geoplot_data(data::FIPSTable{G,L,N}) where {G<:AbstractGeo,L<:Integer,N}
    tmp = FIPSTable{G,Float64,N}(data.fips_index, data.var_index,
        Array{Float64,N}(data.data))
    return geoplot_data(tmp)
end
# ============================================================================ #
function draw_shapes!(::Type{Polygon}, ax, data::GeoplotData, var::AbstractString,
    cm::ColorMap, norm::AbstractNorm, frame::Integer=1)

    polys = Vector{Matrix{Float64}}(undef, n_shape(data))
    for k in 1:n_shape(data)
        ks, ke = largest_part(data.shps[k])
        pts = view(data.shps[k].points, ks:ke)
        polys[k] = [long2merc.(getfield.(pts, :x)) lat2merc.(getfield.(pts, :y))]
    end

    c = map(x -> data.data[frame, x, var], data.fips)
    hp = matplotlib.collections.PolyCollection(polys, facecolor=cm(scale(norm, c)))
    ax.add_collection(hp)
    return hp
end
# ============================================================================ #
function draw_shapes!(::Type{Point}, ax, data::GeoplotData, var::AbstractString,
    cm::ColorMap, ext::Tuple{Float64,Float64}, norm::AbstractNorm, frame::Integer=1)
    xy = centroid.(data.shps)
    c = map(x -> data.data[frame, x, var], data.fips)
    return ax.scatter(getindex.(xy, 1), getindex.(xy, 2), 4.0, c=cm(scale(norm, c)))
end
# ============================================================================ #
quantile_threshold(::Type{County}) = 0.999
quantile_threshold(::Type{<:AbstractGeo}) = 0.99
# ============================================================================ #
function nearest_shape(data::GeoplotData, idx::AbstractArray{<:Integer}, mouseevt)
    n = length(idx)
    kt = 1
    if n > 1
        x = mouseevt.xdata
        y = mouseevt.ydata
        kt = argmin(1:n) do k
            c = centroid(data.shps[idx[k]+1])
            return sqrt((c[1]-x)^2 + (c[2]-y)^2)
        end
    end
    return kt
end
# ============================================================================ #
function map_figure(data::GeoplotData{T}, var::AbstractString, frame::Integer;
    maxq::Real=quantile_threshold(T), cmap::AbstractString="viridis",
    norm::Type{<:AbstractNorm}=ExtremaNorm, title::AbstractString="",
    cb_label::AbstractString="", shape::Type{<:AbstractShape}=default_shape(T),
    outline::Type{<:AbstractGeo}=default_outline(T),
    vmax::Real=NaN) where T<:AbstractGeo

    h, ax = subplots(1,1)
    h.set_size_inches((10,6))
    norm, cm, hp = add_map!(ax, data, var, frame, maxq=maxq, cmap=cmap,
        norm=norm, shape=shape, outline=outline, mx=vmax)

    ax.set_title(title, fontsize=18)

    cb = h.colorbar(
        PyPlot.matplotlib.cm.ScalarMappable(norm=mpl_norm(norm), cmap=cm),
        ax = ax
    )

    if !isempty(cb_label)
        cb.set_label(cb_label, fontsize=14)
    end

    onpick(evt) = begin
        if evt.mouseevent.button == 1
            b, l = hp.contains(evt.mouseevent)
            if b
                n = length(l["ind"])
                kt = nearest_shape(data, l["ind"], evt.mouseevent)
                idx = l["ind"][kt] + 1
                fips_code = data.fips[idx]
                str = @sprintf("%s %d: %.03f", titlecase(geo_name(T)),
                    fips_code, data.data[frame,fips_code,var])
                str = !isempty(title) ? title * " | " * str : str
                ax.set_title(str, fontsize=18)
            end
        end
    end

    h.canvas.mpl_connect("pick_event", onpick)
    h.tight_layout()

    return h, ax, cb
end
# ============================================================================ #
function add_map!(ax::PyCall.PyObject, data::GeoplotData{T}, var::AbstractString,
    frame::Integer; mn::Real=NaN, mx::Real=NaN, maxq::Real=quantile_threshold(T),
    cmap::AbstractString="viridis", norm::Type{<:AbstractNorm}=ExtremaNorm,
    shape::Type{<:AbstractShape}=default_shape(T),
    outline::Type{<:AbstractGeo}=default_outline(T)) where T<:AbstractGeo

    outline_shp, outline_fips = outline_shapes(outline, data)

    data_mat = data_matrix(data, var)
    mn = isnan(mn) ? minimum(data_mat) : mn
    mx = isnan(mx) ? quantile(vec(data_mat), maxq) : mx

    cm = PyPlot.get_cmap(cmap)

    norm = norm(mn, mx)

    hp = draw_shapes!(shape, ax, data, var, cm, norm, frame)

    state_outlines!(ax, outline_shp, outline_fips, state_color(shape))

    ax.set_title("Day $(frame-1)", fontsize=18)

    foreach(ax.spines) do k
        ax.spines[k].set_visible(false)
    end

    ax.set_yticks([])
    ax.set_xticks([])

    return norm, cm, hp
end
# ============================================================================ #
function ticks_to_dates(ax, nt::Integer, start_date::AbstractString,
    month_fmt::AbstractString="m/yy")

    if nt > 30
        use_week = nt < 60
        start = Dates.Date(start_date)
        ref = Dates.Date(Dates.year(start), Dates.month(start), 1)
        xt = Int[]
        xtl = String[]
        while getfield(ref - start, :value) <= nt
            if ref >= start
                push!(xtl, use_week ?
                    Dates.format(ref,"u d") :
                    Dates.format(ref, month_fmt))
                push!(xt, getfield(ref - start, :value))
            end
            ref += use_week ? Dates.Day(7) : Dates.Month(1)
        end
    else
        xt = filter!(x -> 0 <= x < nt, ax.get_xticks())
        xtl = map(x -> start + Dates.Day(x), xt)
    end

    return xt, xtl
end
# ============================================================================ #
function add_state_timeseries!(ax, data::GeoplotData, var::AbstractString,
    frame::Integer=1, vertical::Bool=false, start_date::AbstractString="",
    top::Integer=typemax(Int), AT::Type{<:AbstractGeo}=State)

    # if data are already normalized (cases-per-100k) then simply averaging
    # will maintain the proper units
    state_data = EpicastTables.aggregate(AT, data.data, true)
    states = EpicastTables.all_geo(AT, state_data)
    nstate = length(states)

    colors = map(col -> (red(col), green(col), blue(col)),
        distinguishable_colors(nstate, [RGB(1,1,1), RGB(0,0,0)],
            dropseed=true)
    )

    states_tmp = sort!(collect(states))
    if top < length(states_tmp)
        mx = [maximum(state_data[k, var]) for k in states_tmp]
        states_use = states_tmp[sortperm(mx, rev=true)[1:top]]
    else
        states_use = states_tmp
    end

    nt = n_timepoint(data)
    j = 1
    mx2 = -Inf
    for k in states_use
        v = state_data[k, var]
        ax.plot(0:(nt-1), v, linewidth=2, color=colors[j], label=get(STATE_FIPS, k, string(k)),
            picker=true)
        mx2 = max(mx2, maximum(v))
        j += 1
    end

    ax.set_xlim(-1, nt)

    ax.spines["right"].set_visible(false)
    ax.spines["top"].set_visible(false)
    ax.set_ylabel("New cases per 100k residents", fontsize=14)

    if !isempty(start_date)
        xt, xtl = ticks_to_dates(ax, nt, start_date)

        ax.set_xticks(xt, xtl, rotation=45)
        ax.set_xlabel("Date", fontsize=14)
    else
        ax.set_xlabel("Simulation day", fontsize=14)
    end

    ncol = length(states_use) > 25 ? 2 : 1

    if vertical
        ax.legend(frameon=true, loc="lower left", bbox_to_anchor=(1.02, 0.0), ncol=ncol)
    else
        ax.legend(frameon=true, loc="upper left", bbox_to_anchor=(1.02, 1.0), ncol=ncol)
    end

    ax.set_title(" ", fontsize=18)
    mx2 *= 1.05
    time_idc = ax.plot([frame-1,frame-1],[0,mx2], "--", color="darkgray", linewidth=2)[1]

    return mx2, time_idc
end
# ============================================================================ #
function make_figure(data::GeoplotData{T}; ofile::AbstractString="",
    style_geo::Function=identity, style_line::Function=identity,
    maxq::Real=quantile_threshold(T), frame::Integer=1, vertical::Bool=false,
    norm::Type{<:AbstractNorm}=ExtremaNorm,
    shape::Type{<:AbstractShape}=default_shape(T),
    outline::AbstractGeo=default_outline(T),
    agg_level::Type{<:AbstractGeo}=State) where T<:AbstractGeo

    var = first(keys(data.data.var_index))

    return make_figure(data, var, ofile=ofile,
                       style_geo=style_geo,
                       style_line=style_line,
                       maxq=maxq, frame=frame,
                       vertical=vertical,
                       norm=norm, shape=shape,
                       outline=outline,
                       agg_level=agg_level)
end
# ---------------------------------------------------------------------------- #
function make_figure(data::GeoplotData{T}, var::AbstractString;
    ofile::AbstractString="", style_geo::Function=identity,
    style_line::Function=identity,
    maxq::Real=quantile_threshold(T), frame::Integer=1,
    vertical::Bool=false,
    norm::Type{<:AbstractNorm}=ExtremaNorm,
    shape::Type{<:AbstractShape}=default_shape(T),
    outline::AbstractGeo=default_outline(T),
    agg_level::Type{<:AbstractGeo}=State) where T<:AbstractGeo

    nt = n_timepoint(data)

    @assert(0 < frame <= nt, "given frame $(frame) is out-of-bounds")

    nstate = n_state(data)
    width = nstate > 25 ? 15 : 14
    w_ratio = [1.0, 1.0]
    if 25 < nstate <= 40
        w_ratio .= [1.3, 1.0]
    elseif nstate > 40
        w_ratio .= [1.7, 1.0]
    end

    if vertical
        h, ax = subplots(2, 1)
        h.set_size_inches((width/2,10))
    else
        h, ax = subplots(1, 2, width_ratios=w_ratio)
        h.set_size_inches((width,6.5))
    end

    norm, cm, hp = add_map!(ax[1], data, var, frame, maxq=maxq,
        norm=norm, shape=shape, outline=outline)

    mx2, time_idc = add_state_timeseries!(ax[2], data, var, frame, vertical, "",
        typemax(Int), agg_level)

    county_line = nothing

    IDX = frame

    update_figure(idx::Integer) = begin
        c = map(x -> data.data[idx, x, var], data.fips)
        hp.set_facecolors(cm(scale(norm, c)))
        ax[1].set_title("Day " * string(idx-1), fontsize=18)
        time_idc.set_xdata([idx-1, idx-1])
    end

    onscroll(evt) = begin
        tmp = evt.button == "up" ? IDX + 1 : IDX - 1
        tmp = max(1, min(nt, tmp))
        if tmp != IDX
            IDX = tmp
            update_figure(IDX)
        end
    end

    onpick(evt) = begin
        if evt.mouseevent.button == 1
            b, l = hp.contains(evt.mouseevent)
            if b
                kt = nearest_shape(data, l["ind"], evt.mouseevent)
                idx = l["ind"][kt] + 1
                fips_code = data.fips[idx]
                v = data.data[fips_code, var]
                mx_use = max(mx2, maximum(v))
                if county_line == nothing
                    county_line = ax[2].plot(0:(nt-1), v,
                        color="black", linewidth=2.5)[1]
                else
                    county_line.set_ydata(v)
                end
                ax[2].set_ylim(0, mx_use)
                time_idc.set_ydata([0, mx_use])
                gn = titlecase(geo_name(T)) * " "
                ax[2].set_title(gn * string(fips_code), fontsize=18)
            elseif evt.artist.axes == ax[2]
                lab = evt.artist.get_label()
                ax[2].set_title("State " * lab, fontsize=18)
            end
        end
    end

    h.tight_layout()

    if vertical
        h.subplots_adjust(hspace=0)
        bb = ax[1].get_position()
        ax[1].set_position([0.02, bb.y0, 0.9, bb.height])
    end

    if !isempty(ofile)
        if endswith(ofile, ".mp4")
            anim = animation.FuncAnimation(h, update_figure, frames=2:nt)
            anim.save(ofile, fps=3)
        else
            h.savefig(ofile, dpi=200)
        end
    else
        h.canvas.mpl_connect("scroll_event", onscroll)
        h.canvas.mpl_connect("pick_event", onpick)
    end

    return h, ax
end
# ============================================================================ #
function load_polygons(ifile::AbstractString, geo::AbstractSet{<:Integer}, level::Integer=1)
    tbl = Shapefile.Table(ifile)
    fips = parse.(Int, tbl.GEOID)
    idx = findall(x->in(div(x, level), geo), fips)
    out = filter!(!ismissing, Shapefile.shapes(tbl)[idx])
    return convert(Vector{Shapefile.Polygon}, out), fips[idx]
end
# ============================================================================ #
# /Users/palexander/Documents/emerge+radium/geo-data/cb_2019_us_state_500k/cb_2019_us_state_500k.shp
function state_outlines!(ax, shpfile::AbstractString, states::AbstractSet{<:Integer}, color::AbstractString="white")
    shps, _ = load_polygons(shpfile, states, 1)
    for shp in shps
        add_state_outline!(ax, shp, edgecolor=color, color="none")
    end
end
# ============================================================================ #
function add_state_outline!(ax, shp::Shapefile.Polygon; color=nothing,
    edgecolor=nothing, label=nothing)

    h = Vector{PyPlot.PyCall.PyObject}(undef, length(shp.parts))
    for k in 1:length(shp.parts)
        ks = shp.parts[k] + 1
        ke = k < length(shp.parts) ? shp.parts[k+1] : length(shp.points)
        h[k] = ax.fill(
            long2merc.(getfield.(shp.points[ks:ke], :x)),
            lat2merc.(getfield.(shp.points[ks:ke], :y)),
            facecolor=color,
            edgecolor=edgecolor,
            label=label,
            picker=label!=nothing,
            linewidth=STATE_OUTLINE_WIDTH)[1]
    end
    return h
end
# ============================================================================ #
function smooth_data!(d::GeoplotData)
    for slice in eachslice(d.data.data, dims=3)
        for v in eachcol(slice)
            for k = 8:length(v)
                v[k] = mean(v[k-7:k])
            end
        end
    end
end
# ============================================================================ #
function epidemic_overview(data::GeoplotData, var::AbstractString,
    frames::AbstractVector{<:Integer}, ax::Vector{PyCall.PyObject},
    start_date::AbstractString; labs = ["A.", "B."], laby::Real=0.99,
    lab_loc::AbstractString="left", timeseries_geo::Type{<:AbstractGeo}=State,
    n_geo::Integer=15, vpad::Real=0.05)

    h = ax[1].figure

    add_frame!(h, ax[1], ax[end-1], data, var, frames[1])

    for k in 2:length(frames)
        add_frame!(h, ax[k], nothing, data, var, frames[k])
    end

    mx2, time_idc = add_state_timeseries!(ax[end], data, var,
        frames[1], false, start_date, n_geo, timeseries_geo)

    ax[end].text(frames[1]-1,time_idc.get_ydata()[2], string(frames[1]-1),
        fontsize=12, va="bottom", ha="center")

    for t in frames[2:end]
        ax[end].plot([t-1, t-1], time_idc.get_ydata(), "--", color="darkgray",
            linewidth=2)
        ax[end].text(t-1,time_idc.get_ydata()[2], string(t-1), fontsize=12,
            va="bottom", ha="center")
    end

    bbox = axes_perimiter(ax[end])

    position_epi_frames(ax[1:end-1], bbox, 0.0, vpad)

    pos = axes_perimiter(ax[1]).bounds
    h.text(pos[1],laby,labs[1], fontsize=30, va="top", ha=lab_loc)

    pos = axes_perimiter(ax[end]).bounds
    h.text(pos[1],laby,labs[2], fontsize=30, va="top", ha=lab_loc)

    return h, ax
end
# ============================================================================ #
function add_frame!(h, ax, cbax, data::GeoplotData, var::AbstractString,
    frame::Integer)

    norm, cm, hp = add_map!(ax, data, var, frame)

    if cbax != nothing
        cb = h.colorbar(
            PyPlot.matplotlib.cm.ScalarMappable(norm=mpl_norm(norm),
                cmap=cm),
            cax = cbax
        )

        cb.set_label("New cases per 100k residents", fontsize=12)
    else
        cb = nothing
    end

    return cb
end
# ============================================================================ #
function position_epi_frames(ax, bbox, bpad, vpad)
    bt = bbox.y0 + bpad
    ht = bbox.y1 - bbox.y0

    ncol = Int((length(ax) - 1) / 2)

    for k = 1:ncol
        pos = ax[k].get_position().bounds
        ax[k].set_position([pos[1], bt + (ht/2) + bpad, pos[3], ht/2 - vpad])
    end

    for k = (ncol+1):(ncol*2)
        pos = ax[k].get_position().bounds
        ax[k].set_position([pos[1], bt, pos[3], ht/2 - vpad])
    end

    k = ncol*2 + 1

    pos = ax[k].get_position().bounds
    ax[k].set_position([pos[1], bt, pos[3], ht - vpad])
end
# ============================================================================ #
function axes_perimiter(ax)
    h = ax.figure
    h.draw_without_rendering()
    bbpx = ax.get_tightbbox(h.canvas.get_renderer())
    return bbpx.transformed(h.transFigure.inverted())
end
# ============================================================================ #
const STATE_FIPS = Dict(
    1 => "AL", 2 => "AK", 4 => "AZ", 5 => "AR", 6 => "CA", 8 => "CO", 9 => "CT",
    10 => "DE", 11 => "DC", 12 => "FL", 13 => "GA", 15 => "HI", 16 => "ID",
    17 => "IL", 18 => "IN", 19 => "IA", 20 => "KS", 21 => "KY", 22 => "LA",
    23 => "ME", 24 => "MD", 25 => "MA", 26 => "MI", 27 => "MN", 28 => "MS",
    29 => "MO", 30 => "MT", 31 => "NE", 32 => "NV", 33 => "NH", 34 => "NJ",
    35 => "NM", 36 => "NY", 37 => "NC", 38 => "ND", 39 => "OH", 40 => "OK",
    41 => "OR", 42 => "PA", 44 => "RI", 45 => "SC", 46 => "SD", 47 => "TN",
    48 => "TX", 49 => "UT", 50 => "VT", 51 => "VA", 53 => "WA", 54 => "WV",
    55 => "WI", 56 => "WY"
)
# ============================================================================ #
# From: https://www2.census.gov/geo/pdfs/maps-data/maps/reference/us_regdiv.pdf
const DIVISION_FIPS = Dict(
    [9,23,25,33,44,50] => "NewEng", #"New England",
    [34,36,42] => "MidAtl", #"Middle Atlantic",
    [18,17,26,39,55] => "ENCen", #"East North Central",
    [19, 20,27,29,31,38,46] => "WNCen", #"West North Central",
    [10,11,12,13,24,37,45,51,54] => "SAtl", #"South Atlantic",
    [1,21,28,47] => "ESCen", #"East South Central",
    [5,22,40,48] => "WSCen", #"West South Central",
    [4,8,12,35,30,49,32,56] => "Mtn", #"Mountain",
    [2,6,15,41,53] => "Pac", #"Pacific"
)
# ============================================================================ #
end # module EpicastGeoplot
