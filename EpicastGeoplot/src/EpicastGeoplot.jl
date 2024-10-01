module EpicastGeoplot

using Shapefile, PyPlot, Colors, Printf, Statistics

import EpicastTables

using Epicast, EpicastTables

using PyCall

export Polygon, Point, County, Tract, State

include("./normalize.jl")

const animation = PyNULL()

function __init__()
    copy!(animation, pyimport("matplotlib.animation"))
end

const DATADIR = joinpath(@__DIR__, "..", "assets", "geo-data", "data")

# ============================================================================ #
function case_count!(out::AbstractVector, x::Epicast.RunData,
    var::AbstractString, idx::AbstractVector{<:Integer})

    in = view(Epicast.rundata(x, var), :, idx)

    if Epicast.has_demographic(x, var)
        # number of new cases relative to total size of that demographic
        out .= Epicast.new_cases(in)
        out ./= sum(view(Epicast.demographics(x, var), idx))
        out .*= 1e5

    elseif endswith(var, r"_age\d")
        # number of new [hospitialized, icu, ventiated, dead] for given age
        # group relative to that age group's total population (for each
        # location)
        age = "age_" * var[end]

        out .= Epicast.new_cases(in)
        out ./= sum(view(Epicast.demographics(x, age), idx))
        out .*= 1e5

    elseif startswith(var, "infection-src")
        # % of infections attributable to given context
        out .= Epicast.new_cases(in)
        out ./= Epicast.new_cases(view(Epicast.rundata(x, "total"), idx, :))

    elseif startswith(var, "status_")
        # % of total agents w/in geography that are in given state
        out .= Epicast.total_cases(in)
        out ./= sum(view(Epicast.demographics(x, "total"), idx))
        out .*= 1e5

    else
        # default: total counts per location
        out .= Epicast.total_cases(in)
    end

    return out
end
# ============================================================================ #
function data_dict(data::Epicast.RunData, names::Vector{<:AbstractString},
    conv::Integer=EpicastTables.TRACT2COUNTY)

    if conv > 1
        # dat is: time x fips x var
        dat, grp = Epicast.group_by(data, names, conv, case_count!)
    else
        # time x fips x var
        dat = Array{Float64,3}(undef, Epicast.n_timepoint(data),
            length(data.fips), length(names))
        
        # this is pretty slow, and kinda inefficient, but it's nice to reuse
        # case_count!()
        for k in 1:size(dat, 2)
            for j in 1:size(dat, 3)
                case_count!(view(dat, :, k, j), data, names[j], [k])
            end
        end

        grp = data.fips
    end

    replace!(dat, NaN=>0)

    idx = Dict{UInt64,Int}(grp[k] => k for k in eachindex(grp))
    var = Dict{String,Int}(names[k] => k for k in eachindex(names))
    return FIPSTable(dat, idx, var)
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
@inline point2vec(p::Shapefile.Point) = [p.x, p.y]
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

abstract type AbstractGeo end
struct Tract <: AbstractGeo end
struct County <: AbstractGeo end
struct State <: AbstractGeo end

# ============================================================================ #
shape_file(::Type{County}) = joinpath(DATADIR, "cb_2019_us_county_5m.shp")
shape_file(::Type{Tract}) = joinpath(DATADIR, "cb_2019_us_tract_500k.shp")
shape_file(::Type{State}) = joinpath(DATADIR, "cb_2019_us_state_500k.shp")
to_state(::Type{County}) = EpicastTables.COUNTY2STATE
to_state(::Type{Tract}) = EpicastTables.TRACT2STATE
to_geo(::Type{County}) = EpicastTables.TRACT2COUNTY
to_geo(::Type{Tract}) = 1
geo_name(::Type{County}) = "county"
geo_name(::Type{Tract}) = "tract"
default_shape(::Type{County}) = Polygon
default_shape(::Type{Tract}) = Point
default_outline(::Type{Tract}) = County
default_outline(::Type{County}) = State
# ============================================================================ #
struct GeoplotData{T<:AbstractGeo,N}
    shps::Vector{Shapefile.Polygon}
    fips::Vector{UInt64} # FIPS codes for shapes in <shps>
    data::FIPSTable{Float64,N}
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
outline_shapes(::Type{County}, d::GeoplotData{Tract}) = shape_file(County), Set(all_counties(d.data))
outline_shapes(::Type{State}, d::GeoplotData{Tract}) = shape_file(State), Set(all_states(d.data))
outline_shapes(::Type{County}, d::GeoplotData{County}) = shape_file(County), Set(all_counties(d.data))
outline_shapes(::Type{State}, d::GeoplotData{County}) = shape_file(State), Set(all_states(d.data))
# ============================================================================ #
function geoplot_data(::Type{T}, data_file::AbstractString,
    prefix::AbstractString="") where T<:AbstractGeo

    if endswith(data_file, ".events.bin")
        data = Epicast.read_eventfile(Epicast.RunData, data_file)
    else
        data = Epicast.read_runfile(data_file)
    end

    if isempty(prefix)
        names = sort!(Epicast.column_names(data))
    elseif EpicastTables.has_var(data.data, prefix)
        names = [prefix]
    else
        names = Epicast.filter_columns(x -> startswith(x, prefix), data.data)
    end

    return geoplot_data(T, data, names)
end
# ---------------------------------------------------------------------------- #
function geoplot_data(::Type{T}, data_file::AbstractString, pat::Regex) where T<:AbstractGeo

    if endswith(data_file, ".events.bin")
        data = Epicast.read_eventfile(Epicast.RunData, data_file)
    else
        data = Epicast.read_runfile(data_file)
    end

    names = Epicast.filter_columns(x -> match(pat, x) != nothing, data.data)

    return geoplot_data(T, data, names)
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
    names::AbstractVector{<:AbstractString}) where T<:AbstractGeo

    data = data_dict(all_data, names, to_geo(T))

    states = Set(all_states(all_data.data))

    shps, fips = load_polygons(shape_file(T), states, to_state(T))

    filter_shapes!(shps, fips, data)

    return GeoplotData{T,3}(shps, fips, data)
end
# ---------------------------------------------------------------------------- #
aggregate_to(::Type{County}, data::FIPSTable) = aggregate_county(data, true)
aggregate_to(::Type{Tract}, data::FIPSTable) = data
# ---------------------------------------------------------------------------- #
function geoplot_data(::Type{T}, data::FIPSTable{Float64,N}) where {T<:AbstractGeo,N}
    states = Set(all_states(data))
    shps, fips = load_polygons(shape_file(T), states, to_state(T))
    tmp = aggregate_to(T, data)
    filter_shapes!(shps, fips, tmp)
    return GeoplotData{T,N}(shps, fips, tmp)
end
# ---------------------------------------------------------------------------- #
function geoplot_data(::Type{T}, data::FIPSTable{L,N}) where {T<:AbstractGeo,L<:Integer,N}
    tmp = FIPSTable{Float64,N}(data.fips_index, data.var_index,
        Array{Float64,N}(data.data))
    return geoplot_data(T, tmp)
end
# ============================================================================ #
function draw_shapes!(::Type{Polygon}, ax, data::GeoplotData, var::AbstractString,
    cm::ColorMap, norm::AbstractNorm, frame::Integer=1)

    polys = Vector{Matrix{Float64}}(undef, n_shape(data))
    for k in 1:n_shape(data)
        ks, ke = largest_part(data.shps[k])
        pts = view(data.shps[k].points, ks:ke)
        polys[k] = [getfield.(pts, :x) getfield.(pts, :y)]
    end

    c = map(x -> data.data[frame, x, var], data.fips)
    hp = matplotlib.collections.PolyCollection(polys, facecolor=cm(scale(norm, c)))
    ax.add_collection(hp)
    return hp
end
# ============================================================================ #
function draw_shapes!(::Type{Point}, ax, data::GeoplotData, var::AbstractString,
    cm::ColorMap, norm::AbstractNorm, frame::Integer=1)
    
    xy = centroid.(data.shps)
    c = map(x -> data.data[frame, x, var], data.fips)
    return ax.scatter(getindex.(xy, 1), getindex.(xy, 2), 4.0, c=cm(scale(norm, c)))
end
# ============================================================================ #
quantile_threshold(::Type{County}) = 0.999
quantile_threshold(::Type{Tract}) = 0.99
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
    outline::Type{<:AbstractGeo}=default_outline(T)) where T<:AbstractGeo
    
    h, ax = subplots(1,1)
    h.set_size_inches((10,6))
    norm, cm, hp = add_map!(ax, data, var, frame, maxq=maxq, cmap=cmap,
        norm=norm, shape=shape, outline=outline)
    
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
function add_state_timeseries!(ax, data::GeoplotData, var::AbstractString,
    frame::Integer=1, smooth::Bool=false, vertical::Bool=false)

    # if data are already normalized (cases-per-100k) then simply averaging
    # will maintain the proper units
    state_data = aggregate_state(data.data, true)
    states = all_states(state_data)
    nstate = length(states)

    colors = map(col -> (red(col), green(col), blue(col)), 
        distinguishable_colors(nstate, [RGB(1,1,1), RGB(0,0,0)],
            dropseed=true)
    )

    smooth && @info("smoothing...")

    nt = n_timepoint(data)
    j = 1
    mx2 = -Inf
    for k in sort!(collect(states))
        v = state_data[k, var]
        if smooth
            v = copy(v)
            for k = 8:length(v)
                v[k] = mean(v[k-7:k])
            end
        end
        ax.plot(0:(nt-1), v, linewidth=2, label=STATE_FIPS[k], color=colors[j],
            picker=true)
        mx2 = max(mx2, maximum(v))
        j += 1
    end

    ax.spines["right"].set_visible(false)
    ax.spines["top"].set_visible(false)
    ax.set_ylabel("New cases per 100k residents", fontsize=14)
    ax.set_xlabel("Simulation day", fontsize=14)

    ncol = length(states) > 25 ? 2 : 1

    if vertical
        ax.legend(frameon=true, loc="lower left", bbox_to_anchor=(1.02, 0.0), ncol=ncol)
    else
        ax.legend(frameon=true, loc="upper left", bbox_to_anchor=(1.02, 1.0), ncol=ncol)
    end

    ax.set_title(" ", fontsize=18)

    time_idc = ax.plot([frame-1,frame-1],[0,mx2], "--", color="darkgray", linewidth=2)[1]

    return mx2, time_idc
end
# ============================================================================ #
function make_figure(data::GeoplotData{T}; ofile::AbstractString="",
    maxq::Real=quantile_threshold(T), frame::Integer=1, vertical::Bool=false,
    smooth::Bool=false, norm::Type{<:AbstractNorm}=ExtremaNorm,
    shape::Type{<:AbstractShape}=default_shape(T),
    outline::AbstractGeo=default_outline(T)) where T<:AbstractGeo

    var = first(keys(data.data.var_index))

    return make_figure(data, var, ofile=ofile, maxq=maxq, frame=frame,
        vertical=vertical, smooth=smooth, norm=norm, shape=shape, outline=outline)
end
# ---------------------------------------------------------------------------- #
function make_figure(data::GeoplotData{T}, var::AbstractString;
    ofile::AbstractString="", maxq::Real=quantile_threshold(T),
    frame::Integer=1, vertical::Bool=false, smooth::Bool=false,
    norm::Type{<:AbstractNorm}=ExtremaNorm,
    shape::Type{<:AbstractShape}=default_shape(T),
    outline::Type{<:AbstractGeo}=default_outline(T)) where T<:AbstractGeo

    nt = n_timepoint(data)

    @assert(0 < frame <= nt, "give frame $(frame) is out-of-bounds")

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

    mx2, time_idc = add_state_timeseries!(ax[2], data, var, frame, smooth,
        vertical)

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
                mx_use = max(mx2, maximum(data.data[fips_code, var]))
                if county_line == nothing
                    county_line = ax[2].plot(0:(nt-1), data.data[fips_code, var],
                        color="black", linewidth=2.5)[1]
                else
                    county_line.set_ydata(data.data[fips_code, var])
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
function load_polygons(ifile::AbstractString, geo::Set{<:Integer}, level::Integer=1)
    tbl = Shapefile.Table(ifile)
    fips = parse.(Int, tbl.GEOID)
    idx = findall(x->in(div(x, level), geo), fips)
    out = filter!(!ismissing, Shapefile.shapes(tbl)[idx])
    return convert(Vector{Shapefile.Polygon}, out), fips[idx]
end
# ============================================================================ #
# /Users/palexander/Documents/emerge+radium/geo-data/cb_2019_us_state_500k/cb_2019_us_state_500k.shp
function state_outlines!(ax, shpfile::AbstractString, states::Set{<:Integer}, color::AbstractString="white")
    shps, _ = load_polygons(shpfile, states, 1)
    for shp in shps
        add_state!(ax, shp, edgecolor=color, color="none")
    end
end
# ============================================================================ #
function add_state!(ax, shp::Shapefile.Polygon; color=nothing,
    edgecolor=nothing, label=nothing)

    h = Vector{PyPlot.PyCall.PyObject}(undef, length(shp.parts))
    for k in 1:length(shp.parts)
        ks = shp.parts[k] + 1
        ke = k < length(shp.parts) ? shp.parts[k+1] : length(shp.points)
        h[k] = ax.fill(getfield.(shp.points[ks:ke], :x),
            getfield.(shp.points[ks:ke], :y), facecolor=color,
            edgecolor=edgecolor, label=label, picker=label!=nothing)[1]
    end
    return h
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
end # module EpicastGeoplot
