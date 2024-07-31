module EpicastGeoplot

using Shapefile, PyPlot, Colors, Printf, Statistics

using Epicast

using PyCall

export CountyPolygon, TractPoint

const animation = PyNULL()

function __init__()
    copy!(animation, pyimport("matplotlib.animation"))
end

const DATADIR = joinpath(@__DIR__, "..", "assets", "geo-data", "data")

# ============================================================================ #
struct FIPSTable{N}
    fips_index::Dict{UInt64,Int}
    var_index::Dict{String,Int}
    data::Array{Float64,N}# time x fips x var
end
all_fips(f::FIPSTable) = collect(keys(f.fips_index))
has_fips(f::FIPSTable, fips::Integer) = haskey(f.fips_index, fips)
Base.getindex(f::FIPSTable{2}, fips::Integer) = view(f.data, :, f.fips_index[fips]) #f.data[:,f.fips_index[k]]
Base.getindex(f::FIPSTable{3}, fips::Integer, var::AbstractString) = view(f.data, :, f.fips_index[fips], f.var_index[var])
Base.getindex(f::FIPSTable{3}, var::AbstractString) = view(f.data, :, :, f.var_index[var])

struct FIPSTableIterator
    data::FIPSTable{3}
    var::String
end

function Base.iterate(f::FIPSTable{2}, k::Integer=1)
    tmp = iterate(f.fips_index, k)
    tmp == nothing && return tmp
    return (tmp[1].first, f[tmp[1].second]), tmp[2]
end
Base.IteratorSize(::FIPSTable) = Base.HasLength()
Base.IteratorEltype(::FIPSTable) = Base.HasEltype()
Base.length(g::FIPSTable) = length(g.fips_index)

# Base.eltype(::FIPSTable{2}) = Tuple{UInt64,SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}}
Base.eltype(::FIPSTable{N}) where N = Tuple{UInt64,SubArray{Float64, 1, Array{Float64,N}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}}

eachfips(f::FIPSTable{2}) = f

eachfips(f::FIPSTable{3}, var::AbstractString) = FIPSTableIterator(f, var)

function Base.iterate(f::FIPSTableIterator, k::Integer=1)
    tmp = iterate(f.fips_index, k)
    tmp == nothing && return tmp
    return (tmp[1].first, f.data[tmp[1].second, f.var]), tmp[2]
end
# ============================================================================ #
function case_count!(out::AbstractVector, x::Epicast.RunData,
    var::AbstractString, idx::AbstractVector{<:Integer}, norm::AbstractString="")

    in = view(Epicast.rundata(x, var), idx, :)

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
        #out .*= 1e5

    elseif length(norm) == 2 && norm[1] == "total" && Epicast.has_demographic(x, norm[2])
        # number of total cases relative to total size of that demographic
        out .= Epicast.total_cases(in)
        out ./= sum(view(Epicast.demographics(x, norm[2]), idx))
        #out .*= 1e5

    elseif endswith(norm, r"age_\d")
        # number of new [hospitialized, icu, ventiated, dead] for given age
        # group relative to that age group's total population (for each
        # location)
        age = "age_" * norm[end]

        out .= Epicast.new_cases(in)
        out ./= sum(view(Epicast.demographics(x, age), idx))
        #out .*= 1e5

    elseif norm == "new_cases_prop"
        # % of infections attributable to given context
        out .= Epicast.new_cases(in)
        out ./= Epicast.new_cases(view(Epicast.rundata(x, "total"), idx, :))

    elseif norm == "change_prop"
        # change in state as a % of total agents w/in geography
        out .= Epicast.new_cases(in)
        out ./= sum(view(Epicast.demographics(x, "total"), idx, :))
        #out .*= 1e5

    elseif norm == "prop"
        # % of total agents w/in geography that are in given state
        out .= Epicast.total_cases(in)
        out ./= sum(view(Epicast.demographics(x, "total"), idx))
        #out .*= 1e5

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
function data_dict(data::Epicast.RunData, names::Vector{<:AbstractString},
    conv::Integer=TRACT2COUNTY, norms::Dict{String,String}=Dict{String,String}())

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
    return FIPSTable{3}(idx, var, dat)
end
# ============================================================================ #
function largest_part(poly::Shapefile.Polygon)
    ks = poly.parts[end] + 1
    ke = length(poly.points)
    mx = ke - poly.parts[end]
    for k in 2:length(poly.parts)
        tmp = poly.parts[k] - poly.parts[k-1]
        if tmp > mx
            mx = tmp
            ks = poly.parts[k-1] + 1
            ke = poly.parts[k]
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

    #@inbounds
    for k in ks:(ke-1)
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
struct CountyPolygon <: AbstractShape end
struct TractPoint <: AbstractShape end

shape_file(::Type{CountyPolygon}) = joinpath(DATADIR, "cb_2019_us_county_5m.shp")
shape_file(::Type{TractPoint}) = joinpath(DATADIR, "cb_2019_us_tract_500k.shp")
to_state(::Type{CountyPolygon}) = COUNTY2STATE
to_state(::Type{TractPoint}) = TRACT2STATE
to_geo(::Type{CountyPolygon}) = TRACT2COUNTY
to_geo(::Type{TractPoint}) = 1
state_color(::Type{CountyPolygon}) = "white"
state_color(::Type{TractPoint}) = "black"
# ============================================================================ #
struct GeoplotData{T<:AbstractShape}
    shps::Vector{Shapefile.Polygon}
    fips::Vector{UInt64}
    states::Set{UInt64}
    data::FIPSTable{3}
    state_data::FIPSTable{3}
end
data_matrix(g::GeoplotData, var::AbstractString) = g.data[var]
data_matrix(g::GeoplotData) = data_matrix(g, first(keys(g.data.var_index)))
state_data_matrix(g::GeoplotData, var::AbstractString) = g.state_data[var]
n_geo(g::GeoplotData) = size(data_matrix(g), 2)
n_shape(g::GeoplotData) = length(g.shps)
n_timepoint(g::GeoplotData) = size(data_matrix(g), 1)
n_state(g::GeoplotData) = length(g.states)
column_names(g::GeoplotData) = collect(keys(g.data.var_index))
# ============================================================================ #
function geoplot_data(::Type{T}, data_file::AbstractString,
    prefix::AbstractString="", ::Type{S}=UInt32;
    norms::Dict{String,String}=Dict{String,String}()
    ) where{T <:AbstractShape, S <: Integer}

    if endswith(data_file, ".events.bin")
        data = Epicast.read_eventfile(Epicast.RunData, data_file)
    else
        data = Epicast.read_runfile(data_file, S)
    end

    if isempty(prefix)
        names = sort!(Epicast.column_names(data))
    elseif haskey(data.data.index, prefix)
        names = [prefix]
    else
        names = Epicast.filter_columns(x -> startswith(x, prefix), data.data)
    end

    return geoplot_data(T, data, names; norms=norms)
end
# ---------------------------------------------------------------------------- #
function geoplot_data(::Type{T}, data_file::AbstractString, pat::Regex,
        ::Type{S}=UInt32; norms::Dict{String,String}=Dict{String,String}()
    ) where{T <:AbstractShape, S <: Integer}

    if endswith(data_file, ".events.bin")
        data = Epicast.read_eventfile(Epicast.RunData, data_file)
    else
        data = Epicast.read_runfile(data_file, S)
    end

    names = Epicast.filter_columns(x -> match(pat, x) != nothing, data.data)

    return geoplot_data(T, data, names; norms=norms)
end
# ---------------------------------------------------------------------------- #
function geoplot_data(::Type{T}, all_data::Epicast.RunData,
    names::AbstractVector{<:AbstractString}, ::Type{S}=UInt32;
    norms::Dict{String,String}=Dict{String,String}(),
    ) where{T <:AbstractShape, S <: Integer}

    data = data_dict(all_data, names, to_geo(T), norms)
    state_data = data_dict(all_data, names, TRACT2STATE, norms)

    states = Set(all_fips(state_data))

    shps, fips = load_polygons(shape_file(T), states, to_state(T))

    if T == TractPoint
        idx = findall(!in(all_fips(data)), fips)
        deleteat!(shps, idx)
        deleteat!(fips, idx)
    end

    return GeoplotData{T}(shps, fips, states, data, state_data)
end
# ============================================================================ #
function draw_shapes!(ax, data::GeoplotData{CountyPolygon}, var::AbstractString,
    cm::ColorMap, ext::Tuple{Float64,Float64})

    polys = Vector{Matrix{Float64}}(undef, n_shape(data))
    for k in 1:n_shape(data)
        ks, ke = largest_part(data.shps[k])
        pts = view(data.shps[k].points, ks:ke)
        polys[k] = [getfield.(pts, :x) getfield.(pts, :y)]
    end
    c = map(x -> (data.data[x, var][1] - ext[1]) / ext[2], data.fips)
    hp = matplotlib.collections.PolyCollection(polys, facecolor=cm(c))
    ax.add_collection(hp)
    return hp
end
# ============================================================================ #
function draw_shapes!(ax, data::GeoplotData{TractPoint}, var::AbstractString,
    cm::ColorMap, ext::Tuple{Float64,Float64})

    xy = centroid.(data.shps)
    c = map(x -> (data.data[x, var][1] - ext[1]) / ext[2], data.fips)
    return ax.scatter(getindex.(xy, 1), getindex.(xy, 2), 3.0, c=cm(c))
end
quantile_threshold(::Type{CountyPolygon}) = 0.999
quantile_threshold(::Type{TractPoint}) = 0.99
# ============================================================================ #
function make_figure(data::GeoplotData{T}; ofile::AbstractString="",
    style_geo::Function=identity, style_line::Function=identity,
    maxq::Real=quantile_threshold(T), frame::Integer=1, vertical::Bool=false) where T<:AbstractShape

    var = first(keys(data.data.var_index))

    return make_figure(data, var, ofile=ofile, style_geo=style_geo, style_line=style_line,
        maxq=maxq, frame=frame, vertical=vertical)
end
# ---------------------------------------------------------------------------- #
function make_figure(data::GeoplotData{T}, var::AbstractString;
    ofile::AbstractString="", style_geo::Function=identity, style_line::Function=identity,
    maxq::Real=quantile_threshold(T), frame::Integer=1, vertical::Bool=false) where T<:AbstractShape

    state_shp = joinpath(DATADIR, "cb_2019_us_state_500k.shp")

    data_mat = data_matrix(data, var)
    mn = minimum(data_mat)
    mx = quantile(vec(data_mat), maxq)

    cm = PyPlot.get_cmap("viridis")

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

    hp = draw_shapes!(ax[1], data, var, cm, (mn, mx))

    state_outlines!(ax[1], state_shp, data.states, state_color(T))

    ax[1].set_title("Day 0", fontsize=18)

    foreach(ax[1].spines) do k
        ax[1].spines[k].set_visible(false)
    end

    ax[1].set_yticks([])
    ax[1].set_xticks([])
    style_geo(ax)

    colors = map(col -> (red(col), green(col), blue(col)),
        distinguishable_colors(nstate, [RGB(1,1,1), RGB(0,0,0)],
            dropseed=true)
    )

    nt = n_timepoint(data)
    j = 1
    for k in sort!(collect(keys(data.state_data.fips_index)))
        v = data.state_data[k, var]
        ax[2].plot(0:(nt-1), v, linewidth=2, label=STATE_FIPS[k], color=colors[j],
            picker=true)
        j += 1
    end

    ax[2].spines["right"].set_visible(false)
    ax[2].spines["top"].set_visible(false)
    ax[2].set_ylabel("New cases per 100k residents", fontsize=14)
    ax[2].set_xlabel("Simulation day", fontsize=14)
    style_line(ax[2])

    ncol = length(data.states) > 25 ? 2 : 1

    if vertical
        ax[2].legend(frameon=true, loc="lower left", bbox_to_anchor=(1.02, 0.0), ncol=ncol)
    else
        ax[2].legend(frameon=true, loc="upper left", bbox_to_anchor=(1.02, 1.0), ncol=ncol)
    end
    ax[2].set_title(" ", fontsize=18)
    mx2 = maximum(maximum, state_data_matrix(data, var))

    time_idc = ax[2].plot([0,0],[0,mx2], "--", color="darkgray", linewidth=2)[1]

    county_line = nothing

    IDX = 1

    update_figure(idx::Integer) = begin
        # t0 = time()
        c = map(x -> (data.data[x, var][idx] - mn) / mx, data.fips)
        hp.set_facecolors(cm(c))
        ax[1].set_title("Day " * string(idx-1), fontsize=18)
        time_idc.set_xdata([idx-1, idx-1])
        # println("Frame $(idx): ", time() - t0)
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
                idx = l["ind"][1] + 1
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
                ax[2].set_title("County " * string(fips_code), fontsize=18)
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

    (1 < frame <= nt) && update_figure(frame)

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
# ============================================================================ #
function load_polygons(ifile::AbstractString, geo::AbstractSet{<:Integer}, level::Integer=1)
    tbl = Shapefile.Table(ifile)
    fips = parse.(Int, tbl.GEOID)
    idx = findall(x->in(div(x, level), geo), fips)
    return Shapefile.shapes(tbl)[idx], fips[idx]
end
# ============================================================================ #
# /Users/palexander/Documents/emerge+radium/geo-data/cb_2019_us_state_500k/cb_2019_us_state_500k.shp
function state_outlines!(ax, shpfile::AbstractString, states::AbstractSet{<:Integer}, color::AbstractString="white")
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
