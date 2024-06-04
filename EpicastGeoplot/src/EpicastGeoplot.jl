module EpicastGeoplot

using Shapefile, PyPlot, Colors, Printf

using Epicast

using PyCall
const animation = PyNULL()

function __init__()
    copy!(animation, pyimport("matplotlib.animation"))
end

const TRACT2STATE = 10^9
const TRACT2COUNTY = 10^6
const COUNTY2STATE = 10^3

# ============================================================================ #
# function county_data(x::Epicast.RunData, county_fips, var::AbstractString)
#     idx = findall(x -> div(x, TRACT2COUNTY) == county_fips, x.fips)
#     return Epicast.total_cases(view(Epicast.rundata(x, var), idx, :))
# end
function data_dict(data::Epicast.RunData, conv::Integer=TRACT2COUNTY)
    tmp, grp = Epicast.group_by(data, "total", x->div(x, conv),
        Epicast.new_cases, normalize=true)

    return Dict{Int,Vector{Float64}}(grp[k] => tmp[:,1,k] .* 1e5 for k in eachindex(grp))
end
# ============================================================================ #
struct GeoplotData
    shps::Vector{Shapefile.Polygon}
    fips::Vector{UInt64}
    states::Set{UInt64}
    data::Dict{Int,Vector{Float64}}
    state_data::Dict{Int,Vector{Float64}}
end
# ============================================================================ #
function geoplot_data(data_file::AbstractString)

    all_data = Epicast.read_runfile(data_file)

    states = Set(div.(all_data.fips, TRACT2STATE))

    county_shp = "/Users/palexander/Documents/emerge+radium/geo-data/cb_2019_us_all_5m/cb_2019_us_county_5m/cb_2019_us_county_5m.shp"

    shps, fips = load_polygons(county_shp, states, COUNTY2STATE)

    # tmp, grp = Epicast.group_by(all_data, "total", x->div(x, TRACT2COUNTY),
    #     Epicast.new_cases, normalize=true)
    # data = Dict{Int,Vector{Float64}}(x => log.(1.0 .+ county_data(all_data, x, "total")) for x in fips)
    # data = Dict{Int,Vector{Float64}}(grp[k] => tmp[:,1,k] .* 1e5 for k in eachindex(grp))

    county_data = data_dict(all_data, TRACT2COUNTY)

    state_data = data_dict(all_data, TRACT2STATE)


    return GeoplotData(shps, fips, states, county_data, state_data)
end
# ============================================================================ #
function make_figure(data::GeoplotData, ofile::AbstractString="")

    # data_file = "/Users/palexander/Documents/emerge+radium/testing_results/15844097_run_003/15844097_run_003.bin"
    state_shp = "/Users/palexander/Documents/emerge+radium/geo-data/cb_2019_us_state_500k/cb_2019_us_state_500k.shp"

    mn = minimum(minimum, values(data.data))
    mx = maximum(maximum, values(data.data))

    cm = PyPlot.get_cmap("viridis")
    h, ax = subplots(1,2)
    h.set_size_inches((14,7))

    hp = Vector{Vector{PyPlot.PyCall.PyObject}}(undef, length(data.shps))
    for k in 1:length(data.shps)
        d = (data.data[data.fips[k]][1] - mn) / mx
        hp[k] = fill_shape!(ax[1], data.shps[k], facecolor=cm(d), label=string(Int(data.fips[k])))
    end

    state_outlines!(ax[1], state_shp, data.states)

    ax[1].set_title("Day 0", fontsize=18)

    foreach(ax[1].spines) do k
        ax[1].spines[k].set_visible(false)
    end

    ax[1].set_yticks([])
    ax[1].set_xticks([])

    colors = map(col -> (red(col), green(col), blue(col)), 
        distinguishable_colors(length(data.state_data), [RGB(1,1,1), RGB(0,0,0)],
            dropseed=true)
    )

    nt = length(first(data.data).second)
    j = 1
    for (k,v) in data.state_data
        ax[2].plot(0:(nt-1), v, linewidth=2, label=STATE_FIPS[k])
    end

    ax[2].spines["right"].set_visible(false)
    ax[2].spines["top"].set_visible(false)
    ax[2].set_ylabel("New cases per 100k residents", fontsize=14)
    ax[2].set_xlabel("Simulation day", fontsize=14)

    ax[2].legend(frameon=true, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    ax[2].set_title(" ", fontsize=18)
    mx2 = maximum(maximum, values(data.state_data))

    time_idc = ax[2].plot([0,0],[0,mx2], "--", color="darkgray", linewidth=2)[1]

    county_line = nothing

    IDX = 1

    update_figure(idx::Integer) = begin
        for k in 1:length(data.shps)
            d = (data.data[data.fips[k]][idx] - mn) / mx
            update_color!(hp[k], cm(d))
        end
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
        lab = evt.artist.get_label()
        if lab != nothing
            fips_code = parse(Int, lab)
            mx_use = max(mx2, maximum(data.data[fips_code]))
            if county_line == nothing
                county_line = ax[2].plot(0:(nt-1), data.data[fips_code],
                    color="black", linewidth=2.5)[1]
            else
                county_line.set_ydata(data.data[fips_code])
            end
            ax[2].set_ylim(0, mx_use)
            time_idc.set_ydata([0, mx_use])
            ax[2].set_title("County " * lab, fontsize=18)
        end
    end

    h.tight_layout()

    if !isempty(ofile)
        anim = animation.FuncAnimation(h, update_figure, frames=2:nt)
        anim.save(ofile)
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
    return Shapefile.shapes(tbl)[idx], fips[idx]
end
# ============================================================================ #
function update_color!(patches, color::NTuple{4,Float64})
    for patch in patches
        patch.set(facecolor=color)
    end
end
# ============================================================================ #
# /Users/palexander/Documents/emerge+radium/geo-data/cb_2019_us_state_500k/cb_2019_us_state_500k.shp
function state_outlines!(ax, shpfile::AbstractString, states::Set{<:Integer})
    shps, _ = load_polygons(shpfile, states, 1)
    for shp in shps
        fill_shape!(ax, shp, edgecolor="white", facecolor="none")
    end
end
# ============================================================================ #
function fill_shape!(ax, shp::Shapefile.Polygon; facecolor=nothing,
    edgecolor=nothing, label=nothing)

    h = Vector{PyPlot.PyCall.PyObject}(undef, length(shp.parts))
    for k in 1:length(shp.parts)
        ks = shp.parts[k] + 1
        ke = k < length(shp.parts) ? shp.parts[k+1] : length(shp.points)
        h[k] = ax.fill(getfield.(shp.points[ks:ke], :x),
            getfield.(shp.points[ks:ke], :y), facecolor=facecolor,
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
