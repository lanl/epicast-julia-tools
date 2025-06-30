module PlotHelpers

using PyPlot, PyCall, ColorTypes

export default_axes, axes_layout, axes_label, merge_axes

const RealVec = AbstractVector{<:Real}
# ============================================================================ #
function default_axes(ax=nothing, width=3.0)
    if ax == nothing
        ax =PyPlot.axes()
    end
    #set ticks to face out
    ax.tick_params(direction="out", length=width*2, width=width)

    #turn off top and right axes
    ax.spines["right"].set_visible(false)
    ax.spines["top"].set_visible(false)

    #remove top and right tick marks
    tmp = ax.get_xaxis()
    tmp.tick_bottom()

    tmp = ax.get_yaxis()
    tmp.tick_left()

    ax.spines["left"].set_linewidth(width)
    ax.spines["bottom"].set_linewidth(width)

    return ax
end
# ============================================================================ #
function axes_layout(h::PyPlot.Figure; row_height::RealVec=Float64[], row_spacing::RealVec=Float64[],
    col_width::RealVec=Float64[], col_spacing::RealVec=Float64[])

    nrow = length(row_height)
    ncol = length(col_width)
    ax = Vector{PyCall.PyObject}(undef, nrow*ncol)

    height, btm = axes_position(row_height, spacing=row_spacing, vertical=true)

    for k in eachindex(height)
        width, left = axes_position(col_width, spacing=col_spacing)
        for j in eachindex(width)
            ax[(k-1)*ncol + j] = h.add_axes([left[j], btm[k], width[j], height[k]])
        end
    end
    return ax
end
# ============================================================================ #
function axes_position(width::AbstractVector{<:Real}; pad::Real=0.05, spacing::AbstractVector{<:Real}=fill(pad, length(width)+1), vertical::Bool=false)
    @assert(length(spacing) > length(width), "spacing array does not contain enough items!")
    widthsc = (width ./ sum(width)) .* ((1.0-spacing[end]) - sum(spacing[1:end-1]))
    if vertical
        left = 1 .- cumsum(widthsc .+ spacing[1:end-1])
    else
        left = cumsum(spacing[1:end-1]) .+ vcat(0.0, cumsum(widthsc[1:end-1]))
    end
    return widthsc, left
end
# ============================================================================ #
function axes_label(h, ax, lab::String, x::Real=NaN)

    yl = ax.get_ylim()

    px = matplotlib.transforms.Bbox.from_bounds(10, 10, 60, 10)
    transf = ax.transData.inverted()
    tmp = px.transformed(transf)

    if isnan(x)
        field = :x1
        ytl = ax.get_yticklabels()
        if !isempty(ytl)
            ht = ax.get_yticklabels()[end]
        else
            ht = ax.yaxis.get_label()
            if isempty(ht.get_text())
                ht = ax
                field = :x0
            end
        end
        h.draw_without_rendering()
        bbpx = ht.get_window_extent(renderer=h.canvas.get_renderer())
        bbdata = bbpx.transformed(transf)
        x = getproperty(bbdata, field) - tmp.width
    end

    ax.text(x, yl[2] + tmp.height, lab, fontsize=30, verticalalignment="bottom", horizontalalignment="left")
    return x
end
# ============================================================================ #
function axes_label2(ax, lab::String)

    h = ax.figure
    h.draw_without_rendering()
    bbpx = ax.get_tightbbox(h.canvas.get_renderer())
    bb = bbpx.transformed(h.transFigure.inverted())

    x = bb.x0
    y = bb.y1

    return h.text(x, y, lab, fontsize=30, verticalalignment="bottom", horizontalalignment="left")
end
# ============================================================================ #
function axes_perimiter(ax)
    h = ax.figure
    h.draw_without_rendering()
    bbpx = ax.get_tightbbox(h.canvas.get_renderer())
    return bbpx.transformed(h.transFigure.inverted())
end
# ============================================================================ #
function merge_axes(ax1, ax2)

    pos1 = ax1.get_position()
    pos2 = ax2.get_position()
    l = min(pos1.x0, pos2.x0)
    b = min(pos1.y0, pos2.y0)
    x1 = max(pos1.x1, pos2.x1)
    y1 = max(pos1.y1, pos2.y1)
    pos = [l, b, x1-l, y1-b]

    ax2.remove()
    ax1.set_position(pos)

    return ax1
end
# ============================================================================ #
function plot_with_error(x::RealVec, y::RealVec, yerr::RealVec,
    col::AbstractString, ax=nothing; args...)
    return plot_with_error(x, y, yerr, parse(Colorant, col), ax; args...)
end
# ---------------------------------------------------------------------------- #
function plot_with_error(x::RealVec, y::RealVec, yerr::RealVec,
    col::ColorTypes.RGB, ax=nothing; args...)
    return plot_with_error(x, y, y.-yerr, y.+yerr, col, ax; args...)
end
# ---------------------------------------------------------------------------- #
function plot_with_error(x::RealVec, y::RealVec, yerr::RealVec,
    col::RealVec, ax=nothing; args...)
    return plot_with_error(x, y, y.-yerr, y.+yerr, ColorTypes.RGB(col[1], col[2], col[3]), ax; args...)
end
# ---------------------------------------------------------------------------- #
function plot_with_error(x::RealVec, y::RealVec, ylo::RealVec, yhi::RealVec,
    col::ColorTypes.RGB, ax=nothing; linewidth=4.0, alpha=0.5, label="", ferr=6.0, args...)
    if ax == nothing
        ax = default_axes()
    end
    col_array = [col.r, col.g, col.b]
    fcol, ecol = shading_color(col_array, ferr)
    ax.fill_between(x, ylo, yhi, facecolor=fcol, edgecolor=ecol, alpha=alpha)
    return ax.plot(x, y, "-", color=col_array, linewidth=linewidth, label=label, args...)
end
# ---------------------------------------------------------------------------- #
function shading_color(col::RealVec, ferr::Real=6.0)
    fedge = 0.25
    orig = convert(HSV, RGB(col...))
    hsv = HSV(orig.h, orig.s/ferr, 1.0 - (abs(1.0 - orig.s)^ferr/ferr))
    col_err = hsv2rgb(hsv)
    col_edge = (1.0 - fedge) * col_err + fedge * col
    return col_err, col_edge
end
# ---------------------------------------------------------------------------- #
function hsv2rgb(x::ColorTypes.HSV{T}) where {T<:Number}
    rgb = convert(RGB, x)
    return T[rgb.r, rgb.g, rgb.b]
end
# ============================================================================ #
end # module PlotHelpers
