# ============================================================================ #
# NOTE: PyPlot is assumed to be available w/in the scope of this file
# ============================================================================ #
@inline function scale(mn::Real, rn::Real, v::Real)
    tmp = (v - mn) / rn
    return min(max(0.0, tmp), 1.0)
end
# ============================================================================ #
abstract type AbstractNorm end
@inline minmax(m::AbstractNorm) = m.mn, m.mx
# ============================================================================ #
struct ExtremaNorm <: AbstractNorm
    mn::Float64
    mx::Float64
end
# ---------------------------------------------------------------------------- #
ExtremaNorm(mn::Real, mx::Real, ::Real) = ExtremaNorm(mn, mx)
ExtremaNorm(x::AbstractVector{<:Real}, ::Real=0) = ExtremaNorm(extrema(x)...)
# ---------------------------------------------------------------------------- #
function scale(m::AbstractNorm, v::AbstractVector{<:Real})
    out = similar(v)
    rn = abs(m.mx - m.mn)
    @inbounds for k in eachindex(out)
        out[k] = scale(m.mn, rn, v[k])
    end
    return out
end
mpl_norm(m::AbstractNorm) = PyPlot.matplotlib.colors.Normalize(vmin=m.mn, vmax=m.mx, clip=true)
# ============================================================================ #
struct TwoSlopeNorm <: AbstractNorm
    mn::Float64
    mx::Float64
    mid::Float64
end
# ---------------------------------------------------------------------------- #
function TwoSlopeNorm(mn::Real, mx::Real)
    mid = mn < 0.0 ? 0.0 : mn + (mx-mn)/2
    return TwoSlopeNorm(mn, mx, mid)
end
# ---------------------------------------------------------------------------- #
function TwoSlopeNorm(x::AbstractVector{<:Real}, mid::Real=NaN)
    mn, mx = extrema(x)
    mid = isnan(mid) ? (mn < 0 ? 0.0 : mn + (mx-mn)/2) : mid
    return TwoSlopeNorm(mn, mid, mx)
end
# ---------------------------------------------------------------------------- #
function scale(m::TwoSlopeNorm, v::AbstractVector{<:Real})
    out = similar(v)
    r1 = abs(m.mid - m.mn)
    r2 = abs(m.mx - m.mid)
    @inbounds for k in eachindex(out)
        out[k] = v[k] < m.mid ? scale(m.mn, r1, v[k]) * 0.5 :
            (scale(m.mid, r2, v[k]) * 0.5) + 0.5
    end
    return out
end
mpl_norm(m::TwoSlopeNorm) = PyPlot.matplotlib.colors.TwoSlopeNorm(m.mid, vmin=m.mn, vmax=m.mx)
# ============================================================================ #
struct SymmetricNorm <: AbstractNorm
    mn::Float64
    mx::Float64
    function SymmetricNorm(mn::Real, mx::Real)
        mnt = Float64(mn)
        mxt = Float64(mx)
        if sign(mn) != sign(mx)
            v = max(abs(mn), abs(mx))
            mnt = -v
            mxt = v
        end
        return new(mnt, mxt)
    end
end
function SymmetricNorm(x::AbstractVector{<:Real}, ::Real=0)
    mn, mx = extrema(x)
    return SymmetricNorm(mn, mx)
end
# ============================================================================ #