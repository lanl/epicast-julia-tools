# Dev notes for epicast-utils repo

## Normalizing / preprocessing Epicast run data

```julia
using EpicastGeoplot, Epicast, EpicastTables

ifile = "<path_to_count_or_event_file>"

# load data from <ifile> as a count-based Epicast.RunData object and
# aggregate both the run data and demographic data at the county level
run_data = Epicast.aggregate(County, Epicast.read_rundata(ifile))

# function signature to get the denominator for nomalization, 3 note:
# 1) the full RunData object is passed, so demographics or case data (etc.) can
#    be used for normalization
# 2) if smoothing and/or differencing is requested (see preprocess!() below),
#    that will occur ***before*** get_denom() is called
# 3) the second argument <var> is redundant for the sinlge column syntax, but
#    important for the alternate methods listed below
function get_denom(x::Epicast.RunData, var::AbstractString)
    return Epicast.demographics(x, var) .* 1e-5
end

# preprocess the "total" column by smoothing, differencing (so "new cases"),
# then normalizing by the vector returned by the get_denom() function
# note that only the run data (not demographic the data) is preprocessed
Epicast.preprocess!(run_data, "total", smooth=true, diff=true,
    get_denom=get_denom)

gd = EpicastGeoplot.geoplot_data(run_data)

```

```julia
# preprocess!() also has other methods:

# ---------------------------------------------------------------------------- #
# preprocess column1 and column2 in the same way
Epicast.preprocess!(run_data, ["column1","column2"], smooth=true, diff=true,
    get_denom=(x,_) -> Epicast.demographics(x, "total") .* 1e-5)

# ---------------------------------------------------------------------------- #
# preprocess all columns matching a regex the same way
Epicast.preprocess!(run_data, r"column\d", smooth=true, diff=true,
    get_denom=(x,_) -> Epicast.demographics(x, "total") .* 1e-5)

# ---------------------------------------------------------------------------- #
# preprocess all columns for which the given function return true
fmatch(var::AbstractString) = endswith(var, r"_age\d")

# slightly more complex example of a function to calulate the denominator for
# normalization
function get_age_denom(x::Epicast.RunData, var::AbstractString)
    denom_name = "age_" * var[end]
    return Epicast.demographics(x, denom_name)
end

Epicast.preprocess!(run_data, fmatch, smooth=true, diff=true,
    get_denom=get_age_denom)

# ---------------------------------------------------------------------------- #
# one final example: % of new / daily infections attributable to a given context

function get_src_denom(x::Epicast.RunData, ::AbstractString)
    # NOTE: we have to be careful as to whether the "total" column has been
    # smoothed, diff'd, or normalized, here I'm assuming the "total" column has
    # NOT been preprocessed at all, so we have to diff it to get new cases
    tmp = Epicast.rundata(x, "total")
    tmp[2:end] .= diff(tmp)
    return tmp
end

# NOTE: in this case we are not smoothing the numerator (or denominator)
Epicast.preprocess!(run_data, var->startswith(var, "infection-src"),
    smooth=false, diff=true, get_denom=get_src_denom)

```