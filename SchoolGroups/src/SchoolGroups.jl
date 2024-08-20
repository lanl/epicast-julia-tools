module SchoolGroups

using Statistics, CSV
using UrbanPop
# ============================================================================ #
function parse_number(::Type{T}, x::AbstractString, rep::T=typemax(T)) where T<:Integer
    val = tryparse(T, x)
    return val == nothing ? rep : val
end
# ============================================================================ #
function parse_number(::Type{T}, x::AbstractString, rep::T=T(NaN)) where T<:AbstractFloat
    val = tryparse(T, x)
    return val == nothing ? rep : val
end
# ============================================================================ #
function column_map(x::AbstractString)
    if startswith(x, "Prekinder")
        return 1
    elseif startswith(x, "Kinder")
        return 2
    else
        m = match(r"^Grade\s*(\d+)\s*.*", x)
        m == nothing && error("failed to parse grade level \"$(x)\"")
        n = parse(Int, m[1])
        if n <= 6
            return 3
        elseif n <= 12
            return 4
        else
            return 0
        end
    end
end
# ============================================================================ #
nanadd(a::Float64, b::Float64) = isnan(a) ? b : a + b
# ============================================================================ #
function urbanpop_counties(idir::AbstractString)
    tracts, _ = UrbanPop.all_tract_data(idir)
    return sort!(unique(div.(tracts, 10^6)))
end
function nanmean(data::AbstractVector{T}) where T<:Real
    return sum(x -> isnan(x) ? T(0) : x, data) / count(!isnan, data)
end
# ============================================================================ #
function generate_schoolgroup_file(counties::AbstractVector{<:Integer},
    ofile::AbstractString; max_size::Integer=50)

    idir = "/Users/palexander/Documents/emerge+radium/input-data/school_data"

    student_file = joinpath(idir, "2019_ELSI_student_counts_district.csv")
    teacher_file = joinpath(idir, "2019_ELSI_teacher_counts_district.csv")

    student_data = CSV.File(student_file, header=7, footerskip=7)
    teacher_data = CSV.File(teacher_file, header=7, footerskip=7)

    all_counties = parse_number.(Int, student_data["County Number [District] 2019-20"], 0)
    # counties = sort!(filter!(>(0), unique(all_counties)))

    s_cols = [
        "Prekindergarten Students [District] 2019-20",
        "Kindergarten Students [District] 2019-20",
        "Grade 1 Students [District] 2019-20",
        "Grade 2 Students [District] 2019-20",
        "Grade 3 Students [District] 2019-20",
        "Grade 4 Students [District] 2019-20",
        "Grade 5 Students [District] 2019-20",
        "Grade 6 Students [District] 2019-20",
        "Grade 7 Students [District] 2019-20",
        "Grade 8 Students [District] 2019-20",
        "Grade 9 Students [District] 2019-20",
        "Grade 11 Students [District] 2019-20",
        "Grade 10 Students [District] 2019-20",
        "Grade 12 Students [District] 2019-20"
    ]

    t_cols = [
        "Prekindergarten Teachers [District] 2019-20",
        "Kindergarten Teachers [District] 2019-20",
        "Elementary Teachers [District] 2019-20",
        "Secondary Teachers [District] 2019-20"
    ]

    col_idx = map(column_map, s_cols)

    pt_ratio = parse_number.(Float64,
        teacher_data["Pupil/Teacher Ratio [District] 2019-20"], NaN)

    # county x grade matrix of school group sizes
    out = fill(NaN, length(counties), 4)

    for k in eachindex(counties)
        idx = findall(isequal(counties[k]), all_counties)
        if length(idx) > 0

            pt_ratio_cnty = nanmean(pt_ratio[idx])
            
            for j in eachindex(s_cols)
                # total # of student in each county of each age / grade group
                tmp = sum(x -> parse_number(Float64, x, 0.0),
                    student_data[s_cols[j]][idx])
                out[k, col_idx[j]] = nanadd(out[k, col_idx[j]], tmp)
            end
            
            for j in eachindex(t_cols)
                # # of students per teacher
                nt = sum(x -> parse_number(Float64, x, 0.0), teacher_data[t_cols[j]][idx])

                out[k, j] /= nt

                # if student/teacher is not valid, or > max_size, attempt to
                # replace that data with the county average student/teacher
                if (!isfinite(out[k, j]) || out[k, j] > max_size) && isfinite(pt_ratio_cnty)
                    out[k,j] = pt_ratio_cnty
                end
            end
            
        end
    end

    # interpolate as needed between states via simple mean over all counties
    # for that school type
    if any(!isfinite, out)
        for k in 1:size(out, 2)
            idx = findall(!isfinite, out[:,k])
            other = setdiff(1:size(out,1), idx)
            out[idx,k] .= mean(out[other,k])
        end
    end

    data = permutedims(round.(UInt16, out), (2,1))
    replace!(x -> min(x, max_size), data)

    open(ofile, "w") do io
        write(io, UInt64(size(data,2)))
        write(io, UInt64(size(data,1)))
        write(io, UInt16.(counties))
        write(io, data)
    end

    return data
end
# ============================================================================ #
"""
`counties, data = read_schoolgroup_data(ifile::AbstractString)`

* `counties` - vector of county FIPS codes
* `data` - a county x school type (i.e. 4) matrix of school group sizes,
           order is [pre_k, kindergarten, elementary, secondary]
"""
function read_schoolgroup_data(ifile::AbstractString)
    return open(ifile, "r") do io
        nc = read(io, UInt64)
        nr = read(io, UInt64)
        counties = Vector{UInt16}(undef, nc)
        read!(io, counties)
        data = Matrix{UInt16}(undef, nr, nc)
        read!(io, data)
        return counties, permutedims(data, (2,1))
    end
end
# ============================================================================ #
end # module SchoolGroups
