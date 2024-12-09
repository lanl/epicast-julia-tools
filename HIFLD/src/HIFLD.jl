module HIFLD

using CSV

function get_location(ifile::AbstractString, state::AbstractString, field::Symbol=:State)

    lon = field == :STATE ? :LONGITUDE : :Longitude
    lat = field == :STATE ? :LATITUDE : :Latitude

    csv = filter(CSV.File(ifile)) do row
        row[field] == state
    end

    return getproperty.(csv, lon), getproperty.(csv, lat)
end

end # module HIFLD
