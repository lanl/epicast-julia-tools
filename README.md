# epicast-utils
A collection of Julia packages for wrangeling data related to the [EpiCast](https://gitlab.lanl.gov/palexander/epicast) epidemic simulator. Released under O4913.

## UrbanPop
Utilities for converting UrbanPop data stored as Apache feather files into the EpiCast internal tract and agent DB formats.

## Workerflow
Utilities for converting ASCII workerflow sparse matricies based on tract "ids" to a binary representation based on FIPS codes.

## Checkpoint
Utilities for working with and verifying EpiCast checkpoint files.

## EpicastGeoplot

### Loading data from a runfile
```julia
using EpicastGeoplot; const EG = EpicastGeoplot

ifile = "<path_to_a_count_based_runfile>"

# ===== load all columns
data = EG.geoplot_data(CountyPolygon, ifile)

# ===== load a single, specific column by name
data = EG.geoplot_data(CountyPolygon, ifile, "total")
data = EG.geoplot_data(CountyPolygon, ifile, "age_4")

# ===== load a group of columns by prefix
data = EG.geoplot_data(CountyPolygon, ifile, "age_")
data = EG.geoplot_data(CountyPolygon, ifile, "status_")

# ===== load a group of columns by regex
# all columns that end in _age4
data = EG.geoplot_data(CountyPolygon, ifile, r".*_age4")
EG.column_names(data)

# all columns that are a age-based breakdown
data = EG.geoplot_data(CountyPolygon, ifile, r"(?:age_\d|.*_age\d)")
EG.column_names(data)
```

***NOTE***: see [EpicastGeoplot.case_count!()](./EpicastGeoplot/src/EpicastGeoplot.jl#L58) for default normalization scheme.

### Plotting geo-animation
```julia
ifile = "<path_to_a_count_based_runfile>"

# load all columns
data = EG.geoplot_data(CountyPolygon, ifile)

# plot total infecetions per 100k residents
h, ax = EG.make_figure(data, "total")

# save animation as an mp4: note <ofile> is a kwarg
h, ax = EG.make_figure(data, "hospitalized_age4",
    ofile = splitext(ifile)[1] * "-county.mp4")
```

# Usage

Basic docs go here...
