# epicast-julia-tools
A collection of Julia packages for wrangeling data related to the [EpiCast](https://arxiv.org/abs/2504.03604) epidemic simulator.

## Install

A smoother installation expierence is stil a wip. For now you need to install each sub-module that you would like to use *AFTER* installing it's dependencies. For example:

```julia
import Pkg

# just install everything into the current env
all_submodules = ["EpicastTables", "Epicast", "EpicastPlot", "EpicastGeoplot",
    "PlotHelpers"]

for submod in all_submodules
    Pkg.add("https://github.com/lanl/epicast-julia-tools:" * submod)
end
```

## Sub-modules

### EpicastTables

Basic table-like data structure for storing `time x location x variable` count data as well as 1d and 2d variants. 

### Epicast

Basic I/O and preprocessing routines for data in either count or event format.

### EpicastPlot

Basic plotting of timeseries data.

### PlotHelpers

Basic plotting helper function that are occationally useful. Function provided herein are not specific to Epicast or used by other modules in this repository, but are useful for constructing more complex figures such as those appearing in recent Epicast-related publications.

### EpicastGeoplot

More feature-rich plotting of timeseries and geographic (i.e., map) data.

## Example usage
```julia
using Epicast, EpicastGeoplot; const EG = EpicastGeoplot

ifile = "<path_to_a_count_or_event_file>"

# name of the column to preprocess (we're only going to do one here)
col_name = "total"

# gist: load data, aggregate to County level, apply 7-day moving average
# convert to new counts per day, normalize counts by resident population of each
# county / 100k
data = EG.geoplot_data(
    Epicast.preprocess!(
        # aggregate count data at the county level
        Epicast.aggregate(County, Epicast.read_rundata(ifile)),
        col_name,
        smooth=true, # apply 7-day moving average
        diff=true,   # convert cumulative counts to "new per day"

        # function to "get denominator" for normalizing data
        # x is an Epicast.Rundata object, var is a column name
        get_denom = (x,var) -> Epicast.demographics(x, var) .* 1e-5
    )
)

# make sure the units given by the y-label reflect the normalization performed
# above
EG.make_figure(data, "total", ylab="New cases per 100k residents")

# or to save the animation to an mp4 iff ffmpeg is on the PATH:
# EG.make_figure(data, "total", ofile = "./test.mp4")

```

# Usage

More docs go here...

# Release

This software has been approved for open source release and has been assigned O4913.

# Copyright

© 2025. Triad National Security, LLC. All rights reserved.

This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.

# License

This program is Open-Source under the BSD-3 License.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

