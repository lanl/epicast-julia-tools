# UrbanPop -> EpiCast convertion utility

Converts a directory of `.feather` Apache Arrow files to the new EpiCast AgentDB format.

## Asumptions

The convertion utility assumes that all feather files are named as:

`syp_<FIPS code>.feather`

as they are in the UrbanPop v1 data. Modification will likely be needed in the future.

```julia
using UrbanPop

idir = "/home/user/NM_arrow_data"
odir = "/home/user/epicast-data"

# output file will be named or the state FIPS code of the input files, an error
# is thrown if not all files appear to correspond to the same state
ofile = UrbanPop.convert_feather_dir(idir, odir)

# memmap the Agent data
data = UrbanPop.memmap(ofile)

# first agent...
data[1]

```