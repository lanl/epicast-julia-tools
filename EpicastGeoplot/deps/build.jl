using LibGit2, Tar, Inflate

asset_dir = joinpath(@__DIR__, "..", "assets")
!isdir(asset_dir) && mkpath(asset_dir)

repo_dir = joinpath(asset_dir, "geo-data")
archive = joinpath(repo_dir, "us_shapefile.tar.gz")
datadir = joinpath(repo_dir, "data")

if !isdir(repo_dir) || !isfile(archive)
    LibGit2.clone("https://gitlab.lanl.gov/palexander/epicast-geodata.git", repo_dir)
end

if !isdir(datadir)
    Tar.extract(IOBuffer(inflate_gzip(archive)), datadir)
end

for name in ["county_5m", "state_500k", "tract_500k"]
    dbf_file = joinpath(datadir, "cb_2019_us_" * name * ".dbf")
    shp_file = joinpath(datadir, "cb_2019_us_" * name * ".shp")
    @assert(isfile(dbf_file), "failed to locate dbf file at $(dbf_file)")
    @assert(isfile(shp_file), "failed to locate shp file at $(shp_file)")
end