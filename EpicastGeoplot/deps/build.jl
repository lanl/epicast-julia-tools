using LibGit2, Tar, Inflate

repo_dir = joinpath(@__DIR__, "..", "geo-data")
archive = joinpath(repo_dir, "us_shapefiles.tar.gz")

if !isdir(repo_dir) || !isfile(archive)
    @info("Cloning geodata repository...")
    LibGit2.clone("https://git.lanl.gov/epicast2/geodata.git", repo_dir)
end

datadir = joinpath(repo_dir, "assets")

if !isdir(datadir)
    @info("Extracting data...")
    Tar.extract(IOBuffer(inflate_gzip(archive)), datadir)
end

for name in ["county_5m", "state_5m", "tract_500k", "blockgroup_500k"]
    dbf_file = joinpath(datadir, "cb_2019_us_" * name * ".dbf")
    shp_file = joinpath(datadir, "cb_2019_us_" * name * ".shp")
    @assert(isfile(dbf_file), "failed to locate dbf file at $(dbf_file)")
    @assert(isfile(shp_file), "failed to locate shp file at $(shp_file)")
end