# Example files

`stoney_depth.tif` has the same grid and depths as `depth_raster` of
positionR 1.7.0, a 100 m depth raster of Stoney Lake, Ontario,
derived from the bathymetry made before 2024, with the CRS set to EPSG:32617
(UTM zone 17N on WGS 84). Depth is zero or negative metres and land is `NA`.
4,818 cells are water, and 56 of them hold 0. The file is the same file as
`stoney_depth.tif` in detfield and arraydesignR.

Every other file is simulated or synthetic. The script that writes them,
<https://github.com/jakebrownscombe/arraydesign/blob/main/data-raw/extdata.R>,
is in the repository and not in the installed package. `example_lake.tif` is the toy lake of
`example_lake()`. `planned_stations.csv` is the 20-receiver grid of
`vignette("arraydesign")`, as latitude and longitude. `range_test.csv` is a
simulated range test at the shallowest of those receivers. The receivers,
stations and transmitters in them are invented.
