# ---- Build the example files in inst/extdata ----
# The vignettes start from files that look like a user's own: a depth raster,
# a range test and a table of planned stations. All three describe the Stoney
# example of the suite: Stoney Lake, Ontario, at 100 m, a grid of 20
# receivers, and the range model of the arraydesignR app. Run from the
# package root with arraydesign installed:
#   Rscript data-raw/extdata.R

# ---- Libraries ----
library(arraydesign)

# ---- Stoney Lake depth raster ----
# stoney_depth.tif is not written here. It is a byte-identical copy of
# arraydesignR/inst/extdata/stoney_depth.tif, the raster the arraydesignR app
# loads at start. It has the same grid and depths as depth_raster of
# positionR 1.7.0, a 100 m depth raster of Stoney Lake, derived from
# the bathymetry made before 2024, with the CRS set to EPSG:32617. Depth is
# zero or negative metres (0 to -32.5), land NA. Every depth below is taken
# with abs().
stoney_file <- "inst/extdata/stoney_depth.tif"
stopifnot(unname(tools::md5sum(stoney_file)) == "21ba8e62d4f1fe7050b155a8fa6c59b3")
stoney <- terra::rast(stoney_file)
water <- !is.na(terra::values(stoney, mat = FALSE))
sum(water)  # water cells: 4,818
sum(terra::values(stoney, mat = FALSE) == 0, na.rm = TRUE)  # water cells at 0 m: 56
range(terra::values(stoney, mat = FALSE), na.rm = TRUE)  # -32.5 to 0 m

# ---- Toy lake raster ----
# The raster of example_lake(), kept as a GeoTIFF because scripts outside the
# package may read it. Depth in metres, positive down, NA on land.
lake <- example_lake()
tif <- "inst/extdata/example_lake.tif"
terra::writeRaster(lake, tif, overwrite = TRUE, datatype = "FLT4S",
                   gdal = c("COMPRESS=DEFLATE"))
stopifnot(isTRUE(all.equal(terra::values(terra::rast(tif)), terra::values(lake),
                           tolerance = 1e-6)))
file.size(tif)

# ---- Planned stations ----
# The planned array: 20 receivers on a regular grid over the water at least
# 1 m deep, as the app's deployable_grid() and place_receivers() build it.
deployable <- terra::ifel(abs(stoney) < 1, NA, stoney)
grid_20 <- place_receivers(deployable, n = 20, method = "grid", seed = 1)
nrow(grid_20)  # receivers: 20
range(grid_20$deploy_z)  # receiver depth, m: 1.0 to 24.8

## Write and read back
# Latitude and longitude to 7 decimals (about 1 cm). Grid receivers sit at
# cell centres, 50 m from any cell edge, so the rounding keeps every station
# in its cell.
planned <- data.frame(receiver_serial = grid_20$receiver_serial,
                      station_id = grid_20$station_id,
                      deploy_lat = round(grid_20$deploy_lat, 7),
                      deploy_lon = round(grid_20$deploy_lon, 7))
back <- place_receivers(stoney, method = "manual", coords = planned)
max(sqrt((back$x - grid_20$x)^2 + (back$y - grid_20$y)^2))  # largest shift, m
stopifnot(max(abs(back$x - grid_20$x)) < 0.05, max(abs(back$y - grid_20$y)) < 0.05,
          isTRUE(all.equal(back$deploy_z, grid_20$deploy_z)))
write.csv(planned, "inst/extdata/planned_stations.csv", row.names = FALSE)
planned

# ---- Range test ----
# A range test at the shallowest planned receiver, R006 (station S006), in
# 1.0 m of water (deploy_z 1.004 m). At that depth the range model, d50
# 250 m and d05 500 m at 1 m and d50 750 m and d05 1,500 m at 35 m, passes
# through 0.5 at 250.1 m and 0.05 at 500.1 m. One test transmitter at each of
# eight distances for 48 hours, each with the app's delay: uniform between
# 120 and 180 s from the start of one transmission to the start of the next,
# so about 24 transmissions an hour. Detections are binomial, with no spread
# from hour to hour (the range model is static).
#
# A logistic regression of detections on distance fitted to this file gives
# d50 247.9 m (95% CI 242.5 to 253.3) and d05 492.3 m (481.3 to 503.2), from
# 9,205 transmissions (R 4.2.1, seed 20261102). Both intervals hold the
# curve's 250.1 m and 500.1 m at 1.0 m.
# The Pearson dispersion is 1.01.
model <- detfield::de_logistic(d50 = c(250, 750), d05 = c(500, 1500), depth = c(1, 35))
test_rx <- grid_20[which.min(grid_20$deploy_z), ]
test_rx[c("receiver_serial", "station_id", "deploy_z")]
set.seed(20261102)
distances <- c(50, 150, 250, 350, 450, 550, 650, 800)
start <- as.POSIXct("2024-05-20 00:00:00", tz = "UTC")
n_hours <- 48
one_tag <- function(distance_m, k) {
  # Transmission times: a random first delay, then uniform 120 to 180 s.
  gaps <- stats::runif(ceiling(n_hours * 3600 / 120) + 10, 120, 180)
  gaps[1] <- stats::runif(1, 0, gaps[1])
  at <- cumsum(gaps)
  at <- at[at < n_hours * 3600]
  hour <- floor(at / 3600)
  sent <- tabulate(hour + 1, nbins = n_hours)
  p <- stats::plogis(model$predict(distance_m, depth_m = test_rx$deploy_z))
  data.frame(receiver_serial = test_rx$receiver_serial, station_id = test_rx$station_id,
             deploy_z = round(test_rx$deploy_z, 2),
             transmitter = sprintf("A69-1604-%d", 100 + k), distance_m = distance_m,
             hour = format(start + (seq_len(n_hours) - 1) * 3600, "%Y-%m-%d %H:%M:%S"),
             transmissions = sent, detections = stats::rbinom(n_hours, sent, p))
}
range_test <- do.call(rbind, Map(one_tag, distances, seq_along(distances)))
nrow(range_test)  # rows: 8 distances x 48 hours = 384
sum(range_test$transmissions)  # transmissions sent
mean(range_test$transmissions)  # transmissions per hour, about 24

## Fitted curve
fit <- stats::glm(cbind(detections, transmissions - detections) ~ distance_m,
                  family = stats::binomial, data = range_test)
b <- stats::coef(fit)
c(d50 = -b[[1]] / b[[2]], d05 = (stats::qlogis(0.05) - b[[1]]) / b[[2]])  # 247.9, 492.3 m
sum(stats::residuals(fit, type = "pearson")^2) / stats::df.residual(fit)  # dispersion, 1.01
write.csv(range_test, "inst/extdata/range_test.csv", row.names = FALSE)
