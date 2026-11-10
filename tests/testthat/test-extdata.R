# ---- Example files in inst/extdata ----
# The Stoney example of the package: the depth raster of the arraydesignR
# app, the 20-receiver grid of vignette("arraydesign") and a range test at
# its shallowest receiver. Written by data-raw/extdata.R.

stoney_file <- function() system.file("extdata", "stoney_depth.tif", package = "arraydesign")

test_that("the Stoney raster is the app's copy, byte for byte", {
  tif <- stoney_file()
  expect_true(nzchar(tif))
  expect_identical(unname(tools::md5sum(tif)), "21ba8e62d4f1fe7050b155a8fa6c59b3")
  v <- terra::values(terra::rast(tif), mat = FALSE)
  expect_identical(sum(!is.na(v)), 4818L)
  # Depth is zero or negative metres; 56 water cells hold 0.
  expect_true(all(v[!is.na(v)] <= 0))
  expect_identical(sum(v == 0, na.rm = TRUE), 56L)
})

test_that("the toy lake raster is the raster of example_lake()", {
  tif <- system.file("extdata", "example_lake.tif", package = "arraydesign")
  expect_true(nzchar(tif))
  lake <- example_lake()
  from_file <- terra::rast(tif)
  expect_equal(dim(from_file), dim(lake))
  expect_equal(as.vector(terra::ext(from_file)), as.vector(terra::ext(lake)))
  expect_equal(terra::values(from_file, mat = FALSE), terra::values(lake, mat = FALSE),
               tolerance = 1e-6)
})

test_that("the planned stations are the 20-receiver grid", {
  stoney <- terra::rast(stoney_file())
  planned <- utils::read.csv(system.file("extdata", "planned_stations.csv",
                                         package = "arraydesign"))
  rx <- expect_silent(place_receivers(stoney, method = "manual", coords = planned))
  grid <- place_receivers(terra::ifel(abs(stoney) < 1, NA, stoney), n = 20, method = "grid",
                          seed = 1)
  expect_identical(rx$station_id, grid$station_id)
  expect_identical(rx$receiver_serial, grid$receiver_serial)
  # Grid receivers sit at cell centres; the largest shift after the round
  # trip through 7-decimal latitude and longitude is 0.0064 m, as
  # data-raw/extdata.R prints it.
  expect_lt(max(sqrt((rx$x - grid$x)^2 + (rx$y - grid$y)^2)), 0.05)
  expect_equal(rx$deploy_z, grid$deploy_z)
  expect_true(all(rx$deploy_z >= 1))
})

test_that("the range test reproduces the depth-interpolated curve at its receiver's depth", {
  rt <- utils::read.csv(system.file("extdata", "range_test.csv", package = "arraydesign"))
  expect_named(rt, c("receiver_serial", "station_id", "deploy_z", "transmitter", "distance_m",
                     "hour", "transmissions", "detections"))
  expect_true(all(rt$detections >= 0 & rt$detections <= rt$transmissions))
  expect_identical(anyDuplicated(rt[c("distance_m", "hour")]), 0L)
  expect_identical(unique(rt$station_id), "S006")
  fit <- stats::glm(cbind(detections, transmissions - detections) ~ distance_m,
                    family = stats::binomial, data = rt)
  b <- stats::coef(fit)
  v <- stats::vcov(fit)
  at <- function(q) {
    est <- (q - b[[1]]) / b[[2]]
    g <- c(-1 / b[[2]], -(q - b[[1]]) / b[[2]]^2)
    c(est = est, se = sqrt(drop(t(g) %*% v %*% g)))
  }
  d50 <- at(0)
  d05 <- at(stats::qlogis(0.05))
  # 247.9 m and 492.3 m from the file written by data-raw/extdata.R.
  expect_lt(abs(d50[["est"]] - 247.9), 0.05)
  expect_lt(abs(d05[["est"]] - 492.3), 0.05)
  # The curve de_logistic(d50 = c(250, 750), d05 = c(500, 1500),
  # depth = c(1, 35)) at the test depth lies inside both 95% intervals.
  z <- unique(rt$deploy_z)
  truth <- c(250 + (z - 1) * 500 / 34, 500 + (z - 1) * 1000 / 34)
  expect_lt(abs(d50[["est"]] - truth[1]), 1.96 * d50[["se"]])
  expect_lt(abs(d05[["est"]] - truth[2]), 1.96 * d05[["se"]])
})
