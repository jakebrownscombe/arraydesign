# Cell of each receiver by terra, independent of the package's own arithmetic.
on_water <- function(grid, rx) {
  !is.na(terra::extract(grid, cbind(rx$x, rx$y))[, 1])
}

test_that("grid and random placement return exactly n, all on water, on the toy lake", {
  lake <- example_lake()
  for (n in c(5, 48, 50, 200)) {
    g <- place_receivers(lake, n = n)
    r <- place_receivers(lake, n = n, method = "random", seed = n)
    expect_equal(nrow(g), n)
    expect_equal(nrow(r), n)
    expect_true(all(on_water(lake, g)))
    expect_true(all(on_water(lake, r)))
    expect_equal(anyDuplicated(r[c("x", "y")]), 0L)
  }
})

test_that("grid and random placement return exactly n on a Stoney-sized raster", {
  big <- stoney_sized()
  expect_equal(sum(!is.na(terra::values(big))), 45974)
  for (n in c(5, 48, 50, 200)) {
    g <- place_receivers(big, n = n)
    r <- place_receivers(big, n = n, method = "random", seed = 1)
    expect_equal(nrow(g), n)
    expect_equal(nrow(r), n)
    expect_true(all(on_water(big, g)))
    expect_true(all(on_water(big, r)))
  }
})

test_that("the output follows the table contract", {
  rx <- place_receivers(example_lake(), n = 12)
  expect_s3_class(rx, "data.frame")
  expect_s3_class(rx, "receiver_layout")
  expect_named(rx, c("receiver_serial", "station_id", "x", "y", "deploy_lat", "deploy_lon",
                     "deploy_z"))
  expect_type(rx$receiver_serial, "character")
  expect_type(rx$station_id, "character")
  expect_equal(rx$receiver_serial[1:2], c("R001", "R002"))
  expect_equal(attr(rx, "epsg"), 32617L)
  expect_true(all(rx$deploy_lat > 44.5 & rx$deploy_lat < 44.6))
  expect_true(all(rx$deploy_lon > -78.2 & rx$deploy_lon < -78.1))
  # deploy_z is the depth of the cell, positive down.
  expect_equal(rx$deploy_z, terra::extract(example_lake(), cbind(rx$x, rx$y))[, 1])
  expect_true(all(rx$deploy_z > 0))
  wide <- place_receivers(example_lake(), n = 1200, method = "random", seed = 1)
  expect_equal(wide$receiver_serial[1], "R0001")
})

test_that("a grid at a given spacing is a square lattice on water", {
  lake <- example_lake()
  rx <- place_receivers(lake, spacing = 300)
  expect_true(all(on_water(lake, rx)))
  d <- as.matrix(stats::dist(rx[c("x", "y")]))
  diag(d) <- Inf
  expect_equal(unname(apply(d, 1, min)), rep(300, nrow(rx)), tolerance = 1e-8)
  expect_equal(attr(rx, "spacing"), 300)
  expect_error(place_receivers(lake, spacing = 1e5), "No lattice node")
})

test_that("a grid with n drops surplus nodes only when no lattice fits exactly", {
  lake <- example_lake()
  trimmed <- vapply(1:40, function(n) attr(place_receivers(lake, n = n), "trimmed"), integer(1))
  expect_true(all(trimmed >= 0))
  expect_true(any(trimmed == 0))
  expect_true(any(trimmed > 0))
  # Every water cell holds a receiver at the 30 m lattice, one per cell.
  all_cells <- place_receivers(lake, n = 2152)
  expect_equal(anyDuplicated(terra::cellFromXY(lake, cbind(all_cells$x, all_cells$y))), 0L)
  # At spacings of a cell or more, no lattice holds more nodes than there are
  # water cells.
  wg <- arraydesign:::water_grid(lake)
  expect_error(arraydesign:::grid_exact_n(wg, 2153), "No square lattice")
})

test_that("random placement keeps a minimum spacing, and says when n do not fit", {
  lake <- example_lake()
  rx <- place_receivers(lake, n = 20, spacing = 250, method = "random", seed = 3)
  expect_equal(nrow(rx), 20)
  expect_gte(min(stats::dist(rx[c("x", "y")])), 250)
  expect_error(place_receivers(lake, n = 50, spacing = 600, method = "random", seed = 3),
               "receivers fit at a spacing of 600")
})

test_that("manual placement keeps water points, snaps shore points and refuses others", {
  lake <- example_lake()
  rx <- suppressMessages(place_receivers(
    lake, method = "manual",
    coords = data.frame(x = c(725300, 725880), y = c(4935600, 4935600), station_id = c("A", "B"))
  ))
  expect_equal(rx$station_id, c("A", "B"))
  expect_equal(rx$receiver_serial, c("R001", "R002"))
  expect_equal(attr(rx, "snapped"), 1L)
  expect_true(all(on_water(lake, rx)))
  expect_equal(rx$x[1], 725300)
  # x = 725880 sits on the peninsula (x 725870 to 725930); it moves 29 m.
  expect_equal(rx$x[2], 725855)
  expect_message(place_receivers(lake, method = "manual",
                                 coords = data.frame(x = 725880, y = 4935600)),
                 "1 of 1 points moved off land")
  # The outer corner of the north-west land cell is 62 m from the nearest
  # water cell centre but 41 m from its edge, within the 42 m cell diagonal
  # that detfield::water_distances() allows, so it moves.
  corner <- data.frame(x = c(725300, 725001), y = c(4935600, 4936199))
  expect_message(place_receivers(lake, method = "manual", coords = corner), "moved off land")
  # Beyond the lake, 100 m past its edge, a point is refused.
  wide <- terra::extend(lake, 5)
  far <- data.frame(x = c(725300, terra::xmin(lake) - 100), y = c(4935600, 4936199))
  expect_error(place_receivers(wide, method = "manual", coords = far), "row 2 is on land")
  # Every point far on land is counted, and at most 10 are named
  far12 <- data.frame(x = c(725300, rep(terra::xmin(lake) - 100, 12)),
                      y = c(4935600, rep(4936199, 12)))
  expect_error(place_receivers(wide, method = "manual", coords = far12),
               paste("`coords` has 12 rows on land more than one cell diagonal (42 m) from the",
                     "edge of the nearest water cell: 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 and 2 more."),
               fixed = TRUE, class = "arraydesign_error_placement")
  expect_error(place_receivers(lake, method = "manual", coords = data.frame(x = 1, y = 1)),
               "has 1 row outside the raster: 1.", fixed = TRUE,
               class = "arraydesign_error_placement")
  # Hundreds of rows off the raster are counted, and 10 are named
  off <- data.frame(x = c(725300, rep(1, 316)), y = c(4935600, rep(1, 316)))
  expect_error(place_receivers(lake, method = "manual", coords = off),
               paste("`coords` has 316 rows outside the raster:",
                     "2, 3, 4, 5, 6, 7, 8, 9, 10, 11 and 306 more."),
               fixed = TRUE, class = "arraydesign_error_placement")
  # Points off the raster are refused first: with one far on land as well,
  # the error names the point off the raster only.
  mixed <- data.frame(x = c(725300, terra::xmin(lake) - 100, 1),
                      y = c(4935600, 4936199, 1))
  expect_error(place_receivers(wide, method = "manual", coords = mixed),
               "`coords` has 1 row outside the raster: 3.", fixed = TRUE,
               class = "arraydesign_error_placement")
  expect_error(place_receivers(wide, method = "manual", coords = mixed[1:2, ]),
               "`coords` row 2 is on land", fixed = TRUE,
               class = "arraydesign_error_placement")
  # Both placement refusals name the user's call
  err <- tryCatch(place_receivers(wide, method = "manual", coords = mixed), error = identity)
  expect_identical(conditionCall(err),
                   quote(place_receivers(wide, method = "manual", coords = mixed)))
  err <- tryCatch(place_receivers(wide, method = "manual", coords = far12), error = identity)
  expect_identical(conditionCall(err),
                   quote(place_receivers(wide, method = "manual", coords = far12)))
  ll <- data.frame(deploy_lon = rx$deploy_lon, deploy_lat = rx$deploy_lat)
  back <- place_receivers(lake, method = "manual", coords = ll)
  expect_equal(back$x, rx$x, tolerance = 1e-6)
  expect_error(place_receivers(lake, method = "manual", coords = data.frame(a = 1)),
               "needs columns")
  expect_error(place_receivers(lake, method = "manual", coords = data.frame(x = NA, y = 1)),
               "missing coordinates, the first in row 1")
  typo <- data.frame(x = c("725900", "7259OO"), y = c("4935600", "4935600"))
  expect_error(place_receivers(lake, method = "manual", coords = typo),
               paste("`coords\\$x` holds values that are not numbers: 1 of 2, the first",
                     "\"7259OO\" in row 2"),
               class = "arraydesign_error_input")
  typo_ll <- data.frame(deploy_lon = c(-78.1, -78.1), deploy_lat = c("44.5", "44,5"))
  expect_error(place_receivers(lake, method = "manual", coords = typo_ll),
               "`coords\\$deploy_lat` holds values that are not numbers")
})

test_that("manual placement warns on receivers sharing a cell and refuses repeated serials", {
  lake <- example_lake()
  # Two peninsula points that snap to the same water cell centre.
  snapped <- data.frame(x = c(725880, 725882), y = c(4935600, 4935601))
  expect_warning(rx <- suppressMessages(place_receivers(lake, method = "manual",
                                                        coords = snapped)),
                 "`coords` has 2 rows that share a grid cell with another row", fixed = TRUE)
  expect_equal(anyDuplicated(terra::cellFromXY(lake, cbind(rx$x, rx$y))), 2L)
  twice <- data.frame(x = c(725300, 725300, 725400), y = c(4935600, 4935600, 4935600))
  expect_warning(place_receivers(lake, method = "manual", coords = twice),
                 "share every detection probability: 1, 2.", fixed = TRUE)
  # A long list of shared rows is counted, and at most 10 are named
  eleven <- data.frame(x = rep(725300, 11), y = rep(4935600, 11))
  expect_warning(place_receivers(lake, method = "manual", coords = eleven),
                 "11 rows that share a grid cell .*: 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 and 1 more\\.$")
  expect_silent(place_receivers(lake, method = "manual",
                                coords = data.frame(x = c(725300, 725330), y = 4935600)))
  dup <- data.frame(x = c(725300, 725400), y = 4935600, receiver_serial = c("A", "A"))
  expect_error(place_receivers(lake, method = "manual", coords = dup), "repeats: A")
})

test_that("a depth-masked raster keeps receivers out of shallow water", {
  lake <- example_lake()
  deep <- terra::ifel(lake < 5, NA, lake)
  for (method in c("grid", "random")) {
    rx <- place_receivers(deep, n = 20, method = method, seed = 1)
    expect_equal(nrow(rx), 20L)
    expect_gte(min(rx$deploy_z), 5)
  }
  expect_lt(min(place_receivers(lake, n = 20)$deploy_z), 5)
})

test_that("bad arguments stop with a message", {
  lake <- example_lake()
  expect_error(place_receivers(lake, n = 3000), "only 2152 water cells")
  expect_error(place_receivers(lake, n = 2.5), "whole number")
  expect_error(place_receivers(lake, n = 0), "whole number")
  expect_error(place_receivers(lake), "needs `n` or `spacing`")
  expect_error(place_receivers(lake, method = "random", spacing = 100), "needs `n`")
  expect_error(place_receivers(lake, spacing = -1), "at least the cell size")
  expect_error(place_receivers(lake, spacing = 20), "at least the cell size")
  expect_error(place_receivers(lake, n = 3, method = "manual",
                               coords = data.frame(x = 1, y = 1)), "not used with")
  expect_error(place_receivers(lake, method = "manual"), "needs `coords`")
  expect_error(place_receivers(lake, n = 3, method = "nonsense"), "should be one of",
               class = "arraydesign_error_input")
  for (bad in list(c("grid", "random"), NA, 1)) {
    expect_error(place_receivers(lake, n = 3, method = bad), "must be one string",
                 class = "arraydesign_error_input")
  }
  expect_identical(nrow(place_receivers(lake, n = 3, method = "rand", seed = 1)), 3L)
})

test_that("rasters that cannot hold receivers are refused", {
  lake <- example_lake()
  expect_error(place_receivers(matrix(1, 2, 2), n = 1), "SpatRaster")
  expect_error(place_receivers(c(lake, lake), n = 1), "one layer")
  ll <- terra::rast(nrows = 10, ncols = 10, xmin = -78, xmax = -77, ymin = 44, ymax = 45)
  terra::values(ll) <- 1
  expect_error(place_receivers(ll, n = 1), "`grid` is in longitude and latitude")
  # The refusal is classed and names the user's call
  err <- tryCatch(place_receivers(ll, n = 1), error = identity)
  expect_s3_class(err, "arraydesign_error_grid")
  expect_identical(conditionCall(err), quote(place_receivers(ll, n = 1)))
  expect_match(conditionMessage(err), "Project it to a CRS in metres first")
  flat <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 100, ymin = 0, ymax = 200,
                      crs = "EPSG:32617")
  terra::values(flat) <- 1
  expect_error(place_receivers(flat, n = 1), "square")
  dry <- lake
  terra::values(dry) <- NA
  expect_error(place_receivers(dry, n = 1), "no water cells")
  nocrs <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 100, ymin = 0, ymax = 100,
                       crs = "")
  terra::values(nocrs) <- 1
  expect_error(place_receivers(nocrs, n = 1), "no CRS")
})

test_that("a raster whose CRS unit is not metres is refused, naming the call", {
  lake <- example_lake()
  feet <- terra::project(lake, "EPSG:2263", method = "near")
  err <- tryCatch(place_receivers(feet, spacing = 500), error = identity)
  expect_s3_class(err, "arraydesign_error_grid")
  expect_identical(conditionCall(err), quote(place_receivers(feet, spacing = 500)))
  expect_match(conditionMessage(err), "unit is the US survey foot (0.3048006 m)", fixed = TRUE)
  expect_match(conditionMessage(err), "Project it to a CRS in metres first")
  km <- terra::project(lake, "+proj=utm +zone=17 +datum=WGS84 +units=km", method = "near")
  expect_error(place_receivers(km, n = 3), "the kilometre (1000 m)", fixed = TRUE,
               class = "arraydesign_error_grid")
  expect_error(simulate_tracks(feet, crw_model(30, 15, 40), n_animals = 1, n_steps = 2),
               "US survey foot", class = "arraydesign_error_grid")
  # The same lake in metres is placed
  metres <- terra::project(lake, "EPSG:32617", method = "near")
  expect_identical(nrow(place_receivers(metres, n = 5, seed = 1)), 5L)
})

test_that("each non-metre CRS form is refused with the name of its axis unit", {
  lake <- example_lake()
  utm <- "+proj=utm +zone=17 +datum=WGS84 +units="
  tmerc <- paste("+proj=tmerc +lat_0=40.1666666666667 +lon_0=-74 +k=0.999933333 +x_0=300000",
                 "+y_0=0 +ellps=GRS80")
  forms <- list(
    c("EPSG:2263", "the US survey foot (0.3048006 m)"),
    c("EPSG:2240", "the US survey foot (0.3048006 m)"),
    c(paste0(utm, "km"), "the kilometre (1000 m)"),
    c(paste0(utm, "ft"), "the foot (0.3048 m)"),
    c(paste0(utm, "us-ft"), "the US survey foot (0.3048006 m)"),
    c(paste0(utm, "cm"), "the centimetre (0.01 m)"),
    # A BOUNDCRS also holds the metre of its target ellipsoid
    c(paste(tmerc, "+towgs84=0,0,0,0,0,0,0 +units=us-ft"), "the US survey foot (0.3048006 m)"),
    c("+proj=utm +zone=17 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=us-ft",
      "the US survey foot (0.3048006 m)")
  )
  for (f in forms) {
    g <- terra::project(lake, f[1], method = "near")
    expect_error(place_receivers(g, n = 1), paste0("unit is ", f[2], "."), fixed = TRUE,
                 class = "arraydesign_error_grid", info = f[1])
  }
})

test_that("a CRS whose scale at the raster centre is far from 1 is refused, naming the UTM zone", {
  stoney <- terra::rast(system.file("extdata", "stoney_depth.tif", package = "arraydesign"))
  lake <- example_lake()
  for (g in list(stoney, lake)) {
    mercator <- terra::project(g, "EPSG:3857", method = "near")
    expect_gt(min(crs_scale(mercator)$scale), 1.39)
    expect_lt(max(crs_scale(mercator)$scale), 1.41)
    err <- tryCatch(place_receivers(mercator, spacing = 500), error = identity)
    expect_s3_class(err, "arraydesign_error_grid")
    expect_identical(conditionCall(err), quote(place_receivers(mercator, spacing = 500)))
    expect_match(conditionMessage(err), "scale at the raster centre is 1.401 east-west, 1.405",
                 fixed = TRUE)
    expect_match(conditionMessage(err), paste("1.405 north-south, 1.403 south-west to north-east",
                                              "and 1.403 north-west to south-east, more than 1",
                                              "percent from 1"), fixed = TRUE)
    expect_match(conditionMessage(err), "UTM zone of the centre, 17N (EPSG:32617).", fixed = TRUE)
    expect_error(simulate_tracks(mercator, crw_model(30, 15, 40), n_animals = 1, n_steps = 2),
                 "scale at the raster centre", class = "arraydesign_error_grid")
    # Feet are refused for their unit before their scale is measured
    feet <- terra::project(g, "EPSG:2263", method = "near")
    expect_error(place_receivers(feet, n = 1), "unit is the US survey foot", fixed = TRUE,
                 class = "arraydesign_error_grid")
    # UTM 17N, the zone the refusal names, is placed
    utm <- terra::project(g, "EPSG:32617", method = "near")
    expect_identical(nrow(place_receivers(utm, n = 5, seed = 1)), 5L)
  }
  expect_identical(nrow(place_receivers(stoney, spacing = 500)), 196L)
})

# A 10 km raster of 100 m water cells in `crs`, centred on a longitude and
# latitude.
site_grid <- function(crs, lon, lat) {
  centre <- terra::project(cbind(lon, lat), from = "EPSG:4326", to = crs)
  g <- terra::rast(terra::ext(centre[1] - 5000, centre[1] + 5000, centre[2] - 5000,
                              centre[2] + 5000), resolution = 100, crs = crs)
  terra::values(g) <- 10
  g
}

test_that("the scale is measured on both axes and both diagonals", {
  # Sinusoidal off its central meridian: both axes within 1 percent of 1, the
  # grid sheared, so the diagonals are not.
  sinusoidal <- site_grid("ESRI:54008", 11, 45.6)
  s <- crs_scale(sinusoidal)$scale
  expect_named(s, c("east-west", "north-south", "south-west to north-east",
                    "north-west to south-east"))
  expect_true(all(abs(s[1:2] - 1) <= 0.01))
  expect_equal(unname(s[3:4]), c(0.934, 1.071), tolerance = 1e-3)
  err <- tryCatch(place_receivers(sinusoidal, spacing = 500), error = identity)
  expect_s3_class(err, "arraydesign_error_grid")
  expect_match(conditionMessage(err), "0.934 south-west to north-east and 1.071", fixed = TRUE)
  expect_match(conditionMessage(err), "UTM zone of the centre, 32N (EPSG:32632).", fixed = TRUE)
  expect_error(place_receivers(site_grid("ESRI:54008", 6.5, 46.4), spacing = 500),
               "0.960 south-west to north-east and 1.042", class = "arraydesign_error_grid")
  # On its own central meridian the same projection keeps all four at 1
  centred <- site_grid("+proj=sinu +lon_0=11 +datum=WGS84 +units=m", 11, 45.6)
  expect_equal(unname(crs_scale(centred)$scale), rep(1, 4), tolerance = 1e-4)
  expect_identical(nrow(place_receivers(centred, spacing = 500)), 400L)
  # A conformal CRS in its own zone, and an equal-area CRS within 1 percent, pass
  expect_identical(nrow(place_receivers(site_grid("EPSG:32632", 11, 45.6), spacing = 500)), 400L)
  expect_identical(nrow(place_receivers(site_grid("EPSG:3035", 11, 45.6), spacing = 500)), 400L)
  # Web Mercator at the equator: 1.000 east-west and 1.007 north-south
  equator <- crs_scale(site_grid("EPSG:3857", 33, 0))$scale
  expect_equal(unname(equator[1:2]), c(1, 1.0067), tolerance = 1e-4)
  # Near the equator Web Mercator is within 1 percent, so it is placed
  expect_identical(nrow(place_receivers(site_grid("EPSG:3857", 33, -1), spacing = 500)), 400L)
  # A scale that prints alike on all four is given once
  expect_error(place_receivers(terra::project(example_lake(), "EPSG:3395", method = "near"), n = 1),
               "scale at the raster centre is 1.401 on its two axes and its diagonals, more",
               fixed = TRUE, class = "arraydesign_error_grid")
})

test_that("a CRS whose scale cannot be measured is refused with no other condition", {
  eng_wkt <- paste0('ENGCRS["local",EDATUM["site"],CS[Cartesian,2],',
                    'AXIS["x",east,LENGTHUNIT["metre",1]],',
                    'AXIS["y",north,LENGTHUNIT["metre",1]]]')
  stoney <- terra::rast(system.file("extdata", "stoney_depth.tif", package = "arraydesign"))
  for (g in list(example_lake(), stoney)) {
    eng <- g
    terra::crs(eng) <- eng_wkt
    expect_true(all(is.na(crs_scale(eng)$scale)))
    warned <- 0
    err <- withCallingHandlers(
      tryCatch(place_receivers(eng, n = 2, seed = 1), error = identity),
      warning = function(w) {
        warned <<- warned + 1
        invokeRestart("muffleWarning")
      }
    )
    expect_identical(warned, 0)
    expect_s3_class(err, "arraydesign_error_grid")
    expect_identical(conditionCall(err), quote(place_receivers(eng, n = 2, seed = 1)))
    expect_identical(conditionMessage(err), paste(
      "`grid` is refused. The scale of this CRS at the raster centre cannot be measured, so",
      "give the raster in a projected CRS with a geographic datum, such as its UTM zone."
    ))
    err <- tryCatch(simulate_tracks(eng, crw_model(30, 15, 40), n_animals = 1, n_steps = 2),
                    error = identity)
    expect_s3_class(err, "arraydesign_error_grid")
    expect_match(conditionMessage(err), "cannot be measured", fixed = TRUE)
    expect_identical(conditionCall(err), quote(simulate_tracks(eng, crw_model(30, 15, 40),
                                                               n_animals = 1, n_steps = 2)))
  }
})

test_that("the UTM zone of a point covers both hemispheres and both ends of longitude", {
  expect_identical(utm_epsg(-180, 1), 32601L)
  expect_identical(utm_epsg(179.9, 1), 32660L)
  expect_identical(utm_epsg(180, -1), 32760L)
  expect_identical(utm_epsg(0, 0), 32631L)
})

test_that("depth stored as negative numbers is returned positive down", {
  lake <- example_lake()
  neg <- -lake
  rx <- place_receivers(neg, n = 5)
  expect_equal(rx$deploy_z, place_receivers(lake, n = 5)$deploy_z)
  mixed <- lake
  terra::values(mixed)[100] <- -5
  expect_error(place_receivers(mixed, n = 5), "`grid` values mix positive and negative")
  expect_error(place_receivers(mixed, n = 5), "x[x >= 0] <- NA, then x <- -x", fixed = TRUE)
})

test_that("a seed reproduces a random layout and leaves the caller's stream alone", {
  lake <- example_lake()
  set.seed(42)
  before <- .Random.seed
  a <- place_receivers(lake, n = 10, method = "random", seed = 7)
  expect_identical(.Random.seed, before)
  b <- place_receivers(lake, n = 10, method = "random", seed = 7)
  expect_identical(a, b)
  expect_false(identical(a$x, place_receivers(lake, n = 10, method = "random", seed = 8)$x))
  expect_error(place_receivers(lake, n = 3, method = "random", seed = "a"), "one number")
})

test_that("packed and raster package layers are accepted", {
  lake <- example_lake()
  expect_equal(nrow(place_receivers(terra::wrap(lake), n = 5)), 5)
  skip_if_not_installed("raster")
  expect_equal(nrow(place_receivers(raster::raster(lake), n = 5)), 5)
})

test_that("manual placement takes and refuses the points detfield does", {
  lake <- terra::extend(example_lake(), 5)
  res <- terra::res(lake)[1]
  water <- which(!is.na(terra::values(lake, mat = FALSE)))
  xy <- terra::xyFromCell(lake, water)
  west <- water[which.min(xy[, 1])]
  c0 <- terra::xyFromCell(lake, west)
  # Points west of the westmost water cell centre, from 1 to 4 cells away.
  offs <- c(30, 50, 60, 100)
  pts <- data.frame(station_id = paste0("P", seq_along(offs)), x = unname(c0[1, 1]) - offs,
                    y = unname(c0[1, 2]))
  ad <- vapply(seq_len(nrow(pts)), function(i) {
    !inherits(try(suppressMessages(place_receivers(lake, method = "manual",
                                                   coords = pts[i, c("x", "y")])),
                  silent = TRUE), "try-error")
  }, logical(1))
  df <- vapply(seq_len(nrow(pts)), function(i) {
    !inherits(try(suppressWarnings(detfield::water_distances(lake, pts[i, ])), silent = TRUE),
              "try-error")
  }, logical(1))
  expect_identical(ad, df)
  expect_true(any(ad) && !all(ad))
  expect_equal(res, 30)
})

test_that("a layout written to CSV and read back keeps its cells and depths", {
  # Grid and random receivers sit at cell centres. A receiver on a cell edge
  # can move to the next cell when its coordinates are rounded on the way
  # through a CSV.
  stoney <- terra::rast(system.file("extdata", "stoney_depth.tif", package = "arraydesign"))
  deployable <- terra::ifel(abs(stoney) < 1, NA, stoney)
  wg <- arraydesign:::water_grid(stoney)
  layouts <- list(grid = place_receivers(deployable, n = 20, method = "grid", seed = 1),
                  random = place_receivers(deployable, n = 20, method = "random", seed = 1))
  for (nm in names(layouts)) {
    rx <- layouts[[nm]]
    centre <- terra::xyFromCell(stoney, terra::cellFromXY(stoney, cbind(rx$x, rx$y)))
    expect_equal(unname(cbind(rx$x, rx$y)), unname(centre), info = nm)
    f <- withr::local_tempfile(fileext = ".csv")
    utils::write.csv(as.data.frame(rx)[c("station_id", "deploy_lat", "deploy_lon")], f,
                     row.names = FALSE)
    back <- place_receivers(stoney, method = "manual", coords = utils::read.csv(f))
    expect_identical(arraydesign:::cell_at(wg, back$x, back$y),
                     arraydesign:::cell_at(wg, rx$x, rx$y), info = nm)
    expect_equal(back$deploy_z, rx$deploy_z, info = nm)
  }
})
