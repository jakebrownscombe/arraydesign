test_that("no position and no step crosses a 30 m land wall in 20,000 steps per step size", {
  skip_if_not_installed("sf")
  wall <- wall_lake()
  land <- land_polygons(wall)
  wall_x <- 700000 + 49 * 30
  start <- data.frame(x = 700000 + 25 * 30, y = 4900000 + 30 * 30)
  for (mean_step in c(30, 60, 90)) {
    tr <- simulate_tracks(wall, crw_model(mean_step, mean_step / 2, 45), n_animals = 20,
                          n_steps = 1000, control = track_control(start = start),
                          seed = mean_step)
    expect_equal(sum(attr(tr, "escapes")[, -1]), 20000)
    on_land <- is.na(terra::extract(wall, cbind(tr$x, tr$y))[, 1])
    expect_equal(sum(on_land), 0)
    crossing <- lengths(sf::st_intersects(step_segments(tr), land)) > 0
    expect_equal(sum(crossing), 0)
    expect_equal(sum(tr$x > wall_x), 0)
  }
})

test_that("the segment test catches a corner clip and passes open water", {
  wg <- arraydesign:::water_grid(wall_lake())
  x_wall <- 700000 + 49 * 30
  y_mid <- 4900000 + 30 * 30
  expect_false(arraydesign:::segment_wet(wg, x_wall - 1, y_mid, x_wall + 31, y_mid))
  expect_true(arraydesign:::segment_wet(wg, x_wall - 50, y_mid, x_wall - 1, y_mid + 20))
  expect_true(arraydesign:::segment_wet(wg, x_wall - 5, y_mid, x_wall - 5, y_mid))
  # The north-west corner of the peninsula tip in example_lake() is at
  # (725870, 4935810). Both segments join the same two water cells; the first
  # passes north of the corner, the second cuts 5.7 m through the land cell.
  wl <- arraydesign:::water_grid(example_lake())
  cx <- 725870
  cy <- 4935810
  expect_true(arraydesign:::segment_wet(wl, cx - 5, cy - 1, cx + 1, cy + 5))
  expect_false(arraydesign:::segment_wet(wl, cx - 1, cy - 5, cx + 5, cy + 1))
  # Cell arithmetic agrees with terra.
  pts <- withr::with_seed(1, cbind(stats::runif(500, 699000, 704000),
                                   stats::runif(500, 4899000, 4902500)))
  expect_equal(arraydesign:::cell_at(wg, pts[, 1], pts[, 2]),
               terra::cellFromXY(wg$raster, pts))
})

test_that("realised steps and turns away from the shore match the model (KS)", {
  lake <- open_lake()
  m <- crw_model(30, 10, 30)
  start <- data.frame(x = 700000 + 101 * 30, y = 4900000 + 101 * 30)
  tr <- simulate_tracks(lake, m, n_animals = 20, n_steps = 1000,
                        control = track_control(start = start), seed = 11)
  # Distance from each start of a step to the nearest land cell edge.
  shore <- pmin(tr$x - 700030, 706030 - tr$x, tr$y - 4900030, 4906030 - tr$y)
  n <- nrow(tr)
  same <- tr$path_id[-1] == tr$path_id[-n]
  len <- sqrt(diff(tr$x)^2 + diff(tr$y)^2)
  bearing <- atan2(diff(tr$y), diff(tr$x))
  away <- same & shore[-n] > 200
  steps <- len[away]
  expect_gt(length(steps), 15000)
  ks_step <- suppressWarnings(stats::ks.test(steps, withr::with_seed(2, m$step(1e5))))
  expect_lt(ks_step$statistic, 0.015)
  k <- which(away[-1] & away[-length(away)] & same[-1])
  turns <- wrap_angle(bearing[k + 1] - bearing[k])
  ks_turn <- suppressWarnings(stats::ks.test(turns, withr::with_seed(3, m$turn(1e5))))
  expect_lt(ks_turn$statistic, 0.015)
})

test_that("the move and rest switch gives the configured dwell times", {
  p <- matrix(c(0.95, 0.10, 0.05, 0.90), 2)
  tr <- simulate_tracks(open_lake(), crw_model(30, 10, 30, transition = p), n_animals = 1,
                        n_steps = 20000, seed = 5)
  s <- tr$state[-1]
  runs <- rle(s)
  move <- runs$lengths[runs$values == 1]
  rest <- runs$lengths[runs$values == 2]
  # Geometric dwell: mean 1 / (1 - p_ii); 4 standard errors either way.
  expect_lt(abs(mean(move) - 20), 4 * sqrt(0.95) / 0.05 / sqrt(length(move)))
  expect_lt(abs(mean(rest) - 10), 4 * sqrt(0.90) / 0.10 / sqrt(length(rest)))
  moved <- sqrt(diff(tr$x)^2 + diff(tr$y)^2) > 0
  expect_equal(sum(moved[s == 2]), 0)
  expect_equal(attr(tr, "escapes")$still, sum(s == 2))
})

test_that("a goal pulls paths in and holds them there", {
  lake <- example_lake()
  goal <- data.frame(x = 726500, y = 4935600)
  tr <- simulate_tracks(lake, crw_model(30, 10, 30), n_animals = 3, n_steps = 400,
                        control = track_control(start = data.frame(x = 726100, y = 4935400),
                                                goal = goal, goal_bias = 0.6),
                        seed = 4)
  arrived <- attr(tr, "arrived")
  expect_false(anyNA(arrived))
  for (i in 1:3) {
    p <- tr[tr$path_id == i, ]
    after <- p[p$step >= arrived[i], ]
    expect_lte(sqrt((after$x[1] - 726500)^2 + (after$y[1] - 4935600)^2), 30)
    expect_equal(length(unique(after$x)), 1L)
  }
  expect_warning(
    free <- simulate_tracks(lake, crw_model(30, 10, 30), n_animals = 2, n_steps = 50,
                            control = track_control(goal = c(x = 726500, y = 4935600),
                                                    goal_bias = 0, goal_tolerance = 0),
                            seed = 4),
    "No animal"
  )
  expect_true(all(is.na(attr(free, "arrived"))))
})

test_that("a step that meets a wall is reflected off it, as off a mirror", {
  wg <- arraydesign:::water_grid(wall_lake())
  x_wall <- 700000 + 49 * 30
  y0 <- 4900000 + 30 * 30
  # Head on: 10 m to the wall, 20 m back.
  mv <- arraydesign:::take_step(wg, x_wall - 10, y0, 0, 30, 50)
  expect_equal(c(mv$x, mv$y, mv$bearing), c(x_wall - 20, y0, pi), tolerance = 1e-9)
  expect_equal(mv$outcome, "reflected")
  # At 45 degrees: x mirrored at the wall, y carries on.
  mv <- arraydesign:::take_step(wg, x_wall - 10, y0, pi / 4, 30 * sqrt(2), 50)
  expect_equal(c(mv$x, mv$y, mv$bearing), c(x_wall - 20, y0 + 30, 3 * pi / 4), tolerance = 1e-9)
  # Into the north-west corner of open_lake(): off the west shore, then the north.
  wo <- arraydesign:::water_grid(open_lake())
  mv <- arraydesign:::take_step(wo, 700040, 4906015, 3 * pi / 4, 30 * sqrt(2), 50)
  expect_equal(c(mv$x, mv$y, mv$bearing), c(700050, 4906015, -pi / 4), tolerance = 1e-9)
  # With no reflections allowed the step holds.
  expect_equal(arraydesign:::take_step(wg, x_wall - 10, y0, 0, 30, 0)$outcome, "held")
  # A leg that starts on land has no edge to reflect off; the step holds.
  expect_equal(arraydesign:::first_land(wg, x_wall + 15, y0, x_wall + 45, y0),
               list(t = 0, axis = NA_integer_))
  expect_equal(arraydesign:::take_step(wg, x_wall + 15, y0, 0, 30, 50)$outcome, "held")
})

test_that("a reflected step whose straight line would cross land holds", {
  # The path runs from A north-east to a wall at x = 300 m and back
  # north-west to B. The straight line A-B cuts one land cell that the path
  # itself passes beside. Coordinates below are offsets from the lower-left.
  w <- matrix(TRUE, 20, 12)
  w[, 11] <- FALSE
  wg <- arraydesign:::water_grid(lake_raster(w))
  len <- 435 * sqrt(2)
  mv <- arraydesign:::take_step(wg, 700015, 4900015, pi / 4, len, 50)
  expect_equal(c(mv$x, mv$y), c(700150, 4900450), tolerance = 1e-6)
  w[12, 3] <- FALSE
  wp <- arraydesign:::water_grid(lake_raster(w))
  expect_false(arraydesign:::segment_wet(wp, 700015, 4900015, 700150, 4900450))
  mv <- arraydesign:::take_step(wp, 700015, 4900015, pi / 4, len, 50)
  expect_equal(c(mv$x, mv$y, mv$outcome), c(700015, 4900015, "held"))
})

test_that("time near the shore matches its share of the water, uniform and persistent turns", {
  # A walker that redraws or holds a blocked step can over- or under-fill the
  # ring of cells within one cell diagonal of the shore. The ring's share of
  # positions must match its share of the water, for uniform and for
  # persistent turns.
  lake <- example_lake()
  wg <- arraydesign:::water_grid(lake)
  land <- terra::rast(lake)
  terra::values(land) <- ifelse(wg$wet, NA, 1)
  ring <- wg$wet & terra::values(terra::distance(land), mat = FALSE) <= 30 * sqrt(2) + 1e-6
  share <- mean(ring[wg$wet])
  for (turn_sd in c(Inf, 30)) {
    tr <- simulate_tracks(lake, crw_model(30, 10, turn_sd), n_animals = 8, n_steps = 6000,
                          seed = 40)
    keep <- tr$step > 1000
    by_path <- tapply(ring[tr$cell[keep]], tr$path_id[keep], mean)
    z <- (mean(by_path) - share) / (stats::sd(by_path) / sqrt(8))
    expect_lt(abs(z), 3)
    expect_lt(abs(mean(by_path) / share - 1), 0.1)
  }
})

test_that("a lone water cell keeps the path inside it, and degenerate walks warn", {
  w <- matrix(FALSE, 5, 5)
  w[3, 3] <- TRUE
  pond <- lake_raster(w)
  expect_warning(tr <- simulate_tracks(pond, crw_model(60, 10, 30), n_steps = 25, seed = 1),
                 "were reflected")
  expect_true(all(tr$cell == 13L))
  expect_equal(attr(tr, "escapes")$reflected + attr(tr, "escapes")$held, 25)
  expect_warning(
    held <- simulate_tracks(pond, crw_model(60, 10, 30), n_steps = 25,
                            control = track_control(max_reflections = 0), seed = 1),
    "25 of 25 moving steps \\(100.0%\\) held"
  )
  expect_equal(length(unique(held$x)), 1L)
  # A step in kilometres on a lake 1.8 km across: the shore shapes every step.
  expect_warning(simulate_tracks(example_lake(), crw_model(5000, 100, 30), n_animals = 2,
                                 n_steps = 50, seed = 1), "moving steps")
  expect_silent(simulate_tracks(example_lake(), crw_model(30, 10, 30), n_animals = 2,
                                n_steps = 500, seed = 1))
})

test_that("a goal must connect to every start by water, and a missed goal warns", {
  wall <- wall_lake()
  east <- data.frame(x = 700000 + 75 * 30, y = 4900000 + 30 * 30)
  west <- data.frame(x = 700000 + 25 * 30, y = 4900000 + 30 * 30)
  expect_error(simulate_tracks(wall, crw_model(30, 10, 30), n_steps = 5,
                               control = track_control(start = west, goal = east)),
               "does not connect to the goal")
  # Starts drawn at random come from the goal's basin only.
  tr <- suppressWarnings(simulate_tracks(wall, crw_model(30, 10, 30), n_animals = 30, n_steps = 1,
                                         control = track_control(goal = east, goal_bias = 0),
                                         seed = 1))
  expect_true(all(tr$x > 700000 + 50 * 30))
  near <- data.frame(x = west$x + 30, y = west$y)
  expect_warning(simulate_tracks(wall, crw_model(30, 10, 30), n_animals = 2, n_steps = 5,
                                 control = track_control(start = near, goal = west,
                                                         goal_bias = 0, goal_tolerance = 0),
                                 seed = 1),
                 "No animal of 2 reached the goal in 5 steps")
})

test_that("results record the generator kind, which a seed alone does not fix", {
  lake <- example_lake()
  m <- crw_model(30, 10, 30)
  old <- RNGkind()
  withr::defer(do.call(RNGkind, as.list(old)))
  RNGkind("L'Ecuyer-CMRG")
  a <- simulate_tracks(lake, m, n_animals = 2, n_steps = 30, seed = 3)
  tx_a <- simulate_transmissions(a, seed = 3)
  expect_equal(attr(a, "rng_kind")[1], "L'Ecuyer-CMRG")
  expect_equal(attr(tx_a, "rng_kind")[1], "L'Ecuyer-CMRG")
  RNGkind("Mersenne-Twister")
  b <- simulate_tracks(lake, m, n_animals = 2, n_steps = 30, seed = 3)
  expect_equal(attr(b, "rng_kind")[1], "Mersenne-Twister")
  expect_false(identical(a$x, b$x))
  do.call(RNGkind, as.list(attr(a, "rng_kind")))
  again <- simulate_tracks(lake, m, n_animals = 2, n_steps = 30, seed = 3)
  expect_identical(again$x, a$x)
  expect_identical(simulate_transmissions(again, seed = 3)$x, tx_a$x)
})

test_that("the table has positions, cells and times in the stated form", {
  lake <- example_lake()
  tr <- simulate_tracks(lake, crw_model(30, 10, 30, time_step_s = 120), n_animals = 2,
                        n_steps = 50, seed = 1)
  expect_s3_class(tr, "sim_tracks")
  expect_named(tr, c("path_id", "step", "time", "x", "y", "cell", "state"))
  expect_equal(nrow(tr), 102)
  expect_equal(attr(tr$time, "tzone"), "UTC")
  expect_equal(unique(diff(as.numeric(tr$time[tr$path_id == 1]))), 120)
  expect_equal(tr$cell, terra::cellFromXY(lake, cbind(tr$x, tr$y)))
  expect_equal(attr(tr, "epsg"), 32617L)
  expect_s3_class(attr(tr, "model"), "movement_model")
  expect_equal(attr(tr, "seed"), 1)
})

test_that("starts can be given, one for all or one per path, and must be on water", {
  lake <- example_lake()
  m <- crw_model(30, 10, 30)
  one <- simulate_tracks(lake, m, n_animals = 3, n_steps = 5,
                         control = track_control(start = c(x = 725300, y = 4935600)), seed = 1)
  expect_equal(one$x[one$step == 0], rep(725300, 3))
  each <- data.frame(x = c(725300, 725400), y = c(4935600, 4935600))
  two <- simulate_tracks(lake, m, n_animals = 2, n_steps = 5,
                         control = track_control(start = each), seed = 1)
  expect_equal(two$x[two$step == 0], each$x)
  expect_error(simulate_tracks(lake, m, n_animals = 3, n_steps = 5,
                               control = track_control(start = each)), "1 row or n_animals")
  expect_error(simulate_tracks(lake, m, n_steps = 5,
                               control = track_control(start = c(x = 725900, y = 4935600))),
               "on land")
  expect_error(simulate_tracks(lake, m, control = track_control(start = data.frame(a = 1))),
               "columns `x` and `y`")
  expect_error(simulate_tracks(lake, m, control = track_control(start = c(x = NA, y = 1))),
               "missing coordinates")
  east <- data.frame(x = "east", y = 4935600)
  expect_error(simulate_tracks(lake, m, control = track_control(start = east)),
               "`start\\$x` holds values that are not numbers: 1 of 1, the first \"east\" in row 1",
               class = "arraydesign_error_input")
  expect_error(simulate_tracks(lake, m, control = track_control(goal = each)), "one point")
  expect_error(simulate_tracks(lake, m, control = track_control(goal = c(x = 1, y = 1))),
               "goal` is on land")
})

test_that("bad arguments stop with a message", {
  lake <- example_lake()
  m <- crw_model(30, 10, 30)
  expect_error(simulate_tracks(lake, list()), "movement_model")
  expect_error(simulate_tracks(lake, m, n_animals = 0), "n_animals")
  expect_error(simulate_tracks(lake, m, n_steps = 1.5), "n_steps")
  expect_error(simulate_tracks(lake, m, control = list()), "track_control")
  bad <- m
  bad$step <- function(n) rep(1, max(n, 5) - 1)
  expect_error(simulate_tracks(lake, bad, n_steps = 10), "not valid")
  expect_error(track_control(goal_bias = 2), "goal_bias")
  expect_error(track_control(goal_tolerance = -1), "goal_tolerance")
  expect_error(track_control(start_time = Sys.Date()), "POSIXct")
  expect_error(track_control(start_time = as.POSIXct("2025-01-01", tz = "America/Toronto")),
               "UTC")
  expect_error(track_control(start_time = as.POSIXct(c("2025-01-01", "2025-01-02"), tz = "UTC")),
               "one time")
  expect_error(track_control(max_reflections = -1), "max_reflections")
  expect_error(track_control(max_reflections = 1.5), "max_reflections")
})

test_that("a raster with no CRS is refused, naming the call", {
  nocrs <- example_lake()
  terra::crs(nocrs) <- ""
  m <- crw_model(30, 10, 30)
  err <- tryCatch(simulate_tracks(nocrs, m, n_steps = 2), error = identity)
  expect_s3_class(err, "arraydesign_error_grid")
  expect_identical(conditionMessage(err), paste(
    "`grid` has no CRS, so latitude and longitude cannot be computed.",
    "Set one with terra::crs()."
  ))
  expect_identical(conditionCall(err), quote(simulate_tracks(nocrs, m, n_steps = 2)))
})

test_that("a seed reproduces tracks exactly and leaves the caller's stream alone", {
  lake <- example_lake()
  m <- crw_model(30, 10, 30, transition = matrix(c(0.9, 0.2, 0.1, 0.8), 2))
  set.seed(99)
  before <- .Random.seed
  a <- simulate_tracks(lake, m, n_animals = 4, n_steps = 300, seed = 123)
  expect_identical(.Random.seed, before)
  b <- simulate_tracks(lake, m, n_animals = 4, n_steps = 300, seed = 123)
  expect_identical(as.data.frame(a), as.data.frame(b))
  expect_identical(attr(a, "escapes"), attr(b, "escapes"))
  c2 <- simulate_tracks(lake, m, n_animals = 4, n_steps = 300, seed = 124)
  expect_false(identical(a$x, c2$x))
  set.seed(5)
  d1 <- simulate_tracks(lake, m, n_steps = 20)
  set.seed(5)
  d2 <- simulate_tracks(lake, m, n_steps = 20)
  expect_identical(d1$x, d2$x)
})
