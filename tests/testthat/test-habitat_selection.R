# ---- habitat_selection(): a known preference in the walker ----
# The walker accepts each candidate step with probability w(end) / W, so the
# selection coefficient `beta` is the truth a conditional logit on the true
# steps reads. These tests check the constructor, that a model without a
# preference takes the walk without rejection, the sampler's bookkeeping, and
# that the true steps recover `beta` on Stoney Lake.

stoney_depth <- function() {
  lake <- terra::rast(system.file("extdata", "stoney_depth.tif", package = "arraydesign"))
  depth <- abs(lake)
  names(depth) <- "depth_m"
  list(lake = lake, depth = depth)
}

test_that("the constructor takes exactly one of beta and preference", {
  depth <- example_lake()
  walk <- crw_model(30, 15, 40)
  expect_error(habitat_selection(walk, depth), "exactly one of `beta` and `preference`")
  expect_error(habitat_selection(walk, depth, beta = 0.1, preference = function(z) z),
               "exactly one of `beta` and `preference`")
  expect_error(habitat_selection(walk, depth, beta = c(0.1, 0.2)), "one finite number")
  expect_error(habitat_selection(walk, depth, beta = Inf), "one finite number")
  expect_error(habitat_selection(list(), depth, beta = 0.1), "must be a movement_model")
  expect_error(habitat_selection(walk, terra::values(depth), beta = 0.1), "SpatRaster")
  expect_error(habitat_selection(walk, c(depth, depth), beta = 0.1), "one layer")
  once <- habitat_selection(walk, depth, beta = 0.1)
  expect_error(habitat_selection(once, depth, beta = 0.1), "already has a selection")
})

test_that("the constructor refuses a layer or a preference that is not finite", {
  depth <- example_lake()
  walk <- crw_model(30, 15, 40)
  infinite <- depth
  infinite[100] <- Inf
  expect_error(habitat_selection(walk, infinite, beta = 0.1), "1 infinite values")
  missing_below_5 <- function(z) ifelse(z > 5, 0, NA)
  expect_error(habitat_selection(walk, depth, preference = missing_below_5), "not finite")
  expect_error(habitat_selection(walk, depth, preference = function(z) 1), "one number per")
  expect_error(habitat_selection(walk, depth, preference = "z"), "must be a function")
})

test_that("the model stores beta as given, its source and a print line", {
  depth <- example_lake()
  m <- habitat_selection(crw_model(30, 15, 40), depth, beta = -0.3)
  expect_identical(m$selection$beta, -0.3)
  expect_identical(m$selection$layer_name, "depth_m")
  expect_match(m$source, "selection on depth_m: beta -0.3 per unit$")
  expect_output(print(m), "selection: log weight -0.3 per unit of depth_m")
  q <- habitat_selection(crw_model(30, 15, 40), depth, preference = function(z) -z^2)
  expect_null(q$selection$beta)
  expect_output(print(q), "log weight from a function of depth_m")
})

test_that("simulate_tracks() refuses a layer off the grid or missing on water", {
  lake <- example_lake()
  walk <- crw_model(30, 15, 40)
  other <- terra::extend(lake, 2)
  m <- habitat_selection(walk, other, beta = 0.1)
  expect_error(simulate_tracks(lake, m, n_animals = 1, n_steps = 5, seed = 1),
               "does not match the geometry of `grid`")
  holed <- lake
  holed[which(!is.na(terra::values(lake, mat = FALSE)))[1]] <- NA
  m <- habitat_selection(walk, holed, beta = 0.1)
  expect_error(simulate_tracks(lake, m, n_animals = 1, n_steps = 5, seed = 1),
               "no finite log weight on 1 of 2152 water cells")
})

test_that("a model without a preference gives the tracks of detfield's stoney_track", {
  skip_if_not(exists("stoney_track", envir = asNamespace("detfield")) ||
                "stoney_track" %in% utils::data(package = "detfield")$results[, "Item"],
              "detfield has no stoney_track")
  ref <- detfield::stoney_track$track
  lake <- stoney_depth()$lake
  tracks <- simulate_tracks(
    lake, crw_model(step_mean = 30, step_sd = 15, turn_sd = 40, time_step_s = 60),
    n_animals = 5, n_steps = 48 * 60,
    control = track_control(start_time = as.POSIXct("2024-06-01", tz = "UTC")), seed = 1L
  )
  expect_identical(nrow(tracks), nrow(ref))
  # Equal to floating-point tolerance: on the macOS CI runner the coordinates
  # differed from the stored track from step 777 on, by 1e-10 to 1e-9 m.
  expect_equal(tracks$x, ref$x)
  expect_equal(tracks$y, ref$y)
  expect_identical(tracks$time, ref$datetime)
  expect_false("heading" %in% names(tracks))
  expect_null(attr(tracks, "selection"))
})

test_that("the sampler's counts add up, and resting steps draw no candidate", {
  st <- stoney_depth()
  p <- matrix(c(0.95, 0.10, 0.05, 0.90), 2)
  m <- habitat_selection(crw_model(30, 15, 40, transition = p), st$depth, beta = -0.3)
  tracks <- simulate_tracks(st$lake, m, n_animals = 5, n_steps = 2000, seed = 3)
  sel <- attr(tracks, "selection")
  esc <- attr(tracks, "escapes")
  expect_identical(names(sel), c("path_id", "proposals", "accepted", "over_reach",
                                 "above_bound", "exhausted"))
  expect_true(all(sel$proposals >= sel$accepted))
  expect_equal(rowSums(esc[, -1]), rep(2000, 5))
  # Every moving step is accepted or exhausted; a resting step draws nothing
  moving <- tapply(tracks$state[tracks$step > 0] == 1L, tracks$path_id[tracks$step > 0], sum)
  expect_equal(sel$accepted + sel$exhausted, as.vector(moving))
  expect_equal(esc$still, as.vector(2000 - moving))
  expect_true(all(is.na(tracks$heading[tracks$state == 2L])))
  expect_equal(sum(!is.na(tracks$heading)), sum(sel$accepted))
  # A candidate over the bound must be longer than the reach
  expect_lte(sum(sel$above_bound), sum(sel$over_reach))
  expect_lte(sum(sel$over_reach), sum(sel$accepted) / 1000)
  expect_equal(attr(tracks, "model")$selection$window, 21L)
  expect_gt(attr(tracks, "model")$selection$reach, 110)
  expect_lt(attr(tracks, "model")$selection$reach, 130)
  expect_output(print(tracks), "candidates accepted")
})

test_that("the bound covers every cell a step within the reach can end in", {
  wg <- arraydesign:::water_grid(example_lake())
  m <- habitat_selection(crw_model(30, 15, 40), example_lake(), beta = 0.2)
  sel <- arraydesign:::selection_setup(wg, m)
  # Ends of 20,000 steps of length up to the reach from random water points;
  # no end may weigh more than the bound read at its start
  pts <- withr::with_seed(7, {
    cells <- sample(which(wg$wet), 20000, replace = TRUE)
    xy <- arraydesign:::cell_xy(wg, cells) +
      matrix(stats::runif(40000, -14.9, 14.9), ncol = 2)
    len <- stats::runif(20000, 0, sel$reach)
    b <- stats::runif(20000, -pi, pi)
    cbind(xy, xy[, 1] + len * cos(b), xy[, 2] + len * sin(b))
  })
  end_cell <- arraydesign:::cell_at(wg, pts[, 3], pts[, 4])
  ok <- !is.na(end_cell) & wg$wet[ifelse(is.na(end_cell), 1L, end_cell)]
  bound <- vapply(which(ok), function(i) {
    arraydesign:::point_bound(wg, sel, pts[i, 1], pts[i, 2])
  }, numeric(1))
  expect_gt(sum(ok), 15000)
  expect_true(all(sel$log_w[end_cell[ok]] <= bound + 1e-12))
})

test_that("a goal and a preference together still bring the animals in", {
  lake <- example_lake()
  m <- habitat_selection(crw_model(30, 10, 30), lake, beta = 0.1)
  tracks <- simulate_tracks(lake, m, n_animals = 3, n_steps = 400,
                            control = track_control(start = data.frame(x = 726100, y = 4935400),
                                                    goal = data.frame(x = 726500, y = 4935600),
                                                    goal_bias = 0.6),
                            seed = 4)
  expect_false(anyNA(attr(tracks, "arrived")))
  expect_equal(sum(attr(tracks, "selection")$exhausted), 0)
})

test_that("a quadratic preference puts the median used depth near its optimum", {
  # example_lake() deepens by about 1 m per 30 m cell, so a gentle curvature
  # keeps the acceptance near 0.17 and no step runs out of candidates
  lake <- example_lake()
  start <- data.frame(x = 725500, y = 4935600)
  walk <- crw_model(30, 15, 40)
  m <- habitat_selection(walk, lake, preference = function(z) -0.05 * (z - 15)^2)
  chosen <- simulate_tracks(lake, m, n_animals = 4, n_steps = 1000,
                            control = track_control(start = start), seed = 5)
  plain <- simulate_tracks(lake, walk, n_animals = 4, n_steps = 1000,
                           control = track_control(start = start), seed = 5)
  depth <- terra::values(lake, mat = FALSE)
  used <- depth[chosen$cell[chosen$step > 200]]
  free <- depth[plain$cell[plain$step > 200]]
  expect_lt(abs(stats::median(used) - 15), 1.5)
  expect_gt(abs(stats::median(free) - 15), 3)
})

# ---- Step-scale recovery on Stoney Lake ----
# 20 animals for 48 h from the same 20 starts at each beta.
# The replicate is the animal; the seed only makes the run reproducible.
recovery_runs <- local({
  runs <- NULL
  function() {
    if (!is.null(runs)) return(runs)
    st <- stoney_depth()
    wg <- arraydesign:::water_grid(st$lake)
    starts <- withr::with_seed(20, arraydesign:::cell_xy(wg, sample(which(wg$wet), 20)))
    walk <- crw_model(30, 15, 40)
    runs <<- lapply(c(neg = -0.3, zero = 0, pos = 0.3), function(beta) {
      tracks <- simulate_tracks(st$lake, habitat_selection(walk, st$depth, beta = beta),
                                n_animals = 20, n_steps = 2880,
                                control = track_control(start = starts), seed = 1)
      steps <- available_steps(tracks, n = 10, seed = 2)
      list(tracks = tracks, fits = clogit_per_animal(steps))
    })
    runs
  }
})

clogit_per_animal <- function(steps) {
  do.call(rbind, lapply(split(steps, steps$path_id), function(one) {
    fit <- survival::clogit(case ~ depth_m + strata(stratum), data = one)
    data.frame(path_id = one$path_id[1], estimate = unname(stats::coef(fit)),
               se = sqrt(unname(diag(stats::vcov(fit)))))
  }))
}

test_that("a conditional logit on the true steps recovers beta at -0.3, 0 and 0.3", {
  skip_on_cran()
  skip_if_not_installed("survival")
  # clogit() calls coxph() by name, so survival must be attached
  suppressPackageStartupMessages(library(survival))
  runs <- recovery_runs()
  for (beta in names(runs)) {
    truth <- c(neg = -0.3, zero = 0, pos = 0.3)[[beta]]
    fits <- runs[[beta]]$fits
    expect_identical(nrow(fits), 20L)
    between_se <- stats::sd(fits$estimate) / sqrt(20)
    expect_lt(abs(mean(fits$estimate) - truth), 2 * between_se)
    expect_lt(abs(mean(fits$estimate) - truth), 0.03)
    covered <- sum(abs(fits$estimate - truth) <= stats::qnorm(0.975) * fits$se)
    expect_gte(covered, 17L)
  }
})

test_that("median bottom depth used rises from beta = -0.3 to 0.3 in paired animals", {
  skip_on_cran()
  skip_if_not_installed("survival")
  suppressPackageStartupMessages(library(survival))
  runs <- recovery_runs()
  depth <- terra::values(stoney_depth()$depth, mat = FALSE)
  median_depth <- function(tracks) tapply(depth[tracks$cell], tracks$path_id, stats::median)
  rises <- median_depth(runs$pos$tracks) > median_depth(runs$neg$tracks)
  expect_gte(sum(rises), 16L)
})
