# ---- available_steps(): controls for a step selection check ----

selection_tracks <- function(beta = -0.2, n_animals = 3, n_steps = 300, transition = NULL,
                             seed = 1) {
  lake <- example_lake()
  m <- habitat_selection(crw_model(30, 15, 40, transition = transition), lake, beta = beta)
  simulate_tracks(lake, m, n_animals = n_animals, n_steps = n_steps, seed = seed)
}

test_that("each moving step has its used end and n controls", {
  tracks <- selection_tracks()
  steps <- available_steps(tracks, n = 7, seed = 2)
  moving <- sum(!is.na(tracks$heading))
  expect_identical(nrow(steps), moving * 8L)
  expect_identical(names(steps), c("path_id", "step", "stratum", "case", "x", "y", "cell",
                                   "depth_m", "step_length"))
  expect_true(all(tapply(steps$case, steps$stratum, sum) == 1L))
  expect_true(all(table(steps$stratum) == 8L))
  expect_identical(attr(steps, "n"), 7L)
  expect_identical(attr(steps, "epsg"), 32617L)
  # The used rows are the track's own positions, in order
  used <- steps[steps$case == 1L, ]
  rows <- tracks[!is.na(tracks$heading), ]
  expect_equal(used$x, rows$x)
  expect_equal(used$y, rows$y)
  expect_equal(used$step, rows$step)
  # Step length of a used row is the distance from the previous position
  prev <- which(!is.na(tracks$heading)) - 1L
  expect_equal(used$step_length,
               sqrt((rows$x - tracks$x[prev])^2 + (rows$y - tracks$y[prev])^2))
})

test_that("every end is on water and the layer is read at its cell", {
  lake <- example_lake()
  steps <- available_steps(selection_tracks(), n = 10, seed = 2)
  values <- terra::values(lake, mat = FALSE)
  expect_false(anyNA(steps$cell))
  expect_false(anyNA(values[steps$cell]))
  expect_equal(steps$depth_m, values[steps$cell])
  expect_equal(steps$cell, terra::cellFromXY(lake, cbind(steps$x, steps$y)))
})

test_that("with beta = 0 the used steps and the controls come from one kernel", {
  tracks <- selection_tracks(beta = 0, n_animals = 10, n_steps = 1000)
  expect_equal(attr(tracks, "selection")$accepted, attr(tracks, "selection")$proposals)
  steps <- available_steps(tracks, n = 5, seed = 3)
  ks <- suppressWarnings(stats::ks.test(steps$step_length[steps$case == 1L],
                                        steps$step_length[steps$case == 0L]))
  expect_gt(ks$p.value, 0.001)
})

test_that("a seed gives the same controls and restores the caller's stream", {
  tracks <- selection_tracks()
  set.seed(9)
  before <- stats::runif(1)
  set.seed(9)
  a <- available_steps(tracks, n = 3, seed = 4)
  after <- stats::runif(1)
  b <- available_steps(tracks, n = 3, seed = 4)
  expect_identical(a, b)
  expect_identical(before, after)
})

test_that("resting steps have no controls", {
  p <- matrix(c(0.9, 0.2, 0.1, 0.8), 2)
  tracks <- selection_tracks(transition = p, n_steps = 500)
  steps <- available_steps(tracks, n = 4, seed = 2)
  resting <- tracks$step[tracks$state == 2L & tracks$step > 0]
  expect_gt(length(resting), 0)
  key <- paste(tracks$path_id, tracks$step)[tracks$state == 2L & tracks$step > 0]
  expect_false(any(paste(steps$path_id, steps$step) %in% key))
})

test_that("a layer of the user's gives one column per layer", {
  lake <- example_lake()
  extra <- c(lake, lake * 2)
  names(extra) <- c("depth_m", "double_depth")
  steps <- available_steps(selection_tracks(), layer = extra, n = 2, seed = 1)
  expect_equal(steps$double_depth, 2 * steps$depth_m)
})

test_that("available_steps() refuses what it cannot use", {
  lake <- example_lake()
  plain <- simulate_tracks(lake, crw_model(30, 15, 40), n_animals = 1, n_steps = 20, seed = 1)
  expect_error(available_steps(plain), "hold no `heading`")
  tracks <- selection_tracks()
  expect_error(available_steps(as.data.frame(tracks)), "sim_tracks object")
  expect_error(available_steps(tracks, n = 0), "`n` must be one whole number")
  expect_error(available_steps(tracks, layer = terra::extend(lake, 1)), "geometry")
  named <- lake
  names(named) <- "step"
  expect_error(available_steps(tracks, layer = named), "not one of the output columns")
})
