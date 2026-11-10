# ---- Shared settings ----
# Animals that stay put at a cell centre and a tag that sends exactly every
# 120 s: the conditions under which location_surface()'s likelihood is the
# model that made the data, so its regions are calibrated by construction.
still <- crw_model(1e-6, 1e-6, 40)
exact <- list(movement = still, tag_delay = 120)
toy <- example_lake()
toy_rx <- place_receivers(toy, n = 8)
logistic <- detfield::de_logistic(250, 600)

# ---- evaluate_array(): the object ----

test_that("one call places, simulates, positions and scores", {
  ev <- suppressMessages(evaluate_array(toy_rx, toy, logistic, n_animals = 3,
                                        seed = 1, control = list(hours = 4)))
  expect_s3_class(ev, "array_evaluation")
  expect_s3_class(ev$score, "surface_score")
  expect_s3_class(ev$tracks, "sim_tracks")
  expect_s3_class(ev$detections, "sim_detections")
  expect_s3_class(ev$surfaces, "loc_surface")
  expect_equal(ev$score$overall$n_bins, nrow(attr(ev$surfaces, "bins")))
  expect_equal(length(unique(ev$tracks$path_id)), 3)
  expect_equal(max(ev$tracks$time) - min(ev$tracks$time), as.difftime(4, units = "hours"))
  expect_identical(names(ev$timing), c("field", "simulate", "position", "score"))
  expect_identical(ev$settings$n_animals, 3)
  expect_identical(ev$settings$position_model, logistic$source)
  out <- capture.output(print(ev))
  expect_match(out[1], "8 receivers; 3 animals x 4 hours")
  expect_false(any(grepl("positioned with", out)))
  expect_false(any(grepl("too narrow", out)))
})

test_that("seeds reproduce and the caller's stream is restored", {
  run <- function(seed) {
    suppressMessages(evaluate_array(toy_rx, toy, logistic, n_animals = 2,
                                    seed = seed, control = list(hours = 3)))
  }
  set.seed(99)
  before <- stats::runif(1)
  set.seed(99)
  a <- run(7)
  after <- stats::runif(1)
  b <- run(7)
  c2 <- run(8)
  expect_identical(before, after)
  expect_identical(a$detections, b$detections)
  expect_identical(a$score$bins, b$score$bins)
  expect_false(identical(a$detections$detections, c2$detections$detections))
})

test_that("seeds reproduce with the station-hour effect too", {
  m <- detfield::de_model(logistic, sigma = list(station_hour = 1))
  run <- function() {
    suppressMessages(evaluate_array(toy_rx, toy, m, n_animals = 2,
                                    seed = 3, control = list(hours = 2)))
  }
  a <- run()
  b <- run()
  expect_identical(a$detections, b$detections)
  expect_identical(a$score$overall, b$score$overall)
  expect_output(print(a), "station-hour effect and 2 hours the hit-rate interval is too narrow",
                fixed = TRUE)
  a$settings$hours <- 24
  expect_false(any(grepl("too narrow", capture.output(print(a)))))
  expect_identical(attr(a$detections, "station_hour"), 1)
  expect_identical(sum(a$detections$transmissions$n_detections),
                   nrow(a$detections$detections))
  expect_identical(a$detections$transmissions$transmission_id,
                   seq_len(nrow(a$detections$transmissions)))
})

test_that("a data frame of coordinates is placed on the grid", {
  coords <- data.frame(x = toy_rx$x[1:4], y = toy_rx$y[1:4])
  ev <- suppressMessages(evaluate_array(coords, toy, logistic, n_animals = 2,
                                        seed = 1, control = list(hours = 2)))
  expect_s3_class(ev$layout, "receiver_layout")
  expect_equal(ev$layout$x, coords$x)
  expect_equal(ncol(ev$field$eta), 4)
  # A column subset of a layout keeps its class; it is placed like coordinates.
  sub <- toy_rx[1:4, c("x", "y")]
  expect_s3_class(sub, "receiver_layout")
  ev_sub <- suppressMessages(evaluate_array(sub, toy, logistic, n_animals = 2,
                                            seed = 1, control = list(hours = 2)))
  expect_identical(ev_sub$layout, ev$layout)
  expect_identical(ev_sub$score$bins, ev$score$bins)
})

# ---- evaluate_array(): coverage of the water ----

test_that("k = 1 is coverage_summary() and k = 2 counts cells two receivers cover", {
  ev1 <- suppressMessages(evaluate_array(toy_rx, toy, logistic, n_animals = 1,
                                         seed = 1,
                                         control = list(hours = 1, thresholds = c(0.3, 0.5))))
  ev2 <- suppressMessages(evaluate_array(toy_rx, toy, logistic, n_animals = 1,
                                         seed = 1,
                                         control = list(k = 2, hours = 1,
                                                        thresholds = c(0.3, 0.5))))
  base <- coverage_summary(ev1$field, thresholds = c(0.3, 0.5))
  expect_equal(ev1$coverage[names(base)], base)
  expect_identical(ev1$coverage$k, c(1L, 1L))
  p <- stats::plogis(ev2$field$eta)
  expect_equal(ev2$coverage$n_cells, c(sum(rowSums(p >= 0.3) >= 2), sum(rowSums(p >= 0.5) >= 2)))
  expect_equal(ev2$coverage$share, ev2$coverage$n_cells / 2152)
  expect_true(all(ev2$coverage$n_cells <= ev1$coverage$n_cells))
  expect_output(print(ev2), "cells at p >= 0.3 under 2\\+ receivers")
})

test_that("a k above the number of receivers warns", {
  expect_warning(
    ev <- suppressMessages(evaluate_array(toy_rx, toy, logistic, n_animals = 1,
                                          seed = 1, control = list(k = 9, hours = 1))),
    "`control\\$k` is 9 but the layout has 8 receivers"
  )
  expect_identical(ev$coverage$n_cells, c(0L, 0L))
  expect_warning(suppressMessages(evaluate_array(toy_rx[1, ], toy, logistic, n_animals = 1,
                                                 seed = 1, control = list(k = 2, hours = 1))),
                 "has 1 receiver, so")
  expect_silent(suppressMessages(evaluate_array(toy_rx, toy, logistic, n_animals = 1,
                                                seed = 1, control = list(k = 8, hours = 1))))
})

# ---- evaluate_array(): what was scored ----

test_that("the score states how many animal-hours were scored", {
  one <- place_receivers(toy, n = 1)
  ev <- suppressMessages(evaluate_array(one, toy, logistic, n_animals = 10,
                                        seed = 5, control = list(hours = 3)))
  o <- ev$score$overall
  tx <- ev$detections$transmissions
  hours <- unique(data.frame(tx$animal_id, floor(as.numeric(tx$datetime) / 3600)))
  expect_identical(o$n_animal_hours, nrow(hours))
  expect_identical(o$n_animal_hours, 30L)
  expect_lt(o$n_bins, o$n_animal_hours)
  expect_output(print(ev), "<array_evaluation> 1 receiver; 10 animals", fixed = TRUE)
  expect_equal(o$scored_share, o$n_bins / 30)
  expect_output(print(ev), sprintf("scored %d of 30 animal-hours with transmissions \\(%.1f%%\\)",
                                   o$n_bins, 100 * o$n_bins / 30))
  full <- suppressMessages(evaluate_array(toy_rx, toy, logistic, n_animals = 3,
                                          seed = 1, control = list(hours = 2)))
  expect_identical(full$score$overall$scored_share, 1)
})


# ---- evaluate_array(): calibration under the likelihood's own model ----

test_that("95% regions hold the truth at nominal under the binomial model", {
  ev <- suppressMessages(evaluate_array(toy_rx, toy, logistic, n_animals = 100,
                                        seed = 3, control = c(exact, list(hours = 24))))
  o <- ev$score$overall
  expect_equal(o$n_bins, 2400)
  expect_lt(abs(o$hit_rate - 0.95), 3 * o$hit_rate_se)
})

test_that("95% regions hold the truth at nominal with a station-hour effect", {
  skip_on_cran()
  m <- detfield::de_model(logistic, sigma = list(station_hour = 1))
  ev <- suppressMessages(evaluate_array(toy_rx, toy, m, n_animals = 40,
                                        seed = 4, control = c(exact, list(hours = 12))))
  o <- ev$score$overall
  expect_gt(o$n_bins, 400)
  expect_lt(abs(o$hit_rate - 0.95), 3 * o$hit_rate_se)
})

# ---- evaluate_array(): misspecified range ----

test_that("error rises as the positioning d50 moves from the true d50", {
  scales <- c(0.5, 0.71, 1, 1.41, 2)
  err <- vapply(scales, function(s) {
    ev <- suppressMessages(evaluate_array(
      toy_rx, toy, logistic, n_animals = 40,
      seed = 5,
      control = list(movement = still, hours = 8,
                     position_model = perturb_de_model(logistic, d50_scale = s))
    ))
    unlist(ev$score$overall[c("median_error_mean_m", "median_error_mode_m", "n_bins")])
  }, numeric(3))
  expect_true(all(err[3, ] == err[3, 1]))
  mean_err <- err[1, ]
  expect_true(all(diff(mean_err[1:3]) < 0))
  expect_true(all(diff(mean_err[3:5]) > 0))
  mode_err <- err[2, ]
  expect_true(all(diff(mode_err[1:3]) <= 0) && all(diff(mode_err[3:5]) >= 0))
  expect_gt(mode_err[1], mode_err[3])
  expect_gt(mode_err[5], mode_err[3])
})

test_that("print names the positioning model when it differs", {
  ev <- suppressMessages(evaluate_array(
    toy_rx, toy, logistic, n_animals = 2,
    seed = 1, control = list(hours = 2, position_model = perturb_de_model(logistic, 2, 2))
  ))
  expect_output(print(ev), "positioned with: perturb_de_model\\(d50 x 2, d05 x 2")
  ev$settings$skipped <- 2L
  expect_output(print(ev), "2 degenerate bins left out")
})

# ---- evaluate_array(): speed ----

test_that("the defaults score the toy lake in a median under 10 s at Stoney's sigma", {
  skip_on_cran()
  # A median of three runs: a single run can come close to 10 s on a loaded machine.
  m <- detfield::de_model(logistic, sigma = list(station_hour = 2.55))
  elapsed <- vapply(1:3, function(i) {
    t0 <- proc.time()[["elapsed"]]
    ev <- suppressMessages(evaluate_array(toy_rx, toy, m, seed = i, control = list()))
    expect_equal(ev$score$overall$n_animals, 20)
    proc.time()[["elapsed"]] - t0
  }, numeric(1))
  expect_lt(stats::median(elapsed), 10)
})

# ---- evaluate_array(): nothing detected, and plots ----

test_that("an array that hears nothing scores no bins", {
  deaf <- detfield::de_model(function(distance_m) -40 + 0 * distance_m)
  ev <- evaluate_array(toy_rx, toy, deaf, n_animals = 2, seed = 1, control = list(hours = 2))
  expect_null(ev$surfaces)
  expect_identical(ev$score$overall$n_bins, 0L)
  expect_identical(ev$settings$skipped, 0L)
  expect_identical(ev$coverage$n_cells, c(0L, 0L))
  out <- capture.output(print(ev))
  expect_match(out, "0 animal-bins of 0 animals", all = FALSE)
  p <- plot(ev)
  expect_s3_class(p, "ggplot")
  expect_equal(nrow(ggplot2::layer_data(p, 3)), 0)
  # With no detection every receiver is a white ring.
  expect_equal(nrow(ggplot2::layer_data(p, 4)), nrow(toy_rx))
  expect_identical(toupper(unique(ggplot2::layer_data(p, 5)$colour)), "#FFFFFF")
  expect_equal(nrow(ggplot2::layer_data(p, 7)), 0)
})

test_that("plot draws a map of the errors", {
  ev <- suppressMessages(evaluate_array(toy_rx, toy, logistic, n_animals = 2,
                                        seed = 1, control = list(hours = 2)))
  p <- plot(ev)
  expect_true(withVisible(plot(ev))$visible)
  pts <- ggplot2::layer_data(p, 3)
  expect_equal(nrow(pts), nrow(ev$score$bins))
  expect_equal(pts$x, ev$score$bins$truth_x)
  cls <- findInterval(ev$score$bins$error_mean_m, c(0, 200, 400, 800))
  expect_equal(as.integer(p$data$error_class[p$data$what == "animal-hour"]), cls)
  mode <- plot(ev, error = "mode", breaks = c(0, 100))
  expect_identical(levels(mode$data$error_class), c("0-100", "> 100"))
  # The error classes are on the colour scale; fill is the depth.
  expect_identical(ggplot2::get_labs(mode)$colour, "Error of\nmode (m)")
  png_file <- tempfile(fileext = ".png")
  ggplot2::ggsave(png_file, p, width = 6, height = 4.5, dpi = 100)
  expect_gt(file.size(png_file), 1000)
  expect_error(plot(ev, error = "median"), "should be one of",
               class = "arraydesign_error_input")
})

test_that("plot draws every animal-hour in the class that holds its error", {
  ev <- suppressMessages(evaluate_array(toy_rx, toy, logistic, n_animals = 2,
                                        seed = 1, control = list(hours = 2)))
  e <- ev$score$bins$error_mean_m
  # A first limit between the smallest and the largest error leaves errors
  # below it; they get the class "< limit".
  first <- round(stats::median(e))
  expect_true(any(e < first))
  breaks <- c(first, first + 200, first + 400)
  p <- plot(ev, breaks = breaks)
  shown <- p$data$error_class[p$data$what == "animal-hour"]
  expect_identical(levels(shown), c(sprintf("< %g", first), sprintf("%g-%g", first, first + 200),
                                    sprintf("%g-%g", first + 200, first + 400),
                                    sprintf("> %g", first + 400)))
  expect_length(shown, length(e))
  expect_equal(as.integer(shown), findInterval(e, c(-Inf, breaks)))
  # A first limit of 0 has no class below it
  expect_identical(levels(plot(ev, breaks = c(0, 100))$data$error_class), c("0-100", "> 100"))
  bad <- list("x", c(0, -5), 200, c(0, NA, 400), c(0, Inf), c(400, 200), c(0, 0), numeric())
  for (b in bad) {
    err <- tryCatch(plot(ev, breaks = b), error = identity)
    expect_s3_class(err, "arraydesign_error_input")
    expect_identical(conditionMessage(err),
                     "`breaks` must be at least 2 finite numbers in increasing order.")
    expect_identical(conditionCall(err), quote(plot.array_evaluation(ev, breaks = b)))
  }
})

test_that("plot refusals of `error` and `legend` name the user's call", {
  ev <- suppressMessages(evaluate_array(toy_rx, toy, logistic, n_animals = 1,
                                        seed = 1, control = list(hours = 1)))
  refusals <- list(
    list(quote(plot.array_evaluation(ev, error = "median")),
         "`error` should be one of \"mean\", \"mode\"."),
    list(quote(plot.array_evaluation(ev, error = 1)),
         "`error` must be one string, one of \"mean\", \"mode\"."),
    list(quote(plot.array_evaluation(ev, legend = "yes")),
         "`legend` must be TRUE or FALSE.")
  )
  for (r in refusals) {
    # Dispatch through plot() so the call is the method's, as a user meets it.
    generic <- r[[1]]
    generic[[1]] <- quote(plot)
    err <- tryCatch(eval(generic), error = identity)
    expect_s3_class(err, "arraydesign_error_input")
    expect_identical(conditionMessage(err), r[[2]])
    expect_identical(conditionCall(err), r[[1]])
  }
})

test_that("plot draws every receiver on a wide, short device", {
  # At 874 by 460 pixels, the size of the arraydesignR Score tab, the legend
  # can clip the receivers nearest the top and bottom of the map. Each
  # receiver must change the image as much as the others do.
  skip_if_not_installed("png")
  lake <- terra::unwrap(toy_lake_fixture$grid)
  rx <- place_receivers(lake, n = 5, method = "grid", seed = 1)
  ev <- suppressMessages(evaluate_array(rx, lake, logistic, n_animals = 2,
                                        seed = 1, control = list(hours = 2)))
  render <- function(layout) {
    e <- ev
    e$layout <- layout
    f <- tempfile(fileext = ".png")
    ggplot2::ggsave(f, plot(e), width = 874 / 96, height = 460 / 96, dpi = 96, bg = "white")
    png::readPNG(f)[, , 1:3]
  }
  base <- render(ev$layout[0, ])
  changed <- vapply(seq_len(nrow(ev$layout)), function(i) {
    sum(apply(abs(render(ev$layout[i, ]) - base) > 0.02, c(1, 2), any))
  }, numeric(1))
  expect_gt(length(changed), 2)
  expect_true(all(changed >= 0.6 * max(changed)), label = paste(changed, collapse = ", "))
})

test_that("a positioning model flat in distance scores no bins instead of stopping", {
  flat <- detfield::de_model(function(distance_m) 0 * distance_m + 1)
  expect_warning(
    ev <- evaluate_array(toy_rx, toy, logistic, n_animals = 2,
                         seed = 1, control = list(hours = 3, position_model = flat)),
    "All 6 detected animal-hours are degenerate"
  )
  expect_null(ev$surfaces)
  expect_identical(ev$score$overall$n_bins, 0L)
  expect_identical(ev$settings$skipped, 6L)
  expect_identical(ev$score$overall$n_animal_hours, 6L)
  expect_output(print(ev), "6 degenerate bins left out")
  expect_error(all_degenerate(simpleError("other"), data.frame()), "other")
})

# ---- evaluate_array(): movement and density ----

test_that("movement lowers the hit rate more in a dense array than in a sparse one", {
  skip_on_cran()
  # Still animals hold the hit rate at nominal at 1 and 16 receivers; with the
  # default movement the shortfall grows with density, as regions shrink.
  hit <- function(n_rx, movement) {
    ctrl <- list(hours = 12)
    if (!is.null(movement)) ctrl$movement <- movement
    ev <- suppressMessages(evaluate_array(place_receivers(toy, n = n_rx), toy, logistic,
                                          n_animals = 60, seed = 11, control = ctrl))
    ev$score$overall
  }
  s1 <- hit(1, still)
  s16 <- hit(16, still)
  m1 <- hit(1, NULL)
  m16 <- hit(16, NULL)
  expect_lt(abs(s1$hit_rate - 0.95), 3 * s1$hit_rate_se)
  expect_lt(abs(s16$hit_rate - 0.95), 3 * s16$hit_rate_se)
  expect_gt(s16$hit_rate - m16$hit_rate, s1$hit_rate - m1$hit_rate)
  expect_lt(m16$hit_rate, m1$hit_rate)
})

# ---- evaluate_array(): refusals ----

test_that("bad arguments and settings are refused", {
  ea <- function(...) evaluate_array(toy_rx, toy, logistic, ...)
  expect_error(evaluate_array(toy_rx, toy, list()), "must be a de_model")
  expect_error(ea(n_animals = 0), "`n_animals`")
  expect_error(ea(control = list(k = 1.5)), "`control\\$k`")
  expect_error(ea(control = list(seed = 1)), "`seed` is an argument of evaluate_array")
  expect_error(ea(control = list(seed = 1)), "seed = 1)", fixed = TRUE)
  expect_error(ea(control = list(seed = 123456789)), "seed = 123456789)", fixed = TRUE)
  expect_error(ea(seed = "a"), "`seed` must be one number")
  expect_error(ea(control = 1), "must be a list")
  expect_error(ea(control = list(list(hours = 2))), "needs a name")
  expect_error(ea(control = list(hours = 2, 3)), "needs a name")
  expect_error(ea(control = list(speed = 1)), "Unknown setting in `control`: speed")
  expect_error(ea(control = list(movement = 1)), "movement_model")
  expect_error(ea(control = list(position_model = 1)), "a de_model")
  expect_error(ea(control = list(hours = 0)), "`control\\$hours`")
  expect_error(ea(control = list(level = 1)), "`control\\$level`")
  expect_error(ea(control = list(cores = 0)), "`control\\$cores`")
  expect_error(ea(control = list(start_time = as.POSIXct("2025-01-01 00:30:00", tz = "UTC"))),
               "on the hour")
  expect_error(ea(control = list(movement = crw_model(5, 2, 30, time_step_s = 7))),
               "must divide an hour")
  expect_error(evaluate_array(1:3, toy, logistic), "receiver_layout or a data frame")
})

test_that("every grid refusal through evaluate_array() names the user's call", {
  xy <- as.data.frame(toy_rx)[c("x", "y")]
  square <- function(crs) {
    r <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 100, ymin = 0, ymax = 100,
                     crs = crs)
    terra::values(r) <- 1
    r
  }
  dry <- toy
  terra::values(dry) <- NA
  flat <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 100, ymin = 0, ymax = 200,
                      crs = "EPSG:32617")
  terra::values(flat) <- 1
  ll <- terra::rast(nrows = 10, ncols = 10, xmin = -78, xmax = -77, ymin = 44, ymax = 45)
  terra::values(ll) <- 1
  grids <- list(
    "SpatRaster" = list(matrix(1, 2, 2), "arraydesign_error_grid"),
    "one layer" = list(c(toy, toy), "arraydesign_error_grid"),
    "longitude and latitude" = list(ll, "arraydesign_error_grid"),
    "US survey foot" = list(terra::project(toy, "EPSG:2263", method = "near"),
                            "arraydesign_error_grid"),
    "square" = list(flat, "arraydesign_error_grid"),
    "no CRS" = list(square(""), "arraydesign_error_grid"),
    "cannot be measured" = list(square(paste0(
      'ENGCRS["local",EDATUM["site"],CS[Cartesian,2],',
      'AXIS["x",east,LENGTHUNIT["metre",1]],AXIS["y",north,LENGTHUNIT["metre",1]]]'
    )), "arraydesign_error_grid"),
    "no water cells" = list(dry, "arraydesign_error_grid")
  )
  for (msg in names(grids)) {
    g <- grids[[msg]][[1]]
    err <- tryCatch(evaluate_array(xy, g, logistic), error = identity)
    expect_s3_class(err, grids[[msg]][[2]])
    expect_match(conditionMessage(err), msg, fixed = TRUE, info = msg)
    expect_identical(conditionCall(err), quote(evaluate_array(xy, g, logistic)), info = msg)
  }
  # A whole layout skips placement, and a grid with no CRS is still refused
  nocrs <- toy
  terra::crs(nocrs) <- ""
  err <- tryCatch(evaluate_array(toy_rx, nocrs, logistic), error = identity)
  expect_s3_class(err, "arraydesign_error_grid")
  expect_identical(conditionCall(err), quote(evaluate_array(toy_rx, nocrs, logistic)))
})

# ---- Seeds ----
# The 18 values: vectors, a full .Random.seed, other types, whole numbers, a
# fraction, and values set.seed() refuses (a logical NA, NaN, Inf, 2^31).
seed_values <- function() {
  full <- withr::with_seed(1, get(".Random.seed", envir = globalenv()))
  list(`16 numbers` = c(10403L, 624L, 1:14), `1:16` = 1:16, `.Random.seed` = full,
       `NULL` = NULL, `7` = 7, `"a"` = "a", `c(1, 2)` = c(1, 2), `NA` = NA, `-1` = -1,
       `1.5` = 1.5, `list(1)` = list(1), `1e9` = 1e9, `7L` = 7L, `TRUE` = TRUE,
       `NaN` = NaN, `Inf` = Inf, `2^31` = 2^31, `-0` = -0)
}

test_that("a seed in control stops with advice that set.seed() takes", {
  ea <- function(...) evaluate_array(toy_rx, toy, logistic, ...)
  vals <- seed_values()
  expect_length(vals, 18)
  # A value the argument takes is advised by its whole part, which set.seed()
  # uses.
  advised_same <- c(`7` = "7", `-1` = "-1", `1e9` = "1000000000", `7L` = "7", `-0` = "0",
                    `1.5` = "1")
  held <- c(`16 numbers` = "16 numbers", `1:16` = "16 numbers",
            `.Random.seed` = paste(length(vals$.Random.seed), "numbers"), `NULL` = "NULL",
            `"a"` = "a character of length 1", `c(1, 2)` = "2 numbers",
            `NA` = "a logical of length 1",
            `list(1)` = "a list of length 1", `TRUE` = "a logical of length 1",
            `NaN` = "NaN", `Inf` = "the number Inf",
            `2^31` = "the number 2147483648, outside -2147483647 to 2147483647")
  expect_setequal(c(names(advised_same), names(held)), names(vals))
  for (nm in names(vals)) {
    e <- tryCatch(ea(control = stats::setNames(vals[nm], "seed")), error = function(e) e)
    msg <- conditionMessage(e)
    expect_s3_class(e, "arraydesign_error_input")
    expect_length(strsplit(msg, "\n")[[1]], 1)
    expect_lt(nchar(msg), 200)
    advised <- sub("\\).*$", "", sub("^.*seed = ", "", msg))
    if (nm %in% names(advised_same)) {
      expect_identical(advised, advised_same[[nm]], label = nm)
    } else {
      expect_identical(advised, "1", label = nm)
      expect_match(msg, paste0("`control$seed` holds ", held[[nm]], "."), fixed = TRUE)
    }
    # Advice a user follows must not fail in its turn.
    v <- as.numeric(advised)
    expect_true(arraydesign:::is_seed(v), label = nm)
    expect_no_error(withr::with_seed(v, stats::runif(1)))
  }
  expect_match(conditionMessage(tryCatch(ea(control = list(seed = NA_real_)),
                                         error = function(e) e)),
               "`control$seed` holds NA.", fixed = TRUE)
})

test_that("the advised seed keeps every digit and gives the same stream", {
  advised <- function(s) {
    sub("\\).*$", "", sub("^.*seed = ", "", arraydesign:::seed_advice(s)))
  }
  expect_identical(advised(1234567.5), "1234567")
  expect_identical(advised(12345678.6), "12345678")
  expect_identical(advised(-0.5), "0")
  # Just below a whole number, a 15-digit print would round up to it.
  expect_identical(advised(2147483646.9999995), "2147483646")
  expect_identical(advised(-2147483647), "-2147483647")
  stream <- function(s) withr::with_seed(s, stats::runif(1))
  vals <- withr::with_seed(42, {
    whole <- floor(10^stats::runif(100, 5, log10(2147483646)))
    (whole + stats::runif(100, 0.01, 0.99)) * sample(c(-1, 1), 100, TRUE)
  })
  same <- vapply(vals, function(s) stream(s) == stream(as.numeric(advised(s))), logical(1))
  expect_identical(sum(!same), 0L)
})

test_that("a seed set.seed() would refuse stops classed before set.seed()", {
  ea <- function(...) evaluate_array(toy_rx, toy, logistic, ...)
  vals <- seed_values()
  refused <- c(vals[c("16 numbers", "1:16", ".Random.seed", "\"a\"", "c(1, 2)", "NA",
                      "list(1)", "TRUE", "NaN", "Inf", "2^31")],
               list(`NA_real_` = NA_real_, `-Inf` = -Inf, `-2^31` = -2^31,
                    `2147483647.5` = 2147483647.5))
  stream <- withr::with_seed(3, {
    for (nm in names(refused)) {
      e <- tryCatch(ea(seed = refused[[nm]]), error = function(e) e)
      expect_s3_class(e, "arraydesign_error_input")
      expect_match(conditionMessage(e), "`seed` must be one number from -2147483647 to 2147483647",
                   fixed = TRUE)
      expect_length(strsplit(conditionMessage(e), "\n")[[1]], 1)
      expect_lt(nchar(conditionMessage(e)), 200)
      e2 <- tryCatch(place_receivers(toy, n = 4, method = "random", seed = refused[[nm]]),
                     error = function(e) e)
      expect_s3_class(e2, "arraydesign_error_input")
    }
    stats::runif(1)
  })
  # A refusal comes before set.seed(), so the caller's stream is untouched.
  expect_identical(stream, withr::with_seed(3, stats::runif(1)))
  expect_match(conditionMessage(tryCatch(ea(seed = Inf), error = function(e) e)),
               "`seed` holds the number Inf.", fixed = TRUE)
  expect_match(conditionMessage(tryCatch(ea(seed = NaN), error = function(e) e)),
               "`seed` holds NaN.", fixed = TRUE)
  # The refusal prints the value given with every digit.
  expect_match(conditionMessage(tryCatch(ea(seed = 2147483647.5), error = function(e) e)),
               "`seed` holds the number 2147483647.5, outside", fixed = TRUE)
  expect_match(conditionMessage(tryCatch(ea(seed = -2147483647.5), error = function(e) e)),
               "`seed` holds the number -2147483647.5, outside", fixed = TRUE)
  # The check must not refuse a value set.seed() takes, at either end of its
  # range.
  accepted <- c(vals[c("7", "-1", "1.5", "1e9", "7L", "-0")],
                list(2147483647, -2147483647))
  for (s in accepted) {
    expect_s3_class(place_receivers(toy, n = 4, method = "random", seed = s), "receiver_layout")
  }
})

test_that("a fractional seed gives the result of its whole part, as the help says", {
  run <- function(s) {
    suppressMessages(evaluate_array(toy_rx, toy, logistic, n_animals = 2, seed = s,
                                    control = list(hours = 1)))
  }
  one <- run(1)
  # The results store the seed as given, so that attribute differs.
  expect_equal(run(1.5)$tracks, one$tracks, ignore_attr = "seed", tolerance = 0)
  expect_equal(run(1.9)$detections, one$detections, ignore_attr = "seed", tolerance = 0)
})

test_that("control$tag_model reads one curve of a model with several", {
  m <- detfield::de_model(data.frame(distance_m = rep(c(0, 500, 1000), 2),
                                     eta = c(2, 0, -2, 1, -1, -3),
                                     tag_model = rep(c("A", "B"), each = 3)))
  expect_error(suppressMessages(evaluate_array(toy_rx, toy, m, n_animals = 2,
                                               seed = 1, control = list(hours = 2))),
               "give one tag model")
  ev <- suppressMessages(evaluate_array(toy_rx, toy, m, n_animals = 2,
                                        seed = 1, control = list(hours = 2, tag_model = "B")))
  expect_s3_class(ev, "array_evaluation")
  expect_identical(ev$settings$tag_model, "B")
  one <- detfield::de_model(data.frame(distance_m = c(0, 500, 1000), eta = c(1, -1, -3)))
  ref <- suppressMessages(evaluate_array(toy_rx, toy, one, n_animals = 2,
                                         seed = 1, control = list(hours = 2)))
  expect_equal(ev$field$eta, ref$field$eta)
  expect_error(evaluate_array(toy_rx, toy, m, control = list(tag_model = c("A", "B"))),
               "one tag model name")
})

# ---- The reading of the curve ----

test_that("control$likelihood reaches location_surface() and is recorded", {
  # The second reading is named by detfield; read it from its formals.
  readings <- eval(formals(detfield::location_surface)$likelihood)
  other <- setdiff(readings, "binomial")[1]
  ev <- function(lik) {
    ctrl <- list(hours = 3)
    if (!is.null(lik)) ctrl$likelihood <- lik
    suppressMessages(evaluate_array(toy_rx, toy, logistic, n_animals = 3, seed = 1,
                                    control = ctrl))
  }
  default <- ev(NULL)
  binomial <- ev("binomial")
  weight <- ev(other)
  expect_identical(default$settings$likelihood, "binomial")
  expect_identical(weight$settings$likelihood, other)
  # The binomial default is the same run as asking for it.
  expect_identical(default$surfaces, binomial$surfaces)
  # The same animals and detections, positioned two ways.
  expect_identical(weight$detections, binomial$detections)
  listening <- data.frame(receiver_serial = toy_rx$receiver_serial,
                          station_id = toy_rx$station_id,
                          deploy_datetime = as.POSIXct("2024-12-31", tz = "UTC"),
                          recover_datetime = as.POSIXct("2025-01-03", tz = "UTC"))
  counts <- detfield::bin_detections(binomial$detections$detections, listening,
                                     tag_delay = c(65, 185), bin = "hour")
  direct <- detfield::location_surface(counts, weight$field, overdispersion = "none",
                                       likelihood = other,
                                       control = list(on_degenerate = "skip"))
  expect_equal(weight$surfaces$prob, direct$prob)
  expect_false(isTRUE(all.equal(weight$surfaces$prob, binomial$surfaces$prob)))
  expect_output(print(binomial), "positioned by: binomial grid localization")
  expect_output(print(weight), paste("positioned by:", sub("_", "-", other)))
  expect_error(ev("poisson"), "`control\\$likelihood` must be one of")
})

# ---- Error classes ----

test_that("errors carry arraydesign classes by kind", {
  expect_error(evaluate_array(toy_rx, toy, list()), class = "arraydesign_error_input")
  expect_error(evaluate_array(toy_rx, toy, list()), class = "arraydesign_error")
  expect_error(place_receivers(1:3, n = 2), class = "arraydesign_error_grid")
  expect_error(place_receivers(toy, n = 2, spacing = 1e6, method = "random"),
               class = "arraydesign_error_placement")
})

test_that("a degenerate stop from detfield is caught by its class", {
  counts <- data.frame(animal_id = c("a", "a", "b"), bin_start = c(1, 2, 1))
  by_class <- structure(class = c("detfield_error_degenerate", "error", "condition"),
                        list(message = "every cell has the same likelihood", call = NULL))
  expect_warning(out <- arraydesign:::all_degenerate(by_class, counts),
                 "All 3 detected animal-hours are degenerate")
  expect_null(out)
  expect_error(arraydesign:::all_degenerate(simpleError("another failure"), counts),
               "another failure")
})
