# The oracle: glatos 0.9.8 simulates the same walk, tag and curve with its own
# code (crw_in_polygon, transmit_along_path, detect_transmissions). The two
# chains must give the same detection rate in every 100 m distance band,
# within a binomial 95% interval for the difference, Bonferroni-adjusted over
# the bands. arraydesign runs on a 10 m copy of the lake so that taking the
# probability at the cell centre moves a transmission at most 7 m.
# The test is written against glatos 0.9.8, the floor in DESCRIPTION; with an
# older glatos it skips, and the skip names the installed and the needed
# versions.

test_that("detection rates per distance band agree with the glatos chain", {
  skip_on_cran()
  skip_if_not_installed("glatos", "0.9.8")
  skip_if_not_installed("sf")
  lake <- example_lake()
  poly <- sf::st_cast(sf::st_as_sf(terra::as.polygons(terra::ifel(is.na(lake), NA, 1),
                                                      dissolve = TRUE)), "POLYGON")
  rx <- place_receivers(lake, n = 6)
  curve <- function(d) stats::plogis(2.5 - d / 100)
  n_animals <- 6
  n_steps <- 1000
  start <- c(725300, 4935600)

  glatos_runs <- withr::with_seed(31, lapply(seq_len(n_animals), function(i) {
    path <- glatos::crw_in_polygon(poly, theta = c(0, 10), stepLen = 30, initPos = start,
                                   nsteps = n_steps, sp_out = FALSE, show_progress = FALSE)
    trns <- glatos::transmit_along_path(path, vel = 0.5, delayRng = c(60, 180), burstDur = 5,
                                        sp_out = FALSE)
    dtc <- glatos::detect_transmissions(trns, rx[c("x", "y")], detRngFun = curve,
                                        inputCRS = 32617, sp_out = FALSE, show_progress = FALSE)
    heard <- matrix(FALSE, nrow(trns), nrow(rx))
    heard[cbind(dtc$trns_id, dtc$rec_id)] <- TRUE
    list(d = sqrt(outer(trns$x, rx$x, "-")^2 + outer(trns$y, rx$y, "-")^2), heard = heard)
  }))
  g_d <- do.call(rbind, lapply(glatos_runs, `[[`, "d"))
  g_heard <- do.call(rbind, lapply(glatos_runs, `[[`, "heard"))

  fine <- terra::disagg(lake, 3)
  cells <- which(!is.na(terra::values(fine, mat = FALSE)))
  xy <- terra::xyFromCell(fine, cells)
  dc <- sqrt(outer(xy[, 1], rx$x, "-")^2 + outer(xy[, 2], rx$y, "-")^2)
  field <- list(eta = 2.5 - dc / 100, cells = cells, receivers = rx$receiver_serial, grid = fine)
  m <- movement_model(function(n) rep(30, n), function(n) stats::rnorm(n, 0, 10 * pi / 180),
                      time_step_s = 60)
  tr <- simulate_tracks(fine, m, n_animals, n_steps,
                        track_control(start = data.frame(x = start[1], y = start[2])), seed = 32)
  tx <- simulate_transmissions(tr, seed = 33)
  det <- simulate_detections(tx, field, seed = 34)
  a_d <- sqrt(outer(tx$x, rx$x, "-")^2 + outer(tx$y, rx$y, "-")^2)
  a_heard <- matrix(FALSE, nrow(tx), nrow(rx))
  a_heard[cbind(det$detections$transmission_id,
                match(det$detections$receiver_serial, rx$receiver_serial))] <- TRUE

  # Both tags send about 28.8 transmissions an hour.
  hours <- n_animals * n_steps * 60 / 3600
  expect_equal(nrow(g_d) / hours, 28.8, tolerance = 0.05)
  expect_equal(nrow(a_d) / hours, 28.8, tolerance = 0.05)

  breaks <- c(seq(0, 800, 100), Inf)
  gb <- cut(g_d, breaks, right = FALSE)
  ab <- cut(a_d, breaks, right = FALSE)
  k <- nlevels(gb)
  agree <- vapply(levels(gb), function(b) {
    x <- c(sum(g_heard[gb == b]), sum(a_heard[ab == b]))
    n <- c(sum(gb == b), sum(ab == b))
    ci <- suppressWarnings(stats::prop.test(x, n, conf.level = 1 - 0.05 / k,
                                            correct = FALSE)$conf.int)
    ci[1] <= 0 && ci[2] >= 0
  }, logical(1))
  expect_true(all(agree))
})
