test_that("detection rates per distance band match a known curve (binomial 95%)", {
  # One set of transmissions, detected under 50 seeds. In each 100 m band the
  # exact binomial 95% interval of the observed rate should hold the mean
  # true probability in about 95% of seeds; the trials in a band have unequal
  # probabilities, which makes the interval slightly conservative.
  lake <- example_lake()
  rx <- place_receivers(lake, n = 6)
  field <- make_field(lake, rx, a = 2.5, s = 100)
  tr <- simulate_tracks(lake, crw_model(30, 10, 30), n_animals = 20, n_steps = 2000, seed = 21)
  tx <- simulate_transmissions(tr, seed = 22)
  row <- match(tx$cell, field$cells)
  dist <- matrix(field$distances$distance_m, ncol = 6)[row, ]
  p <- stats::plogis(2.5 - dist / 100)
  band <- cut(dist, c(seq(0, 900, 100), Inf), right = FALSE)
  covered <- vapply(1:50, function(seed) {
    det <- simulate_detections(tx, field, seed = seed)
    heard <- matrix(FALSE, nrow(tx), 6)
    heard[cbind(det$detections$transmission_id,
                match(det$detections$receiver_serial, rx$receiver_serial))] <- TRUE
    if (seed == 1) {
      expect_equal(sum(heard), nrow(det$detections))
      expect_equal(det$transmissions$n_detections, rowSums(heard))
    }
    vapply(levels(band), function(b) {
      sel <- band == b
      ci <- stats::binom.test(sum(heard[sel]), sum(sel))$conf.int
      mean(p[sel]) >= ci[1] && mean(p[sel]) <= ci[2]
    }, logical(1))
  }, logical(nlevels(band)))
  expect_gte(mean(covered), 0.93)
  expect_true(all(rowMeans(covered) >= 0.86))
})

test_that("detections carry the contract columns, ordered by time", {
  lake <- example_lake()
  rx <- place_receivers(lake, n = 4)
  tr <- simulate_tracks(lake, crw_model(30, 10, 30), n_animals = 2, n_steps = 200, seed = 1)
  det <- simulate_detections(simulate_transmissions(tr, seed = 2), make_field(lake, rx),
                             seed = 3)
  expect_s3_class(det, "sim_detections")
  d <- det$detections
  expect_named(d, c("datetime", "receiver_serial", "transmitter", "animal_id",
                    "transmission_id"))
  expect_s3_class(d$datetime, "POSIXct")
  expect_equal(attr(d$datetime, "tzone"), "UTC")
  expect_type(d$receiver_serial, "character")
  expect_type(d$transmitter, "character")
  expect_type(d$animal_id, "character")
  expect_false(is.unsorted(as.numeric(d$datetime)))
  expect_output(print(det), "detections of")
})

test_that("probabilities of 0 and 1 give no detections and every detection", {
  lake <- example_lake()
  rx <- place_receivers(lake, n = 3)
  tr <- simulate_tracks(lake, crw_model(30, 10, 30), n_steps = 100, seed = 1)
  tx <- simulate_transmissions(tr, seed = 2)
  never <- make_field(lake, rx)
  never$eta[] <- -Inf
  always <- never
  always$eta[] <- Inf
  expect_equal(nrow(simulate_detections(tx, never, seed = 1)$detections), 0L)
  expect_equal(nrow(simulate_detections(tx, always, seed = 1)$detections), 3L * nrow(tx))
})

test_that("a missing probability or an unknown cell stops; nothing is assumed", {
  lake <- example_lake()
  rx <- place_receivers(lake, n = 3)
  tr <- simulate_tracks(lake, crw_model(30, 10, 30), n_steps = 100, seed = 1)
  tx <- simulate_transmissions(tr, seed = 2)
  field <- make_field(lake, rx)
  holed <- field
  holed$eta[match(tx$cell[1], field$cells), 2] <- NA
  expect_error(simulate_detections(tx, holed), "no probability is assumed")
  short <- field
  keep <- field$cells != tx$cell[1]
  short$eta <- field$eta[keep, ]
  short$cells <- field$cells[keep]
  expect_error(simulate_detections(tx, short), "cells the field does not hold")
})

test_that("a malformed field or a different grid is refused", {
  lake <- example_lake()
  rx <- place_receivers(lake, n = 3)
  tr <- simulate_tracks(lake, crw_model(30, 10, 30), n_steps = 20, seed = 1)
  tx <- simulate_transmissions(tr, seed = 2)
  field <- make_field(lake, rx)
  expect_error(simulate_detections(tr, field), "simulate_transmissions")
  expect_error(simulate_detections(tx, list(eta = 1)), "must be a de_field")
  bad <- field
  bad$eta <- as.data.frame(bad$eta)
  expect_error(simulate_detections(tx, bad), "numeric matrix")
  bad <- field
  bad$receivers <- bad$receivers[1:2]
  expect_error(simulate_detections(tx, bad), "3 receivers|2 receivers")
  bad <- field
  bad$receivers[2] <- bad$receivers[1]
  expect_error(simulate_detections(tx, bad), "must not repeat")
  moved <- field
  moved$grid <- terra::shift(lake, dx = 15)
  expect_error(simulate_detections(tx, moved), "different grids")
  # Same numbers, other CRS: the cell ids would name other places.
  other <- field
  other$grid <- lake
  terra::crs(other$grid) <- "EPSG:32618"
  expect_error(simulate_detections(tx, other), "different CRSs")
})

test_that("a seed reproduces detections and leaves the caller's stream alone", {
  lake <- example_lake()
  rx <- place_receivers(lake, n = 4)
  tr <- simulate_tracks(lake, crw_model(30, 10, 30), n_animals = 2, n_steps = 200, seed = 1)
  tx <- simulate_transmissions(tr, seed = 2)
  field <- make_field(lake, rx)
  set.seed(3)
  before <- .Random.seed
  a <- simulate_detections(tx, field, seed = 4)
  expect_identical(.Random.seed, before)
  expect_identical(a, simulate_detections(tx, field, seed = 4))
  expect_false(identical(a$detections, simulate_detections(tx, field, seed = 5)$detections))
})

test_that("draws in blocks give the same successes as one block", {
  p <- matrix(0, 1001, 5000)
  p[cbind(c(1, 1000, 1001), c(1, 4999, 5000))] <- 1
  hits <- arraydesign:::bernoulli_hits(p)
  expect_equal(unname(hits[order(hits[, 1]), , drop = FALSE]),
               cbind(c(1, 1000, 1001), c(1, 4999, 5000)))
})

test_that("a transmission's draws do not depend on how many transmissions follow it", {
  lake <- example_lake()
  rx <- place_receivers(lake, n = 6)
  field <- make_field(lake, rx)
  tr <- simulate_tracks(lake, crw_model(30, 10, 30), n_animals = 4, n_steps = 300, seed = 1)
  tx <- simulate_transmissions(tr, seed = 2)
  first <- tx[seq_len(200), ]
  attributes(first) <- c(attributes(first)[c("names", "row.names", "class")],
                         attributes(tx)[c("tag_delay", "seed", "grid", "epsg")])
  alone <- simulate_detections(first, field, seed = 5)$detections
  full <- simulate_detections(tx, field, seed = 5)$detections
  inside <- full[full$transmission_id <= 200, ]
  rownames(inside) <- NULL
  expect_identical(alone, inside)
  # Row by row against a plain draw of n x m uniforms.
  p <- stats::plogis(field$eta[match(tx$cell, field$cells), ])
  u <- withr::with_seed(5, matrix(stats::runif(length(p)), nrow(p), byrow = TRUE))
  expect_equal(nrow(full), sum(u < p))
})

# ---- A field keyed by station, built by detfield from deployments ----
# On a station-keyed field, simulate_detections() writes the receiver's serial
# as receiver_serial. A station id in that column would make bin_detections()
# mark every detection not listening.

toy_start <- as.POSIXct("2024-06-01 06:00:00", tz = "UTC")

toy_field <- function(station_hour = 0, deployments = toy_lake_fixture$receivers) {
  grid <- terra::unwrap(toy_lake_fixture$grid)
  model <- detfield::de_model(detfield::de_logistic(d50 = 250, d05 = 600),
                              sigma = list(station_hour = station_hour))
  dist <- suppressWarnings(detfield::water_distances(grid, deployments))
  list(grid = grid, field = suppressMessages(detfield::de_field(dist, model)))
}

toy_transmissions <- function(grid, n_animals = 10, hours = 6, movement = crw_model(3, 1, 40)) {
  tr <- simulate_tracks(grid, movement, n_animals = n_animals, n_steps = hours * 60,
                        control = track_control(start_time = toy_start), seed = 1)
  simulate_transmissions(tr, tag_delay = c(65, 185), seed = 2)
}

test_that("the round trip through detfield loses no detection and positions well", {
  dep <- toy_lake_fixture$receivers
  toy <- toy_field()
  expect_identical(toy$field$key, "station_id")
  tx <- toy_transmissions(toy$grid)
  det <- simulate_detections(tx, toy$field, seed = 3)
  d <- det$detections
  expect_named(d, c("datetime", "receiver_serial", "station_id", "transmitter", "animal_id",
                    "transmission_id"))
  expect_identical(d$receiver_serial, dep$receiver_serial[match(d$station_id, dep$station_id)])
  counts <- detfield::bin_detections(d, dep, tag_delay = c(65, 185))
  expect_true(all(counts$listening))
  expect_equal(attr(counts, "not_listening"), 0)
  expect_equal(sum(counts$y), nrow(d))
  surfaces <- detfield::location_surface(counts, toy$field, overdispersion = "none")
  expect_equal(NROW(attr(surfaces, "not_listening_bins")), 0L)
  score <- score_surfaces(surfaces, det$transmissions)
  expect_equal(score$overall$n_bins, nrow(unique(counts[c("animal_id", "bin_start")])))
  expect_lt(score$overall$median_error_mode_m, 200)
})

test_that("a station's detections go to the receiver deployed there at the time", {
  dep <- toy_lake_fixture$receivers
  # Station S1 swaps receivers at 09:00 with a 30 min overlap.
  swap <- dep[dep$station_id == "S1", ]
  swap$receiver_serial <- "480099"
  swap$deploy_datetime <- toy_start + 3 * 3600 - 1800
  first <- dep
  first$recover_datetime[first$station_id == "S1"] <- toy_start + 3 * 3600
  deployments <- rbind(first, swap)
  toy <- toy_field(deployments = deployments)
  tx <- toy_transmissions(toy$grid)
  d <- simulate_detections(tx, toy$field, seed = 3)$detections
  s1 <- d[d$station_id == "S1", ]
  expect_gt(nrow(s1), 0)
  later <- s1$datetime >= toy_start + 3 * 3600 - 1800
  expect_true(all(s1$receiver_serial[later] == "480099"))
  expect_true(all(s1$receiver_serial[!later] == "480001"))
  counts <- detfield::bin_detections(d, deployments, tag_delay = c(65, 185))
  expect_true(all(counts$listening))
  expect_equal(sum(counts$y), nrow(d))
  # A track outside every deployment makes no trial rather than invent a
  # receiver: 5 stations times every transmission.
  late <- tx
  late$datetime <- late$datetime + 5 * 86400
  expect_message(none <- simulate_detections(late, toy$field, seed = 3),
                 paste(5 * nrow(tx), "of", 5 * nrow(tx), "transmission and station trials"))
  expect_equal(nrow(none$detections), 0L)
  expect_equal(attr(none, "trials_not_listening"), 5 * nrow(tx))
  expect_true(all(none$transmissions$n_detections == 0L))
})

# ---- Stations deployed late or recovered early ----
# A station with no receiver in the water makes no trial, by
# bin_detections()'s rule: a receiver in the water, and min_cover of the UTC
# hour covered. Each case below moves one end of S1: deployment at 06:40 or
# 06:20, recovery at 08:20 or 11:50.

toy_at <- function(h, m = 0) toy_start + h * 3600 + m * 60

s1_changed <- function(deploy = NULL, recover = NULL) {
  dep <- toy_lake_fixture$receivers
  s1 <- dep$station_id == "S1"
  if (!is.null(deploy)) dep$deploy_datetime[s1] <- deploy
  if (!is.null(recover)) dep$recover_datetime[s1] <- recover
  dep
}

test_that("a station deployed late or recovered early makes no trial while not listening", {
  full <- toy_field()
  tx <- toy_transmissions(full$grid)
  base <- simulate_detections(tx, full$field, seed = 3)$detections
  # 1,364 detections, with distance along detfield's exact path through water.
  expect_equal(nrow(base), 1364L)
  t <- tx$datetime
  # The trials S1 should not make, from the times alone (min_cover 0.5): late
  # by 40 min loses the whole first hour; late by 20 min loses only the 20 min
  # with no receiver; recovered at 08:20 loses 08:00 onwards (20 min of the
  # hour is under half); recovered at 11:50 keeps the hour and loses 10 min.
  cases <- list(
    list(dep = s1_changed(deploy = toy_at(0, 40)), off = t < toy_at(1)),
    list(dep = s1_changed(deploy = toy_at(0, 20)), off = t < toy_at(0, 20)),
    list(dep = s1_changed(recover = toy_at(2, 20)), off = t >= toy_at(2)),
    list(dep = s1_changed(recover = toy_at(5, 50)), off = t >= toy_at(5, 50))
  )
  for (cs in cases) {
    toy <- toy_field(deployments = cs$dep)
    expect_message(det <- simulate_detections(tx, toy$field, seed = 3),
                   paste(sum(cs$off), "of", 5 * nrow(tx)))
    expect_equal(attr(det, "trials_not_listening"), sum(cs$off))
    d <- det$detections
    # The draws of the other trials do not move: the result is the full
    # array's detections less those of S1 at the times it is off.
    gone <- base$station_id == "S1" & base$transmission_id %in% tx$transmission_id[cs$off]
    expect_gt(sum(gone), 0)
    kept <- base[!gone, ]
    rownames(kept) <- NULL
    expect_identical(d, kept)
    expect_equal(det$transmissions$n_detections,
                 tabulate(d$transmission_id, nbins = nrow(tx)))
    # Round trip: every detection in a listening row, none lost.
    counts <- detfield::bin_detections(d, cs$dep, tag_delay = c(65, 185))
    expect_true(all(counts$listening))
    expect_equal(attr(counts, "not_listening"), 0)
    expect_equal(sum(counts$y), nrow(d))
    surfaces <- detfield::location_surface(counts, toy$field, overdispersion = "none")
    expect_equal(NROW(attr(surfaces, "not_listening_bins")), 0L)
  }
})

test_that("min_cover is the share of the hour bin_detections() asks for", {
  late <- s1_changed(deploy = toy_at(0, 40))
  toy <- toy_field(deployments = late)
  tx <- toy_transmissions(toy$grid, n_animals = 4, hours = 3)
  t <- tx$datetime
  # At 0 the 20 min in the water count; at 1 the first hour of a receiver
  # deployed at 06:20 does not.
  expect_message(zero <- simulate_detections(tx, toy$field, min_cover = 0, seed = 3),
                 "with no receiver deployed")
  expect_equal(attr(zero, "trials_not_listening"), sum(t < toy_at(0, 40)))
  counts <- detfield::bin_detections(zero$detections, late, tag_delay = c(65, 185),
                                     min_cover = 0)
  expect_true(all(counts$listening))
  expect_equal(sum(counts$y), nrow(zero$detections))
  early <- s1_changed(deploy = toy_at(0, 20))
  toy20 <- toy_field(deployments = early)
  expect_message(one <- simulate_detections(tx, toy20$field, min_cover = 1, seed = 3),
                 "less than `min_cover` \\(1\\)")
  expect_equal(attr(one, "trials_not_listening"), sum(t < toy_at(1)))
  counts1 <- detfield::bin_detections(one$detections, early, tag_delay = c(65, 185),
                                      min_cover = 1)
  expect_true(all(counts1$listening))
  expect_equal(sum(counts1$y), nrow(one$detections))
  for (bad in list(-0.1, 1.5, c(0.2, 0.3), NA_real_, "half")) {
    expect_error(simulate_detections(tx, toy$field, min_cover = bad), "`min_cover` must be")
  }
})

test_that("a station with no deployment in the lookup never listens", {
  toy <- toy_field()
  tx <- toy_transmissions(toy$grid, n_animals = 4, hours = 3)
  field <- toy$field
  field$lookup <- field$lookup[field$lookup$station_id != "S2", ]
  expect_message(det <- simulate_detections(tx, field, seed = 3), "Stations: S2\\.")
  expect_equal(attr(det, "trials_not_listening"), nrow(tx))
  expect_false("S2" %in% det$detections$station_id)
  # A full lookup leaves nothing out and asks for no message.
  expect_silent(all_in <- simulate_detections(tx, toy$field, seed = 3))
  expect_equal(attr(all_in, "trials_not_listening"), 0L)
})

test_that("with no lookup the serial is NA and detfield matches by station", {
  dep <- toy_lake_fixture$receivers
  bare <- dep[c("station_id", "x", "y", "deploy_z")]
  toy <- toy_field(deployments = bare)
  expect_null(toy$field$lookup)
  tx <- toy_transmissions(toy$grid, n_animals = 4, hours = 3)
  d <- simulate_detections(tx, toy$field, seed = 3)$detections
  expect_true(all(is.na(d$receiver_serial)))
  listening <- dep
  listening$receiver_serial <- NA_character_
  counts <- suppressMessages(detfield::bin_detections(d, listening, tag_delay = c(65, 185)))
  expect_true(all(counts$listening))
  expect_equal(sum(counts$y), nrow(d))
  # A station with two receivers and no times has no one receiver.
  two <- dep[c("receiver_serial", "station_id", "x", "y", "deploy_z")]
  two <- rbind(two, transform(two[1, ], receiver_serial = "480099"))
  toy2 <- toy_field(deployments = two)
  expect_error(simulate_detections(tx, toy2$field, seed = 3), "no deployment times")
})

# ---- Overdispersion ----
# The field's station-hour spread, field$model$sigma, changes the detections.

test_that("the field's station-hour spread changes the detections", {
  toy0 <- toy_field(0)
  toy2 <- toy_field(2.55)
  tx <- toy_transmissions(toy0$grid)
  a <- simulate_detections(tx, toy0$field, seed = 3)
  expect_message(b <- simulate_detections(tx, toy2$field, seed = 3),
                 "standard deviation\\s+2.55")
  expect_false(identical(a$detections, b$detections))
  expect_identical(attr(a, "overdispersion"), "none")
  expect_identical(attr(b, "overdispersion"), "station_hour")
  expect_equal(attr(b, "station_hour"), 2.55)
  # "none" is the binomial whatever the field holds, and asks for no message.
  expect_silent(n <- simulate_detections(tx, toy2$field, overdispersion = "none", seed = 3))
  expect_identical(n$detections, a$detections)
  expect_silent(simulate_detections(tx, toy2$field, overdispersion = "station_hour", seed = 3))
  expect_error(simulate_detections(tx, toy2$field, overdispersion = "day"),
               "should be one of", class = "arraydesign_error_input")
})

test_that("a field with no station-hour spread gives binomial draws without a message", {
  # ?simulate_detections, Details: a model with no such value (NA, absent or
  # 0) gives independent binomial draws.
  toy0 <- toy_field(0)
  tx <- toy_transmissions(toy0$grid)
  a <- simulate_detections(tx, toy0$field, seed = 3)
  no_spread <- toy0$field
  no_spread$model$sigma$station_hour <- NA_real_
  absent <- toy0$field
  absent$model$sigma$station_hour <- NULL
  for (f in list(no_spread, absent)) {
    expect_silent(d <- simulate_detections(tx, f, seed = 3))
    expect_identical(d$detections, a$detections)
    expect_identical(attr(d, "overdispersion"), "none")
    expect_identical(attr(d, "station_hour"), 0)
  }
})

test_that("the station-hour effect is shared by every animal within an hour", {
  # At a logit of 0 and a standard deviation of 500, a receiver-hour falls
  # between p = 0.01 and 0.99 with probability 0.007, so nearly all
  # transmissions in one are heard or none are.
  lake <- example_lake()
  rx <- place_receivers(lake, n = 4)
  field <- make_field(lake, rx)
  field$eta[] <- 0
  tr <- simulate_tracks(lake, crw_model(30, 10, 30), n_animals = 6, n_steps = 24 * 60, seed = 1)
  tx <- simulate_transmissions(tr, seed = 2)
  hour <- floor(as.numeric(tx$datetime) / 3600)
  share <- function(det) {
    heard <- matrix(FALSE, nrow(tx), nrow(rx))
    heard[cbind(det$detections$transmission_id,
                match(det$detections$receiver_serial, rx$receiver_serial))] <- TRUE
    as.vector(apply(heard, 2, function(h) tapply(h, hour, mean)))
  }
  field$model <- list(sigma = list(station_hour = 500))
  wide <- share(suppressMessages(simulate_detections(tx, field, seed = 3)))
  expect_gte(mean(wide < 0.02 | wide > 0.98), 0.95)
  expect_true(any(wide < 0.02) && any(wide > 0.98))
  # Without the effect every receiver-hour sits near 0.5.
  flat <- share(simulate_detections(tx, field, overdispersion = "none", seed = 3))
  expect_true(all(flat > 0.2 & flat < 0.8))
  # A spread of 1 gives receiver-hour shares more variable than binomial.
  field$model$sigma$station_hour <- 1
  one <- share(suppressMessages(simulate_detections(tx, field, seed = 3)))
  expect_gt(stats::var(one), 3 * stats::var(flat))
})

test_that("one receiver per station needs no times, and an open deployment covers", {
  dep <- toy_lake_fixture$receivers
  untimed <- toy_field(deployments = dep[c("receiver_serial", "station_id", "x", "y",
                                           "deploy_z")])
  expect_false("deploy_datetime" %in% names(untimed$field$lookup))
  tx <- toy_transmissions(untimed$grid, n_animals = 4, hours = 3)
  d <- simulate_detections(tx, untimed$field, seed = 3)$detections
  expect_identical(d$receiver_serial, dep$receiver_serial[match(d$station_id, dep$station_id)])
  open <- dep
  open$recover_datetime <- as.POSIXct(NA, tz = "UTC")
  still_in <- toy_field(deployments = open)
  d2 <- simulate_detections(tx, still_in$field, seed = 3)$detections
  expect_identical(d2$receiver_serial, d$receiver_serial)
})
