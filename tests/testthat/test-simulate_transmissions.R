# Expected transmissions per hour for a start-to-start delay range: one per
# mean interval, 3600 / mean(tag_delay).
per_hour <- function(tag_delay) 3600 / mean(tag_delay)

test_that("transmissions per hour follow the tag delay, not the time step", {
  tr <- simulate_tracks(example_lake(), crw_model(30, 10, 30, time_step_s = 600), n_animals = 10,
                        n_steps = 600, seed = 1)
  hours <- 10 * 600 * 600 / 3600
  for (d in list(c(65, 185), c(125, 245), c(35, 95))) {
    tx <- simulate_transmissions(tr, tag_delay = d, seed = 2)
    mu <- mean(d)
    # Renewal process: the count over T seconds has variance T var / mu^3.
    sd_count <- sqrt(hours * 3600 * diff(d)^2 / 12 / mu^3)
    expect_lt(abs(nrow(tx) - hours * per_hour(d)), 4 * sd_count)
  }
  expect_equal(per_hour(c(65, 185)), 28.8)
  # The same tag on a 60 s walk sends the same number per hour.
  fast <- simulate_tracks(example_lake(), crw_model(30, 10, 30, time_step_s = 60), n_animals = 10,
                          n_steps = 6000, seed = 1)
  n_fast <- nrow(simulate_transmissions(fast, seed = 2))
  expect_lt(abs(n_fast / hours - 28.8), 4 * sqrt(hours * 3600 * 120^2 / 12 / 125^3) / hours)
})

test_that("each transmission lies on the step that holds it, on water", {
  lake <- example_lake()
  tr <- simulate_tracks(lake, crw_model(60, 20, 40), n_animals = 3, n_steps = 300, seed = 3)
  tx <- simulate_transmissions(tr, seed = 4)
  expect_false(anyNA(terra::extract(lake, cbind(tx$x, tx$y))[, 1]))
  for (i in sample(nrow(tx), 50)) {
    p <- tr[tr$path_id == tx$path_id[i], ]
    k <- findInterval(as.numeric(tx$datetime[i]), as.numeric(p$time))
    a <- c(p$x[k], p$y[k])
    b <- c(p$x[k + 1], p$y[k + 1])
    f <- as.numeric(difftime(tx$datetime[i], p$time[k], units = "secs")) / 60
    expect_equal(c(tx$x[i], tx$y[i]), a + f * (b - a), tolerance = 1e-9)
  }
  expect_equal(tx$cell, terra::cellFromXY(lake, cbind(tx$x, tx$y)))
})

test_that("the table carries contract columns and the state of the step", {
  p <- matrix(c(0.9, 0.2, 0.1, 0.8), 2)
  tr <- simulate_tracks(example_lake(), crw_model(30, 10, 30, transition = p), n_animals = 12,
                        n_steps = 100, seed = 1)
  tx <- simulate_transmissions(tr, seed = 2)
  expect_s3_class(tx, "sim_transmissions")
  expect_named(tx, c("transmission_id", "path_id", "animal_id", "transmitter", "datetime",
                     "x", "y", "cell", "state", "min_delay", "max_delay"))
  expect_equal(tx$transmission_id, seq_len(nrow(tx)))
  expect_equal(attr(tx$datetime, "tzone"), "UTC")
  expect_type(tx$animal_id, "character")
  expect_equal(unique(tx$transmitter)[1], "SIM-0001")
  expect_setequal(unique(tx$state), 1:2)
  resting <- tx[tx$state == 2, ]
  k <- vapply(seq_len(nrow(resting)), function(i) {
    p <- tr[tr$path_id == resting$path_id[i], ]
    p$state[findInterval(as.numeric(resting$datetime[i]), as.numeric(p$time)) + 1]
  }, integer(1))
  expect_true(all(k == 2L))
  expect_equal(attr(tx, "tag_delay"), c(65, 185))
  expect_true(all(tx$min_delay == 65 & tx$max_delay == 185))
  expect_equal(attr(tx, "epsg"), 32617L)
})

test_that("a seed reproduces transmissions, and bad arguments stop", {
  tr <- simulate_tracks(example_lake(), crw_model(30, 10, 30), n_animals = 2, n_steps = 100,
                        seed = 1)
  set.seed(1)
  before <- .Random.seed
  a <- simulate_transmissions(tr, seed = 9)
  expect_identical(.Random.seed, before)
  expect_identical(a, simulate_transmissions(tr, seed = 9))
  expect_false(identical(a$datetime, simulate_transmissions(tr, seed = 10)$datetime))
  expect_error(simulate_transmissions(data.frame()), "simulate_tracks")
  expect_error(simulate_transmissions(tr, tag_delay = c(60, 120, 180)), "one or two numbers")
  expect_error(simulate_transmissions(tr, tag_delay = c(180, 60)), "shortest delay first")
  expect_error(simulate_transmissions(tr, tag_delay = c(0, 10)), "above 0")
  expect_error(simulate_transmissions(tr, tag_delay = NA_real_), "above 0")
  expect_error(simulate_transmissions(tr, tag_delay = "120"), "one or two numbers")
})

test_that("a short smallest interval does not size a huge vector", {
  # Draws sized by the smallest interval would ask for 8.6e11 numbers at a
  # smallest interval of 1e-6 s.
  tr <- simulate_tracks(example_lake(), crw_model(30, 10, 30, time_step_s = 600), n_animals = 2,
                        n_steps = 60, seed = 1)
  hours <- 2 * 60 * 600 / 3600
  tx <- simulate_transmissions(tr, tag_delay = c(1e-6, 10), seed = 2)
  sd_count <- sqrt(hours * 3600 * 10^2 / 12 / 5^3)
  expect_lt(abs(nrow(tx) - hours * 720), 4 * sd_count)
  long <- simulate_tracks(example_lake(), crw_model(30, 10, 30, time_step_s = 86400), n_steps = 10,
                          seed = 1)
  tx <- simulate_transmissions(long, tag_delay = c(1e-6, 10), seed = 2)
  expect_lt(abs(nrow(tx) / 172800 - 1), 0.01)
})

test_that("one number is a fixed delay", {
  tr <- simulate_tracks(example_lake(), crw_model(30, 10, 30), n_animals = 1, n_steps = 120,
                        seed = 1)
  tx <- simulate_transmissions(tr, tag_delay = 120, seed = 2)
  expect_equal(unique(round(diff(as.numeric(tx$datetime)), 6)), 120)
  expect_equal(attr(tx, "tag_delay"), c(120, 120))
})

test_that("the hourly count matches detfield's n for the same tag_delay", {
  # Both read tag_delay as start-to-start, so a delay of 120 s gives 30
  # transmissions an hour.
  hours <- 200
  tr <- simulate_tracks(example_lake(), crw_model(30, 10, 30, time_step_s = 600), n_animals = 1,
                        n_steps = hours * 6, seed = 1)
  for (d in list(c(60, 180), c(90, 150), 120)) {
    tx <- simulate_transmissions(tr, tag_delay = d, seed = 2)
    one <- data.frame(datetime = tx$datetime[1], receiver_serial = "R1",
                      transmitter = tx$transmitter[1], animal_id = tx$animal_id[1])
    listening <- data.frame(receiver_serial = "R1", date = as.Date(tx$datetime[1]))
    n_num <- detfield::bin_detections(one, listening, tag_delay = d)$n
    tags <- unique(tx[c("transmitter", "min_delay", "max_delay")])
    n_tab <- detfield::bin_detections(one, listening, tag_delay = tags)$n
    expect_equal(n_num, 30L)
    expect_equal(n_tab, 30L)
    dd <- range(d)
    # Renewal process: hourly rate SE = sqrt(3600 var / mu^3 / hours).
    se <- sqrt(3600 * diff(dd)^2 / 12 / mean(dd)^3 / hours)
    rate <- nrow(tx) / (hours * 6 * 600 / 3600)
    expect_lt(abs(rate - n_num), max(4 * se, 0.2))
  }
})
