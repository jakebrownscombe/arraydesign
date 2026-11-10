test_that("coverage counts cells above each threshold over all water cells", {
  lake <- example_lake()
  rx <- place_receivers(lake, n = 2)
  field <- make_field(lake, rx)
  # Known answer: receiver 1 gives p = 0.9 on the first 100 cells, receiver 2
  # gives p = 0.6 on the next 50; every other cell is at p = 0.01.
  field$eta[] <- stats::qlogis(0.01)
  field$eta[1:100, 1] <- stats::qlogis(0.9)
  field$eta[101:150, 2] <- stats::qlogis(0.6)
  cs <- coverage_summary(field, thresholds = c(0.5, 0.8, 0.95))
  expect_equal(cs$n_cells, c(150L, 100L, 0L))
  expect_equal(cs$n_water_cells, rep(2152L, 3))
  expect_equal(cs$share, c(150, 100, 0) / 2152)
  expect_equal(cs$n_receivers, rep(2L, 3))
  one <- coverage_summary(field, listening = "R002")
  expect_equal(one$n_cells, c(50L, 0L))
  expect_equal(one$n_receivers, c(1L, 1L))
  none <- coverage_summary(field, listening = character())
  expect_equal(none$n_cells, c(0L, 0L))
})

test_that("cells the field leaves out count as not covered", {
  lake <- example_lake()
  field <- make_field(lake, place_receivers(lake, n = 3))
  half <- field
  keep <- seq_len(1000)
  half$eta <- field$eta[keep, ]
  half$cells <- field$cells[keep]
  full <- coverage_summary(field, thresholds = 1e-12)
  part <- coverage_summary(half, thresholds = 1e-12)
  expect_equal(full$n_cells, 2152L)
  expect_equal(part$n_cells, 1000L)
  expect_equal(part$n_water_cells, 2152L)
})

test_that("bad inputs stop; a missing value is never filled", {
  lake <- example_lake()
  field <- make_field(lake, place_receivers(lake, n = 3))
  expect_error(coverage_summary(field, thresholds = 1.5), "above 0 and at most 1")
  # A threshold of 0 would count every cell as covered, even with no receivers.
  expect_error(coverage_summary(field, thresholds = 0), "above 0")
  none <- coverage_summary(field, listening = character(0), thresholds = 1e-12)
  expect_equal(c(none$n_cells, none$n_receivers), c(0L, 0L))
  expect_error(coverage_summary(field, thresholds = numeric()), "probabilities")
  expect_error(coverage_summary(field, listening = "R999"), "R999")
  holed <- field
  holed$eta[5, 1] <- NA
  expect_error(coverage_summary(holed), "no probability is assumed")
  expect_equal(coverage_summary(holed, listening = c("R002", "R003"))$n_receivers, c(2L, 2L))
  dry <- field
  dry$cells[1] <- 1L
  expect_error(coverage_summary(dry), "land or outside")
})

test_that("a station-keyed field takes stations or serials in `listening`", {
  grid <- terra::unwrap(toy_lake_fixture$grid)
  dep <- toy_lake_fixture$receivers
  field <- suppressMessages(detfield::de_field(detfield::water_distances(grid, dep),
                                               detfield::de_logistic(d50 = 250, d05 = 600)))
  by_station <- coverage_summary(field, listening = c("S1", "S2"))
  by_serial <- coverage_summary(field, listening = c("480001", "480002"))
  expect_identical(by_station, by_serial)
  expect_equal(by_station$n_receivers, c(2L, 2L))
  expect_error(coverage_summary(field, listening = "R999"), "does not hold: R999")
})
