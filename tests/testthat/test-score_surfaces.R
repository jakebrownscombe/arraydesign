# ---- Fixture: surfaces of simulated animals on the example lake ----
# Built once; every test reads it.
score_fixture <- local({
  lake <- example_lake()
  rx <- place_receivers(lake, n = 6)
  tr <- simulate_tracks(lake, crw_model(20, 10, 40), n_animals = 3, n_steps = 180,
                        control = track_control(start_time = as.POSIXct("2025-01-01", tz = "UTC")),
                        seed = 11)
  tx <- simulate_transmissions(tr, seed = 12)
  dist <- detfield::water_distances(lake, rx[c("receiver_serial", "x", "y")])
  field <- detfield::de_field(dist, detfield::de_logistic(250, 600))
  det <- simulate_detections(tx, field, seed = 13)
  listening <- data.frame(receiver_serial = rx$receiver_serial, date = as.Date("2025-01-01"))
  counts <- suppressMessages(detfield::bin_detections(det$detections, listening,
                                                      tag_delay = 125))
  surfaces <- detfield::location_surface(counts, field, overdispersion = "none")
  list(lake = lake, truth = det$transmissions, surfaces = surfaces)
})

# ---- score_surfaces(): known answers from an independent oracle ----

test_that("errors, coverage and areas match a recomputation with terra", {
  fx <- score_fixture
  sc <- score_surfaces(fx$surfaces, fx$truth, level = 0.9)
  s <- as.data.frame(unclass(fx$surfaces)[c("animal_id", "bin_start", "cell", "prob")])
  region <- detfield::credible_region(fx$surfaces, level = 0.9)
  expect_equal(nrow(sc$bins), nrow(attr(fx$surfaces, "bins")))
  for (i in seq_len(nrow(sc$bins))) {
    a <- sc$bins$animal_id[i]
    b0 <- sc$bins$bin_start[i]
    si <- s[s$animal_id == a & s$bin_start == b0, ]
    xy <- terra::xyFromCell(fx$lake, si$cell)
    mode_xy <- xy[which.max(si$prob), ]
    mean_xy <- colSums(xy * si$prob) / sum(si$prob)
    ti <- fx$truth[fx$truth$animal_id == a & fx$truth$datetime >= b0 &
                     fx$truth$datetime < b0 + 3600, ]
    cen <- c(mean(ti$x), mean(ti$y))
    cells <- region$cell[region$animal_id == a & region$bin_start == b0]
    inside <- terra::cellFromXY(fx$lake, cbind(ti$x, ti$y)) %in% cells
    expect_equal(sc$bins$n_truth[i], nrow(ti))
    expect_equal(sc$bins$error_mode_m[i], sqrt(sum((mode_xy - cen)^2)), tolerance = 1e-9)
    expect_equal(sc$bins$error_mean_m[i], sqrt(sum((mean_xy - cen)^2)), tolerance = 1e-9)
    expect_equal(sc$bins$in_region[i], mean(inside))
    expect_equal(sc$bins$area_m2[i], length(cells) * 900)
    expect_equal(sc$bins$region_mass[i],
                 sum(region$prob[region$animal_id == a & region$bin_start == b0]))
    expect_gte(sc$bins$region_mass[i], 0.9)
  }
})

test_that("the summaries are the medians, 90th percentiles and means of the bins", {
  sc <- score_surfaces(score_fixture$surfaces, score_fixture$truth)
  b <- sc$bins
  o <- sc$overall
  expect_equal(o$n_bins, nrow(b))
  expect_equal(o$n_animals, 3)
  expect_equal(o$median_error_mode_m, stats::median(b$error_mode_m))
  expect_equal(o$p90_error_mean_m, unname(stats::quantile(b$error_mean_m, 0.9)))
  expect_equal(o$hit_rate, mean(b$in_region))
  expect_equal(o$median_area_m2, stats::median(b$area_m2))
  expect_equal(o$mean_region_mass, mean(b$region_mass))
  a2 <- b[b$animal_id == "2", ]
  row <- sc$animals[sc$animals$animal_id == "2", ]
  expect_equal(row$n_bins, nrow(a2))
  expect_equal(row$p90_error_mode_m, unname(stats::quantile(a2$error_mode_m, 0.9)))
  expect_equal(row$hit_rate, mean(a2$in_region))
})

test_that("the hit-rate SE is clustered by animal and by bin start", {
  t0 <- as.POSIXct("2025-01-01", tz = "UTC")
  scored <- data.frame(animal_id = c("a", "a", "a", "b", "b", "c"),
                       bin_start = t0 + 3600 * c(0, 1, 2, 0, 1, 0),
                       in_region = c(1, 1, 0, 1, 0.5, 1), error_mode_m = 1:6,
                       error_mean_m = 1:6, area_m2 = 1:6, region_mass = 0.95)
  o <- overall_summary(scored)
  cbar <- 4.5 / 6
  e <- scored$in_region - cbar
  r_a <- c(sum(e[1:3]), sum(e[4:5]), e[6])
  r_t <- c(e[1] + e[4] + e[6], e[2] + e[5], e[3])
  v_a <- 3 / 2 * sum(r_a^2) / 36
  v_t <- 3 / 2 * sum(r_t^2) / 36
  v_b <- 6 / 5 * sum(e^2) / 36
  expect_equal(o$hit_rate_se_animal, sqrt(v_a))
  expect_equal(o$hit_rate_se, sqrt(max(v_a + v_t - v_b, v_a, v_t)))
  expect_equal(o$hit_rate_se_bins, stats::sd(scored$in_region) / sqrt(6))
  expect_identical(o$n_bin_starts, 3L)
  # One animal: the hour clusters alone. One animal and one hour: NA.
  one <- overall_summary(scored[1:3, ])
  e1 <- scored$in_region[1:3] - 2 / 3
  expect_equal(one$hit_rate_se, sqrt(3 / 2 * sum(e1^2) / 9))
  expect_identical(one$hit_rate_se_from, "bin start")
  expect_true(is.na(one$hit_rate_se_animal))
  expect_true(is.na(overall_summary(scored[1, ])$hit_rate_se))
  expect_identical(overall_summary(scored[0, ])$n_bins, 0L)
  expect_identical(nrow(animal_summaries(scored[0, ])), 0L)
})

test_that("the SE says when the one-way floor replaced the two-way value", {
  # Residuals that cancel within animals: V_a is 0, V_t is not, and the
  # two-way value falls below V_t, so the floor is the by-hour SE.
  e <- c(1, -1, 1, -1)
  animal <- c(1, 1, 2, 2)
  time <- c(1, 2, 1, 2)
  v_t <- 2 * sum(rowsum(e, time)^2) / 16
  got <- cluster_se(e, animal, time)
  expect_equal(got$se, sqrt(v_t))
  expect_identical(got$from, "bin start (floor)")
  expect_identical(cluster_se(e, time, animal)$from, "animal (floor)")
  t0 <- as.POSIXct("2025-01-01", tz = "UTC")
  scored <- data.frame(animal_id = rep(c("a", "b", "c"), each = 4),
                       bin_start = t0 + 3600 * rep(0:3, 3),
                       in_region = c(1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1),
                       error_mode_m = 1, error_mean_m = 1, area_m2 = 1, region_mass = 0.95)
  o <- overall_summary(scored)
  expect_identical(o$hit_rate_se_from, cluster_se(scored$in_region - o$hit_rate,
                                                  scored$animal_id,
                                                  as.numeric(scored$bin_start))$from)
  how <- c("animal and bin start" = "clustered two ways",
           "animal (floor)" = "clustered by animal, the one-way floor",
           "bin start (floor)" = "clustered by bin start, the one-way floor")
  # A wide console keeps each summary on one line for the match.
  out <- withr::with_options(list(width = 200), capture.output(print_score_lines(o, 0.95)))
  expect_match(paste(out, collapse = "\n"),
               sprintf("%s; 3 animals and 4 bin starts", how[[o$hit_rate_se_from]]),
               fixed = TRUE)
  o$hit_rate_se_from <- "animal"
  expect_output(print_score_lines(o, 0.95), "clustered by animal; 3 animals", fixed = TRUE)
  o$hit_rate_se <- NA_real_
  expect_output(print_score_lines(o, 0.95), "SE not estimable, one cluster each way",
                fixed = TRUE)
})

test_that("the 95% interval uses t with the smaller cluster count less one", {
  t0 <- as.POSIXct("2025-01-01", tz = "UTC")
  scored <- data.frame(animal_id = rep(c("a", "b", "c", "d"), each = 3),
                       bin_start = t0 + 3600 * rep(0:2, 4),
                       in_region = c(1, 0, 1, 1, 1, 0.5, 0, 1, 1, 1, 1, 1),
                       error_mode_m = 1, error_mean_m = 1, area_m2 = 1, region_mass = 0.95)
  o <- overall_summary(scored)
  expect_identical(o$hit_rate_df, 2L)
  half <- stats::qt(0.975, 2) * o$hit_rate_se
  expect_equal(o$hit_rate_lo, max(0, o$hit_rate - half))
  expect_equal(o$hit_rate_hi, min(1, o$hit_rate + half))
  expect_output(print_score_lines(o, 0.95), "t with 2 df", fixed = TRUE)
})

test_that("a hit rate of 0 or 1 gives no SE and an exact bound", {
  t0 <- as.POSIXct("2025-01-01", tz = "UTC")
  scored <- data.frame(animal_id = rep(c("a", "b", "c", "d", "e", "f"), each = 6),
                       bin_start = t0 + 3600 * rep(0:5, 6), in_region = 1,
                       error_mode_m = 1, error_mean_m = 1, area_m2 = 1, region_mass = 0.95)
  o <- overall_summary(scored)
  expect_true(is.na(o$hit_rate_se))
  expect_identical(o$hit_rate_se_from, "not estimable")
  expect_equal(c(o$hit_rate_lo, o$hit_rate_hi), c(0.025^(1 / 6), 1))
  out <- capture.output(print_score_lines(o, 0.95))
  expect_match(out, "SE not estimable, every bin scores the same", all = FALSE)
  expect_false(any(grepl("SE 0.000", out, fixed = TRUE)))
  scored$in_region <- 0
  o0 <- overall_summary(scored)
  expect_equal(c(o0$hit_rate_lo, o0$hit_rate_hi), c(0, 1 - 0.025^(1 / 6)))
  scored$in_region <- 0.5
  expect_true(is.na(overall_summary(scored)$hit_rate_lo))
})

test_that("a larger level gives larger regions", {
  a50 <- score_surfaces(score_fixture$surfaces, score_fixture$truth, level = 0.5)
  a99 <- score_surfaces(score_fixture$surfaces, score_fixture$truth, level = 0.99)
  expect_true(all(a99$bins$area_m2 >= a50$bins$area_m2))
  expect_gte(a99$overall$hit_rate, a50$overall$hit_rate)
  expect_identical(attr(a99, "level"), 0.99)
})

test_that("bins without truth are counted and left out", {
  truth <- score_fixture$truth
  first <- min(attr(score_fixture$surfaces, "bins")$bin_start)
  drop <- truth$animal_id == "1" & truth$datetime < first + 3600
  expect_warning(sc <- score_surfaces(score_fixture$surfaces, truth[!drop, ]),
                 "1 of 9 bins hold no transmission")
  expect_identical(attr(sc, "no_truth"), 1L)
  expect_equal(nrow(sc$bins), 8)
  expect_output(print(sc), "1 bins had no truth")
  # No truth in any bin: every bin is left out, and the grid check has
  # nothing to test, so it does not stop.
  later <- truth
  later$datetime <- later$datetime + 10 * 86400
  expect_warning(none <- score_surfaces(score_fixture$surfaces, later), "9 of 9 bins")
  expect_identical(none$overall$n_bins, 0L)
  expect_output(print(none), "0 animal-bins of 0 animals")
})

test_that("print shows the summaries", {
  sc <- score_surfaces(score_fixture$surfaces, score_fixture$truth)
  # A wide console keeps each summary on one line for the matches.
  out <- withr::with_options(list(width = 200), capture.output(print(sc)))
  expect_match(out[1], "9 animal-bins of 3 animals; 95% credible regions")
  expect_match(out[4], sprintf("hit rate %.3f: share of the truth inside the 95%% regions",
                               sc$overall$hit_rate), fixed = TRUE)
  expect_match(out[4], "regions (a calibration check)", fixed = TRUE)
  expect_match(out[5], "95% interval", fixed = TRUE)
  expect_match(out[6], "median region")
  expect_false(any(grepl("animal-hours with transmissions", out)))
})

test_that("print wraps to the console width", {
  sc <- score_surfaces(score_fixture$surfaces, score_fixture$truth)
  out <- withr::with_options(list(width = 60), capture.output(print(sc)))
  expect_lte(max(nchar(out)), 60)
  expect_match(paste(trimws(out), collapse = " "),
               "regions (a calibration check)", fixed = TRUE)
})

# ---- score_surfaces(): refusals ----

test_that("bad inputs are refused", {
  fx <- score_fixture
  expect_error(score_surfaces(data.frame(), fx$truth), "location_surface")
  expect_error(score_surfaces(fx$surfaces, fx$truth[c("x", "y")]), "columns animal_id")
  bad <- fx$truth
  bad$datetime <- as.POSIXct(format(bad$datetime), tz = "America/Toronto")
  expect_error(score_surfaces(fx$surfaces, bad), "UTC")
  bad <- fx$truth
  bad$x[1] <- NA
  expect_error(score_surfaces(fx$surfaces, bad), "no NA")
  expect_error(score_surfaces(fx$surfaces, fx$truth, level = 1), "kept mass")
  bad <- fx$truth
  bad$datetime[1] <- NA
  expect_error(score_surfaces(fx$surfaces, bad), "`truth\\$datetime` must have no NA; 1 of")
})

test_that("truth off the surfaces' grid is refused or counted", {
  fx <- score_fixture
  far <- fx$truth
  far$x <- far$x + 1e6
  expect_error(score_surfaces(fx$surfaces, far), "None of the .* Is `truth` in the CRS")
  lonlat <- fx$truth
  lonlat$x <- -78.1
  lonlat$y <- 44.5
  expect_error(score_surfaces(fx$surfaces, lonlat), "None of the")
  part <- fx$truth
  moved <- part$animal_id == "1"
  part$x[moved] <- part$x[moved] + 1e6
  joined <- join_truth(check_truth(part), attr(fx$surfaces, "bins"))$truth
  expect_warning(score_surfaces(fx$surfaces, part),
                 sprintf("^%d of %d transmissions .* lie outside the grid and 0 on land",
                         sum(joined$animal_id == "1"), nrow(joined)))
  land_xy <- terra::xyFromCell(fx$lake, which(is.na(terra::values(fx$lake)))[1])
  onland <- fx$truth
  two <- onland$animal_id == "2"
  onland$x[two] <- land_xy[1, 1]
  onland$y[two] <- land_xy[1, 2]
  joined <- join_truth(check_truth(onland), attr(fx$surfaces, "bins"))$truth
  expect_warning(score_surfaces(fx$surfaces, onland),
                 sprintf("^0 of %d transmissions .* outside the grid and %d on land",
                         nrow(joined), sum(joined$animal_id == "2")))
  # Simulated truth is always on water: no warning.
  expect_silent(score_surfaces(fx$surfaces, fx$truth))
})
