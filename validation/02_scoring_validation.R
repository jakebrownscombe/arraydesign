# ---- Array scoring: validation runs ----
# Runs the checks of the scoring module at a larger size than the tests and
# prints every number with its denominator: calibration of the 95% regions
# (the hit rate) when one model simulates and positions, what movement within
# the hour and a random tag delay do to it, how the movement shortfall grows
# with array density, how well the hit-rate interval holds at the default six
# hours, the error under a mis-scaled d50 and a too-small sigma, and the time
# evaluate_array() takes at its defaults. Writes validation/outputs/02_*.
# Run from the package root, with the package and detfield installed:
#   Rscript validation/02_scoring_validation.R

# ---- Libraries ----
library(arraydesign)
library(detfield)
out <- "validation/outputs"
cat("arraydesign", as.character(packageVersion("arraydesign")), "; detfield",
    as.character(packageVersion("detfield")), ";", R.version.string, "\n")

# ---- Setting ----
# example_lake(): 2,152 water cells of 30 m. Eight receivers on a lattice. A
# logistic range curve with d50 250 m and d05 600 m, read at the distance
# through water.
lake <- example_lake()
rx <- place_receivers(lake, n = 8)
logistic <- de_logistic(250, 600)
with_sigma <- function(s) de_model(logistic, sigma = list(station_hour = s))
still <- crw_model(1e-6, 1e-6, 40)
exact <- list(movement = still, tag_delay = 120)
cat("water cells:", sum(!is.na(terra::values(lake))), "; receivers:", nrow(rx), "\n")

# The seed is an argument of evaluate_array(); the runs below list it with
# their settings, so it is taken out of the list here.
run <- function(model, n_animals, control) {
  seed <- control$seed
  control$seed <- NULL
  suppressMessages(evaluate_array(rx, lake, model, n_animals = n_animals, seed = seed,
                                  control = control))
}
row_of <- function(label, ev) {
  o <- ev$score$overall
  data.frame(case = label, n_animals = o$n_animals, n_bin_starts = o$n_bin_starts,
             n_bins = o$n_bins, n_animal_hours = o$n_animal_hours, hit_rate = o$hit_rate,
             hit_rate_se = o$hit_rate_se, hit_rate_se_from = o$hit_rate_se_from,
             hit_rate_lo = o$hit_rate_lo, hit_rate_hi = o$hit_rate_hi,
             hit_rate_se_animal = o$hit_rate_se_animal,
             z = (o$hit_rate - 0.95) / o$hit_rate_se, mean_region_mass = o$mean_region_mass,
             median_error_mean_m = o$median_error_mean_m,
             median_error_mode_m = o$median_error_mode_m,
             p90_error_mean_m = o$p90_error_mean_m, median_area_m2 = o$median_area_m2)
}

# ---- 1. Calibration under the likelihood's own model ----
# Animals still at a cell centre, a transmission every 120 s, so n is exactly
# 30 an hour: the data come from the model location_surface() computes.
cal <- rbind(
  row_of("binomial, exact", run(logistic, 200, c(exact, list(hours = 24, seed = 101)))),
  row_of("station_hour 1, exact", run(with_sigma(1), 60, c(exact, list(hours = 24,
                                                                       seed = 102)))),
  row_of("station_hour 2.55, exact", run(with_sigma(2.55), 60, c(exact, list(hours = 24,
                                                                             seed = 103))))
)
print(cal, digits = 4)
write.csv(cal, file.path(out, "02_calibration.csv"), row.names = FALSE)

# ---- 2. Departures from the likelihood's model ----
# Random delay (60 to 180 s plus a 5 s burst, the default) with still
# animals; then the default movement (a random walk of 30 m a minute) and a slower
# walk of 6 m a minute, both with the random delay.
dep <- rbind(
  row_of("binomial, still, random delay",
         run(logistic, 200, list(movement = still, hours = 24, seed = 201))),
  row_of("binomial, default walk, random delay", run(logistic, 200, list(hours = 24,
                                                                          seed = 202))),
  row_of("binomial, 6 m/min, random delay",
         run(logistic, 200, list(movement = crw_model(6, 4, 40), hours = 24, seed = 203))),
  row_of("station_hour 1, default walk, random delay",
         run(with_sigma(1), 60, list(hours = 24, seed = 204)))
)
print(dep, digits = 4)
write.csv(dep, file.path(out, "02_departures.csv"), row.names = FALSE)

# ---- 2b. Movement shortfall against array density ----
# Still animals and the default movement at 1, 4, 8 and 16 receivers, 60
# animals x 12 hours, binomial, one seed. The shortfall is still minus moving.
dens <- do.call(rbind, lapply(c(1, 4, 8, 16), function(n_rx) {
  rx_n <- place_receivers(lake, n = n_rx)
  one <- function(label, ctrl) {
    ev <- suppressMessages(evaluate_array(rx_n, lake, logistic, n_animals = 60, seed = 11,
                                          control = ctrl))
    cbind(receivers = n_rx, row_of(label, ev))
  }
  rbind(one("still", list(movement = still, hours = 12)),
        one("default movement", list(hours = 12)))
}))
print(dens, digits = 4)
write.csv(dens, file.path(out, "02_density.csv"), row.names = FALSE)
shortfall <- with(dens, hit_rate[case == "still"] - hit_rate[case == "default movement"])
names(shortfall) <- unique(dens$receivers)
print(round(shortfall, 4))

# ---- 2c. The hit-rate interval at the defaults ----
# Replicate runs at the defaults (20 animals x 6 hours, so 6 bin starts), 40
# seeds each, binomial and station-hour sigma 1; then sigma 1 at 24 hours, 20
# seeds. The replicate mean stands in for the true hit rate. A run "misses"
# when its 95% interval (t with min(animals, hours) - 1 df) leaves out the
# replicate mean; nominal is 5%. Runs in parallel: each seed fixes its run.
se_rep <- function(model, seeds, ctrl = list()) {
  rows <- parallel::mclapply(seeds, function(sd) {
    o <- suppressMessages(evaluate_array(rx, lake, model, seed = sd, control = ctrl))
    o <- o$score$overall
    data.frame(seed = sd, hit_rate = o$hit_rate, se = o$hit_rate_se, lo = o$hit_rate_lo,
               hi = o$hit_rate_hi, df = o$hit_rate_df, n_bins = o$n_bins)
  }, mc.cores = 8)
  do.call(rbind, rows)
}
se_cases <- list(list("binomial, 6 h", logistic, 1001:1040, list()),
                 list("sigma 1, 6 h", with_sigma(1), 1001:1040, list()),
                 list("sigma 1, 24 h", with_sigma(1), 2001:2020, list(hours = 24)))
se_runs <- do.call(rbind, lapply(se_cases, function(cs) {
  cbind(case = cs[[1]], se_rep(cs[[2]], cs[[3]], cs[[4]]))
}))
write.csv(se_runs, file.path(out, "02_se_replicates.csv"), row.names = FALSE)
se_summary <- do.call(rbind, lapply(split(se_runs, factor(se_runs$case, unique(se_runs$case))),
                                    function(r) {
  m <- mean(r$hit_rate)
  data.frame(case = r$case[1], runs = nrow(r), bins_per_run = stats::median(r$n_bins),
             df = stats::median(r$df), replicate_mean = m, replicate_sd = stats::sd(r$hit_rate),
             mean_se = mean(r$se, na.rm = TRUE), se_na = sum(is.na(r$se)),
             miss_t = sum(r$lo > m | r$hi < m, na.rm = TRUE),
             miss_196 = sum(abs(r$hit_rate - m) > 1.96 * r$se, na.rm = TRUE))
}))
rownames(se_summary) <- NULL
print(se_summary, digits = 4)
write.csv(se_summary, file.path(out, "02_se_summary.csv"), row.names = FALSE)

# ---- 3. Misspecified d50 and sigma ----
# One simulation per sweep: the seed fixes tracks, transmissions and
# detections, so only the positioning model changes between rows.
scales <- c(0.5, 0.71, 1, 1.41, 2)
sweep_ev <- lapply(scales, function(s) {
  run(logistic, 100, c(exact, list(hours = 24, seed = 301,
                                   position_model = perturb_de_model(logistic, s))))
})
sweep_d50 <- do.call(rbind, Map(function(s, ev) {
  cbind(d50_scale = s, row_of(sprintf("d50 x %.2f", s), ev))
}, scales, sweep_ev))
# Same detections in every row: the transmissions heard must not change.
stopifnot(length(unique(lapply(sweep_ev, function(ev) ev$detections$detections))) == 1)
print(sweep_d50, digits = 4)
write.csv(sweep_d50, file.path(out, "02_d50_sweep.csv"), row.names = FALSE)
# detfield's own simulator gave coverage 0.6745 when positioning at sigma
# 1.362 against a true 2.55; here the walker and simulate_detections() make
# the data.
sweep_sigma <- do.call(rbind, lapply(c(1.362 / 2.55, 1), function(s) {
  ev <- run(with_sigma(2.55), 60, c(exact, list(
    hours = 24, seed = 302, position_model = perturb_de_model(with_sigma(2.55), sigma_scale = s)
  )))
  cbind(sigma_position = 2.55 * s, row_of(sprintf("sigma %.3f", 2.55 * s), ev))
}))
print(sweep_sigma, digits = 4)
write.csv(sweep_sigma, file.path(out, "02_sigma_sweep.csv"), row.names = FALSE)

# ---- 4. Time at the defaults ----
# 20 animals x 6 hours on example_lake(), ten runs per model, one R process
# (cores = 1), run one after another. The host is printed, because other work
# on the machine slows the runs.
host <- paste(Sys.info()[["sysname"]], Sys.info()[["machine"]], parallel::detectCores(),
              "cores")
timing <- do.call(rbind, lapply(c(NA, 1, 2.55), function(s) {
  m <- if (is.na(s)) logistic else with_sigma(s)
  secs <- vapply(1:10, function(i) {
    t0 <- proc.time()[["elapsed"]]
    ev <- run(m, 20, list(seed = 400 + i))
    proc.time()[["elapsed"]] - t0
  }, numeric(1))
  data.frame(station_hour = s, runs = 10, median_s = stats::median(secs), max_s = max(secs),
             over_10_s = sum(secs >= 10))
}))
cat("host:", host, "\n")
print(timing)
write.csv(timing, file.path(out, "02_timing.csv"), row.names = FALSE)

# ---- 5. Figures ----
# (a) Map of the error at the defaults with Stoney's sigma.
ev_map <- run(with_sigma(2.55), 20, list(seed = 501))
fig_map <- file.path(out, "02_error_map.png")
png(fig_map, width = 7, height = 5, units = "in", res = 300)
plot(ev_map)
dev.off()
cat("map:", ev_map$score$overall$n_bins, "animal-hours of", ev_map$score$overall$n_animals,
    "animals\n")

# (b) Median error of the surface mean, with a 95% interval from 2,000
# bootstrap resamples of animals, and the hit rate with its 95% interval,
# against the d50 scale of the positioning model.
boot_median <- function(b, reps = 2000) {
  by <- split(b$error_mean_m, b$animal_id)
  withr::with_seed(502, stats::quantile(replicate(reps, {
    stats::median(unlist(by[sample(length(by), replace = TRUE)]))
  }), c(0.025, 0.975)))
}
ci <- t(vapply(sweep_ev, function(ev) boot_median(ev$score$bins), numeric(2)))
sweep_d50$median_error_mean_lo <- ci[, 1]
sweep_d50$median_error_mean_hi <- ci[, 2]
write.csv(sweep_d50, file.path(out, "02_d50_sweep.csv"), row.names = FALSE)
fig_sweep <- file.path(out, "02_d50_sweep.png")
png(fig_sweep, width = 7, height = 3.4, units = "in", res = 300)
op <- par(mfrow = c(1, 2), mar = c(4.2, 4.4, 1.2, 0.8), cex = 0.9)
x <- log2(sweep_d50$d50_scale)
plot(x, sweep_d50$median_error_mean_m, pch = 19, xaxt = "n", ylim = c(0, max(ci[, 2]) * 1.05),
     xlab = "d50 of the positioning model / true d50", ylab = "Median error of mean (m)")
axis(1, at = x, labels = format(sweep_d50$d50_scale))
arrows(x, ci[, 1], x, ci[, 2], angle = 90, code = 3, length = 0.04)
mtext("(a)", side = 3, adj = 0, line = 0.1)
plot(x, sweep_d50$hit_rate, pch = 19, xaxt = "n", ylim = c(0, 1),
     xlab = "d50 of the positioning model / true d50", ylab = "Hit rate of 95% regions")
axis(1, at = x, labels = format(sweep_d50$d50_scale))
ok <- sweep_d50$hit_rate_hi > sweep_d50$hit_rate_lo
arrows(x[ok], sweep_d50$hit_rate_lo[ok], x[ok], sweep_d50$hit_rate_hi[ok], angle = 90, code = 3,
       length = 0.04)
abline(h = 0.95, lty = 2)
mtext("(b)", side = 3, adj = 0, line = 0.1)
par(op)
dev.off()

sessionInfo()
