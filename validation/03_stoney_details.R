# ---- Array scoring on the Stoney example: numbers for the help pages ----
# Computes on the shipped Stoney Lake raster every number in the details of
# ?evaluate_array: the scored share of a reference run, the hit rate under
# departures the likelihood does not model, the hit rate against array
# density, how well the hit-rate interval holds, and the time a call takes.
# The setting is that of vignette("arraydesign"): the 100 m Stoney raster, 20
# receivers on the grid over water at least 1 m deep (seed 1), the
# depth-interpolated logistic curve, crw_model(30, 15, 40, time_step_s = 60)
# and a delay of 120 to 180 s. Writes validation/outputs/03_*.csv.
# Run from the package root, with the package and detfield installed:
#   Rscript validation/03_stoney_details.R

# ---- Libraries ----
library(arraydesign)
library(detfield)
out <- "validation/outputs"
cat("arraydesign", as.character(packageVersion("arraydesign")), "; detfield",
    as.character(packageVersion("detfield")), ";", R.version.string, "\n")

# ---- Setting ----
stoney <- terra::rast(system.file("extdata", "stoney_depth.tif", package = "arraydesign"))
deployable <- terra::ifel(abs(stoney) < 1, NA, stoney)
rx <- place_receivers(deployable, n = 20, method = "grid", seed = 1)
curve <- de_logistic(d50 = c(250, 750), d05 = c(500, 1500), depth = c(1, 35))
with_sigma <- function(s) de_model(curve, sigma = list(station_hour = s))
delay <- c(120, 180)
still <- crw_model(1e-6, 1e-6, 40)
cat("water cells:", sum(!is.na(terra::values(stoney))), "; receivers:", nrow(rx), "\n")

# The seed is an argument of evaluate_array(); the runs below list it with
# their settings, so it is taken out of the list here.
run <- function(model, n_animals, control, layout = rx) {
  seed <- control$seed
  control$seed <- NULL
  suppressMessages(evaluate_array(layout, stoney, model, n_animals = n_animals, seed = seed,
                                  control = c(list(tag_delay = delay), control)))
}
row_of <- function(label, ev) {
  o <- ev$score$overall
  data.frame(case = label, n_animals = o$n_animals, n_animal_hours = o$n_animal_hours,
             n_bins = o$n_bins, scored_share = o$scored_share, hit_rate = o$hit_rate,
             hit_rate_se = o$hit_rate_se, median_error_mean_m = o$median_error_mean_m)
}

# ---- 1. Scored share of the reference run ----
# The app's reference run: 20 animals x 6 hours, seed 1.
ref <- row_of("reference run", run(curve, 20, list(hours = 6, seed = 1)))
print(ref, digits = 4)
write.csv(ref, file.path(out, "03_reference.csv"), row.names = FALSE)

# ---- 2. Departures from the likelihood's model ----
# 200 animals x 24 hours, binomial: still animals, a walk of 6 m a minute and
# the default walk of 30 m a minute, all with the 120 to 180 s delay.
dep <- rbind(
  row_of("still", run(curve, 200, list(movement = still, hours = 24, seed = 201))),
  row_of("6 m/min", run(curve, 200, list(movement = crw_model(6, 4, 40), hours = 24,
                                         seed = 203))),
  row_of("default walk", run(curve, 200, list(hours = 24, seed = 202)))
)
print(dep, digits = 4)
write.csv(dep, file.path(out, "03_departures.csv"), row.names = FALSE)

# ---- 3. Hit rate against array density ----
# Grids of 1, 5, 10 and 20 receivers on the deployable water, 60 animals x 12
# hours, still and the default walk, one seed.
dens <- do.call(rbind, lapply(c(1, 5, 10, 20), function(n_rx) {
  lay <- place_receivers(deployable, n = n_rx, method = "grid", seed = 1)
  rbind(cbind(receivers = n_rx, row_of("still", run(curve, 60, list(movement = still,
                                                     hours = 12, seed = 11), lay))),
        cbind(receivers = n_rx, row_of("default walk", run(curve, 60, list(hours = 12,
                                                            seed = 11), lay))))
}))
print(dens, digits = 4)
write.csv(dens, file.path(out, "03_density.csv"), row.names = FALSE)

# ---- 4. The hit-rate interval at the defaults ----
# 40 seeds at 20 animals x 6 hours (binomial, then station-hour sigma 1) and
# 20 seeds at sigma 1 x 24 hours. A run misses when its 95% interval leaves
# out the mean of its replicates; nominal 5%.
se_rep <- function(model, seeds, hours) {
  rows <- parallel::mclapply(seeds, function(sd) {
    o <- run(model, 20, list(hours = hours, seed = sd))$score$overall
    data.frame(seed = sd, hit_rate = o$hit_rate, se = o$hit_rate_se, lo = o$hit_rate_lo,
               hi = o$hit_rate_hi, n_bins = o$n_bins)
  }, mc.cores = 8)
  do.call(rbind, rows)
}
cases <- list(list("binomial, 6 h", curve, 1001:1040, 6),
              list("sigma 1, 6 h", with_sigma(1), 1001:1040, 6),
              list("sigma 1, 24 h", with_sigma(1), 2001:2020, 24))
se_summary <- do.call(rbind, lapply(cases, function(cs) {
  r <- se_rep(cs[[2]], cs[[3]], cs[[4]])
  m <- mean(r$hit_rate)
  data.frame(case = cs[[1]], runs = nrow(r), replicate_mean = m,
             replicate_sd = stats::sd(r$hit_rate), mean_se = mean(r$se, na.rm = TRUE),
             miss = sum(r$lo > m | r$hi < m, na.rm = TRUE))
}))
print(se_summary, digits = 4)
write.csv(se_summary, file.path(out, "03_se_summary.csv"), row.names = FALSE)

# ---- 5. Time at the defaults ----
# The reference call, 20 animals x 6 hours, ten runs per model, one R process.
host <- paste(Sys.info()[["sysname"]], Sys.info()[["machine"]], parallel::detectCores(),
              "cores")
timing <- do.call(rbind, lapply(c(NA, 2.55), function(s) {
  m <- if (is.na(s)) curve else with_sigma(s)
  secs <- vapply(1:10, function(i) {
    t0 <- proc.time()[["elapsed"]]
    run(m, 20, list(hours = 6, seed = 400 + i))
    proc.time()[["elapsed"]] - t0
  }, numeric(1))
  data.frame(station_hour = s, runs = 10, median_s = stats::median(secs), max_s = max(secs))
}))
cat("host:", host, "\n")
print(timing)
write.csv(timing, file.path(out, "03_timing.csv"), row.names = FALSE)

sessionInfo()
