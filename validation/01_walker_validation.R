# ---- Walker and detection simulator: validation runs ----
# Runs the checks of the walker and the detection simulator at a larger size
# than the tests and prints every number with its denominator. Writes
# validation/outputs/. Run from the package root, with the package installed:
#   STONEY_30M=<.rds file of a 30 m depth raster of Stoney Lake> \
#     GLATOS_LIB=<library with glatos 0.9.8> Rscript validation/01_walker_validation.R
# Without STONEY_30M the Stoney section is skipped. The occupancy table
# then prints, and the table and its figure (08_shore_occupancy.csv and .png)
# go to tempdir(), because the files in validation/outputs/ hold the two
# Stoney walks too.
# Without glatos 0.9.8 the oracle section is skipped.

# ---- Libraries ----
if (nzchar(Sys.getenv("GLATOS_LIB"))) .libPaths(c(Sys.getenv("GLATOS_LIB"), .libPaths()))
library(arraydesign)
library(terra)
library(sf)
# The Stoney raster is a raster package layer. raster fails to load its Rcpp
# module once this script has filled the global environment, so it loads here.
if (nzchar(Sys.getenv("STONEY_30M"))) invisible(loadNamespace("raster"))
# The fixtures read their files through testthat's test_path(). This script
# runs without testthat attached, so the same path is built here.
test_path <- function(...) file.path("tests", "testthat", ...)
source("tests/testthat/helper-fixtures.R")
out <- "validation/outputs"
cat("arraydesign", as.character(packageVersion("arraydesign")), ";", R.version.string, "\n")

# ---- 1. Land wall: 20 paths x 1,000 steps per mean step ----
# Positive control: a walker with the same samplers that tests end points
# only, as positionR 1.x did, and holds when an end point is on land. It must
# show crossings, or the oracle cannot see them.
wall <- wall_lake()
land <- land_polygons(wall)
wall_x <- 700000 + 49 * 30
start <- data.frame(x = 700000 + 25 * 30, y = 4900000 + 30 * 30)
end_only_walk <- function(grid, model, n_animals, n_steps, start, seed) {
  wg <- arraydesign:::water_grid(grid)
  walks <- withr::with_seed(seed, lapply(seq_len(n_animals), function(i) {
    len <- model$step(n_steps)
    turn <- model$turn(n_steps)
    x <- y <- numeric(n_steps + 1)
    x[1] <- start$x
    y[1] <- start$y
    b <- runif(1, -pi, pi)
    for (k in seq_len(n_steps)) {
      x1 <- x[k] + len[k] * cos(b + turn[k])
      y1 <- y[k] + len[k] * sin(b + turn[k])
      ok <- arraydesign:::is_wet(wg, x1, y1)
      x[k + 1] <- if (ok) x1 else x[k]
      y[k + 1] <- if (ok) y1 else y[k]
      if (ok) b <- b + turn[k]
    }
    data.frame(path_id = i, step = seq_len(n_steps + 1) - 1L, x = x, y = y)
  }))
  structure(do.call(rbind, walks), class = c("sim_tracks", "data.frame"),
            grid = terra::wrap(grid))
}
wall_run <- function(mean_step, test) {
  m <- crw_model(mean_step, mean_step / 2, 45)
  tr <- if (test == "whole step") {
    simulate_tracks(wall, m, n_animals = 20, n_steps = 1000, control = track_control(start = start),
                    seed = mean_step)
  } else {
    end_only_walk(wall, m, 20, 1000, start, mean_step)
  }
  segs <- step_segments(tr)
  list(tracks = tr, steps = 20000L, moving = length(segs),
       on_land = sum(is.na(terra::extract(wall, cbind(tr$x, tr$y))[, 1])),
       crossing = sum(lengths(sf::st_intersects(segs, land)) > 0),
       paths_across = length(unique(tr$path_id[tr$x > wall_x])))
}
wall_tab <- do.call(rbind, lapply(c(30, 60, 90), function(s) {
  do.call(rbind, lapply(c("whole step", "end point only"), function(m) {
    r <- wall_run(s, m)
    data.frame(mean_step_m = s, land_test = m, steps = r$steps, moving_segments = r$moving,
               positions_on_land = r$on_land, segments_crossing_land = r$crossing,
               paths_in_far_basin = r$paths_across, of_paths = 20L)
  }))
}))
print(wall_tab)
write.csv(wall_tab, file.path(out, "01_wall_test.csv"), row.names = FALSE)
## Map of the 60 m runs, both land tests; the land wall is at x = 701,470 to 701,500
## Counts in the titles are over all 20 paths x 1,000 steps.
# plot() of tracks returns a ggplot, so each panel is printed into its own
# viewport of one device.
png(file.path(out, "01_wall_test_map.png"), width = 1600, height = 1900, res = 200)
grid::grid.newpage()
grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 1)))
for (i in 1:2) {
  m <- c("whole step", "end point only")[i]
  r <- wall_run(60, m)
  # The first 300 steps of 4 paths, in the colours that read on dark water;
  # all 1,000 steps of 20 paths hide the scale bar.
  shown <- r$tracks[r$tracks$step <= 300, ]
  attr(shown, "grid") <- attr(r$tracks, "grid")
  p <- plot(shown, paths = 1:4) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::ggtitle(sprintf("(%s) %s: %d of %d steps cross", letters[i],
                             c("Whole-step test", "End points only")[i], r$crossing, r$steps))
  print(p, vp = grid::viewport(layout.pos.row = i, layout.pos.col = 1))
}
dev.off()

# ---- 2. Steps and turns: offshore against the model, near shore as reflected ----
lake <- open_lake()
m <- crw_model(30, 10, 30)
tr <- simulate_tracks(lake, m, n_animals = 20, n_steps = 1000,
                      control = track_control(start = data.frame(x = 703030, y = 4903030)),
                      seed = 11)
shore <- pmin(tr$x - 700030, 706030 - tr$x, tr$y - 4900030, 4906030 - tr$y)
n <- nrow(tr)
same <- tr$path_id[-1] == tr$path_id[-n]
len <- sqrt(diff(tr$x)^2 + diff(tr$y)^2)
bearing <- atan2(diff(tr$y), diff(tr$x))
away <- same & shore[-n] > 200
k <- which(away[-1] & away[-length(away)] & same[-1])
ks_step <- suppressWarnings(ks.test(len[away], withr::with_seed(2, m$step(1e5))))
ks_turn <- suppressWarnings(ks.test(wrap_angle(bearing[k + 1] - bearing[k]),
                                    withr::with_seed(3, m$turn(1e5))))
cat(sprintf("KS steps: D = %.4f, p = %.3f, %d steps of %d starting > 200 m from shore\n",
            ks_step$statistic, ks_step$p.value, sum(away), sum(same)))
cat(sprintf("KS turns: D = %.4f, p = %.3f, %d turns\n", ks_turn$statistic, ks_turn$p.value,
            length(k)))
## Turns by distance from shore on example_lake(), where the shore is near
toy <- example_lake()
twg <- arraydesign:::water_grid(toy)
tland <- rast(toy)
values(tland) <- ifelse(twg$wet, NA, 1)
tdist <- values(distance(tland), mat = FALSE)
tt <- simulate_tracks(toy, m, n_animals = 20, n_steps = 20000, seed = 12)
tn <- nrow(tt)
tsame <- tt$path_id[-1] == tt$path_id[-tn]
tb <- atan2(diff(tt$y), diff(tt$x))
tmoved <- tsame & (diff(tt$x) != 0 | diff(tt$y) != 0)
kk <- which(tmoved[-1] & tmoved[-length(tmoved)] & tsame[-1])
turn <- abs(wrap_angle(tb[kk + 1] - tb[kk]))
d_at <- tdist[tt$cell[kk + 1]]
model_wide <- mean(abs(withr::with_seed(4, m$turn(1e6))) > pi / 2)
turn_tab <- data.frame(
  where = c("all", "within 60 m of land", "more than 200 m from land"),
  turns = c(length(turn), sum(d_at <= 60), sum(d_at > 200)),
  share_over_90 = c(mean(turn > pi / 2), mean(turn[d_at <= 60] > pi / 2),
                    mean(turn[d_at > 200] > pi / 2)))
print(turn_tab)
cat(sprintf("model: %.2f%% of 1e6 turn draws exceed 90 degrees\n", 100 * model_wide))

# ---- 2b. Time spent near the shore against its share of the water ----
# Per path, the share of positions in each band of distance from land, after
# a 2,000-step burn-in, against the band's share of the water cells. The
# standard error is across the 28 paths.
occupancy <- function(grid, model, dist, label, n_animals = 28, n_steps = 20000, burn = 2000) {
  wet <- !is.na(values(grid, mat = FALSE))
  bands <- cut(dist, c(0, 30 * sqrt(2) + 1e-6, 90, 150, Inf), right = TRUE,
               labels = c("<= 42 m", "42-90 m", "90-150 m", "> 150 m"))
  water_share <- as.vector(table(bands[wet])) / sum(wet)
  per_path <- do.call(rbind, parallel::mclapply(seq_len(n_animals), function(i) {
    tr <- simulate_tracks(grid, model, n_steps = n_steps + burn, seed = 9000 + i)
    keep <- tr$step > burn
    c(as.vector(table(bands[tr$cell[keep]])) / sum(keep),
      held = sum(attr(tr, "escapes")$held), reflected = sum(attr(tr, "escapes")$reflected))
  }, mc.cores = min(28L, parallel::detectCores())))
  share <- colMeans(per_path[, 1:4])
  se <- apply(per_path[, 1:4], 2, sd) / sqrt(n_animals)
  data.frame(walk = label, band = levels(bands), water_cells = as.vector(table(bands[wet])),
             water_share = water_share, position_share = share, se = se,
             ratio = share / water_share, z = (share - water_share) / se,
             positions = n_animals * n_steps, reflected = sum(per_path[, "reflected"]),
             held = sum(per_path[, "held"]), steps = n_animals * (n_steps + burn))
}
occ <- rbind(occupancy(toy, crw_model(30, 10, Inf), tdist, "example_lake, crw(30, 10, Inf)"),
             occupancy(toy, crw_model(30, 10, 30), tdist, "example_lake, crw(30, 10, 30)"))
# ---- 3. Move and rest dwell times ----
p <- matrix(c(0.95, 0.10, 0.05, 0.90), 2)
tr2 <- simulate_tracks(open_lake(), crw_model(30, 10, 30, transition = p), n_animals = 1,
                       n_steps = 20000, seed = 5)
runs <- rle(tr2$state[-1])
dwell <- tapply(runs$lengths, runs$values, mean)
n_runs <- tapply(runs$lengths, runs$values, length)
cat(sprintf("dwell: move %.2f steps (%d runs; configured 20), rest %.2f (%d runs; configured 10)\n",
            dwell[1], n_runs[1], dwell[2], n_runs[2]))

# ---- 4. Transmissions per hour against the tag delay ----
tr3 <- simulate_tracks(example_lake(), crw_model(30, 10, 30, time_step_s = 600), n_animals = 10,
                       n_steps = 600, seed = 1)
hours <- 10 * 600 * 600 / 3600
tx_tab <- do.call(rbind, lapply(list(c(60, 180), c(120, 240), c(30, 90)), function(d) {
  # tag_delay is start-to-start: the programmed delay plus the 5 s transmission.
  tx <- simulate_transmissions(tr3, tag_delay = d + 5, seed = 2)
  data.frame(delay = paste(d, collapse = " to "), transmissions = nrow(tx), track_hours = hours,
             per_hour = nrow(tx) / hours, expected = 3600 / (mean(d) + 5))
}))
print(tx_tab)

# ---- 5. Detection rate per distance band against the known curve ----
lake <- example_lake()
rx <- place_receivers(lake, n = 6)
field <- make_field(lake, rx, a = 2.5, s = 100)
tr4 <- simulate_tracks(lake, crw_model(30, 10, 30), n_animals = 20, n_steps = 2000, seed = 21)
tx4 <- simulate_transmissions(tr4, seed = 22)
dist <- matrix(field$distances$distance_m, ncol = 6)[match(tx4$cell, field$cells), ]
pr <- plogis(2.5 - dist / 100)
band <- cut(dist, c(seq(0, 900, 100), Inf), right = FALSE)
cover <- sapply(1:200, function(seed) {
  det <- simulate_detections(tx4, field, seed = seed)
  heard <- matrix(FALSE, nrow(tx4), 6)
  heard[cbind(det$detections$transmission_id,
              match(det$detections$receiver_serial, rx$receiver_serial))] <- TRUE
  vapply(levels(band), function(b) {
    ci <- binom.test(sum(heard[band == b]), sum(band == b))$conf.int
    mean(pr[band == b]) >= ci[1] && mean(pr[band == b]) <= ci[2]
  }, logical(1))
})
band_tab <- data.frame(band = levels(band), pairs = as.vector(table(band)),
                       mean_p = as.vector(tapply(pr, band, mean)),
                       covered_of_200 = rowSums(cover))
print(band_tab)
cat(sprintf("binomial 95%% interval held the true rate in %d of %d band-seeds (%.1f%%)\n",
            sum(cover), length(cover), 100 * mean(cover)))

# ---- 6. Oracle: glatos 0.9.8 on the same lake, tag and curve ----
if (requireNamespace("glatos", quietly = TRUE) && packageVersion("glatos") >= "0.9.8") {
  poly <- st_cast(st_as_sf(as.polygons(ifel(is.na(lake), NA, 1), dissolve = TRUE)), "POLYGON")
  det_curve <- function(d) plogis(2.5 - d / 100)
  g <- withr::with_seed(31, lapply(1:10, function(i) {
    path <- glatos::crw_in_polygon(poly, theta = c(0, 10), stepLen = 30,
                                   initPos = c(725300, 4935600), nsteps = 1000,
                                   sp_out = FALSE, show_progress = FALSE)
    trns <- glatos::transmit_along_path(path, vel = 0.5, delayRng = c(60, 180), burstDur = 5,
                                        sp_out = FALSE)
    dtc <- glatos::detect_transmissions(trns, rx[c("x", "y")], detRngFun = det_curve,
                                        inputCRS = 32617, sp_out = FALSE, show_progress = FALSE)
    heard <- matrix(FALSE, nrow(trns), 6)
    heard[cbind(dtc$trns_id, dtc$rec_id)] <- TRUE
    list(d = sqrt(outer(trns$x, rx$x, "-")^2 + outer(trns$y, rx$y, "-")^2), heard = heard)
  }))
  fine <- disagg(lake, 3)
  cells <- which(!is.na(values(fine, mat = FALSE)))
  xy <- xyFromCell(fine, cells)
  dc <- sqrt(outer(xy[, 1], rx$x, "-")^2 + outer(xy[, 2], rx$y, "-")^2)
  ff <- list(eta = 2.5 - dc / 100, cells = cells, receivers = rx$receiver_serial, grid = fine)
  mm <- movement_model(function(n) rep(30, n), function(n) rnorm(n, 0, 10 * pi / 180),
                       time_step_s = 60)
  tr5 <- simulate_tracks(fine, mm, 10, 1000,
                         track_control(start = data.frame(x = 725300, y = 4935600)), seed = 32)
  tx5 <- simulate_transmissions(tr5, seed = 33)
  det5 <- simulate_detections(tx5, ff, seed = 34)
  a_d <- sqrt(outer(tx5$x, rx$x, "-")^2 + outer(tx5$y, rx$y, "-")^2)
  a_h <- matrix(FALSE, nrow(tx5), 6)
  a_h[cbind(det5$detections$transmission_id,
            match(det5$detections$receiver_serial, rx$receiver_serial))] <- TRUE
  g_d <- do.call(rbind, lapply(g, `[[`, "d"))
  g_h <- do.call(rbind, lapply(g, `[[`, "heard"))
  br <- c(seq(0, 800, 100), Inf)
  gb <- cut(g_d, br, right = FALSE)
  ab <- cut(a_d, br, right = FALSE)
  oracle <- do.call(rbind, lapply(levels(gb), function(b) {
    x <- c(sum(g_h[gb == b]), sum(a_h[ab == b]))
    nn <- c(sum(gb == b), sum(ab == b))
    ci <- suppressWarnings(prop.test(x, nn, conf.level = 1 - 0.05 / nlevels(gb),
                                     correct = FALSE)$conf.int)
    data.frame(band = b, pairs_glatos = nn[1], rate_glatos = x[1] / nn[1],
               pairs_arraydesign = nn[2], rate_arraydesign = x[2] / nn[2],
               diff_lo = ci[1], diff_hi = ci[2], agree = ci[1] <= 0 && ci[2] >= 0)
  }))
  print(oracle, digits = 3)
  cat(sprintf("transmissions: glatos %d, arraydesign %d over %g track-hours each\n",
              nrow(g_d), nrow(a_d), 10 * 1000 * 60 / 3600))
  write.csv(oracle, file.path(out, "06_glatos_oracle.csv"), row.names = FALSE)
  png(file.path(out, "06_glatos_oracle.png"), width = 1600, height = 1100, res = 200)
  mid <- c(seq(50, 750, 100), 850)
  plot(mid, oracle$rate_glatos, pch = 1, ylim = c(0, 1), xlab = "Distance band (m)",
       ylab = "Detection rate", main = "Detection rate per 100 m band")
  points(mid + 15, oracle$rate_arraydesign, pch = 19)
  curve(plogis(2.5 - x / 100), 0, 900, add = TRUE, lty = 2)
  legend("topright", c("glatos 0.9.8", "arraydesign", "Curve"), pch = c(1, 19, NA),
         lty = c(NA, NA, 2), bty = "n")
  invisible(dev.off())
} else {
  cat(if (requireNamespace("glatos", quietly = TRUE)) {
    sprintf("oracle skipped: glatos %s is older than 0.9.8\n", packageVersion("glatos"))
  } else if (nzchar(system.file(package = "glatos"))) {
    "oracle skipped: glatos is installed but does not load; a dependency is missing\n"
  } else {
    "oracle skipped: glatos 0.9.8 not on the library path\n"
  })
}

# ---- 7. Stoney Lake raster: placement and the land test ----
stoney_file <- Sys.getenv("STONEY_30M")
if (nzchar(stoney_file)) {
  st <- rast(readRDS(stoney_file))
  cat("Stoney grid:", nrow(st), "x", ncol(st), "cells,", sum(!is.na(values(st))), "water\n")
  cat("Stoney raster values (stored):", range(values(st), na.rm = TRUE), "\n")
  st_z <- arraydesign:::depth_down(arraydesign:::water_grid(st),
                                   which(!is.na(values(st, mat = FALSE))))
  cat(sprintf("Stoney depth: %d of %d water cells (%.1f%%) are 1 m deep or less\n",
              sum(st_z <= 1), length(st_z), 100 * mean(st_z <= 1)))
  st_deep <- ifel(st < 1, NA, st)
  place_tab <- do.call(rbind, lapply(c(5, 48, 50, 200), function(n) {
    t0 <- Sys.time()
    g <- place_receivers(st, n = n)
    r <- place_receivers(st, n = n, method = "random", seed = 1)
    gd <- place_receivers(st_deep, n = n)
    data.frame(n = n, grid = nrow(g), grid_trimmed = attr(g, "trimmed"),
               grid_on_land = sum(is.na(extract(st, cbind(g$x, g$y))[, 1])), random = nrow(r),
               random_on_land = sum(is.na(extract(st, cbind(r$x, r$y))[, 1])),
               grid_at_1m_or_less = sum(g$deploy_z <= 1),
               random_at_1m_or_less = sum(r$deploy_z <= 1),
               masked_grid = nrow(gd), masked_under_1m = sum(gd$deploy_z < 1),
               seconds = round(as.numeric(Sys.time() - t0, units = "secs"), 2))
  }))
  print(place_tab)
  st_land <- land_polygons(st)
  st_tab <- do.call(rbind, lapply(c(50, 90), function(s) {
    t0 <- Sys.time()
    trs <- simulate_tracks(st, crw_model(s, s / 2, 45), n_animals = 20, n_steps = 1000, seed = s)
    secs <- as.numeric(Sys.time() - t0, units = "secs")
    segs <- step_segments(trs, crs = st_crs(st_land))
    esc <- attr(trs, "escapes")
    data.frame(mean_step_m = s, steps = sum(esc[, -1]), seconds = round(secs, 1),
               positions_on_land = sum(is.na(extract(st, cbind(trs$x, trs$y))[, 1])),
               segments_crossing_land = sum(lengths(st_intersects(segs, st_land)) > 0),
               reflected = sum(esc$reflected), held = sum(esc$held))
  }))
  print(st_tab)
  trs <- simulate_tracks(st, crw_model(50, 25, 45), n_animals = 8, n_steps = 3000, seed = 7)
  # The tracks figure is drawn for a read on screen and is not kept in the
  # repository, so it goes to tempdir().
  png(file.path(tempdir(), "07_stoney_tracks.png"), width = 1800, height = 1350, res = 150)
  # Inside if(), a ggplot is not printed unless asked, and png() then writes no file.
  print(plot(trs) + ggplot2::ggtitle("Stoney Lake, 8 paths x 3,000 steps of 50 m mean"))
  dev.off()
  st_wg <- arraydesign:::water_grid(st)
  st_landr <- rast(st)
  values(st_landr) <- ifelse(st_wg$wet, NA, 1)
  st_dist <- values(distance(st_landr), mat = FALSE)
  occ <- rbind(occ,
               occupancy(st, crw_model(30, 10, 30), st_dist, "Stoney, crw(30, 10, 30)"),
               occupancy(st, crw_model(90, 30, 30), st_dist, "Stoney, crw(90, 30, 30)"))
}

# ---- 8. Occupancy by distance from land ----
print(occ, digits = 3)
# The table and figure in validation/outputs/ hold the Stoney walks; a run
# without them prints the example-lake rows and leaves both files as they are.
occ_out <- if (nzchar(stoney_file)) out else tempdir()
if (!nzchar(stoney_file)) {
  cat("08_shore_occupancy.csv and .png go to tempdir(): STONEY_30M is unset\n")
}
write.csv(occ, file.path(occ_out, "08_shore_occupancy.csv"), row.names = FALSE)
walks <- unique(occ$walk)
png(file.path(occ_out, "08_shore_occupancy.png"), width = 1600, height = 1100, res = 200)
op <- par(mar = c(4, 4.2, 1, 1))
xs <- seq_len(4)
plot(NA, xlim = c(0.6, 4.4), ylim = c(0.7, 1.3), xaxt = "n", xlab = "Distance from land",
     ylab = "Share of positions / share of water")
axis(1, at = xs, labels = levels(factor(occ$band, levels = unique(occ$band))))
abline(h = 1, lty = 2)
pal <- c("#E69F00", "#D55E00", "#009E73", "#000000")
for (i in seq_along(walks)) {
  d <- occ[occ$walk == walks[i], ]
  xi <- xs + (i - (length(walks) + 1) / 2) * 0.12
  arrows(xi, (d$position_share - 2 * d$se) / d$water_share, xi,
         (d$position_share + 2 * d$se) / d$water_share, angle = 90, code = 3,
         length = 0.03, col = pal[i])
  points(xi, d$ratio, pch = 19, col = pal[i])
}
legend("topright", walks, col = pal[seq_along(walks)], pch = 19, bty = "n", cex = 0.8)
par(op)
dev.off()
sessionInfo()
