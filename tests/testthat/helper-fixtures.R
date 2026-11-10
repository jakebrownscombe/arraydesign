# ---- Test fixtures ----
# Rasters built in code, and a detection field with the shape of a
# detfield::de_field: eta (cells by receivers, logit scale), cells, receivers,
# grid, model and distances. The detection probability is a logistic decline
# over straight-line distance from the cell centre to the receiver.

lake_raster <- function(water, res = 30, xmin = 700000, ymin = 4900000) {
  r <- terra::rast(nrows = nrow(water), ncols = ncol(water), xmin = xmin,
                   xmax = xmin + res * ncol(water), ymin = ymin,
                   ymax = ymin + res * nrow(water), crs = "EPSG:32617")
  terra::values(r) <- ifelse(as.vector(t(water)), 10, NA)
  r
}

# Two basins of 48 by 58 cells separated by a land wall one cell (30 m) wide,
# with a one-cell land margin all round. No water path joins the basins.
wall_lake <- function() {
  w <- matrix(TRUE, 60, 100)
  w[c(1, 60), ] <- FALSE
  w[, c(1, 100)] <- FALSE
  w[, 50] <- FALSE
  lake_raster(w)
}

# Open water of 200 by 200 cells inside a one-cell land margin.
open_lake <- function() {
  w <- matrix(TRUE, 202, 202)
  w[c(1, 202), ] <- FALSE
  w[, c(1, 202)] <- FALSE
  lake_raster(w)
}

# A lake on the 421 by 560 grid of the Stoney Lake depth raster at 30 m, with
# about as many water cells (Stoney has 46,097): an irregular ellipse with two
# islands.
stoney_sized <- function() {
  rows <- matrix(seq_len(421), 421, 560)
  cols <- matrix(seq_len(560), 421, 560, byrow = TRUE)
  ang <- atan2(rows - 211, cols - 280)
  radius <- 1 + 0.12 * sin(5 * ang) + 0.05 * cos(11 * ang)
  w <- ((cols - 280) / 205)^2 + ((rows - 211) / 72)^2 < radius^2
  w[((cols - 200) / 12)^2 + ((rows - 200) / 9)^2 < 1] <- FALSE
  w[((cols - 360) / 15)^2 + ((rows - 230) / 10)^2 < 1] <- FALSE
  lake_raster(w)
}

# A de_field-shaped list with a logistic curve: p = plogis(a - d / s).
make_field <- function(grid, receivers, a = 2.5, s = 100) {
  cells <- which(!is.na(terra::values(grid, mat = FALSE)))
  xy <- terra::xyFromCell(grid, cells)
  d <- sqrt(outer(xy[, 1], receivers$x, "-")^2 + outer(xy[, 2], receivers$y, "-")^2)
  structure(
    list(eta = a - d / s, cells = cells, receivers = receivers$receiver_serial, grid = grid,
         model = list(a = a, s = s),
         distances = data.frame(cell = rep(cells, ncol(d)),
                                receiver_serial = rep(receivers$receiver_serial, each = nrow(d)),
                                distance_m = as.vector(d))),
    class = "de_field"
  )
}

# Land cells of a grid as sf polygons, for an oracle that does not share code
# with the walker.
land_polygons <- function(grid) {
  land <- terra::ifel(is.na(grid), 1, NA)
  sf::st_as_sf(terra::as.polygons(land, dissolve = TRUE))
}

# Step segments of a track as sf linestrings, one per moving step.
step_segments <- function(tracks, crs = 32617) {
  n <- nrow(tracks)
  keep <- which(tracks$path_id[-1] == tracks$path_id[-n] &
                  (tracks$x[-1] != tracks$x[-n] | tracks$y[-1] != tracks$y[-n]))
  lines <- lapply(keep, function(i) {
    sf::st_linestring(rbind(c(tracks$x[i], tracks$y[i]), c(tracks$x[i + 1], tracks$y[i + 1])))
  })
  sf::st_sfc(lines, crs = crs)
}

# Wraps angles to [-pi, pi).
wrap_angle <- function(a) (a + pi) %% (2 * pi) - pi

# A toy lake: a grid of 1,690 water cells of 30 m and five deployments. The
# tests keep their own copy so that their numbers do not depend on detfield's
# example data.
toy_lake_fixture <- readRDS(test_path("fixtures", "toy_lake.rds"))
