#' Place receivers on the water cells of a raster
#'
#' Places a candidate receiver array on a depth raster in which NA cells are
#' land. Three methods:
#'
#' * `"grid"`: a square lattice. With `spacing`, the lattice is centred on the
#'   water and every node on a water cell becomes a receiver. With `n`, the
#'   spacing and the lattice offset are searched until exactly `n` nodes fall
#'   on water; spacings below the cell size are not tried. If no lattice
#'   gives exactly `n`, the nearest lattice with more nodes is taken and the
#'   surplus nodes with the fewest water neighbours, which are the nodes
#'   nearest the shore, are dropped.
#' * `"random"`: `n` water cells drawn at random. With `spacing` as well, no
#'   two receivers are closer than
#'   `spacing`, and the function stops if `n` such receivers do not fit.
#' * `"manual"`: the points in `coords`. A point on a land cell moves to the
#'   centre of the nearest water cell when the edge of that cell is within
#'   one cell diagonal, the rule [detfield::water_distances()] applies, so
#'   both packages take or refuse the same points. It covers points that
#'   rasterisation of the shoreline put on land. A point farther from water,
#'   or outside the raster, is refused with an error that names it, or the
#'   first 10 and the count. Points outside the raster are refused first, so
#'   on an input that has both kinds the error names only those; the points
#'   on land are named once every point is on the raster. Two receivers in
#'   one cell share every cell-based detection probability, so the function
#'   warns when that happens, naming the rows, or the first 10 and the count;
#'   repeated `receiver_serial` values are refused.
#'
#' The `"grid"` and `"random"` methods return exactly `n` receivers when `n`
#' is given. They place each receiver at the centre of its cell, so
#' coordinates written to a file and read back with `"manual"` fall in the
#' same cells.
#'
#' Every water cell is a candidate, however shallow. To keep receivers out of
#' water too shallow to deploy in, set those cells to NA first, for example
#' `place_receivers(terra::ifel(grid < 1, NA, grid), n = 48)` keeps receivers
#' in 1 m of water or more when depth is stored as positive metres.
#'
#' @param grid A [terra::SpatRaster] of water depth in metres, with NA on
#'   land, in a projected CRS whose unit is the metre. A raster in longitude
#'   and latitude or in another unit, such as feet, is refused. The scale of
#'   the CRS is measured at the raster centre on its two axes and its
#'   diagonals. A CRS is refused when any of the four is more than 1 percent
#'   from 1, as in Web Mercator (EPSG:3857) away from the equator, whose
#'   metres are not metres on the ground. The message names the scale and the
#'   UTM zone of the centre. A CRS whose scale cannot be measured, such as a
#'   local engineering CRS, is refused too. Depth
#'   may be stored as positive or as negative numbers; `deploy_z` is returned
#'   positive down.
#' @param n Number of receivers, for `"grid"` and `"random"`.
#' @param spacing Distance between receivers in metres, at least the cell
#'   size. For `"grid"`, the lattice spacing, used when `n` is `NULL`. For
#'   `"random"`, the smallest allowed distance between two receivers.
#' @param method One string: `"grid"` (the default), `"random"` or `"manual"`.
#' @param coords For `"manual"`, a data frame with columns `x` and `y` in the
#'   CRS of `grid`, or `deploy_lon` and `deploy_lat` in WGS 84. Columns
#'   `receiver_serial` and `station_id`, when present, are kept.
#' @param seed Seed for `"random"` and for the tie-breaks of `"grid"`; `NULL`
#'   uses the current random number stream. The caller's stream is restored
#'   after a seeded call. The draws depend on the generator kind as well as
#'   the seed; the kind is stored with the result.
#'
#' @return A data frame of class `receiver_layout`, one row per receiver, with
#'   the columns of the positionR table contract: `receiver_serial` (`"R001"`,
#'   `"R002"`, ...), `station_id` (`"S001"`, ...), `x` and `y` in the CRS of
#'   `grid`, `deploy_lat` and `deploy_lon` in WGS 84, and `deploy_z`, the
#'   depth of the cell in metres. Attributes: `epsg`, the EPSG code of `grid`
#'   (NA when terra cannot name one), `method`, `spacing` (the lattice spacing
#'   used, for `"grid"`), `trimmed` (surplus lattice nodes dropped, for
#'   `"grid"` with `n`), `snapped` (points moved off land, for `"manual"`)
#'   and `rng_kind` (the [RNGkind()] in force).
#' @source The example raster `stoney_depth.tif` in
#'   `system.file("extdata", package = "arraydesign")` has the same grid and
#'   depths as `depth_raster` of positionR 1.7.0, a 100 m depth raster
#'   of Stoney Lake, derived from the bathymetry made before 2024,
#'   with the CRS set to EPSG:32617. Depth is zero or negative metres and land
#'   is `NA`; 56 of its 4,818 water cells hold 0.
#' @export
#' @examples
#' # Stoney Lake, Ontario: 100 m cells, depth zero or negative metres, land NA.
#' stoney <- terra::rast(system.file("extdata", "stoney_depth.tif",
#'                                   package = "arraydesign"))
#' # The 20 planned stations of vignette("arraydesign"), as latitude and longitude.
#' planned <- read.csv(system.file("extdata", "planned_stations.csv",
#'                                 package = "arraydesign"))
#' rx <- place_receivers(stoney, method = "manual", coords = planned)
#' nrow(rx)
#' head(rx)
#' plot(rx, stoney)
#'
#' # No receiver can be moored in water under 1 m, so drop it before a grid or
#' # random layout. The grid of 20 is the planned array.
#' deployable <- terra::ifel(abs(stoney) < 1, NA, stoney)
#' place_receivers(deployable, n = 20, method = "grid", seed = 1)
#' place_receivers(deployable, n = 20, method = "random", spacing = 1000, seed = 1)
#' nrow(place_receivers(deployable, spacing = 2000))
place_receivers <- function(grid, n = NULL, spacing = NULL,
                            method = c("grid", "random", "manual"),
                            coords = NULL, seed = NULL) {
  place_layout(grid, n, spacing, method, coords, seed, call = sys.call())
}

# The body of place_receivers(). `call` is the user's call that a grid or
# placement refusal names, so evaluate_array() can pass its own.
place_layout <- function(grid, n, spacing, method, coords, seed, call) {
  method <- match_choice(method, c("grid", "random", "manual"), "method")
  wg <- water_grid(grid, require_crs = TRUE, call = call)
  check_placement_args(wg, n, spacing, method, coords)
  placed <- with_seed_or_not(seed, switch(
    method,
    grid = place_grid(wg, n, spacing),
    random = list(xy = random_cells(wg, n, spacing)),
    manual = manual_points(wg, coords, call)
  ))
  # A grid node can fall exactly on the edge between two cells, and a
  # coordinate written to a file and read back can then land in the other
  # cell, with another depth. Grid and random receivers sit at the centre of
  # their cell, so a round trip always finds the same cell.
  if (method != "manual") placed$xy <- cell_xy(wg, cell_at(wg, placed$xy[, 1], placed$xy[, 2]))
  out <- receiver_table(wg, placed$xy, placed$ids)
  attr(out, "method") <- method
  attr(out, "spacing") <- placed$spacing
  attr(out, "trimmed") <- placed$trimmed
  attr(out, "snapped") <- placed$snapped
  attr(out, "rng_kind") <- RNGkind()
  out
}

check_placement_args <- function(wg, n, spacing, method, coords) {
  if (method == "manual") {
    stop_unless(is.null(n) && is.null(spacing),
                "`n` and `spacing` are not used with method = \"manual\"; give `coords` only.")
    stop_unless(is.data.frame(coords) && nrow(coords) > 0,
                "method = \"manual\" needs `coords`, a data frame with at least one row.")
    stop_unless(!anyDuplicated(coords$receiver_serial),
                "`coords$receiver_serial` repeats: ",
                paste(unique(coords$receiver_serial[duplicated(coords$receiver_serial)]),
                      collapse = ", "), ".")
    return(invisible())
  }
  stop_unless(is.null(n) || is_scalar(n, 1, whole = TRUE),
              "`n` must be one whole number of at least 1.")
  stop_unless(is.null(n) || n <= sum(wg$wet),
              "`n` is ", n, " but `grid` has only ", sum(wg$wet), " water cells.")
  stop_unless(is.null(spacing) || is_scalar(spacing, wg$res),
              "`spacing` must be one number of metres, at least the cell size (", wg$res, " m).")
  needs <- if (method == "grid") "`n` or `spacing`" else "`n`"
  stop_unless(!is.null(n) || (method == "grid" && !is.null(spacing)),
              "method = \"", method, "\" needs ", needs, ".")
}

# ---- Grid placement ----
place_grid <- function(wg, n, spacing) {
  if (!is.null(n)) return(grid_exact_n(wg, n))
  xy <- lattice_nodes(wg, spacing, 0, 0)
  stop_unless(nrow(xy) > 0, "No lattice node falls on water at spacing ", spacing, " m.")
  list(xy = xy, spacing = spacing, trimmed = 0L)
}

# Nodes of a square lattice at spacing `s`, centred on the water cells and
# shifted by fractions `ox` and `oy` of the spacing, kept where they fall on water.
lattice_nodes <- function(wg, s, ox, oy, box = water_box(wg)) {
  cx <- mean(box[1:2]) + ox * s
  cy <- mean(box[3:4]) + oy * s
  kx <- seq(floor((box[1] - cx) / s), ceiling((box[2] - cx) / s))
  ky <- seq(floor((box[3] - cy) / s), ceiling((box[4] - cy) / s))
  xy <- cbind(x = rep(cx + kx * s, times = length(ky)), y = rep(cy + ky * s, each = length(kx)))
  xy[is_wet(wg, xy[, 1], xy[, 2]), , drop = FALSE]
}

# Bounding box of the water cell centres: xmin, xmax, ymin, ymax.
water_box <- function(wg) {
  wxy <- cell_xy(wg, which(wg$wet))
  c(range(wxy[, 1]), range(wxy[, 2]))
}

# Searches spacing and offset for a lattice with exactly n nodes on water,
# trying spacings nearest the area-based guess first.
grid_exact_n <- function(wg, n) {
  s0 <- wg$res * sqrt(sum(wg$wet) / n)
  # Spacings below the cell size are not tried: two receivers in one cell
  # would share every cell-based detection probability.
  f <- seq(0.8, 1.25, by = 0.002)
  f <- f[order(abs(log(f)))]
  f <- f[s0 * f >= wg$res]
  box <- water_box(wg)
  best <- NULL
  for (fi in f) {
    found <- lattice_offsets(wg, s0 * fi, box, n)
    if (nrow(found$xy) == n) return(list(xy = found$xy, spacing = s0 * fi, trimmed = 0L))
    if (nrow(found$xy) > n && (is.null(best) || nrow(found$xy) < nrow(best$xy))) {
      best <- list(xy = found$xy, spacing = s0 * fi)
    }
  }
  stop_unless(!is.null(best), "No square lattice puts ", n, " or more receivers on water; ",
              "use method = \"random\".")
  best$trimmed <- nrow(best$xy) - as.integer(n)
  best$xy <- trim_surplus(wg, best$xy, n)
  best
}

# The lattice at spacing s with exactly n nodes on water among 16 offsets,
# or else the one with the fewest nodes above n (none above n gives 0 rows).
lattice_offsets <- function(wg, s, box, n) {
  offsets <- expand.grid(ox = c(0, 0.5, 0.25, 0.75), oy = c(0, 0.5, 0.25, 0.75))
  best <- matrix(numeric(), 0, 2)
  for (j in seq_len(nrow(offsets))) {
    xy <- lattice_nodes(wg, s, offsets$ox[j], offsets$oy[j], box)
    if (nrow(xy) == n) return(list(xy = xy))
    if (nrow(xy) > n && (!nrow(best) || nrow(xy) < nrow(best))) best <- xy
  }
  list(xy = best)
}

# Drops the nodes nearest the shore until n remain: those with the fewest
# water cells among the eight around them, ties broken at random.
trim_surplus <- function(wg, xy, n) {
  wet_around <- vapply(seq_len(nrow(xy)), function(i) {
    dx <- rep(c(-1, 0, 1), 3) * wg$res
    dy <- rep(c(-1, 0, 1), each = 3) * wg$res
    sum(is_wet(wg, xy[i, 1] + dx, xy[i, 2] + dy))
  }, numeric(1))
  keep <- order(-wet_around, stats::runif(nrow(xy)))[seq_len(n)]
  xy[sort(keep), , drop = FALSE]
}

# ---- Random placement ----
random_cells <- function(wg, n, spacing) {
  wet_cells <- which(wg$wet)
  if (is.null(spacing)) return(cell_xy(wg, wet_cells[sample.int(length(wet_cells), n)]))
  candidates <- cell_xy(wg, wet_cells[sample.int(length(wet_cells))])
  chosen <- matrix(numeric(), 0, 2, dimnames = list(NULL, c("x", "y")))
  for (i in seq_len(nrow(candidates))) {
    d2 <- (chosen[, 1] - candidates[i, 1])^2 + (chosen[, 2] - candidates[i, 2])^2
    if (all(d2 >= spacing^2)) chosen <- rbind(chosen, candidates[i, ])
    if (nrow(chosen) == n) return(chosen)
  }
  stop_classed("placement", "Only ", nrow(chosen), " receivers fit at a spacing of ", spacing,
               " m; asked for ", n, ".")
}

# ---- Manual placement ----
manual_points <- function(wg, coords, call = NULL) {
  xy <- manual_xy(wg, coords)
  cell <- cell_at(wg, xy[, 1], xy[, 2])
  outside <- which(is.na(cell))
  if (length(outside)) {
    rt <- rows_text(outside)
    stop_classed("placement", "`coords` has ", rt$count, " outside the raster: ", rt$shown, ".",
                 call = call)
  }
  dry <- which(!wg$wet[cell])
  wxy <- cell_xy(wg, which(wg$wet))
  # For each point on land: the gap to the edge of the nearest water cell,
  # measured as detfield's snap measures it, so a station that
  # water_distances() accepts is accepted here too; and the water cell whose
  # centre is nearest, where the point moves.
  gap <- numeric(length(dry))
  nearest <- integer(length(dry))
  for (j in seq_along(dry)) {
    i <- dry[j]
    to_edge <- sqrt(pmax(abs(wxy[, 1] - xy[i, 1]) - wg$res / 2, 0)^2 +
                      pmax(abs(wxy[, 2] - xy[i, 2]) - wg$res / 2, 0)^2)
    gap[j] <- min(to_edge)
    nearest[j] <- which.min((wxy[, 1] - xy[i, 1])^2 + (wxy[, 2] - xy[i, 2])^2)
  }
  refuse_far_on_land(dry, gap, sqrt(2) * wg$res, call)
  xy[dry, ] <- wxy[nearest, , drop = FALSE]
  if (length(dry)) message(length(dry), " of ", nrow(xy), " points moved off land.")
  warn_shared_cells(wg, xy)
  list(xy = xy, ids = coords[intersect(c("receiver_serial", "station_id"), names(coords))],
       snapped = length(dry))
}

# Refuses the points on land farther than `limit` from the edge of water.
# One such point is named with its gap; more are counted and at most 10
# named. `call` is the user's call that the refusal names.
refuse_far_on_land <- function(dry, gap, limit, call = NULL) {
  far <- which(gap > limit)
  if (length(far) == 1L) {
    stop_classed("placement", "`coords` row ", dry[far], " is on land ", round(gap[far]),
                 " m from the edge of the nearest water cell, more than one cell diagonal (",
                 round(limit), " m).", call = call)
  }
  if (length(far)) {
    rt <- rows_text(dry[far])
    stop_classed("placement", "`coords` has ", rt$count, " on land more than one cell ",
                 "diagonal (", round(limit), " m) from the edge of the nearest water cell: ",
                 rt$shown, ".", call = call)
  }
  invisible()
}

# Receivers in one cell cannot be told apart by a cell-based field.
warn_shared_cells <- function(wg, xy) {
  cell <- cell_at(wg, xy[, 1], xy[, 2])
  shared <- which(cell %in% cell[duplicated(cell)])
  if (length(shared)) {
    rt <- rows_text(shared)
    warning("`coords` has ", rt$count, " that share a grid cell with another row, so they ",
            "share every detection probability: ", rt$shown, ".", call. = FALSE)
  }
  invisible()
}

manual_xy <- function(wg, coords) {
  if (all(c("x", "y") %in% names(coords))) {
    return(numeric_xy(coords$x, coords$y, "coords"))
  }
  stop_unless(all(c("deploy_lon", "deploy_lat") %in% names(coords)),
              "`coords` needs columns `x` and `y`, or `deploy_lon` and `deploy_lat`.")
  ll <- numeric_xy(coords$deploy_lon, coords$deploy_lat, "coords",
                   columns = c("deploy_lon", "deploy_lat"))
  xy <- terra::project(unname(ll), from = "EPSG:4326", to = terra::crs(wg$raster))
  colnames(xy) <- c("x", "y")
  xy
}

# ---- Output table ----
receiver_table <- function(wg, xy, ids = NULL) {
  n <- nrow(xy)
  width <- max(3L, nchar(n))
  ll <- terra::project(unname(xy[, 1:2, drop = FALSE]), from = terra::crs(wg$raster),
                       to = "EPSG:4326")
  out <- data.frame(
    receiver_serial = sprintf("R%0*d", width, seq_len(n)),
    station_id = sprintf("S%0*d", width, seq_len(n)),
    x = unname(xy[, 1]), y = unname(xy[, 2]),
    deploy_lat = ll[, 2], deploy_lon = ll[, 1],
    deploy_z = depth_down(wg, cell_at(wg, xy[, 1], xy[, 2])),
    stringsAsFactors = FALSE
  )
  for (nm in names(ids)) out[[nm]] <- as.character(ids[[nm]])
  attr(out, "epsg") <- grid_epsg(wg$raster)
  class(out) <- c("receiver_layout", "data.frame")
  out
}
