# ---- Water grid helpers ----
# The walker, the transmission model and the receiver placement all read the
# raster through one small list, so cell arithmetic is done in R without a
# terra call per step. Cells are numbered as terra numbers them: row by row
# from the top left, starting at 1.

# Turns a SpatRaster (or a raster package layer) into the list the internal
# functions use. NA cells are land. `arg` is the name of the caller's
# argument, so that a message names what the user passed. `call` is the
# user's call that a refusal names; NULL names none.
water_grid <- function(grid, require_crs = FALSE, arg = "grid", call = NULL) {
  a <- sprintf("`%s`", arg)
  if (inherits(grid, "PackedSpatRaster")) grid <- terra::unwrap(grid)
  if (inherits(grid, c("RasterLayer", "RasterBrick", "RasterStack"))) grid <- terra::rast(grid)
  if (!inherits(grid, "SpatRaster")) {
    stop_classed("grid", a, " must be a terra SpatRaster in which NA cells are land.",
                 call = call)
  }
  if (terra::nlyr(grid) != 1L) {
    stop_classed("grid", a, " must have one layer; it has ", terra::nlyr(grid), ".",
                 call = call)
  }
  check_metres(grid, arg, call)
  res <- terra::res(grid)
  if (abs(res[1] - res[2]) > 1e-6 * res[1]) {
    stop_classed("grid", a, " cells must be square; resolution is ", res[1], " by ", res[2], ".",
                 call = call)
  }
  if (require_crs && !nzchar(terra::crs(grid))) {
    stop_classed("grid", a, " has no CRS, so latitude and longitude cannot be computed. ",
                 "Set one with terra::crs().", call = call)
  }
  vals <- terra::values(grid, mat = FALSE)
  wet <- !is.na(vals)
  if (!any(wet)) stop_classed("grid", a, " has no water cells: every cell is NA.", call = call)
  e <- as.vector(terra::ext(grid))
  list(
    raster = grid, values = vals, wet = wet, res = res[1], arg = arg,
    xmin = e[["xmin"]], xmax = e[["xmax"]], ymin = e[["ymin"]], ymax = e[["ymax"]],
    nrow = terra::nrow(grid), ncol = terra::ncol(grid)
  )
}

# Refuses a grid whose CRS is in longitude and latitude, in a unit other
# than the metre, or at a scale at the raster centre that is far from 1 or
# cannot be measured. Spacing, step lengths and ranges are read as metres on
# the ground, so a grid in feet or in Web Mercator would run without a sign
# of the error and place more receivers at a spacing. A grid with no CRS
# passes here; water_grid() refuses it when the caller sets require_crs.
check_metres <- function(grid, arg, call) {
  a <- sprintf("`%s`", arg)
  if (!nzchar(terra::crs(grid))) return(invisible())
  if (isTRUE(terra::is.lonlat(grid))) {
    stop_classed("grid", a, " is in longitude and latitude.\nProject it to a CRS in metres first.",
                 call = call)
  }
  if (!isTRUE(all.equal(terra::linearUnits(grid), 1))) {
    stop_classed("grid", a, " has a CRS whose unit is ", crs_unit(grid),
                 ".\nProject it to a CRS in metres first.", call = call)
  }
  check_scale(grid, call, arg)
}

# ---- Scale of the CRS at the raster centre ----
# A projected CRS in metres can still give lengths that differ from lengths
# on the ground. Web Mercator has a scale of about 1 / cos(latitude), 1.41
# at 45 degrees; at the equator it is 1.000 east-west and 1.007
# north-south. A sinusoidal CRS away from its central meridian has a scale
# near 1 on both axes and is sheared, so its scale on the diagonals runs
# from 0.93 to 1.07 at 11 E 45.6 N. The check therefore measures four
# segments. UTM stays within 0.1 percent across a zone.

# The length of each segment, in metres of the grid, and the largest
# departure from 1 that any of the four scales may have.
scale_segment_m <- 1000
max_scale_error <- 0.01

# The scale of a grid's CRS at the raster centre on four 1 km segments
# through it: east-west, north-south and the two diagonals. Each scale is
# the segment's length in the grid over its geodesic length on WGS 84. Also
# returns the centre in longitude and latitude. The four scales are NA when
# the points cannot be projected to finite longitudes and latitudes, as in
# a local engineering CRS. The warnings of a failed projection are dropped,
# so that the refusal of check_scale() is the only condition the user sees.
crs_scale <- function(grid) {
  e <- as.vector(terra::ext(grid))
  cx <- (e[["xmin"]] + e[["xmax"]]) / 2
  cy <- (e[["ymin"]] + e[["ymax"]]) / 2
  half <- scale_segment_m / 2
  # Half of each segment as an x and a y offset from the centre.
  dx <- c(half, 0, half / sqrt(2), half / sqrt(2))
  dy <- c(0, half, half / sqrt(2), -half / sqrt(2))
  pts <- cbind(c(cx, cx - dx, cx + dx), c(cy, cy - dy, cy + dy))
  ll <- tryCatch(
    suppressWarnings(terra::project(pts, from = terra::crs(grid), to = "EPSG:4326")),
    error = function(e) matrix(NA_real_, nrow(pts), 2)
  )
  scale <- rep(NA_real_, 4)
  if (all(is.finite(ll))) {
    geo <- terra::distance(ll[2:5, , drop = FALSE], ll[6:9, , drop = FALSE],
                           lonlat = TRUE, pairwise = TRUE)
    scale <- scale_segment_m / geo
  }
  names(scale) <- c("east-west", "north-south", "south-west to north-east",
                    "north-west to south-east")
  list(scale = scale, lon = ll[1, 1], lat = ll[1, 2])
}

# Refuses a grid whose CRS has a scale that cannot be measured, or a scale
# more than max_scale_error from 1 on any of the four segments of
# crs_scale(). Without a measured scale, metres of the grid cannot be
# compared with metres on the ground. The second message names the four
# scales, or one when all four print alike, and the UTM zone of the centre.
# `call` is the user's call that the refusal names; `arg` is the argument.
check_scale <- function(grid, call = NULL, arg = "grid") {
  s <- crs_scale(grid)
  if (anyNA(s$scale)) {
    stop_classed("grid", "`", arg, "` is refused. The scale of this CRS at the raster centre ",
                 "cannot be measured, so give the raster in a projected CRS with a geographic ",
                 "datum, such as its UTM zone.", call = call)
  }
  if (all(abs(s$scale - 1) <= max_scale_error)) return(invisible())
  shown <- sprintf("%.3f", s$scale)
  scale_text <- if (all(shown == shown[1])) {
    paste(shown[1], "on its two axes and its diagonals")
  } else {
    paste0(shown[1], " east-west, ", shown[2], " north-south, ", shown[3], " ",
           names(s$scale)[3], " and ", shown[4], " ", names(s$scale)[4])
  }
  epsg <- utm_epsg(s$lon, s$lat)
  stop_classed("grid", "`", arg, "` has a CRS whose scale at the raster centre is ", scale_text,
               ", more than ", 100 * max_scale_error,
               " percent from 1, so its metres are not metres on the ground.\n",
               "Project it to the UTM zone of the centre, ", epsg %% 100,
               if (epsg < 32700) "N" else "S", " (EPSG:", epsg, ").", call = call)
}

# EPSG code of the WGS 84 UTM zone that holds a point: 326zz north of the
# equator, 327zz south. A longitude of exactly 180 falls in zone 60.
utm_epsg <- function(lon, lat) {
  zone <- min(floor((lon + 180) / 6) + 1, 60)
  as.integer(if (lat >= 0) 32600 + zone else 32700 + zone)
}

# The name of a projected CRS's axis unit and its size in metres, such as
# "the US survey foot (0.3048006 m)". The name is that of the LENGTHUNIT in
# the WKT whose size equals terra::linearUnits(), which is the size of the
# axis unit. A WKT holds other LENGTHUNITs: the datum's, and in a BOUNDCRS
# the metre of the target ellipsoid.
crs_unit <- function(grid) {
  wkt <- terra::crs(grid)
  size <- terra::linearUnits(grid)
  units <- regmatches(wkt, gregexpr('LENGTHUNIT\\["[^"]+",[-0-9.eE+]+', wkt))[[1]]
  unit_names <- sub('^LENGTHUNIT\\["([^"]+)".*$', "\\1", units)
  unit_sizes <- as.numeric(sub("^.*,", "", units))
  axis <- unit_names[abs(unit_sizes - size) <= 1e-9 * size]
  name <- if (length(axis)) paste("the", axis[1]) else "a unit with no name"
  paste0(name, " (", format(size, digits = 7), " m)")
}

# Cell number of each point, NA outside the extent.
cell_at <- function(wg, x, y) {
  col <- floor((x - wg$xmin) / wg$res)
  row <- floor((wg$ymax - y) / wg$res)
  inside <- col >= 0 & col < wg$ncol & row >= 0 & row < wg$nrow
  out <- row * wg$ncol + col + 1
  out[!inside] <- NA
  out
}

# Centre coordinates of cells.
cell_xy <- function(wg, cell) {
  row <- (cell - 1) %/% wg$ncol
  col <- (cell - 1) %% wg$ncol
  cbind(x = wg$xmin + (col + 0.5) * wg$res, y = wg$ymax - (row + 0.5) * wg$res)
}

is_wet <- function(wg, x, y) {
  cell <- cell_at(wg, x, y)
  !is.na(cell) & wg$wet[ifelse(is.na(cell), 1L, cell)]
}

# TRUE when the straight segment from (x0, y0) to (x1, y1) passes through
# water cells only. The test is exact: it finds every grid line the segment
# crosses and checks one point inside each cell between two crossings, so a
# step cannot jump a land wall however thin, and a corner clip is caught.
segment_wet <- function(wg, x0, y0, x1, y1) {
  is.null(first_land(wg, x0, y0, x1, y1))
}

# Where a segment first enters land: `t`, the fraction along the segment of
# the grid line it crosses into the land cell, and `axis`, 1 when that line is
# a column edge (x is constant) and 2 when it is a row edge. NULL when the
# whole segment is on water. An `axis` of NA means the segment starts on land.
first_land <- function(wg, x0, y0, x1, y1) {
  tu <- line_crossings((x0 - wg$xmin) / wg$res, (x1 - wg$xmin) / wg$res)
  tv <- line_crossings((wg$ymax - y0) / wg$res, (wg$ymax - y1) / wg$res)
  t_cross <- c(tu, tv)
  axis <- rep(1:2, c(length(tu), length(tv)))
  keep <- t_cross > 0 & t_cross < 1
  o <- order(t_cross[keep])
  t_cross <- t_cross[keep][o]
  axis <- axis[keep][o]
  t_all <- c(0, t_cross, 1)
  t_mid <- (t_all[-1] + t_all[-length(t_all)]) / 2
  wet <- is_wet(wg, x0 + t_mid * (x1 - x0), y0 + t_mid * (y1 - y0))
  if (all(wet)) return(NULL)
  j <- which.min(wet)
  if (j == 1L) return(list(t = 0, axis = NA_integer_))
  list(t = t_all[j], axis = axis[j - 1L])
}

# Fractions along a segment at which it crosses integer grid lines of one axis.
line_crossings <- function(a0, a1) {
  lo <- ceiling(min(a0, a1))
  hi <- floor(max(a0, a1))
  if (a0 == a1 || lo > hi) return(numeric())
  (seq(lo, hi) - a0) / (a1 - a0)
}

# EPSG code of a raster as an integer, NA when terra cannot name one. The
# body is that of detfield's grid_epsg(), so both packages store the same
# code for the same raster. Callers pass the raster of water_grid().
grid_epsg <- function(grid) {
  code <- terra::crs(grid, describe = TRUE)$code
  if (length(code) && !is.na(code) && grepl("^[0-9]+$", code)) as.integer(code) else NA_integer_
}

# Depth of water cells as positive metres down. A raster that stores depth as
# negative numbers is read as such; a raster that mixes signs is refused,
# because neither reading is safe.
depth_down <- function(wg, cell) {
  v <- wg$values[wg$wet]
  z <- wg$values[cell]
  if (all(v >= 0)) return(z)
  if (all(v <= 0)) return(-z)
  arg <- if (is.null(wg$arg)) "grid" else wg$arg
  stop_classed("grid", "`", arg, "` values mix positive and negative depths (", sum(v < 0),
               " negative and ", sum(v > 0), " positive water cells). Give depth in metres, ",
               "positive down, with NA on land. For elevation, with water below 0 and land ",
               "above, run x[x >= 0] <- NA, then x <- -x.")
}

# Runs `code` under `seed` and restores the caller's random number state; with
# no seed it uses the caller's stream as it is. The generator kind is the
# caller's, as with set.seed(), so results store it beside the seed.
with_seed_or_not <- function(seed, code) {
  if (is.null(seed)) return(force(code))
  check_seed(seed)
  withr::with_seed(seed, code)
}

# set.seed() takes one finite number whose whole part fits in an R integer,
# and truncates a fraction toward zero. It stops with an unclassed error on
# Inf, NaN, NA or a whole part beyond the integer range, so a seed is checked
# here first and refused with the package's classed error. The check keeps to
# the range the help gives, so it also refuses a fraction beyond either end,
# such as 2147483647.5, which set.seed() truncates to that end.
is_seed <- function(x) is_scalar(x, -.Machine$integer.max, .Machine$integer.max)

check_seed <- function(seed) {
  stop_unless(is.null(seed) || is_seed(seed), "`seed` must be one number from -",
              .Machine$integer.max, " to ", .Machine$integer.max,
              ", the whole numbers set.seed() takes, or NULL; `seed` holds ",
              describe_seed(seed), ".")
}

# Describes a seed value in a few words. A full .Random.seed printed in a
# message runs to thousands of characters, so a vector is given by its length.
describe_seed <- function(s) {
  if (is.null(s)) return("NULL")
  if (is.numeric(s) && length(s) == 1L) return(describe_number(s))
  if (is.numeric(s)) return(paste(length(s), "numbers"))
  paste0("a ", class(s)[1], " of length ", length(s))
}

describe_number <- function(s) {
  if (is.nan(s)) return("NaN")
  if (is.na(s)) return("NA")
  shown <- paste("the number", number_text(s))
  if (is.finite(s) && abs(s) > .Machine$integer.max) {
    shown <- paste0(shown, ", outside -", .Machine$integer.max, " to ", .Machine$integer.max)
  }
  shown
}

# Prints the shortest decimal, up to 17 significant digits, that reads back as
# the same double. Where no count reads back, it prints 17 digits: R's parser
# cannot read some doubles back at any count. A value typed with 17 digits may
# print as a double other than the one typed, differing in its last digits. The
# rule is the same for both signs, and a value of 1e15 or more in absolute size
# prints in scientific form.
number_text <- function(s) {
  if (!is.finite(s)) return(format(s))
  for (d in 15:17) {
    out <- format(s, digits = d, scientific = abs(s) >= 1e15)
    if (as.numeric(out) == s) break
  }
  out
}

# ---- Argument checks ----
# TRUE when x is one non-missing number from lower to upper, and whole if asked.
is_scalar <- function(x, lower = -Inf, upper = Inf, whole = FALSE) {
  ok <- is.numeric(x) && length(x) == 1L && !is.na(x) && x >= lower && x <= upper
  ok && (!whole || x == round(x))
}

# Rows named in a refusal: the count, then at most `max_shown` row numbers,
# so a table of hundreds of rows gives a message of one line.
rows_text <- function(rows, max_shown = 10L) {
  shown <- paste(utils::head(rows, max_shown), collapse = ", ")
  more <- length(rows) - max_shown
  if (more > 0) shown <- paste0(shown, " and ", more, " more")
  list(count = paste(length(rows), if (length(rows) == 1L) "row" else "rows"), shown = shown)
}

stop_unless <- function(ok, ..., kind = "input") {
  if (!isTRUE(ok)) stop_classed(kind, ...)
  invisible()
}

# ---- Errors ----
# Every error arraydesign raises with stop_classed() has the classes
# `arraydesign_error_<kind>`, `arraydesign_error`, `error` and `condition`, so
# a caller can catch one kind by class with tryCatch(). The kinds: `input`, an
# argument of the wrong type or out of range; `grid`, a raster that cannot be
# read as water or whose CRS is not in metres on the ground; `placement`,
# receivers that cannot be placed; `track`, a
# start or goal the walker cannot use; `field`, a detection field that does
# not match the tracks or the receivers. `call` is the user's call that the
# error names; NULL names none.
stop_classed <- function(kind, ..., call = NULL) {
  cond <- structure(
    class = c(paste0("arraydesign_error_", kind), "arraydesign_error", "error", "condition"),
    # collapse = "" keeps the message one string when a piece is a vector.
    list(message = paste0(..., collapse = ""), call = call)
  )
  stop(cond)
}

# match.arg() with the same matching, whose refusal is an error of kind
# `input`. A value of the wrong type or length gets its own message, since
# "should be one of" would point at the spelling. `call` is the user's call
# that the refusal names; NULL names none.
match_choice <- function(arg, choices, name, call = NULL) {
  if (identical(arg, choices)) return(choices[1])
  one_of <- paste0("\"", choices, "\"", collapse = ", ")
  stop_unless(is.character(arg) && length(arg) == 1 && !is.na(arg),
              "`", name, "` must be one string, one of ", one_of, ".", call = call)
  tryCatch(match.arg(arg, choices), error = function(e) {
    stop_classed("input", "`", name, "` should be one of ", one_of, ".", call = call)
  })
}

# The x and y coordinates as a two-column numeric matrix. A missing value and
# a value that does not read as a number get different messages; the second
# names the first such value and its row.
numeric_xy <- function(x, y, what, columns = c("x", "y")) {
  xy <- cbind(x = as_number(x, what, columns[1]), y = as_number(y, what, columns[2]))
  bad <- which(is.na(xy), arr.ind = TRUE)
  stop_unless(!nrow(bad), "`", what, "` holds missing coordinates, the first in row ",
              min(bad[, 1]), ".")
  xy
}

as_number <- function(v, what, column) {
  if (is.numeric(v)) return(as.numeric(v))
  out <- suppressWarnings(as.numeric(as.character(v)))
  not_number <- which(is.na(out) & !is.na(v))
  stop_unless(!length(not_number), "`", what, "$", column, "` holds values that are not ",
              "numbers: ", length(not_number), " of ", length(v), ", the first \"",
              as.character(v[not_number[1]]), "\" in row ", not_number[1], ".")
  out
}

# ---- Console text ----
# Text wrapped to the console width, so prints and messages fit a narrow
# console or a knitted page. `indent` is the first line's indent, `exdent`
# the indent of the lines after it.
wrap_text <- function(text, indent = 0, exdent = indent + 2) {
  paste(strwrap(text, width = getOption("width", 80), indent = indent, exdent = exdent),
        collapse = "\n")
}

cat_wrapped <- function(text, indent = 0, exdent = indent + 2) {
  cat(wrap_text(text, indent, exdent), "\n", sep = "")
}
