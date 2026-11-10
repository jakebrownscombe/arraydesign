#' Simulate tag transmissions along tracks
#'
#' Places the transmissions of a tag with a random delay along each simulated
#' track. `tag_delay` is defined as in [detfield::bin_detections()]: the time
#' from the start of one transmission to the start of the next, the
#' transmission itself included. Each interval is drawn uniformly between the
#' shortest and the longest delay, as in the transmission model of
#' `glatos::transmit_along_path()`, and the first transmission falls at a
#' uniform random time within the first interval. The transmissions do not
#' depend on the step of the walk: a tag with a 120 s mean delay sends
#' 3600 / 120 = 30 transmissions an hour whatever the time step, the `n` that
#' `detfield::bin_detections()` gives for the same `tag_delay`. Each
#' transmission's position is interpolated along the straight step that holds
#' it, so it lies on water.
#'
#' A tag sheet usually gives the programmed delay without the transmission, a
#' few seconds: add it. The default, 65 to 185 s, is a tag programmed at 60 to
#' 180 s with a 5 s transmission, 28.8 transmissions an hour.
#'
#' @param tracks A `sim_tracks` object from [simulate_tracks()].
#' @param tag_delay Delay from the start of one transmission to the start of
#'   the next, in seconds: one number for a fixed delay, or two, the shortest
#'   and the longest, both above 0.
#' @param seed Seed; `NULL` uses the current random number stream. The
#'   caller's stream is restored after a seeded call. The draws depend on the
#'   generator kind as well as the seed; the kind is stored with the result.
#'
#' @return A data frame of class `sim_transmissions`, one row per
#'   transmission, ordered by path and time: `transmission_id` (integer),
#'   `path_id`, `animal_id` (the path as character), `transmitter`
#'   (`"SIM-0001"`, ...), `datetime` (POSIXct, UTC), `x`, `y`, `cell`,
#'   `state`, and `min_delay` and `max_delay`, the delay range used, in
#'   seconds, under the suite's column names (`transmitter`, `min_delay`,
#'   `max_delay`, as [detfield::bin_detections()] reads them), so
#'   `unique(tx[c("transmitter", "min_delay", "max_delay")])` is a `tag_delay`
#'   table for [detfield::bin_detections()]. Attributes: `tag_delay` (the two
#'   numbers), `seed`, `rng_kind` (the [RNGkind()] in force), and `epsg` and
#'   `grid`, taken from the tracks.
#' @seealso [simulate_tracks()], [simulate_detections()]
#' @export
#' @examples
#' stoney <- terra::rast(system.file("extdata", "stoney_depth.tif",
#'                                   package = "arraydesign"))
#' # Two animals for 2 hours.
#' walk <- crw_model(30, 15, 40, time_step_s = 60)
#' tr <- simulate_tracks(stoney, walk, n_animals = 2, n_steps = 120, seed = 1)
#' tx <- simulate_transmissions(tr, tag_delay = c(120, 180), seed = 2)
#' head(tx)
#' # Transmissions per hour of track: about 3600 / 150 = 24.
#' nrow(tx) / (2 * 120 * 60 / 3600)
simulate_transmissions <- function(tracks, tag_delay = c(65, 185), seed = NULL) {
  stop_unless(inherits(tracks, "sim_tracks"), "`tracks` must come from simulate_tracks().")
  tag_delay <- delay_range(tag_delay)
  wg <- water_grid(attr(tracks, "grid"))
  by_path <- split(tracks, tracks$path_id)
  tx <- with_seed_or_not(seed, lapply(by_path, transmit_one, tag_delay = tag_delay))
  out <- as.data.frame(data.table::rbindlist(tx))
  width <- max(4L, nchar(max(tracks$path_id)))
  out$transmission_id <- seq_len(nrow(out))
  out$animal_id <- as.character(out$path_id)
  out$transmitter <- sprintf("SIM-%0*d", width, out$path_id)
  out$cell <- as.integer(cell_at(wg, out$x, out$y))
  out$min_delay <- rep(tag_delay[1], nrow(out))
  out$max_delay <- rep(tag_delay[2], nrow(out))
  out <- out[c("transmission_id", "path_id", "animal_id", "transmitter", "datetime",
               "x", "y", "cell", "state", "min_delay", "max_delay")]
  structure(out, class = c("sim_transmissions", "data.frame"), tag_delay = tag_delay,
            seed = seed, rng_kind = RNGkind(), epsg = attr(tracks, "epsg"),
            grid = attr(tracks, "grid"))
}

# One number is a fixed delay; two are the shortest and longest. Both above
# 0, as detfield::bin_detections() reads them.
delay_range <- function(tag_delay) {
  ok <- is.numeric(tag_delay) && length(tag_delay) %in% 1:2 && all(is.finite(tag_delay)) &&
    all(tag_delay > 0)
  stop_unless(ok, "`tag_delay` must be one or two numbers of seconds, each above 0.")
  if (length(tag_delay) == 1L) tag_delay <- c(tag_delay, tag_delay)
  stop_unless(tag_delay[1] <= tag_delay[2], "`tag_delay` must give the shortest delay first.")
  as.numeric(tag_delay)
}

# Transmission times along one path, with positions on the steps that hold them.
transmit_one <- function(path, tag_delay) {
  t0 <- path$time[1]
  elapsed <- as.numeric(difftime(path$time, t0, units = "secs"))
  span <- max(elapsed)
  at <- cumsum(transmission_gaps(span, tag_delay))
  at <- at[at <= span]
  # The state of a transmission is the state of the step it falls on: the
  # step that ends at the first position later than it.
  k <- pmin(findInterval(at, elapsed) + 1L, nrow(path))
  data.frame(
    path_id = rep(path$path_id[1], length(at)),
    datetime = t0 + at,
    x = stats::approx(elapsed, path$x, xout = at)$y,
    y = stats::approx(elapsed, path$y, xout = at)$y,
    state = path$state[k]
  )
}

# Intervals between transmissions until they pass `span`, the first one cut
# to a uniform share of itself. Draws come in chunks sized by the mean
# interval, so a short smallest interval does not size a huge vector.
transmission_gaps <- function(span, tag_delay) {
  lo <- tag_delay[1]
  hi <- tag_delay[2]
  chunk <- function(left) stats::runif(ceiling(1.05 * left / ((lo + hi) / 2)) + 10L, lo, hi)
  gaps <- chunk(span)
  gaps[1] <- stats::runif(1, 0, gaps[1])
  while (sum(gaps) <= span) gaps <- c(gaps, chunk(span - sum(gaps)))
  gaps
}
