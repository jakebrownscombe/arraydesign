#' Simulate detections of transmissions by a receiver array
#'
#' Draws, for every transmission and every receiver in a detection field, one
#' Bernoulli trial with the detection probability the field gives that
#' receiver at the cell holding the transmission. The field is a `de_field`
#' from [detfield::de_field()]: `eta`, a matrix of cells by receivers on the
#' logit scale; `cells`, the water-cell ids of its rows; `receivers`, the
#' receiver keys of its columns; `key`, the column the keys come from; and
#' `grid`, the raster template the cell ids refer to. A list of the same shape
#' with no `key` is read as keyed by `receiver_serial`.
#'
#' **Receivers and stations.** A field built from a table with a `station_id`
#' column has stations as its columns. The detections then carry `station_id`,
#' and `receiver_serial` is the receiver deployed at that station at the time
#' of the detection, read from the field's `lookup` (the deployments given to
#' [detfield::water_distances()]). Where two deployments of a station overlap,
#' the receiver deployed later takes the detection, as in
#' [detfield::bin_detections()]. With no `lookup`, `receiver_serial` is `NA`,
#' and [detfield::bin_detections()] matches the detections by station to
#' listening rows that have no serial.
#'
#' **Stations not listening.** When the `lookup` has deployment times, a
#' station makes a trial only when it listens, by the rule of
#' [detfield::bin_detections()] with hour bins: a receiver is in the water at
#' the station at the time of the transmission, and the station's deployments
#' together cover at least `min_cover` of that UTC hour. A station deployed
#' late, recovered early or between receivers therefore makes no detection at
#' those times, and the detections that remain all fall in listening rows of
#' [detfield::bin_detections()] given the same deployments and `min_cover`.
#' A message gives the number of trials not made. A `lookup` with no times
#' takes every station as listening throughout.
#'
#' **Overdispersion.** With `overdispersion = "station_hour"`, the default,
#' each receiver and UTC hour gets one logit-scale effect, normal with standard
#' deviation `field$model$sigma$station_hour` and shared by every animal in the
#' hour: the form [detfield::location_surface()] integrates. When the field's
#' model has no such value (`NA`, absent, or 0), or with `"none"`, the trials
#' are independent binomial draws. Left at its default with a spread above 0,
#' the call says so in a message.
#'
#' The function stops when a transmission falls on a cell the field does not
#' hold, or when any detection probability it needs is missing. It never
#' fills a gap with a guessed probability, because a guess would be scored as
#' if the range model had predicted it.
#'
#' @param transmissions A `sim_transmissions` object from
#'   [simulate_transmissions()].
#' @param field A `de_field`: see Details.
#' @param overdispersion `"station_hour"` (default) or `"none"`; see Details.
#' @param min_cover Share of a UTC hour, from 0 to 1, that a station's
#'   deployments must cover for it to listen in the hour; the default is the
#'   default of [detfield::bin_detections()]. Used only when the field's
#'   `lookup` has deployment times; see Details.
#' @param seed Seed; `NULL` uses the current random number stream. The
#'   caller's stream is restored after a seeded call. The draws depend on the
#'   generator kind as well as the seed; the kind is stored with the result.
#'   Without overdispersion, transmission `i` always takes the `i`-th run of
#'   uniform numbers, one per receiver, so the first `m` transmissions get the
#'   same draws whether they are simulated alone or with more. With it, the
#'   hourly effects are drawn first, one run per hour in time order, and the
#'   uniform numbers follow, hour by hour. Trials of a station not listening
#'   are drawn and then discarded, so they do not change the other draws.
#'
#' @return A list of class `sim_detections` with two data frames:
#'   `detections`, one row per detection, with the columns of the positionR
#'   table contract (`datetime`, `receiver_serial`, `transmitter`, `animal_id`,
#'   and `station_id` when the field is keyed by station) and
#'   `transmission_id`, ordered by time and receiver key; and `transmissions`,
#'   the truth: every transmission with its position, heard or not, and
#'   `n_detections`, the receivers that heard it. Attributes: `seed`,
#'   `rng_kind` (the [RNGkind()] in force), `overdispersion`,
#'   `station_hour`, the standard deviation used (0 for none), and
#'   `trials_not_listening`, the transmission and station trials not made
#'   because the station was not listening (0 when none were left out).
#' @seealso [simulate_transmissions()], [coverage_summary()]
#' @export
#' @examples
#' stoney <- terra::rast(system.file("extdata", "stoney_depth.tif",
#'                                   package = "arraydesign"))
#' planned <- read.csv(system.file("extdata", "planned_stations.csv",
#'                                 package = "arraydesign"))
#' rx <- place_receivers(stoney, method = "manual", coords = planned)
#' # A deployments table: one receiver per station for two days.
#' deployments <- data.frame(
#'   receiver_serial = rx$receiver_serial, station_id = rx$station_id,
#'   x = rx$x, y = rx$y, deploy_z = rx$deploy_z,
#'   deploy_datetime = as.POSIXct("2024-12-31", tz = "UTC"),
#'   recover_datetime = as.POSIXct("2025-01-02", tz = "UTC")
#' )
#' range_model <- detfield::de_logistic(d50 = c(250, 750), d05 = c(500, 1500),
#'                                      depth = c(1, 35))
#' # Pairs more than 4 km apart through water are left out. At the deepest
#' # receiver the curve gives them a detection probability below 1e-7, and
#' # the distances take a third of the time.
#' field <- detfield::de_field(detfield::water_distances(stoney, deployments,
#'                                                     max_distance = 4000),
#'                             range_model)
#'
#' # Five tracks of 2 hours.
#' tr <- simulate_tracks(stoney, crw_model(30, 15, 40, time_step_s = 60),
#'                       n_animals = 5, n_steps = 120, seed = 1)
#' tx <- simulate_transmissions(tr, tag_delay = c(120, 180), seed = 2)
#' # This range model has no station-hour effect. Each detection probability
#' # comes from the range curve alone, with no random effect per receiver or
#' # hour.
#' det <- simulate_detections(tx, field, overdispersion = "none", seed = 3)
#' head(det$detections)
#' table(det$transmissions$n_detections)
#'
#' # The detections go straight into detfield.
#' counts <- detfield::bin_detections(det$detections, deployments,
#'                                    tag_delay = c(120, 180))
#' table(counts$listening)
simulate_detections <- function(transmissions, field,
                                overdispersion = c("station_hour", "none"), min_cover = 0.5,
                                seed = NULL) {
  if (!inherits(transmissions, "sim_transmissions")) {
    stop_classed("input", "`transmissions` must come from simulate_transmissions().")
  }
  fg <- check_field(field)
  same_grid(water_grid(attr(transmissions, "grid")), fg)
  sigma <- detection_spread(field, overdispersion, missing(overdispersion))
  stop_unless(is.numeric(min_cover) && length(min_cover) == 1L && isTRUE(min_cover >= 0) &&
                isTRUE(min_cover <= 1), "`min_cover` must be one number from 0 to 1.")
  row <- match(transmissions$cell, field$cells)
  if (anyNA(row)) {
    stop_classed("field", sum(is.na(row)), " of ", nrow(transmissions),
                 " transmissions fall on cells the field does not hold.")
  }
  eta <- field$eta[row, , drop = FALSE]
  if (anyNA(eta)) {
    bad <- which(is.na(eta), arr.ind = TRUE)
    stop_classed("field", "The field gives no detection probability for ", nrow(bad),
                 " transmission and receiver pairs (", length(unique(bad[, 2])), " receivers, ",
                 length(unique(bad[, 1])), " transmissions). Fix the field; no probability is ",
                 "assumed.")
  }
  hits <- with_seed_or_not(seed, draw_hits(eta, transmissions$datetime, sigma))
  trials <- station_trials(field, transmissions$datetime, min_cover)
  if (!is.null(trials$row)) hits <- hits[!is.na(trials$row[hits]), , drop = FALSE]
  det <- detection_rows(transmissions, field, hits, trials)
  truth <- as.data.frame(transmissions)
  class(truth) <- "data.frame"
  truth$n_detections <- tabulate(hits[, 1], nbins = nrow(truth))
  structure(list(detections = det, transmissions = truth), class = "sim_detections",
            seed = seed, rng_kind = RNGkind(),
            overdispersion = if (sigma > 0) "station_hour" else "none", station_hour = sigma,
            trials_not_listening = if (is.null(trials)) 0L else trials$masked)
}

# ---- Overdispersion ----
# The station-hour standard deviation to draw with: the field model's value,
# or 0 for independent binomial draws.
detection_spread <- function(field, overdispersion, by_default) {
  overdispersion <- match_choice(overdispersion, c("station_hour", "none"), "overdispersion")
  if (overdispersion == "none") return(0)
  sigma <- field$model$sigma$station_hour
  if (is.null(sigma) || length(sigma) != 1L || is.na(sigma)) return(0)
  if (sigma > 0 && by_default) {
    message(wrap_text(paste0(
      "Detections share one logit effect per receiver and UTC hour, standard deviation ",
      format(sigma), " (the field's sigma$station_hour); `overdispersion = \"none\"` ",
      "leaves it out."
    ), exdent = 0))
  }
  as.numeric(sigma)
}

# Rows (transmissions) and columns (receivers) of the successes. With a
# spread, one effect per hour and receiver is drawn first, by hour in time
# order, and the trials follow hour by hour, so all animals in an hour share
# the hour's effects.
draw_hits <- function(eta, datetime, sigma) {
  if (sigma == 0) return(bernoulli_hits(stats::plogis(eta)))
  hour <- floor(as.numeric(datetime) / 3600)
  hours <- sort(unique(hour))
  eff <- matrix(stats::rnorm(length(hours) * ncol(eta), 0, sigma), length(hours), ncol(eta),
                byrow = TRUE)
  at <- match(hour, hours)
  ord <- order(at)
  hits <- bernoulli_hits(stats::plogis(eta[ord, , drop = FALSE] + eff[at[ord], , drop = FALSE]))
  hits[, 1] <- ord[hits[, 1]]
  hits
}

# ---- Detection table ----
detection_rows <- function(tx, field, hits, trials) {
  keys <- as.character(field$receivers)[hits[, 2]]
  det <- data.frame(
    datetime = tx$datetime[hits[, 1]],
    receiver_serial = keys,
    transmitter = tx$transmitter[hits[, 1]],
    animal_id = tx$animal_id[hits[, 1]],
    transmission_id = tx$transmission_id[hits[, 1]],
    stringsAsFactors = FALSE
  )
  if (identical(field$key, "station_id")) {
    at <- if (is.null(trials)) {
      rep(NA_integer_, nrow(hits))
    } else if (!is.null(trials$row)) {
      trials$row[hits]
    } else {
      trials$station_row[hits[, 2]]
    }
    det$receiver_serial <- if (is.null(trials)) {
      rep(NA_character_, nrow(hits))
    } else {
      as.character(field$lookup$receiver_serial)[at]
    }
    det$station_id <- keys
    det <- det[c("datetime", "receiver_serial", "station_id", "transmitter", "animal_id",
                 "transmission_id")]
  }
  det <- det[order(det$datetime, keys), ]
  rownames(det) <- NULL
  det
}

# ---- Stations and their receivers ----
# For a field keyed by station with a lookup: the lookup row of the receiver
# that makes each trial. With deployment times, `row` is a matrix of
# transmissions by stations, NA where the station is not listening; without
# them, `station_row` gives one row per station and every station listens.
# NULL for other fields: bin_detections() then matches by station.
station_trials <- function(field, datetime, min_cover) {
  lookup <- field$lookup
  if (!identical(field$key, "station_id") || is.null(lookup)) return(NULL)
  stations <- as.character(field$receivers)
  if (!all(c("deploy_datetime", "recover_datetime") %in% names(lookup))) {
    n <- tabulate(match(lookup$station_id, stations), nbins = length(stations))
    many <- stations[n > 1L]
    stop_unless(!length(many), "Station ", many[1], " has ", max(n), " receivers in the ",
                "field's lookup and no deployment times, so a detection there has no one ",
                "receiver. Give water_distances() deploy_datetime and recover_datetime.")
    return(list(station_row = match(stations, lookup$station_id), masked = 0L))
  }
  t <- as.numeric(datetime)
  row <- matrix(NA_integer_, length(t), length(stations))
  absent <- 0L
  short <- 0L
  for (j in seq_along(stations)) {
    at <- which(lookup$station_id %in% stations[j])
    at <- at[order(as.numeric(lookup$deploy_datetime[at]))]
    on <- trial_deployment(lookup[at, , drop = FALSE], t, min_cover)
    row[, j] <- at[on$row]
    absent <- absent + sum(is.na(on$deployed))
    short <- short + sum(!is.na(on$deployed) & is.na(on$row))
  }
  masked <- absent + short
  if (masked > 0L) {
    off <- stations[colSums(is.na(row)) > 0L]
    message(wrap_text(paste0(
      masked, " of ", length(row), " transmission and station trials were not made ",
      "because the station was not listening: ", absent, " with no receiver deployed ",
      "at the time, and ", short, " in an hour the station's deployments cover for less ",
      "than `min_cover` (", format(min_cover), "). Stations: ", paste(off, collapse = ", "), "."
    ), exdent = 0))
  }
  list(row = row, masked = masked)
}

# The trials one station makes at times `t` (seconds), from its deployments
# ordered by deploy time. `deployed` is the deployment in the water at each
# time (the later one where two overlap, as in bin_detections()), and `row`
# is the same where the deployments also cover min_cover of the UTC hour,
# else NA. Cover sums over the deployments, at most 1, as in bin_detections().
trial_deployment <- function(dep, t, min_cover) {
  none <- rep(NA_integer_, length(t))
  if (!nrow(dep) || !length(t)) return(list(deployed = none, row = none))
  start <- as.numeric(dep$deploy_datetime)
  end <- as.numeric(dep$recover_datetime)
  end[is.na(end)] <- Inf
  inside <- outer(t, start, ">=") & outer(t, end, "<")
  deployed <- max.col(inside * rep(seq_len(nrow(dep)), each = length(t)), ties.method = "last")
  deployed[rowSums(inside) == 0] <- NA_integer_
  hour <- floor(t / 3600) * 3600
  hours <- unique(hour)
  share <- vapply(hours, function(h) {
    min(1, sum(pmax(0, pmin(end, h + 3600) - pmax(start, h))) / 3600)
  }, numeric(1))
  listens <- share > 0 & share >= min_cover
  row <- deployed
  row[!listens[match(hour, hours)]] <- NA_integer_
  list(deployed = deployed, row = row)
}

# Rows and columns of the successes of one Bernoulli trial per cell of p.
# Draws are made in blocks of rows so a long simulation does not hold every
# uniform number at once. Each block is filled row by row, so row i takes
# uniforms (i - 1) * ncol + 1 to i * ncol of the stream whatever the row count.
bernoulli_hits <- function(p) {
  block <- max(1L, floor(5e6 / ncol(p)))
  starts <- seq(1L, nrow(p), by = block)
  hits <- lapply(starts, function(s) {
    rows <- s:min(s + block - 1L, nrow(p))
    u <- matrix(stats::runif(length(rows) * ncol(p)), length(rows), byrow = TRUE)
    h <- which(u < p[rows, , drop = FALSE], arr.ind = TRUE)
    h[, 1] <- h[, 1] + s - 1L
    h
  })
  do.call(rbind, c(list(matrix(integer(), 0, 2)), hits))
}

# ---- Detection field checks ----
# A de_field from detfield, or a list of the same shape.
check_field <- function(field) {
  need <- c("eta", "cells", "receivers", "grid")
  if (!is.list(field) || !all(need %in% names(field))) {
    stop_classed("input", "`field` must be a de_field, or a list with elements ",
                 paste0("`", need, "`", collapse = ", "), ".")
  }
  if (!is.matrix(field$eta) || !is.numeric(field$eta)) {
    stop_classed("input", "`field$eta` must be a numeric matrix of cells by receivers.")
  }
  if (nrow(field$eta) != length(field$cells) || ncol(field$eta) != length(field$receivers)) {
    stop_classed("input", "`field$eta` is ", nrow(field$eta), " by ", ncol(field$eta),
                 " but the field has ", length(field$cells), " cells and ", length(field$receivers),
                 " receivers.")
  }
  if (anyDuplicated(field$cells) || anyDuplicated(field$receivers)) {
    stop_classed("input", "`field$cells` and `field$receivers` must not repeat.")
  }
  water_grid(field$grid)
}

# Cell ids mean the same place only when the two grids share geometry and CRS.
same_grid <- function(a, b) {
  crs_a <- terra::crs(a$raster)
  crs_b <- terra::crs(b$raster)
  if (nzchar(crs_a) && nzchar(crs_b) && !terra::same.crs(crs_a, crs_b)) {
    stop_classed("field",
                 "The tracks and the field are in different CRSs, so their cell ids do not mean ",
                 "the same places.")
  }
  ga <- unlist(a[c("res", "xmin", "ymax", "nrow", "ncol")])
  gb <- unlist(b[c("res", "xmin", "ymax", "nrow", "ncol")])
  if (any(abs(ga - gb) > 1e-6 * pmax(1, abs(ga)))) {
    stop_classed("field",
                 "The tracks and the field are on different grids, so their cell ids do not match.")
  }
  invisible()
}
