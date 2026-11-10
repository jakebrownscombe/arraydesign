#' Score a receiver array by simulation
#'
#' Scores a candidate array in one call. It builds the array's detection
#' field with [detfield::water_distances()] and [detfield::de_field()],
#' simulates animals with [simulate_tracks()], their transmissions with
#' [simulate_transmissions()] and the detections with
#' [simulate_detections()], counts the detections per hour with
#' [detfield::bin_detections()], computes a location surface for every
#' detected animal-hour with [detfield::location_surface()], and scores the
#' surfaces against the true positions with [score_surfaces()]. The share of
#' the water the array covers is returned beside the score.
#'
#' **What is scored.** Only animal-hours with at least one detection get a
#' surface, so the errors and the hit rate describe the detected hours. The
#' score's `overall` row adds `n_animal_hours`, the animal-hours that hold a
#' transmission, and `scored_share`, the scored bins over them. On the
#' 20-receiver grid of `vignette("arraydesign")` (20 animals x 6 hours,
#' seed 1), 77 of 120 animal-hours were scored; read the error of a sparse
#' array beside this share.
#'
#' **The detection model.** Detections are simulated from `model`, and
#' positions are computed with `control$position_model`, which is `model`
#' unless given; pass a [perturb_de_model()] copy to see what a wrong range
#' model does. When a model has `sigma$station_hour`, the simulation draws one
#' logit-scale effect per receiver and UTC hour, normal with that standard
#' deviation and shared by every animal in the hour, and the positioning
#' integrates over it (`overdispersion = "station_hour"`). Without it the
#' simulation has no such effect and the positioning uses independent
#' binomials (`overdispersion = "none"`). With the same model on both sides
#' the 95% regions are calibrated by construction, apart from the movement
#' of an animal within the hour, which the likelihood does not model, and
#' the spread of the random delay around its mean, since the counts take the
#' mean number of transmissions. Both lower the hit rate. On the same
#' 20-receiver grid with the binomial and a delay of 120 to 180 s, 200 animals over
#' 24 hours gave 0.958 with still animals, 0.886 with a walk of 6 m a
#' minute, and 0.546 with the default movement, 30 m a minute, which carries
#' an animal well outside a 95% region within the hour. The script,
#' <https://github.com/jakebrownscombe/arraydesign/blob/main/validation/03_stoney_details.R>,
#' is in the repository and not in the installed package.
#'
#' **The hit rate measures movement within the hour.** On Stoney, with 60
#' animals x 12 hours under the binomial, grids of 5, 10 and 20 receivers
#' gave 0.499, 0.550 and 0.561 with the default movement, and 0.979, 0.941
#' and 0.966 with still animals. One receiver scored 7 of 720 animal-hours,
#' too few to read. Compare designs on the error, the scored share and the
#' water covered. Use the hit rate to check that the
#' regions are calibrated, with still animals (`crw_model(1e-6, 1e-6, 40)`)
#' when the design itself is the question.
#'
#' **The hit-rate interval needs hours.** Its SE is clustered two ways, by
#' animal and by hour (bin start), and is never taken below the larger
#' one-way SE; [score_surfaces()] gives the formula, and the print names the
#' estimate used. The interval's degrees of freedom are the smaller cluster
#' count less one, so the default six hours give 5. On the same grid the 95%
#' interval left out the mean of 40 replicate runs at six hours in 0 of 40
#' under the binomial and 1 of 40 at a station-hour standard deviation of 1,
#' and the mean of 20 runs at 24 hours
#' in 0 of 20 at 1. The print notes when a station-hour effect meets fewer
#' than 20 hours. Use 24 hours or more for a calibration check.
#'
#' **Transmissions per hour.** The simulation and the counts read the same
#' `control$tag_delay`, the start-to-start delay of [simulate_transmissions()]
#' and [detfield::bin_detections()], so `n` matches the simulated
#' transmissions on average.
#'
#' **Coverage.** With `control$k = 1`, [coverage_summary()]: water cells
#' where at least one receiver on its own detects a transmission with
#' probability at or above each threshold. The print adds the cells where the
#' array does, any receiver hearing: one less the product over the receivers
#' of the probability of a miss, the coverage of Pedersen et al. (2014). With
#' `k` above 1, cells where at least `k` receivers do, the case for
#' positioning.
#' Coverage is read from the field of `model`, the simulated truth.
#'
#' @param layout The receivers: a `receiver_layout` from [place_receivers()],
#'   or a data frame of coordinates that `place_receivers(method = "manual")`
#'   takes. A layout missing any of its columns is placed as coordinates.
#' @param grid A [terra::SpatRaster] of water depth in metres, `NA` on land,
#'   in a projected CRS whose unit is the metre and whose scale at the raster
#'   centre, on its two axes and its diagonals, is within 1 percent of 1, as
#'   in [place_receivers()].
#' @param model A `de_model` from [detfield::de_model()] or
#'   [detfield::de_logistic()]: the detection model of the simulation.
#' @param n_animals Number of simulated animals.
#' @param seed A seed, one number from -2147483647 to 2147483647 as
#'   [set.seed()] takes, or `NULL` (the default) for the current random number
#'   stream. The seed is passed to [simulate_tracks()],
#'   [simulate_transmissions()] and [simulate_detections()] in turn, as the
#'   arraydesignR app does, so the same seed gives the app's animals and
#'   detections. The caller's stream is restored after a seeded call. A
#'   fraction is truncated toward zero, as [set.seed()] does, so `seed = 1.5`
#'   gives the result of `seed = 1`.
#' @param control A list of settings; any not given take their defaults:
#'   * `movement`: a [movement_model()]; default
#'     `crw_model(30, 15, 40, time_step_s = 60)`, a correlated random walk with
#'     gamma steps of mean 30 m and standard deviation 15 m per 60 s (0.5 m/s)
#'     and turns of standard deviation 40 degrees, the walk of
#'     `vignette("movement-models")`. No species was measured to set these
#'     values; set a walk for your species.
#'   * `hours`: hours simulated per animal, from the hour; default 6. With the
#'     20-receiver grid of `vignette("arraydesign")` on the 100 m raster, and
#'     the range curve and tag delay of the example below, 20 animals over 6
#'     hours took a median of 3.3 s under the binomial and 11.2 s at a
#'     station-hour standard deviation of 2.55 (10 runs each). The time
#'     depends on the machine and its load. Under the binomial most of
#'     the time goes to the distances through water of the detection field;
#'     with a station-hour standard deviation, to positioning, in proportion
#'     to the animal-hours scored.
#'   * `k`: number of receivers a cell needs at a threshold to count as
#'     covered; default 1. A `k` above the number of receivers warns, and
#'     every count is 0.
#'   * `tag_delay`: the delay from the start of one transmission to the start
#'     of the next, in seconds, one number or the shortest and longest,
#'     passed to [simulate_transmissions()] and [detfield::bin_detections()];
#'     default `c(65, 185)`, a tag programmed at 60 to 180 s with a 5 s
#'     transmission.
#'   * `level`: level of the credible regions; default 0.95.
#'   * `tag_model`: the tag model to read from `model` and `position_model`
#'     when they hold several, one string passed to [detfield::de_field()];
#'     default `NULL`, for a model with one curve.
#'   * `position_model`: the `de_model` used to position; default `model`.
#'     A model that gives every cell the same likelihood, such as one flat in
#'     distance, makes every bin degenerate: the call warns and scores none.
#'   * `thresholds`: detection probabilities for the coverage table; default
#'     `c(0.5, 0.8)`.
#'   * `start_time`: POSIXct in UTC on the hour; default 2025-01-01 00:00.
#'   * `cores`: processes for [detfield::location_surface()]; default 1.
#'   * `likelihood`: how [detfield::location_surface()] reads the curve of
#'     `position_model`. `"binomial"` (the default), binomial grid
#'     localization: the curve is the probability that each transmission is
#'     heard, so misses and silent receivers count. `"presence_only"`,
#'     presence-only grid localization: the product of the curve over the
#'     detections, with misses and silent receivers dropped. Its credible
#'     regions are not calibrated; they sharpen with the count of detections.
#'
#' @return A list of class `array_evaluation`: `layout`, the receivers;
#'   `coverage`, the coverage table with `k`; `score`, a `surface_score`
#'   whose `overall` row also holds `n_animal_hours` and `scored_share`;
#'   `field`, the `de_field` of `model`; `tracks`; `detections`, the
#'   `sim_detections`; `surfaces`, the location surfaces (`NULL` when nothing
#'   was detected); `settings`, the control values used plus `n_animals`,
#'   `seed`, the model sources and `skipped`, the degenerate bins left out; and
#'   `timing`, seconds per step.
#' @references
#' Pedersen, M. W., Burgess, G. and Weng, K. C. (2014). A quantitative
#' approach to static sensor network design. *Methods in Ecology and
#' Evolution*, 5, 1043-1051. \doi{10.1111/2041-210X.12255}
#' (eqn 5: a cell's coverage is the chance that one or more receivers detect a
#' transmission from it.)
#' @seealso [score_surfaces()], [perturb_de_model()], [coverage_summary()]
#' @export
#' @examples
#' # The 20 planned stations, the depth-dependent range curve, and 20 animals
#' # for 6 hours, as in the app's reference run. The depth raster is averaged
#' # to 200 m cells (1,345 water cells instead of 4,818), so the example runs
#' # in about a second; `hours` under `control` gives the time on 100 m cells.
#' stoney <- terra::rast(system.file("extdata", "stoney_depth.tif",
#'                                   package = "arraydesign"))
#' stoney_200 <- terra::aggregate(stoney, 2, na.rm = TRUE)
#' planned <- read.csv(system.file("extdata", "planned_stations.csv",
#'                                 package = "arraydesign"))
#' rx <- place_receivers(stoney_200, method = "manual", coords = planned)
#' range_model <- detfield::de_logistic(d50 = c(250, 750), d05 = c(500, 1500),
#'                                      depth = c(1, 35))
#' ev <- evaluate_array(rx, stoney_200, range_model, n_animals = 20, seed = 1,
#'                      control = list(hours = 6, tag_delay = c(120, 180)))
#' ev
#' plot(ev)
#'
#' # Positioning with a shorter range than the truth: see ?perturb_de_model.
evaluate_array <- function(layout, grid, model, n_animals = 20, seed = NULL,
                           control = list()) {
  stop_unless(inherits(model, "de_model"),
              "`model` must be a de_model from detfield::de_model() or de_logistic().")
  stop_unless(is_scalar(n_animals, 1, whole = TRUE),
              "`n_animals` must be one whole number, 1 or more.")
  check_seed(seed)
  ctrl <- evaluation_control(control, model)
  ctrl["seed"] <- list(seed)
  k <- ctrl$k
  timing <- c(field = NA_real_, simulate = NA_real_, position = NA_real_, score = NA_real_)
  t0 <- proc.time()[["elapsed"]]
  call <- sys.call()
  water_grid(grid, require_crs = TRUE, call = call)
  receivers <- layout_receivers(layout, grid, call)
  check_k(k, nrow(receivers))
  dist <- detfield::water_distances(grid, receivers[c("receiver_serial", "x", "y", "deploy_z")])
  # A range model fitted over several tag models needs one named.
  field <- detfield::de_field(dist, model, tag_model = ctrl$tag_model)
  pos_field <- if (identical(ctrl$position_model, model)) field else
    detfield::de_field(dist, ctrl$position_model, tag_model = ctrl$tag_model)
  timing[["field"]] <- since(t0)
  t0 <- proc.time()[["elapsed"]]
  sims <- simulate_array(grid, field, n_animals, ctrl)
  timing[["simulate"]] <- since(t0)
  t0 <- proc.time()[["elapsed"]]
  pos <- position_animals(sims$detections, receivers, pos_field, ctrl)
  timing[["position"]] <- since(t0)
  t0 <- proc.time()[["elapsed"]]
  score <- if (is.null(pos$surfaces)) empty_score(ctrl$level) else
    score_surfaces(pos$surfaces, sims$detections$transmissions, ctrl$level)
  score <- add_denominator(score, sims$detections$transmissions)
  coverage <- array_coverage(field, k, ctrl$thresholds)
  timing[["score"]] <- since(t0)
  settings <- c(ctrl[setdiff(names(ctrl), "position_model")],
                list(n_animals = n_animals, model = model$source,
                     position_model = ctrl$position_model$source, skipped = pos$skipped))
  structure(list(layout = receivers, coverage = coverage, score = score, field = field,
                 tracks = sims$tracks, detections = sims$detections, surfaces = pos$surfaces,
                 settings = settings, timing = timing),
            class = "array_evaluation")
}

since <- function(t0) proc.time()[["elapsed"]] - t0

# A seed in `control` always stops the call. When the argument takes the value,
# the message advises its whole part, every digit printed: set.seed() truncates
# a fraction, so the advised seed gives the same stream. A fraction printed to
# 15 significant digits can round up to the next whole number. Any other value
# is described, and the advice is seed = 1.
seed_advice <- function(s) {
  if (is_seed(s)) {
    return(paste0("Give it as evaluate_array(..., seed = ", sprintf("%.0f", trunc(s) + 0), ")."))
  }
  paste0("Give one whole number, such as evaluate_array(..., seed = 1); `control$seed` holds ",
         describe_seed(s), ".")
}

# ---- Settings ----
evaluation_defaults <- function() {
  list(movement = crw_model(30, 15, 40, time_step_s = 60), hours = 6, tag_delay = c(65, 185),
       level = 0.95, position_model = NULL, tag_model = NULL, thresholds = c(0.5, 0.8),
       k = 1, start_time = as.POSIXct("2025-01-01 00:00:00", tz = "UTC"), cores = 1,
       likelihood = "binomial")
}

evaluation_control <- function(control, model) {
  stop_unless(is.list(control), "`control` must be a list.")
  stop_unless(!length(control) || (!is.null(names(control)) && all(nzchar(names(control)))),
              "Every element of `control` needs a name, such as `list(hours = 4)`.")
  # Building a movement model draws from its sampler once as a check, so the
  # defaults are built without moving the caller's random number stream.
  defaults <- withr::with_preserve_seed(evaluation_defaults())
  stop_unless(!"seed" %in% names(control), "`seed` is an argument of evaluate_array(). ",
              seed_advice(control$seed))
  unknown <- setdiff(names(control), names(defaults))
  stop_unless(!length(unknown), "Unknown setting in `control`: ", paste(unknown, collapse = ", "),
              ". Known: ", paste(names(defaults), collapse = ", "), ".")
  ctrl <- utils::modifyList(defaults, control, keep.null = TRUE)
  if (is.null(ctrl$position_model)) ctrl$position_model <- model
  stop_unless(inherits(ctrl$movement, "movement_model"), "`control$movement` must be a ",
              "movement_model.")
  stop_unless(inherits(ctrl$position_model, "de_model"), "`control$position_model` must be ",
              "a de_model.")
  stop_unless(is_scalar(ctrl$hours, 1, whole = TRUE), "`control$hours` must be one whole ",
              "number, 1 or more.")
  stop_unless(is_scalar(ctrl$k, 1, whole = TRUE), "`control$k` must be one whole number, ",
              "1 or more.")
  stop_unless(is_scalar(ctrl$level, 1e-9, 0.999), "`control$level` must be above 0 and at ",
              "most 0.999, the kept mass of the surfaces.")
  stop_unless(is.null(ctrl$tag_model) || (is.character(ctrl$tag_model) &&
                                            length(ctrl$tag_model) == 1 &&
                                            !is.na(ctrl$tag_model)),
              "`control$tag_model` must be NULL or one tag model name.")
  stop_unless(is_scalar(ctrl$cores, 1, whole = TRUE), "`control$cores` must be one whole ",
              "number, 1 or more.")
  # The readings detfield offers, taken from location_surface() itself.
  readings <- eval(formals(detfield::location_surface)$likelihood)
  stop_unless(is.character(ctrl$likelihood) && length(ctrl$likelihood) == 1L &&
                ctrl$likelihood %in% readings, "`control$likelihood` must be one of ",
              paste0("\"", readings, "\"", collapse = " or "), ".")
  ctrl$tag_delay <- delay_range(ctrl$tag_delay)
  check_utc(ctrl$start_time, "control$start_time")
  stop_unless(as.numeric(ctrl$start_time) %% 3600 == 0, "`control$start_time` must be on ",
              "the hour, so the simulated hours are the counted hours.")
  steps <- ctrl$hours * 3600 / ctrl$movement$time_step
  stop_unless(abs(steps - round(steps)) < 1e-9, "The movement model's time step must divide ",
              "an hour.")
  ctrl
}

# A whole receiver_layout as given; anything else, such as coordinates or a
# column subset of a layout, which keeps its class, is placed on the grid.
# `call` is evaluate_array()'s call, which a grid refusal names.
layout_receivers <- function(layout, grid, call) {
  need <- c("receiver_serial", "station_id", "x", "y", "deploy_z")
  if (inherits(layout, "receiver_layout") && all(need %in% names(layout))) return(layout)
  stop_unless(is.data.frame(layout), "`layout` must be a receiver_layout or a data frame of ",
              "coordinates.")
  place_layout(grid, NULL, NULL, "manual", as.data.frame(layout), NULL, call = call)
}

# A k above the receivers counts no cell; say so rather than return zeros
# that look like a poor array.
check_k <- function(k, n_rx) {
  if (k > n_rx) {
    warning("`control$k` is ", k, " but the layout has ", n_rx, " receiver", if (n_rx > 1) "s",
            ", so no cell can be covered by ", k, " receivers and every count is 0.",
            call. = FALSE)
  }
  invisible()
}

# ---- Simulation ----
# A seed seeds each of the three steps on its own, as the arraydesignR app
# does, so one seed gives the app's tracks, transmissions and detections.
simulate_array <- function(grid, field, n_animals, ctrl) {
  steps <- round(ctrl$hours * 3600 / ctrl$movement$time_step)
  tracks <- simulate_tracks(grid, ctrl$movement, n_animals = n_animals, n_steps = steps,
                            control = track_control(start_time = ctrl$start_time),
                            seed = ctrl$seed)
  tx <- simulate_transmissions(tracks, tag_delay = ctrl$tag_delay, seed = ctrl$seed)
  # The field's station-hour spread, when it has one, as location_surface()
  # integrates it. Without one the draws are those of overdispersion = "none".
  det <- simulate_detections(tx, field, overdispersion = "station_hour", seed = ctrl$seed)
  list(tracks = tracks, detections = det)
}

# ---- Positioning ----
# The surfaces, or NULL, and the number of degenerate bins left out. A
# positioning model that gives every cell the same likelihood makes every bin
# degenerate; that is reported as skipped bins with a warning.
position_animals <- function(sim, receivers, field, ctrl) {
  if (!nrow(sim$detections)) return(list(surfaces = NULL, skipped = 0L))
  span <- c(ctrl$start_time - 86400, ctrl$start_time + (ctrl$hours + 24) * 3600)
  listening <- data.frame(receiver_serial = receivers$receiver_serial,
                          station_id = receivers$station_id,
                          deploy_datetime = span[1], recover_datetime = span[2],
                          stringsAsFactors = FALSE)
  counts <- detfield::bin_detections(sim$detections, listening,
                                     tag_delay = ctrl$tag_delay, bin = "hour")
  sigma <- field$model$sigma$station_hour
  od <- if (is.null(sigma) || is.na(sigma) || sigma == 0) "none" else "station_hour"
  surfaces <- tryCatch(
    detfield::location_surface(counts, field, overdispersion = od,
                               likelihood = ctrl$likelihood,
                               control = list(cores = ctrl$cores, on_degenerate = "skip")),
    error = function(e) all_degenerate(e, counts)
  )
  if (is.null(surfaces)) {
    return(list(surfaces = NULL, skipped = nrow(unique(counts[c("animal_id", "bin_start")]))))
  }
  list(surfaces = surfaces, skipped = NROW(attr(surfaces, "skipped")))
}

# NULL, with a warning, when location_surface() stopped because every bin is
# degenerate, the condition of class `detfield_error_degenerate`; any other
# error is passed on.
all_degenerate <- function(e, counts) {
  if (!inherits(e, "detfield_error_degenerate")) stop(e)
  n <- nrow(unique(counts[c("animal_id", "bin_start")]))
  warning("All ", n, " detected animal-hours are degenerate: the positioning model gives ",
          "every cell the same likelihood, so none is scored. Does its curve fall with ",
          "distance?", call. = FALSE)
  NULL
}

empty_score <- function(level) {
  structure(list(bins = data.frame(), animals = data.frame(),
                 overall = data.frame(n_animals = 0L, n_bin_starts = 0L, n_bins = 0L)),
            class = "surface_score", level = level, no_truth = 0L)
}

# The scored share: bins scored over the animal-hours that hold a
# transmission. Hours with no detection get no surface, so the error and the
# hit rate describe the detected hours only.
add_denominator <- function(score, tx) {
  hours <- unique(data.frame(a = tx$animal_id, h = floor(as.numeric(tx$datetime) / 3600)))
  score$overall$n_animal_hours <- nrow(hours)
  score$overall$scored_share <- if (nrow(hours)) score$overall$n_bins / nrow(hours) else NA_real_
  score
}

# ---- Coverage by k receivers ----
array_coverage <- function(field, k, thresholds) {
  if (k == 1) return(cbind(coverage_summary(field, thresholds = thresholds), k = 1L))
  base <- coverage_summary(field, thresholds = thresholds)
  p <- stats::plogis(field$eta)
  base$n_cells <- vapply(thresholds, function(t) sum(rowSums(p >= t) >= k), integer(1))
  base$share <- base$n_cells / base$n_water_cells
  cbind(base, k = as.integer(k))
}

# ---- Methods ----
#' @export
print.array_evaluation <- function(x, ...) {
  s <- x$settings
  cat(sprintf("<array_evaluation> %d receiver%s; %d animals x %d hours\n", nrow(x$layout),
              if (nrow(x$layout) == 1) "" else "s", s$n_animals, s$hours))
  cat_wrapped(paste("simulated with:", s$model), indent = 2, exdent = 4)
  if (!identical(s$position_model, s$model)) {
    cat_wrapped(paste("positioned with:", s$position_model), indent = 2, exdent = 4)
  }
  cat("  positioned by:", sub("_", "-", s$likelihood), "grid localization\n")
  cv <- x$coverage
  # With k = 1, coverage_summary() counts cells that at least one receiver
  # hears at the threshold on its own; the array line counts cells that any
  # receiver hears, 1 - prod(1 - p) over the receivers.
  n_water <- format(cv$n_water_cells[1], big.mark = ",")
  if (cv$k[1] == 1) {
    cat_wrapped(sprintf("cells at p >= %s under one receiver alone: %s of %s (%.1f%%)",
                        format(cv$threshold), format(cv$n_cells, big.mark = ","), n_water,
                        100 * cv$share), indent = 2, exdent = 4)
    p_any <- 1 - apply(1 - x$field$p, 1, prod)
    n_any <- vapply(cv$threshold, function(t) sum(p_any >= t), integer(1))
    cat_wrapped(sprintf("cells at p >= %s under the array, any receiver: %s of %s (%.1f%%)",
                        format(cv$threshold), format(n_any, big.mark = ","), n_water,
                        100 * n_any / cv$n_water_cells[1]), indent = 2, exdent = 4)
  } else {
    cat_wrapped(sprintf("cells at p >= %s under %d+ receivers: %s of %s (%.1f%%)",
                        format(cv$threshold), cv$k, format(cv$n_cells, big.mark = ","),
                        n_water, 100 * cv$share), indent = 2, exdent = 4)
  }
  print(x$score)
  few_hours_note(x)
  if (s$skipped > 0) cat(sprintf("  %d degenerate bins left out\n", s$skipped))
  cat_wrapped(sprintf("%.1f s: %s", sum(x$timing),
                      paste(sprintf("%s %.1f", names(x$timing), x$timing), collapse = ", ")),
              indent = 2, exdent = 4)
  invisible(x)
}

# Under a station-hour effect the hit-rate SE rests on the hours, and fewer
# than 20 give an interval that is too narrow (validation/02_scoring_validation.R
# measures it).
few_hours_note <- function(x) {
  sigma <- x$field$model$sigma$station_hour
  has_sigma <- !is.null(sigma) && !is.na(sigma) && sigma > 0
  if (has_sigma && x$settings$hours < 20 && x$score$overall$n_bins > 0) {
    cat_wrapped(sprintf("with a station-hour effect and %d hours %s; see ?evaluate_array",
                        x$settings$hours, "the hit-rate interval is too narrow"),
                indent = 2, exdent = 4)
  }
  invisible()
}

#' Map the positioning error of an array evaluation
#'
#' Draws the water by depth on [detfield::scale_fill_depth()], land blank,
#' the receivers that heard the simulated animals as yellow rings and those
#' that heard none as white rings, each over a dark ring, and one point per scored
#' animal-hour at the centroid of the animal's true positions in the hour,
#' coloured by the error class of the location surface's mean (or mode) and
#' edged in black. The keys give the error classes and the depth in metres.
#'
#' The plot's data hold the scored animal-hours and the receivers, one row
#' each, told apart by the column `what`; `n_detections` holds the detections
#' of each receiver over the run. To map several evaluations of one
#' lake side by side, bind their data with a column that names each, replace
#' the data with [ggplot2::%+%], and facet by that column.
#'
#' @param x An `array_evaluation` from [evaluate_array()].
#' @param error `"mean"` or `"mode"`: which point estimate's error to draw.
#' @param breaks Class limits of the error in metres: at least 2 finite
#'   numbers in increasing order. The last class is open. When the first
#'   limit is above 0, the errors below it form a class of their own, so
#'   every animal-hour is drawn in the class that holds its error.
#' @param legend Whether to draw the keys of the error classes and of the
#'   depth to the right of the map.
#' @param ... Not used. Add ggplot2 layers, scales or [ggplot2::theme()]
#'   settings to the plot this method returns instead.
#'
#' @return A ggplot object, drawn when printed.
#' @seealso [evaluate_array()]
#' @export
#' @examples
#' # 20 animals for 6 hours on the 20 planned stations, with the depth raster
#' # averaged to 200 m cells to keep the example short.
#' stoney <- terra::rast(system.file("extdata", "stoney_depth.tif",
#'                                   package = "arraydesign"))
#' stoney_200 <- terra::aggregate(stoney, 2, na.rm = TRUE)
#' planned <- read.csv(system.file("extdata", "planned_stations.csv",
#'                                 package = "arraydesign"))
#' rx <- place_receivers(stoney_200, method = "manual", coords = planned)
#' range_model <- detfield::de_logistic(d50 = c(250, 750), d05 = c(500, 1500),
#'                                      depth = c(1, 35))
#' ev <- evaluate_array(rx, stoney_200, range_model, n_animals = 20, seed = 1,
#'                      control = list(hours = 6, tag_delay = c(120, 180)))
#' plot(ev)
#' # The error of the surface mode instead of the mean.
#' plot(ev, error = "mode")
plot.array_evaluation <- function(x, error = c("mean", "mode"), breaks = c(0, 200, 400, 800),
                                  legend = TRUE, ...) {
  call <- sys.call()
  error <- match_choice(error, c("mean", "mode"), "error", call = call)
  stop_unless(isTRUE(legend) || isFALSE(legend), "`legend` must be TRUE or FALSE.", call = call)
  stop_unless(is.numeric(breaks) && length(breaks) >= 2 && all(is.finite(breaks)) &&
                all(diff(breaks) > 0),
              "`breaks` must be at least 2 finite numbers in increasing order.",
              call = call)
  unused_dots(...)
  b <- x$score$bins
  # The classes run from -Inf to Inf, so every error falls in one. An error
  # is a distance, so the class below a first limit of 0 or less is empty
  # and is dropped from the key.
  labels <- c(sprintf("< %g", breaks[1]),
              sprintf("%g-%g", utils::head(breaks, -1), breaks[-1]),
              sprintf("> %g", breaks[length(breaks)]))
  # An evaluation that heard nothing has a bins table with no columns.
  if (!nrow(b)) b <- data.frame(truth_x = numeric(), truth_y = numeric(), error = numeric())
  e <- if (nrow(b)) b[[paste0("error_", error, "_m")]] else b$error
  cls <- findInterval(e, c(-Inf, breaks), left.open = FALSE)
  if (breaks[1] <= 0 && !any(cls == 1L, na.rm = TRUE)) {
    labels <- labels[-1]
    cls <- cls - 1L
  }
  # Detections of each receiver over the run, 0 for a receiver that heard none
  heard <- table(factor(x$detections$detections$receiver_serial,
                        levels = x$layout$receiver_serial))
  d <- rbind(
    data.frame(what = rep("animal-hour", length(e)), x = b$truth_x, y = b$truth_y,
               error_class = factor(labels[cls], levels = labels),
               n_detections = rep(NA_integer_, length(e))),
    data.frame(what = rep("receiver", nrow(x$layout)), x = x$layout$x, y = x$layout$y,
               error_class = factor(rep(NA_character_, nrow(x$layout)), levels = labels),
               n_detections = as.integer(heard))
  )
  grid <- water_grid(x$field$grid)$raster
  hours <- function(d) d[d$what == "animal-hour", ]
  ggplot2::ggplot(d, ggplot2::aes(x = .data$x, y = .data$y)) +
    depth_layers(grid, legend = legend) +
    # A black disc under each point gives it an edge on light and dark water,
    # since the fill scale belongs to the depth.
    ggplot2::geom_point(data = hours, shape = 16, colour = "black", size = 2.3) +
    ggplot2::geom_point(data = hours, ggplot2::aes(colour = .data$error_class), shape = 16,
                        size = 1.6, show.legend = TRUE) +
    receiver_rings(function(d) d[d$what == "receiver" & d$n_detections == 0, ],
                   colour = ring_white) +
    receiver_rings(function(d) d[d$what == "receiver" & d$n_detections > 0, ]) +
    # Reversed, so that a darker point is a larger error. Plasma ends short of
    # yellow, the colour of the receivers, and its steps stay apart at the
    # dark end, where viridis's last two classes look alike.
    ggplot2::scale_colour_viridis_d(name = sprintf("Error of\n%s (m)", error), direction = -1,
                                    option = "plasma", end = 0.8,
                                    drop = FALSE,
                                    guide = if (legend) ggplot2::guide_legend(
                                      override.aes = list(size = 2.6), order = 1
                                    ) else "none") +
    map_frame()
}
