#' Available steps for a step selection check on simulated tracks
#'
#' For each moving step of tracks simulated with a [habitat_selection()]
#' model, returns the step the animal took and `n` control steps drawn from
#' the same movement kernel: from the step's start and `heading`, one step
#' length and one turn from the model's samplers, any goal pull, and
#' reflection at the shore. A conditional logit of `case` on the layer, with
#' one stratum per step, then estimates the step-scale selection coefficient,
#' which for a log-linear preference is the `beta` the tracks were simulated
#' with:
#' `clogit(case ~ depth_m + strata(stratum), data = steps)`, with the
#' **survival** package attached.
#'
#' The estimator is consistent here because the controls come from the
#' kernel that proposed each step, so the ratio of used to available is the
#' preference alone. Steps of a resting animal, of an animal at its goal and
#' steps held after 1,000 rejected candidates have no `heading` and are left
#' out.
#'
#' Matching each used step with control steps drawn from a distribution of
#' steps, then fitting a conditional logistic regression, is step selection
#' analysis as Avgar et al. (2016) describe it. Their integrated form also
#' estimates the movement kernel from the data. Here the kernel is the
#' simulation's own and is known exactly.
#'
#' @param tracks A `sim_tracks` object from [simulate_tracks()] with a model
#'   from [habitat_selection()], as returned, with its attributes. Use
#'   `beta = 0` for a walk with no preference.
#' @param layer `NULL` to read the selection's own layer, or a
#'   [terra::SpatRaster] with the geometry of the tracks' grid, one column per
#'   layer, named by the layer names.
#' @param n Number of control steps per step.
#' @param seed Seed for the control steps; `NULL` uses the current random
#'   number stream. The caller's stream is restored after a seeded call.
#'
#' @return A data frame, `n + 1` rows per moving step: `path_id`, `step` (the
#'   step that ends at the used position, as in `tracks`), `stratum` (one per
#'   path and step, for `strata()`), `case` (1 for the step taken, 0 for a
#'   control), `x`, `y` and `cell` of the end of the step, one column per
#'   layer with its value at the cell, and `step_length`, the straight-line
#'   distance in metres from the start of the step to its end. Attributes `n`
#'   and `epsg`.
#' @references
#' Avgar, T., Potts, J. R., Lewis, M. A. and Boyce, M. S. (2016). Integrated
#' step selection analysis: bridging the gap between resource selection and
#' animal movement. *Methods in Ecology and Evolution*, 7, 619-630.
#' \doi{10.1111/2041-210X.12528}
#' (Step selection analysis, p. 620: each used step matched with available
#' steps drawn from a distribution of steps.)
#' @seealso [habitat_selection()], [simulate_tracks()]
#' @export
#' @examples
#' lake <- terra::rast(system.file("extdata", "stoney_depth.tif",
#'                                 package = "arraydesign"))
#' depth <- abs(lake)
#' names(depth) <- "depth_m"
#' shallow <- habitat_selection(crw_model(30, 15, 40), depth, beta = -0.3)
#'
#' # Three animals for 6 hours, fewer than a check needs, to keep the
#' # example short
#' tracks <- simulate_tracks(lake, shallow, n_animals = 3, n_steps = 360,
#'                           seed = 1)
#' steps <- available_steps(tracks, n = 10, seed = 2)
#' head(steps)
#'
#' # A conditional logit on the three animals, to set against beta = -0.3
#' # (clogit() reads a bare strata() term, so survival is attached)
#' if (requireNamespace("survival", quietly = TRUE)) {
#'   library(survival)
#'   fit <- clogit(case ~ depth_m + strata(stratum), data = steps)
#'   coef(fit)
#' }
#' # Three animals are too few for the estimate to sit near -0.3; section 5 of
#' # vignette("movement-models") fits 20 animals over 48 hours, and the 95%
#' # interval of their mean holds -0.3.
available_steps <- function(tracks, layer = NULL, n = 10, seed = NULL) {
  stop_unless(inherits(tracks, "sim_tracks") && !is.null(attr(tracks, "grid")) &&
                !is.null(attr(tracks, "model")),
              "`tracks` must be a sim_tracks object from simulate_tracks(), with its ",
              "attributes.")
  stop_unless("heading" %in% names(tracks),
              "`tracks` hold no `heading`: simulate them with a model from ",
              "habitat_selection(); `beta = 0` gives a walk with no preference.")
  stop_unless(is_scalar(n, 1, whole = TRUE), "`n` must be one whole number of at least 1.")
  model <- attr(tracks, "model")
  control <- attr(tracks, "control")
  wg <- water_grid(attr(tracks, "grid"))
  layer <- step_layers(layer, model, wg)
  d <- unclass_df(tracks)
  used <- which(!is.na(d$heading))
  stop_unless(length(used) > 0, "`tracks` have no moving steps.")
  stop_unless(all(d$step[used] >= 1L) && all(d$path_id[used - 1L] == d$path_id[used]),
              "`tracks` are not in path and step order, as simulate_tracks() returns them.")
  goal <- goal_point(wg, control)
  ends <- with_seed_or_not(seed, control_ends(wg, model, d, used, goal,
                                              control$max_reflections, n))
  steps_table(d, used, ends, layer, wg, n)
}

# The layer to read at the step ends: the selection's own, or the user's.
step_layers <- function(layer, model, wg) {
  if (is.null(layer)) {
    stop_unless(!is.null(model$selection), "The tracks' model has no selection; give `layer`.")
    return(terra::unwrap(model$selection$layer))
  }
  if (inherits(layer, "PackedSpatRaster")) layer <- terra::unwrap(layer)
  stop_unless(inherits(layer, "SpatRaster"), "`layer` must be a terra SpatRaster or NULL.")
  same <- isTRUE(terra::compareGeom(wg$raster, layer, crs = TRUE, ext = TRUE, rowcol = TRUE,
                                    stopOnError = FALSE))
  stop_unless(same, "`layer` does not match the geometry of the tracks' grid: give a ",
              "layer with the grid's extent, cells and CRS.")
  reserved <- c("path_id", "step", "stratum", "case", "x", "y", "cell", "step_length")
  clash <- intersect(names(layer), reserved)
  stop_unless(!anyDuplicated(names(layer)) && !length(clash),
              "`layer` names must be unique and not one of the output columns",
              if (length(clash)) paste0(" (", paste(clash, collapse = ", "), ")"), ".")
  layer
}

# Ends of `n` control steps for each used step, drawn from the kernel at the
# step's start and heading. Rows run step by step, n controls per step.
control_ends <- function(wg, model, d, used, goal, max_reflections, n) {
  from <- rep(used - 1L, each = n)
  heading <- rep(d$heading[used], each = n)
  m <- length(from)
  len <- model$step(m)
  turn <- model$turn(m)
  stop_unless(valid_draws(len, m, 0) && valid_draws(turn, m, -Inf),
              "The model's samplers returned draws that are not valid: ",
              "n finite numbers, and steps of at least 0.")
  x <- y <- numeric(m)
  for (i in seq_len(m)) {
    x0 <- d$x[from[i]]
    y0 <- d$y[from[i]]
    b <- pull_to_goal(heading[i] + turn[i], x0, y0, goal)
    mv <- take_step(wg, x0, y0, b, len[i], max_reflections)
    x[i] <- mv$x
    y[i] <- mv$y
  }
  cbind(x = x, y = y)
}

# One block of n + 1 rows per used step: the step taken, then its controls.
steps_table <- function(d, used, ends, layer, wg, n) {
  k <- length(used)
  block <- rep(seq_len(k), each = n + 1L)
  is_used <- rep(c(TRUE, rep(FALSE, n)), k)
  x <- y <- numeric(k * (n + 1L))
  x[is_used] <- d$x[used]
  y[is_used] <- d$y[used]
  x[!is_used] <- ends[, "x"]
  y[!is_used] <- ends[, "y"]
  start <- used[block] - 1L
  cell <- as.integer(cell_at(wg, x, y))
  out <- data.frame(path_id = d$path_id[used][block], step = d$step[used][block],
                    stratum = block, case = as.integer(is_used), x = x, y = y, cell = cell)
  values <- terra::values(layer, mat = TRUE)[cell, , drop = FALSE]
  for (j in seq_len(ncol(values))) out[[names(layer)[j]]] <- values[, j]
  out$step_length <- sqrt((x - d$x[start])^2 + (y - d$y[start])^2)
  attr(out, "n") <- as.integer(n)
  attr(out, "epsg") <- grid_epsg(wg$raster)
  out
}
