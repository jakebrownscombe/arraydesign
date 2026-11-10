#' Give a movement model a preference for a habitat layer
#'
#' Attaches a habitat preference to a [movement_model()], so that
#' [simulate_tracks()] walks animals that select for, or against, the values
#' of a raster layer with a known strength. The preference is a log weight per
#' cell, `beta * layer` for a log-linear preference, or `preference(layer)` for
#' any other shape, such as an optimum depth.
#'
#' The walker draws each moving step from the model without the preference: a
#' step length and a turn, any goal pull on the bearing, and reflection at the
#' shore. It accepts the step with probability `w(end) / W`, where `w` is the
#' weight of the cell the step ends in and `W` the largest weight the step
#' could reach from its start. A rejected step is drawn again from the same
#' start and heading. This is rejection sampling, so the steps have exactly
#' the density `K(x, y) w(y) / integral of K w`, where `K` is the walk without
#' the preference. That density is the normalised product of a movement kernel
#' and a habitat-selection function, Eqn 1 of Avgar et al. (2016). `beta` is
#' therefore the step-scale selection coefficient that a conditional logit on
#' the true steps estimates (see [available_steps()]).
#'
#' `W` is the largest weight among the cells within the reach `R` of the start
#' of the step, where `R` is the 0.9999 quantile of 100,000 step lengths from
#' the model. A step longer than `R`, about 1 in 10,000, can end beyond those
#' cells; it is accepted with probability `min(1, w / W)` and counted as
#' `over_reach`. A step whose weight exceeds `W` is counted as `above_bound`.
#' A step with 1,000 rejected candidates holds its position and is counted as
#' `exhausted`. Each moving step therefore costs about `1 / acceptance`
#' candidates. On Stoney Lake with bottom depth, 20 animals walked for 48
#' hours from the same starts kept 0.51 of their candidates at `beta = -0.3`
#' per m and 0.39 at `beta = 0.3`. The script that gives these rates,
#' <https://github.com/jakebrownscombe/arraydesign/blob/main/validation/04_selection_acceptance.R>,
#' and its log in `validation/outputs/` are in the repository and not in the
#' installed package.
#'
#' The weight is read at the cell that holds the end of the step, without
#' interpolation, which is the value an analysis of the positions reads.
#'
#' @param model A [movement_model()] without a preference.
#' @param layer A one-layer [terra::SpatRaster] of the habitat variable. At
#'   simulation it must have the geometry of the grid
#'   ([terra::compareGeom()]) and a finite value on every water cell. For
#'   bottom depth from a raster of negative depths, pass `abs(lake)` and name
#'   the layer, for example `depth_m`.
#' @param beta One number: the log selection strength per unit of `layer` at
#'   the scale of one step. Negative values select low values of the layer.
#'   It is stored as given, so a fit to the tracks can be compared with it.
#' @param preference A function of a vector of layer values that returns one
#'   finite log weight per value, for a preference that is not log-linear.
#'   Give exactly one of `beta` and `preference`.
#'
#' @return The `movement_model`, with an element `selection` (the layer, its
#'   name, `beta` and `preference`) and the source line extended. Tracks
#'   simulated with it carry a column `heading` and an attribute `selection`
#'   that counts, per animal, the candidates drawn (`proposals`), the steps
#'   accepted (`accepted`), `over_reach`, `above_bound` and `exhausted`.
#' @references
#' Avgar, T., Potts, J. R., Lewis, M. A. and Boyce, M. S. (2016). Integrated
#' step selection analysis: bridging the gap between resource selection and
#' animal movement. *Methods in Ecology and Evolution*, 7, 619-630.
#' \doi{10.1111/2041-210X.12528}
#' (Eqn 1, p. 622: the density of the next position as the normalised product
#' of a movement kernel and a habitat-selection function.)
#' @seealso [available_steps()], [simulate_tracks()]
#' @export
#' @examples
#' lake <- terra::rast(system.file("extdata", "stoney_depth.tif",
#'                                 package = "arraydesign"))
#' depth <- abs(lake)
#' names(depth) <- "depth_m"
#' walk <- crw_model(30, 15, 40, time_step_s = 60)
#'
#' # Select shallow water: log weight -0.3 per m of bottom depth
#' shallow <- habitat_selection(walk, depth, beta = -0.3)
#' shallow
#'
#' # Three animals for 6 hours. The print counts candidates and accepted
#' # steps.
#' tracks <- simulate_tracks(lake, shallow, n_animals = 3, n_steps = 360,
#'                           seed = 1)
#' tracks
#' attr(tracks, "selection")
#'
#' # An optimum at 4 m bottom depth instead of a straight-line preference
#' habitat_selection(walk, depth,
#'                   preference = function(z) -0.1 * (z - 4)^2)
habitat_selection <- function(model, layer, beta = NULL, preference = NULL) {
  stop_unless(inherits(model, "movement_model"),
              "`model` must be a movement_model; see movement_model().")
  stop_unless(is.null(model$selection),
              "`model` already has a selection; start from the model without one.")
  layer <- one_layer(layer, "layer")
  stop_unless(xor(is.null(beta), is.null(preference)),
              "Give exactly one of `beta` and `preference`.")
  values <- terra::values(layer, mat = FALSE)
  values <- values[!is.na(values)]
  stop_unless(length(values) > 0, "`layer` has no values: every cell is NA.")
  stop_unless(all(is.finite(values)), "`layer` holds ", sum(!is.finite(values)),
              " infinite values; give finite numbers, with NA on land.")
  if (is.null(preference)) {
    stop_unless(is_scalar(beta) && is.finite(beta), "`beta` must be one finite number.")
    what <- sprintf("beta %g per unit", beta)
  } else {
    check_preference(preference, values)
    what <- "a log-preference function"
  }
  name <- names(layer)
  model$selection <- list(layer = terra::wrap(layer), layer_name = name, beta = beta,
                          preference = preference)
  model$source <- paste0(model$source, if (nzchar(model$source)) "; ",
                         "selection on ", name, ": ", what)
  model
}

# A one-layer SpatRaster, or a stop naming the argument.
one_layer <- function(layer, arg) {
  if (inherits(layer, "PackedSpatRaster")) layer <- terra::unwrap(layer)
  stop_unless(inherits(layer, "SpatRaster"), "`", arg, "` must be a terra SpatRaster.")
  stop_unless(terra::nlyr(layer) == 1L, "`", arg, "` must have one layer; it has ",
              terra::nlyr(layer), ".")
  layer
}

# The function must return one finite log weight per value it is given.
check_preference <- function(preference, values) {
  stop_unless(is.function(preference), "`preference` must be a function of layer values.")
  out <- preference(values)
  ok <- is.numeric(out) && length(out) == length(values)
  stop_unless(ok, "`preference` must return one number per layer value; given ",
              length(values), " values it returned ", length(out), ".")
  stop_unless(all(is.finite(out)), "`preference` returned ", sum(!is.finite(out)),
              " values that are not finite numbers; it must return finite log weights.")
  invisible()
}

# ---- Selection at simulation ----
# Log weight of every cell (-Inf on land), the reach R of a step, and the
# window of cell offsets that can hold a point within R of a point in the
# centre cell. The bound W of a step is read from that window at its start.
selection_setup <- function(wg, model) {
  sel <- model$selection
  layer <- terra::unwrap(sel$layer)
  same <- isTRUE(terra::compareGeom(wg$raster, layer, crs = TRUE, ext = TRUE, rowcol = TRUE,
                                    stopOnError = FALSE))
  stop_unless(same, "The selection layer (", sel$layer_name, ") does not match the ",
              "geometry of `grid`: give a layer with the grid's extent, cells and CRS.")
  z <- terra::values(layer, mat = FALSE)[wg$wet]
  log_w <- rep(-Inf, length(wg$wet))
  log_w[wg$wet] <- if (is.null(sel$preference)) sel$beta * z else sel$preference(z)
  bad <- sum(!is.finite(log_w[wg$wet]))
  stop_unless(bad == 0, "The selection layer (", sel$layer_name, ") gives no finite log ",
              "weight on ", bad, " of ", sum(wg$wet), " water cells; it must have a ",
              "finite value on every water cell of `grid`.")
  reach <- step_reach(model)
  offsets <- window_offsets(reach, wg$res)
  list(log_w = log_w, offsets = offsets, reach = reach, window = nrow(offsets),
       max_tries = 1000L)
}

# The 0.9999 quantile of 100,000 step lengths. It is drawn under its own seed
# so that it is a property of the model, and the caller's stream is restored.
step_reach <- function(model) {
  draws <- withr::with_seed(1L, model$step(1e5))
  stats::quantile(draws, 0.9999, names = FALSE)
}

# Row and column offsets of the cells that a step of length `reach` or less
# can end in from some point of the centre cell: those whose square lies
# within `reach` of the centre cell's square.
window_offsets <- function(reach, res) {
  k <- ceiling(reach / res) + 1L
  off <- expand.grid(di = -k:k, dj = -k:k)
  gap <- res * sqrt(pmax(abs(off$di) - 1, 0)^2 + pmax(abs(off$dj) - 1, 0)^2)
  as.matrix(off[gap <= reach, ])
}

# Log bound W of a step from (x, y): the largest log weight among the cells
# whose square lies within the reach of the point. Every step of length
# `reach` or less ends in one of them, so w(end) / W is at most 1.
point_bound <- function(wg, sel, x, y) {
  cell <- cell_at(wg, x, y)
  row <- (cell - 1) %/% wg$ncol + sel$offsets[, 1]
  col <- (cell - 1) %% wg$ncol + sel$offsets[, 2]
  inside <- row >= 0 & row < wg$nrow & col >= 0 & col < wg$ncol
  row <- row[inside]
  col <- col[inside]
  left <- wg$xmin + col * wg$res
  top <- wg$ymax - row * wg$res
  dx <- pmax(left - x, 0, x - left - wg$res)
  dy <- pmax(top - wg$res - y, 0, y - top)
  near <- dx^2 + dy^2 <= sel$reach^2
  max(sel$log_w[row[near] * wg$ncol + col[near] + 1])
}

# One line for the print of a movement model.
selection_text <- function(selection) {
  what <- if (is.null(selection$preference)) {
    sprintf("log weight %g per unit of %s", selection$beta, selection$layer_name)
  } else {
    sprintf("log weight from a function of %s", selection$layer_name)
  }
  paste("selection:", what, "(rejection sampling)")
}
