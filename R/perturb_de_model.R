#' A detection model with its range shifted, for misspecification tests
#'
#' Returns a copy of a detection model whose curve reaches 0.5 at
#' `d50_scale` times the original 50 percent distance and 0.05 at `d05_scale`
#' times the original 5 percent distance, and whose overdispersion standard
#' deviations are multiplied by `sigma_scale`. Simulating detections with one
#' model and positioning them with a perturbed copy shows what a wrong range
#' model does to the positions.
#'
#' The new curve reads the original one at a warped distance. With equal
#' scales the warp is a plain rescaling, `distance / scale`, which is exact
#' for any curve: every distance at which the curve takes a given value is
#' multiplied by the scale. With unequal scales the warp is the straight line
#' through the two anchor points, so a curve that is linear in distance on the
#' logit scale, such as [detfield::de_logistic()] or a logistic `glm` on
#' distance, stays logistic: `perturb_de_model(de_logistic(250, 600), 0.5)`
#' predicts what `de_logistic(125, 600)` predicts. The anchors are found by
#' bisection for each combination of depth, transmitter model and modulators
#' the curve is asked about. Near the receiver the line can give a negative
#' distance; a logistic curve or a `glm` extends its slope there, and a table
#' of curve values holds its first value.
#'
#' @param model A `de_model` from [detfield::de_model()] or
#'   [detfield::de_logistic()].
#' @param d50_scale Factor for the distance at which half of the
#'   transmissions are detected; above 0.
#' @param d05_scale Factor for the distance at which 5 percent are detected;
#'   above 0. With unequal scales the scaled 5 percent distance must stay
#'   beyond the scaled 50 percent distance. This is checked when the model is
#'   made, at the curve's default conditions; a curve that needs a depth or a
#'   transmitter model to be read is checked at each condition when it is
#'   asked.
#' @param sigma_scale Factor for `sigma$station_hour` and `sigma$station_day`;
#'   0 or more. A standard deviation that is `NA` stays `NA`.
#'
#' @return A `de_model` with the same modulators, `blocked` offset and
#'   transmitter models as `model`, a new `predict` function, the scaled
#'   `sigma`, and a `source` line that names the scales. The scales are also
#'   stored in the attribute `perturbation`.
#' @seealso [evaluate_array()], which takes a perturbed model as
#'   `control$position_model`.
#' @export
#' @examples
#' m <- detfield::de_logistic(d50 = 250, d05 = 600)
#' short <- perturb_de_model(m, d50_scale = 0.5)
#' short
#' plogis(short$predict(c(125, 600)))
#'
#' # Twice the range everywhere, and twice the per-hour spread.
#' m_sd <- detfield::de_model(m, sigma = list(station_hour = 1))
#' long <- perturb_de_model(m_sd, d50_scale = 2, d05_scale = 2, sigma_scale = 2)
#' long$sigma$station_hour
#' plogis(long$predict(c(500, 1200)))
#'
#' # Score the 20 planned stations when animals are detected under the true
#' # curve and positioned with one whose 50 and 5 percent distances are half
#' # of it. The depth raster is averaged to 200 m cells to keep the example
#' # short.
#' stoney <- terra::rast(system.file("extdata", "stoney_depth.tif",
#'                                   package = "arraydesign"))
#' stoney_200 <- terra::aggregate(stoney, 2, na.rm = TRUE)
#' planned <- read.csv(system.file("extdata", "planned_stations.csv",
#'                                 package = "arraydesign"))
#' rx <- place_receivers(stoney_200, method = "manual", coords = planned)
#' range_model <- detfield::de_logistic(d50 = c(250, 750), d05 = c(500, 1500),
#'                                      depth = c(1, 35))
#' short <- perturb_de_model(range_model, d50_scale = 0.5, d05_scale = 0.5)
#' ev_short <- evaluate_array(rx, stoney_200, range_model, n_animals = 20, seed = 1,
#'                            control = list(hours = 6, tag_delay = c(120, 180),
#'                                           position_model = short))
#' ev_short$score$overall[c("median_error_mean_m", "hit_rate")]
perturb_de_model <- function(model, d50_scale = 1, d05_scale = 1, sigma_scale = 1) {
  stop_unless(inherits(model, "de_model") && is.function(model$predict),
              "`model` must be a de_model from detfield::de_model() or de_logistic().")
  stop_unless(is_scalar(d50_scale, 1e-12) && is_scalar(d05_scale, 1e-12),
              "`d50_scale` and `d05_scale` must be single numbers above 0.")
  stop_unless(is_scalar(sigma_scale, 0), "`sigma_scale` must be one number, 0 or more.")
  sigma <- lapply(model$sigma, function(s) s * sigma_scale)
  out <- detfield::de_model(model, sigma = sigma, modulators = model$modulators,
                            blocked = if (is.null(model$blocked)) 0 else model$blocked)
  out$predict <- warped_curve(model$predict, d50_scale, d05_scale)
  if (abs(d50_scale - d05_scale) > 1e-12 * max(d50_scale, d05_scale)) {
    check_default_anchors(model, d50_scale, d05_scale)
  }
  out$source <- sprintf("perturb_de_model(d50 x %g, d05 x %g, sigma x %g) of: %s", d50_scale,
                        d05_scale, sigma_scale, model$source)
  attr(out, "perturbation") <- c(d50_scale = d50_scale, d05_scale = d05_scale,
                                 sigma_scale = sigma_scale)
  out
}

# The original curve read at a warped distance. The predict function keeps the
# signature detfield calls: distance_m, depth_m, tag_model and modulators.
warped_curve <- function(f, s50, s05) {
  force(f)
  if (abs(s50 - s05) <= 1e-12 * max(s50, s05)) {
    return(function(distance_m, depth_m = NA_real_, tag_model = NA_character_, ...) {
      f(distance_m / s50, depth_m = depth_m, tag_model = tag_model, ...)
    })
  }
  function(distance_m, depth_m = NA_real_, tag_model = NA_character_, ...) {
    stop_unless(is.numeric(distance_m), "`distance_m` must be numeric.")
    n <- length(distance_m)
    cond <- curve_conditions(n, depth_m, tag_model, list(...))
    key <- do.call(paste, c(unname(as.list(cond)), sep = "\r"))
    first <- !duplicated(key)
    anchors <- curve_anchors(f, cond[first, , drop = FALSE], s50, s05)
    at <- match(key, key[first])
    a <- anchors$d50[at]
    b <- anchors$d05[at]
    w <- a + (distance_m - s50 * a) * (b - a) / (s05 * b - s50 * a)
    args <- c(list(w), as.list(cond))
    names(args)[1] <- "distance_m"
    do.call(f, args)
  }
}

# The anchors at the curve's default conditions: no depth, no transmitter
# model and the modulators' default values. A swap, a curve at or below 0.5
# at the receiver, or one that never falls to 0.05 stops here, when the model
# is made. A curve that cannot be read without a depth or a transmitter model
# is left to the check in warped_curve().
check_default_anchors <- function(model, s50, s05) {
  cond <- curve_conditions(1L, NA_real_, NA_character_, model$modulators)
  at0 <- tryCatch(do.call(model$predict, c(list(distance_m = 0), as.list(cond))),
                  error = function(e) NA_real_)
  if (length(at0) == 1L && is.finite(at0)) curve_anchors(model$predict, cond, s50, s05)
  invisible()
}

# One row per distance asked about, holding every argument of the curve other
# than distance, so rows that share a row of this table share their anchors.
curve_conditions <- function(n, depth_m, tag_model, mods) {
  cond <- data.frame(depth_m = rep_len(as.numeric(depth_m), n),
                     tag_model = rep_len(as.character(tag_model), n), stringsAsFactors = FALSE)
  for (nm in names(mods)) cond[[nm]] <- rep_len(mods[[nm]], n)
  cond
}

# The 50 and 5 percent distances of the curve under each condition, by bisection on
# all conditions at once. Stops when a curve is not above 0.5 at the receiver,
# never falls to 0.05, or when the scaled anchors would swap order.
curve_anchors <- function(f, cond, s50, s05) {
  at <- function(d) do.call(f, c(list(distance_m = d), as.list(cond)))
  stop_unless(all(at(rep(0, nrow(cond))) > 0),
              "The curve is at or below 0.5 at the receiver, so it has no 50 percent distance to ",
              "scale. Use equal `d50_scale` and `d05_scale`.")
  hi <- rep(1, nrow(cond))
  for (i in seq_len(60)) {
    low <- at(hi) > stats::qlogis(0.05)
    if (!any(low)) break
    hi[low] <- hi[low] * 2
  }
  stop_unless(!any(low), "The curve does not fall to 0.05 within 1e18 m.")
  d50 <- bisect_distance(at, 0, rep(0, nrow(cond)), hi)
  d05 <- bisect_distance(at, stats::qlogis(0.05), d50, hi)
  stop_unless(all(s05 * d05 > s50 * d50),
              "The scaled 5 percent distance must stay beyond the scaled 50 percent ",
              "distance; it does not for d50_scale ", s50, " and d05_scale ", s05, ".")
  list(d50 = d50, d05 = d05)
}

# Distances at which a decreasing curve crosses `target`, between lo and hi,
# to a relative precision of 2^-50.
bisect_distance <- function(at, target, lo, hi) {
  for (i in seq_len(50)) {
    mid <- (lo + hi) / 2
    above <- at(mid) > target
    lo[above] <- mid[above]
    hi[!above] <- mid[!above]
  }
  (lo + hi) / 2
}
