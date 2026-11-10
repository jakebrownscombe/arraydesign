#' Score location surfaces against known positions
#'
#' For simulated animals, whose true positions are known, measures how far
#' each location surface lies from the truth and whether its credible region
#' holds the truth. The truth of an animal-bin is every transmission the
#' animal sent in the bin: the transmissions are joined to the bins by
#' `animal_id` and time, `bin_start <= datetime < bin_end`.
#'
#' Per animal-bin:
#' * `error_mode_m` and `error_mean_m`: the distance from the centre of the
#'   surface's most probable cell, and from the probability-weighted mean of
#'   its cell centres, to the centroid of the bin's transmission positions.
#' * `in_region`: the share of the bin's transmissions that fall in a cell of
#'   the credible region at `level`. For an animal that stays in one cell for
#'   the bin, as the likelihood assumes, this is 0 or 1; an animal that moves
#'   across the region's edge in the bin scores in between.
#' * `area_m2` and `region_mass`: the area of the credible region and the
#'   probability it holds. A region is whole cells, so it holds at least
#'   `level`; compare the hit rate with the mean mass as well as with the level.
#'
#' **Hit rate.** The hit rate is the mean of `in_region` over bins, so each
#' bin counts once: the share of the truth that the credible regions hold. It
#' checks whether the regions are as wide as their level says. Compare designs
#' on the error and the scored share. When animals move within the hour it falls as
#' the regions shrink, so a denser array scores lower; see [evaluate_array()].
#'
#' Bins share an animal, and animals share an hour: a receiver that hears
#' badly in one hour does so for every animal. `hit_rate_se` is therefore
#' clustered two ways, by animal and by bin start (Cameron, Gelbach and
#' Miller 2011, \doi{10.1198/jbes.2010.07136}): `sqrt(V_a + V_t - V_b)`,
#' where `V_a` and `V_t` are the one-way cluster variances, each
#' `G / (G - 1) * sum_g r_g^2 / N^2` with `r_g` the summed residuals of
#' cluster `g`, `G` the clusters and `N` the bins, and `V_b` the same over
#' single bins. A way with one cluster drops out. The result is never taken
#' below the larger one-way SE, a conservative floor chosen here: the
#' subtraction can pull the two-way value under both when one way carries no
#' correlation, and the paper's own fallback for a negative value is the
#' larger one-way SE. `hit_rate_se_from` says which estimate was used:
#' `"animal and bin start"`, `"animal (floor)"` or `"bin start (floor)"`;
#' `"animal"` or `"bin start"` when the other way has one cluster;
#' `"not estimable"` when every bin scores the same; `NA` for one bin.
#' `hit_rate_se_animal` clusters by animal only and `hit_rate_se_bins` treats
#' bins as independent; both are too small when bins share an hour's
#' detection conditions.
#'
#' The 95% interval `hit_rate_lo` to `hit_rate_hi` uses the t distribution
#' with `hit_rate_df` degrees of freedom, the smaller cluster count less one,
#' the small-sample adjustment Cameron, Gelbach and Miller (2011) use because
#' the SE rests on few clusters. With few bin starts, such as the six hours of
#' [evaluate_array()]'s defaults, the SE itself is noisy and the interval is
#' rough. The script that measures how often it holds,
#' <https://github.com/jakebrownscombe/arraydesign/blob/main/validation/02_scoring_validation.R>,
#' is in the repository and not in the installed package. When every bin
#' scores the same, the SE cannot be estimated: `hit_rate_se` is `NA` and,
#' at a hit rate of 0 or 1, the
#' interval is the exact binomial bound for that many successes or failures
#' in as many trials as the smaller cluster count.
#'
#' @param surfaces Location surfaces from [detfield::location_surface()].
#' @param truth The true transmissions: a data frame with `animal_id`,
#'   `datetime` (POSIXct, UTC, no `NA`), `x` and `y` in the CRS of the
#'   surfaces' grid, such as [simulate_transmissions()] returns or the `transmissions` element
#'   of [simulate_detections()].
#' @param level Share of the probability the credible region holds; above 0
#'   and at most the kept mass of the surfaces.
#'
#' @return A list of class `surface_score` with three data frames:
#'   * `bins`: one row per animal-bin, with `animal_id`, `bin_start`,
#'     `n_truth` (transmissions in the bin), `truth_x`, `truth_y`, `mode_x`,
#'     `mode_y`, `mean_x`, `mean_y`, `error_mode_m`, `error_mean_m`,
#'     `in_region`, `area_m2` and `region_mass`.
#'   * `animals`: one row per animal, with `n_bins` and the summaries below.
#'   * `overall`: one row, with `n_animals`, `n_bin_starts`, `n_bins`, the
#'     summaries below, `mean_region_mass`, `hit_rate_se`,
#'     `hit_rate_se_from`, `hit_rate_df`, `hit_rate_lo`, `hit_rate_hi`,
#'     `hit_rate_se_animal` and `hit_rate_se_bins`.
#'
#'   The summaries are `median_error_mode_m`, `p90_error_mode_m`,
#'   `median_error_mean_m`, `p90_error_mean_m`, `hit_rate` and
#'   `median_area_m2`. Attributes: `level`, and `no_truth`, the number of bins
#'   left out because no transmission of the animal fell in them.
#' @seealso [evaluate_array()], which simulates, positions and scores in one
#'   call.
#' @export
#' @examples
#' stoney <- terra::rast(system.file("extdata", "stoney_depth.tif",
#'                                   package = "arraydesign"))
#' planned <- read.csv(system.file("extdata", "planned_stations.csv",
#'                                 package = "arraydesign"))
#' rx <- place_receivers(stoney, method = "manual", coords = planned)
#' range_model <- detfield::de_logistic(d50 = c(250, 750), d05 = c(500, 1500),
#'                                      depth = c(1, 35))
#' # Pairs more than 4 km apart through water are left out. At the deepest
#' # receiver the curve gives them a detection probability below 1e-7, and
#' # the distances take a third of the time.
#' field <- detfield::de_field(detfield::water_distances(stoney, rx, max_distance = 4000),
#'                             range_model)
#' # Five animals for 3 hours each.
#' tr <- simulate_tracks(stoney, crw_model(30, 15, 40, time_step_s = 60),
#'                       n_animals = 5, n_steps = 180, seed = 1)
#' tx <- simulate_transmissions(tr, tag_delay = c(120, 180), seed = 2)
#' det <- simulate_detections(tx, field, overdispersion = "none", seed = 3)
#' listening <- data.frame(receiver_serial = rx$receiver_serial,
#'                         date = as.Date("2025-01-01"))
#' counts <- detfield::bin_detections(det$detections, listening,
#'                                    tag_delay = c(120, 180))
#' surfaces <- detfield::location_surface(counts, field, overdispersion = "none")
#' score <- score_surfaces(surfaces, det$transmissions)
#' score
#' head(score$bins)
score_surfaces <- function(surfaces, truth, level = 0.95) {
  stop_unless(inherits(surfaces, "loc_surface"),
              "`surfaces` must come from detfield::location_surface().")
  truth <- check_truth(truth)
  wg <- water_grid(attr(surfaces, "grid"))
  bins <- attr(surfaces, "bins")[c("animal_id", "bin_start", "bin_end")]
  region <- detfield::credible_region(surfaces, level = level)
  joined <- join_truth(truth, bins)
  check_truth_on_grid(joined$truth, wg)
  scored <- score_bins(surfaces, region, joined, wg)
  missing <- is.na(scored$n_truth)
  if (any(missing)) {
    warning(sum(missing), " of ", nrow(scored), " bins hold no transmission of their animal ",
            "in `truth` and are left out.", call. = FALSE)
  }
  scored <- scored[!missing, , drop = FALSE]
  rownames(scored) <- NULL
  structure(list(bins = scored, animals = animal_summaries(scored),
                 overall = overall_summary(scored)),
            class = "surface_score", level = level, no_truth = sum(missing))
}

check_truth <- function(truth) {
  need <- c("animal_id", "datetime", "x", "y")
  stop_unless(is.data.frame(truth) && all(need %in% names(truth)),
              "`truth` must be a data frame with columns ", paste(need, collapse = ", "), ".")
  check_utc(truth$datetime, "truth$datetime")
  stop_unless(!anyNA(truth$datetime), "`truth$datetime` must have no NA; ",
              sum(is.na(truth$datetime)), " of ", nrow(truth), " rows have one.")
  stop_unless(is.numeric(truth$x) && is.numeric(truth$y) && !anyNA(truth$x) && !anyNA(truth$y),
              "`truth$x` and `truth$y` must be numbers with no NA.")
  data.frame(animal_id = as.character(truth$animal_id), datetime = truth$datetime,
             x = truth$x, y = truth$y, stringsAsFactors = FALSE)
}

# Truth in the wrong CRS scores silently as large errors, so truth off the
# grid stops when none of it is on, and warns otherwise.
check_truth_on_grid <- function(t, wg) {
  if (!nrow(t)) return(invisible())
  cell <- cell_at(wg, t$x, t$y)
  off <- is.na(cell)
  stop_unless(!all(off), "None of the ", nrow(t), " transmissions in the surfaces' bins lies ",
              "on the surfaces' grid. Is `truth` in the CRS of the grid?")
  land <- !off & !wg$wet[ifelse(off, 1L, cell)]
  if (any(off | land)) {
    warning(sum(off), " of ", nrow(t), " transmissions in the surfaces' bins lie outside the ",
            "grid and ", sum(land), " on land. Is `truth` in the CRS of the grid?", call. = FALSE)
  }
  invisible()
}

# Each transmission gets the row of the bin that holds it, or NA.
join_truth <- function(truth, bins) {
  bins$row <- seq_len(nrow(bins))
  by_animal <- split(bins, bins$animal_id)
  truth$row <- NA_integer_
  for (a in intersect(names(by_animal), unique(truth$animal_id))) {
    b <- by_animal[[a]][order(by_animal[[a]]$bin_start), ]
    i <- which(truth$animal_id == a)
    k <- findInterval(as.numeric(truth$datetime[i]), as.numeric(b$bin_start))
    inside <- k > 0
    inside[inside] <- truth$datetime[i][inside] < b$bin_end[k[inside]]
    truth$row[i[inside]] <- b$row[k[inside]]
  }
  list(truth = truth[!is.na(truth$row), ], bins = bins)
}

# One row per bin of the surfaces: truth centroid, mode, mean, region.
score_bins <- function(surfaces, region, joined, wg) {
  bins <- joined$bins
  s <- data.frame(row = match(paste(surfaces$animal_id, as.numeric(surfaces$bin_start)),
                              paste(bins$animal_id, as.numeric(bins$bin_start))),
                  cell = surfaces$cell, prob = surfaces$prob)
  xy <- cell_xy(wg, s$cell)
  w <- rowsum(s$prob, s$row)
  mean_xy <- rowsum(xy * s$prob, s$row) / as.vector(w)
  s <- s[order(s$row, -s$prob), ]
  top <- s[!duplicated(s$row), ]
  mode_xy <- cell_xy(wg, top$cell)
  t <- joined$truth
  n_truth <- tabulate(t$row, nbins = nrow(bins))
  cen <- by_bin(cbind(t$x, t$y), t$row, nrow(bins)) / n_truth
  summ <- attr(region, "summary")
  r_row <- match(paste(region$animal_id, as.numeric(region$bin_start)),
                 paste(bins$animal_id, as.numeric(bins$bin_start)))
  in_reg <- paste(t$row, cell_at(wg, t$x, t$y)) %in% paste(r_row, region$cell)
  out <- data.frame(animal_id = bins$animal_id, bin_start = bins$bin_start,
                    n_truth = ifelse(n_truth > 0, n_truth, NA_integer_),
                    truth_x = cen[, 1], truth_y = cen[, 2], stringsAsFactors = FALSE)
  out$mode_x <- mode_xy[match(seq_len(nrow(bins)), top$row), 1]
  out$mode_y <- mode_xy[match(seq_len(nrow(bins)), top$row), 2]
  out$mean_x <- mean_xy[match(seq_len(nrow(bins)), rownames(w)), 1]
  out$mean_y <- mean_xy[match(seq_len(nrow(bins)), rownames(w)), 2]
  out$error_mode_m <- sqrt((out$mode_x - out$truth_x)^2 + (out$mode_y - out$truth_y)^2)
  out$error_mean_m <- sqrt((out$mean_x - out$truth_x)^2 + (out$mean_y - out$truth_y)^2)
  out$in_region <- as.vector(by_bin(as.numeric(in_reg), t$row, nrow(bins))) / n_truth
  at <- match(paste(bins$animal_id, as.numeric(bins$bin_start)),
              paste(summ$animal_id, as.numeric(summ$bin_start)))
  out$area_m2 <- summ$area_m2[at]
  out$region_mass <- summ$mass[at]
  out
}

# Column sums of x by bin row, one row for every bin, 0 where a bin has none.
by_bin <- function(x, row, n_bins) {
  x <- as.matrix(x)
  out <- matrix(0, n_bins, ncol(x))
  if (length(row)) {
    rs <- rowsum(x, row)
    out[as.integer(rownames(rs)), ] <- rs
  }
  out
}

# The six summaries of a set of scored bins.
score_summaries <- function(b) {
  q <- function(v, p) unname(stats::quantile(v, p, names = FALSE))
  data.frame(median_error_mode_m = q(b$error_mode_m, 0.5),
             p90_error_mode_m = q(b$error_mode_m, 0.9),
             median_error_mean_m = q(b$error_mean_m, 0.5),
             p90_error_mean_m = q(b$error_mean_m, 0.9),
             hit_rate = mean(b$in_region), median_area_m2 = q(b$area_m2, 0.5))
}

animal_summaries <- function(scored) {
  if (!nrow(scored)) return(data.frame())
  by <- split(scored, factor(scored$animal_id, unique(scored$animal_id)))
  out <- do.call(rbind, lapply(by, score_summaries))
  out <- cbind(data.frame(animal_id = names(by), n_bins = vapply(by, nrow, integer(1)),
                          stringsAsFactors = FALSE), out)
  rownames(out) <- NULL
  out
}

overall_summary <- function(scored) {
  n <- nrow(scored)
  if (!n) return(data.frame(n_animals = 0L, n_bin_starts = 0L, n_bins = 0L))
  out <- cbind(data.frame(n_animals = length(unique(scored$animal_id)),
                          n_bin_starts = length(unique(as.numeric(scored$bin_start))),
                          n_bins = n),
               score_summaries(scored))
  out$mean_region_mass <- mean(scored$region_mass)
  cbind(out, hit_rate_interval(scored, out$hit_rate))
}

# SE, its source, degrees of freedom and 95% interval of the hit rate.
hit_rate_interval <- function(scored, hr) {
  e <- scored$in_region - hr
  n <- length(e)
  time <- as.numeric(scored$bin_start)
  cs <- cluster_se(e, scored$animal_id, time)
  g <- c(length(unique(scored$animal_id)), length(unique(time)))
  m <- if (any(g >= 2)) min(g[g >= 2]) else 1L
  out <- data.frame(hit_rate_se = cs$se, hit_rate_se_from = cs$from, hit_rate_df = m - 1L,
                    hit_rate_lo = NA_real_, hit_rate_hi = NA_real_, stringsAsFactors = FALSE)
  if (all(e == 0)) {
    out$hit_rate_se <- NA_real_
    out$hit_rate_se_from <- "not estimable"
    # Exact bound with the smaller cluster count as the number of trials.
    if (hr == 1) out[c("hit_rate_lo", "hit_rate_hi")] <- c(0.025^(1 / m), 1)
    if (hr == 0) out[c("hit_rate_lo", "hit_rate_hi")] <- c(0, 1 - 0.025^(1 / m))
  } else if (!is.na(cs$se) && m >= 2) {
    half <- stats::qt(0.975, m - 1) * cs$se
    out$hit_rate_lo <- max(0, hr - half)
    out$hit_rate_hi <- min(1, hr + half)
  }
  out$hit_rate_se_animal <- cluster_se(e, scored$animal_id, rep(0, n))$se
  out$hit_rate_se_bins <- if (n > 1) stats::sd(scored$in_region) / sqrt(n) else NA_real_
  out
}

# Standard error of a mean from its residuals `e`, clustered two ways: by
# animal and by bin start (Cameron, Gelbach and Miller 2011). Each animal-bin
# is one row, so the intersection term is the rows themselves. A way with one
# cluster drops out; with neither way clustered the SE is NA. The result is
# floored at the larger one-way value. `from` names the estimate returned.
cluster_se <- function(e, animal, time) {
  n <- length(e)
  part <- function(g) {
    k <- length(unique(g))
    if (k < 2) return(NA_real_)
    k / (k - 1) * sum(rowsum(e, g)^2) / n^2
  }
  v_a <- part(animal)
  v_t <- part(time)
  if (is.na(v_a) && is.na(v_t)) return(list(se = NA_real_, from = NA_character_))
  if (is.na(v_t)) return(list(se = sqrt(v_a), from = "animal"))
  if (is.na(v_a)) return(list(se = sqrt(v_t), from = "bin start"))
  v <- v_a + v_t - n / (n - 1) * sum(e^2) / n^2
  from <- if (v >= max(v_a, v_t)) "animal and bin start" else if (v_a >= v_t)
    "animal (floor)" else "bin start (floor)"
  list(se = sqrt(max(v, v_a, v_t)), from = from)
}

#' @export
print.surface_score <- function(x, ...) {
  o <- x$overall
  cat_wrapped(sprintf("<surface_score> %d animal-bins of %d animals; %g%% credible regions",
                      o$n_bins, o$n_animals, 100 * attr(x, "level")))
  if (!is.null(o$n_animal_hours)) {
    cat_wrapped(sprintf("scored %d of %d animal-hours with transmissions (%.1f%%); %s",
                        o$n_bins, o$n_animal_hours, 100 * o$scored_share,
                        "hours with no detection are not scored"), indent = 2, exdent = 4)
  }
  if (o$n_bins > 0) print_score_lines(o, attr(x, "level"))
  if (attr(x, "no_truth") > 0) cat(sprintf("  %d bins had no truth\n", attr(x, "no_truth")))
  invisible(x)
}

print_score_lines <- function(o, level) {
  cat(sprintf("  error of the mode: median %.0f m, 90th percentile %.0f m\n",
              o$median_error_mode_m, o$p90_error_mode_m))
  cat(sprintf("  error of the mean: median %.0f m, 90th percentile %.0f m\n",
              o$median_error_mean_m, o$p90_error_mean_m))
  cat_wrapped(sprintf("hit rate %.3f: share of the truth inside the %g%% regions %s",
                      o$hit_rate, 100 * level, "(a calibration check)"),
              indent = 2, exdent = 4)
  cat_wrapped(sprintf("%s; regions hold %.3f on average", hit_rate_ci_text(o),
                      o$mean_region_mass), indent = 4, exdent = 6)
  cat(sprintf("  median region %s m^2\n", format(round(o$median_area_m2), big.mark = ",")))
}

hit_rate_ci_text <- function(o) {
  ci <- if (is.na(o$hit_rate_lo)) "95% interval not available" else
    sprintf("95%% interval %.3f to %.3f", o$hit_rate_lo, o$hit_rate_hi)
  clusters <- sprintf("%d animals and %d bin starts", o$n_animals, o$n_bin_starts)
  if (is.na(o$hit_rate_se)) {
    why <- if (identical(o$hit_rate_se_from, "not estimable")) "every bin scores the same" else
      "one cluster each way"
    return(sprintf("%s (SE not estimable, %s; %s)", ci, why, clusters))
  }
  how <- switch(o$hit_rate_se_from,
                "animal (floor)" = "clustered by animal, the one-way floor",
                "bin start (floor)" = "clustered by bin start, the one-way floor",
                "animal and bin start" = "clustered two ways",
                paste("clustered by", o$hit_rate_se_from))
  sprintf("%s (SE %.3f %s; %s; t with %d df)", ci, o$hit_rate_se, how, clusters, o$hit_rate_df)
}
