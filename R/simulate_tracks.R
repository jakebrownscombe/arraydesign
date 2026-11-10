#' Settings for the walker
#'
#' Collects the settings of [simulate_tracks()] that most runs leave at their
#' defaults: where the animals start, an optional goal, when the clock starts,
#' and how many reflections a step may take at the shore.
#'
#' @param start `NULL` to start each animal at the centre of a random water
#'   cell, or a data frame (or matrix) with columns `x` and `y`: one row for
#'   all animals, or one row per animal. Every start must be on water. With a
#'   goal and `start = NULL`, starts are drawn from the water connected to
#'   the goal.
#' @param goal `NULL`, or one point with `x` and `y` on water that pulls every
#'   animal towards it. An animal that comes within `goal_tolerance` of the goal
#'   stays there for the rest of the simulation. The goal must lie in water
#'   that every start can reach.
#' @param goal_bias Weight of the bearing to the goal, 0 to 1, against the
#'   bearing the walk would take; 0 ignores the goal and 1 heads straight for
#'   it. The pull acts on the heading only, so a goal behind a peninsula can
#'   hold an animal against the shore.
#' @param goal_tolerance Distance in metres at which an animal has reached the
#'   goal; `NULL` uses the cell size of the grid.
#' @param start_time Time of the first position, POSIXct in UTC.
#' @param max_reflections Reflections off the shore allowed within one step
#'   before the animal holds its position for that step.
#'
#' @return A list of class `track_control`.
#' @seealso [simulate_tracks()]
#' @export
#' @examples
#' # Start in the south-west arm of Stoney Lake and head for a bay on the
#' # north shore, 8.3 km away; an animal within 300 m has arrived.
#' track_control(start = data.frame(x = 721950, y = 4930450),
#'               goal = data.frame(x = 725500, y = 4938000), goal_bias = 0.2,
#'               goal_tolerance = 300)
track_control <- function(start = NULL, goal = NULL, goal_bias = 0.3, goal_tolerance = NULL,
                          start_time = as.POSIXct("2025-01-01 00:00:00", tz = "UTC"),
                          max_reflections = 50) {
  stop_unless(is_scalar(goal_bias, 0, 1), "`goal_bias` must be one number from 0 to 1.")
  stop_unless(is.null(goal_tolerance) || is_scalar(goal_tolerance, 0),
              "`goal_tolerance` must be one number of metres, 0 or more.")
  check_utc(start_time, "start_time")
  stop_unless(length(start_time) == 1L && !is.na(start_time), "`start_time` must be one time.")
  stop_unless(is_scalar(max_reflections, 0, whole = TRUE),
              "`max_reflections` must be one whole number, 0 or more.")
  structure(list(start = start, goal = goal, goal_bias = goal_bias,
                 goal_tolerance = goal_tolerance, start_time = start_time,
                 max_reflections = as.integer(max_reflections)), class = "track_control")
}

check_utc <- function(x, what) {
  if (!inherits(x, "POSIXct") || !identical(attr(x, "tzone"), "UTC")) {
    stop_classed("input", "`", what, "` must be POSIXct with time zone \"UTC\".")
  }
  invisible()
}

#' Simulate animal tracks on a water raster
#'
#' Walks animals on a raster in which NA cells are land, with the step and
#' turn samplers of a [movement_model()]. A step that meets land is
#' reflected off the edge of the land cell, as light off a mirror, for the
#' rest of its length, so the animal follows the shore back into the water
#' and no part of the step crosses land. Reflection keeps the time spent in
#' each water cell in proportion to its area, near the shore as well as
#' offshore, for uniform and for persistent turns. A step holds its position
#' only when it needs more than `max_reflections` reflections, or when the
#' straight line from its start to its end would cross land; positions are
#' joined by straight lines downstream, so that line must be on water.
#'
#' Bearings are in radians, anticlockwise from east. The first bearing of a
#' track is uniform; after a reflection the walk keeps the reflected bearing.
#' With a two-state model the animal does not move while it rests, and keeps
#' its bearing for when it moves again.
#'
#' The function warns when more than 1 percent of moving steps held, when more than
#' half were reflected (often a step length in the wrong units, or a goal
#' pulling an animal against the shore), when more than 0.1 percent held after 1,000
#' rejected candidates under a [habitat_selection()] preference, or when a
#' goal was given and no animal reached it. With a goal, every start must be in
#' water that connects to the goal through cells that share an edge; the
#' function stops otherwise.
#'
#' @param grid A [terra::SpatRaster] in a projected CRS whose unit is the
#'   metre and whose scale at the raster centre, on its two axes and its
#'   diagonals, is within 1 percent of 1, as in [place_receivers()], with
#'   square cells, NA on land.
#' @param model A [movement_model()], with or without a habitat preference
#'   from [habitat_selection()].
#' @param n_animals Number of animals, each walked on its own track.
#' @param n_steps Steps per animal; each track has `n_steps + 1` positions.
#' @param control A [track_control()] list: start, goal, clock, reflections.
#' @param seed Seed for the simulation; `NULL` uses the current random number
#'   stream. The caller's stream is restored after a seeded call. The draws
#'   depend on the generator kind as well as the seed, as with [set.seed()];
#'   the kind is stored with the result.
#'
#' @return A data frame of class `sim_tracks`, one row per position, with
#'   columns `path_id` (integer), `step` (0 to `n_steps`), `time` (POSIXct,
#'   UTC), `x`, `y`, `cell` (the grid cell holding the position) and `state`
#'   (1 moving, 2 resting; the state of the step that ends at the position).
#'   Attributes: `model`, `control`, `seed`, `rng_kind` (the [RNGkind()] in
#'   force), `grid` (the raster, packed with [terra::wrap()]), `epsg`,
#'   `arrived` (step at which each animal reached the goal, NA if it did not)
#'   and `escapes`, a data frame per animal counting the steps taken as drawn
#'   (`clear`), reflected at the shore (`reflected`), held (`held`), and
#'   without movement because the animal rested or had reached the goal
#'   (`still`). The four counts of an animal sum to `n_steps`. A model from
#'   [habitat_selection()] adds a column `heading`, the bearing in radians a
#'   moving step started from before its turn (NA at step 0, at resting
#'   steps and at steps held after 1,000 rejected candidates), and an
#'   attribute `selection`, a data frame per animal counting the candidate
#'   steps drawn (`proposals`), the steps accepted (`accepted`), candidates
#'   longer than the reach (`over_reach`), candidates whose weight exceeded
#'   the bound (`above_bound`) and steps held after 1,000 rejections
#'   (`exhausted`, counted in `held` as well).
#' @seealso [simulate_transmissions()], [movement_model()]
#' @export
#' @examples
#' stoney <- terra::rast(system.file("extdata", "stoney_depth.tif",
#'                                   package = "arraydesign"))
#' # Three animals for 6 hours. evaluate_array() walks its animals with this
#' # model unless `control$movement` says otherwise.
#' walk <- crw_model(30, 15, 40, time_step_s = 60)
#' tr <- simulate_tracks(stoney, walk, n_animals = 3, n_steps = 360, seed = 1)
#' tr
#' attr(tr, "escapes")
#' plot(tr)
#'
#' # The same walk from the south-west arm, pulled towards a bay on the north
#' # shore for up to 15 hours; `arrived` gives the step of arrival.
#' to_bay <- track_control(start = data.frame(x = 721950, y = 4930450),
#'                         goal = data.frame(x = 725500,
#'                                           y = 4938000), goal_bias = 0.2,
#'                         goal_tolerance = 300)
#' tr_bay <- simulate_tracks(stoney, walk, n_animals = 3, n_steps = 900,
#'                           control = to_bay,
#'                           seed = 1)
#' attr(tr_bay, "arrived")
#' plot(tr_bay)
simulate_tracks <- function(grid, model, n_animals = 1, n_steps = 100, control = track_control(),
                            seed = NULL) {
  wg <- water_grid(grid, require_crs = TRUE, call = sys.call())
  stop_unless(inherits(model, "movement_model"),
              "`model` must be a movement_model; see movement_model().")
  stop_unless(is_scalar(n_animals, 1, whole = TRUE),
              "`n_animals` must be one whole number of at least 1.")
  stop_unless(is_scalar(n_steps, 1, whole = TRUE),
              "`n_steps` must be one whole number of at least 1.")
  stop_unless(inherits(control, "track_control"), "`control` must come from track_control().")
  goal <- goal_point(wg, control)
  reach <- reachable(wg, goal)
  # A model with a habitat preference walks by rejection. A model without one
  # skips the sampler, so its tracks for a given seed do not depend on it.
  sel <- if (is.null(model$selection)) NULL else selection_setup(wg, model)
  walks <- with_seed_or_not(seed, {
    starts <- start_points(wg, control$start, n_animals, reach)
    lapply(seq_len(n_animals), function(i) {
      if (!is.null(sel)) {
        return(walk_path_selection(wg, model, n_steps, starts[i, ], goal,
                                   control$max_reflections, sel))
      }
      walk_path(wg, model, n_steps, starts[i, ], goal, control$max_reflections)
    })
  })
  if (!is.null(sel)) model$selection[c("reach", "window")] <- sel[c("reach", "window")]
  out <- assemble_tracks(wg, walks, model, control, seed)
  warn_degenerate(out, goal, n_steps)
  out
}

# ---- Starts and goal ----
start_points <- function(wg, start, n_animals, reach) {
  if (is.null(start)) {
    cells <- which(reach)
    return(cell_xy(wg, cells[sample.int(length(cells), n_animals, replace = TRUE)]))
  }
  xy <- as_xy(start, "start")
  if (!nrow(xy) %in% c(1L, n_animals)) {
    stop_classed("track", "`start` must have 1 row or n_animals (", n_animals, ") rows; it has ",
                 nrow(xy), ".")
  }
  dry <- which(!is_wet(wg, xy[, 1], xy[, 2]))
  if (length(dry)) {
    rt <- rows_text(dry)
    stop_classed("track", "`start` has ", rt$count, " on land or outside the grid: ", rt$shown, ".")
  }
  cut_off <- which(!reach[cell_at(wg, xy[, 1], xy[, 2])])
  if (length(cut_off)) {
    rt <- rows_text(cut_off)
    stop_classed("track", "`start` has ", rt$count, " in water that does not connect to the goal: ",
                 rt$shown, ".")
  }
  xy[rep_len(seq_len(nrow(xy)), n_animals), , drop = FALSE]
}

goal_point <- function(wg, control) {
  if (is.null(control$goal)) return(NULL)
  xy <- as_xy(control$goal, "goal")
  if (nrow(xy) != 1L) stop_classed("track", "`goal` must be one point.")
  if (!is_wet(wg, xy[1, 1], xy[1, 2])) {
    stop_classed("track", "`goal` is on land or outside the grid.")
  }
  tol <- if (is.null(control$goal_tolerance)) wg$res else control$goal_tolerance
  list(x = xy[1, 1], y = xy[1, 2], bias = control$goal_bias, tolerance = tol)
}

# Water cells a path can share with the goal. A step cannot pass between two
# cells that touch only at a corner, so water connects through shared edges.
reachable <- function(wg, goal) {
  if (is.null(goal)) return(wg$wet)
  r <- terra::rast(wg$raster)
  terra::values(r) <- ifelse(wg$wet, 1L, NA_integer_)
  part <- terra::values(terra::patches(r, directions = 4), mat = FALSE)
  !is.na(part) & part == part[cell_at(wg, goal$x, goal$y)]
}

as_xy <- function(p, what) {
  if (is.numeric(p) && is.null(dim(p)) && length(p) == 2L) {
    p <- matrix(p, 1, dimnames = list(NULL, c("x", "y")))
  }
  stop_unless((is.data.frame(p) || is.matrix(p)) && all(c("x", "y") %in% colnames(p)),
              "`", what, "` must have columns `x` and `y`.")
  numeric_xy(p[, "x"], p[, "y"], what)
}

# ---- One path ----
walk_path <- function(wg, model, n_steps, start, goal, max_reflections) {
  len <- model$step(n_steps)
  turn <- model$turn(n_steps)
  stop_unless(valid_draws(len, n_steps, 0) && valid_draws(turn, n_steps, -Inf),
              "The model's samplers returned draws that are not valid: ",
              "n finite numbers, and steps of at least 0.")
  state <- state_sequence(model, n_steps)
  x <- y <- numeric(n_steps + 1L)
  x[1] <- start[1]
  y[1] <- start[2]
  bearing <- stats::runif(1, -pi, pi)
  outcome <- rep("still", n_steps)
  arrived <- NA_integer_
  for (k in seq_len(n_steps)) {
    x[k + 1] <- x[k]
    y[k + 1] <- y[k]
    if (state[k] == 2L || !is.na(arrived)) next
    b <- pull_to_goal(bearing + turn[k], x[k], y[k], goal)
    mv <- take_step(wg, x[k], y[k], b, len[k], max_reflections)
    x[k + 1] <- mv$x
    y[k + 1] <- mv$y
    bearing <- mv$bearing
    outcome[k] <- mv$outcome
    if (at_goal(goal, mv$x, mv$y)) arrived <- k
  }
  list(x = x, y = y, state = c(state[1], state), outcome = outcome, arrived = arrived)
}

# ---- One path with a habitat preference ----
# As walk_path(), but each moving step is drawn by rejection: candidates come
# from the same kernel (step, turn, goal pull, reflection) and are accepted
# with probability w(end) / W(start). `heading` is the bearing a step started
# from, before its turn; NA where the step drew no candidate or none was
# accepted.
walk_path_selection <- function(wg, model, n_steps, start, goal, max_reflections, sel) {
  state <- state_sequence(model, n_steps)
  x <- y <- numeric(n_steps + 1L)
  x[1] <- start[1]
  y[1] <- start[2]
  heading <- rep(NA_real_, n_steps + 1L)
  bearing <- stats::runif(1, -pi, pi)
  outcome <- rep("still", n_steps)
  tally <- selection_tally()
  arrived <- NA_integer_
  for (k in seq_len(n_steps)) {
    x[k + 1] <- x[k]
    y[k + 1] <- y[k]
    if (state[k] == 2L || !is.na(arrived)) next
    st <- selected_step(wg, model, c(x[k], y[k]), bearing, goal, max_reflections, sel)
    tally <- tally + st$tally
    if (is.null(st$move)) {
      outcome[k] <- "held"
      next
    }
    heading[k + 1] <- bearing
    x[k + 1] <- st$move$x
    y[k + 1] <- st$move$y
    bearing <- st$move$bearing
    outcome[k] <- st$move$outcome
    if (at_goal(goal, st$move$x, st$move$y)) arrived <- k
  }
  list(x = x, y = y, state = c(state[1], state), outcome = outcome, arrived = arrived,
       heading = heading, tally = tally)
}

selection_tally <- function() {
  c(proposals = 0, accepted = 0, over_reach = 0, above_bound = 0, exhausted = 0)
}

# One accepted step from `from`, or `move = NULL` after `sel$max_tries`
# rejections. Candidates are drawn in batches; each is one step length, one
# turn and one uniform, so the batch size does not change the distribution.
selected_step <- function(wg, model, from, bearing, goal, max_reflections, sel) {
  bound <- point_bound(wg, sel, from[1], from[2])
  tally <- selection_tally()
  batch <- 16L
  repeat {
    len <- model$step(batch)
    turn <- model$turn(batch)
    stop_unless(valid_draws(len, batch, 0) && valid_draws(turn, batch, -Inf),
                "The model's samplers returned draws that are not valid: ",
                "n finite numbers, and steps of at least 0.")
    u <- stats::runif(batch)
    for (i in seq_len(batch)) {
      tally[["proposals"]] <- tally[["proposals"]] + 1
      b <- pull_to_goal(bearing + turn[i], from[1], from[2], goal)
      mv <- take_step(wg, from[1], from[2], b, len[i], max_reflections)
      ratio <- exp(sel$log_w[cell_at(wg, mv$x, mv$y)] - bound)
      tally[["over_reach"]] <- tally[["over_reach"]] + (len[i] > sel$reach)
      tally[["above_bound"]] <- tally[["above_bound"]] + (ratio > 1)
      if (u[i] < ratio) {
        tally[["accepted"]] <- 1
        return(list(move = mv, tally = tally))
      }
      if (tally[["proposals"]] >= sel$max_tries) {
        tally[["exhausted"]] <- 1
        return(list(move = NULL, tally = tally))
      }
    }
  }
}

at_goal <- function(goal, x, y) {
  !is.null(goal) && sqrt((goal$x - x)^2 + (goal$y - y)^2) <= goal$tolerance
}

# Bearing pulled towards the goal by a weighted mean of the two unit vectors.
pull_to_goal <- function(b, x, y, goal) {
  if (is.null(goal)) return(b)
  g <- atan2(goal$y - y, goal$x - x)
  atan2((1 - goal$bias) * sin(b) + goal$bias * sin(g),
        (1 - goal$bias) * cos(b) + goal$bias * cos(g))
}

# One step from (x0, y0). A step that meets land is reflected off the cell
# edge it crosses, as light off a mirror, for the rest of its length, and
# again at each further edge. Reflection keeps the time spent in each water
# cell in proportion to its area for both uniform and persistent turns; the
# rejection and replacement rules of positionR 1.x left the shore under-used.
# The straight line from start to end must also be on water, because
# positions are joined by straight lines downstream; a step that fails that
# test, or needs more than `max_reflections`, holds and is counted.
take_step <- function(wg, x0, y0, b, len, max_reflections) {
  d <- c(len * cos(b), len * sin(b))
  p <- c(x0, y0)
  n <- 0L
  repeat {
    hit <- first_land(wg, p[1], p[2], p[1] + d[1], p[2] + d[2])
    if (is.null(hit)) break
    if (n >= max_reflections || is.na(hit$axis)) return(held_step(x0, y0, d))
    # Stop 1 micrometre short of the edge, so the next leg starts on water
    # at any coordinate size.
    t <- max(hit$t - 1e-6 / sqrt(sum(d^2)), 0)
    p <- p + t * d
    d <- (1 - t) * d
    d[hit$axis] <- -d[hit$axis]
    n <- n + 1L
  }
  end <- p + d
  if (n == 0L) return(list(x = end[1], y = end[2], bearing = b, outcome = "clear"))
  if (!segment_wet(wg, x0, y0, end[1], end[2])) return(held_step(x0, y0, d))
  list(x = end[1], y = end[2], bearing = atan2(d[2], d[1]), outcome = "reflected")
}

# A held step keeps the heading the reflections ended on, which points back
# into the water, so the next step does not meet the same edge head on.
held_step <- function(x0, y0, d) {
  list(x = x0, y = y0, bearing = atan2(d[2], d[1]), outcome = "held")
}

# State of each step: 1 moving, 2 resting. The first state is drawn from the
# stationary distribution of the chain, so a path does not start biased.
state_sequence <- function(model, n_steps) {
  if (model$states == 1L) return(rep(1L, n_steps))
  p <- model$transition
  s <- integer(n_steps)
  s[1] <- if (stats::runif(1) < p[2, 1] / (p[1, 2] + p[2, 1])) 1L else 2L
  u <- stats::runif(n_steps)
  for (k in seq_len(n_steps)[-1]) s[k] <- if (u[k] < p[s[k - 1], 1]) 1L else 2L
  s
}

# ---- Output ----
assemble_tracks <- function(wg, walks, model, control, seed) {
  n1 <- length(walks[[1]]$x)
  out <- data.frame(
    path_id = rep(seq_along(walks), each = n1),
    step = rep(seq_len(n1) - 1L, times = length(walks)),
    x = unlist(lapply(walks, `[[`, "x")), y = unlist(lapply(walks, `[[`, "y"))
  )
  out$time <- control$start_time + out$step * model$time_step
  out$cell <- as.integer(cell_at(wg, out$x, out$y))
  out$state <- unlist(lapply(walks, `[[`, "state"))
  out <- out[c("path_id", "step", "time", "x", "y", "cell", "state")]
  if (!is.null(walks[[1]]$heading)) out$heading <- unlist(lapply(walks, `[[`, "heading"))
  levels <- c("clear", "reflected", "held", "still")
  esc <- t(vapply(walks, function(w) table(factor(w$outcome, levels))[levels], numeric(4)))
  escapes <- data.frame(path_id = seq_along(walks), esc)
  out <- structure(
    out, class = c("sim_tracks", "data.frame"), model = model, control = control, seed = seed,
    rng_kind = RNGkind(), grid = terra::wrap(wg$raster), epsg = grid_epsg(wg$raster),
    arrived = vapply(walks, `[[`, integer(1), "arrived"), escapes = escapes
  )
  if (!is.null(walks[[1]]$tally)) {
    tally <- do.call(rbind, lapply(walks, `[[`, "tally"))
    attr(out, "selection") <- data.frame(path_id = seq_along(walks), tally)
  }
  out
}

# A warning for walks the shoreline or the goal shaped more than the model.
warn_degenerate <- function(tracks, goal, n_steps) {
  esc <- attr(tracks, "escapes")
  moving <- sum(esc$clear, esc$reflected, esc$held)
  held <- sum(esc$held)
  reflected <- sum(esc$reflected)
  msg <- character()
  if (held > 0.01 * moving) {
    msg <- c(msg, sprintf("%d of %d moving steps (%.1f%%) held at the shore; %s", held, moving,
                          100 * held / moving,
                          "steps may be long for this shoreline at this cell size."))
  }
  if (reflected > 0.5 * moving) {
    msg <- c(msg, sprintf("%d of %d moving steps (%.1f%%) were reflected at the shore; %s",
                          reflected, moving, 100 * reflected / moving,
                          "check the step length units and any goal."))
  }
  exhausted <- sum(attr(tracks, "selection")$exhausted)
  if (exhausted > 0.001 * moving) {
    msg <- c(msg, sprintf("%d of %d moving steps held after %s; %s", exhausted, moving,
                          "1,000 rejected candidates",
                          "the preference is strong for this step length and layer."))
  }
  if (!is.null(goal) && all(is.na(attr(tracks, "arrived")))) {
    msg <- c(msg, sprintf("No animal of %d reached the goal in %d steps.", nrow(esc), n_steps))
  }
  if (length(msg)) warning(paste(msg, collapse = "\n"), call. = FALSE)
  invisible()
}
