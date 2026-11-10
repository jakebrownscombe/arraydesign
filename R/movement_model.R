#' A movement model for the walker
#'
#' Holds what [simulate_tracks()] needs to move an animal: a sampler of step
#' lengths, a sampler of turning angles, the time one step takes, and,
#' optionally, a switch between a moving state and a resting state. Any pair
#' of samplers can be supplied, so a model fitted elsewhere, such as the step
#' and turn distributions of a hidden Markov model, can drive the walker.
#' [crw_model()] builds the common case, a correlated random walk.
#'
#' The walker moves by a sequence of step lengths and turning angles, each
#' drawn from its own probability distribution. The move and rest switch below
#' is a two-state Markov chain between steps, a switch between movement states
#' by a transition matrix as in Morales et al. (2004).
#'
#' With `states = 2`, state 1 moves with the two samplers and state 2 rests:
#' the animal stays where it is. The state changes between steps by the
#' transition matrix, so the number of steps spent in state `i` before a
#' change is geometric with mean `1 / (1 - transition[i, i])`.
#'
#' @param step A function of `n` that returns `n` step lengths in metres, each
#'   the straight-line distance covered in one time step.
#' @param turn A function of `n` that returns `n` turning angles in radians,
#'   the change of bearing between one step and the next.
#' @param states 1 for one moving state, or 2 for a move and rest switch.
#' @param transition For `states = 2`, a 2 by 2 matrix of transition
#'   probabilities between steps. Rows are the state at one step and columns
#'   the state at the next, so `transition[i, j]` is the probability of state
#'   `j` after a step in state `i`, and each row sums to 1. State 1 is moving
#'   and state 2 resting; `transition[2, 2]` is the probability of resting
#'   again after a rest.
#' @param time_step_s Seconds per step.
#' @param source Where the samplers come from, in words. It is printed with the
#'   model and kept with every simulation made from it.
#'
#' @return An object of class `movement_model`: a list with elements `step`,
#'   `turn`, `states`, `transition`, `time_step` (seconds per step) and
#'   `source`.
#' @references
#' Morales, J. M., Haydon, D. T., Frair, J., Holsinger, K. E. and Fryxell,
#' J. M. (2004). Extracting more out of relocation data: building movement
#' models as mixtures of random walks. *Ecology*, 85, 2436-2445.
#' \doi{10.1890/03-0269}
#' (Switching between behavioural states by a state transition matrix.)
#' @seealso [crw_model()], [simulate_tracks()]
#' @export
#' @examples
#' # Steps of 20 to 40 m per minute and turns of up to 30 degrees either way.
#' m <- movement_model(
#'   step = function(n) stats::runif(n, 20, 40),
#'   turn = function(n) stats::runif(n, -pi / 6, pi / 6),
#'   time_step_s = 60, source = "uniform steps and turns, for illustration"
#' )
#' m
#'
#' # The same walk, resting for a mean of 10 steps after a mean of 20 moving.
#' movement_model(
#'   step = function(n) stats::runif(n, 20, 40),
#'   turn = function(n) stats::runif(n, -pi / 6, pi / 6),
#'   states = 2, transition = matrix(c(0.95, 0.10, 0.05, 0.90), 2)
#' )
movement_model <- function(step, turn, states = 1, transition = NULL, time_step_s = 60,
                           source = "") {
  check_sampler(step, "step", lower = 0)
  check_sampler(turn, "turn", lower = -Inf)
  stop_unless(is_scalar(states) && states %in% 1:2, "`states` must be 1 or 2.")
  check_transition(transition, states)
  stop_unless(is_scalar(time_step_s, 1e-9, .Machine$double.xmax),
              "`time_step_s` must be one positive number of seconds.")
  stop_unless(is.character(source) && length(source) == 1L,
              "`source` must be one character string.")
  structure(
    list(step = step, turn = turn, states = as.integer(states), transition = transition,
         time_step = time_step_s, source = source),
    class = "movement_model"
  )
}

# A sampler must return n finite numbers at or above `lower`. Calling it once
# here catches a wrong sampler before a long simulation starts.
check_sampler <- function(f, what, lower) {
  stop_unless(is.function(f), "`", what, "` must be a function of n.")
  stop_unless(valid_draws(f(5L), 5L, lower), "`", what, "`(5) must return 5 finite numbers",
              if (lower == 0) " of at least 0", ".")
}

# TRUE when `draw` holds n finite numbers, none below `lower`.
valid_draws <- function(draw, n, lower) {
  is.numeric(draw) && length(draw) == n && all(is.finite(draw)) && all(draw >= lower)
}

check_transition <- function(transition, states) {
  if (states == 1) {
    stop_unless(is.null(transition), "`transition` is for `states = 2` only.")
    return(invisible())
  }
  shape_ok <- is.matrix(transition) && is.numeric(transition) && all(dim(transition) == 2L)
  stop_unless(shape_ok && valid_draws(transition, 4L, 0) &&
                all(abs(rowSums(transition) - 1) < 1e-8),
              "`transition` must be a 2 by 2 matrix of probabilities whose rows sum to 1.")
  stop_unless(all(diag(transition) < 1), "`transition` has an absorbing state (a diagonal of 1).")
}

#' A correlated random walk with gamma steps
#'
#' Builds a [movement_model()] with gamma-distributed step lengths and turning
#' angles from a wrapped normal distribution centred on zero, so the animal
#' tends to keep its heading. The turning spread is the input to vary in a
#' sensitivity analysis: `turn_sd = Inf` gives a walk with no memory of its
#' heading, and small values give straight paths.
#'
#' @param step_mean,step_sd Mean and standard deviation of the step length in
#'   metres per time step.
#' @param turn_sd Standard deviation of the turning angle in degrees, before
#'   wrapping to the circle; `Inf` for uniform turns.
#' @param time_step_s Seconds per step.
#' @param transition `NULL` for one moving state, or a 2 by 2 transition matrix
#'   for a move and rest switch, as in [movement_model()].
#'
#' @return A `movement_model`.
#' @seealso [movement_model()], [simulate_tracks()]
#' @export
#' @examples
#' crw_model(step_mean = 30, step_sd = 10, turn_sd = 30)
#' crw_model(30, 10, 30, transition = matrix(c(0.95, 0.10, 0.05, 0.90), 2))
crw_model <- function(step_mean, step_sd, turn_sd, time_step_s = 60, transition = NULL) {
  stop_unless(is_scalar(step_mean, 1e-9, .Machine$double.xmax) &&
                is_scalar(step_sd, 1e-9, .Machine$double.xmax),
              "`step_mean` and `step_sd` must be positive numbers.")
  shape <- (step_mean / step_sd)^2
  rate <- step_mean / step_sd^2
  movement_model(
    step = function(n) stats::rgamma(n, shape = shape, rate = rate),
    turn = wrapped_normal(turn_sd),
    states = if (is.null(transition)) 1 else 2, transition = transition,
    time_step_s = time_step_s,
    source = sprintf(paste(
      "crw_model(): gamma steps, mean %g m and sd %g m per %g s;",
      "wrapped normal turns, sd %g degrees"
    ), step_mean, step_sd, time_step_s, turn_sd)
  )
}

# Sampler of turning angles from a normal distribution with sd in degrees,
# wrapped to [-pi, pi).
wrapped_normal <- function(turn_sd) {
  stop_unless(is_scalar(turn_sd, 0),
              "`turn_sd` must be one number of degrees, 0 or more, or Inf.")
  if (is.infinite(turn_sd)) return(function(n) stats::runif(n, -pi, pi))
  sd_rad <- turn_sd * pi / 180
  function(n) (stats::rnorm(n, 0, sd_rad) + pi) %% (2 * pi) - pi
}

#' @export
print.movement_model <- function(x, ...) {
  cat("<movement_model>", x$states, if (x$states == 1) "state" else "states (move, rest)", "\n")
  cat("time step:", x$time_step, "s\n")
  if (x$states == 2) {
    dwell <- 1 / (1 - diag(x$transition))
    cat(sprintf("mean steps in a state: move %.1f, rest %.1f\n", dwell[1], dwell[2]))
  }
  if (!is.null(x$selection)) cat_wrapped(selection_text(x$selection))
  cat_wrapped(paste("source:", if (nzchar(x$source)) x$source else "not stated"))
  invisible(x)
}
