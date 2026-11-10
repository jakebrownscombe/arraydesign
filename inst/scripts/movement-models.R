# Code of the vignette "Specifying and testing movement models".
# Extracted from vignettes/movement-models.Rmd by tools/purl_vignettes.R, a repository file
# the installed package does not carry. The vignette is the source; this file follows it.

# ---- Read the lake ----
# Load the package, and ggplot2 for the maps
library(arraydesign)
library(ggplot2)

# Read the depth raster of Stoney Lake
lake <- terra::rast(system.file("extdata", "stoney_depth.tif",
                                package = "arraydesign"))

# ---- Depth map helper ----
# Take the water cells, with depth as positive metres
water <- terra::as.data.frame(lake, xy = TRUE, na.rm = TRUE)
names(water)[3] <- "depth_m"
water$depth_m <- abs(water$depth_m)

# Draw water by depth, land blank; `zoom` takes the x and y limits of a
# part of the lake, and NULL draws all of it
depth_map <- function(zoom = NULL) {
  ggplot(water, aes(x, y)) +
    geom_raster(aes(fill = depth_m)) +
    detfield::scale_fill_depth() +
    coord_equal(xlim = zoom$x, ylim = zoom$y, expand = FALSE) +
    labs(x = "Easting (m)", y = "Northing (m)") +
    detfield::annotation_scale_bar(location = "br") +
    detfield::theme_positionr()
}

# ---- Correlated random walk ----
# Build a walk: gamma steps of mean 30 m and sd 15 m per 60 s, turns of
# sd 40 degrees
walk <- crw_model(step_mean = 30, step_sd = 15, turn_sd = 40, time_step_s = 60)
walk

# ---- Explore step lengths ----
# Draw 10,000 step lengths from the walk, m per 60 s
set.seed(1)
step_draws <- walk$step(10000)
summary(step_draws)

# ---- Tracks by turning spread ----
# Set a release point in the middle of the lake
middle <- data.frame(x = 727850, y = 4937250)

# Set one turning spread per walk, in degrees; Inf gives uniform turns
turn_sds <- c(ten = 10, forty = 40, uniform = Inf)

# Walk five animals for 240 steps of 60 s under each walk
turn_tracks <- lapply(turn_sds, function(turn_sd) {
  spread_walk <- crw_model(step_mean = 30, step_sd = 15, turn_sd = turn_sd,
                           time_step_s = 60)
  simulate_tracks(lake, spread_walk, n_animals = 5, n_steps = 240,
                  control = track_control(start = middle), seed = 1)
})

# ---- Map turning spread ----
# Name the panels, in the order of turn_sds
panel_names <- c(ten = "(a) 10 degrees", forty = "(b) 40 degrees",
                 uniform = "(c) uniform turns")

# Stack the three runs in one table, with the panel each is drawn in
spread_paths <- do.call(rbind, lapply(names(turn_tracks), function(spread) {
  data.frame(path_id = turn_tracks[[spread]]$path_id,
             x = turn_tracks[[spread]]$x, y = turn_tracks[[spread]]$y,
             panel = panel_names[[spread]])
}))

# Draw one panel per turning spread over the middle of the lake, each
# track over a white casing
depth_map(zoom = list(x = c(724500, 731500), y = c(4934800, 4939800))) +
  geom_path(data = spread_paths, aes(group = path_id), colour = "white",
            linewidth = 1.2) +
  geom_path(data = spread_paths, aes(group = path_id, colour = factor(path_id)),
            linewidth = 0.5) +
  geom_point(data = middle, shape = 21, fill = "white", colour = "black",
             size = 2, stroke = 0.6) +
  scale_colour_tracks(name = "Animal") +
  facet_wrap(~panel, ncol = 2) +
  theme(legend.position = "inside", legend.position.inside = c(0.75, 0.25),
        legend.box = "horizontal")

# ---- Count shore reflections ----
# Count the steps per animal: as drawn, reflected, held and still
attr(turn_tracks$uniform, "escapes")

# ---- Two state samplers ----
# Set the moving-state estimates: mean and sd of step length (m per 60 s)
# and the concentration of the turning angles
hmm_step_mean <- 25
hmm_step_sd <- 15
hmm_rho <- 0.6

# Draw gamma step lengths with that mean and sd
gamma_steps <- function(n_draws) {
  stats::rgamma(n_draws, shape = (hmm_step_mean / hmm_step_sd)^2,
                rate = hmm_step_mean / hmm_step_sd^2)
}

# Draw wrapped Cauchy turning angles with concentration rho
wrapped_cauchy <- function(n_draws) {
  (stats::rcauchy(n_draws, 0, -log(hmm_rho)) + pi) %% (2 * pi) - pi
}

# ---- Two state model ----
# Set the transition matrix: rows are the state at one step, columns the
# state at the next (1 moving, 2 resting)
move_rest <- matrix(c(0.95, 0.05,    # moving: keep moving 0.95, rest 0.05
                      0.10, 0.90),   # resting: move 0.10, keep resting 0.90
                    nrow = 2, byrow = TRUE)

# Move with the samplers (state 1) or rest (state 2)
hmm_model <- movement_model(
  step = gamma_steps, turn = wrapped_cauchy, time_step_s = 60,
  states = 2, transition = move_rest,
  source = paste("HMM, illustrative values: gamma steps 25 m (sd 15)",
                 "per 60 s, wrapped Cauchy turns rho 0.6, and rest")
)
hmm_model

# ---- Two state track ----
# Walk one animal for 480 steps of 60 s from the middle of the lake
two_state <- simulate_tracks(lake, hmm_model, n_animals = 1, n_steps = 480,
                             control = track_control(start = middle), seed = 4)

# Count the positions in each state
table(two_state$state)

# ---- Resting places ----
# Take every position of the track, and those at which the animal rested
positions <- data.frame(x = two_state$x, y = two_state$y)
resting <- positions[two_state$state == 2, ]

# Count the resting positions, and the places at which the animal rested
nrow(resting)
nrow(unique(resting))

# ---- Map resting positions ----
# Map the track over a white casing in the moving colour, the resting
# positions on top as squares in the resting colour, and the start, zoomed
# to the track
depth_map(zoom = list(x = range(two_state$x) + c(-800, 800),
                      y = range(two_state$y) + c(-800, 800))) +
  geom_path(data = positions, colour = "white", linewidth = 1.4) +
  geom_path(data = positions, aes(colour = "Moving"), linewidth = 0.6) +
  geom_point(data = resting, aes(colour = "Resting"), shape = 15, size = 1.8) +
  geom_point(data = positions[1, ], shape = 21, fill = "white",
             colour = "black", size = 2.2, stroke = 0.6) +
  scale_colour_tracks(name = "State",
                      colours = c(Moving = "bluish_green", Resting = "orange")) +
  guides(colour = guide_legend(override.aes = list(linetype = c(1, 0),
                                                   shape = c(NA, 15))))

# ---- Own steps sampler ----
# Enter twelve measured steps, m per 180 s
my_steps <- c(12, 34, 8, 56, 120, 44, 19, 73, 33, 251, 62, 28)

# Resample the steps; turn up to 90 degrees either way
own_model <- movement_model(
  step = function(n_draws) sample(my_steps, n_draws, replace = TRUE),
  turn = function(n_draws) stats::runif(n_draws, -pi / 2, pi / 2),
  time_step_s = 180,
  source = "twelve example steps per 180 s, resampled; turns up to 90 degrees"
)
own_model

# ---- Sampler error ----
# Try a step sampler that ignores n_draws and draws one value
try(movement_model(
  step = function(n_draws) stats::rgamma(1, shape = 2, rate = 0.1),
  turn = function(n_draws) stats::runif(n_draws, -pi, pi)
))

# ---- Grid array ----
# Lay the 20-receiver grid of Get started over water at least 1 m deep
deployable <- terra::ifel(abs(lake) < 1, NA, lake)
receivers <- place_receivers(deployable, n = 20, method = "grid", seed = 1)

# Set the assumed range model: d50 and d05 for receivers 1 m and 35 m deep
range_model <- detfield::de_logistic(d50 = c(250, 750), d05 = c(500, 1500),
                                     depth = c(1, 35))

# ---- Release and goal ----
# Set a release point in the south-west arm and a goal on the north shore
release <- data.frame(x = 721950, y = 4930450)
goal <- data.frame(x = 725500, y = 4938000)
start_time <- as.POSIXct("2025-06-01 06:00:00", tz = "UTC")

# Compute the straight-line distance from the release to the goal, m
sqrt((goal$x - release$x)^2 + (goal$y - release$y)^2)

# ---- Head for a goal ----
# Set the goal biases to compare
biases <- c(a = 0.2, b = 0.3, c = 0.5)

# Walk ten animals for 1,440 steps of 60 s at each bias
goal_tracks <- lapply(biases, function(bias) {
  control <- track_control(start = release, goal = goal, goal_bias = bias,
                           goal_tolerance = 300, start_time = start_time)
  simulate_tracks(lake, walk, n_animals = 10, n_steps = 1440, control = control,
                  seed = 1)
})

# Print the first and last time of the run at a bias of 0.2
range(goal_tracks$a$time)

# ---- Arrival steps ----
# Take the step of arrival of each animal, at each bias
arrived <- lapply(goal_tracks, attr, "arrived")
arrived

# ---- Reflected share ----
# Compute the share of moving steps reflected at the shore, per bias
reflected_share <- sapply(goal_tracks, function(tracks) {
  escapes <- attr(tracks, "escapes")
  sum(escapes$reflected) / sum(escapes$clear, escapes$reflected, escapes$held)
})
round(reflected_share, 3)

# ---- Map goal tracks ----
# Draw a circle of 300 m round the goal, the arrival distance
angle <- seq(0, 2 * pi, length.out = 100)
arrival_ring <- data.frame(x = goal$x + 300 * cos(angle),
                           y = goal$y + 300 * sin(angle))

# Stack the tracks of the three biases, with the panel each is drawn in
bias_paths <- do.call(rbind, lapply(names(biases), function(run) {
  data.frame(as.data.frame(goal_tracks[[run]])[c("path_id", "x", "y")],
             bias = as.character(biases[[run]]),
             panel = sprintf("(%s) goal bias %s", run, biases[[run]]))
}))

# Map the tracks over the west half of the lake, one panel and one colour
# per bias, each track over a white casing
goal_map <- depth_map(zoom = list(x = c(720300, 729500),
                                  y = c(4929240, 4939200))) +
  geom_path(data = bias_paths, aes(group = path_id), colour = "white",
            linewidth = 0.8) +
  geom_path(data = bias_paths, aes(group = path_id, colour = bias),
            linewidth = 0.4) +
  scale_colour_tracks(name = "Goal bias",
                      colours = c("0.2" = "bluish_green", "0.3" = "orange",
                                  "0.5" = "reddish_purple")) +
  geom_point(data = release, shape = 21, fill = "white", colour = "black",
             size = 2, stroke = 0.6) +
  geom_path(data = arrival_ring, linetype = "dashed", colour = "black") +
  geom_point(data = goal, shape = 23, size = 2, stroke = 0.8,
             colour = "black", fill = "white") +
  facet_wrap(~panel, ncol = 3) +
  theme(legend.position = "bottom")
goal_map

# ---- Depth preference ----
# Make a layer of bottom depth in positive metres, named for the analysis
depth <- abs(lake)
names(depth) <- "depth_m"

# Give the walk of section 1 a preference for shallow water: a log weight
# of -0.3 per m of bottom depth
shallow_walk <- habitat_selection(walk, depth, beta = -0.3)
shallow_walk

# ---- Walk with preference ----
# Draw 20 release cells on water, used for every run of this section
set.seed(20)
releases <- water[sample(nrow(water), 20), c("x", "y")]

# Walk 20 animals for 2,880 steps of 60 s (48 hours) from those cells
shallow_tracks <- simulate_tracks(lake, shallow_walk, n_animals = 20,
                                  n_steps = 2880,
                                  control = track_control(start = releases),
                                  seed = 1)
shallow_tracks

# ---- Walk three preferences ----
# Take the first five releases in the eastern part of the lake, the area
# of Figure 5, which holds the north-east basin and the deepest water
basin <- list(x = c(727000, 737400), y = c(4935600, 4941900))
basin_releases <- releases[releases$x > basin$x[1] &
                             releases$y > basin$y[1], ][1:5, ]

# Set three preferences for bottom depth, per m
betas <- c(shallow = -0.3, none = 0, deep = 0.3)

# Walk the five animals for 48 hours under each preference
preference_tracks <- lapply(betas, function(beta) {
  simulate_tracks(lake, habitat_selection(walk, depth, beta = beta),
                  n_animals = 5, n_steps = 2880,
                  control = track_control(start = basin_releases), seed = 1)
})

# Median bottom depth of the positions of each run, m
depth_values <- terra::values(depth, mat = FALSE)
used_depth <- sapply(preference_tracks, function(tracks) {
  median(depth_values[tracks$cell])
})
round(used_depth, 2)

# ---- Map three preferences ----
# Name the panels, in the order of betas
preference_names <- c(shallow = "(a) beta -0.3 per m, shallow",
                      none = "(b) beta 0, no preference",
                      deep = "(c) beta 0.3 per m, deep")

# Stack the three runs, with the panel each is drawn in
preference_paths <- do.call(rbind, lapply(
  names(preference_tracks),
  function(run) {
    data.frame(path_id = preference_tracks[[run]]$path_id,
               x = preference_tracks[[run]]$x, y = preference_tracks[[run]]$y,
               preference = run, panel = preference_names[[run]])
  }
))

# Draw one panel and one colour per preference over the eastern part of
# the lake, each track over a white casing
depth_map(zoom = basin) +
  geom_path(data = preference_paths, aes(group = path_id), colour = "white",
            linewidth = 0.6) +
  geom_path(data = preference_paths,
            aes(group = path_id, colour = preference), linewidth = 0.3) +
  geom_point(data = basin_releases, shape = 21, fill = "white",
             colour = "black", size = 1.6, stroke = 0.5) +
  scale_colour_tracks(name = "Preference",
                      colours = c(shallow = "orange", none = "bluish_green",
                                  deep = "reddish_purple"),
                      breaks = names(betas),
                      labels = c("prefers shallow", "no preference",
                                 "prefers deep")) +
  facet_wrap(~panel, ncol = 2) +
  theme(legend.position = "inside", legend.position.inside = c(0.75, 0.25),
        legend.box = "horizontal")

# ---- Draw available steps ----
# Draw 10 control steps for each step of the 20 animals
shallow_steps <- available_steps(shallow_tracks, n = 10, seed = 2)
head(shallow_steps)

# ---- Fit per animal ----
# Attach survival; clogit() reads a bare strata() term
library(survival)

# Fit a conditional logit per animal: step taken (case 1) against controls
per_animal <- do.call(rbind, lapply(
  split(shallow_steps, shallow_steps$path_id),
  function(steps) {
    fit <- clogit(case ~ depth_m + strata(stratum), data = steps)
    data.frame(path_id = steps$path_id[1], estimate = unname(coef(fit)),
               se = sqrt(unname(diag(vcov(fit)))))
  }
))
head(per_animal)

# ---- Compare with beta ----
# Mean across animals and its 95% interval, t on 19 degrees of freedom
mean_estimate <- mean(per_animal$estimate)
half_width <- qt(0.975, 19) * sd(per_animal$estimate) / sqrt(20)
round(c(beta = -0.3, mean = mean_estimate,
        lower = mean_estimate - half_width,
        upper = mean_estimate + half_width), 3)

# Count the animals whose 95% interval holds beta
covered <- sum(abs(per_animal$estimate + 0.3) <=
                 qnorm(0.975) * per_animal$se)
covered

# ---- Score by movement model ----
# Build a slow walk: 6 m a minute
slow_walk <- crw_model(step_mean = 6, step_sd = 4, turn_sd = 40,
                       time_step_s = 60)

# Name the four movement models
models <- list(slow = slow_walk, crw = walk, two_state = hmm_model,
               own_steps = own_model)

# Evaluate the grid under each
evaluations <- lapply(models, function(model) {
  evaluate_array(receivers, lake, range_model, n_animals = 20,
                 seed = 1,
                 control = list(hours = 6, tag_delay = c(120, 180),
                                movement = model))
})

# ---- Rate per animal ----
# Tabulate the detection rate, one row per model and animal
rates <- do.call(rbind, lapply(names(evaluations), function(model) {
  sent <- evaluations[[model]]$detections$transmissions
  rate <- tapply(sent$n_detections > 0, sent$path_id, mean)
  data.frame(movement = model, animal_id = names(rate),
             rate = as.vector(rate))
}))
head(rates)

# ---- Compare the models ----
# Tabulate detection, positioning error and hit rate, one row per model
by_model <- do.call(rbind, lapply(names(evaluations), function(model) {
  overall <- evaluations[[model]]$score$overall
  model_rates <- rates$rate[rates$movement == model]
  data.frame(movement = model, mean_rate = round(mean(model_rates), 3),
             never_heard = sum(model_rates == 0),
             median_error_m = round(overall$median_error_mean_m),
             hit_rate = round(overall$hit_rate, 3))
}))
by_model

# ---- Plot rate by movement model ----
# Name each model in words on the axis
model_names <- c(slow = "walk, 6 m a minute", crw = "walk, 30 m a minute",
                 two_state = "move and rest", own_steps = "resampled steps")
rates$movement <- factor(model_names[rates$movement], levels = model_names)

# Compute the mean rate per model, to draw as a bar
means <- aggregate(rate ~ movement, data = rates, FUN = mean)

# Draw one point per animal, spread sideways, and the mean per model
ggplot(rates, aes(movement, rate)) +
  geom_point(position = position_jitter(width = 0.15, height = 0, seed = 1),
             colour = "#21908C", alpha = 0.7, size = 1.6) +
  geom_errorbar(data = means, aes(ymin = rate, ymax = rate), width = 0.4,
                colour = "grey20", linewidth = 0.8) +
  scale_y_continuous(labels = function(share) sprintf("%g%%", 100 * share),
                     limits = c(0, 1)) +
  labs(x = "Movement model", y = "Share of transmissions heard") +
  detfield::theme_positionr()
