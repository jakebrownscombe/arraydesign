# Code of the vignette "Design an array for your own water body: compare layouts".
# Extracted from vignettes/design-an-array.Rmd by tools/purl_vignettes.R, a repository file
# the installed package does not carry. The vignette is the source; this file follows it.

# ---- Read your files ----
# Load the package, and ggplot2 for the figures
library(arraydesign)
library(ggplot2)

# Read your depth raster and candidate stations
lake <- terra::rast(system.file("extdata", "stoney_depth.tif",
                                package = "arraydesign"))
planned <- read.csv(system.file("extdata", "planned_stations.csv",
                                package = "arraydesign"))

# Take the range test: one row per receiver, tag and hour
range_test <- detfield::stoney_range_test

# Check the raster's coordinate system and cell size, m
terra::crs(lake, describe = TRUE)$code
terra::res(lake)

# Check the stations, and the receiver depths of the range test, m
head(planned, 3)
table(range_test$depth_m)

# ---- Your own range model ----
# Fit a logistic regression of detections on distance and receiver depth
range_fit <- glm(cbind(heard, sent - heard) ~ distance_m * depth_m,
                 family = binomial, data = range_test)

# Pass the fit to detfield as the range model of this page
range_model <- detfield::de_model(range_fit)
range_model

# ---- Compare the two models ----
# Set the stand-in of Get started: d50 and d05 for receivers 1 m and 35 m deep
stand_in <- detfield::de_logistic(d50 = c(250, 750), d05 = c(500, 1500),
                                  depth = c(1, 35))

# Read the detection probability of each model at 250, 500 and 750 m
distances <- c(250, 500, 750)
probability <- do.call(rbind, lapply(c(1, 10, 25), function(depth) {
  data.frame(depth_m = depth, distance_m = distances,
             fitted = round(plogis(range_model$predict(distances,
                                                       depth_m = depth)), 3),
             stand_in = round(plogis(stand_in$predict(distances,
                                                      depth_m = depth)), 3))
}))
probability

# ---- Largest difference ----
# Find the largest difference in detection probability, and where, per depth
fine_distances <- seq(0, 1500, by = 10)
largest <- do.call(rbind, lapply(c(1, 10, 25), function(depth) {
  gap <- abs(plogis(range_model$predict(fine_distances, depth_m = depth)) -
               plogis(stand_in$predict(fine_distances, depth_m = depth)))
  data.frame(depth_m = depth, largest = round(max(gap), 3),
             at_m = fine_distances[which.max(gap)])
}))
largest

# ---- Deployable water ----
# Set water shallower than 1 m to NA
deployable <- terra::ifel(abs(lake) < 1, NA, lake)

# ---- Lay out planned array ----
# Place the planned stations on the raster
planned_receivers <- place_receivers(lake, method = "manual",
                                     coords = planned)
head(planned_receivers, 3)

# ---- Station off the raster ----
# Copy the planned stations with a slip in the longitude of row 3
stations_with_slip <- planned
stations_with_slip$deploy_lon[3] <- -78.29706

# Try to place them; place_receivers() refuses the row and names it
try(place_receivers(lake, method = "manual", coords = stations_with_slip))

# ---- Lay out other arrays ----
# Name the three candidate arrays
designs <- list(
  planned = planned_receivers,
  random_20 = place_receivers(deployable, n = 20, method = "random",
                              spacing = 1000, seed = 1),
  grid_30 = place_receivers(deployable, n = 30, method = "grid", seed = 1)
)

# Count the receivers in each
sapply(designs, nrow)

# ---- Deepest receivers ----
# Give the depth of the deepest receiver in each design, m
round(sapply(designs, function(layout) max(layout$deploy_z)), 1)

# ---- Grid by spacing ----
# Count the receivers in a grid 2 km apart
nrow(place_receivers(deployable, spacing = 2000))

# ---- Build detection fields ----
# Build one detection field per design
fields <- lapply(designs, function(layout) {
  detfield::de_field(detfield::water_distances(lake, layout), range_model)
})

# ---- Coverage table ----
# Count the coverage of each design, one row per design and threshold
coverage <- do.call(rbind, lapply(names(fields), function(design) {
  cbind(design = design, coverage_summary(fields[[design]]))
}))
coverage

# ---- Array coverage ----
# Name the panels, in the order of `designs`
panel_names <- c(planned = "(a) planned, 20", random_20 = "(b) 20 at random",
                 grid_30 = "(c) grid of 30")

# Compute the probability that any receiver hears, per design and cell
array_cover <- do.call(rbind, lapply(names(fields), function(design) {
  data.frame(design = panel_names[[design]],
             p = 1 - apply(1 - fields[[design]]$p, 1, prod),
             terra::xyFromCell(lake, fields[[design]]$cells))
}))

# Give the share of water cells at p >= 0.5, per design
tapply(array_cover$p >= 0.5, array_cover$design, mean)

# ---- Map candidate coverage ----
# Stack the receivers of the designs, with the panel each is drawn in
design_rings <- do.call(rbind, lapply(names(designs), function(design) {
  data.frame(x = designs[[design]]$x, y = designs[[design]]$y,
             design = panel_names[[design]])
}))

# Draw one panel per design: the array's probability, and the receivers on
# top as a dark ring under a coloured one, keyed by the shape scale
ring <- detfield::positionr_ring_colours()
coverage_map <- ggplot(array_cover, aes(x, y)) +
  geom_raster(aes(fill = p)) +
  detfield::scale_fill_probability(name = "Detection\nprobability") +
  geom_point(data = design_rings, aes(shape = "Receiver"),
             colour = ring[["edge"]], size = 0.9, stroke = 0.7) +
  geom_point(data = design_rings, aes(shape = "Receiver"),
             colour = ring[["heard"]], size = 0.9, stroke = 0.35) +
  scale_shape_manual(name = NULL, values = c(Receiver = 21)) +
  facet_wrap(~design, ncol = 2) +
  coord_equal(expand = FALSE) +
  labs(x = "Easting (m)", y = "Northing (m)") +
  detfield::annotation_scale_bar(location = "br") +
  detfield::theme_positionr() +
  theme(legend.position = "inside", legend.position.inside = c(0.75, 0.25))
coverage_map

# ---- Coverage while two receivers are out ----
# Count the coverage of the planned array with S006 and S016 pulled
pulled <- c("S006", "S016")
coverage_summary(fields$planned,
                 listening = setdiff(planned_receivers$station_id, pulled))

# ---- Score the designs ----
# Evaluate each design on the same 20 animals, 6 hours each
evaluations <- lapply(designs, function(layout) {
  evaluate_array(layout, lake, range_model, n_animals = 20,
                 seed = 1, control = list(hours = 6, tag_delay = c(120, 180)))
})

# ---- Rate per animal ----
# Tabulate the detection rate, one row per design and animal
rates <- do.call(rbind, lapply(names(evaluations), function(design) {
  sent <- evaluations[[design]]$detections$transmissions
  rate <- tapply(sent$n_detections > 0, sent$path_id, mean)
  data.frame(design = panel_names[[design]], animal_id = names(rate),
             rate = as.vector(rate))
}))
head(rates, 3)

# ---- Compare the designs ----
# Tabulate coverage, detection and positioning, one row per design
scores <- do.call(rbind, lapply(names(evaluations), function(design) {
  overall <- evaluations[[design]]$score$overall
  design_rates <- rates$rate[rates$design == panel_names[[design]]]
  data.frame(design = design, receivers = nrow(designs[[design]]),
             covered_p50 = round(evaluations[[design]]$coverage$share[1], 3),
             mean_rate = round(mean(design_rates), 3),
             never_heard = sum(design_rates == 0),
             scored_share = round(overall$scored_share, 3),
             median_error_m = round(overall$median_error_mean_m))
}))
scores

# ---- Plot rate by design ----
# Compute the mean rate, and the animals never heard, per design
rate_means <- aggregate(rate ~ design, data = rates, FUN = mean)
never <- aggregate(rate ~ design, data = rates,
                   FUN = function(rate) sum(rate == 0))
never$label <- sprintf("%d at 0%%", never$rate)

# Draw each animal across the designs, its points joined, the means, and
# the count of animals never heard, which overlap at 0
ggplot(rates, aes(design, rate)) +
  geom_line(aes(group = animal_id), colour = "grey70", linewidth = 0.4) +
  geom_point(colour = "#21908C", alpha = 0.8, size = 1.8) +
  geom_errorbar(data = rate_means, aes(ymin = rate, ymax = rate), width = 0.4,
                colour = "grey20", linewidth = 0.8) +
  geom_text(data = never, aes(y = -0.07, label = label), size = 3.2,
            colour = "grey20") +
  scale_y_continuous(labels = function(share) sprintf("%g%%", 100 * share),
                     breaks = seq(0, 1, 0.25), limits = c(-0.1, 1)) +
  labs(x = "Candidate array", y = "Share of transmissions heard") +
  detfield::theme_positionr()

# ---- Planned under both models ----
# Score the grid of Get started and the planned array under the stand-in
reference_settings <- list(hours = 6, tag_delay = c(120, 180))
grid_stand_in <- evaluate_array(
  place_receivers(deployable, n = 20, method = "grid", seed = 1), lake,
  stand_in, n_animals = 20, seed = 1, control = reference_settings
)
planned_stand_in <- evaluate_array(designs$planned, lake, stand_in, n_animals = 20,
                                   seed = 1, control = reference_settings)

# Compare the hit rates, with the planned array under the fitted model
round(c(grid_stand_in = grid_stand_in$score$overall$hit_rate,
        planned_stand_in = planned_stand_in$score$overall$hit_rate,
        planned_fitted = evaluations$planned$score$overall$hit_rate), 3)

# ---- Station file against grid ----
# Compare the station file and the grid, detections and each animal-hour
same_detections <- identical(planned_stand_in$detections$detections,
                             grid_stand_in$detections$detections)
by_hour <- merge(planned_stand_in$score$bins, grid_stand_in$score$bins,
                 by = c("animal_id", "bin_start"), suffixes = c("_file", "_grid"))
differs <- by_hour[by_hour$in_region_file != by_hour$in_region_grid,
                   c("animal_id", "bin_start", "in_region_file", "in_region_grid")]
same_detections
differs

# ---- Position with half the range ----
# Halve d50 and d05 of the range model
half_range_model <- perturb_de_model(range_model, d50_scale = 0.5,
                                     d05_scale = 0.5)

# Detect under the true curve; position with the halved one
half_range_eval <- evaluate_array(
  designs$planned, lake, range_model, n_animals = 20, seed = 1,
  control = list(hours = 6, tag_delay = c(120, 180),
                 position_model = half_range_model)
)

# ---- Compare the ranges ----
# Compare the median error and the hit rate, true range against half
score_columns <- c("median_error_mean_m", "hit_rate")
rbind(true_range = evaluations$planned$score$overall[score_columns],
      half_range = half_range_eval$score$overall[score_columns])

# ---- Write the chosen stations ----
# Keep station, position and depth of the chosen design, depth to 0.1 m
chosen <- designs$planned
field_sheet <- chosen[c("station_id", "deploy_lat", "deploy_lon", "deploy_z")]
field_sheet$deploy_z <- round(field_sheet$deploy_z, 1)

# Write the sheet and read it back
sheet_file <- file.path(tempdir(), "stations_for_deployment.csv")
write.csv(field_sheet, sheet_file, row.names = FALSE)
head(read.csv(sheet_file), 3)

# Count the stations under 2 m deep
sum(chosen$deploy_z < 2)
