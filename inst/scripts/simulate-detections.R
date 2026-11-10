# Code of the vignette "Simulate and score an array step by step".
# Extracted from vignettes/simulate-detections.Rmd by tools/purl_vignettes.R, a repository file
# the installed package does not carry. The vignette is the source; this file follows it.

# ---- Set up the array ----
# Load the package, and ggplot2 for the maps
library(arraydesign)
library(ggplot2)

# Read the depth raster of Stoney Lake
lake <- terra::rast(system.file("extdata", "stoney_depth.tif",
                                package = "arraydesign"))

# Lay a grid of 20 receivers over water at least 1 m deep
deployable <- terra::ifel(abs(lake) < 1, NA, lake)
receivers <- place_receivers(deployable, n = 20, method = "grid", seed = 1)
nrow(receivers)

# ---- Detection field ----
# Set the assumed d50 and d05 for receivers 1 m and 35 m deep
range_model <- detfield::de_logistic(d50 = c(250, 750), d05 = c(500, 1500),
                                     depth = c(1, 35))

# Build the field: distances through water, then the curve at each
# distance and depth
field <- detfield::de_field(detfield::water_distances(lake, receivers),
                            range_model)

# Print the column that names the receivers of the field
field$key

# ---- Simulate tracks ----
# Build a walk: gamma steps of mean 30 m a minute, turns of sd 40 degrees
walk <- crw_model(step_mean = 30, step_sd = 15, turn_sd = 40, time_step_s = 60)

# Walk 20 animals for 360 steps of 60 s from midnight UTC
start <- as.POSIXct("2025-01-01 00:00:00", tz = "UTC")
tracks <- simulate_tracks(lake, walk, n_animals = 20, n_steps = 360,
                          control = track_control(start_time = start),
                          seed = 1)
tracks

# ---- Simulate transmissions ----
# Send a transmission every 120 to 180 s along each track
transmissions <- simulate_transmissions(tracks, tag_delay = c(120, 180),
                                        seed = 1)
head(transmissions)

# Count the transmissions in all, and per animal-hour
nrow(transmissions)
nrow(transmissions) / (20 * 6)

# ---- Simulate detections ----
# Give each receiver one trial per transmission
sim <- simulate_detections(transmissions, field, overdispersion = "none",
                           seed = 1)
sim

# ---- Detection rate per animal ----
# Take every transmission, with the number of receivers that heard it
truth <- sim$transmissions

# Compute the share of each animal's transmissions heard by any receiver
rate <- tapply(truth$n_detections > 0, truth$path_id, mean)
round(rate, 3)

# Give the mean, standard error and median across animals, and the
# animals never heard
c(mean = mean(rate), se = sd(rate) / sqrt(length(rate)), median = median(rate))
sum(rate == 0)

# ---- Animal hours heard ----
# Find the hour of the run in which each transmission was sent, 0 to 5
hour_of_run <- floor(as.numeric(truth$datetime - start, units = "hours"))

# Find the animal-hours with at least one detection, as a share of all 120
heard_hours <- unique(data.frame(path_id = truth$path_id,
                                 hour = hour_of_run)[truth$n_detections > 0, ])
nrow(heard_hours) / (20 * 6)

# ---- Map detected transmissions ----
# Take the water cells, with depth as positive metres
water <- terra::as.data.frame(lake, xy = TRUE, na.rm = TRUE)
names(water)[3] <- "depth_m"
water$depth_m <- abs(water$depth_m)

# Convert the tracks to a plain table
paths <- as.data.frame(tracks)[c("path_id", "x", "y")]

# Map the tracks, misses, heard transmissions, the starts (the first
# position of each track) and receivers over the basin. The colour scale keys
# the tracks and transmissions; the shape scale keys the receivers, a dark
# ring under a coloured one
ring <- detfield::positionr_ring_colours()
ggplot(water, aes(x, y)) +
  geom_raster(aes(fill = depth_m)) +
  detfield::scale_fill_depth() +
  geom_path(data = paths, aes(group = path_id, colour = "Track"),
            linewidth = 1) +
  geom_point(data = truth[truth$n_detections == 0, ], aes(colour = "Not heard"),
             shape = 4, size = 0.8, stroke = 0.4, alpha = 0.8) +
  geom_point(data = truth[truth$n_detections > 0, ], aes(colour = "Heard"),
             shape = 21, fill = "#2ECC40", size = 1.8, stroke = 0.3) +
  geom_point(data = paths[!duplicated(paths$path_id), ], aes(colour = "Start"),
             shape = 21, fill = "white", size = 2, stroke = 0.5) +
  scale_colour_manual(
    name = NULL, breaks = c("Track", "Start", "Heard", "Not heard"),
    values = c(Track = "grey85", Start = "black", Heard = "black",
               `Not heard` = "red3"),
    # The key draws the cross larger than the map does, so that it reads
    guide = guide_legend(position = "bottom",
                         override.aes = list(size = c(1, 2, 1.8, 1.6),
                                             stroke = c(0.5, 0.5, 0.3, 0.6),
                                             alpha = 1))
  ) +
  geom_point(data = receivers, aes(shape = "Receiver"), colour = ring[["edge"]],
             size = 4, stroke = 2.2) +
  geom_point(data = receivers, aes(shape = "Receiver"), colour = ring[["heard"]],
             size = 4, stroke = 1.1) +
  scale_shape_manual(name = NULL, values = c(Receiver = 21)) +
  coord_equal(xlim = c(730300, 737100), ylim = c(4937700, 4941940),
              expand = FALSE) +
  labs(x = "Easting (m)", y = "Northing (m)") +
  detfield::annotation_scale_bar(location = "br") +
  detfield::theme_positionr()

# ---- Deployments table ----
# Put every receiver in the water from a day before the run to a day after
deployments <- data.frame(receiver_serial = receivers$receiver_serial,
                          station_id = receivers$station_id,
                          x = receivers$x, y = receivers$y,
                          deploy_z = receivers$deploy_z,
                          deploy_datetime = start - 86400,
                          recover_datetime = start + 86400)

# Recover the receiver at S016 three hours into the run
recovered <- deployments$station_id == "S016"
deployments$recover_datetime[recovered] <- start + 3 * 3600
deployments[recovered, c("station_id", "deploy_z", "deploy_datetime",
                         "recover_datetime")]

# ---- Receiver recovered early ----
# Build the field from the deployments, then simulate with the same
# transmissions and seed
field_recovered <- detfield::de_field(
  detfield::water_distances(lake, deployments), range_model
)
sim_recovered <- simulate_detections(transmissions, field_recovered,
                                     overdispersion = "none", seed = 1)

# ---- Detections at s016 ----
# Count the detections at S016, listening for 6 hours and recovered at
# 3 hours
c(listening_6_hours = sum(sim$detections$station_id == "S016"),
  recovered_at_3_hours = sum(sim_recovered$detections$station_id == "S016"))

# ---- Hourly spread ----
# Add a station-hour sd of 1 on the logit scale to the range model
spread_model <- detfield::de_model(range_model, sigma = list(station_hour = 1))

# Build its detection field
field_spread <- detfield::de_field(detfield::water_distances(lake, receivers),
                                   spread_model)

# ---- Detections with spread ----
# Simulate with the same transmissions and seed, with the spread drawn
sim_spread <- simulate_detections(transmissions, field_spread, seed = 1)

# Compare the detection rate per animal, without and with the spread
rate_spread <- tapply(sim_spread$transmissions$n_detections > 0,
                      sim_spread$transmissions$path_id, mean)
round(rbind(static = summary(rate), with_spread = summary(rate_spread)), 3)

# ---- Listening table ----
# Set every receiver listening from a day before the run to a day after
listening <- data.frame(receiver_serial = receivers$receiver_serial,
                        station_id = receivers$station_id,
                        deploy_datetime = start - 86400,
                        recover_datetime = start + 86400)

# ---- Bin detections ----
# Count the detections per animal, receiver and hour
counts <- detfield::bin_detections(sim$detections, listening,
                                   tag_delay = c(120, 180))
head(counts)

# ---- Location surfaces ----
# Compute one surface per animal-hour with a detection
surfaces <- detfield::location_surface(counts, field, overdispersion = "none")

# ---- Score surfaces ----
# Score the error and hit rate of each surface against the truth
score <- score_surfaces(surfaces, truth)
score
head(score$bins[c("animal_id", "bin_start", "n_truth", "error_mean_m",
                  "in_region", "area_m2")])

# ---- Pick the mapped hour ----
# Count the receivers that heard each animal-hour
heard_by <- aggregate(y ~ animal_id + bin_start, data = counts,
                      FUN = function(detections) sum(detections > 0))
names(heard_by)[3] <- "receivers_heard"

# Find the animal-hours heard by the most receivers, with their scores
most_heard <- merge(
  heard_by[heard_by$receivers_heard == max(heard_by$receivers_heard), ],
  score$bins[c("animal_id", "bin_start", "n_truth", "error_mean_m",
               "in_region")]
)
most_heard

# Pick the one nearest the run's hit rate
mapped_hour <- most_heard[which.min(abs(most_heard$in_region -
                                          score$overall$hit_rate)), ]
mapped_hour

# ---- Surface map helper ----
# Map one animal-hour: water by depth, the surface in magma from 0 to its
# largest cell, the 95% region outlined, receivers that heard the animal as
# yellow rings sized by their detections, silent receivers as white rings,
# and the animal's true transmissions in the hour, keyed with the outline
surface_map <- function(surfaces, animal, hour, lake, receivers, detections,
                        truth, margin = 500) {
  # Take the cells of the one surface, with their coordinates
  hour_cells <- surfaces[surfaces$animal_id == animal &
                           surfaces$bin_start == hour, ]
  hour_cells <- cbind(hour_cells, terra::xyFromCell(lake, hour_cells$cell))

  # Outline its 95% region
  region <- detfield::credible_region(surfaces, level = 0.95)
  region <- region[region$animal_id == animal & region$bin_start == hour, ]
  outline <- terra::rast(lake)
  outline[region$cell] <- 1
  outline <- terra::geom(terra::as.polygons(outline), df = TRUE)

  # Take the true transmissions and detections of the animal in the hour
  in_hour <- function(times) times >= hour & times < hour + 3600
  hour_truth <- truth[truth$animal_id == animal & in_hour(truth$datetime), ]
  hour_detections <- detections[detections$animal_id == animal &
                                  in_hour(detections$datetime), ]

  # Count the detections per receiver in the hour
  rings <- data.frame(receiver_serial = receivers$receiver_serial,
                      x = receivers$x, y = receivers$y)
  rings$n <- as.vector(table(factor(hour_detections$receiver_serial,
                                    levels = rings$receiver_serial)))
  heard <- rings[rings$n > 0, ]

  # Zoom to the surface, the receivers that heard and the truth
  in_view <- rbind(hour_cells[c("x", "y")], heard[c("x", "y")],
                   hour_truth[c("x", "y")])
  xlim <- range(in_view$x) + c(-margin, margin)
  ylim <- range(in_view$y) + c(-margin, margin)

  # Take the water in view, coloured by depth with the ramp of
  # detfield::scale_fill_depth(). The fill scale holds the probability, so
  # the depth key is an alpha scale held at 1 whose boxes take the ramp
  depth_ramp <- scales::seq_gradient_pal("cornflowerblue", "blue4")
  water <- terra::as.data.frame(abs(lake), xy = TRUE, na.rm = TRUE)
  names(water)[3] <- "depth_m"
  deepest <- max(water$depth_m)
  water$colour <- depth_ramp(water$depth_m / deepest)
  depth_breaks <- pretty(c(0, deepest), n = 3)
  depth_breaks <- depth_breaks[depth_breaks <= deepest]
  water <- water[water$x > xlim[1] - 100 & water$x < xlim[2] + 100 &
                   water$y > ylim[1] - 100 & water$y < ylim[2] + 100, ]

  # Draw the layers
  top_prob <- max(hour_cells$prob)
  cell_size <- terra::res(lake)[1]
  ggplot(hour_cells, aes(x, y)) +
    geom_tile(data = water, aes(fill = I(.data$colour), alpha = .data$depth_m),
              width = cell_size, height = cell_size) +
    scale_alpha_continuous(
      name = "Depth (m)", range = c(1, 1), limits = c(0, deepest),
      breaks = depth_breaks,
      guide = guide_legend(position = "bottom", nrow = 1, override.aes = list(
        fill = depth_ramp(depth_breaks / deepest)
      ))
    ) +
    geom_tile(aes(fill = .data$prob), width = cell_size, height = cell_size) +
    detfield::scale_fill_probability(
      name = "Probability\nper cell", limits = c(0, top_prob),
      breaks = c(0, top_prob), labels = signif(c(0, top_prob), 2)
    ) +
    geom_path(data = outline, aes(group = interaction(.data$part, .data$hole),
                                  colour = "95% region"), linewidth = 0.9) +
    detfield::layer_receivers(rings, count = "n", size = 3, max_size = 9) +
    geom_point(data = hour_truth[hour_truth$n_detections == 0, ],
               aes(colour = "Not heard"), shape = 4, size = 1.6, stroke = 0.7) +
    geom_point(data = hour_truth[hour_truth$n_detections > 0, ],
               aes(colour = "Heard"), shape = 21, fill = "#2ECC40", size = 2,
               stroke = 0.3) +
    scale_colour_manual(
      name = NULL, breaks = c("95% region", "Heard", "Not heard"),
      values = c(`95% region` = "#00A087", Heard = "black", `Not heard` = "red3")
    ) +
    coord_equal(xlim = xlim, ylim = ylim, expand = FALSE) +
    labs(title = sprintf("Animal %s, %s to %s UTC", animal,
                         format(hour, "%H:%M", tz = "UTC"),
                         format(hour + 3600, "%H:%M", tz = "UTC")),
         x = "Easting (m)", y = "Northing (m)") +
    detfield::annotation_scale_bar(location = "br") +
    detfield::theme_positionr()
}

# ---- Pick a poor hour ----
# Take the most-heard animal-hour with the lowest share inside its region
poor_hour <- most_heard[which.min(most_heard$in_region), ]
poor_hour

# ---- Map a poor hour ----
# Map the poor hour from the surfaces of this section
surface_map(surfaces, poor_hour$animal_id, poor_hour$bin_start, lake,
            receivers, sim$detections, truth)
