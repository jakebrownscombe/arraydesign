# Code of the vignette "Design and test an array for Stoney Lake".
# Extracted from vignettes/arraydesign.Rmd by tools/purl_vignettes.R, a repository file
# the installed package does not carry. The vignette is the source; this file follows it.

# ---- Read lake ----
# Load the package, and ggplot2 for the maps
library(arraydesign)
library(ggplot2)

# Read the depth raster included with the package
lake <- terra::rast(system.file("extdata", "stoney_depth.tif",
                                package = "arraydesign"))
lake

# ---- Place on raster ----
# Set water shallower than 1 m to NA, so that no receiver goes there
deployable <- terra::ifel(abs(lake) < 1, NA, lake)

# Lay a grid of 20 receivers over the deployable water
receivers <- place_receivers(deployable, n = 20, method = "grid", seed = 1)
head(receivers)

# Print the depth range under the receivers, m
range(receivers$deploy_z)

# ---- Place receivers ----
# Map the layout over the lake
plot(receivers, lake)

# ---- Range model ----
# Set the assumed range model: d50 and d05 for receivers 1 m and 35 m deep
range_model <- detfield::de_logistic(d50 = c(250, 750), d05 = c(500, 1500),
                                     depth = c(1, 35))
range_model

# ---- Detection field build ----
# Measure the distance through water from each receiver to each water cell
distances <- detfield::water_distances(lake, receivers)

# Read the curve at each distance: one probability per receiver and cell
field <- detfield::de_field(distances, range_model)

# ---- Coverage summary ----
# Count the cells one receiver alone hears at p >= 0.5 and p >= 0.8
coverage_summary(field)

# ---- Array probability ----
# Compute the chance that one or more receivers hear, per water cell
array_p <- data.frame(p = 1 - apply(1 - field$p, 1, prod),
                      terra::xyFromCell(lake, field$cells))

# Count the cells at p >= 0.5 under the array, and all water cells
c(cells_at_0.5 = sum(array_p$p >= 0.5), water_cells = nrow(array_p))

# ---- Detection field ----
# Take the ring colours of the positionR maps
ring <- detfield::positionr_ring_colours()

# Map the array's detection probability, with the receivers on top as a dark
# ring under a coloured one; the shape scale keys them
ggplot(array_p, aes(x, y)) +
  geom_raster(aes(fill = p)) +
  detfield::scale_fill_probability(name = "Detection\nprobability") +
  geom_point(data = receivers, aes(shape = "Receiver"), colour = ring[["edge"]],
             size = 1.3, stroke = 0.9) +
  geom_point(data = receivers, aes(shape = "Receiver"), colour = ring[["heard"]],
             size = 1.3, stroke = 0.45) +
  scale_shape_manual(name = NULL, values = c(Receiver = 21)) +
  coord_equal(expand = FALSE) +
  labs(x = "Easting (m)", y = "Northing (m)") +
  detfield::annotation_scale_bar(location = "br") +
  detfield::theme_positionr()

# ---- Movement walk ----
# Build a correlated random walk: gamma steps, wrapped normal turns,
# 60 s a step
walk <- crw_model(step_mean = 30, step_sd = 15, turn_sd = 40, time_step_s = 60)
walk

# ---- Evaluate ----
# Simulate 20 animals for 6 hours, detect, position and score them
evaluation <- evaluate_array(receivers, lake, range_model, n_animals = 20,
                             seed = 1,
                             control = list(hours = 6, tag_delay = c(120, 180),
                                            movement = walk))
evaluation

# ---- Detection rate ----
# Take the transmissions, with the number of receivers that heard each
sent <- evaluation$detections$transmissions

# Compute the share of each animal's transmissions heard by any receiver
rate_per_animal <- tapply(sent$n_detections > 0, sent$path_id, mean)
c(mean = mean(rate_per_animal), median = median(rate_per_animal),
  never_heard = sum(rate_per_animal == 0))

# ---- Receiver detections ----
# Count the detections of each receiver; a receiver that heard none counts 0
heard_at <- evaluation$detections$detections$receiver_serial
receiver_counts <- data.frame(
  station_id = receivers$station_id, x = receivers$x, y = receivers$y,
  detections = as.vector(table(factor(heard_at,
                                      levels = receivers$receiver_serial)))
)

# List the stations that heard none of the animals
receiver_counts$station_id[receiver_counts$detections == 0]

# ---- Simulate animals ----
# Take the water cells, with depth as positive metres
water <- terra::as.data.frame(abs(lake), xy = TRUE, na.rm = TRUE)
names(water)[3] <- "depth_m"

# Take the first position of each track, its start
tracks_table <- as.data.frame(evaluation$tracks)
starts <- tracks_table[!duplicated(tracks_table$path_id), ]

# Map the water, the tracks, the misses, the heard transmissions, the starts
# and the receivers, in that order; the colour scale keys the marks
ggplot(tracks_table, aes(x, y)) +
  geom_raster(data = water, aes(fill = depth_m)) +
  detfield::scale_fill_depth() +
  geom_path(aes(group = path_id, colour = "Track"), linewidth = 1.1) +
  geom_point(data = sent[sent$n_detections == 0, ], aes(colour = "Not heard"),
             shape = 4, size = 0.35, stroke = 0.25, alpha = 0.7) +
  geom_point(data = sent[sent$n_detections > 0, ], aes(colour = "Heard"),
             shape = 21, fill = "#2ECC40", size = 1.7, stroke = 0.3) +
  geom_point(data = starts, aes(colour = "Start"), shape = 21, fill = "white",
             size = 1.8, stroke = 0.5) +
  scale_colour_manual(
    name = NULL, breaks = c("Track", "Start", "Heard", "Not heard"),
    values = c(Track = "grey85", Start = "black", Heard = "black",
               `Not heard` = "red3"),
    # The key draws the cross larger than the map does, so that it reads
    guide = guide_legend(position = "bottom",
                         override.aes = list(size = c(1, 1.8, 1.7, 1.6),
                                             stroke = c(0.5, 0.5, 0.3, 0.6),
                                             alpha = 1))
  ) +
  detfield::layer_receivers(receiver_counts, count = "detections",
                            size = 3.4) +
  coord_equal(expand = FALSE) +
  labs(x = "Easting (m)", y = "Northing (m)") +
  detfield::annotation_scale_bar(location = "br") +
  detfield::theme_positionr()

# ---- Evaluate array ----
# Map the positioning error of each scored animal-hour
plot(evaluation)
