# ---- Habitat selection: acceptance rates for the help page ----
# Computes the two acceptance rates in the details of ?habitat_selection: the
# share of candidate steps the walker keeps under a preference for bottom
# depth of beta = -0.3 and beta = 0.3 per m. The setting is that of section 5
# of vignette("movement-models"): the 100 m Stoney raster, crw_model(30, 15,
# 40, time_step_s = 60), 20 release cells drawn once on water with
# set.seed(20), 2,880 steps of 60 s (48 hours) per animal and seed = 1.
# Writes validation/outputs/04_selection_acceptance.csv.
# Run from the package root, with the package installed:
#   Rscript validation/04_selection_acceptance.R \
#     > validation/outputs/04_selection_acceptance.log 2>&1

# ---- Libraries ----
library(arraydesign)
out <- "validation/outputs"
cat("arraydesign", as.character(packageVersion("arraydesign")), ";", R.version.string, "\n")

# ---- Setting ----
lake <- terra::rast(system.file("extdata", "stoney_depth.tif", package = "arraydesign"))
water <- terra::as.data.frame(lake, xy = TRUE, na.rm = TRUE)
walk <- crw_model(step_mean = 30, step_sd = 15, turn_sd = 40, time_step_s = 60)
depth <- abs(lake)
names(depth) <- "depth_m"
set.seed(20)
releases <- water[sample(nrow(water), 20), c("x", "y")]
cat("water cells:", nrow(water), "; release cells:", nrow(unique(releases)), "\n")

# ---- Walk under each preference ----
rates <- do.call(rbind, lapply(c(-0.3, 0.3), function(b) {
  tracks <- simulate_tracks(lake, habitat_selection(walk, depth, beta = b), n_animals = 20,
                            n_steps = 2880, control = track_control(start = releases),
                            seed = 1)
  sel <- attr(tracks, "selection")
  data.frame(beta = b, n_animals = nrow(sel), steps_per_animal = 2880,
             proposals = sum(sel$proposals), accepted = sum(sel$accepted),
             acceptance = round(sum(sel$accepted) / sum(sel$proposals), 4),
             over_reach = sum(sel$over_reach), above_bound = sum(sel$above_bound),
             exhausted = sum(sel$exhausted))
}))
print(rates, row.names = FALSE)
write.csv(rates, file.path(out, "04_selection_acceptance.csv"), row.names = FALSE)
