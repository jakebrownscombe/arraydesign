# arraydesign 0.1.1

First release.

## Arrays

- `place_receivers()` puts receivers on a depth raster on a regular grid, at random water cells
  or at the coordinates you give. Grid and random methods return exactly `n` receivers, each at
  the centre of its cell. Coordinates in latitude and longitude are projected to the raster; a
  point just on land, within one cell diagonal of water, moves to the nearest water cell with a
  message. Any other point on land or off the raster is refused with an error that names it,
  or the first 10 and the count. Points off the raster are refused first; the points on land
  are named once every point is on the raster. A raster in longitude and latitude or in a
  projected CRS whose unit is not the metre is refused. The scale of the CRS is measured at
  the raster centre on its two axes and its diagonals. A raster is refused when any of the
  four is more than 1 percent from 1, as in Web Mercator away from the equator. The message
  names the scale and the UTM zone of the centre. A raster whose scale cannot be measured,
  such as one in a local engineering CRS, is refused as well.
- `coverage_summary()` counts, at each threshold, the water cells that at least one receiver
  hears on its own, optionally for the receivers in the water on a given day (`listening`). The
  print of `evaluate_array()` also gives the cells the array hears, any receiver hearing.

## Simulation

- `movement_model()` holds step and turn samplers, a time step and an optional switch between
  moving and resting. `crw_model()` builds a correlated random walk with gamma steps and
  wrapped-normal turns.
- `simulate_tracks()` walks animals through water only: a step that meets the shore is
  reflected off it, so no track crosses land. `track_control()` sets start positions, the
  start time and a goal that pulls the animals towards it, with the step at which each animal
  arrived kept in `attr(tracks, "arrived")`.
- `habitat_selection()` gives a movement model a preference for a habitat layer, as a
  coefficient per unit of the layer or as any log-preference function. The walker keeps each
  step with a probability set by the preference at its end, so the coefficient is known
  exactly at the scale of one step; tracks record each step's heading and the candidates drawn
  per animal. A model without a preference skips the sampler.
- `available_steps()` draws control steps from the walk for each step an animal took, for a
  conditional logit that checks the preference on the true tracks. Its `case` and `stratum`
  columns carry the names of **detfield**'s used and available points, so one model call fits
  both: `clogit(case ~ depth_m + strata(stratum))`.
- `simulate_transmissions()` sends transmissions at a random start-to-start delay
  (`tag_delay`), at positions interpolated along the track.
- `simulate_detections()` gives each receiver one chance to hear each transmission, with the
  probability of a `detfield` detection field. It draws an effect per receiver and hour when
  the range model has one (`overdispersion = "station_hour"`), and makes trials only while a
  receiver is deployed when the field carries deployment times.

## Scoring

- `evaluate_array()` scores a layout in one call: it builds the detection field, simulates
  animals, transmissions and detections, positions every detected animal-hour with
  `detfield::location_surface()`, and scores the positions against the true tracks. Its
  `seed` is passed to each of the three simulation steps in turn, so `seed = 1` gives the same
  animals and detections as the arraydesignR app's run at seed 1. The count of receivers a
  covered cell needs is `control$k`.
- `evaluate_array(control = list(likelihood = ...))` positions by binomial (the default) or
  presence-only grid localization; `control$position_model` positions with a range model other
  than the one that simulated the detections.
- `score_surfaces()` measures the error of the mode and mean of each location surface, the
  share of the truth inside its credible region (`hit_rate`, with a clustered standard error and
  95% interval) and the region's area.
- `perturb_de_model()` scales the 50% and 5% distances of a range model, to see what a wrong
  range model does to the score.

## Maps

- `plot()` draws tracks, a receiver layout and an array evaluation on the positionR colours:
  water by depth, land blank, tracks in grey, receivers as yellow rings and a scale bar. The
  error map of an evaluation colours each scored animal-hour by error class, by default 0 to
  200, 200 to 400, 400 to 800 and over 800 m, and draws the receivers that heard none of the
  animals as white rings. Errors below a first limit above 0 get a class of their own, so
  every animal-hour is drawn in the class that holds its error.
- `scale_colour_tracks()` colours tracks that are the subject of a map, one colour per animal
  or group, from six colours of the Okabe-Ito palette without its pale yellow: orange, bluish
  green, reddish purple, blue, vermillion and sky blue, in that order. `colours` ties a colour
  to each level by name, such as bluish green for moving and orange for resting.

## Interface

- The functions share argument names: the raster is `grid` in every function, including
  `place_receivers()`; the count of simulated animals is `n_animals` in `simulate_tracks()` and
  `evaluate_array()`; and the time one step takes is `time_step_s`, in seconds, in
  `movement_model()` and `crw_model()`.
- Every error that arraydesign raises itself has the classes `arraydesign_error_<kind>` and
  `arraydesign_error`, so it can be caught by class. The exception is an error about `grid` in
  `plot()`, which reaches the caller as a plain error because S4 method dispatch drops the
  class. The kinds are `input`, `grid`, `placement`, `track` and `field`.
- Coordinates that do not read as numbers are refused with the column, the count and the first
  such value and its row.

## Example data and documentation

- Example files on Stoney Lake, Ontario, in `inst/extdata`: the 100 m depth raster, a 48-hour
  range test and 20 planned stations. Every help example and article runs on them.
- Get started, `vignette("arraydesign")`, tests an array on Stoney Lake: place receivers, map
  coverage, set a movement model, simulate detections and score positions.
- Three articles: "Design an array for your own water body: compare layouts" (from
  your own raster, station lists and range model to a station sheet), "Specifying and testing
  movement models", which includes a preference for bottom depth checked on the true tracks, and
  "Simulate and score an array step by step". Each article is also a plain R script in
  `system.file("scripts", package = "arraydesign")`.
- Every page and example measures distance along **detfield**'s exact shortest path through
  water, which bends at the shore instead of following the grid. The planned stations are the
  grid written in latitude and longitude; the round trip through degrees moves each receiver by
  a few millimetres, so the two give hit rates that differ in the third decimal.
- The help of `movement_model()` and `crw_model()` cites no source for the correlated random
  walk itself.
