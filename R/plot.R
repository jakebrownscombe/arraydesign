#' Map simulated tracks
#'
#' Draws the water of the raster the tracks were simulated on by depth, on
#' [detfield::scale_fill_depth()], with land left blank, one grey45 line per
#' path and a scale bar. Paths have no start or end marker.
#'
#' The plot's data are the positions drawn, one row per position with the
#' columns of `x`. To draw several simulations on one lake side by side, bind
#' their positions with a column that names each, replace the data with
#' [ggplot2::%+%], and facet by that column.
#'
#' @param x A `sim_tracks` object from [simulate_tracks()].
#' @param paths Path ids to draw; `NULL` draws every path.
#' @param col Line colours, recycled over the paths drawn; `NULL` draws every
#'   path in grey45, the track colour of the positionR maps.
#' @param ... Not used. Add ggplot2 layers, scales or [ggplot2::theme()]
#'   settings to the plot this method returns instead.
#'
#' @return A ggplot object, drawn when printed.
#' @seealso [simulate_tracks()]
#' @export
#' @examples
#' stoney <- terra::rast(system.file("extdata", "stoney_depth.tif",
#'                                   package = "arraydesign"))
#' # Three animals for 6 hours.
#' walk <- crw_model(30, 15, 40, time_step_s = 60)
#' tr <- simulate_tracks(stoney, walk, n_animals = 3, n_steps = 360, seed = 1)
#' plot(tr)
#' # One animal, zoomed to its extent with a 500 m margin.
#' one <- tr[tr$path_id == 2, ]
#' plot(tr, paths = 2) +
#'   ggplot2::coord_equal(xlim = range(one$x) + c(-500, 500),
#'                        ylim = range(one$y) + c(-500, 500), expand = FALSE)
plot.sim_tracks <- function(x, paths = NULL, col = NULL, ...) {
  unused_dots(...)
  grid <- terra::unwrap(attr(x, "grid"))
  ids <- if (is.null(paths)) unique(x$path_id) else paths
  missing <- setdiff(ids, x$path_id)
  if (length(missing)) {
    stop_classed("input", "`paths` names paths not in `x`: ", paste(missing, collapse = ", "), ".")
  }
  d <- unclass_df(x)
  d <- d[d$path_id %in% ids, , drop = FALSE]
  d$path <- factor(d$path_id, levels = ids)
  if (is.null(col)) col <- track_grey
  ggplot2::ggplot(d, ggplot2::aes(x = .data$x, y = .data$y, colour = .data$path)) +
    depth_layers(grid) +
    ggplot2::geom_path(ggplot2::aes(group = .data$path), linewidth = 0.4) +
    ggplot2::scale_colour_manual(values = rep_len(col, length(ids)), guide = "none") +
    map_frame()
}

#' Colours for simulated tracks
#'
#' A discrete ggplot2 colour scale for tracks drawn over water by depth. Its
#' six colours are the Okabe-Ito palette of [grDevices::palette.colors()]
#' less its black, grey and pale yellow, in this order: orange `#E69F00`,
#' bluish green `#009E73`, reddish purple `#CC79A7`, blue `#0072B2`,
#' vermillion `#D55E00` and sky blue `#56B4E9`. With more than six animals
#' or groups the colours repeat. Draw each track over a wider white line, so
#' that the blues read on deep water.
#'
#' Give `colours` to tie a colour to each level by name. The suite uses
#' bluish green for moving and orange for resting, with resting positions
#' drawn as squares so that the state is not carried by colour alone.
#'
#' @param ... Passed to [ggplot2::discrete_scale()], such as `name`,
#'   `labels` or `guide`.
#' @param colours `NULL` to give the levels the six colours in order, or a
#'   character vector named by level whose values name colours of the set:
#'   `"orange"`, `"bluish_green"`, `"reddish_purple"`, `"blue"`,
#'   `"vermillion"` or `"sky_blue"`.
#'
#' @return A ggplot2 scale for the `colour` aesthetic.
#' @seealso [plot.sim_tracks()]
#' @export
#' @examples
#' stoney <- terra::rast(system.file("extdata", "stoney_depth.tif",
#'                                   package = "arraydesign"))
#' # A walk that rests: keep moving 0.95, rest 0.05; keep resting 0.90
#' rests <- matrix(c(0.95, 0.05, 0.10, 0.90), nrow = 2, byrow = TRUE)
#' walk <- crw_model(30, 15, 40, time_step_s = 60, transition = rests)
#' # Three animals for 4 hours from one point in the middle of the lake
#' middle <- data.frame(x = 727850, y = 4937250)
#' tracks <- as.data.frame(simulate_tracks(
#'   stoney, walk, n_animals = 3, n_steps = 240,
#'   control = track_control(start = middle), seed = 1
#' ))
#' # One colour per animal, each track over a white casing
#' ggplot2::ggplot(tracks, ggplot2::aes(x, y, group = path_id)) +
#'   ggplot2::geom_path(colour = "white", linewidth = 1.2) +
#'   ggplot2::geom_path(ggplot2::aes(colour = factor(path_id))) +
#'   scale_colour_tracks(name = "Animal") +
#'   ggplot2::coord_equal() +
#'   detfield::theme_positionr()
#' # One colour per state, named; resting positions as squares
#' resting <- tracks[tracks$state == 2, ]
#' ggplot2::ggplot(tracks, ggplot2::aes(x, y)) +
#'   ggplot2::geom_path(ggplot2::aes(group = path_id, colour = "Moving")) +
#'   ggplot2::geom_point(data = resting, ggplot2::aes(colour = "Resting"),
#'                       shape = 15) +
#'   scale_colour_tracks(name = "State",
#'                       colours = c(Moving = "bluish_green",
#'                                   Resting = "orange")) +
#'   ggplot2::coord_equal() +
#'   detfield::theme_positionr()
scale_colour_tracks <- function(..., colours = NULL) {
  if (is.null(colours)) {
    return(ggplot2::discrete_scale(
      "colour", palette = function(n) rep_len(unname(track_colours), n), ...
    ))
  }
  stop_unless(named_track_colours(colours),
              paste0("`colours` must be named by level, with values from: ",
                     paste(names(track_colours), collapse = ", "), "."))
  ggplot2::scale_colour_manual(
    values = stats::setNames(unname(track_colours[colours]), names(colours)), ...
  )
}

#' Map a receiver layout
#'
#' Draws the water of the depth raster by depth on
#' [detfield::scale_fill_depth()], with land left blank, each receiver as a
#' yellow ring over a dark ring labelled with its station name, and a scale
#' bar.
#'
#' @param x A `receiver_layout` from [place_receivers()].
#' @param grid The depth raster to draw under the receivers, or `NULL` for the
#'   receivers alone.
#' @param labels The column of `x` to label each receiver with:
#'   `"station_id"` (the default, also used for `TRUE`), `"receiver_serial"`,
#'   or `FALSE` for no labels.
#' @param ... Not used. Add ggplot2 layers, scales or [ggplot2::theme()]
#'   settings to the plot this method returns instead.
#'
#' @return A ggplot object, drawn when printed.
#' @seealso [place_receivers()]
#' @export
#' @examples
#' stoney <- terra::rast(system.file("extdata", "stoney_depth.tif",
#'                                   package = "arraydesign"))
#' planned <- read.csv(system.file("extdata", "planned_stations.csv",
#'                                 package = "arraydesign"))
#' rx <- place_receivers(stoney, method = "manual", coords = planned)
#' plot(rx, stoney)
#' # Labelled by receiver serial instead of station.
#' plot(rx, stoney, labels = "receiver_serial")
plot.receiver_layout <- function(x, grid = NULL, labels = "station_id", ...) {
  unused_dots(...)
  if (isTRUE(labels)) labels <- "station_id"
  stop_unless(isFALSE(labels) || (is.character(labels) && length(labels) == 1L &&
                                    labels %in% names(x)),
              "`labels` must be FALSE or the name of a column of `x`, such as \"station_id\".")
  d <- unclass_df(x)
  if (!isFALSE(labels)) d$label <- as.character(d[[labels]])
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$x, y = .data$y))
  if (!is.null(grid)) p <- p + depth_layers(water_grid(grid)$raster)
  p <- p + receiver_rings()
  if (!isFALSE(labels)) {
    p <- p + ggplot2::geom_label(
      ggplot2::aes(label = .data$label), size = 2.4, vjust = -0.7,
      label.padding = ggplot2::unit(0.1, "lines"), fill = grDevices::adjustcolor("white", 0.8),
      linewidth = 0
    )
  }
  # With no raster, the bottom margin holds the scale bar below the lowest
  # receivers, so that the bar does not cover a receiver or its label.
  if (is.null(grid)) {
    return(p + ggplot2::coord_equal(clip = "off") +
             ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.3, 0.1))) +
             ggplot2::labs(x = "Easting (m)", y = "Northing (m)") +
             detfield::annotation_scale_bar(location = "br") +
             detfield::theme_positionr())
  }
  p + map_frame()
}

# TRUE when `colours` names levels and its values name colours of the set.
named_track_colours <- function(colours) {
  levels <- names(colours)
  is.character(colours) && length(colours) > 0L &&
    all(colours %in% names(track_colours)) &&
    length(levels) == length(colours) && all(nzchar(levels), na.rm = TRUE)
}

# ---- Map layers ----
# The colours of the positionR maps: water by
# depth on detfield's cornflowerblue to blue4 ramp, land left blank on the
# white panel, tracks grey45, and receivers as a yellow ring over a dark ring
# one step wider, so a receiver reads on light and dark water alike. A
# receiver that heard none of the animals is a white ring over the dark ring.
track_grey <- "grey45"
# Tracks that are the subject of a figure are drawn in colour, one per animal
# or group: the Okabe-Ito palette of grDevices::palette.colors() less black
# and grey, which read as the grey track, and the pale yellow, which washes
# out on the light end of the depth ramp. The first two are the two-state
# colours, moving bluish green and resting orange, the pair furthest apart.
track_colours <- c(orange = "#E69F00", bluish_green = "#009E73",
                   reddish_purple = "#CC79A7", blue = "#0072B2",
                   vermillion = "#D55E00", sky_blue = "#56B4E9")
# The ring colours are detfield's, so every positionR map draws the same.
ring_yellow <- detfield::positionr_ring_colours()[["heard"]]
ring_white <- detfield::positionr_ring_colours()[["silent"]]
ring_edge <- detfield::positionr_ring_colours()[["edge"]]

# Water cells by depth. Depth is drawn as its absolute value, so a raster
# that stores depth as negative numbers (the bundled Stoney Lake raster runs
# -32.5 to 0 m) still draws its deepest cell blue4. Land cells are not drawn.
depth_layers <- function(grid, legend = TRUE) {
  d <- terra::as.data.frame(grid, xy = TRUE, na.rm = TRUE)
  names(d)[3] <- "depth"
  d$depth <- abs(d$depth)
  list(
    ggplot2::geom_raster(data = d, ggplot2::aes(x = .data$x, y = .data$y, fill = .data$depth),
                         inherit.aes = FALSE),
    detfield::scale_fill_depth(guide = if (legend) "colourbar" else "none")
  )
}

# Receivers as rings: the dark ring first, the coloured one over it. `data`
# is passed to both layers, as a data frame or a function of the plot's data.
receiver_rings <- function(data = NULL, size = 3, colour = ring_yellow) {
  list(
    ggplot2::geom_point(data = data, shape = 1, colour = ring_edge, size = size, stroke = 1.8),
    ggplot2::geom_point(data = data, shape = 1, colour = colour, size = size, stroke = 0.7)
  )
}

# Equal coordinates with no margin round the raster, axis titles in metres,
# a scale bar, and the suite theme. Drawing is not clipped to the panel, so
# a receiver on the edge of the raster, and its label, are drawn whole.
map_frame <- function() {
  list(
    ggplot2::coord_equal(expand = FALSE, clip = "off"),
    ggplot2::labs(x = "Easting (m)", y = "Northing (m)"),
    detfield::annotation_scale_bar(location = "br"),
    detfield::theme_positionr()
  )
}

unused_dots <- function(...) {
  if (...length()) {
    warning("Arguments in `...` are not used: add ggplot2 layers, scales or theme() ",
            "settings to the plot this method returns.", call. = FALSE)
  }
  invisible()
}

#' @export
print.sim_tracks <- function(x, ...) {
  esc <- attr(x, "escapes")
  steps <- sum(esc[, -1])
  n_animals <- length(unique(x$path_id))
  cat(sprintf("<sim_tracks> %d animal%s, %d positions, time step %g s\n",
              n_animals, if (n_animals == 1) "" else "s", nrow(x), attr(x, "model")$time_step))
  cat(sprintf("steps: %d; as drawn %d, reflected at the shore %d, held %d, still %d\n",
              steps, sum(esc$clear), sum(esc$reflected), sum(esc$held), sum(esc$still)))
  sel <- attr(x, "selection")
  if (!is.null(sel)) {
    cat_wrapped(sprintf(paste("selection: %d of %d candidates accepted (%.3f); longer than",
                              "the reach %d, above the bound %d, exhausted %d"),
                        sum(sel$accepted), sum(sel$proposals),
                        sum(sel$accepted) / max(sum(sel$proposals), 1), sum(sel$over_reach),
                        sum(sel$above_bound), sum(sel$exhausted)))
  }
  print(utils::head(as.data.frame(unclass_df(x))), ...)
  invisible(x)
}

#' @export
print.sim_detections <- function(x, ...) {
  tx <- x$transmissions
  cat(sprintf("<sim_detections> %d detections of %d transmissions; %d heard (%.1f%%)\n",
              nrow(x$detections), nrow(tx), sum(tx$n_detections > 0),
              100 * mean(tx$n_detections > 0)))
  print(utils::head(x$detections), ...)
  invisible(x)
}

unclass_df <- function(x) {
  attributes(x) <- attributes(x)[c("names", "row.names")]
  class(x) <- "data.frame"
  x
}
