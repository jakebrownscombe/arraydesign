# ---- Plot methods ----
# Each method returns a ggplot object. The tests read its layers, mappings,
# scales and labels, and build it, rather than draw it on a device.

geoms_of <- function(p) unname(vapply(p$layers, function(l) class(l$geom)[1], character(1)))
scale_call <- function(p, aes) paste(deparse(p$scales$get_scales(aes)$call), collapse = "")

# The receiver ring colours of the positionR maps.
ring_colours <- function(p, at) {
  toupper(vapply(at, function(i) unique(ggplot2::layer_data(p, i)$colour), character(1)))
}

test_that("depth stored as negative numbers still draws the deepest cell blue4", {
  # The bundled Stoney raster stores depth as -32.5 to 0 m.
  neg <- -example_lake()
  rx <- place_receivers(example_lake(), n = 3)
  p <- plot(rx, neg)
  water <- ggplot2::layer_data(p, 1)
  depth <- terra::values(neg)[!is.na(terra::values(neg))]
  expect_equal(toupper(water$fill[which.min(depth)]), "#00008B")
  expect_equal(toupper(water$fill[which.max(depth)]), "#6495ED")
  pos <- ggplot2::layer_data(plot(rx, example_lake()), 1)
  expect_identical(water$fill, pos$fill)
})

test_that("the track map draws depth, grey paths with no markers, and a scale bar", {
  lake <- example_lake()
  tr <- simulate_tracks(lake, crw_model(30, 10, 30), n_animals = 3, n_steps = 50, seed = 1)
  p <- plot(tr)
  expect_s3_class(p, "ggplot")
  expect_true(withVisible(plot(tr))$visible)
  # Land is blank: one raster layer, of the water. No start or end markers.
  expect_identical(geoms_of(p), c("GeomRaster", "GeomPath", "ScaleBarRect", "ScaleBarText"))
  expect_match(scale_call(p, "fill"), "scale_fill_gradient")
  expect_identical(unique(ggplot2::layer_data(p, 2)$colour), "grey45")
  # Deeper water is darker: the deepest cell takes blue4, the shallowest
  # cornflowerblue.
  water <- ggplot2::layer_data(p, 1)
  depth <- terra::values(lake)[!is.na(terra::values(lake))]
  expect_equal(toupper(water$fill[which.max(depth)]), "#00008B")
  expect_equal(toupper(water$fill[which.min(depth)]), "#6495ED")
  labs <- ggplot2::get_labs(p)
  expect_identical(c(labs$x, labs$y, labs$fill), c("Easting (m)", "Northing (m)", "Depth (m)"))
  expect_equal(p$coordinates$ratio, 1)
  # 51 positions per path; the water layer holds the 2152 water cells only.
  expect_equal(nrow(ggplot2::layer_data(p, 2)), 153)
  expect_equal(nrow(ggplot2::layer_data(p, 1)), 2152)

  one <- plot(tr, paths = 2, col = "red")
  expect_equal(unique(ggplot2::layer_data(one, 2)$colour), "red")
  expect_equal(nrow(ggplot2::layer_data(one, 2)), 51)
  expect_error(plot(tr, paths = 9), "not in `x`")
  expect_warning(plot(tr, legend = FALSE), "not used")
})

test_that("the layout map draws receivers as yellow rings and labels them", {
  lake <- example_lake()
  rx <- place_receivers(lake, n = 8)
  p <- plot(rx, lake)
  expect_s3_class(p, "ggplot")
  expect_identical(geoms_of(p), c("GeomRaster", "GeomPoint", "GeomPoint", "GeomLabel",
                                  "ScaleBarRect", "ScaleBarText"))
  expect_identical(ring_colours(p, 2:3), c("#1A1A1A", "#FDE725"))
  expect_true(all(ggplot2::layer_data(p, 2)$shape == 1))
  expect_identical(ggplot2::layer_data(p, 4)$label, rx$station_id)
  expect_identical(ggplot2::layer_data(plot(rx, lake, labels = "receiver_serial"), 4)$label,
                   rx$receiver_serial)
  bare <- plot(rx, labels = FALSE)
  expect_identical(geoms_of(bare), c("GeomPoint", "GeomPoint", "ScaleBarRect", "ScaleBarText"))
  expect_identical(ggplot2::get_labs(bare)$x, "Easting (m)")
  expect_error(plot(rx, labels = "name"), "`labels` must be FALSE or the name of a column")
})

test_that("the error map draws depth, edged error points and receiver rings", {
  lake <- example_lake()
  ev <- suppressMessages(evaluate_array(place_receivers(lake, n = 8), lake,
                                        detfield::de_logistic(250, 600), n_animals = 2,
                                        seed = 1, control = list(hours = 2)))
  p <- plot(ev)
  expect_s3_class(p, "ggplot")
  expect_identical(geoms_of(p), c("GeomRaster", "GeomPoint", "GeomPoint", "GeomPoint",
                                  "GeomPoint", "GeomPoint", "GeomPoint", "ScaleBarRect",
                                  "ScaleBarText"))
  expect_identical(unique(ggplot2::layer_data(p, 2)$colour), "black")
  expect_equal(nrow(ggplot2::layer_data(p, 3)), nrow(ev$score$bins))
  # Receivers that heard none are white rings (layers 4 and 5), those that
  # heard the animals yellow rings (layers 6 and 7), each over a dark ring.
  heard <- table(factor(ev$detections$detections$receiver_serial,
                        levels = ev$layout$receiver_serial))
  expect_equal(nrow(ggplot2::layer_data(p, 5)), sum(heard == 0))
  expect_equal(nrow(ggplot2::layer_data(p, 7)), sum(heard > 0))
  expect_identical(p$data$n_detections[p$data$what == "receiver"], as.integer(heard))
  expect_identical(ring_colours(p, 6:7), c("#1A1A1A", "#FDE725"))
  if (any(heard == 0)) expect_identical(ring_colours(p, 4:5), c("#1A1A1A", "#FFFFFF"))
  expect_match(scale_call(p, "colour"), "scale_colour_viridis_d")
  expect_s3_class(p$scales$get_scales("colour")$guide, "GuideLegend")
  expect_identical(p$scales$get_scales("fill")$name, "Depth (m)")
  quiet <- plot(ev, legend = FALSE)
  expect_identical(quiet$scales$get_scales("colour")$guide, "none")
  expect_identical(quiet$scales$get_scales("fill")$guide, "none")
  expect_error(plot(ev, legend = NA), "`legend` must be TRUE or FALSE")
  withr::local_pdf(NULL)
  expect_s3_class(ggplot2::ggplotGrob(p), "gtable")
})

test_that("print methods summarise", {
  tr <- simulate_tracks(example_lake(), crw_model(30, 10, 30), n_animals = 2, n_steps = 50,
                        seed = 1)
  expect_output(print(tr), "<sim_tracks> 2 animals, 102 positions, time step 60 s")
  expect_output(print(tr), "steps: 100")
})

test_that("example_lake has the stated size, CRS and peninsula", {
  lake <- example_lake()
  expect_equal(dim(lake)[1:2], c(40, 60))
  expect_equal(terra::res(lake), c(30, 30))
  expect_equal(sum(!is.na(terra::values(lake))), 2152)
  expect_equal(arraydesign:::grid_epsg(lake), 32617L)
  # The peninsula: columns 30 and 31 are land from row 14 down, water above.
  v <- matrix(terra::values(lake), 40, 60, byrow = TRUE)
  expect_true(all(is.na(v[14:40, 30:31])))
  expect_false(anyNA(v[2:13, 30:31]))
  expect_equal(range(v, na.rm = TRUE), c(1.18, 19.97), tolerance = 0.01)
})

test_that("receivers at the shore are drawn whole", {
  skip_if_not_installed("png")
  # These 4 receivers sit 45 to 50 m inside the raster edge, where a clipped
  # panel cuts their symbols. Each receiver must change the rendered map as
  # much as the others do.
  lake <- example_lake()
  rx <- place_receivers(lake, method = "manual",
                        coords = data.frame(x = c(725045, 725450, 726750, 726300),
                                            y = c(4935600, 4935045, 4935600, 4936155)))
  p <- plot(rx, lake, labels = FALSE)
  expect_identical(p$coordinates$clip, "off")
  expect_equal(ggplot2::layer_data(p, 3)[c("x", "y")], as.data.frame(rx[c("x", "y")]),
               ignore_attr = TRUE)
  render <- function(layout) {
    f <- tempfile(fileext = ".png")
    ggplot2::ggsave(f, plot(layout, lake, labels = FALSE), width = 7, height = 5, dpi = 100,
                    bg = "white")
    png::readPNG(f)
  }
  base <- render(rx[0, ])
  changed <- vapply(seq_len(nrow(rx)), function(i) {
    sum(apply(abs(render(rx[i, ])[, , 1:3] - base[, , 1:3]) > 0.02, c(1, 2), any))
  }, numeric(1))
  expect_true(all(changed >= 0.6 * max(changed)), label = paste(changed, collapse = ", "))
})

test_that("scale_colour_tracks() gives the six track colours in order and repeats them", {
  lake <- example_lake()
  tr <- simulate_tracks(lake, crw_model(30, 10, 30), n_animals = 8, n_steps = 20, seed = 1)
  p <- ggplot2::ggplot(as.data.frame(tr),
                       ggplot2::aes(x, y, colour = factor(path_id))) +
    ggplot2::geom_path() +
    scale_colour_tracks(name = "Animal")
  colours <- unique(ggplot2::layer_data(p, 1)[c("group", "colour")])
  colours <- toupper(colours$colour[order(colours$group)])
  # Okabe-Ito less black, grey and pale yellow, in this order.
  six <- c("#E69F00", "#009E73", "#CC79A7", "#0072B2", "#D55E00", "#56B4E9")
  expect_identical(colours, c(six, six[1:2]))
  expect_false("#F0E442" %in% colours)
  # The six are entries of R's own Okabe-Ito palette.
  expect_true(all(six %in% toupper(grDevices::palette.colors(NULL, "Okabe-Ito"))))
  expect_identical(p$scales$get_scales("colour")$name, "Animal")
})

test_that("scale_colour_tracks() ties named colours to levels", {
  states <- data.frame(x = 1:4, y = 1:4,
                       state = c("Resting", "Moving", "Moving", "Resting"))
  p <- ggplot2::ggplot(states, ggplot2::aes(x, y, colour = state)) +
    ggplot2::geom_point() +
    scale_colour_tracks(colours = c(Moving = "bluish_green", Resting = "orange"))
  drawn <- toupper(ggplot2::layer_data(p, 1)$colour)
  expect_identical(drawn, c("#E69F00", "#009E73", "#009E73", "#E69F00"))
  expect_error(scale_colour_tracks(colours = c(Moving = "yellow")),
               "named by level")
  expect_error(scale_colour_tracks(colours = "orange"), "named by level")
})
