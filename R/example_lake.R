#' A small synthetic lake for unit tests
#'
#' A small synthetic lake for unit tests; examples and articles use the Stoney
#' Lake raster in `inst/extdata/stoney_depth.tif`.
#'
#' @return A [terra::SpatRaster] of depth in metres, `NA` on land.
#' @keywords internal
#' @export
example_lake <- function() {
  nc <- 60L
  nr <- 40L
  col <- rep(seq_len(nc), times = nr)
  row <- rep(seq_len(nr), each = nc)
  # Two bowls, one each side of the peninsula, so depth falls towards every shore.
  half_centre <- ifelse(col <= 30, 15.5, 45.5)
  depth <- 1 + 19 * (1 - ((col - half_centre) / 15)^2) * (1 - ((row - 20.5) / 20)^2)
  depth <- pmax(depth, 1)
  land <- col == 1 | col == nc | row == 1 | row == nr | (col %in% 30:31 & row >= 14)
  depth[land] <- NA
  r <- terra::rast(
    nrows = nr, ncols = nc, xmin = 725000, xmax = 725000 + 30 * nc,
    ymin = 4935000, ymax = 4935000 + 30 * nr, crs = "EPSG:32617"
  )
  terra::values(r) <- depth
  names(r) <- "depth_m"
  r
}
