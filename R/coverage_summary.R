#' Water cells that at least one receiver hears on its own
#'
#' For each threshold, counts the water cells where at least one receiver
#' detects a transmission with probability at or above the threshold, and
#' divides by the number of water cells in the field's grid. A water cell
#' that the field does not hold counts as not covered.
#'
#' @param field A `de_field` from [detfield::de_field()], or a list of the
#'   same shape (see [simulate_detections()]).
#' @param listening `NULL` for every receiver in the field, or the receivers
#'   to count, such as those deployed and listening on one day: the field's
#'   receiver keys, which are `station_id` values when the field is keyed by
#'   station, or `receiver_serial` values, which a station-keyed field maps to
#'   their stations through its `lookup`.
#' @param thresholds Detection probabilities, each above 0 and at most 1. A
#'   threshold of 0 would count every cell as covered whatever the receivers.
#'
#' @return A data frame, one row per threshold: `threshold`, `n_cells` (water
#'   cells covered), `n_water_cells` (all water cells in the grid, the
#'   denominator), `share` (`n_cells / n_water_cells`) and `n_receivers` (the
#'   receivers counted).
#' @seealso [simulate_detections()]
#' @export
#' @examples
#' stoney <- terra::rast(system.file("extdata", "stoney_depth.tif",
#'                                   package = "arraydesign"))
#' planned <- read.csv(system.file("extdata", "planned_stations.csv",
#'                                 package = "arraydesign"))
#' rx <- place_receivers(stoney, method = "manual", coords = planned)
#' # A deployments table: one receiver per station for two days.
#' deployments <- data.frame(
#'   receiver_serial = rx$receiver_serial, station_id = rx$station_id,
#'   x = rx$x, y = rx$y, deploy_z = rx$deploy_z,
#'   deploy_datetime = as.POSIXct("2024-12-31", tz = "UTC"),
#'   recover_datetime = as.POSIXct("2025-01-02", tz = "UTC")
#' )
#' range_model <- detfield::de_logistic(d50 = c(250, 750), d05 = c(500, 1500),
#'                                      depth = c(1, 35))
#' # Pairs more than 4 km apart through water are left out. At the deepest
#' # receiver the curve gives them a detection probability below 1e-7, and
#' # the distances take a third of the time.
#' field <- detfield::de_field(detfield::water_distances(stoney, deployments,
#'                                                     max_distance = 4000),
#'                             range_model)
#' coverage_summary(field)
#' # Two stations, named by station or by the receivers deployed there.
#' coverage_summary(field, listening = c("S001", "S002"))
#' coverage_summary(field, listening = c("R001", "R002"))
coverage_summary <- function(field, listening = NULL, thresholds = c(0.5, 0.8)) {
  fg <- check_field(field)
  stop_unless(is.numeric(thresholds) && length(thresholds) && !anyNA(thresholds) &&
                all(thresholds > 0 & thresholds <= 1),
              "`thresholds` must be probabilities above 0 and at most 1.")
  cols <- listening_columns(field, listening)
  stop_unless(all(fg$wet[field$cells] %in% TRUE),
              "`field$cells` holds cells that are land or outside the grid.")
  eta <- field$eta[, cols, drop = FALSE]
  if (anyNA(eta)) {
    stop_classed("field", "The field has ", sum(is.na(eta)),
                 " missing values among the receivers counted; no probability is assumed.")
  }
  best <- rep(0, nrow(eta))
  if (length(cols)) best <- stats::plogis(do.call(pmax, as.data.frame(eta)))
  n_water <- sum(fg$wet)
  n_cells <- vapply(thresholds, function(t) sum(best >= t), integer(1))
  data.frame(threshold = thresholds, n_cells = n_cells, n_water_cells = n_water,
             share = n_cells / n_water, n_receivers = length(cols))
}

# Columns of the field for the receivers counted. A station-keyed field also
# takes serials, mapped to their stations through its lookup.
listening_columns <- function(field, listening) {
  if (is.null(listening)) return(seq_along(field$receivers))
  lk <- field$lookup
  if (identical(field$key, "station_id") && is.data.frame(lk)) {
    as_station <- lk$station_id[match(listening, lk$receiver_serial)]
    listening <- ifelse(listening %in% field$receivers | is.na(as_station), listening,
                        as_station)
  }
  unknown <- setdiff(listening, field$receivers)
  stop_unless(!length(unknown), "`listening` names receivers the field does not hold: ",
              paste(unknown, collapse = ", "), ".")
  which(field$receivers %in% listening)
}
