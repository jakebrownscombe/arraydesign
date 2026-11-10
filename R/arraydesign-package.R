# The detfield import loads detfield's namespace at install, so the version
# floor in DESCRIPTION is checked then.
#' @keywords internal
#' @importFrom ggplot2 .data
#' @importFrom detfield positionr_ring_colours
"_PACKAGE"

# data.table syntax works inside a package only when its namespace imports
# data.table or sets this flag. The packages call data.table through `::`,
# so they set the flag.
.datatable.aware <- TRUE # nolint: object_name_linter.

# Symbols that data.table and the table contract (CONTRACT.md) use inside
# `[.data.table` calls. R CMD check would report them as undefined globals.
utils::globalVariables(c(
  ".I", ".N", ":=", "datetime", "receiver_serial", "transmitter",
  "animal_id", "station_id", "sensor_raw", "sensor_value", "sensor_unit",
  "deploy_datetime", "recover_datetime", "deploy_lat", "deploy_lon",
  "deploy_z", "x", "y"
))
