# The hard-dependencies-only CI job runs R CMD check without testthat
# installed, so the test run is skipped when testthat is absent.
if (requireNamespace("testthat", quietly = TRUE)) {
  library(testthat)
  library(arraydesign)
  test_check("arraydesign")
}
