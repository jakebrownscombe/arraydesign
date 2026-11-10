# ---- Vignette scripts ----
# Each vignette's code is also a plain R script, vignettes/<name>.R, written
# by tools/purl_vignettes.R and committed, with a copy in inst/scripts/ that
# installs. A script must match its vignette and run on the files the
# package installs. In a source checkout the scripts are read from
# vignettes/; in R CMD check, from the scripts folder of the installed
# package.

slow_scripts <- "design-an-array.R"

source_vignettes <- function() {
  testthat::test_path("..", "..", "vignettes")
}

vignette_scripts <- function() {
  src <- source_vignettes()
  dir <- if (dir.exists(src)) src else system.file("scripts", package = "arraydesign")
  normalizePath(list.files(dir, pattern = "\\.R$", full.names = TRUE))
}

test_that("every vignette has a script", {
  scripts <- vignette_scripts()
  skip_if(length(scripts) == 0, "no vignette scripts beside the tests or in the installed package")
  expect_gt(length(scripts), 0)
  src <- source_vignettes()
  skip_if_not(dir.exists(src), "vignette sources are not beside the tests")
  rmds <- list.files(src, pattern = "\\.Rmd$")
  expect_setequal(sub("\\.Rmd$", ".R", rmds), basename(scripts))
})

test_that("each committed script matches its vignette", {
  skip_if_not_installed("knitr")
  skip_if_not_installed("rmarkdown")
  src <- source_vignettes()
  tool <- test_path("..", "..", "tools", "purl_vignettes.R")
  skip_if_not(dir.exists(src) && file.exists(tool), "vignette sources are not beside the tests")
  env <- new.env()
  sys.source(tool, envir = env)
  out_dir <- withr::local_tempdir()
  fresh <- env$purl_vignettes(dirname(src), out_dir = out_dir, install_dir = NULL)
  installed <- file.path(dirname(src), "inst", "scripts")
  for (f in fresh) {
    for (committed in file.path(c(src, installed), basename(f))) {
      expect_true(file.exists(committed), info = committed)
      expect_identical(readLines(committed), readLines(f), info = committed)
    }
  }
})

for (script in vignette_scripts()) {
  test_that(paste("vignette script", basename(script), "runs without error or warning"), {
    if (basename(script) %in% slow_scripts) skip_on_cran()
    withr::local_dir(withr::local_tempdir())
    grDevices::pdf(NULL)
    withr::defer(grDevices::dev.off())
    # The scripts show the errors a reader will meet inside try(); their text
    # goes to a null connection, not the test log.
    sink_to <- file(nullfile(), open = "w")
    withr::defer(close(sink_to))
    withr::local_options(try.outFile = sink_to)
    warned <- character()
    # print.eval = TRUE prints each figure as a reader's session would, so a
    # ggplot that fails or warns only when drawn fails here too.
    expect_error(
      withCallingHandlers(
        utils::capture.output(source(script, local = new.env(), print.eval = TRUE)),
        warning = function(w) {
          warned <<- c(warned, conditionMessage(w))
          invokeRestart("muffleWarning")
        },
        message = function(m) invokeRestart("muffleMessage")
      ),
      NA
    )
    expect_identical(warned, character())
  })
}

# A reader who starts from the vignettes meets every export there. An export
# counts as shown when a vignette script calls it, not when its name appears
# in a file name.
test_that("every export is called in a vignette script", {
  code <- unlist(lapply(vignette_scripts(), readLines, warn = FALSE))
  code <- code[!grepl("^\\s*#", code)]
  # example_lake() is exported for the tests only. An export used only by
  # tests is not shown to users.
  exports <- setdiff(sort(getNamespaceExports("arraydesign")), "example_lake")
  called <- vapply(exports, function(f) any(grepl(paste0("\\b", f, "\\("), code)), logical(1))
  expect_identical(exports[!called], character(0))
})
