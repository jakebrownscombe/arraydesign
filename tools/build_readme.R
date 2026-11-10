#!/usr/bin/env Rscript
# ---- Build the README ----
# Knits README.Rmd to README.md on the installed package, times the knit, and
# prints the size of the README's figures and the coverage counts its example
# map shows. The one figure is made by the README's own chunk,
# man/figures/README-coverage-map-1.png. The walkthrough is the Get started
# vignette; tools/build_get_started.R knits it.
# Usage, from the package root after installing the package:
#   Rscript tools/build_readme.R

# ---- Knit ----
stopifnot(file.exists("README.Rmd"))
cat("arraydesign", as.character(packageVersion("arraydesign")), "; detfield",
    as.character(packageVersion("detfield")), "\n")
unlink(Sys.glob("man/figures/README-*.png"))
env <- new.env(parent = globalenv())
t0 <- proc.time()[["elapsed"]]
rmarkdown::render("README.Rmd", output_file = "README.md", envir = env, quiet = TRUE)
secs <- proc.time()[["elapsed"]] - t0
cat(sprintf("knit took %.1f s\n", secs))
md <- readLines("README.md", warn = FALSE)
cat("README.md lines:", length(md), "; chunk errors:", sum(grepl("^#> Error", md)),
    "; chunk warnings:", sum(grepl("^#> Warning", md)), "\n")

# ---- Figures ----
array_p <- 1 - apply(1 - env$field$p, 1, prod)
cat(sprintf("coverage map: %d water cells, %d (%.1f%%) heard with probability 0.5 or more\n",
            length(array_p), sum(array_p >= 0.5), 100 * mean(array_p >= 0.5)))
for (f in c("man/figures/README-coverage-map-1.png", "man/figures/concept_suite.png")) {
  d <- dim(png::readPNG(f))
  cat(sprintf("%-46s %4d x %4d px  %5.0f KB\n", f, d[2], d[1], file.size(f) / 1024))
}
