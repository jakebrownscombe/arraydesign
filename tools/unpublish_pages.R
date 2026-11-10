#!/usr/bin/env Rscript
# ---- Remove maintainer pages from the built site ----
# Usage: Rscript tools/unpublish_pages.R [package-dir]
# CONTRACT.md (and DEPLOY.md where a package has one) are written for the
# maintainers. pkgdown 2.0.9 renders every top-level .md file as a page and
# does not read .Rbuildignore, so run this as the last step after
# pkgdown::build_site(). It deletes the pages, then rebuilds the search index
# and sitemap from the pages left. The same step as detfield's
# tools/unpublish_contract.R.

args <- commandArgs(trailingOnly = TRUE)
pkg_dir <- normalizePath(if (length(args)) args[1] else ".")
docs <- file.path(pkg_dir, "docs")
pages <- c("CONTRACT.html", "DEPLOY.html")

present <- file.exists(file.path(docs, pages))
cat("present before:", paste(pages, present, sep = " ", collapse = "; "), "\n")
unlink(file.path(docs, pages))
# Point rdrr.io links to suite functions at the suite's own reference pages.
system2("python3", c("-I", file.path(dirname(pkg_dir), "tools", "relink_reference.py"), docs))
pkgdown::build_search(pkg_dir)
# The sitemap still names the deleted pages after build_search(), which does
# not rewrite it.
pkgdown:::build_sitemap(pkgdown::as_pkgdown(pkg_dir))

mentions <- function(file) {
  lines <- readLines(file, warn = FALSE)
  sum(vapply(pages, function(p) sum(grepl(p, lines, fixed = TRUE)), numeric(1)))
}
cat("sitemap.xml lines naming the pages:", mentions(file.path(docs, "sitemap.xml")), "\n")
cat("search.json lines naming the pages:", mentions(file.path(docs, "search.json")), "\n")
html <- list.files(docs, pattern = "\\.html$", recursive = TRUE, full.names = TRUE)
linking <- html[vapply(html, mentions, numeric(1)) > 0]
cat("pages linking them:", length(linking), "of", length(html), "\n")
