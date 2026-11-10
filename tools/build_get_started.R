#!/usr/bin/env Rscript
# ---- Knit the Get started vignette ----
# Knits vignettes/arraydesign.Rmd on the installed package into a temporary
# folder, times the knit (CRAN builds every vignette; it must take under a
# minute), and prints the numbers its text and figures quote.
# Needs arraydesign and detfield installed in a library that R finds (set
# R_LIBS for another library), and rmarkdown with pandoc.
# Usage, from the package root after installing the package:
#   R_LIBS=/path/to/library Rscript tools/build_get_started.R

# ---- Knit ----
rmd <- "vignettes/arraydesign.Rmd"
stopifnot(file.exists(rmd))
out_dir <- tempfile("get_started_")
dir.create(out_dir)
env <- new.env(parent = globalenv())
t0 <- proc.time()[["elapsed"]]
rmarkdown::render(rmd, output_dir = out_dir, intermediates_dir = out_dir,
                  envir = env, quiet = TRUE)
secs <- proc.time()[["elapsed"]] - t0
cat(sprintf("Get started knit took %.1f s\n", secs))
html <- readLines(file.path(out_dir, "arraydesign.html"), warn = FALSE)
n_embedded <- sum(lengths(regmatches(html, gregexpr("data:image/png;base64", html))))
cat("figures embedded:", n_embedded, "\n")

# ---- Numbers the page quotes ----
o <- env$evaluation$score$overall
sent <- env$sent
array_p <- env$array_p
rate <- env$rate_per_animal
cat(sprintf("step 2: %d of %d water cells (%.1f%%) heard with probability 0.5 or more\n",
            sum(array_p$p >= 0.5), nrow(array_p), 100 * mean(array_p$p >= 0.5)))
cat(sprintf(paste("step 4: %d transmissions, %d heard; %d receivers heard, %d none;",
                  "mean detection rate per animal %.1f%%, %d animals never heard\n"),
            nrow(sent), sum(sent$n_detections > 0), sum(env$receiver_counts$detections > 0),
            sum(env$receiver_counts$detections == 0), 100 * mean(rate), sum(rate == 0)))
cat(sprintf(paste("step 5: %d of %d animal-hours scored, median error of the mean %.0f m,",
                  "hit rate %.3f\n"),
            o$n_bins, o$n_animal_hours, o$median_error_mean_m, o$hit_rate))
