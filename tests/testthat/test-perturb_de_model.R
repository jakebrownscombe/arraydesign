# ---- perturb_de_model(): known answers ----

test_that("scales of 1 leave the curve unchanged", {
  m <- detfield::de_logistic(c(200, 400), c(600, 1000), depth = c(2, 20))
  p <- perturb_de_model(m)
  d <- c(0, 150, 400, 800, 1500)
  expect_equal(p$predict(d, depth_m = 7), m$predict(d, depth_m = 7), tolerance = 1e-12)
  expect_identical(attr(p, "perturbation"),
                   c(d50_scale = 1, d05_scale = 1, sigma_scale = 1))
})

test_that("a scaled logistic is the logistic through the scaled distances", {
  d <- c(0, 60, 125, 250, 400, 600, 900, 1400)
  m <- detfield::de_logistic(250, 600)
  expect_equal(perturb_de_model(m, d50_scale = 0.5)$predict(d),
               detfield::de_logistic(125, 600)$predict(d), tolerance = 1e-9)
  expect_equal(perturb_de_model(m, d05_scale = 1.5)$predict(d),
               detfield::de_logistic(250, 900)$predict(d), tolerance = 1e-9)
  expect_equal(perturb_de_model(m, 2, 2)$predict(d),
               detfield::de_logistic(500, 1200)$predict(d), tolerance = 1e-12)
})

test_that("a depth-varying logistic is scaled at every depth", {
  m <- detfield::de_logistic(c(200, 400), c(600, 1000), depth = c(2, 20))
  target <- detfield::de_logistic(c(100, 200), c(1200, 2000), depth = c(2, 20))
  p <- perturb_de_model(m, d50_scale = 0.5, d05_scale = 2)
  d <- rep(c(0, 90, 300, 700, 1600), 4)
  z <- rep(c(1, 2, 11, 25), each = 5)
  expect_equal(p$predict(d, depth_m = z), target$predict(d, depth_m = z), tolerance = 1e-8)
})

test_that("equal scales rescale distance for a curve that is not logistic", {
  tab <- data.frame(distance_m = c(0, 100, 300, 700), p = c(0.95, 0.8, 0.3, 0.02))
  m <- detfield::de_model(tab)
  p <- perturb_de_model(m, 1.7, 1.7)
  d <- c(10, 170, 400, 900, 1300)
  expect_equal(p$predict(d), m$predict(d / 1.7), tolerance = 1e-12)
})

test_that("a glm's scaled 50 and 5 percent distances fall where the scales put them", {
  range_test <- data.frame(distance_m = rep(c(50, 200, 400, 600, 800), 2),
                           depth_m = rep(c(4, 10), each = 5),
                           heard = c(92, 71, 40, 14, 4, 95, 80, 52, 22, 7))
  fit <- stats::glm(cbind(heard, 100 - heard) ~ distance_m + depth_m, family = binomial,
                    data = range_test)
  m <- detfield::de_model(fit)
  p <- perturb_de_model(m, d50_scale = 0.8, d05_scale = 1.3)
  for (z in c(4, 10)) {
    b <- stats::coef(fit)
    d50 <- (b[1] + b[3] * z) / -b[2]
    d05 <- (stats::qlogis(0.05) - b[1] - b[3] * z) / b[2]
    expect_equal(stats::plogis(p$predict(c(0.8 * d50, 1.3 * d05), depth_m = z)),
                 c(0.5, 0.05), tolerance = 1e-8, ignore_attr = TRUE)
  }
})

test_that("modulators reach the curve and set their own anchors", {
  f <- function(distance_m, noise, ...) 3 - distance_m / 100 - noise
  m <- detfield::de_model(f, modulators = list(noise = 0))
  p <- perturb_de_model(m, d50_scale = 0.5, d05_scale = 1)
  d05 <- function(noise) 100 * (3 - noise - stats::qlogis(0.05))
  got <- stats::plogis(p$predict(c(150, 100, d05(0), d05(1)), noise = c(0, 1, 0, 1)))
  expect_equal(got, c(0.5, 0.5, 0.05, 0.05), tolerance = 1e-8)
  expect_identical(p$modulators, m$modulators)
})

test_that("sigma is scaled and NA stays NA", {
  m <- detfield::de_model(detfield::de_logistic(250, 600),
                          sigma = list(station_hour = 1.2, station_day = NA))
  p <- perturb_de_model(m, sigma_scale = 2)
  expect_equal(p$sigma$station_hour, 2.4)
  expect_true(is.na(p$sigma$station_day))
  expect_match(p$source, "sigma x 2) of: de_logistic, d50 250 m", fixed = TRUE)
  expect_s3_class(p, "de_model")
})

test_that("blocked, uses and transmitter models carry over", {
  m <- detfield::de_model(detfield::de_logistic(c(200, 400), c(600, 1000), depth = c(2, 20)),
                          blocked = -3)
  p <- perturb_de_model(m, 0.5)
  expect_identical(p$blocked, -3)
  expect_identical(p$uses, m$uses)
  expect_identical(p$tag_models, m$tag_models)
})

test_that("a field from a perturbed model equals the field of the scaled logistic", {
  lake <- example_lake()
  rx <- place_receivers(lake, n = 4)
  dist <- detfield::water_distances(lake, rx[c("receiver_serial", "x", "y")])
  f1 <- detfield::de_field(dist, perturb_de_model(detfield::de_logistic(250, 600), 0.5))
  f2 <- detfield::de_field(dist, detfield::de_logistic(125, 600))
  expect_equal(f1$eta, f2$eta, tolerance = 1e-9)
})

# ---- perturb_de_model(): refusals ----

test_that("bad arguments are refused", {
  m <- detfield::de_logistic(250, 600)
  expect_error(perturb_de_model(list()), "must be a de_model")
  expect_error(perturb_de_model(m, d50_scale = 0), "above 0")
  expect_error(perturb_de_model(m, d05_scale = c(1, 2)), "above 0")
  expect_error(perturb_de_model(m, sigma_scale = -1), "0 or more")
})

test_that("scaled anchors that swap order are refused when the model is made", {
  m <- detfield::de_logistic(250, 600)
  expect_error(perturb_de_model(m, d50_scale = 3),
               "must stay beyond the scaled 50 percent distance")
  expect_error(perturb_de_model(m, d50_scale = 2, d05_scale = 0.5),
               "must stay beyond the scaled 50 percent distance")
  expect_s3_class(perturb_de_model(m, d50_scale = 2, d05_scale = 1), "de_model")
})

test_that("a depth curve is checked when it is asked, at each depth", {
  m <- detfield::de_logistic(c(200, 400), c(600, 1000), depth = c(2, 20))
  # At 2 m the scaled 50 percent distance (750 m) passes the 5 percent distance (600 m).
  p <- perturb_de_model(m, d50_scale = 3.75)
  expect_error(p$predict(100, depth_m = 2), "must stay beyond the scaled 50 percent distance")
})

test_that("a curve without a 50 or 5 percent distance cannot take unequal scales", {
  low <- detfield::de_model(function(distance_m) -1 - distance_m / 100)
  expect_error(perturb_de_model(low, 0.5), "no 50 percent distance")
  expect_equal(perturb_de_model(low, 2, 2)$predict(200), -2)
  flat <- detfield::de_model(function(distance_m) 2 + 0 * distance_m)
  expect_error(perturb_de_model(flat, 0.5), "does not fall to 0.05")
  expect_error(perturb_de_model(detfield::de_logistic(250, 600), 0.5)$predict("a"),
               "must be numeric")
})
