# The estimation log, the result and the printed results state the integration of a fit:
# the quadrature points per factor, and with adaptive quadrature the threshold.
#
# Design: one factor measured by three linear indicators, then a second stage with one linear
# outcome, estimated with the factor scores of the first stage (adaptive quadrature).

.qio_data <- function(n = 300, seed = 11) {
  set.seed(seed)
  f <- rnorm(n)
  data.frame(intercept = 1,
             y1 = f + rnorm(n, 0, 0.6), y2 = 0.8 * f + rnorm(n, 0, 0.7), y3 = 0.6 * f + rnorm(n, 0, 0.5),
             w = 0.5 + 0.7 * f + rnorm(n, 0, 0.8))
}

.qio_stage1 <- function(dat, fm) {
  ind <- function(nm, norm) define_model_component(nm, dat, nm, fm, covariates = "intercept",
                                                   model_type = "linear", loading_normalization = norm)
  define_model_system(components = list(ind("y1", 1), ind("y2", NA_real_), ind("y3", NA_real_)), factor = fm)
}

test_that("full quadrature: the log, the result and the printed results give the points per factor", {
  skip_on_cran()
  dat <- .qio_data()
  fm <- define_factor_model(n_factors = 1)
  ms <- .qio_stage1(dat, fm)
  ctrl <- define_estimation_control(n_quad_points = 6, num_cores = 1)
  log <- character(0)
  out <- capture.output(
    fit <- withCallingHandlers(
      estimate_model_rcpp(ms, dat, control = ctrl, parallel = FALSE, verbose = TRUE),
      message = function(m) { log <<- c(log, conditionMessage(m)); invokeRestart("muffleMessage") }))
  # in the log, above the benchmark of the evaluation times
  at <- grep("^Quadrature: 6 points per factor", log)
  expect_length(at, 1)
  expect_lt(at, grep("Benchmarking computation times", log)[1])
  expect_identical(fit$quadrature, list(n_points = 6L, adaptive = FALSE, threshold = NA_real_))
  expect_output(print(fit), "Quadrature: +6 points per factor")
  expect_output(print(summary(fit)), "Quadrature: +6 points per factor")
})

test_that("adaptive quadrature: the output gives the most points per factor and the threshold", {
  skip_on_cran()
  dat <- .qio_data()
  fm <- define_factor_model(n_factors = 1)
  ctrl <- define_estimation_control(n_quad_points = 8, num_cores = 1)
  fit1 <- estimate_model_rcpp(.qio_stage1(dat, fm), dat, control = ctrl, parallel = FALSE, verbose = FALSE)
  fs <- estimate_factorscores_rcpp(result = fit1, data = dat, control = ctrl, verbose = FALSE)
  ms2 <- define_model_system(
    components = list(define_model_component("w", dat, "w", fm, covariates = "intercept", model_type = "linear",
                                             loading_normalization = NA_real_)),
    factor = fm, previous_stage = fit1)
  ctrl2 <- define_estimation_control(n_quad_points = 8, num_cores = 1, adaptive_integration = TRUE, adapt_int_thresh = 0.25)
  fit2 <- estimate_model_rcpp(ms2, dat, control = ctrl2, parallel = FALSE, verbose = FALSE,
                              factor_scores = as.matrix(fs[, "factor_1", drop = FALSE]),
                              factor_ses = as.matrix(fs[, "se_factor_1", drop = FALSE]),
                              factor_vars = as.numeric(fit1$estimates["factor_var_1"]))
  expect_identical(fit2$quadrature, list(n_points = 8L, adaptive = TRUE, threshold = 0.25))
  expect_output(print(fit2), "Quadrature: +adaptive, at most 8 points per factor \\(threshold 0.25\\)")
})

test_that("a result saved by an earlier version prints without a quadrature line", {
  skip_on_cran()
  dat <- .qio_data()
  fm <- define_factor_model(n_factors = 1)
  fit <- estimate_model_rcpp(.qio_stage1(dat, fm), dat, control = define_estimation_control(n_quad_points = 6, num_cores = 1),
                             parallel = FALSE, verbose = FALSE)
  fit$quadrature <- NULL
  out <- capture.output(print(fit))
  expect_false(any(grepl("Quadrature", out)))
  expect_true(any(grepl("^Parameters:", out)))
})
