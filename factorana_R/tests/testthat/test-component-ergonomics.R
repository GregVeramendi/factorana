# Ergonomics of define_model_component(): the covariates default and the
# factor_index convenience for assigning an item to a single factor.

test_that("covariates defaults to NULL (no boilerplate covariates = NULL)", {
  dat <- data.frame(intercept = 1, yo = sample(1:4, 80, TRUE))
  fm <- define_factor_model(n_factors = 1)
  # oprobit measurement item with no covariates: omitting `covariates` must work
  mc <- define_model_component("mo", dat, "yo", fm, model_type = "oprobit",
                               num_choices = 4, loading_normalization = 1)
  expect_equal(length(mc$covariates), 0L)
})

test_that("factor_index assigns an item to a single factor", {
  dat <- data.frame(y1 = rnorm(50), y2 = rnorm(50))
  fm <- define_factor_model(n_factors = 2)

  a <- define_model_component("a", dat, "y1", fm, factor_index = 2,
                              loading_normalization = NA_real_)
  expect_equal(a$loading_normalization, c(0, NA_real_))

  b <- define_model_component("b", dat, "y2", fm, factor_index = 1,
                              loading_normalization = 1)
  expect_equal(b$loading_normalization, c(1, 0))

  # default loading_normalization (NULL) under factor_index = free loading there
  d <- define_model_component("d", dat, "y1", fm, factor_index = 1, intercept = FALSE)
  expect_equal(d$loading_normalization, c(NA_real_, 0))

  # identical to the explicit length-n_factors vector
  v <- define_model_component("a", dat, "y1", fm, loading_normalization = c(0, NA_real_), intercept = FALSE)
  expect_identical(a$loading_normalization, v$loading_normalization)
})

test_that("factor_index validates its inputs", {
  dat <- data.frame(y1 = rnorm(50))
  fm <- define_factor_model(n_factors = 2)
  expect_error(define_model_component("x", dat, "y1", fm, factor_index = 3,
                                      loading_normalization = NA_real_), "1\\.\\.2")
  expect_error(define_model_component("x", dat, "y1", fm, factor_index = 1,
                                      loading_normalization = c(1, 0)), "single value")
})

test_that("factor_index gives the same fit as the explicit vector", {
  skip_on_cran()
  set.seed(1); n <- 400; f1 <- rnorm(n); f2 <- rnorm(n)
  dat <- data.frame(
    y1 = 1.0 * f1 + rnorm(n, 0, 0.5), y2 = 0.8 * f1 + rnorm(n, 0, 0.5),
    y3 = 1.0 * f2 + rnorm(n, 0, 0.5), y4 = 0.9 * f2 + rnorm(n, 0, 0.5))
  fm <- define_factor_model(n_factors = 2)
  ctrl <- define_estimation_control(n_quad_points = 8, num_cores = 1)

  # explicit vectors
  ms_vec <- define_model_system(factor = fm, components = list(
    define_model_component("y1", dat, "y1", fm, loading_normalization = c(1, 0), intercept = FALSE),
    define_model_component("y2", dat, "y2", fm, loading_normalization = c(NA_real_, 0), intercept = FALSE),
    define_model_component("y3", dat, "y3", fm, loading_normalization = c(0, 1), intercept = FALSE),
    define_model_component("y4", dat, "y4", fm, loading_normalization = c(0, NA_real_), intercept = FALSE)))
  # factor_index shorthand
  ms_idx <- define_model_system(factor = fm, components = list(
    define_model_component("y1", dat, "y1", fm, factor_index = 1, loading_normalization = 1, intercept = FALSE),
    define_model_component("y2", dat, "y2", fm, factor_index = 1, loading_normalization = NA_real_, intercept = FALSE),
    define_model_component("y3", dat, "y3", fm, factor_index = 2, loading_normalization = 1, intercept = FALSE),
    define_model_component("y4", dat, "y4", fm, factor_index = 2, loading_normalization = NA_real_, intercept = FALSE)))

  f_vec <- estimate_model_rcpp(ms_vec, dat, control = ctrl, optimizer = "nlminb",
                               parallel = FALSE, verbose = FALSE)
  f_idx <- estimate_model_rcpp(ms_idx, dat, control = ctrl, optimizer = "nlminb",
                               parallel = FALSE, verbose = FALSE)
  expect_equal(f_vec$convergence, 0)
  expect_equal(f_idx$convergence, 0)
  expect_equal(unname(f_idx$estimates), unname(f_vec$estimates), tolerance = 1e-6)
})
