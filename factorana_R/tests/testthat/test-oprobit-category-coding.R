# Guard tests for ordered-probit outcome coding.
#
# The C++ ordered-probit evaluator uses the raw outcome value as the category
# index (1..num_choices) and indexes threshold parameters by it. The component
# does not store its data and the estimation pipeline does not recode the
# outcome, so a non-contiguous or 0-indexed coding would reach C++ unchanged
# and either read threshold parameters out of bounds (value > num_choices) or
# be off by one (0-indexed), silently corrupting the likelihood.
# define_model_component() must reject these at construction with a clear
# message. A contiguous 1..num_choices coding must continue to work.

make_oprobit_df <- function(n = 600, seed = 11, cuts = c(-0.8, 0.0, 0.8)) {
  set.seed(seed)
  f <- rnorm(n)
  mk <- function(load) findInterval(load * f + rnorm(n), cuts) + 1L  # 1..(length(cuts)+1)
  data.frame(intercept = 1, f = f,
             y1 = mk(1.0), y2 = mk(0.9), y3 = mk(1.1))
}

test_that("contiguous 1..K oprobit outcome builds and estimates (no guard error)", {
  skip_on_cran()
  dat <- make_oprobit_df()
  fm <- define_factor_model(n_factors = 1)
  mc1 <- define_model_component("m1", dat, "y1", fm, model_type = "oprobit",
                                covariates = NULL, num_choices = 4, loading_normalization = 1)
  mc2 <- define_model_component("m2", dat, "y2", fm, model_type = "oprobit",
                                covariates = NULL, num_choices = 4, loading_normalization = NA_real_)
  mc3 <- define_model_component("m3", dat, "y3", fm, model_type = "oprobit",
                                covariates = NULL, num_choices = 4, loading_normalization = NA_real_)
  ms <- define_model_system(components = list(mc1, mc2, mc3), factor = fm)
  ctrl <- define_estimation_control(n_quad_points = 12, num_cores = 1)
  res <- estimate_model_rcpp(ms, dat, control = ctrl, parallel = FALSE,
                             optimizer = "nlminb", verbose = FALSE)
  expect_equal(res$convergence, 0)
})

test_that("oprobit value exceeding num_choices is rejected (out-of-bounds guard)", {
  dat <- make_oprobit_df()
  # Collapse category 3 into 4 so observed values are {1,2,4} but declare 3 cats.
  dat$y1[dat$y1 == 3L] <- 4L
  fm <- define_factor_model(n_factors = 1)
  expect_error(
    define_model_component("m1", dat, "y1", fm, model_type = "oprobit",
                           covariates = NULL, num_choices = 3, loading_normalization = 1),
    "exceeding num_choices"
  )
})

test_that("0-indexed oprobit outcome is rejected with a recode hint", {
  dat <- make_oprobit_df()
  dat$y1 <- dat$y1 - 1L  # now 0..3
  fm <- define_factor_model(n_factors = 1)
  expect_error(
    define_model_component("m1", dat, "y1", fm, model_type = "oprobit",
                           covariates = NULL, num_choices = 4, loading_normalization = 1),
    "1-indexed|0-indexed|minimum value"
  )
})

test_that("empty interior oprobit category is rejected and names the category", {
  dat <- make_oprobit_df()
  # Remove category 3 (shift 4->3 leaves a gap: observed {1,2,3} mapped from {1,2,4}).
  # Build values {1,2,4,5} with a true gap at 3 by widening to 5 categories.
  set.seed(99)
  f <- rnorm(800)
  y <- findInterval(1.0 * f + rnorm(800), c(-1.0, -0.3, 0.3, 1.0)) + 1L  # 1..5
  y[y == 3L] <- 2L  # empty out category 3 -> observed {1,2,4,5}
  dat2 <- data.frame(intercept = 1, y1 = y)
  fm <- define_factor_model(n_factors = 1)
  expect_error(
    define_model_component("m1", dat2, "y1", fm, model_type = "oprobit",
                           covariates = NULL, num_choices = 5, loading_normalization = 1),
    "no observations in category 3"
  )
})

test_that("data declared with too few categories (more observed than num_choices) errors", {
  dat <- make_oprobit_df()  # y1 in 1..4
  fm <- define_factor_model(n_factors = 1)
  expect_error(
    define_model_component("m1", dat, "y1", fm, model_type = "oprobit",
                           covariates = NULL, num_choices = 3, loading_normalization = 1),
    "exceeding num_choices"
  )
})
