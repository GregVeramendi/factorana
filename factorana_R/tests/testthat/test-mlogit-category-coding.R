# Regression guard tests for multinomial-logit outcome coding.
#
# Like ordered probit, the C++ multinomial-logit evaluator uses the outcome
# value as a 1-indexed choice index and indexes per-choice parameter blocks by
# it, and the component does not store its (validated) data. Standard logit
# validation in define_model_component() already requires an exact contiguous
# 1..num_choices coding on the evaluation subset, which prevents the
# out-of-bounds / off-by-one indexing that the ordered-probit guard was added
# to close. These tests lock that behavior in so it cannot silently regress.

make_mlogit_df <- function(n = 600, seed = 7) {
  set.seed(seed)
  f <- rnorm(n)
  U <- cbind(0, 0.5 * f + rnorm(n), 0.2 * f + rnorm(n), -0.3 * f + rnorm(n))
  data.frame(intercept = 1, x = rnorm(n), f = f, y = max.col(U))  # y in 1..4
}

test_that("contiguous 1..K multinomial-logit outcome builds and estimates", {
  skip_on_cran()
  dat <- make_mlogit_df()
  fm <- define_factor_model(n_factors = 1)
  # Identify the single factor with two linear indicators, plus the mlogit item.
  dat$z1 <- 1.0 * dat$f + rnorm(nrow(dat), 0, 0.5)
  dat$z2 <- 0.8 * dat$f + rnorm(nrow(dat), 0, 0.5)
  mc1 <- define_model_component("z1", dat, "z1", fm, covariates = "intercept",
                                model_type = "linear", loading_normalization = 1)
  mc2 <- define_model_component("z2", dat, "z2", fm, covariates = "intercept",
                                model_type = "linear", loading_normalization = NA_real_)
  mc3 <- define_model_component("m", dat, "y", fm, covariates = c("intercept", "x"),
                                model_type = "logit", num_choices = 4,
                                loading_normalization = NA_real_)
  ms <- define_model_system(components = list(mc1, mc2, mc3), factor = fm)
  ctrl <- define_estimation_control(n_quad_points = 10, num_cores = 1)
  res <- estimate_model_rcpp(ms, dat, control = ctrl, parallel = FALSE,
                             optimizer = "nlminb", verbose = FALSE)
  expect_equal(res$convergence, 0)
})

test_that("non-contiguous mlogit outcome is rejected", {
  dat <- make_mlogit_df()
  dat$y[dat$y == 3L] <- 4L   # observed {1,2,4}
  fm <- define_factor_model(n_factors = 1)
  expect_error(
    define_model_component("m", dat, "y", fm, covariates = c("intercept", "x"),
                           model_type = "logit", num_choices = 3,
                           loading_normalization = NA_real_),
    "contiguous"
  )
})

test_that("0-indexed mlogit outcome is rejected with a recode hint", {
  dat <- make_mlogit_df()
  dat$y <- dat$y - 1L        # 0..3
  fm <- define_factor_model(n_factors = 1)
  expect_error(
    define_model_component("m", dat, "y", fm, covariates = c("intercept", "x"),
                           model_type = "logit", num_choices = 4,
                           loading_normalization = NA_real_),
    "1-indexed|minimum value"
  )
})

test_that("empty mlogit category is rejected (count or contiguity mismatch)", {
  # Interior empty -> breaks contiguity.
  dat <- make_mlogit_df()
  dat$y[dat$y == 3L] <- 4L   # {1,2,4}, declare 4 -> not contiguous
  fm <- define_factor_model(n_factors = 1)
  expect_error(
    define_model_component("m", dat, "y", fm, covariates = c("intercept", "x"),
                           model_type = "logit", num_choices = 4,
                           loading_normalization = NA_real_),
    "contiguous"
  )
  # Top category empty -> count mismatch.
  dat2 <- make_mlogit_df()
  dat2$y[dat2$y == 4L] <- 3L  # {1,2,3}, declare 4
  expect_error(
    define_model_component("m", dat2, "y", fm, covariates = c("intercept", "x"),
                           model_type = "logit", num_choices = 4,
                           loading_normalization = NA_real_),
    "does not match detected unique outcome values"
  )
})
