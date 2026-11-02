# A fit produced with previous_stage must itself be usable as previous_stage. define_model_system()
# stores bookkeeping fields in the factor model it returns, and the factor-model identity check
# used to compare those too, so a third stage failed against a second-stage fit (estimate_hjv.R
# Stage 4). The setup mirrors HJV: variances fixed at 1, free loadings with exclusions.

.chain_data <- function(n = 400, seed = 3) {
  set.seed(seed)
  f1 <- rnorm(n); f2 <- rnorm(n)
  data.frame(intercept = 1, x = rnorm(n),
             a1 = 1.0 * f1 + rnorm(n, 0, 0.5), a2 = 0.8 * f1 + rnorm(n, 0, 0.5), a3 = 1.2 * f1 + rnorm(n, 0, 0.5),
             b1 = 1.0 * f2 + rnorm(n, 0, 0.5), b2 = 0.9 * f2 + rnorm(n, 0, 0.5), b3 = 1.1 * f2 + rnorm(n, 0, 0.5),
             y2 = 0.3 + 0.5 * f1 + 0.4 * f2 + rnorm(n, 0, 0.6),
             y3 = -0.2 + 0.6 * f1 - 0.3 * f2 + rnorm(n, 0, 0.6))
}

test_that("previous_stage chains across three stages", {
  skip_on_cran()
  dat <- .chain_data()
  fm <- fix_factor_param(define_factor_model(n_factors = 2), c(factor_var_1 = 1, factor_var_2 = 1))
  mk <- function(nm, y, norm) define_model_component(nm, dat, y, fm, covariates = "intercept",
                                                     model_type = "linear", loading_normalization = norm)
  f1only <- c(NA_real_, 0); f2only <- c(0, NA_real_); both <- c(NA_real_, NA_real_)
  meas <- list(mk("a1", "a1", f1only), mk("a2", "a2", f1only), mk("a3", "a3", f1only),
               mk("b1", "b1", f2only), mk("b2", "b2", f2only), mk("b3", "b3", f2only))
  ctrl <- define_estimation_control(n_quad_points = 8, num_cores = 1)

  ms1 <- define_model_system(components = meas, factor = fm)
  fit1 <- estimate_model_rcpp(ms1, dat, control = ctrl, optimizer = "nlminb", parallel = FALSE, verbose = FALSE)
  expect_equal(fit1$convergence, 0)

  ms2 <- define_model_system(components = list(mk("y2", "y2", both)), factor = fm, previous_stage = fit1)
  # bookkeeping covers every factor, not only the first
  expect_true(isTRUE(ms2$factor$variance_fixed))
  expect_equal(names(ms2$factor$variance_value), c("factor_var_1", "factor_var_2"))
  fit2 <- estimate_model_rcpp(ms2, dat, control = ctrl, optimizer = "nlminb", parallel = FALSE, verbose = FALSE)
  expect_equal(fit2$convergence, 0)

  # the third stage takes the second-stage fit as previous_stage (this used to error)
  ms3 <- define_model_system(components = list(mk("y3", "y3", both)), factor = fm, previous_stage = fit2)
  expect_s3_class(ms3, "model_system")
  fixed <- ms3$previous_stage_info$fixed_param_names
  expect_true(all(grepl("^(a[123]|b[123]|y2)_|^factor_var_", fixed)))
  expect_true(all(c("y2_intercept", "y2_loading_1", "y2_sigma") %in% fixed))
  fit3 <- estimate_model_rcpp(ms3, dat, control = ctrl, optimizer = "nlminb", parallel = FALSE, verbose = FALSE)
  expect_equal(fit3$convergence, 0)
  expect_true(all(is.finite(fit3$std_errors[c("y3_intercept", "y3_loading_1", "y3_loading_2", "y3_sigma")])))
  # stage-2 values carried into stage 3 unchanged
  expect_equal(unname(fit3$estimates["y2_loading_1"]), unname(fit2$estimates["y2_loading_1"]))
  expect_equal(unname(fit3$estimates["a2_loading_1"]), unname(fit1$estimates["a2_loading_1"]))
})

test_that("a fresh factor model that differs in specification is still rejected", {
  dat <- .chain_data(n = 200)
  fm <- fix_factor_param(define_factor_model(n_factors = 2), c(factor_var_1 = 1, factor_var_2 = 1))
  mc <- define_model_component("a1", dat, "a1", fm, covariates = "intercept", model_type = "linear",
                               loading_normalization = c(NA_real_, 0))
  fake_fit <- structure(list(model_system = define_model_system(components = list(mc), factor = fm),
                             estimates = c(factor_var_1 = 1, factor_var_2 = 1, a1_intercept = 0, a1_loading_1 = 1, a1_sigma = 0.5),
                             std_errors = c(factor_var_1 = 0, factor_var_2 = 0, a1_intercept = 0.1, a1_loading_1 = 0.1, a1_sigma = 0.1)),
                        class = "factorana_result")
  fm3 <- define_factor_model(n_factors = 3)
  mc3 <- define_model_component("y2", dat, "y2", fm3, covariates = "intercept", model_type = "linear",
                                loading_normalization = c(NA_real_, NA_real_, 1))
  expect_error(define_model_system(components = list(mc3), factor = fm3, previous_stage = fake_fit), "must be identical")
})
