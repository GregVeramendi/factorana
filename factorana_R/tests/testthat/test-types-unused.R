# A factor model with n_types > 1 but no component that uses types must estimate like the n_types = 1 model.
# Before 1.8.0 initialize_factor_model_cpp() allocated the type-probability parameters whenever n_types > 1 while
# the R side added them only when some component had use_types = TRUE, so the first stage of a typed pipeline
# (measurement system without types, outcome stage with types) failed with "Fixed values size mismatch".

.types_unused_data <- function(n = 400, seed = 11) {
  set.seed(seed)
  f <- rnorm(n)
  data.frame(intercept = 1, x1 = rnorm(n),
             y1 = 0.2 + 1.0 * f + rnorm(n, 0, 0.6),
             y2 = -0.1 + 0.8 * f + rnorm(n, 0, 0.7),
             y3 = 0.3 + 0.6 * f + rnorm(n, 0, 0.5))
}

.types_unused_system <- function(dat, n_types) {
  fm <- define_factor_model(n_factors = 1, n_types = n_types)
  mk <- function(nm, norm) define_model_component(nm, dat, nm, fm, covariates = c("intercept", "x1"),
                                                   model_type = "linear", loading_normalization = norm)
  define_model_system(components = list(mk("y1", 1), mk("y2", NA_real_), mk("y3", NA_real_)), factor = fm)
}

test_that("n_types > 1 without typed components estimates like n_types = 1", {
  dat <- .types_unused_data()
  ctrl <- define_estimation_control(n_quad_points = 8, num_cores = 1)
  fit1 <- estimate_model_rcpp(.types_unused_system(dat, 1), dat, control = ctrl, optimizer = "nlminb",
                              parallel = FALSE, verbose = FALSE)
  fit3 <- estimate_model_rcpp(.types_unused_system(dat, 3), dat, control = ctrl, optimizer = "nlminb",
                              parallel = FALSE, verbose = FALSE)
  expect_equal(fit1$convergence, 0)
  expect_equal(fit3$convergence, 0)
  expect_false(any(grepl("type", names(fit3$estimates))))
  expect_equal(names(fit3$estimates), names(fit1$estimates))
  expect_equal(fit3$estimates, fit1$estimates, tolerance = 1e-6)
  expect_equal(fit3$loglik, fit1$loglik, tolerance = 1e-8)
})

test_that("a typed outcome component still brings the type probability parameters", {
  dat <- .types_unused_data()
  dat$y4 <- 0.5 + 0.4 * (dat$y1 - 0.2) + rnorm(nrow(dat), 0, 0.8)
  fm <- define_factor_model(n_factors = 1, n_types = 2)
  mk <- function(nm, norm, types = FALSE) define_model_component(nm, dat, nm, fm, covariates = c("intercept", "x1"),
                                                                  model_type = "linear", loading_normalization = norm,
                                                                  use_types = types)
  ms <- define_model_system(components = list(mk("y1", 1), mk("y2", NA_real_), mk("y3", NA_real_), mk("y4", NA_real_, TRUE)),
                            factor = fm)
  init <- initialize_parameters(ms, dat, verbose = FALSE)
  expect_true(all(c("typeprob_2_intercept", "type_2_loading_1", "y4_type_2_intercept") %in% init$param_names))
  fm_ptr <- initialize_factor_model_cpp(ms, as.matrix(dat), n_quad = 8L, init_params = init$init_params)
  expect_true(is.finite(evaluate_loglik_only_cpp(fm_ptr, init$init_params)))
})
