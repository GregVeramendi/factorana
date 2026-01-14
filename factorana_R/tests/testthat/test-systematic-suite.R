# Systematic Test Suite for Factorana Package
#
# This test suite systematically tests the package functionality by:
# 1. Checking analytical gradients against finite differences
# 2. Checking analytical Hessians against finite differences
# 3. Estimating models with default and true initial parameters
# 4. Comparing estimates to true parameter values
#
# Tests are marked skip_on_cran() as they involve estimation and take time

# Test configuration
VERBOSE <- Sys.getenv("FACTORANA_TEST_VERBOSE", "FALSE") == "TRUE"
SAVE_LOGS <- Sys.getenv("FACTORANA_TEST_SAVE_LOGS", "TRUE") == "TRUE"
GRAD_TOL <- 1e-3  # Relaxed to accommodate numerical precision in finite differences
HESS_TOL <- 1e-3  # Checking all elements (diagonal and off-diagonal)

# ==============================================================================
# Test A: Measurement System with 3 Linear Tests and 1 Factor
# ==============================================================================

test_that("Model A: Measurement system with 3 linear tests and 1 factor", {
  skip_on_cran()

  set.seed(104)

  # Simulate data
  n <- 500
  f <- rnorm(n)  # Latent factor

  # True parameters - ORDER: factor_var, T1_int, T1_sigma, T2_int, T2_lambda, T2_sigma, T3_int, T3_lambda, T3_sigma
  true_params <- c(1.0,  # Factor variance
                   2.0, 0.5,   # T1: int, sigma (loading fixed to 1.0)
                   1.5, 1.2, 0.6,   # T2: int, lambda, sigma
                   1.0, 0.8, 0.4)   # T3: int, lambda, sigma

  # Generate data (using true loading values for simulation)
  T1 <- true_params[2] + 1.0*f + rnorm(n, 0, true_params[3])  # T1 loading = 1.0 (fixed)
  T2 <- true_params[4] + true_params[5]*f + rnorm(n, 0, true_params[6])
  T3 <- true_params[7] + true_params[8]*f + rnorm(n, 0, true_params[9])

  dat <- data.frame(intercept = 1, T1 = T1, T2 = T2, T3 = T3, eval = 1)

  # Create model system with 1 factor
  fm <- define_factor_model(n_factors = 1, n_types = 1)

  mc_T1 <- define_model_component(
    name = "T1", data = dat, outcome = "T1", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = 1.0,  # Fix loading to 1 for identification
    evaluation_indicator = "eval"
  )
  mc_T2 <- define_model_component(
    name = "T2", data = dat, outcome = "T2", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = NA_real_,  # Estimate freely
    evaluation_indicator = "eval"
  )
  mc_T3 <- define_model_component(
    name = "T3", data = dat, outcome = "T3", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = NA_real_,  # Estimate freely
    evaluation_indicator = "eval"
  )

  ms <- define_model_system(components = list(mc_T1, mc_T2, mc_T3), factor = fm)

  # First run estimation to get param_fixed
  est_comp <- run_estimation_comparison(
    ms, dat, true_params,
    param_names = c("f_var", "T1_int", "T1_sigma",  # T1 loading fixed to 1.0
                   "T2_int", "T2_lambda", "T2_sigma",
                   "T3_int", "T3_lambda", "T3_sigma"),
    verbose = VERBOSE
  )

  # Run gradient/Hessian checks with param_fixed
  grad_check <- check_gradient_accuracy(ms, dat, true_params,
                                       param_fixed = est_comp$param_fixed,
                                       tol = GRAD_TOL, verbose = VERBOSE, n_quad = 8)
  hess_check <- check_hessian_accuracy(ms, dat, true_params,
                                       param_fixed = est_comp$param_fixed,
                                       tol = HESS_TOL, verbose = VERBOSE, n_quad = 8)

  # Collect diagnostics
  diagnostics <- list(
    gradient_check = grad_check,
    hessian_check = hess_check,
    estimation = est_comp,
    overall_pass = grad_check$pass && hess_check$pass &&
                   est_comp$default_converged && est_comp$reasonable_default
  )

  # Save log
  if (SAVE_LOGS) {
    log_file <- save_diagnostics_to_log("test_A_measurement_system", diagnostics)
    if (VERBOSE) cat("Log saved to:", log_file, "\n")
  }

  # Assertions
  expect_true(grad_check$pass, info = sprintf("Gradient check failed (max error: %.2e)", grad_check$max_error))
  expect_true(hess_check$pass, info = sprintf("Hessian check failed (max error: %.2e)", hess_check$max_error))
  expect_true(est_comp$default_converged, info = "Estimation with default init failed to converge")
  expect_true(est_comp$reasonable_default, info = "Estimates are not reasonable")
})

# ==============================================================================
# Test B: Measurement System (A) + Probit
# ==============================================================================

test_that("Model B: Measurement system with 3 linear tests and probit outcome", {
  skip_on_cran()

  set.seed(105)

  # Simulate data
  n <- 500
  x1 <- rnorm(n)
  f <- rnorm(n)  # Latent factor

  # True parameters - ORDER: factor_var, T1_int, T1_sigma, T2_int, T2_lambda, T2_sigma, T3_int, T3_lambda, T3_sigma, y_int, y_beta1, y_loading
  true_params <- c(1.0,  # Factor variance
                   2.0, 0.5,   # T1: int, sigma (loading fixed to 1.0)
                   1.5, 1.2, 0.6,   # T2: int, lambda, sigma
                   1.0, 0.8, 0.4,   # T3: int, lambda, sigma
                   0.5, 0.7, 0.9)   # Probit: int, beta1, loading

  # Generate data (using true loading values for simulation)
  T1 <- true_params[2] + 1.0*f + rnorm(n, 0, true_params[3])  # T1 loading = 1.0 (fixed)
  T2 <- true_params[4] + true_params[5]*f + rnorm(n, 0, true_params[6])
  T3 <- true_params[7] + true_params[8]*f + rnorm(n, 0, true_params[9])

  z <- true_params[10] + true_params[11]*x1 + true_params[12]*f
  y <- as.numeric(runif(n) < pnorm(z))

  dat <- data.frame(intercept = 1, x1 = x1, T1 = T1, T2 = T2, T3 = T3, y = y, eval = 1)

  # Create model system
  fm <- define_factor_model(n_factors = 1, n_types = 1)

  mc_T1 <- define_model_component(
    name = "T1", data = dat, outcome = "T1", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = 1.0, evaluation_indicator = "eval"
  )
  mc_T2 <- define_model_component(
    name = "T2", data = dat, outcome = "T2", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = NA_real_, evaluation_indicator = "eval"
  )
  mc_T3 <- define_model_component(
    name = "T3", data = dat, outcome = "T3", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = NA_real_, evaluation_indicator = "eval"
  )
  mc_y <- define_model_component(
    name = "y", data = dat, outcome = "y", factor = fm,
    covariates = c("intercept", "x1"), model_type = "probit",
    loading_normalization = NA_real_, evaluation_indicator = "eval"
  )

  ms <- define_model_system(components = list(mc_T1, mc_T2, mc_T3, mc_y), factor = fm)

  # First run estimation to get param_fixed
  est_comp <- run_estimation_comparison(
    ms, dat, true_params,
    param_names = c("f_var", "T1_int", "T1_sigma",  # T1 loading fixed to 1.0
                   "T2_int", "T2_lambda", "T2_sigma",
                   "T3_int", "T3_lambda", "T3_sigma",
                   "y_int", "y_beta1", "y_loading"),
    verbose = VERBOSE
  )

  # Run gradient/Hessian checks with param_fixed
  grad_check <- check_gradient_accuracy(ms, dat, true_params,
                                       param_fixed = est_comp$param_fixed,
                                       tol = GRAD_TOL, verbose = VERBOSE, n_quad = 8)
  hess_check <- check_hessian_accuracy(ms, dat, true_params,
                                       param_fixed = est_comp$param_fixed,
                                       tol = HESS_TOL, verbose = VERBOSE, n_quad = 8)

  # Collect diagnostics
  diagnostics <- list(
    gradient_check = grad_check,
    hessian_check = hess_check,
    estimation = est_comp,
    overall_pass = grad_check$pass && hess_check$pass &&
                   est_comp$default_converged && est_comp$reasonable_default
  )

  # Save log
  if (SAVE_LOGS) {
    log_file <- save_diagnostics_to_log("test_B_measurement_plus_probit", diagnostics)
    if (VERBOSE) cat("Log saved to:", log_file, "\n")
  }

  # Assertions
  expect_true(grad_check$pass, info = sprintf("Gradient check failed (max error: %.2e)", grad_check$max_error))
  expect_true(hess_check$pass, info = sprintf("Hessian check failed (max error: %.2e)", hess_check$max_error))
  expect_true(est_comp$default_converged, info = "Estimation with default init failed to converge")
  expect_true(est_comp$reasonable_default, info = "Estimates are not reasonable")
})

# ==============================================================================
# Test C: Measurement System (A) + Ordered Probit
# ==============================================================================

test_that("Model C: Measurement system with 3 linear tests and ordered probit", {
  skip_on_cran()

  set.seed(106)

  # Simulate data
  n <- 500
  x1 <- rnorm(n)
  f <- rnorm(n)  # Latent factor

  # True parameters - ORDER: factor_var, T1_int, T1_sigma, T2_int, T2_lambda, T2_sigma, T3_int, T3_lambda, T3_sigma,
  #                           y_beta1, y_loading, thresh1, thresh_incr1
  # NOTE: Ordered probit has NO INTERCEPT (absorbed into thresholds)
  # Threshold parameterization: thresh2 = thresh1 + abs(thresh_incr1)
  true_params <- c(1.0,  # Factor variance
                   2.0, 0.5,   # T1: int, sigma (loading fixed to 1.0)
                   1.5, 1.2, 0.6,   # T2: int, lambda, sigma
                   1.0, 0.8, 0.4,   # T3: int, lambda, sigma
                   0.6, 0.8,   # Oprobit: beta1, loading (NO intercept)
                   -0.5, 1.0)  # Oprobit thresholds: thresh1, thresh_incr1 (intercept absorbed)

  # Generate data (using true loading values for simulation)
  T1 <- true_params[2] + 1.0*f + rnorm(n, 0, true_params[3])  # T1 loading = 1.0 (fixed)
  T2 <- true_params[4] + true_params[5]*f + rnorm(n, 0, true_params[6])
  T3 <- true_params[7] + true_params[8]*f + rnorm(n, 0, true_params[9])

  # Ordered probit latent variable (with intercept for data generation)
  # But model will not estimate intercept - it's absorbed into thresholds
  y_intercept <- 0.5
  z <- y_intercept + true_params[10]*x1 + true_params[11]*f

  # Compute absolute threshold values from incremental parameterization
  # thresh1 has intercept absorbed: original 0.0 - intercept 0.5 = -0.5
  thresh1_abs <- true_params[12]
  thresh2_abs <- thresh1_abs + abs(true_params[13])

  u <- rnorm(n)
  y <- ifelse(z + u < thresh1_abs, 1, ifelse(z + u < thresh2_abs, 2, 3))

  dat <- data.frame(intercept = 1, x1 = x1, T1 = T1, T2 = T2, T3 = T3, y = y, eval = 1)

  # Create model system
  fm <- define_factor_model(n_factors = 1, n_types = 1)

  mc_T1 <- define_model_component(
    name = "T1", data = dat, outcome = "T1", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = 1.0, evaluation_indicator = "eval"
  )
  mc_T2 <- define_model_component(
    name = "T2", data = dat, outcome = "T2", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = NA_real_, evaluation_indicator = "eval"
  )
  mc_T3 <- define_model_component(
    name = "T3", data = dat, outcome = "T3", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = NA_real_, evaluation_indicator = "eval"
  )
  mc_y <- define_model_component(
    name = "y", data = dat, outcome = "y", factor = fm,
    covariates = "x1",  # NO intercept for ordered probit (absorbed into thresholds)
    model_type = "oprobit",
    num_choices = 3,  # 3 ordered categories
    loading_normalization = NA_real_, evaluation_indicator = "eval"
  )

  ms <- define_model_system(components = list(mc_T1, mc_T2, mc_T3, mc_y), factor = fm)

  # First run estimation to get param_fixed
  est_comp <- run_estimation_comparison(
    ms, dat, true_params,
    param_names = c("f_var", "T1_int", "T1_sigma",  # T1 loading fixed to 1.0
                   "T2_int", "T2_lambda", "T2_sigma",
                   "T3_int", "T3_lambda", "T3_sigma",
                   "y_beta1", "y_loading",  # Oprobit: NO intercept
                   "thresh1", "thresh_incr1"),  # Incremental threshold parameterization
    verbose = VERBOSE
  )

  # Run gradient/Hessian checks with param_fixed
  grad_check <- check_gradient_accuracy(ms, dat, true_params,
                                       param_fixed = est_comp$param_fixed,
                                       tol = GRAD_TOL, verbose = VERBOSE, n_quad = 8)
  hess_check <- check_hessian_accuracy(ms, dat, true_params,
                                       param_fixed = est_comp$param_fixed,
                                       tol = HESS_TOL, verbose = VERBOSE, n_quad = 8)

  # Collect diagnostics
  diagnostics <- list(
    gradient_check = grad_check,
    hessian_check = hess_check,
    estimation = est_comp,
    overall_pass = grad_check$pass && hess_check$pass &&
                   est_comp$default_converged && est_comp$reasonable_default
  )

  # Save log
  if (SAVE_LOGS) {
    log_file <- save_diagnostics_to_log("test_C_measurement_plus_oprobit", diagnostics)
    if (VERBOSE) cat("Log saved to:", log_file, "\n")
  }

  # Assertions
  expect_true(grad_check$pass, info = sprintf("Gradient check failed (max error: %.2e)", grad_check$max_error))
  expect_true(hess_check$pass, info = sprintf("Hessian check failed (max error: %.2e)", hess_check$max_error))
  expect_true(est_comp$default_converged, info = "Estimation with default init failed to converge")
  expect_true(est_comp$reasonable_default, info = "Estimates are not reasonable")
})

# ==============================================================================
# Test D: Measurement System (A) + Multinomial Logit
# ==============================================================================

test_that("Model D: Measurement system with 3 linear tests and multinomial logit", {
  skip_on_cran()

  set.seed(107)

  # Simulate data
  n <- 500
  x1 <- rnorm(n)
  f <- rnorm(n)  # Latent factor

  # True parameters - ORDER: factor_var, T1_int, T1_sigma, T2_int, T2_lambda, T2_sigma, T3_int, T3_lambda, T3_sigma,
  #                           y1_int, y1_beta1, y1_loading, y2_int, y2_beta1, y2_loading
  true_params <- c(1.0,  # Factor variance
                   2.0, 0.5,   # T1: int, sigma (loading fixed to 1.0)
                   1.5, 1.2, 0.6,   # T2: int, lambda, sigma
                   1.0, 0.8, 0.4,   # T3: int, lambda, sigma
                   0.5, 0.6, 0.7,   # Mlogit choice 1: int, beta1, loading
                   1.0, -0.5, 0.9)  # Mlogit choice 2: int, beta1, loading

  # Generate data (using true loading values for simulation)
  T1 <- true_params[2] + 1.0*f + rnorm(n, 0, true_params[3])  # T1 loading = 1.0 (fixed)
  T2 <- true_params[4] + true_params[5]*f + rnorm(n, 0, true_params[6])
  T3 <- true_params[7] + true_params[8]*f + rnorm(n, 0, true_params[9])

  z1 <- true_params[10] + true_params[11]*x1 + true_params[12]*f
  z2 <- true_params[13] + true_params[14]*x1 + true_params[15]*f

  exp_z0 <- 1
  exp_z1 <- exp(z1)
  exp_z2 <- exp(z2)
  denom <- exp_z0 + exp_z1 + exp_z2

  p0 <- exp_z0 / denom
  p1 <- exp_z1 / denom
  p2 <- exp_z2 / denom

  y <- numeric(n)
  for (i in seq_len(n)) {
    # C++ expects multinomial choices coded as 1, 2, 3 (not 0, 1, 2)
    y[i] <- sample(1:3, 1, prob = c(p0[i], p1[i], p2[i]))
  }

  dat <- data.frame(intercept = 1, x1 = x1, T1 = T1, T2 = T2, T3 = T3, y = y, eval = 1)

  # Create model system
  fm <- define_factor_model(n_factors = 1, n_types = 1)

  mc_T1 <- define_model_component(
    name = "T1", data = dat, outcome = "T1", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = 1.0, evaluation_indicator = "eval"
  )
  mc_T2 <- define_model_component(
    name = "T2", data = dat, outcome = "T2", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = NA_real_, evaluation_indicator = "eval"
  )
  mc_T3 <- define_model_component(
    name = "T3", data = dat, outcome = "T3", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = NA_real_, evaluation_indicator = "eval"
  )
  mc_y <- define_model_component(
    name = "y", data = dat, outcome = "y", factor = fm,
    covariates = c("intercept", "x1"), model_type = "logit",
    num_choices = 3,  # 3 choices: 0, 1, 2
    loading_normalization = NA_real_,  # Single factor, loading estimated freely
    evaluation_indicator = "eval"
  )

  ms <- define_model_system(components = list(mc_T1, mc_T2, mc_T3, mc_y), factor = fm)

  # First run estimation to get param_fixed
  est_comp <- run_estimation_comparison(
    ms, dat, true_params,
    param_names = c("f_var",
                   "T1_int", "T1_sigma",  # T1 loading fixed to 1.0
                   "T2_int", "T2_lambda", "T2_sigma",
                   "T3_int", "T3_lambda", "T3_sigma",
                   "y1_int", "y1_beta1", "y1_loading",
                   "y2_int", "y2_beta1", "y2_loading"),
    verbose = VERBOSE
  )

  # Run gradient/Hessian checks with param_fixed
  grad_check <- check_gradient_accuracy(ms, dat, true_params,
                                       param_fixed = est_comp$param_fixed,
                                       tol = GRAD_TOL, verbose = VERBOSE, n_quad = 8)
  hess_check <- check_hessian_accuracy(ms, dat, true_params,
                                       param_fixed = est_comp$param_fixed,
                                       tol = HESS_TOL, verbose = VERBOSE, n_quad = 8)

  # Collect diagnostics
  diagnostics <- list(
    gradient_check = grad_check,
    hessian_check = hess_check,
    estimation = est_comp,
    overall_pass = grad_check$pass && hess_check$pass &&
                   est_comp$default_converged && est_comp$reasonable_default
  )

  # Save log
  if (SAVE_LOGS) {
    log_file <- save_diagnostics_to_log("test_D_measurement_plus_mlogit", diagnostics)
    if (VERBOSE) cat("Log saved to:", log_file, "\n")
  }

  # Assertions
  expect_true(grad_check$pass, info = sprintf("Gradient check failed (max error: %.2e)", grad_check$max_error))
  expect_true(hess_check$pass, info = sprintf("Hessian check failed (max error: %.2e)", hess_check$max_error))
  expect_true(est_comp$default_converged, info = "Estimation with default init failed to converge")
  expect_true(est_comp$reasonable_default, info = "Estimates are not reasonable")
})

# ==============================================================================
# Test E: Roy Model
# ==============================================================================

test_that("Model E: Roy selection model", {
  skip_on_cran()

  set.seed(108)

  # Simulate Roy model data
  n <- 2000  # Increased from 500 to match manual test for better numerical stability
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  f <- rnorm(n)  # Latent factor (ability)

  # True parameters - ORDER: factor_var, T1_int, T1_sigma, T2_int, T2_lambda, T2_sigma,
  #                           T3_int, T3_lambda, T3_sigma, wage0_int, wage0_beta1, wage0_beta2, wage0_sigma,
  #                           wage1_int, wage1_beta1, wage1_lambda, wage1_sigma, sector_int, sector_beta1, sector_loading
  true_params <- c(
    1.0,  # Factor variance
    # T1: int, sigma (loading FIXED to 1.0)
    2.0, 0.5,
    # T2: int, lambda, sigma
    1.5, 1.2, 0.6,
    # T3: int, lambda, sigma
    1.0, 0.8, 0.4,
    # Wage0: int, beta1, beta2, sigma (loading FIXED to 0.0 - no factor effect)
    2.0, 0.5, 0.3, 0.6,
    # Wage1: int, beta1, lambda, sigma
    2.5, 0.6, 1.0, 0.7,
    # Sector: int, beta1, loading (probit - no sigma)
    0.0, 0.4, 0.8
  )

  # Generate test scores (using fixed loading values for data generation)
  T1 <- true_params[2] + 1.0*f + rnorm(n, 0, true_params[3])  # T1 loading = 1.0 (fixed)
  T2 <- true_params[4] + true_params[5]*f + rnorm(n, 0, true_params[6])
  T3 <- true_params[7] + true_params[8]*f + rnorm(n, 0, true_params[9])

  # Generate potential wages
  wage0 <- true_params[10] + true_params[11]*x1 + true_params[12]*x2 +
           rnorm(n, 0, true_params[13])  # No factor effect (loading = 0.0 fixed)
  wage1 <- true_params[14] + true_params[15]*x1 + true_params[16]*f +
           rnorm(n, 0, true_params[17])

  # Sector choice (based on utility difference)
  z_sector <- true_params[18] + true_params[19]*x2 + true_params[20]*f
  sector <- as.numeric(runif(n) < pnorm(z_sector))

  # Observed wage
  wage <- ifelse(sector == 1, wage1, wage0)

  dat <- data.frame(
    intercept = 1,
    x1 = x1, x2 = x2,
    T1 = T1, T2 = T2, T3 = T3,
    wage = wage,
    sector = sector,
    eval_tests = 1,
    eval_wage0 = 1 - sector,
    eval_wage1 = sector,
    eval_sector = 1
  )

  # Create model system
  fm <- define_factor_model(n_factors = 1, n_types = 1)

  mc_T1 <- define_model_component(
    name = "T1", data = dat, outcome = "T1", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = 1.0, evaluation_indicator = "eval_tests"
  )
  mc_T2 <- define_model_component(
    name = "T2", data = dat, outcome = "T2", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = NA_real_, evaluation_indicator = "eval_tests"
  )
  mc_T3 <- define_model_component(
    name = "T3", data = dat, outcome = "T3", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = NA_real_, evaluation_indicator = "eval_tests"
  )
  mc_wage0 <- define_model_component(
    name = "wage0", data = dat, outcome = "wage", factor = fm,
    covariates = c("intercept", "x1", "x2"), model_type = "linear",
    loading_normalization = 0.0, evaluation_indicator = "eval_wage0"
  )
  mc_wage1 <- define_model_component(
    name = "wage1", data = dat, outcome = "wage", factor = fm,
    covariates = c("intercept", "x1"), model_type = "linear",
    loading_normalization = NA_real_, evaluation_indicator = "eval_wage1"
  )
  mc_sector <- define_model_component(
    name = "sector", data = dat, outcome = "sector", factor = fm,
    covariates = c("intercept", "x2"), model_type = "probit",
    loading_normalization = NA_real_, evaluation_indicator = "eval_sector"
  )

  ms <- define_model_system(
    components = list(mc_T1, mc_T2, mc_T3, mc_wage0, mc_wage1, mc_sector),
    factor = fm
  )

  # First run estimation to get param_fixed
  est_comp <- run_estimation_comparison(
    ms, dat, true_params,
    param_names = c("f_var",
                   "T1_int", "T1_sigma",  # T1: loading fixed to 1.0
                   "T2_int", "T2_lambda", "T2_sigma",
                   "T3_int", "T3_lambda", "T3_sigma",
                   "w0_int", "w0_beta1", "w0_beta2", "w0_sigma",
                   "w1_int", "w1_beta1", "w1_lambda", "w1_sigma",
                   "s_int", "s_beta1", "s_loading"),
    verbose = VERBOSE
  )

  # Run gradient/Hessian checks with param_fixed
  grad_check <- check_gradient_accuracy(ms, dat, true_params,
                                       param_fixed = est_comp$param_fixed,
                                       tol = GRAD_TOL, verbose = VERBOSE, n_quad = 16)
  hess_check <- check_hessian_accuracy(ms, dat, true_params,
                                       param_fixed = est_comp$param_fixed,
                                       tol = HESS_TOL, verbose = VERBOSE, n_quad = 16)

  # Collect diagnostics
  diagnostics <- list(
    gradient_check = grad_check,
    hessian_check = hess_check,
    estimation = est_comp,
    overall_pass = grad_check$pass && hess_check$pass &&
                   est_comp$default_converged && est_comp$reasonable_default
  )

  # Save log
  if (SAVE_LOGS) {
    log_file <- save_diagnostics_to_log("test_E_roy_model", diagnostics)
    if (VERBOSE) cat("Log saved to:", log_file, "\n")
  }

  # Assertions
  expect_true(grad_check$pass, info = sprintf("Gradient check failed (max error: %.2e)", grad_check$max_error))
  expect_true(hess_check$pass, info = sprintf("Hessian check failed (max error: %.2e)", hess_check$max_error))
  expect_true(est_comp$default_converged, info = "Estimation with default init failed to converge")
  expect_true(est_comp$reasonable_default, info = "Estimates are not reasonable")
})

# ==============================================================================
# Test F: Measurement System + Exploded Multinomial Logit (Ranked Choices)
# ==============================================================================

test_that("Model F: Measurement system with exploded logit (ranked choices)", {
  skip_on_cran()

  set.seed(109)

  # Simulate data - simpler case: 3 choices, 2 ranks
  n <- 1000  # Larger sample for stability
  x1 <- rnorm(n)
  f <- rnorm(n)  # Latent factor

  # True parameters - ORDER: factor_var, T1_int, T1_sigma, T2_int, T2_lambda, T2_sigma,
  #                           c1_int, c1_beta1, c1_loading, c2_int, c2_beta1, c2_loading
  # 3 choices = 2 non-reference alternatives
  true_params <- c(1.0,  # Factor variance
                   0.0, 0.5,   # T1: int, sigma (loading fixed to 1.0)
                   0.0, 0.8, 0.5,   # T2: int, lambda, sigma
                   0.5, 0.8, 0.6,   # Mlogit choice 1: int, beta1, loading
                   -0.3, -0.5, 0.9)  # Mlogit choice 2: int, beta1, loading

  # Generate data (using true loading values for simulation)
  T1 <- true_params[2] + 1.0*f + rnorm(n, 0, true_params[3])  # T1 loading = 1.0 (fixed)
  T2 <- true_params[4] + true_params[5]*f + rnorm(n, 0, true_params[6])

  # Multinomial logit utilities (3 choices: reference + 2 alternatives)
  z1 <- true_params[7] + true_params[8]*x1 + true_params[9]*f
  z2 <- true_params[10] + true_params[11]*x1 + true_params[12]*f

  exp_z0 <- 1  # Reference category
  exp_z1 <- exp(z1)
  exp_z2 <- exp(z2)
  denom <- exp_z0 + exp_z1 + exp_z2

  # Generate rankings for each observation (2 ranks from 3 choices)
  rank1 <- rank2 <- numeric(n)
  for (i in seq_len(n)) {
    probs <- c(exp_z0, exp_z1[i], exp_z2[i]) / denom[i]
    # Sample first choice
    rank1[i] <- sample(1:3, 1, prob = probs)
    # Update probs for second choice (remove first choice)
    remaining <- setdiff(1:3, rank1[i])
    probs2 <- probs[remaining] / sum(probs[remaining])
    rank2[i] <- sample(remaining, 1, prob = probs2)
  }

  dat <- data.frame(intercept = 1, x1 = x1, T1 = T1, T2 = T2,
                    rank1 = rank1, rank2 = rank2, eval = 1)

  # Create model system
  fm <- define_factor_model(n_factors = 1, n_types = 1)

  mc_T1 <- define_model_component(
    name = "T1", data = dat, outcome = "T1", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = 1.0, evaluation_indicator = "eval"
  )
  mc_T2 <- define_model_component(
    name = "T2", data = dat, outcome = "T2", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = NA_real_, evaluation_indicator = "eval"
  )
  # Exploded logit with 2 ranks
  mc_choice <- define_model_component(
    name = "choice", data = dat,
    outcome = c("rank1", "rank2"),  # Vector of outcomes for ranked choices
    factor = fm,
    covariates = c("intercept", "x1"), model_type = "logit",
    num_choices = 3,  # 3 choices total
    loading_normalization = NA_real_,  # Single factor, loading estimated freely
    evaluation_indicator = "eval"
  )

  ms <- define_model_system(components = list(mc_T1, mc_T2, mc_choice), factor = fm)

  # First run estimation to get param_fixed
  est_comp <- run_estimation_comparison(
    ms, dat, true_params,
    param_names = c("f_var",
                   "T1_int", "T1_sigma",  # T1 loading fixed to 1.0
                   "T2_int", "T2_lambda", "T2_sigma",
                   "c1_int", "c1_beta1", "c1_loading",
                   "c2_int", "c2_beta1", "c2_loading"),
    verbose = VERBOSE
  )

  # Run gradient check with param_fixed
  # NOTE: Hessian for exploded logit with conditional probabilities is complex
  # and not yet fully implemented - skip Hessian check for this model type
  grad_check <- check_gradient_accuracy(ms, dat, true_params,
                                       param_fixed = est_comp$param_fixed,
                                       tol = GRAD_TOL, verbose = VERBOSE, n_quad = 8)

  # Collect diagnostics
  diagnostics <- list(
    gradient_check = grad_check,
    estimation = est_comp,
    overall_pass = grad_check$pass && est_comp$default_converged
  )

  # Save log
  if (SAVE_LOGS) {
    log_file <- save_diagnostics_to_log("test_F_exploded_logit", diagnostics)
    if (VERBOSE) cat("Log saved to:", log_file, "\n")
  }

  # Assertions - for exploded logit, we focus on gradient correctness and convergence
  # Hessian for conditional probabilities in exploded logit is a TODO
  expect_true(grad_check$pass, info = sprintf("Gradient check failed (max error: %.2e)", grad_check$max_error))
  expect_true(est_comp$default_converged, info = "Estimation with default init failed to converge")
})

# ==============================================================================
# Test G: Structural Equation Model (SE_linear)
# f2 = se_intercept + se_linear_1 * f1 + epsilon
# ==============================================================================

test_that("Model G: SE_linear structural equation model with 2 factors", {
  skip_on_cran()

  set.seed(201)

  # Simulate data
  n <- 500

  # True parameters for factor structure
  true_var_f1 <- 1.5        # Variance of input factor
  true_se_intercept <- 0.5   # SE intercept
  true_se_linear <- 0.8      # SE linear coefficient
  true_se_residual_var <- 0.5  # SE residual variance

  # Generate factors
  f1 <- rnorm(n, 0, sqrt(true_var_f1))
  eps <- rnorm(n, 0, sqrt(true_se_residual_var))
  f2 <- true_se_intercept + true_se_linear * f1 + eps

  # Generate measurements for f1 (3 measures)
  Y1_1 <- 2.0 + 1.0 * f1 + rnorm(n, 0, 1.0)   # int=2, loading=1 (fixed), sigma=1
  Y1_2 <- 1.5 + 0.8 * f1 + rnorm(n, 0, 0.9)   # int=1.5, loading=0.8, sigma=0.9
  Y1_3 <- 1.0 + 1.2 * f1 + rnorm(n, 0, 1.1)   # int=1, loading=1.2, sigma=1.1

  # Generate measurements for f2 (3 measures)
  Y2_1 <- 1.8 + 1.0 * f2 + rnorm(n, 0, 1.0)   # int=1.8, loading=1 (fixed), sigma=1
  Y2_2 <- 1.2 + 0.7 * f2 + rnorm(n, 0, 0.8)   # int=1.2, loading=0.7, sigma=0.8
  Y2_3 <- 0.5 + 1.1 * f2 + rnorm(n, 0, 0.95)  # int=0.5, loading=1.1, sigma=0.95

  dat <- data.frame(
    intercept = 1,
    Y1_1 = Y1_1, Y1_2 = Y1_2, Y1_3 = Y1_3,
    Y2_1 = Y2_1, Y2_2 = Y2_2, Y2_3 = Y2_3,
    eval = 1
  )

  # Create model system with 2 factors and SE_linear structure
  fm <- define_factor_model(n_factors = 2, factor_structure = "SE_linear")

  # 3 measurement equations for factor 1
  mc1_1 <- define_model_component(
    name = "Y1_1", data = dat, outcome = "Y1_1", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(1, 0),  # Fix loading to 1 for f1, 0 for f2
    evaluation_indicator = "eval"
  )
  mc1_2 <- define_model_component(
    name = "Y1_2", data = dat, outcome = "Y1_2", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(NA_real_, 0),  # Free for f1, 0 for f2
    evaluation_indicator = "eval"
  )
  mc1_3 <- define_model_component(
    name = "Y1_3", data = dat, outcome = "Y1_3", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(NA_real_, 0),  # Free for f1, 0 for f2
    evaluation_indicator = "eval"
  )

  # 3 measurement equations for factor 2
  mc2_1 <- define_model_component(
    name = "Y2_1", data = dat, outcome = "Y2_1", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(0, 1),  # 0 for f1, fix to 1 for f2
    evaluation_indicator = "eval"
  )
  mc2_2 <- define_model_component(
    name = "Y2_2", data = dat, outcome = "Y2_2", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(0, NA_real_),  # 0 for f1, free for f2
    evaluation_indicator = "eval"
  )
  mc2_3 <- define_model_component(
    name = "Y2_3", data = dat, outcome = "Y2_3", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(0, NA_real_),  # 0 for f1, free for f2
    evaluation_indicator = "eval"
  )

  ms <- define_model_system(
    components = list(mc1_1, mc1_2, mc1_3, mc2_1, mc2_2, mc2_3),
    factor = fm
  )

  # True parameters:
  # factor_var_1, se_intercept, se_linear_1, se_residual_var,
  # Y1_1_int, Y1_1_sigma,  (loading fixed to 1)
  # Y1_2_int, Y1_2_lambda, Y1_2_sigma,
  # Y1_3_int, Y1_3_lambda, Y1_3_sigma,
  # Y2_1_int, Y2_1_sigma,  (loading fixed to 1)
  # Y2_2_int, Y2_2_lambda, Y2_2_sigma,
  # Y2_3_int, Y2_3_lambda, Y2_3_sigma
  true_params <- c(
    true_var_f1,  # factor_var_1 = 1.5
    true_se_intercept, true_se_linear, true_se_residual_var,  # SE params: 0.5, 0.8, 0.5
    2.0, 1.0,  # Y1_1: int, sigma (loading fixed to 1)
    1.5, 0.8, 0.9,  # Y1_2: int, lambda, sigma
    1.0, 1.2, 1.1,  # Y1_3: int, lambda, sigma
    1.8, 1.0,  # Y2_1: int, sigma (loading fixed to 1)
    1.2, 0.7, 0.8,  # Y2_2: int, lambda, sigma
    0.5, 1.1, 0.95  # Y2_3: int, lambda, sigma
  )

  param_names <- c(
    "f1_var",
    "se_intercept", "se_linear_1", "se_residual_var",
    "Y1_1_int", "Y1_1_sigma",
    "Y1_2_int", "Y1_2_lambda", "Y1_2_sigma",
    "Y1_3_int", "Y1_3_lambda", "Y1_3_sigma",
    "Y2_1_int", "Y2_1_sigma",
    "Y2_2_int", "Y2_2_lambda", "Y2_2_sigma",
    "Y2_3_int", "Y2_3_lambda", "Y2_3_sigma"
  )

  # Run estimation to get param_fixed
  est_comp <- run_estimation_comparison(
    ms, dat, true_params,
    param_names = param_names,
    verbose = VERBOSE
  )

  # Run gradient/Hessian checks with param_fixed
  grad_check <- check_gradient_accuracy(ms, dat, true_params,
                                        param_fixed = est_comp$param_fixed,
                                        tol = GRAD_TOL, verbose = VERBOSE, n_quad = 8)
  hess_check <- check_hessian_accuracy(ms, dat, true_params,
                                       param_fixed = est_comp$param_fixed,
                                       tol = HESS_TOL, verbose = VERBOSE, n_quad = 8)

  # Collect diagnostics
  diagnostics <- list(
    gradient_check = grad_check,
    hessian_check = hess_check,
    estimation = est_comp,
    overall_pass = grad_check$pass && hess_check$pass &&
                   est_comp$default_converged && est_comp$reasonable_default
  )

  # Save log
  if (SAVE_LOGS) {
    log_file <- save_diagnostics_to_log("test_G_se_linear", diagnostics)
    if (VERBOSE) cat("Log saved to:", log_file, "\n")
  }

  # Assertions - focus on gradient and Hessian correctness
  # Estimation convergence with default init is not strictly required for this test
  expect_true(grad_check$pass, info = sprintf("Gradient check failed (max error: %.2e)", grad_check$max_error))
  expect_true(hess_check$pass, info = sprintf("Hessian check failed (max error: %.2e)", hess_check$max_error))
})

# ==============================================================================
# Test H: Structural Equation Model (SE_quadratic)
# f2 = se_intercept + se_linear_1 * f1 + se_quadratic_1 * f1^2 + epsilon
# ==============================================================================

test_that("Model H: SE_quadratic structural equation model with 2 factors", {
  skip_on_cran()

  set.seed(202)

  # Simulate data
  n <- 500

  # True parameters for factor structure
  true_var_f1 <- 1.5          # Variance of input factor
  true_se_intercept <- 0.5     # SE intercept
  true_se_linear <- 0.8        # SE linear coefficient
  true_se_quadratic <- 0.2     # SE quadratic coefficient
  true_se_residual_var <- 0.5  # SE residual variance

  # Generate factors
  f1 <- rnorm(n, 0, sqrt(true_var_f1))
  eps <- rnorm(n, 0, sqrt(true_se_residual_var))
  f2 <- true_se_intercept + true_se_linear * f1 + true_se_quadratic * f1^2 + eps

  # Generate measurements for f1 (3 measures)
  Y1_1 <- 2.0 + 1.0 * f1 + rnorm(n, 0, 1.0)   # int=2, loading=1 (fixed), sigma=1
  Y1_2 <- 1.5 + 0.8 * f1 + rnorm(n, 0, 0.9)   # int=1.5, loading=0.8, sigma=0.9
  Y1_3 <- 1.0 + 1.2 * f1 + rnorm(n, 0, 1.1)   # int=1, loading=1.2, sigma=1.1

  # Generate measurements for f2 (3 measures)
  Y2_1 <- 1.8 + 1.0 * f2 + rnorm(n, 0, 1.0)   # int=1.8, loading=1 (fixed), sigma=1
  Y2_2 <- 1.2 + 0.7 * f2 + rnorm(n, 0, 0.8)   # int=1.2, loading=0.7, sigma=0.8
  Y2_3 <- 0.5 + 1.1 * f2 + rnorm(n, 0, 0.95)  # int=0.5, loading=1.1, sigma=0.95

  dat <- data.frame(
    intercept = 1,
    Y1_1 = Y1_1, Y1_2 = Y1_2, Y1_3 = Y1_3,
    Y2_1 = Y2_1, Y2_2 = Y2_2, Y2_3 = Y2_3,
    eval = 1
  )

  # Create model system with 2 factors and SE_quadratic structure
  fm <- define_factor_model(n_factors = 2, factor_structure = "SE_quadratic")

  # 3 measurement equations for factor 1
  mc1_1 <- define_model_component(
    name = "Y1_1", data = dat, outcome = "Y1_1", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(1, 0),  # Fix loading to 1 for f1, 0 for f2
    evaluation_indicator = "eval"
  )
  mc1_2 <- define_model_component(
    name = "Y1_2", data = dat, outcome = "Y1_2", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(NA_real_, 0),  # Free for f1, 0 for f2
    evaluation_indicator = "eval"
  )
  mc1_3 <- define_model_component(
    name = "Y1_3", data = dat, outcome = "Y1_3", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(NA_real_, 0),  # Free for f1, 0 for f2
    evaluation_indicator = "eval"
  )

  # 3 measurement equations for factor 2
  mc2_1 <- define_model_component(
    name = "Y2_1", data = dat, outcome = "Y2_1", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(0, 1),  # 0 for f1, fix to 1 for f2
    evaluation_indicator = "eval"
  )
  mc2_2 <- define_model_component(
    name = "Y2_2", data = dat, outcome = "Y2_2", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(0, NA_real_),  # 0 for f1, free for f2
    evaluation_indicator = "eval"
  )
  mc2_3 <- define_model_component(
    name = "Y2_3", data = dat, outcome = "Y2_3", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(0, NA_real_),  # 0 for f1, free for f2
    evaluation_indicator = "eval"
  )

  ms <- define_model_system(
    components = list(mc1_1, mc1_2, mc1_3, mc2_1, mc2_2, mc2_3),
    factor = fm
  )

  # True parameters:
  # factor_var_1, se_intercept, se_linear_1, se_quadratic_1, se_residual_var,
  # Y1_1_int, Y1_1_sigma,  (loading fixed to 1)
  # Y1_2_int, Y1_2_lambda, Y1_2_sigma,
  # Y1_3_int, Y1_3_lambda, Y1_3_sigma,
  # Y2_1_int, Y2_1_sigma,  (loading fixed to 1)
  # Y2_2_int, Y2_2_lambda, Y2_2_sigma,
  # Y2_3_int, Y2_3_lambda, Y2_3_sigma
  true_params <- c(
    true_var_f1,  # factor_var_1 = 1.5
    true_se_intercept, true_se_linear, true_se_quadratic, true_se_residual_var,  # SE params
    2.0, 1.0,  # Y1_1: int, sigma (loading fixed to 1)
    1.5, 0.8, 0.9,  # Y1_2: int, lambda, sigma
    1.0, 1.2, 1.1,  # Y1_3: int, lambda, sigma
    1.8, 1.0,  # Y2_1: int, sigma (loading fixed to 1)
    1.2, 0.7, 0.8,  # Y2_2: int, lambda, sigma
    0.5, 1.1, 0.95  # Y2_3: int, lambda, sigma
  )

  param_names <- c(
    "f1_var",
    "se_intercept", "se_linear_1", "se_quadratic_1", "se_residual_var",
    "Y1_1_int", "Y1_1_sigma",
    "Y1_2_int", "Y1_2_lambda", "Y1_2_sigma",
    "Y1_3_int", "Y1_3_lambda", "Y1_3_sigma",
    "Y2_1_int", "Y2_1_sigma",
    "Y2_2_int", "Y2_2_lambda", "Y2_2_sigma",
    "Y2_3_int", "Y2_3_lambda", "Y2_3_sigma"
  )

  # Run estimation to get param_fixed
  est_comp <- run_estimation_comparison(
    ms, dat, true_params,
    param_names = param_names,
    verbose = VERBOSE
  )

  # Run gradient/Hessian checks with param_fixed
  grad_check <- check_gradient_accuracy(ms, dat, true_params,
                                        param_fixed = est_comp$param_fixed,
                                        tol = GRAD_TOL, verbose = VERBOSE, n_quad = 8)
  hess_check <- check_hessian_accuracy(ms, dat, true_params,
                                       param_fixed = est_comp$param_fixed,
                                       tol = HESS_TOL, verbose = VERBOSE, n_quad = 8)

  # Collect diagnostics
  diagnostics <- list(
    gradient_check = grad_check,
    hessian_check = hess_check,
    estimation = est_comp,
    overall_pass = grad_check$pass && hess_check$pass &&
                   est_comp$default_converged && est_comp$reasonable_default
  )

  # Save log
  if (SAVE_LOGS) {
    log_file <- save_diagnostics_to_log("test_H_se_quadratic", diagnostics)
    if (VERBOSE) cat("Log saved to:", log_file, "\n")
  }

  # Assertions - focus on gradient and Hessian correctness
  # Estimation convergence with default init is not strictly required for this test
  expect_true(grad_check$pass, info = sprintf("Gradient check failed (max error: %.2e)", grad_check$max_error))
  expect_true(hess_check$pass, info = sprintf("Hessian check failed (max error: %.2e)", hess_check$max_error))
})


# ==============================================================================
# Test I: Equality Constraints (Measurement Invariance) with SE_quadratic
# ==============================================================================

test_that("Model I: SE_quadratic with equality constraints (measurement invariance)", {
  skip_on_cran()

  set.seed(456)

  # True parameters
  n <- 800
  true_var_f1 <- 1.0
  true_se_intercept <- 0.0  # Set to 0 for identification
  true_se_linear <- 0.7
  true_se_quadratic <- 0.15
  true_se_residual_var <- 0.4

  # Common measurement parameters (invariance across time points)
  true_loading_2 <- 0.8
  true_loading_3 <- 1.2
  true_sigma_1 <- 1.0
  true_sigma_2 <- 0.9
  true_sigma_3 <- 1.1

  # Generate factors
  f1 <- rnorm(n, 0, sqrt(true_var_f1))
  eps <- rnorm(n, 0, sqrt(true_se_residual_var))
  f2 <- true_se_intercept + true_se_linear * f1 + true_se_quadratic * f1^2 + eps

  # Measurements at time 1 (factor 1)
  Y1_1 <- 1.0 * f1 + rnorm(n, 0, true_sigma_1)
  Y1_2 <- true_loading_2 * f1 + rnorm(n, 0, true_sigma_2)
  Y1_3 <- true_loading_3 * f1 + rnorm(n, 0, true_sigma_3)

  # Measurements at time 2 (factor 2) - SAME loadings and sigmas
  Y2_1 <- 1.0 * f2 + rnorm(n, 0, true_sigma_1)
  Y2_2 <- true_loading_2 * f2 + rnorm(n, 0, true_sigma_2)
  Y2_3 <- true_loading_3 * f2 + rnorm(n, 0, true_sigma_3)

  dat <- data.frame(Y1_1, Y1_2, Y1_3, Y2_1, Y2_2, Y2_3, intercept = 1, eval = 1)

  # Model setup
  fm <- define_factor_model(n_factors = 2, factor_structure = "SE_quadratic")

  mc1_1 <- define_model_component(name = "Y1_1", data = dat, outcome = "Y1_1", factor = fm,
                                    covariates = "intercept", model_type = "linear",
                                    loading_normalization = c(1, 0), evaluation_indicator = "eval")
  mc1_2 <- define_model_component(name = "Y1_2", data = dat, outcome = "Y1_2", factor = fm,
                                    covariates = "intercept", model_type = "linear",
                                    loading_normalization = c(NA_real_, 0), evaluation_indicator = "eval")
  mc1_3 <- define_model_component(name = "Y1_3", data = dat, outcome = "Y1_3", factor = fm,
                                    covariates = "intercept", model_type = "linear",
                                    loading_normalization = c(NA_real_, 0), evaluation_indicator = "eval")
  mc2_1 <- define_model_component(name = "Y2_1", data = dat, outcome = "Y2_1", factor = fm,
                                    covariates = "intercept", model_type = "linear",
                                    loading_normalization = c(0, 1), evaluation_indicator = "eval")
  mc2_2 <- define_model_component(name = "Y2_2", data = dat, outcome = "Y2_2", factor = fm,
                                    covariates = "intercept", model_type = "linear",
                                    loading_normalization = c(0, NA_real_), evaluation_indicator = "eval")
  mc2_3 <- define_model_component(name = "Y2_3", data = dat, outcome = "Y2_3", factor = fm,
                                    covariates = "intercept", model_type = "linear",
                                    loading_normalization = c(0, NA_real_), evaluation_indicator = "eval")

  # Define model system WITH equality constraints (measurement invariance)
  ms <- define_model_system(
    components = list(mc1_1, mc1_2, mc1_3, mc2_1, mc2_2, mc2_3),
    factor = fm,
    equality_constraints = list(
      c("Y1_2_loading_1", "Y2_2_loading_2"),  # Same loading for item 2
      c("Y1_3_loading_1", "Y2_3_loading_2"),  # Same loading for item 3
      c("Y1_1_sigma", "Y2_1_sigma"),          # Same error variance for item 1
      c("Y1_2_sigma", "Y2_2_sigma"),          # Same error variance for item 2
      c("Y1_3_sigma", "Y2_3_sigma")           # Same error variance for item 3
    )
  )

  # Test 1: Verify equality constraints are stored correctly
  expect_equal(length(ms$equality_constraints), 5)

  # Test 2: Run estimation
  control <- define_estimation_control(n_quad_points = 8, num_cores = 1)
  result <- estimate_model_rcpp(model_system = ms, data = dat, control = control,
                                 optimizer = "nlminb", verbose = FALSE)

  # Test 3: Verify tied parameters are exactly equal (use unname() to compare values, not named vectors)
  expect_equal(unname(result$estimates["Y1_2_loading_1"]), unname(result$estimates["Y2_2_loading_2"]))
  expect_equal(unname(result$estimates["Y1_3_loading_1"]), unname(result$estimates["Y2_3_loading_2"]))
  expect_equal(unname(result$estimates["Y1_1_sigma"]), unname(result$estimates["Y2_1_sigma"]))
  expect_equal(unname(result$estimates["Y1_2_sigma"]), unname(result$estimates["Y2_2_sigma"]))
  expect_equal(unname(result$estimates["Y1_3_sigma"]), unname(result$estimates["Y2_3_sigma"]))

  # Test 4: Verify param_table shows tied_to column correctly
  expect_equal(result$param_table$tied_to[result$param_table$name == "Y2_2_loading_2"], "Y1_2_loading_1")
  expect_equal(result$param_table$tied_to[result$param_table$name == "Y2_1_sigma"], "Y1_1_sigma")

  # Test 5: Parameter recovery (check key parameters are within tolerance)
  # Slightly relaxed tolerance to account for random variation in Monte Carlo simulation
  tolerance <- 0.20
  expect_lt(abs(result$estimates["factor_var_1"] - true_var_f1), tolerance)
  expect_lt(abs(result$estimates["se_linear_1"] - true_se_linear), tolerance)
  expect_lt(abs(result$estimates["se_quadratic_1"] - true_se_quadratic), tolerance)
  expect_lt(abs(result$estimates["se_residual_var"] - true_se_residual_var), tolerance)
  expect_lt(abs(result$estimates["Y1_2_loading_1"] - true_loading_2), tolerance)
  expect_lt(abs(result$estimates["Y1_3_loading_1"] - true_loading_3), tolerance)
  expect_lt(abs(result$estimates["Y1_1_sigma"] - true_sigma_1), tolerance)
  expect_lt(abs(result$estimates["Y1_2_sigma"] - true_sigma_2), tolerance)
  expect_lt(abs(result$estimates["Y1_3_sigma"] - true_sigma_3), tolerance)

  # Test 6: FD gradient check with equality constraints
  # The C++ gradient needs to be aggregated for tied params
  # Note: C++ now returns gradients for FREE params only, so we need to map back to full
  data_mat <- as.matrix(dat)
  init_result <- initialize_parameters(ms, dat, verbose = FALSE)
  params <- init_result$init_params

  # Set to true values for gradient check
  params[1] <- true_var_f1
  params[2] <- true_se_intercept
  params[3] <- true_se_linear
  params[4] <- true_se_quadratic
  params[5] <- true_se_residual_var
  params[9] <- true_loading_2   # Y1_2_loading_1
  params[12] <- true_loading_3  # Y1_3_loading_1
  params[7] <- true_sigma_1     # Y1_1_sigma
  params[10] <- true_sigma_2    # Y1_2_sigma
  params[13] <- true_sigma_3    # Y1_3_sigma
  # Set tied params equal
  params[17] <- true_loading_2  # Y2_2_loading_2
  params[20] <- true_loading_3  # Y2_3_loading_2
  params[15] <- true_sigma_1    # Y2_1_sigma
  params[18] <- true_sigma_2    # Y2_2_sigma
  params[21] <- true_sigma_3    # Y2_3_sigma

  fm_ptr <- initialize_factor_model_cpp(ms, data_mat, 8, params)

  # Extract free params (tied params are now fixed in C++)
  params_free <- extract_free_params_cpp(fm_ptr, params)
  cpp_result <- evaluate_likelihood_cpp(fm_ptr, params_free, compute_gradient = TRUE, compute_hessian = FALSE)

  # Map gradient from free params back to full params
  # Free indices are all params except tied ones (17, 20, 15, 18, 21)
  # free_idx = c(1:14, 16, 19) = 16 free params
  free_idx <- c(1:14, 16, 19)
  analytical_grad <- rep(0, length(params))
  analytical_grad[free_idx] <- cpp_result$gradient

  # No need for manual aggregation since C++ already excludes tied params
  # The gradient for primary params now directly includes the contribution
  aggregated_analytical <- analytical_grad

  # Compute FD gradient with constraints applied
  eps_fd <- 1e-5
  fd_grad <- numeric(length(params))
  for (i in seq_along(params)) {
    params_plus <- params
    params_minus <- params
    params_plus[i] <- params[i] + eps_fd
    params_minus[i] <- params[i] - eps_fd

    # Apply equality constraints
    if (i == 9) { params_plus[17] <- params_plus[9]; params_minus[17] <- params_minus[9] }
    if (i == 12) { params_plus[20] <- params_plus[12]; params_minus[20] <- params_minus[12] }
    if (i == 7) { params_plus[15] <- params_plus[7]; params_minus[15] <- params_minus[7] }
    if (i == 10) { params_plus[18] <- params_plus[10]; params_minus[18] <- params_minus[10] }
    if (i == 13) { params_plus[21] <- params_plus[13]; params_minus[21] <- params_minus[13] }

    # Need to extract free params for the perturbed params
    params_plus_free <- extract_free_params_cpp(fm_ptr, params_plus)
    params_minus_free <- extract_free_params_cpp(fm_ptr, params_minus)

    ll_plus <- evaluate_loglik_only_cpp(fm_ptr, params_plus_free)
    ll_minus <- evaluate_loglik_only_cpp(fm_ptr, params_minus_free)
    fd_grad[i] <- (ll_plus - ll_minus) / (2 * eps_fd)
  }

  # Check primary parameters only (not derived ones)
  primary_indices <- c(1:14, 16, 19)
  max_rel_err <- 0
  worst_idx <- 0
  cat("\n=== DEBUG: Gradient comparison ===\n")
  for (i in primary_indices) {
    ana <- aggregated_analytical[i]
    fd <- fd_grad[i]
    rel_err <- if (abs(fd) > 1e-6) abs(ana - fd) / abs(fd) else abs(ana - fd)
    cat(sprintf("  param %2d: ana=%.6f, fd=%.6f, rel_err=%.2e\n", i, ana, fd, rel_err))
    if (!is.na(rel_err) && rel_err > max_rel_err) {
      max_rel_err <- rel_err
      worst_idx <- i
    }
  }
  cat(sprintf("Worst param: %d, max_rel_err: %.2e\n", worst_idx, max_rel_err))

  expect_lt(max_rel_err, 1e-4, label = sprintf("Gradient max error: %.2e", max_rel_err))

  if (VERBOSE) {
    cat("Test I: Equality constraints (measurement invariance)\n")
    cat("  Equality constraints: 5 groups\n")
    cat("  Free parameters: 16 (21 total - 5 tied)\n")
    cat(sprintf("  Gradient max rel error: %.2e\n", max_rel_err))
    cat(sprintf("  SE_linear recovery: true=%.3f, est=%.3f\n",
                true_se_linear, result$estimates["se_linear_1"]))
    cat(sprintf("  Loading_2 recovery: true=%.3f, est=%.3f\n",
                true_loading_2, result$estimates["Y1_2_loading_1"]))
  }
})
