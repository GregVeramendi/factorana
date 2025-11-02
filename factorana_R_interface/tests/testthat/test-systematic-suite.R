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
GRAD_TOL <- 1e-5
HESS_TOL <- 1e-4

# ==============================================================================
# Test A: Linear Model with Zero Factors
# ==============================================================================

test_that("Model A: Linear model with zero factors", {
  skip_on_cran()
  skip("Zero-factor models not yet supported in C++ - requires handling n_factors=0 throughout codebase")

  set.seed(100)

  # Simulate data
  n <- 500
  x1 <- rnorm(n)
  x2 <- rnorm(n)

  # True parameters: intercept, beta1, beta2, sigma
  true_params <- c(2.0, 1.5, -0.8, 0.5)

  y <- true_params[1] + true_params[2]*x1 + true_params[3]*x2 +
       rnorm(n, 0, true_params[4])

  dat <- data.frame(intercept = 1, x1 = x1, x2 = x2, y = y, eval = 1)

  # Create model system with 0 factors
  fm <- define_factor_model(n_factors = 0, n_types = 1, n_quad = 8)
  mc <- define_model_component(
    name = "y",
    data = dat,
    outcome = "y",
    factor = fm,
    covariates = c("intercept", "x1", "x2"),
    model_type = "linear",
    loading_normalization = numeric(0),
    evaluation_indicator = "eval"
  )
  ms <- define_model_system(components = list(mc), factor = fm)

  # Run checks
  grad_check <- check_gradient_accuracy(ms, dat, true_params, tol = GRAD_TOL, verbose = VERBOSE)
  hess_check <- check_hessian_accuracy(ms, dat, true_params, tol = HESS_TOL, verbose = VERBOSE)
  est_comp <- run_estimation_comparison(
    ms, dat, true_params,
    param_names = c("intercept", "beta1", "beta2", "sigma"),
    verbose = VERBOSE
  )

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
    log_file <- save_diagnostics_to_log("test_A_linear_zero_factors", diagnostics)
    if (VERBOSE) cat("Log saved to:", log_file, "\n")
  }

  # Assertions
  expect_true(grad_check$pass, info = sprintf("Gradient check failed (max error: %.2e)", grad_check$max_error))
  expect_true(hess_check$pass, info = sprintf("Hessian check failed (max error: %.2e)", hess_check$max_error))
  expect_true(est_comp$default_converged, info = "Estimation with default init failed to converge")
  expect_true(est_comp$reasonable_default, info = "Estimates are not reasonable")
})

# ==============================================================================
# Test B: Probit Model with Zero Factors
# ==============================================================================

test_that("Model B: Probit model with zero factors", {
  skip_on_cran()
  skip("Zero-factor models not yet supported in C++ - requires handling n_factors=0 throughout codebase")

  set.seed(101)

  # Simulate data
  n <- 500
  x1 <- rnorm(n)
  x2 <- rnorm(n)

  # True parameters: intercept, beta1, beta2
  true_params <- c(0.5, 1.0, -0.7)

  z <- true_params[1] + true_params[2]*x1 + true_params[3]*x2
  y <- as.numeric(runif(n) < pnorm(z))

  dat <- data.frame(intercept = 1, x1 = x1, x2 = x2, y = y, eval = 1)

  # Create model system with 0 factors
  fm <- define_factor_model(n_factors = 0, n_types = 1, n_quad = 8)
  mc <- define_model_component(
    name = "y",
    data = dat,
    outcome = "y",
    factor = fm,
    covariates = c("intercept", "x1", "x2"),
    model_type = "probit",
    loading_normalization = numeric(0),
    evaluation_indicator = "eval"
  )
  ms <- define_model_system(components = list(mc), factor = fm)

  # Run checks
  grad_check <- check_gradient_accuracy(ms, dat, true_params, tol = GRAD_TOL, verbose = VERBOSE)
  hess_check <- check_hessian_accuracy(ms, dat, true_params, tol = HESS_TOL, verbose = VERBOSE)
  est_comp <- run_estimation_comparison(
    ms, dat, true_params,
    param_names = c("intercept", "beta1", "beta2"),
    verbose = VERBOSE
  )

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
    log_file <- save_diagnostics_to_log("test_B_probit_zero_factors", diagnostics)
    if (VERBOSE) cat("Log saved to:", log_file, "\n")
  }

  # Assertions
  expect_true(grad_check$pass, info = sprintf("Gradient check failed (max error: %.2e)", grad_check$max_error))
  expect_true(hess_check$pass, info = sprintf("Hessian check failed (max error: %.2e)", hess_check$max_error))
  expect_true(est_comp$default_converged, info = "Estimation with default init failed to converge")
  expect_true(est_comp$reasonable_default, info = "Estimates are not reasonable")
})

# ==============================================================================
# Test C: Ordered Probit with 3 Choices and Zero Factors
# ==============================================================================

test_that("Model C: Ordered probit with 3 choices and zero factors", {
  skip_on_cran()
  skip("Zero-factor models not yet supported in C++ - requires handling n_factors=0 throughout codebase")

  set.seed(102)

  # Simulate data
  n <- 500
  x1 <- rnorm(n)
  x2 <- rnorm(n)

  # True parameters: intercept, beta1, beta2, tau1, tau2
  # (3 choices means 2 thresholds)
  true_params <- c(0.5, 0.8, -0.6, 0.0, 1.0)

  z <- true_params[1] + true_params[2]*x1 + true_params[3]*x2
  tau1 <- true_params[4]
  tau2 <- tau1 + abs(true_params[5])  # Ensure tau2 > tau1

  # Generate ordered outcome (1, 2, or 3)
  u <- rnorm(n)
  y <- ifelse(z + u < tau1, 1,
              ifelse(z + u < tau2, 2, 3))

  dat <- data.frame(intercept = 1, x1 = x1, x2 = x2, y = y, eval = 1)

  # Create model system with 0 factors
  fm <- define_factor_model(n_factors = 0, n_types = 1, n_quad = 8)
  mc <- define_model_component(
    name = "y",
    data = dat,
    outcome = "y",
    factor = fm,
    covariates = c("intercept", "x1", "x2"),
    model_type = "oprobit",
    loading_normalization = numeric(0),
    evaluation_indicator = "eval"
  )
  ms <- define_model_system(components = list(mc), factor = fm)

  # Run checks
  grad_check <- check_gradient_accuracy(ms, dat, true_params, tol = GRAD_TOL, verbose = VERBOSE)
  hess_check <- check_hessian_accuracy(ms, dat, true_params, tol = HESS_TOL, verbose = VERBOSE)
  est_comp <- run_estimation_comparison(
    ms, dat, true_params,
    param_names = c("intercept", "beta1", "beta2", "tau1", "tau2"),
    verbose = VERBOSE
  )

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
    log_file <- save_diagnostics_to_log("test_C_oprobit_zero_factors", diagnostics)
    if (VERBOSE) cat("Log saved to:", log_file, "\n")
  }

  # Assertions
  expect_true(grad_check$pass, info = sprintf("Gradient check failed (max error: %.2e)", grad_check$max_error))
  expect_true(hess_check$pass, info = sprintf("Hessian check failed (max error: %.2e)", hess_check$max_error))
  expect_true(est_comp$default_converged, info = "Estimation with default init failed to converge")
  expect_true(est_comp$reasonable_default, info = "Estimates are not reasonable")
})

# ==============================================================================
# Test D: Multinomial Logit with 3 Choices and Zero Factors
# ==============================================================================

test_that("Model D: Multinomial logit with 3 choices and zero factors", {
  skip_on_cran()
  skip("Zero-factor models not yet supported in C++ - requires handling n_factors=0 throughout codebase")

  set.seed(103)

  # Simulate data
  n <- 500
  x1 <- rnorm(n)
  x2 <- rnorm(n)

  # True parameters: intercept1, beta1_1, beta2_1, intercept2, beta1_2, beta2_2
  # (3 choices with choice 0 as reference means 2 sets of parameters)
  true_params <- c(0.5, 0.8, -0.4, 1.0, -0.6, 0.7)

  # Linear predictors for choices 1 and 2 (choice 0 is reference)
  z1 <- true_params[1] + true_params[2]*x1 + true_params[3]*x2
  z2 <- true_params[4] + true_params[5]*x1 + true_params[6]*x2

  # Multinomial logit probabilities
  exp_z0 <- 1
  exp_z1 <- exp(z1)
  exp_z2 <- exp(z2)
  denom <- exp_z0 + exp_z1 + exp_z2

  p0 <- exp_z0 / denom
  p1 <- exp_z1 / denom
  p2 <- exp_z2 / denom

  # Generate choice outcome
  y <- numeric(n)
  for (i in seq_len(n)) {
    y[i] <- sample(0:2, 1, prob = c(p0[i], p1[i], p2[i]))
  }

  dat <- data.frame(intercept = 1, x1 = x1, x2 = x2, y = y, eval = 1)

  # Create model system with 0 factors
  fm <- define_factor_model(n_factors = 0, n_types = 1, n_quad = 8)
  mc <- define_model_component(
    name = "y",
    data = dat,
    outcome = "y",
    factor = fm,
    covariates = c("intercept", "x1", "x2"),
    model_type = "mlogit",
    loading_normalization = numeric(0),
    evaluation_indicator = "eval"
  )
  ms <- define_model_system(components = list(mc), factor = fm)

  # Run checks
  grad_check <- check_gradient_accuracy(ms, dat, true_params, tol = GRAD_TOL, verbose = VERBOSE)
  hess_check <- check_hessian_accuracy(ms, dat, true_params, tol = HESS_TOL, verbose = VERBOSE)
  est_comp <- run_estimation_comparison(
    ms, dat, true_params,
    param_names = c("int1", "beta1_1", "beta2_1", "int2", "beta1_2", "beta2_2"),
    verbose = VERBOSE
  )

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
    log_file <- save_diagnostics_to_log("test_D_mlogit_zero_factors", diagnostics)
    if (VERBOSE) cat("Log saved to:", log_file, "\n")
  }

  # Assertions
  expect_true(grad_check$pass, info = sprintf("Gradient check failed (max error: %.2e)", grad_check$max_error))
  expect_true(hess_check$pass, info = sprintf("Hessian check failed (max error: %.2e)", hess_check$max_error))
  expect_true(est_comp$default_converged, info = "Estimation with default init failed to converge")
  expect_true(est_comp$reasonable_default, info = "Estimates are not reasonable")
})

# ==============================================================================
# Test E: Measurement System with 3 Linear Tests and 1 Factor
# ==============================================================================

test_that("Model E: Measurement system with 3 linear tests and 1 factor", {
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
  fm <- define_factor_model(n_factors = 1, n_types = 1, n_quad = 8)

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

  # Run checks
  grad_check <- check_gradient_accuracy(ms, dat, true_params, tol = GRAD_TOL, verbose = VERBOSE)
  hess_check <- check_hessian_accuracy(ms, dat, true_params, tol = HESS_TOL, verbose = VERBOSE)
  est_comp <- run_estimation_comparison(
    ms, dat, true_params,
    param_names = c("f_var", "T1_int", "T1_sigma",  # T1 loading fixed to 1.0
                   "T2_int", "T2_lambda", "T2_sigma",
                   "T3_int", "T3_lambda", "T3_sigma"),
    verbose = VERBOSE
  )

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
    log_file <- save_diagnostics_to_log("test_E_measurement_system", diagnostics)
    if (VERBOSE) cat("Log saved to:", log_file, "\n")
  }

  # Assertions
  expect_true(grad_check$pass, info = sprintf("Gradient check failed (max error: %.2e)", grad_check$max_error))
  expect_true(hess_check$pass, info = sprintf("Hessian check failed (max error: %.2e)", hess_check$max_error))
  expect_true(est_comp$default_converged, info = "Estimation with default init failed to converge")
  expect_true(est_comp$reasonable_default, info = "Estimates are not reasonable")
})

# ==============================================================================
# Test F: Measurement System (E) + Probit (B)
# ==============================================================================

test_that("Model F: Measurement system with 3 linear tests and probit outcome", {
  skip_on_cran()

  set.seed(105)

  # Simulate data
  n <- 500
  x1 <- rnorm(n)
  f <- rnorm(n)  # Latent factor

  # True parameters: factor_var, T1-T3 params (9), probit params (4)
  true_params <- c(1.0,  # Factor variance (fixed)
                   2.0, 1.0, 0.5,   # T1: int, loading (fixed to 1), sigma
                   1.5, 1.2, 0.6,   # T2: int, loading, sigma
                   1.0, 0.8, 0.4,   # T3: int, loading, sigma
                   0.5, 0.7, 0.9)   # Probit: int, beta1, loading

  T1 <- true_params[2] + true_params[3]*f + rnorm(n, 0, true_params[4])
  T2 <- true_params[5] + true_params[6]*f + rnorm(n, 0, true_params[7])
  T3 <- true_params[8] + true_params[9]*f + rnorm(n, 0, true_params[10])

  z <- true_params[11] + true_params[12]*x1 + true_params[13]*f
  y <- as.numeric(runif(n) < pnorm(z))

  dat <- data.frame(intercept = 1, x1 = x1, T1 = T1, T2 = T2, T3 = T3, y = y, eval = 1)

  # Create model system
  fm <- define_factor_model(n_factors = 1, n_types = 1, n_quad = 8)

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

  # Run checks
  grad_check <- check_gradient_accuracy(ms, dat, true_params, tol = GRAD_TOL, verbose = VERBOSE)
  hess_check <- check_hessian_accuracy(ms, dat, true_params, tol = HESS_TOL, verbose = VERBOSE)
  est_comp <- run_estimation_comparison(
    ms, dat, true_params,
    param_names = c("f_var",
                   "T1_int", "T1_sig",  # T1: loading fixed to 1.0
                   "T2_int", "T2_sig", "T2_load",
                   "T3_int", "T3_sig", "T3_load",
                   "y_int", "y_beta1", "y_load"),
    verbose = VERBOSE
  )

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
    log_file <- save_diagnostics_to_log("test_F_measurement_plus_probit", diagnostics)
    if (VERBOSE) cat("Log saved to:", log_file, "\n")
  }

  # Assertions
  expect_true(grad_check$pass, info = sprintf("Gradient check failed (max error: %.2e)", grad_check$max_error))
  expect_true(hess_check$pass, info = sprintf("Hessian check failed (max error: %.2e)", hess_check$max_error))
  expect_true(est_comp$default_converged, info = "Estimation with default init failed to converge")
  expect_true(est_comp$reasonable_default, info = "Estimates are not reasonable")
})

# ==============================================================================
# Test G: Measurement System (E) + Ordered Probit (C)
# ==============================================================================

test_that("Model G: Measurement system with 3 linear tests and ordered probit", {
  skip_on_cran()

  set.seed(106)

  # Simulate data
  n <- 500
  x1 <- rnorm(n)
  f <- rnorm(n)  # Latent factor

  # True parameters: factor_var, T1-T3 params (9), oprobit params (5)
  true_params <- c(1.0,  # Factor variance (fixed)
                   2.0, 1.0, 0.5,   # T1: int, loading (fixed to 1), sigma
                   1.5, 1.2, 0.6,   # T2: int, loading, sigma
                   1.0, 0.8, 0.4,   # T3: int, loading, sigma
                   0.5, 0.6, 0.8, 0.0, 1.0)   # Oprobit: int, beta1, loading, tau1, tau2

  T1 <- true_params[2] + true_params[3]*f + rnorm(n, 0, true_params[4])
  T2 <- true_params[5] + true_params[6]*f + rnorm(n, 0, true_params[7])
  T3 <- true_params[8] + true_params[9]*f + rnorm(n, 0, true_params[10])

  z <- true_params[11] + true_params[12]*x1 + true_params[13]*f
  tau1 <- true_params[14]
  tau2 <- tau1 + abs(true_params[15])

  u <- rnorm(n)
  y <- ifelse(z + u < tau1, 1, ifelse(z + u < tau2, 2, 3))

  dat <- data.frame(intercept = 1, x1 = x1, T1 = T1, T2 = T2, T3 = T3, y = y, eval = 1)

  # Create model system
  fm <- define_factor_model(n_factors = 1, n_types = 1, n_quad = 8)

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
    covariates = c("intercept", "x1"), model_type = "oprobit",
    loading_normalization = NA_real_, evaluation_indicator = "eval"
  )

  ms <- define_model_system(components = list(mc_T1, mc_T2, mc_T3, mc_y), factor = fm)

  # Run checks
  grad_check <- check_gradient_accuracy(ms, dat, true_params, tol = GRAD_TOL, verbose = VERBOSE)
  hess_check <- check_hessian_accuracy(ms, dat, true_params, tol = HESS_TOL, verbose = VERBOSE)
  est_comp <- run_estimation_comparison(
    ms, dat, true_params,
    param_names = c("f_var",
                   "T1_int", "T1_sig",  # T1: loading fixed to 1.0
                   "T2_int", "T2_sig", "T2_load",
                   "T3_int", "T3_sig", "T3_load",
                   "y_int", "y_beta1", "y_load", "tau2"),  # tau1 might be fixed to 0
    verbose = VERBOSE
  )

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
    log_file <- save_diagnostics_to_log("test_G_measurement_plus_oprobit", diagnostics)
    if (VERBOSE) cat("Log saved to:", log_file, "\n")
  }

  # Assertions
  expect_true(grad_check$pass, info = sprintf("Gradient check failed (max error: %.2e)", grad_check$max_error))
  expect_true(hess_check$pass, info = sprintf("Hessian check failed (max error: %.2e)", hess_check$max_error))
  expect_true(est_comp$default_converged, info = "Estimation with default init failed to converge")
  expect_true(est_comp$reasonable_default, info = "Estimates are not reasonable")
})

# ==============================================================================
# Test H: Measurement System (E) + Multinomial Logit (D)
# ==============================================================================

test_that("Model H: Measurement system with 3 linear tests and multinomial logit", {
  skip_on_cran()

  set.seed(107)

  # Simulate data
  n <- 500
  x1 <- rnorm(n)
  f <- rnorm(n)  # Latent factor

  # True parameters: factor_var, T1-T3 params (9), mlogit params (6)
  true_params <- c(1.0,  # Factor variance (fixed)
                   2.0, 1.0, 0.5,   # T1: int, loading (fixed to 1), sigma
                   1.5, 1.2, 0.6,   # T2: int, loading, sigma
                   1.0, 0.8, 0.4,   # T3: int, loading, sigma
                   0.5, 0.6, 0.7,   # Mlogit choice 1: int, beta1, loading
                   1.0, -0.5, 0.9)  # Mlogit choice 2: int, beta1, loading

  T1 <- true_params[2] + true_params[3]*f + rnorm(n, 0, true_params[4])
  T2 <- true_params[5] + true_params[6]*f + rnorm(n, 0, true_params[7])
  T3 <- true_params[8] + true_params[9]*f + rnorm(n, 0, true_params[10])

  z1 <- true_params[11] + true_params[12]*x1 + true_params[13]*f
  z2 <- true_params[14] + true_params[15]*x1 + true_params[16]*f

  exp_z0 <- 1
  exp_z1 <- exp(z1)
  exp_z2 <- exp(z2)
  denom <- exp_z0 + exp_z1 + exp_z2

  p0 <- exp_z0 / denom
  p1 <- exp_z1 / denom
  p2 <- exp_z2 / denom

  y <- numeric(n)
  for (i in seq_len(n)) {
    y[i] <- sample(0:2, 1, prob = c(p0[i], p1[i], p2[i]))
  }

  dat <- data.frame(intercept = 1, x1 = x1, T1 = T1, T2 = T2, T3 = T3, y = y, eval = 1)

  # Create model system
  fm <- define_factor_model(n_factors = 1, n_types = 1, n_quad = 8)

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

  # Run checks
  grad_check <- check_gradient_accuracy(ms, dat, true_params, tol = GRAD_TOL, verbose = VERBOSE)
  hess_check <- check_hessian_accuracy(ms, dat, true_params, tol = HESS_TOL, verbose = VERBOSE)
  est_comp <- run_estimation_comparison(
    ms, dat, true_params,
    param_names = c("f_var", "T1_int", "T1_load", "T1_sig",
                   "T2_int", "T2_load", "T2_sig",
                   "T3_int", "T3_load", "T3_sig",
                   "y1_int", "y1_beta1", "y1_load",
                   "y2_int", "y2_beta1", "y2_load"),
    verbose = VERBOSE
  )

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
    log_file <- save_diagnostics_to_log("test_H_measurement_plus_mlogit", diagnostics)
    if (VERBOSE) cat("Log saved to:", log_file, "\n")
  }

  # Assertions
  expect_true(grad_check$pass, info = sprintf("Gradient check failed (max error: %.2e)", grad_check$max_error))
  expect_true(hess_check$pass, info = sprintf("Hessian check failed (max error: %.2e)", hess_check$max_error))
  expect_true(est_comp$default_converged, info = "Estimation with default init failed to converge")
  expect_true(est_comp$reasonable_default, info = "Estimates are not reasonable")
})

# ==============================================================================
# Test I: Roy Model
# ==============================================================================

test_that("Model I: Roy selection model", {
  skip_on_cran()

  set.seed(108)

  # Simulate Roy model data
  n <- 500
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  f <- rnorm(n)  # Latent factor (ability)

  # True parameters (based on our previous Roy model tests)
  # factor_var, T1-T3 params (9), wage0 params (4), wage1 params (4), sector params (3)
  true_params <- c(
    1.0,  # Factor variance (fixed)
    # T1: Test 1
    0.0, 1.0, 0.5,   # int, loading (fixed to 1), sigma
    # T2: Test 2
    0.0, 1.2, 0.6,   # int, loading, sigma
    # T3: Test 3
    0.0, 0.8, 0.4,   # int, loading, sigma
    # Wage0 (sector 0)
    2.0, 0.5, 0.3, 0.6,   # int, beta1, beta2, sigma
    # Wage1 (sector 1)
    2.5, 0.6, 1.0, 0.7,   # int, beta1, loading, sigma
    # Sector choice (probit)
    0.0, 0.4, 0.8   # int, beta1, loading
  )

  # Generate test scores
  T1 <- true_params[2] + true_params[3]*f + rnorm(n, 0, true_params[4])
  T2 <- true_params[5] + true_params[6]*f + rnorm(n, 0, true_params[7])
  T3 <- true_params[8] + true_params[9]*f + rnorm(n, 0, true_params[10])

  # Generate potential wages
  wage0 <- true_params[11] + true_params[12]*x1 + true_params[13]*x2 +
           rnorm(n, 0, true_params[14])
  wage1 <- true_params[15] + true_params[16]*x1 + true_params[17]*f +
           rnorm(n, 0, true_params[18])

  # Sector choice (based on utility difference)
  z_sector <- true_params[19] + true_params[20]*x2 + true_params[21]*f
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
  fm <- define_factor_model(n_factors = 1, n_types = 1, n_quad = 16)

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

  # Run checks
  grad_check <- check_gradient_accuracy(ms, dat, true_params, tol = GRAD_TOL, verbose = VERBOSE)
  hess_check <- check_hessian_accuracy(ms, dat, true_params, tol = HESS_TOL, verbose = VERBOSE)
  est_comp <- run_estimation_comparison(
    ms, dat, true_params,
    param_names = c("f_var",
                   "T1_int", "T1_sig",  # T1: loading fixed to 1.0
                   "T2_int", "T2_sig", "T2_load",
                   "T3_int", "T3_sig", "T3_load",
                   "w0_int", "w0_beta1", "w0_beta2", "w0_sig",
                   "w1_int", "w1_beta1", "w1_load", "w1_sig",
                   "s_int", "s_beta1", "s_load"),
    verbose = VERBOSE
  )

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
    log_file <- save_diagnostics_to_log("test_I_roy_model", diagnostics)
    if (VERBOSE) cat("Log saved to:", log_file, "\n")
  }

  # Assertions
  expect_true(grad_check$pass, info = sprintf("Gradient check failed (max error: %.2e)", grad_check$max_error))
  expect_true(hess_check$pass, info = sprintf("Hessian check failed (max error: %.2e)", hess_check$max_error))
  expect_true(est_comp$default_converged, info = "Estimation with default init failed to converge")
  expect_true(est_comp$reasonable_default, info = "Estimates are not reasonable")
})
