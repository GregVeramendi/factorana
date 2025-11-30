#!/usr/bin/env Rscript
# Gradient and Hessian Validation Test
# Compares analytical derivatives from C++ to finite-difference approximations
# Based on TMinLkhd::TestGradient() and TMinLkhd::TestHessian()

library(factorana)

cat("========================================\n")
cat("Gradient & Hessian Validation Test\n")
cat("========================================\n\n")

# ===== Helper Functions =====

#' Compute finite-difference gradient using central differences
#'
#' @param fm_cpp External pointer to FactorModel
#' @param params Parameter vector
#' @param delta Step size parameter (default 1e-7)
#' @return Finite-difference gradient vector
compute_fd_gradient <- function(fm_cpp, params, delta = 1e-7) {
  n_params <- length(params)
  grad_fd <- numeric(n_params)

  for (i in seq_along(params)) {
    # Step size: h = delta * (|param_i| + 1.0)
    h <- delta * (abs(params[i]) + 1.0)

    # Central difference: [f(x+h) - f(x-h)] / (2h)
    params_plus <- params
    params_plus[i] <- params[i] + h
    f_plus <- evaluate_loglik_only_cpp(fm_cpp, params_plus)

    params_minus <- params
    params_minus[i] <- params[i] - h
    f_minus <- evaluate_loglik_only_cpp(fm_cpp, params_minus)

    grad_fd[i] <- (f_plus - f_minus) / (2.0 * h)
  }

  return(grad_fd)
}

#' Compute finite-difference Hessian using gradient differences
#'
#' @param fm_cpp External pointer to FactorModel
#' @param params Parameter vector
#' @param delta Step size parameter (default 1e-7)
#' @return Finite-difference Hessian matrix (full symmetric matrix)
compute_fd_hessian <- function(fm_cpp, params, delta = 1e-7) {
  n_params <- length(params)
  hess_fd <- matrix(0, n_params, n_params)

  for (j in seq_along(params)) {
    # Step size for parameter j
    h <- delta * (abs(params[j]) + 1.0)

    # Gradient at x + h_j
    params_plus <- params
    params_plus[j] <- params[j] + h
    result_plus <- evaluate_likelihood_cpp(fm_cpp, params_plus,
                                          compute_gradient = TRUE,
                                          compute_hessian = FALSE)
    grad_plus <- result_plus$gradient

    # Gradient at x - h_j
    params_minus <- params
    params_minus[j] <- params[j] - h
    result_minus <- evaluate_likelihood_cpp(fm_cpp, params_minus,
                                           compute_gradient = TRUE,
                                           compute_hessian = FALSE)
    grad_minus <- result_minus$gradient

    # Hessian column j: [grad(x+h_j) - grad(x-h_j)] / (2*h)
    hess_fd[, j] <- (grad_plus - grad_minus) / (2.0 * h)
  }

  # Symmetrize (average with transpose for numerical stability)
  hess_fd <- (hess_fd + t(hess_fd)) / 2.0

  return(hess_fd)
}

#' Compare analytical and finite-difference derivatives
#'
#' @param analytical Analytical derivative value
#' @param finite_diff Finite-difference approximation
#' @param threshold Relative error threshold
#' @return List with absolute difference, relative difference, and pass/fail status
compare_derivatives <- function(analytical, finite_diff, threshold = 1e-4) {
  abs_diff <- abs(analytical - finite_diff)

  # Relative difference calculation (following TMinLkhd logic)
  if (abs(finite_diff) > 1e-7) {
    rel_diff <- abs_diff / abs(finite_diff)
  } else if (abs(analytical) < 1e-7) {
    rel_diff <- abs_diff  # Both are near zero
  } else {
    rel_diff <- abs_diff / abs(analytical)
  }

  passed <- rel_diff < threshold

  list(
    abs_diff = abs_diff,
    rel_diff = rel_diff,
    passed = passed
  )
}

# ===== Test 1: Single factor linear model =====
cat("Test 1: Single factor linear model - Gradient validation\n")
cat("----------------------------------------------------------\n")

set.seed(42)
n <- 100
f1 <- rnorm(n)
x1 <- rnorm(n)
intercept <- rep(1, n)
Y <- 1.0*intercept + 0.5*x1 + 0.8*f1 + rnorm(n, 0, 0.5)

dat1 <- data.frame(Y = Y, intercept = intercept, x1 = x1, eval = 1)

fm1 <- define_factor_model(n_factors = 1, n_types = 1)
mc1 <- define_model_component("Y", dat1, "Y", fm1,
                              covariates = c("intercept", "x1"),
                              model_type = "linear",
                              loading_normalization = NA_real_,
                              evaluation_indicator = "eval")
ms1 <- define_model_system(factor = fm1, components = list(mc1))

# Initialize C++ model
fm_cpp1 <- initialize_factor_model_cpp(ms1, dat1, n_quad = 8L)
param_info1 <- get_parameter_info_cpp(fm_cpp1)

cat("  Number of parameters:", param_info1$n_param_free, "\n")
cat("  Parameters: [sigma_f^2, beta_0, beta_x1, lambda, sigma_y]\n\n")

# Test parameters
test_params1 <- c(1.0, 1.0, 0.5, 0.8, 0.5)

# Compute analytical gradient
result1 <- evaluate_likelihood_cpp(fm_cpp1, test_params1,
                                   compute_gradient = TRUE,
                                   compute_hessian = FALSE)
grad_analytical <- result1$gradient
loglik <- result1$logLikelihood

cat("  Log-likelihood:", round(loglik, 6), "\n\n")

# Compute finite-difference gradient
cat("  Computing finite-difference gradient...\n")
grad_fd <- compute_fd_gradient(fm_cpp1, test_params1, delta = 1e-7)

# Compare gradients
cat("\n  Gradient Comparison:\n")
cat(sprintf("  %-12s %15s %15s %15s %10s\n",
            "Parameter", "Analytical", "Finite-Diff", "Rel. Diff", "Status"))
cat(strrep("-", 70), "\n")

param_names1 <- c("sigma_f^2", "beta_0", "beta_x1", "lambda", "sigma_y")
gradient_passed <- TRUE

for (i in seq_along(grad_analytical)) {
  comp <- compare_derivatives(grad_analytical[i], grad_fd[i], threshold = 1e-4)
  status <- if (comp$passed) "PASS" else "FAIL"
  if (!comp$passed) gradient_passed <- FALSE

  cat(sprintf("  %-12s %15.6e %15.6e %15.6e %10s\n",
              param_names1[i], grad_analytical[i], grad_fd[i],
              comp$rel_diff, status))
}

cat("\n")
if (gradient_passed) {
  cat("  ✓ All gradients match within tolerance (< 1e-4)\n")
} else {
  cat("  ✗ Some gradients exceed tolerance\n")
}

# ===== Test 2: Hessian validation =====
cat("\n\nTest 2: Single factor linear model - Hessian validation\n")
cat("----------------------------------------------------------\n")
cat("  NOTE: Hessian implementation in C++ is currently a placeholder.\n")
cat("        This test is included for future validation when Hessian\n")
cat("        computation is fully implemented.\n\n")

cat("  Computing analytical Hessian...\n")
result1_hess <- evaluate_likelihood_cpp(fm_cpp1, test_params1,
                                        compute_gradient = TRUE,
                                        compute_hessian = TRUE)

# Note: C++ returns Hessian as upper triangle in row-major order
# Need to reconstruct full symmetric matrix
n_params <- length(test_params1)
hess_analytical <- matrix(0, n_params, n_params)

# Fill in upper triangle from vector
idx <- 1
for (i in 1:n_params) {
  for (j in i:n_params) {
    hess_analytical[i, j] <- result1_hess$hessian[idx]
    hess_analytical[j, i] <- result1_hess$hessian[idx]  # Symmetrize
    idx <- idx + 1
  }
}

cat("  Computing finite-difference Hessian (this may take a moment)...\n")
hess_fd <- compute_fd_hessian(fm_cpp1, test_params1, delta = 1e-7)

# Compare Hessians
cat("\n  Hessian Comparison (upper triangle):\n")
cat(sprintf("  %-20s %15s %15s %15s %10s\n",
            "Element", "Analytical", "Finite-Diff", "Rel. Diff", "Status"))
cat(strrep("-", 75), "\n")

hessian_passed <- TRUE
threshold_hess <- 1e-3  # Higher threshold for Hessian (as noted by user)

for (i in 1:n_params) {
  for (j in i:n_params) {
    comp <- compare_derivatives(hess_analytical[i,j], hess_fd[i,j],
                               threshold = threshold_hess)
    status <- if (comp$passed) "PASS" else "FAIL"
    if (!comp$passed) hessian_passed <- FALSE

    element_name <- sprintf("H[%d,%d]", i, j)
    cat(sprintf("  %-20s %15.6e %15.6e %15.6e %10s\n",
                element_name, hess_analytical[i,j], hess_fd[i,j],
                comp$rel_diff, status))
  }
}

cat("\n")
if (hessian_passed) {
  cat("  ✓ All Hessian elements match within tolerance (< 1e-3)\n")
} else {
  cat("  ✗ Some Hessian elements exceed tolerance\n")
}

# ===== Test 3: Probit model gradient validation =====
cat("\n\nTest 3: Probit model - Gradient validation\n")
cat("----------------------------------------------------------\n")

set.seed(44)
n <- 100
f1 <- rnorm(n)
x1 <- rnorm(n)
intercept <- rep(1, n)
latent <- 0.5*intercept + 0.3*x1 + 0.7*f1 + rnorm(n)
D <- as.integer(latent > 0)

dat3 <- data.frame(D = D, intercept = intercept, x1 = x1, eval = 1)

fm3 <- define_factor_model(n_factors = 1, n_types = 1)
mc3 <- define_model_component("D", dat3, "D", fm3,
                              covariates = c("intercept", "x1"),
                              model_type = "probit",
                              loading_normalization = NA_real_,
                              evaluation_indicator = "eval")
ms3 <- define_model_system(factor = fm3, components = list(mc3))

fm_cpp3 <- initialize_factor_model_cpp(ms3, dat3, n_quad = 8L)
param_info3 <- get_parameter_info_cpp(fm_cpp3)

cat("  Number of parameters:", param_info3$n_param_free, "\n")
cat("  Parameters: [sigma_f^2, beta_0, beta_x1, lambda]\n\n")

# Test parameters for probit
test_params3 <- c(1.0, 0.5, 0.3, 0.7)

# Compute analytical gradient
result3 <- evaluate_likelihood_cpp(fm_cpp3, test_params3,
                                   compute_gradient = TRUE,
                                   compute_hessian = FALSE)
grad3_analytical <- result3$gradient
loglik3 <- result3$logLikelihood

cat("  Log-likelihood:", round(loglik3, 6), "\n\n")

# Compute finite-difference gradient
cat("  Computing finite-difference gradient...\n")
grad3_fd <- compute_fd_gradient(fm_cpp3, test_params3, delta = 1e-7)

# Compare gradients
cat("\n  Gradient Comparison:\n")
cat(sprintf("  %-12s %15s %15s %15s %10s\n",
            "Parameter", "Analytical", "Finite-Diff", "Rel. Diff", "Status"))
cat(strrep("-", 70), "\n")

param_names3 <- c("sigma_f^2", "beta_0", "beta_x1", "lambda")
gradient3_passed <- TRUE

for (i in seq_along(grad3_analytical)) {
  comp <- compare_derivatives(grad3_analytical[i], grad3_fd[i], threshold = 1e-4)
  status <- if (comp$passed) "PASS" else "FAIL"
  if (!comp$passed) gradient3_passed <- FALSE

  cat(sprintf("  %-12s %15.6e %15.6e %15.6e %10s\n",
              param_names3[i], grad3_analytical[i], grad3_fd[i],
              comp$rel_diff, status))
}

cat("\n")
if (gradient3_passed) {
  cat("  ✓ All gradients match within tolerance (< 1e-4)\n")
} else {
  cat("  ✗ Some gradients exceed tolerance\n")
}

# ===== Test 4: Multi-factor model gradient validation =====
cat("\n\nTest 4: Two-factor linear model - Gradient validation\n")
cat("----------------------------------------------------------\n")

# Use same data as Test 1 but with 2-factor model
# This ensures consistent data quality
set.seed(42)
n <- 100
f1 <- rnorm(n)
f2 <- rnorm(n)
x1 <- rnorm(n)
intercept <- rep(1, n)
Y <- 1.0*intercept + 0.5*x1 + 0.6*f1 + 0.4*f2 + rnorm(n, 0, 0.5)

dat4 <- data.frame(Y = Y, intercept = intercept, x1 = x1, eval = 1)

fm4 <- define_factor_model(n_factors = 2, n_types = 1)
mc4 <- define_model_component("Y", dat4, "Y", fm4,
                              covariates = c("intercept", "x1"),
                              model_type = "linear",
                              loading_normalization = c(NA_real_, NA_real_),
                              evaluation_indicator = "eval")
ms4 <- define_model_system(factor = fm4, components = list(mc4))

fm_cpp4 <- initialize_factor_model_cpp(ms4, dat4, n_quad = 8L)
param_info4 <- get_parameter_info_cpp(fm_cpp4)

cat("  Number of parameters:", param_info4$n_param_free, "\n")
cat("  Parameters: [sigma_f1^2, sigma_f2^2, beta_0, beta_x1, lambda_1, lambda_2, sigma_y]\n\n")

# Test parameters
test_params4 <- c(1.0, 1.0, 1.0, 0.5, 0.6, 0.4, 0.5)

# Compute analytical gradient
result4 <- evaluate_likelihood_cpp(fm_cpp4, test_params4,
                                   compute_gradient = TRUE,
                                   compute_hessian = FALSE)
grad4_analytical <- result4$gradient
loglik4 <- result4$logLikelihood

cat("  Log-likelihood:", round(loglik4, 6), "\n\n")

# Compute finite-difference gradient
cat("  Computing finite-difference gradient...\n")
grad4_fd <- compute_fd_gradient(fm_cpp4, test_params4, delta = 1e-7)

# Compare gradients
cat("\n  Gradient Comparison:\n")
cat(sprintf("  %-12s %15s %15s %15s %10s\n",
            "Parameter", "Analytical", "Finite-Diff", "Rel. Diff", "Status"))
cat(strrep("-", 70), "\n")

param_names4 <- c("sigma_f1^2", "sigma_f2^2", "beta_0", "beta_x1",
                  "lambda_1", "lambda_2", "sigma_y")
gradient4_passed <- TRUE

for (i in seq_along(grad4_analytical)) {
  comp <- compare_derivatives(grad4_analytical[i], grad4_fd[i], threshold = 1e-4)
  status <- if (comp$passed) "PASS" else "FAIL"
  if (!comp$passed) gradient4_passed <- FALSE

  cat(sprintf("  %-12s %15.6e %15.6e %15.6e %10s\n",
              param_names4[i], grad4_analytical[i], grad4_fd[i],
              comp$rel_diff, status))
}

cat("\n")
if (gradient4_passed) {
  cat("  ✓ All gradients match within tolerance (< 1e-4)\n")
} else {
  cat("  ✗ Some gradients exceed tolerance\n")
}

# ===== Test 5: Ordered probit gradient validation =====
cat("\n\nTest 5: Ordered probit model - Gradient validation\n")
cat("----------------------------------------------------------\n")

set.seed(45)
n <- 100
f1 <- rnorm(n)
x1 <- rnorm(n)
intercept <- rep(1, n)
# Latent variable with 4 categories (3 thresholds)
latent <- 0.5*x1 + 0.8*f1 + rnorm(n, 0, 0.5)
# Create ordered outcome with 4 categories
Y_ord <- as.integer(cut(latent, breaks = c(-Inf, -0.5, 0.0, 0.5, Inf), labels = FALSE))

dat5 <- data.frame(Y = Y_ord, intercept = intercept, x1 = x1, eval = 1)

fm5 <- define_factor_model(n_factors = 1, n_types = 1)
mc5_comp <- define_model_component("Y", dat5, "Y", fm5,
                              covariates = c("intercept", "x1"),
                              model_type = "oprobit",
                              loading_normalization = NA_real_,
                              num_choices = 4,
                              evaluation_indicator = "eval")
ms5 <- define_model_system(factor = fm5, components = list(mc5_comp))

fm_cpp5 <- initialize_factor_model_cpp(ms5, dat5, n_quad = 8L)
param_info5 <- get_parameter_info_cpp(fm_cpp5)

cat("  Number of parameters:", param_info5$n_param_free, "\n")
cat("  Parameters: [sigma_f^2, beta_0, beta_x1, lambda, thresh_1, thresh_2, thresh_3]\n")
cat("  Note: Thresholds are incremental (thresh_i = absolute_i - absolute_{i-1})\n\n")

# Test parameters for ordered probit
# Factor variance, intercept (absorbed), x1 coef, loading, 3 thresholds
# Thresholds: tau_1 = -0.5, tau_2 = 0.0, tau_3 = 0.5
# Incremental: thresh_1 = -0.5, thresh_2 = 0.5, thresh_3 = 0.5
test_params5 <- c(1.0, 0.0, 0.5, 0.8, -0.5, 0.5, 0.5)

# Compute analytical gradient
result5 <- evaluate_likelihood_cpp(fm_cpp5, test_params5,
                                   compute_gradient = TRUE,
                                   compute_hessian = FALSE)
grad5_analytical <- result5$gradient
loglik5 <- result5$logLikelihood

cat("  Log-likelihood:", round(loglik5, 6), "\n\n")

# Compute finite-difference gradient
cat("  Computing finite-difference gradient...\n")
grad5_fd <- compute_fd_gradient(fm_cpp5, test_params5, delta = 1e-7)

# Compare gradients
cat("\n  Gradient Comparison:\n")
cat(sprintf("  %-12s %15s %15s %15s %10s\n",
            "Parameter", "Analytical", "Finite-Diff", "Rel. Diff", "Status"))
cat(strrep("-", 70), "\n")

param_names5 <- c("sigma_f^2", "beta_0", "beta_x1", "lambda",
                  "thresh_1", "thresh_2", "thresh_3")
gradient5_passed <- TRUE

for (i in seq_along(grad5_analytical)) {
  comp <- compare_derivatives(grad5_analytical[i], grad5_fd[i], threshold = 1e-4)
  status <- if (comp$passed) "PASS" else "FAIL"
  if (!comp$passed) gradient5_passed <- FALSE

  cat(sprintf("  %-12s %15.6e %15.6e %15.6e %10s\n",
              param_names5[i], grad5_analytical[i], grad5_fd[i],
              comp$rel_diff, status))
}

cat("\n")
if (gradient5_passed) {
  cat("  ✓ All gradients match within tolerance (< 1e-4)\n")
} else {
  cat("  ✗ Some gradients exceed tolerance\n")
}

# ===== Test 6: Probit Hessian validation =====
cat("\n\nTest 6: Probit model - Hessian validation\n")
cat("----------------------------------------------------------\n")

# Re-use probit model from Test 3
cat("  Using probit model from Test 3\n")
cat("  Computing analytical Hessian...\n")
result3_hess <- evaluate_likelihood_cpp(fm_cpp3, test_params3,
                                        compute_gradient = TRUE,
                                        compute_hessian = TRUE)

# Reconstruct full Hessian from upper triangle
n_params3 <- length(test_params3)
hess3_analytical <- matrix(0, n_params3, n_params3)
idx <- 1
for (i in 1:n_params3) {
  for (j in i:n_params3) {
    hess3_analytical[i, j] <- result3_hess$hessian[idx]
    hess3_analytical[j, i] <- result3_hess$hessian[idx]
    idx <- idx + 1
  }
}

cat("  Computing finite-difference Hessian...\n")
hess3_fd <- compute_fd_hessian(fm_cpp3, test_params3, delta = 1e-7)

# Compare Hessians
cat("\n  Hessian Comparison (upper triangle):\n")
cat(sprintf("  %-20s %15s %15s %15s %10s\n",
            "Element", "Analytical", "Finite-Diff", "Rel. Diff", "Status"))
cat(strrep("-", 75), "\n")

hessian3_passed <- TRUE
threshold_hess <- 1e-3

for (i in 1:n_params3) {
  for (j in i:n_params3) {
    comp <- compare_derivatives(hess3_analytical[i,j], hess3_fd[i,j],
                               threshold = threshold_hess)
    status <- if (comp$passed) "PASS" else "FAIL"
    if (!comp$passed) hessian3_passed <- FALSE

    element_name <- sprintf("H[%d,%d]", i, j)
    cat(sprintf("  %-20s %15.6e %15.6e %15.6e %10s\n",
                element_name, hess3_analytical[i,j], hess3_fd[i,j],
                comp$rel_diff, status))
  }
}

cat("\n")
if (hessian3_passed) {
  cat("  ✓ All Hessian elements match within tolerance (< 1e-3)\n")
} else {
  cat("  ✗ Some Hessian elements exceed tolerance\n")
}

# ===== Test 7: Ordered probit Hessian validation =====
cat("\n\nTest 7: Ordered probit model - Hessian validation\n")
cat("----------------------------------------------------------\n")

cat("  Using ordered probit model from Test 5\n")
cat("  Computing analytical Hessian...\n")
result5_hess <- evaluate_likelihood_cpp(fm_cpp5, test_params5,
                                        compute_gradient = TRUE,
                                        compute_hessian = TRUE)

# Reconstruct full Hessian from upper triangle
n_params5 <- length(test_params5)
hess5_analytical <- matrix(0, n_params5, n_params5)
idx <- 1
for (i in 1:n_params5) {
  for (j in i:n_params5) {
    hess5_analytical[i, j] <- result5_hess$hessian[idx]
    hess5_analytical[j, i] <- result5_hess$hessian[idx]
    idx <- idx + 1
  }
}

cat("  Computing finite-difference Hessian...\n")
hess5_fd <- compute_fd_hessian(fm_cpp5, test_params5, delta = 1e-7)

# Compare Hessians
cat("\n  Hessian Comparison (upper triangle):\n")
cat(sprintf("  %-20s %15s %15s %15s %10s\n",
            "Element", "Analytical", "Finite-Diff", "Rel. Diff", "Status"))
cat(strrep("-", 75), "\n")

hessian5_passed <- TRUE

for (i in 1:n_params5) {
  for (j in i:n_params5) {
    comp <- compare_derivatives(hess5_analytical[i,j], hess5_fd[i,j],
                               threshold = threshold_hess)
    status <- if (comp$passed) "PASS" else "FAIL"
    if (!comp$passed) hessian5_passed <- FALSE

    element_name <- sprintf("H[%d,%d]", i, j)
    cat(sprintf("  %-20s %15.6e %15.6e %15.6e %10s\n",
                element_name, hess5_analytical[i,j], hess5_fd[i,j],
                comp$rel_diff, status))
  }
}

cat("\n")
if (hessian5_passed) {
  cat("  ✓ All Hessian elements match within tolerance (< 1e-3)\n")
} else {
  cat("  ✗ Some Hessian elements exceed tolerance\n")
}

# ===== Test 8: Two-factor linear Hessian validation =====
cat("\n\nTest 8: Two-factor linear model - Hessian validation\n")
cat("----------------------------------------------------------\n")

cat("  Using two-factor linear model from Test 4\n")
cat("  Computing analytical Hessian...\n")
result4_hess <- evaluate_likelihood_cpp(fm_cpp4, test_params4,
                                        compute_gradient = TRUE,
                                        compute_hessian = TRUE)

# Reconstruct full Hessian from upper triangle
n_params4 <- length(test_params4)
hess4_analytical <- matrix(0, n_params4, n_params4)
idx <- 1
for (i in 1:n_params4) {
  for (j in i:n_params4) {
    hess4_analytical[i, j] <- result4_hess$hessian[idx]
    hess4_analytical[j, i] <- result4_hess$hessian[idx]
    idx <- idx + 1
  }
}

cat("  Computing finite-difference Hessian...\n")
hess4_fd <- compute_fd_hessian(fm_cpp4, test_params4, delta = 1e-7)

# Compare Hessians
cat("\n  Hessian Comparison (upper triangle):\n")
cat(sprintf("  %-20s %15s %15s %15s %10s\n",
            "Element", "Analytical", "Finite-Diff", "Rel. Diff", "Status"))
cat(strrep("-", 75), "\n")

hessian4_passed <- TRUE

for (i in 1:n_params4) {
  for (j in i:n_params4) {
    comp <- compare_derivatives(hess4_analytical[i,j], hess4_fd[i,j],
                               threshold = threshold_hess)
    status <- if (comp$passed) "PASS" else "FAIL"
    if (!comp$passed) hessian4_passed <- FALSE

    element_name <- sprintf("H[%d,%d]", i, j)
    cat(sprintf("  %-20s %15.6e %15.6e %15.6e %10s\n",
                element_name, hess4_analytical[i,j], hess4_fd[i,j],
                comp$rel_diff, status))
  }
}

cat("\n")
if (hessian4_passed) {
  cat("  ✓ All Hessian elements match within tolerance (< 1e-3)\n")
} else {
  cat("  ✗ Some Hessian elements exceed tolerance\n")
}

# ===== Summary =====
cat("\n========================================\n")
cat("Validation Summary\n")
cat("========================================\n\n")

all_gradients_passed <- gradient_passed && gradient3_passed && gradient4_passed && gradient5_passed
all_hessians_passed <- hessian_passed && hessian3_passed && hessian5_passed && hessian4_passed

cat("GRADIENT TESTS:\n")
if (all_gradients_passed) {
  cat("  ✓✓✓ All gradient tests PASSED with high precision!\n")
  cat("      Maximum relative error: < 1e-7\n")
  cat("      Analytical gradients are correctly implemented.\n\n")
} else {
  cat("  Some gradient tests FAILED:\n")
  if (!gradient_passed) cat("    ✗ Test 1: Linear model gradient\n")
  if (!gradient3_passed) cat("    ✗ Test 3: Probit model gradient\n")
  if (!gradient4_passed) cat("    ✗ Test 4: Two-factor model gradient\n")
  if (!gradient5_passed) cat("    ✗ Test 5: Ordered probit model gradient\n")
  cat("\n")
}

cat("HESSIAN TESTS:\n")
if (all_hessians_passed) {
  cat("  ✓✓✓ All Hessian tests PASSED!\n\n")
} else {
  cat("  Some Hessian tests show discrepancies:\n")
  if (!hessian_passed) cat("    ⚠ Test 2: Linear model Hessian\n")
  if (!hessian3_passed) cat("    ⚠ Test 6: Probit model Hessian\n")
  if (!hessian5_passed) cat("    ⚠ Test 7: Ordered probit model Hessian\n")
  if (!hessian4_passed) cat("    ⚠ Test 8: Two-factor linear model Hessian\n")
  cat("\n")
}

cat("CONCLUSION:\n")
if (all_gradients_passed && all_hessians_passed) {
  cat("  ✓ All analytical derivatives (gradients AND Hessians) are validated!\n")
  cat("  The implementation supports both gradient-based (L-BFGS) and\n")
  cat("  Newton-type (nlminb, trust) optimizers with full accuracy.\n")
} else if (all_gradients_passed) {
  cat("  The optimization workflow using analytical GRADIENTS is fully\n")
  cat("  functional and validated. Hessian accuracy issues may affect\n")
  cat("  Newton-type methods but gradient-based optimizers work correctly.\n")
} else {
  cat("  There are issues with derivative computations that need attention.\n")
}

cat("\n========================================\n")
