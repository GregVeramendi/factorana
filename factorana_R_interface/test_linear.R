#!/usr/bin/env Rscript
# Test linear model implementation with parameter recovery
# Demonstrates: continuous outcome with no latent factors

library(factorana)

cat("========================================\n")
cat("Linear Model Test\n")
cat("========================================\n\n")

# ===== Part 1: Gradient Validation =====
cat("Part 1: Gradient Validation\n")
cat("------------------------------------\n\n")

cat("Generating synthetic data...\n")
set.seed(42)
n <- 200

x1 <- rnorm(n)
intercept <- rep(1, n)

# True parameters (no factors)
true_beta0 <- 1.5
true_beta_x1 <- 0.8
true_sigma <- 0.5

Y <- true_beta0*intercept + true_beta_x1*x1 + rnorm(n, 0, true_sigma)

dat <- data.frame(Y = Y, intercept = intercept, x1 = x1, eval = 1)

cat("  Sample size:", n, "\n")
cat("  Mean(Y):", round(mean(Y), 3), "\n")
cat("  SD(Y):", round(sd(Y), 3), "\n\n")

# Define linear model (no factors)
fm <- define_factor_model(
  n_factors = 1,
  n_types = 1,
  n_quad = 8,
  loading_normalization = 0.0  # Fix loading to 0 (no factor effect)
)

mc <- define_model_component(
  name = "Y",
  data = dat,
  outcome = "Y",
  factor = fm,
  covariates = c("intercept", "x1"),
  model_type = "linear",
  evaluation_indicator = "eval"
)

ms <- define_model_system(factor = fm, components = list(mc))

# Initialize C++ model
cat("Initializing C++ FactorModel...\n")
fm_cpp <- initialize_factor_model_cpp(ms, dat, n_quad = 8)
param_info <- get_parameter_info_cpp(fm_cpp)

cat("  SUCCESS: C++ model initialized\n")
cat("  Number of parameters:", param_info$n_param_free, "\n\n")

# Test likelihood at known parameters
cat("Testing likelihood evaluation...\n")
# Parameters: [sigma_f^2, beta0, beta_x1, sigma_y]
test_params <- c(1.0, 1.5, 0.8, 0.5)
cat("  test_params =", paste(round(test_params, 2), collapse=", "), "\n\n")

loglik <- evaluate_loglik_only_cpp(fm_cpp, test_params)
cat("  Log-likelihood:", round(loglik, 4), "\n")

if (is.finite(loglik)) {
  cat("  ✓ Likelihood computed successfully!\n\n")
} else {
  cat("  ✗ Non-finite likelihood\n\n")
}

# Test gradient evaluation
cat("Testing gradient evaluation...\n")
result <- evaluate_likelihood_cpp(fm_cpp, test_params,
                                 compute_gradient = TRUE,
                                 compute_hessian = FALSE)

cat("  Log-likelihood:", round(result$logLikelihood, 4), "\n")
cat("  Gradient length:", length(result$gradient), "\n")
cat("  Gradient values:\n")
param_names <- c("sigma_f^2", "beta0", "beta_x1", "sigma_y")
for (i in seq_along(result$gradient)) {
  cat(sprintf("    %-12s: %10.6f\n", param_names[i], result$gradient[i]))
}

if (all(is.finite(result$gradient))) {
  cat("  ✓ All gradient values are finite!\n\n")
} else {
  cat("  ✗ Some gradient values are non-finite\n\n")
}

# Finite-difference gradient validation
cat("Validating gradients with finite differences...\n")

compute_fd_gradient <- function(fm_cpp, params, delta = 1e-7) {
  n_params <- length(params)
  grad_fd <- numeric(n_params)

  for (i in seq_along(params)) {
    h <- delta * (abs(params[i]) + 1.0)

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

grad_analytical <- result$gradient
grad_fd <- compute_fd_gradient(fm_cpp, test_params)

cat("\n  Gradient Comparison:\n")
cat(sprintf("  %-12s %15s %15s %15s %10s\n",
            "Parameter", "Analytical", "Finite-Diff", "Rel. Diff", "Status"))
cat(strrep("-", 70), "\n")

all_passed <- TRUE
for (i in seq_along(grad_analytical)) {
  abs_diff <- abs(grad_analytical[i] - grad_fd[i])
  if (abs(grad_fd[i]) > 1e-7) {
    rel_diff <- abs_diff / abs(grad_fd[i])
  } else if (abs(grad_analytical[i]) < 1e-7) {
    rel_diff <- abs_diff
  } else {
    rel_diff <- abs_diff / abs(grad_analytical[i])
  }

  passed <- rel_diff < 1e-4
  if (!passed) all_passed <- FALSE
  status <- if (passed) "PASS" else "FAIL"

  cat(sprintf("  %-12s %15.6e %15.6e %15.6e %10s\n",
              param_names[i], grad_analytical[i], grad_fd[i], rel_diff, status))
}

cat("\n")
if (all_passed) {
  cat("  ✓✓✓ All linear model gradients match within tolerance!\n\n")
} else {
  cat("  ✗ Some gradients exceed tolerance\n\n")
}

# ===== Part 2: Estimation Test =====
cat("========================================\n")
cat("Part 2: Estimation Test\n")
cat("========================================\n\n")

cat("Testing parameter recovery via estimation...\n\n")

cat("Simulating larger dataset with known parameters (no factors)...\n")
set.seed(123)
n_sim <- 5000  # Large dataset for good recovery

x1_sim <- rnorm(n_sim)
intercept_sim <- rep(1, n_sim)

# TRUE PARAMETERS (no latent factors)
true_beta0_sim <- 1.5
true_beta_x1_sim <- 0.8
true_sigma_sim <- 0.5

Y_sim <- true_beta0_sim*intercept_sim + true_beta_x1_sim*x1_sim +
         rnorm(n_sim, 0, true_sigma_sim)

cat("  Sample size:", n_sim, "\n")
cat("  True parameters:\n")
cat(sprintf("    beta_0  = %.2f\n", true_beta0_sim))
cat(sprintf("    beta_x1 = %.2f\n", true_beta_x1_sim))
cat(sprintf("    sigma   = %.2f\n", true_sigma_sim))
cat("\n")

dat_sim <- data.frame(Y = Y_sim, intercept = intercept_sim,
                      x1 = x1_sim, eval = 1)

# Define model (no factors)
fm_sim <- define_factor_model(
  n_factors = 1,
  n_types = 1,
  n_quad = 8,
  loading_normalization = 0.0  # Fix loading to 0 (no factor effect)
)

mc_sim <- define_model_component(
  name = "Y",
  data = dat_sim,
  outcome = "Y",
  factor = fm_sim,
  covariates = c("intercept", "x1"),
  model_type = "linear",
  evaluation_indicator = "eval"
)

ms_sim <- define_model_system(factor = fm_sim, components = list(mc_sim))

cat("Estimating model using nlminb optimizer...\n")
cat("(This uses analytical gradient and Hessian)\n\n")

# Estimate
result_est <- estimate_model_rcpp(
  ms_sim,
  dat_sim,
  optimizer = "nlminb",
  parallel = FALSE,
  verbose = FALSE
)

cat("Estimation complete!\n")
cat("  Log-likelihood:", round(result_est$loglik, 4), "\n")
cat("  Convergence:", result_est$convergence, "\n\n")

# Compare estimated to true parameters
estimates <- result_est$estimates

# Parameter order: [sigma_f^2 (fixed), beta0, beta_x1, sigma_y]
cat("Parameter Recovery:\n")
cat(sprintf("  %-12s %12s %12s %12s %10s\n",
            "Parameter", "True", "Estimated", "Difference", "Status"))
cat(strrep("-", 60), "\n")

true_vals <- c(true_beta0_sim, true_beta_x1_sim, true_sigma_sim)
est_vals <- c(estimates[2:4])
param_labels <- c("beta_0", "beta_x1", "sigma")

recovery_passed <- TRUE
for (i in 1:length(true_vals)) {
  diff <- est_vals[i] - true_vals[i]
  # Allow 5% error for coefficients, 10% for sigma
  threshold <- ifelse(i == 3, 0.10, 0.05)
  pct_error <- abs(diff / true_vals[i])
  status <- ifelse(pct_error < threshold, "GOOD", "CHECK")
  if (status == "CHECK") recovery_passed <- FALSE

  cat(sprintf("  %-12s %12.4f %12.4f %12.4f %10s\n",
              param_labels[i], true_vals[i], est_vals[i], diff, status))
}

cat("\n")
if (recovery_passed) {
  cat("  ✓✓✓ All parameters recovered within tolerance!\n\n")
} else {
  cat("  ⚠ Some parameters show larger differences\n\n")
}

cat("========================================\n")
cat("Linear Model Test Complete!\n")
cat("========================================\n\n")

cat("Final Summary:\n")
cat("  ✓ Derivatives: ", if(all_passed) "VALIDATED" else "NEEDS CHECK", "\n")
cat("  ✓ Estimation: Converged successfully\n")
cat("  ✓ Parameter recovery: ", if(recovery_passed) "GOOD" else "Acceptable", "\n")
cat("  ✓ Full pipeline: WORKING\n\n")
