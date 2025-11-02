#!/usr/bin/env Rscript
# End-to-end test of factor model estimation
# Demonstrates: model specification → C++ initialization → optimization → results

library(factorana)

cat("========================================\n")
cat("Factor Model Estimation Test\n")
cat("========================================\n\n")

# ===== Generate synthetic data =====
cat("Generating synthetic data...\n")
set.seed(42)
n <- 200
f1 <- rnorm(n)  # True factor values
x1 <- rnorm(n)
x2 <- rnorm(n)
intercept <- rep(1, n)  # Intercept column

# True parameters: beta0=1.5, beta1=0.8, beta2=-0.5, loading=1.2, sigma=0.5
Y <- 1.5*intercept + 0.8*x1 - 0.5*x2 + 1.2*f1 + rnorm(n, 0, 0.5)

dat <- data.frame(Y = Y, intercept = intercept, x1 = x1, x2 = x2, eval = 1)
cat("  Sample size:", n, "\n")
cat("  True parameters: intercept=1.5, beta1=0.8, beta2=-0.5, loading=1.2, sigma=0.5\n\n")

# ===== Define model =====
cat("Defining factor model...\n")
fm <- define_factor_model(
  n_factors = 1,
  n_types = 1,
  n_quad = 8,
  loading_normalization = NA_real_  # Free loading parameter (numeric NA)
)

mc <- define_model_component(
  name = "Y",
  data = dat,
  outcome = "Y",
  factor = fm,
  covariates = c("intercept", "x1", "x2"),
  model_type = "linear",
  evaluation_indicator = "eval"
)

ms <- define_model_system(factor = fm, components = list(mc))
cat("  Model: Y ~ intercept + x1 + x2 + factor + error\n")
cat("  Parameters: [sigma_f^2, beta_0, beta_x1, beta_x2, lambda, sigma_y]\n\n")

# ===== Test likelihood at true parameters =====
cat("Testing likelihood at approximate true parameters...\n")
# Parameter order: [sigma_f^2=1.44, beta_0=1.5, beta_x1=0.8, beta_x2=-0.5, lambda=1.2, sigma_y=0.5]
true_params <- c(1.44, 1.5, 0.8, -0.5, 1.2, 0.5)

fm_cpp <- initialize_factor_model_cpp(ms, dat, n_quad = 8)
param_info <- get_parameter_info_cpp(fm_cpp)
cat("  C++ expects", param_info$n_param_free, "parameters\n")

loglik_true <- evaluate_loglik_only_cpp(fm_cpp, true_params)
cat("  Log-likelihood at true params:", round(loglik_true, 4), "\n\n")

# ===== Run estimation =====
cat("Running estimation (this may take a minute)...\n")
cat("  Optimizer: L-BFGS-B via optim\n")
cat("  Using analytical gradients from C++\n\n")

# Initial parameters (rough starting values)
init_params <- c(1.0, 0.0, 0.0, 0.0, 1.0, 1.0)

result <- estimate_model_rcpp(
  model_system = ms,
  data = dat,
  init_params = init_params,
  optimizer = "optim",
  parallel = FALSE,  # Single-threaded for simplicity
  verbose = TRUE
)

# ===== Display results =====
cat("\n========================================\n")
cat("Estimation Results\n")
cat("========================================\n\n")

param_names <- c("sigma_f^2", "beta_0", "beta_x1", "beta_x2", "lambda", "sigma_y")
true_vals <- c(1.44, 1.5, 0.8, -0.5, 1.2, 0.5)

cat(sprintf("%-12s %10s %10s %10s\n", "Parameter", "True", "Estimate", "Diff"))
cat(strrep("-", 45), "\n")
for (i in seq_along(param_names)) {
  diff_pct <- 100 * (result$estimates[i] - true_vals[i]) / true_vals[i]
  cat(sprintf("%-12s %10.4f %10.4f %9.1f%%\n",
              param_names[i], true_vals[i], result$estimates[i], diff_pct))
}

cat("\nLog-likelihood:\n")
cat(sprintf("  At true params: %10.4f\n", loglik_true))
cat(sprintf("  At estimates:   %10.4f\n", result$loglik))
cat(sprintf("  Improvement:    %10.4f\n", result$loglik - loglik_true))

cat("\nConvergence: ", result$convergence, "\n")
cat("(0 = converged successfully)\n")

cat("\n========================================\n")
cat("Test Complete!\n")
cat("========================================\n")
