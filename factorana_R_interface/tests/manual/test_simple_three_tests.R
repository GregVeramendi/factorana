#!/usr/bin/env Rscript
# Simplest possible test: 3 test scores, 1 factor, all linear

library(factorana)

cat("=============================================================\n")
cat("Simple Factor Model: 3 Test Scores\n")
cat("=============================================================\n\n")

set.seed(123)
n <- 500

# Generate data with known parameters
f <- rnorm(n, 0, 1)  # factor with variance = 1

# Test scores with known parameters
T1 <- 1.0 + 0.5 * f + rnorm(n, 0, 0.3)  # intercept=1, lambda=0.5, sigma=0.3
T2 <- 2.0 + 0.8 * f + rnorm(n, 0, 0.4)  # intercept=2, lambda=0.8, sigma=0.4
T3 <- 0.5 + 1.0 * f + rnorm(n, 0, 0.5)  # intercept=0.5, lambda=1.0, sigma=0.5

dat <- data.frame(
  intercept = 1,
  T1 = T1,
  T2 = T2,
  T3 = T3,
  eval = 1
)

cat("Generated data:\n")
cat(sprintf("  n = %d observations\n", n))
cat(sprintf("  True factor variance = 1.0\n\n"))

# Define model
fm <- define_factor_model(n_factors = 1, n_types = 1, n_quad = 8)

# T1: fixed loading at 1.0 (to identify factor variance)
mc_T1 <- define_model_component(
  name = "T1", data = dat, outcome = "T1", factor = fm,
  covariates = "intercept", model_type = "linear",
  loading_normalization = 1.0, evaluation_indicator = "eval"
)

# T2: free loading
mc_T2 <- define_model_component(
  name = "T2", data = dat, outcome = "T2", factor = fm,
  covariates = "intercept", model_type = "linear",
  loading_normalization = NA_real_, evaluation_indicator = "eval"
)

# T3: free loading
mc_T3 <- define_model_component(
  name = "T3", data = dat, outcome = "T3", factor = fm,
  covariates = "intercept", model_type = "linear",
  loading_normalization = NA_real_, evaluation_indicator = "eval"
)

ms <- define_model_system(components = list(mc_T1, mc_T2, mc_T3), factor = fm)
fm_cpp <- initialize_factor_model_cpp(ms, as.matrix(dat), n_quad = 8)

# Parameters: factor_var, T1_int, T1_sigma, T2_int, T2_lambda, T2_sigma, T3_int, T3_lambda, T3_sigma
true_params <- c(1.0, 1.0, 0.3, 2.0, 0.8, 0.4, 0.5, 1.0, 0.5)
param_names <- c("factor_var", "T1_int", "T1_sigma", 
                 "T2_int", "T2_lambda", "T2_sigma",
                 "T3_int", "T3_lambda", "T3_sigma")

cat("=== PART 1: Gradient Validation ===\n\n")

result <- evaluate_likelihood_cpp(fm_cpp, true_params,
                                  compute_gradient = TRUE,
                                  compute_hessian = FALSE)

cat(sprintf("Log-likelihood at true params: %.4f\n\n", result$loglik))

# Finite difference gradient check
cat("Gradient validation:\n")
delta <- 1e-7
max_diff <- 0

for (i in seq_along(true_params)) {
  h <- delta * (abs(true_params[i]) + 1.0)
  
  params_plus <- true_params
  params_plus[i] <- true_params[i] + h
  f_plus <- evaluate_loglik_only_cpp(fm_cpp, params_plus)
  
  params_minus <- true_params
  params_minus[i] <- true_params[i] - h
  f_minus <- evaluate_loglik_only_cpp(fm_cpp, params_minus)
  
  grad_fd <- (f_plus - f_minus) / (2.0 * h)
  grad_analytical <- result$gradient[i]
  
  diff <- abs(grad_fd - grad_analytical)
  rel_diff <- diff / (abs(grad_analytical) + 1e-10)
  max_diff <- max(max_diff, rel_diff)
  
  status <- if (rel_diff < 1e-4) "✓" else "FAIL"
  
  cat(sprintf("[%d] %-12s: FD=%10.6f, Analytical=%10.6f, RelDiff=%.2e %s\n",
              i, param_names[i], grad_fd, grad_analytical, rel_diff, status))
}

cat(sprintf("\nMax gradient relative difference: %.2e\n", max_diff))

cat("\n=== PART 2: Hessian Validation ===\n\n")

# Get analytical Hessian
result_hess <- evaluate_likelihood_cpp(fm_cpp, true_params,
                                       compute_gradient = TRUE,
                                       compute_hessian = TRUE)

# Convert to matrix
n_params <- length(true_params)
hess_analytical <- matrix(0, n_params, n_params)
idx <- 1
for (i in 1:n_params) {
  for (j in i:n_params) {
    hess_analytical[i, j] <- result_hess$hessian[idx]
    hess_analytical[j, i] <- result_hess$hessian[idx]
    idx <- idx + 1
  }
}

# Compute FD Hessian
cat("Computing Hessian via finite differences...\n")
hess_fd <- matrix(0, n_params, n_params)
delta_hess <- 1e-6

for (i in seq_along(true_params)) {
  h <- delta_hess * (abs(true_params[i]) + 1.0)
  
  params_plus <- true_params
  params_plus[i] <- true_params[i] + h
  grad_plus <- evaluate_likelihood_cpp(fm_cpp, params_plus,
                                       compute_gradient = TRUE,
                                       compute_hessian = FALSE)$gradient
  
  params_minus <- true_params
  params_minus[i] <- true_params[i] - h
  grad_minus <- evaluate_likelihood_cpp(fm_cpp, params_minus,
                                        compute_gradient = TRUE,
                                        compute_hessian = FALSE)$gradient
  
  hess_fd[, i] <- (grad_plus - grad_minus) / (2.0 * h)
}

# Compare ONLY sigma parameters (indices 3, 6, 9)
sigma_indices <- c(3, 6, 9)

cat("\nHessian DIAGONAL comparison (sigma parameters only):\n")
cat(sprintf("%-12s %15s %15s %12s\n", "Parameter", "Analytical", "Finite Diff", "RelDiff"))
cat(strrep("-", 65), "\n")

for (i in sigma_indices) {
  analytical_val <- hess_analytical[i, i]
  fd_val <- hess_fd[i, i]
  diff <- abs(analytical_val - fd_val)
  rel_diff <- diff / (abs(analytical_val) + 1e-6)
  
  status <- if (rel_diff < 1e-3) "✓" else "FAIL"
  
  cat(sprintf("%-12s %15.6f %15.6f %11.2e %s\n",
              param_names[i], analytical_val, fd_val, rel_diff, status))
}

cat("\nAll other parameters:\n")
other_indices <- setdiff(2:n_params, sigma_indices)  # Skip factor_var (fixed)
for (i in other_indices) {
  analytical_val <- hess_analytical[i, i]
  fd_val <- hess_fd[i, i]
  diff <- abs(analytical_val - fd_val)
  rel_diff <- diff / (abs(analytical_val) + 1e-6)
  
  status <- if (rel_diff < 1e-3) "✓" else "FAIL"
  
  cat(sprintf("%-12s %15.6f %15.6f %11.2e %s\n",
              param_names[i], analytical_val, fd_val, rel_diff, status))
}

cat("\n=============================================================\n")
