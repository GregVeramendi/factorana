#!/usr/bin/env Rscript
# Test probit Hessian with simple binary outcome model

library(factorana)

set.seed(789)
n <- 500

# Generate simple probit data: D = 1(intercept + beta*X + lambda*f > 0)
f <- rnorm(n, 0, 1)
X <- rnorm(n)
latent <- 0.5 + 1.0 * X + 1.0 * f + rnorm(n, 0, 1)
D <- as.numeric(latent > 0)

dat <- data.frame(intercept = 1, D = D, X = X, eval = 1)

cat("Probit model test\n")
cat(sprintf("  n = %d observations\n", n))
cat(sprintf("  Selection rate: %.1f%%\n", 100 * mean(D)))

# Define model
fm <- define_factor_model(n_factors = 1, n_types = 1, n_quad = 8)

mc_D <- define_model_component(
  name = "D", data = dat, outcome = "D", factor = fm,
  covariates = c("intercept", "X"), model_type = "probit",
  loading_normalization = NA_real_, evaluation_indicator = "eval"
)

ms <- define_model_system(components = list(mc_D), factor = fm)
fm_cpp <- initialize_factor_model_cpp(ms, as.matrix(dat), n_quad = 8)

# Parameters: factor_var, intercept, beta_X, lambda
true_params <- c(1.0, 0.5, 1.0, 1.0)
param_names <- c("factor_var", "intercept", "beta_X", "lambda")

cat("\n=== Gradient Validation ===\n\n")

result <- evaluate_likelihood_cpp(fm_cpp, true_params, TRUE, FALSE)

cat(sprintf("Log-likelihood: %.4f\n\n", result$loglik))

# Finite difference gradient check
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

cat("\n=== Hessian Validation ===\n\n")

# Get analytical Hessian
result_hess <- evaluate_likelihood_cpp(fm_cpp, true_params, TRUE, TRUE)

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
  grad_plus <- evaluate_likelihood_cpp(fm_cpp, params_plus, TRUE, FALSE)$gradient

  params_minus <- true_params
  params_minus[i] <- true_params[i] - h
  grad_minus <- evaluate_likelihood_cpp(fm_cpp, params_minus, TRUE, FALSE)$gradient

  hess_fd[, i] <- (grad_plus - grad_minus) / (2.0 * h)
}

cat("\nHessian DIAGONAL comparison:\n")
cat(sprintf("%-12s %15s %15s %12s\n", "Parameter", "Analytical", "Finite Diff", "RelDiff"))
cat(strrep("-", 65), "\n")

for (i in 1:n_params) {
  analytical_val <- hess_analytical[i, i]
  fd_val <- hess_fd[i, i]
  diff <- abs(analytical_val - fd_val)
  rel_diff <- diff / (abs(analytical_val) + 1e-6)

  status <- if (rel_diff < 1e-3) "✓" else "FAIL"

  cat(sprintf("%-12s %15.6f %15.6f %11.2e %s\n",
              param_names[i], analytical_val, fd_val, rel_diff, status))
}

cat("\nFull analytical Hessian:\n")
print(round(hess_analytical, 4))

cat("\nFull FD Hessian:\n")
print(round(hess_fd, 4))
