#!/usr/bin/env Rscript
# Test if sigma Hessian issue is due to fabs() or quadrature

library(factorana)

set.seed(456)
n <- 100

# Simple data - one linear outcome, one test score, no missing data
dat <- data.frame(
  intercept = 1,
  Y = rnorm(n, mean = 2, sd = 0.5),
  T = rnorm(n, mean = 1, sd = 0.5),
  x1 = rnorm(n),
  q1 = rnorm(n),
  eval = 1
)

fm <- define_factor_model(n_factors = 1, n_types = 1, n_quad = 8)

mc_Y <- define_model_component(
  name = "Y", data = dat, outcome = "Y", factor = fm,
  covariates = c("intercept", "x1"), model_type = "linear",
  loading_normalization = NA_real_, evaluation_indicator = "eval"
)

mc_T <- define_model_component(
  name = "T", data = dat, outcome = "T", factor = fm,
  covariates = c("intercept", "q1"), model_type = "linear",
  loading_normalization = NA_real_, evaluation_indicator = "eval"
)

ms <- define_model_system(components = list(mc_Y, mc_T), factor = fm)
fm_cpp <- initialize_factor_model_cpp(ms, as.matrix(dat), n_quad = 8)

# Test parameters
# factor_var, Y: intercept, beta_x1, lambda, sigma, T: intercept, beta_q1, lambda, sigma
test_params <- c(1.0, 2.0, 0.5, 1.0, 0.5, 1.0, 0.5, 1.0, 0.5)
param_names <- c("factor_var", "Y_int", "Y_beta", "Y_lambda", "Y_sigma", 
                 "T_int", "T_beta", "T_lambda", "T_sigma")

# Get analytical Hessian
result <- evaluate_likelihood_cpp(fm_cpp, test_params, 
                                 compute_gradient = TRUE,
                                 compute_hessian = TRUE)

# Convert to matrix
n_params <- length(test_params)
hess_analytical <- matrix(0, n_params, n_params)
idx <- 1
for (i in 1:n_params) {
  for (j in i:n_params) {
    hess_analytical[i, j] <- result$hessian[idx]
    if (i != j) hess_analytical[j, i] <- result$hessian[idx]
    idx <- idx + 1
  }
}

# Compute FD Hessian
delta <- 1e-6
hess_fd <- matrix(0, n_params, n_params)

for (i in seq_along(test_params)) {
  h <- delta * (abs(test_params[i]) + 1.0)
  
  params_plus <- test_params
  params_plus[i] <- test_params[i] + h
  grad_plus <- evaluate_likelihood_cpp(fm_cpp, params_plus,
                                       compute_gradient = TRUE,
                                       compute_hessian = FALSE)$gradient
  
  params_minus <- test_params
  params_minus[i] <- test_params[i] - h
  grad_minus <- evaluate_likelihood_cpp(fm_cpp, params_minus,
                                        compute_gradient = TRUE,
                                        compute_hessian = FALSE)$gradient
  
  hess_fd[, i] <- (grad_plus - grad_minus) / (2.0 * h)
}

# Compare sigma diagonals only
sigma_indices <- c(5, 9)  # Y_sigma, T_sigma

cat("Sigma Hessian Diagonal Comparison:\n")
cat(sprintf("%-10s %15s %15s %12s\n", "Parameter", "Analytical", "Finite Diff", "RelDiff"))
cat(strrep("-", 60), "\n")

for (i in sigma_indices) {
  analytical_val <- hess_analytical[i, i]
  fd_val <- hess_fd[i, i]
  diff <- abs(analytical_val - fd_val)
  rel_diff <- diff / (abs(analytical_val) + 1e-6)
  
  status <- if (rel_diff < 1e-3) "✓" else "FAIL"
  
  cat(sprintf("%-10s %15.6f %15.6f %11.2e %s\n",
              param_names[i], analytical_val, fd_val, rel_diff, status))
}

# Check gradient sign too
cat("\n\nGradient check:\n")
cat(sprintf("Gradient at test params: %s\n", paste(round(result$gradient, 4), collapse=", ")))

# Is the gradient negative of what we expect?
# For minimizing -loglik, gradient should be negative of d(loglik)/dθ
# So if loglik increases when sigma increases, gradient for minimization should be negative

cat("\n\nCheck: is C++ returning Hessian of loglik or -loglik?\n")
cat("If Hessian is of loglik, it should be negated for minimization context.\n")
