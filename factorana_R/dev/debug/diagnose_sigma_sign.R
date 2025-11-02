#!/usr/bin/env Rscript
library(factorana)

set.seed(456)
n <- 50

# Simple case: one outcome with FIXED loading to identify factor variance
dat <- data.frame(
  intercept = 1,
  Y = rnorm(n, mean = 2, sd = 0.5),
  x1 = rnorm(n),
  eval = 1
)

fm <- define_factor_model(n_factors = 1, n_types = 1, n_quad = 8)

mc_Y <- define_model_component(
  name = "Y", data = dat, outcome = "Y", factor = fm,
  covariates = c("intercept", "x1"), model_type = "linear",
  loading_normalization = 1.0,  # FIXED loading
  evaluation_indicator = "eval"
)

ms <- define_model_system(components = list(mc_Y), factor = fm)
fm_cpp <- initialize_factor_model_cpp(ms, as.matrix(dat), n_quad = 8)

# Parameters: factor_var, intercept, beta_x1, sigma (lambda is fixed at 1.0)
test_params <- c(1.0, 2.0, 0.5, 0.5)

cat("=== Testing sigma Hessian with FIXED loading ===\n\n")

# Analytical
result <- evaluate_likelihood_cpp(fm_cpp, test_params, 
                                 compute_gradient = TRUE,
                                 compute_hessian = TRUE)

cat(sprintf("Log-likelihood: %.6f\n", result$loglik))
cat(sprintf("Gradient: %s\n\n", paste(round(result$gradient, 4), collapse=", ")))

# Convert hessian
n_params <- 4
hess_mat <- matrix(0, n_params, n_params)
idx <- 1
for (i in 1:n_params) {
  for (j in i:n_params) {
    hess_mat[i, j] <- result$hessian[idx]
    hess_mat[j, i] <- result$hessian[idx]
    idx <- idx + 1
  }
}

cat("Hessian diagonal:\n")
for (i in 1:n_params) {
  cat(sprintf("  [%d] %10.4f\n", i, hess_mat[i,i]))
}

# Finite difference for sigma (parameter 4)
delta <- 1e-6
h <- delta * (abs(test_params[4]) + 1.0)

params_plus <- test_params
params_plus[4] <- test_params[4] + h
grad_plus <- evaluate_likelihood_cpp(fm_cpp, params_plus,
                                     compute_gradient = TRUE,
                                     compute_hessian = FALSE)$gradient[4]

params_minus <- test_params
params_minus[4] <- test_params[4] - h
grad_minus <- evaluate_likelihood_cpp(fm_cpp, params_minus,
                                      compute_gradient = TRUE,
                                      compute_hessian = FALSE)$gradient[4]

hess_fd <- (grad_plus - grad_minus) / (2.0 * h)

cat("\n\nSigma Hessian comparison:\n")
cat(sprintf("  Analytical: %10.4f\n", hess_mat[4,4]))
cat(sprintf("  Finite Diff: %10.4f\n", hess_fd))
cat(sprintf("  Difference: %10.4f\n", hess_mat[4,4] - hess_fd))
cat(sprintf("  Ratio (should be ~1): %10.4f\n", hess_mat[4,4] / hess_fd))

if (abs(hess_mat[4,4] / hess_fd + 1.0) < 0.01) {
  cat("\n⚠️  FOUND IT: Ratio is ~-1, meaning OPPOSITE SIGNS!\n")
}
