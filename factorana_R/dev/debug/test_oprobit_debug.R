#!/usr/bin/env Rscript
# Debug ordered probit Hessian issue

library(factorana)

set.seed(42)
n <- 50  # Smaller sample for debugging
f1 <- rnorm(n)
x1 <- rnorm(n)
intercept <- rep(1, n)

latent <- 0.5*intercept + 0.3*x1 + 0.7*f1 + rnorm(n)
true_thresholds <- c(-1, 0, 1)
Y_ordered <- cut(latent, breaks = c(-Inf, true_thresholds, Inf), labels = FALSE)

dat <- data.frame(Y = Y_ordered, intercept = intercept, x1 = x1, eval = 1)

fm <- define_factor_model(
  n_factors = 1,
  n_types = 1,
  n_quad = 8,
  loading_normalization = NA_real_
)

mc <- define_model_component(
  name = "Y",
  data = dat,
  outcome = "Y",
  factor = fm,
  covariates = c("intercept", "x1"),
  model_type = "oprobit",
  evaluation_indicator = "eval",
  num_choices = 4
)

ms <- define_model_system(factor = fm, components = list(mc))
fm_cpp <- initialize_factor_model_cpp(ms, dat, n_quad = 8)

test_params <- c(1.0, 0.5, 0.3, 0.7, -1.0, 1.0, 1.0)

# Compute analytical Hessian
result <- evaluate_likelihood_cpp(fm_cpp, test_params,
                                  compute_gradient = TRUE,
                                  compute_hessian = TRUE)

# Convert to matrix
hess_mat <- matrix(0, 7, 7)
idx <- 1
for (i in 1:7) {
  for (j in i:7) {
    hess_mat[i, j] <- result$hessian[idx]
    hess_mat[j, i] <- result$hessian[idx]
    idx <- idx + 1
  }
}

cat("Analytical Hessian (full matrix):\n")
print(round(hess_mat, 4))

cat("\n\nFocus on factor variance row/col (index 1):\n")
cat("H[1,1] (sigma, sigma):", hess_mat[1,1], "\n")
cat("H[1,2] (sigma, beta0):", hess_mat[1,2], "\n")
cat("H[1,3] (sigma, beta1):", hess_mat[1,3], "\n")
cat("H[1,4] (sigma, lambda):", hess_mat[1,4], "\n")

# Compute finite-difference for just the diagonal
delta <- 1e-5
h <- delta * 2.0  # For sigma_f^2 = 1.0

params_plus <- test_params
params_plus[1] <- test_params[1] + h
f_plus <- evaluate_loglik_only_cpp(fm_cpp, params_plus)

params_minus <- test_params
params_minus[1] <- test_params[1] - h
f_minus <- evaluate_loglik_only_cpp(fm_cpp, params_minus)

f_center <- evaluate_loglik_only_cpp(fm_cpp, test_params)

hess_fd_11 <- (f_plus - 2*f_center + f_minus) / (h^2)

cat("\nFinite-difference H[1,1]:", hess_fd_11, "\n")
cat("Analytical H[1,1]:", hess_mat[1,1], "\n")
cat("Ratio:", hess_mat[1,1] / hess_fd_11, "\n")
