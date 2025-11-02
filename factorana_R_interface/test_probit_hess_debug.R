#!/usr/bin/env Rscript
# Check if probit also has factor variance Hessian issues

library(factorana)

set.seed(42)
n <- 50
f1 <- rnorm(n)
x1 <- rnorm(n)
intercept <- rep(1, n)

latent <- 0.5*intercept + 0.3*x1 + 0.7*f1 + rnorm(n)
Y_binary <- as.numeric(latent > 0)

dat <- data.frame(Y = Y_binary, intercept = intercept, x1 = x1, eval = 1)

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
  model_type = "probit",
  evaluation_indicator = "eval"
)

ms <- define_model_system(factor = fm, components = list(mc))
fm_cpp <- initialize_factor_model_cpp(ms, dat, n_quad = 8)

test_params <- c(1.0, 0.5, 0.3, 0.7)

# Compute analytical Hessian
result <- evaluate_likelihood_cpp(fm_cpp, test_params,
                                  compute_gradient = TRUE,
                                  compute_hessian = TRUE)

if (length(result$hessian) > 0) {
  # Convert to matrix
  n_par <- length(test_params)
  hess_mat <- matrix(0, n_par, n_par)
  idx <- 1
  for (i in 1:n_par) {
    for (j in i:n_par) {
      hess_mat[i, j] <- result$hessian[idx]
      hess_mat[j, i] <- result$hessian[idx]
      idx <- idx + 1
    }
  }

  cat("Probit Hessian (full matrix):\n")
  print(round(hess_mat, 4))

  # Finite-difference for diagonal
  delta <- 1e-5
  h <- delta * 2.0

  params_plus <- test_params
  params_plus[1] <- test_params[1] + h
  f_plus <- evaluate_loglik_only_cpp(fm_cpp, params_plus)

  params_minus <- test_params
  params_minus[1] <- test_params[1] - h
  f_minus <- evaluate_loglik_only_cpp(fm_cpp, params_minus)

  f_center <- evaluate_loglik_only_cpp(fm_cpp, test_params)

  hess_fd_11 <- (f_plus - 2*f_center + f_minus) / (h^2)

  cat("\nProbit H[1,1] (sigma_f^2):\n")
  cat("  Finite-difference:", hess_fd_11, "\n")
  cat("  Analytical:", hess_mat[1,1], "\n")
  cat("  Ratio:", hess_mat[1,1] / hess_fd_11, "\n")
} else {
  cat("Probit Hessian not implemented\n")
}
