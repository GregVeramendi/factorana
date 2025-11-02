#!/usr/bin/env Rscript
# Test component-specific loading normalizations

library(factorana)

set.seed(123)
n <- 500

# Simulate data with one factor
f1 <- rnorm(n)
x1 <- rnorm(n)
intercept <- rep(1, n)

# Outcome Y1: fixed loading at 1.0
Y1 <- 0.5*intercept + 0.3*x1 + 1.0*f1 + rnorm(n, sd = 0.5)

# Outcome Y2: free loading (true value = 0.7)
Y2 <- -0.3*intercept + 0.5*x1 + 0.7*f1 + rnorm(n, sd = 0.5)

dat <- data.frame(
  Y1 = Y1,
  Y2 = Y2,
  intercept = intercept,
  x1 = x1,
  eval = 1
)

cat("=== Testing Component-Specific Loading Normalizations ===\n\n")

# Define factor model
fm <- define_factor_model(
  n_factors = 1,
  n_types = 1,
  n_quad = 8
)

# Component 1: FIXED loading at 1.0
cat("Component 1: Fixed loading at 1.0\n")
mc1 <- define_model_component(
  name = "Y1",
  data = dat,
  outcome = "Y1",
  factor = fm,
  covariates = c("intercept", "x1"),
  model_type = "linear",
  evaluation_indicator = "eval",
  loading_normalization = 1.0  # FIXED at 1.0
)

cat("  Number of parameters:", mc1$nparam_model, "\n")
cat("  Loading normalization:", mc1$loading_normalization, "\n")
print(mc1)
cat("\n")

# Component 2: FREE loading
cat("Component 2: Free loading\n")
mc2 <- define_model_component(
  name = "Y2",
  data = dat,
  outcome = "Y2",
  factor = fm,
  covariates = c("intercept", "x1"),
  model_type = "linear",
  evaluation_indicator = "eval",
  loading_normalization = NA_real_  # FREE
)

cat("  Number of parameters:", mc2$nparam_model, "\n")
cat("  Loading normalization:", mc2$loading_normalization, "\n")
print(mc2)
cat("\n")

# Create model system
ms <- define_model_system(factor = fm, components = list(mc1, mc2))

# Initialize C++ object
fm_cpp <- initialize_factor_model_cpp(ms, as.matrix(dat), n_quad = 8)

# Get parameter info
param_info <- get_parameter_info_cpp(fm_cpp)
cat("Total number of parameters:", param_info$n_param_free, "\n")
cat("Expected: sigma_f^2(1) + mc1_params(3: beta0, beta_x1, sigma) + mc2_params(4: beta0, beta_x1, lambda, sigma) = 8\n\n")

# Test with known parameters
# sigma_f^2, mc1: beta0, beta_x1, sigma, mc2: beta0, beta_x1, lambda, sigma
test_params <- c(1.0, 0.5, 0.3, 0.5, -0.3, 0.5, 0.7, 0.5)
param_names <- c("sigma_f^2", "Y1_beta0", "Y1_beta_x1", "Y1_sigma",
                 "Y2_beta0", "Y2_beta_x1", "Y2_lambda", "Y2_sigma")

cat("=== Testing Likelihood Evaluation ===\n")
result <- evaluate_likelihood_cpp(fm_cpp, test_params,
                                 compute_gradient = TRUE,
                                 compute_hessian = FALSE)

cat("Log-likelihood:", result$loglik, "\n\n")

cat("=== Gradient Check (Finite Differences) ===\n")
delta <- 1e-7

max_diff <- 0
for (i in seq_along(test_params)) {
  h <- delta * (abs(test_params[i]) + 1.0)

  params_plus <- test_params
  params_plus[i] <- test_params[i] + h
  f_plus <- evaluate_loglik_only_cpp(fm_cpp, params_plus)

  params_minus <- test_params
  params_minus[i] <- test_params[i] - h
  f_minus <- evaluate_loglik_only_cpp(fm_cpp, params_minus)

  grad_fd <- (f_plus - f_minus) / (2.0 * h)

  diff <- abs(grad_fd - result$gradient[i])
  max_diff <- max(max_diff, diff)

  status <- if (diff < 1e-4) "✓" else "FAIL"
  cat(sprintf("[%d] %-12s: FD=%10.6f, Analytical=%10.6f, Diff=%10.6f %s\n",
              i, param_names[i], grad_fd, result$gradient[i], diff, status))
}

cat("\nMax gradient difference:", sprintf("%.8f", max_diff), "\n")
if (max_diff < 1e-4) {
  cat("✓ All gradients match!\n")
} else {
  cat("⚠️  Some gradients don't match\n")
}

cat("\n=== Parameter Estimation ===\n")
init_params <- c(1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)

result <- estimate_model_rcpp(
  ms, dat,
  init_params = init_params,
  optimizer = "nlminb",
  parallel = FALSE,
  verbose = TRUE
)

cat("\n=== Parameter Recovery ===\n")
true_params <- c(1.0, 0.5, 0.3, 0.5, -0.3, 0.5, 0.7, 0.5)

cat(sprintf("%-15s %10s %10s %10s %10s\n",
            "Parameter", "True", "Estimated", "Diff", "% Error"))
cat(strrep("-", 65), "\n")

for (i in seq_along(true_params)) {
  diff <- result$estimates[i] - true_params[i]
  pct_error <- 100 * diff / true_params[i]

  status <- if (abs(pct_error) < 10) "✓" else "⚠"

  cat(sprintf("%-15s %10.4f %10.4f %10.4f %9.1f%% %s\n",
              param_names[i], true_params[i], result$estimates[i],
              diff, pct_error, status))
}
