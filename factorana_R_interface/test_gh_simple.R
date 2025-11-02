#!/usr/bin/env Rscript
# Simple test: verify GH quadrature with one linear outcome
# Model: Y = intercept + lambda*f + epsilon
# where f ~ N(0,1) and epsilon ~ N(0, sigma^2)

library(factorana)

set.seed(456)
n <- 100

# True parameters
true_intercept <- 1.0
true_lambda <- 1.0
true_sigma <- 0.5
true_factor_var <- 1.0

# Generate data
f <- rnorm(n, 0, sqrt(true_factor_var))
Y <- true_intercept + true_lambda * f + rnorm(n, 0, true_sigma)

dat <- data.frame(intercept = 1, Y = Y, eval = 1)

# Define model
fm <- define_factor_model(n_factors = 1, n_types = 1, n_quad = 16)

mc_Y <- define_model_component(
  name = "Y", data = dat, outcome = "Y", factor = fm,
  covariates = "intercept", model_type = "linear",
  loading_normalization = NA_real_, evaluation_indicator = "eval"
)

ms <- define_model_system(components = list(mc_Y), factor = fm)
fm_cpp <- initialize_factor_model_cpp(ms, as.matrix(dat), n_quad = 16)

# Parameters: factor_var (fixed to 1), intercept, lambda, sigma
true_params <- c(true_factor_var, true_intercept, true_lambda, true_sigma)

cat("=== Testing Gauss-Hermite Integration ===\n\n")
cat("True parameters:\n")
cat(sprintf("  factor_var = %.2f (will be fixed)\n", true_factor_var))
cat(sprintf("  intercept  = %.2f\n", true_intercept))
cat(sprintf("  lambda     = %.2f\n", true_lambda))
cat(sprintf("  sigma      = %.2f\n", true_sigma))
cat("\n")

# Compute likelihood
result <- evaluate_likelihood_cpp(fm_cpp, true_params, TRUE, FALSE)
cat(sprintf("Log-likelihood at true params: %.6f\n\n", result$loglik))

# Now compute "by hand" for first observation to verify
cat("=== Manual Calculation for First Observation ===\n\n")
y1 <- dat$Y[1]

# Get quadrature nodes and weights
quad <- gauss_hermite_quadrature(16)
nodes <- quad$nodes
weights <- quad$weights

cat(sprintf("Y[1] = %.4f\n", y1))
cat(sprintf("Number of quadrature points: %d\n", length(nodes)))
cat(sprintf("Sum of weights: %.6f (should be sqrt(pi) = %.6f)\n",
    sum(weights), sqrt(pi)))

# Compute likelihood by integrating over f
manual_prob <- 0
for (i in seq_along(nodes)) {
  # Transform node: f = sqrt(2) * sqrt(factor_var) * node
  f_val <- sqrt(2.0) * sqrt(true_factor_var) * nodes[i]

  # Compute predicted Y
  y_pred <- true_intercept + true_lambda * f_val

  # Compute normal density for Y|f
  residual <- y1 - y_pred
  density <- (1/sqrt(2*pi*true_sigma^2)) * exp(-residual^2 / (2*true_sigma^2))

  # Add weighted contribution
  manual_prob <- manual_prob + weights[i] * density
}

cat(sprintf("\nManual calculation:\n"))
cat(sprintf("  Integrated probability: %.6e\n", manual_prob))
cat(sprintf("  Log-likelihood: %.6f\n", log(manual_prob)))

# Compare to C++ result for this observation
cat(sprintf("\nC++ calculation:\n"))
cat(sprintf("  Log-likelihood (all %d obs): %.6f\n", n, result$loglik))
cat(sprintf("  Average log-lik per obs: %.6f\n", result$loglik / n))

cat("\n=== Gradient Check ===\n\n")
delta <- 1e-7
param_names <- c("factor_var", "intercept", "lambda", "sigma")

for (i in 2:4) {  # Skip factor_var (will be fixed)
  h <- delta * (abs(true_params[i]) + 1.0)

  params_plus <- true_params
  params_plus[i] <- true_params[i] + h
  f_plus <- evaluate_loglik_only_cpp(fm_cpp, params_plus)

  params_minus <- true_params
  params_minus[i] <- true_params[i] - h
  f_minus <- evaluate_loglik_only_cpp(fm_cpp, params_minus)

  grad_fd <- (f_plus - f_minus) / (2.0 * h)
  grad_analytical <- result$gradient[i]

  rel_diff <- abs(grad_fd - grad_analytical) / (abs(grad_analytical) + 1e-10)
  status <- if (rel_diff < 1e-5) "✓" else "FAIL"

  cat(sprintf("[%d] %-12s: FD=%10.6f, Analytical=%10.6f, RelDiff=%.2e %s\n",
              i, param_names[i], grad_fd, grad_analytical, rel_diff, status))
}
