#!/usr/bin/env Rscript
# Test: Do we recover factor loadings correctly in a simple three-test model?

library(factorana)

set.seed(999)
n <- 1000

# True parameters
true_T1_intercept <- 1.0
true_T1_lambda <- 1.0  # This will be normalized/fixed
true_T1_sigma <- 0.5

true_T2_intercept <- 0.5
true_T2_lambda <- 1.2
true_T2_sigma <- 0.6

true_T3_intercept <- -0.3
true_T3_lambda <- 0.8
true_T3_sigma <- 0.4

true_factor_var <- 1.0

# Generate data
f <- rnorm(n, 0, sqrt(true_factor_var))
T1 <- true_T1_intercept + true_T1_lambda * f + rnorm(n, 0, true_T1_sigma)
T2 <- true_T2_intercept + true_T2_lambda * f + rnorm(n, 0, true_T2_sigma)
T3 <- true_T3_intercept + true_T3_lambda * f + rnorm(n, 0, true_T3_sigma)

dat <- data.frame(intercept = 1, T1 = T1, T2 = T2, T3 = T3, eval = 1)

cat("=== Three Test Scores - Loading Recovery Test ===\n\n")
cat(sprintf("n = %d observations\n", n))
cat(sprintf("True parameters:\n"))
cat(sprintf("  T1: intercept=%.2f, lambda=%.2f (fixed), sigma=%.2f\n",
    true_T1_intercept, true_T1_lambda, true_T1_sigma))
cat(sprintf("  T2: intercept=%.2f, lambda=%.2f, sigma=%.2f\n",
    true_T2_intercept, true_T2_lambda, true_T2_sigma))
cat(sprintf("  T3: intercept=%.2f, lambda=%.2f, sigma=%.2f\n",
    true_T3_intercept, true_T3_lambda, true_T3_sigma))
cat(sprintf("  Factor variance: %.2f (fixed)\n\n", true_factor_var))

# Define model - fix T1 loading to 1.0 for identification
fm <- define_factor_model(n_factors = 1, n_types = 1, n_quad = 16)

mc_T1 <- define_model_component(
  name = "T1", data = dat, outcome = "T1", factor = fm,
  covariates = "intercept", model_type = "linear",
  loading_normalization = 1.0, evaluation_indicator = "eval"
)

mc_T2 <- define_model_component(
  name = "T2", data = dat, outcome = "T2", factor = fm,
  covariates = "intercept", model_type = "linear",
  loading_normalization = NA_real_, evaluation_indicator = "eval"
)

mc_T3 <- define_model_component(
  name = "T3", data = dat, outcome = "T3", factor = fm,
  covariates = "intercept", model_type = "linear",
  loading_normalization = NA_real_, evaluation_indicator = "eval"
)

ms <- define_model_system(components = list(mc_T1, mc_T2, mc_T3), factor = fm)
fm_cpp <- initialize_factor_model_cpp(ms, as.matrix(dat), n_quad = 16)

# Parameters: factor_var (fixed), T1_int, T1_sigma, T2_int, T2_lambda, T2_sigma, T3_int, T3_lambda, T3_sigma
true_params <- c(true_factor_var,
                 true_T1_intercept, true_T1_sigma,
                 true_T2_intercept, true_T2_lambda, true_T2_sigma,
                 true_T3_intercept, true_T3_lambda, true_T3_sigma)

cat("=== Gradient Check at True Parameters ===\n\n")
result <- evaluate_likelihood_cpp(fm_cpp, true_params, TRUE, FALSE)
cat(sprintf("Log-likelihood: %.6f\n\n", result$loglik))

delta <- 1e-7
param_names <- c("factor_var", "T1_int", "T1_sigma", "T2_int", "T2_lambda", "T2_sigma",
                 "T3_int", "T3_lambda", "T3_sigma")

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

  rel_diff <- abs(grad_fd - grad_analytical) / (abs(grad_analytical) + 1e-10)
  status <- if (rel_diff < 1e-5) "✓" else "FAIL"

  cat(sprintf("[%d] %-12s: FD=%10.6f, Analytical=%10.6f, RelDiff=%.2e %s\n",
              i, param_names[i], grad_fd, grad_analytical, rel_diff, status))
}

cat("\n=== Parameter Estimation ===\n\n")

# Initialize parameters
init_params <- c(
  1.0,  # factor_var (will be fixed)
  mean(dat$T1), sd(dat$T1),  # T1: intercept, sigma
  mean(dat$T2), 0.8, sd(dat$T2),  # T2: intercept, lambda, sigma
  mean(dat$T3), 0.8, sd(dat$T3)  # T3: intercept, lambda, sigma
)

# Fix factor variance
param_fixed <- c(TRUE, rep(FALSE, 8))

cat("Starting parameter estimation...\n")
opt_result <- nlminb(
  init_params[!param_fixed],
  function(p) {
    full_p <- init_params
    full_p[!param_fixed] <- p
    -evaluate_loglik_only_cpp(fm_cpp, full_p)
  },
  gradient = function(p) {
    full_p <- init_params
    full_p[!param_fixed] <- p
    grad_full <- evaluate_likelihood_cpp(fm_cpp, full_p, TRUE, FALSE)$gradient
    -grad_full[!param_fixed]
  },
  control = list(eval.max = 2000, iter.max = 2000)
)

# Extract final parameters
final_params <- init_params
final_params[!param_fixed] <- opt_result$par

cat(sprintf("\nConverged in %d iterations\n", opt_result$iterations))
cat(sprintf("Final log-likelihood: %.6f\n\n", -opt_result$objective))

cat("=== Parameter Recovery ===\n\n")
cat(sprintf("%-12s %8s %10s %10s %10s\n",
            "Parameter", "True", "Estimated", "Diff", "% Error"))
cat(strrep("-", 60), "\n")

errors <- numeric(length(param_names))
for (i in seq_along(param_names)) {
  if (param_fixed[i]) {
    cat(sprintf("%-12s %8.4f %10.4f %10s %10s (fixed)\n",
                param_names[i], true_params[i], final_params[i], "-", "-"))
  } else {
    diff <- final_params[i] - true_params[i]
    pct_error <- 100 * diff / true_params[i]
    errors[i] <- abs(pct_error)
    status <- if (abs(pct_error) < 10) "✓" else "⚠"
    cat(sprintf("%-12s %8.4f %10.4f %10.4f %9.1f%% %s\n",
                param_names[i], true_params[i], final_params[i], diff, pct_error, status))
  }
}

cat(strrep("-", 60), "\n")
cat(sprintf("\nMax error: %.1f%%\n", max(errors)))

# Evaluate at true params
ll_true <- evaluate_loglik_only_cpp(fm_cpp, true_params)
cat(sprintf("\nLog-likelihood at true params: %.6f\n", ll_true))
cat(sprintf("Log-likelihood at estimates: %.6f\n", -opt_result$objective))
cat(sprintf("Difference: %.6f ", -opt_result$objective - ll_true))
if (-opt_result$objective > ll_true) {
  cat("✓ (estimates better, as expected)\n")
} else {
  cat("⚠ (true params better - problem!)\n")
}
