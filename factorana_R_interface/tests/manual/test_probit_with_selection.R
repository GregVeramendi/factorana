#!/usr/bin/env Rscript
# Test: Probit selection + one outcome observed only when D=1 + test scores
# Mimics Roy model structure but simpler

library(factorana)

set.seed(666)
n <- 1000

# True parameters
true_factor_var <- 1.0  # FIXED

# Probit selection
true_D_intercept <- 0.5
true_D_beta_Z <- 1.0
true_D_lambda <- 1.0  # free

# Y1 outcome (only observed when D=1)
true_Y1_intercept <- 2.0
true_Y1_beta_X <- 1.0
true_Y1_lambda <- 2.0  # free
true_Y1_sigma <- 0.5

# Test scores (observed for everyone)
true_T1_intercept <- 1.0
true_T1_lambda <- 1.0  # free
true_T1_sigma <- 0.5

true_T2_intercept <- 0.5
true_T2_lambda <- 1.2  # free
true_T2_sigma <- 0.6

# Generate data
f <- rnorm(n, 0, sqrt(true_factor_var))
Z <- rnorm(n, 0, 1)
X <- rnorm(n, 0, 1)

# Selection
D_latent <- true_D_intercept + true_D_beta_Z * Z + true_D_lambda * f + rnorm(n, 0, 1)
D <- as.integer(D_latent > 0)

# Y1 (generated for all but only "observed" when D=1)
Y1 <- true_Y1_intercept + true_Y1_beta_X * X + true_Y1_lambda * f + rnorm(n, 0, true_Y1_sigma)

# Test scores
T1 <- true_T1_intercept + true_T1_lambda * f + rnorm(n, 0, true_T1_sigma)
T2 <- true_T2_intercept + true_T2_lambda * f + rnorm(n, 0, true_T2_sigma)

# Evaluation indicators
dat <- data.frame(
  intercept = 1, D = D, Z = Z, X = X,
  Y1 = Y1, T1 = T1, T2 = T2,
  eval_all = 1,      # Everyone
  eval_Y1 = D        # Only D=1
)

cat("=== Selection Model with Missing Outcome ===\n\n")
cat(sprintf("n = %d observations\n", n))
cat(sprintf("Selection rate (D=1): %.1f%%\n", 100 * mean(D)))
cat(sprintf("Y1 observations: %d (only when D=1)\n", sum(D)))
cat(sprintf("\nTrue parameters:\n"))
cat(sprintf("  Factor variance: %.2f (FIXED)\n", true_factor_var))
cat(sprintf("  Selection: intercept=%.2f, beta_Z=%.2f, lambda=%.2f\n",
    true_D_intercept, true_D_beta_Z, true_D_lambda))
cat(sprintf("  Y1 (D=1 only): intercept=%.2f, beta_X=%.2f, lambda=%.2f, sigma=%.2f\n",
    true_Y1_intercept, true_Y1_beta_X, true_Y1_lambda, true_Y1_sigma))
cat(sprintf("  T1: intercept=%.2f, lambda=%.2f, sigma=%.2f\n",
    true_T1_intercept, true_T1_lambda, true_T1_sigma))
cat(sprintf("  T2: intercept=%.2f, lambda=%.2f, sigma=%.2f\n\n",
    true_T2_intercept, true_T2_lambda, true_T2_sigma))

# Define model
fm <- define_factor_model(n_factors = 1, n_types = 1, n_quad = 16)

mc_D <- define_model_component(
  name = "D", data = dat, outcome = "D", factor = fm,
  covariates = c("intercept", "Z"), model_type = "probit",
  loading_normalization = NA_real_, evaluation_indicator = "eval_all"
)

mc_Y1 <- define_model_component(
  name = "Y1", data = dat, outcome = "Y1", factor = fm,
  covariates = c("intercept", "X"), model_type = "linear",
  loading_normalization = NA_real_, evaluation_indicator = "eval_Y1"  # KEY: Only D=1
)

mc_T1 <- define_model_component(
  name = "T1", data = dat, outcome = "T1", factor = fm,
  covariates = "intercept", model_type = "linear",
  loading_normalization = NA_real_, evaluation_indicator = "eval_all"
)

mc_T2 <- define_model_component(
  name = "T2", data = dat, outcome = "T2", factor = fm,
  covariates = "intercept", model_type = "linear",
  loading_normalization = NA_real_, evaluation_indicator = "eval_all"
)

ms <- define_model_system(components = list(mc_D, mc_Y1, mc_T1, mc_T2), factor = fm)
fm_cpp <- initialize_factor_model_cpp(ms, as.matrix(dat), n_quad = 16)

# Parameters
true_params <- c(true_factor_var,
                 true_D_intercept, true_D_beta_Z, true_D_lambda,
                 true_Y1_intercept, true_Y1_beta_X, true_Y1_lambda, true_Y1_sigma,
                 true_T1_intercept, true_T1_lambda, true_T1_sigma,
                 true_T2_intercept, true_T2_lambda, true_T2_sigma)

cat("=== Gradient Check at True Parameters ===\n\n")
result <- evaluate_likelihood_cpp(fm_cpp, true_params, TRUE, FALSE)
cat(sprintf("Log-likelihood: %.6f\n\n", result$loglik))

delta <- 1e-7
param_names <- c("factor_var", "D_int", "D_beta_Z", "D_lambda",
                 "Y1_int", "Y1_beta_X", "Y1_lambda", "Y1_sigma",
                 "T1_int", "T1_lambda", "T1_sigma",
                 "T2_int", "T2_lambda", "T2_sigma")

max_rel_diff <- 0
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
  max_rel_diff <- max(max_rel_diff, rel_diff)
  status <- if (rel_diff < 1e-4) "✓" else "FAIL"

  cat(sprintf("[%2d] %-12s: FD=%11.6f, Analytical=%11.6f, RelDiff=%.2e %s\n",
              i, param_names[i], grad_fd, grad_analytical, rel_diff, status))
}

cat(sprintf("\nMax gradient relative difference: %.2e\n", max_rel_diff))

cat("\n=== Parameter Estimation ===\n\n")

init_params <- c(
  1.0,  # factor_var (FIXED)
  0.0, 0.5, 0.8,  # D
  1.5, 0.8, 1.5, 1.0,  # Y1
  mean(dat$T1), 0.8, sd(dat$T1),  # T1
  mean(dat$T2), 0.8, sd(dat$T2)   # T2
)

param_fixed <- c(TRUE, rep(FALSE, 13))

cat("Starting estimation...\n\n")

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

final_params <- init_params
final_params[!param_fixed] <- opt_result$par

cat(sprintf("Converged in %d iterations\n", opt_result$iterations))
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
    errors[i] <- 0
  } else {
    diff <- final_params[i] - true_params[i]
    pct_error <- 100 * diff / true_params[i]
    errors[i] <- abs(pct_error)
    status <- if (abs(pct_error) < 15) "✓" else "⚠"
    cat(sprintf("%-12s %8.4f %10.4f %10.4f %9.1f%% %s\n",
                param_names[i], true_params[i], final_params[i], diff, pct_error, status))
  }
}

cat(strrep("-", 60), "\n")
cat(sprintf("\nMax error: %.1f%%\n", max(errors)))

ll_true <- evaluate_loglik_only_cpp(fm_cpp, true_params)
cat(sprintf("\nLog-likelihood at true params: %.6f\n", ll_true))
cat(sprintf("Log-likelihood at estimates: %.6f\n", -opt_result$objective))
cat(sprintf("Difference: %.6f ", -opt_result$objective - ll_true))
if (-opt_result$objective > ll_true) {
  cat("✓\n")
} else {
  cat("⚠ (TRUE PARAMS BETTER!)\n")
}
