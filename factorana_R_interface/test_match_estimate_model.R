#!/usr/bin/env Rscript
# Test: Exactly match estimate_model_rcpp's approach
# Use FULL parameter vector and fix factor variance via box constraints

library(factorana)

set.seed(456)
n <- 2000

# Generate same data as Roy model
f <- rnorm(n, 0, 1)
U1 <- rnorm(n, 0, 1) * 0.5
U0 <- rnorm(n, 0, 1) * 0.5
UV <- rnorm(n, 0, 1)
UT1 <- rnorm(n, 0, 1) * 0.5
UT2 <- rnorm(n, 0, 1) * 0.6
UT3 <- rnorm(n, 0, 1) * 0.4

intercept <- rep(1, n)
Z1 <- rnorm(n, 0, 1)
X1 <- rnorm(n, 0, 1)
Q1 <- rnorm(n, 0, 1)
Q2 <- rnorm(n, 0, 1)
Q3 <- rnorm(n, 0, 1)

Y1 <- 2.0 * intercept + 1.0 * X1 + 2.0 * f + U1
Y0 <- 1.0 * intercept + 1.0 * X1 + 1.0 * f + U0
I_latent <- 0.5 * intercept + 1.0 * Z1 + 1.0 * f + UV
D <- as.integer(I_latent > 0)
T1 <- 1.0 * intercept + 1.0 * Q1 + 1.0 * f + UT1
T2 <- 0.5 * intercept + 0.8 * Q2 + 1.2 * f + UT2
T3 <- -0.5 * intercept + 1.1 * Q3 + 0.9 * f + UT3

dat <- data.frame(
  intercept = intercept, D = D, Z1 = Z1, X1 = X1,
  Y1 = Y1, Y0 = Y0, Q1 = Q1, Q2 = Q2, Q3 = Q3,
  T1 = T1, T2 = T2, T3 = T3,
  eval_Y1 = D, eval_Y0 = 1L - D, eval_all = 1L
)

# Setup model (same as test_roy_model.R)
fm <- define_factor_model(n_factors = 1, n_types = 1, n_quad = 16)
mc_sel <- define_model_component(name = "Selection", data = dat, outcome = "D", factor = fm,
  evaluation_indicator = "eval_all", covariates = c("intercept", "Z1"),
  model_type = "probit", loading_normalization = NA_real_)
mc_Y1 <- define_model_component(name = "Y1", data = dat, outcome = "Y1", factor = fm,
  evaluation_indicator = "eval_Y1", covariates = c("intercept", "X1"),
  model_type = "linear", loading_normalization = NA_real_)
mc_Y0 <- define_model_component(name = "Y0", data = dat, outcome = "Y0", factor = fm,
  evaluation_indicator = "eval_Y0", covariates = c("intercept", "X1"),
  model_type = "linear", loading_normalization = NA_real_)
mc_T1 <- define_model_component(name = "T1", data = dat, outcome = "T1", factor = fm,
  evaluation_indicator = "eval_all", covariates = c("intercept", "Q1"),
  model_type = "linear", loading_normalization = NA_real_)
mc_T2 <- define_model_component(name = "T2", data = dat, outcome = "T2", factor = fm,
  evaluation_indicator = "eval_all", covariates = c("intercept", "Q2"),
  model_type = "linear", loading_normalization = NA_real_)
mc_T3 <- define_model_component(name = "T3", data = dat, outcome = "T3", factor = fm,
  evaluation_indicator = "eval_all", covariates = c("intercept", "Q3"),
  model_type = "linear", loading_normalization = NA_real_)

ms <- define_model_system(components = list(mc_sel, mc_Y1, mc_Y0, mc_T1, mc_T2, mc_T3), factor = fm)
fm_cpp <- initialize_factor_model_cpp(ms, as.matrix(dat), n_quad = 16)

true_params <- c(1.0, 0.5, 1.0, 1.0, 2.0, 1.0, 2.0, 0.5, 1.0, 1.0, 1.0, 0.5,
                 1.0, 1.0, 1.0, 0.5, 0.5, 0.8, 1.2, 0.6, -0.5, 1.1, 0.9, 0.4)

# Setup box constraints EXACTLY like estimate_model_rcpp
n_params <- length(true_params)
lower_bounds <- rep(-Inf, n_params)
upper_bounds <- rep(Inf, n_params)

# Fix factor variance (parameter 1) to 1.0
lower_bounds[1] <- 1.0
upper_bounds[1] <- 1.0

cat("=== TEST 1: nlminb WITHOUT Hessian (like estimate_model_rcpp but no Hessian) ===\n\n")

opt_no_hess <- nlminb(
  true_params,  # FULL parameter vector
  function(params) {
    -evaluate_loglik_only_cpp(fm_cpp, params)
  },
  gradient = function(params) {
    -evaluate_likelihood_cpp(fm_cpp, params, TRUE, FALSE)$gradient
  },
  lower = lower_bounds,
  upper = upper_bounds,
  control = list(eval.max = 2000, iter.max = 2000)
)

cat(sprintf("Converged in %d iterations\n", opt_no_hess$iterations))
cat(sprintf("Final LL: %.6f\n", -opt_no_hess$objective))
cat(sprintf("Convergence: %d (%s)\n\n", opt_no_hess$convergence, opt_no_hess$message))

cat("=== TEST 2: nlminb WITH Hessian (exactly like estimate_model_rcpp) ===\n\n")

opt_with_hess <- nlminb(
  true_params,  # FULL parameter vector
  function(params) {
    -evaluate_loglik_only_cpp(fm_cpp, params)
  },
  gradient = function(params) {
    -evaluate_likelihood_cpp(fm_cpp, params, TRUE, FALSE)$gradient
  },
  hessian = function(params) {
    result <- evaluate_likelihood_cpp(fm_cpp, params, FALSE, TRUE)
    # Convert upper triangle to full symmetric matrix
    hess_mat <- matrix(0, n_params, n_params)
    idx <- 1
    for (i in 1:n_params) {
      for (j in i:n_params) {
        hess_mat[i,j] <- result$hessian[idx]
        hess_mat[j,i] <- result$hessian[idx]
        idx <- idx + 1
      }
    }
    -hess_mat  # Negative for minimization
  },
  lower = lower_bounds,
  upper = upper_bounds,
  control = list(eval.max = 2000, iter.max = 2000)
)

cat(sprintf("Converged in %d iterations\n", opt_with_hess$iterations))
cat(sprintf("Final LL: %.6f\n", -opt_with_hess$objective))
cat(sprintf("Convergence: %d (%s)\n\n", opt_with_hess$convergence, opt_with_hess$message))

cat("=== TEST 3: Using estimate_model_rcpp (the wrapper function) ===\n\n")

result_est <- estimate_model_rcpp(
  ms, dat,
  init_params = true_params,
  optimizer = "nlminb",
  parallel = FALSE,
  verbose = FALSE
)

cat(sprintf("Final LL: %.6f\n", result_est$loglik))
cat("opt_result structure:\n")
print(str(result_est$opt_result))
cat("\n")

cat("=== COMPARISON ===\n\n")
cat(sprintf("LL without Hessian (manual):  %.6f\n", -opt_no_hess$objective))
cat(sprintf("LL with Hessian (manual):     %.6f\n", -opt_with_hess$objective))
cat(sprintf("LL estimate_model_rcpp:       %.6f\n", result_est$loglik))
cat(sprintf("\nDifference (manual - wrapper): %.6f\n\n", -opt_with_hess$objective - result_est$loglik))

# Compare parameters
cat("Parameter comparison (first 10 params):\n")
cat(sprintf("%-10s %12s %12s %12s\n", "Param", "True", "Manual", "Wrapper"))
cat(strrep("-", 50), "\n")
for (i in 1:min(10, length(true_params))) {
  cat(sprintf("%-10d %12.6f %12.6f %12.6f\n",
              i, true_params[i], opt_with_hess$par[i], result_est$estimates[i]))
}

if (abs(-opt_with_hess$objective - result_est$loglik) < 0.01) {
  cat("✓ estimate_model_rcpp matches manual nlminb\n")
} else {
  cat("⚠ estimate_model_rcpp finds DIFFERENT optimum!\n")
  if (result_est$loglik < -opt_with_hess$objective) {
    cat("  → Wrapper found WORSE optimum\n")
  } else {
    cat("  → Wrapper found BETTER optimum\n")
  }
}
