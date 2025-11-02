#!/usr/bin/env Rscript
# Test: Does using analytical Hessian lead to worse optimum?

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

# Setup model (abbreviated for speed)
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
param_fixed <- c(TRUE, rep(FALSE, 23))

# Setup box constraints for sigma parameters
# Sigma parameters are at positions 8, 12, 16, 20, 24 in full params
# Which are positions 7, 11, 15, 19, 23 in free params
n_free <- sum(!param_fixed)
lower_bounds <- rep(-Inf, n_free)
upper_bounds <- rep(Inf, n_free)
lower_bounds[c(7, 11, 15, 19, 23)] <- 0.01  # Sigma parameters must be positive

cat("=== TEST 1: nlminb WITHOUT Hessian (quasi-Newton) ===\n\n")

opt_no_hess <- nlminb(
  true_params[!param_fixed],
  function(p) {
    full_p <- true_params
    full_p[!param_fixed] <- p
    -evaluate_loglik_only_cpp(fm_cpp, full_p)
  },
  gradient = function(p) {
    full_p <- true_params
    full_p[!param_fixed] <- p
    -evaluate_likelihood_cpp(fm_cpp, full_p, TRUE, FALSE)$gradient[!param_fixed]
  },
  lower = lower_bounds,
  upper = upper_bounds,
  control = list(eval.max = 2000, iter.max = 2000)
)

cat(sprintf("Converged in %d iterations\n", opt_no_hess$iterations))
cat(sprintf("Final LL: %.6f\n", -opt_no_hess$objective))
cat(sprintf("Convergence: %d (%s)\n\n", opt_no_hess$convergence, opt_no_hess$message))

cat("=== TEST 2: nlminb WITH Hessian (Newton method) ===\n\n")

opt_with_hess <- nlminb(
  true_params[!param_fixed],
  function(p) {
    full_p <- true_params
    full_p[!param_fixed] <- p
    -evaluate_loglik_only_cpp(fm_cpp, full_p)
  },
  gradient = function(p) {
    full_p <- true_params
    full_p[!param_fixed] <- p
    -evaluate_likelihood_cpp(fm_cpp, full_p, TRUE, FALSE)$gradient[!param_fixed]
  },
  hessian = function(p) {
    full_p <- true_params
    full_p[!param_fixed] <- p
    result <- evaluate_likelihood_cpp(fm_cpp, full_p, TRUE, TRUE)
    # Extract free parameter Hessian (upper triangle)
    n_free <- sum(!param_fixed)
    hess_mat <- matrix(0, n_free, n_free)
    idx <- 1
    for (i in 1:n_free) {
      for (j in i:n_free) {
        hess_mat[i,j] <- result$hessian[idx]
        hess_mat[j,i] <- result$hessian[idx]
        idx <- idx + 1
      }
    }
    -hess_mat  # Negative because we're minimizing negative log-likelihood
  },
  lower = lower_bounds,
  upper = upper_bounds,
  control = list(eval.max = 2000, iter.max = 2000)
)

cat(sprintf("Converged in %d iterations\n", opt_with_hess$iterations))
cat(sprintf("Final LL: %.6f\n", -opt_with_hess$objective))
cat(sprintf("Convergence: %d (%s)\n\n", opt_with_hess$convergence, opt_with_hess$message))

cat("=== COMPARISON ===\n\n")
cat(sprintf("LL without Hessian: %.6f\n", -opt_no_hess$objective))
cat(sprintf("LL with Hessian:    %.6f\n", -opt_with_hess$objective))
cat(sprintf("Difference:         %.6f\n\n", -opt_with_hess$objective - (-opt_no_hess$objective)))

if (-opt_no_hess$objective > -opt_with_hess$objective) {
  cat("✓ WITHOUT Hessian found BETTER optimum!\n")
  cat("  → Analytical Hessian may be leading optimizer to suboptimal solution\n")
} else if (-opt_with_hess$objective > -opt_no_hess$objective) {
  cat("✓ WITH Hessian found BETTER optimum!\n")
  cat("  → Newton method is superior (as expected)\n")
} else {
  cat("≈ Both found same optimum\n")
}
