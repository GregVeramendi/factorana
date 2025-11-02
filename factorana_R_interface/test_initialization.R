#!/usr/bin/env Rscript
# Test parameter initialization and factor identification

library(factorana)

set.seed(456)
n <- 500

# Simulate data with one factor
f1 <- rnorm(n)
x1 <- rnorm(n)
intercept <- rep(1, n)

# Two outcomes with different loading values
Y1 <- 0.5*intercept + 0.3*x1 + 1.0*f1 + rnorm(n, sd = 0.5)
Y2 <- -0.3*intercept + 0.5*x1 + 0.7*f1 + rnorm(n, sd = 0.5)

dat <- data.frame(
  Y1 = Y1,
  Y2 = Y2,
  intercept = intercept,
  x1 = x1,
  eval = 1
)

cat("=============================================================\n")
cat("TEST 1: Factor Variance IDENTIFIED (Component 1 has fixed loading)\n")
cat("=============================================================\n\n")

# Define factor model
fm1 <- define_factor_model(
  n_factors = 1,
  n_types = 1,
  n_quad = 8
)

# Component 1: FIXED loading at 1.0 (identifies factor variance)
mc1_identified <- define_model_component(
  name = "Y1",
  data = dat,
  outcome = "Y1",
  factor = fm1,
  covariates = c("intercept", "x1"),
  model_type = "linear",
  evaluation_indicator = "eval",
  loading_normalization = 1.0  # FIXED - identifies factor variance
)

# Component 2: FREE loading
mc2_identified <- define_model_component(
  name = "Y2",
  data = dat,
  outcome = "Y2",
  factor = fm1,
  covariates = c("intercept", "x1"),
  model_type = "linear",
  evaluation_indicator = "eval",
  loading_normalization = NA_real_  # FREE
)

ms_identified <- define_model_system(factor = fm1, components = list(mc1_identified, mc2_identified))

# Test initialization
cat("Testing initialization...\n")
init_result1 <- initialize_parameters(ms_identified, dat, verbose = TRUE)

cat("\nInitialized parameters:\n")
cat(sprintf("  Total: %d parameters\n", length(init_result1$init_params)))
cat(sprintf("  Factor variance fixed: %s\n", init_result1$factor_variance_fixed))
cat(sprintf("  Factor variance initial value: %.4f\n\n", init_result1$init_params[1]))

# Run estimation
cat("Running estimation with automatic initialization...\n")
result1 <- estimate_model_rcpp(
  ms_identified, dat,
  init_params = NULL,  # Will use automatic initialization
  optimizer = "nlminb",
  parallel = FALSE,
  verbose = TRUE
)

cat("\n=== RESULTS (Identified Case) ===\n")
cat(sprintf("Final log-likelihood: %.4f\n", result1$loglik))
cat(sprintf("Factor variance estimate: %.4f (should be estimated, not fixed)\n\n", result1$estimates[1]))


cat("\n=============================================================\n")
cat("TEST 2: Factor Variance NOT IDENTIFIED (All loadings free)\n")
cat("=============================================================\n\n")

# Define factor model with all free loadings
fm2 <- define_factor_model(
  n_factors = 1,
  n_types = 1,
  n_quad = 8
)

# Component 1: FREE loading (does NOT identify factor variance)
mc1_unidentified <- define_model_component(
  name = "Y1",
  data = dat,
  outcome = "Y1",
  factor = fm2,
  covariates = c("intercept", "x1"),
  model_type = "linear",
  evaluation_indicator = "eval",
  loading_normalization = NA_real_  # FREE
)

# Component 2: FREE loading
mc2_unidentified <- define_model_component(
  name = "Y2",
  data = dat,
  outcome = "Y2",
  factor = fm2,
  covariates = c("intercept", "x1"),
  model_type = "linear",
  evaluation_indicator = "eval",
  loading_normalization = NA_real_  # FREE
)

ms_unidentified <- define_model_system(factor = fm2, components = list(mc1_unidentified, mc2_unidentified))

# Test initialization
cat("Testing initialization...\n")
init_result2 <- initialize_parameters(ms_unidentified, dat, verbose = TRUE)

cat("\nInitialized parameters:\n")
cat(sprintf("  Total: %d parameters\n", length(init_result2$init_params)))
cat(sprintf("  Factor variance fixed: %s\n", init_result2$factor_variance_fixed))
cat(sprintf("  Factor variance initial value: %.4f\n\n", init_result2$init_params[1]))

# Run estimation
cat("Running estimation with automatic initialization...\n")
result2 <- estimate_model_rcpp(
  ms_unidentified, dat,
  init_params = NULL,  # Will use automatic initialization
  optimizer = "nlminb",
  parallel = FALSE,
  verbose = TRUE
)

cat("\n=== RESULTS (Unidentified Case) ===\n")
cat(sprintf("Final log-likelihood: %.4f\n", result2$loglik))
cat(sprintf("Factor variance estimate: %.4f (should be FIXED at 1.0)\n\n", result2$estimates[1]))


cat("\n=============================================================\n")
cat("SUMMARY\n")
cat("=============================================================\n\n")

cat("Test 1 (Identified):\n")
cat(sprintf("  Factor variance fixed: %s\n", init_result1$factor_variance_fixed))
cat(sprintf("  Factor variance estimate: %.4f\n", result1$estimates[1]))
if (!init_result1$factor_variance_fixed && abs(result1$estimates[1] - 1.0) > 0.1) {
  cat("  ✓ Factor variance was estimated (not fixed)\n\n")
} else {
  cat("  ⚠️  Warning: Factor variance may not have been estimated correctly\n\n")
}

cat("Test 2 (Unidentified):\n")
cat(sprintf("  Factor variance fixed: %s\n", init_result2$factor_variance_fixed))
cat(sprintf("  Factor variance estimate: %.4f\n", result2$estimates[1]))
if (init_result2$factor_variance_fixed && abs(result2$estimates[1] - 1.0) < 1e-6) {
  cat("  ✓ Factor variance was correctly fixed at 1.0\n\n")
} else {
  cat("  ⚠️  Warning: Factor variance should be fixed at 1.0\n\n")
}

cat("=============================================================\n")
cat("TEST 3: Verify initialization improves starting values\n")
cat("=============================================================\n\n")

# Get simple zero initialization
n_params <- length(init_result1$init_params)
zero_params <- c(1.0, rep(0.0, n_params - 1))

# Evaluate likelihood at zero initialization vs smart initialization
fm_cpp <- initialize_factor_model_cpp(ms_identified, as.matrix(dat), n_quad = 8)

ll_zero <- evaluate_loglik_only_cpp(fm_cpp, zero_params)
ll_init <- evaluate_loglik_only_cpp(fm_cpp, init_result1$init_params)
ll_final <- result1$loglik

cat(sprintf("Log-likelihood at zero initialization: %.4f\n", ll_zero))
cat(sprintf("Log-likelihood at smart initialization: %.4f\n", ll_init))
cat(sprintf("Log-likelihood at final estimates: %.4f\n", ll_final))
cat(sprintf("\nImprovement from zero to smart init: %.4f\n", ll_init - ll_zero))
cat(sprintf("Improvement from smart init to final: %.4f\n", ll_final - ll_init))

if (ll_init > ll_zero) {
  cat("\n✓ Smart initialization provides better starting point than zero initialization!\n")
} else {
  cat("\n⚠️  Smart initialization did not improve over zero initialization\n")
}
