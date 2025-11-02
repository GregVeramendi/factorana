#!/usr/bin/env Rscript
# Test: Verify parameter constraints and free parameter handling

library(factorana)

set.seed(456)
n <- 200  # Small sample for quick test

# Generate simple data
f <- rnorm(n)
T1 <- 1.0 + 1.0*f + rnorm(n, 0, 0.5)
T2 <- 0.5 + 1.2*f + rnorm(n, 0, 0.6)

dat <- data.frame(intercept=1, T1=T1, T2=T2, eval=1)

# Define model with factor variance fixed
fm <- define_factor_model(n_factors=1, n_types=1, n_quad=8)
mc_T1 <- define_model_component(name="T1", data=dat, outcome="T1", factor=fm,
  covariates="intercept", model_type="linear",
  loading_normalization=NA_real_, evaluation_indicator="eval")
mc_T2 <- define_model_component(name="T2", data=dat, outcome="T2", factor=fm,
  covariates="intercept", model_type="linear",
  loading_normalization=NA_real_, evaluation_indicator="eval")

ms <- define_model_system(components=list(mc_T1, mc_T2), factor=fm)

# Test estimation with the new parameter handling
cat("Testing estimate_model_rcpp with new parameter handling...\n")
cat("Expected behavior:\n")
cat("  - Factor variance (param 1) should be fixed to 1.0\n")
cat("  - Sigma parameters (params 4, 7) should have lower bound = 0.01\n")
cat("  - Only 6 free parameters should be optimized\n\n")

init_params <- c(1.0, 1.0, 1.0, 0.5, 0.5, 1.2, 0.6)
result <- estimate_model_rcpp(
  ms, dat,
  init_params = init_params,
  optimizer = "nlminb",
  parallel = FALSE,
  verbose = TRUE
)

cat(sprintf("\n✓ Estimation completed: LL = %.4f\n", result$loglik))
cat(sprintf("  Estimated parameters: %s\n", paste(round(result$estimates, 3), collapse=", ")))
