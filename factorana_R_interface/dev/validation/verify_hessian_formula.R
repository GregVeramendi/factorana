#!/usr/bin/env Rscript
# Verify the Hessian formula by checking the actual values being returned

library(factorana)

set.seed(456)
n <- 20  # Small for clarity

dat <- data.frame(
  intercept = 1,
  Y = c(1.8, 2.1, 2.3, 1.9, 2.0, 2.2, 1.7, 2.4, 2.1, 1.9,
        2.0, 2.2, 1.8, 2.1, 2.3, 2.0, 1.9, 2.1, 2.2, 2.0),
  x1 = rnorm(n),
  eval = 1
)

fm <- define_factor_model(n_factors = 1, n_types = 1, n_quad = 8)

mc_Y <- define_model_component(
  name = "Y", data = dat, outcome = "Y", factor = fm,
  covariates = c("intercept", "x1"), model_type = "linear",
  loading_normalization = 1.0,  # FIXED
  evaluation_indicator = "eval"
)

ms <- define_model_system(components = list(mc_Y), factor = fm)
fm_cpp <- initialize_factor_model_cpp(ms, as.matrix(dat), n_quad = 8)

# params: factor_var, intercept, beta_x1, sigma
params <- c(1.0, 2.0, 0.5, 0.5)

result <- evaluate_likelihood_cpp(fm_cpp, params, TRUE, TRUE)

cat("=== Checking if C++ returns Hessian for ALL parameters or only FREE ===\n\n")
cat("Number of input parameters:", length(params), "\n")
cat("Number of gradient elements returned:", length(result$gradient), "\n")
cat("Number of Hessian elements returned:", length(result$hessian), "\n")
cat("Expected Hessian size for 4x4 upper triangle:", 4*5/2, "\n\n")

if (length(result$hessian) == 10) {
  cat("✓ Hessian is for ALL 4 parameters\n")
} else {
  cat("⚠️  Hessian size doesn't match!\n")
}

# Check which parameters are supposed to be fixed
cat("\n=== Checking parameter fixing ===\n")
cat("Factor_var should be FIXED (all loadings are NA or fixed)\n")
cat("Sigma loading is FIXED at 1.0\n")
cat("So: factor_var should be fixed, sigma_loading is fixed\n")
cat("Free parameters should be: intercept, beta_x1, sigma (3 total)\n\n")

cat("But we're getting Hessian for 4 parameters!\n")
cat("This suggests evaluate_likelihood_cpp doesn't know about parameter fixing\n")
