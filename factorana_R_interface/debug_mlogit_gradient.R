#!/usr/bin/env Rscript
# Debug: Check if choice 3 gradients are computed

library(factorana)

set.seed(123)
n <- 100

# Create simple 3-choice data
x1 <- rnorm(n)
intercept <- rep(1, n)

# Make sure we have choice 3 observations
Y_choice <- sample(1:3, n, replace = TRUE)
cat("Choice distribution:\n")
print(table(Y_choice))
cat("\n")

dat <- data.frame(Y = Y_choice, intercept = intercept, x1 = x1, eval = 1)

# Define model (no factors)
fm <- define_factor_model(
  n_factors = 1,
  n_types = 1,
  n_quad = 8,
  loading_normalization = 0.0
)

mc <- define_model_component(
  name = "Y",
  data = dat,
  outcome = "Y",
  factor = fm,
  covariates = c("intercept", "x1"),
  model_type = "logit",
  num_choices = 3,
  evaluation_indicator = "eval"
)

ms <- define_model_system(factor = fm, components = list(mc))

# Initialize C++ model
fm_cpp <- initialize_factor_model_cpp(ms, dat, n_quad = 8)

# Test parameters (all non-zero to see clear gradients)
test_params <- c(1.0, 0.5, 0.3, -0.5, -0.3)
cat("Test parameters:", paste(test_params, collapse=", "), "\n\n")

# Compute gradient
result <- evaluate_likelihood_cpp(fm_cpp, test_params,
                                 compute_gradient = TRUE,
                                 compute_hessian = FALSE)

cat("Log-likelihood:", result$logLikelihood, "\n")
cat("Gradients:\n")
param_names <- c("sigma_f^2", "beta0_c2", "beta_x1_c2", "beta0_c3", "beta_x1_c3")
for (i in seq_along(result$gradient)) {
  cat(sprintf("  [%d] %-12s: %12.6f\n", i, param_names[i], result$gradient[i]))
}

# Check if choice 3 gradients are zero
if (abs(result$gradient[4]) < 1e-10 && abs(result$gradient[5]) < 1e-10) {
  cat("\n⚠️  Choice 3 gradients are ZERO - this is the bug!\n")
} else {
  cat("\n✓  Choice 3 gradients are NON-ZERO - working correctly!\n")
}
