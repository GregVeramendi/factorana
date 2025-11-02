#!/usr/bin/env Rscript
# Test if multiple C++ objects give consistent results

library(factorana)

set.seed(123)
n <- 100

Y_choice <- sample(1:3, n, replace = TRUE)
intercept <- rep(1, n)
x1 <- rnorm(n)

dat_df <- data.frame(Y = Y_choice, intercept = intercept, x1 = x1, eval = 1)
dat_mat <- as.matrix(dat_df)

fm <- define_factor_model(
  n_factors = 1,
  n_types = 1,
  n_quad = 8,
  loading_normalization = 0.0
)

mc <- define_model_component(
  name = "Y",
  data = dat_df,
  outcome = "Y",
  factor = fm,
  covariates = c("intercept", "x1"),
  model_type = "logit",
  num_choices = 3,
  evaluation_indicator = "eval"
)

ms <- define_model_system(factor = fm, components = list(mc))

test_params <- c(1.0, 0.5, 0.3, -0.5, -0.3)

cat("=== CONSISTENCY TEST ===\n\n")

# Create 3 C++ objects with the same matrix data
fm_cpp_1 <- initialize_factor_model_cpp(ms, dat_mat, n_quad = 8)
fm_cpp_2 <- initialize_factor_model_cpp(ms, dat_mat, n_quad = 8)
fm_cpp_3 <- initialize_factor_model_cpp(ms, dat_mat, n_quad = 8)

ll_1 <- evaluate_loglik_only_cpp(fm_cpp_1, test_params)
ll_2 <- evaluate_loglik_only_cpp(fm_cpp_2, test_params)
ll_3 <- evaluate_loglik_only_cpp(fm_cpp_3, test_params)

cat("C++ object 1:", ll_1, "\n")
cat("C++ object 2:", ll_2, "\n")
cat("C++ object 3:", ll_3, "\n\n")

if (identical(ll_1, ll_2) && identical(ll_2, ll_3)) {
  cat("✓ All C++ objects give identical results\n")
} else {
  cat("⚠️  C++ objects give DIFFERENT results - serious bug!\n")
}

# Now test after evaluating multiple times
cat("\n=== TESTING MULTIPLE EVALUATIONS ===\n\n")

ll_1a <- evaluate_loglik_only_cpp(fm_cpp_1, test_params)
ll_1b <- evaluate_loglik_only_cpp(fm_cpp_1, test_params)
ll_1c <- evaluate_loglik_only_cpp(fm_cpp_1, test_params)

cat("Same object, eval 1:", ll_1a, "\n")
cat("Same object, eval 2:", ll_1b, "\n")
cat("Same object, eval 3:", ll_1c, "\n\n")

if (identical(ll_1a, ll_1b) && identical(ll_1b, ll_1c)) {
  cat("✓ Same C++ object gives consistent results across calls\n")
} else {
  cat("⚠️  Same C++ object gives DIFFERENT results - state corruption!\n")
}

# Test with gradient evaluation
cat("\n=== TESTING WITH GRADIENT COMPUTATION ===\n\n")

result_grad <- evaluate_likelihood_cpp(fm_cpp_1, test_params,
                                       compute_gradient = TRUE,
                                       compute_hessian = FALSE)
ll_after_grad <- evaluate_loglik_only_cpp(fm_cpp_1, test_params)

cat("Likelihood from gradient call:", result_grad$logLikelihood, "\n")
cat("Likelihood after gradient call:", ll_after_grad, "\n\n")

if (identical(result_grad$logLikelihood, ll_after_grad)) {
  cat("✓ Gradient computation doesn't affect state\n")
} else {
  cat("⚠️  Gradient computation CHANGES state!\n")
}
