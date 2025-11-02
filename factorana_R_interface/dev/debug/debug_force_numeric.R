#!/usr/bin/env Rscript
# Test if forcing Y to be numeric fixes the issue

library(factorana)

set.seed(123)
n <- 100

# Create data
Y_choice <- sample(1:3, n, replace = TRUE)
intercept <- rep(1, n)
x1 <- rnorm(n)

# Three versions of data
dat_int <- data.frame(Y = as.integer(Y_choice), intercept = intercept, x1 = x1, eval = 1)
dat_num <- data.frame(Y = as.numeric(Y_choice), intercept = intercept, x1 = x1, eval = 1)
dat_mat <- as.matrix(dat_int)

cat("=== DATA TYPE COMPARISON ===\n")
cat("Y type (integer df):", class(dat_int$Y), "\n")
cat("Y type (numeric df):", class(dat_num$Y), "\n")
cat("Y type (matrix):    ", class(dat_mat[, "Y"]), "\n\n")

# Define model
fm <- define_factor_model(
  n_factors = 1,
  n_types = 1,
  n_quad = 8,
  loading_normalization = 0.0
)

mc_int <- define_model_component(
  name = "Y",
  data = dat_int,
  outcome = "Y",
  factor = fm,
  covariates = c("intercept", "x1"),
  model_type = "logit",
  num_choices = 3,
  evaluation_indicator = "eval"
)

mc_num <- define_model_component(
  name = "Y",
  data = dat_num,
  outcome = "Y",
  factor = fm,
  covariates = c("intercept", "x1"),
  model_type = "logit",
  num_choices = 3,
  evaluation_indicator = "eval"
)

ms_int <- define_model_system(factor = fm, components = list(mc_int))
ms_num <- define_model_system(factor = fm, components = list(mc_num))

# Initialize
fm_cpp_int <- initialize_factor_model_cpp(ms_int, dat_int, n_quad = 8)
fm_cpp_num <- initialize_factor_model_cpp(ms_num, dat_num, n_quad = 8)
fm_cpp_mat <- initialize_factor_model_cpp(ms_int, dat_mat, n_quad = 8)

test_params <- c(1.0, 0.5, 0.3, -0.5, -0.3)

ll_int <- evaluate_loglik_only_cpp(fm_cpp_int, test_params)
ll_num <- evaluate_loglik_only_cpp(fm_cpp_num, test_params)
ll_mat <- evaluate_loglik_only_cpp(fm_cpp_mat, test_params)

cat("=== LIKELIHOOD COMPARISON ===\n")
cat("Log-likelihood (integer Y)::", ll_int, "\n")
cat("Log-likelihood (numeric Y)::", ll_num, "\n")
cat("Log-likelihood (matrix):   ", ll_mat, "\n\n")

cat("Difference (int vs num):", abs(ll_int - ll_num), "\n")
cat("Difference (int vs mat):", abs(ll_int - ll_mat), "\n")
cat("Difference (num vs mat):", abs(ll_num - ll_mat), "\n\n")

if (abs(ll_num - ll_mat) < 1e-10) {
  cat("✓ Using numeric Y in data.frame matches matrix!\n")
  cat("  => The issue is integer vs numeric type conversion\n")
} else {
  cat("⚠️ Still different - the issue is elsewhere\n")
}
