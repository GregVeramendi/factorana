#!/usr/bin/env Rscript
# Debug: Check if data.frame and matrix give same likelihood

library(factorana)

set.seed(123)
n <- 100

# Create simple data
Y_choice <- sample(1:3, n, replace = TRUE)
intercept <- rep(1, n)
x1 <- rnorm(n)

# Data frame version
dat_df <- data.frame(Y = Y_choice, intercept = intercept, x1 = x1, eval = 1)

# Matrix version
dat_mat <- as.matrix(dat_df)

cat("Data frame columns:", paste(names(dat_df), collapse=", "), "\n")
cat("Matrix columns:", paste(colnames(dat_mat), collapse=", "), "\n")
cat("\n")

# Check if values match
cat("First 5 rows of data.frame:\n")
print(head(dat_df, 5))
cat("\nFirst 5 rows of matrix:\n")
print(head(dat_mat, 5))
cat("\n")

# Define model (no factors)
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

# Initialize with data.frame
fm_cpp_df <- initialize_factor_model_cpp(ms, dat_df, n_quad = 8)

# Initialize with matrix
fm_cpp_mat <- initialize_factor_model_cpp(ms, dat_mat, n_quad = 8)

# Test same parameters on both
test_params <- c(1.0, 0.5, 0.3, -0.5, -0.3)

loglik_df <- evaluate_loglik_only_cpp(fm_cpp_df, test_params)
loglik_mat <- evaluate_loglik_only_cpp(fm_cpp_mat, test_params)

cat("Log-likelihood (data.frame):", loglik_df, "\n")
cat("Log-likelihood (matrix):    ", loglik_mat, "\n")
cat("Difference:                 ", abs(loglik_df - loglik_mat), "\n\n")

if (abs(loglik_df - loglik_mat) < 1e-10) {
  cat("✓ Likelihoods match - no data format issue\n")
} else {
  cat("⚠️  Likelihoods DON'T match - data format issue!\n")
}
