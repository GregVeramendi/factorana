#!/usr/bin/env Rscript
# Detailed investigation of data.frame vs matrix issue

library(factorana)

set.seed(123)
n <- 20  # Small sample for detailed inspection

# Create simple data
Y_choice <- sample(1:3, n, replace = TRUE)
intercept <- rep(1, n)
x1 <- rnorm(n)

# Data frame version
dat_df <- data.frame(Y = Y_choice, intercept = intercept, x1 = x1, eval = 1)

# Matrix version
dat_mat <- as.matrix(dat_df)

cat("=== DATA COMPARISON ===\n\n")
cat("Data frame class:", class(dat_df), "\n")
cat("Matrix class:", class(dat_mat), "\n\n")

cat("Data frame structure:\n")
str(dat_df)
cat("\n")

cat("Matrix structure:\n")
str(dat_mat)
cat("\n")

# Check if conversion preserves values
cat("First 10 Y values (data.frame):", dat_df$Y[1:10], "\n")
cat("First 10 Y values (matrix):    ", dat_mat[1:10, "Y"], "\n\n")

cat("First 10 x1 values (data.frame):", dat_df$x1[1:10], "\n")
cat("First 10 x1 values (matrix):    ", dat_mat[1:10, "x1"], "\n\n")

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
cat("=== INITIALIZING WITH DATA.FRAME ===\n")
fm_cpp_df <- initialize_factor_model_cpp(ms, dat_df, n_quad = 8)

# Initialize with matrix
cat("=== INITIALIZING WITH MATRIX ===\n")
fm_cpp_mat <- initialize_factor_model_cpp(ms, dat_mat, n_quad = 8)

# Test same parameters on both
test_params <- c(1.0, 0.5, 0.3, -0.5, -0.3)

cat("\n=== LIKELIHOOD EVALUATION ===\n")
cat("Test parameters:", paste(test_params, collapse=", "), "\n\n")

# Evaluate with data.frame version
result_df <- evaluate_likelihood_cpp(fm_cpp_df, test_params,
                                     compute_gradient = TRUE,
                                     compute_hessian = FALSE)

# Evaluate with matrix version
result_mat <- evaluate_likelihood_cpp(fm_cpp_mat, test_params,
                                      compute_gradient = TRUE,
                                      compute_hessian = FALSE)

cat("Log-likelihood (data.frame):", result_df$logLikelihood, "\n")
cat("Log-likelihood (matrix):    ", result_mat$logLikelihood, "\n")
cat("Difference:                 ", abs(result_df$logLikelihood - result_mat$logLikelihood), "\n\n")

cat("Gradients (data.frame):\n")
for (i in seq_along(result_df$gradient)) {
  cat(sprintf("  [%d] = %.6f\n", i, result_df$gradient[i]))
}

cat("\nGradients (matrix):\n")
for (i in seq_along(result_mat$gradient)) {
  cat(sprintf("  [%d] = %.6f\n", i, result_mat$gradient[i]))
}

cat("\nGradient differences:\n")
for (i in seq_along(result_df$gradient)) {
  diff <- abs(result_df$gradient[i] - result_mat$gradient[i])
  cat(sprintf("  [%d] = %.6f\n", i, diff))
}

# Check individual observation contributions
cat("\n=== PER-OBSERVATION CHECK ===\n")
cat("Checking if specific observations give different results...\n\n")

# Create single-observation datasets
for (i in 1:min(5, n)) {
  dat_single_df <- dat_df[i, , drop = FALSE]
  dat_single_mat <- as.matrix(dat_single_df)

  mc_single <- define_model_component(
    name = "Y",
    data = dat_single_df,
    outcome = "Y",
    factor = fm,
    covariates = c("intercept", "x1"),
    model_type = "logit",
    num_choices = 3,
    evaluation_indicator = "eval"
  )

  ms_single <- define_model_system(factor = fm, components = list(mc_single))

  fm_single_df <- initialize_factor_model_cpp(ms_single, dat_single_df, n_quad = 8)
  fm_single_mat <- initialize_factor_model_cpp(ms_single, dat_single_mat, n_quad = 8)

  ll_df <- evaluate_loglik_only_cpp(fm_single_df, test_params)
  ll_mat <- evaluate_loglik_only_cpp(fm_single_mat, test_params)

  cat(sprintf("Obs %d (Y=%d, x1=%.3f): LL_df=%.6f, LL_mat=%.6f, diff=%.6f\n",
              i, dat_df$Y[i], dat_df$x1[i], ll_df, ll_mat, abs(ll_df - ll_mat)))
}
