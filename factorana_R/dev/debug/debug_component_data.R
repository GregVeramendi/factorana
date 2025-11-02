#!/usr/bin/env Rscript
# Check if model component stores data differently for df vs matrix

library(factorana)

set.seed(123)
n <- 20

# Create simple data
Y_choice <- sample(1:3, n, replace = TRUE)
intercept <- rep(1, n)
x1 <- rnorm(n)

# Data frame version
dat_df <- data.frame(Y = Y_choice, intercept = intercept, x1 = x1, eval = 1)

# Matrix version - convert BEFORE passing to define_model_component
dat_mat <- as.matrix(dat_df)
dat_mat_df <- as.data.frame(dat_mat)  # Convert back to df to pass to model component

cat("=== CHECKING MODEL COMPONENT DATA STORAGE ===\n\n")

# Define model
fm <- define_factor_model(
  n_factors = 1,
  n_types = 1,
  n_quad = 8,
  loading_normalization = 0.0
)

# Component with original data.frame
mc_df <- define_model_component(
  name = "Y",
  data = dat_df,
  outcome = "Y",
  factor = fm,
  covariates = c("intercept", "x1"),
  model_type = "logit",
  num_choices = 3,
  evaluation_indicator = "eval"
)

# Component with matrix-converted data.frame
mc_mat <- define_model_component(
  name = "Y",
  data = dat_mat_df,
  outcome = "Y",
  factor = fm,
  covariates = c("intercept", "x1"),
  model_type = "logit",
  num_choices = 3,
  evaluation_indicator = "eval"
)

cat("Data stored in mc_df (first 5 rows):\n")
print(head(mc_df$data, 5))
cat("\n")

cat("Data stored in mc_mat (first 5 rows):\n")
print(head(mc_mat$data, 5))
cat("\n")

cat("Data types in mc_df:\n")
str(mc_df$data)
cat("\n")

cat("Data types in mc_mat:\n")
str(mc_mat$data)
cat("\n")

# Check if Y column is different type
cat("Class of Y in mc_df:", class(mc_df$data$Y), "\n")
cat("Class of Y in mc_mat:", class(mc_mat$data$Y), "\n")
cat("\n")

# Check exact values
cat("First 10 Y values in mc_df:", mc_df$data$Y[1:10], "\n")
cat("First 10 Y values in mc_mat:", mc_mat$data$Y[1:10], "\n")
cat("Identical?", identical(mc_df$data$Y, mc_mat$data$Y), "\n\n")

cat("First 10 x1 values in mc_df:", mc_df$data$x1[1:10], "\n")
cat("First 10 x1 values in mc_mat:", mc_mat$data$x1[1:10], "\n")
cat("Identical?", identical(mc_df$data$x1, mc_mat$data$x1), "\n\n")

# Now test if using the model component data directly makes a difference
ms_df <- define_model_system(factor = fm, components = list(mc_df))
ms_mat <- define_model_system(factor = fm, components = list(mc_mat))

# Initialize using the actual data stored in the model component
fm_cpp_df <- initialize_factor_model_cpp(ms_df, mc_df$data, n_quad = 8)
fm_cpp_mat <- initialize_factor_model_cpp(ms_mat, mc_mat$data, n_quad = 8)

test_params <- c(1.0, 0.5, 0.3, -0.5, -0.3)

ll_df <- evaluate_loglik_only_cpp(fm_cpp_df, test_params)
ll_mat <- evaluate_loglik_only_cpp(fm_cpp_mat, test_params)

cat("=== LIKELIHOOD TEST USING MODEL COMPONENT DATA ===\n")
cat("Log-likelihood (from mc_df$data):", ll_df, "\n")
cat("Log-likelihood (from mc_mat$data):", ll_mat, "\n")
cat("Difference:", abs(ll_df - ll_mat), "\n")
