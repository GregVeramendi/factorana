#!/usr/bin/env Rscript
# Test what happens when Rcpp converts data.frame to matrix

library(factorana)

set.seed(123)
n <- 10

Y_int <- sample(1:3, n, replace = TRUE)
intercept <- rep(1, n)
x1 <- rnorm(n)

# Data frame with integer Y
dat_df_int <- data.frame(Y = as.integer(Y_int), intercept = intercept, x1 = x1, eval = 1)

# Data frame with numeric Y
dat_df_num <- data.frame(Y = as.numeric(Y_int), intercept = intercept, x1 = x1, eval = 1)

cat("=== TESTING R's as.matrix() CONVERSION ===\n\n")

cat("Original data.frame (integer Y):\n")
str(dat_df_int)
cat("\n")

mat_from_int <- as.matrix(dat_df_int)
cat("Result of as.matrix() on integer df:\n")
str(mat_from_int)
cat("Y column class:", class(mat_from_int[, "Y"]), "\n")
cat("Y values:", mat_from_int[, "Y"], "\n\n")

cat("Original data.frame (numeric Y):\n")
str(dat_df_num)
cat("\n")

mat_from_num <- as.matrix(dat_df_num)
cat("Result of as.matrix() on numeric df:\n")
str(mat_from_num)
cat("Y column class:", class(mat_from_num[, "Y"]), "\n")
cat("Y values:", mat_from_num[, "Y"], "\n\n")

# Direct matrix creation from integer df
cat("Matrix created directly from integer df:\n")
mat_direct <- as.matrix(dat_df_int)
str(mat_direct)

# Key question: Are they identical?
cat("\n=== COMPARISON ===\n")
cat("mat_from_int == mat_from_num:", identical(mat_from_int, mat_from_num), "\n")
cat("Y columns equal:", all(mat_from_int[, "Y"] == mat_from_num[, "Y"]), "\n")
cat("Y columns identical:", identical(mat_from_int[, "Y"], mat_from_num[, "Y"]), "\n")
