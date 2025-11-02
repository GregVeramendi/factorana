#!/usr/bin/env Rscript
# Check if Hessian storage order might be the issue

library(factorana)

set.seed(456)
n <- 20
dat <- data.frame(intercept = 1, Y = rnorm(n, 2, 0.5), x1 = rnorm(n), eval = 1)
fm <- define_factor_model(n_factors = 1, n_types = 1, n_quad = 8)
mc_Y <- define_model_component(
  name = "Y", data = dat, outcome = "Y", factor = fm,
  covariates = c("intercept", "x1"), model_type = "linear",
  loading_normalization = 1.0, evaluation_indicator = "eval"
)
ms <- define_model_system(components = list(mc_Y), factor = fm)
fm_cpp <- initialize_factor_model_cpp(ms, as.matrix(dat), n_quad = 8)

params <- c(1.0, 2.0, 0.5, 0.5)

result <- evaluate_likelihood_cpp(fm_cpp, params, TRUE, TRUE)

cat("Hessian as returned (upper triangle, length 10):\n")
print(round(result$hessian, 4))

cat("\n\nConverting to 4x4 matrix (assuming ROW-MAJOR upper triangle):\n")
n_p <- 4
hess_row <- matrix(0, n_p, n_p)
idx <- 1
for (i in 1:n_p) {
  for (j in i:n_p) {
    hess_row[i, j] <- result$hessian[idx]
    hess_row[j, i] <- result$hessian[idx]
    idx <- idx + 1
  }
}
print(round(hess_row, 4))
cat("Diagonal:", round(diag(hess_row), 4), "\n")

cat("\n\nWhat if it's COLUMN-MAJOR upper triangle?\n")
hess_col <- matrix(0, n_p, n_p)
idx <- 1
for (j in 1:n_p) {
  for (i in 1:j) {
    hess_col[i, j] <- result$hessian[idx]
    hess_col[j, i] <- result$hessian[idx]
    idx <- idx + 1
  }
}
print(round(hess_col, 4))
cat("Diagonal:", round(diag(hess_col), 4), "\n")
