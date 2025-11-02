#!/usr/bin/env Rscript
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

# Test with POSITIVE sigma
params_pos <- c(1.0, 2.0, 0.5, 0.5)

# Test with NEGATIVE sigma
params_neg <- c(1.0, 2.0, 0.5, -0.5)

cat("Testing if fabs() is handled correctly:\n\n")

result_pos <- evaluate_likelihood_cpp(fm_cpp, params_pos, TRUE, FALSE)
result_neg <- evaluate_likelihood_cpp(fm_cpp, params_neg, TRUE, FALSE)

cat("With sigma = +0.5:\n")
cat(sprintf("  loglik = %.6f, grad[sigma] = %.6f\n", result_pos$loglik, result_pos$gradient[4]))

cat("With sigma = -0.5:\n")
cat(sprintf("  loglik = %.6f, grad[sigma] = %.6f\n", result_neg$loglik, result_neg$gradient[4]))

cat("\nExpected: same loglik, opposite sign gradients\n")
cat(sprintf("Loglik diff: %.10f\n", abs(result_pos$loglik - result_neg$loglik)))
cat(sprintf("Grad sum: %.10f (should be ~0)\n", result_pos$gradient[4] + result_neg$gradient[4]))
