#!/usr/bin/env Rscript
# Test if deep copy in estimate_model_rcpp works

library(factorana)

set.seed(123)
n <- 100

Y_choice <- sample(1:3, n, replace = TRUE)
intercept <- rep(1, n)
x1 <- rnorm(n)

dat <- data.frame(Y = Y_choice, intercept = intercept, x1 = x1, eval = 1)

fm <- define_factor_model(n_factors = 1, n_types = 1, n_quad = 8, loading_normalization = 0.0)
mc <- define_model_component(name = "Y", data = dat, outcome = "Y",
                             factor = fm, covariates = c("intercept", "x1"),
                             model_type = "logit", num_choices = 3, evaluation_indicator = "eval")
ms <- define_model_system(factor = fm, components = list(mc))

cat("=== Testing if estimate_model_rcpp uses fresh model_system ===\n\n")

# estimate_model_rcpp uses matrix internally, so test with matrix too
dat_mat <- as.matrix(dat)

# Create first C++ object (use matrix like estimate_model_rcpp does)
fm_cpp_1 <- initialize_factor_model_cpp(ms, dat_mat, n_quad = 8)

# Call estimate_model_rcpp (which should make a deep copy)
result <- estimate_model_rcpp(
  ms, dat,
  init_params = c(1.0, 0.5, 0.3, -0.5, -0.3),
  optimizer = "nlminb",
  parallel = FALSE,
  verbose = FALSE
)

cat("Estimation complete\n")
cat("  Final LL from estimate_model_rcpp:", result$loglik, "\n")

# Check with the first C++ object
ll_check <- evaluate_loglik_only_cpp(fm_cpp_1, result$estimates)
cat("  LL check with fm_cpp_1:", ll_check, "\n\n")

# Create a FRESH C++ object after estimation (use matrix)
fm2 <- define_factor_model(n_factors = 1, n_types = 1, n_quad = 8, loading_normalization = 0.0)
mc2 <- define_model_component(name = "Y", data = dat, outcome = "Y",
                              factor = fm2, covariates = c("intercept", "x1"),
                              model_type = "logit", num_choices = 3, evaluation_indicator = "eval")
ms2 <- define_model_system(factor = fm2, components = list(mc2))
fm_cpp_2 <- initialize_factor_model_cpp(ms2, dat_mat, n_quad = 8)

ll_fresh <- evaluate_loglik_only_cpp(fm_cpp_2, result$estimates)
cat("  LL check with fresh C++ object:", ll_fresh, "\n\n")

if (abs(result$loglik - ll_check) < 1e-6 && abs(result$loglik - ll_fresh) < 1e-6) {
  cat("✓ All likelihoods match - deep copy is working!\n")
} else {
  cat("⚠️  Likelihoods differ:\n")
  cat("   estimate_model_rcpp:", result$loglik, "\n")
  cat("   fm_cpp_1:", ll_check, "\n")
  cat("   fm_cpp_fresh:", ll_fresh, "\n")
}
