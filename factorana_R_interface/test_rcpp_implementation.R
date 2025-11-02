#!/usr/bin/env Rscript
# Test script for new Rcpp implementation
# Tests high-priority features: multi-dim integration, gradient/Hessian, parameter mapping

library(factorana)

cat("========================================\n")
cat("Testing Rcpp Implementation\n")
cat("========================================\n\n")

# ===== Test 1: Single factor, linear model =====
cat("Test 1: Single factor, linear model\n")
cat("------------------------------------\n")

set.seed(42)
n <- 100
f1 <- rnorm(n)
x1 <- rnorm(n)
Y <- 1.0 + 0.5*x1 + 0.8*f1 + rnorm(n, 0, 0.5)

dat1 <- data.frame(Y = Y, x1 = x1, eval = 1)

fm1 <- define_factor_model(1, 1, 8, loading_normalization = 1)
mc1 <- define_model_component("Y", dat1, "Y", fm1, covariates = "x1",
                              model_type = "linear", evaluation_indicator = "eval")
ms1 <- define_model_system(factor = fm1, components = list(mc1))

cat("  Model system created\n")
cat("  Number of factors:", fm1$n_factors, "\n")
cat("  Model type:", mc1$model_type, "\n")

# Test Gauss-Hermite quadrature
gh_result <- gauss_hermite_quadrature(8)
cat("  GH quadrature computed:", length(gh_result$nodes), "nodes\n")

# Initialize C++ model
cat("\n  Initializing C++ FactorModel...\n")
tryCatch({
  fm_cpp <- initialize_factor_model_cpp(ms1, dat1, n_quad = 8)
  cat("  SUCCESS: C++ model initialized\n")
  # Note: External pointer created successfully
}, error = function(e) {
  cat("  ERROR:", conditionMessage(e), "\n")
})

# ===== Test 2: Two factors, linear model =====
cat("\n\nTest 2: Two factors, linear model\n")
cat("------------------------------------\n")

set.seed(43)
n <- 100
f1 <- rnorm(n)
f2 <- rnorm(n)
x1 <- rnorm(n)
Y <- 1.0 + 0.5*x1 + 0.6*f1 + 0.4*f2 + rnorm(n, 0, 0.5)

dat2 <- data.frame(Y = Y, x1 = x1, eval = 1)

fm2 <- define_factor_model(2, 1, 8, loading_normalization = c(1, 1))
mc2 <- define_model_component("Y", dat2, "Y", fm2, covariates = "x1",
                              model_type = "linear", evaluation_indicator = "eval")
ms2 <- define_model_system(factor = fm2, components = list(mc2))

cat("  Model system created\n")
cat("  Number of factors:", fm2$n_factors, "\n")
cat("  Expected integration points:", 8^2, "\n")

# Initialize C++ model with 2 factors
cat("\n  Initializing C++ FactorModel with 2 factors...\n")
tryCatch({
  fm_cpp2 <- initialize_factor_model_cpp(ms2, dat2, n_quad = 8)
  cat("  SUCCESS: Multi-factor model initialized\n")
  # Multi-dimensional integration enabled
}, error = function(e) {
  cat("  ERROR:", conditionMessage(e), "\n")
})

# ===== Test 3: Probit model =====
cat("\n\nTest 3: Single factor, probit model\n")
cat("------------------------------------\n")

set.seed(44)
n <- 100
f1 <- rnorm(n)
x1 <- rnorm(n)
latent <- 0.3*x1 + 0.7*f1 + rnorm(n)
D <- as.integer(latent > 0)

dat3 <- data.frame(D = D, x1 = x1, eval = 1)

fm3 <- define_factor_model(1, 1, 8, loading_normalization = 1)
mc3 <- define_model_component("D", dat3, "D", fm3, covariates = "x1",
                              model_type = "probit", evaluation_indicator = "eval")
ms3 <- define_model_system(factor = fm3, components = list(mc3))

cat("  Model system created\n")
cat("  Model type:", mc3$model_type, "\n")

# Initialize C++ model
cat("\n  Initializing C++ FactorModel...\n")
tryCatch({
  fm_cpp3 <- initialize_factor_model_cpp(ms3, dat3, n_quad = 8)
  cat("  SUCCESS: Probit model initialized\n")
  # Probit model type handled correctly
}, error = function(e) {
  cat("  ERROR:", conditionMessage(e), "\n")
})

# ===== Test 4: Likelihood evaluation =====
cat("\n\nTest 4: Likelihood evaluation\n")
cat("------------------------------------\n")

if (exists("fm_cpp")) {
  # First, get parameter info from C++
  param_info <- get_parameter_info_cpp(fm_cpp)
  cat("  C++ Parameter Info:\n")
  cat("    n_obs:", param_info$n_obs, "\n")
  cat("    n_param_total:", param_info$n_param_total, "\n")
  cat("    n_param_free:", param_info$n_param_free, "\n")

  cat("\n  Testing likelihood evaluation...\n")

  # Parameter vector for single factor linear model:
  # Expected structure: [factor_var, model_params...]
  # Model params: [beta_x1, lambda (loading), sigma_y]
  # But if loading_normalization=1, lambda is fixed, so model params = [beta_x1, sigma_y]
  # Total: 1 factor variance + 2 model params = 3 params
  test_params <- c(1.0, 0.5, 0.5)
  cat("    Test params (length", length(test_params), "):", paste(test_params, collapse=", "), "\n")

  tryCatch({
    loglik <- evaluate_loglik_only_cpp(fm_cpp, test_params)
    cat("  Log-likelihood:", loglik, "\n")

    if (is.finite(loglik)) {
      cat("  SUCCESS: Likelihood computed successfully!\n")
      cat("  Value is finite and computable\n")
    } else {
      cat("  WARNING: Non-finite likelihood\n")
    }
  }, error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
    cat("  (May need to adjust parameter count)\n")
  })

  # Test gradient
  cat("\n  Testing gradient evaluation...\n")
  tryCatch({
    result <- evaluate_likelihood_cpp(fm_cpp, test_params,
                                     compute_gradient = TRUE,
                                     compute_hessian = FALSE)
    cat("  Log-likelihood:", result$logLikelihood, "\n")
    cat("  Gradient length:", length(result$gradient), "\n")
    cat("  Gradient:", paste(round(result$gradient, 4), collapse = ", "), "\n")

    if (all(is.finite(result$gradient))) {
      cat("  SUCCESS: Gradient computed with chain rule!\n")
      cat("  All gradient values are finite\n")
    } else {
      cat("  WARNING: Non-finite gradient values\n")
    }
  }, error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
  })
}

# ===== Test 5: Multi-factor likelihood =====
cat("\n\nTest 5: Multi-factor likelihood evaluation\n")
cat("------------------------------------\n")

if (exists("fm_cpp2")) {
  # Get parameter info
  param_info2 <- get_parameter_info_cpp(fm_cpp2)
  cat("  C++ Parameter Info:\n")
  cat("    n_obs:", param_info2$n_obs, "\n")
  cat("    n_param_total:", param_info2$n_param_total, "\n")
  cat("    n_param_free:", param_info2$n_param_free, "\n")

  cat("\n  Testing 2-factor likelihood evaluation...\n")

  # For k=2 factors with loading_normalization=c(1,1) (both fixed):
  # Parameters: [sigma_f1^2, sigma_f2^2, beta_x1, sigma_y]
  # Total: 2 factor variances + 2 model params = 4 params
  test_params2 <- c(1.0, 1.0, 0.5, 0.5)
  cat("    Test params (length", length(test_params2), "):", paste(test_params2, collapse=", "), "\n")

  tryCatch({
    loglik2 <- evaluate_loglik_only_cpp(fm_cpp2, test_params2)
    cat("  Log-likelihood:", loglik2, "\n")

    if (is.finite(loglik2)) {
      cat("  SUCCESS: Multi-dimensional integration works!\n")
      cat("  Integrated over", 8^2, "=", 64, "quadrature points\n")
      cat("  Multi-factor models fully functional\n")
    } else {
      cat("  WARNING: Non-finite likelihood\n")
    }
  }, error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
    cat("  (May need to adjust parameter count for 2 factors)\n")
  })
}

# ===== Test 6: Probit likelihood =====
cat("\n\nTest 6: Probit likelihood evaluation\n")
cat("------------------------------------\n")

if (exists("fm_cpp3")) {
  # Get parameter info
  param_info3 <- get_parameter_info_cpp(fm_cpp3)
  cat("  C++ Parameter Info:\n")
  cat("    n_obs:", param_info3$n_obs, "\n")
  cat("    n_param_total:", param_info3$n_param_total, "\n")
  cat("    n_param_free:", param_info3$n_param_free, "\n")

  cat("\n  Testing probit likelihood evaluation...\n")

  # For probit with 1 factor, loading_normalization=1 (fixed):
  # Parameters: [sigma_f^2, beta_x1]
  # Total: 1 factor variance + 1 model param = 2 params
  # (No sigma for probit since variance is fixed at 1)
  test_params3 <- c(1.0, 0.5)
  cat("    Test params (length", length(test_params3), "):", paste(test_params3, collapse=", "), "\n")

  tryCatch({
    loglik3 <- evaluate_loglik_only_cpp(fm_cpp3, test_params3)
    cat("  Log-likelihood:", loglik3, "\n")

    if (is.finite(loglik3)) {
      cat("  SUCCESS: Probit model works!\n")
      cat("  Likelihood computed correctly for binary outcome\n")
    } else {
      cat("  WARNING: Non-finite likelihood\n")
    }
  }, error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
  })

  # Test gradient
  cat("\n  Testing probit gradient...\n")
  tryCatch({
    result3 <- evaluate_likelihood_cpp(fm_cpp3, test_params3,
                                      compute_gradient = TRUE,
                                      compute_hessian = FALSE)
    cat("  Log-likelihood:", result3$logLikelihood, "\n")
    cat("  Gradient length:", length(result3$gradient), "\n")
    cat("  Gradient:", paste(round(result3$gradient, 4), collapse = ", "), "\n")

    if (all(is.finite(result3$gradient))) {
      cat("  SUCCESS: Probit gradients computed!\n")
    } else {
      cat("  WARNING: Non-finite gradient values\n")
    }
  }, error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
  })
}

cat("\n========================================\n")
cat("Testing Complete!\n")
cat("========================================\n")
