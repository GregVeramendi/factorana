#!/usr/bin/env Rscript
# Test if using fresh model_system for each init fixes the issue

library(factorana)

set.seed(123)
n <- 100

Y_choice <- sample(1:3, n, replace = TRUE)
intercept <- rep(1, n)
x1 <- rnorm(n)

dat_mat <- as.matrix(data.frame(Y = Y_choice, intercept = intercept, x1 = x1, eval = 1))

test_params <- c(1.0, 0.5, 0.3, -0.5, -0.3)

cat("=== TEST 1: Reuse same model_system ===\n")

fm <- define_factor_model(n_factors = 1, n_types = 1, n_quad = 8, loading_normalization = 0.0)
mc <- define_model_component(name = "Y", data = as.data.frame(dat_mat), outcome = "Y",
                             factor = fm, covariates = c("intercept", "x1"),
                             model_type = "logit", num_choices = 3, evaluation_indicator = "eval")
ms <- define_model_system(factor = fm, components = list(mc))

fm_cpp_1 <- initialize_factor_model_cpp(ms, dat_mat, n_quad = 8)
fm_cpp_2 <- initialize_factor_model_cpp(ms, dat_mat, n_quad = 8)
fm_cpp_3 <- initialize_factor_model_cpp(ms, dat_mat, n_quad = 8)

cat("Object 1:", evaluate_loglik_only_cpp(fm_cpp_1, test_params), "\n")
cat("Object 2:", evaluate_loglik_only_cpp(fm_cpp_2, test_params), "\n")
cat("Object 3:", evaluate_loglik_only_cpp(fm_cpp_3, test_params), "\n\n")

cat("=== TEST 2: Create fresh model_system for each ===\n")

fm1 <- define_factor_model(n_factors = 1, n_types = 1, n_quad = 8, loading_normalization = 0.0)
mc1 <- define_model_component(name = "Y", data = as.data.frame(dat_mat), outcome = "Y",
                              factor = fm1, covariates = c("intercept", "x1"),
                              model_type = "logit", num_choices = 3, evaluation_indicator = "eval")
ms1 <- define_model_system(factor = fm1, components = list(mc1))

fm2 <- define_factor_model(n_factors = 1, n_types = 1, n_quad = 8, loading_normalization = 0.0)
mc2 <- define_model_component(name = "Y", data = as.data.frame(dat_mat), outcome = "Y",
                              factor = fm2, covariates = c("intercept", "x1"),
                              model_type = "logit", num_choices = 3, evaluation_indicator = "eval")
ms2 <- define_model_system(factor = fm2, components = list(mc2))

fm3 <- define_factor_model(n_factors = 1, n_types = 1, n_quad = 8, loading_normalization = 0.0)
mc3 <- define_model_component(name = "Y", data = as.data.frame(dat_mat), outcome = "Y",
                              factor = fm3, covariates = c("intercept", "x1"),
                              model_type = "logit", num_choices = 3, evaluation_indicator = "eval")
ms3 <- define_model_system(factor = fm3, components = list(mc3))

fm_cpp_fresh1 <- initialize_factor_model_cpp(ms1, dat_mat, n_quad = 8)
fm_cpp_fresh2 <- initialize_factor_model_cpp(ms2, dat_mat, n_quad = 8)
fm_cpp_fresh3 <- initialize_factor_model_cpp(ms3, dat_mat, n_quad = 8)

cat("Fresh 1:", evaluate_loglik_only_cpp(fm_cpp_fresh1, test_params), "\n")
cat("Fresh 2:", evaluate_loglik_only_cpp(fm_cpp_fresh2, test_params), "\n")
cat("Fresh 3:", evaluate_loglik_only_cpp(fm_cpp_fresh3, test_params), "\n")
