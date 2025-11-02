#!/usr/bin/env Rscript
# Check if loading_normalization gets modified

library(factorana)

set.seed(123)
n <- 20
dat_mat <- as.matrix(data.frame(Y = sample(1:3, n, replace = TRUE),
                                intercept = rep(1, n),  x1 = rnorm(n), eval = 1))

fm <- define_factor_model(n_factors = 1, n_types = 1, n_quad = 8, loading_normalization = 0.0)
mc <- define_model_component(name = "Y", data = as.data.frame(dat_mat), outcome = "Y",
                             factor = fm, covariates = c("intercept", "x1"),
                             model_type = "logit", num_choices = 3, evaluation_indicator = "eval")
ms <- define_model_system(factor = fm, components = list(mc))

cat("=== TRACKING loading_normalization ===\n\n")

cat("BEFORE init 1:\n")
cat("  Value:", ms$factor$loading_normalization, "\n")
cat("  Class:", class(ms$factor$loading_normalization), "\n")
cat("  is.na?:", is.na(ms$factor$loading_normalization), "\n\n")

fm_cpp_1 <- initialize_factor_model_cpp(ms, dat_mat, n_quad = 8)

cat("AFTER init 1:\n")
cat("  Value:", ms$factor$loading_normalization, "\n")
cat("  Class:", class(ms$factor$loading_normalization), "\n")
cat("  is.na?:", is.na(ms$factor$loading_normalization), "\n\n")

fm_cpp_2 <- initialize_factor_model_cpp(ms, dat_mat, n_quad = 8)

cat("AFTER init 2:\n")
cat("  Value:", ms$factor$loading_normalization, "\n")
cat("  Class:", class(ms$factor$loading_normalization), "\n")
cat("  is.na?:", is.na(ms$factor$loading_normalization), "\n\n")

test_params <- c(1.0, 0.5, 0.3, -0.5, -0.3)
cat("Likelihood 1:", evaluate_loglik_only_cpp(fm_cpp_1, test_params), "\n")
cat("Likelihood 2:", evaluate_loglik_only_cpp(fm_cpp_2, test_params), "\n")
