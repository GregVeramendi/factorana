#!/usr/bin/env Rscript
# Check what gets modified in model_system

library(factorana)

set.seed(123)
n <- 20

Y_choice <- sample(1:3, n, replace = TRUE)
intercept <- rep(1, n)
x1 <- rnorm(n)

dat_mat <- as.matrix(data.frame(Y = Y_choice, intercept = intercept, x1 = x1, eval = 1))

fm <- define_factor_model(n_factors = 1, n_types = 1, n_quad = 8, loading_normalization = 0.0)
mc <- define_model_component(name = "Y", data = as.data.frame(dat_mat), outcome = "Y",
                             factor = fm, covariates = c("intercept", "x1"),
                             model_type = "logit", num_choices = 3, evaluation_indicator = "eval")
ms <- define_model_system(factor = fm, components = list(mc))

cat("=== BEFORE first initialization ===\n")
cat("Y column class:", class(ms$components[[1]]$data$Y), "\n")
cat("First 10 Y values:", ms$components[[1]]$data$Y[1:10], "\n")
cat("Y storage mode:", storage.mode(ms$components[[1]]$data$Y), "\n\n")

fm_cpp_1 <- initialize_factor_model_cpp(ms, dat_mat, n_quad = 8)

cat("=== AFTER first initialization ===\n")
cat("Y column class:", class(ms$components[[1]]$data$Y), "\n")
cat("First 10 Y values:", ms$components[[1]]$data$Y[1:10], "\n")
cat("Y storage mode:", storage.mode(ms$components[[1]]$data$Y), "\n\n")

fm_cpp_2 <- initialize_factor_model_cpp(ms, dat_mat, n_quad = 8)

cat("=== AFTER second initialization ===\n")
cat("Y column class:", class(ms$components[[1]]$data$Y), "\n")
cat("First 10 Y values:", ms$components[[1]]$data$Y[1:10], "\n")
cat("Y storage mode:", storage.mode(ms$components[[1]]$data$Y), "\n\n")

test_params <- c(1.0, 0.5, 0.3, -0.5, -0.3)
cat("Likelihood from object 1:", evaluate_loglik_only_cpp(fm_cpp_1, test_params), "\n")
cat("Likelihood from object 2:", evaluate_loglik_only_cpp(fm_cpp_2, test_params), "\n")
