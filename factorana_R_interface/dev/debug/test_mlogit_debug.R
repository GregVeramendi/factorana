#!/usr/bin/env Rscript
library(factorana)

set.seed(42)
n <- 50
f1 <- rnorm(n)
x1 <- rnorm(n)
intercept <- rep(1, n)

U1 <- 0.5*intercept + 0.3*x1 + 0.4*f1 + rnorm(n)
U2 <- -0.3*intercept + 0.5*x1 + 0.6*f1 + rnorm(n)

Y_choice <- rep(1, n)
Y_choice[U1 > 0 & U1 > U2] <- 2
Y_choice[U2 > 0 & U2 > U1] <- 3

dat <- data.frame(Y = Y_choice, intercept = intercept, x1 = x1, eval = 1)

fm <- define_factor_model(
  n_factors = 1,
  n_types = 1,
  n_quad = 8,
  loading_normalization = NA_real_
)

mc <- define_model_component(
  name = "Y",
  data = dat,
  outcome = "Y",
  factor = fm,
  covariates = c("intercept", "x1"),
  model_type = "logit",
  evaluation_indicator = "eval",
  num_choices = 3
)

ms <- define_model_system(factor = fm, components = list(mc))
fm_cpp <- initialize_factor_model_cpp(ms, dat, n_quad = 8)
param_info <- get_parameter_info_cpp(fm_cpp)

cat("Expected parameters:\n")
cat("  Factor variance: 1\n")
cat("  Choice 1 (vs ref): intercept + x1 + loading = 3\n")
cat("  Choice 2 (vs ref): intercept + x1 + loading = 3\n")
cat("  Total: 1 + 3 + 3 = 7\n\n")

cat("Actual from C++:\n")
cat("  n_param_free:", param_info$n_param_free, "\n")
cat("  n_param_total:", param_info$n_param_total, "\n")
