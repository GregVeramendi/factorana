#!/usr/bin/env Rscript
library(factorana)

set.seed(456)
n <- 2000

true_Y1_sigma <- 0.5
true_Y0_sigma <- 0.5
true_T1_sigma <- 0.5
true_T2_sigma <- 0.6
true_T3_sigma <- 0.4

f <- rnorm(n, 0, 1)
U1 <- rnorm(n, 0, 1) * true_Y1_sigma
U0 <- rnorm(n, 0, 1) * true_Y0_sigma
UV <- rnorm(n, 0, 1)
UT1 <- rnorm(n, 0, 1) * true_T1_sigma
UT2 <- rnorm(n, 0, 1) * true_T2_sigma
UT3 <- rnorm(n, 0, 1) * true_T3_sigma

intercept <- rep(1, n)
Z1 <- rnorm(n, 0, 1)
X1 <- rnorm(n, 0, 1)
Q1 <- rnorm(n, 0, 1)
Q2 <- rnorm(n, 0, 1)
Q3 <- rnorm(n, 0, 1)

Y1 <- 2.0 * intercept + 1.0 * X1 + 2.0 * f + U1
Y0 <- 1.0 * intercept + 1.0 * X1 + 1.0 * f + U0
I_latent <- 0.5 * intercept + 1.0 * Z1 + 1.0 * f + UV
D <- as.integer(I_latent > 0)
T1_score <- 1.0 * intercept + 1.0 * Q1 + 1.0 * f + UT1
T2_score <- 0.5 * intercept + 0.8 * Q2 + 1.2 * f + UT2
T3_score <- -0.5 * intercept + 1.1 * Q3 + 0.9 * f + UT3
Y <- ifelse(D == 1, Y1, Y0)

dat <- data.frame(
  intercept = intercept, Y = Y, D = D, Z1 = Z1, X1 = X1,
  Q1 = Q1, Q2 = Q2, Q3 = Q3,
  T1 = T1_score, T2 = T2_score, T3 = T3_score,
  eval_Y1 = D, eval_Y0 = 1L - D
)

fm <- define_factor_model(n_factors = 1, n_types = 1, n_quad = 20)  # More quad points

mc_sel <- define_model_component(
  name = "Selection", data = dat, outcome = "D", factor = fm,
  covariates = c("intercept", "Z1"), model_type = "probit",
  loading_normalization = NA_real_
)
mc_Y1 <- define_model_component(
  name = "Y1", data = dat, outcome = "Y", factor = fm,
  covariates = c("intercept", "X1"), model_type = "linear",
  loading_normalization = NA_real_, evaluation_indicator = "eval_Y1"
)
mc_Y0 <- define_model_component(
  name = "Y0", data = dat, outcome = "Y", factor = fm,
  covariates = c("intercept", "X1"), model_type = "linear",
  loading_normalization = NA_real_, evaluation_indicator = "eval_Y0"
)
mc_T1 <- define_model_component(
  name = "TestScore1", data = dat, outcome = "T1", factor = fm,
  covariates = c("intercept", "Q1"), model_type = "linear",
  loading_normalization = NA_real_
)
mc_T2 <- define_model_component(
  name = "TestScore2", data = dat, outcome = "T2", factor = fm,
  covariates = c("intercept", "Q2"), model_type = "linear",
  loading_normalization = NA_real_
)
mc_T3 <- define_model_component(
  name = "TestScore3", data = dat, outcome = "T3", factor = fm,
  covariates = c("intercept", "Q3"), model_type = "linear",
  loading_normalization = NA_real_
)

ms <- define_model_system(components = list(mc_sel, mc_Y1, mc_Y0, mc_T1, mc_T2, mc_T3), factor = fm)

true_params <- c(
  1.0,  0.5, 1.0, 1.0,  2.0, 1.0, 2.0, 0.5,  1.0, 1.0, 1.0, 0.5,
  1.0, 1.0, 1.0, 0.5,  0.5, 0.8, 1.2, 0.6,  -0.5, 1.1, 0.9, 0.4
)

cat("Testing with different quadrature points:\n\n")

for (nq in c(8, 15, 20)) {
  fm_cpp <- initialize_factor_model_cpp(ms, as.matrix(dat), n_quad = nq)
  result <- evaluate_likelihood_cpp(fm_cpp, true_params, compute_gradient = TRUE, compute_hessian = FALSE)
  
  cat(sprintf("n_quad = %2d: loglik = %10.4f, max|grad| = %10.4f\n", 
              nq, result$loglik, max(abs(result$gradient))))
}
