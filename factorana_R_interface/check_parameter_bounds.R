library(factorana)

set.seed(456)
n <- 100

dat <- data.frame(
  intercept = 1,
  Y = rnorm(n),
  D = sample(0:1, n, replace = TRUE),
  Z1 = rnorm(n),
  X1 = rnorm(n),
  Q1 = rnorm(n),
  Q2 = rnorm(n),
  Q3 = rnorm(n),
  T1 = rnorm(n),
  T2 = rnorm(n),
  T3 = rnorm(n),
  eval_Y1 = sample(0:1, n, replace = TRUE),
  eval_Y0 = sample(0:1, n, replace = TRUE)
)

fm <- define_factor_model(n_factors = 1, n_types = 1, n_quad = 8)

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
fm_cpp <- initialize_factor_model_cpp(ms, as.matrix(dat), n_quad = 8)

param_info <- get_parameter_info_cpp(fm_cpp)

cat("Parameter bounds:\n")
cat(sprintf("Total parameters: %d\n", param_info$n_param_total))
cat(sprintf("Free parameters: %d\n", param_info$n_param_free))
cat(sprintf("Fixed parameters: %d\n\n", param_info$n_param_total - param_info$n_param_free))

cat("Lower bounds:\n")
print(param_info$lower_bounds)

cat("\nUpper bounds:\n")
print(param_info$upper_bounds)

cat("\nFixed parameters (lower == upper):\n")
fixed_idx <- which(abs(param_info$lower_bounds - param_info$upper_bounds) < 1e-10)
cat(sprintf("Indices: %s\n", paste(fixed_idx, collapse = ", ")))
if (length(fixed_idx) > 0) {
  cat(sprintf("Fixed at values: %s\n", paste(param_info$lower_bounds[fixed_idx], collapse = ", ")))
}
