# Tests for simulate_factor_model(). The strongest validation is parameter
# recovery: simulate from known parameters, re-estimate, and check the
# estimates land near the truth. We also check the output structure, seed
# reproducibility, and (for mixtures, which are hard to re-estimate because of
# label switching) that the drawn factor distribution matches the specification.

.sim_base <- function(nb = 500, seed = 1) {
  set.seed(seed)
  data.frame(intercept = 1, x = rnorm(nb),
             y1 = rnorm(nb), y2 = rnorm(nb), y3 = rnorm(nb), y4 = rnorm(nb),
             yp = rbinom(nb, 1, 0.5), yo = sample(1:4, nb, TRUE),
             yl = sample(1:3, nb, TRUE), yb = sample(1:2, nb, TRUE))
}

test_that("output structure and seed reproducibility", {
  base <- .sim_base()
  fm <- define_factor_model(n_factors = 1)
  m1 <- define_model_component("m1", base, "y1", fm, covariates = "intercept",
                               model_type = "linear", loading_normalization = 1)
  m2 <- define_model_component("m2", base, "y2", fm, covariates = c("intercept", "x"),
                               model_type = "linear", loading_normalization = NA_real_)
  ms <- define_model_system(components = list(m1, m2), factor = fm)
  truth <- c(factor_var_1 = 1, m1_intercept = 0.2, m1_sigma = 0.5,
             m2_intercept = -0.1, m2_x = 0.3, m2_loading_1 = 0.9, m2_sigma = 0.5)

  s1 <- simulate_factor_model(ms, base, n = 200, params = truth, seed = 42)
  expect_equal(nrow(s1), 200)
  expect_true(all(c("xid", "intercept", "x", "factor_1", "y1", "y2") %in% names(s1)))
  expect_true(all(c("m1_Vobs", "m1_eps", "m2_Vobs") %in% names(s1)))  # detail columns
  # same seed -> identical data
  s2 <- simulate_factor_model(ms, base, n = 200, params = truth, seed = 42)
  expect_equal(s1, s2)
  # detail = FALSE drops diagnostic columns but keeps outcomes
  s3 <- simulate_factor_model(ms, base, n = 50, params = truth, seed = 1, detail = FALSE)
  expect_false("m1_Vobs" %in% names(s3))
  expect_true(all(c("y1", "y2") %in% names(s3)))
})

test_that("recovers parameters for all four model types (1 factor)", {
  skip_on_cran()
  base <- .sim_base()
  fm <- define_factor_model(n_factors = 1)
  comps <- list(
    define_model_component("m1", base, "y1", fm, covariates = "intercept",
                           model_type = "linear", loading_normalization = 1),
    define_model_component("m2", base, "y2", fm, covariates = "intercept",
                           model_type = "linear", loading_normalization = NA_real_),
    define_model_component("m3", base, "y3", fm, covariates = "intercept",
                           model_type = "linear", loading_normalization = NA_real_),
    define_model_component("mp", base, "yp", fm, covariates = c("intercept", "x"),
                           model_type = "probit", loading_normalization = NA_real_),
    define_model_component("mo", base, "yo", fm, covariates = NULL,
                           model_type = "oprobit", num_choices = 4, loading_normalization = NA_real_),
    define_model_component("ml", base, "yl", fm, covariates = "intercept",
                           model_type = "logit", num_choices = 3, loading_normalization = NA_real_))
  ms <- define_model_system(components = comps, factor = fm)
  truth <- c(factor_var_1 = 1,
             m1_intercept = 0.5, m1_sigma = 0.7,
             m2_intercept = -0.2, m2_loading_1 = 0.8, m2_sigma = 0.6,
             m3_intercept = 0.1, m3_loading_1 = 1.2, m3_sigma = 0.65,
             mp_intercept = 0.3, mp_x = 0.5, mp_loading_1 = 0.9,
             mo_loading_1 = 1.0, mo_thresh_1 = -0.8, mo_thresh_2 = 0.9, mo_thresh_3 = 0.9,
             ml_c1_intercept = 0.4, ml_c1_loading_1 = 0.8,
             ml_c2_intercept = -0.3, ml_c2_loading_1 = 1.2)
  sim <- simulate_factor_model(ms, base, n = 6000, params = truth, seed = 42)
  expect_setequal(sort(unique(sim$yo)), 1:4)
  expect_setequal(sort(unique(sim$yl)), 1:3)

  ctrl <- define_estimation_control(n_quad_points = 12, num_cores = 1)
  fit <- estimate_model_rcpp(ms, sim, control = ctrl, optimizer = "nlminb",
                             parallel = FALSE, verbose = FALSE)
  expect_equal(fit$convergence, 0)
  z <- (fit$estimates[names(truth)] - truth) / fit$std_errors[names(truth)]
  expect_true(all(abs(z) < 3.5),
              info = paste("max |z| =", round(max(abs(z)), 2),
                           "at", names(truth)[which.max(abs(z))]))
})

test_that("recovers a correlated two-factor model", {
  skip_on_cran()
  base <- .sim_base()
  fm <- define_factor_model(n_factors = 2, factor_structure = "correlation")
  comps <- list(
    define_model_component("a1", base, "y1", fm, covariates = "intercept",
                           model_type = "linear", loading_normalization = c(1, 0)),
    define_model_component("a2", base, "y2", fm, covariates = "intercept",
                           model_type = "linear", loading_normalization = c(NA_real_, 0)),
    define_model_component("a3", base, "y3", fm, covariates = "intercept",
                           model_type = "linear", loading_normalization = c(0, 1)),
    define_model_component("a4", base, "y4", fm, covariates = "intercept",
                           model_type = "linear", loading_normalization = c(0, NA_real_)))
  ms <- define_model_system(components = comps, factor = fm)
  truth <- c(factor_var_1 = 1, factor_var_2 = 1, factor_corr_1_2 = 0.4,
             a1_intercept = 0, a1_sigma = 0.6, a2_intercept = 0, a2_loading_1 = 0.9, a2_sigma = 0.6,
             a3_intercept = 0, a3_sigma = 0.6, a4_intercept = 0, a4_loading_2 = 1.1, a4_sigma = 0.6)
  sim <- simulate_factor_model(ms, base, n = 6000, params = truth, seed = 5)
  ctrl <- define_estimation_control(n_quad_points = 12, num_cores = 1)
  fit <- estimate_model_rcpp(ms, sim, control = ctrl, optimizer = "nlminb",
                             parallel = FALSE, verbose = FALSE)
  expect_equal(fit$convergence, 0)
  key <- c("factor_var_1", "factor_var_2", "factor_corr_1_2", "a2_loading_1", "a4_loading_2")
  z <- (fit$estimates[key] - truth[key]) / fit$std_errors[key]
  expect_true(all(abs(z) < 3.5))
})

test_that("recovers an SE_linear structural model", {
  skip_on_cran()
  base <- .sim_base()
  fm <- define_factor_model(n_factors = 2, factor_structure = "SE_linear")
  comps <- list(
    define_model_component("d1", base, "y1", fm, covariates = "intercept",
                           model_type = "linear", loading_normalization = c(1, 0)),
    define_model_component("d2", base, "y2", fm, covariates = "intercept",
                           model_type = "linear", loading_normalization = c(NA_real_, 0)),
    fix_coefficient(define_model_component("d3", base, "y3", fm, covariates = "intercept",
                           model_type = "linear", loading_normalization = c(0, 1)), "intercept", 0),
    fix_coefficient(define_model_component("d4", base, "y4", fm, covariates = "intercept",
                           model_type = "linear", loading_normalization = c(0, NA_real_)), "intercept", 0))
  ms <- define_model_system(components = comps, factor = fm)
  truth <- c(factor_var_1 = 1, se_intercept = 0.4, se_linear_1 = 0.6, se_residual_var = 0.5,
             d1_intercept = 0, d1_sigma = 0.6, d2_intercept = 0, d2_loading_1 = 0.9, d2_sigma = 0.6,
             d3_intercept = 0, d3_sigma = 0.6, d4_intercept = 0, d4_loading_2 = 1.1, d4_sigma = 0.6)
  sim <- simulate_factor_model(ms, base, n = 6000, params = truth, seed = 6)
  ctrl <- define_estimation_control(n_quad_points = 12, num_cores = 1)
  fit <- estimate_model_rcpp(ms, sim, control = ctrl, optimizer = "nlminb",
                             parallel = FALSE, verbose = FALSE)
  expect_equal(fit$convergence, 0)
  key <- c("factor_var_1", "se_intercept", "se_linear_1", "se_residual_var")
  z <- (fit$estimates[key] - truth[key]) / fit$std_errors[key]
  expect_true(all(abs(z) < 3.5))
})

test_that("recovers latent-type parameters and draws the right type fractions", {
  skip_on_cran()
  base <- .sim_base()
  fm <- define_factor_model(n_factors = 1, n_types = 2)
  comps <- list(
    define_model_component("m1", base, "y1", fm, covariates = "intercept",
                           model_type = "linear", loading_normalization = 1),
    define_model_component("m2", base, "y2", fm, covariates = "intercept",
                           model_type = "linear", loading_normalization = NA_real_),
    define_model_component("m3", base, "y3", fm, covariates = "intercept",
                           model_type = "linear", loading_normalization = NA_real_, use_types = TRUE),
    define_model_component("m4", base, "y4", fm, covariates = "intercept",
                           model_type = "linear", loading_normalization = NA_real_, use_types = TRUE))
  ms <- define_model_system(components = comps, factor = fm)
  truth <- c(factor_var_1 = 1, typeprob_2_intercept = 0.2, type_2_loading_1 = 0.0,
             m1_intercept = 0, m1_sigma = 0.5, m2_intercept = 0, m2_loading_1 = 0.9, m2_sigma = 0.5,
             m3_intercept = 0, m3_loading_1 = 1.0, m3_sigma = 0.5, m3_type_2_intercept = 1.5,
             m4_intercept = 0, m4_loading_1 = 1.1, m4_sigma = 0.5, m4_type_2_intercept = -1.2)
  sim <- simulate_factor_model(ms, base, n = 10000, params = truth, seed = 8)
  # type fraction matches the type-model intercept
  expect_equal(mean(sim$type == 2), plogis(0.2), tolerance = 0.03)
  ctrl <- define_estimation_control(n_quad_points = 12, num_cores = 1)
  fit <- estimate_model_rcpp(ms, sim, control = ctrl, optimizer = "nlminb",
                             parallel = FALSE, verbose = FALSE)
  expect_equal(fit$convergence, 0)
  key <- c("typeprob_2_intercept", "m3_type_2_intercept", "m4_type_2_intercept")
  z <- (fit$estimates[key] - truth[key]) / fit$std_errors[key]
  expect_true(all(abs(z) < 3.5))
})

test_that("mixture-of-normals factor is drawn from the specified distribution", {
  skip_on_cran()
  base <- .sim_base()
  fm <- define_factor_model(n_factors = 1, n_mixtures = 2)
  comps <- lapply(1:4, function(i) define_model_component(
    paste0("m", i), base, paste0("y", i), fm, covariates = "intercept",
    model_type = "linear", loading_normalization = if (i == 1) 1 else NA_real_))
  ms <- define_model_system(components = comps, factor = fm)
  lw <- log(0.6 / 0.4)  # weights ~ (0.6, 0.4)
  truth <- c(mix1_factor_var_1 = 0.5, mix2_factor_var_1 = 1.2,
             mix1_factor_mean_1 = 1.0, mix1_logweight = lw,
             m1_intercept = 0, m1_sigma = 0.4, m2_intercept = 0, m2_loading_1 = 0.9, m2_sigma = 0.4,
             m3_intercept = 0, m3_loading_1 = 1.1, m3_sigma = 0.4,
             m4_intercept = 0, m4_loading_1 = 1.0, m4_sigma = 0.4)
  sim <- simulate_factor_model(ms, base, n = 12000, params = truth, seed = 9)
  # mixture weight, overall mean (constrained to 0), and per-mixture moments
  expect_equal(mean(sim$mixture == 1), 0.6, tolerance = 0.03)
  expect_equal(mean(sim$factor_1), 0, tolerance = 0.05)
  expect_equal(mean(sim$factor_1[sim$mixture == 1]), 1.0, tolerance = 0.05)
  expect_equal(sd(sim$factor_1[sim$mixture == 1]), sqrt(0.5), tolerance = 0.05)
  expect_equal(mean(sim$factor_1[sim$mixture == 2]), -1.5, tolerance = 0.05)
  expect_equal(sd(sim$factor_1[sim$mixture == 2]), sqrt(1.2), tolerance = 0.06)
})

test_that("errors on missing parameters and unsupported outcomes", {
  base <- .sim_base()
  fm <- define_factor_model(n_factors = 1)
  m1 <- define_model_component("m1", base, "y1", fm, covariates = "intercept",
                               model_type = "linear", loading_normalization = 1)
  ms <- define_model_system(components = list(m1), factor = fm)
  expect_error(simulate_factor_model(ms, base, n = 10, params = c(factor_var_1 = 1)),
               "not found")
  expect_error(simulate_factor_model(ms, base), "params")  # no params for model_system
})
