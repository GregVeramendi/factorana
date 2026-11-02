# Tests for endogenous regressors in simulate_factor_model(): a covariate in one
# component is the *simulated* outcome of another component (ported from the legacy
# TModel::SetEndogenousRegs / Sim mechanism). A linear/probit producer feeds its
# continuous/binary value; a discrete (logit/oprobit) producer feeds one-hot
# <outcome><choice> dummies. Producers are simulated before consumers via a
# topological sort, so the listing order does not matter.
#
# Recovery on these systems has two requirements. First, a consumer that loads on
# its producer's factor needs at least four indicators on that factor. Three linear
# indicators give 6 moments for 7 parameters, so the coefficient is not identified.
# Second, the quadrature must resolve the factor posterior. With precise indicators,
# n_quad = 12 put one sigma 4.1 SE from truth, so these fits use n_quad = 48.

.endog_base <- function(nb = 600, seed = 1) {
  set.seed(seed)
  catv <- sample(1:3, nb, TRUE)
  # cat1/cat2/cat3 carry realistic one-hot values (as a real data set would) so
  # the component-definition collinearity check passes; simulation overwrites them.
  data.frame(intercept = 1, x = rnorm(nb),
             yA = rnorm(nb), yB = rnorm(nb),
             m1 = rnorm(nb), m2 = rnorm(nb), m3 = rnorm(nb),
             cat = catv,
             cat1 = as.integer(catv == 1L),
             cat2 = as.integer(catv == 2L),
             cat3 = as.integer(catv == 3L),
             pbin = rbinom(nb, 1, 0.5))
}

# Simulate from `truth`, re-estimate on the simulated data, and check that every
# parameter lies within 3.5 SE of truth (the bound used in test-simulate-factor-model.R).
.endog_recover <- function(ms, base, truth, n_quad = 48) {
  sim <- simulate_factor_model(ms, base, n = nrow(base), params = truth, seed = 42)
  ctrl <- define_estimation_control(n_quad_points = n_quad, num_cores = 1)
  fit <- estimate_model_rcpp(ms, sim, control = ctrl, optimizer = "nlminb",
                             parallel = FALSE, verbose = FALSE)
  expect_equal(fit$convergence, 0)
  expect_true(all(names(truth) %in% names(fit$estimates)))
  se <- fit$std_errors[names(truth)]
  expect_true(all(is.finite(se) & se < 5))
  z <- (fit$estimates[names(truth)] - truth) / se
  expect_true(all(abs(z) < 3.5),
              info = paste("max |z| =", round(max(abs(z)), 2),
                           "at", names(truth)[which.max(abs(z))]))
}

test_that("linear producer feeds its simulated outcome into a consumer's predictor", {
  base <- .endog_base()
  fm <- define_factor_model(n_factors = 1)
  mA <- define_model_component("A", base, "yA", fm, covariates = "intercept",
                               model_type = "linear", loading_normalization = 1)
  mB <- define_model_component("B", base, "yB", fm, covariates = c("intercept", "yA"),
                               model_type = "linear", loading_normalization = NA_real_,
                               endogenous = "yA")
  ms <- define_model_system(components = list(mA, mB), factor = fm)
  truth <- c(factor_var_1 = 1, A_intercept = 0.3, A_sigma = 0.6,
             B_intercept = -0.2, B_yA = 0.8, B_loading_1 = 1.1, B_sigma = 0.5)
  sim <- simulate_factor_model(ms, base, n = nrow(base), params = truth, seed = 42)

  # B's linear predictor must be built from the SIMULATED yA, not the resampled one
  pred <- truth["B_intercept"] + truth["B_yA"] * sim$yA +
    truth["B_loading_1"] * sim$factor_1
  expect_equal(unname(pred), sim$B_Vobs, tolerance = 1e-10)
  # the returned yA is A's simulated outcome, not the resampled data value
  expect_equal(sim$yA, sim$A_Vobs + sim$A_eps)
  expect_false(isTRUE(all.equal(sim$yA, base$yA[sim$xid])))
})

test_that("recovers the endogenous coefficient when the consumer skips the factor", {
  skip_on_cran()
  base <- .endog_base(nb = 6000, seed = 11)
  fm <- define_factor_model(n_factors = 1)
  # factor measured by A, m1, m2; consumer B depends on yA but NOT on the factor
  mA <- define_model_component("A", base, "yA", fm, covariates = "intercept",
                               model_type = "linear", loading_normalization = 1)
  m1 <- define_model_component("M1", base, "m1", fm, covariates = "intercept",
                               model_type = "linear", loading_normalization = NA_real_)
  m2 <- define_model_component("M2", base, "m2", fm, covariates = "intercept",
                               model_type = "linear", loading_normalization = NA_real_)
  mB <- define_model_component("B", base, "yB", fm, covariates = c("intercept", "yA"),
                               model_type = "linear", loading_normalization = 0,
                               endogenous = "yA")
  ms <- define_model_system(components = list(mA, m1, m2, mB), factor = fm)
  truth <- c(factor_var_1 = 1,
             A_intercept = 0.3, A_sigma = 0.6,
             M1_intercept = 0.1, M1_loading_1 = 0.9, M1_sigma = 0.7,
             M2_intercept = 0.0, M2_loading_1 = 1.2, M2_sigma = 0.5,
             B_intercept = -0.2, B_yA = 0.8, B_sigma = 0.5)
  .endog_recover(ms, base, truth)
})

test_that("recovers the endogenous coefficient when the consumer shares the factor", {
  skip_on_cran()
  base <- .endog_base(nb = 6000, seed = 12)
  fm <- define_factor_model(n_factors = 1)
  mA <- define_model_component("A", base, "yA", fm, covariates = "intercept",
                               model_type = "linear", loading_normalization = 1)
  mB <- define_model_component("B", base, "yB", fm, covariates = c("intercept", "yA"),
                               model_type = "linear", loading_normalization = NA_real_,
                               endogenous = "yA")
  m1 <- define_model_component("M1", base, "m1", fm, covariates = "intercept",
                               model_type = "linear", loading_normalization = NA_real_)
  m2 <- define_model_component("M2", base, "m2", fm, covariates = "intercept",
                               model_type = "linear", loading_normalization = NA_real_)
  ms <- define_model_system(components = list(mA, mB, m1, m2), factor = fm)
  truth <- c(factor_var_1 = 1,
             A_intercept = 0.3, A_sigma = 0.7,
             B_intercept = -0.2, B_yA = 0.8, B_loading_1 = 0.9, B_sigma = 0.8,
             M1_intercept = 0.1, M1_loading_1 = 0.9, M1_sigma = 0.8,
             M2_intercept = 0, M2_loading_1 = 1.1, M2_sigma = 0.8)
  .endog_recover(ms, base, truth)
})

test_that("recovers the coefficient on a probit producer's outcome", {
  skip_on_cran()
  base <- .endog_base(nb = 6000, seed = 13)
  fm <- define_factor_model(n_factors = 1)
  m1 <- define_model_component("M1", base, "m1", fm, covariates = "intercept",
                               model_type = "linear", loading_normalization = 1)
  m2 <- define_model_component("M2", base, "m2", fm, covariates = "intercept",
                               model_type = "linear", loading_normalization = NA_real_)
  mP <- define_model_component("P", base, "pbin", fm, covariates = "intercept",
                               model_type = "probit", loading_normalization = NA_real_)
  mC <- define_model_component("C", base, "yB", fm, covariates = c("intercept", "pbin"),
                               model_type = "linear", loading_normalization = NA_real_,
                               endogenous = "pbin")
  ms <- define_model_system(components = list(m1, m2, mP, mC), factor = fm)
  truth <- c(factor_var_1 = 1,
             M1_intercept = 0.1, M1_sigma = 0.8,
             M2_intercept = 0, M2_loading_1 = 1.1, M2_sigma = 0.8,
             P_intercept = 0.2, P_loading_1 = 0.9,
             C_intercept = 0.1, C_pbin = 1.3, C_loading_1 = 0.7, C_sigma = 0.8)
  .endog_recover(ms, base, truth)
})

test_that("recovers the dummy coefficients on a discrete producer's outcome", {
  skip_on_cran()
  base <- .endog_base(nb = 6000, seed = 14)
  fm <- define_factor_model(n_factors = 1)
  m1 <- define_model_component("M1", base, "m1", fm, covariates = "intercept",
                               model_type = "linear", loading_normalization = 1)
  m2 <- define_model_component("M2", base, "m2", fm, covariates = "intercept",
                               model_type = "linear", loading_normalization = NA_real_)
  mP <- define_model_component("P", base, "cat", fm, covariates = NULL,
                               model_type = "oprobit", num_choices = 3,
                               loading_normalization = NA_real_)
  mC <- define_model_component("C", base, "yB", fm,
                               covariates = c("intercept", "cat2", "cat3"),
                               model_type = "linear", loading_normalization = NA_real_,
                               endogenous = "cat")
  ms <- define_model_system(components = list(m1, m2, mP, mC), factor = fm)
  truth <- c(factor_var_1 = 1,
             M1_intercept = 0.1, M1_sigma = 0.8,
             M2_intercept = 0, M2_loading_1 = 1.1, M2_sigma = 0.8,
             P_loading_1 = 0.9, P_thresh_1 = -0.5, P_thresh_2 = 1.0,
             C_intercept = 0.1, C_cat2 = 0.7, C_cat3 = -0.4,
             C_loading_1 = 0.7, C_sigma = 0.8)
  .endog_recover(ms, base, truth)
})

test_that("a missing producer outcome propagates to its consumer", {
  base <- .endog_base()
  base$ev <- rep(c(0L, 1L), length.out = nrow(base))
  fm <- define_factor_model(n_factors = 1)
  mA <- define_model_component("A", base, "yA", fm, covariates = "intercept",
                               evaluation_indicator = "ev", model_type = "linear",
                               loading_normalization = 1)
  mB <- define_model_component("B", base, "yB", fm, covariates = c("intercept", "yA"),
                               model_type = "linear", loading_normalization = NA_real_,
                               endogenous = "yA")
  ms <- define_model_system(components = list(mA, mB), factor = fm)
  truth <- c(factor_var_1 = 1, A_intercept = 0.3, A_sigma = 0.6,
             B_intercept = -0.2, B_yA = 0.8, B_loading_1 = 1.1, B_sigma = 0.5)
  sim <- simulate_factor_model(ms, base, n = nrow(base), params = truth, seed = 3)
  expect_true(any(sim$A_miss == 1L))
  expect_true(all(is.na(sim$yA[sim$A_miss == 1L])))
  expect_true(all(is.na(sim$yB[sim$A_miss == 1L])))
  expect_false(anyNA(sim$yB[sim$A_miss == 0L]))
})

test_that("discrete producer feeds one-hot dummies (legacy parity)", {
  base <- .endog_base()
  fm <- define_factor_model(n_factors = 1)
  # oprobit producer 'cat' with 3 categories; consumer uses dummies cat2, cat3
  mP <- define_model_component("P", base, "cat", fm, covariates = NULL,
                               model_type = "oprobit", num_choices = 3,
                               loading_normalization = 1)
  mC <- define_model_component("C", base, "yB", fm,
                               covariates = c("intercept", "cat2", "cat3"),
                               model_type = "linear", loading_normalization = NA_real_,
                               endogenous = "cat")
  ms <- define_model_system(components = list(mP, mC), factor = fm)
  truth <- c(factor_var_1 = 1,
             P_loading_1 = 1.0, P_thresh_1 = -0.5, P_thresh_2 = 1.0,
             C_intercept = 0.1, C_cat2 = 0.7, C_cat3 = -0.4,
             C_loading_1 = 0.9, C_sigma = 0.5)
  sim <- simulate_factor_model(ms, base, n = nrow(base), params = truth, seed = 7)

  expect_true(all(c("cat2", "cat3") %in% names(sim)))
  expect_setequal(unique(sim$cat2), c(0L, 1L))
  expect_equal(sim$cat2, as.integer(sim$cat == 2))
  expect_equal(sim$cat3, as.integer(sim$cat == 3))
  # the dummies (not the resampled real cat2/cat3) drive the consumer predictor
  pred <- truth["C_intercept"] + truth["C_cat2"] * sim$cat2 +
    truth["C_cat3"] * sim$cat3 + truth["C_loading_1"] * sim$factor_1
  expect_equal(unname(pred), sim$C_Vobs, tolerance = 1e-10)
})

test_that("multinomial logit producer feeds one-hot dummies", {
  base <- .endog_base()
  fm <- define_factor_model(n_factors = 1)
  mP <- define_model_component("P", base, "cat", fm, covariates = "intercept",
                               model_type = "logit", num_choices = 3,
                               loading_normalization = NA_real_)
  mC <- define_model_component("C", base, "yB", fm,
                               covariates = c("intercept", "cat2", "cat3"),
                               model_type = "linear", loading_normalization = 0,
                               endogenous = "cat")
  ms <- define_model_system(components = list(mP, mC), factor = fm)
  truth <- c(factor_var_1 = 1,
             P_c1_intercept = 0.3, P_c1_loading_1 = 0.8,
             P_c2_intercept = -0.2, P_c2_loading_1 = 1.1,
             C_intercept = 0.1, C_cat2 = 0.7, C_cat3 = -0.4, C_sigma = 0.5)
  sim <- simulate_factor_model(ms, base, n = nrow(base), params = truth, seed = 7)
  expect_setequal(sort(unique(sim$cat)), 1:3)
  expect_equal(sim$cat2, as.integer(sim$cat == 2))
  expect_equal(sim$cat3, as.integer(sim$cat == 3))
})

test_that("probit producer feeds its 0/1 value", {
  base <- .endog_base()
  fm <- define_factor_model(n_factors = 1)
  mP <- define_model_component("P", base, "pbin", fm, covariates = "intercept",
                               model_type = "probit", loading_normalization = 1)
  mC <- define_model_component("C", base, "yB", fm, covariates = c("intercept", "pbin"),
                               model_type = "linear", loading_normalization = 0,
                               endogenous = "pbin")
  ms <- define_model_system(components = list(mP, mC), factor = fm)
  truth <- c(factor_var_1 = 1, P_intercept = 0.2, P_loading_1 = 0.9,
             C_intercept = 0.1, C_pbin = 1.3, C_sigma = 0.5)
  sim <- simulate_factor_model(ms, base, n = nrow(base), params = truth, seed = 7)
  expect_setequal(unique(sim$pbin), c(0L, 1L))
  pred <- truth["C_intercept"] + truth["C_pbin"] * sim$pbin
  expect_equal(unname(pred), sim$C_Vobs, tolerance = 1e-10)
})

test_that("components are auto-reordered so producers run before consumers", {
  base <- .endog_base()
  fm <- define_factor_model(n_factors = 1)
  mA <- define_model_component("A", base, "yA", fm, covariates = "intercept",
                               model_type = "linear", loading_normalization = 1)
  mB <- define_model_component("B", base, "yB", fm, covariates = c("intercept", "yA"),
                               model_type = "linear", loading_normalization = NA_real_,
                               endogenous = "yA")
  truth <- c(factor_var_1 = 1, A_intercept = 0.3, A_sigma = 0.6,
             B_intercept = -0.2, B_yA = 0.8, B_loading_1 = 1.1, B_sigma = 0.5)
  # producer-first vs consumer-first must give identical draws for the same seed
  ms_pf <- define_model_system(components = list(mA, mB), factor = fm)
  ms_cf <- define_model_system(components = list(mB, mA), factor = fm)
  s_pf <- simulate_factor_model(ms_pf, base, n = nrow(base), params = truth, seed = 99)
  s_cf <- simulate_factor_model(ms_cf, base, n = nrow(base), params = truth, seed = 99)
  expect_equal(s_pf$yA, s_cf$yA)
  expect_equal(s_pf$yB, s_cf$yB)
})

test_that("cyclic endogenous dependencies are rejected", {
  base <- .endog_base()
  fm <- define_factor_model(n_factors = 1)
  mA <- define_model_component("A", base, "yA", fm, covariates = c("intercept", "yB"),
                               model_type = "linear", loading_normalization = 1,
                               endogenous = "yB")
  mB <- define_model_component("B", base, "yB", fm, covariates = c("intercept", "yA"),
                               model_type = "linear", loading_normalization = NA_real_,
                               endogenous = "yA")
  ms <- define_model_system(components = list(mA, mB), factor = fm)
  truth <- c(factor_var_1 = 1, A_intercept = 0, A_yB = 0.3, A_sigma = 0.5,
             B_intercept = 0, B_yA = 0.3, B_loading_1 = 1, B_sigma = 0.5)
  expect_error(simulate_factor_model(ms, base, n = nrow(base), params = truth, seed = 1),
               "cyclic")
})

test_that("invalid endogenous declarations are caught", {
  base <- .endog_base()
  fm <- define_factor_model(n_factors = 1)
  # define-time: declared name drives none of the component's covariates
  expect_error(
    define_model_component("Z", base, "yA", fm, covariates = "intercept",
                           model_type = "linear", loading_normalization = 1,
                           endogenous = "yB"),
    "drives none")
  # sim-time: declared name matches no component's outcome
  mX <- define_model_component("X", base, "yA", fm, covariates = c("intercept", "yB"),
                               model_type = "linear", loading_normalization = 1,
                               endogenous = "yB")
  msX <- define_model_system(components = list(mX), factor = fm)
  expect_error(
    simulate_factor_model(msX, base, n = nrow(base),
                          params = c(factor_var_1 = 1, X_intercept = 0,
                                     X_yB = 0.3, X_sigma = 0.5), seed = 1),
    "no component produces")
})

test_that("seed reproducibility holds with an endogenous link", {
  base <- .endog_base()
  fm <- define_factor_model(n_factors = 1)
  mA <- define_model_component("A", base, "yA", fm, covariates = "intercept",
                               model_type = "linear", loading_normalization = 1)
  mB <- define_model_component("B", base, "yB", fm, covariates = c("intercept", "yA"),
                               model_type = "linear", loading_normalization = NA_real_,
                               endogenous = "yA")
  ms <- define_model_system(components = list(mA, mB), factor = fm)
  truth <- c(factor_var_1 = 1, A_intercept = 0.3, A_sigma = 0.6,
             B_intercept = -0.2, B_yA = 0.8, B_loading_1 = 1.1, B_sigma = 0.5)
  s1 <- simulate_factor_model(ms, base, n = 200, params = truth, seed = 5)
  s2 <- simulate_factor_model(ms, base, n = 200, params = truth, seed = 5)
  expect_equal(s1, s2)
})
