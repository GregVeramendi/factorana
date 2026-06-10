# Tests for robust / cluster-robust sandwich covariance (vcov_factorana).
#
# Strategy: with all factor loadings fixed at 0 the latent factor drops out and
# the model reduces to a plain GLM, so factorana's per-observation scores and
# inverse-information bread coincide with glm()'s. Binary logit is the canonical
# link, so observed information equals expected information and the sandwich
# matches sandwich::vcovHC / vcovCL essentially to machine precision. We also
# check the type='hessian' passthrough, the one-obs-per-cluster reduction to the
# robust estimator, and a genuine latent-factor model for sanity.

make_logit_glm_data <- function(n = 1500, seed = 2024, n_clusters = 40) {
  set.seed(seed)
  x1 <- rnorm(n); x2 <- rnorm(n)
  eta <- -0.3 + 0.8 * x1 - 0.5 * x2
  y01 <- rbinom(n, 1, 1 / (1 + exp(-eta)))
  data.frame(intercept = 1, x1 = x1, x2 = x2,
             y = y01 + 1L, y01 = y01,
             cl = sample(seq_len(n_clusters), n, replace = TRUE))
}

fit_glm_equiv <- function(dat) {
  fm <- define_factor_model(n_factors = 1)
  fm <- fix_factor_param(fm, c(factor_var_1 = 1))
  mc <- define_model_component("y", dat, "y", fm,
                               covariates = c("intercept", "x1", "x2"),
                               model_type = "logit", num_choices = 2,
                               loading_normalization = 0)
  ms <- define_model_system(components = list(mc), factor = fm)
  estimate_model_rcpp(ms, dat,
                      control = define_estimation_control(n_quad_points = 8, num_cores = 1),
                      parallel = FALSE, optimizer = "nlminb", verbose = FALSE)
}

test_that("robust and cluster SEs match sandwich on a GLM-equivalent logit", {
  skip_on_cran()
  skip_if_not_installed("sandwich")
  dat <- make_logit_glm_data()
  res <- fit_glm_equiv(dat)
  expect_equal(res$convergence, 0)

  bn <- c("y_intercept", "y_x1", "y_x2")
  g <- glm(y01 ~ x1 + x2, family = binomial("logit"), data = dat)

  # Point estimates match glm.
  expect_equal(unname(res$estimates[bn]), unname(coef(g)), tolerance = 1e-4)

  # Model-based (inverse information) SEs match glm vcov.
  se_model_fa  <- res$std_errors[bn]
  se_model_glm <- sqrt(diag(vcov(g)))
  expect_equal(unname(se_model_fa), unname(se_model_glm), tolerance = 1e-5)

  # Robust matches vcovHC(HC0).
  se_rob_fa  <- sqrt(diag(vcov_factorana(res, dat, type = "robust", finite_sample = FALSE)))[bn]
  se_rob_glm <- sqrt(diag(sandwich::vcovHC(g, type = "HC0")))
  expect_equal(unname(se_rob_fa), unname(se_rob_glm), tolerance = 1e-5)

  # Cluster matches vcovCL(HC0, cadjust = FALSE).
  se_clu_fa  <- sqrt(diag(vcov_factorana(res, dat, type = "cluster", cluster = "cl",
                                         finite_sample = FALSE)))[bn]
  se_clu_glm <- sqrt(diag(sandwich::vcovCL(g, cluster = dat$cl, type = "HC0", cadjust = FALSE)))
  expect_equal(unname(se_clu_fa), unname(se_clu_glm), tolerance = 1e-5)
})

test_that("type='hessian' reproduces the fit's standard errors", {
  skip_on_cran()
  dat <- make_logit_glm_data(n = 600)
  res <- fit_glm_equiv(dat)
  expect_equal(res$convergence, 0)
  V <- vcov_factorana(res, dat, type = "hessian")
  free <- res$free_idx
  se_V <- sqrt(diag(V))[free]
  expect_equal(unname(se_V), unname(res$std_errors[free]), tolerance = 1e-8)
})

test_that("cluster with one observation per cluster equals the robust estimator", {
  skip_on_cran()
  dat <- make_logit_glm_data(n = 500)
  dat$rowid <- seq_len(nrow(dat))  # each row its own cluster
  res <- fit_glm_equiv(dat)
  expect_equal(res$convergence, 0)
  V_rob <- vcov_factorana(res, dat, type = "robust", finite_sample = FALSE)
  V_clu <- suppressWarnings(
    vcov_factorana(res, dat, type = "cluster", cluster = "rowid", finite_sample = FALSE))
  expect_equal(V_clu, V_rob, tolerance = 1e-8)
})

test_that("vcov_factorana works on a genuine latent-factor model", {
  skip_on_cran()
  skip_if_not_installed("sandwich")
  set.seed(7)
  n <- 1200
  f <- rnorm(n)
  dat <- data.frame(
    intercept = 1, x = rnorm(n),
    y1 = 1.0 * f + rnorm(n, 0, 0.6),
    y2 = 0.8 * f + rnorm(n, 0, 0.6),
    y3 = 0.9 * f + 0.5 * rnorm(n) + rnorm(n, 0, 0.6),
    cl = sample(1:50, n, replace = TRUE)
  )
  fm <- define_factor_model(n_factors = 1)
  mc1 <- define_model_component("m1", dat, "y1", fm, covariates = "intercept",
                                model_type = "linear", loading_normalization = 1)
  mc2 <- define_model_component("m2", dat, "y2", fm, covariates = "intercept",
                                model_type = "linear", loading_normalization = NA_real_)
  mc3 <- define_model_component("m3", dat, "y3", fm, covariates = c("intercept", "x"),
                                model_type = "linear", loading_normalization = NA_real_)
  ms <- define_model_system(components = list(mc1, mc2, mc3), factor = fm)
  res <- estimate_model_rcpp(ms, dat,
                             control = define_estimation_control(n_quad_points = 12, num_cores = 1),
                             parallel = FALSE, optimizer = "nlminb", verbose = FALSE)
  expect_equal(res$convergence, 0)

  V_rob <- vcov_factorana(res, dat, type = "robust")
  V_clu <- vcov_factorana(res, dat, type = "cluster", cluster = "cl")
  free <- res$free_idx

  # Covariances are finite, symmetric, with positive diagonal on free params.
  expect_true(all(is.finite(V_rob)) && all(is.finite(V_clu)))
  expect_equal(V_rob, t(V_rob), tolerance = 1e-8)
  expect_true(all(diag(V_rob)[free] > 0))
  expect_true(all(diag(V_clu)[free] > 0))

  # Robust SEs are in a sane neighbourhood of the model SEs (same order of magnitude).
  se_rob <- sqrt(diag(V_rob))[free]
  se_mod <- res$std_errors[free]
  expect_true(all(se_rob / se_mod > 0.25 & se_rob / se_mod < 4))
})

test_that("robust_se returns a named vector (name-based indexing works)", {
  skip_on_cran()
  dat <- make_logit_glm_data(n = 600)
  res <- fit_glm_equiv(dat)
  expect_equal(res$convergence, 0)

  se <- robust_se(res, dat, type = "robust")
  # the naming contract: one named entry per model parameter
  expect_named(se, names(res$estimates))
  # name-based indexing must not silently return NA for free parameters
  free <- res$free_idx
  expect_false(anyNA(se[names(res$estimates)[free]]))

  # cluster path keeps names too
  se_cl <- robust_se(res, dat, type = "cluster", cluster = "cl")
  expect_named(se_cl, names(res$estimates))
})

test_that("vcov_factorana rejects mismatched data (fingerprint guard)", {
  skip_on_cran()
  dat <- make_logit_glm_data(n = 500)
  res <- fit_glm_equiv(dat)
  expect_equal(res$convergence, 0)

  # same data: fine
  expect_silent(invisible(vcov_factorana(res, dat, type = "robust")))
  # same shape, different row order: caught
  expect_error(vcov_factorana(res, dat[sample(nrow(dat)), ], type = "robust"),
               "does not match")
  # same shape, different rows: caught
  expect_error(vcov_factorana(res, dat[sample(nrow(dat), replace = TRUE), ], type = "robust"),
               "does not match")
})
