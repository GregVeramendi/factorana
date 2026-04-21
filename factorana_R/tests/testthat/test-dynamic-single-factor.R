# Dynamic single-factor test: one latent construct measured at two time
# points, coupled by a linear structural equation.
#
# This test uses the define_dynamic_measurement() / build_dynamic_previous_stage()
# wrapper functions, which encapsulate:
#   - Stage 1: 2-factor independent measurement system with
#     equality_constraints tying loadings and sigmas across periods and
#     period-specific intercepts.
#   - Stage 2 bridge: build a dummy previous_stage that carries the
#     anchor-period (wave 1) intercepts into every factor slot.
#
# Then Stage 2 SE_linear recovers alpha, beta, sigma_eps^2, and Var(f_1).

VERBOSE <- Sys.getenv("FACTORANA_TEST_VERBOSE", "FALSE") == "TRUE"


# ---- DGP helper -------------------------------------------------------------

.simulate_dynamic_single_factor_dgp <- function(
    n            = 1500,
    seed         = 41,
    true_var_f1  = 1.0,
    true_alpha   = 0.4,
    true_beta    = 0.6,
    true_sigma_e = sqrt(0.5),
    item_int     = c(1.5, 1.0, 0.8),
    item_load    = c(1.0, 0.9, 1.1),
    item_sigma   = c(0.7, 0.75, 0.65)) {

  set.seed(seed)

  f1  <- rnorm(n, 0, sqrt(true_var_f1))
  eps <- rnorm(n, 0, true_sigma_e)
  f2  <- true_alpha + true_beta * f1 + eps

  gen_Y <- function(f, i) {
    item_int[i] + item_load[i] * f + rnorm(length(f), 0, item_sigma[i])
  }

  dat_wide <- data.frame(
    id        = seq_len(n),
    intercept = 1,
    eval      = 1L,
    Y_t1_m1 = gen_Y(f1, 1), Y_t1_m2 = gen_Y(f1, 2), Y_t1_m3 = gen_Y(f1, 3),
    Y_t2_m1 = gen_Y(f2, 1), Y_t2_m2 = gen_Y(f2, 2), Y_t2_m3 = gen_Y(f2, 3)
  )

  list(
    wide = dat_wide,
    true = list(
      var_f1     = true_var_f1,
      alpha      = true_alpha,
      beta       = true_beta,
      sigma_e2   = true_sigma_e^2,
      item_int   = item_int,
      item_load  = item_load,
      item_sigma = item_sigma
    )
  )
}


# =============================================================================
# TEST 1: Stage 1 tied measurement recovers DGP intercepts on wave 1
# =============================================================================
test_that("define_dynamic_measurement: Stage 1 recovers wave-1 tau_m, shared loadings/sigmas", {
  skip_on_cran()

  sim   <- .simulate_dynamic_single_factor_dgp(n = 1500, seed = 41)
  truth <- sim$true

  dyn <- define_dynamic_measurement(
    data                 = sim$wide,
    items                = c("m1", "m2", "m3"),
    period_prefixes      = c("Y_t1_", "Y_t2_"),
    model_type           = "linear",
    evaluation_indicator = "eval"
  )
  expect_s3_class(dyn, "dynamic_measurement")

  # Equality constraints: 2 loadings (items m2, m3) + 3 sigmas = 5 groups
  expect_equal(length(dyn$equality_constraints), 5L)

  ctrl <- define_estimation_control(n_quad_points = 8, num_cores = 1)
  s1 <- estimate_model_rcpp(
    dyn$model_system, sim$wide, control = ctrl,
    optimizer = "nlminb", parallel = FALSE, verbose = FALSE
  )
  expect_equal(s1$convergence, 0,
               info = "Stage 1 must converge strictly")

  est <- s1$estimates

  # Wave-1 intercepts recover DGP (E[f_1] = 0 by convention)
  for (i in 1:3) {
    nm <- paste0("Y_t1_m", i, "_intercept")
    expect_equal(unname(est[nm]), truth$item_int[i], tolerance = 0.10,
                 info = sprintf("%s: true=%.3f est=%.3f",
                                nm, truth$item_int[i], est[nm]))
  }

  # Tied loadings and sigmas recover DGP
  expect_equal(unname(est["Y_t1_m2_loading_1"]), truth$item_load[2], tolerance = 0.05)
  expect_equal(unname(est["Y_t1_m3_loading_1"]), truth$item_load[3], tolerance = 0.05)
  for (i in 1:3) {
    nm <- paste0("Y_t1_m", i, "_sigma")
    expect_equal(unname(est[nm]), truth$item_sigma[i], tolerance = 0.05)
  }

  # Equality constraints hold exactly in the estimates
  expect_equal(unname(est["Y_t1_m2_loading_1"]),
               unname(est["Y_t2_m2_loading_2"]), tolerance = 1e-10)
  expect_equal(unname(est["Y_t1_m1_sigma"]),
               unname(est["Y_t2_m1_sigma"]),     tolerance = 1e-10)
})


# =============================================================================
# TEST 2: End-to-end: wrapper plus Stage 2 SE_linear recovers structural params
# =============================================================================
test_that("define_dynamic_measurement + Stage 2 SE_linear recovers alpha, beta, sigma_eps, var_f1", {
  skip_on_cran()

  sim   <- .simulate_dynamic_single_factor_dgp(n = 1500, seed = 41)
  truth <- sim$true
  ctrl  <- define_estimation_control(n_quad_points = 8, num_cores = 1)

  # Stage 1 via the wrapper
  dyn <- define_dynamic_measurement(
    data                 = sim$wide,
    items                = c("m1", "m2", "m3"),
    period_prefixes      = c("Y_t1_", "Y_t2_"),
    model_type           = "linear",
    evaluation_indicator = "eval"
  )
  s1 <- estimate_model_rcpp(
    dyn$model_system, sim$wide, control = ctrl,
    optimizer = "nlminb", parallel = FALSE, verbose = FALSE
  )
  expect_equal(s1$convergence, 0)

  # Build the Stage 2 previous_stage bridge via the wrapper
  dummy <- build_dynamic_previous_stage(dyn, s1, sim$wide, anchor_period = 1L)

  # Stage 2: SE_linear
  fm_s2 <- define_factor_model(n_factors = 2, n_types = 1,
                                factor_structure = "SE_linear")
  ms_s2 <- define_model_system(components = list(), factor = fm_s2,
                                previous_stage = dummy)

  init_s2 <- initialize_parameters(ms_s2, sim$wide, verbose = FALSE)
  init_s2$init_params["factor_var_1"]    <- unname(dummy$estimates["factor_var_1"])
  init_s2$init_params["se_intercept"]    <- 0.0
  init_s2$init_params["se_linear_1"]     <- 0.5
  init_s2$init_params["se_residual_var"] <- 0.5

  r2 <- estimate_model_rcpp(
    ms_s2, sim$wide, init_params = init_s2$init_params, control = ctrl,
    optimizer = "nlminb", parallel = FALSE, verbose = FALSE
  )
  expect_equal(r2$convergence, 0)

  est <- r2$estimates
  se  <- r2$std_errors

  # The load-bearing assertion: alpha recovers.
  expect_equal(unname(est["se_intercept"]), truth$alpha, tolerance = 0.12,
               info = sprintf("se_intercept (alpha): true=%.3f est=%.3f se=%.3f",
                              truth$alpha, est["se_intercept"], se["se_intercept"]))

  # Other structural params
  expect_equal(unname(est["se_linear_1"]),     truth$beta,     tolerance = 0.10)
  expect_equal(unname(est["se_residual_var"]), truth$sigma_e2, tolerance = 0.15)
  expect_equal(unname(est["factor_var_1"]),    truth$var_f1,   tolerance = 0.15)

  new_free <- c("factor_var_1", "se_intercept", "se_linear_1", "se_residual_var")
  expect_true(all(se[new_free] > 0))

  if (VERBOSE) {
    cat("\n=== Stage 2 SE_linear recovery (via wrapper) ===\n")
    cat(sprintf("%-22s %10s %10s %10s\n", "param", "true", "est", "se"))
    for (p in new_free) {
      tr <- switch(p,
                   factor_var_1    = truth$var_f1,
                   se_intercept    = truth$alpha,
                   se_linear_1     = truth$beta,
                   se_residual_var = truth$sigma_e2)
      cat(sprintf("  %-22s %+10.4f %+10.4f %10.4f\n",
                  p, tr, est[p], se[p]))
    }
  }
})
