# Dynamic single-factor test: one latent construct measured at two time
# points, coupled by a linear structural equation.
#
# WORKFLOW (intended standard workflow for SE_linear two-stage estimation
# with a single latent construct observed longitudinally):
#
#   Stage 1. Fit a 2-factor independent measurement model on wide data,
#            with equality constraints tying factor loadings and residual
#            sigmas across periods (measurement invariance), but leaving
#            measurement intercepts period-specific. Under the DGP with
#            E[f_1] = 0 by convention, the wave-1 intercepts identify the
#            true DGP intercepts tau_m.
#
#   Stage 2. Build a dummy 2-factor `previous_stage` object that carries
#            the WAVE-1 intercepts into both factor slots (discarding the
#            wave-2 intercepts, which absorb the factor-mean drift). Fit
#            SE_linear on top, recovering the structural intercept alpha,
#            slope beta, residual variance sigma_eps^2, and f_1 variance.
#
# The "wave-1 intercepts for both slots" move is the critical bit: pooling
# measurement intercepts naively across periods biases tau_m by
# lambda_m * E[f_2] / 2, which propagates into an under-estimate of alpha
# in Stage 2. Using the wave-1 intercepts for both slots anchors the
# measurement level under the factor-identification convention
# E[f_1] = 0, leaving the mean shift between observed wave-1 and wave-2 Y
# to identify alpha in Stage 2.

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
      var_f1    = true_var_f1,
      alpha     = true_alpha,
      beta      = true_beta,
      sigma_e2  = true_sigma_e^2,
      item_int  = item_int,
      item_load = item_load,
      item_sigma = item_sigma
    )
  )
}


# ---- Stage 1: 2-factor independent + tied loadings and sigmas --------------
# Period-specific intercepts. Loadings and sigmas tied via equality
# constraints (measurement invariance for slopes and noise). Returns the
# result AND the constructed model_system (needed to build the dummy
# previous_stage object for Stage 2).

.stage1_2factor_tied_measurement <- function(dat_wide, n_quad = 8) {

  fm <- define_factor_model(n_factors = 2, n_types = 1,
                            factor_structure = "independent")

  mk <- function(name, out, norm) {
    define_model_component(
      name = name, data = dat_wide, outcome = out, factor = fm,
      covariates = "intercept", model_type = "linear",
      loading_normalization = norm,
      use_types = FALSE, evaluation_indicator = "eval"
    )
  }

  comps <- list(
    mk("Y_t1_m1", "Y_t1_m1", c(1, 0)),          # loading fixed to 1 on f_1
    mk("Y_t1_m2", "Y_t1_m2", c(NA_real_, 0)),
    mk("Y_t1_m3", "Y_t1_m3", c(NA_real_, 0)),
    mk("Y_t2_m1", "Y_t2_m1", c(0, 1)),          # loading fixed to 1 on f_2
    mk("Y_t2_m2", "Y_t2_m2", c(0, NA_real_)),
    mk("Y_t2_m3", "Y_t2_m3", c(0, NA_real_))
  )

  # Tie free loadings (m2, m3) and ALL sigmas across the two factors.
  # Intercepts stay free and period-specific. m1 loadings are both fixed to
  # 1 (no free parameter), so no tie is needed there.
  eq_constraints <- list(
    c("Y_t1_m2_loading_1", "Y_t2_m2_loading_2"),
    c("Y_t1_m3_loading_1", "Y_t2_m3_loading_2"),
    c("Y_t1_m1_sigma",     "Y_t2_m1_sigma"),
    c("Y_t1_m2_sigma",     "Y_t2_m2_sigma"),
    c("Y_t1_m3_sigma",     "Y_t2_m3_sigma")
  )

  ms <- define_model_system(
    components           = comps,
    factor               = fm,
    equality_constraints = eq_constraints
  )

  ctrl <- define_estimation_control(n_quad_points = n_quad, num_cores = 1)
  res  <- estimate_model_rcpp(
    model_system = ms, data = dat_wide, control = ctrl,
    optimizer = "nlminb", parallel = FALSE, verbose = FALSE
  )
  list(result = res, model_system = ms)
}


# ---- Build dummy 2-factor previous_stage from Stage 1 ----------------------
# Stage 1 gives period-specific intercepts for each item. We carry the
# WAVE-1 intercepts into BOTH factor slots of the dummy (discarding the
# wave-2 intercepts). This anchors the measurement level under
# E[f_1] = 0, leaving the mean shift between observed Y_t1 and Y_t2 for
# Stage 2's se_intercept to pick up.

.build_dummy_2factor_from_stage1 <- function(stage1, dat_wide) {

  s1_est <- stage1$result$estimates

  fm <- define_factor_model(n_factors = 2, n_types = 1,
                            factor_structure = "independent")
  mk <- function(name, out, norm) {
    define_model_component(
      name = name, data = dat_wide, outcome = out, factor = fm,
      covariates = "intercept", model_type = "linear",
      loading_normalization = norm,
      use_types = FALSE, evaluation_indicator = "eval"
    )
  }
  comps <- list(
    mk("Y_t1_m1", "Y_t1_m1", c(1, 0)),
    mk("Y_t1_m2", "Y_t1_m2", c(NA_real_, 0)),
    mk("Y_t1_m3", "Y_t1_m3", c(NA_real_, 0)),
    mk("Y_t2_m1", "Y_t2_m1", c(0, 1)),
    mk("Y_t2_m2", "Y_t2_m2", c(0, NA_real_)),
    mk("Y_t2_m3", "Y_t2_m3", c(0, NA_real_))
  )
  ms_dummy <- define_model_system(components = comps, factor = fm)

  # Pull the wave-1 intercepts and the tied (shared) loadings / sigmas.
  get <- function(nm) unname(s1_est[nm])
  tau_wave1 <- c(
    m1 = get("Y_t1_m1_intercept"),
    m2 = get("Y_t1_m2_intercept"),
    m3 = get("Y_t1_m3_intercept")
  )
  lambda <- c(
    m2 = get("Y_t1_m2_loading_1"),
    m3 = get("Y_t1_m3_loading_1")
  )
  sigma <- c(
    m1 = get("Y_t1_m1_sigma"),
    m2 = get("Y_t1_m2_sigma"),
    m3 = get("Y_t1_m3_sigma")
  )

  # Build combined parameter vector in the canonical order expected by
  # build_parameter_metadata(): factor_var_1, factor_var_2, then per
  # component (intercept, [loading], sigma).
  vals  <- c(get("factor_var_1"), get("factor_var_2"))
  names_ <- c("factor_var_1", "factor_var_2")

  # Factor-1 components (Y_t1_m1 has no free loading)
  for (i in 1:3) {
    comp <- paste0("Y_t1_m", i)
    vals   <- c(vals,   tau_wave1[[paste0("m", i)]])
    names_ <- c(names_, paste0(comp, "_intercept"))
    if (i > 1) {
      vals   <- c(vals,   lambda[[paste0("m", i)]])
      names_ <- c(names_, paste0(comp, "_loading_1"))
    }
    vals   <- c(vals,   sigma[[paste0("m", i)]])
    names_ <- c(names_, paste0(comp, "_sigma"))
  }
  # Factor-2 components — USE THE WAVE-1 INTERCEPTS HERE TOO
  for (i in 1:3) {
    comp <- paste0("Y_t2_m", i)
    vals   <- c(vals,   tau_wave1[[paste0("m", i)]])   # wave-1, not wave-2
    names_ <- c(names_, paste0(comp, "_intercept"))
    if (i > 1) {
      vals   <- c(vals,   lambda[[paste0("m", i)]])
      names_ <- c(names_, paste0(comp, "_loading_2"))
    }
    vals   <- c(vals,   sigma[[paste0("m", i)]])
    names_ <- c(names_, paste0(comp, "_sigma"))
  }
  names(vals) <- names_

  se_dummy <- setNames(rep(0, length(vals)), names_)

  list(
    model_system = ms_dummy,
    estimates    = vals,
    std_errors   = se_dummy,
    convergence  = 0L,
    loglik       = 0.0
  )
}


# =============================================================================
# TEST 1 — Stage 1 tied measurement recovers DGP intercepts on wave 1
# =============================================================================
test_that("dynamic single-factor Stage 1 (tied loadings/sigmas, period-specific tau) recovers measurement params", {
  skip_on_cran()

  sim <- .simulate_dynamic_single_factor_dgp(n = 1500, seed = 41)
  s1  <- .stage1_2factor_tied_measurement(sim$wide, n_quad = 8)
  expect_equal(s1$result$convergence, 0,
               info = "Stage 1 2-factor tied-measurement must converge strictly")

  est <- s1$result$estimates
  truth <- sim$true

  # Wave-1 intercepts should recover the DGP intercepts (E[f_1] = 0 by
  # convention, so E[Y_t1,m] = tau_m^DGP identifies tau).
  for (i in 1:3) {
    nm <- paste0("Y_t1_m", i, "_intercept")
    expect_equal(unname(est[nm]), truth$item_int[i], tolerance = 0.10,
                 info = sprintf("%s: true=%.3f est=%.3f",
                                nm, truth$item_int[i], est[nm]))
  }

  # Wave-2 intercepts are BIASED by lambda_m * E[f_2] in expectation; we
  # don't assert anything about them — they will be discarded in Stage 2.

  # Tied loadings should match DGP (they were tied across periods and
  # estimated on 2n observations' worth of information).
  expect_equal(unname(est["Y_t1_m2_loading_1"]), truth$item_load[2], tolerance = 0.05)
  expect_equal(unname(est["Y_t1_m3_loading_1"]), truth$item_load[3], tolerance = 0.05)

  # Tied sigmas should match DGP.
  for (i in 1:3) {
    nm <- paste0("Y_t1_m", i, "_sigma")
    expect_equal(unname(est[nm]), truth$item_sigma[i], tolerance = 0.05,
                 info = sprintf("%s: true=%.3f est=%.3f",
                                nm, truth$item_sigma[i], est[nm]))
  }

  if (VERBOSE) {
    cat("\n=== Stage 1 recovery (tied loadings/sigmas, period-specific tau) ===\n")
    cat("DGP tau_m vs wave-1 tau_m (should match):\n")
    for (i in 1:3) {
      cat(sprintf("  m%d: true=%.3f  wave1_est=%.4f  wave2_est=%.4f\n",
                  i, truth$item_int[i],
                  est[paste0("Y_t1_m", i, "_intercept")],
                  est[paste0("Y_t2_m", i, "_intercept")]))
    }
  }
})


# =============================================================================
# TEST 2 — Full workflow: Stage 1 tied -> dummy -> Stage 2 SE_linear
# =============================================================================
test_that("dynamic single-factor: Stage 2 SE_linear recovers alpha, beta, sigma_eps, var_f1", {
  skip_on_cran()

  sim   <- .simulate_dynamic_single_factor_dgp(n = 1500, seed = 41)
  truth <- sim$true

  # Stage 1
  s1 <- .stage1_2factor_tied_measurement(sim$wide, n_quad = 8)
  expect_equal(s1$result$convergence, 0)

  # Dummy previous_stage with wave-1 intercepts in both factor slots
  dummy <- .build_dummy_2factor_from_stage1(s1, sim$wide)

  # Stage 2: SE_linear on top
  fm_s2 <- define_factor_model(n_factors = 2, n_types = 1,
                                factor_structure = "SE_linear")
  ms_s2 <- define_model_system(components = list(), factor = fm_s2,
                                previous_stage = dummy)

  init_s2 <- initialize_parameters(ms_s2, sim$wide, verbose = FALSE)
  init_s2$init_params["factor_var_1"]    <- unname(dummy$estimates["factor_var_1"])
  init_s2$init_params["se_intercept"]    <- 0.0
  init_s2$init_params["se_linear_1"]     <- 0.5
  init_s2$init_params["se_residual_var"] <- 0.5

  ctrl <- define_estimation_control(n_quad_points = 8, num_cores = 1)
  r2 <- estimate_model_rcpp(
    model_system = ms_s2, data = sim$wide,
    init_params = init_s2$init_params, control = ctrl,
    optimizer = "nlminb", parallel = FALSE, verbose = FALSE
  )
  expect_equal(r2$convergence, 0,
               info = "Stage 2 SE_linear must converge strictly")

  est <- r2$estimates
  se  <- r2$std_errors

  # The load-bearing assertion: alpha is recovered. This is the whole
  # reason for using the wave-1 intercepts in the dummy — pooling
  # naively would bias alpha towards zero by roughly half of E[f_2].
  expect_equal(unname(est["se_intercept"]), truth$alpha, tolerance = 0.12,
               info = sprintf("se_intercept (alpha): true=%.3f est=%.3f se=%.3f",
                              truth$alpha, est["se_intercept"], se["se_intercept"]))

  # Other structural params
  expect_equal(unname(est["se_linear_1"]),     truth$beta,     tolerance = 0.10,
               info = sprintf("se_linear_1 (beta): true=%.3f est=%.3f",
                              truth$beta, est["se_linear_1"]))
  expect_equal(unname(est["se_residual_var"]), truth$sigma_e2, tolerance = 0.15,
               info = sprintf("se_residual_var: true=%.3f est=%.3f",
                              truth$sigma_e2, est["se_residual_var"]))
  expect_equal(unname(est["factor_var_1"]),    truth$var_f1,   tolerance = 0.15,
               info = sprintf("factor_var_1: true=%.3f est=%.3f",
                              truth$var_f1, est["factor_var_1"]))

  # All new Stage 2 free params must have positive SEs
  new_free <- c("factor_var_1", "se_intercept", "se_linear_1", "se_residual_var")
  expect_true(all(se[new_free] > 0),
              info = "All Stage 2 free params must have positive std_errors")

  if (VERBOSE) {
    cat("\n=== Stage 2 SE_linear recovery ===\n")
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
