# Tests for the bootstrap primitives (generate_bootstrap_samples,
# bootstrap_fit_sample, collect_bootstrap, bootstrap_factorana).
#
# The bootstrap resamples clusters (or rows) with replacement and applies the
# multiplicities as integer frequency weights, which is numerically identical
# to estimating on the physically expanded data. Each (stage, sample) fit is an
# idempotent file in a directory, so the workflow is restartable and
# distributable, and multi-stage fits chain by pointing each replicate's later
# stage at its own earlier-stage result.

make_boot_data <- function(n = 250, seed = 3, n_clusters = 30) {
  set.seed(seed)
  f <- rnorm(n)
  data.frame(
    intercept = 1,
    school = sample(seq_len(n_clusters), n, replace = TRUE),
    y1 = 1.0 * f + rnorm(n, 0, 0.6),
    y2 = 0.8 * f + rnorm(n, 0, 0.6),
    y3 = 1.1 * f + rnorm(n, 0, 0.6)
  )
}

build_1f_system <- function(dat, weights = NULL) {
  fm <- define_factor_model(n_factors = 1)
  mc1 <- define_model_component("m1", dat, "y1", fm, covariates = "intercept",
                                model_type = "linear", loading_normalization = 1)
  mc2 <- define_model_component("m2", dat, "y2", fm, covariates = "intercept",
                                model_type = "linear", loading_normalization = NA_real_)
  mc3 <- define_model_component("m3", dat, "y3", fm, covariates = "intercept",
                                model_type = "linear", loading_normalization = NA_real_)
  define_model_system(components = list(mc1, mc2, mc3), factor = fm, weights = weights)
}

test_that("generate_bootstrap_samples produces valid frequency weights", {
  dat <- make_boot_data()
  # Row bootstrap: each replicate's weights sum to nobs (n draws).
  s_row <- generate_bootstrap_samples(dat, R = 20, cluster = NULL, seed = 1)
  expect_equal(dim(s_row$weights), c(nrow(dat), 20))
  expect_true(all(colSums(s_row$weights) == nrow(dat)))
  expect_true(all(s_row$weights >= 0))

  # Cluster bootstrap: each row's weight equals its cluster's draw count, so
  # within a cluster all rows share the same weight.
  s_cl <- generate_bootstrap_samples(dat, R = 20, cluster = "school", seed = 1)
  expect_equal(s_cl$n_clusters, length(unique(dat$school)))
  for (b in 1:3) {
    w <- s_cl$weights[, b]
    by_cluster <- tapply(w, dat$school, function(z) length(unique(z)))
    expect_true(all(by_cluster == 1))  # constant within cluster
  }
})

test_that("frequency-weighted estimation equals estimation on expanded data", {
  skip_on_cran()
  dat <- make_boot_data(n = 300, seed = 5)
  set.seed(9)
  freq <- as.integer(tabulate(sample(nrow(dat), nrow(dat), replace = TRUE),
                              nbins = nrow(dat)))
  ctrl <- define_estimation_control(n_quad_points = 10, num_cores = 1)

  keep <- freq > 0
  dw <- dat[keep, ]; dw$w <- freq[keep]
  fit_w <- estimate_model_rcpp(build_1f_system(dw, weights = "w"), dw, control = ctrl,
                               optimizer = "nlminb", parallel = FALSE, verbose = FALSE)

  de <- dat[rep(seq_len(nrow(dat)), freq), ]
  fit_e <- estimate_model_rcpp(build_1f_system(de), de, control = ctrl,
                               optimizer = "nlminb", parallel = FALSE, verbose = FALSE)

  expect_equal(fit_w$convergence, 0)
  expect_equal(fit_e$convergence, 0)
  expect_equal(unname(fit_w$estimates[names(fit_e$estimates)]),
               unname(fit_e$estimates), tolerance = 1e-5)
})

test_that("cluster bootstrap runs, is restartable, and agrees with the cluster sandwich", {
  skip_on_cran()
  dat <- make_boot_data(n = 250, seed = 3)
  ms <- build_1f_system(dat)
  ctrl <- define_estimation_control(n_quad_points = 8, num_cores = 1)
  fit <- estimate_model_rcpp(ms, dat, control = ctrl, optimizer = "nlminb",
                             parallel = FALSE, verbose = FALSE)
  expect_equal(fit$convergence, 0)

  bdir <- file.path(tempdir(), paste0("boot_", as.integer(Sys.time()) %% 100000))
  on.exit(unlink(bdir, recursive = TRUE), add = TRUE)

  boot <- bootstrap_factorana(ms, dat, R = 40, cluster = "school", stage_dir = bdir,
                              control = ctrl, seed = 123, optimizer = "nlminb",
                              parallel = FALSE, verbose = FALSE)
  expect_s3_class(boot, "factorana_bootstrap")
  expect_equal(boot$n_total, 40)
  expect_gte(boot$n_converged, 38)  # essentially all converge
  expect_true(all(is.finite(boot$boot_se)))

  free <- fit$free_idx
  expect_true(all(boot$boot_se[names(fit$estimates)[free]] > 0))
  # Percentile interval should bracket the point estimate for most parameters.
  ci <- boot$ci
  rownames(ci) <- ci$parameter
  brackets <- ci[names(fit$estimates)[free], "lower"] <= fit$estimates[free] &
              fit$estimates[free] <= ci[names(fit$estimates)[free], "upper"]
  expect_gte(mean(brackets), 0.8)

  # Bootstrap SE should be the same order of magnitude as the cluster sandwich.
  se_clu <- robust_se(fit, dat, type = "cluster", cluster = "school",
                      finite_sample = FALSE)
  ratio <- boot$boot_se[names(fit$estimates)[free]] / se_clu[free]
  expect_true(all(ratio > 0.4 & ratio < 2.5))

  # Restartability: a second call re-uses the saved samples (instant, identical).
  t0 <- proc.time()[3]
  boot2 <- bootstrap_factorana(ms, dat, R = 40, cluster = "school", stage_dir = bdir,
                               control = ctrl, seed = 123, optimizer = "nlminb",
                               parallel = FALSE, verbose = FALSE)
  expect_lt(proc.time()[3] - t0, 5)  # all skipped, no re-estimation
  expect_equal(boot2$boot_se, boot$boot_se)
})

test_that("two-stage bootstrap chains each replicate on its own stage-1 fit", {
  skip_on_cran()
  # Stage 1: dynamic measurement (2 periods, 3 items). Stage 2: SE_linear.
  set.seed(41)
  n <- 600
  f1 <- rnorm(n)
  f2 <- 0.4 + 0.6 * f1 + rnorm(n, 0, sqrt(0.5))
  gy <- function(f, i) c(1.5, 1.0, 0.8)[i] + c(1.0, 0.9, 1.1)[i] * f +
    rnorm(length(f), 0, c(0.7, 0.75, 0.65)[i])
  dat <- data.frame(
    intercept = 1, eval = 1L, person = seq_len(n),
    Y_t1_m1 = gy(f1, 1), Y_t1_m2 = gy(f1, 2), Y_t1_m3 = gy(f1, 3),
    Y_t2_m1 = gy(f2, 1), Y_t2_m2 = gy(f2, 2), Y_t2_m3 = gy(f2, 3)
  )
  dyn <- define_dynamic_measurement(
    data = dat, items = c("m1", "m2", "m3"),
    period_prefixes = c("Y_t1_", "Y_t2_"),
    model_type = "linear", evaluation_indicator = "eval")
  ctrl <- define_estimation_control(n_quad_points = 6, num_cores = 1)

  R <- 6
  samples <- generate_bootstrap_samples(dat, R = R, cluster = "person", seed = 7)
  s1dir <- file.path(tempdir(), paste0("b1_", as.integer(Sys.time()) %% 100000))
  s2dir <- file.path(tempdir(), paste0("b2_", as.integer(Sys.time()) %% 100000))
  on.exit(unlink(c(s1dir, s2dir), recursive = TRUE), add = TRUE)

  for (b in seq_len(R)) {
    # Stage 1 for sample b.
    f1file <- bootstrap_fit_sample(dyn$model_system, dat, samples, sample_id = b,
                                   stage_dir = s1dir, control = ctrl,
                                   optimizer = "nlminb", parallel = FALSE, verbose = FALSE)
    res1 <- readRDS(f1file)

    # Build Stage 2 for sample b, chaining on THIS sample's Stage 1 fit.
    prev <- build_dynamic_previous_stage(dyn = dyn, stage1_result = res1,
                                         data = dat, anchor_period = 1L)
    fm2 <- define_factor_model(n_factors = 2, n_types = 1, factor_structure = "SE_linear")
    ms2 <- define_model_system(components = list(), factor = fm2, previous_stage = prev)

    bootstrap_fit_sample(ms2, dat, samples, sample_id = b,
                         stage_dir = s2dir, control = ctrl,
                         optimizer = "nlminb", parallel = FALSE, verbose = FALSE)
  }

  boot2 <- collect_bootstrap(s2dir, require_convergence = FALSE)
  expect_equal(boot2$n_total, R)
  expect_true("se_linear_1" %in% colnames(boot2$estimates))
  expect_true("se_intercept" %in% colnames(boot2$estimates))
  # Stage-2 structural estimates are finite across replicates.
  expect_true(all(is.finite(boot2$estimates[, "se_linear_1"])))
})

test_that("bootstrap_factorana_multistage drives a two-stage bootstrap in one call", {
  skip_on_cran()
  set.seed(41)
  n <- 600
  f1 <- rnorm(n)
  f2 <- 0.4 + 0.6 * f1 + rnorm(n, 0, sqrt(0.5))
  gy <- function(f, i) c(1.5, 1.0, 0.8)[i] + c(1.0, 0.9, 1.1)[i] * f +
    rnorm(length(f), 0, c(0.7, 0.75, 0.65)[i])
  dat <- data.frame(
    intercept = 1, eval = 1L, person = seq_len(n),
    Y_t1_m1 = gy(f1, 1), Y_t1_m2 = gy(f1, 2), Y_t1_m3 = gy(f1, 3),
    Y_t2_m1 = gy(f2, 1), Y_t2_m2 = gy(f2, 2), Y_t2_m3 = gy(f2, 3)
  )
  dyn <- define_dynamic_measurement(
    data = dat, items = c("m1", "m2", "m3"),
    period_prefixes = c("Y_t1_", "Y_t2_"),
    model_type = "linear", evaluation_indicator = "eval")
  ctrl <- define_estimation_control(n_quad_points = 6, num_cores = 1)

  bdir <- file.path(tempdir(), paste0("bms_", as.integer(Sys.time()) %% 100000))
  on.exit(unlink(bdir, recursive = TRUE), add = TRUE)

  stage_builders <- list(
    # Stage 1: measurement system (prev_fits is empty here).
    function(prev_fits, data) dyn$model_system,
    # Stage 2: SE_linear, chaining on this replicate's Stage 1 fit.
    function(prev_fits, data) {
      prev <- build_dynamic_previous_stage(dyn = dyn, stage1_result = prev_fits[[1]],
                                           data = data, anchor_period = 1L)
      fm2 <- define_factor_model(n_factors = 2, n_types = 1, factor_structure = "SE_linear")
      define_model_system(components = list(), factor = fm2, previous_stage = prev)
    }
  )

  boot <- bootstrap_factorana_multistage(
    stage_builders, dat, R = 6, cluster = "person", dir = bdir,
    control = ctrl, seed = 7, optimizer = "nlminb", parallel = FALSE, verbose = FALSE)

  expect_s3_class(boot, "factorana_bootstrap_multistage")
  expect_equal(boot$n_stages, 2)
  expect_false(is.null(boot$final))
  expect_true("se_linear_1" %in% boot$final$ci$parameter)

  # Restartable: a second call skips all fits (fast) and returns the same SEs.
  t0 <- proc.time()[3]
  boot_b <- bootstrap_factorana_multistage(
    stage_builders, dat, R = 6, cluster = "person", dir = bdir,
    control = ctrl, seed = 7, optimizer = "nlminb", parallel = FALSE, verbose = FALSE)
  expect_lt(proc.time()[3] - t0, 8)
  expect_equal(boot_b$final$boot_se, boot$final$boot_se)
})

test_that("bootstrap fits are identical in parallel and serial mode", {
  skip_on_cran()
  skip_on_os("windows")  # FORK cluster required for load_all + parallel
  dat <- make_boot_data(n = 300, seed = 3)
  ms <- build_1f_system(dat)
  samples <- generate_bootstrap_samples(dat, R = 3, cluster = "school", seed = 123)
  ds <- file.path(tempdir(), paste0("bps_", as.integer(Sys.time()) %% 100000))
  dp <- file.path(tempdir(), paste0("bpp_", as.integer(Sys.time()) %% 100000))
  on.exit(unlink(c(ds, dp), recursive = TRUE), add = TRUE)
  ctrl1 <- define_estimation_control(n_quad_points = 8, num_cores = 1)
  ctrl2 <- define_estimation_control(n_quad_points = 8, num_cores = 2)

  f_ser <- readRDS(bootstrap_fit_sample(ms, dat, samples, 1, ds, control = ctrl1,
                                        optimizer = "nlminb", parallel = FALSE, verbose = FALSE))
  f_par <- readRDS(bootstrap_fit_sample(ms, dat, samples, 1, dp, control = ctrl2,
                                        optimizer = "nlminb", parallel = TRUE, verbose = FALSE))
  expect_equal(f_ser$convergence, 0)
  expect_equal(f_par$convergence, 0)
  # Weighted estimation must give the same answer whether the (weighted) data is
  # split across workers or not.
  expect_equal(unname(f_par$estimates[names(f_ser$estimates)]),
               unname(f_ser$estimates), tolerance = 1e-6)
})

test_that("bootstrap_factorana runs end-to-end with parallel estimation", {
  skip_on_cran()
  skip_on_os("windows")
  dat <- make_boot_data(n = 250, seed = 3)
  ms <- build_1f_system(dat)
  ctrl <- define_estimation_control(n_quad_points = 8, num_cores = 2)
  bdir <- file.path(tempdir(), paste0("bpe_", as.integer(Sys.time()) %% 100000))
  on.exit(unlink(bdir, recursive = TRUE), add = TRUE)
  boot <- bootstrap_factorana(ms, dat, R = 8, cluster = "school", stage_dir = bdir,
                              control = ctrl, seed = 1, optimizer = "nlminb",
                              parallel = TRUE, verbose = FALSE)
  expect_equal(boot$n_total, 8)
  expect_gte(boot$n_converged, 7)
  expect_true(all(is.finite(boot$boot_se)))
})
