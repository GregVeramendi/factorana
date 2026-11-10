# Tests for ranked (exploded) multinomial logit in simulate_factor_model(). The
# simulation mirrors the rank loop of Model::EvalLogit() in src/Model.cpp: one
# draw per rank column, a choice set that shrinks by rank when exclude_chosen is
# TRUE and stays fixed when it is FALSE, and rank-share corrections added to the
# utilities of the non-reference alternatives. Validation is by parameter
# recovery, as in test-simulate-factor-model.R: simulate from known parameters,
# re-estimate with nlminb, and require every estimate within 3.5 SE of truth.
#
# The recovery fits use n_quad = 48, as test-simulate-endogenous.R does. Three
# precise indicators plus three ranks of choice information give a tight factor
# posterior. At 12 nodes the exclude_chosen fit converged with m3_sigma 5.7 SE
# from truth and a log-likelihood 51 units below the 48-node value, and at 24
# nodes every estimate was within 1.7 SE (measured November 2026).

K_EXP <- 5L
NRANK_EXP <- 3L

# Base rows: covariates, three linear indicators, rank columns and a rank-share
# block of (K - 1) * nrank = 12 contiguous columns. Rank 1 is always filled, and
# ranks 2 and 3 are filled for 80% and 60% of rows, so the "base" fill rule has a
# pattern to reproduce. One correction in ten is missing (NA).
.exp_base <- function(nb = 4000, seed = 1) {
  set.seed(seed)
  d <- data.frame(intercept = 1, x = rnorm(nb),
                  y1 = rnorm(nb), y2 = rnorm(nb), y3 = rnorm(nb))
  u <- runif(nb)
  d$r1 <- sample.int(K_EXP, nb, TRUE)
  d$r2 <- ifelse(u < 0.8, sample.int(K_EXP, nb, TRUE), NA_integer_)
  d$r3 <- ifelse(u < 0.6, sample.int(K_EXP, nb, TRUE), NA_integer_)
  for (r in 0:(NRANK_EXP - 1L)) for (c in 0:(K_EXP - 2L)) {
    v <- rnorm(nb, 0, 0.5)
    v[runif(nb) < 0.1] <- NA
    d[[sprintf("rs_r%d_c%d", r, c)]] <- v
  }
  d
}

# Rank-share block column names in data order: rank 0 alternatives 1..K-1, then
# rank 1, and so on, as Model::EvalLogit() reads them by position.
.exp_rs_cols <- unlist(lapply(0:(NRANK_EXP - 1L), function(r)
  sprintf("rs_r%d_c%d", r, 0:(K_EXP - 2L))))

.exp_system <- function(base, exclude_chosen, rankshare = FALSE) {
  fm <- define_factor_model(n_factors = 1)
  comps <- list(
    define_model_component("m1", base, "y1", fm, covariates = "intercept",
                           model_type = "linear", loading_normalization = 1),
    define_model_component("m2", base, "y2", fm, covariates = "intercept",
                           model_type = "linear", loading_normalization = NA_real_),
    define_model_component("m3", base, "y3", fm, covariates = "intercept",
                           model_type = "linear", loading_normalization = NA_real_),
    define_model_component("app", base, c("r1", "r2", "r3"), fm,
                           covariates = c("intercept", "x"), model_type = "logit",
                           num_choices = K_EXP, exclude_chosen = exclude_chosen,
                           rankshare_var = if (rankshare) "rs_r0_c0" else NULL,
                           loading_normalization = NA_real_))
  define_model_system(components = comps, factor = fm)
}

.exp_truth <- c(factor_var_1 = 1,
                m1_intercept = 0.2, m1_sigma = 0.6,
                m2_intercept = -0.1, m2_loading_1 = 0.9, m2_sigma = 0.6,
                m3_intercept = 0.3, m3_loading_1 = 1.1, m3_sigma = 0.5,
                app_c1_intercept = 0.4, app_c1_x = 0.5, app_c1_loading_1 = 0.8,
                app_c2_intercept = -0.2, app_c2_x = -0.4, app_c2_loading_1 = 0.3,
                app_c3_intercept = 0.1, app_c3_x = 0.2, app_c3_loading_1 = -0.5,
                app_c4_intercept = -0.3, app_c4_x = 0.6, app_c4_loading_1 = 1.0)

# Rank-invariant utility of non-reference alternative a (1..K-1) for each
# simulated row, from the truth and the simulated covariate and factor.
.exp_utility <- function(sim, a) {
  unname(.exp_truth[sprintf("app_c%d_intercept", a)] +
           .exp_truth[sprintf("app_c%d_x", a)] * sim$x +
           .exp_truth[sprintf("app_c%d_loading_1", a)] * sim$factor_1)
}

.exp_rank_mat <- function(sim) unname(as.matrix(sim[, c("r1", "r2", "r3")]))
.exp_prob_mat <- function(sim, r) {
  unname(as.matrix(sim[, sprintf("app_r%d_prob_C%d", r, seq_len(K_EXP))]))
}

# The draws at rank r must follow the reported probabilities: for each alternative
# the number drawn is within 3.5 binomial SEs of the summed probability.
.exp_check_draws <- function(sim, r) {
  y <- .exp_rank_mat(sim)[, r]
  p <- .exp_prob_mat(sim, r)
  ok <- !is.na(y)
  for (j in seq_len(K_EXP)) {
    pj <- p[ok, j]
    z <- (sum(y[ok] == j) - sum(pj)) / sqrt(sum(pj * (1 - pj)))
    expect_true(abs(z) < 3.5, info = sprintf("rank %d alternative %d: z = %.2f", r, j, z))
  }
}

# Re-estimate on the simulated data and require strict convergence, finite and
# moderate SEs, and every estimate within 3.5 SE of truth.
.exp_recover <- function(ms, sim, truth = .exp_truth, n_quad = 48) {
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
  invisible(fit)
}

test_that("exclude_chosen = TRUE: no repeats within a row, base fill pattern, recovery", {
  skip_on_cran()
  base <- .exp_base()
  ms <- .exp_system(base, exclude_chosen = TRUE)
  sim <- simulate_factor_model(ms, base, n = 6000, params = .exp_truth, seed = 11)
  R <- .exp_rank_mat(sim)

  # every filled rank holds a valid alternative and no alternative repeats within a row
  expect_true(all(R[!is.na(R)] %in% seq_len(K_EXP)))
  dup <- apply(R, 1, function(v) anyDuplicated(v[!is.na(v)]) > 0)
  expect_false(any(dup))

  # rank 1 is always filled, and ranks 2 and 3 are filled exactly where the base row was
  expect_false(anyNA(R[, 1]))
  expect_equal(is.na(R[, 2]), is.na(base$r2[sim$xid]))
  expect_equal(is.na(R[, 3]), is.na(base$r3[sim$xid]))
  expect_true(mean(is.na(R[, 3])) > 0.3)

  # utilities equal the index from the truth at every rank (no rank-share correction)
  for (a in seq_len(K_EXP - 1L)) {
    expect_equal(sim[[sprintf("app_r1_Vobs_alt%d", a)]], .exp_utility(sim, a))
    expect_equal(sim[[sprintf("app_r3_Vobs_alt%d", a)]], .exp_utility(sim, a))
  }
  # probabilities sum to one over the available set, all positive at rank 1, and
  # zero at rank 2 for the alternative drawn at rank 1
  p1 <- .exp_prob_mat(sim, 1)
  p2 <- .exp_prob_mat(sim, 2)
  expect_equal(rowSums(p1), rep(1, nrow(sim)))
  expect_equal(rowSums(p2), rep(1, nrow(sim)))
  expect_true(all(p1 > 0))
  expect_true(all(p2[cbind(seq_len(nrow(sim)), R[, 1])] == 0))
  for (r in seq_len(NRANK_EXP)) .exp_check_draws(sim, r)

  .exp_recover(ms, sim)
})

test_that("exclude_chosen = FALSE: alternatives repeat across ranks, recovery", {
  skip_on_cran()
  base <- .exp_base()
  ms <- .exp_system(base, exclude_chosen = FALSE)
  sim <- simulate_factor_model(ms, base, n = 6000, params = .exp_truth, seed = 12)
  R <- .exp_rank_mat(sim)
  expect_true(all(R[!is.na(R)] %in% seq_len(K_EXP)))

  # with five alternatives and three ranks a repeat is common
  full <- !is.na(R[, 2]) & !is.na(R[, 3])
  dup <- apply(R[full, , drop = FALSE], 1, function(v) anyDuplicated(v) > 0)
  expect_true(mean(dup) > 0.3)

  # the rank-1 draw stays available at rank 2, and the same set with the same
  # utilities gives the same probabilities at every rank
  p1 <- .exp_prob_mat(sim, 1)
  p2 <- .exp_prob_mat(sim, 2)
  expect_true(all(p2[cbind(seq_len(nrow(sim)), R[, 1])] > 0))
  expect_equal(p2, p1)
  expect_equal(.exp_prob_mat(sim, 3), p1)
  for (r in seq_len(NRANK_EXP)) .exp_check_draws(sim, r)

  .exp_recover(ms, sim)
})

test_that("rankshare_var: corrections enter the rank utilities by position, recovery", {
  skip_on_cran()
  base <- .exp_base()
  ms <- .exp_system(base, exclude_chosen = FALSE, rankshare = TRUE)
  sim <- simulate_factor_model(ms, base, n = 6000, params = .exp_truth, seed = 13)

  # the 12-column block is carried over, contiguous and in data order
  expect_true(all(.exp_rs_cols %in% names(sim)))
  expect_equal(diff(match(.exp_rs_cols, names(sim))), rep(1L, length(.exp_rs_cols) - 1L))
  for (col in .exp_rs_cols) expect_equal(sim[[col]], base[[col]][sim$xid])

  # rank r utility of alternative a is the index plus the correction, NA read as 0
  for (r in seq_len(NRANK_EXP)) for (a in seq_len(K_EXP - 1L)) {
    corr <- sim[[sprintf("rs_r%d_c%d", r - 1L, a - 1L)]]
    corr[is.na(corr)] <- 0
    expect_equal(sim[[sprintf("app_r%d_Vobs_alt%d", r, a)]], .exp_utility(sim, a) + corr)
  }
  # corrections differ by rank, so the rank probabilities differ
  p1 <- .exp_prob_mat(sim, 1)
  p2 <- .exp_prob_mat(sim, 2)
  expect_false(isTRUE(all.equal(p1, p2)))
  V2 <- cbind(0, sapply(seq_len(K_EXP - 1L), function(a) sim[[sprintf("app_r2_Vobs_alt%d", a)]]))
  expect_equal(p2, exp(V2) / rowSums(exp(V2)))
  for (r in seq_len(NRANK_EXP)) .exp_check_draws(sim, r)

  .exp_recover(ms, sim)
})

test_that("fill_ranks = \"all\" fills every rank, \"base\" needs the rank columns, detail = FALSE", {
  base <- .exp_base(nb = 300)
  ms <- .exp_system(base, exclude_chosen = TRUE)
  sim <- simulate_factor_model(ms, base, n = 500, params = .exp_truth, seed = 14,
                               fill_ranks = "all")
  R <- .exp_rank_mat(sim)
  expect_false(anyNA(R))
  expect_false(any(apply(R, 1, function(v) anyDuplicated(v) > 0)))

  # detail = FALSE keeps the outcomes and drops the per-rank diagnostics
  s2 <- simulate_factor_model(ms, base, n = 50, params = .exp_truth, seed = 1, detail = FALSE)
  expect_true(all(c("r1", "r2", "r3") %in% names(s2)))
  expect_false(any(grepl("^app_r[0-9]+_", names(s2))))

  # without the rank columns in `data`, "base" cannot know which ranks to fill
  base_norank <- base[, setdiff(names(base), c("r1", "r2", "r3"))]
  expect_error(simulate_factor_model(ms, base_norank, n = 50, params = .exp_truth, seed = 1),
               "fill_ranks")
  expect_no_error(simulate_factor_model(ms, base_norank, n = 50, params = .exp_truth, seed = 1,
                                        fill_ranks = "all"))
})

test_that("fill_ranks = \"first\" draws rank 1 only, as the legacy TModel::Sim() does", {
  base <- .exp_base()
  ms <- .exp_system(base, exclude_chosen = TRUE)
  sim <- simulate_factor_model(ms, base, n = 6000, params = .exp_truth, seed = 15,
                               fill_ranks = "first")
  R <- .exp_rank_mat(sim)
  expect_false(anyNA(R[, 1]))
  expect_true(all(is.na(R[, 2:3])))
  expect_true(all(R[, 1] %in% seq_len(K_EXP)))
  # only rank 1 has diagnostic columns
  expect_true(all(sprintf("app_r1_prob_C%d", seq_len(K_EXP)) %in% names(sim)))
  expect_false(any(grepl("^app_r[23]_", names(sim))))
  # the empirical choice frequencies match the multinomial logit probabilities at
  # the true parameters, computed here from the truth, within 3.5 binomial SEs
  V <- cbind(0, sapply(seq_len(K_EXP - 1L), function(a) .exp_utility(sim, a)))
  P <- exp(V) / rowSums(exp(V))
  for (j in seq_len(K_EXP)) {
    z <- (sum(R[, 1] == j) - sum(P[, j])) / sqrt(sum(P[, j] * (1 - P[, j])))
    expect_true(abs(z) < 3.5, info = sprintf("alternative %d: z = %.2f", j, z))
  }
  expect_equal(.exp_prob_mat(sim, 1), P)
  # a re-estimation on these data sees rank 1 only: the first rank column is the
  # outcome a plain multinomial logit would use
  expect_equal(sum(!is.na(R)), nrow(sim))
})

test_that("an evaluation indicator of 0 leaves every rank NA", {
  base <- .exp_base(nb = 300)
  base$ev <- as.integer(runif(300) < 0.7)
  fm <- define_factor_model(n_factors = 1)
  m1 <- define_model_component("m1", base, "y1", fm, covariates = "intercept",
                               model_type = "linear", loading_normalization = 1)
  app <- define_model_component("app", base, c("r1", "r2", "r3"), fm,
                                covariates = c("intercept", "x"), model_type = "logit",
                                num_choices = K_EXP, loading_normalization = NA_real_,
                                evaluation_indicator = "ev")
  ms <- define_model_system(components = list(m1, app), factor = fm)
  truth <- .exp_truth[!grepl("^m[23]_", names(.exp_truth))]
  sim <- simulate_factor_model(ms, base, n = 400, params = truth, seed = 2, fill_ranks = "all")
  off <- sim$app_miss == 1L
  expect_true(any(off) && any(!off))
  expect_true(all(is.na(.exp_rank_mat(sim)[off, ])))
  expect_false(anyNA(.exp_rank_mat(sim)[!off, ]))
})

test_that("a rank column cannot be an endogenous regressor", {
  base <- .exp_base(nb = 300)
  fm <- define_factor_model(n_factors = 1)
  m1 <- define_model_component("m1", base, "y1", fm, covariates = "intercept",
                               model_type = "linear", loading_normalization = 1)
  app <- define_model_component("app", base, c("r1", "r2", "r3"), fm,
                                covariates = c("intercept", "x"), model_type = "logit",
                                num_choices = K_EXP, loading_normalization = NA_real_)
  cons <- define_model_component("cons", base, "y2", fm, covariates = c("intercept", "r1"),
                                 model_type = "linear", loading_normalization = NA_real_,
                                 endogenous = "r1")
  ms <- define_model_system(components = list(m1, app, cons), factor = fm)
  expect_error(simulate_factor_model(ms, base, n = 20, params = .exp_truth, seed = 1),
               "rank column")
})

test_that("a missing rank-share correction is NA and counts as 0 in the likelihood", {
  base <- .exp_base(nb = 400, seed = 5)
  expect_true(anyNA(base[.exp_rs_cols]))
  zero <- base
  for (col in .exp_rs_cols) zero[[col]][is.na(zero[[col]])] <- 0
  eval_at <- function(d) {
    ms <- .exp_system(d, exclude_chosen = FALSE, rankshare = TRUE)
    p <- initialize_parameters(ms, d, verbose = FALSE)$init_params
    p[names(.exp_truth)] <- .exp_truth
    evaluate_likelihood_cpp(initialize_factor_model_cpp(ms, d, 8), p, compute_gradient = TRUE)
  }
  na_coded <- eval_at(base)
  zero_coded <- eval_at(zero)
  expect_true(is.finite(na_coded$logLikelihood))
  expect_identical(na_coded$logLikelihood, zero_coded$logLikelihood)
  expect_identical(na_coded$gradient, zero_coded$gradient)

  # the -9999 code of the legacy data files is not read as missing, and is rejected
  old <- base
  old$rs_r1_c2[3] <- -9999
  expect_error(.exp_system(old, exclude_chosen = FALSE, rankshare = TRUE), "Recode those values to NA")
})
