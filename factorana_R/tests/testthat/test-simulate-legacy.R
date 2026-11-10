# Tests for the legacy layout of simulate_factor_model(): the table the C++
# TMinLkhd::Simulate() / TModel::Sim() wrote to simulation_data.csv, its file
# writer write_legacy_simulation_csv(), and legacy_param_table().
#
# The column order and the file format asserted here were checked against the
# C++ program itself in November 2026: on the mini-HJV model with three types
# (28 components, 344 parameters, 100,000 simulated units) the header of the R
# file equals the header of the C++ file, all 363 columns have the same field
# format, and no column mean differs by more than 2.3 standard errors.
#
# These tests simulate only, so no optimizer runs here.

.leg_base <- function(nb = 500, seed = 3) {
  set.seed(seed)
  d <- data.frame(id = 5000L + seq_len(nb), intercept = 1, x = rnorm(nb))
  d$m1 <- rnorm(nb)
  d$d <- rbinom(nb, 1, 0.5)
  d$ev_d <- as.integer(runif(nb) < 0.8)
  d$d_miss <- 1L - d$ev_d
  d$ch <- sample.int(3, nb, TRUE)
  d$ch2 <- as.integer(d$ch == 2L)
  d$ch3 <- as.integer(d$ch == 3L)
  d$w <- rnorm(nb)
  d$ev_w <- as.integer(runif(nb) < 0.6)
  d$op <- sample.int(3, nb, TRUE)
  d
}

# One factor, two types. Listing order m1, w, ch, d, op: the consumers w and ch
# are listed before their producers, so the simulation order is m1, d, ch, w, op.
.leg_system <- function(base, n_types = 2) {
  fm <- define_factor_model(n_factors = 1, n_types = n_types)
  ut <- n_types > 1
  define_model_system(factor = fm, components = list(
    define_model_component("m1", base, "m1", fm, covariates = "intercept",
                           model_type = "linear", loading_normalization = 1),
    define_model_component("w", base, "w", fm,
                           covariates = c("intercept", "ch2", "ch3", "d", "d_miss"),
                           model_type = "linear", loading_normalization = NA_real_,
                           evaluation_indicator = "ev_w", endogenous = c("ch", "d"),
                           use_types = ut),
    define_model_component("ch", base, "ch", fm, covariates = c("intercept", "x", "d"),
                           model_type = "logit", num_choices = 3,
                           loading_normalization = NA_real_, endogenous = "d",
                           use_types = ut),
    define_model_component("d", base, "d", fm, covariates = c("intercept", "x"),
                           model_type = "probit", loading_normalization = NA_real_,
                           evaluation_indicator = "ev_d", use_types = ut),
    define_model_component("op", base, "op", fm, covariates = "x",
                           model_type = "oprobit", num_choices = 3,
                           loading_normalization = NA_real_)))
}

.leg_truth <- c(
  factor_var_1 = 0.8, typeprob_2_intercept = 0.3, type_2_loading_1 = 0.7,
  m1_intercept = 0.2, m1_sigma = 0.5,
  w_intercept = 1.0, w_ch2 = 0.5, w_ch3 = -0.4, w_d = 0.3, w_d_miss = 2.5,
  w_loading_1 = 0.6, w_sigma = 0.7, w_type_2_intercept = 0.9,
  ch_c1_intercept = 0.2, ch_c1_x = 0.4, ch_c1_d = 0.6, ch_c1_loading_1 = 0.5,
  ch_c1_type_2_intercept = -0.8,
  ch_c2_intercept = -0.3, ch_c2_x = -0.2, ch_c2_d = 0.3, ch_c2_loading_1 = -0.6,
  ch_c2_type_2_intercept = 0.5,
  d_intercept = 0.1, d_x = 0.5, d_loading_1 = 0.8, d_type_2_intercept = -0.4,
  op_x = 0.3, op_loading_1 = 0.5, op_thresh_1 = -0.5, op_thresh_2 = 1.1)

.leg_header <- c(
  "id",
  "sim_TypeProb_miss", "sim_TypeProb_Vobs", "sim_TypeProb_Vfac0", "sim_TypeProb_eps",
  "sim_TypeProb_prob", "sim_TypeProb",
  "sim_m1_miss", "sim_m1_Vobs", "sim_m1_Vfac0", "sim_m1_Vtyp0", "sim_m1_eps", "sim_m1",
  "sim_d_miss", "sim_d_Vobs", "sim_d_Vfac0", "sim_d_Vtyp0", "sim_d_eps", "sim_d_prob", "sim_d",
  "sim_ch_miss",
  "sim_ch_Vobs_C2", "sim_ch_Vend_C2", "sim_ch_Vfac0_C2", "sim_ch_Vtyp0_C2",
  "sim_ch_Vobs_C3", "sim_ch_Vend_C3", "sim_ch_Vfac0_C3", "sim_ch_Vtyp0_C3",
  "sim_ch_eps_C2", "sim_ch_eps_C3", "sim_ch_prob_C2", "sim_ch_prob_C3", "sim_ch",
  "sim_w_miss", "sim_w_Vobs", "sim_w_Vend", "sim_w_Vfac0", "sim_w_Vtyp0", "sim_w_eps", "sim_w",
  "sim_op_miss", "sim_op_Vobs", "sim_op_Vfac0", "sim_op_Vtyp0", "sim_op_eps",
  "sim_op_prob_C2", "sim_op_prob_C3", "sim_op",
  "sim_fac0", "sim_type1", "sim_flag")

.leg_sim <- function(base, ms, n = 3000, seed = 21, ...) {
  as.data.frame(simulate_factor_model(ms, base, n = n, params = .leg_truth, seed = seed,
                                      layout = "legacy", index_var = "id", ...))
}

test_that("the truth vector names every parameter of the test system", {
  base <- .leg_base()
  ms <- .leg_system(base)
  expect_setequal(names(.leg_truth), names(initialize_parameters(ms, base, verbose = FALSE)$init_params))
})

test_that("legacy layout has the columns of the C++ header, in order", {
  base <- .leg_base()
  L <- .leg_sim(base, .leg_system(base))
  expect_identical(names(L), .leg_header)
  expect_equal(nrow(L), 3000)
  # integer columns are those the C++ code writes with %10d
  int_cols <- c("id", grep("_miss$", .leg_header, value = TRUE),
                "sim_TypeProb", "sim_d", "sim_ch", "sim_op", "sim_flag")
  is_int <- vapply(L, is.integer, logical(1))
  expect_setequal(names(L)[is_int], int_cols)
  expect_true(all(L$sim_flag == 1L))
  expect_true(all(L$id %in% base$id))
  expect_equal(L$sim_type1, as.numeric(L$sim_TypeProb == 2L))
})

test_that("each outcome is rebuilt from its detail columns", {
  base <- .leg_base()
  L <- .leg_sim(base, .leg_system(base))
  index <- function(pre, sfx = "") {
    v <- L[[paste0(pre, "_Vobs", sfx)]] + L[[paste0(pre, "_Vfac0", sfx)]] + L[[paste0(pre, "_Vtyp0", sfx)]]
    if (paste0(pre, "_Vend", sfx) %in% names(L)) v <- v + L[[paste0(pre, "_Vend", sfx)]]
    v
  }
  expect_equal(index("sim_m1") + L$sim_m1_eps, L$sim_m1, tolerance = 1e-12)
  expect_equal(index("sim_w") + L$sim_w_eps, L$sim_w, tolerance = 1e-12)
  expect_identical(as.integer(index("sim_d") + L$sim_d_eps > 0), L$sim_d)
  u <- cbind(0, index("sim_ch", "_C2") + L$sim_ch_eps_C2, index("sim_ch", "_C3") + L$sim_ch_eps_C3)
  expect_identical(max.col(u, ties.method = "first"), L$sim_ch)
  # ordered probit: category from the index plus the shock against the thresholds
  tau <- cumsum(c(.leg_truth[["op_thresh_1"]], abs(.leg_truth[["op_thresh_2"]])))
  expect_identical(findInterval(index("sim_op") + L$sim_op_eps, tau) + 1L, L$sim_op)
  # type model: the shock is drawn conditional on the type and reproduces it
  ut <- L$sim_TypeProb_Vobs + L$sim_TypeProb_Vfac0 + L$sim_TypeProb_eps
  expect_identical(as.integer(ut > 0) + 1L, L$sim_TypeProb)
  expect_true(all(L$sim_TypeProb_Vobs == .leg_truth[["typeprob_2_intercept"]]))

  # the pieces are what the C++ code names them
  expect_equal(L$sim_m1_Vobs, rep(.leg_truth[["m1_intercept"]], nrow(L)))
  expect_equal(L$sim_m1_Vfac0, L$sim_fac0)                 # loading fixed at 1
  expect_true(all(L$sim_m1_Vtyp0 == 0))                    # component without types
  expect_equal(L$sim_w_Vtyp0, .leg_truth[["w_type_2_intercept"]] * L$sim_type1)
  expect_equal(L$sim_ch_Vtyp0_C3, .leg_truth[["ch_c2_type_2_intercept"]] * L$sim_type1)
  expect_equal(L$sim_ch_Vend_C2, .leg_truth[["ch_c1_d"]] * L$sim_d)
  expect_equal(L$sim_w_Vend, .leg_truth[["w_ch2"]] * (L$sim_ch == 2L) +
                 .leg_truth[["w_ch3"]] * (L$sim_ch == 3L) + .leg_truth[["w_d"]] * L$sim_d)
  expect_equal(L$sim_w_Vobs, rep(.leg_truth[["w_intercept"]], nrow(L)), tolerance = 1e-12)
})

test_that("probability columns follow the conventions of the C++ file", {
  base <- .leg_base()
  L <- .leg_sim(base, .leg_system(base))
  expect_equal(L$sim_d_prob, pnorm(L$sim_d_Vobs + L$sim_d_Vfac0 + L$sim_d_Vtyp0))
  v2 <- L$sim_ch_Vobs_C2 + L$sim_ch_Vend_C2 + L$sim_ch_Vfac0_C2 + L$sim_ch_Vtyp0_C2
  v3 <- L$sim_ch_Vobs_C3 + L$sim_ch_Vend_C3 + L$sim_ch_Vfac0_C3 + L$sim_ch_Vtyp0_C3
  den <- 1 + exp(v2) + exp(v3)
  expect_equal(L$sim_ch_prob_C2, exp(v2) / den)
  expect_equal(L$sim_ch_prob_C3, exp(v3) / den)
  # ordered probit: the columns C2 and C3 hold the probabilities of categories 1 and 2
  v <- L$sim_op_Vobs + L$sim_op_Vfac0
  tau <- cumsum(c(.leg_truth[["op_thresh_1"]], abs(.leg_truth[["op_thresh_2"]])))
  expect_equal(L$sim_op_prob_C2, pnorm(tau[1] - v))
  expect_equal(L$sim_op_prob_C3, pnorm(tau[2] - v) - pnorm(tau[1] - v))
  # the realized share of category 1 is within 4 standard errors of its mean probability
  p1 <- L$sim_op_prob_C2
  expect_lt(abs(mean(L$sim_op == 1L) - mean(p1)), 4 * sqrt(mean(p1 * (1 - p1)) / nrow(L)))
})

test_that("every component is simulated for every unit, with the legacy missing codes", {
  base <- .leg_base()
  ms <- .leg_system(base)
  L <- .leg_sim(base, ms)
  # code 1: the base row does not evaluate the component, the outcome is simulated anyway
  expect_setequal(unique(L$sim_w_miss), c(0L, 1L))
  expect_identical(L$sim_w_miss, 1L - base$ev_w[match(L$id, base$id)])
  expect_false(anyNA(L$sim_w))
  expect_false(anyNA(L$sim_d))
  # the missing indicator of an endogenous producer contributes nothing: with
  # w_d_miss = 2.5 a data value of 1 would move sim_w_Vobs away from the intercept
  expect_equal(L$sim_w_Vobs, rep(.leg_truth[["w_intercept"]], nrow(L)), tolerance = 1e-12)

  # code 2: a covariate is missing (NA), the block is missing, and so are the
  # consumers of the missing outcome
  bad <- base
  bad$x[seq(1, nrow(bad), by = 5)] <- NA
  B <- .leg_sim(bad, ms)
  lost <- !(B$id %in% bad$id[!is.na(bad$x)])
  expect_true(any(lost))
  for (cmp in c("d", "ch", "w", "op")) {
    expect_true(all(B[[paste0("sim_", cmp, "_miss")]][lost] == 2L), info = cmp)
    block <- grep(sprintf("^sim_%s(_|$)", cmp), names(B), value = TRUE)
    block <- setdiff(block, paste0("sim_", cmp, "_miss"))
    expect_true(all(is.na(as.matrix(B[lost, block]))), info = cmp)
    expect_false(anyNA(as.matrix(B[!lost, block])), info = cmp)
  }
  expect_true(all(B$sim_m1_miss == 0L))      # m1 does not use x
  expect_false(anyNA(B$sim_m1))

  # only NA is missing: the -9999 code of the legacy data files would enter the
  # index as a number, so it stops the simulation
  old <- base
  old$x[7] <- -9999
  expect_error(.leg_sim(old, ms), "Recode those values to NA")
})

test_that("the legacy layout holds the same draws as the standard layout", {
  base <- .leg_base()
  ms <- .leg_system(base)
  L <- .leg_sim(base, ms, seed = 77)
  S <- simulate_factor_model(ms, base, n = 3000, params = .leg_truth, seed = 77,
                             simulate_unevaluated = TRUE)
  expect_identical(S$m1, L$sim_m1)
  expect_identical(S$d, L$sim_d)
  expect_identical(S$ch, L$sim_ch)
  expect_identical(S$w, L$sim_w)
  expect_identical(S$op, L$sim_op)
  expect_identical(S$type, L$sim_TypeProb)
  expect_identical(S$factor_1, L$sim_fac0)
  expect_identical(base$id[S$xid], L$id)
  # the standard layout with simulate_unevaluated = TRUE reports the same codes
  expect_identical(S$w_miss, L$sim_w_miss)
})

test_that("a ranked logit component is written as one draw over all alternatives", {
  base <- .leg_base()
  base$r1 <- sample.int(4, nrow(base), TRUE)
  base$r2 <- sample.int(4, nrow(base), TRUE)
  fm <- define_factor_model(n_factors = 1)
  ms <- define_model_system(factor = fm, components = list(
    define_model_component("m1", base, "m1", fm, covariates = "intercept",
                           model_type = "linear", loading_normalization = 1),
    define_model_component("app", base, c("r1", "r2"), fm, covariates = c("intercept", "x"),
                           model_type = "logit", num_choices = 4,
                           loading_normalization = NA_real_)))
  p <- c(factor_var_1 = 1, m1_intercept = 0, m1_sigma = 0.5,
         app_c1_intercept = 0.3, app_c1_x = 0.5, app_c1_loading_1 = 0.4,
         app_c2_intercept = -0.2, app_c2_x = -0.3, app_c2_loading_1 = 0.8,
         app_c3_intercept = 0.1, app_c3_x = 0.2, app_c3_loading_1 = -0.5)
  L <- as.data.frame(simulate_factor_model(ms, base, n = 4000, params = p, seed = 5,
                                           layout = "legacy"))
  expect_identical(names(L)[1], "xid")
  app_cols <- grep("^sim_app", names(L), value = TRUE)
  expect_identical(app_cols, c(
    "sim_app_miss",
    "sim_app_Vobs_C2", "sim_app_Vfac0_C2", "sim_app_Vobs_C3", "sim_app_Vfac0_C3",
    "sim_app_Vobs_C4", "sim_app_Vfac0_C4",
    "sim_app_eps_C2", "sim_app_eps_C3", "sim_app_eps_C4",
    "sim_app_prob_C2", "sim_app_prob_C3", "sim_app_prob_C4", "sim_app"))
  u <- cbind(0, sapply(2:4, function(j)
    L[[sprintf("sim_app_Vobs_C%d", j)]] + L[[sprintf("sim_app_Vfac0_C%d", j)]] +
      L[[sprintf("sim_app_eps_C%d", j)]]))
  expect_identical(max.col(u, ties.method = "first"), L$sim_app)
  p3 <- L$sim_app_prob_C3
  expect_lt(abs(mean(L$sim_app == 3L) - mean(p3)), 4 * sqrt(mean(p3 * (1 - p3)) / nrow(L)))
  # the other fill rules write several rank columns, which the C++ file does not have
  expect_error(simulate_factor_model(ms, base, params = p, layout = "legacy", fill_ranks = "all"),
               "fill_ranks")
})

test_that("write_legacy_simulation_csv() writes the C++ file format", {
  tab <- data.frame(id = c(7L, 12L), sim_a_miss = c(0L, 2L), sim_a_Vobs = c(1.5, NA),
                    sim_a = c(-0.123456, NA), sim_b = c(1L, NA), sim_flag = c(1L, 1L))
  f <- tempfile(fileext = ".csv")
  on.exit(unlink(f))
  write_legacy_simulation_csv(tab, f)
  expect_identical(readLines(f), c(
    "id, sim_a_miss, sim_a_Vobs, sim_a, sim_b, sim_flag",
    "         7,          0,    1.50000,   -0.12346,          1,          1 ",
    "        12,          2,         NA,         NA,         NA,          1 "))
  # chunking does not change the file
  f2 <- tempfile(fileext = ".csv")
  on.exit(unlink(f2), add = TRUE)
  write_legacy_simulation_csv(tab, f2, chunk_rows = 1L)
  expect_identical(readLines(f2), readLines(f))
  expect_error(write_legacy_simulation_csv(data.frame(a = "x"), f), "non-numeric")
})

test_that("csv_file and rds_file save the simulated table", {
  base <- .leg_base()
  base$x[seq(1, nrow(base), by = 10)] <- NA
  ms <- .leg_system(.leg_base())
  csv <- tempfile(fileext = ".csv")
  rds <- tempfile(fileext = ".rds")
  on.exit(unlink(c(csv, rds)))
  L <- .leg_sim(base, ms, n = 800, csv_file = csv, rds_file = rds)

  expect_identical(as.data.frame(readRDS(rds)), L)

  lines <- readLines(csv)
  expect_length(lines, 801)
  expect_identical(lines[1], paste(.leg_header, collapse = ", "))
  expect_true(all(grepl(" $", lines[-1])))
  expect_length(unique(nchar(lines[-1])), 1L)          # fixed-width fields
  R <- utils::read.csv(csv, strip.white = TRUE, check.names = FALSE)
  expect_identical(names(R), .leg_header)
  expect_true(anyNA(L))
  expect_identical(is.na(as.matrix(L)), is.na(as.matrix(R)), ignore_attr = TRUE)
  expect_equal(as.matrix(R), round(as.matrix(L), 5), tolerance = 1e-9, ignore_attr = TRUE)

  # the standard layout can be saved as a compressed R file, but not as the C++ file
  rds2 <- tempfile(fileext = ".rds")
  on.exit(unlink(rds2), add = TRUE)
  S <- simulate_factor_model(ms, .leg_base(), n = 200, params = .leg_truth, seed = 1, rds_file = rds2)
  expect_identical(readRDS(rds2), S)
  expect_error(simulate_factor_model(ms, .leg_base(), params = .leg_truth, csv_file = csv), "csv_file")
})

test_that("detail and component_order select and order the blocks", {
  base <- .leg_base()
  ms <- .leg_system(base)
  L <- .leg_sim(base, ms, n = 300, detail = c("d", "TypeProb"))
  expect_true(all(c("sim_d_Vobs", "sim_d_eps", "sim_TypeProb_eps") %in% names(L)))
  expect_identical(grep("^sim_w", names(L), value = TRUE), c("sim_w_miss", "sim_w"))
  expect_identical(grep("^sim_ch", names(L), value = TRUE),
                   c("sim_ch_miss", "sim_ch_prob_C2", "sim_ch_prob_C3", "sim_ch"))
  N <- .leg_sim(base, ms, n = 300, detail = FALSE)
  expect_identical(names(N), c(
    "id", "sim_TypeProb_miss", "sim_TypeProb_prob", "sim_TypeProb",
    "sim_m1_miss", "sim_m1", "sim_d_miss", "sim_d_prob", "sim_d",
    "sim_ch_miss", "sim_ch_prob_C2", "sim_ch_prob_C3", "sim_ch", "sim_w_miss", "sim_w",
    "sim_op_miss", "sim_op_prob_C2", "sim_op_prob_C3", "sim_op",
    "sim_fac0", "sim_type1", "sim_flag"))
  expect_error(.leg_sim(base, ms, detail = "nosuch"), "names no component")

  O <- .leg_sim(base, ms, n = 300, detail = FALSE, component_order = c("op", "d"))
  blocks <- unique(sub("^sim_([A-Za-z0-9]+?)(_miss|_prob.*)?$", "\\1", grep("^sim_", names(O), value = TRUE), perl = TRUE))
  expect_identical(blocks, c("TypeProb", "op", "d", "m1", "ch", "w", "fac0", "type1", "flag"))
  expect_error(.leg_sim(base, ms, component_order = c("w", "d")), "before its endogenous producer")
  expect_error(.leg_sim(base, ms, component_order = "nosuch"), "component_order")
})

test_that("legacy_param_table() uses the model and variable names of the C++ list", {
  base <- .leg_base()
  ms <- .leg_system(base)
  pt <- as.data.frame(legacy_param_table(ms, .leg_truth))
  expect_identical(names(pt), c("V1", "model", "variable", "estimate", "se", "parameter"))
  expect_identical(pt$V1, seq_len(nrow(pt)) - 1L)
  expect_setequal(pt$parameter, names(.leg_truth))
  get <- function(m, v) pt$estimate[pt$model == m & pt$variable == v]
  expect_equal(get("factor0_mix0", "sigma"), sqrt(0.8))
  expect_equal(get("TypeProb", "constant"), 0.3)
  expect_equal(get("TypeProb", "factor1"), 0.7)
  expect_equal(get("m1", "intercept"), 0.2)
  expect_equal(get("m1", "sigma"), 0.5)
  expect_length(get("m1", "factor1"), 0L)               # loading fixed by normalization
  expect_equal(get("w", "d_miss"), 2.5)
  expect_equal(get("w", "factor1"), 0.6)
  expect_equal(get("w", "type1"), 0.9)
  expect_equal(get("d", "x"), 0.5)
  expect_equal(get("op", "threshold1"), -0.5)
  expect_equal(get("op", "threshold2"), 1.1)
  # three alternatives: the C++ list repeats the plain model name, in the order of the alternatives
  expect_equal(get("ch", "x"), c(0.4, -0.2))
  expect_equal(get("ch", "type1"), c(-0.8, 0.5))
  expect_identical(legacy_param_table(ms, .leg_truth, type_intercept_name = "intercept")$variable[2], "intercept")

  # more than three alternatives: <component>_ch<j>, and a fixed coefficient keeps its value
  base$r4 <- sample.int(4, nrow(base), TRUE)
  fm <- define_factor_model(n_factors = 1)
  c4 <- define_model_component("c4", base, "r4", fm, covariates = c("intercept", "x"),
                               model_type = "logit", num_choices = 4, loading_normalization = NA_real_)
  c4 <- fix_coefficient(c4, "x", 0.25, choice = 3)
  ms4 <- define_model_system(factor = fm, components = list(
    define_model_component("m1", base, "m1", fm, covariates = "intercept",
                           model_type = "linear", loading_normalization = 1), c4))
  p4 <- c(factor_var_1 = 4, m1_intercept = 0, m1_sigma = 1,
          c4_c1_intercept = 1, c4_c1_x = 2, c4_c1_loading_1 = 3,
          c4_c2_intercept = 4, c4_c2_x = 5, c4_c2_loading_1 = 6,
          c4_c3_intercept = 7, c4_c3_loading_1 = 9)
  t4 <- as.data.frame(legacy_param_table(ms4, p4))
  expect_equal(t4$estimate[t4$model == "factor0_mix0"], 2)
  expect_equal(t4$estimate[t4$model == "c4_ch3" & t4$variable == "x"], 5)
  expect_equal(t4$estimate[t4$model == "c4_ch4" & t4$variable == "x"], 0.25)
  expect_equal(t4$estimate[t4$model == "c4_ch4" & t4$variable == "factor1"], 9)
})

test_that("a multinomial logit with type intercepts is simulated in the standard layout", {
  base <- .leg_base()
  ms <- .leg_system(base)
  S <- simulate_factor_model(ms, base, n = 2000, params = .leg_truth, seed = 4)
  v <- .leg_truth[["ch_c2_intercept"]] + .leg_truth[["ch_c2_x"]] * S$x +
    .leg_truth[["ch_c2_d"]] * S$d + .leg_truth[["ch_c2_loading_1"]] * S$factor_1 +
    .leg_truth[["ch_c2_type_2_intercept"]] * (S$type == 2L)
  ok <- !is.na(S$d)
  expect_equal(S$ch_Vobs_alt2[ok], v[ok], tolerance = 1e-12)
})

test_that("only a declared consumer reads the simulated value of a producer", {
  base <- .leg_base()
  base$u <- rnorm(nrow(base))
  fm <- define_factor_model(n_factors = 1)
  mk <- function(endogenous) define_model_system(factor = fm, components = list(
    define_model_component("m1", base, "m1", fm, covariates = "intercept",
                           model_type = "linear", loading_normalization = 1),
    define_model_component("d", base, "d", fm, covariates = c("intercept", "x"),
                           model_type = "probit", loading_normalization = NA_real_),
    define_model_component("u", base, "u", fm, covariates = c("intercept", "d"),
                           model_type = "linear", loading_normalization = NA_real_,
                           endogenous = endogenous)))
  p <- c(factor_var_1 = 1, m1_intercept = 0, m1_sigma = 0.5,
         d_intercept = 0.2, d_x = 0.6, d_loading_1 = 0.7,
         u_intercept = 0.1, u_d = 1.5, u_loading_1 = 0.3, u_sigma = 0.4)
  D <- simulate_factor_model(mk("d"), base, n = 1500, params = p, seed = 8)
  U <- simulate_factor_model(mk(NULL), base, n = 1500, params = p, seed = 8)
  expect_identical(D$d, U$d)                               # same draws of the producer
  expect_false(identical(D$d, base$d[D$xid]))
  expect_equal(D$u_Vobs, 0.1 + 1.5 * D$d + 0.3 * D$factor_1, tolerance = 1e-12)
  expect_equal(U$u_Vobs, 0.1 + 1.5 * base$d[U$xid] + 0.3 * U$factor_1, tolerance = 1e-12)
})

test_that("a model declared with types that no component uses is simulated without types", {
  base <- .leg_base()
  fm <- define_factor_model(n_factors = 1, n_types = 3)
  ms <- define_model_system(factor = fm, components = list(
    define_model_component("m1", base, "m1", fm, covariates = "intercept",
                           model_type = "linear", loading_normalization = 1),
    define_model_component("d", base, "d", fm, covariates = c("intercept", "x"),
                           model_type = "probit", loading_normalization = NA_real_)))
  p <- c(factor_var_1 = 1, m1_intercept = 0, m1_sigma = 0.5,
         d_intercept = 0.2, d_x = 0.6, d_loading_1 = 0.7)
  L <- as.data.frame(simulate_factor_model(ms, base, n = 200, params = p, seed = 2, layout = "legacy"))
  expect_identical(names(L), c("xid", "sim_m1_miss", "sim_m1_Vobs", "sim_m1_Vfac0", "sim_m1_eps", "sim_m1",
                               "sim_d_miss", "sim_d_Vobs", "sim_d_Vfac0", "sim_d_eps", "sim_d_prob", "sim_d",
                               "sim_fac0", "sim_flag"))
  S <- simulate_factor_model(ms, base, n = 200, params = p, seed = 2)
  expect_false("type" %in% names(S))
})

test_that("endogenous_bins feeds the threshold bins of a linear producer to its consumers", {
  base <- .leg_base()
  base$m12 <- rbinom(nrow(base), 1, 0.3)
  base$m13 <- (1L - base$m12) * rbinom(nrow(base), 1, 0.4)
  base$u <- rnorm(nrow(base))
  fm <- define_factor_model(n_factors = 1)
  ms <- define_model_system(factor = fm, components = list(
    define_model_component("u", base, "u", fm, covariates = c("intercept", "m12", "m13"),
                           model_type = "linear", loading_normalization = NA_real_,
                           endogenous = "m1"),
    define_model_component("m1", base, "m1", fm, covariates = "intercept",
                           model_type = "linear", loading_normalization = 1)))
  p <- c(factor_var_1 = 1, m1_intercept = 0.1, m1_sigma = 0.6,
         u_intercept = 0.2, u_m12 = 1.5, u_m13 = -2, u_loading_1 = 0.3, u_sigma = 0.4)
  thr <- c(-0.3, 0.6)
  # without thresholds a linear producer drives only a covariate of its own name
  expect_error(simulate_factor_model(ms, base, params = p), "none of its covariates")
  S <- simulate_factor_model(ms, base, n = 2000, params = p, seed = 6, endogenous_bins = list(m1 = thr))
  bin <- 1L + (S$m1 > thr[1]) + (S$m1 > thr[2])        # the C++ rule: strictly above a threshold
  expect_setequal(unique(bin), 1:3)
  expect_identical(S$m12, as.integer(bin == 2L))
  expect_identical(S$m13, as.integer(bin == 3L))
  expect_equal(S$u_Vobs, 0.2 + 1.5 * (bin == 2L) - 2 * (bin == 3L) + 0.3 * S$factor_1, tolerance = 1e-12)
  L <- as.data.frame(simulate_factor_model(ms, base, n = 2000, params = p, seed = 6, layout = "legacy",
                                           endogenous_bins = list(m1 = thr)))
  expect_identical(L$sim_m1, S$m1)
  expect_equal(L$sim_u_Vend, 1.5 * (bin == 2L) - 2 * (bin == 3L))
  expect_equal(L$sim_u_Vobs, rep(0.2, nrow(L)), tolerance = 1e-12)
  expect_error(simulate_factor_model(ms, base, params = p, endogenous_bins = list(m1 = c(1, 0))),
               "endogenous_bins")
  expect_error(simulate_factor_model(ms, base, params = p, endogenous_bins = list(nosuch = 0)),
               "endogenous_bins")
})

test_that("the first component with an outcome variable is its producer", {
  base <- .leg_base()
  base$u <- rnorm(nrow(base))
  fm <- define_factor_model(n_factors = 1)
  ms <- define_model_system(factor = fm, components = list(
    define_model_component("m1", base, "m1", fm, covariates = "intercept",
                           model_type = "linear", loading_normalization = 1),
    define_model_component("da", base, "d", fm, covariates = "intercept",
                           model_type = "probit", loading_normalization = NA_real_),
    define_model_component("db", base, "d", fm, covariates = c("intercept", "x"),
                           model_type = "probit", loading_normalization = NA_real_),
    define_model_component("u", base, "u", fm, covariates = c("intercept", "d"),
                           model_type = "linear", loading_normalization = NA_real_,
                           endogenous = "d")))
  p <- c(factor_var_1 = 1, m1_intercept = 0, m1_sigma = 0.5,
         da_intercept = 0.3, da_loading_1 = 0.5, db_intercept = -0.8, db_x = 1, db_loading_1 = -0.5,
         u_intercept = 0, u_d = 1, u_loading_1 = 0, u_sigma = 0.1)
  L <- as.data.frame(simulate_factor_model(ms, base, n = 2000, params = p, seed = 12, layout = "legacy"))
  expect_false(identical(L$sim_da, L$sim_db))
  expect_equal(L$sim_u_Vend, as.numeric(L$sim_da))     # the first listed component, as in the C++ code
})
