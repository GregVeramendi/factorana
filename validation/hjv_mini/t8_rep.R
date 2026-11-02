# T8 replication helper: cost (points per person, seconds) against accuracy of adaptive rules, posterior-mode scores.
# Reuses the T6 replications (same seeds, legacy stage-1 estimates, full-quadrature reference).
# Usage: Rscript t8_rep.R <N> <rep> <step> [variant]
#  step "postse":  stage-1 fit from the replication's legacy meas_par + posterior SEs (include_prior = TRUE)
#                  -> results/t8/rep<r>/r_scores_post.rds (scores must equal the T6 obs-only run's scores)
#  step "adapt":   R adaptive education fit for a variant "<se>_t<thresh>_m<max>", se in {obs, post}
#                  -> results/t8/rep<r>/r_<variant>.csv with est, se, conv, optimizer seconds, wall seconds
#  step "full":    R full-quadrature education fit (timing reference) -> results/t8/rep<r>/r_full.csv
args <- commandArgs(trailingOnly = TRUE); N <- as.integer(args[1]); rep <- as.integer(args[2]); step <- args[3]
suppressMessages(pkgload::load_all("../../factorana_R", compile = FALSE, quiet = TRUE)); source("model_r.R"); source("legacy_names.R")
dat <- read.csv(sprintf("data/mini_%d.csv", N))
t6 <- sprintf("results/t6/rep%d", rep); out <- sprintf("results/t8/rep%d", rep); dir.create(out, recursive = TRUE, showWarnings = FALSE)
tr_saved <- readRDS(file.path(t6, sprintf("truth_%d.rds", N))); tr_now <- readRDS(sprintf("data/truth_%d.rds", N))
stopifnot(identical(tr_saved$f, tr_now$f))                       # regenerated data equal the T6 replication's data
ctrl <- function(adapt = FALSE, thresh = 0.5, mx = 8L)
  define_estimation_control(n_quad_points = mx, num_cores = 8, adaptive_integration = adapt, adapt_int_thresh = thresh)

stage1_fit <- function() {
  fit <- readRDS(sprintf("results/t6_template_s1_%d.rds", N))
  tab <- legacy_param_table(file.path(t6, sprintf("legacy_s1_%d", N), "run.log"))
  par <- read.table(file.path(t6, sprintf("legacy_s1_%d", N), "meas_par.txt"), col.names = c("idx", "est", "se"))
  leg <- merge(tab[, c("idx", "param")], par, by = "idx")
  stopifnot(all(leg$param %in% names(fit$estimates)))
  fit$estimates[leg$param] <- leg$est; fit$std_errors[leg$param] <- leg$se
  fit
}

if (step == "postse") {
  fit <- stage1_fit(); saveRDS(fit, file.path(out, "s1_fit.rds"))
  t0 <- Sys.time()
  fs <- estimate_factorscores_rcpp(fit, dat, control = ctrl(), parallel = TRUE, verbose = FALSE, include_prior = TRUE, id_var = "id")
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  obs <- readRDS(file.path(t6, "r_scores.rds"))
  gap <- max(abs(as.matrix(fs[, c("factor_1", "factor_2", "factor_3")]) - as.matrix(obs[, c("factor_1", "factor_2", "factor_3")])))
  attr(fs, "seconds") <- secs; saveRDS(fs, file.path(out, "r_scores_post.rds"))
  cat(sprintf("posterior SEs done in %.1f s, converged %.3f, max |score gap| vs T6 obs-only run %.2e\n", secs, mean(fs$converged), gap))
} else if (step %in% c("adapt", "full")) {
  fit <- readRDS(file.path(out, "s1_fit.rds"))
  if (step == "adapt") {
    v <- args[4]; m <- regmatches(v, regexec("^(obs|obsP|post)_t([0-9.]+)_m([0-9]+)$", v))[[1]]; stopifnot(length(m) == 4)
    se_src <- m[2]; thresh <- as.numeric(m[3]); mx <- as.integer(m[4])
    fs <- if (se_src == "post") readRDS(file.path(out, "r_scores_post.rds")) else readRDS(file.path(t6, "r_scores.rds"))
    sc <- as.matrix(fs[, c("factor_1", "factor_2", "factor_3")]); se <- as.matrix(fs[, c("se_factor_1", "se_factor_2", "se_factor_3")])
    # obsP: person-level fallback. A person with any factor SE above the prior SD gets full quadrature on every factor
    # (SE set far above the SD, so the C++ per-factor rule centers all three factors at the prior with max points);
    # no other person then triggers the per-factor fallback, so the rest follow the shared formula at the scores.
    if (se_src == "obsP") { fb <- apply(se > 1, 1, any); se[fb, ] <- 10; cat(sprintf("obsP: full quadrature for %d of %d persons\n", sum(fb), nrow(se))) }
    t0 <- Sys.time()
    f3 <- estimate_model_rcpp(build_mini_system(dat, 3, fit), dat, control = ctrl(TRUE, thresh, mx), optimizer = "nlminb", parallel = TRUE,
                              verbose = FALSE, factor_scores = sc, factor_ses = se, factor_vars = c(factor_var_1 = 1, factor_var_2 = 1, factor_var_3 = 1))
    tag <- v
  } else {
    t0 <- Sys.time()
    f3 <- estimate_model_rcpp(build_mini_system(dat, 3, fit), dat, control = ctrl(), optimizer = "nlminb", parallel = TRUE, verbose = FALSE)
    tag <- "full"
  }
  wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  write.csv(data.frame(param = names(f3$estimates), est = unname(f3$estimates), se = unname(f3$std_errors), conv = f3$convergence,
                       opt_secs = f3$optimization_time, wall_secs = wall, iterations = if (!is.null(f3$iterations)) f3$iterations else NA),
            file.path(out, sprintf("r_%s.csv", tag)), row.names = FALSE)
  cat(sprintf("R %s: convergence=%d, optimizer %.1f s, wall %.1f s\n", tag, f3$convergence, f3$optimization_time, wall))
}
