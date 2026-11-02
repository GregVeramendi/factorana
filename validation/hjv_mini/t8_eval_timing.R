# T8 addendum: single-core evaluation times of the education-stage system (stage 1 fixed) in full quadrature and in
# adaptive mode, to separate the integration cost from the fixed per-observation cost. Replication 1 data (N = 2001).
# Usage: Rscript t8_eval_timing.R [repeats]
args <- commandArgs(trailingOnly = TRUE); reps <- if (length(args) >= 1) as.integer(args[1]) else 3L
suppressMessages(pkgload::load_all("../../factorana_R", compile = FALSE, quiet = TRUE)); source("model_r.R")
dat <- read.csv("data/mini_2001.csv"); s1 <- readRDS("results/t8/rep1/s1_fit.rds"); ms <- build_mini_system(dat, 3, s1)
init <- factorana:::initialize_parameters(ms, dat, verbose = FALSE); meta <- factorana:::build_parameter_metadata(ms)
cons <- factorana:::setup_parameter_constraints(ms, init$init_params, meta, init$factor_variance_fixed, verbose = FALSE)
p <- init$init_params; pf <- p[cons$free_idx]
post <- readRDS("results/t8/rep1/r_scores_post.rds"); obs <- readRDS("results/t6/rep1/r_scores.rds")
fp <- read.table("results/t6/rep1/legacy_s5_2001/factor_predictions.txt", col.names = c("id", "f1", "se1", "f2", "se2", "f3", "se3")); fp <- fp[match(dat$id, fp$id), ]
sc <- function(fs) as.matrix(fs[, c("factor_1", "factor_2", "factor_3")]); se <- function(fs) as.matrix(fs[, c("se_factor_1", "se_factor_2", "se_factor_3")])
tm <- function(f) median(vapply(seq_len(reps), function(i) system.time(f())[["elapsed"]], numeric(1)))
pts <- function(S, thr) { nq <- pmin(1 + 2 * floor(S / thr), 8); nq[S > 1] <- 8; mean(apply(nq, 1, prod)) }
run <- function(label, scores = NULL, ses = NULL, thr = NA) {
  fm_ptr <- initialize_factor_model_cpp(ms, as.matrix(dat), n_quad = 8L, init_params = p)
  if (!is.null(scores)) set_adaptive_quadrature_cpp(fm_ptr, scores, ses, c(1, 1, 1), thr, 8L, FALSE)
  t_l <- tm(function() evaluate_loglik_only_cpp(fm_ptr, pf)); t_g <- tm(function() evaluate_likelihood_cpp(fm_ptr, pf, TRUE, FALSE)); t_h <- tm(function() evaluate_likelihood_cpp(fm_ptr, pf, TRUE, TRUE))
  cat(sprintf("R %-34s points/person %6.1f  loglik %6.3f s  +grad %6.3f s  +hess %7.3f s\n", label, if (is.null(scores)) 512 else pts(ses, thr), t_l, t_g, t_h))
}
run("full quadrature (8 points/factor)")
run("adaptive, posterior SE, t = 1.0", sc(post), se(post), 1.0)
run("adaptive, posterior SE, t = 0.5", sc(post), se(post), 0.5)
run("adaptive, ML scores (legacy SE), t = 0.5", as.matrix(fp[, c("f1", "f2", "f3")]), as.matrix(fp[, c("se1", "se2", "se3")]), 0.5)
cat("free parameters:", length(pf), " components:", length(ms$components), " N:", nrow(dat), "\n")
