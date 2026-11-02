# T8 at scale, R side. Usage: Rscript t8_big_r.R <N> <step> [args]
#  scores                      posterior-mode scores with posterior SE from results/r_s1_<N>.rds -> results/t8big/r_scores_post_<N>.rds
#                              and factor_predictions_P_<N>.txt (legacy format) in results/t8big/
#  adapt <se> <thresh> <cores> education stage (stage 1 fixed) with adaptive integration; se = post | ml (ml reads
#                              results/legacy_s5_<N>/factor_predictions.txt); appends a line to results/t8big/r_fits_<N>.csv
#  evaltime <repeats>          single-core evaluation times: full, post t 0.5, post t 1.0, ml t 0.5
args <- commandArgs(trailingOnly = TRUE); N <- as.integer(args[1]); step <- args[2]
suppressMessages(pkgload::load_all("../../factorana_R", compile = FALSE, quiet = TRUE)); source("model_r.R")
dat <- read.csv(sprintf("data/mini_%d.csv", N)); dir.create("results/t8big", showWarnings = FALSE)
s1 <- readRDS(sprintf("results/r_s1_%d.rds", N))$fit
sc <- function(fs) as.matrix(fs[, c("factor_1", "factor_2", "factor_3")]); se <- function(fs) as.matrix(fs[, c("se_factor_1", "se_factor_2", "se_factor_3")])
ml <- function() { fp <- read.table(sprintf("results/legacy_s5_%d/factor_predictions.txt", N), col.names = c("id", "f1", "se1", "f2", "se2", "f3", "se3"))
  fp <- fp[match(dat$id, fp$id), ]; list(sc = as.matrix(fp[, c("f1", "f2", "f3")]), se = as.matrix(fp[, c("se1", "se2", "se3")])) }
if (step == "scores") {
  t0 <- Sys.time(); fs <- estimate_factorscores_rcpp(s1, dat, control = define_estimation_control(n_quad_points = 8, num_cores = 8), parallel = TRUE, verbose = FALSE, include_prior = TRUE, id_var = "id")
  cat(sprintf("R posterior scores N=%d: %.1f s, converged %.3f\n", N, as.numeric(difftime(Sys.time(), t0, units = "secs")), mean(fs$converged)))
  saveRDS(fs, sprintf("results/t8big/r_scores_post_%d.rds", N))
  write.table(data.frame(fs$id, fs$factor_1, fs$se_factor_1, fs$factor_2, fs$se_factor_2, fs$factor_3, fs$se_factor_3), sprintf("results/t8big/factor_predictions_P_%d.txt", N), row.names = FALSE, col.names = FALSE, quote = FALSE)
} else if (step == "adapt") {
  src <- args[3]; thr <- as.numeric(args[4]); cores <- as.integer(args[5])
  S <- if (src == "post") { fs <- readRDS(sprintf("results/t8big/r_scores_post_%d.rds", N)); list(sc = sc(fs), se = se(fs)) } else ml()
  ctrl <- define_estimation_control(n_quad_points = 8, num_cores = cores, adaptive_integration = TRUE, adapt_int_thresh = thr)
  t0 <- Sys.time(); f <- estimate_model_rcpp(build_mini_system(dat, 3, s1), dat, control = ctrl, optimizer = "nlminb", parallel = cores > 1, verbose = FALSE,
                                             factor_scores = S$sc, factor_ses = S$se, factor_vars = c(factor_var_1 = 1, factor_var_2 = 1, factor_var_3 = 1))
  wall <- as.numeric(difftime(Sys.time(), t0, units = "secs")); nq <- pmin(1 + 2 * floor(S$se / thr), 8); nq[S$se > 1] <- 8; pts <- mean(apply(nq, 1, prod))
  line <- sprintf("%s,%g,%d,%d,%d,%.2f,%.1f,%.2f,%.6f", src, thr, cores, f$convergence, f$iterations, f$optimization_time, wall, pts, f$loglik)
  cat("R adaptive", line, "\n"); cat(line, "\n", file = sprintf("results/t8big/r_fits_%d.csv", N), append = TRUE)
} else if (step == "evaltime") {
  reps <- if (length(args) >= 3) as.integer(args[3]) else 1L
  ms <- build_mini_system(dat, 3, s1); init <- factorana:::initialize_parameters(ms, dat, verbose = FALSE); meta <- factorana:::build_parameter_metadata(ms)
  cons <- factorana:::setup_parameter_constraints(ms, init$init_params, meta, init$factor_variance_fixed, verbose = FALSE); p <- init$init_params; pf <- p[cons$free_idx]
  tm <- function(f) median(vapply(seq_len(reps), function(i) system.time(f())[["elapsed"]], numeric(1)))
  run <- function(label, S = NULL, thr = NA) { fm_ptr <- initialize_factor_model_cpp(ms, as.matrix(dat), n_quad = 8L, init_params = p)
    if (!is.null(S)) set_adaptive_quadrature_cpp(fm_ptr, S$sc, S$se, c(1, 1, 1), thr, 8L, FALSE)
    pts <- if (is.null(S)) 512 else { nq <- pmin(1 + 2 * floor(S$se / thr), 8); nq[S$se > 1] <- 8; mean(apply(nq, 1, prod)) }
    cat(sprintf("R eval N=%d %-30s points/person %6.1f  loglik %7.3f s  +grad %7.3f s  +hess %8.3f s\n", N, label, pts, tm(function() evaluate_loglik_only_cpp(fm_ptr, pf)),
                tm(function() evaluate_likelihood_cpp(fm_ptr, pf, TRUE, FALSE)), tm(function() evaluate_likelihood_cpp(fm_ptr, pf, TRUE, TRUE)))) }
  fs <- readRDS(sprintf("results/t8big/r_scores_post_%d.rds", N)); P <- list(sc = sc(fs), se = se(fs))
  run("adaptive, posterior SE, t = 1.0", P, 1.0); run("adaptive, posterior SE, t = 0.5", P, 0.5)
  if (file.exists(sprintf("results/legacy_s5_%d/factor_predictions.txt", N))) run("adaptive, ML scores, t = 0.5", ml(), 0.5)
  run("full quadrature")
}
