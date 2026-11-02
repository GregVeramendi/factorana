# T5, R side: education stage with adaptive quadrature. Usage: Rscript t5_r_adapt.R <N> <src: L|R> [cores]
# src L = legacy ML factor scores (factor_predictions.txt), src R = R posterior-mode scores (obs-only SE, as R's pipeline).
args <- commandArgs(trailingOnly = TRUE); N <- as.integer(args[1]); src <- args[2]; cores <- if (length(args) >= 3) as.integer(args[3]) else 8L
suppressMessages(pkgload::load_all("../../factorana_R", compile = FALSE, quiet = TRUE)); source("model_r.R")
dat <- read.csv(sprintf("data/mini_%d.csv", N)); s1 <- readRDS(sprintf("results/r_s1_%d.rds", N))$fit
if (src == "L") {
  fp <- read.table(sprintf("results/legacy_s5_%d/factor_predictions.txt", N), col.names = c("id","f1","se1","f2","se2","f3","se3"))
  fp <- fp[match(dat$id, fp$id), ]; scores <- as.matrix(fp[, c("f1","f2","f3")]); ses <- as.matrix(fp[, c("se1","se2","se3")])
} else {
  fs <- readRDS(sprintf("results/r_scores_%d.rds", N))$obs_se$scores; stopifnot(all(fs$id == dat$id))
  scores <- as.matrix(fs[, c("factor_1","factor_2","factor_3")]); ses <- as.matrix(fs[, c("se_factor_1","se_factor_2","se_factor_3")])
}
ms <- build_mini_system(dat, 3, s1)
ctrl <- define_estimation_control(n_quad_points = 8, num_cores = cores, adaptive_integration = TRUE, adapt_int_thresh = 0.5)
t0 <- Sys.time()
fit <- estimate_model_rcpp(ms, dat, control = ctrl, optimizer = "nlminb", parallel = cores > 1, verbose = FALSE,
                           factor_scores = scores, factor_ses = ses, factor_vars = c(factor_var_1 = 1, factor_var_2 = 1, factor_var_3 = 1))
secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("R stage 3 adaptive (scores %s) N=%d: convergence=%d loglik=%.6f secs=%.1f\n", src, N, fit$convergence, fit$loglik, secs))
out <- data.frame(param = names(fit$estimates), est = unname(fit$estimates), se = unname(fit$std_errors))
saveRDS(list(fit = fit, secs = secs, cores = cores), sprintf("results/r_s3a%s_%d.rds", src, N)); write.csv(out, sprintf("results/r_s3a%s_%d.csv", src, N), row.names = FALSE)
