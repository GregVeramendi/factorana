# T4, R side: factor scores from the Stage 1 fit. Usage: Rscript t4_r_scores.R <N> <cores>
# Writes results/r_scores_<N>.rds with both SE definitions (include_prior FALSE and TRUE).
args <- commandArgs(trailingOnly = TRUE); N <- as.integer(args[1]); cores <- if (length(args) >= 2) as.integer(args[2]) else 8L
suppressMessages(pkgload::load_all("../../factorana_R", compile = FALSE, quiet = TRUE))
dat <- read.csv(sprintf("data/mini_%d.csv", N)); fit <- readRDS(sprintf("results/r_s1_%d.rds", N))$fit
ctrl <- define_estimation_control(n_quad_points = 8, num_cores = cores)
out <- list()
for (ip in c(FALSE, TRUE)) {
  t0 <- Sys.time()
  fs <- estimate_factorscores_rcpp(fit, dat, control = ctrl, parallel = cores > 1, verbose = FALSE, include_prior = ip, id_var = "id")
  out[[if (ip) "prior_se" else "obs_se"]] <- list(scores = fs, secs = as.numeric(difftime(Sys.time(), t0, units = "secs")))
  cat(sprintf("R scores include_prior=%s: %d rows, %.1f s, converged %.3f\n", ip, nrow(fs), out[[length(out)]]$secs, mean(fs$converged)))
}
saveRDS(out, sprintf("results/r_scores_%d.rds", N))
