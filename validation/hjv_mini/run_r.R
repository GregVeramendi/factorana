# Usage: Rscript run_r.R <stage> <N> [cores] [prev_stage] [ntypes] [init: default|truth]. Saves results/r_<tag>_<N>.rds/.csv.
# Stage 5 = education and labor market jointly. ntypes > 1 (T9) adds latent types; tags then end in t<T>, data suffix _t<T>.
# init = "truth" starts the free parameters of this stage at the true values (legacy init_par.txt equivalent); tag suffix "T".
args <- commandArgs(trailingOnly = TRUE)
stage <- as.integer(args[1]); N <- as.integer(args[2]); cores <- if (length(args) >= 3) as.integer(args[3]) else 8L
prev_stage <- if (length(args) >= 4) as.integer(args[4]) else if (stage == 4) 1L else 1L
ntypes <- if (length(args) >= 5) as.integer(args[5]) else 1L; init <- if (length(args) >= 6) args[6] else "default"
ttag <- if (ntypes > 1) sprintf("t%d", ntypes) else ""; dtag <- if (ntypes > 1) sprintf("_t%d", ntypes) else ""
tag <- paste0(if (stage == 4 && prev_stage == 3) "s4p3" else sprintf("s%d", stage), ttag, if (init == "truth") "T" else "")
suppressMessages(pkgload::load_all("../../factorana_R", compile = FALSE, quiet = TRUE))
source("model_r.R"); source("truth_params.R")
dat <- read.csv(sprintf("data/mini_%d%s.csv", N, dtag))
previous <- if (stage >= 3) readRDS(sprintf("results/r_s%d%s_%d.rds", prev_stage, ttag, N))$fit
ms <- build_mini_system(dat, stage, previous, ntypes = ntypes)
ctrl <- define_estimation_control(n_quad_points = 8, num_cores = cores)
init_params <- NULL
if (init == "truth") {   # previous-stage estimates where fixed, true values for this stage's free parameters
  ip <- factorana:::initialize_parameters(ms, dat, verbose = FALSE)
  tr <- readRDS(sprintf("data/truth_%d%s.rds", N, dtag))$truth
  init_params <- truth_r_params(tr, ip$param_names, subset = TRUE)
  if (!is.null(previous)) { keep <- intersect(names(init_params), names(previous$estimates)); init_params[keep] <- previous$estimates[keep] }
}
t0 <- Sys.time()
fit <- estimate_model_rcpp(ms, dat, control = ctrl, optimizer = "nlminb", parallel = cores > 1, verbose = FALSE, init_params = init_params)
secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("R stage %d (previous fixed: stage %d, types %d, init %s) N=%d: convergence=%d loglik=%.6f iterations=%s secs=%.1f cores=%d\n",
            stage, prev_stage, ntypes, init, N, fit$convergence, fit$loglik, paste(fit$iterations, collapse = "/"), secs, cores))
out <- data.frame(param = names(fit$estimates), est = unname(fit$estimates), se = unname(fit$std_errors))
dir.create("results", showWarnings = FALSE)
saveRDS(list(fit = fit, secs = secs, cores = cores, prev_stage = prev_stage), sprintf("results/r_%s_%d.rds", tag, N))
write.csv(out, sprintf("results/r_%s_%d.csv", tag, N), row.names = FALSE)
