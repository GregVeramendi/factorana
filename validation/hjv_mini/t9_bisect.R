# T9 bisection: where do the typed log-likelihoods differ? Evaluates R and legacy at variants of the true parameter
# vector in which subsets of the type-related parameters are set to zero. Legacy runs TestCode with LKHD_ONLY=1 in
# results/legacy_s0t<T>_<N> (init_par.txt rewritten per variant; the truth file is restored at the end).
# Usage: Rscript t9_bisect.R <N> <T>
args <- commandArgs(trailingOnly = TRUE); N <- as.integer(args[1]); T <- as.integer(args[2])
suppressMessages(pkgload::load_all("../../factorana_R", compile = FALSE, quiet = TRUE)); source("model_r.R"); source("truth_params.R")
dtag <- sprintf("_t%d", T); d0 <- sprintf("results/legacy_s0t%d_%d", T, N)
dat <- read.csv(sprintf("data/mini_%d%s.csv", N, dtag)); tr <- readRDS(sprintf("data/truth_%d%s.rds", N, dtag))$truth
fm <- mini_factor_model(T)
ms <- define_model_system(components = c(build_mini_system(dat, 1, ntypes = T)$components, build_mini_system(dat, 5, ntypes = T)$components), factor = fm)
init <- factorana:::initialize_parameters(ms, dat, verbose = FALSE); p0 <- truth_r_params(tr, init$param_names, subset = TRUE)
meta <- factorana:::build_parameter_metadata(ms); cons <- factorana:::setup_parameter_constraints(ms, p0, meta, init$factor_variance_fixed, verbose = FALSE)
map <- read.csv(file.path(d0, "legacy_index_map.csv")); stopifnot(all(map$param %in% names(p0)))
file.copy(file.path(d0, "init_par.txt"), file.path(d0, "init_par_truth.txt"), overwrite = TRUE)
const <- N * 3 * 0.5 * log(pi)
is_tp <- grepl("^typeprob_[0-9]+_intercept$", names(p0)); is_tl <- grepl("^type_[0-9]+_loading_[0-9]+$", names(p0))
is_ci <- grepl("_type_[0-9]+_intercept$", names(p0))
comp_of <- function(nm) sub("_(c[0-9]+_)?type_[0-9]+_intercept$", "", nm)
lin <- c("swesat", "wage1", "wage2"); lgt <- c("HS_track", "application", "switch_major")
variants <- list(
  truth = rep(TRUE, length(p0)),
  no_type_params = !(is_tp | is_tl | is_ci),
  typeprob_intercepts_only = !(is_tl | is_ci),
  type_loadings_only = !(is_tp | is_ci),
  component_type_intercepts_only = !(is_tp | is_tl),
  linear_type_intercepts_only = !(is_tp | is_tl) & !(is_ci & !(comp_of(names(p0)) %in% lin)),
  logit_type_intercepts_only = !(is_tp | is_tl) & !(is_ci & !(comp_of(names(p0)) %in% lgt)),
  probit_type_intercepts_only = !(is_tp | is_tl) & !(is_ci & (comp_of(names(p0)) %in% c(lin, lgt))),
  advmath_type_intercepts_only = !(is_tp | is_tl) & !(is_ci & comp_of(names(p0)) != "advmath"))
if (length(args) >= 3 && args[3] == "logit") variants <- list(   # one logit component at a time, plus first-type-only
  HS_track_only = !(is_tp | is_tl) & !(is_ci & comp_of(names(p0)) != "HS_track"),
  application_only = !(is_tp | is_tl) & !(is_ci & comp_of(names(p0)) != "application"),
  switch_major_only = !(is_tp | is_tl) & !(is_ci & comp_of(names(p0)) != "switch_major"),
  HS_track_c1_type2_only = !(is_tp | is_tl) & !(is_ci & names(p0) != "HS_track_c1_type_2_intercept"),
  HS_track_c3_type3_only = !(is_tp | is_tl) & !(is_ci & names(p0) != "HS_track_c3_type_3_intercept"))
res <- list()
for (v in names(variants)) {
  p <- p0; p[!variants[[v]]] <- 0
  fm_ptr <- initialize_factor_model_cpp(ms, as.matrix(dat), n_quad = 8L, init_params = p)
  ll_r <- evaluate_loglik_only_cpp(fm_ptr, p[cons$free_idx])
  write.table(data.frame(map$idx, sprintf("%.12f", p[map$param]), 0.1), file.path(d0, "init_par.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
  rc <- system(sprintf("LKHD_ONLY=1 ./prob_t9_s0t%d > %s/bisect_%s.log 2>&1", T, d0, v))
  obj_l <- scan(file.path(d0, "t1_dump.txt"), quiet = TRUE)[1]; ll_l <- -obj_l - const
  res[[v]] <- data.frame(variant = v, nonzero_type_params = sum(variants[[v]] & (is_tp | is_tl | is_ci)), loglik_R = ll_r, loglik_legacy = ll_l, gap = ll_r - ll_l, legacy_rc = rc)
  cat(sprintf("%-32s R %16.6f  legacy %16.6f  gap %12.6f\n", v, ll_r, ll_l, ll_r - ll_l))
}
file.copy(file.path(d0, "init_par_truth.txt"), file.path(d0, "init_par.txt"), overwrite = TRUE)
out <- do.call(rbind, res); write.csv(out, sprintf("results/t9_bisect_%d_t%d%s.csv", N, T, if (length(args) >= 3) paste0("_", args[3]) else ""), row.names = FALSE); print(out, row.names = FALSE)
