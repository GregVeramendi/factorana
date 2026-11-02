# T8 summary: cost (integration points per person, seconds) against accuracy for each adaptive estimator, over the
# T6/T8 replications. Accuracy in units of the full-quadrature SE, as in t6_summary.R. Usage: Rscript t8_summary.R [N]
args <- commandArgs(trailingOnly = TRUE); N <- if (length(args) >= 1) as.integer(args[1]) else 2001L
source("legacy_names.R"); source("truth_params.R")
reps <- sort(as.integer(sub("rep", "", list.files("results/t8", "^rep[0-9]+$"))))
read_leg <- function(dir, col = "est") { tab <- legacy_param_table(file.path(dir, "run.log"))
  par <- read.table(file.path(dir, "fullmodel_par.txt"), col.names = c("idx", "est", "se"))
  m <- merge(tab[!tab$fixed, c("idx", "param")], par, by = "idx"); setNames(m[[col]], m$param) }
ipopt_secs <- function(dir) { x <- grep("Total seconds in IPOPT", readLines(file.path(dir, "run.log")), value = TRUE)
  as.numeric(sub(".*=\\s*", "", x[1])) }                                     # first block = education group
read_r <- function(f) { d <- read.csv(f); list(est = setNames(d$est, d$param), conv = d$conv[1], secs = d$opt_secs[1], wall = d$wall_secs[1]) }
# Points per person from the factor-score SEs (prior SD = 1). rule: "factor" (R: SE > SD sends that factor to mx points
# at the prior), "person" (any SE > SD sends the person to full quadrature), "legacy" (mean SE/SD > 1.5 does the same).
points <- function(se, rule, thr, mx) { nq <- pmin(1 + 2 * floor(se / thr), mx)
  if (rule == "factor") nq[se > 1] <- mx
  if (rule == "person") nq[apply(se > 1, 1, any), ] <- mx
  if (rule == "legacy") nq[rowMeans(se) > 1.5, ] <- mx
  mean(apply(nq, 1, prod)) }
se_of <- function(fs) as.matrix(fs[, c("se_factor_1", "se_factor_2", "se_factor_3")])
variants <- c("obs_t0.5_m8", "obsP_t0.5_m8", "post_t0.5_m8", "post_t0.3_m8", "obs_t1.0_m8", "post_t1.0_m8",
              if (length(args) >= 2) args[-1])      # extra variants, e.g. post_t0.7_m8 post_t0.85_m8
vspec <- function(v) { m <- regmatches(v, regexec("^(obs|obsP|post)_t([0-9.]+)_m([0-9]+)$", v))[[1]]
  list(se = m[2], thr = as.numeric(m[3]), mx = as.integer(m[4]), rule = if (m[2] == "obsP") "person" else "factor") }
rows <- list(); cost <- list()
for (r in reps) {
  d6 <- sprintf("results/t6/rep%d", r); d8 <- sprintf("results/t8/rep%d", r)
  if (!all(file.exists(file.path(d8, sprintf("r_%s.csv", variants))))) next
  tr <- readRDS(file.path(d6, sprintf("truth_%d.rds", N)))$truth
  full <- read_leg(file.path(d6, sprintf("legacy_s2_%d", N))); se_full <- read_leg(file.path(d6, sprintf("legacy_s2_%d", N)), "se")
  ed <- names(full)[!grepl("^(g_|cog|misslead|lead|stresstol|wage|factor_var)", names(full))]
  truth <- truth_r_params(tr, ed, subset = TRUE)
  fp <- read.table(file.path(d6, sprintf("legacy_s5_%d", N), "factor_predictions.txt"), col.names = c("id", "f1", "se1", "f2", "se2", "f3", "se3"))
  se_ml <- as.matrix(fp[, c("se1", "se2", "se3")]); se_obs <- se_of(readRDS(file.path(d6, "r_scores.rds"))); se_post <- se_of(readRDS(file.path(d8, "r_scores_post.rds")))
  est <- list(); cst <- list()
  add <- function(name, e, pts, secs, code, conv = 0) { est[[name]] <<- e[ed]; cst[[name]] <<- data.frame(estimator = name, code = code, points = pts, secs = secs, conv = conv) }
  add("full quadrature (legacy)", full, 512, ipopt_secs(file.path(d6, sprintf("legacy_s2_%d", N))), "legacy")
  add("ML scores, legacy rule [production]", read_leg(file.path(d6, sprintf("legacy_s2aL_%d", N))), points(se_ml, "legacy", 0.5, 8),
      ipopt_secs(file.path(d6, sprintf("legacy_s2aL_%d", N))), "legacy")
  add("posterior mode, obs SE, legacy rule (no fallback fires)", read_leg(file.path(d6, sprintf("legacy_s2aR_%d", N))), points(se_obs, "legacy", 0.5, 8),
      ipopt_secs(file.path(d6, sprintf("legacy_s2aR_%d", N))), "legacy")
  for (v in variants) { s <- vspec(v); x <- read_r(file.path(d8, sprintf("r_%s.csv", v)))
    se <- switch(s$se, obs = se_obs, obsP = se_obs, post = se_post)
    lab <- sprintf("posterior mode, %s SE, %s, t=%.1f, max=%d", if (s$se == "post") "posterior" else "obs", c(factor = "per-factor fallback", person = "per-person fallback")[s$rule], s$thr, s$mx)
    add(lab, x$est, points(se, s$rule, s$thr, s$mx), x$secs, "R", x$conv != 0) }
  if (file.exists(file.path(d8, "r_full.csv"))) { x <- read_r(file.path(d8, "r_full.csv")); cst[["full quadrature (R)"]] <- data.frame(estimator = "full quadrature (R)", code = "R", points = 512, secs = x$secs, conv = x$conv != 0) }
  for (nm in names(est)) rows[[length(rows) + 1]] <- data.frame(rep = r, estimator = nm, param = ed, z = (est[[nm]] - truth) / se_full[ed], dev_full = (est[[nm]] - full[ed]) / se_full[ed])
  cost[[length(cost) + 1]] <- cbind(rep = r, do.call(rbind, cst))
}
d <- do.call(rbind, rows); cs <- do.call(rbind, cost)
cat("replications used:", length(unique(d$rep)), " education parameters:", length(unique(d$param)), "\n")
acc <- do.call(rbind, lapply(split(d, d$estimator), function(x) { pm <- tapply(x$z, x$param, mean)
  data.frame(estimator = x$estimator[1], mean_abs_bias_SE = round(mean(abs(pm)), 3), RMSE_SE = round(sqrt(mean(x$z^2)), 3),
             mean_dev_full_SE = round(mean(abs(x$dev_full)), 3), max_dev_full_SE = round(max(abs(x$dev_full)), 2)) }))
cc <- do.call(rbind, lapply(split(cs, cs$estimator), function(x) data.frame(estimator = x$estimator[1], code = x$code[1], points_per_person = round(mean(x$points), 1),
             secs = round(mean(x$secs, na.rm = TRUE), 1), nonconv = sum(x$conv))))
full_secs <- c(legacy = cc$secs[cc$estimator == "full quadrature (legacy)"], R = cc$secs[cc$estimator == "full quadrature (R)"])
cc$time_vs_full_same_code <- round(cc$secs / full_secs[cc$code], 3)
res <- merge(cc, acc, by = "estimator", all.x = TRUE)
ord <- c("full quadrature (legacy)", "full quadrature (R)", "ML scores, legacy rule [production]", "posterior mode, obs SE, legacy rule (no fallback fires)",
         grep("^posterior mode, obs SE, per-factor", res$estimator, value = TRUE), grep("per-person", res$estimator, value = TRUE), grep("posterior SE", res$estimator, value = TRUE))
res <- res[order(match(res$estimator, ord)), ]; print(res, row.names = FALSE, right = FALSE)
cat("\nPer-parameter bias (mean z) above 0.5 SE for any adaptive estimator:\n")
pm <- tapply(d$z, list(d$param, d$estimator), mean); flag <- apply(abs(pm[, colnames(pm) != "full quadrature (legacy)", drop = FALSE]), 1, max) > 0.5
if (any(flag)) print(round(pm[flag, , drop = FALSE], 2)) else cat("(none)\n")
write.csv(d, sprintf("results/t8_long_%d.csv", N), row.names = FALSE); write.csv(res, sprintf("results/t8_summary_%d.csv", N), row.names = FALSE)
