# T6 replication helper. Usage: Rscript t6_rep.R <N> <rep> <step>
#  step "fit1":   R stage 1 (first replication only) -> results/r_s1_<N>.rds (template for later replications)
#  step "scores": build an R stage-1 fit from legacy meas_par (estimates agree to 1e-8), R posterior-mode scores,
#                 and factor_predictions.txt from them -> results/t6_tmp/
#  step "adapt":  R adaptive education stage with scores "L" or "R" (4th arg) -> results/t6_tmp/r_adapt_<src>.csv
args <- commandArgs(trailingOnly = TRUE); N <- as.integer(args[1]); rep <- as.integer(args[2]); step <- args[3]
suppressMessages(pkgload::load_all("../../factorana_R", compile = FALSE, quiet = TRUE)); source("model_r.R"); source("legacy_names.R")
dat <- read.csv(sprintf("data/mini_%d.csv", N)); dir.create("results/t6_tmp", showWarnings = FALSE)
ctrl <- function(adapt = FALSE) define_estimation_control(n_quad_points = 8, num_cores = 8, adaptive_integration = adapt, adapt_int_thresh = 0.5)
if (step == "fit1") {
  fit <- estimate_model_rcpp(build_mini_system(dat, 1), dat, control = ctrl(), optimizer = "nlminb", parallel = TRUE, verbose = FALSE)
  stopifnot(fit$convergence == 0); saveRDS(fit, sprintf("results/t6_template_s1_%d.rds", N)); cat("template fit saved\n")
} else if (step == "scores") {
  fit <- readRDS(sprintf("results/t6_template_s1_%d.rds", N))
  tab <- legacy_param_table(sprintf("results/legacy_s1_%d/run.log", N)); par <- read.table(sprintf("results/legacy_s1_%d/meas_par.txt", N), col.names = c("idx","est","se"))
  leg <- merge(tab[, c("idx", "param")], par, by = "idx")   # keep one est/se pair (tab also has est/se columns)
  stopifnot(all(leg$param %in% names(fit$estimates)), length(leg$est) == nrow(leg))
  fit$estimates[leg$param] <- leg$est; fit$std_errors[leg$param] <- leg$se   # legacy stage-1 values into the R fit object
  saveRDS(fit, "results/t6_tmp/s1_fit.rds")
  fs <- estimate_factorscores_rcpp(fit, dat, control = ctrl(), parallel = TRUE, verbose = FALSE, include_prior = FALSE, id_var = "id")
  saveRDS(fs, "results/t6_tmp/r_scores.rds")
  write.table(data.frame(fs$id, fs$factor_1, fs$se_factor_1, fs$factor_2, fs$se_factor_2, fs$factor_3, fs$se_factor_3),
              "results/t6_tmp/factor_predictions_R.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
  cat("R scores done, converged", mean(fs$converged), "\n")
} else if (step == "adapt") {
  src <- args[4]; fit <- readRDS("results/t6_tmp/s1_fit.rds")
  if (src == "L") { fp <- read.table(sprintf("results/legacy_s5_%d/factor_predictions.txt", N), col.names = c("id","f1","se1","f2","se2","f3","se3"))
    fp <- fp[match(dat$id, fp$id), ]; sc <- as.matrix(fp[, c("f1","f2","f3")]); se <- as.matrix(fp[, c("se1","se2","se3")])
  } else { fs <- readRDS("results/t6_tmp/r_scores.rds"); sc <- as.matrix(fs[, c("factor_1","factor_2","factor_3")]); se <- as.matrix(fs[, c("se_factor_1","se_factor_2","se_factor_3")]) }
  f3 <- estimate_model_rcpp(build_mini_system(dat, 3, fit), dat, control = ctrl(TRUE), optimizer = "nlminb", parallel = TRUE, verbose = FALSE,
                            factor_scores = sc, factor_ses = se, factor_vars = c(factor_var_1 = 1, factor_var_2 = 1, factor_var_3 = 1))
  write.csv(data.frame(param = names(f3$estimates), est = unname(f3$estimates), se = unname(f3$std_errors), conv = f3$convergence),
            sprintf("results/t6_tmp/r_adapt_%s.csv", src), row.names = FALSE)
  cat(sprintf("R adaptive (%s) convergence=%d\n", src, f3$convergence))
}
