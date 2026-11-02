# T4: compare factor scores. Legacy default (ML, no prior), legacy with prior, R (posterior mode), and the truth.
args <- commandArgs(trailingOnly = TRUE); N <- as.integer(args[1])
tr <- readRDS(sprintf("data/truth_%d.rds", N)); f <- tr$f
rd <- function(d) { x <- read.table(d, col.names = c("id", "f1", "se1", "f2", "se2", "f3", "se3")); x[order(x$id), ] }
L0 <- rd(sprintf("results/legacy_s5_%d/factor_predictions.txt", N))
f1 <- sprintf("results/legacy_s5p_%d/factor_predictions.txt", N); has_prior <- file.exists(f1) && length(readLines(f1)) == nrow(L0)
L1 <- if (has_prior) rd(f1) else NULL
R <- readRDS(sprintf("results/r_scores_%d.rds", N))
Rs <- R$obs_se$scores; stopifnot(all(Rs$id == L0$id))
cmp <- function(a, b, lab) cat(sprintf("%-38s max|diff| f: %.2e %.2e %.2e   se: %.2e %.2e %.2e\n", lab,
  max(abs(a$f1 - b$f1)), max(abs(a$f2 - b$f2)), max(abs(a$f3 - b$f3)), max(abs(a$se1 - b$se1)), max(abs(a$se2 - b$se2)), max(abs(a$se3 - b$se3))))
Rd <- data.frame(f1 = Rs$factor_1, f2 = Rs$factor_2, f3 = Rs$factor_3, se1 = Rs$se_factor_1, se2 = Rs$se_factor_2, se3 = Rs$se_factor_3)
Rp <- R$prior_se$scores; Rdp <- data.frame(f1 = Rp$factor_1, f2 = Rp$factor_2, f3 = Rp$factor_3, se1 = Rp$se_factor_1, se2 = Rp$se_factor_2, se3 = Rp$se_factor_3)
if (has_prior) { cmp(L1, Rdp, "legacy+prior vs R (prior SE)"); cmp(L1, Rd, "legacy+prior vs R (obs-only SE)") } else cat("legacy+prior scores: not available\n")
cmp(L0, Rd, "legacy default (ML) vs R (obs-only SE)")
cat(sprintf("legacy ML vs R: share of obs with |score diff| > 0.01: %.3f %.3f %.3f; mean |diff|: %.3f %.3f %.3f\n", mean(abs(L0$f1-Rd$f1)>.01), mean(abs(L0$f2-Rd$f2)>.01), mean(abs(L0$f3-Rd$f3)>.01), mean(abs(L0$f1-Rd$f1)), mean(abs(L0$f2-Rd$f2)), mean(abs(L0$f3-Rd$f3))))
rmse <- function(a) sprintf("%.4f %.4f %.4f", sqrt(mean((a$f1 - f[, 1])^2)), sqrt(mean((a$f2 - f[, 2])^2)), sqrt(mean((a$f3 - f[, 3])^2)))
cat("RMSE vs truth  legacy ML:", rmse(L0), " legacy+prior:", if (has_prior) rmse(L1) else "n/a", " R:", rmse(Rd), "\n")
if (has_prior) cat("legacy SE guard (0.8 shrink): share of obs with legacy+prior SE differing from R prior SE by > 1e-3:",
    sprintf("%.3f %.3f %.3f", mean(abs(L1$se1 - Rdp$se1) > 1e-3), mean(abs(L1$se2 - Rdp$se2) > 1e-3), mean(abs(L1$se3 - Rdp$se3) > 1e-3)), "\n")
cat(sprintf("SE comparison legacy ML vs R obs-only: mean |diff| %.4f %.4f %.4f; mean legacy SE %.3f %.3f %.3f; mean R SE %.3f %.3f %.3f\n", mean(abs(L0$se1-Rd$se1)), mean(abs(L0$se2-Rd$se2)), mean(abs(L0$se3-Rd$se3)), mean(L0$se1), mean(L0$se2), mean(L0$se3), mean(Rd$se1), mean(Rd$se2), mean(Rd$se3)))
