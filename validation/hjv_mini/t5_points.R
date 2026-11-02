# T5: points per observation implied by each code's adaptive rule, from the same factor scores (threshold 0.5, max 8).
# Usage: Rscript t5_points.R <N> <src L|R>
args <- commandArgs(trailingOnly = TRUE); N <- as.integer(args[1]); src <- args[2]
if (src == "L") { fp <- read.table(sprintf("results/legacy_s5_%d/factor_predictions.txt", N), col.names = c("id","f1","se1","f2","se2","f3","se3"))
  se <- as.matrix(fp[, c("se1","se2","se3")]) } else { fs <- readRDS(sprintf("results/r_scores_%d.rds", N))$obs_se$scores
  se <- as.matrix(fs[, c("se_factor_1","se_factor_2","se_factor_3")]) }
sd <- 1; thr <- 0.5; mx <- 8
nq <- pmin(1 + 2 * floor(se / sd / thr), mx)                     # shared formula
nq_r <- nq; nq_r[se > sd] <- mx                                   # R: per-factor fallback when SE > SD
nq_l <- nq; nq_l[nq_l > mx] <- mx; fb <- rowMeans(se / sd) > 1.5; nq_l[fb, ] <- mx   # legacy: observation-level fallback
tot_r <- apply(nq_r, 1, prod); tot_l <- apply(nq_l, 1, prod)
cat(sprintf("scores %s: obs where the two rules give different point counts: %d of %d (%.1f%%); legacy fallback fires for %d obs, R per-factor fallback for %d factor-obs\n",
            src, sum(tot_r != tot_l), nrow(se), 100 * mean(tot_r != tot_l), sum(fb), sum(se > sd)))
cat(sprintf("mean total points per obs: R rule %.1f, legacy rule %.1f (full quadrature = 512)\n", mean(tot_r), mean(tot_l)))
print(table(nq_r[,1]) ); print(table(nq_l[,1]))
