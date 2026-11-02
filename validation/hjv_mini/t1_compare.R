# T1: compare log-likelihood, gradient and Hessian at the true parameters. Usage: Rscript t1_compare.R <N> [tag] [outcome]
# With "outcome" the gradient and Hessian comparison is restricted to the outcome-model parameters: legacy stage 0
# returns zero derivatives for the measurement models when outcome models are present, and those parameters are fixed
# in the final stage anyway. The log-likelihood comparison always covers the whole system.
args <- commandArgs(trailingOnly = TRUE); N <- as.integer(args[1]); tag <- if (length(args) >= 2) args[2] else ""
outcome_only <- length(args) >= 3 && args[3] == "outcome"
source("legacy_names.R")
r <- readRDS(sprintf("results/t1_r%s_%d.rds", tag, N))
ldir <- sprintf("results/legacy_s0%s_%d", tag, N)
dump <- scan(file.path(ldir, "t1_dump.txt"), quiet = TRUE)
free_l <- scan(file.path(ldir, "t1_freeparlist.txt"), quiet = TRUE)
map <- read.csv(file.path(ldir, "legacy_index_map.csv"))
nf <- length(free_l); stopifnot(length(dump) == 1 + nf + nf * (nf + 1) / 2)
obj_l <- dump[1]; grad_l <- dump[2:(1 + nf)]; hess_l <- dump[-(1:(1 + nf))]
names_l <- map$param[match(free_l, map$idx)]
# R: free parameters are those with a finite non-NA gradient slot; names in R order
free_r <- r$free_names
stopifnot(length(free_r) == length(r$gradient), setequal(free_r, names_l))
unpack <- function(h, n) { H <- matrix(0, n, n); k <- 1; for (i in 1:n) for (j in i:n) { H[i, j] <- H[j, i] <- h[k]; k <- k + 1 }; H }
H_r <- unpack(r$hessian, nf); dimnames(H_r) <- list(free_r, free_r)
H_l <- unpack(hess_l, nf); dimnames(H_l) <- list(names_l, names_l)
g_r <- setNames(r$gradient, free_r); g_l <- setNames(grad_l, names_l)[free_r]; H_l <- H_l[free_r, free_r]
if (outcome_only) {
  keep <- !grepl("^(g_|cog|misslead|lead|stresstol|factor_var)", free_r)
  cat(sprintf("restricting derivatives to %d outcome-model parameters (of %d free)\n", sum(keep), length(free_r)))
  g_r <- g_r[keep]; g_l <- g_l[keep]; H_r <- H_r[keep, keep]; H_l <- H_l[keep, keep]
}
const <- N * 3 * 0.5 * log(pi)
d_ll <- (-r$loglik) - obj_l - const
rel <- function(a, b) abs(a - b) / pmax(abs(b), 1e-8)
g_rel <- rel(g_r, -g_l); h_rel <- rel(H_r, -H_l)
cat(sprintf("loglik: R = %.10f  legacy = %.10f  gap after sqrt(pi) constant = %.3e (rel %.1e)\n", r$loglik, -obj_l - const, d_ll, abs(d_ll / r$loglik)))
cat(sprintf("gradient: max |g_r + g_l| = %.3e, max relative gap = %.3e at %s; max |g| = %.3e\n",
            max(abs(g_r + g_l)), max(g_rel), names(g_r)[which.max(g_rel)], max(abs(g_r))))
big <- abs(H_r) > 1e-6 * max(abs(H_r))      # relative gaps on elements that are not structurally ~0
cat(sprintf("hessian: %d elements; max |H_r + H_l| = %.3e; max relative gap on non-tiny elements = %.3e\n",
            length(H_r), max(abs(H_r + H_l)), max(h_rel[big])))
pass <- abs(d_ll / r$loglik) <= 1e-8 && max(g_rel) <= 1e-6 && max(h_rel[big]) <= 1e-6
cat("T1 pass (loglik rel 1e-8, gradient and Hessian elements rel 1e-6):", pass, "\n")
saveRDS(list(d_ll = d_ll, g_rel = g_rel, h_rel = h_rel, big = big), sprintf("results/t1_compare%s_%d.rds", tag, N))
