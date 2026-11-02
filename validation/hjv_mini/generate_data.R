# Mini-HJV data generating process. Writes data/mini_<N>.{csv,raw}, data/mini_varlist.txt, data/truth_<N>.rds.
# Usage: Rscript generate_data.R <N> <seed> [ntypes]
# ntypes > 1 adds latent types (T9): type probabilities follow a logit in the factors, and every education and
# labor-market component gets type-specific intercepts (per alternative for the logits). Files then carry the suffix
# _t<ntypes>. With ntypes = 1 nothing is drawn for types, so the untyped data are unchanged for a given seed.
args <- commandArgs(trailingOnly = TRUE)
N <- if (length(args) >= 1) as.integer(args[1]) else 2000L
seed <- if (length(args) >= 2) as.integer(args[2]) else 1L
NT <- if (length(args) >= 3) as.integer(args[3]) else 1L
tag <- if (NT > 1) sprintf("_t%d", NT) else ""
set.seed(seed)
MISS <- -9999

# ---- factors and covariates --------------------------------------------------
f <- matrix(rnorm(3 * N), N, 3, dimnames = list(NULL, c("f1", "f2", "f3")))
x1 <- rnorm(N); x2 <- rbinom(N, 1, 0.4)
X <- cbind(intercept = 1, x1 = x1, x2 = x2)

lin <- function(b, lam, sig) as.numeric(X %*% b + f %*% lam + rnorm(N, 0, sig))
pro <- function(b, lam, extra = 0) as.integer(X %*% b + f %*% lam + extra + rnorm(N) > 0)
# multinomial logit draw, K alternatives, alternative 1 is the reference (utility 0)
mnl <- function(B, L, extra = NULL) {        # B: (K-1) x 3, L: (K-1) x 3, extra: N x (K-1)
  U <- X %*% t(B) + f %*% t(L); if (!is.null(extra)) U <- U + extra
  U <- cbind(0, U); U <- U - log(-log(matrix(runif(length(U)), nrow(U))))
  max.col(U, ties.method = "first")
}
truth <- list()

# ---- latent types (T9 only) ----------------------------------------------------
# P(type = t | f) = exp(a_t + L_t f) / (1 + sum_s exp(a_s + L_s f)), type 1 the reference. Type effects are
# deterministic patterns: dvec(s) for a single-index component, dmat(s) per alternative for a logit.
if (NT > 1) {
  tt <- seq_len(NT - 1)
  truth$types <- list(a = seq(0.3, -0.3, length.out = NT - 1),
                      L = t(sapply(tt, function(t) c(0.5, -0.3, 0.2) * (1 - 0.3 * (t - 1)) * (-1)^(t - 1))))
  eta <- cbind(0, matrix(truth$types$a, N, NT - 1, byrow = TRUE) + f %*% t(truth$types$L))
  type <- max.col(eta - log(-log(matrix(runif(length(eta)), nrow(eta)))), ties.method = "first")
} else type <- rep(1L, N)
tdum <- if (NT > 1) outer(type, 2:NT, "==") * 1 else matrix(0, N, 0)
dvec <- function(s) s * seq(-1, 1, length.out = max(NT - 1, 1))[seq_len(NT - 1)]          # length NT-1
dmat <- function(s) outer(s, seq(-1, 1, length.out = max(NT - 1, 1))[seq_len(NT - 1)])     # (K-1) x (NT-1)
dT  <- function(d) if (NT > 1) as.numeric(tdum %*% d) else 0
dTm <- function(D) if (NT > 1) tdum %*% t(D) else 0

# ---- Stage 1: measurement system -------------------------------------------------
# grades load on f1, f2, f3; sports mostly on f2
truth$g_eng9    <- list(b = c(0.3, 0.2, -0.1), lam = c(0.8, 0.3, 0.4), sig = 0.7)
truth$g_math9   <- list(b = c(0.1, 0.3,  0.0), lam = c(0.9, 0.1, 0.5), sig = 0.6)
truth$g_swe9    <- list(b = c(0.2, 0.1, -0.2), lam = c(0.7, 0.3, 0.4), sig = 0.7)
truth$g_sports9 <- list(b = c(0.0, 0.1,  0.2), lam = c(0.1, 0.8, 0.2), sig = 0.7)
truth$g_gpa9    <- list(b = c(0.0, 0.2,  0.0), lam = c(0.6, 0.3, 0.6), sig = 0.3)
for (m in c("g_eng9", "g_math9", "g_swe9", "g_sports9", "g_gpa9")) assign(m, with(truth[[m]], lin(b, lam, sig)))
# cognitive (f1 only), observed for 70 percent
truth$cog1 <- list(b = c(0.0, 0.1, 0.0), lam = c(1.0, 0, 0), sig = 0.5)
truth$cog2 <- list(b = c(0.1, 0.0, 0.1), lam = c(0.9, 0, 0), sig = 0.6)
truth$cog3 <- list(b = c(-0.1, 0.1, 0.0), lam = c(1.1, 0, 0), sig = 0.7)
ind_cog <- rbinom(N, 1, 0.7)
for (m in c("cog1", "cog2", "cog3")) assign(m, ifelse(ind_cog == 1, with(truth[[m]], lin(b, lam, sig)), MISS))
# leadership and stress tolerance (f1, f2), missing when misslead = 1 (probit on f1)
truth$misslead  <- list(b = c(-0.3, 0.2, -0.1), lam = c(-0.5, 0, 0))
misslead <- with(truth$misslead, pro(b, lam)); ind_lead <- 1L - misslead
truth$lead      <- list(b = c(0.0, 0.1, 0.0), lam = c(0.3, 1.0, 0), sig = 0.6)
truth$stresstol <- list(b = c(0.1, 0.0, 0.1), lam = c(0.2, 0.9, 0), sig = 0.7)
lead      <- ifelse(ind_lead == 1, with(truth$lead, lin(b, lam, sig)), MISS)
stresstol <- ifelse(ind_lead == 1, with(truth$stresstol, lin(b, lam, sig)), MISS)

# ---- Stage 3: education (needed before the high-school grades) -------------------
truth$advmath <- list(b = c(-0.2, 0.3, -0.1), lam = c(0.9, 0.2, 0.3), d = dvec(0.5))
truth$adveng  <- list(b = c(-0.1, 0.2,  0.0), lam = c(0.8, 0.3, 0.2), d = dvec(-0.4))
advmath <- with(truth$advmath, pro(b, lam, dT(d))); adveng <- with(truth$adveng, pro(b, lam, dT(d)))
# HS track: 4 choices; covariates intercept, x1, x2, advmath, adveng. x2 coefficient of choice 4 fixed at 0.
truth$HS_track <- list(B = rbind(c(-0.5, 0.2, 0.1), c(-1.0, 0.3, 0.2), c(-1.5, 0.4, 0.0)),
                       g_adv = rbind(c(0.4, 0.3), c(0.8, 0.5), c(1.2, 0.7)),
                       L = rbind(c(0.5, 0.2, 0.3), c(0.9, 0.3, 0.4), c(1.3, 0.4, 0.5)), D = dmat(c(0.3, -0.5, 0.4)))
HS_track <- with(truth$HS_track, mnl(B, L, cbind(advmath, adveng) %*% t(g_adv) + dTm(D)))
HS_track2 <- as.integer(HS_track == 2); HS_track3 <- as.integer(HS_track == 3); HS_track4 <- as.integer(HS_track == 4)
ind_hs <- rbinom(N, 1, 0.85)                       # in high-school sample
HS <- cbind(HS_track2, HS_track3, HS_track4)
# high-school grades, observed if ind_hs
truth$g_swe10    <- list(b = c(0.2, 0.1, -0.1), g_hs = c(0.3, 0.5, 0.8), lam = c(0.7, 0.3, 0.5), sig = 0.7)
truth$g_sports10 <- list(b = c(0.0, 0.1,  0.2), g_hs = c(0.1, 0.2, 0.3), lam = c(0.1, 0.8, 0.3), sig = 0.7)
truth$g_gpahs    <- list(b = c(150, 8, -5),     g_hs = c(10, 20, 30),    lam = c(20, 8, 25),     sig = 30)
for (m in c("g_swe10", "g_sports10", "g_gpahs")) assign(m, ifelse(ind_hs == 1,
  with(truth[[m]], lin(b, lam, sig) + as.numeric(HS %*% g_hs)), MISS))
# apply to college, then ranked application (5 majors, up to 5 ranks, same nest may repeat)
truth$apply <- list(b = c(-0.3, 0.3, 0.1), g = c(0.5, 0.4, 0.3, 0.6), lam = c(0.6, 0.4, 0.3), d = dvec(0.6))
Z <- cbind(advmath, adveng, HS_track3, HS_track4)
apply <- with(truth$apply, pro(b, lam, as.numeric(Z %*% g) + dT(d)))
truth$app <- list(B = rbind(c(-0.2, 0.3, 0.1), c(-0.4, 0.1, 0.2), c(-0.6, -0.2, 0.0), c(-0.8, 0.2, -0.1)),
                  G = rbind(c(0.3, 0.2, 0.1, 0.2), c(0.5, 0.1, 0.2, 0.3), c(0.2, 0.4, 0.3, 0.1), c(0.6, 0.3, 0.0, 0.4)),
                  L = rbind(c(0.4, 0.1, 0.2), c(0.6, 0.2, 0.1), c(0.3, 0.5, 0.2), c(0.7, 0.3, 0.3)), D = dmat(c(0.4, -0.3, 0.2, -0.5)))
nranks <- sample(1:5, N, replace = TRUE)
app <- sapply(1:5, function(r) with(truth$app, mnl(B, L, Z %*% t(G) + dTm(D))))
app[apply == 0, ] <- MISS
for (r in 1:5) app[apply == 1 & nranks < r, r] <- 0   # unused ranks are skipped
colnames(app) <- paste0("app", 1:5)
# SweSAT: take it, then score
truth$Dswesat <- list(b = c(-0.2, 0.2, 0.0), g = c(0.3, 0.2, 0.4, 0.5), lam = c(0.5, 0.2, 0.3), d = dvec(-0.5))
Dswesat <- with(truth$Dswesat, pro(b, lam, as.numeric(Z %*% g) + dT(d)))
truth$swesat  <- list(b = c(0.1, 0.3, -0.1), g = c(0.4, 0.3, 0.2, 0.3), lam = c(0.9, 0.2, 0.3), sig = 0.8, d = dvec(0.5))
swesat <- ifelse(Dswesat == 1, with(truth$swesat, lin(b, lam, sig) + as.numeric(Z %*% g) + dT(d)), MISS)
# enroll (if applied), switch major (if enrolled), graduation probits (if enrolled)
truth$enroll <- list(b = c(0.4, 0.1, 0.0), g = c(0.2, 0.1, 0.1, 0.2), lam = c(0.3, 0.2, 0.1), d = dvec(0.4))
enroll <- ifelse(apply == 1, with(truth$enroll, pro(b, lam, as.numeric(Z %*% g) + dT(d))), MISS)
ind_enroll <- as.integer(enroll == 1)
truth$switch <- list(B = rbind(c(-1.0, 0.1, 0.0), c(-1.2, 0.2, 0.1), c(-1.4, 0.0, 0.2)),
                     G = rbind(c(0.2, 0.1, 0.1, 0.2), c(0.1, 0.3, 0.2, 0.1), c(0.3, 0.2, 0.0, 0.3)),
                     L = rbind(c(0.3, 0.2, 0.1), c(0.4, 0.1, 0.2), c(0.2, 0.3, 0.3)), D = dmat(c(-0.3, 0.4, 0.2)))
switch_major <- ifelse(ind_enroll == 1, with(truth$switch, mnl(B, L, Z %*% t(G) + dTm(D))), MISS)
grad <- matrix(MISS, N, 3, dimnames = list(NULL, paste0("grad", 1:3)))
truth$grad <- list(b = rbind(c(0.2, 0.1, 0.0), c(0.0, 0.2, 0.1), c(-0.2, 0.1, 0.2)),
                   g = rbind(c(0.3, 0.2, 0.1, 0.2), c(0.2, 0.3, 0.2, 0.1), c(0.1, 0.1, 0.3, 0.3)),
                   lam = rbind(c(0.5, 0.3, 0.2), c(0.4, 0.4, 0.1), c(0.6, 0.2, 0.3)), d = rbind(dvec(0.5), dvec(-0.3), dvec(0.4)))
for (j in 1:3) grad[ind_enroll == 1, j] <- with(truth$grad, pro(b[j, ], lam[j, ], as.numeric(Z %*% g[j, ]) + dT(d[j, ])))[ind_enroll == 1]

# ---- Stage 4: labor market -------------------------------------------------------
truth$wage1 <- list(b = c(2.0, 0.1, -0.1), g = c(0.1, 0.2), lam = c(0.3, 0.2, 0.1), sig = 0.5, d = dvec(0.3))
truth$wage2 <- list(b = c(2.2, 0.1,  0.0), g = c(0.2, 0.3), lam = c(0.4, 0.1, 0.2), sig = 0.6, d = dvec(-0.4))
W <- cbind(HS_track3, HS_track4)
wage1 <- with(truth$wage1, lin(b, lam, sig) + as.numeric(W %*% g) + dT(d))
wage2 <- with(truth$wage2, lin(b, lam, sig) + as.numeric(W %*% g) + dT(d))

dat <- data.frame(id = seq_len(N), intercept = 1, x1, x2,
                  g_eng9, g_math9, g_swe9, g_sports9, g_gpa9, cog1, cog2, cog3, ind_cog,
                  misslead, ind_lead, lead, stresstol,
                  advmath, adveng, HS_track, HS_track2, HS_track3, HS_track4, ind_hs,
                  g_swe10, g_sports10, g_gpahs, apply, app, Dswesat, swesat, enroll, ind_enroll,
                  switch_major, grad, wage1, wage2)
if (NT > 1) dat$true_type <- type
dir.create("data", showWarnings = FALSE)
write.csv(dat, sprintf("data/mini_%d%s.csv", N, tag), row.names = FALSE)
write.table(dat, sprintf("data/mini_%d%s.raw", N, tag), row.names = FALSE, col.names = FALSE)
writeLines(c(ncol(dat), names(dat)), sprintf("data/mini_varlist%s.txt", tag))
saveRDS(list(truth = truth, f = f, N = N, seed = seed, ntypes = NT, type = type), sprintf("data/truth_%d%s.rds", N, tag))
if (NT > 1) cat("types:", round(prop.table(table(type)), 2), "\n")
cat(sprintf("N=%d  cols=%d  shares: cog %.2f lead %.2f hs %.2f apply %.2f enroll %.2f\n", N, ncol(dat),
            mean(ind_cog), mean(ind_lead), mean(ind_hs), mean(apply), mean(ind_enroll)))
cat("HS_track:", round(prop.table(table(HS_track)), 2), " app1:", round(prop.table(table(app[apply == 1, 1])), 2), "\n")
