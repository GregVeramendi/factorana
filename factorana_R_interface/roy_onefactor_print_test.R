#!/usr/bin/env Rscript

cat("=== Roy one-factor demo ===\n")

# ---- Source your package code ----
r_files <- c("define_model_component.R",
             "define_model_system.R",
             "define_estimation_control.R",
             "define_factor_model.R",
             "estimate_model.R",
             "initialize_parameters.R")
for (f in r_files) {
  p <- file.path("R", f)
  if (file.exists(p)) {
    cat("Sourcing:", p, "\n")
    source(p, chdir = TRUE)
  }
}

# ---- Generate Fake Data ----
set.seed(123)
n  <- 5000
f  <- rnorm(n, 0, 1)
U1 <- rnorm(n, 0, 1)
U0 <- rnorm(n, 0, 1)
UV <- rnorm(n, 0, 1)
UT <- rnorm(n, 0, 1)

Z0 <- X0 <- Q0 <- 1
Z1 <- rnorm(n, 0, 2)
X1 <- rnorm(n, 0, 2)
Q1 <- rnorm(n, 0, 1)

Y1 <- 2*Z0 + X1 + 2*f + U1
Y0 <- 1*Z0 + X1 + 1*f + U0
I  <- 0.5*Z0 + Z1 + f + UV
D  <- as.integer(I > 0)
Tscore <- Q0 + Q1 + f + UT
Y  <- ifelse(D == 1L, Y1, Y0)

dat <- data.frame(Y, D, Z1, X1, Q1, T = Tscore)
dat$Z0 <- 1; dat$X0 <- 1; dat$Q0 <- 1
dat$eval_y0 <- 1L - dat$D

cat("\n--- Simulated data ---\n")
print(head(dat))

# ---- Factor model ----
factor_1 <- define_factor_model(1, 1, 8)
cat("\n--- Factor model ---\n")
print(factor_1)

# ---- Components ----
mc_sel <- define_model_component("selection", dat, "D", factor_1,
                                 covariates = c("Z1"), model_type = "probit")
mc_y1  <- define_model_component("Y1", dat, "Y", factor_1,
                                 evaluation_indicator = "D",
                                 covariates = c("X1"), model_type = "linear")
mc_y0  <- define_model_component("Y0", dat, "Y", factor_1,
                                 evaluation_indicator = "eval_y0",
                                 covariates = c("X1"), model_type = "linear")
mc_T   <- define_model_component("T",  dat, "T", factor_1,
                                 covariates = c("Q1"), model_type = "linear")

cat("\n--- Components ---\n")
print(mc_sel)
print(mc_y1)
print(mc_y0)
print(mc_T)

# ---- System & control ----
ms   <- define_model_system(factor = factor_1, components = list(mc_sel, mc_y1, mc_y0, mc_T))
ctrl <- define_estimation_control(num_cores = 1)

cat("\n--- Model system ---\n")
print(ms)
cat("\n--- Estimation control ---\n")
print(ctrl)

# ---- Initialize parameters ----
cat("\n--- Initialization results ---\n")
inits <- lapply(ms$components, initialize_parameters)
print(inits)

# ---- Estimate model ----
cat("\n--- Estimate model ---\n")
out <- estimate_model(ms, ctrl)
print(out)

# ---- LONG FORMAT ----
init_long <- do.call(rbind, lapply(seq_along(inits), function(i) {
  comp <- ms$components[[i]]
  ini  <- inits[[i]]

  rows <- list(
    data.frame(component=comp$name, param="intercept", value=unname(ini$intercept)),
    if (length(ini$betas)) data.frame(
      component=comp$name,
      param=paste0("beta[", names(ini$betas), "]"),
      value=unname(ini$betas)
    ),
    data.frame(component=comp$name, param="loading[f1]", value=unname(ini$loading))
  )
  do.call(rbind, rows)
}))

#prepend factor variance row first
fv_row <- data.frame(component = "factor",
                     param = "factor_var[f1]",
                     value = 1.0)
init_long <- rbind(fv_row, init_long)
rownames(init_long) <- NULL

# # add factor variance
# init_long <- rbind(init_long,
#                    data.frame(component="factor", param="factor_var[f1]", value=1.0)
# )
# rownames(init_long) <- NULL

write.csv(init_long, "results/system_inits_long.csv", row.names = FALSE)
cat("Wrote results/system_inits_long.csv\n")



# ---- Sanity stats ----
cat("\n--- Sanity stats ---\n")
m1 <- mean(dat$Y[dat$D == 1])
m0 <- mean(dat$Y[dat$D == 0])
sel <- mean(dat$D)
cat(sprintf("Treated mean: %.3f\n", m1))
cat(sprintf("Control mean: %.3f\n", m0))
cat(sprintf("Selection rate (D=1): %.3f\n", sel))

cat("\nDone\n")
