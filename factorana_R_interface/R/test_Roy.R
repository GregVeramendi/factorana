#Load package / source files
devtools::load_all(".")  # if in a package
source("R/define_model_component.R")
source("R/estimate_model.R")
source("R/define_model_system.R")
source ("R/define_factor_model.R")
source("R/define_estimation_control.R")

set.seed(123)

#Generate Fake Data
n <- 5000

f <- rnorm(n, 0, 1)
U1 <- rnorm(n, 0 ,1)
U0 <- rnorm(n, 0 ,1)
UV <- rnorm(n, 0 ,1)
UT <- rnorm(n, 0 ,1)

Z0 <- X0 <- Q0 <- 1
Z1 <- rnorm(n, 0, 2); X1 <- rnorm(n, 0, 2); Q1 <- rnorm(n, 0, 1)

Y1 <- 2*X0 + X1 + 2*f + U1
Y0 <- 1*X0 + X1 + 1*f + U0
I  <- 0.5*Z0 + Z1 + f + UV
T  <- Q0 + Q1 + f + UT
D  <- as.integer(I > 0)

Y  <- ifelse(D == 1L, Y1, Y0)

dat <- data.frame(Y, D, Z1, X1, Q1, T)
dat$Z0 <- 1; dat$X0 <- 1; dat$Q0 <- 1
dat$eval_y0 <- 1L - dat$D


# Fake factor model (simplified placeholder)

factor_1 <- define_factor_model(1, 1, 8)

# Components
mc_sel <- define_model_component("selection", dat, "D", factor_1, covariates=c("Z0","Z1"), model_type="probit")
mc_y1  <- define_model_component("Y1", dat, "Y", factor_1, evaluation_indicator="D",       covariates=c("X0","X1"), model_type="linear")
mc_y0  <- define_model_component("Y0", dat, "Y", factor_1, evaluation_indicator="eval_y0", covariates=c("X0","X1"), model_type="linear")
mc_T   <- define_model_component("T",  dat, "T", factor_1, covariates=c("Q0","Q1"), model_type="linear")

ms <- define_model_system(factor=factor_1, components=list(mc_sel, mc_y1, mc_y0, mc_T))
ctrl <- define_estimation_control(num_cores = 1)

# Initialize each component
inits <- lapply(ms$components, initialize_parameters)

for (ini in inits) {
  expect_true(all(c("intercept","betas","loading","factor_var","factor_cor") %in% names(ini)))
  expect_false(any(is.na(c(ini$intercept, ini$loading, ini$factor_var, ini$factor_cor))))
}

# Estimator smoke test (+ optional CSV)
tmp <- tempdir(); owd <- setwd(tmp); on.exit(setwd(owd), add=TRUE)
if (!dir.exists("results")) dir.create("results")
out <- estimate_model(ms, ctrl)
expect_s3_class(out, "data.frame")
# If your estimator writes a CSV, check it exists:
csvs <- list.files("results", pattern="system_inits_.*\\.csv$", full.names=TRUE)
if (length(csvs) > 0) {
  expect_true(file.exists(tail(csvs, 1)))
}
