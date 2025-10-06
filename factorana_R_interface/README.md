
# factorana

<!-- badges: start -->
<!-- badges: end -->

usethis::use_readme_md()

# factorana (R interface)

R front-end for specifying model components (linear/logit/probit/ordered probit), defining latent factor structures (with any number of factors), and producing initialization values before passing to C++ for estimation.


---

## Features
- Define any number of latent factors with flexible loading normalization
- Specify model components independently (linear, logit, probit, oprobit)
- Automatically initialize parameters and compute approximate standard errors
- Export files for C++ back-end estimation.

## Quick start

```r
# devtools::load_all() inside factorana_R_interface
# Quick demo: Simulated system with continuous, binary, and ordered outcomes

library(factorana)

set.seed(42)
n <- 500
f  <- rnorm(n)
x1 <- rnorm(n); x2 <- rnorm(n)
z1 <- rnorm(n)
q1 <- rnorm(n)

# latent variable generating an ordered outcome
latent <- 0.5*x1 - 0.3*x2 + 0.8*f + rnorm(n)
y_ord <- cut(latent, breaks = quantile(latent, seq(0,1,0.25)), include.lowest = TRUE, ordered_result = TRUE)

Y  <- 1 + x1 + 0.8*f + rnorm(n)
D  <- as.integer(0.5*z1 + f + rnorm(n) > 0)
T  <- q1 + 0.5*f + rnorm(n)
dat <- data.frame(Y, D, T, x1, x2, z1, q1, y_ord)
dat$eval_y0 <- 1L - dat$D

# ---- Define system ----
fm  <- define_factor_model(2, 1, 8, loading_normalization = c(NA, 1))
mc_sel <- define_model_component("selection", dat, "D", fm, covariates = "z1", model_type = "probit")
mc_y1  <- define_model_component("Y1", dat, "Y", fm, evaluation_indicator = "D", covariates = c("x1","x2"), model_type = "linear")
mc_y0  <- define_model_component("Y0", dat, "Y", fm, evaluation_indicator = "eval_y0", covariates = c("x1","x2"), model_type = "linear")
mc_T   <- define_model_component("Tscore", dat, "T", fm, covariates = "q1", model_type = "linear")
mc_op  <- define_model_component("Satisfaction", dat, "y_ord", fm, covariates = "x1", model_type = "oprobit")

ms   <- define_model_system(factor = fm, components = list(mc_sel, mc_y1, mc_y0, mc_T, mc_op))
ctrl <- define_estimation_control(num_cores = 1)

# ---- Estimate and export ----
out <- estimate_and_write(ms, fm, ctrl, data = dat)

```

Typical printed output:
```
=== Running estimation ===
     component   intercept                               betas loading factor_var factor_cor
1    selection -0.03187048                     0.4646297187091   0.3;1          1          0
2           Y1  1.33174200 1.0760519470925;-0.0357340595369585   0.3;1          1          0
3           Y0  0.58041749 0.990587130317505;0.093872170432524   0.3;1          1          0
4       Tscore -0.04325684                   0.997703682546533   0.3;1          1          0
5 Satisfaction  0.00000000                                   0   0.3;1          1          0

Wrote all output files to /Users/.../factorana/factorana_R_interface/results 
```

- Number of latent factors = 2 -> length-2 `loading` for each component 
- Normalization `c(NA, 1)` fixes the second loading to 1; the first is initialized at default `0.3`

---

## Concepts

### Latent factors & normalization

- Choose `k = n_factors` in `define_factor_model()`.  
- Provide `loading_normalization` as a numeric vector of length `k`:
  - `NA`  → loading is **free** (initialized at a small base value).
  - value → loading is **fixed** to that numeric value (e.g., `0` or `1`).

Examples:
```r
# k = 1; fix loading at 1
define_factor_model(1, 1, 8, loading_normalization = 1)

# k = 3; mix of free/fixed
define_factor_model(3, 1, 8, loading_normalization = c(NA, 0, 1))
```

### Ordered probit (OProbit)

- Outcome **must be** an **ordered factor**.  
  If you pass integers, `define_model_component()` will coerce to an ordered factor with contiguous labels.
- With \(J\) categories, initialization returns \(J-1\) strictly increasing **thresholds**.

```r
# 4 categories -> 3 thresholds
stopifnot(is.numeric(ini$thresholds), length(ini$thresholds) == 3,
          all(diff(ini$thresholds) > 0))
```

---

## API (R)
More detailed explanations within functions.

### `define_factor_model(n_factors, n_types, n_quad_points, correlation = FALSE, n_mixtures = 1, loading_normalization = NULL)`
- `n_factors` (int ≥1): number of latent factors \(k\)
- `loading_normalization` (length-k numeric): `NA` = free, number = fixed  
- `correlation` logical, whether factors are correlated
- `n_quad_points` number of quadrature points used in integration
- Returns an object of class `"factor_model"`. Stores `n_factors`, `loading_normalization`, and derived counts (`nfac_param`, …).

### `define_model_component(name, data, outcome, factor, evaluation_indicator = NULL, covariates, model_type, intercept = TRUE, num_choices = 2, nrank = NULL)`
- Validates data (no missing in eval subset), coerces outcome for `oprobit` to ordered factors if needed.
- `model_type`: `"linear"`, `"probit"`, `"logit"`, `"oprobit"`.
- Returns a `"model_component"` with pointers to `factor` (includes normalization).

### `define_model_system(components, factor)`
- Bundles components and the shared factor model into a `"model_system"`.

### `define_estimation_control(num_cores = 1)`
- Simple container for runtime controls such as the number of parallel cores.

### `initialize_parameters(mc)`
- (Internal helper) Produces initial values for a **single** model component:
  - `intercept`, `betas`
  - `loading` (length-k; applies normalization)
  - For `oprobit`: `thresholds` (length \(J-1\))
  - `factor_var`, `factor_cor` placeholders

### `estimate_model(ms, control)`
- Current prototype: initializes each component and returns a small summary `data.frame`.  
  (The C++ estimation step will consume these later.)
  
### `estimate_and_write(ms, factor_model, control, data, outdir = "results")
- wrapper that calls estimate_model() and exports 4 files:
1. system_inits_long.csv : Long format list of component level parameters
2. meas_par.csv : Numeric parameter + SE table for the backend
3. model_config.csv : Key-value configuration table (system, factor, control)
4. simulated_data.csv : the dataset used in estimation 

---

## Examples

### k = 3 with mixed constraints
```r
set.seed(1)
n  <- 400
x1 <- rnorm(n); x2 <- rnorm(n)
p  <- pnorm(0.5*x1 - 0.4*x2)
y  <- as.integer(runif(n) < p)
dat <- data.frame(y=y, x1=x1, x2=x2, eval=1L)

fm3 <- define_factor_model(3, 1, 8, loading_normalization = c(NA, 0, 1))
mc3 <- define_model_component("prb3", dat, "y", fm3, "eval", c("x1","x2"), "probit")
ini3 <- initialize_parameters(mc3)

# Expect free/fixed: [*, 0, 1]
ini3$loading
```

### OProbit with 5 categories
```r
set.seed(2)
n  <- 350
x1 <- rnorm(n); e <- rnorm(n)
z  <- 0.8*x1 + e
Y5 <- cut(z, breaks = quantile(z, probs = seq(0,1,by=0.2)),
          include.lowest = TRUE, ordered_result = TRUE)
df5 <- data.frame(Y=Y5, x1=x1, eval=1L)

fm1 <- define_factor_model(1, 1, 8, loading_normalization = 1)
mc5 <- define_model_component("op5", df5, "Y", fm1, "eval", "x1", "oprobit")
ini5 <- initialize_parameters(mc5)

length(ini5$thresholds)  # 4
all(diff(ini5$thresholds) > 0)  # TRUE
ini5$loading   # 1
```

---

## Testing

Run all tests:
```r
devtools::test()
```

Run a subset (while iterating):
```r
devtools::test(filter = "modeltypes|multifactor|oprobit")
```

What’s covered:
- Validation of inputs & conditioning on evaluation subset
- Ordered probit thresholds (J-1, strictly increasing)
- Multi-factor loading normalization (length-k; fixed entries works)

---
