
# factorana

<!-- badges: start -->
<!-- badges: end -->

# factorana (R interface)

R front-end for specifying model components (linear/logit/probit/ordered probit), defining latent factor structures (with any number of factors), and producing initialization values before passing to C++ for estimation.


---

## Installation

Install directly from GitHub using `devtools` or `remotes`:

```r
# Install devtools if you don't have it
install.packages("devtools")

# Install factorana from GitHub
devtools::install_github("GregVeramendi/factorana", subdir = "factorana_R")

# Load the package
library(factorana)
```

**Requirements**:
- **R** ≥ 4.0.0
- **C++ compiler**:
  - **Windows**: Install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) (matches your R version)
  - **macOS**: Install Xcode Command Line Tools via `xcode-select --install`
  - **Linux**: Install `build-essential` (Ubuntu/Debian) or `gcc-c++` (Fedora/RHEL)
- **R packages**: Rcpp, RcppEigen (automatically installed as dependencies)

The package will automatically compile the C++ code during installation.

### Troubleshooting Installation

**Windows**: If you get compilation errors, make sure Rtools is installed and on your PATH:
```r
# Check if Rtools is found
Sys.which("make")
```

**macOS**: If you get "clang: error: unsupported option '-fopenmp'", this is expected and can be ignored (OpenMP is optional).

**All platforms**: If installation fails, try:
```r
# Install with verbose output to see errors
devtools::install_github("GregVeramendi/factorana",
                        subdir = "factorana_R",
                        build_vignettes = FALSE,
                        force = TRUE)
```

## Features
- Define any number of latent factors with flexible loading normalization
- Specify model components independently (linear, logit, probit, oprobit)
- **Dynamic factor models**: Structural equations between latent factors via `define_dyn_model_component()`
- **Nonlinear factor effects**: Model quadratic (`f²`) and interaction (`f_j × f_k`) factor terms via `factor_spec`
- Fix regression coefficients to specific values via `fix_coefficient()`
- Fix type-specific intercepts to zero via `fix_type_intercepts()` (for multi-type models)
- Multi-stage/sequential estimation with fixed early-stage parameters
- Automatically initialize parameters using component-by-component estimation
- Fast C++ backend with R-level parallelization for large datasets
- Analytical Hessian for fast convergence and accurate standard errors

## Quick start: Roy model example

This example demonstrates a Roy selection model with unobserved ability (latent factor), sector choice, test scores, and wages.

<details>
<summary><b>▶ Click to expand full example code</b></summary>

```r
library(factorana)

# Generate Roy model data
set.seed(108)
n <- 10000

# Covariates
x1 <- rnorm(n)  # Affects wages
x2 <- rnorm(n)  # Affects wages and sector choice
f <- rnorm(n)   # Latent ability (unobserved)

# Test scores (measure ability with error)
T1 <- 2.0 + 1.0*f + rnorm(n, 0, 0.5)
T2 <- 1.5 + 1.2*f + rnorm(n, 0, 0.6)
T3 <- 1.0 + 0.8*f + rnorm(n, 0, 0.4)

# Potential wages in each sector
wage0 <- 2.0 + 0.5*x1 + 0.3*x2 + 0.5*f + rnorm(n, 0, 0.6)     # Ability affects wage0
wage1 <- 2.5 + 0.6*x1 + 1.0*f + rnorm(n, 0, 0.7)              # Ability affects wage1

# Sector choice (high ability → more likely sector 1)
z_sector <- 0.0 + 0.4*x2 + 0.8*f
sector <- as.numeric(runif(n) < pnorm(z_sector))

# Observed wage (only see wage in chosen sector)
wage <- ifelse(sector == 1, wage1, wage0)

# Create dataset with evaluation indicators
dat <- data.frame(
  intercept = 1,
  x1 = x1, x2 = x2,
  T1 = T1, T2 = T2, T3 = T3,
  wage = wage,
  sector = sector,
  eval_tests = 1,                    # Always observe test scores
  eval_wage0 = 1 - sector,           # Observe wage0 when sector=0
  eval_wage1 = sector,               # Observe wage1 when sector=1
  eval_sector = 1                    # Always observe sector choice
)

# Define factor model (1 latent ability factor)
fm <- define_factor_model(n_factors = 1, n_types = 1)

# Define model components
# Test 1: Normalize loading to 1.0 for identification
mc_T1 <- define_model_component(
  name = "T1", data = dat, outcome = "T1", factor = fm,
  covariates = "intercept", model_type = "linear",
  loading_normalization = 1.0,
  evaluation_indicator = "eval_tests"
)

# Tests 2 and 3: Free loadings
mc_T2 <- define_model_component(
  name = "T2", data = dat, outcome = "T2", factor = fm,
  covariates = "intercept", model_type = "linear",
  loading_normalization = NA_real_,
  evaluation_indicator = "eval_tests"
)
mc_T3 <- define_model_component(
  name = "T3", data = dat, outcome = "T3", factor = fm,
  covariates = "intercept", model_type = "linear",
  loading_normalization = NA_real_,
  evaluation_indicator = "eval_tests"
)

# Wage in sector 0: Ability effect (free loading)
mc_wage0 <- define_model_component(
  name = "wage0", data = dat, outcome = "wage", factor = fm,
  covariates = c("intercept", "x1", "x2"), model_type = "linear",
  loading_normalization = NA_real_,
  evaluation_indicator = "eval_wage0"
)

# Wage in sector 1: Ability matters (free loading)
mc_wage1 <- define_model_component(
  name = "wage1", data = dat, outcome = "wage", factor = fm,
  covariates = c("intercept", "x1"), model_type = "linear",
  loading_normalization = NA_real_,
  evaluation_indicator = "eval_wage1"
)

# Sector choice: Probit model
mc_sector <- define_model_component(
  name = "sector", data = dat, outcome = "sector", factor = fm,
  covariates = c("intercept", "x2"), model_type = "probit",
  loading_normalization = NA_real_,
  evaluation_indicator = "eval_sector"
)

# Define model system
ms <- define_model_system(
  components = list(mc_T1, mc_T2, mc_T3, mc_wage0, mc_wage1, mc_sector),
  factor = fm
)

# Single-core estimation
# Note: init_params = NULL triggers automatic initialization
ctrl_single <- define_estimation_control(n_quad_points = 16, num_cores = 1)
time_start_single <- Sys.time()
result_single <- estimate_model_rcpp(
  model_system = ms,
  data = dat,
  init_params = NULL,       # Automatic initialization (recommended)
  control = ctrl_single,
  parallel = FALSE,
  optimizer = "nlminb",     # Default: fast with analytical Hessian
  verbose = TRUE
)
time_end_single <- Sys.time()
time_single <- as.numeric(difftime(time_end_single, time_start_single, units = "secs"))

# Parallel estimation with 4 cores (3x speedup on large datasets)
ctrl_parallel <- define_estimation_control(n_quad_points = 16, num_cores = 4)
time_start_parallel <- Sys.time()
result_parallel <- estimate_model_rcpp(
  model_system = ms,
  data = dat,
  init_params = NULL,       # Automatic initialization
  control = ctrl_parallel,
  parallel = TRUE,
  optimizer = "nlminb",
  verbose = TRUE
)
time_end_parallel <- Sys.time()
time_parallel <- as.numeric(difftime(time_end_parallel, time_start_parallel, units = "secs"))

# Create formatted results tables
# Define true parameter values
true_params <- c(
  1.0,      # factor_var (fixed to estimates)
  2.0, 0.5, # T1: intercept, sigma
  1.5, 1.2, 0.6, # T2: intercept, loading, sigma
  1.0, 0.8, 0.4, # T3: intercept, loading, sigma
  2.0, 0.5, 0.3, 0.5, 0.6, # wage0: intercept, x1, x2, loading, sigma
  2.5, 0.6, 1.0, 0.7, # wage1: intercept, x1, loading, sigma
  0.0, 0.4, 0.8  # sector: intercept, x2, loading
)

# Update factor variance to match estimate
true_params[1] <- result_parallel$estimates[1]

# Parameter names
param_names <- c(
  "factor_var",
  "T1_intercept", "T1_sigma",
  "T2_intercept", "T2_loading", "T2_sigma",
  "T3_intercept", "T3_loading", "T3_sigma",
  "wage0_intercept", "wage0_x1", "wage0_x2", "wage0_loading", "wage0_sigma",
  "wage1_intercept", "wage1_x1", "wage1_loading", "wage1_sigma",
  "sector_intercept", "sector_x2", "sector_loading"
)

# Component labels for grouping
components <- c(
  "Factor",
  "T1", "T1",
  "T2", "T2", "T2",
  "T3", "T3", "T3",
  "wage0", "wage0", "wage0", "wage0", "wage0",
  "wage1", "wage1", "wage1", "wage1",
  "sector", "sector", "sector"
)

# Table 1: Parameter estimates
results_table <- data.frame(
  Component = components,
  Parameter = param_names,
  True = sprintf("%.3f", true_params),
  Estimate = sprintf("%.3f", result_parallel$estimates),
  Std_Error = sprintf("%.3f", result_parallel$std_errors)
)

cat("\n=== Parameter Estimates ===\n")
print(results_table, row.names = FALSE, right = FALSE)

# Table 2: Estimation diagnostics
speedup <- time_single / time_parallel
diagnostics <- data.frame(
  Method = c("Single-core", "Parallel (4 cores)"),
  Log_Likelihood = sprintf("%.2f", c(result_single$loglik, result_parallel$loglik)),
  Time_sec = sprintf("%.2f", c(time_single, time_parallel)),
  Speedup = c("1.0x", sprintf("%.2fx", speedup)),
  Convergence = c(result_single$convergence, result_parallel$convergence),
  N_Parameters = c(length(result_single$estimates), length(result_parallel$estimates))
)

cat("\n=== Estimation Diagnostics ===\n")
print(diagnostics, row.names = FALSE, right = FALSE)
```

</details>

**Key features demonstrated:**

- **Latent factor**: Unobserved ability (`f`) affects test scores, wages, and sector choice
- **Loading normalization**:
  - `T1`: Fixed to 1.0 for identification (first test score normalizes scale)
  - `T2, T3, wage0, wage1, sector`: Free parameters (estimated)
- **Evaluation indicators**: Only observe wages in chosen sector (partial observability)
- **Automatic initialization**: `init_params = NULL` uses smart component-by-component initialization
- **Parallelization**: 4 cores provides ~3x speedup on large datasets (n=10,000)
- **Optimizer**: Default `nlminb` uses analytical Hessian for fast convergence

---

## Two-Stage (Sequential) Estimation

For complex models, you can estimate in multiple stages to improve convergence and interpretability. In the first stage, estimate a subset of components (e.g., measurement system). In subsequent stages, fix those components and add new ones.

<details>
<summary><b>▶ Click to expand two-stage estimation example</b></summary>

**Example**: Estimate test scores first, then add wage and sector equations:

```r
# Using the same Roy model data from above...

# ======================================================================
# STAGE 1: Estimate measurement system (test scores only)
# ======================================================================

fm <- define_factor_model(n_factors = 1, n_types = 1)
ctrl <- define_estimation_control(n_quad_points = 16, num_cores = 1)

# Define three test score components
mc_T1 <- define_model_component(
  name = "T1", data = dat, outcome = "T1", factor = fm,
  covariates = "intercept", model_type = "linear",
  loading_normalization = 1.0,
  evaluation_indicator = "eval_tests"
)

mc_T2 <- define_model_component(
  name = "T2", data = dat, outcome = "T2", factor = fm,
  covariates = "intercept", model_type = "linear",
  loading_normalization = NA_real_,
  evaluation_indicator = "eval_tests"
)

mc_T3 <- define_model_component(
  name = "T3", data = dat, outcome = "T3", factor = fm,
  covariates = "intercept", model_type = "linear",
  loading_normalization = NA_real_,
  evaluation_indicator = "eval_tests"
)

# Create model system for stage 1
ms_stage1 <- define_model_system(
  components = list(mc_T1, mc_T2, mc_T3),
  factor = fm
)

# Estimate stage 1
result_stage1 <- estimate_model_rcpp(
  model_system = ms_stage1,
  data = dat,
  init_params = NULL,
  control = ctrl,
  optimizer = "nlminb",
  verbose = TRUE
)

# ======================================================================
# STAGE 2: Add wage/sector equations, fixing stage 1 parameters
# ======================================================================

# Define wage and sector components (same as before)
mc_wage0 <- define_model_component(
  name = "wage0", data = dat, outcome = "wage", factor = fm,
  covariates = c("intercept", "x1", "x2"), model_type = "linear",
  loading_normalization = NA_real_,
  evaluation_indicator = "eval_wage0"
)

mc_wage1 <- define_model_component(
  name = "wage1", data = dat, outcome = "wage", factor = fm,
  covariates = c("intercept", "x1"), model_type = "linear",
  loading_normalization = NA_real_,
  evaluation_indicator = "eval_wage1"
)

mc_sector <- define_model_component(
  name = "sector", data = dat, outcome = "sector", factor = fm,
  covariates = c("intercept", "x2"), model_type = "probit",
  loading_normalization = NA_real_,
  evaluation_indicator = "eval_sector"
)

# Create model system for stage 2, passing stage 1 results
ms_stage2 <- define_model_system(
  components = list(mc_wage0, mc_wage1, mc_sector),
  factor = fm,
  previous_stage = result_stage1  # Fix stage 1 parameters
)

# Estimate stage 2
result_stage2 <- estimate_model_rcpp(
  model_system = ms_stage2,
  data = dat,
  init_params = NULL,
  control = ctrl,
  optimizer = "nlminb",
  verbose = TRUE
)

# View combined results
print(result_stage2$estimates)    # Includes both stage 1 and stage 2 parameters
print(result_stage2$std_errors)   # Standard errors preserved from stage 1
```

</details>

**How it works:**
- Stage 1 estimates only the measurement system (9 parameters)
- Stage 2 fixes those 9 parameters and estimates wage/sector parameters (12 new parameters)
- `result_stage2` contains all 21 parameters with correct standard errors
- Factor variance is fixed at stage 1 estimate
- Computational efficiency: gradients/Hessians not computed for fixed components

**Benefits:**
- Improved convergence for complex models
- Interpretable intermediate results
- Faster optimization (skip derivatives for fixed parameters)
- Can chain multiple stages: `stage3 <- define_model_system(..., previous_stage = result_stage2)`

---

## Concepts

### Latent factors & normalization

- Choose `k = n_factors` in `define_factor_model()`.
- Specify `loading_normalization` at the **component level** via `define_model_component()`:
  - `NA` or `NA_real_` → loading is **free** (estimated parameter)
  - Numeric value → loading is **fixed** to that value (e.g., `0` or `1.0`)

**Identification**: At least one component must have a fixed loading to normalize the factor scale.

Examples:
```r
# Single factor model
fm <- define_factor_model(n_factors = 1, n_types = 1)
ctrl <- define_estimation_control(n_quad_points = 16, num_cores = 1)

# First component: fix loading to 1.0 for identification
mc1 <- define_model_component(
  name = "test1", data = dat, outcome = "y1", factor = fm,
  covariates = "intercept", model_type = "linear",
  loading_normalization = 1.0
)

# Other components: free loadings
mc2 <- define_model_component(
  name = "test2", data = dat, outcome = "y2", factor = fm,
  covariates = "intercept", model_type = "linear",
  loading_normalization = NA_real_  # Free parameter
)
```

### Ordered probit (OProbit)

- Outcome **must be** an **ordered factor**.
  If you pass integers, `define_model_component()` will coerce to an ordered factor with contiguous labels.
- With \(J\) categories, the model estimates \(J-1\) **threshold increments**.
- The first threshold can be any value, subsequent thresholds are positive increments.
- During likelihood evaluation, thresholds are accumulated: thresh_k = thresh_1 + |incr_2| + ... + |incr_k|

---

## API (R)
More detailed explanations within functions.

### `define_factor_model(n_factors, n_types, correlation = FALSE, n_mixtures = 1)`
- `n_factors` (int ≥0): number of latent factors (use 0 for models without factors)
- `n_types` (int ≥1): number of types
- `correlation` (logical): whether factors are correlated (default: FALSE)
- `n_mixtures` (int 1-3): number of discrete mixtures (default: 1)
- Returns an object of class `"factor_model"`
- **Note**: Loading normalization is specified at the component level via `define_model_component()`. Quadrature points are specified in `define_estimation_control()`.

### `define_model_component(name, data, outcome, factor, evaluation_indicator = NULL, covariates, model_type, loading_normalization = NA_real_, factor_spec = "linear", intercept = TRUE, num_choices = 2, nrank = NULL)`
- Validates data (no missing in eval subset), coerces outcome for `oprobit` to ordered factors if needed.
- `model_type`: `"linear"`, `"probit"`, `"logit"`, `"oprobit"`.
- `loading_normalization`: Normalization for factor loadings (NA or NA_real_ = free, numeric = fixed)
- `factor_spec`: Factor specification for nonlinear effects (default: `"linear"`)
  - `"linear"`: Standard linear factor terms only (λ × f)
  - `"quadratic"`: Include quadratic terms (λ × f + λ_quad × f²)
  - `"interactions"`: Include interaction terms (λ × f + λ_inter × f_j × f_k) for multi-factor models
  - `"full"`: Include both quadratic and interaction terms
- Returns a `"model_component"` with pointers to `factor`.

### `define_dyn_model_component(name, data, outcome_factor, factor, covariates = NULL, dyn_type = "linear", factor_spec = "linear", intercept = TRUE, evaluation_indicator = NULL)`
- Creates a **dynamic factor model component** representing a structural equation between latent factors.
- The model estimates: `f_outcome = X'β + Σ_{k≠outcome} λ_k × f_k + [higher-order terms] + ε`
- **Parameters**:
  - `name`: Name for the component
  - `data`: Data frame (only used for covariates and evaluation_indicator)
  - `outcome_factor`: Integer index (1-based) of the outcome factor
  - `factor`: Factor model from `define_factor_model()`
  - `covariates`: Optional character vector of covariate names
  - `dyn_type`: Currently only `"linear"` supported (future: `"cobb_douglas"`, `"ces"`)
  - `factor_spec`: Factor specification for non-outcome factors
    - `"linear"`: Linear factor terms only
    - `"quadratic"`: Include f² terms for non-outcome factors
    - `"interactions"`: Include f_j × f_k terms (excluding outcome factor)
    - `"full"`: Both quadratic and interaction terms
  - `intercept`: Whether to include an intercept (default: TRUE)
  - `evaluation_indicator`: Optional indicator variable for which observations to include
- **Returns**: Object of class `c("dyn_model_component", "model_component")`
- **Key differences from `define_model_component()`**:
  - No `outcome` variable needed (outcome is implicitly zero)
  - Loading on outcome factor is fixed to -1 (not estimated)
  - Factor normalization is automatic: `-1` for outcome factor, `NA` (free) for others
  - Quadratic/interaction terms only apply to non-outcome factors
- **Example**:
  ```r
  # 2-factor model: f1 = intercept + lambda*f2 + epsilon
  fm <- define_factor_model(n_factors = 2)
  dyn <- define_dyn_model_component(
    name = "structural", data = dat, outcome_factor = 1, factor = fm,
    intercept = TRUE
  )
  ```

### `fix_coefficient(component, covariate, value, choice = NULL)`
- Fixes a regression coefficient to a specified value during estimation.
- **Parameters**:
  - `component`: Model component from `define_model_component()`
  - `covariate`: Name of the covariate whose coefficient should be fixed
  - `value`: Numeric value to fix the coefficient to
  - `choice`: For multinomial logit only - which choice (1-indexed, excluding reference)
- **Returns**: Modified model component with the constraint added
- **Example**:
  ```r
  # Fix intercept to 0
  mc <- fix_coefficient(mc, "intercept", 0)

  # Fix x1 coefficient to a specific value
  mc <- fix_coefficient(mc, "x1", 1.5)

  # For multinomial logit: fix x1 for first non-reference choice
  mc <- fix_coefficient(mc, "x1", 0, choice = 1)
  ```
- **Notes**:
  - Multiple coefficients can be fixed by chaining calls
  - Only regression coefficients (betas) can be fixed, not sigma, thresholds, or loadings
  - Fixed coefficients are held constant during optimization

### `fix_type_intercepts(component, types = NULL, choice = NULL)`
- Fixes type-specific intercepts to zero for models with `n_types > 1`.
- **Parameters**:
  - `component`: Model component from `define_model_component()`
  - `types`: Integer vector of which types to fix (default: all non-reference types, i.e., types 2 through n_types)
  - `choice`: For multinomial logit only - which choice's type intercepts to fix (default: all choices)
- **Returns**: Modified model component with the constraints added
- **Example**:
  ```r
  # Create 2-type factor model
  fm <- define_factor_model(n_factors = 1, n_types = 2)
  mc <- define_model_component(name = "Y", data = dat, ...)

  # Fix all type intercepts to 0
  mc <- fix_type_intercepts(mc)

  # Or fix specific type(s) in a 3-type model
  mc <- fix_type_intercepts(mc, types = 2)  # Only fix type 2

  # For multinomial logit: fix type intercepts for specific choice
  mc <- fix_type_intercepts(mc, choice = 1)
  ```
- **Notes**:
  - Type 1 is the reference type (no intercept parameter)
  - Fixed type intercepts are excluded from the parameter vector entirely
  - Useful when type effects should operate only through factor loadings, not direct intercept shifts

### `define_model_system(components, factor, previous_stage = NULL)`
- Bundles components and the shared factor model into a `"model_system"`.
- `previous_stage` (optional): Result object from a previous `estimate_model_rcpp()` call
  - Enables multi-stage/sequential estimation
  - Fixes all previous-stage component parameters and factor variance
  - Previous-stage components are prepended to new components
  - Standard errors are preserved from previous stage
  - See "Two-Stage Estimation" section for usage example

### `define_estimation_control(n_quad_points = 16, num_cores = 1)`
- Container for estimation settings including numerical integration and parallelization.
- `n_quad_points` (int ≥1): Number of Gauss-Hermite quadrature points for numerical integration (default: 16)
- `num_cores`: Number of CPU cores to use for parallel estimation (default: 1)

### `estimate_model_rcpp(model_system, data, init_params = NULL, control = NULL, optimizer = "nlminb", parallel = TRUE, verbose = TRUE)`
- Main estimation function using C++ backend with R-level parallelization
- **Parameters**:
  - `model_system`: Output from `define_model_system()`
  - `data`: Data frame containing all variables
  - `init_params`: Optional initial parameter values (auto-initialized if NULL)
  - `control`: Output from `define_estimation_control()` (uses 1 core if NULL)
  - `optimizer`: Optimization algorithm (default: `"nlminb"`)
    - `"nlminb"`: Fast, uses analytical Hessian (recommended)
    - `"nloptr"`: L-BFGS with bounds (slower, ~4.6x)
    - `"optim"`: L-BFGS-B (similar to nloptr)
    - `"trust"`: Trust region (experimental)
  - `parallel`: Enable parallelization when `num_cores > 1` (default: TRUE)
  - `verbose`: Print progress messages (default: TRUE)
- **Returns**: List with
  - `estimates`: Parameter estimates
  - `std_errors`: Standard errors (from Hessian)
  - `loglik`: Log-likelihood at optimum
  - `convergence`: Convergence code (0 = success)
  - `model_system`: Original model system
  - `optimizer`: Optimizer used

**Parallelization notes**:
- Splits observations across workers (not parameters)
- Each worker evaluates likelihood on its data subset
- Results are aggregated across workers
- Achieves ~3x speedup with 4 cores on large datasets
- Uses `doParallel` for Windows compatibility
- Single evaluation indicator per observation (one worker evaluates it)

### `initialize_parameters(model_system, data, verbose = TRUE)`
- Produces initial values for all parameters in the model system
- Estimates each component separately (ignoring factors) to get starting values
- **Returns**: List with
  - `init_params`: Initial parameter vector
  - `factor_variance_fixed`: Whether factor variance is identified

---

## Displaying Results

The package provides functions to display estimation results in formatted tables, either to the screen or as LaTeX output.

### `components_table(result, components = NULL)`
- Creates a table with model components as columns and parameter types as rows
- Similar to SEM/CFA output format for easy comparison across components
- `result`: Output from `estimate_model_rcpp()`
- `components`: Optional character vector of component names to include (default: all)
- **Returns**: Object of class `"components_table"` that prints nicely

### `components_to_latex(result, components = NULL, caption = NULL, label = NULL, digits = 3)`
- Exports the component table to LaTeX format (booktabs style)
- Same parameters as `components_table()`, plus:
  - `caption`: Optional table caption
  - `label`: Optional LaTeX label for cross-referencing
  - `digits`: Number of decimal places (default: 3)
- **Returns**: Character string with LaTeX code

### Example: Displaying All Components

After estimating a Roy model (see Quick Start example):

```r
# Display all components in one table
components_table(result_parallel)
```

Output:
```
Factor Model Results by Component
=========================================================================================================

Parameter          sector         wage1         wage0            T1            T2            T3
---------------------------------------------------------------------------------------------
beta_intercept                   2.531***      2.013***      1.998***      1.500***      0.994***
                                  (0.043)       (0.035)       (0.012)       (0.016)       (0.011)
beta_x1                          0.594***      0.498***
                                  (0.048)       (0.035)
beta_x2             0.416***                   0.299***
                     (0.044)                    (0.030)
Loading 1           0.807***      0.983***      0.521***      1.000         1.191***      0.797***
                     (0.058)       (0.037)       (0.038)         (-)        (0.021)       (0.016)
Sigma                             0.697***      0.609***      0.501***      0.600***      0.399***
                                   (0.027)       (0.024)       (0.011)       (0.013)       (0.009)
---------------------------------------------------------------------------------------------------------
N                    10000          5012          4988         10000         10000         10000
Log-likelihood: -42156.78

Standard errors in parentheses
Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1
```

### Example: Splitting Components into Separate Tables

For complex models, you may want to display the measurement system separately from the structural equations:

```r
# Table 1: Measurement system (test scores)
components_table(result_parallel, components = c("T1", "T2", "T3"))

# Table 2: Selection and outcome equations
components_table(result_parallel, components = c("sector", "wage0", "wage1"))
```

**Measurement System Table:**
```
Factor Model Results by Component
=============================================================================

Parameter               T1            T2            T3
-------------------------------------------------------------
beta_intercept      1.998***      1.500***      0.994***
                     (0.012)       (0.016)       (0.011)
Loading 1           1.000         1.191***      0.797***
                       (-)         (0.021)       (0.016)
Sigma               0.501***      0.600***      0.399***
                     (0.011)       (0.013)       (0.009)
-----------------------------------------------------------------------------
N                    10000         10000         10000
Log-likelihood: -42156.78

Standard errors in parentheses
Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1
```

**Selection and Outcomes Table:**
```
Factor Model Results by Component
=============================================================================

Parameter          sector         wage1         wage0
-------------------------------------------------------------
beta_intercept                   2.531***      2.013***
                                  (0.043)       (0.035)
beta_x1                          0.594***      0.498***
                                  (0.048)       (0.035)
beta_x2             0.416***                   0.299***
                     (0.044)                    (0.030)
Loading 1           0.807***      0.983***      0.521***
                     (0.058)       (0.037)       (0.038)
Sigma                             0.697***      0.609***
                                   (0.027)       (0.024)
-----------------------------------------------------------------------------
N                    10000          5012          4988
Log-likelihood: -42156.78

Standard errors in parentheses
Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1
```

### Example: LaTeX Export

```r
# Export measurement system to LaTeX
latex_code <- components_to_latex(
  result_parallel,
  components = c("T1", "T2", "T3"),
  caption = "Measurement System Estimates",
  label = "tab:measurement"
)
cat(latex_code)

# Export selection/outcome to LaTeX
latex_code2 <- components_to_latex(
  result_parallel,
  components = c("sector", "wage0", "wage1"),
  caption = "Selection and Outcome Equations",
  label = "tab:structural"
)
writeLines(latex_code2, "structural_table.tex")
```

### Other Display Functions

- `print(result)` - Basic summary of estimation results
- `summary(result)` - Detailed summary with z-values and p-values
- `results_table(result1, result2, ...)` - Compare multiple model results side by side
- `results_to_latex(result1, result2, ...)` - Export multi-model comparison to LaTeX

---

## Examples

### Dynamic Factor Model (Structural Equation between Factors)

Model structural relationships between latent factors, e.g., `f1 = β₀ + λ₂×f2 + ε`:

<details>
<summary><b>▶ Click to expand dynamic factor model example</b></summary>

```r
set.seed(111)
n <- 500

# Two latent factors with correlation
f1 <- rnorm(n)
f2 <- rnorm(n)

# Measurement system (test scores to identify factors)
# Factor 1 indicators
T1 <- 2.0 + 1.0*f1 + rnorm(n, 0, 0.5)  # Loading fixed to 1
T2 <- 1.5 + 1.2*f1 + rnorm(n, 0, 0.6)

# Factor 2 indicators
T3 <- 1.0 + 1.0*f2 + rnorm(n, 0, 0.4)  # Loading fixed to 1
T4 <- 0.8 + 0.9*f2 + rnorm(n, 0, 0.5)

dat <- data.frame(intercept = 1, T1 = T1, T2 = T2, T3 = T3, T4 = T4, eval = 1)

# Define 2-factor model
fm <- define_factor_model(n_factors = 2, n_types = 1)
ctrl <- define_estimation_control(n_quad_points = 8, num_cores = 1)

# Factor 1 measurement equations (identify factor 1)
mc_T1 <- define_model_component(
  name = "T1", data = dat, outcome = "T1", factor = fm,
  covariates = "intercept", model_type = "linear",
  loading_normalization = c(1.0, 0),  # f1 loading=1, f2 loading=0
  evaluation_indicator = "eval"
)

mc_T2 <- define_model_component(
  name = "T2", data = dat, outcome = "T2", factor = fm,
  covariates = "intercept", model_type = "linear",
  loading_normalization = c(NA_real_, 0),  # f1 free, f2=0
  evaluation_indicator = "eval"
)

# Factor 2 measurement equations (identify factor 2)
mc_T3 <- define_model_component(
  name = "T3", data = dat, outcome = "T3", factor = fm,
  covariates = "intercept", model_type = "linear",
  loading_normalization = c(0, 1.0),  # f1=0, f2 loading=1
  evaluation_indicator = "eval"
)

mc_T4 <- define_model_component(
  name = "T4", data = dat, outcome = "T4", factor = fm,
  covariates = "intercept", model_type = "linear",
  loading_normalization = c(0, NA_real_),  # f1=0, f2 free
  evaluation_indicator = "eval"
)

# Dynamic structural equation: f1 = intercept + lambda*f2 + epsilon
# This models how factor 1 depends on factor 2
dyn <- define_dyn_model_component(
  name = "structural",
  data = dat,
  outcome_factor = 1,     # f1 is the "outcome" factor
  factor = fm,
  intercept = TRUE,       # Include intercept
  evaluation_indicator = "eval"
)

# Estimate
ms <- define_model_system(components = list(mc_T1, mc_T2, mc_T3, mc_T4, dyn), factor = fm)
result <- estimate_model_rcpp(ms, dat, init_params = NULL, control = ctrl, verbose = TRUE)

# View results
# structural_loading_2: effect of f2 on f1
# structural_intercept: intercept in structural equation
# structural_sigma: error variance in structural equation
print(result$estimates)
print(result$std_errors)
```

</details>

**Key points:**
- `define_dyn_model_component()` creates structural equations between factors
- `outcome_factor = 1` means f1 is the dependent variable
- Loading on f1 is automatically fixed to -1 (identity for outcome)
- Loading on f2 (λ₂) is estimated as `structural_loading_2`
- Works with `factor_spec = "quadratic"` or `"interactions"` for nonlinear effects between factors
- Measurement equations are needed to identify the factors

---

### Measurement System with Three Tests

<details>
<summary><b>▶ Click to expand example</b></summary>

A standard factor analysis setup with three test scores measuring a latent ability factor:

```r
set.seed(104)
n <- 500
f <- rnorm(n)  # Latent ability

# Three test scores (T1 loading fixed to 1 for identification)
T1 <- 2.0 + 1.0*f + rnorm(n, 0, 0.5)
T2 <- 1.5 + 1.2*f + rnorm(n, 0, 0.6)
T3 <- 1.0 + 0.8*f + rnorm(n, 0, 0.4)

dat <- data.frame(intercept = 1, T1 = T1, T2 = T2, T3 = T3, eval = 1)

# Define factor model and estimation control
fm <- define_factor_model(n_factors = 1, n_types = 1)
ctrl <- define_estimation_control(n_quad_points = 8, num_cores = 1)

# Define three test components
mc_T1 <- define_model_component(
  name = "T1", data = dat, outcome = "T1", factor = fm,
  covariates = "intercept", model_type = "linear",
  loading_normalization = 1.0,  # Fix for identification
  evaluation_indicator = "eval"
)

mc_T2 <- define_model_component(
  name = "T2", data = dat, outcome = "T2", factor = fm,
  covariates = "intercept", model_type = "linear",
  loading_normalization = NA_real_,  # Free parameter
  evaluation_indicator = "eval"
)

mc_T3 <- define_model_component(
  name = "T3", data = dat, outcome = "T3", factor = fm,
  covariates = "intercept", model_type = "linear",
  loading_normalization = NA_real_,  # Free parameter
  evaluation_indicator = "eval"
)

# Estimate
ms <- define_model_system(components = list(mc_T1, mc_T2, mc_T3), factor = fm)
result <- estimate_model_rcpp(ms, dat, init_params = NULL, control = ctrl, verbose = TRUE)

print(result$estimates)
print(result$std_errors)
```

</details>

### Ordered Probit with Measurement System

<details>
<summary><b>▶ Click to expand example</b></summary>

Combining test scores with an ordered outcome (e.g., educational attainment):

```r
set.seed(106)
n <- 500
x1 <- rnorm(n)
f <- rnorm(n)  # Latent ability

# Three test scores (same as above)
T1 <- 2.0 + 1.0*f + rnorm(n, 0, 0.5)
T2 <- 1.5 + 1.2*f + rnorm(n, 0, 0.6)
T3 <- 1.0 + 0.8*f + rnorm(n, 0, 0.4)

# Ordered outcome with 3 categories
# Note: Latent z = 0.5 (intercept) + 0.6*x1 + 0.8*f + error
# Intercept absorbed into thresholds, so model estimates only beta and loading
z <- 0.5 + 0.6*x1 + 0.8*f + rnorm(n)
y <- cut(z, breaks = c(-Inf, -0.5, 0.5, Inf), labels = FALSE)

dat <- data.frame(intercept = 1, x1 = x1, T1 = T1, T2 = T2, T3 = T3, y = y, eval = 1)

# Define factor model and estimation control
fm <- define_factor_model(n_factors = 1, n_types = 1)
ctrl <- define_estimation_control(n_quad_points = 8, num_cores = 1)

mc_T1 <- define_model_component(
  name = "T1", data = dat, outcome = "T1", factor = fm,
  covariates = "intercept", model_type = "linear",
  loading_normalization = 1.0, evaluation_indicator = "eval"
)

mc_T2 <- define_model_component(
  name = "T2", data = dat, outcome = "T2", factor = fm,
  covariates = "intercept", model_type = "linear",
  loading_normalization = NA_real_, evaluation_indicator = "eval"
)

mc_T3 <- define_model_component(
  name = "T3", data = dat, outcome = "T3", factor = fm,
  covariates = "intercept", model_type = "linear",
  loading_normalization = NA_real_, evaluation_indicator = "eval"
)

# Ordered probit component (NO intercept - absorbed into thresholds)
mc_y <- define_model_component(
  name = "y", data = dat, outcome = "y", factor = fm,
  covariates = "x1",  # No intercept for oprobit
  model_type = "oprobit",
  num_choices = 3,  # 3 ordered categories
  loading_normalization = NA_real_,
  evaluation_indicator = "eval"
)

# Estimate
ms <- define_model_system(components = list(mc_T1, mc_T2, mc_T3, mc_y), factor = fm)
result <- estimate_model_rcpp(ms, dat, init_params = NULL, control = ctrl, verbose = TRUE)

print(result$estimates)
```

</details>

### Multinomial Logit with Measurement System

<details>
<summary><b>▶ Click to expand example</b></summary>

Combining test scores with a discrete choice outcome (3 alternatives):

```r
set.seed(107)
n <- 500
x1 <- rnorm(n)
f <- rnorm(n)  # Latent ability

# Three test scores (same as above)
T1 <- 2.0 + 1.0*f + rnorm(n, 0, 0.5)
T2 <- 1.5 + 1.2*f + rnorm(n, 0, 0.6)
T3 <- 1.0 + 0.8*f + rnorm(n, 0, 0.4)

# Multinomial choice with 3 alternatives (choice 0 is reference)
z1 <- 0.5 + 0.6*x1 + 0.7*f
z2 <- 1.0 - 0.5*x1 + 0.9*f

exp_z0 <- 1
exp_z1 <- exp(z1)
exp_z2 <- exp(z2)
denom <- exp_z0 + exp_z1 + exp_z2

p0 <- exp_z0 / denom
p1 <- exp_z1 / denom
p2 <- exp_z2 / denom

y <- numeric(n)
for (i in seq_len(n)) {
  # C++ expects choices coded as 1, 2, 3 (not 0, 1, 2)
  y[i] <- sample(1:3, 1, prob = c(p0[i], p1[i], p2[i]))
}

dat <- data.frame(intercept = 1, x1 = x1, T1 = T1, T2 = T2, T3 = T3, y = y, eval = 1)

# Define factor model and estimation control
fm <- define_factor_model(n_factors = 1, n_types = 1)
ctrl <- define_estimation_control(n_quad_points = 8, num_cores = 1)

mc_T1 <- define_model_component(
  name = "T1", data = dat, outcome = "T1", factor = fm,
  covariates = "intercept", model_type = "linear",
  loading_normalization = 1.0, evaluation_indicator = "eval"
)

mc_T2 <- define_model_component(
  name = "T2", data = dat, outcome = "T2", factor = fm,
  covariates = "intercept", model_type = "linear",
  loading_normalization = NA_real_, evaluation_indicator = "eval"
)

mc_T3 <- define_model_component(
  name = "T3", data = dat, outcome = "T3", factor = fm,
  covariates = "intercept", model_type = "linear",
  loading_normalization = NA_real_, evaluation_indicator = "eval"
)

# Multinomial logit component
mc_y <- define_model_component(
  name = "y", data = dat, outcome = "y", factor = fm,
  covariates = c("intercept", "x1"),
  model_type = "logit",
  num_choices = 3,  # 3 alternatives
  loading_normalization = NA_real_,
  evaluation_indicator = "eval"
)

# Estimate
ms <- define_model_system(components = list(mc_T1, mc_T2, mc_T3, mc_y), factor = fm)
result <- estimate_model_rcpp(ms, dat, init_params = NULL, control = ctrl, verbose = TRUE)

print(result$estimates)
```

</details>

### Quadratic Factor Effects

<details>
<summary><b>▶ Click to expand example</b></summary>

When factor effects are nonlinear, use `factor_spec = "quadratic"` to include f² terms:

```r
set.seed(109)
n <- 500
f <- rnorm(n)  # Latent ability
x1 <- rnorm(n)

# Three test scores (standard measurement system)
T1 <- 2.0 + 1.0*f + rnorm(n, 0, 0.5)
T2 <- 1.5 + 1.2*f + rnorm(n, 0, 0.6)
T3 <- 1.0 + 0.8*f + rnorm(n, 0, 0.4)

# Outcome with quadratic factor effect
# Y = intercept + beta*x1 + lambda*f + lambda_quad*f^2 + error
Y <- 3.0 + 0.5*x1 + 0.8*f + 0.3*f^2 + rnorm(n, 0, 0.5)

dat <- data.frame(intercept = 1, x1 = x1, T1 = T1, T2 = T2, T3 = T3, Y = Y, eval = 1)

# Define factor model and estimation control
fm <- define_factor_model(n_factors = 1, n_types = 1)
ctrl <- define_estimation_control(n_quad_points = 8, num_cores = 1)

# Standard measurement equations (linear factor terms)
mc_T1 <- define_model_component(
  name = "T1", data = dat, outcome = "T1", factor = fm,
  covariates = "intercept", model_type = "linear",
  loading_normalization = 1.0, evaluation_indicator = "eval"
)

mc_T2 <- define_model_component(
  name = "T2", data = dat, outcome = "T2", factor = fm,
  covariates = "intercept", model_type = "linear",
  loading_normalization = NA_real_, evaluation_indicator = "eval"
)

mc_T3 <- define_model_component(
  name = "T3", data = dat, outcome = "T3", factor = fm,
  covariates = "intercept", model_type = "linear",
  loading_normalization = NA_real_, evaluation_indicator = "eval"
)

# Outcome with QUADRATIC factor effect
mc_Y <- define_model_component(
  name = "Y", data = dat, outcome = "Y", factor = fm,
  covariates = c("intercept", "x1"), model_type = "linear",
  loading_normalization = NA_real_,
  factor_spec = "quadratic",  # Enable f^2 term
  evaluation_indicator = "eval"
)

# Estimate
ms <- define_model_system(components = list(mc_T1, mc_T2, mc_T3, mc_Y), factor = fm)
result <- estimate_model_rcpp(ms, dat, init_params = NULL, control = ctrl, verbose = TRUE)

# View results - includes Y_loading_quad_1 parameter
print(result$estimates)
```

</details>

**Key points:**
- `factor_spec = "quadratic"` adds quadratic loading parameters (λ_quad)
- Linear predictor becomes: xβ + λf + λ_quad f²
- Works with all model types: linear, probit, logit, oprobit
- Quadratic loadings are always free (estimated), never fixed
- Parameter naming: `{component}_loading_quad_{factor_index}`

### Factor Interaction Terms (Multi-Factor)

<details>
<summary><b>▶ Click to expand example</b></summary>

For multi-factor models, use `factor_spec = "interactions"` to include cross-product terms (f_j × f_k):

```r
set.seed(110)
n <- 500
f1 <- rnorm(n)  # Latent factor 1
f2 <- rnorm(n)  # Latent factor 2
x1 <- rnorm(n)

# Measurement system for 2-factor model (4 test scores, 2 per factor)
T1 <- 2.0 + 1.0*f1 + rnorm(n, 0, 0.5)  # Factor 1 indicator (loading fixed to 1)
T2 <- 1.5 + 1.2*f1 + rnorm(n, 0, 0.6)  # Factor 1 indicator
T3 <- 1.0 + 1.0*f2 + rnorm(n, 0, 0.4)  # Factor 2 indicator (loading fixed to 1)
T4 <- 0.8 + 0.9*f2 + rnorm(n, 0, 0.5)  # Factor 2 indicator

# Outcome with factor interaction effect
# Y = intercept + beta*x1 + lambda1*f1 + lambda2*f2 + lambda_inter*f1*f2 + error
Y <- 3.0 + 0.5*x1 + 0.8*f1 + 0.6*f2 + 0.4*f1*f2 + rnorm(n, 0, 0.5)

dat <- data.frame(intercept = 1, x1 = x1, T1 = T1, T2 = T2, T3 = T3, T4 = T4, Y = Y, eval = 1)

# Define 2-factor model
fm <- define_factor_model(n_factors = 2, n_types = 1)
ctrl <- define_estimation_control(n_quad_points = 8, num_cores = 1)

# Measurement equations for factor identification
mc_T1 <- define_model_component(
  name = "T1", data = dat, outcome = "T1", factor = fm,
  covariates = "intercept", model_type = "linear",
  loading_normalization = c(1.0, 0),  # Factor 1 loading=1, Factor 2 loading=0
  evaluation_indicator = "eval"
)

mc_T2 <- define_model_component(
  name = "T2", data = dat, outcome = "T2", factor = fm,
  covariates = "intercept", model_type = "linear",
  loading_normalization = c(NA_real_, 0),  # Factor 1 free, Factor 2 = 0
  evaluation_indicator = "eval"
)

mc_T3 <- define_model_component(
  name = "T3", data = dat, outcome = "T3", factor = fm,
  covariates = "intercept", model_type = "linear",
  loading_normalization = c(0, 1.0),  # Factor 1 = 0, Factor 2 loading=1
  evaluation_indicator = "eval"
)

mc_T4 <- define_model_component(
  name = "T4", data = dat, outcome = "T4", factor = fm,
  covariates = "intercept", model_type = "linear",
  loading_normalization = c(0, NA_real_),  # Factor 1 = 0, Factor 2 free
  evaluation_indicator = "eval"
)

# Outcome with INTERACTION factor effect
mc_Y <- define_model_component(
  name = "Y", data = dat, outcome = "Y", factor = fm,
  covariates = c("intercept", "x1"), model_type = "linear",
  loading_normalization = c(NA_real_, NA_real_),  # Both linear loadings free
  factor_spec = "interactions",  # Enable f1*f2 interaction term
  evaluation_indicator = "eval"
)

# Estimate
ms <- define_model_system(components = list(mc_T1, mc_T2, mc_T3, mc_T4, mc_Y), factor = fm)
result <- estimate_model_rcpp(ms, dat, init_params = NULL, control = ctrl, verbose = TRUE)

# View results - includes Y_loading_inter_1_2 parameter
print(result$estimates)
```

</details>

**Key points:**
- `factor_spec = "interactions"` adds interaction loading parameters for all factor pairs (f_j × f_k, j < k)
- For k factors, adds k(k-1)/2 interaction terms (e.g., 2 factors → 1 interaction term)
- Linear predictor becomes: xβ + Σ λ_j f_j + Σ_{j<k} λ_inter_{jk} f_j f_k
- Works with all model types: linear, probit, logit, oprobit
- Requires k ≥ 2 factors (silently ignored for single-factor models)
- Parameter naming: `{component}_loading_inter_{j}_{k}` (e.g., `Y_loading_inter_1_2`)
- Use `factor_spec = "full"` to include both quadratic and interaction terms

---

## Performance and optimization

### Choosing an optimizer

**nlminb (default, recommended)**:
- Uses analytical Hessian for fast convergence
- ~4.6x faster than nloptr on typical problems
- Supports box constraints (sigma >= 0.01, cutpoint increments >= 0.01)
- Best choice for most applications

**nloptr**:
- L-BFGS with gradient only (approximates Hessian)
- Slower but more conservative with constraints
- Use if nlminb has convergence issues

**optim (L-BFGS-B)**:
- Similar to nloptr but uses stats::optim
- Good alternative to nloptr

**trust**:
- Trust region method with Hessian
- Experimental, doesn't support fixed parameters well

### Parallelization guidelines

**When to use parallelization**:
- Large datasets (n > 5,000)
- Complex models with many components
- Multiple local optimizations (e.g., testing initial values)

**Typical speedups** (n=10,000 Roy model on laptop):
- 2 cores: 1.7x speedup (85% efficiency)
- 4 cores: 3.0x speedup (75% efficiency)

**Note**: Speedup depends on several factors:
- CPU architecture and speed
- System load (other processes running)
- Number of observations and model complexity
- Memory bandwidth

For complex models with large datasets, higher core counts can be beneficial. Running with 32 cores on a server has been effective for very complex models.

**Best practices**:
```r
# Use 1 core for small datasets or quick tests
ctrl <- define_estimation_control(num_cores = 1)

# On a shared server, be considerate of other users
n_cores <- min(8, parallel::detectCores() / 2)
ctrl <- define_estimation_control(num_cores = n_cores)

# On a dedicated machine with no other users, use all available cores
n_cores <- parallel::detectCores()
ctrl <- define_estimation_control(num_cores = n_cores)

# Enable parallel mode
result <- estimate_model_rcpp(ms, dat, control = ctrl, parallel = TRUE)
```

**Note**: Parallelization splits observations across workers. Each worker gets a subset of the data and evaluates the likelihood independently. Results are aggregated to compute the total log-likelihood and gradients.

---

## Testing

### Quick Tests

Run all automated tests (should complete in ~15 seconds):
```r
devtools::test()
```

Run a subset while developing:
```r
devtools::test(filter = "modeltypes|multifactor|oprobit")
```

### Comprehensive Tests

The systematic test suite performs extensive validation including:
- Analytical gradient verification vs. finite differences
- Analytical Hessian verification vs. finite differences
- Parameter recovery from simulated data
- Convergence from default and true initial values

Run the comprehensive suite with verbose output:
```r
# Enable verbose output and log saving
Sys.setenv(FACTORANA_TEST_VERBOSE = "TRUE")
Sys.setenv(FACTORANA_TEST_SAVE_LOGS = "TRUE")

# Run systematic tests (takes ~2-3 minutes)
devtools::test(filter = "systematic")
```

This will test:
- **Model A**: Measurement system (3 linear tests)
- **Model B**: Measurement + probit outcome
- **Model C**: Measurement + ordered probit (3 categories)
- **Model D**: Measurement + multinomial logit (3 choices)
- **Model E**: Roy selection model (wages + sector choice)

Logs are saved to `tests/testthat/test_logs/` with detailed diagnostics.

### What's Covered

- Validation of inputs & conditioning on evaluation subset
- Ordered probit thresholds (J-1, strictly increasing)
- Multi-factor loading normalization (length-k; fixed entries works)
- Parallelization correctness and performance (Roy model, n=10,000)
- Two-stage/sequential estimation (measurement system → full Roy model)
- Quadrature accuracy with properly identified three-measurement system
- Parameter recovery from simulated data
- Gradient and Hessian accuracy checks
- Factor interaction terms (quadratic and cross-product effects for all model types)

---
