
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

## Quick start: Roy model example

This example demonstrates a Roy selection model with unobserved ability (latent factor), sector choice, test scores, and wages.

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
fm <- define_factor_model(n_factors = 1, n_types = 1, n_quad = 16)

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
ctrl_single <- define_estimation_control(num_cores = 1)
result_single <- estimate_model_rcpp(
  model_system = ms,
  data = dat,
  init_params = NULL,       # Automatic initialization (recommended)
  control = ctrl_single,
  parallel = FALSE,
  optimizer = "nlminb",     # Default: fast with analytical Hessian
  verbose = TRUE
)

# Parallel estimation with 4 cores (3x speedup on large datasets)
ctrl_parallel <- define_estimation_control(num_cores = 4)
result_parallel <- estimate_model_rcpp(
  model_system = ms,
  data = dat,
  init_params = NULL,       # Automatic initialization
  control = ctrl_parallel,
  parallel = TRUE,
  optimizer = "nlminb",
  verbose = TRUE
)

# View results
print(result_parallel$estimates)
print(result_parallel$std_errors)
print(result_parallel$loglik)
```

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
fm <- define_factor_model(n_factors = 1, n_types = 1, n_quad = 16)

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
- With \(J\) categories, initialization returns \(J-1\) strictly increasing **thresholds**.

```r
# 4 categories -> 3 thresholds
stopifnot(is.numeric(ini$thresholds), length(ini$thresholds) == 3,
          all(diff(ini$thresholds) > 0))
```

---

## API (R)
More detailed explanations within functions.

### `define_factor_model(n_factors, n_types, n_quad_points, correlation = FALSE, n_mixtures = 1)`
- `n_factors` (int ≥0): number of latent factors (use 0 for models without factors)
- `n_types` (int ≥1): number of types
- `n_quad_points` (int ≥1): number of quadrature points used in integration
- `correlation` (logical): whether factors are correlated (default: FALSE)
- `n_mixtures` (int 1-3): number of discrete mixtures (default: 1)
- Returns an object of class `"factor_model"`
- **Note**: Loading normalization is now specified at the component level via `define_model_component()`

### `define_model_component(name, data, outcome, factor, evaluation_indicator = NULL, covariates, model_type, intercept = TRUE, num_choices = 2, nrank = NULL)`
- Validates data (no missing in eval subset), coerces outcome for `oprobit` to ordered factors if needed.
- `model_type`: `"linear"`, `"probit"`, `"logit"`, `"oprobit"`.
- Returns a `"model_component"` with pointers to `factor` (includes normalization).

### `define_model_system(components, factor)`
- Bundles components and the shared factor model into a `"model_system"`.

### `define_estimation_control(num_cores = 1)`
- Simple container for runtime controls such as the number of parallel cores.
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

## Examples

### Simple probit with factor
```r
set.seed(1)
n <- 400
x1 <- rnorm(n); x2 <- rnorm(n)
f <- rnorm(n)  # Latent factor
z <- 0.5*x1 - 0.4*x2 + 0.6*f
y <- as.integer(runif(n) < pnorm(z))
dat <- data.frame(y=y, x1=x1, x2=x2, eval=1L)

# Define factor and component
fm <- define_factor_model(n_factors = 1, n_types = 1, n_quad = 8)
mc_probit <- define_model_component(
  name = "choice", data = dat, outcome = "y", factor = fm,
  evaluation_indicator = "eval", covariates = c("x1", "x2"),
  model_type = "probit",
  loading_normalization = 1.0  # Fix for identification
)

# Initialize and estimate
ms <- define_model_system(components = list(mc_probit), factor = fm)
result <- estimate_model_rcpp(ms, dat, init_params = NULL)
print(result$estimates)
```

### Ordered probit with 5 categories
```r
set.seed(2)
n <- 350
x1 <- rnorm(n); e <- rnorm(n)
z <- 0.8*x1 + e
Y5 <- cut(z, breaks = quantile(z, probs = seq(0,1,by=0.2)),
          include.lowest = TRUE, ordered_result = TRUE)
df5 <- data.frame(Y=Y5, x1=x1, eval=1L)

# Define factor and component
fm1 <- define_factor_model(n_factors = 1, n_types = 1, n_quad = 8)
mc5 <- define_model_component(
  name = "satisfaction", data = df5, outcome = "Y", factor = fm1,
  evaluation_indicator = "eval", covariates = "x1",
  model_type = "oprobit",
  loading_normalization = 1.0  # Fix for identification
)

# Initialize (automatic in estimate_model_rcpp, but can call separately)
ini5 <- initialize_parameters(define_model_system(list(mc5), fm1), df5)
length(ini5$init_params)  # Includes covariates, loading, thresholds
```

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

**Typical speedups** (n=10,000 Roy model):
- 2 cores: 1.7x speedup (85% efficiency)
- 4 cores: 3.0x speedup (75% efficiency)
- 8 cores: Diminishing returns due to overhead

**Best practices**:
```r
# Use 1 core for small datasets
ctrl <- define_estimation_control(num_cores = 1)

# Use multiple cores for large datasets
n_cores <- min(4, parallel::detectCores() - 1)  # Leave 1 core free
ctrl <- define_estimation_control(num_cores = n_cores)

# Enable parallel mode
result <- estimate_model_rcpp(ms, dat, control = ctrl, parallel = TRUE)
```

**Note**: Parallelization splits observations across workers. Each worker gets a subset of the data and evaluates the likelihood independently. Results are aggregated to compute the total log-likelihood and gradients.

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

What's covered:
- Validation of inputs & conditioning on evaluation subset
- Ordered probit thresholds (J-1, strictly increasing)
- Multi-factor loading normalization (length-k; fixed entries works)
- Parallelization correctness and performance (Roy model, n=10,000)
- Parameter recovery from simulated data

---
