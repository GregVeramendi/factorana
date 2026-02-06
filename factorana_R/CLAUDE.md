# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

**Note:** Any working notes, debug logs, or analysis documents created during development should be named `CLAUDE_<name>.md` (e.g., `CLAUDE_debug_notes.md`, `CLAUDE_optimization_analysis.md`). These files are gitignored via the `CLAUDE*.md` pattern and won't be committed to the repository.

## Project Overview

**factorana** is a hybrid R/C++ system for estimating latent factor models with flexible measurement systems. There are two implementations:

1. **Rcpp-based R Package** (Primary Implementation)
   - R handles model specification, optimization, and parallelization
   - C++ (via Rcpp) provides fast likelihood/gradient/Hessian evaluation
   - Works on Windows, Mac, Linux without external dependencies
   - Uses doParallel for Windows-compatible parallelization

2. **Legacy: Linux C++ Backend** (For comparison/validation)
   - Standalone C++ implementation using ROOT, KNITRO, Ipopt
   - Linux-only, requires specialized libraries
   - See "Legacy C++ Backend" section below

## Architecture (Rcpp-based)

### Overview

The R package is self-contained and handles the entire estimation workflow:

1. **Model Specification** (R)
   - Define factor models with flexible loading normalization
   - Specify model components (linear/probit/logit/oprobit)
   - Bundle into model_system objects

2. **Likelihood Evaluation** (C++ via Rcpp)
   - Fast evaluation of log-likelihood, gradient, and Hessian
   - Uses Eigen for linear algebra
   - Gauss-Hermite quadrature for integration over latent factors
   - No ROOT or external optimizer dependencies

3. **Optimization** (R)
   - Uses nloptr or optim with analytical gradients
   - doParallel splits observations across workers (Windows-compatible)
   - Each worker has C++ likelihood evaluator for its data subset
   - Aggregates likelihood/gradients across workers

### Key Components

**R Package Structure** (`factorana_R/`):
- `R/define_factor_model.R` - Defines k factors with flexible loading normalization
- `R/define_model_component.R` - Specifies individual measurement equations (linear/probit/logit/oprobit)
- `R/define_model_system.R` - Bundles components and shared factor model
- `R/fix_coefficient.R` - Fix regression coefficients to specified values
- `R/estimate_model.R` - Initializes parameters for all components (existing method)
- `R/optimize_model.R` - **NEW**: Rcpp-based optimization with parallelization
- `R/write_files.R` - Exports CSV files for legacy C++ backend
- `tests/testthat/` - Validation tests for inputs, conditioning, model types, multi-factor systems

**Rcpp C++ Source** (`factorana_R/src/`):
- `Model.{h,cpp}` - Individual model likelihood evaluation (linear/probit/logit/oprobit)
  - Stripped-down version of legacy TModel
  - Only includes variables needed for likelihood calculation
  - Computes analytical gradients and Hessians
- `FactorModel.{h,cpp}` - System-level likelihood aggregation
  - Adapted from legacy TMinLkhd::CalcLkhd()
  - Integrates over latent factors using Gauss-Hermite quadrature
  - Manages multiple Model objects and parameter indexing
- `distributions.h` - Normal PDF/CDF wrappers using R's built-in functions
- `gauss_hermite.{h,cpp}` - Ported from legacy gaushermite.C (no ROOT dependencies)
- `rcpp_interface.cpp` - Rcpp exports to R:
  - `gauss_hermite_quadrature()` - Compute GH nodes/weights
  - `initialize_factor_model_cpp()` - Create FactorModel from model_system
  - `evaluate_likelihood_cpp()` - Compute likelihood/gradient/Hessian
  - `evaluate_loglik_only_cpp()` - Fast likelihood-only evaluation
  - `set_adaptive_quadrature_cpp()` - Enable adaptive integration with factor scores
  - `disable_adaptive_quadrature_cpp()` - Revert to standard quadrature
  - `set_observation_weights_cpp()` - Set per-observation likelihood weights
- `Makevars` / `Makevars.win` - Compilation settings (Eigen, C++11)

**Implementation Status**:

✅ **COMPLETED**:
- Multi-dimensional GH quadrature using "odometer" method for any number of factors
- Complete gradient/Hessian accumulation with proper chain rule application
- Parameter mapping: factor variances → model-specific parameters
- R to C++ interface for model_component objects (variable name → index mapping)
- Full end-to-end optimization workflow with analytical gradients and Hessians
- Linear model: full implementation with gradients/Hessians
- Probit model: full implementation with gradients/Hessians
- Ordered probit: full implementation with gradients/Hessians
- Logit model: full implementation with gradients/Hessians (binary and multinomial)
- Integration with R optimizers (optim, nloptr, nlminb)
- Parallel computation support via doParallel
- Standard error computation from Hessian
- Hessian aggregation in parallel mode
- Fixed coefficients: ability to fix regression coefficients to specified values via `fix_coefficient()`
- Factor interaction terms: quadratic (`f_k²`) and interaction (`f_j * f_k`) factor terms via `factor_spec` parameter
- Adaptive integration: per-observation quadrature based on factor score uncertainty for efficient two-stage estimation
- Observation weights: per-observation likelihood weighting for importance sampling
- Exploded multinomial logit: ranked choice models where individuals rank multiple alternatives

✅ **NEWLY COMPLETED**:
- Structural equation models between factors: `factor_structure = "SE_linear"` and `"SE_quadratic"`
- Equality constraints for measurement invariance via `equality_constraints` in `define_model_system()`
- Exploded nested logit: `exclude_chosen = FALSE` allows same nest to be chosen multiple times
- Rank-share corrections via `rankshare_var` parameter for rank-and-choice-specific utility adjustments
- Component-level type control: `use_types` parameter allows selective type intercepts per component
- Mixture models (nmix > 1): `n_mixtures` parameter in `define_factor_model()` for mixture of normals factor distribution

## Development Commands

### Versioning

The package uses semantic versioning (MAJOR.MINOR.PATCH). **Increment the version with each commit/push:**

- Current version is in `DESCRIPTION` file
- Increment PATCH version (e.g., 0.2.1 → 0.2.2) for each push
- Update the `Version:` field in `DESCRIPTION` before committing

Example workflow:
```bash
# Before committing, update version in DESCRIPTION
# Then commit with version bump included
git add DESCRIPTION src/*.cpp
git commit -m "Your commit message"
git push
```

### Rcpp Package Development (Primary Workflow)

```r
# Load package for development (compiles C++ automatically)
setwd("factorana_R")
devtools::load_all()

# Run all tests
devtools::test()

# Run specific test files
devtools::test(filter = "modeltypes|multifactor|oprobit")

# After modifying C++ code, recompile
devtools::clean_dll()
devtools::load_all()

# Regenerate Rcpp exports after adding/modifying exported functions
Rcpp::compileAttributes()

# Full package install from R console
install.packages(".", repos = NULL, type = "source")
```

**Build from command line:**
```bash
cd factorana_R
R CMD INSTALL --no-multiarch --with-keep.source .
```

### Using the Rcpp-based Estimator

```r
library(factorana)

# Define model system
fm <- define_factor_model(n_factors = 1, n_types = 1)
mc <- define_model_component("Y", data, "outcome", fm, covariates = "x1", model_type = "linear",
                             loading_normalization = 1)  # loading constraint at component level
ms <- define_model_system(factor = fm, components = list(mc))
control <- define_estimation_control(n_quad_points = 8, num_cores = 4)

# Estimate using Rcpp backend with parallelization
result <- estimate_model_rcpp(
  model_system = ms,
  data = data,
  control = control,
  optimizer = "nlminb",  # Uses analytical gradient + Hessian (most efficient!)
  # Other options:
  #   "nloptr" - L-BFGS (gradient only, default)
  #   "optim"  - L-BFGS-B (gradient only)
  #   "trust"  - Trust region (gradient + Hessian, requires trustOptim package)
  parallel = TRUE,
  verbose = TRUE
)

# Legacy: Export CSVs for comparison with old C++ backend
estimate_and_write(ms, fm, control, data)
```

**Optimizer Comparison:**

| Optimizer | Uses Gradient | Uses Hessian | Speed | Best For |
|-----------|--------------|--------------|-------|----------|
| `nloptr` (L-BFGS) | ✅ | ❌ | Good | Default, robust |
| `optim` (L-BFGS-B) | ✅ | ❌ | Good | Bounded parameters |
| `nlminb` | ✅ | ✅ | **Faster** | **Most efficient** |
| `trust` | ✅ | ✅ | **Fastest** | Ill-conditioned problems |

**Note:** `nlminb` and `trust` exploit analytical Hessian for 30-50% faster convergence with fewer iterations.

### Running Tests

Run all automated tests:
```r
devtools::test()  # 516 tests covering all model types and configurations
```

Run specific test subsets:
```r
devtools::test(filter = "modeltypes")    # Model type tests
devtools::test(filter = "multifactor")   # Multi-factor CFA tests
devtools::test(filter = "systematic")    # Comprehensive gradient/Hessian validation
```

The test suite validates:
- Single and multi-factor models with all model types (linear, probit, logit, oprobit)
- Analytical gradients vs finite-difference approximations (relative error < 1e-7)
- Analytical Hessians vs finite-difference approximations
- Parameter recovery from simulated data
- Parallelization correctness
- Two-stage/sequential estimation
- Fixed coefficient constraints
- Factor interaction terms (quadratic and cross-product terms)

**Important: Hessian FD tests must check EVERY element of the upper triangle, not just diagonal elements.** Use `check_hessian_accuracy()` from `helper-test-suite.R` which validates all elements.

**CRITICAL: All tests that use the optimizer MUST check convergence.** Every test that calls `estimate_model_rcpp()` must verify `expect_equal(result$convergence, 0)` (0 = success in nlminb). A test that doesn't check convergence can pass spuriously even when the model fails to converge, giving misleading results. Additionally:
- Tests with "within 2 SEs" parameter checks must first validate that SEs are reasonable (not absurdly large, e.g., < 10)
- Never accept `convergence == 1` (relative convergence) as success - only `convergence == 0` is acceptable

### Debugging Speed Issues / Benchmarking

When debugging performance issues or benchmarking code changes, run only this subset of tests for faster iteration:

```r
# 1. First run finite-difference gradient/Hessian tests (validates correctness)
devtools::test(filter = "systematic")

# 2. Then run two-stage tests (validates fixed parameter handling)
devtools::test(filter = "two-stage")
```

These tests are sufficient because:
- `systematic` tests validate gradients and Hessians against finite differences for all model types
- `two-stage` tests exercise the fixed parameter code paths (Stage 1 params fixed in Stage 2)

Only run the full test suite (`devtools::test()`) after confirming these subset tests pass.

### Modifying C++ Likelihood Code

When modifying the C++ likelihood evaluation:

1. Edit `src/Model.cpp` for model-specific likelihood/gradient/Hessian
2. Edit `src/FactorModel.cpp` for overall likelihood aggregation
3. Recompile: `devtools::clean_dll(); devtools::load_all()`
4. Test with simple example before full estimation

Key files:
- `src/Model.cpp`: Lines for Linear (~150-250), Probit (~300-400)
- `src/FactorModel.cpp`: CalcLkhd() method (~140-260)

---

## Legacy C++ Backend (Linux Only)

The original C++ backend uses ROOT framework and requires environment variables:
- `$ROOTSYS` - ROOT installation directory
- `$FACTORDIR` - This repository root (for includes)
- `$IPOPTDIR` - Ipopt optimizer installation

**Build libraries:**
```bash
cd libraries/
make -f Makefile
```

**Build problem executables:**
```bash
cd bin/
# Edit makeprob script to set correct paths for your system
./makeprob
# Produces: prob (single-threaded) and prob_mpi (MPI-parallel)
```

Note: `makeprob` scripts are system-specific and reference:
- KNITRO library paths
- ROOT libraries
- MPI compiler (mpicxx)
- Ipopt libraries

**Legacy C++ Components** (`libraries/`, `include/`):
- `TMinLkhd.{hh,cc}` - Maximum likelihood estimation engine
- `TModel.{hh,cc}` - Model specification class
- `minlkhd_nlp.{hh,cc}` - NLP interface for Ipopt/KNITRO
- `gaushermite.C` - Gauss-Hermite quadrature
- `probmain.cc` - Single-threaded entry point
- `probmain_mpi.cc` - MPI-parallel entry point

**Problem Files** (`bin/`):
- `prob.cc` - Generated problem definition file
- `makeprob*` - Scripts that generate `prob.cc`
- `knitro.opt` - KNITRO solver options

## Model Type Reference

When working with TModel or model components:
- `modtype = 1` → Linear regression
- `modtype = 2` → Probit
- `modtype = 3` → Logit
- `modtype = 4` → Ordered probit

## Factor Model Normalization

Factor loadings require normalization for identification. In R API:
- `loading_normalization = c(NA, 1, 0)` for k=3 factors:
  - `NA` → free parameter (initialized at 0.3 by default)
  - `1` → fixed at 1
  - `0` → fixed at 0

This is critical for multi-factor models to be identifiable.

## Common Workflows

### Fixing Regression Coefficients

Use `fix_coefficient()` to constrain regression coefficients to specific values during estimation:

```r
# Define a model component
mc <- define_model_component(
  name = "Y", data = dat, outcome = "Y", factor = fm,
  covariates = c("intercept", "x1", "x2"), model_type = "linear",
  loading_normalization = NA_real_, evaluation_indicator = "eval"
)

# Fix intercept to 0
mc <- fix_coefficient(mc, "intercept", 0)

# Fix x1 coefficient to a specific value
mc <- fix_coefficient(mc, "x1", 1.5)

# For multinomial logit with multiple choices, specify which choice (1-indexed, excluding reference)
mc_mlogit <- fix_coefficient(mc_mlogit, "x1", 0, choice = 1)  # Fix x1 for first non-reference choice

# Build model system and estimate as usual
ms <- define_model_system(components = list(mc), factor = fm)
result <- estimate_model_rcpp(ms, dat, control)
```

**Notes:**
- Fixed coefficients are held constant during optimization
- Only regression coefficients (betas) can be fixed, not sigma, thresholds, or loadings
- Multiple coefficients can be fixed by chaining `fix_coefficient()` calls
- The function validates that the covariate exists in the component

### Component-Level Type Usage (use_types)

When using latent types (`n_types > 1`), you can control which components have type-specific intercepts using the `use_types` parameter. This is useful for two-stage estimation where:
- Stage 1 (measurement system) doesn't need type intercepts
- Stage 2 (outcome model) uses type intercepts to capture type-specific effects

```r
# Define a model with latent types
fm <- define_factor_model(n_factors = 1, n_types = 2)

# Measurement equations WITHOUT type intercepts
mc1 <- define_model_component(
  name = "m1", data = dat, outcome = "Y1", factor = fm,
  covariates = "intercept", model_type = "linear",
  loading_normalization = 1.0,
  use_types = FALSE  # No type intercepts (default)
)
mc2 <- define_model_component(
  name = "m2", data = dat, outcome = "Y2", factor = fm,
  covariates = "intercept", model_type = "linear",
  use_types = FALSE
)

# Outcome equation WITH type intercepts
mc3 <- define_model_component(
  name = "outcome", data = dat, outcome = "Y_outcome", factor = fm,
  covariates = c("intercept", "x1"), model_type = "linear",
  use_types = TRUE  # Add type-specific intercepts
)

# Build model system and estimate
ms <- define_model_system(components = list(mc1, mc2, mc3), factor = fm)
result <- estimate_model_rcpp(ms, dat, control)
```

**Key Points:**
- `use_types = FALSE` (default): Component has same linear predictor for all types
- `use_types = TRUE`: Component gets type-specific intercepts (one per non-reference type)
- Type loadings (for computing type probabilities) are only included if at least one component uses types
- For multinomial logit with `use_types = TRUE`, each non-reference choice gets type intercepts
- Parameter naming: `{component}_type_{t}_intercept` (e.g., `outcome_type_2_intercept`)
- For multinomial logit: `{component}_c{choice}_type_{t}_intercept` (e.g., `choice_c1_type_2_intercept`)

### Using Factor Interaction Terms

The `factor_spec` parameter in `define_model_component()` enables quadratic and interaction terms for factor effects:

```r
# Define a model with factor interactions (f1 * f2 term)
mc <- define_model_component(
  name = "Y", data = dat, outcome = "Y", factor = fm,
  covariates = c("intercept", "x1"), model_type = "linear",
  loading_normalization = c(NA_real_, NA_real_),  # 2 factors
  factor_spec = "interactions",  # Add f1*f2 interaction term
  evaluation_indicator = "eval"
)
```

**Available `factor_spec` options:**
- `"linear"` (default) - Standard linear factor terms: `λ₁f₁ + λ₂f₂`
- `"quadratic"` - Add quadratic terms: `λ₁f₁ + λ₂f₂ + λ_q₁f₁² + λ_q₂f₂²`
- `"interactions"` - Add interaction terms only: `λ₁f₁ + λ₂f₂ + λ_inter*f₁*f₂`
- `"full"` - Both quadratic and interactions: `λ₁f₁ + λ₂f₂ + λ_q₁f₁² + λ_q₂f₂² + λ_inter*f₁*f₂`

**Parameter naming:**
- Quadratic loadings: `{comp}_loading_quad_{k}` (e.g., `Y_loading_quad_1`)
- Interaction loadings: `{comp}_loading_inter_{j}_{k}` (e.g., `Y_loading_inter_1_2`)
- For multinomial logit: `{comp}_c{choice}_loading_inter_{j}_{k}` (e.g., `Y_c1_loading_inter_1_2`)

**Notes:**
- Interaction terms require at least 2 factors (`n_factors >= 2`)
- Works with all model types: linear, probit, logit, and oprobit
- Factor identification via linear measurement equations is recommended when using interactions

### Exploded Multinomial Logit (Ranked Choices)

Exploded logit models are used when individuals rank multiple alternatives rather than choosing just one. For example, students ranking their preferred college majors on an application.

**Usage:**

```r
# Define model with multiple ranked outcome columns
mc <- define_model_component(
  name = "major_application",
  data = dat,
  outcome = c("major_rank1", "major_rank2", "major_rank3", "major_rank4", "major_rank5"),
  factor = fm,
  covariates = c("intercept", "ability"),
  model_type = "logit",
  n_choices = 10,  # 10 possible majors to choose from
  loading_normalization = NA_real_,
  evaluation_indicator = "eval"
)

# nrank is automatically inferred from length(outcome) = 5
```

**Key Points:**
- Specify `outcome` as a character vector of column names (one per rank)
- `n_choices` is the total number of alternatives (e.g., majors)
- `nrank` is inferred from `length(outcome)`, or can be set explicitly
- Missing ranks are allowed: values outside 1..n_choices in any rank column are skipped
- Likelihood is the product of per-rank likelihoods: L = ∏ P(choice_r | X, f)
- Same regression coefficients and factor loadings apply to all ranks

**Data Format:**
- Each rank column contains integers 1..n_choices indicating which alternative was chosen at that rank
- Missing/unused ranks can be indicated by 0, NA, or any value outside 1..n_choices
- Example: If a student only ranked 3 of 5 ranks, the last 2 columns would have NA or 0

### Exploded Nested Logit (Same Nest Can Be Chosen Multiple Times)

For nested logit models where choices represent **nests** (groups of similar items) rather than individual alternatives, the same nest can be chosen multiple times across ranks. Use `exclude_chosen = FALSE` to enable this behavior.

**Usage:**

```r
# Exploded nested logit - same nest can be chosen again
mc <- define_model_component(
  name = "nest_choice",
  data = dat,
  outcome = c("rank1", "rank2", "rank3"),
  factor = fm,
  covariates = c("intercept", "x1"),
  model_type = "logit",
  num_choices = 5,  # 5 nests
  exclude_chosen = FALSE,  # Allow same nest to be chosen again
  loading_normalization = NA_real_,
  evaluation_indicator = "eval"
)
```

**Key Differences from Standard Exploded Logit:**
- `exclude_chosen = TRUE` (default): Standard exploded logit - each alternative can only be chosen once
- `exclude_chosen = FALSE`: Nested logit - alternatives remain available for subsequent ranks

### Rank-Share Corrections (rankshare_var)

The `rankshare_var` parameter provides rank-and-choice-specific corrections to the linear predictor. This is useful when an attractive item within a nest affects the probability of choosing that nest again.

**Usage:**

```r
# Create rankshare correction columns in data
# Layout: (num_choices-1) * nrank columns
# Access pattern: rankshare_var + (num_choices-1)*irank + icat
# where irank = 0..nrank-1 and icat = 0..num_choices-2

# Example for 4 choices, 3 ranks: need (4-1)*3 = 9 columns
dat$rs_r0_c0 <- correction_values_rank0_choice1
dat$rs_r0_c1 <- correction_values_rank0_choice2
dat$rs_r0_c2 <- correction_values_rank0_choice3
dat$rs_r1_c0 <- correction_values_rank1_choice1
# ... and so on

mc <- define_model_component(
  name = "nested_choice",
  data = dat,
  outcome = c("rank1", "rank2", "rank3"),
  factor = fm,
  covariates = c("intercept", "x1"),
  model_type = "logit",
  num_choices = 4,
  exclude_chosen = FALSE,
  rankshare_var = "rs_r0_c0",  # First rankshare column
  loading_normalization = NA_real_,
  evaluation_indicator = "eval"
)
```

**Notes:**
- Corrections are added to the linear predictor: `V_k = X*beta + lambda*f + rankshare_correction_k`
- The reference category (choice 0) has no correction
- Missing corrections (values <= -9998) are treated as 0
- Rankshare corrections don't add parameters - they're data values, so gradients w.r.t. model parameters are unaffected

### Adaptive Integration (Two-Stage Estimation)

Adaptive integration reduces computational cost in second-stage estimation by using factor scores and their standard errors from Stage 1 to determine how many integration points each observation needs:

- **Small SE** (factor well-identified): Use 1 point centered at factor score
- **Medium SE**: Use fewer points centered at factor score with importance sampling
- **Large SE** (factor poorly identified): Use full quadrature centered at prior mean (0)

**Formula:** `n_quad_obs = 1 + 2 × floor(factor_se / factor_sd / threshold)`

where `factor_sd = sqrt(factor_variance)` is the factor standard deviation.

**Typical Two-Stage Workflow:**

```r
library(factorana)

# === STAGE 1: Estimate measurement model ===
# Define Stage 1 model (measurement equations only)
fm <- define_factor_model(n_factors = 2)
mc1 <- define_model_component("m1", data, "Y1", fm, model_type = "linear",
                               loading_normalization = c(1, 0))
mc2 <- define_model_component("m2", data, "Y2", fm, model_type = "linear",
                               loading_normalization = c(NA, 1))
ms_stage1 <- define_model_system(factor = fm, components = list(mc1, mc2))
control <- define_estimation_control(n_quad_points = 16)

# Estimate Stage 1
result_s1 <- estimate_model_rcpp(ms_stage1, data, control)

# Extract factor scores and SEs
fscores <- estimate_factorscores_rcpp(ms_stage1, data, result_s1$par, control)
factor_scores <- fscores$factor_scores   # matrix [nobs x nfac]
factor_ses <- fscores$factor_ses         # matrix [nobs x nfac]

# === STAGE 2: Estimate structural model with adaptive integration ===
# Add structural component
mc3 <- define_model_component("outcome", data, "Y_outcome", fm,
                               covariates = c("intercept", "x1"),
                               model_type = "linear",
                               loading_normalization = c(NA, NA))
ms_stage2 <- define_model_system(factor = fm, components = list(mc1, mc2, mc3))

# Initialize FactorModel for Stage 2
fm_ptr <- initialize_factor_model_cpp(ms_stage2, data, n_quad = 16)

# Extract factor variances from Stage 1 results
factor_vars <- c(result_s1$par["factor_var_1"], result_s1$par["factor_var_2"])

# Enable adaptive integration (prints diagnostic summary)
set_adaptive_quadrature_cpp(
  fm_ptr,
  factor_scores,   # From Stage 1
  factor_ses,      # From Stage 1
  factor_vars,     # Factor variances from Stage 1
  threshold = 0.3, # Smaller = more integration points
  max_quad = 16,   # Maximum points per factor
  verbose = TRUE   # Print summary table
)

# Output example:
# Adaptive Integration Summary
# ----------------------------
# Threshold: 0.3, Max quad points: 16
#
# Integration points per observation:
#   Points   Observations   Percent
#        1            400      80.0%
#        3            100      20.0%
#
# Average integration points: 1.4 (vs 16 standard)
# Computational reduction: 91.2%

# Now run estimation with the adaptive FactorModel
# ... (use evaluate_likelihood_cpp with fm_ptr)
```

**Control Parameters:**

In `define_estimation_control()`:
- `adaptive_integration = FALSE` (default): Standard quadrature
- `adaptive_integration = TRUE`: Enable adaptive mode
- `adapt_int_thresh = 0.3` (default): Threshold for determining points

**Low-Level Functions:**

| Function | Description |
|----------|-------------|
| `set_adaptive_quadrature_cpp(fm_ptr, scores, ses, vars, threshold, max_quad, verbose)` | Enable adaptive integration |
| `disable_adaptive_quadrature_cpp(fm_ptr)` | Revert to standard quadrature |
| `set_observation_weights_cpp(fm_ptr, weights)` | Set per-observation likelihood weights |

**Notes:**
- Importance sampling correction is applied automatically to account for non-centered integration
- The IS formula is: `IS = (SE/σ) × exp(z²/2 - f²/(2σ²))` where:
  - `SE` is the factor score standard error from Stage 1
  - `σ` is the factor standard deviation (sqrt of variance parameter)
  - `z` is the Gauss-Hermite quadrature node
  - `f = factor_score + SE × z` is the factor value at this integration point
- When factor SE is very large (> sqrt(variance)), full quadrature is used at prior center
- **Expected approximation error**: Adaptive quadrature with mostly 1-point integration (well-identified factors) typically gives ~15-20% likelihood discrepancy compared to full quadrature. This is inherent to the mode approximation. Parameter estimates remain close (typically <5% relative difference).
- Disabling adaptive mode clears both adaptive settings and observation weights

### Observation Weights

Observation weights allow different observations to have different influence on the likelihood. This is useful for:
- Survey weights (to make estimates representative of a population)
- Importance sampling (as used internally by adaptive integration)
- Down-weighting outliers or problematic observations

**Usage:**

```r
# Add weights variable to your data
data$my_weights <- some_weight_vector

# Specify weights in define_model_system()
ms <- define_model_system(
  components = list(mc1, mc2),
  factor = fm,
  weights = "my_weights"  # Name of weights column in data
)

# Estimate as usual - weights are applied automatically
result <- estimate_model_rcpp(ms, data, control)
```

**Requirements:**
- Weights must be a numeric column in the data passed to `estimate_model_rcpp()`
- Weights must be positive (warning issued if any are <= 0)
- Weights cannot contain NA values (will error)
- Each observation's log-likelihood contribution is multiplied by its weight

### Checkpointing (Long-Running Estimations)

For long-running estimations that may take hours or days, the `checkpoint_file` parameter saves parameters periodically so you can restart from the last checkpoint if needed.

**Usage:**

```r
# Run estimation with checkpointing enabled
result <- estimate_model_rcpp(
  model_system = ms,
  data = dat,
  control = control,
  checkpoint_file = "my_checkpoint.txt"  # Parameters saved here
)
```

**How it works:**
- Parameters are saved each time the Hessian is evaluated at a point with improved likelihood
- This "smart" approach only saves when the likelihood actually improves, not every iteration
- The checkpoint file contains:
  - Header with timestamp, log-likelihood, and iteration number
  - Parameter names and values in CSV format

**Example checkpoint file:**
```
# Checkpoint saved: 2026-01-22 16:18:33
# Log-likelihood: -502.64076359
# Iteration: 8

"parameter","value"
"factor_var_1",0.718778794640373
"m1_intercept",-0.00617248432231106
"m1_sigma",0.607418329076786
...
```

**Restarting from a checkpoint:**

```r
# Read checkpoint file
checkpoint <- read.csv("my_checkpoint.txt", comment.char = "#")
init_params <- setNames(checkpoint$value, checkpoint$parameter)

# Resume estimation from checkpoint
result <- estimate_model_rcpp(
  model_system = ms,
  data = dat,
  init_params = init_params,  # Start from checkpoint parameters
  control = control,
  checkpoint_file = "my_checkpoint.txt"  # Continue saving checkpoints
)
```

**Notes:**
- The checkpoint file is overwritten each time a better point is found
- If estimation is interrupted, the file contains the last best parameters
- Reading the checkpoint with `comment.char = "#"` skips the metadata header lines

### Equality Constraints (Measurement Invariance)

Equality constraints allow constraining parameters to be equal during estimation. This is essential for longitudinal models where measurement invariance is assumed (e.g., same factor loadings for measures at different ages).

**Usage:**

```r
# Define model with two time periods using the same measures
fm <- define_factor_model(n_factors = 2, factor_structure = "SE_linear")

# Measurement equations for age 10 (factor 1)
mc_t1_age10 <- define_model_component("t1_age10", data, "test1_age10", fm,
                                       model_type = "linear", loading_normalization = c(1, 0))
mc_t2_age10 <- define_model_component("t2_age10", data, "test2_age10", fm,
                                       model_type = "linear", loading_normalization = c(NA, 0))

# Measurement equations for age 15 (factor 2)
mc_t1_age15 <- define_model_component("t1_age15", data, "test1_age15", fm,
                                       model_type = "linear", loading_normalization = c(0, 1))
mc_t2_age15 <- define_model_component("t2_age15", data, "test2_age15", fm,
                                       model_type = "linear", loading_normalization = c(0, NA))

# Specify equality constraints: same loadings and sigmas across time
ms <- define_model_system(
  components = list(mc_t1_age10, mc_t2_age10, mc_t1_age15, mc_t2_age15),
  factor = fm,
  equality_constraints = list(
    c("t2_age10_loading_1", "t2_age15_loading_2"),  # Equal loadings
    c("t1_age10_sigma", "t1_age15_sigma"),          # Equal sigmas for test 1
    c("t2_age10_sigma", "t2_age15_sigma")           # Equal sigmas for test 2
  )
)

result <- estimate_model_rcpp(ms, data, control)

# Tied parameters will be exactly equal in result$estimates
```

**Notes:**
- First parameter in each group is the "primary" (freely estimated)
- All other parameters in the group are set equal to the primary
- Gradients are correctly aggregated for tied parameters
- Parameter table in verbose output shows `tied_to` column for constrained parameters

### Structural Equation Models (SE_linear, SE_quadratic)

Structural equation models allow specifying causal relationships between factors, replacing the correlation-based approach for multi-factor models.

**Key concept:** Instead of integrating over all factors directly, we integrate over *primitive* random variables: input factors (f₁) and residuals (ε). The outcome factor is deterministic given these:
- `SE_linear`: f₂ = α + α₁*f₁ + ε
- `SE_quadratic`: f₂ = α + α₁*f₁ + α₂*f₁² + ε

**Usage:**

```r
# Define 2-factor SE model
fm <- define_factor_model(
  n_factors = 2,
  factor_structure = "SE_quadratic"  # f₂ = α + α₁*f₁ + α₂*f₁² + ε
)

# Measurement equations for factor 1 (input factor)
mc1 <- define_model_component("Y1", data, "Y1", fm, model_type = "linear",
                               loading_normalization = c(1, 0))
mc2 <- define_model_component("Y2", data, "Y2", fm, model_type = "linear",
                               loading_normalization = c(NA, 0))

# Measurement equations for factor 2 (outcome factor)
mc3 <- define_model_component("Y3", data, "Y3", fm, model_type = "linear",
                               loading_normalization = c(0, 1))
mc4 <- define_model_component("Y4", data, "Y4", fm, model_type = "linear",
                               loading_normalization = c(0, NA))

ms <- define_model_system(components = list(mc1, mc2, mc3, mc4), factor = fm)
result <- estimate_model_rcpp(ms, data, control)

# SE parameters in result$estimates:
# - se_intercept: α (intercept in structural equation)
# - se_linear_1: α₁ (coefficient of f₁)
# - se_quadratic_1: α₂ (coefficient of f₁², for SE_quadratic only)
# - se_residual_var: σ²_ε (residual variance)
```

**Parameter names:**
- `factor_var_1`: Variance of input factor f₁
- `se_intercept`: Intercept α in structural equation
- `se_linear_k`: Coefficient of f_k in structural equation (for each input factor)
- `se_quadratic_k`: Coefficient of f_k² (SE_quadratic only)
- `se_residual_var`: Variance of residual ε

**Identification:**
- Input factor variance (f₁) is identified via measurement equations
- Residual variance (σ²_ε) is a free parameter in the SE specification
- Outcome factor variance is derived: Var(f₂) = α₁²·Var(f₁) + σ²_ε

### Adding a New Model Type to Rcpp Backend

1. Add model type to `ModelType` enum in `src/Model.h`
2. Implement `Eval[NewType]()` method in `src/Model.cpp` with likelihood/gradient/Hessian
3. Add dispatch case in `Model::Eval()` main method
4. Update `src/rcpp_interface.cpp` to handle new type in initialization
5. Add test cases in `factorana_R/tests/testthat/`
6. Recompile and test: `devtools::clean_dll(); devtools::load_all()`

### Multi-Dimensional Integration

`FactorModel::CalcLkhd()` supports multi-factor models using the "odometer" method:
- Uses `nquad^nfactors` total integration points (e.g., 8 quadrature points × 2 factors = 64 points)
- Factor indices cycle like an odometer: `facint[0]` increments first, then `facint[1]`, etc.
- Each integration point computes: `fac[i] = sigma[i] * node[facint[i]] + mean[i]`
- Weight is product of individual weights: `w = prod(weights[facint[i]])`
- See `src/FactorModel.cpp:CalcLkhd()` lines 191-327 for implementation

### Comparing Rcpp vs Legacy Results

1. Run Rcpp estimation: `result_rcpp <- estimate_model_rcpp(ms, data, control)`
2. Export CSVs for legacy: `estimate_and_write(ms, fm, control, data)`
3. Copy CSVs to `bin/`, compile and run legacy C++ backend
4. Compare likelihoods, parameters, gradients between implementations

### Performance Profiling

```r
# Profile R code
profvis::profvis({
  result <- estimate_model_rcpp(ms, data, control)
})

# Benchmark parallelization
library(microbenchmark)
microbenchmark(
  sequential = estimate_model_rcpp(ms, data, control, parallel = FALSE),
  parallel_2 = estimate_model_rcpp(ms, data, control_2cores, parallel = TRUE),
  parallel_4 = estimate_model_rcpp(ms, data, control_4cores, parallel = TRUE),
  times = 5
)
```

## Important Details

- **Optimizer Choice**: Always use `optimizer = "nlminb"` for estimation. This is the only optimizer that uses the analytical Hessian, which provides faster convergence and more accurate standard errors. Other optimizers (`optim`, `nloptr`) only use gradients and are slower. **DO NOT switch to a different optimizer without explicit user permission.**
- **Convergence Settings**: All optimization routines (model estimation, factor score estimation, etc.) use R default convergence settings. These are well-scaled, smooth optimization problems. **DO NOT CHANGE convergence settings (tolerances, iteration limits, etc.) without explicit user permission.**
- **Test Criteria**: When a test fails, fix the underlying code, not the test. **DO NOT relax test tolerances, thresholds, or acceptance criteria without explicit user permission.** Changing test criteria to make a failing test pass is not acceptable.
- **Ordered Probit**: Outcomes MUST be ordered factors in R. With J categories, initialization produces J-1 strictly increasing thresholds.
- **Evaluation Indicators**: Components can be evaluated on subsets via `evaluation_indicator` parameter (e.g., for treatment effect models with potential outcomes).
- **File Formats**: The C++ backend historically read `.txt` files; recent R interface exports `.csv`. Update `makeprob` scripts accordingly if interfacing the two systems.
- **KNITRO Configuration**: `bin/knitro.opt` contains solver parameters. Algorithm options: 0=auto, 1=interior direct, 2=interior CG, 3=active set.
- **macOS Compatibility**: Do NOT use the `timeout` command in bash scripts - it is not available on macOS by default. Instead, use R's `setTimeLimit()` or run scripts directly without timeout wrappers.
- **Legacy Code Naming**: In the legacy C++ code (`TMinLkhd.cc`), `Getfvar(imix, ifac)` returns the factor **standard deviation** (σ), NOT the variance (σ²). This is confusing because "fvar" suggests variance. The Rcpp code uses `sigma_fac` for SD and `factor_var` for variance to be clearer.
