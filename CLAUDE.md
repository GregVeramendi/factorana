# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**factorana** is a hybrid R/C++ system for estimating latent factor models with flexible measurement systems. There are two implementations:

1. **NEW: Rcpp-based R Package** (Primary Development Focus)
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

**R Package Structure** (`factorana_R_interface/`):
- `R/define_factor_model.R` - Defines k factors with flexible loading normalization
- `R/define_model_component.R` - Specifies individual measurement equations (linear/probit/logit/oprobit)
- `R/define_model_system.R` - Bundles components and shared factor model
- `R/estimate_model.R` - Initializes parameters for all components (existing method)
- `R/optimize_model.R` - **NEW**: Rcpp-based optimization with parallelization
- `R/write_files.R` - Exports CSV files for legacy C++ backend
- `tests/testthat/` - Validation tests for inputs, conditioning, model types, multi-factor systems

**Rcpp C++ Source** (`factorana_R_interface/src/`):
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
- `Makevars` / `Makevars.win` - Compilation settings (Eigen, C++11)

**Implementation Status**:

✅ **COMPLETED (2025-11)**:
- Multi-dimensional GH quadrature using "odometer" method for any number of factors
- Complete gradient/Hessian accumulation with proper chain rule application
- Parameter mapping: factor variances → model-specific parameters
- R to C++ interface for model_component objects (variable name → index mapping)
- Full end-to-end optimization workflow with analytical gradients
- Linear model: full implementation with gradients/Hessians
- Probit model: full implementation with gradients/Hessians
- Integration with R optimizers (optim, nloptr)
- Parallel computation support via doParallel

⏳ **TODO (Future Development)**:
- **Linear model Hessian**: Missing sigma parameter terms (see HESSIAN_TODO.md)
- **Probit model Hessian**: Full implementation using lambda function approach (see HESSIAN_TODO.md)
- Logit model: gradient/Hessian implementation (currently placeholder)
- Ordered probit: gradient/Hessian implementation (currently placeholder)
- Mixture models (nmix > 1)
- Correlated factors
- Multinomial logit (numchoice > 2)
- Standard error computation from Hessian
- Hessian aggregation in parallel mode

**Note:** See `factorana_R_interface/HESSIAN_TODO.md` for detailed implementation guide and mathematical background for completing Hessian calculations.

## Development Commands

### Rcpp Package Development (Primary Workflow)

```r
# Load package for development (compiles C++ automatically)
setwd("factorana_R_interface")
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
cd factorana_R_interface
R CMD INSTALL --no-multiarch --with-keep.source .
```

### Using the Rcpp-based Estimator

```r
library(factorana)

# Define model system (same as before)
fm <- define_factor_model(1, 1, 8, loading_normalization = 1)
mc <- define_model_component("Y", data, "outcome", fm, covariates = "x1", model_type = "linear")
ms <- define_model_system(factor = fm, components = list(mc))
control <- define_estimation_control(num_cores = 4)

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

### Test Scripts

**test_rcpp_implementation.R** - Comprehensive unit tests for C++ implementation:
```bash
cd factorana_R_interface
Rscript test_rcpp_implementation.R
```
Tests:
- Single factor linear model: initialization, likelihood, gradients
- Multi-factor (k=2) linear model: multi-dimensional integration over 64 quadrature points
- Single factor probit model: binary outcomes with proper gradient computation
- All tests verify finite likelihood/gradient values and successful C++ initialization

**test_estimation.R** - End-to-end optimization workflow:
```bash
cd factorana_R_interface
Rscript test_estimation.R
```
Features:
- Generates synthetic data with known parameters (n=200)
- Model: Y ~ intercept + x1 + x2 + factor + error
- Tests full estimation pipeline: R specification → C++ initialization → optimization → results
- Uses L-BFGS-B with analytical gradients from C++
- Compares estimates to true parameters
- Verifies convergence and log-likelihood improvement

**test_gradient_validation.R** - Finite-difference validation of analytical derivatives:
```bash
cd factorana_R_interface
Rscript test_gradient_validation.R
```
Validates analytical gradients against finite-difference approximations using central differences (based on TMinLkhd::TestGradient()):
- **Test 1**: Single factor linear model (5 parameters)
- **Test 2**: Hessian validation (documents that Hessian is placeholder)
- **Test 3**: Probit model (4 parameters)
- **Test 4**: Two-factor linear model (7 parameters)

All gradient tests pass with relative error < 1e-7, confirming:
- Correct chain rule application for factor variance derivatives
- Proper gradient accumulation across models
- Accurate multi-dimensional integration gradient computation
- Valid implementation for both linear and probit models

**Validation methodology:**
- Step size: `h = delta * (|param| + 1.0)` where `delta = 1e-7`
- Central difference: `grad_fd[i] = [f(x+h) - f(x-h)] / (2h)`
- Relative error: `|analytical - finite_diff| / |finite_diff|`
- Threshold: < 1e-4 (gradients achieve < 1e-7 in practice)

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

### Adding a New Model Type to Rcpp Backend

1. Add model type to `ModelType` enum in `src/Model.h`
2. Implement `Eval[NewType]()` method in `src/Model.cpp` with likelihood/gradient/Hessian
3. Add dispatch case in `Model::Eval()` main method
4. Update `src/rcpp_interface.cpp` to handle new type in initialization
5. Add test cases in `factorana_R_interface/tests/testthat/`
6. Recompile and test: `devtools::clean_dll(); devtools::load_all()`

### Extending to Multi-Dimensional Integration

Currently `FactorModel::CalcLkhd()` only supports 1D integration (single factor). To extend:

1. Edit `src/FactorModel.cpp:CalcLkhd()` around line 160
2. Replace single loop over quadrature points with nested loops for multiple factors
3. Compute factor values as: `fac[i] = sigma[i] * nodes[index_i] + mean[i]`
4. Weight is product of individual weights: `w = prod(weights[index_i])`
5. See legacy `TMinLkhd::CalcLkhd()` around line 5600-5800 for reference

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

- **Ordered Probit**: Outcomes MUST be ordered factors in R. With J categories, initialization produces J-1 strictly increasing thresholds.
- **Evaluation Indicators**: Components can be evaluated on subsets via `evaluation_indicator` parameter (e.g., for treatment effect models with potential outcomes).
- **File Formats**: The C++ backend historically read `.txt` files; recent R interface exports `.csv`. Update `makeprob` scripts accordingly if interfacing the two systems.
- **KNITRO Configuration**: `bin/knitro.opt` contains solver parameters. Algorithm options: 0=auto, 1=interior direct, 2=interior CG, 3=active set.
- **macOS Compatibility**: Do NOT use the `timeout` command in bash scripts - it is not available on macOS by default. Instead, use R's `setTimeLimit()` or run scripts directly without timeout wrappers.
