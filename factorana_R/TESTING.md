# Testing Organization

This document describes the testing structure for the factorana package.

## Directory Structure

```
factorana_R/
├── tests/
│   ├── testthat/           # Automated tests (run via devtools::test())
│   │   ├── test-*.R        # Test files
│   │   └── helper-*.R      # Helper functions for tests
│   ├── testthat.R          # Test runner configuration
│   └── manual/             # Manual/exploratory tests (not automated)
│       ├── README.md       # Documentation
│       ├── run_all_tests.R # Script to run all manual tests
│       └── test_*.R        # Manual test files
├── dev/
│   ├── debug/              # Debug scripts from development
│   │   ├── README.md
│   │   └── debug_*.R       # Debug scripts
│   └── validation/         # Validation scripts
│       ├── README.md
│       └── check_*.R       # Validation scripts
```

## Running Tests

### Automated Tests
Run the automated test suite using devtools:
```r
devtools::test()
```

Or from command line:
```bash
R CMD check .
```

### Manual Tests
Run all manual tests:
```r
source("tests/manual/run_all_tests.R")
```

Run individual manual test:
```r
source("tests/manual/test_roy_model.R")
```

### Validation Scripts
Run validation checks:
```r
source("dev/validation/check_gradient_at_true.R")
```

## Test Categories

### Automated Tests (`tests/testthat/`)
- **test-components.R** - Model component functionality
- **test-correlated-factors.R** - Correlated two-factor models
- **test-factorscores.R** - Factor score estimation (MAP estimates)
- **test-fixed-coefficients.R** - Fixed coefficient constraints
- **test-parallelization.R** - Parallel estimation correctness and performance
- **test-parameter-constraints.R** - Parameter constraint system
- **test-quadrature.R** - Gauss-Hermite quadrature accuracy with three-measurement system
- **test-systematic-suite.R** - Systematic model validation
- **test-two-stage.R** - Multi-stage/sequential estimation with Roy model
- **test-types.R** - Latent types models
- **test-test-01-validation.R** - Basic validation
- **test-test-02-conditioning.R** - Conditioning tests
- **test-test-03-modeltypes.R** - Different model types (linear, probit, ordered probit)
- **test-test-04-initialize.R** - Parameter initialization
- **test-test-05-multifactor.R** - Multi-factor models with component-level normalization

### Manual Tests (`tests/manual/`)
**Model Types:**
- test_linear.R - Linear regression components
- test_probit.R - Probit model components
- test_mlogit.R - Multinomial logit components
- test_oprobit.R - Ordered probit components

**Complex Models:**
- test_roy_model.R - Full Roy selection model with missing data
- test_probit_with_selection.R - Probit with selection
- test_three_tests_plus_probit.R - Multiple components
- test_simple_three_tests.R - Three test scores model
- test_three_tests_all_free.R - All free loadings
- test_three_tests_from_truth.R - Optimization from true values

**Validation:**
- test_gradient_validation.R - Comprehensive gradient checks
- test_rcpp_implementation.R - C++ implementation verification
- test_match_estimate_model.R - Wrapper consistency
- test_probit_hessian.R - Probit Hessian validation

**Investigation:**
- test_roy_two_starts.R - Multiple local optima investigation
- test_with_vs_without_hessian.R - Hessian impact on convergence
- test_roy_init_at_truth.R - Optimization from true parameters
- test_check_optimum.R - Optimum verification
- test_simple_selection_from_truth.R - Selection model optimization

**Specialized Tests:**
- test_loading_recovery.R - Factor loading parameter recovery
- test_loading_only_fixed.R - Fixed loading specifications
- test_component_normalizations.R - Component-specific normalization strategies
- test_initialization.R - Parameter initialization strategies
- test_parameter_constraints.R - Parameter constraint handling
- test_estimation.R - General estimation workflow
- test_gh_simple.R - Gauss-Hermite integration

### Debug Scripts (`dev/debug/`)
Scripts used during development to diagnose issues:
- Data handling: debug_data_*.R
- C++ interface: debug_cpp_*.R, debug_fresh_ms.R
- Gradients: debug_gradient_*.R
- Diagnostics: diagnose_*.R

### Validation Scripts (`dev/validation/`)
Mathematical correctness validation:
- check_gradient_at_true.R
- check_gradient_quad.R
- check_parameter_bounds.R
- verify_hessian_formula.R

## Build Configuration

The `.Rbuildignore` excludes development files from package builds:
- `^dev$` - All development scripts
- `^tests/manual$` - Manual tests
- Documentation markdown files

## Adding New Tests

### For Automated Testing
1. Create `tests/testthat/test-feature-name.R`
2. Use `test_that()` blocks with `expect_*()` assertions
3. Mark slow tests with `skip_on_cran()`

Example:
```r
test_that("feature works correctly", {
  skip_on_cran()

  # Setup
  result <- my_function(params)

  # Assertions
  expect_true(!is.null(result))
  expect_equal(result$value, expected_value)
})
```

### For Manual Testing
1. Create `tests/manual/test_feature.R`
2. Use plain R script with informative output
3. Document in `tests/manual/README.md`

### For Debugging
1. Create `dev/debug/debug_issue_name.R`
2. Include detailed diagnostics
3. Document what issue it helped solve

## Test Results

Current test suite status:
- **Automated tests**: 17 test files covering core functionality
  - **254 tests passing**, 0 failures, 0 warnings (as of last run)
  - Includes parallelization test with Roy model (n=10,000)
  - Includes two-stage estimation test with Roy model
  - Quadrature test uses three-measurement system for proper identification
  - Fixed coefficient constraints tests
  - Factor score estimation tests
  - Correlated factor model tests
  - Latent types model tests
- **Manual tests**: 26 comprehensive tests
  - Model types: linear, probit, logit, oprobit
  - Complex models: Roy selection, multiple components
  - Validation: gradients, Hessians, parameter recovery
- **Debug scripts**: Historical diagnostic scripts in `dev/debug/`
- **Validation scripts**: 4 mathematical validation scripts in `dev/validation/`

## Notes

- Manual tests are preserved but not run automatically
- Debug scripts are kept for historical reference
- Validation scripts provide mathematical verification
- All development files are excluded from CRAN builds
