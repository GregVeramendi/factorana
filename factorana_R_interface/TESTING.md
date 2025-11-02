# Testing Organization

This document describes the testing structure for the factorana package.

## Directory Structure

```
factorana_R_interface/
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
- **test-parameter-constraints.R** - Parameter constraint system
- **test-quadrature.R** - Gauss-Hermite quadrature accuracy
- **test-components.R** - Model component functionality
- **test-01-validation.R** - Basic validation
- **test-02-conditioning.R** - Conditioning tests
- **test-03-modeltypes.R** - Different model types
- **test-04-initialize.R** - Parameter initialization
- **test-05-multifactor.R** - Multi-factor models

### Manual Tests (`tests/manual/`)
**Model Types:**
- test_linear.R
- test_probit.R
- test_mlogit.R
- test_oprobit.R

**Complex Models:**
- test_roy_model.R - Full Roy selection model
- test_probit_with_selection.R
- test_three_tests_plus_probit.R

**Validation:**
- test_gradient_validation.R - Comprehensive gradient checks
- test_rcpp_implementation.R - C++ verification
- test_match_estimate_model.R - Wrapper consistency

**Investigation:**
- test_roy_two_starts.R - Local optima
- test_with_vs_without_hessian.R - Hessian impact
- test_roy_init_at_truth.R - Optimization behavior

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
- **Automated tests**: 30 passing, 9 failing (pre-existing issues)
- **Manual tests**: 29 comprehensive tests
- **Debug scripts**: 22 diagnostic scripts
- **Validation scripts**: 4 mathematical validation scripts

## Notes

- Manual tests are preserved but not run automatically
- Debug scripts are kept for historical reference
- Validation scripts provide mathematical verification
- All development files are excluded from CRAN builds
