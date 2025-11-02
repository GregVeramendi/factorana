# Manual Tests

This directory contains comprehensive manual test scripts that are not part of the automated `devtools::test()` suite.

## Purpose

These tests are more exploratory and comprehensive than unit tests. They:
- Validate full model estimation workflows
- Test parameter recovery with simulated data
- Investigate optimization behavior
- Compare different estimation strategies

## Running Manual Tests

Individual tests can be run directly:
```r
source("tests/manual/test_roy_model.R")
```

Or run all manual tests:
```r
source("tests/manual/run_all_tests.R")
```

## Test Categories

### Model Type Tests
- `test_linear.R` - Linear regression components
- `test_probit.R` - Probit model components
- `test_mlogit.R` - Multinomial logit components
- `test_oprobit.R` - Ordered probit components

### Complex Models
- `test_roy_model.R` - Full Roy selection model with missing data
- `test_probit_with_selection.R` - Probit with selection
- `test_three_tests_plus_probit.R` - Multiple components

### Validation Tests
- `test_gradient_validation.R` - Gradient validation via finite differences
- `test_rcpp_implementation.R` - C++ implementation verification
- `test_match_estimate_model.R` - Wrapper vs direct optimization comparison

### Investigation Tests
- `test_roy_two_starts.R` - Multiple local optima investigation
- `test_with_vs_without_hessian.R` - Hessian impact on convergence
- `test_roy_init_at_truth.R` - Optimization from true parameters

### Specialized Tests
- `test_loading_recovery.R` - Factor loading parameter recovery
- `test_loading_only_fixed.R` - Fixed loading specifications
- `test_component_normalizations.R` - Normalization strategies
- `test_initialization.R` - Parameter initialization strategies

## Note

These tests are intentionally excluded from automated testing because they:
- Take longer to run
- Are exploratory in nature
- Test edge cases and unusual scenarios
- Generate detailed diagnostic output
