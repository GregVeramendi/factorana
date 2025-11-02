# Validation Scripts

This directory contains scripts for validating the correctness of the implementation.

## Purpose

These scripts verify:
- Gradient accuracy via finite differences
- Hessian accuracy via finite differences
- Quadrature integration accuracy
- Parameter bounds handling
- Optimization convergence

## Validation Scripts

### Gradient Validation
- `check_gradient_at_true.R` - Check gradients at true parameter values
- `check_gradient_quad.R` - Gradient accuracy with different quadrature points

### Parameter Validation
- `check_parameter_bounds.R` - Verify parameter bound handling

### Formula Verification
- `verify_hessian_formula.R` - Verify Hessian mathematical formulas

## Running Validations

Run individual validation:
```r
source("dev/validation/check_gradient_at_true.R")
```

## Validation Criteria

- **Gradients**: Relative difference < 1e-4 vs finite differences
- **Hessians**: Relative difference < 1e-10 vs finite differences
- **Quadrature**: Results stable across different n_quad values

## Note

These validation scripts complement the automated tests in `tests/testthat/` by providing:
- More detailed diagnostic output
- Validation at specific parameter values
- Comprehensive checking of mathematical correctness
