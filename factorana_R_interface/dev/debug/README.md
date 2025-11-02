# Debug Scripts

This directory contains scripts used during development to debug specific issues.

## Purpose

These scripts helped diagnose and fix:
- C++ implementation bugs
- Data handling issues
- Gradient calculation errors
- Memory management problems
- Model component interactions

## Debug Script Categories

### Data Handling
- `debug_component_data.R` - Component data structure issues
- `debug_data_detail.R` - Detailed data inspection
- `debug_data_mismatch.R` - Data type mismatches
- `debug_force_numeric.R` - Numeric conversion issues

### C++ Interface
- `debug_cpp_consistency.R` - C++ object consistency
- `debug_fresh_ms.R` - Model system reuse bugs
- `debug_ms_modification.R` - Model system modification
- `debug_norm_modification.R` - Normalization modifications
- `debug_rcpp_conversion.R` - R/C++ conversion

### Gradient/Hessian
- `debug_gradient_indexing.R` - Gradient indexing bugs
- `debug_mlogit_gradient.R` - Multinomial logit gradients
- `debug_modhess.R` - Hessian computation
- `diagnose_sigma_sign.R` - Sigma sign issues

### Model Diagnostics
- `diagnose_roy_model.R` - Roy model diagnostics
- `test_deep_copy.R` - Deep copy behavior
- `test_fabs_derivative.R` - Absolute value derivatives
- `test_hessian_storage.R` - Hessian storage format

### Test Debugging
- `test_mlogit_debug.R` - Multinomial logit debugging
- `test_oprobit_debug.R` - Ordered probit debugging
- `test_probit_hess_debug.R` - Probit Hessian debugging
- `test_no_factor_normalization.R` - Normalization issues
- `test_sigma_hessian.R` - Sigma parameter Hessian

## Note

These scripts are preserved for:
- Historical reference
- Similar future debugging needs
- Understanding past issues
- Development documentation

They are excluded from the package build but kept in version control.
