# Known Issues and Debugging Notes

## Issue 1: Correlated Factor Model Convergence (nlminb code 1)

**Status:** RESOLVED
**Date:** 2024-12-01
**Resolved:** 2025-12-01
**Tests affected:** `test-correlated-factors.R` (lines 133 and 272)

### Problem

The correlated two-factor model tests sometimes fail with `nlminb` returning convergence code 1 ("singular convergence") instead of 0 ("success"). This indicates the optimizer is having difficulty finding the optimum properly.

From nlminb documentation:
- Code 0: "relative convergence" (success)
- Code 1: "singular convergence" (step size became too small - indicates a problem)

### Symptoms

1. Test "correlated two-factor linear model converges and recovers correlation" (seed=42, rho=0.5) returns convergence code 1
2. Test "negative correlation is recovered correctly" (seed=456, rho=-0.4) returns convergence code 1

### Observations

- The parameter estimates appear close to the true values despite code 1
- The estimated correlation is within tolerance of the true value
- Log-likelihood improves during optimization but improvement rate slows significantly
- The issue may be related to flat likelihood surface near the optimum for correlated factor models

### Possible Causes

1. **Identification issues**: The correlation parameter may interact poorly with factor variances
2. **Cholesky parameterization**: May need reparameterization for better numerical properties
3. **Starting values**: Initial correlation of 0.0 may not be optimal
4. **Bounds**: Correlation bounded to (-0.99, 0.99) - may need adjustment
5. **Hessian issues**: The analytical Hessian for the correlation parameter may have numerical issues

### Temporary Workaround

Tests were modified to accept both codes 0 and 1 as "passing" (this is NOT a proper fix):
```r
expect_true(result$convergence %in% c(0, 1), ...)
```

### ROOT CAUSE IDENTIFIED (2025-12-01)

**The Hessian computation in `FactorModel.cpp` does NOT handle the correlated factor case!**

Investigation results:
- **Gradient check**: PASSED - All gradients match finite differences (max relative error < 1e-4)
- **Hessian check**: FAILED - Most diagonal elements have ~100% relative error

The gradient code (lines 269-312 in `FactorModel.cpp`) correctly handles the Cholesky transformation for the correlation parameter:
```cpp
// Gradient w.r.t. ρ (correlation parameter at index 2)
// df2/dρ = σ2 * (z1 - ρ*z2/sqrt_1_minus_rho2)
double df2_drho = sigma2 * (z1 - rho * z2 / sqrt_1_minus_rho2);
gradilk[2] += dL_df2 * df2_drho;
```

However, the Hessian code (lines 315-364) only handles independent factors and:
1. Does NOT compute Hessian terms for the correlation parameter (index 2)
2. Uses incorrect chain rule for the Cholesky transformation
3. Missing cross-derivatives between ρ, σ₁², and σ₂²

This explains why nlminb returns code 1 (singular convergence) - the incorrect Hessian provides bad curvature information to the optimizer, causing step size to become too small.

### RESOLUTION (2025-12-01)

The Hessian computation for correlated factor models was implemented in `FactorModel.cpp`. The fix:
1. Added proper Hessian terms for the correlation parameter (index 2)
2. Implemented cross-derivatives between ρ, σ₁², and σ₂²
3. Applied correct chain rule for the Cholesky transformation

After the fix:
- Both correlated factor tests now converge with code 0
- Gradient and Hessian pass finite-difference validation
- Test file updated to require convergence code 0 only

### TODO: Fix (COMPLETED)

1. [x] Check gradient accuracy for correlation parameter using finite differences - PASSED
2. [x] Check Hessian accuracy for correlation parameter - FAILED (bug found)
3. [x] **FIX: Implement proper Hessian for correlated factor models in `FactorModel.cpp`**
4. [x] Verify fix with finite-difference validation - PASSED
5. [x] Verify correlated factor tests pass with convergence code 0 - PASSED
6. [x] Update test file to require convergence code 0 only

### Workaround: Use Gradient-Only Optimizer

Until the Hessian bug is fixed, use a gradient-only optimizer for correlated factor models:
```r
result <- estimate_model_rcpp(..., optimizer = "nloptr")  # Uses L-BFGS (gradient only)
# OR
result <- estimate_model_rcpp(..., optimizer = "optim")   # Uses L-BFGS-B (gradient only)
```

These optimizers don't use the analytical Hessian and should converge properly.

### Files Involved

- `tests/testthat/test-correlated-factors.R` - Test definitions
- `R/optimize_model.R` - Optimizer settings (lines 598-609)
- `R/initialize_parameters.R` - Correlation parameter initialization (lines 80-85)
- `src/FactorModel.cpp` - C++ likelihood evaluation with Cholesky decomposition

### Related Code

Cholesky decomposition for correlated factors in C++:
```cpp
// f1 = sigma1 * z1
// f2 = rho * sigma2 * z1 + sigma2 * sqrt(1-rho^2) * z2
```

Parameter initialization in R:
```r
# Add correlation parameters if correlation = TRUE
if (isTRUE(model_system$factor$correlation) && n_factors == 2) {
  init_params <- c(init_params, 0.0)  # Initialize correlation to 0
  param_names <- c(param_names, "factor_corr_1_2")
}
```
