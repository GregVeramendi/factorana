# Optimizer Convergence Criteria

This document describes convergence codes and criteria for optimizers used in this package.

## nlminb (Default Optimizer)

`nlminb` is the default optimizer used in `estimate_model_rcpp()`. It uses gradient and Hessian information for efficient optimization.

### Convergence Codes

| Code | Status | Meaning |
|------|--------|---------|
| 0 | **SUCCESS** | Relative convergence - optimizer found an optimum |
| 1 | **FAILURE** | Singular convergence - step size became too small |
| 2 | **FAILURE** | False convergence |
| 3 | **FAILURE** | Function evaluation limit reached |
| 4 | **FAILURE** | Iteration limit reached |

**IMPORTANT:** Only convergence code 0 indicates successful convergence. Code 1 ("singular convergence") is NOT acceptable - it indicates the optimizer encountered numerical difficulties.

### Message Field

The `message` field in nlminb output provides a text description:
- `"relative convergence (4)"` - Success
- `"singular convergence (7)"` - Failure (numerical issues)
- `"false convergence (8)"` - Failure (gradient appears zero but optimum not reached)

### Evaluations Field

The `evaluations` list contains:
- `function` - Number of objective function evaluations
- `gradient` - Number of gradient evaluations

### Diagnosing Convergence Issues

When nlminb returns convergence code 1 (singular convergence):

1. **Check final gradient magnitudes**: Small gradients (< 1e-4) suggest you're near a stationary point
2. **Check Hessian condition number**: Poor conditioning (> 1e10) indicates numerical difficulties
3. **Check parameter bounds**: Parameters near bounds can cause convergence issues
4. **Try different starting values**: Poor initialization can lead to bad local optima
5. **Try different parameterization**: Transform parameters to improve numerical properties

## optim (L-BFGS-B)

### Convergence Codes

| Code | Status | Meaning |
|------|--------|---------|
| 0 | **SUCCESS** | Converged successfully |
| 1 | **FAILURE** | Iteration limit reached |
| 10 | **FAILURE** | Degeneracy in Nelder-Mead simplex |
| 51 | **WARNING** | Non-finite function value |
| 52 | **FAILURE** | Non-finite gradient or Hessian |

## nloptr (L-BFGS)

nloptr uses different status codes. See nloptr documentation for details.

## trust (Trust Region)

The `trustOptim` package returns a `converged` field (TRUE/FALSE).

## General Guidelines

### Tests Should Check for Code 0 Only

When writing tests for model convergence:

```r
# CORRECT: Only accept success
expect_equal(result$convergence, 0)

# WRONG: Don't accept failure codes as "passing"
expect_true(result$convergence %in% c(0, 1))  # BAD - code 1 is failure!
```

### Investigating Failed Convergence

1. **Print the full result object** to see all diagnostic information
2. **Check gradient at final point** - should be near zero for true optimum
3. **Check Hessian eigenvalues** - should all be positive (negative for maximization)
4. **Try alternative optimizers** - different algorithms may succeed where others fail
5. **Increase iteration/evaluation limits** if needed
6. **Improve starting values** via better initialization

### Common Causes of Convergence Failure

1. **Model identification issues** - parameters not uniquely determined by data
2. **Poor parameterization** - consider transforms (log, atanh) for bounded parameters
3. **Flat likelihood surface** - weak identification makes optimization difficult
4. **Numerical precision** - very small/large values can cause issues
5. **Bad starting values** - poor initialization leads to bad regions of parameter space
