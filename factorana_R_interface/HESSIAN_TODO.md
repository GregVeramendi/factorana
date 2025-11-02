# Hessian Implementation Status and TODO

## Current Status (as of 2025-11)

### ✅ What's Already Implemented

**R Optimization Infrastructure (optimize_model.R)**:
- ✅ Hessian function that converts C++ upper triangle to full matrix
- ✅ Support for `nlminb` optimizer (uses analytical gradient + Hessian)
- ✅ Support for `trust` optimizer (trust region with analytical Hessian)
- ✅ Parallel Hessian aggregation (TODO: verify correctness)
- ⚠️ L-BFGS methods available but DON'T use Hessian (gradient only)

### ✅ What's Already Implemented in C++

**Linear Model (Model.cpp:143-235)**:
- Hessian matrix initialization: `H[i,j] = -1/(σ²)` for all upper triangle elements
- Factor variance parameter Hessian terms: multiply rows/columns by loading parameters
- Factor loading parameter Hessian terms: multiply rows/columns by factor values
- Cross-derivative term: `H[θᵢ, αᵢ] += Z/σ²` (factor variance × loading)
- Regression coefficient Hessian terms: multiply rows/columns by regressor data values

**Probit Model (Model.cpp:237+)**:
- Basic structure exists but needs completion

## ⏳ What Remains to Be Done

### 1. **Linear Model - Complete Sigma Parameter Hessian** (PRIORITY: HIGH)

Currently missing the Hessian terms for the error variance parameter (σ).

**From legacy TModel.cc (lines ~2100-2110):**
```cpp
// sigma_i index
Int_t index = numfac+nregressors+ifreefac;
// d/dsigma_i row
for (int j = index; j < npar ; j++)
    hess[index*npar+j] *= 2.0*Z/sigma;
// d/dsigma_i col
for (int j = 0; j < index ; j++)
    hess[j*npar+index] *= 2.0*Z/sigma;
// separately handle diagonal element
hess[index*npar+index] = -3.0*Z*Z/pow(sigma,4) + 1.0/(sigma*sigma);
```

**What this does:**
- Off-diagonal elements of sigma's row/column: multiply by `2Z/σ`
- Diagonal element: `d²L/dσ² = -3Z²/σ⁴ + 1/σ²`

**Implementation location:** Add after line 234 in Model.cpp:EvalLinear()

---

### 2. **Probit Model - Complete Hessian Implementation** (PRIORITY: MEDIUM)

The probit Hessian uses a different approach based on the "lambda function": λ(z) = φ(z)/Φ(z)

**From legacy TModel.cc (lines ~2200-2350):**

**Step 1:** Initialize Hessian to 1.0 (not -1/σ² like linear):
```cpp
hess.resize(npar*npar,0.0);
for (int i = 0; i < npar ; i++) {
    for (int j = i ; j < npar ; j++) {
        hess[i*npar+j] = 1.0;
    }
}
```

**Step 2:** For factor variance parameters (θ):
```cpp
// lambda^L(theta) - lambda^Prob(theta) (row)
for (int j = i; j < npar ; j++)
    hess[i*npar+j] *= -Z*obsSign*param[...] - modEval[i+1];
// dZ/dtheta_i (col)
for (int j = 0; j <= i ; j++)
    hess[j*npar+i] *= obsSign*param[...];
```

Where:
- `Z = expres[0]` (linear predictor)
- `obsSign = 1.0` if Y=1, `-1.0` if Y=0
- `lambda^Prob(θ) = modEval[i+1]` (gradient value)

**Step 3:** For factor loadings (α):
```cpp
Int_t index = numfac+nregressors+ifreefac;
// lambda^L(alpha) - lambda^Prob(alpha) (row)
for (int j = index; j < npar ; j++)
    hess[index*npar+j] *= -Z*obsSign*fac[i] - modEval[...];
// dZ/dalpha (col)
for (int j = 0; j <= index ; j++)
    hess[j*npar+index] *= obsSign*fac[i];
```

**Step 4:** For regression coefficients (β):
```cpp
Int_t index = numfac+ireg;
// lambda^L(X) - lambda^Prob(X) (row)
for (int j = index; j < npar ; j++)
    hess[index*npar+j] *= -Z*obsSign*data[...] - modEval[...];
// dZ/dX_i (col)
for (int j = 0; j <= index ; j++)
    hess[j*npar+index] *= obsSign*data[...];
```

**Step 5:** Add cross-derivative terms dZ/dθᵢ dαⱼ:
```cpp
ifreefac = 0;
for (int i = 0 ; i < numfac; i++) {
    if (facnorm.size()==0 || facnorm[i]<-9998.0) {
        Int_t index = numfac+nregressors+ifreefac;
        hess[i*npar+index] += obsSign;
        ifreefac++;
    }
}
```

**Step 6:** Final scaling - multiply all by pdf/cdf:
```cpp
for (int i = 0; i < npar ; i++) {
    for (int j = i ; j < npar ; j++) {
        hess[i*npar+j] *= pdf/cdf;
    }
}
```

**Key mathematical insight:**
The probit Hessian uses the identity:
```
d²log(Φ(z))/dxᵢdxⱼ = -λ(z)[z(dz/dxᵢ)(dz/dxⱼ) + d²z/dxᵢdxⱼ] - λ(z)²(dz/dxᵢ)(dz/dxⱼ)
                     = (pdf/cdf)[-(z·dz/dxᵢ + λ)(dz/dxⱼ) + d²z/dxᵢdxⱼ]
```

**Implementation location:** Replace placeholder in Model.cpp:EvalProbit()

---

### 3. **Ordered Probit (Oprobit) - Full Implementation** (PRIORITY: LOW)

Currently has placeholder code. Needs complete gradient and Hessian.

**Complexity:** High - involves multiple threshold parameters and cumulative differences.

**From legacy TModel.cc (lines ~2500-2800):**
- Multiple PDF/CDF evaluations at different thresholds
- Gradient involves ratios of PDFs to CDF differences
- Hessian has complex cross terms between thresholds and other parameters

**Defer until:** Linear and probit Hessians are validated.

---

### 4. **Multinomial Logit - Full Implementation** (PRIORITY: LOW)

Currently has placeholder code.

**From legacy TModel.cc (lines ~2900-3200):**
- Uses softmax probabilities
- Hessian involves Jacobian of softmax transformation
- Very complex cross terms between categories

**Defer until:** Other models are complete.

---

## Implementation Strategy

### Phase 1: Complete Linear Model Hessian (IMMEDIATE)
1. Add sigma parameter Hessian terms to EvalLinear()
2. Validate with finite-difference test
3. Expected effort: 30 minutes

### Phase 2: Implement Probit Model Hessian (NEXT)
1. Study legacy implementation thoroughly
2. Implement 6-step algorithm above
3. Add helper functions for lambda calculations if needed
4. Validate with finite-difference test
5. Expected effort: 2-3 hours

### Phase 3: Testing and Validation
1. Extend test_gradient_validation.R to properly test Hessians
2. Lower threshold to 1e-3 as user suggested (Hessian is noisier)
3. Test on multiple parameter configurations
4. Expected effort: 1 hour

### Phase 4: Advanced Models (FUTURE)
- Ordered probit: 4-6 hours
- Multinomial logit: 4-6 hours

---

## Why Analytical Hessian Matters

### Computational Efficiency

**Current L-BFGS (gradient only):**
- Approximates Hessian from gradient history
- Requires storing past gradients (memory overhead)
- Approximation quality depends on curvature
- Typically requires more iterations to converge

**With Analytical Hessian (nlminb/trust):**
- Exact second-order information at each iteration
- Newton steps directly target optimum
- Typically **2-5x fewer iterations** to converge
- Each iteration is slightly more expensive (Hessian evaluation)
- **Overall: Usually 30-50% faster for well-conditioned problems**

### When Analytical Hessian Helps Most

1. **Near the optimum**: Newton methods have quadratic convergence
2. **Moderate dimensionality** (10-100 parameters): Hessian computation is still cheap
3. **Ill-conditioned problems**: Exact curvature helps navigate narrow valleys
4. **High-precision requirements**: Gets closer to true optimum

### When L-BFGS Is Fine

1. **Very high dimensions** (>100 parameters): Hessian becomes expensive
2. **Rough landscapes**: Far from optimum, quasi-Newton is robust
3. **Quick approximations**: If you just need "good enough"

### For Inference

- Standard errors = sqrt(diag(inv(-Hessian)))
- Confidence intervals require standard errors
- Currently SE calculation is incomplete (line 257+ in optimize_model.R)

### For Diagnostics

- Check positive definiteness → local minimum
- Eigenvalues → curvature of likelihood surface
- Conditioning number → numerical stability

---

## Current Workaround

The optimization workflow currently works well with **gradients only** using L-BFGS:
- L-BFGS approximates Hessian from gradient history
- Works very well in practice
- All validation tests pass

**Users can run optimizations NOW without waiting for full Hessian implementation.**

The Hessian is needed for:
- More advanced optimizers (trust region, Newton-CG)
- Accurate standard error calculation
- Model diagnostics

---

## Testing Methodology

When implementing Hessian terms, follow this validation process:

1. **Unit test each component:**
   ```r
   # Test sigma Hessian term
   test_params <- c(1.0, 1.0, 0.5, 0.8, 0.5)

   h_analytical <- evaluate_likelihood_cpp(fm_cpp, test_params,
                                          compute_hessian = TRUE)$hessian
   h_fd <- compute_fd_hessian(fm_cpp, test_params, delta = 1e-7)

   # Check relative error < 1e-3
   rel_error <- abs(h_analytical - h_fd) / abs(h_fd)
   ```

2. **Test multiple configurations:**
   - Different numbers of factors (1, 2, 3)
   - Different numbers of regressors
   - Mixed fixed/free loadings
   - Edge cases (parameters near zero, large values)

3. **Compare to legacy C++:**
   - Run same model in legacy backend
   - Extract Hessian from TModel output
   - Compare element-by-element

---

## Files to Modify

1. **src/Model.cpp**
   - `EvalLinear()` - Add sigma Hessian (lines 234+)
   - `EvalProbit()` - Replace placeholder (lines 237+)

2. **test_gradient_validation.R**
   - Update Test 2 to properly validate linear Hessian
   - Add Test 5 for probit Hessian
   - Adjust thresholds to 1e-3 for Hessian tests

3. **R/optimize_model.R** (LATER)
   - Complete standard error calculation (lines 177-179)
   - Convert upper triangle Hessian to full matrix
   - Invert and extract diagonal

---

## Mathematical Background

### Linear Model Second Derivatives

For Y ~ N(μ, σ²) where μ = β₀ + Σβᵢxᵢ + Σλⱼfⱼ:

```
L = log(φ((Y-μ)/σ)) = -log(σ) - (Y-μ)²/(2σ²)
```

**First derivatives:**
```
∂L/∂β = (Y-μ)·x/σ²
∂L/∂λ = (Y-μ)·f/σ²
∂L/∂σ = (Y-μ)²/σ³ - 1/σ
```

**Second derivatives:**
```
∂²L/∂β∂β' = -x·x'/σ²
∂²L/∂β∂λ = -x·f/σ²
∂²L/∂β∂σ = -2(Y-μ)·x/σ³
∂²L/∂λ∂λ' = -f·f'/σ²
∂²L/∂λ∂σ = -2(Y-μ)·f/σ³
∂²L/∂σ² = -3(Y-μ)²/σ⁴ + 1/σ²
```

For factor variances (θ = σ_f²), use chain rule:
```
∂/∂θ = ∂f/∂θ · ∂/∂f  where ∂f/∂σ_f² = x/(2σ_f)
```

### Probit Model Second Derivatives

For Φ(z) where z = β₀ + Σβᵢxᵢ + Σλⱼfⱼ, define λ(z) = φ(z)/Φ(z):

```
∂log(Φ)/∂β = λ(z)·x
∂²log(Φ)/∂β∂β' = -[z·λ(z) + λ(z)²]·x·x'
```

Similar expressions for other parameters with appropriate chain rule application.

---

## Next Steps

1. **Immediate:** Complete linear model sigma Hessian (~30 min)
2. **This week:** Implement probit Hessian (~3 hours)
3. **Next week:** Comprehensive validation testing
4. **Future:** Advanced models (oprobit, logit) as needed

The gradient implementation is **complete and validated**. The Hessian work is **optional enhancement** for advanced features.
