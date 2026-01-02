# Hessian Performance Analysis: Rcpp vs Legacy C++

## Problem Statement
The Rcpp implementation is approximately **8x slower** for Hessian calculation compared to the legacy C++ code.

## Recent Optimization (commit 993095a)
Changed Hessian accumulation loops from iterating over `nparam × nparam` to only iterating over `nparam_free × nparam_free`:

```cpp
// BEFORE (slow):
for (int i = 0; i < nparam; i++) {
    for (int j = i; j < nparam; j++) { ... }
}

// AFTER (faster):
for (int fi = 0; fi < nparam_free; fi++) {
    int i = freeparlist[fi];
    for (int fj = fi; fj < nparam_free; fj++) {
        int j = freeparlist[fj];
        ...
    }
}
```

This was applied to 4 locations in `CalcLkhd()`:
- Line 775-783: type_weighted_hess accumulation
- Line 1005-1013: hessilk computation from type_weighted_hess
- Line 1027-1035: totalhess accumulation from hessilk
- Line 1064-1072: full_hessL accumulation from totalhess

---

## Remaining Performance Gaps to Investigate

### 1. Vector Initialization (std::fill) - HIGH PRIORITY
Check if `std::fill` is called on full `nparam × nparam` arrays every iteration:

| Location | Current Behavior | Legacy Behavior |
|----------|------------------|-----------------|
| Line 393-394: `totalgrad`, `totalhess` | `std::fill` on full arrays | Only zeroes free param entries |
| Line 424-425: `gradilk`, `hessilk` | `std::fill` on full arrays | Only zeroes free param entries |
| Line 498-499: type_weighted_grad/hess | `std::fill` on full arrays | N/A (no type mixture in legacy?) |
| Line 509-510: grad_this_type/hess_this_type | `std::fill` on full arrays | N/A |

**Legacy code (lines 5545-5559):**
```cpp
if (stochflag==3) {
  if (predicting==0) {
    for (ifreepar = 0 ; ifreepar < nfreeparam; ifreepar++) {
      for (UInt_t jfreepar = ifreepar ; jfreepar < nfreeparam; jfreepar++) {
        Int_t fullHessindex = freeparlist[ifreepar]*nparam+freeparlist[jfreepar];
        totalhess[fullHessindex] = 0.0;  // Only zeroes free param entries!
      }
    }
  }
}
```

### 2. Gradient Loops Still Using nparam
Several loops still iterate over all `nparam` instead of just free parameters:

| Line | Loop | Notes |
|------|------|-------|
| 726 | type_weighted_grad accumulation | Can skip fixed params |
| 831 | Type loading cross-terms | Complex, may need full loop |
| 867 | Factor variance cross-terms | Complex, may need full loop |
| 996-997 | gradilk from type_weighted_grad | Can skip fixed params |
| 1020-1021 | totalgrad accumulation | Can skip fixed params |
| 1056-1057 | full_gradL accumulation | Can skip fixed params |

### 3. Storage Allocation
Check if these are allocated once vs. per-iteration:

| Array | Size | Allocation Pattern |
|-------|------|-------------------|
| `hessilk` | `nparam * nparam` | Pre-allocated line 376 |
| `totalhess` | `nparam * nparam` | Pre-allocated line 369 |
| `type_weighted_hess` | `nparam * nparam` | Pre-allocated line 383 |
| `hess_this_type` | `nparam * nparam` | Pre-allocated line 386 |
| `full_hessL` | `nparam * nparam` | Pre-allocated line 354 |

✅ All pre-allocated outside observation loop - this is correct.

### 4. Model.cpp Hessian Computation
Need to check if `Model::Eval()` has similar inefficiencies.

---

## Detailed Comparison: Legacy vs Rcpp Structure

### Legacy CalcLkhd (TMinLkhd.cc)
```
Loop over observations (i)
  Zero totalhess ONLY for free param pairs
  Loop over mixtures (imix)
    Loop over types (itype)
      Loop over integration points (intpt)
        Zero hessilk ONLY for free param pairs
        Loop over models (imod)
          Model.Eval() -> get modhess
          Accumulate into hessilk
        Accumulate hessilk into totalhess (only free params)
    Convert totalhess to hessL (only free params)
```

### Rcpp CalcLkhd (FactorModel.cpp)
```
Loop over observations (iobs)
  std::fill totalhess (ALL entries)  <-- INEFFICIENT
  Loop over integration points (intpt)
    std::fill hessilk (ALL entries)   <-- INEFFICIENT
    std::fill type_weighted_hess (ALL entries)  <-- INEFFICIENT
    Loop over types (ityp)
      std::fill hess_this_type (ALL entries)  <-- INEFFICIENT
      Loop over models (imod)
        Model.Eval() -> get modHess
        Accumulate into hess_this_type
      Accumulate hess_this_type into type_weighted_hess (only free)  ✓
    Convert type_weighted_hess to hessilk (only free)  ✓
    Accumulate hessilk into totalhess (only free)  ✓
  Convert totalhess to full_hessL (only free)  ✓
```

---

## Proposed Fixes

### Fix 1: Selective Zeroing of Hessian Arrays
Replace `std::fill` with selective zeroing only for free parameter entries:

```cpp
// Instead of:
if (iflag == 3) std::fill(totalhess.begin(), totalhess.end(), 0.0);

// Use:
if (iflag == 3) {
    for (int fi = 0; fi < nparam_free; fi++) {
        int i = freeparlist[fi];
        for (int fj = fi; fj < nparam_free; fj++) {
            int j = freeparlist[fj];
            totalhess[i * nparam + j] = 0.0;
        }
    }
}
```

Locations to fix:
- Line 394: totalhess
- Line 425: hessilk
- Line 499: type_weighted_hess
- Line 510: hess_this_type

### Fix 2: Optimize Gradient Loops
For gradient accumulation, consider skipping fixed parameters.

---

## Performance Testing

### Test Setup (Server Benchmark)
- Model: Production model
- Number of parameters: **304**
- Data: 5k and 10k observations

### Benchmark Results

**Legacy C++ (100k observations, 304 parameters):**
| Operation | Time |
|-----------|------|
| Log-likelihood | 2.20 sec |
| Gradient | 4.93 sec |
| Hessian | 264.06 sec |

**Rcpp (10k and 5k observations, 304 parameters):**
| Observations | Log-likelihood | Gradient | Hessian | Hessian/Gradient Ratio |
|--------------|----------------|----------|---------|------------------------|
| 10,000 | 1.561 sec | 2.375 sec | 186.013 sec | **78×** |
| 5,000 | 0.882 sec | 1.474 sec | 94.794 sec | **64×** |

**Rcpp (20k observations):**
| Operation | Time |
|-----------|------|
| Log-likelihood | 1.669 sec |
| Gradient | 3.815 sec |
| Hessian | 385.801 sec |

**Direct Comparison (100k observations, 304 parameters):**
| Operation | Legacy C++ | Rcpp | Rcpp / Legacy |
|-----------|------------|------|---------------|
| Log-likelihood | 2.20 sec | 7.14 sec | **3.2×** slower |
| Gradient | 4.93 sec | 15.98 sec | **3.2×** slower |
| Hessian | 264.06 sec | 2120.81 sec | **8.0×** slower |

**Key Findings:**
1. Rcpp is uniformly ~3× slower for likelihood and gradient
2. Rcpp is **8× slower** for Hessian specifically
3. The Hessian has **additional** performance issues beyond the general slowdown

**Key Observations:**
1. Hessian time scales **linearly with observations** (186/94.8 ≈ 1.96× for 2× obs)
2. This confirms the overhead is NOT the issue - it's per-observation computation
3. Hessian takes **~78× longer than gradient** - this is very high
4. Expected ratio should be ~O(nparam) ≈ 304× in worst case, but we're doing work proportional to nparam² per observation

### Expected Complexity Analysis

For 304 parameters:
- Gradient: O(nobs × nint × nmodels × nparam)
- Hessian: O(nobs × nint × nmodels × nparam²)
- Expected ratio: ~nparam = 304×, but we see only 78×

This suggests either:
1. Gradient is slower than optimal, OR
2. Hessian is actually faster than naive O(nparam²), OR
3. The 78× ratio indicates significant unnecessary work in Hessian

### std::fill Impact Estimate

With nparam = 304:
- Full matrix: 304² = 92,416 elements
- Upper triangle: 304×305/2 = 46,360 elements
- Ratio: only **2×** improvement from selective zeroing

**Conclusion: std::fill fix alone will NOT explain 8× slowdown. Must find other issues.**

---

## Investigation Log

### 2026-01-02: Initial Analysis
1. Reviewed git diff showing recent optimization (commit 993095a)
2. Identified 4 Hessian accumulation loops were optimized to use nparam_free
3. Found remaining issues:
   - `std::fill` still zeroes entire nparam×nparam arrays
   - Some gradient loops still use nparam
4. Next: Profile to confirm std::fill is the bottleneck

### 2026-01-02: Root Cause Identified - std::fill on Full Arrays

**The main bottleneck is `std::fill` being called on full `nparam × nparam` arrays.**

Found 6 `std::fill` calls in the hot path (FactorModel.cpp):

| Line | Array | Size | Frequency |
|------|-------|------|-----------|
| 393 | totalgrad | nparam | 1× per obs |
| **394** | **totalhess** | **nparam²** | **1× per obs** |
| 424 | gradilk | nparam | 1× per intpt |
| **425** | **hessilk** | **nparam²** | **1× per intpt** |
| 498 | type_weighted_grad | nparam | 1× per intpt |
| **499** | **type_weighted_hess** | **nparam²** | **1× per intpt** |
| 509 | grad_this_type | nparam | 1× per type×intpt |
| **510** | **hess_this_type** | **nparam²** | **1× per type×intpt** |

**Complexity Analysis:**

For typical model: nparam = 100, nparam_free = 50, nobs = 500, nint = 8, ntyp = 1

| Array | Current (std::fill) | Optimal (selective) | Ratio |
|-------|---------------------|---------------------|-------|
| hessilk per intpt | 10,000 ops | 1,275 ops | 7.8× |
| Total hessilk | 500×8×10,000 = 40M | 500×8×1,275 = 5.1M | 7.8× |

**This explains the 8× slowdown!**

**Legacy code only zeroes free parameter entries:**
```cpp
// TMinLkhd.cc lines 5545-5559
for (ifreepar = 0 ; ifreepar < nfreeparam; ifreepar++) {
    for (UInt_t jfreepar = ifreepar ; jfreepar < nfreeparam; jfreepar++) {
        Int_t fullHessindex = freeparlist[ifreepar]*nparam+freeparlist[jfreepar];
        totalhess[fullHessindex] = 0.0;  // Only free param entries!
    }
}
```

---

## Potential Optimization Areas Identified

### 1. Repeated sqrt() Calls in Hot Loops ⚠️ MAJOR ISSUE

In `FactorModel::CalcLkhd()` lines 680-708 (modHess mapping to global Hessian):
```cpp
for (int i = 0; i < nDimModHess; i++) {
    for (int j = i; j < nDimModHess; j++) {
        if (i < nfac) {
            double sigma_i = std::sqrt(factor_var[i]);  // sqrt() in INNER LOOP!
            double x_node_i = quad_nodes[facint[i]];
            chain_factor *= x_node_i / (2.0 * sigma_i);
        }
        if (j < nfac) {
            double sigma_j = std::sqrt(factor_var[j]);  // sqrt() in INNER LOOP!
            // ...
        }
        if (i < nfac && i == j) {
            double sigma = std::sqrt(factor_var[i]);    // THIRD sqrt() call!
            // ...
        }
    }
}
```

This is nested inside: model loop → type loop → integration point loop → observation loop.

**Impact calculation:**
- 100k observations × 8 integration points × 10 models × 50 Hessian entries = 400 million sqrt() calls!
- sqrt() is expensive (~20-30 cycles)

**Additional locations with repeated sqrt():**
- Lines 576-594: Gradient chain rule inside type loop
- Lines 627-628: Correlated factor Hessian
- Lines 787-802: Type probability gradient
- Lines 887-898: Type probability Hessian

**Fix:** Pre-compute `sigma[k] = sqrt(factor_var[k])` and `x_nodes[k] = quad_nodes[facint[k]]`
at the start of each integration point, then reuse throughout.

### 2. Array Resizing Inside Model::Eval()

In `Model::Eval()` and `Model::EvalLinear()`:
```cpp
modEval.resize(ngrad);
std::fill(modEval.begin(), modEval.end(), 0.0);
hess.resize(npar * npar);
```

This happens once per model per integration point per observation.
For 100k obs × 8 quad points × 10 models = 8 million resize/fill operations.

**Fix:** Pre-allocate modEval and hess vectors and reuse them.

### 3. Vector Allocation in Legacy vs Rcpp

Legacy code creates `vector<double> expres` on stack (line 366):
```cpp
vector<double> expres(1,0.0);
```

Rcpp code uses `static thread_local` (line 89-90):
```cpp
static thread_local std::vector<double> expres;
expres.assign(1, 0.0);
```

This is actually an optimization, but the `.assign()` still does work.

### 4. Hessian Cross-Term Computations

The Hessian cross-term loops (factor variance × model params) iterate over ALL model parameters even when many are fixed. This could be optimized to skip fixed parameters.

### 5. Distribution Functions

Rcpp uses R's distribution functions via wrappers:
```cpp
modEval[0] = normal_pdf(Z, 0.0, sigma);  // Calls R::dnorm
```

Legacy uses ROOT:
```cpp
modEval[0] = ROOT::Math::normal_pdf(Z, sigma);
```

R's dnorm may have additional overhead (bounds checking, special case handling).

---

## Questions to Answer
1. What percentage of time is spent in std::fill vs actual computation?
2. Does Model::Eval() have similar inefficiencies?
3. Is the type model (ntyp > 1) case being tested?
4. Are there any algorithmic differences between legacy and Rcpp?
