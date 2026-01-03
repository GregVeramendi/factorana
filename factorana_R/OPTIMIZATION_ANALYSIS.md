# Deep Performance Analysis: Remaining Optimization Opportunities

## Current Performance (v0.2.10, 10k obs, 32 cores, 304 params)
- Log-likelihood: 0.457s
- Gradient: 1.039s
- Hessian: 77.4s

## Comparison with Legacy C++ (scaled to 100k obs)
| Operation | Legacy | Rcpp (est) | Ratio |
|-----------|--------|------------|-------|
| Log-likelihood | 2.20s | ~4.5s | 2.0× slower |
| Gradient | 4.93s | ~10s | 2.0× slower |
| Hessian | 264s | ~770s | 2.9× slower |

---

## Optimization Opportunity 1: Custom Distribution Functions (HIGH IMPACT)

### Problem
Rcpp uses R's C API (`R::dnorm`, `R::pnorm`) which include:
- Bounds checking and NA handling
- Special case handling for edge values
- Function call overhead (not inlined)
- No SIMD vectorization

Legacy uses ROOT's highly optimized `ROOT::Math::normal_pdf/cdf` which:
- Use SIMD instructions when available
- Are fully inlined
- Have minimal overhead

### Proposed Solution
Replace R's distribution functions with fast inline implementations:

```cpp
// Fast normal PDF - assumes valid inputs, no NA checking
inline double fast_normal_pdf(double x, double sigma = 1.0) {
    constexpr double inv_sqrt_2pi = 0.3989422804014327;
    double z = x / sigma;
    return (inv_sqrt_2pi / sigma) * std::exp(-0.5 * z * z);
}

// Fast normal CDF using Abramowitz & Stegun approximation
// Maximum error < 7.5e-8 for all x
inline double fast_normal_cdf(double x) {
    constexpr double a1 =  0.254829592;
    constexpr double a2 = -0.284496736;
    constexpr double a3 =  1.421413741;
    constexpr double a4 = -1.453152027;
    constexpr double a5 =  1.061405429;
    constexpr double p  =  0.3275911;

    int sign = (x < 0) ? -1 : 1;
    x = std::fabs(x) / 1.4142135623730951;  // sqrt(2)

    double t = 1.0 / (1.0 + p * x);
    double t2 = t * t;
    double t3 = t2 * t;
    double t4 = t3 * t;
    double t5 = t4 * t;
    double y = 1.0 - (a1*t + a2*t2 + a3*t3 + a4*t4 + a5*t5) * std::exp(-x*x);

    return 0.5 * (1.0 + sign * y);
}
```

### Expected Impact
- Distribution functions are called once per model per observation per integration point
- For 10k obs × 8 quad × 10 models = 800,000 calls per likelihood evaluation
- Could save 20-40% of likelihood/gradient time

### Validation Required
- Verify accuracy vs R::dnorm/pnorm across range [-10, 10]
- Ensure no numerical issues at extreme values

---

## Optimization Opportunity 2: Hessian Loop Restructuring (MEDIUM IMPACT)

### Problem
The inner Hessian loop (FactorModel.cpp lines 704-728) has:
- Conditional branches inside inner loop (`if (i < nfac)`)
- Non-contiguous memory access due to index mapping

```cpp
for (int i = 0; i < nDimModHess; i++) {
    for (int j = i; j < nDimModHess; j++) {
        // ... conditionals and scattered memory access
    }
}
```

### Proposed Solution
Split into separate loops for different regions:

```cpp
// 1. Factor-factor block (no chain rule)
for (int i = 0; i < nfac; i++) {
    for (int j = i; j < nfac; j++) {
        int idx = i * nparam + j;
        double chain = df_dsigma2[i] * df_dsigma2[j];
        hess[idx] += modHess[i * nDimModHess + j] * chain;
        if (i == j) hess[idx] += modEval[i+1] * d2f_dsigma2_sq[i];
    }
}

// 2. Factor-param cross block
for (int i = 0; i < nfac; i++) {
    for (int j = nfac; j < nDimModHess; j++) {
        int full_j = firstpar + j - nfac;
        int idx = i * nparam + full_j;
        hess[idx] += modHess[i * nDimModHess + j] * df_dsigma2[i];
    }
}

// 3. Param-param block (no chain rule)
for (int i = nfac; i < nDimModHess; i++) {
    for (int j = i; j < nDimModHess; j++) {
        int full_i = firstpar + i - nfac;
        int full_j = firstpar + j - nfac;
        hess[full_i * nparam + full_j] += modHess[i * nDimModHess + j];
    }
}
```

### Expected Impact
- Eliminates branch mispredictions in hot loop
- Better instruction cache utilization
- Estimated 10-15% improvement in Hessian time

---

## Optimization Opportunity 3: Cache-Friendly Memory Layout (MEDIUM IMPACT)

### Problem
With nparam = 304:
- Full Hessian: 304 × 304 × 8 bytes = 740 KB
- L2 cache: typically 256-512 KB
- L3 cache: typically 8-16 MB (shared)

Random access to hess_this_type causes cache misses when nparam is large.

### Proposed Solution
Use blocked/tiled access pattern:

```cpp
constexpr int BLOCK_SIZE = 64;  // Fits in L1 cache

for (int bi = 0; bi < nparam; bi += BLOCK_SIZE) {
    for (int bj = bi; bj < nparam; bj += BLOCK_SIZE) {
        // Process block [bi:bi+BLOCK_SIZE, bj:bj+BLOCK_SIZE]
        for (int i = bi; i < std::min(bi + BLOCK_SIZE, nparam); i++) {
            for (int j = std::max(bj, i); j < std::min(bj + BLOCK_SIZE, nparam); j++) {
                // ... accumulate
            }
        }
    }
}
```

### Expected Impact
- Reduces L2/L3 cache misses
- Most beneficial for large parameter counts (>200)
- Estimated 5-10% improvement for Hessian

---

## Optimization Opportunity 4: Vectorization with Eigen (MEDIUM IMPACT)

### Problem
Current code uses scalar operations; modern CPUs can do 4-8 doubles per cycle with AVX/AVX2.

### Proposed Solution
Use Eigen's vectorized operations for bulk array operations:

```cpp
#include <Eigen/Dense>

// Instead of loop over nparam entries:
Eigen::Map<Eigen::VectorXd> grad(totalgrad.data(), nparam);
Eigen::Map<Eigen::VectorXd> ilk(gradilk.data(), nparam);
grad += probilk * ilk;  // Vectorized

// For Hessian outer product:
Eigen::Map<Eigen::MatrixXd> hess(totalhess.data(), nparam, nparam);
Eigen::Map<Eigen::VectorXd> g(gradilk.data(), nparam);
hess.noalias() += probilk * (hess_ilk + g * g.transpose());
```

### Expected Impact
- Could provide 2-4× speedup for large array operations
- Most beneficial for gradient accumulation
- Requires careful integration to avoid overhead

---

## Optimization Opportunity 5: Reduce Type Model Overhead (LOW IMPACT for ntyp=1)

### Problem
Lines 842-1000+ in CalcLkhd handle type model Hessian terms. When ntyp=1, these should be completely skipped, but there may be overhead from checking conditions.

### Current Code
```cpp
if (ntyp > 1) {
    // Complex nested loops for type loading Hessians
}
```

### Verification Needed
- Confirm that ntyp=1 case truly skips all type-related Hessian work
- Profile to verify no hidden overhead

---

## Implementation Priority

| Optimization | Impact | Effort | Risk |
|--------------|--------|--------|------|
| 1. Fast distribution functions | HIGH | LOW | LOW |
| 2. Loop restructuring | MEDIUM | MEDIUM | LOW |
| 3. Cache blocking | MEDIUM | MEDIUM | LOW |
| 4. Eigen vectorization | MEDIUM | HIGH | MEDIUM |
| 5. Type model skip | LOW | LOW | LOW |

### Recommended Order
1. **Fast distribution functions** - Quick win, easy to implement and test
2. **Loop restructuring** - Clear improvement, moderate effort
3. **Type model verification** - Quick check
4. **Cache blocking** - If still needed after 1-2
5. **Eigen vectorization** - More invasive, consider last

---

## Validation Strategy

For each optimization:
1. Run FD tests to verify gradient/Hessian accuracy
2. Run full test suite (513+ tests)
3. Benchmark on 10k observations
4. Compare to legacy C++ results if possible

---

## Theoretical Performance Limit

Given that:
- Rcpp has inherent overhead vs standalone C++
- R's memory management differs from ROOT
- Compiler optimizations may differ

Realistic target: **~1.5-2× of legacy performance** (vs current 2.9×)

With all optimizations: Hessian could potentially reach ~500-600s for 100k obs (vs current ~770s estimated).
