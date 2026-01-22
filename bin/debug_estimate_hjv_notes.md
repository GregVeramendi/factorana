# Debug Notes: estimate_hjv.R Performance Investigation

## Problem Summary
The R package (factorana) Hessian computation is ~2x slower than the legacy C++ code when the switching model is included.

## Benchmark Results (50k observations)

| Metric | R (v0.2.38) | R (v0.2.35) | Legacy |
|--------|-------------|-------------|--------|
| Log-likelihood | 0.448s | 0.457s | 0.200s |
| Gradient | 0.979s | 0.893s | 0.630s |
| Hessian | 1285.7s | 1220.2s | 600.4s |

**Ratio (R/Legacy):** LL ~2.2x, Grad ~1.5x, Hess ~2.1x slower

### Model Configuration Tests

| Configuration | R vs Legacy Speed |
|---------------|-------------------|
| Full model without switching component | **Same speed** |
| First 4 components (adveng, advmath, HS_Track, application) | **Same speed** |
| Full model with switching component | ~2x slower |

**Key observation:** The exploded logit (application model with 13 ranks × 13 choices) runs at legacy speed. Only the regular multinomial logit (switching model with 1 rank × 12 choices) causes the slowdown. This is counterintuitive since the application model has ~15x more work per observation.

## Parameter Count Discrepancy

| | R Code | Legacy |
|--|--------|--------|
| Total params | 1266 | 1266 |
| Free params | 458 | 458 |
| Fixed params | 802 | 808 |
| **Unaccounted** | **6** | 0 |

R code: 458 + 802 = 1260, but reports 1266 total (6 unaccounted)
Legacy: 458 + 808 = 1266 total

**Key observation:** When excluding application and switching models, parameter counts match between R and legacy.

## Questions to Investigate

1. **Where do the 6 extra parameters come from?**
   - Tied parameters from equality constraints?
   - Type intercepts counted differently?
   - Factor loadings handled differently?

2. **Switching model configuration:**
   - Number of choices (n_choices)?
   - Number of factors affecting it?
   - Loading normalizations?
   - Fixed coefficients via fix_coefficient()?

3. **Is the fast path being used?**
   - Check if `use_fast_path` evaluates to true for the switching model
   - Conditions: all_loadings_free && no quadratic && no interaction

4. **Application model vs Switching model:**
   - Application (exploded logit) runs at legacy speed
   - Switching (regular logit) is 2x slower
   - Both use EvalLogit - what's different?

## Code Changes Made (v0.2.38)

1. Added fast path in EvalLogit Hessian (Model.cpp:1327-1427)
2. Cached regressor values to avoid getRegValue() function call overhead
3. Removed `use_free_opt` from fast path condition
4. Fast path matches legacy TModel.cc structure

## Investigation Findings (2026-01-22)

### Root Cause of 6 Parameter Discrepancy

The discrepancy comes from `build_parameter_metadata()` in R generating FEWER parameters than `initialize_parameters()` for multinomial logit models.

**Naming Convention Mismatch:**
- `initialize_parameters.R` uses: `{comp}_c{choice}_loading_{k}` (e.g., `educDSwitchtoMajor_c1_loading_1`)
- `build_parameter_metadata` uses: `{comp}_loading_{k}_alt{alt}` (e.g., `educDSwitchtoMajor_loading_1_alt1`)

**Parameter Count Sources:**
1. `n_params_total` comes from C++ (`get_parameter_info_cpp`)
2. `param_metadata$n_params` comes from R's `build_parameter_metadata()`
3. These are NOT validated to match!

**Key Issue:**
- C++: `regressor_idx.size() + n_free_loadings` per alternative
- R metadata: `1 (intercept) + length(covariates excluding "intercept") + n_free_loadings` per alternative

For multinomial logit where `comp$intercept = TRUE` but "intercept" is not in covariates:
- C++ counts: N covariates + 3 loadings = N+3 per alternative
- R metadata counts: 1 intercept + N covariates + 3 loadings = N+4 per alternative

This would actually make R have MORE params, but the user sees FEWER. Need to investigate if intercept is being added differently.

### loadings Listed as Fixed

The `educDSwitchtoMajor_loading_*_alt11` loadings at the END of the fixed list are suspicious:
- Model has `loading_normalization = c(NA_real_, NA_real_, NA_real_)` → all loadings FREE
- Yet they appear in "Fixed parameters (802)"
- This suggests parameter indexing is misaligned

### Next Investigation Steps

1. **Add validation** that `param_metadata$n_params == n_params_total`
2. **Debug print** both parameter lists to find where they diverge
3. **Check intercept handling** for multinomial logit in both:
   - `initialize_parameters.R` (line ~488)
   - `build_parameter_metadata` (line ~207)

## Memory Allocation Optimization (v0.2.39)

### Issue Identified

Inside `Model::EvalLogit()`, several vectors were being allocated fresh on EVERY call:
- `free_loading_cumcount` (numfac+1 elements)
- `chosen` (numchoice elements)
- `logitgrad` (numchoice × npar elements for Hessian)
- `hess_idx_free` (npar elements)
- `rankedChoiceCorr` (numchoice-1 elements)
- `pdf` (numchoice elements)
- `reg_values` (nregressors elements) in fast path

For a model with 50k observations × 64 quadrature points (8² for 2 factors) = 3.2M calls to EvalLogit per likelihood evaluation, these allocations cause significant overhead.

**Example: `logitgrad` allocation**
For switching model: 12 choices × ~696 params = 8,352 doubles × 8 bytes = 67KB per call
For 3.2M calls: cumulative allocation requests of ~200TB (obviously reused by allocator, but still overhead)

### Fix Implemented

Changed all hot-path vectors to `static thread_local`:
```cpp
// Before:
std::vector<bool> chosen(numchoice, false);
std::vector<double> logitgrad;

// After:
static thread_local std::vector<bool> chosen;
if (chosen.size() < static_cast<size_t>(numchoice)) {
    chosen.resize(numchoice);
}
std::fill(chosen.begin(), chosen.begin() + numchoice, false);

static thread_local std::vector<double> logitgrad;
// ...resize only if needed, then fill
```

This matches the pattern already used for `expres` at Model.cpp:119.

### Expected Impact

The 2x slowdown is likely a combination of:
1. Memory allocation overhead (addressed by this fix)
2. Possible differences in model structure/evaluation count (still investigating)

The allocation overhead affects ALL evaluation levels (LL, Grad, Hess) proportionally, which matches the benchmark pattern where all three are ~2x slower.

## Next Steps

- [x] Identify that param_metadata count doesn't match C++ count
- [x] Add validation to warn on param_metadata/C++ count mismatch
- [x] Optimize EvalLogit memory allocation with thread_local vectors
- [x] Fix "Free parameter size mismatch" error in optimization (v0.2.40)
- [ ] Fix param_metadata to match C++ structure exactly
- [ ] Add debug output to verify fast path is being taken
- [ ] Re-run benchmarks with v0.2.40 to measure improvement

## Bug Fix: Free Parameter Size Mismatch (v0.2.40)

**Issue:** Optimization failed with "Free parameter size mismatch" error after benchmark completed successfully.

**Root Cause:** Line 1310 in optimize_model.R extracted free parameters using R's `param_constraints$free_idx`:
```r
current_start <- full_init_params[param_constraints$free_idx]  # WRONG
```

But `param_constraints$free_idx` comes from `param_metadata` which may have different parameter indexing than C++. The benchmark used `init_params` (extracted via C++'s `extract_free_params_cpp()`) and worked correctly.

**Fix:** Use `init_params` directly instead of re-extracting from R metadata:
```r
current_start <- init_params  # CORRECT - uses C++ extracted free params
```

## Files
- Model file: `/Users/mendi/Dropbox (Personal)/Greg/Work/Research/Project-code/factorana/bin/estimate_hjv.R`
- R package: `/Users/mendi/Dropbox (Personal)/Greg/Work/Research/Project-code/factorana/factorana_R/`
- Legacy code: `/Users/mendi/Dropbox (Personal)/Greg/Work/Research/Project-code/factorana/libraries/`
