#!/usr/bin/env Rscript
# Verify that define_factor_model no longer has loading_normalization parameter

library(factorana)

cat("=== Testing that factor model no longer has loading_normalization ===\n\n")

# Should work: define factor model without loading_normalization
fm <- define_factor_model(
  n_factors = 1,
  n_types = 1,
  n_quad = 8
)

cat("✓ Factor model created successfully\n")
print(fm)

# Verify factor model does not have loading_normalization field
if ("loading_normalization" %in% names(fm)) {
  stop("❌ Factor model still has loading_normalization field!")
} else {
  cat("\n✓ Factor model does NOT have loading_normalization field (correct)\n")
}

# Should error: trying to pass loading_normalization
cat("\n=== Testing that loading_normalization parameter is rejected ===\n")
tryCatch({
  fm_bad <- define_factor_model(
    n_factors = 1,
    n_types = 1,
    n_quad = 8,
    loading_normalization = 1.0
  )
  stop("❌ Should have rejected loading_normalization parameter!")
}, error = function(e) {
  if (grepl("unused argument", e$message)) {
    cat("✓ Correctly rejected loading_normalization parameter\n")
  } else {
    stop("❌ Wrong error: ", e$message)
  }
})

cat("\n=== Testing component-specific normalizations work ===\n")

# Create simple data
set.seed(100)
n <- 100
Y <- rnorm(n)
x1 <- rnorm(n)
dat <- data.frame(Y = Y, x1 = x1, intercept = 1, eval = 1)

# Component with fixed loading
mc_fixed <- define_model_component(
  name = "Y",
  data = dat,
  outcome = "Y",
  factor = fm,
  covariates = c("intercept", "x1"),
  model_type = "linear",
  evaluation_indicator = "eval",
  loading_normalization = 1.0  # Fixed at 1.0
)

cat("✓ Component with fixed loading created\n")
cat("  Loading normalization:", mc_fixed$loading_normalization, "\n")

# Component with free loading
mc_free <- define_model_component(
  name = "Y",
  data = dat,
  outcome = "Y",
  factor = fm,
  covariates = c("intercept", "x1"),
  model_type = "linear",
  evaluation_indicator = "eval",
  loading_normalization = NA_real_  # Free
)

cat("✓ Component with free loading created\n")
cat("  Loading normalization:", mc_free$loading_normalization, "\n")

# Component with default (all free)
mc_default <- define_model_component(
  name = "Y",
  data = dat,
  outcome = "Y",
  factor = fm,
  covariates = c("intercept", "x1"),
  model_type = "linear",
  evaluation_indicator = "eval"
  # No loading_normalization specified
)

cat("✓ Component with default normalization created\n")
cat("  Loading normalization:", mc_default$loading_normalization, "\n")

if (is.na(mc_default$loading_normalization)) {
  cat("  ✓ Default is NA (free) as expected\n")
} else {
  stop("❌ Default should be NA (free)")
}

cat("\n=== All tests passed! ===\n")
cat("\nSummary:\n")
cat("  - Factor model no longer has loading_normalization parameter ✓\n")
cat("  - Loading normalizations are now component-specific ✓\n")
cat("  - Components can have fixed, free, or default (free) loadings ✓\n")
