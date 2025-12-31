# =============================================================================
# Stage 1 Estimation Script for Major Choice Project
# Translated from bin/prob_hjv.cc (lines 80-750)
# =============================================================================
#
# 3-factor model with 8 quadrature points, no correlation
# - Factor 1: Academic/cognitive ability
# - Factor 2: Non-cognitive (emotional stability, leadership)
# - Factor 3: General ability
#
# Data: /home/Projects/MajorChoice/intdata/MajorChoice/majorchoice_wfac_v14f.dta

library(factorana)
library(haven)

# =============================================================================
# Load Data
# =============================================================================
data_path <- "/home/Projects/MajorChoice/intdata/MajorChoice/majorchoice_wfac_v14f.dta"
dat <- haven::read_dta(data_path)

# Convert to data.frame (factorana expects data.frame, not tibble)
dat <- as.data.frame(dat)

# Add intercept column
dat$constant <- 1

# =============================================================================
# Define Covariate Sets
# =============================================================================

# Background variables (common across most models)
bkgnd1990 <- c(
  "college_mom",
  "hs_mom", "hs_mom_miss",
  "college_dad",
  "hs_dad", "hs_dad_miss",
  "fam_inc73_mom",
  "health1", "health2", "health_missing",
  "schoolave_faminc",
  "constant"
)

# Background without intercept (for models that add it separately)
bkgnd1990_nocons <- setdiff(bkgnd1990, "constant")

# =============================================================================
# Define Factor Model
# =============================================================================

# 3-factor model, no types (n_types = 1 means single type = no discrete heterogeneity)
fm <- define_factor_model(
  n_factors = 3,
  n_types = 1,  # Single type = no unobserved discrete heterogeneity
  corr = FALSE
)

# =============================================================================
# Define Model Components
# =============================================================================

# Loading normalization conventions in factorana_R:
#   NA_real_ = free parameter (to be estimated)
#   0        = loading fixed to 0 (factor does not load on this outcome)
#   1        = loading fixed to 1 (scale normalization)

components <- list()

# -----------------------------------------------------------------------------
# 9th Grade Measurements (all 3 factors load freely)
# -----------------------------------------------------------------------------

norm_all_free <- c(NA_real_, NA_real_, NA_real_)

# seng9: 9th Grade English Grade
components$seng9 <- define_model_component(
  name = "seng9",
  data = dat,
  outcome = "seng9",
  factor = fm,
  evaluation_indicator = "ind_seng9",
  covariates = c(bkgnd1990, "adveng"),
  model_type = "linear",
  loading_normalization = norm_all_free
)

# smath9: 9th Grade Math Grade
components$smath9 <- define_model_component(
  name = "smath9",
  data = dat,
  outcome = "smath9",
  factor = fm,
  evaluation_indicator = "ind_smath9",
  covariates = c(bkgnd1990, "advmath"),
  model_type = "linear",
  loading_normalization = norm_all_free
)

# ssports9: 9th Grade Sports Grade
components$ssports9 <- define_model_component(
  name = "ssports9",
  data = dat,
  outcome = "ssports9",
  factor = fm,
  evaluation_indicator = "ind_ssports9",
  covariates = bkgnd1990,
  model_type = "linear",
  loading_normalization = norm_all_free
)

# sswe9: 9th Grade Swedish Grade
components$sswe9 <- define_model_component(
  name = "sswe9",
  data = dat,
  outcome = "sswe9",
  factor = fm,
  evaluation_indicator = "ind_sswe9",
  covariates = bkgnd1990,
  model_type = "linear",
  loading_normalization = norm_all_free
)

# sgpa9: 9th Grade GPA (includes other 9th grade scores as covariates)
components$sgpa9 <- define_model_component(
  name = "sgpa9",
  data = dat,
  outcome = "sgpa9",
  factor = fm,
  evaluation_indicator = "ind_sgpa9",
  covariates = c(bkgnd1990,
                 "smath9", "seng9",
                 "ssports9", "ssports9_miss",
                 "sswe9", "sswe9_miss"),
  model_type = "linear",
  loading_normalization = norm_all_free
)

# -----------------------------------------------------------------------------
# 10th Grade / High School Measurements
# -----------------------------------------------------------------------------

# ssports10: 10th Grade Sports Grade
components$ssports10 <- define_model_component(
  name = "ssports10",
  data = dat,
  outcome = "ssports10",
  factor = fm,
  evaluation_indicator = "ind_ssports10",
  covariates = c(bkgnd1990, "HS_track3", "HS_track4"),
  model_type = "linear",
  loading_normalization = norm_all_free
)

# smath10: HS Math Grade
components$smath10 <- define_model_component(
  name = "smath10",
  data = dat,
  outcome = "smath10",
  factor = fm,
  evaluation_indicator = "ind_smath10",
  covariates = c(bkgnd1990, "advmath", "HS_track3", "HS_track4"),
  model_type = "linear",
  loading_normalization = norm_all_free
)

# sgpahs_voc: HS GPA (Vocational track)
components$sgpahs_voc <- define_model_component(
  name = "sgpahs_voc",
  data = dat,
  outcome = "gpahs",
  factor = fm,
  evaluation_indicator = "ind_gpahs_voc",
  covariates = c(bkgnd1990,
                 "adveng", "advmath",
                 "ssports10", "ssports10_miss",
                 "smath10", "smath10_miss"),
  model_type = "linear",
  loading_normalization = norm_all_free
)

# sgpahs_aca: HS GPA (Academic track)
components$sgpahs_aca <- define_model_component(
  name = "sgpahs_aca",
  data = dat,
  outcome = "gpahs",
  factor = fm,
  evaluation_indicator = "ind_gpahs_aca",
  covariates = c(bkgnd1990,
                 "adveng", "advmath",
                 "ssports10", "ssports10_miss",
                 "smath10", "smath10_miss",
                 "HS_track4"),
  model_type = "linear",
  loading_normalization = norm_all_free
)

# -----------------------------------------------------------------------------
# Cognitive Measures (only factor 1 loads)
# -----------------------------------------------------------------------------

norm_factor1_only <- c(NA_real_, 0, 0)

# cog1: Cognitive measure 1
components$cog1 <- define_model_component(
  name = "cog1",
  data = dat,
  outcome = "cog1",
  factor = fm,
  evaluation_indicator = "ind_cog1",
  covariates = c(bkgnd1990,
                 "adveng", "advmath",
                 "HS_track2", "HS_track3", "HS_track4"),
  model_type = "linear",
  loading_normalization = norm_factor1_only
)

# cog2: Cognitive measure 2
components$cog2 <- define_model_component(
  name = "cog2",
  data = dat,
  outcome = "cog2",
  factor = fm,
  evaluation_indicator = "ind_cog2",
  covariates = c(bkgnd1990,
                 "adveng", "advmath",
                 "HS_track2", "HS_track3", "HS_track4"),
  model_type = "linear",
  loading_normalization = norm_factor1_only
)

# cog3: Cognitive measure 3
components$cog3 <- define_model_component(
  name = "cog3",
  data = dat,
  outcome = "cog3",
  factor = fm,
  evaluation_indicator = "ind_cog3",
  covariates = c(bkgnd1990,
                 "adveng", "advmath",
                 "HS_track2", "HS_track3", "HS_track4"),
  model_type = "linear",
  loading_normalization = norm_factor1_only
)

# cog4: Cognitive measure 4
components$cog4 <- define_model_component(
  name = "cog4",
  data = dat,
  outcome = "cog4",
  factor = fm,
  evaluation_indicator = "ind_cog4",
  covariates = c(bkgnd1990,
                 "adveng", "advmath",
                 "HS_track2", "HS_track3", "HS_track4"),
  model_type = "linear",
  loading_normalization = norm_factor1_only
)

# -----------------------------------------------------------------------------
# Probit Model for Missing Leadership Data
# Loading normalization: [NA, 0, 0] (only factor 1)
# -----------------------------------------------------------------------------

# misslead: Probability of missing leadership measure (probit)
# Note: indicator="none" in C++ means all observations are used
components$misslead <- define_model_component(
  name = "misslead",
  data = dat,
  outcome = "miss_lead",
  factor = fm,
  evaluation_indicator = NULL,  # All observations evaluated
  covariates = c(bkgnd1990,
                 "adveng", "advmath",
                 "HS_track2", "HS_track3", "HS_track4"),
  model_type = "probit",
  loading_normalization = norm_factor1_only
)

# -----------------------------------------------------------------------------
# Non-cognitive Measures (factors 1 and 2 load)
# -----------------------------------------------------------------------------

norm_factors12 <- c(NA_real_, NA_real_, 0)

# stresstol: Emotional Stability (measured by pf)
components$stresstol <- define_model_component(
  name = "stresstol",
  data = dat,
  outcome = "pf",
  factor = fm,
  evaluation_indicator = "ind_pf",
  covariates = c(bkgnd1990,
                 "adveng", "advmath",
                 "HS_track2", "HS_track3", "HS_track4"),
  model_type = "linear",
  loading_normalization = norm_factors12
)

# lead: Leadership
components$lead <- define_model_component(
  name = "lead",
  data = dat,
  outcome = "lead",
  factor = fm,
  evaluation_indicator = "ind_lead",
  covariates = c(bkgnd1990,
                 "adveng", "advmath",
                 "HS_track2", "HS_track3", "HS_track4"),
  model_type = "linear",
  loading_normalization = norm_factors12
)

# =============================================================================
# Define Model System
# =============================================================================

ms <- define_model_system(
  components = components,
  factor = fm
)

# Print model system summary
print(ms)

# =============================================================================
# Define Estimation Control
# =============================================================================

control <- define_estimation_control(
  n_quad_points = 8,  # 8 quadrature points as specified
  num_cores = 4       # Adjust based on available cores
)

print(control)

# =============================================================================
# Estimate Model
# =============================================================================

cat("\n", rep("=", 70), "\n", sep = "")
cat("Starting Stage 1 Estimation\n")
cat(rep("=", 70), "\n\n", sep = "")

# Estimate using nlminb optimizer (uses analytical gradient + Hessian)
result <- estimate_model_rcpp(
  model_system = ms,
  data = dat,
  control = control,
  optimizer = "nlminb",
  parallel = TRUE,
  verbose = TRUE
)

# =============================================================================
# Print Results
# =============================================================================

cat("\n", rep("=", 70), "\n", sep = "")
cat("Stage 1 Estimation Results\n")
cat(rep("=", 70), "\n\n", sep = "")

print(result)

# =============================================================================
# Save Results
# =============================================================================

output_path <- "/home/Projects/MajorChoice/intdata/MajorChoice/stage1_results.rds"
saveRDS(result, output_path)
cat("\nResults saved to:", output_path, "\n")

# Save estimates in CSV format (similar to meas_par.txt but with names)
params_df <- data.frame(
  index = seq_along(result$estimates) - 1,  # 0-based index
  model = result$param_table$component,
  coefficient = result$param_table$parameter,
  estimate = result$estimates,
  std_error = result$std_errors,
  row.names = NULL
)

params_csv_path <- "/home/Projects/MajorChoice/intdata/MajorChoice/stage1_params.csv"
write.csv(params_df, file = params_csv_path, row.names = FALSE)
cat("Parameter estimates saved to:", params_csv_path, "\n")

# =============================================================================
# Factor Score Estimation (for Stage 2)
# =============================================================================

cat("\n", rep("=", 70), "\n", sep = "")
cat("Estimating Factor Scores\n")
cat(rep("=", 70), "\n\n", sep = "")

# Estimate factor scores for use in stage 2
factor_scores <- estimate_factorscores_rcpp(
  model_system = ms,
  data = dat,
  estimates = result$estimates,
  control = control,
  parallel = TRUE
)

# Save factor scores
fs_output_path <- "/home/Projects/MajorChoice/intdata/MajorChoice/stage1_factorscores.rds"
saveRDS(factor_scores, fs_output_path)
cat("Factor scores saved to:", fs_output_path, "\n")

# Summary of factor scores
cat("\nFactor Score Summary:\n")
cat("Factor 1 (Cognitive):\n")
cat("  Mean:", round(mean(factor_scores$posterior_means[, 1]), 4), "\n")
cat("  SD:  ", round(sd(factor_scores$posterior_means[, 1]), 4), "\n")
cat("  Avg SE:", round(mean(factor_scores$posterior_ses[, 1]), 4), "\n")

cat("Factor 2 (Non-cognitive):\n")
cat("  Mean:", round(mean(factor_scores$posterior_means[, 2]), 4), "\n")
cat("  SD:  ", round(sd(factor_scores$posterior_means[, 2]), 4), "\n")
cat("  Avg SE:", round(mean(factor_scores$posterior_ses[, 2]), 4), "\n")

cat("Factor 3 (General):\n")
cat("  Mean:", round(mean(factor_scores$posterior_means[, 3]), 4), "\n")
cat("  SD:  ", round(sd(factor_scores$posterior_means[, 3]), 4), "\n")
cat("  Avg SE:", round(mean(factor_scores$posterior_ses[, 3]), 4), "\n")
