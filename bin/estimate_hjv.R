# =============================================================================
# Multi-Stage Estimation Script for Major Choice Project
# Translated from bin/prob_hjv.cc
# =============================================================================
#
# Stage 1: Measurement system (16 linear/probit models)
# Stage 2: Estimate factor scores (using Stage 1 results)
# Stage 3: Education choices (fixed Stage 1 + new education models)
# Stage 4: Labor market outcomes (fixed Stages 1+3 + new labor market models)
#
# 3-factor model with 8 quadrature points, no correlation, single type

# Make sure we find the library in Rscript:
.libPaths(c("/opt/R/4.4.2/lib/R/library", .libPaths()))

library(factorana)
library(haven)

# =============================================================================
# STAGE SELECTION - Modify this to run different stages
# =============================================================================

# Which stage to estimate (1, 2, 3, or 4)
STAGE_TO_RUN <- 1

# Directory for saving/loading results
RESULTS_DIR <- "/home/Projects/MajorChoice/intdata/MajorChoice"

# =============================================================================
# ADAPTIVE QUADRATURE OPTIONS (for Stages 3 and 4 only)
# =============================================================================
# Adaptive quadrature uses factor scores from Stage 2 to reduce
# integration points for observations where factors are well-identified.
# This can dramatically speed up estimation without loss of accuracy.

USE_ADAPTIVE_QUADRATURE <- FALSE  # Set to TRUE to enable for stages 3 and 4
ADAPTIVE_THRESHOLD <- 0.5         # Threshold for determining n_quad per obs (smaller = more points)
                                  # Formula: n_quad = 1 + 2*floor(factor_se/factor_var/threshold)
                                  # Legacy default is 0.5
ADAPTIVE_MAX_QUAD <- 16           # Maximum quadrature points per factor

# File paths for each stage's results
stage_result_files <- list(
  stage1 = file.path(RESULTS_DIR, "hjv_stage1_results.rds"),
  stage2 = file.path(RESULTS_DIR, "hjv_stage2_factorscores.rds"),
  stage3 = file.path(RESULTS_DIR, "hjv_stage3_results.rds"),
  stage4 = file.path(RESULTS_DIR, "hjv_stage4_results.rds")
)

cat("\n", rep("=", 70), "\n", sep = "")
cat("HJV Major Choice Model - Running Stage", STAGE_TO_RUN, "\n")
cat(rep("=", 70), "\n\n", sep = "")

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

# 3-factor model, single type (n_types = 1)
fm <- define_factor_model(
  n_factors = 3,
  n_types = 1  # Single type = no unobserved discrete heterogeneity
)

# =============================================================================
# Loading Normalization Conventions
# =============================================================================

# NA_real_ = free parameter (to be estimated)
# 0        = loading fixed to 0 (factor does not load on this outcome)
# 1        = loading fixed to 1 (scale normalization)

norm_all_free <- c(NA_real_, NA_real_, NA_real_)
norm_factor1_only <- c(NA_real_, 0, 0)
norm_factors12 <- c(NA_real_, NA_real_, 0)

# =============================================================================
# Define Estimation Control
# =============================================================================

control <- define_estimation_control(
  n_quad_points = if (USE_ADAPTIVE_QUADRATURE) ADAPTIVE_MAX_QUAD else 8,  # Use max_quad for adaptive
  num_cores = 32,                                                          # Adjust based on available cores
  adaptive_integration = USE_ADAPTIVE_QUADRATURE,
  adapt_int_thresh = ADAPTIVE_THRESHOLD
)

print(control)

# =============================================================================
# STAGE 1: Measurement System
# =============================================================================

if (STAGE_TO_RUN == 1) {

cat("\n--- Building Stage 1: Measurement System ---\n")

components <- list()

# -----------------------------------------------------------------------------
# 9th Grade Measurements (all 3 factors load freely)
# -----------------------------------------------------------------------------

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
# -----------------------------------------------------------------------------

# misslead: Probability of missing leadership measure (probit)
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

cat("Stage 1 components:", length(components), "\n")

# Define model system for Stage 1 (no previous stage)
ms <- define_model_system(
  components = components,
  factor = fm
)

} # End of Stage 1 definition

# =============================================================================
# STAGE 2: Estimate Factor Scores (using Stage 1 results)
# =============================================================================

if (STAGE_TO_RUN == 2) {

cat("\n--- Stage 2: Estimating Factor Scores ---\n")

# Load Stage 1 results
stage1_file <- stage_result_files$stage1
if (!file.exists(stage1_file)) {
  stop("Stage 1 results not found at: ", stage1_file, "\n",
       "Please run Stage 1 first (set STAGE_TO_RUN <- 1)")
}
cat("Loading Stage 1 results from:", stage1_file, "\n")
stage1_result <- readRDS(stage1_file)
cat("Loaded", length(stage1_result$estimates), "parameters from Stage 1\n\n")

# Use high quadrature for accurate factor scores
cat("Computing factor scores with 16 quadrature points...\n")
fscore_control <- define_estimation_control(n_quad_points = 16)
fscores <- estimate_factorscores_rcpp(
  result = stage1_result,
  data = dat,
  control = fscore_control,
  verbose = TRUE
)

# Extract factor score columns
n_factors <- stage1_result$model_system$factor$n_factors
factor_cols <- paste0("factor_", seq_len(n_factors))
se_cols <- paste0("se_factor_", seq_len(n_factors))

# Create result object with factor scores
result <- list(
  factor_scores = as.matrix(fscores[, factor_cols, drop = FALSE]),
  factor_ses = as.matrix(fscores[, se_cols, drop = FALSE]),
  factor_vars = stage1_result$estimates[grep("^factor_var_", names(stage1_result$estimates))],
  stage1_result = stage1_result,
  n_obs = nrow(fscores),
  converged = fscores$converged
)

cat("\nFactor scores computed for", result$n_obs, "observations\n")
cat("Convergence rate:", sprintf("%.1f%%", 100 * mean(result$converged)), "\n")

# Summary statistics
cat("\nFactor Score Summary:\n")
for (k in seq_len(n_factors)) {
  fs <- result$factor_scores[, k]
  se <- result$factor_ses[, k]
  cat(sprintf("  Factor %d: mean=%.3f, sd=%.3f, mean_se=%.3f\n",
              k, mean(fs), sd(fs), mean(se)))
}

# Save factor scores
output_path <- stage_result_files$stage2
saveRDS(result, output_path)
cat("\nFactor scores saved to:", output_path, "\n")

cat("\n", rep("=", 70), "\n", sep = "")
cat("Stage 2 Complete - Factor scores ready for Stages 3 and 4\n")
cat(rep("=", 70), "\n", sep = "")

# Skip the rest of the estimation code
q(save = "no")

} # End of Stage 2 definition

# =============================================================================
# STAGE 3: Education Choices (with fixed Stage 1)
# =============================================================================

if (STAGE_TO_RUN == 3) {

cat("\n--- Building Stage 3: Education Choices ---\n")

# Load Stage 1 results
stage1_file <- stage_result_files$stage1
if (!file.exists(stage1_file)) {
  stop("Stage 1 results not found at: ", stage1_file, "\n",
       "Please run Stage 1 first (set STAGE_TO_RUN <- 1)")
}
cat("Loading Stage 1 results from:", stage1_file, "\n")
stage1_result <- readRDS(stage1_file)
cat("Loaded", length(stage1_result$estimates), "fixed parameters from Stage 1\n\n")

# Load factor scores from Stage 2
stage2_file <- stage_result_files$stage2
if (!file.exists(stage2_file)) {
  stop("Stage 2 results (factor scores) not found at: ", stage2_file, "\n",
       "Please run Stage 2 first (set STAGE_TO_RUN <- 2)")
}
cat("Loading factor scores from Stage 2:", stage2_file, "\n")
stage2_result <- readRDS(stage2_file)
factor_scores_s1 <- stage2_result$factor_scores
factor_ses_s1 <- stage2_result$factor_ses
factor_vars_s1 <- stage2_result$factor_vars
cat("  Factor scores for", nrow(factor_scores_s1), "observations\n\n")

# Use factor model from previous stage (must be identical object for define_model_system)
fm <- stage1_result$model_system$factor

# Define only Stage 3 components (education choices)
components <- list()

# -----------------------------------------------------------------------------
# 9th Grade Advanced Math (Probit)
# -----------------------------------------------------------------------------

components$advmath <- define_model_component(
  name = "advmath",
  data = dat,
  outcome = "advmath",
  factor = fm,
  evaluation_indicator = "ind_advmath",
  covariates = c(bkgnd1990, "schoolave_advmath", "advmathIV"),
  model_type = "probit",
  loading_normalization = norm_all_free
)

# -----------------------------------------------------------------------------
# 9th Grade Advanced English (Probit)
# -----------------------------------------------------------------------------

components$adveng <- define_model_component(
  name = "adveng",
  data = dat,
  outcome = "adveng",
  factor = fm,
  evaluation_indicator = "ind_adveng",
  covariates = c(bkgnd1990, "schoolave_adveng", "advengIV"),
  model_type = "probit",
  loading_normalization = norm_all_free
)

# -----------------------------------------------------------------------------
# High School Track Choice (Multinomial Logit with 4 choices)
# -----------------------------------------------------------------------------

components$HS_track <- define_model_component(
  name = "HS_track",
  data = dat,
  outcome = "HS_track",
  factor = fm,
  evaluation_indicator = NULL,
  covariates = c(bkgnd1990,
                 "adveng", "advmath",
                 "sgpa92", "sgpa93", "sgpa94",
                 "schoolave_hstrk2", "hstrk2IV",
                 "schoolave_hstrk3", "hstrk3IV",
                 "schoolave_hstrk4", "hstrk4IV"),
  model_type = "logit",
  num_choices = 4,
  loading_normalization = norm_all_free
)

# Apply fix_coefficient constraints
components$HS_track <- fix_coefficient(components$HS_track, "schoolave_hstrk4", 0.0, choice = 1)
components$HS_track <- fix_coefficient(components$HS_track, "hstrk4IV", 0.0, choice = 1)
components$HS_track <- fix_coefficient(components$HS_track, "schoolave_hstrk3", 0.0, choice = 1)
components$HS_track <- fix_coefficient(components$HS_track, "hstrk3IV", 0.0, choice = 1)

components$HS_track <- fix_coefficient(components$HS_track, "schoolave_hstrk4", 0.0, choice = 2)
components$HS_track <- fix_coefficient(components$HS_track, "hstrk4IV", 0.0, choice = 2)
components$HS_track <- fix_coefficient(components$HS_track, "schoolave_hstrk2", 0.0, choice = 2)
components$HS_track <- fix_coefficient(components$HS_track, "hstrk2IV", 0.0, choice = 2)

components$HS_track <- fix_coefficient(components$HS_track, "schoolave_hstrk3", 0.0, choice = 3)
components$HS_track <- fix_coefficient(components$HS_track, "hstrk3IV", 0.0, choice = 3)
components$HS_track <- fix_coefficient(components$HS_track, "schoolave_hstrk2", 0.0, choice = 3)
components$HS_track <- fix_coefficient(components$HS_track, "hstrk2IV", 0.0, choice = 3)

# -----------------------------------------------------------------------------
# Apply to College (Probit)
# -----------------------------------------------------------------------------

components$educDapply <- define_model_component(
  name = "educDapply",
  data = dat,
  outcome = "applyCollege",
  factor = fm,
  evaluation_indicator = "hsgrad",
  covariates = c(bkgnd1990,
                 "adveng", "advmath",
                 "HS_track3", "HS_track4",
                 "gpahs2", "gpahs3", "gpahs4"),
  model_type = "probit",
  loading_normalization = norm_all_free
)

# -----------------------------------------------------------------------------
# Take SweSAT (Probit)
# -----------------------------------------------------------------------------

components$Dswesat <- define_model_component(
  name = "Dswesat",
  data = dat,
  outcome = "Dswesat",
  factor = fm,
  evaluation_indicator = "applyCollege",
  covariates = c(bkgnd1990,
                 "adveng", "advmath",
                 "HS_track3", "HS_track4",
                 "gpahs2", "gpahs3", "gpahs4"),
  model_type = "probit",
  loading_normalization = norm_all_free
)

# -----------------------------------------------------------------------------
# Total SweSAT Score (Linear)
# -----------------------------------------------------------------------------

components$swesat <- define_model_component(
  name = "swesat",
  data = dat,
  outcome = "swesat_total",
  factor = fm,
  evaluation_indicator = "Dswesat",
  covariates = c(bkgnd1990,
                 "adveng", "advmath",
                 "HS_track3", "HS_track4"),
  model_type = "linear",
  loading_normalization = norm_all_free
)

# -----------------------------------------------------------------------------
# Application Major (Exploded Nested Logit with 13 choices, 13 ranks)
# -----------------------------------------------------------------------------

# Build covariate lists
# Note: Major 1 is omitted category, so GPA vars are for majors 2-13
# GPA bins are 1-4 (quintiles)
enrollmaj_gpa <- character(0)
for (imaj in 2:13) {
  for (ibin in 1:4) {
    enrollmaj_gpa <- c(enrollmaj_gpa, paste0("gpahs_maj", imaj, "_", ibin))
  }
}

switchmaj_gpa <- character(0)
for (imaj in 2:12) {
  for (ibin in 1:4) {
    switchmaj_gpa <- c(switchmaj_gpa, paste0("gpahs_maj", imaj, "_", ibin))
  }
}

enrollmaj_iv <- c(
  paste0("maj", 2:13, "IV"),
  paste0("schoolmean_m1_maj", 2:13),
  paste0("logadjadmitshare_m", 2:13),
  "min_logdist_comb"
)

enrollmaj_indic <- paste0("Denroll", 2:13)

major_apply_outcomes <- paste0("major_apply", 1:13)

components$educDapplication <- define_model_component(
  name = "educDapplication",
  data = dat,
  outcome = major_apply_outcomes,
  factor = fm,
  evaluation_indicator = "ind_major_apply_pm",
  covariates = c(bkgnd1990,
                 "adveng", "advmath",
                 "HS_track3", "HS_track4",
                 enrollmaj_gpa,
                 "swesat_total", "swesat_total_miss",
                 enrollmaj_iv),
  model_type = "logit",
  num_choices = 13,
  exclude_chosen = FALSE,
  loading_normalization = norm_all_free,
  skip_collinearity_check = TRUE  # Many coefficients fixed via fix_coefficient() below
)

# Apply coefficient constraints for educDapplication
# For each choice (2-13), fix GPA, IV, school mean, and logadmitshare to 0
# for majors that don't match the choice
# Note: GPA vars are for majors 2-13 (maj1 omitted), bins 1-4
for (ichoice in 2:13) {
  # Fix GPA vars for non-matching majors (majors 2-13)
  for (imaj in 2:13) {
    if (ichoice != imaj) {
      for (ibin in 1:4) {
        gpavar <- paste0("gpahs_maj", imaj, "_", ibin)
        components$educDapplication <- fix_coefficient(
          components$educDapplication, gpavar, 0.0, choice = ichoice - 1
        )
      }
    }
  }
  # Fix IV, school mean, logadmitshare for non-matching majors (2-13)
  for (imaj in 2:13) {
    if (ichoice != imaj) {
      ivvar <- paste0("maj", imaj, "IV")
      components$educDapplication <- fix_coefficient(
        components$educDapplication, ivvar, 0.0, choice = ichoice - 1
      )
      schave <- paste0("schoolmean_m1_maj", imaj)
      components$educDapplication <- fix_coefficient(
        components$educDapplication, schave, 0.0, choice = ichoice - 1
      )
      logadmitshare <- paste0("logadjadmitshare_m", imaj)
      components$educDapplication <- fix_coefficient(
        components$educDapplication, logadmitshare, 0.0, choice = ichoice - 1
      )
    }
  }
}

# -----------------------------------------------------------------------------
# Enroll after Admission (Probit)
# -----------------------------------------------------------------------------

components$educDadmitenroll <- define_model_component(
  name = "educDadmitenroll",
  data = dat,
  outcome = "enroll",
  factor = fm,
  evaluation_indicator = "ind_admitenroll",
  covariates = c(bkgnd1990,
                 "adveng", "advmath",
                 "HS_track3", "HS_track4"),
  model_type = "probit",
  loading_normalization = norm_all_free
)

# -----------------------------------------------------------------------------
# Switch to Major (Multinomial Logit with 12 choices)
# -----------------------------------------------------------------------------

components$educDSwitchtoMajor <- define_model_component(
  name = "educDSwitchtoMajor",
  data = dat,
  outcome = "DSwitchMajor",
  factor = fm,
  evaluation_indicator = "ind_DSwitchMajor",
  covariates = c(bkgnd1990,
                 "adveng", "advmath",
                 "HS_track3", "HS_track4",
                 switchmaj_gpa,
                 "swesat_total", "swesat_total_miss",
                 enrollmaj_indic),
  model_type = "logit",
  num_choices = 12,
  loading_normalization = norm_all_free,
  skip_collinearity_check = TRUE  # Many coefficients fixed via fix_coefficient() below
)

# Apply GPA coefficient constraints for educDSwitchtoMajor
# Note: GPA vars are for majors 2-12 only (maj1 omitted), bins 1-4
for (ichoice in 2:12) {
  for (imaj in 2:12) {
    if (ichoice != imaj) {
      for (ibin in 1:4) {
        gpavar <- paste0("gpahs_maj", imaj, "_", ibin)
        components$educDSwitchtoMajor <- fix_coefficient(
          components$educDSwitchtoMajor, gpavar, 0.0, choice = ichoice - 1
        )
      }
    }
  }
}

# Apply Denroll coefficient constraints
components$educDSwitchtoMajor <- fix_coefficient(components$educDSwitchtoMajor, "Denroll7", 0.0, choice = 1)
components$educDSwitchtoMajor <- fix_coefficient(components$educDSwitchtoMajor, "Denroll11", 0.0, choice = 1)
components$educDSwitchtoMajor <- fix_coefficient(components$educDSwitchtoMajor, "Denroll13", 0.0, choice = 1)

components$educDSwitchtoMajor <- fix_coefficient(components$educDSwitchtoMajor, "Denroll2", 0.0, choice = 2)
components$educDSwitchtoMajor <- fix_coefficient(components$educDSwitchtoMajor, "Denroll5", 0.0, choice = 2)
components$educDSwitchtoMajor <- fix_coefficient(components$educDSwitchtoMajor, "Denroll6", 0.0, choice = 2)
components$educDSwitchtoMajor <- fix_coefficient(components$educDSwitchtoMajor, "Denroll7", 0.0, choice = 2)
components$educDSwitchtoMajor <- fix_coefficient(components$educDSwitchtoMajor, "Denroll8", 0.0, choice = 2)
components$educDSwitchtoMajor <- fix_coefficient(components$educDSwitchtoMajor, "Denroll10", 0.0, choice = 2)
components$educDSwitchtoMajor <- fix_coefficient(components$educDSwitchtoMajor, "Denroll11", 0.0, choice = 2)
components$educDSwitchtoMajor <- fix_coefficient(components$educDSwitchtoMajor, "Denroll13", 0.0, choice = 2)

components$educDSwitchtoMajor <- fix_coefficient(components$educDSwitchtoMajor, "Denroll4", 0.0, choice = 3)
components$educDSwitchtoMajor <- fix_coefficient(components$educDSwitchtoMajor, "Denroll7", 0.0, choice = 3)
components$educDSwitchtoMajor <- fix_coefficient(components$educDSwitchtoMajor, "Denroll13", 0.0, choice = 3)

components$educDSwitchtoMajor <- fix_coefficient(components$educDSwitchtoMajor, "Denroll11", 0.0, choice = 4)
components$educDSwitchtoMajor <- fix_coefficient(components$educDSwitchtoMajor, "Denroll13", 0.0, choice = 4)

components$educDSwitchtoMajor <- fix_coefficient(components$educDSwitchtoMajor, "Denroll4", 0.0, choice = 5)
components$educDSwitchtoMajor <- fix_coefficient(components$educDSwitchtoMajor, "Denroll5", 0.0, choice = 5)
components$educDSwitchtoMajor <- fix_coefficient(components$educDSwitchtoMajor, "Denroll11", 0.0, choice = 5)
components$educDSwitchtoMajor <- fix_coefficient(components$educDSwitchtoMajor, "Denroll12", 0.0, choice = 5)
components$educDSwitchtoMajor <- fix_coefficient(components$educDSwitchtoMajor, "Denroll13", 0.0, choice = 5)

components$educDSwitchtoMajor <- fix_coefficient(components$educDSwitchtoMajor, "Denroll5", 0.0, choice = 6)
components$educDSwitchtoMajor <- fix_coefficient(components$educDSwitchtoMajor, "Denroll11", 0.0, choice = 6)

components$educDSwitchtoMajor <- fix_coefficient(components$educDSwitchtoMajor, "Denroll2", 0.0, choice = 7)
components$educDSwitchtoMajor <- fix_coefficient(components$educDSwitchtoMajor, "Denroll5", 0.0, choice = 7)
components$educDSwitchtoMajor <- fix_coefficient(components$educDSwitchtoMajor, "Denroll11", 0.0, choice = 7)

components$educDSwitchtoMajor <- fix_coefficient(components$educDSwitchtoMajor, "Denroll2", 0.0, choice = 8)
components$educDSwitchtoMajor <- fix_coefficient(components$educDSwitchtoMajor, "Denroll4", 0.0, choice = 8)
components$educDSwitchtoMajor <- fix_coefficient(components$educDSwitchtoMajor, "Denroll5", 0.0, choice = 8)
components$educDSwitchtoMajor <- fix_coefficient(components$educDSwitchtoMajor, "Denroll8", 0.0, choice = 8)
components$educDSwitchtoMajor <- fix_coefficient(components$educDSwitchtoMajor, "Denroll11", 0.0, choice = 8)
components$educDSwitchtoMajor <- fix_coefficient(components$educDSwitchtoMajor, "Denroll13", 0.0, choice = 8)

components$educDSwitchtoMajor <- fix_coefficient(components$educDSwitchtoMajor, "Denroll2", 0.0, choice = 9)
components$educDSwitchtoMajor <- fix_coefficient(components$educDSwitchtoMajor, "Denroll3", 0.0, choice = 9)
components$educDSwitchtoMajor <- fix_coefficient(components$educDSwitchtoMajor, "Denroll4", 0.0, choice = 9)
components$educDSwitchtoMajor <- fix_coefficient(components$educDSwitchtoMajor, "Denroll6", 0.0, choice = 9)
components$educDSwitchtoMajor <- fix_coefficient(components$educDSwitchtoMajor, "Denroll7", 0.0, choice = 9)
components$educDSwitchtoMajor <- fix_coefficient(components$educDSwitchtoMajor, "Denroll8", 0.0, choice = 9)
components$educDSwitchtoMajor <- fix_coefficient(components$educDSwitchtoMajor, "Denroll12", 0.0, choice = 9)
components$educDSwitchtoMajor <- fix_coefficient(components$educDSwitchtoMajor, "Denroll13", 0.0, choice = 9)

components$educDSwitchtoMajor <- fix_coefficient(components$educDSwitchtoMajor, "Denroll5", 0.0, choice = 10)
components$educDSwitchtoMajor <- fix_coefficient(components$educDSwitchtoMajor, "Denroll6", 0.0, choice = 10)
components$educDSwitchtoMajor <- fix_coefficient(components$educDSwitchtoMajor, "Denroll11", 0.0, choice = 10)

components$educDSwitchtoMajor <- fix_coefficient(components$educDSwitchtoMajor, "Denroll3", 0.0, choice = 11)
components$educDSwitchtoMajor <- fix_coefficient(components$educDSwitchtoMajor, "Denroll4", 0.0, choice = 11)
components$educDSwitchtoMajor <- fix_coefficient(components$educDSwitchtoMajor, "Denroll5", 0.0, choice = 11)
components$educDSwitchtoMajor <- fix_coefficient(components$educDSwitchtoMajor, "Denroll6", 0.0, choice = 11)
components$educDSwitchtoMajor <- fix_coefficient(components$educDSwitchtoMajor, "Denroll7", 0.0, choice = 11)
components$educDSwitchtoMajor <- fix_coefficient(components$educDSwitchtoMajor, "Denroll11", 0.0, choice = 11)

# -----------------------------------------------------------------------------
# Graduate by Major (12 Probit models)
# -----------------------------------------------------------------------------

for (imajor in 1:12) {
  comp_name <- paste0("educDGradCollMaj", imajor)
  outcome_var <- paste0("DGradColl", imajor)
  ind_var <- paste0("ind_DGradColl", imajor)

  components[[comp_name]] <- define_model_component(
    name = comp_name,
    data = dat,
    outcome = outcome_var,
    factor = fm,
    evaluation_indicator = ind_var,
    covariates = c(bkgnd1990,
                   "adveng", "advmath",
                   "HS_track3", "HS_track4"),
    model_type = "probit",
    loading_normalization = norm_all_free
  )
}

cat("Stage 3 new components:", length(components), "\n")

# Define model system with previous_stage (Stage 1 fixed)
ms <- define_model_system(
  components = components,
  factor = fm,
  previous_stage = stage1_result
)

} # End of Stage 3 definition

# =============================================================================
# STAGE 4: Labor Market Outcomes (with fixed Stages 1+3)
# =============================================================================

if (STAGE_TO_RUN == 4) {

cat("\n--- Building Stage 4: Labor Market Outcomes ---\n")

# Load Stage 3 results (which includes Stage 1 as previous_stage)
stage3_file <- stage_result_files$stage3
if (!file.exists(stage3_file)) {
  stop("Stage 3 results not found at: ", stage3_file, "\n",
       "Please run Stage 3 first (set STAGE_TO_RUN <- 3)")
}
cat("Loading Stage 3 results from:", stage3_file, "\n")
stage3_result <- readRDS(stage3_file)
cat("Loaded", length(stage3_result$estimates), "fixed parameters from Stages 1+3\n\n")

# Load factor scores from Stage 2
stage2_file <- stage_result_files$stage2
if (!file.exists(stage2_file)) {
  stop("Stage 2 results (factor scores) not found at: ", stage2_file, "\n",
       "Please run Stage 2 first (set STAGE_TO_RUN <- 2)")
}
cat("Loading factor scores from Stage 2:", stage2_file, "\n")
stage2_result <- readRDS(stage2_file)
factor_scores_s1 <- stage2_result$factor_scores
factor_ses_s1 <- stage2_result$factor_ses
factor_vars_s1 <- stage2_result$factor_vars
cat("  Factor scores for", nrow(factor_scores_s1), "observations\n\n")

# Use factor model from previous stage (must be identical object for define_model_system)
fm <- stage3_result$model_system$factor

# Define only Stage 4 components (labor market outcomes)
components <- list()

# Income measures
income_meas <- c("logpvdispinc", "logwage")

# Education categories
education <- c("all", "HSpDO", "HStrk1", "HStrk2", "HStrk3", "HStrk4",
               "CollDO_low", "CollDO_high", "CollGrad",
               "Maj1", "Maj2", "Maj3", "Maj4", "Maj5", "Maj6",
               "Maj7", "Maj8", "Maj9", "Maj10", "Maj11", "Maj12")

# Loop over income measures
for (imeas in 1:2) {

  cat("\n--- Adding", income_meas[imeas], "models ---\n")

  # Loop over education categories
  # Skip ieduc=1 (HSpDO) and ieduc=8 (CollGrad)
  for (ieduc in 2:20) {

    # Skip CollGrad (legacy ieduc=8, R index 9)
    if (ieduc == 9) next

    # Component naming (matching legacy)
    educ_name <- education[ieduc + 1]  # R is 1-indexed
    outcome_var <- income_meas[imeas]
    ind_var <- paste0("ind_", income_meas[imeas], "_", educ_name)
    comp_name <- paste0(income_meas[imeas], educ_name)

    # Build covariate list
    modelvars <- c(bkgnd1990, "adveng", "advmath")

    # Add HS_track for non-HS-track-specific models (ieduc > 5)
    if (ieduc > 5) {
      modelvars <- c(modelvars, "HS_track3", "HS_track4")
    }

    # Create model component
    components[[comp_name]] <- define_model_component(
      name = comp_name,
      data = dat,
      outcome = outcome_var,
      factor = fm,
      evaluation_indicator = ind_var,
      covariates = modelvars,
      model_type = "linear",
      loading_normalization = norm_all_free
    )
  }
}

cat("Stage 4 new components:", length(components), "\n")

# Define model system with previous_stage (Stages 1+3 fixed)
ms <- define_model_system(
  components = components,
  factor = fm,
  previous_stage = stage3_result
)

} # End of Stage 4 definition

# =============================================================================
# Print Model System Summary
# =============================================================================

cat("\nTotal components in model system:", length(ms$components), "\n")
if (!is.null(ms$previous_stage_info)) {
  cat("  - Fixed from previous stage(s):", ms$previous_stage_info$n_components, "\n")
  cat("  - New components to estimate:", length(ms$components) - ms$previous_stage_info$n_components, "\n")
}
print(ms)

# =============================================================================
# Estimate Model
# =============================================================================

cat("\n", rep("=", 70), "\n", sep = "")
cat("Starting Stage", STAGE_TO_RUN, "Estimation\n")
cat(rep("=", 70), "\n\n", sep = "")

# Estimate using nlminb optimizer (uses analytical gradient + Hessian)
# For Stages 3 and 4, use factor scores for initialization and optionally for adaptive quadrature
if (STAGE_TO_RUN > 2 && !is.null(factor_scores_s1)) {
  if (USE_ADAPTIVE_QUADRATURE) {
    # Use factor scores for BOTH initialization AND adaptive quadrature
    result <- estimate_model_rcpp(
      model_system = ms,
      data = dat,
      control = control,
      optimizer = "nlminb",
      parallel = TRUE,
      verbose = TRUE,
      factor_scores = factor_scores_s1,
      factor_ses = factor_ses_s1,
      factor_vars = factor_vars_s1,
      init_factor_scores = factor_scores_s1  # Use factor scores for loading initialization
    )
  } else {
    # Use factor scores only for initialization (not adaptive quadrature)
    result <- estimate_model_rcpp(
      model_system = ms,
      data = dat,
      control = control,
      optimizer = "nlminb",
      parallel = TRUE,
      verbose = TRUE,
      init_factor_scores = factor_scores_s1  # Use factor scores for loading initialization
    )
  }
} else {
  # Stage 1 or no factor scores available: standard estimation
  result <- estimate_model_rcpp(
    model_system = ms,
    data = dat,
    control = control,
    optimizer = "nlminb",
    parallel = TRUE,
    verbose = TRUE
  )
}

# =============================================================================
# Print Results
# =============================================================================

cat("\n", rep("=", 70), "\n", sep = "")
cat("Stage", STAGE_TO_RUN, "Estimation Results\n")
cat(rep("=", 70), "\n\n", sep = "")

print(result)

# =============================================================================
# Save Results
# =============================================================================

# Save to stage-specific file
output_path <- stage_result_files[[paste0("stage", STAGE_TO_RUN)]]
saveRDS(result, output_path)
cat("\nResults saved to:", output_path, "\n")

# Save estimates in CSV format
parse_param_name <- function(name) {
  if (grepl("^factor_var_", name)) {
    return(c(model = "factor", coefficient = name))
  }
  if (grepl("^factor_corr_", name)) {
    return(c(model = "factor", coefficient = name))
  }
  parts <- strsplit(name, "_")[[1]]
  if (length(parts) >= 2) {
    suffixes <- c("intercept", "loading", "sigma", "thresh", "beta")
    for (i in seq_along(parts)) {
      if (parts[i] %in% suffixes || grepl("^loading$", parts[i])) {
        model_name <- paste(parts[1:(i-1)], collapse = "_")
        coef_name <- paste(parts[i:length(parts)], collapse = "_")
        return(c(model = model_name, coefficient = coef_name))
      }
    }
  }
  return(c(model = name, coefficient = name))
}

parsed <- t(sapply(result$param_names, parse_param_name))
params_df <- data.frame(
  index = seq_along(result$estimates) - 1,
  model = parsed[, "model"],
  coefficient = parsed[, "coefficient"],
  estimate = result$estimates,
  std_error = result$std_errors,
  row.names = NULL
)

params_csv_path <- file.path(RESULTS_DIR, paste0("hjv_stage", STAGE_TO_RUN, "_params.csv"))
write.csv(params_df, file = params_csv_path, row.names = FALSE)
cat("Parameter estimates saved to:", params_csv_path, "\n")

cat("\nStage", STAGE_TO_RUN, "estimation complete!\n")
if (STAGE_TO_RUN < 4) {
  cat("\nTo run the next stage, set STAGE_TO_RUN <-", STAGE_TO_RUN + 1, "and re-run this script.\n")
}
