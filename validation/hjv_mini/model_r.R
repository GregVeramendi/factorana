# Mini-HJV model systems for the R package. Mirrors prob_hjv_mini.cc component by component.
# build_mini_system(dat, stage, previous) returns a model_system for stage 1 (measurement),
# 3 (education, with stage 1 fixed) or 4 (labor market, with stages 1 and 3 fixed).
# ntypes > 1 (T9): latent types with type-specific intercepts in every education and labor component
# (use_types = TRUE), none in the measurement system; stage 5 = education and labor market jointly (production layout).
mini_factor_model <- function(ntypes = 1) {
  fix_factor_param(define_factor_model(n_factors = 3, n_types = ntypes), c(factor_var_1 = 1, factor_var_2 = 1, factor_var_3 = 1))
}
build_mini_system <- function(dat, stage = 1, previous = NULL, ntypes = 1) {
  fm <- mini_factor_model(ntypes)
  free3 <- c(NA_real_, NA_real_, NA_real_); f1 <- c(NA_real_, 0, 0); f12 <- c(NA_real_, NA_real_, 0)
  X <- c("intercept", "x1", "x2"); HS <- c("HS_track2", "HS_track3", "HS_track4")
  Z <- c("advmath", "adveng", "HS_track3", "HS_track4")
  comp <- function(name, outcome = name, covariates = X, type = "linear", norm = free3, ind = NULL, types = FALSE, ...) {
    suppressWarnings(define_model_component(name, dat, outcome, fm, covariates = covariates, model_type = type,
                                            loading_normalization = norm, evaluation_indicator = ind,
                                            use_types = types && ntypes > 1, ...))
  }
  meas <- list(
    comp("g_eng9"), comp("g_math9"), comp("g_swe9"), comp("g_sports9"), comp("g_gpa9"),
    comp("g_swe10", covariates = c(X, HS), ind = "ind_hs"),
    comp("g_sports10", covariates = c(X, HS), ind = "ind_hs"),
    comp("g_gpahs", covariates = c(X, HS), ind = "ind_hs"),
    comp("cog1", norm = f1, ind = "ind_cog"), comp("cog2", norm = f1, ind = "ind_cog"), comp("cog3", norm = f1, ind = "ind_cog"),
    comp("misslead", type = "probit", norm = f1),
    comp("lead", norm = f12, ind = "ind_lead"), comp("stresstol", norm = f12, ind = "ind_lead"))
  if (stage == 1) return(define_model_system(components = meas, factor = fm))

  educ <- list(
    comp("advmath", type = "probit", types = TRUE), comp("adveng", type = "probit", types = TRUE),
    # HS track: 4 choices. x2 coefficient of the third non-reference alternative (legacy choice 4) fixed at 0.
    fix_coefficient(comp("HS_track", covariates = c(X, "advmath", "adveng"), type = "logit", num_choices = 4, types = TRUE),
                    "x2", 0, choice = 3),
    comp("apply", covariates = c(X, Z), type = "probit", types = TRUE),
    comp("application", outcome = paste0("app", 1:5), covariates = c(X, Z), type = "logit", num_choices = 5,
         exclude_chosen = FALSE, ind = "apply", types = TRUE),
    comp("Dswesat", covariates = c(X, Z), type = "probit", types = TRUE),
    comp("swesat", covariates = c(X, Z), ind = "Dswesat", types = TRUE),
    comp("enroll", covariates = c(X, Z), type = "probit", ind = "apply", types = TRUE),
    comp("switch_major", covariates = c(X, Z), type = "logit", num_choices = 4, ind = "ind_enroll", types = TRUE),
    comp("grad1", covariates = c(X, Z), type = "probit", ind = "ind_enroll", types = TRUE),
    comp("grad2", covariates = c(X, Z), type = "probit", ind = "ind_enroll", types = TRUE),
    comp("grad3", covariates = c(X, Z), type = "probit", ind = "ind_enroll", types = TRUE))
  if (stage == 3) return(define_model_system(components = educ, factor = fm, previous_stage = previous))

  labor <- list(comp("wage1", covariates = c(X, "HS_track3", "HS_track4"), types = TRUE),
                comp("wage2", covariates = c(X, "HS_track3", "HS_track4"), types = TRUE))
  if (stage == 5) return(define_model_system(components = c(educ, labor), factor = fm, previous_stage = previous))
  define_model_system(components = labor, factor = fm, previous_stage = previous)
}
