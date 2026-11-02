# Map the generator's truth list to the R package's parameter names.
# truth_r_params(truth, param_names) returns a named vector over param_names; it stops on any unmapped name.
truth_r_params <- function(truth, param_names, subset = FALSE) {
  X <- c("intercept", "x1", "x2"); HS <- c("HS_track2", "HS_track3", "HS_track4")
  Z <- c("advmath", "adveng", "HS_track3", "HS_track4"); W <- c("HS_track3", "HS_track4")
  v <- c(factor_var_1 = 1, factor_var_2 = 1, factor_var_3 = 1)
  add <- function(nm, val) { v[nm] <<- val }
  # latent types (T9): typeprob_<t>_intercept, type_<t>_loading_<k>, and per-component type intercepts
  if (!is.null(truth$types)) for (t in seq_along(truth$types$a)) {
    add(sprintf("typeprob_%d_intercept", t + 1), truth$types$a[t]); add(sprintf("type_%d_loading_%d", t + 1, 1:3), truth$types$L[t, ]) }
  single <- function(c, t, extra = NULL, extra_names = NULL, sigma = FALSE) {
    add(paste0(c, "_", X), t$b)
    if (!is.null(extra)) add(paste0(c, "_", extra_names), extra)
    for (k in 1:3) if (t$lam[k] != 0 || c %in% c("g_eng9")) add(paste0(c, "_loading_", k), t$lam[k])
    if (sigma) add(paste0(c, "_sigma"), t$sig)
    if (length(t$d)) add(sprintf("%s_type_%d_intercept", c, seq_along(t$d) + 1), t$d)
  }
  for (c in c("g_eng9", "g_math9", "g_swe9", "g_sports9", "g_gpa9")) single(c, truth[[c]], sigma = TRUE)
  for (c in c("g_swe10", "g_sports10", "g_gpahs")) single(c, truth[[c]], truth[[c]]$g_hs, HS, sigma = TRUE)
  for (c in c("cog1", "cog2", "cog3")) single(c, truth[[c]], sigma = TRUE)
  single("misslead", truth$misslead); single("lead", truth$lead, sigma = TRUE); single("stresstol", truth$stresstol, sigma = TRUE)
  single("advmath", truth$advmath); single("adveng", truth$adveng)
  multi <- function(c, t, covs, G = NULL) {
    for (a in seq_len(nrow(t$B))) {
      add(paste0(c, "_c", a, "_", X), t$B[a, ])
      if (!is.null(G)) add(paste0(c, "_c", a, "_", covs), G[a, ])
      add(paste0(c, "_c", a, "_loading_", 1:3), t$L[a, ])
      if (!is.null(t$D) && ncol(t$D)) add(sprintf("%s_c%d_type_%d_intercept", c, a, seq_len(ncol(t$D)) + 1), t$D[a, ])
    }
  }
  multi("HS_track", truth$HS_track, c("advmath", "adveng"), truth$HS_track$g_adv)
  single("apply", truth$apply, truth$apply$g, Z)
  multi("application", truth$app, Z, truth$app$G)
  single("Dswesat", truth$Dswesat, truth$Dswesat$g, Z)
  single("swesat", truth$swesat, truth$swesat$g, Z, sigma = TRUE)
  single("enroll", truth$enroll, truth$enroll$g, Z)
  multi("switch_major", truth$switch, Z, truth$switch$G)
  for (j in 1:3) single(paste0("grad", j), list(b = truth$grad$b[j, ], lam = truth$grad$lam[j, ],
                                                 d = if (!is.null(truth$grad$d)) truth$grad$d[j, ] else NULL), truth$grad$g[j, ], Z)
  single("wage1", truth$wage1, truth$wage1$g, W, sigma = TRUE); single("wage2", truth$wage2, truth$wage2$g, W, sigma = TRUE)
  missing <- setdiff(param_names, names(v)); extra <- setdiff(names(v), param_names)
  if (length(missing)) stop("unmapped R parameters: ", paste(missing, collapse = ", "))
  if (length(extra) && !subset) stop("truth names not in model: ", paste(extra, collapse = ", "))
  v[param_names]
}
