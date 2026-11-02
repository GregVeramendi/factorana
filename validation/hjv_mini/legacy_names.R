# Parse the final PrintParam block of a legacy run log into (idx, model, var, est, se, param) with R parameter names.
# Legacy prints "****<model> Model****" or "****<model>(Choice c)  Model****" headers and
# " idx. name: est +/- se" rows; loadings are "Fac<k> loading", the linear sigma is "1/Precision", and
# fixed parameters print "+/- -9999.000". Legacy choice c is R alternative c - 1.
legacy_param_table <- function(logfile) {
  log <- readLines(logfile)
  start <- max(grep("^\\*{8}Parameter values:", log))
  model <- NA_character_; alt <- NA_integer_; rows <- list()
  for (ln in log[start:length(log)]) {
    h <- regmatches(ln, regexec("^\\*{4}([^(*]+?) *(\\(Choice ([0-9]+)\\))? +Model\\*+", ln))[[1]]
    if (length(h)) { model <- trimws(h[2]); alt <- if (nzchar(h[4])) as.integer(h[4]) - 1L else NA_integer_; next }
    if (grepl("^\\*{4}Factors", ln)) { model <- "Factors"; alt <- NA_integer_; next }
    m <- regmatches(ln, regexec("^ *([0-9]+)\\. +(.+?): +(-?[0-9.]+(?:e[-+]?[0-9]+)?) \\+/- +(-?[0-9.]+)", ln, perl = TRUE))[[1]]
    if (length(m)) rows[[length(rows) + 1]] <- data.frame(idx = as.integer(m[2]), model = model, alt = alt,
      var = trimws(m[3]), est = as.numeric(m[4]), se = as.numeric(m[5]), stringsAsFactors = FALSE)
  }
  d <- do.call(rbind, rows)
  v <- sub("^Fac([0-9]) loading$", "loading_\\1", d$var); v <- sub("^1/Precision$", "sigma", v)
  # Type intercepts print as "Typ<j> loading" (legacy type j = R type j + 1, type 1 the reference)
  tj <- regmatches(d$var, regexec("^Typ([0-9]+) loading$", d$var)); has <- lengths(tj) > 0
  v[has] <- sprintf("type_%d_intercept", as.integer(vapply(tj[has], `[`, "", 2)) + 1L)
  # Type probability model ("Type Probability" title, legacy choice c = R type c; with two types no choice header)
  istype <- d$model == "Type Probability"; ttyp <- ifelse(is.na(d$alt), 1L, d$alt) + 1L
  tpar <- ifelse(d$var == "intercept", sprintf("typeprob_%d_intercept", ttyp),
                 sprintf("type_%d_loading_%s", ttyp, sub("^Fac([0-9]) loading$", "\\1", d$var)))
  d$param <- ifelse(d$model == "Factors", sub("^Factor([0-9]) sigma$", "factor_var_\\1", d$var),
             ifelse(istype, tpar,
                    ifelse(is.na(d$alt), paste0(d$model, "_", v), paste0(d$model, "_c", d$alt, "_", v))))
  d$fixed <- d$se < 0 | (d$model == "Factors")
  d
}
