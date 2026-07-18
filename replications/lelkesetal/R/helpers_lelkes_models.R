lelkes_iv_formula <- function() {
  affective_polarization ~ log_providers +
    as.factor(year) + region + percent_black + percent_white +
    percent_male + lowed + unemploymentrate + density + log_HHINC |
    log_Total + as.factor(year) + region + percent_black + percent_white +
    percent_male + lowed + unemploymentrate + density + log_HHINC
}

lelkes_extract_term <- function(fit, term = "log_providers") {
  tab <- as.data.frame(fit$estimates)
  if (!"term" %in% names(tab)) {
    stop("`fit$estimates` does not contain a `term` column.", call. = FALSE)
  }
  row <- tab[tab$term == term, , drop = FALSE]
  if (!nrow(row)) {
    stop("Could not find term `", term, "` in fitted SPLIV object.", call. = FALSE)
  }
  row[1, , drop = FALSE]
}

lelkes_baseline_fit <- function(dat) {
  spliv(
    formula = lelkes_iv_formula(),
    data = dat,
    vcov = "cluster",
    cluster = ~ state,
    method = "ltz",
    delta = 0,
    scale_instrument = "residual_sd"
  )
}

lelkes_controls_formula <- function() {
  ~ as.factor(year) + region + percent_black + percent_white +
    percent_male + lowed + unemploymentrate + density + log_HHINC
}

lelkes_residual_sd <- function(response, dat) {
  fml <- stats::as.formula(paste(response, paste(deparse(lelkes_controls_formula()), collapse = "")))
  fit <- stats::lm(fml, data = dat)
  stats::sd(stats::residuals(fit))
}

lelkes_cluster_vcov_lm <- function(fit, cluster) {
  X <- stats::model.matrix(fit)
  u <- stats::residuals(fit)
  cluster <- as.factor(cluster)
  if (length(cluster) != length(u)) {
    stop("Cluster vector length must match fitted residual length.", call. = FALSE)
  }
  keep <- stats::complete.cases(X, u, cluster)
  X <- X[keep, , drop = FALSE]
  u <- u[keep]
  cluster <- droplevels(cluster[keep])
  coef_names <- colnames(X)
  G <- nlevels(cluster)
  n <- nrow(X)
  k <- ncol(X)
  if (G < 2) {
    stop("Need at least two clusters for clustered first-stage diagnostics.", call. = FALSE)
  }
  bread <- qr.solve(crossprod(X))
  meat <- matrix(0, nrow = k, ncol = k)
  for (g in levels(cluster)) {
    idx <- which(cluster == g)
    xu <- crossprod(X[idx, , drop = FALSE], u[idx])
    meat <- meat + xu %*% t(xu)
  }
  meat <- (G / (G - 1)) * ((n - 1) / max(1, n - k)) * meat
  out <- bread %*% meat %*% bread
  dimnames(out) <- list(coef_names, coef_names)
  out
}

lelkes_first_stage_full <- function(dat) {
  fit <- stats::lm(
    log_providers ~ as.factor(year) + region + percent_black + percent_white +
      percent_male + lowed + unemploymentrate + density + log_HHINC + log_Total,
    data = dat
  )
  V <- lelkes_cluster_vcov_lm(fit, dat$state)
  coef_name <- "log_Total"
  beta <- stats::coef(fit)[[coef_name]]
  se <- sqrt(pmax(0, V[coef_name, coef_name]))
  c(
    coefficient = as.numeric(beta),
    se = as.numeric(se),
    f_statistic = as.numeric((beta / se)^2)
  )
}

lelkes_residual_metrics <- function(dat, fit) {
  y_sd <- lelkes_residual_sd("affective_polarization", dat)
  x_sd <- lelkes_residual_sd("log_providers", dat)
  z_sd <- as.numeric(fit$residualized_instrument_sd[["log_Total"]])
  list(
    residual_y_sd = y_sd,
    residual_x_sd = x_sd,
    residual_z_sd = z_sd,
    nobs = length(fit$internals$y),
    n_clusters = length(unique(as.character(fit$internals$cluster_id))),
    first_stage_full = lelkes_first_stage_full(dat)
  )
}

lelkes_run_sensitivity_path <- function(dat,
                                        pattern,
                                        method = c("uci", "ltz"),
                                        delta_grid_y,
                                        uci_steps = 21L) {
  method <- match.arg(method)
  extra <- if (identical(method, "uci")) {
    list(grid = list(steps = as.integer(uci_steps), level = 0.95))
  } else {
    list(grid = list(level = 0.95))
  }
  path <- do.call(
    spliv_sensitivity_path,
    c(
      list(
        formula = lelkes_iv_formula(),
        data = dat,
        vcov = "cluster",
        cluster = ~ state,
        method = method,
        delta_grid = delta_grid_y,
        violation_pattern = pattern,
        scale_instrument = "residual_sd"
      ),
      extra
    )
  )
  path[path$term == "log_providers", , drop = FALSE]
}

lelkes_run_bpe <- function(dat,
                           design_info,
                           residual_x_sd,
                           min_n_S = 200,
                           min_clusters_S = 30) {
  margin_main <- 0.05 * residual_x_sd
  margin_robust <- 0.10 * residual_x_sd

  validate_one <- function(margin) {
    bpe_validate_design(
      formula = lelkes_iv_formula(),
      data = dat,
      design = design_info$design,
      vcov = "cluster",
      cluster = ~ state,
      z_names = "log_Total",
      bpe_min_n_S = min_n_S,
      bpe_min_clusters_S = min_clusters_S,
      bpe_min_varZ_S = 1e-6,
      bpe_equiv_margin = margin,
      bpe_equiv_level = 0.95,
      bpe_transport = "sampling",
      bpe_transport_kappa = 0,
      scale_instrument = "residual_sd"
    )
  }

  validation_main <- validate_one(margin_main)
  validation_robust <- validate_one(margin_robust)

  bpe_fit <- NULL
  bpe_error <- NA_character_
  bpe_source <- NA_character_
  if (isTRUE(validation_main$eligibility_passed)) {
    bpe_fit <- tryCatch(
      spliv(
        formula = lelkes_iv_formula(),
        data = dat,
        vcov = "cluster",
        cluster = ~ state,
        method = "bpe",
        bpe_design = design_info$design,
        bpe_spec = list(z_names = "log_Total"),
        bpe_min_n_S = min_n_S,
        bpe_min_clusters_S = min_clusters_S,
        bpe_min_varZ_S = 1e-6,
        bpe_equiv_margin = margin_main,
        bpe_equiv_level = 0.95,
        bpe_transport = "sampling",
        bpe_not_applicable = "na",
        scale_instrument = "residual_sd"
      ),
      error = function(e) e
    )
    if (inherits(bpe_fit, "error")) {
      bpe_error <- conditionMessage(bpe_fit)
      bpe_source <- "direct_spliv_method_bpe_failed"
    } else {
      bpe_source <- "spliv_method_bpe"
    }
  }

  list(
    diagnostics = rbind(
      lelkes_bpe_validation_row(
        validation = validation_main,
        design_info = design_info,
        margin_share = 0.05,
        residual_x_sd = residual_x_sd,
        bpe_fit = if (!inherits(bpe_fit, "error")) bpe_fit else NULL,
        fit_error = bpe_error,
        estimation_source = bpe_source
      ),
      lelkes_bpe_validation_row(
        validation = validation_robust,
        design_info = design_info,
        margin_share = 0.10,
        residual_x_sd = residual_x_sd,
        bpe_fit = NULL,
        fit_error = NA_character_,
        estimation_source = "validation_only"
      )
    ),
    validation_main = validation_main,
    validation_robust = validation_robust,
    fit = bpe_fit,
    main_margin = margin_main,
    robust_margin = margin_robust
  )
}
