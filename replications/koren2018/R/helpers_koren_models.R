koren_extract_term <- function(fit, term) {
  tab <- as.data.frame(fit$estimates)
  if (!"term" %in% names(tab)) {
    stop("`fit$estimates` does not contain a `term` column.", call. = FALSE)
  }
  rw <- tab[tab$term == term, , drop = FALSE]
  if (nrow(rw) == 0) {
    stop("Could not find term `", term, "` in fitted SPLIV object.", call. = FALSE)
  }
  rw[1, , drop = FALSE]
}

koren_baseline_fit <- function(dat, treatment) {
  spliv(
    formula = koren_iv_formula(treatment),
    data = dat,
    fe = ~ gid + year,
    vcov = "cluster",
    cluster = ~ gid,
    method = "ltz",
    delta = 0,
    scale_instrument = "residual_sd"
  )
}

koren_full_first_stage <- function(fit, treatment) {
  x <- as.numeric(fit$internals$X[, treatment])
  z <- as.numeric(fit$internals$Z[, "spi6"])
  lm_fit <- stats::lm(x ~ z)
  sm <- summary(lm_fit)
  coef_tab <- sm$coefficients
  z_row <- coef_tab["z", , drop = FALSE]
  fstat <- sm$fstatistic
  f_val <- if (length(fstat)) as.numeric(fstat[["value"]]) else NA_real_
  c(
    coefficient = as.numeric(z_row[1, "Estimate"]),
    se = as.numeric(z_row[1, "Std. Error"]),
    f_statistic = f_val
  )
}

koren_residual_metrics <- function(fit, treatment) {
  y_resid <- as.numeric(fit$internals$y)
  x_resid <- as.numeric(fit$internals$X[, treatment])
  z_sd <- as.numeric(fit$residualized_instrument_sd[["spi6"]])
  list(
    residual_y_sd = stats::sd(y_resid),
    residual_x_sd = stats::sd(x_resid),
    residual_z_sd = z_sd,
    nobs = length(y_resid),
    n_clusters = length(unique(as.character(fit$internals$cluster_id))),
    first_stage_full = koren_full_first_stage(fit, treatment)
  )
}

koren_run_sensitivity_path <- function(dat,
                                       treatment,
                                       pattern,
                                       delta_grid_y,
                                       method = c("uci", "ltz"),
                                       uci_steps = 21L) {
  method <- match.arg(method)
  extra_args <- if (identical(method, "uci")) {
    list(grid = list(steps = as.integer(uci_steps), level = 0.95))
  } else {
    list(grid = list(level = 0.95))
  }
  path <- do.call(
    spliv_sensitivity_path,
    c(
      list(
        formula = koren_iv_formula(treatment),
        data = dat,
        fe = ~ gid + year,
        vcov = "cluster",
        cluster = ~ gid,
        method = method,
        delta_grid = delta_grid_y,
        violation_pattern = pattern,
        scale_instrument = "residual_sd"
      ),
      extra_args
    )
  )
  path[path$term == treatment, , drop = FALSE]
}

koren_scalar <- function(x) {
  if (is.null(x) || length(x) == 0) return(NA_real_)
  out <- suppressWarnings(as.numeric(x[[1]]))
  if (length(out) == 0) NA_real_ else out
}

koren_ci_values <- function(ci, z_name = "spi6") {
  if (is.null(ci) || !is.matrix(ci) || !z_name %in% rownames(ci)) {
    return(c(lower = NA_real_, upper = NA_real_))
  }
  c(lower = as.numeric(ci[z_name, "lower"]), upper = as.numeric(ci[z_name, "upper"]))
}

koren_bpe_subset_index <- function(dat, design_info) {
  idx <- tryCatch(design_info$design$subset(dat), error = function(e) rep(FALSE, nrow(dat)))
  idx <- idx %in% TRUE
  if (length(idx) != nrow(dat)) {
    idx <- rep(FALSE, nrow(dat))
  }
  idx
}

koren_bpe_not_applicable_validation <- function(dat, design_info, margin, message) {
  idx <- koren_bpe_subset_index(dat, design_info)
  n_S <- sum(idx, na.rm = TRUE)
  G_S <- if ("gid" %in% names(dat) && n_S > 0) {
    length(unique(dat$gid[idx]))
  } else {
    0L
  }
  list(
    instrument = "spi6",
    n_S = n_S,
    share_S = if (nrow(dat)) n_S / nrow(dat) else NA_real_,
    G_S = G_S,
    varZ_S = NA_real_,
    residualized_instrument_sd_S = NA_real_,
    residualized_treatment_sd_S = NA_real_,
    first_stage_coefficient = NA_real_,
    first_stage_se = NA_real_,
    first_stage_ci = NULL,
    first_stage_f_statistic = NA_real_,
    equivalence_margin = margin,
    equivalence_level = 0.95,
    equivalence_passed = FALSE,
    eligibility_checks = list(
      minimum_n = n_S >= 200,
      minimum_clusters = G_S >= 30,
      residual_variation = FALSE
    ),
    eligibility_passed = FALSE,
    reduced_form_direct_effect = NA_real_,
    reduced_form_sampling_cov = matrix(NA_real_, 1, 1, dimnames = list("spi6", "spi6")),
    prior_Omega_sub = matrix(NA_real_, 1, 1, dimnames = list("spi6", "spi6")),
    transport_mode = "sampling",
    transport_uncertainty_inflation = NA_real_,
    message = message
  )
}

koren_bpe_validation_row <- function(validation,
                                     crop,
                                     treatment,
                                     design_info,
                                     margin_share,
                                     bpe_fit = NULL,
                                     fit_error = NA_character_,
                                     estimation_source = NA_character_) {
  ci <- koren_ci_values(validation$first_stage_ci, validation$instrument %||% "spi6")
  gamma <- koren_scalar(validation$reduced_form_direct_effect)
  gamma_cov <- validation$reduced_form_sampling_cov
  gamma_se <- if (is.matrix(gamma_cov) && all(dim(gamma_cov) >= 1)) sqrt(pmax(0, gamma_cov[1, 1])) else NA_real_
  final_cov <- validation$prior_Omega_sub
  final_cov_scalar <- if (is.matrix(final_cov) && all(dim(final_cov) >= 1)) final_cov[1, 1] else NA_real_

  beta_est <- beta_se <- beta_lo <- beta_hi <- NA_real_
  if (!is.null(bpe_fit) && !inherits(bpe_fit, "error")) {
    beta_row <- tryCatch(koren_extract_term(bpe_fit, treatment), error = function(e) NULL)
    if (!is.null(beta_row)) {
      beta_est <- koren_scalar(beta_row$estimate)
      beta_se <- koren_scalar(beta_row$std.error)
      beta_lo <- koren_scalar(beta_row$conf.low)
      beta_hi <- koren_scalar(beta_row$conf.high)
    }
  }

  checks <- validation$eligibility_checks %||% list()
  data.frame(
    crop = crop,
    treatment = treatment,
    design_role = as.character(design_info$design_role %||% NA_character_),
    design_id = as.character(design_info$design_id %||% NA_character_),
    design_name = as.character(design_info$design_name %||% design_info$design$name %||% NA_character_),
    margin_share = margin_share,
    subset_definition = as.character(design_info$subset_rule %||% NA_character_),
    design_variable = as.character(design_info$design_variable %||% NA_character_),
    subset_threshold = koren_scalar(design_info$threshold),
    threshold_label = as.character(design_info$threshold_label %||% NA_character_),
    sparsebare_threshold = if (identical(design_info$design_variable, "sparsebare")) koren_scalar(design_info$threshold) else NA_real_,
    n_S = koren_scalar(validation$n_S),
    share_S = koren_scalar(validation$share_S),
    G_S = koren_scalar(validation$G_S),
    varZ_S = koren_scalar(validation$varZ_S),
    residualized_instrument_sd_S = koren_scalar(validation$residualized_instrument_sd_S),
    residualized_treatment_sd_S = koren_scalar(validation$residualized_treatment_sd_S),
    residualized_treatment_sd_full = NA_real_,
    first_stage_coefficient = koren_scalar(validation$first_stage_coefficient),
    first_stage_se = koren_scalar(validation$first_stage_se),
    first_stage_ci_low = ci[["lower"]],
    first_stage_ci_high = ci[["upper"]],
    first_stage_f_statistic_diagnostic_only = koren_scalar(validation$first_stage_f_statistic),
    equivalence_margin = koren_scalar(validation$equivalence_margin),
    equivalence_level = koren_scalar(validation$equivalence_level),
    equivalence_passed = isTRUE(validation$equivalence_passed),
    minimum_n_passed = isTRUE(checks$minimum_n),
    minimum_clusters_passed = isTRUE(checks$minimum_clusters),
    residual_variation_passed = isTRUE(checks$residual_variation),
    bpe_eligibility_passed = isTRUE(validation$eligibility_passed),
    gamma_direct_effect_estimate = gamma,
    gamma_direct_effect_se = gamma_se,
    gamma_sampling_cov_scalar = if (is.matrix(gamma_cov) && all(dim(gamma_cov) >= 1)) gamma_cov[1, 1] else NA_real_,
    final_prior_cov_scalar = final_cov_scalar,
    transport_mode = as.character(validation$transport_mode %||% NA_character_),
    transport_uncertainty_inflation = koren_scalar(validation$transport_uncertainty_inflation),
    bpe_beta_estimate = beta_est,
    bpe_beta_se = beta_se,
    bpe_beta_conf_low = beta_lo,
    bpe_beta_conf_high = beta_hi,
    bpe_estimation_source = as.character(estimation_source %||% NA_character_),
    reason_if_not_applicable = as.character(validation$message %||% ""),
    bpe_fit_error = as.character(fit_error %||% NA_character_),
    stringsAsFactors = FALSE
  )
}

koren_run_bpe <- function(dat,
                          crop,
                          treatment,
                          design_info,
                          residual_x_sd,
                          min_n_S = 200,
                          min_clusters_S = 30) {
  margin_main <- 0.05 * residual_x_sd
  margin_robust <- 0.10 * residual_x_sd

  validate_one <- function(margin) {
    tryCatch(
      bpe_validate_design(
        formula = koren_iv_formula(treatment),
        data = dat,
        design = design_info$design,
        fe = ~ gid + year,
        vcov = "cluster",
        cluster = ~ gid,
        z_names = "spi6",
        bpe_min_n_S = min_n_S,
        bpe_min_clusters_S = min_clusters_S,
        bpe_min_varZ_S = 1e-6,
        bpe_equiv_margin = margin,
        bpe_equiv_level = 0.95,
        bpe_transport = "sampling",
        bpe_transport_kappa = 0,
        scale_instrument = "residual_sd"
      ),
      error = function(e) {
        msg <- paste0(
          "BPE design `", design_info$design_name, "` is not applicable in the crop-specific analysis sample: ",
          conditionMessage(e)
        )
        koren_bpe_not_applicable_validation(dat, design_info, margin, msg)
      }
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
        formula = koren_iv_formula(treatment),
        data = dat,
        fe = ~ gid + year,
        vcov = "cluster",
        cluster = ~ gid,
        method = "bpe",
        bpe_design = design_info$design,
        bpe_spec = list(z_names = "spi6"),
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

  main_row <- koren_bpe_validation_row(
    validation = validation_main,
    crop = crop,
    treatment = treatment,
    design_info = design_info,
    margin_share = 0.05,
    bpe_fit = if (!inherits(bpe_fit, "error")) bpe_fit else NULL,
    fit_error = bpe_error,
    estimation_source = bpe_source
  )
  robust_row <- koren_bpe_validation_row(
    validation = validation_robust,
    crop = crop,
    treatment = treatment,
    design_info = design_info,
    margin_share = 0.10,
    bpe_fit = NULL,
    fit_error = NA_character_,
    estimation_source = "validation_only"
  )
  robust_row$reason_if_not_applicable <- paste(
    robust_row$reason_if_not_applicable,
    "10% margin is reported for appendix robustness only and is not used to rescue main BPE eligibility."
  )

  list(
    diagnostics = rbind(main_row, robust_row),
    validation_main = validation_main,
    validation_robust = validation_robust,
    fit = bpe_fit,
    main_margin = margin_main,
    robust_margin = margin_robust
  )
}
