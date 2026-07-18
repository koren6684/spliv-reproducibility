lelkes_interval_string <- function(lo, hi, digits = 3) {
  if (!is.finite(lo) || !is.finite(hi)) {
    return(NA_character_)
  }
  paste0("[", formatC(lo, digits = digits, format = "f"), ", ", formatC(hi, digits = digits, format = "f"), "]")
}

lelkes_scalar <- function(x) {
  if (is.null(x) || length(x) == 0) {
    return(NA_real_)
  }
  out <- suppressWarnings(as.numeric(x[[1]]))
  if (length(out) == 0) NA_real_ else out
}

lelkes_ci_values <- function(ci, z_name = "log_Total") {
  if (is.null(ci) || !is.matrix(ci) || !z_name %in% rownames(ci)) {
    return(c(lower = NA_real_, upper = NA_real_))
  }
  c(lower = as.numeric(ci[z_name, "lower"]), upper = as.numeric(ci[z_name, "upper"]))
}

lelkes_bpe_validation_row <- function(validation,
                                      design_info,
                                      margin_share,
                                      residual_x_sd,
                                      bpe_fit = NULL,
                                      fit_error = NA_character_,
                                      estimation_source = NA_character_) {
  instrument <- if (is.null(validation$instrument)) "log_Total" else validation$instrument
  ci <- lelkes_ci_values(validation$first_stage_ci, instrument)
  gamma <- lelkes_scalar(validation$reduced_form_direct_effect)
  gamma_cov <- validation$reduced_form_sampling_cov
  gamma_se <- if (is.matrix(gamma_cov) && all(dim(gamma_cov) >= 1)) sqrt(pmax(0, gamma_cov[1, 1])) else NA_real_
  final_cov <- validation$prior_Omega_sub
  final_cov_scalar <- if (is.matrix(final_cov) && all(dim(final_cov) >= 1)) final_cov[1, 1] else NA_real_

  beta_est <- beta_se <- beta_lo <- beta_hi <- NA_real_
  if (!is.null(bpe_fit) && !inherits(bpe_fit, "error")) {
    beta_row <- tryCatch(lelkes_extract_term(bpe_fit, "log_providers"), error = function(e) NULL)
    if (!is.null(beta_row)) {
      beta_est <- lelkes_scalar(beta_row$estimate)
      beta_se <- lelkes_scalar(beta_row$std.error)
      beta_lo <- lelkes_scalar(beta_row$conf.low)
      beta_hi <- lelkes_scalar(beta_row$conf.high)
    }
  }

  checks <- if (is.null(validation$eligibility_checks)) list() else validation$eligibility_checks
  reason <- as.character(if (is.null(validation$message)) "" else validation$message)
  if (identical(margin_share, 0.10)) {
    reason <- paste(
      reason,
      "10% margin is reported for appendix robustness only and is not used to rescue main BPE eligibility."
    )
  }

  data.frame(
    application = "Lelkes, Sood, and Iyengar Table 1",
    margin_share = margin_share,
    design_name = design_info$design$name,
    subset_definition = "density <= 10th percentile",
    density_threshold = as.numeric(design_info$threshold),
    n_S = lelkes_scalar(validation$n_S),
    share_S = lelkes_scalar(validation$share_S),
    G_S = lelkes_scalar(validation$G_S),
    varZ_S = lelkes_scalar(validation$varZ_S),
    residualized_instrument_sd_S = lelkes_scalar(validation$residualized_instrument_sd_S),
    residualized_treatment_sd_S = lelkes_scalar(validation$residualized_treatment_sd_S),
    residualized_treatment_sd_full = residual_x_sd,
    first_stage_coefficient = lelkes_scalar(validation$first_stage_coefficient),
    first_stage_se = lelkes_scalar(validation$first_stage_se),
    first_stage_ci_low = ci[["lower"]],
    first_stage_ci_high = ci[["upper"]],
    first_stage_f_statistic_diagnostic_only = lelkes_scalar(validation$first_stage_f_statistic),
    equivalence_margin = lelkes_scalar(validation$equivalence_margin),
    equivalence_level = lelkes_scalar(validation$equivalence_level),
    equivalence_passed = isTRUE(validation$equivalence_passed),
    minimum_n_passed = isTRUE(checks$minimum_n),
    minimum_clusters_passed = isTRUE(checks$minimum_clusters),
    residual_variation_passed = isTRUE(checks$residual_variation),
    bpe_eligibility_passed = isTRUE(validation$eligibility_passed),
    gamma_direct_effect_estimate = gamma,
    gamma_direct_effect_se = gamma_se,
    gamma_sampling_cov_scalar = if (is.matrix(gamma_cov) && all(dim(gamma_cov) >= 1)) gamma_cov[1, 1] else NA_real_,
    final_prior_cov_scalar = final_cov_scalar,
    transport_mode = as.character(if (is.null(validation$transport_mode)) NA_character_ else validation$transport_mode),
    transport_uncertainty_inflation = lelkes_scalar(validation$transport_uncertainty_inflation),
    bpe_beta_estimate = beta_est,
    bpe_beta_se = beta_se,
    bpe_beta_conf_low = beta_lo,
    bpe_beta_conf_high = beta_hi,
    bpe_estimation_source = as.character(if (is.null(estimation_source)) NA_character_ else estimation_source),
    reason_if_not_applicable = reason,
    bpe_fit_error = as.character(if (is.null(fit_error)) NA_character_ else fit_error),
    stringsAsFactors = FALSE
  )
}

lelkes_bpe_status_from_diagnostics <- function(bpe_diag) {
  main <- bpe_diag[bpe_diag$margin_share == 0.05, , drop = FALSE]
  if (!nrow(main)) {
    return(data.frame(
      application = "Lelkes, Sood, and Iyengar Table 1",
      bpe_status = "Not applicable",
      reason = "No 5% BPE diagnostic row was produced.",
      density_threshold = NA_real_,
      n_S = NA_real_,
      G_S = NA_real_,
      bpe_eligibility_passed = FALSE,
      bpe_estimation_source = NA_character_,
      stringsAsFactors = FALSE
    ))
  }
  caveat <- paste(
    "This is a conditional diagnostic using a pre-specified bottom-decile density subset.",
    "Transportability from very low-density counties to denser counties should be interpreted cautiously."
  )
  status <- if (isTRUE(main$bpe_eligibility_passed[[1]])) {
    "Eligible"
  } else {
    "Not applicable"
  }
  reason <- if (isTRUE(main$bpe_eligibility_passed[[1]])) {
    caveat
  } else {
    paste(main$reason_if_not_applicable[[1]], caveat)
  }
  data.frame(
    application = main$application[[1]],
    bpe_status = status,
    reason = reason,
    density_threshold = main$density_threshold[[1]],
    n_S = main$n_S[[1]],
    G_S = main$G_S[[1]],
    bpe_eligibility_passed = main$bpe_eligibility_passed[[1]],
    bpe_estimation_source = main$bpe_estimation_source[[1]],
    stringsAsFactors = FALSE
  )
}

lelkes_nearest_delta_row <- function(path, target_share) {
  idx <- which.min(abs(path$delta_share_of_residual_y_sd - target_share))
  if (!length(idx) ||
      !is.finite(path$delta_share_of_residual_y_sd[[idx]]) ||
      abs(path$delta_share_of_residual_y_sd[[idx]] - target_share) > 1e-8) {
    return(NULL)
  }
  path[idx, , drop = FALSE]
}

lelkes_tipping_share <- function(path) {
  rows <- path[order(path$delta_share_of_residual_y_sd), , drop = FALSE]
  if (!nrow(rows)) {
    return(NA_real_)
  }
  if (isTRUE(rows$contains_zero[[1]])) {
    return(0)
  }
  hit <- rows[rows$contains_zero %in% TRUE, , drop = FALSE]
  if (!nrow(hit)) {
    return(NA_real_)
  }
  min(hit$delta_share_of_residual_y_sd, na.rm = TRUE)
}

lelkes_sensitivity_summary <- function(paths, baseline) {
  keys <- unique(paths[, c("method", "pattern_name")])
  out <- lapply(seq_len(nrow(keys)), function(i) {
    key <- keys[i, , drop = FALSE]
    p <- paths[
      paths$method == key$method & paths$pattern_name == key$pattern_name,
      ,
      drop = FALSE
    ]
    rows <- lapply(c(0.05, 0.10, 0.20), function(s) lelkes_nearest_delta_row(p, s))
    names(rows) <- c("s05", "s10", "s20")
    tp <- lelkes_tipping_share(p)
    interpretation <- if (!is.finite(tp)) {
      "No zero crossing on the supplied delta grid."
    } else if (tp == 0) {
      "Baseline interval already includes zero."
    } else {
      paste0("Interval first includes zero at delta share ", formatC(tp, digits = 2, format = "f"), ".")
    }
    data.frame(
      outcome = "affective_polarization",
      treatment = "log_providers",
      instrument = "log_Total",
      method = key$method,
      pattern_name = key$pattern_name,
      baseline_iv_estimate = baseline$estimate,
      baseline_conf_low = baseline$conf_low,
      baseline_conf_high = baseline$conf_high,
      first_stage_coefficient_full = baseline$first_stage_coefficient_full,
      first_stage_f_full_diagnostic_only = baseline$first_stage_f_full_diagnostic_only,
      residual_y_sd = baseline$residual_y_sd,
      residual_x_sd = baseline$residual_x_sd,
      residual_z_sd = baseline$residual_z_sd,
      delta_share_0.05_interval = if (is.null(rows$s05)) NA_character_ else lelkes_interval_string(rows$s05$conf_low, rows$s05$conf_high),
      delta_share_0.10_interval = if (is.null(rows$s10)) NA_character_ else lelkes_interval_string(rows$s10$conf_low, rows$s10$conf_high),
      delta_share_0.20_interval = if (is.null(rows$s20)) NA_character_ else lelkes_interval_string(rows$s20$conf_low, rows$s20$conf_high),
      tipping_point_share = tp,
      interpretation = interpretation,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

lelkes_application_summary <- function(baseline, sensitivity_summary, bpe_status, bpe_diag = NULL) {
  wide <- baseline
  for (method_i in c("uci", "ltz")) {
    for (pattern_i in c("Uniform direct effect", "Density direct-effect pattern")) {
      row <- sensitivity_summary[
        sensitivity_summary$method == method_i & sensitivity_summary$pattern_name == pattern_i,
        ,
        drop = FALSE
      ]
      prefix <- paste0(method_i, "_", if (grepl("Density", pattern_i)) "density" else "uniform")
      for (share in c("0.05", "0.10", "0.20")) {
        src <- paste0("delta_share_", share, "_interval")
        wide[[paste0(prefix, "_delta_", share, "_interval")]] <- if (nrow(row)) row[[src]][[1]] else NA_character_
      }
    }
  }
  wide$bpe_status <- bpe_status$bpe_status[[1]]
  wide$bpe_reason <- bpe_status$reason[[1]]
  if (!is.null(bpe_diag) && nrow(bpe_diag)) {
    main_bpe <- bpe_diag[bpe_diag$margin_share == 0.05, , drop = FALSE]
    if (nrow(main_bpe)) {
      wide$bpe_density_threshold <- main_bpe$density_threshold[[1]]
      wide$bpe_n_S <- main_bpe$n_S[[1]]
      wide$bpe_G_S <- main_bpe$G_S[[1]]
      wide$bpe_equivalence_margin <- main_bpe$equivalence_margin[[1]]
      wide$bpe_equivalence_passed <- main_bpe$equivalence_passed[[1]]
      wide$bpe_eligibility_passed <- main_bpe$bpe_eligibility_passed[[1]]
      wide$bpe_beta_estimate <- main_bpe$bpe_beta_estimate[[1]]
      wide$bpe_beta_conf_low <- main_bpe$bpe_beta_conf_low[[1]]
      wide$bpe_beta_conf_high <- main_bpe$bpe_beta_conf_high[[1]]
      wide$bpe_estimation_source <- main_bpe$bpe_estimation_source[[1]]
    }
  }
  wide
}

lelkes_write_csv <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE)
}

lelkes_plot_paths <- function(paths,
                              file,
                              device = c("pdf", "png"),
                              method_label = "UCI") {
  device <- match.arg(device)
  if (device == "pdf") {
    grDevices::pdf(file, width = 10, height = 5)
  } else {
    grDevices::png(file, width = 1500, height = 800, res = 160)
  }
  on.exit(grDevices::dev.off(), add = TRUE)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mfrow = c(1, 2), mar = c(4, 4, 4, 1), oma = c(4, 0, 2.5, 0))

  patterns <- unique(paths$pattern_name)
  for (pat in patterns) {
    d <- paths[paths$pattern_name == pat, , drop = FALSE]
    d <- d[order(d$delta_share_of_residual_y_sd), , drop = FALSE]
    ylim <- range(c(d$conf_low, d$conf_high, d$estimate, 0), na.rm = TRUE)
    graphics::plot(
      d$delta_share_of_residual_y_sd,
      d$estimate,
      type = "l",
      lwd = 2,
      col = "#1f4e79",
      ylim = ylim,
      xlab = "Delta as share of residual outcome SD",
      ylab = paste0("IV estimate / ", method_label, " interval"),
      main = pat
    )
    graphics::lines(d$delta_share_of_residual_y_sd, d$conf_low, col = "#c0392b", lwd = 1.6, lty = 2)
    graphics::lines(d$delta_share_of_residual_y_sd, d$conf_high, col = "#c0392b", lwd = 1.6, lty = 2)
    graphics::abline(h = 0, col = "gray40", lty = 3)
    graphics::grid(col = "gray90")
  }
  graphics::mtext(paste(method_label, "sensitivity paths"), outer = TRUE, cex = 1.2, font = 2)
}

lelkes_plot_combined_paths <- function(uci_paths, ltz_paths, file, device = c("pdf", "png")) {
  device <- match.arg(device)
  if (device == "pdf") {
    grDevices::pdf(file, width = 10, height = 5)
  } else {
    grDevices::png(file, width = 1500, height = 800, res = 160)
  }
  on.exit(grDevices::dev.off(), add = TRUE)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mfrow = c(1, 2), mar = c(4, 4, 4, 1), oma = c(0, 0, 2.5, 0))

  patterns <- unique(uci_paths$pattern_name)
  for (i in seq_along(patterns)) {
    pat <- patterns[[i]]
    u <- uci_paths[uci_paths$pattern_name == pat, , drop = FALSE]
    l <- ltz_paths[ltz_paths$pattern_name == pat, , drop = FALSE]
    u <- u[order(u$delta_share_of_residual_y_sd), , drop = FALSE]
    l <- l[order(l$delta_share_of_residual_y_sd), , drop = FALSE]
    ylim <- range(c(u$conf_low, u$conf_high, l$conf_low, l$conf_high, 0), na.rm = TRUE)
    graphics::plot(
      u$delta_share_of_residual_y_sd,
      u$estimate,
      type = "l",
      lwd = 2,
      col = "#1f4e79",
      ylim = ylim,
      xlab = "Delta as share of residual outcome SD",
      ylab = "IV estimate / interval",
      main = pat
    )
    graphics::lines(u$delta_share_of_residual_y_sd, u$conf_low, col = "#c0392b", lwd = 1.5, lty = 2)
    graphics::lines(u$delta_share_of_residual_y_sd, u$conf_high, col = "#c0392b", lwd = 1.5, lty = 2)
    graphics::lines(l$delta_share_of_residual_y_sd, l$conf_low, col = "#1b7f5a", lwd = 1.5, lty = 3)
    graphics::lines(l$delta_share_of_residual_y_sd, l$conf_high, col = "#1b7f5a", lwd = 1.5, lty = 3)
    graphics::abline(h = 0, col = "gray40", lty = 3)
    if (i == 1L) {
      graphics::legend(
        "topleft",
        legend = c("Estimate", "UCI bounds", "LTZ bounds"),
        col = c("#1f4e79", "#c0392b", "#1b7f5a"),
        lty = c(1, 2, 3),
        lwd = c(2, 1.5, 1.5),
        bty = "n"
      )
    }
    graphics::grid(col = "gray90")
  }
  graphics::mtext("Lelkes Table 1 patterned sensitivity paths", outer = TRUE, cex = 1.2, font = 2)
}

lelkes_plot_bpe <- function(bpe_diag, file, device = c("pdf", "png")) {
  device <- match.arg(device)
  if (device == "pdf") {
    grDevices::pdf(file, width = 8, height = 4.8)
  } else {
    grDevices::png(file, width = 1300, height = 760, res = 160)
  }
  on.exit(grDevices::dev.off(), add = TRUE)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mar = c(5, 11, 4, 1), xpd = FALSE)

  d <- bpe_diag[bpe_diag$margin_share == 0.05, , drop = FALSE]
  if (!nrow(d)) {
    graphics::plot.new()
    graphics::text(0.5, 0.5, "No BPE diagnostics available.")
    return(invisible(NULL))
  }
  denom <- d$equivalence_margin
  coef_norm <- d$first_stage_coefficient / denom
  lo_norm <- d$first_stage_ci_low / denom
  hi_norm <- d$first_stage_ci_high / denom
  passed <- d$bpe_eligibility_passed %in% TRUE
  point_col <- ifelse(passed, "#1b7f5a", "#b23b3b")
  status <- ifelse(passed, "eligible", "not eligible")
  y <- 1
  label <- paste0("Low-density: ", status)
  xlim <- range(c(lo_norm, hi_norm, coef_norm, -1, 0, 1), na.rm = TRUE)
  if (!all(is.finite(xlim))) {
    xlim <- c(-1.5, 1.5)
  }
  xpad <- diff(xlim) * 0.08
  if (!is.finite(xpad) || xpad == 0) {
    xpad <- 0.25
  }
  xlim <- xlim + c(-xpad, xpad)

  graphics::plot(
    NA_real_,
    NA_real_,
    xlim = xlim,
    ylim = c(0.5, 1.5),
    yaxt = "n",
    xlab = "First-stage coefficient divided by 5% equivalence margin",
    ylab = "",
    main = "Lelkes confirmatory BPE normalized equivalence diagnostic"
  )
  graphics::grid(col = "gray90")
  graphics::abline(v = c(-1, 1), col = "#c0392b", lty = 2, lwd = 1.5)
  graphics::abline(v = 0, col = "gray35", lty = 3, lwd = 1.5)
  graphics::axis(2, at = y, labels = label, las = 1)
  graphics::segments(lo_norm, y, hi_norm, y, lwd = 2, col = point_col)
  graphics::points(coef_norm, y, pch = 19, col = point_col, cex = 1.3)
  graphics::legend(
    "topright",
    legend = c("+/-1 equivalence margin", "Zero"),
    col = c("#c0392b", "gray35"),
    lty = c(2, 3),
    lwd = c(1.5, 1.5),
    bty = "n"
  )
}

lelkes_write_log <- function(path, lines) {
  writeLines(lines, con = path)
}
