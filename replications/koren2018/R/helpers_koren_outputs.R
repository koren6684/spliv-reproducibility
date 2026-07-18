koren_interval_string <- function(lo, hi, digits = 3) {
  if (!is.finite(lo) || !is.finite(hi)) return(NA_character_)
  paste0("[", formatC(lo, digits = digits, format = "f"), ", ", formatC(hi, digits = digits, format = "f"), "]")
}

koren_nearest_delta_row <- function(path, target_share) {
  idx <- which.min(abs(path$delta_share_of_residual_y_sd - target_share))
  if (!length(idx) || !is.finite(path$delta_share_of_residual_y_sd[[idx]]) ||
      abs(path$delta_share_of_residual_y_sd[[idx]] - target_share) > 1e-8) {
    return(NULL)
  }
  path[idx, , drop = FALSE]
}

koren_tipping_share <- function(path) {
  rows <- path[order(path$delta_share_of_residual_y_sd), , drop = FALSE]
  if (!nrow(rows)) return(NA_real_)
  if (isTRUE(rows$contains_zero[[1]])) return(0)
  hit <- rows[rows$contains_zero %in% TRUE, , drop = FALSE]
  if (!nrow(hit)) return(NA_real_)
  min(hit$delta_share_of_residual_y_sd, na.rm = TRUE)
}

koren_sensitivity_summary <- function(paths, baseline) {
  keys <- unique(paths[, c("crop", "treatment", "pattern_name")])
  out <- lapply(seq_len(nrow(keys)), function(i) {
    key <- keys[i, , drop = FALSE]
    p <- paths[
      paths$crop == key$crop & paths$treatment == key$treatment & paths$pattern_name == key$pattern_name,
      ,
      drop = FALSE
    ]
    b <- baseline[baseline$crop == key$crop & baseline$treatment == key$treatment, , drop = FALSE][1, , drop = FALSE]
    rows <- lapply(c(0.05, 0.10, 0.20), function(s) koren_nearest_delta_row(p, s))
    names(rows) <- c("s05", "s10", "s20")
    tp <- koren_tipping_share(p)
    interpretation <- if (!is.finite(tp)) {
      "No zero crossing on the supplied delta grid."
    } else if (tp == 0) {
      "Baseline interval already includes zero."
    } else {
      paste0("Interval first includes zero at delta share ", formatC(tp, digits = 2, format = "f"), ".")
    }
    data.frame(
      crop = key$crop,
      treatment = key$treatment,
      pattern_name = key$pattern_name,
      baseline_iv_estimate = b$estimate,
      baseline_conf_low = b$conf_low,
      baseline_conf_high = b$conf_high,
      residual_y_sd = b$residual_y_sd,
      delta_share_0.05_interval = if (is.null(rows$s05)) NA_character_ else koren_interval_string(rows$s05$conf_low, rows$s05$conf_high),
      delta_share_0.10_interval = if (is.null(rows$s10)) NA_character_ else koren_interval_string(rows$s10$conf_low, rows$s10$conf_high),
      delta_share_0.20_interval = if (is.null(rows$s20)) NA_character_ else koren_interval_string(rows$s20$conf_low, rows$s20$conf_high),
      tipping_point_share = tp,
      interpretation = interpretation,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

koren_add_ltz_summary_columns <- function(out, ltz_summary) {
  if (is.null(ltz_summary) || !nrow(ltz_summary)) {
    return(out)
  }
  target_patterns <- list(
    ltz_uniform = "Uniform direct effect",
    ltz_sparsebare = "Sparse/bare direct-effect pattern"
  )
  target_shares <- c("0.05", "0.10", "0.20")
  for (prefix in names(target_patterns)) {
    rows <- ltz_summary[ltz_summary$pattern_name == target_patterns[[prefix]], , drop = FALSE]
    keep_cols <- c("crop", "treatment")
    for (share in target_shares) {
      src <- paste0("delta_share_", share, "_interval")
      dst <- paste0(prefix, "_delta_", share, "_interval")
      rows[[dst]] <- if (src %in% names(rows)) rows[[src]] else NA_character_
      keep_cols <- c(keep_cols, dst)
    }
    out <- merge(out, rows[, keep_cols, drop = FALSE], by = c("crop", "treatment"), all.x = TRUE)
  }
  out
}

koren_application_summary <- function(baseline, sensitivity_summary, bpe_diag, ltz_summary = NULL) {
  main_bpe <- bpe_diag[bpe_diag$margin_share == 0.05, , drop = FALSE]
  out <- merge(
    baseline,
    main_bpe[, c("crop", "treatment", "design_name", "subset_definition", "design_variable",
                 "subset_threshold", "threshold_label", "sparsebare_threshold",
                 "n_S", "G_S", "equivalence_margin",
                 "equivalence_passed", "bpe_eligibility_passed", "bpe_beta_estimate",
                 "bpe_beta_conf_low", "bpe_beta_conf_high", "bpe_estimation_source",
                 "reason_if_not_applicable", "bpe_fit_error")],
    by = c("crop", "treatment"),
    all.x = TRUE
  )
  sparse_tp <- sensitivity_summary[sensitivity_summary$pattern_name == "Sparse/bare direct-effect pattern",
                                   c("crop", "treatment", "tipping_point_share"),
                                   drop = FALSE]
  names(sparse_tp)[names(sparse_tp) == "tipping_point_share"] <- "sparsebare_tipping_point_share"
  uniform_tp <- sensitivity_summary[sensitivity_summary$pattern_name == "Uniform direct effect",
                                    c("crop", "treatment", "tipping_point_share"),
                                    drop = FALSE]
  names(uniform_tp)[names(uniform_tp) == "tipping_point_share"] <- "uniform_tipping_point_share"
  out <- merge(out, sparse_tp, by = c("crop", "treatment"), all.x = TRUE)
  out <- merge(out, uniform_tp, by = c("crop", "treatment"), all.x = TRUE)
  out <- koren_add_ltz_summary_columns(out, ltz_summary)
  out
}

koren_write_csv <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE)
}

koren_plot_sensitivity <- function(paths,
                                   file,
                                   device = c("pdf", "png"),
                                   method_label = "UCI",
                                   overall_title = NULL) {
  device <- match.arg(device)
  if (device == "pdf") {
    grDevices::pdf(file, width = 10, height = 7)
  } else {
    grDevices::png(file, width = 1600, height = 1100, res = 160)
  }
  on.exit(grDevices::dev.off(), add = TRUE)

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(
    mfrow = c(2, 2),
    mar = c(4, 4, 3, 1),
    oma = if (is.null(overall_title)) c(0, 0, 0, 0) else c(0, 0, 3, 0)
  )

  crops <- unique(paths$crop)
  patterns <- unique(paths$pattern_name)
  for (crop in crops) {
    for (pat in patterns) {
      d <- paths[paths$crop == crop & paths$pattern_name == pat, , drop = FALSE]
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
        main = paste0(tools::toTitleCase(crop), ": ", pat)
      )
      graphics::lines(d$delta_share_of_residual_y_sd, d$conf_low, col = "#c0392b", lwd = 1.6, lty = 2)
      graphics::lines(d$delta_share_of_residual_y_sd, d$conf_high, col = "#c0392b", lwd = 1.6, lty = 2)
      graphics::abline(h = 0, col = "gray40", lty = 3)
      graphics::grid(col = "gray90")
    }
  }
  if (!is.null(overall_title)) {
    graphics::mtext(overall_title, outer = TRUE, cex = 1.15, font = 2)
  }
}

koren_plot_bpe <- function(bpe_diag, file, device = c("pdf", "png")) {
  device <- match.arg(device)
  if (device == "pdf") {
    grDevices::pdf(file, width = 9, height = 5.8)
  } else {
    grDevices::png(file, width = 1500, height = 950, res = 160)
  }
  on.exit(grDevices::dev.off(), add = TRUE)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mar = c(8, 9, 4, 1), xpd = FALSE)

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
  pass <- d$bpe_eligibility_passed %in% TRUE
  point_col <- ifelse(pass, "#1b7f5a", "#b23b3b")
  status <- ifelse(pass, "eligible", "not eligible")
  y <- rev(seq_len(nrow(d)))
  labels <- paste0(tools::toTitleCase(d$crop), ": ", status)
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
    ylim = range(y) + c(-0.6, 0.6),
    yaxt = "n",
    xlab = "First-stage coefficient divided by 5% equivalence margin",
    ylab = "",
    main = "Primary BPE normalized equivalence diagnostics"
  )
  graphics::grid(col = "gray90")
  graphics::abline(v = c(-1, 1), col = "#c0392b", lty = 2, lwd = 1.5)
  graphics::abline(v = 0, col = "gray35", lty = 3, lwd = 1.5)
  graphics::axis(2, at = y, labels = labels, las = 1)
  finite_rows <- is.finite(coef_norm) & is.finite(lo_norm) & is.finite(hi_norm)
  graphics::segments(lo_norm[finite_rows], y[finite_rows], hi_norm[finite_rows], y[finite_rows], lwd = 2, col = point_col[finite_rows])
  graphics::points(coef_norm[finite_rows], y[finite_rows], pch = 19, col = point_col[finite_rows], cex = 1.2)
  if (any(!finite_rows)) {
    graphics::text(
      x = 0,
      y = y[!finite_rows],
      labels = paste0("not applicable (n_S=", d$n_S[!finite_rows], ")"),
      col = "#b23b3b",
      cex = 0.9
    )
  }
  graphics::mtext(
    "Red dashed lines mark +/-1 equivalence margins; gray dotted line marks zero.",
    side = 1,
    line = 4,
    cex = 0.8,
    col = "gray35"
  )
}

koren_write_log <- function(path, lines) {
  writeLines(lines, con = path)
}
