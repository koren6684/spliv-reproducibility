spliv_sim_has_ggplot2 <- function() {
  requireNamespace("ggplot2", quietly = TRUE)
}

plot_patterned_paths <- function(path_df, file, title = "Sensitivity Paths") {
  if (nrow(path_df) == 0L) {
    return(invisible(NULL))
  }

  spliv_sim_assert_columns(
    path_df,
    c("scenario_name", "truth_pattern", "pattern_label", "delta", "estimate", "conf_low", "conf_high"),
    df_name = "path_df"
  )

  summary_df <- spliv_sim_grouped_summary(
    path_df,
    group_vars = c("scenario_name", "truth_pattern", "pattern_label", "delta"),
    summarize_fn = function(df) {
      data.frame(
        scenario_name = df$scenario_name[1],
        truth_pattern = df$truth_pattern[1],
        pattern_label = df$pattern_label[1],
        delta = df$delta[1],
        estimate = spliv_sim_safe_mean(df$estimate),
        conf_low = spliv_sim_safe_mean(df$conf_low),
        conf_high = spliv_sim_safe_mean(df$conf_high),
        stringsAsFactors = FALSE
      )
    },
    df_name = "path_df"
  )

  if (spliv_sim_has_ggplot2()) {
    gg <- ggplot2::ggplot(
      summary_df,
      ggplot2::aes(x = delta, y = estimate, color = pattern_label, fill = pattern_label)
    ) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = conf_low, ymax = conf_high),
        alpha = 0.12,
        color = NA
      ) +
      ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::facet_wrap(~ scenario_name, scales = "free_y") +
      ggplot2::geom_hline(yintercept = 0, linetype = 3) +
      ggplot2::labs(
        title = title,
        x = expression(delta),
        y = "Estimate / interval"
      ) +
      ggplot2::theme_minimal()
    ggplot2::ggsave(file, gg, width = 11, height = 7, dpi = 180)
    return(invisible(file))
  }

  grDevices::png(file, width = 1200, height = 800, res = 160)
  on.exit(grDevices::dev.off(), add = TRUE)
  scenarios <- unique(summary_df$scenario_name)
  graphics::par(mfrow = c(length(scenarios), 1), mar = c(4, 4, 2, 1))
  for (scenario_name in scenarios) {
    dat <- summary_df[summary_df$scenario_name == scenario_name, , drop = FALSE]
    finite_vals <- c(dat$estimate, dat$conf_low, dat$conf_high)
    finite_vals <- finite_vals[is.finite(finite_vals)]
    ylim <- if (length(finite_vals) == 0L) c(-1, 1) else range(finite_vals)
    graphics::plot(
      NA,
      xlim = range(dat$delta),
      ylim = ylim,
      xlab = "delta",
      ylab = "Estimate / interval",
      main = scenario_name
    )
    graphics::abline(h = 0, lty = 3)
    labels <- unique(dat$pattern_label)
    cols <- seq_along(labels)
    for (idx in seq_along(labels)) {
      subdat <- dat[dat$pattern_label == labels[idx], , drop = FALSE]
      graphics::lines(subdat$delta, subdat$estimate, col = cols[idx], lwd = 2)
      graphics::lines(subdat$delta, subdat$conf_low, col = cols[idx], lty = 2)
      graphics::lines(subdat$delta, subdat$conf_high, col = cols[idx], lty = 2)
    }
    graphics::legend("topright", legend = labels, col = cols, lty = 1, bty = "n")
  }
  invisible(file)
}

plot_patterned_metric <- function(summary_df, value_col, file, title, ylab) {
  if (nrow(summary_df) == 0L) {
    return(invisible(NULL))
  }
  spliv_sim_assert_columns(summary_df, c("delta", "pattern_label", "truth_pattern", value_col), df_name = "summary_df")

  if (spliv_sim_has_ggplot2()) {
    gg <- ggplot2::ggplot(
      summary_df,
      ggplot2::aes(x = delta, y = .data[[value_col]], color = pattern_label)
    ) +
      ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::facet_wrap(~ truth_pattern) +
      ggplot2::labs(title = title, x = expression(delta), y = ylab) +
      ggplot2::theme_minimal()
    ggplot2::ggsave(file, gg, width = 10, height = 6, dpi = 180)
    return(invisible(file))
  }

  grDevices::png(file, width = 1200, height = 700, res = 160)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::plot(
    summary_df$delta,
    summary_df[[value_col]],
    xlab = "delta",
    ylab = ylab,
    main = title,
    pch = 19
  )
  invisible(file)
}

plot_bpe_eligibility <- function(summary_df, file, title = "BPE Eligibility Rates") {
  if (nrow(summary_df) == 0L) {
    return(invisible(NULL))
  }
  spliv_sim_assert_columns(summary_df, c("scenario_name", "eligibility_rate", "transport_gap"), df_name = "summary_df")

  if (spliv_sim_has_ggplot2()) {
    gg <- ggplot2::ggplot(
      summary_df,
      ggplot2::aes(x = scenario_name, y = eligibility_rate, fill = factor(transport_gap))
    ) +
      ggplot2::geom_col(position = "dodge") +
      ggplot2::coord_cartesian(ylim = c(0, 1)) +
      ggplot2::labs(
        title = title,
        x = "Scenario",
        y = "Eligibility rate",
        fill = "Transport gap"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    ggplot2::ggsave(file, gg, width = 10, height = 6, dpi = 180)
    return(invisible(file))
  }

  grDevices::png(file, width = 1200, height = 700, res = 160)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::barplot(summary_df$eligibility_rate, names.arg = summary_df$scenario_name, las = 2, ylim = c(0, 1), main = title)
  invisible(file)
}

plot_search_false_rates <- function(summary_df, file, title = "False Subgroup-Search Rates") {
  if (nrow(summary_df) == 0L) {
    return(invisible(NULL))
  }
  spliv_sim_assert_columns(
    summary_df,
    c("K", "false_f_rate", "false_ci_rate", "false_equivalence_rate"),
    df_name = "summary_df"
  )

  long_df <- spliv_sim_rbind_fill(list(
    data.frame(K = summary_df$K, rule = "F <= 5", rate = summary_df$false_f_rate, stringsAsFactors = FALSE),
    data.frame(K = summary_df$K, rule = "CI includes 0", rate = summary_df$false_ci_rate, stringsAsFactors = FALSE),
    data.frame(K = summary_df$K, rule = "Equivalence pass", rate = summary_df$false_equivalence_rate, stringsAsFactors = FALSE)
  ))

  if (spliv_sim_has_ggplot2()) {
    gg <- ggplot2::ggplot(long_df, ggplot2::aes(x = factor(K), y = rate, fill = rule)) +
      ggplot2::geom_col(position = "dodge") +
      ggplot2::coord_cartesian(ylim = c(0, 1)) +
      ggplot2::labs(title = title, x = "Number of searched groups (K)", y = "False-selection rate") +
      ggplot2::theme_minimal()
    ggplot2::ggsave(file, gg, width = 9, height = 6, dpi = 180)
    return(invisible(file))
  }

  grDevices::png(file, width = 1200, height = 700, res = 160)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::barplot(long_df$rate, names.arg = paste(long_df$rule, long_df$K, sep = "\n"), las = 2, ylim = c(0, 1), main = title)
  invisible(file)
}

plot_search_stress_false_rates <- function(summary_df, file, title = "False Subgroup-Search Stress Test") {
  if (nrow(summary_df) == 0L) {
    return(invisible(NULL))
  }
  spliv_sim_assert_columns(
    summary_df,
    c(
      "K", "pi", "subgroup_share", "noise_multiplier",
      "false_f_rate", "false_ci_rate", "false_equivalence_rate"
    ),
    df_name = "summary_df"
  )

  long_df <- spliv_sim_rbind_fill(list(
    data.frame(
      K = summary_df$K,
      pi = summary_df$pi,
      subgroup_share = summary_df$subgroup_share,
      noise_multiplier = summary_df$noise_multiplier,
      rule = "F <= 5",
      rate = summary_df$false_f_rate,
      stringsAsFactors = FALSE
    ),
    data.frame(
      K = summary_df$K,
      pi = summary_df$pi,
      subgroup_share = summary_df$subgroup_share,
      noise_multiplier = summary_df$noise_multiplier,
      rule = "CI includes 0",
      rate = summary_df$false_ci_rate,
      stringsAsFactors = FALSE
    ),
    data.frame(
      K = summary_df$K,
      pi = summary_df$pi,
      subgroup_share = summary_df$subgroup_share,
      noise_multiplier = summary_df$noise_multiplier,
      rule = "Equivalence pass",
      rate = summary_df$false_equivalence_rate,
      stringsAsFactors = FALSE
    )
  ))

  if (spliv_sim_has_ggplot2()) {
    gg <- ggplot2::ggplot(
      long_df,
      ggplot2::aes(x = factor(K), y = rate, fill = rule)
    ) +
      ggplot2::geom_col(position = "dodge") +
      ggplot2::coord_cartesian(ylim = c(0, 1)) +
      ggplot2::facet_grid(
        rows = ggplot2::vars(noise_multiplier, subgroup_share),
        cols = ggplot2::vars(pi),
        labeller = ggplot2::label_both
      ) +
      ggplot2::labs(
        title = title,
        x = "Number of searched groups (K)",
        y = "False-selection rate",
        fill = "Rule"
      ) +
      ggplot2::theme_minimal()
    ggplot2::ggsave(file, gg, width = 14, height = 10, dpi = 180)
    return(invisible(file))
  }

  panel_keys <- unique(long_df[, c("noise_multiplier", "subgroup_share", "pi"), drop = FALSE])
  grDevices::png(file, width = 1800, height = 1400, res = 180)
  on.exit(grDevices::dev.off(), add = TRUE)
  n_panels <- nrow(panel_keys)
  n_cols <- max(1L, ceiling(sqrt(n_panels)))
  n_rows <- max(1L, ceiling(n_panels / n_cols))
  graphics::par(mfrow = c(n_rows, n_cols), mar = c(5, 4, 3, 1))
  for (idx in seq_len(n_panels)) {
    key <- panel_keys[idx, , drop = FALSE]
    dat <- long_df[
      long_df$noise_multiplier == key$noise_multiplier &
        long_df$subgroup_share == key$subgroup_share &
        long_df$pi == key$pi,
      ,
      drop = FALSE
    ]
    labels <- paste(dat$rule, dat$K, sep = "\n")
    graphics::barplot(
      dat$rate,
      names.arg = labels,
      las = 2,
      ylim = c(0, 1),
      main = paste(
        "pi =", key$pi,
        "| share =", key$subgroup_share,
        "| noise =", key$noise_multiplier
      )
    )
  }
  invisible(file)
}
