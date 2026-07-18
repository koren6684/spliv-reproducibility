spliv_sim_safe_mean <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) == 0L || all(is.na(x))) {
    return(NA_real_)
  }
  mean(x, na.rm = TRUE)
}

spliv_sim_safe_median <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) == 0L || all(is.na(x))) {
    return(NA_real_)
  }
  stats::median(x, na.rm = TRUE)
}

spliv_sim_assert_columns <- function(df, columns, df_name = "data") {
  missing_cols <- setdiff(columns, names(df))
  if (length(missing_cols) > 0L) {
    stop(
      sprintf(
        "`%s` is missing required column(s): %s",
        df_name,
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(df)
}

spliv_sim_group_keys <- function(df, group_vars, df_name = "data") {
  spliv_sim_assert_columns(df, group_vars, df_name)
  do.call(paste, c(lapply(group_vars, function(var_name) {
    values <- as.character(df[[var_name]])
    values[is.na(values)] <- "<NA>"
    values
  }), sep = "\r"))
}

spliv_sim_grouped_summary <- function(df, group_vars, summarize_fn, df_name = "data") {
  spliv_sim_assert_columns(df, group_vars, df_name)
  if (nrow(df) == 0L) {
    return(data.frame())
  }

  split_keys <- spliv_sim_group_keys(df, group_vars, df_name = df_name)
  spliv_sim_rbind_fill(lapply(split(df, split_keys, drop = TRUE), summarize_fn))
}

summarize_patterned_family <- function(results, paths, config) {
  path_rows <- spliv_sim_rbind_fill(lapply(results, function(x) x$outputs$path_rows %||% NULL))
  if (nrow(path_rows) == 0L) {
    return(list(path_rows = data.frame(), summary = data.frame()))
  }

  spliv_sim_assert_columns(
    path_rows,
    c(
      "profile", "scenario_name", "truth_pattern", "pattern_label", "delta",
      "contains_true_beta", "contains_zero", "interval_width", "is_true_delta",
      "tipping_point", "false_robustness", "false_fragility"
    ),
    df_name = "path_rows"
  )

  path_rows$contains_true_beta_num <- suppressWarnings(as.numeric(path_rows$contains_true_beta))
  path_rows$contains_zero_num <- suppressWarnings(as.numeric(path_rows$contains_zero))

  summary_df <- spliv_sim_grouped_summary(
    path_rows,
    group_vars = c("profile", "scenario_name", "truth_pattern", "pattern_label", "delta"),
    summarize_fn = function(df) {
      data.frame(
        profile = df$profile[1],
        scenario_name = df$scenario_name[1],
        truth_pattern = df$truth_pattern[1],
        pattern_label = df$pattern_label[1],
        delta = df$delta[1],
        coverage_rate = spliv_sim_safe_mean(df$contains_true_beta_num),
        mean_interval_width = spliv_sim_safe_mean(df$interval_width),
        zero_inclusion_rate = spliv_sim_safe_mean(df$contains_zero_num),
        stringsAsFactors = FALSE
      )
    },
    df_name = "path_rows"
  )

  true_delta_rows <- path_rows[!is.na(path_rows$is_true_delta) & path_rows$is_true_delta, , drop = FALSE]
  if (nrow(true_delta_rows) > 0L) {
    spliv_sim_assert_columns(
      true_delta_rows,
      c(
        "profile", "scenario_name", "truth_pattern", "pattern_label",
        "tipping_point", "false_robustness", "false_fragility"
      ),
      df_name = "true_delta_rows"
    )

    true_delta_rows$false_robustness_num <- suppressWarnings(as.numeric(true_delta_rows$false_robustness))
    true_delta_rows$false_fragility_num <- suppressWarnings(as.numeric(true_delta_rows$false_fragility))

    delta_summary <- spliv_sim_grouped_summary(
      true_delta_rows,
      group_vars = c("profile", "scenario_name", "truth_pattern", "pattern_label"),
      summarize_fn = function(df) {
        data.frame(
          profile = df$profile[1],
          scenario_name = df$scenario_name[1],
          truth_pattern = df$truth_pattern[1],
          pattern_label = df$pattern_label[1],
          mean_tipping_point = spliv_sim_safe_mean(df$tipping_point),
          false_robustness_rate = spliv_sim_safe_mean(df$false_robustness_num),
          false_fragility_rate = spliv_sim_safe_mean(df$false_fragility_num),
          stringsAsFactors = FALSE
        )
      },
      df_name = "true_delta_rows"
    )

    summary_df <- merge(
      summary_df,
      delta_summary,
      by = c("profile", "scenario_name", "truth_pattern", "pattern_label"),
      all.x = TRUE
    )
  }

  out_file <- file.path(paths$tables_dir, sprintf("%s_patterned_summary.csv", config$profile))
  utils::write.csv(summary_df, out_file, row.names = FALSE)

  plot_patterned_paths(
    path_rows,
    file = file.path(paths$figures_dir, sprintf("%s_sensitivity_paths.png", config$profile)),
    title = sprintf("Patterned Sensitivity Paths (%s)", config$profile)
  )
  plot_patterned_metric(
    summary_df,
    value_col = "coverage_rate",
    file = file.path(paths$figures_dir, sprintf("%s_coverage_by_delta.png", config$profile)),
    title = sprintf("Coverage by Delta (%s)", config$profile),
    ylab = "Coverage"
  )
  plot_patterned_metric(
    summary_df,
    value_col = "mean_interval_width",
    file = file.path(paths$figures_dir, sprintf("%s_interval_width_by_delta.png", config$profile)),
    title = sprintf("Interval Width by Delta (%s)", config$profile),
    ylab = "Mean interval width"
  )

  list(path_rows = path_rows, summary = summary_df, file = out_file)
}

summarize_bpe_family <- function(results, paths, config) {
  validation_rows <- spliv_sim_rbind_fill(lapply(results, function(x) x$outputs$validation_summary %||% NULL))
  bpe_rows <- spliv_sim_rbind_fill(lapply(results, function(x) x$outputs$bpe_summary %||% NULL))
  if (nrow(validation_rows) == 0L && nrow(bpe_rows) == 0L) {
    return(list(summary = data.frame()))
  }

  if (nrow(validation_rows) > 0L) {
    spliv_sim_assert_columns(
      validation_rows,
      c(
        "scenario_id", "scenario_name", "replicate_id", "beta_true",
        "pi_active", "inactive_share", "theta_inactive", "transport_gap",
        "profile", "eligibility_passed", "equivalence_passed",
        "first_stage_coefficient", "first_stage_f_statistic",
        "gamma_hat", "gamma_var"
      ),
      df_name = "validation_rows"
    )
  }
  if (nrow(bpe_rows) > 0L) {
    spliv_sim_assert_columns(
      bpe_rows,
      c(
        "scenario_id", "scenario_name", "replicate_id", "beta_true",
        "pi_active", "inactive_share", "theta_inactive", "transport_gap",
        "profile", "fit_ok", "contains_true_beta", "interval_width"
      ),
      df_name = "bpe_rows"
    )
    bpe_rows$contains_true_beta_num <- suppressWarnings(as.numeric(bpe_rows$contains_true_beta))
  }

  merged <- merge(
    validation_rows,
    bpe_rows,
    by = c(
      "scenario_id", "scenario_name", "replicate_id", "beta_true",
      "pi_active", "inactive_share", "theta_inactive", "transport_gap",
      "profile"
    ),
    all = TRUE,
    suffixes = c("_validation", "_bpe")
  )

  spliv_sim_assert_columns(
    merged,
    c(
      "profile", "scenario_name", "pi_active", "inactive_share",
      "theta_inactive", "transport_gap", "eligibility_passed_validation",
      "equivalence_passed_validation", "fit_ok", "first_stage_coefficient",
      "first_stage_f_statistic", "gamma_hat_validation", "gamma_var_validation",
      "contains_true_beta_num", "interval_width"
    ),
    df_name = "merged_bpe_rows"
  )

  merged$eligibility_num <- suppressWarnings(as.numeric(merged$eligibility_passed_validation))
  merged$equivalence_num <- suppressWarnings(as.numeric(merged$equivalence_passed_validation))
  merged$fit_ok_num <- suppressWarnings(as.numeric(merged$fit_ok))

  group_vars <- c(
    "profile", "scenario_name", "pi_active", "inactive_share",
    "theta_inactive", "transport_gap"
  )
  summary_df <- spliv_sim_grouped_summary(
    merged,
    group_vars = group_vars,
    summarize_fn = function(df) {
      false_eligibility_failure_rate <- if (isTRUE(df$transport_gap[1] == 0)) {
        spliv_sim_safe_mean(as.numeric(!df$eligibility_passed_validation))
      } else {
        NA_real_
      }

      data.frame(
        profile = df$profile[1],
        scenario_name = df$scenario_name[1],
        pi_active = df$pi_active[1],
        inactive_share = df$inactive_share[1],
        theta_inactive = df$theta_inactive[1],
        transport_gap = df$transport_gap[1],
        eligibility_rate = spliv_sim_safe_mean(df$eligibility_num),
        equivalence_rate = spliv_sim_safe_mean(df$equivalence_num),
        fit_ok_rate = spliv_sim_safe_mean(df$fit_ok_num),
        mean_first_stage = spliv_sim_safe_mean(df$first_stage_coefficient),
        mean_first_stage_f = spliv_sim_safe_mean(df$first_stage_f_statistic),
        mean_gamma_hat = spliv_sim_safe_mean(df$gamma_hat_validation),
        mean_gamma_var = spliv_sim_safe_mean(df$gamma_var_validation),
        coverage_rate = spliv_sim_safe_mean(df$contains_true_beta_num),
        mean_interval_width = spliv_sim_safe_mean(df$interval_width),
        false_eligibility_failure_rate = false_eligibility_failure_rate,
        stringsAsFactors = FALSE
      )
    },
    df_name = "merged_bpe_rows"
  )

  out_file <- file.path(paths$tables_dir, sprintf("%s_bpe_summary.csv", config$profile))
  utils::write.csv(summary_df, out_file, row.names = FALSE)
  plot_bpe_eligibility(
    summary_df,
    file = file.path(paths$figures_dir, sprintf("%s_bpe_eligibility_rates.png", config$profile)),
    title = sprintf("BPE Eligibility Rates (%s)", config$profile)
  )

  list(summary = summary_df, file = out_file, merged = merged)
}

summarize_search_family <- function(results, paths, config) {
  summary_rows <- spliv_sim_rbind_fill(lapply(results, function(x) x$outputs$summary %||% NULL))
  if (nrow(summary_rows) == 0L) {
    return(list(summary = data.frame()))
  }

  spliv_sim_assert_columns(
    summary_rows,
    c(
      "profile", "scenario_name", "K", "subgroup_share", "pi",
      "any_false_f", "any_false_ci", "any_false_equivalence",
      "min_f_stat", "selected_se", "selected_n_S"
    ),
    df_name = "summary_rows"
  )

  summary_rows$any_false_f_num <- suppressWarnings(as.numeric(summary_rows$any_false_f))
  summary_rows$any_false_ci_num <- suppressWarnings(as.numeric(summary_rows$any_false_ci))
  summary_rows$any_false_equivalence_num <- suppressWarnings(as.numeric(summary_rows$any_false_equivalence))

  summary_df <- spliv_sim_grouped_summary(
    summary_rows,
    group_vars = c("profile", "scenario_name", "K", "subgroup_share", "pi"),
    summarize_fn = function(df) {
      data.frame(
        profile = df$profile[1],
        scenario_name = df$scenario_name[1],
        K = df$K[1],
        subgroup_share = df$subgroup_share[1],
        pi = df$pi[1],
        false_f_rate = spliv_sim_safe_mean(df$any_false_f_num),
        false_ci_rate = spliv_sim_safe_mean(df$any_false_ci_num),
        false_equivalence_rate = spliv_sim_safe_mean(df$any_false_equivalence_num),
        mean_min_f_stat = spliv_sim_safe_mean(df$min_f_stat),
        mean_selected_se = spliv_sim_safe_mean(df$selected_se),
        mean_selected_n_S = spliv_sim_safe_mean(df$selected_n_S),
        stringsAsFactors = FALSE
      )
    },
    df_name = "summary_rows"
  )

  out_file <- file.path(paths$tables_dir, sprintf("%s_subgroup_search_failure.csv", config$profile))
  utils::write.csv(summary_df, out_file, row.names = FALSE)
  plot_search_false_rates(
    summary_df,
    file = file.path(paths$figures_dir, sprintf("%s_false_subgroup_search_rates.png", config$profile)),
    title = sprintf("False Subgroup-Search Rates (%s)", config$profile)
  )

  list(summary = summary_df, file = out_file)
}

summarize_search_stress_family <- function(results, paths, config) {
  summary_rows <- spliv_sim_rbind_fill(lapply(results, function(x) x$outputs$summary %||% NULL))
  if (nrow(summary_rows) == 0L) {
    return(list(summary = data.frame()))
  }

  spliv_sim_assert_columns(
    summary_rows,
    c(
      "profile", "scenario_name", "pi", "K", "subgroup_share", "noise_multiplier",
      "replicate_id", "any_false_f", "any_false_ci", "any_false_equivalence",
      "min_f_stat", "selected_group_share", "selected_group_pi_hat",
      "selected_group_pi_se", "selected_group_abs_t"
    ),
    df_name = "summary_rows"
  )

  summary_rows$any_false_f_num <- suppressWarnings(as.numeric(summary_rows$any_false_f))
  summary_rows$any_false_ci_num <- suppressWarnings(as.numeric(summary_rows$any_false_ci))
  summary_rows$any_false_equivalence_num <- suppressWarnings(as.numeric(summary_rows$any_false_equivalence))

  summary_df <- spliv_sim_grouped_summary(
    summary_rows,
    group_vars = c("profile", "scenario_name", "pi", "K", "subgroup_share", "noise_multiplier"),
    summarize_fn = function(df) {
      data.frame(
        profile = df$profile[1],
        scenario_name = df$scenario_name[1],
        pi = df$pi[1],
        K = df$K[1],
        subgroup_share = df$subgroup_share[1],
        noise_multiplier = df$noise_multiplier[1],
        R = nrow(df),
        false_f_rate = spliv_sim_safe_mean(df$any_false_f_num),
        false_ci_rate = spliv_sim_safe_mean(df$any_false_ci_num),
        false_equivalence_rate = spliv_sim_safe_mean(df$any_false_equivalence_num),
        mean_min_F = spliv_sim_safe_mean(df$min_f_stat),
        median_min_F = spliv_sim_safe_median(df$min_f_stat),
        mean_selected_group_share = spliv_sim_safe_mean(df$selected_group_share),
        mean_selected_group_pi_hat = spliv_sim_safe_mean(df$selected_group_pi_hat),
        mean_selected_group_pi_se = spliv_sim_safe_mean(df$selected_group_pi_se),
        mean_selected_group_abs_t = spliv_sim_safe_mean(df$selected_group_abs_t),
        stringsAsFactors = FALSE
      )
    },
    df_name = "summary_rows"
  )

  out_file <- file.path(paths$tables_dir, sprintf("%s_subgroup_search_stress.csv", config$profile))
  utils::write.csv(summary_df, out_file, row.names = FALSE)

  figure_file <- file.path(paths$figures_dir, sprintf("%s_false_subgroup_search_stress.png", config$profile))
  plot_search_stress_false_rates(
    summary_df,
    file = figure_file,
    title = sprintf("Subgroup-Search Stress Test (%s)", config$profile)
  )

  list(summary = summary_df, file = out_file, figure = figure_file)
}
