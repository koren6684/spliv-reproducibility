spliv_sim_formula <- function() {
  stats::as.formula("y ~ x - 1 | z - 1")
}

.spliv_sim_extract_term_row <- function(fit, term = "x") {
  estimates <- fit$estimates
  if (!is.data.frame(estimates) || nrow(estimates) == 0L) {
    stop("`fit$estimates` is unavailable or empty.")
  }
  idx <- match(term, estimates$term)
  if (is.na(idx)) {
    idx <- 1L
  }

  estimate <- if ("estimate" %in% names(estimates)) {
    as.numeric(estimates$estimate[idx])
  } else if (!is.null(fit$beta_hat) && length(fit$beta_hat) >= idx) {
    as.numeric(fit$beta_hat[idx])
  } else {
    NA_real_
  }
  se <- if ("std.error" %in% names(estimates)) {
    as.numeric(estimates$std.error[idx])
  } else if ("std_error" %in% names(estimates)) {
    as.numeric(estimates$std_error[idx])
  } else if ("se" %in% names(estimates)) {
    as.numeric(estimates$se[idx])
  } else {
    NA_real_
  }

  data.frame(
    term = as.character(estimates$term[idx]),
    estimate = estimate,
    conf_low = as.numeric(estimates$conf.low[idx]),
    conf_high = as.numeric(estimates$conf.high[idx]),
    se = se,
    interval_width = as.numeric(estimates$conf.high[idx] - estimates$conf.low[idx]),
    stringsAsFactors = FALSE
  )
}

run_baseline_iv <- function(data,
                            formula = spliv_sim_formula(),
                            vcov = "hc1",
                            cluster = NULL,
                            scale_instrument = "residual_sd") {
  start_time <- proc.time()[["elapsed"]]
  fit <- spliv(
    formula = formula,
    data = data,
    method = "ltz",
    delta = 0,
    vcov = vcov,
    cluster = cluster,
    scale_instrument = scale_instrument
  )
  summary_row <- .spliv_sim_extract_term_row(fit)
  summary_row$runtime_sec <- proc.time()[["elapsed"]] - start_time
  summary_row$method <- "baseline_iv"
  list(fit = fit, summary = summary_row)
}

run_uci_path <- function(data,
                         delta_grid,
                         violation_pattern = NULL,
                         formula = spliv_sim_formula(),
                         vcov = "hc1",
                         cluster = NULL,
                         scale_instrument = "residual_sd") {
  start_time <- proc.time()[["elapsed"]]
  path_obj <- spliv_sensitivity_path(
    formula = formula,
    data = data,
    method = "uci",
    delta_grid = delta_grid,
    violation_pattern = violation_pattern,
    vcov = vcov,
    cluster = cluster,
    scale_instrument = scale_instrument
  )
  path_df <- as.data.frame(path_obj, stringsAsFactors = FALSE)
  tipping <- spliv_tipping_point(path_obj)
  tipping_x <- if ("x" %in% names(tipping)) tipping[["x"]] else as.numeric(tipping[[1]])
  path_df$runtime_sec <- proc.time()[["elapsed"]] - start_time
  path_df$tipping_point <- tipping_x
  list(path = path_df, tipping_point = tipping_x, object = path_obj)
}

run_ltz_path <- function(data,
                         delta_grid,
                         violation_pattern = NULL,
                         formula = spliv_sim_formula(),
                         vcov = "hc1",
                         cluster = NULL,
                         scale_instrument = "residual_sd") {
  start_time <- proc.time()[["elapsed"]]
  path_obj <- spliv_sensitivity_path(
    formula = formula,
    data = data,
    method = "ltz",
    delta_grid = delta_grid,
    violation_pattern = violation_pattern,
    vcov = vcov,
    cluster = cluster,
    scale_instrument = scale_instrument
  )
  path_df <- as.data.frame(path_obj, stringsAsFactors = FALSE)
  tipping <- spliv_tipping_point(path_obj)
  tipping_x <- if ("x" %in% names(tipping)) tipping[["x"]] else as.numeric(tipping[[1]])
  path_df$runtime_sec <- proc.time()[["elapsed"]] - start_time
  path_df$tipping_point <- tipping_x
  list(path = path_df, tipping_point = tipping_x, object = path_obj)
}

run_bpe_validation <- function(data,
                               design,
                               formula = spliv_sim_formula(),
                               vcov = "cluster",
                               cluster = ~ unit_id,
                               bpe_equiv_margin = 0.15,
                               bpe_min_n_S = 20,
                               bpe_min_clusters_S = 4,
                               scale_instrument = "residual_sd",
                               bpe_transport = "sampling",
                               bpe_transport_kappa = 0) {
  start_time <- proc.time()[["elapsed"]]
  out <- tryCatch(
    bpe_validate_design(
      formula = formula,
      data = data,
      design = design,
      vcov = vcov,
      cluster = cluster,
      bpe_min_n_S = bpe_min_n_S,
      bpe_min_clusters_S = bpe_min_clusters_S,
      bpe_equiv_margin = bpe_equiv_margin,
      scale_instrument = scale_instrument,
      bpe_transport = bpe_transport,
      bpe_transport_kappa = bpe_transport_kappa
    ),
    error = function(e) e
  )

  if (inherits(out, "error")) {
    summary_row <- data.frame(
      ok = FALSE,
      eligibility_passed = NA,
      equivalence_passed = NA,
      n_S = NA_real_,
      share_S = NA_real_,
      G_S = NA_real_,
      varZ_S = NA_real_,
      residualized_instrument_sd_S = NA_real_,
      residualized_treatment_sd_S = NA_real_,
      first_stage_coefficient = NA_real_,
      first_stage_se = NA_real_,
      first_stage_ci_low = NA_real_,
      first_stage_ci_high = NA_real_,
      first_stage_f_statistic = NA_real_,
      first_stage_effect_one_residual_sd_Z = NA_real_,
      standardized_first_stage_effect = NA_real_,
      gamma_hat = NA_real_,
      gamma_var = NA_real_,
      transport_mode = NA_character_,
      transport_uncertainty_inflation = NA_real_,
      runtime_sec = proc.time()[["elapsed"]] - start_time,
      error = conditionMessage(out),
      stringsAsFactors = FALSE
    )
    return(list(validation = NULL, summary = summary_row, error = conditionMessage(out)))
  }

  inst_name <- out$instrument %||% names(out$first_stage_coefficient)[1]
  x_name <- out$first_stage_target %||% names(out$residualized_treatment_sd_S)[1]
  ci_low <- out$first_stage_ci[inst_name, "lower"]
  ci_high <- out$first_stage_ci[inst_name, "upper"]
  gamma_cov <- out$reduced_form_direct_effect_cov
  summary_row <- data.frame(
    ok = TRUE,
    eligibility_passed = isTRUE(out$eligibility_passed),
    equivalence_passed = isTRUE(out$equivalence_passed),
    n_S = out$n_S,
    share_S = out$share_S,
    G_S = out$G_S %||% NA_real_,
    varZ_S = unname(out$varZ_S[inst_name]),
    residualized_instrument_sd_S = unname(out$residualized_instrument_sd_S[inst_name]),
    residualized_treatment_sd_S = unname(out$residualized_treatment_sd_S[x_name]),
    first_stage_coefficient = unname(out$first_stage_coefficient[inst_name]),
    first_stage_se = unname(out$first_stage_se[inst_name]),
    first_stage_ci_low = unname(ci_low),
    first_stage_ci_high = unname(ci_high),
    first_stage_f_statistic = unname(out$first_stage_f_statistic[inst_name]),
    first_stage_effect_one_residual_sd_Z = unname(out$first_stage_effect_one_residual_sd_Z[inst_name]),
    standardized_first_stage_effect = unname(out$standardized_first_stage_effect[inst_name]),
    gamma_hat = unname(out$reduced_form_direct_effect[inst_name]),
    gamma_var = if (is.matrix(gamma_cov)) gamma_cov[1, 1] else NA_real_,
    transport_mode = out$transport_mode %||% NA_character_,
    transport_uncertainty_inflation = out$transport_uncertainty_inflation %||% NA_real_,
    runtime_sec = proc.time()[["elapsed"]] - start_time,
    error = NA_character_,
    stringsAsFactors = FALSE
  )
  list(validation = out, summary = summary_row, error = NULL)
}

run_bpe_confirmatory <- function(data,
                                 design,
                                 formula = spliv_sim_formula(),
                                 vcov = "cluster",
                                 cluster = ~ unit_id,
                                 bpe_equiv_margin = 0.15,
                                 bpe_min_n_S = 20,
                                 bpe_min_clusters_S = 4,
                                 scale_instrument = "residual_sd",
                                 bpe_transport = "sampling",
                                 bpe_transport_kappa = 0) {
  start_time <- proc.time()[["elapsed"]]
  fit <- tryCatch(
    spliv(
      formula = formula,
      data = data,
      method = "bpe",
      bpe_design = design,
      vcov = vcov,
      cluster = cluster,
      bpe_equiv_margin = bpe_equiv_margin,
      bpe_min_n_S = bpe_min_n_S,
      bpe_min_clusters_S = bpe_min_clusters_S,
      scale_instrument = scale_instrument,
      bpe_transport = bpe_transport,
      bpe_transport_kappa = bpe_transport_kappa,
      bpe_not_applicable = "na"
    ),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    summary_row <- data.frame(
      fit_ok = FALSE,
      estimate = NA_real_,
      conf_low = NA_real_,
      conf_high = NA_real_,
      interval_width = NA_real_,
      eligibility_passed = NA,
      equivalence_passed = NA,
      n_S = NA_real_,
      G_S = NA_real_,
      gamma_hat = NA_real_,
      gamma_var = NA_real_,
      runtime_sec = proc.time()[["elapsed"]] - start_time,
      error = conditionMessage(fit),
      stringsAsFactors = FALSE
    )
    return(list(fit = NULL, summary = summary_row, error = conditionMessage(fit)))
  }

  diag <- fit$bpe_diagnostics %||% list()
  row <- .spliv_sim_extract_term_row(fit)
  gamma_cov <- diag$reduced_form_direct_effect_cov
  inst_name <- diag$instrument %||% names(diag$reduced_form_direct_effect)[1]

  summary_row <- data.frame(
    fit_ok = all(is.finite(c(row$conf_low, row$conf_high))),
    estimate = row$estimate,
    conf_low = row$conf_low,
    conf_high = row$conf_high,
    interval_width = row$interval_width,
    eligibility_passed = diag$eligibility_passed %||% NA,
    equivalence_passed = diag$equivalence_passed %||% NA,
    n_S = diag$n_S %||% NA_real_,
    G_S = diag$G_S %||% NA_real_,
    gamma_hat = if (!is.null(diag$reduced_form_direct_effect)) {
      unname(diag$reduced_form_direct_effect[inst_name])
    } else {
      NA_real_
    },
    gamma_var = if (is.matrix(gamma_cov)) gamma_cov[1, 1] else NA_real_,
    runtime_sec = proc.time()[["elapsed"]] - start_time,
    error = NA_character_,
    stringsAsFactors = FALSE
  )
  list(fit = fit, summary = summary_row, error = NULL)
}

run_subgroup_search_diagnostic <- function(data,
                                           formula = spliv_sim_formula(),
                                           K = 20L,
                                           subgroup_share = 0.10,
                                           seed = 1L,
                                           vcov = "hc1",
                                           cluster = NULL,
                                           bpe_equiv_margin = 0.15,
                                           bpe_min_n_S = 20) {
  start_time <- proc.time()[["elapsed"]]
  candidate_data <- spliv_sim_candidate_groups(
    data = data,
    K = K,
    subgroup_share = subgroup_share,
    seed = seed
  )

  candidate_rows <- vector("list", length(candidate_data$candidate_names))
  for (idx in seq_along(candidate_data$candidate_names)) {
    candidate_name <- candidate_data$candidate_names[idx]
    design <- bpe_design(
      name = candidate_name,
      subset = candidate_name,
      rationale = "Exploratory diagnostic subgroup only; not valid confirmatory BPE.",
      variables_used = candidate_name,
      subset_type = "exploratory_diagnostic",
      pre_specified = TRUE
    )

    validation <- tryCatch(
      bpe_validate_design(
        formula = formula,
        data = candidate_data$data,
        design = design,
        vcov = vcov,
        cluster = cluster,
        bpe_min_n_S = bpe_min_n_S,
        bpe_equiv_margin = bpe_equiv_margin,
        scale_instrument = "residual_sd"
      ),
      error = function(e) e
    )

    if (inherits(validation, "error")) {
      candidate_rows[[idx]] <- data.frame(
        candidate = candidate_name,
        n_S = NA_real_,
        share_S = NA_real_,
        first_stage_coefficient = NA_real_,
        first_stage_se = NA_real_,
        ci_low = NA_real_,
        ci_high = NA_real_,
        ci_includes_zero = NA,
        first_stage_f_statistic = NA_real_,
        equivalence_passed = NA,
        eligibility_passed = NA,
        error = conditionMessage(validation),
        stringsAsFactors = FALSE
      )
      next
    }

    inst_name <- validation$instrument %||% names(validation$first_stage_coefficient)[1]
    ci_low <- validation$first_stage_ci[inst_name, "lower"]
    ci_high <- validation$first_stage_ci[inst_name, "upper"]

    candidate_rows[[idx]] <- data.frame(
      candidate = candidate_name,
      n_S = validation$n_S,
      share_S = validation$share_S,
      first_stage_coefficient = unname(validation$first_stage_coefficient[inst_name]),
      first_stage_se = unname(validation$first_stage_se[inst_name]),
      ci_low = unname(ci_low),
      ci_high = unname(ci_high),
      ci_includes_zero = isTRUE(ci_low <= 0 && ci_high >= 0),
      first_stage_f_statistic = unname(validation$first_stage_f_statistic[inst_name]),
      equivalence_passed = isTRUE(validation$equivalence_passed),
      eligibility_passed = isTRUE(validation$eligibility_passed),
      error = NA_character_,
      stringsAsFactors = FALSE
    )
  }

  candidate_df <- spliv_sim_rbind_fill(candidate_rows)
  valid_df <- candidate_df[is.na(candidate_df$error) & is.finite(candidate_df$first_stage_f_statistic), , drop = FALSE]
  selected_row <- if (nrow(valid_df) > 0L) {
    valid_df[which.min(valid_df$first_stage_f_statistic), , drop = FALSE]
  } else {
    data.frame(
      candidate = NA_character_,
      n_S = NA_real_,
      share_S = NA_real_,
      first_stage_coefficient = NA_real_,
      first_stage_se = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      ci_includes_zero = NA,
      first_stage_f_statistic = NA_real_,
      equivalence_passed = NA,
      eligibility_passed = NA,
      error = NA_character_,
      stringsAsFactors = FALSE
    )
  }

  summary_row <- data.frame(
    any_false_f = if (nrow(valid_df) > 0L) any(valid_df$first_stage_f_statistic <= 5, na.rm = TRUE) else FALSE,
    any_false_ci = if (nrow(valid_df) > 0L) any(valid_df$ci_includes_zero, na.rm = TRUE) else FALSE,
    any_false_equivalence = if (nrow(valid_df) > 0L) any(valid_df$equivalence_passed, na.rm = TRUE) else FALSE,
    min_f_stat = if (nrow(valid_df) > 0L) min(valid_df$first_stage_f_statistic, na.rm = TRUE) else NA_real_,
    selected_candidate = selected_row$candidate[1],
    selected_pi_hat = selected_row$first_stage_coefficient[1],
    selected_se = selected_row$first_stage_se[1],
    selected_n_S = selected_row$n_S[1],
    selected_share_S = selected_row$share_S[1],
    selected_equivalence = selected_row$equivalence_passed[1],
    n_candidates_valid = nrow(valid_df),
    runtime_sec = proc.time()[["elapsed"]] - start_time,
    stringsAsFactors = FALSE
  )
  list(summary = summary_row, candidates = candidate_df)
}

run_subgroup_search_stress_diagnostic <- function(data,
                                                  K = 20L,
                                                  subgroup_share = 0.10,
                                                  seed = 1L,
                                                  equivalence_margin = 0.25) {
  start_time <- proc.time()[["elapsed"]]
  required_cols <- c("unit_id", "x", "z")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0L) {
    stop(
      "Stress-test subgroup search requires columns: ",
      paste(required_cols, collapse = ", "),
      ". Missing: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  group_info <- spliv_sim_candidate_group_matrix(
    unit_id = data$unit_id,
    K = K,
    subgroup_share = subgroup_share,
    seed = seed
  )
  group_factor <- factor(as.character(data$unit_id), levels = group_info$unit_levels)

  z2_unit <- as.numeric(rowsum(data$z^2, group_factor, reorder = FALSE))
  zx_unit <- as.numeric(rowsum(data$z * data$x, group_factor, reorder = FALSE))
  x2_unit <- as.numeric(rowsum(data$x^2, group_factor, reorder = FALSE))
  n_obs_unit <- as.numeric(rowsum(rep.int(1L, nrow(data)), group_factor, reorder = FALSE))

  membership <- group_info$membership
  group_obs_counts <- as.numeric(crossprod(membership, n_obs_unit))
  group_unit_counts <- as.numeric(colSums(membership))
  szz <- as.numeric(crossprod(membership, z2_unit))
  szx <- as.numeric(crossprod(membership, zx_unit))
  sxx <- as.numeric(crossprod(membership, x2_unit))

  pi_hat <- rep(NA_real_, K)
  pi_se <- rep(NA_real_, K)
  ci_low <- rep(NA_real_, K)
  ci_high <- rep(NA_real_, K)
  abs_t <- rep(NA_real_, K)
  f_stat <- rep(NA_real_, K)

  valid <- is.finite(szz) & szz > 0 & is.finite(sxx) & group_obs_counts > 1
  if (any(valid)) {
    pi_hat[valid] <- szx[valid] / szz[valid]
    sse <- rep(NA_real_, K)
    sse[valid] <- pmax(sxx[valid] - (szx[valid]^2 / szz[valid]), 0)
    df <- pmax(group_obs_counts - 1, 1)
    sigma2 <- rep(NA_real_, K)
    sigma2[valid] <- sse[valid] / df[valid]
    positive_se <- valid & is.finite(sigma2) & sigma2 >= 0
    pi_se[positive_se] <- sqrt(sigma2[positive_se] / szz[positive_se])
    positive_se <- positive_se & is.finite(pi_se) & pi_se > 0
    crit <- rep(NA_real_, K)
    crit[positive_se] <- stats::qt(0.975, df[positive_se])
    ci_low[positive_se] <- pi_hat[positive_se] - crit[positive_se] * pi_se[positive_se]
    ci_high[positive_se] <- pi_hat[positive_se] + crit[positive_se] * pi_se[positive_se]
    abs_t[positive_se] <- abs(pi_hat[positive_se]) / pi_se[positive_se]
    f_stat[positive_se] <- abs_t[positive_se]^2
  }

  group_share <- group_obs_counts / nrow(data)
  ci_includes_zero <- is.finite(ci_low) & is.finite(ci_high) & ci_low <= 0 & ci_high >= 0
  equivalence_pass <- is.finite(ci_low) & is.finite(ci_high) &
    ci_low >= -equivalence_margin & ci_high <= equivalence_margin

  selected_idx <- if (any(is.finite(f_stat))) {
    which.min(ifelse(is.finite(f_stat), f_stat, Inf))
  } else {
    NA_integer_
  }

  selected_group_name <- if (is.na(selected_idx)) NA_character_ else colnames(membership)[selected_idx]
  selected_group_size <- if (is.na(selected_idx)) NA_real_ else group_obs_counts[selected_idx]
  selected_group_units <- if (is.na(selected_idx)) NA_real_ else group_unit_counts[selected_idx]
  selected_group_share <- if (is.na(selected_idx)) NA_real_ else group_share[selected_idx]
  selected_group_pi_hat <- if (is.na(selected_idx)) NA_real_ else pi_hat[selected_idx]
  selected_group_pi_se <- if (is.na(selected_idx)) NA_real_ else pi_se[selected_idx]
  selected_group_abs_t <- if (is.na(selected_idx)) NA_real_ else abs_t[selected_idx]

  summary_row <- data.frame(
    any_false_f = any(f_stat <= 5, na.rm = TRUE),
    any_false_ci = any(ci_includes_zero, na.rm = TRUE),
    any_false_equivalence = any(equivalence_pass, na.rm = TRUE),
    n_groups_passing_f = sum(f_stat <= 5, na.rm = TRUE),
    n_groups_passing_ci = sum(ci_includes_zero, na.rm = TRUE),
    n_groups_passing_equivalence = sum(equivalence_pass, na.rm = TRUE),
    min_f_stat = if (any(is.finite(f_stat))) min(f_stat, na.rm = TRUE) else NA_real_,
    selected_group = selected_group_name,
    selected_group_size = selected_group_size,
    selected_group_units = selected_group_units,
    selected_group_share = selected_group_share,
    selected_group_pi_hat = selected_group_pi_hat,
    selected_group_pi_se = selected_group_pi_se,
    selected_group_abs_t = selected_group_abs_t,
    equivalence_margin = equivalence_margin,
    runtime_sec = proc.time()[["elapsed"]] - start_time,
    stringsAsFactors = FALSE
  )

  list(summary = summary_row)
}

.spliv_sim_result_meta <- function(family,
                                   scenario,
                                   replicate_id,
                                   seed,
                                   config,
                                   package_info,
                                   runtime_sec,
                                   ok,
                                   error = NULL) {
  list(
    family = family,
    scenario = as.list(scenario),
    replicate_id = replicate_id,
    seed = seed,
    profile = config$profile,
    config = list(
      R = config$R,
      nx = config$nx,
      ny = config$ny,
      T_periods = config$T_periods,
      delta_grid = config$delta_grid,
      n_cores = config$n_cores,
      overwrite = config$overwrite,
      base_seed = config$base_seed
    ),
    package_info = package_info,
    runtime_sec = runtime_sec,
    timestamp = spliv_sim_now(),
    ok = ok,
    error = error
  )
}

.spliv_sim_annotate_path <- function(path_df,
                                     pattern_label,
                                     scenario,
                                     truth,
                                     replicate_id,
                                     config) {
  path_df <- path_df[path_df$term == "x", , drop = FALSE]
  true_delta <- spliv_sim_nearest_delta(config$delta_grid, truth$theta_true %||% 0)
  path_df$pattern_label <- pattern_label
  path_df$scenario_id <- scenario$scenario_id
  path_df$scenario_name <- scenario$scenario_name
  path_df$truth_pattern <- scenario$truth_pattern
  path_df$beta_true <- truth$beta_true
  path_df$theta_true <- truth$theta_true
  path_df$replicate_id <- replicate_id
  path_df$profile <- config$profile
  path_df$true_delta <- true_delta
  path_df$is_true_delta <- abs(path_df$delta - true_delta) < 1e-8
  path_df$contains_true_beta <- path_df$conf_low <= path_df$beta_true &
    path_df$conf_high >= path_df$beta_true
  path_df$interval_width <- path_df$conf_high - path_df$conf_low
  path_df$false_robustness <- ifelse(path_df$is_true_delta, !path_df$contains_true_beta, NA)
  path_df$false_fragility <- ifelse(
    path_df$is_true_delta,
    path_df$contains_zero & path_df$beta_true > 0,
    NA
  )
  path_df
}

run_patterned_family <- function(config, paths, package_info, log_file = NULL) {
  scenarios <- make_patterned_scenarios(config)
  tasks <- spliv_sim_make_tasks(scenarios, config$R)

  if (!is.null(log_file)) {
    spliv_sim_log(log_file, "Starting patterned sensitivity family with ", length(tasks), " tasks.")
  }

  worker <- function(task) {
    scenario <- scenarios[task$scenario_idx, , drop = FALSE]
    result_file <- spliv_sim_result_file(
      config = config,
      family = "patterned",
      scenario_name = scenario$scenario_name,
      replicate_id = task$replicate_id
    )
    if (file.exists(result_file) && !isTRUE(config$overwrite)) {
      return(list(status = "skipped", file = result_file))
    }

    seed <- spliv_sim_task_seed(config$base_seed, scenario$scenario_id, task$replicate_id)
    start_time <- proc.time()[["elapsed"]]
    result <- tryCatch({
      sim <- simulate_patterned_panel(
        nx = config$nx,
        ny = config$ny,
        T_periods = config$T_periods,
        beta = 1,
        pi = scenario$pi,
        theta_true = scenario$theta_true,
        truth_pattern = scenario$truth_pattern,
        seed = seed
      )
      patterns <- spliv_sim_pattern_objects(sim$data_est)
      incorrect_key <- spliv_sim_incorrect_pattern_key(scenario$truth_pattern)

      baseline <- run_baseline_iv(sim$data_est, vcov = "hc1")
      default_path <- run_uci_path(sim$data_est, config$delta_grid, violation_pattern = NULL, vcov = "hc1")
      uniform_path <- run_uci_path(sim$data_est, config$delta_grid, violation_pattern = patterns$uniform, vcov = "hc1")
      correct_path <- run_uci_path(sim$data_est, config$delta_grid, violation_pattern = patterns[[scenario$truth_pattern]], vcov = "hc1")
      incorrect_path <- run_uci_path(sim$data_est, config$delta_grid, violation_pattern = patterns[[incorrect_key]], vcov = "hc1")

      path_rows <- spliv_sim_rbind_fill(list(
        .spliv_sim_annotate_path(default_path$path, "default", scenario, sim$truth, task$replicate_id, config),
        .spliv_sim_annotate_path(uniform_path$path, "uniform", scenario, sim$truth, task$replicate_id, config),
        .spliv_sim_annotate_path(correct_path$path, "correct", scenario, sim$truth, task$replicate_id, config),
        .spliv_sim_annotate_path(incorrect_path$path, "incorrect", scenario, sim$truth, task$replicate_id, config)
      ))

      baseline_row <- baseline$summary
      baseline_row$scenario_id <- scenario$scenario_id
      baseline_row$scenario_name <- scenario$scenario_name
      baseline_row$truth_pattern <- scenario$truth_pattern
      baseline_row$beta_true <- sim$truth$beta_true
      baseline_row$theta_true <- sim$truth$theta_true
      baseline_row$replicate_id <- task$replicate_id
      baseline_row$profile <- config$profile

      list(
        meta = .spliv_sim_result_meta(
          family = "patterned",
          scenario = scenario,
          replicate_id = task$replicate_id,
          seed = seed,
          config = config,
          package_info = package_info,
          runtime_sec = proc.time()[["elapsed"]] - start_time,
          ok = TRUE
        ),
        outputs = list(
          baseline_summary = baseline_row,
          path_rows = path_rows
        )
      )
    }, error = function(e) {
      list(
        meta = .spliv_sim_result_meta(
          family = "patterned",
          scenario = scenario,
          replicate_id = task$replicate_id,
          seed = seed,
          config = config,
          package_info = package_info,
          runtime_sec = proc.time()[["elapsed"]] - start_time,
          ok = FALSE,
          error = conditionMessage(e)
        ),
        outputs = list()
      )
    })

    spliv_sim_save_result(result, result_file)
    list(
      status = if (isTRUE(result$meta$ok)) "completed" else "failed",
      file = result_file,
      error = result$meta$error %||% NA_character_
    )
  }

  status_rows <- spliv_sim_apply(tasks, worker, n_cores = config$n_cores)
  counts <- spliv_sim_count_status(status_rows)
  if (!is.null(log_file)) {
    spliv_sim_log(log_file, "Patterned family complete.")
    spliv_sim_log_footer(log_file, counts)
  }

  list(
    family = "patterned",
    scenarios = scenarios,
    status = status_rows,
    counts = counts,
    results = spliv_sim_collect_results(config, "patterned")
  )
}

run_bpe_family <- function(config, paths, package_info, log_file = NULL) {
  scenarios <- make_bpe_scenarios(config)
  tasks <- spliv_sim_make_tasks(scenarios, config$R)

  if (!is.null(log_file)) {
    spliv_sim_log(log_file, "Starting confirmatory BPE family with ", length(tasks), " tasks.")
  }

  worker <- function(task) {
    scenario <- scenarios[task$scenario_idx, , drop = FALSE]
    result_file <- spliv_sim_result_file(
      config = config,
      family = "bpe",
      scenario_name = scenario$scenario_name,
      replicate_id = task$replicate_id
    )
    if (file.exists(result_file) && !isTRUE(config$overwrite)) {
      return(list(status = "skipped", file = result_file))
    }

    seed <- spliv_sim_task_seed(config$base_seed, scenario$scenario_id, task$replicate_id)
    start_time <- proc.time()[["elapsed"]]
    result <- tryCatch({
      sim <- simulate_bpe_panel(
        nx = config$nx,
        ny = config$ny,
        T_periods = config$T_periods,
        beta = 1,
        pi_active = scenario$pi_active,
        theta_inactive = scenario$theta_inactive,
        inactive_share = scenario$inactive_share,
        transport_gap = scenario$transport_gap,
        seed = seed
      )
      design <- spliv_sim_bpe_design()
      min_n_s <- max(20L, floor(0.05 * nrow(sim$data_est)))
      min_clusters <- max(4L, floor(config$nx * config$ny * scenario$inactive_share * 0.4))

      baseline <- run_baseline_iv(
        sim$data_est,
        vcov = "cluster",
        cluster = ~ unit_id
      )
      uci_path <- run_uci_path(
        sim$data_est,
        delta_grid = config$delta_grid,
        vcov = "cluster",
        cluster = ~ unit_id
      )
      validation <- run_bpe_validation(
        sim$data_est,
        design = design,
        vcov = "cluster",
        cluster = ~ unit_id,
        bpe_equiv_margin = 0.15,
        bpe_min_n_S = min_n_s,
        bpe_min_clusters_S = min_clusters,
        bpe_transport = "sampling"
      )
      bpe_fit <- run_bpe_confirmatory(
        sim$data_est,
        design = design,
        vcov = "cluster",
        cluster = ~ unit_id,
        bpe_equiv_margin = 0.15,
        bpe_min_n_S = min_n_s,
        bpe_min_clusters_S = min_clusters,
        bpe_transport = "sampling"
      )

      baseline_row <- baseline$summary
      baseline_row$scenario_id <- scenario$scenario_id
      baseline_row$scenario_name <- scenario$scenario_name
      baseline_row$replicate_id <- task$replicate_id
      baseline_row$beta_true <- sim$truth$beta_true
      baseline_row$profile <- config$profile

      validation_row <- validation$summary
      validation_row$scenario_id <- scenario$scenario_id
      validation_row$scenario_name <- scenario$scenario_name
      validation_row$replicate_id <- task$replicate_id
      validation_row$beta_true <- sim$truth$beta_true
      validation_row$pi_active <- scenario$pi_active
      validation_row$inactive_share <- scenario$inactive_share
      validation_row$theta_inactive <- scenario$theta_inactive
      validation_row$transport_gap <- scenario$transport_gap
      validation_row$profile <- config$profile

      bpe_row <- bpe_fit$summary
      bpe_row$scenario_id <- scenario$scenario_id
      bpe_row$scenario_name <- scenario$scenario_name
      bpe_row$replicate_id <- task$replicate_id
      bpe_row$beta_true <- sim$truth$beta_true
      bpe_row$pi_active <- scenario$pi_active
      bpe_row$inactive_share <- scenario$inactive_share
      bpe_row$theta_inactive <- scenario$theta_inactive
      bpe_row$transport_gap <- scenario$transport_gap
      bpe_row$profile <- config$profile
      bpe_row$contains_true_beta <- with(
        bpe_row,
        is.finite(conf_low) && is.finite(conf_high) &&
          conf_low <= beta_true && conf_high >= beta_true
      )

      uci_rows <- uci_path$path
      uci_rows <- uci_rows[uci_rows$term == "x", , drop = FALSE]
      uci_rows$scenario_id <- scenario$scenario_id
      uci_rows$scenario_name <- scenario$scenario_name
      uci_rows$replicate_id <- task$replicate_id
      uci_rows$beta_true <- sim$truth$beta_true
      uci_rows$theta_inactive <- scenario$theta_inactive
      uci_rows$transport_gap <- scenario$transport_gap
      uci_rows$profile <- config$profile
      uci_rows$contains_true_beta <- uci_rows$conf_low <= uci_rows$beta_true &
        uci_rows$conf_high >= uci_rows$beta_true
      uci_rows$interval_width <- uci_rows$conf_high - uci_rows$conf_low

      list(
        meta = .spliv_sim_result_meta(
          family = "bpe",
          scenario = scenario,
          replicate_id = task$replicate_id,
          seed = seed,
          config = config,
          package_info = package_info,
          runtime_sec = proc.time()[["elapsed"]] - start_time,
          ok = TRUE
        ),
        outputs = list(
          baseline_summary = baseline_row,
          validation_summary = validation_row,
          bpe_summary = bpe_row,
          uci_path_rows = uci_rows
        )
      )
    }, error = function(e) {
      list(
        meta = .spliv_sim_result_meta(
          family = "bpe",
          scenario = scenario,
          replicate_id = task$replicate_id,
          seed = seed,
          config = config,
          package_info = package_info,
          runtime_sec = proc.time()[["elapsed"]] - start_time,
          ok = FALSE,
          error = conditionMessage(e)
        ),
        outputs = list()
      )
    })

    spliv_sim_save_result(result, result_file)
    list(
      status = if (isTRUE(result$meta$ok)) "completed" else "failed",
      file = result_file,
      error = result$meta$error %||% NA_character_
    )
  }

  status_rows <- spliv_sim_apply(tasks, worker, n_cores = config$n_cores)
  counts <- spliv_sim_count_status(status_rows)
  if (!is.null(log_file)) {
    spliv_sim_log(log_file, "BPE family complete.")
    spliv_sim_log_footer(log_file, counts)
  }

  list(
    family = "bpe",
    scenarios = scenarios,
    status = status_rows,
    counts = counts,
    results = spliv_sim_collect_results(config, "bpe")
  )
}

run_search_family <- function(config, paths, package_info, log_file = NULL) {
  scenarios <- make_search_scenarios(config)
  tasks <- spliv_sim_make_tasks(scenarios, config$R)

  if (!is.null(log_file)) {
    spliv_sim_log(log_file, "Starting subgroup-search diagnostic family with ", length(tasks), " tasks.")
  }

  worker <- function(task) {
    scenario <- scenarios[task$scenario_idx, , drop = FALSE]
    result_file <- spliv_sim_result_file(
      config = config,
      family = "search",
      scenario_name = scenario$scenario_name,
      replicate_id = task$replicate_id
    )
    if (file.exists(result_file) && !isTRUE(config$overwrite)) {
      return(list(status = "skipped", file = result_file))
    }

    seed <- spliv_sim_task_seed(config$base_seed, scenario$scenario_id, task$replicate_id)
    start_time <- proc.time()[["elapsed"]]
    result <- tryCatch({
      sim <- simulate_search_panel(
        nx = config$nx,
        ny = config$ny,
        T_periods = config$T_periods,
        beta = 1,
        pi = scenario$pi,
        seed = seed
      )
      diagnostic <- run_subgroup_search_diagnostic(
        sim$data_est,
        K = scenario$K,
        subgroup_share = scenario$subgroup_share,
        seed = seed + 5000L,
        vcov = "hc1",
        bpe_equiv_margin = 0.15,
        bpe_min_n_S = max(20L, floor(0.05 * nrow(sim$data_est)))
      )

      summary_row <- diagnostic$summary
      summary_row$scenario_id <- scenario$scenario_id
      summary_row$scenario_name <- scenario$scenario_name
      summary_row$replicate_id <- task$replicate_id
      summary_row$beta_true <- sim$truth$beta_true
      summary_row$pi <- scenario$pi
      summary_row$K <- scenario$K
      summary_row$subgroup_share <- scenario$subgroup_share
      summary_row$profile <- config$profile

      candidate_rows <- diagnostic$candidates
      candidate_rows$scenario_id <- scenario$scenario_id
      candidate_rows$scenario_name <- scenario$scenario_name
      candidate_rows$replicate_id <- task$replicate_id
      candidate_rows$pi <- scenario$pi
      candidate_rows$K <- scenario$K
      candidate_rows$subgroup_share <- scenario$subgroup_share
      candidate_rows$profile <- config$profile

      list(
        meta = .spliv_sim_result_meta(
          family = "search",
          scenario = scenario,
          replicate_id = task$replicate_id,
          seed = seed,
          config = config,
          package_info = package_info,
          runtime_sec = proc.time()[["elapsed"]] - start_time,
          ok = TRUE
        ),
        outputs = list(
          summary = summary_row,
          candidates = candidate_rows
        )
      )
    }, error = function(e) {
      list(
        meta = .spliv_sim_result_meta(
          family = "search",
          scenario = scenario,
          replicate_id = task$replicate_id,
          seed = seed,
          config = config,
          package_info = package_info,
          runtime_sec = proc.time()[["elapsed"]] - start_time,
          ok = FALSE,
          error = conditionMessage(e)
        ),
        outputs = list()
      )
    })

    spliv_sim_save_result(result, result_file)
    list(
      status = if (isTRUE(result$meta$ok)) "completed" else "failed",
      file = result_file,
      error = result$meta$error %||% NA_character_
    )
  }

  status_rows <- spliv_sim_apply(tasks, worker, n_cores = config$n_cores)
  counts <- spliv_sim_count_status(status_rows)
  if (!is.null(log_file)) {
    spliv_sim_log(log_file, "Subgroup-search diagnostic family complete.")
    spliv_sim_log_footer(log_file, counts)
  }

  list(
    family = "search",
    scenarios = scenarios,
    status = status_rows,
    counts = counts,
    results = spliv_sim_collect_results(config, "search")
  )
}

run_search_stress_family <- function(config, paths, package_info, log_file = NULL) {
  scenarios <- make_search_stress_scenarios(config)
  tasks <- spliv_sim_make_tasks(scenarios, config$R)

  if (!is.null(log_file)) {
    spliv_sim_log(log_file, "Starting subgroup-search stress diagnostic family with ", length(tasks), " tasks.")
  }

  worker <- function(task) {
    scenario <- scenarios[task$scenario_idx, , drop = FALSE]
    result_file <- spliv_sim_result_file(
      config = config,
      family = "subgroup_search_stress",
      scenario_name = scenario$scenario_name,
      replicate_id = task$replicate_id
    )
    if (file.exists(result_file) && !isTRUE(config$overwrite)) {
      return(list(status = "skipped", file = result_file))
    }

    seed <- spliv_sim_task_seed(config$base_seed, scenario$scenario_id, task$replicate_id)
    start_time <- proc.time()[["elapsed"]]
    result <- tryCatch({
      sim <- simulate_search_panel(
        nx = config$nx,
        ny = config$ny,
        T_periods = config$T_periods,
        beta = 1,
        pi = scenario$pi,
        first_stage_noise_multiplier = scenario$noise_multiplier,
        seed = seed
      )
      diagnostic <- run_subgroup_search_stress_diagnostic(
        sim$data_est,
        K = scenario$K,
        subgroup_share = scenario$subgroup_share,
        seed = seed + 7000L,
        equivalence_margin = 0.25
      )

      summary_row <- diagnostic$summary
      summary_row$scenario_id <- scenario$scenario_id
      summary_row$scenario_name <- scenario$scenario_name
      summary_row$replicate_id <- task$replicate_id
      summary_row$beta_true <- sim$truth$beta_true
      summary_row$pi <- scenario$pi
      summary_row$K <- scenario$K
      summary_row$subgroup_share <- scenario$subgroup_share
      summary_row$noise_multiplier <- scenario$noise_multiplier
      summary_row$profile <- config$profile

      list(
        meta = .spliv_sim_result_meta(
          family = "subgroup_search_stress",
          scenario = scenario,
          replicate_id = task$replicate_id,
          seed = seed,
          config = config,
          package_info = package_info,
          runtime_sec = proc.time()[["elapsed"]] - start_time,
          ok = TRUE
        ),
        outputs = list(
          summary = summary_row
        )
      )
    }, error = function(e) {
      list(
        meta = .spliv_sim_result_meta(
          family = "subgroup_search_stress",
          scenario = scenario,
          replicate_id = task$replicate_id,
          seed = seed,
          config = config,
          package_info = package_info,
          runtime_sec = proc.time()[["elapsed"]] - start_time,
          ok = FALSE,
          error = conditionMessage(e)
        ),
        outputs = list()
      )
    })

    spliv_sim_save_result(result, result_file)
    list(
      status = if (isTRUE(result$meta$ok)) "completed" else "failed",
      file = result_file,
      error = result$meta$error %||% NA_character_
    )
  }

  status_rows <- spliv_sim_apply(tasks, worker, n_cores = config$n_cores)
  counts <- spliv_sim_count_status(status_rows)
  if (!is.null(log_file)) {
    spliv_sim_log(log_file, "Subgroup-search stress diagnostic family complete.")
    spliv_sim_log_footer(log_file, counts)
  }

  list(
    family = "subgroup_search_stress",
    scenarios = scenarios,
    status = status_rows,
    counts = counts,
    results = spliv_sim_collect_results(config, "subgroup_search_stress")
  )
}
