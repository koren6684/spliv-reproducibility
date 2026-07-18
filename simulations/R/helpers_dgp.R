spliv_sim_make_grid <- function(nx, ny) {
  grid <- expand.grid(
    x_coord = seq_len(nx),
    y_coord = seq_len(ny)
  )
  grid$unit_id <- seq_len(nrow(grid))
  zone_cut_x <- max(1L, ceiling(nx / 3))
  zone_cut_y <- max(1L, ceiling(ny / 2))
  grid$zone_pattern <- as.numeric(
    grid$x_coord <= zone_cut_x & grid$y_coord <= zone_cut_y
  )
  if (nx <= 1L) {
    grid$exposure_gradient <- 1
  } else {
    grid$exposure_gradient <- (grid$x_coord - 1) / (nx - 1)
  }
  grid
}

spliv_sim_expand_panel <- function(grid_df, T_periods) {
  n_units <- nrow(grid_df)
  out <- grid_df[rep(seq_len(n_units), times = T_periods), , drop = FALSE]
  out$unit_id <- factor(out$unit_id)
  out$period <- factor(rep(seq_len(T_periods), each = n_units))
  rownames(out) <- NULL
  out
}

spliv_sim_make_inactive_condition <- function(panel_df, inactive_share = 0.25) {
  unit_df <- unique(panel_df[, c("unit_id", "exposure_gradient")])
  unit_df <- unit_df[order(unit_df$exposure_gradient, decreasing = TRUE), , drop = FALSE]
  n_inactive <- max(1L, floor(inactive_share * nrow(unit_df)))
  inactive_units <- as.character(unit_df$unit_id[seq_len(n_inactive)])
  as.logical(as.character(panel_df$unit_id) %in% inactive_units)
}

spliv_sim_truth_pattern <- function(panel_df, truth_pattern = c("uniform", "zone", "gradient")) {
  truth_pattern <- match.arg(truth_pattern)
  if (identical(truth_pattern, "uniform")) {
    return(rep(1, nrow(panel_df)))
  }
  if (identical(truth_pattern, "zone")) {
    return(as.numeric(panel_df$zone_pattern))
  }
  as.numeric(panel_df$exposure_gradient)
}

spliv_sim_twoway_demean <- function(x, unit_id, period) {
  x <- as.numeric(x)
  x - ave(x, unit_id, FUN = mean) - ave(x, period, FUN = mean) + mean(x)
}

spliv_sim_residualize_panel <- function(data,
                                        y_var = "y",
                                        x_vars = "x",
                                        z_vars = "z",
                                        unit_var = "unit_id",
                                        period_var = "period",
                                        keep_vars = c(
                                          "unit_id", "period",
                                          "x_coord", "y_coord",
                                          "zone_pattern", "exposure_gradient",
                                          "inactive_condition"
                                        )) {
  keep_vars <- intersect(keep_vars, names(data))
  out <- data[keep_vars]

  out[[y_var]] <- spliv_sim_twoway_demean(data[[y_var]], data[[unit_var]], data[[period_var]])
  for (x_var in x_vars) {
    out[[x_var]] <- spliv_sim_twoway_demean(data[[x_var]], data[[unit_var]], data[[period_var]])
  }
  for (z_var in z_vars) {
    out[[z_var]] <- spliv_sim_twoway_demean(data[[z_var]], data[[unit_var]], data[[period_var]])
  }

  rownames(out) <- NULL
  out
}

.spliv_sim_common_panel <- function(nx, ny, T_periods, seed) {
  set.seed(seed)
  panel <- spliv_sim_expand_panel(spliv_sim_make_grid(nx, ny), T_periods)
  n <- nrow(panel)
  n_units <- nlevels(panel$unit_id)
  n_periods <- nlevels(panel$period)

  panel$inactive_condition <- spliv_sim_make_inactive_condition(panel, inactive_share = 0.25)

  panel$z <- rnorm(n)
  panel$unit_fe_x <- rnorm(n_units, sd = 0.4)[panel$unit_id]
  panel$unit_fe_y <- rnorm(n_units, sd = 0.6)[panel$unit_id]
  panel$time_fe_x <- rnorm(n_periods, sd = 0.3)[panel$period]
  panel$time_fe_y <- rnorm(n_periods, sd = 0.4)[panel$period]

  shared_shock <- rnorm(n)
  panel$v <- 0.6 * shared_shock + sqrt(1 - 0.6^2) * rnorm(n, sd = 0.8)
  panel$e <- shared_shock + rnorm(n, sd = 0.8)

  z_resid <- spliv_sim_twoway_demean(panel$z, panel$unit_id, panel$period)
  z_resid_sd <- stats::sd(z_resid)
  if (!is.finite(z_resid_sd) || z_resid_sd <= 0) {
    z_resid_sd <- 1
  }
  panel$z_scaled_resid <- z_resid / z_resid_sd
  panel
}

simulate_patterned_panel <- function(nx = 8L,
                                     ny = 8L,
                                     T_periods = 8L,
                                     beta = 1,
                                     pi = 1,
                                     theta_true = 0.05,
                                     truth_pattern = c("uniform", "zone", "gradient"),
                                     seed = 1L) {
  truth_pattern <- match.arg(truth_pattern)
  panel <- .spliv_sim_common_panel(nx, ny, T_periods, seed)
  d_i <- spliv_sim_truth_pattern(panel, truth_pattern)

  panel$x <- pi * panel$z + panel$unit_fe_x + panel$time_fe_x + panel$v
  panel$y <- beta * panel$x + theta_true * d_i * panel$z_scaled_resid +
    panel$unit_fe_y + panel$time_fe_y + panel$e

  data_est <- spliv_sim_residualize_panel(panel)
  list(
    raw_data = panel,
    data_est = data_est,
    truth = list(
      beta_true = beta,
      pi = pi,
      theta_true = theta_true,
      truth_pattern = truth_pattern,
      z_resid_sd = stats::sd(data_est$z)
    )
  )
}

simulate_bpe_panel <- function(nx = 8L,
                               ny = 8L,
                               T_periods = 8L,
                               beta = 1,
                               pi_active = 1,
                               theta_inactive = 0.05,
                               inactive_share = 0.25,
                               transport_gap = 0,
                               seed = 1L) {
  panel <- .spliv_sim_common_panel(nx, ny, T_periods, seed)
  panel$inactive_condition <- spliv_sim_make_inactive_condition(panel, inactive_share)

  first_stage_strength <- ifelse(panel$inactive_condition, 0, pi_active)
  theta_active <- theta_inactive + transport_gap
  theta_unit <- ifelse(panel$inactive_condition, theta_inactive, theta_active)

  panel$x <- first_stage_strength * panel$z + panel$unit_fe_x + panel$time_fe_x + panel$v
  panel$y <- beta * panel$x + theta_unit * panel$z_scaled_resid +
    panel$unit_fe_y + panel$time_fe_y + panel$e

  data_est <- spliv_sim_residualize_panel(panel)
  list(
    raw_data = panel,
    data_est = data_est,
    truth = list(
      beta_true = beta,
      pi_active = pi_active,
      theta_inactive = theta_inactive,
      theta_active = theta_active,
      inactive_share = inactive_share,
      transport_gap = transport_gap,
      z_resid_sd = stats::sd(data_est$z)
    )
  )
}

simulate_search_panel <- function(nx = 8L,
                                  ny = 8L,
                                  T_periods = 8L,
                                  beta = 1,
                                  pi = 1,
                                  first_stage_noise_multiplier = 1,
                                  seed = 1L) {
  panel <- .spliv_sim_common_panel(nx, ny, T_periods, seed)
  panel$inactive_condition <- FALSE

  panel$x <- pi * panel$z + panel$unit_fe_x + panel$time_fe_x +
    first_stage_noise_multiplier * panel$v
  panel$y <- beta * panel$x + panel$unit_fe_y + panel$time_fe_y + panel$e

  data_est <- spliv_sim_residualize_panel(panel)
  list(
    raw_data = panel,
    data_est = data_est,
    truth = list(
      beta_true = beta,
      pi = pi,
      first_stage_noise_multiplier = first_stage_noise_multiplier,
      has_true_inactive_subset = FALSE
    )
  )
}

spliv_sim_candidate_groups <- function(data, K, subgroup_share = 0.10, seed = 1L) {
  set.seed(seed)
  unit_levels <- unique(as.character(data$unit_id))
  out <- vector("list", K)

  for (k in seq_len(K)) {
    keep_trying <- TRUE
    tries <- 0L
    while (keep_trying) {
      tries <- tries + 1L
      unit_flags <- stats::rbinom(length(unit_levels), 1L, subgroup_share) == 1L
      keep_trying <- (sum(unit_flags) < 2L || sum(unit_flags) > (length(unit_levels) - 2L)) && tries < 100L
    }
    if (sum(unit_flags) < 2L) {
      unit_flags[seq_len(min(2L, length(unit_flags)))] <- TRUE
    }
    if (sum(unit_flags) >= length(unit_flags)) {
      unit_flags[length(unit_flags)] <- FALSE
    }
    column_name <- sprintf("random_group_%03d", k)
    out[[k]] <- data.frame(
      setNames(
        list(as.logical(unit_flags[match(as.character(data$unit_id), unit_levels)])),
        column_name
      ),
      check.names = FALSE
    )
  }

  extras <- do.call(cbind, out)
  list(
    data = cbind(data, extras),
    candidate_names = names(extras)
  )
}

spliv_sim_candidate_group_matrix <- function(unit_id,
                                             K,
                                             subgroup_share = 0.10,
                                             seed = 1L) {
  set.seed(seed)
  unit_levels <- unique(as.character(unit_id))
  n_units <- length(unit_levels)
  if (n_units < 4L) {
    stop("Need at least 4 units to build subgroup-search candidate groups.")
  }

  group_size <- as.integer(round(subgroup_share * n_units))
  group_size <- max(2L, min(n_units - 2L, group_size))
  membership <- matrix(0L, nrow = n_units, ncol = K)

  for (k in seq_len(K)) {
    membership[sample.int(n_units, size = group_size, replace = FALSE), k] <- 1L
  }

  colnames(membership) <- sprintf("random_group_%03d", seq_len(K))
  rownames(membership) <- unit_levels

  list(
    membership = membership,
    unit_levels = unit_levels,
    group_size_units = group_size,
    subgroup_share_realized = group_size / n_units
  )
}

make_patterned_scenarios <- function(config) {
  scenarios <- expand.grid(
    truth_pattern = c("uniform", "zone", "gradient"),
    pi = c(0.3, 1.0),
    theta_true = c(0.00, 0.05, 0.10),
    stringsAsFactors = FALSE
  )

  if (identical(config$profile, "pilot")) {
    scenarios <- rbind(
      scenarios[scenarios$truth_pattern == "uniform" & scenarios$pi == 1.0 & scenarios$theta_true == 0.00, ],
      scenarios[scenarios$truth_pattern == "zone" & scenarios$pi == 1.0 & scenarios$theta_true == 0.05, ],
      scenarios[scenarios$truth_pattern == "gradient" & scenarios$pi == 0.3 & scenarios$theta_true == 0.10, ]
    )
  } else if (identical(config$profile, "moderate")) {
    keep <- scenarios$theta_true != 0.00 | scenarios$truth_pattern == "uniform"
    scenarios <- scenarios[keep, , drop = FALSE]
  }

  scenarios$scenario_id <- seq_len(nrow(scenarios))
  scenarios$scenario_name <- sprintf(
    "patterned_%02d_%s_pi_%s_theta_%s",
    scenarios$scenario_id,
    scenarios$truth_pattern,
    gsub("\\.", "p", format(scenarios$pi, nsmall = 1)),
    gsub("\\.", "p", format(scenarios$theta_true, nsmall = 2))
  )
  rownames(scenarios) <- NULL
  scenarios
}

make_bpe_scenarios <- function(config) {
  scenarios <- expand.grid(
    pi_active = c(0.3, 1.0),
    inactive_share = c(0.10, 0.25),
    theta_inactive = c(0.00, 0.05, 0.10),
    transport_gap = c(0.00, 0.05),
    stringsAsFactors = FALSE
  )

  if (identical(config$profile, "pilot")) {
    scenarios <- rbind(
      scenarios[scenarios$pi_active == 1.0 & scenarios$inactive_share == 0.25 &
                  scenarios$theta_inactive == 0.00 & scenarios$transport_gap == 0.00, ],
      scenarios[scenarios$pi_active == 0.3 & scenarios$inactive_share == 0.25 &
                  scenarios$theta_inactive == 0.05 & scenarios$transport_gap == 0.00, ],
      scenarios[scenarios$pi_active == 1.0 & scenarios$inactive_share == 0.10 &
                  scenarios$theta_inactive == 0.05 & scenarios$transport_gap == 0.05, ]
    )
  } else if (identical(config$profile, "moderate")) {
    keep <- !(scenarios$theta_inactive == 0.10 & scenarios$transport_gap == 0.05)
    scenarios <- scenarios[keep, , drop = FALSE]
  }

  scenarios$scenario_id <- seq_len(nrow(scenarios))
  scenarios$scenario_name <- sprintf(
    "bpe_%02d_pi_%s_share_%s_theta_%s_gap_%s",
    scenarios$scenario_id,
    gsub("\\.", "p", format(scenarios$pi_active, nsmall = 1)),
    gsub("\\.", "p", format(scenarios$inactive_share, nsmall = 2)),
    gsub("\\.", "p", format(scenarios$theta_inactive, nsmall = 2)),
    gsub("\\.", "p", format(scenarios$transport_gap, nsmall = 2))
  )
  rownames(scenarios) <- NULL
  scenarios
}

make_search_scenarios <- function(config) {
  scenarios <- expand.grid(
    K = c(5L, 20L, 100L),
    subgroup_share = c(0.10, 0.25),
    pi = c(0.3, 1.0),
    stringsAsFactors = FALSE
  )

  if (identical(config$profile, "pilot")) {
    scenarios <- rbind(
      scenarios[scenarios$K == 5L & scenarios$subgroup_share == 0.10 & scenarios$pi == 1.0, ],
      scenarios[scenarios$K == 20L & scenarios$subgroup_share == 0.10 & scenarios$pi == 0.3, ],
      scenarios[scenarios$K == 100L & scenarios$subgroup_share == 0.25 & scenarios$pi == 1.0, ]
    )
  } else if (identical(config$profile, "moderate")) {
    scenarios <- scenarios[!(scenarios$K == 100L & scenarios$subgroup_share == 0.25), , drop = FALSE]
  }

  scenarios$scenario_id <- seq_len(nrow(scenarios))
  scenarios$scenario_name <- sprintf(
    "search_%02d_K_%03d_share_%s_pi_%s",
    scenarios$scenario_id,
    scenarios$K,
    gsub("\\.", "p", format(scenarios$subgroup_share, nsmall = 2)),
    gsub("\\.", "p", format(scenarios$pi, nsmall = 1))
  )
  rownames(scenarios) <- NULL
  scenarios
}

make_search_stress_scenarios <- function(config) {
  scenarios <- expand.grid(
    pi = c(0.05, 0.10, 0.20),
    K = c(20L, 100L, 500L),
    subgroup_share = c(0.01, 0.02, 0.05, 0.10),
    noise_multiplier = c(1, 2),
    stringsAsFactors = FALSE
  )

  scenarios$scenario_id <- seq_len(nrow(scenarios))
  scenarios$scenario_name <- sprintf(
    "search_stress_%02d_pi_%s_K_%03d_share_%s_noise_%s",
    scenarios$scenario_id,
    gsub("\\.", "p", format(scenarios$pi, nsmall = 2)),
    scenarios$K,
    gsub("\\.", "p", format(scenarios$subgroup_share, nsmall = 2)),
    as.character(scenarios$noise_multiplier)
  )
  rownames(scenarios) <- NULL
  scenarios
}
