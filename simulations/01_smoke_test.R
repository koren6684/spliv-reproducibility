#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) > 0L) {
  raw_path <- sub("^--file=", "", file_arg[1])
  raw_path <- gsub("~\\+~", " ", raw_path)
  script_path <- tryCatch(normalizePath(raw_path, mustWork = TRUE), error = function(e) raw_path)
  script_dir <- dirname(script_path)
  if (dir.exists(script_dir)) {
    setwd(script_dir)
  }
}

source("00_paths.R")
source("00_config.R")

paths <- spliv_sim_paths(getwd())
spliv_sim_source_helpers(paths)
package_info <- spliv_sim_load_package(paths)

config <- spliv_sim_resolve_config(paths, profile = "pilot")
config$R <- 2L
config$nx <- 4L
config$ny <- 4L
config$T_periods <- 4L
config$delta_grid <- c(0, 0.05)
config$n_cores <- 1L
dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(config$raw_dir, recursive = TRUE, showWarnings = FALSE)

spliv_sim_print_config(config)
log_file <- spliv_sim_new_log_file(paths, "01_smoke_test", config$profile)
spliv_sim_log_header(log_file, "01_smoke_test", config, package_info)

patterned <- simulate_patterned_panel(
  nx = config$nx,
  ny = config$ny,
  T_periods = config$T_periods,
  beta = 1,
  pi = 1,
  theta_true = 0.05,
  truth_pattern = "zone",
  seed = config$base_seed
)
patterns <- spliv_sim_pattern_objects(patterned$data_est)
baseline <- run_baseline_iv(patterned$data_est)
path <- run_uci_path(
  patterned$data_est,
  delta_grid = config$delta_grid,
  violation_pattern = patterns$zone
)

bpe_sim <- simulate_bpe_panel(
  nx = config$nx,
  ny = config$ny,
  T_periods = config$T_periods,
  beta = 1,
  pi_active = 1,
  theta_inactive = 0.05,
  inactive_share = 0.25,
  transport_gap = 0,
  seed = config$base_seed + 1L
)
inactive_design <- spliv_sim_bpe_design()
bpe_validation <- run_bpe_validation(
  bpe_sim$data_est,
  design = inactive_design,
  vcov = "cluster",
  cluster = ~ unit_id,
  bpe_equiv_margin = 0.20,
  bpe_min_n_S = 8,
  bpe_min_clusters_S = 3
)
bpe_fit <- run_bpe_confirmatory(
  bpe_sim$data_est,
  design = inactive_design,
  vcov = "cluster",
  cluster = ~ unit_id,
  bpe_equiv_margin = 0.20,
  bpe_min_n_S = 8,
  bpe_min_clusters_S = 3
)

bpe_fit_reason <- bpe_fit$error %||% if (!isTRUE(bpe_fit$summary$fit_ok)) {
  if (isTRUE(bpe_validation$summary$eligibility_passed)) {
    "BPE returned NA intervals in the smoke configuration."
  } else {
    "BPE eligibility failed in the smoke configuration."
  }
} else {
  NA_character_
}

summary_row <- data.frame(
  package_loaded = TRUE,
  package_path = package_info$path,
  package_version = package_info$version,
  baseline_ok = is.finite(baseline$summary$estimate),
  sensitivity_path_ok = nrow(path$path) > 0L,
  pattern_name = unique(path$path$pattern_name)[1],
  bpe_validation_ok = isTRUE(bpe_validation$summary$ok),
  bpe_validation_eligibility = bpe_validation$summary$eligibility_passed,
  bpe_fit_ok = isTRUE(bpe_fit$summary$fit_ok),
  bpe_fit_error = bpe_fit_reason,
  output_dir = config$output_dir,
  stringsAsFactors = FALSE
)

smoke_rds <- file.path(config$raw_dir, "smoke_test_summary.rds")
smoke_csv <- file.path(paths$tables_dir, "smoke_test_summary.csv")
saveRDS(
  list(
    summary = summary_row,
    baseline = baseline$summary,
    path = path$path,
    bpe_validation = bpe_validation$summary,
    bpe_fit = bpe_fit$summary
  ),
  smoke_rds
)
utils::write.csv(summary_row, smoke_csv, row.names = FALSE)

stopifnot(file.exists(smoke_rds))
stopifnot(file.exists(smoke_csv))
stopifnot(isTRUE(summary_row$baseline_ok))
stopifnot(isTRUE(summary_row$sensitivity_path_ok))
stopifnot(isTRUE(summary_row$bpe_validation_ok))

spliv_sim_log(log_file, "Smoke-test baseline estimate: ", round(baseline$summary$estimate, 4))
spliv_sim_log(log_file, "Smoke-test path tipping point: ", path$tipping_point)
if (isTRUE(summary_row$bpe_fit_ok)) {
  spliv_sim_log(log_file, "Smoke-test confirmatory BPE fit succeeded.")
} else {
  spliv_sim_log(
    log_file,
    "Smoke-test confirmatory BPE fit did not complete; expected diagnostic: ",
    bpe_fit_reason
  )
}
spliv_sim_log(log_file, "Smoke-test outputs written to: ", smoke_rds, " and ", smoke_csv)

print(summary_row)
