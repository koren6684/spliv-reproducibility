#!/usr/bin/env Rscript

script_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    raw_path <- gsub("~\\+~", " ", sub("^--file=", "", file_arg[[1]]))
    dirname(normalizePath(raw_path, mustWork = TRUE))
  } else {
    normalizePath(getwd(), mustWork = TRUE)
  }
})

source(file.path(script_dir, "R", "helpers_lelkes_paths.R"))
source(file.path(script_dir, "R", "helpers_lelkes_data.R"))
source(file.path(script_dir, "R", "helpers_lelkes_models.R"))
source(file.path(script_dir, "R", "helpers_lelkes_outputs.R"))

run_start <- Sys.time()
out_dirs <- lelkes_prepare_dirs(script_dir)
pkg_info <- lelkes_load_spliv(script_dir)
source(file.path(script_dir, "lelkes_design.R"))

fast_mode <- lelkes_bool_env("SPLIV_REPL_FAST_MODE", default = TRUE)
delta_share_grid <- if (isTRUE(fast_mode)) {
  seq(0, 0.20, by = 0.05)
} else {
  seq(0, 0.20, by = 0.01)
}
uci_steps_default <- if (isTRUE(fast_mode)) 5L else 21L
uci_steps <- as.integer(Sys.getenv("LELKES_UCI_STEPS", as.character(uci_steps_default)))
if (!is.finite(uci_steps) || uci_steps < 2) {
  stop("`LELKES_UCI_STEPS` must be an integer >= 2.", call. = FALSE)
}

message("FAST mode: ", fast_mode)
message("Delta share grid: ", paste(delta_share_grid, collapse = ", "))
message("UCI theta grid steps: ", uci_steps)

data_files <- lelkes_find_data_files(script_dir)
message("Reading Lelkes merged data: ", data_files$merged)
message("Reading Lelkes county data: ", data_files$county)
raw <- lelkes_load_raw_data(data_files)
dat <- lelkes_prepare_analysis_data(raw)

patterns <- make_lelkes_patterns()
bpe_design_info <- make_lelkes_bpe_design(dat)

message("Running Lelkes Table 1 baseline IV.")
fit0 <- lelkes_baseline_fit(dat)
beta0 <- lelkes_extract_term(fit0, "log_providers")
metrics <- lelkes_residual_metrics(dat, fit0)
delta_grid_y <- delta_share_grid * metrics$residual_y_sd

message("Running confirmatory BPE validation for the pre-specified low-density subset.")
bpe <- lelkes_run_bpe(
  dat = dat,
  design_info = bpe_design_info,
  residual_x_sd = metrics$residual_x_sd,
  min_n_S = 200,
  min_clusters_S = 30
)
bpe_dt <- bpe$diagnostics
bpe_status <- lelkes_bpe_status_from_diagnostics(bpe_dt)

baseline_dt <- data.frame(
  application = "Lelkes, Sood, and Iyengar Table 1",
  outcome = "affective_polarization",
  treatment = "log_providers",
  instrument = "log_Total",
  controls = paste(
    "as.factor(year)", "region", "percent_black", "percent_white",
    "percent_male", "lowed", "unemploymentrate", "density", "log_HHINC",
    sep = " + "
  ),
  cluster = "state",
  nobs = metrics$nobs,
  n_clusters = metrics$n_clusters,
  estimate = as.numeric(beta0$estimate),
  std_error = as.numeric(beta0$std.error),
  conf_low = as.numeric(beta0$conf.low),
  conf_high = as.numeric(beta0$conf.high),
  residual_y_sd = metrics$residual_y_sd,
  residual_x_sd = metrics$residual_x_sd,
  residual_z_sd = metrics$residual_z_sd,
  first_stage_coefficient_full = metrics$first_stage_full[["coefficient"]],
  first_stage_se_full = metrics$first_stage_full[["se"]],
  first_stage_f_full_diagnostic_only = metrics$first_stage_full[["f_statistic"]],
  stringsAsFactors = FALSE
)

uci_rows <- list()
ltz_rows <- list()
for (pattern_name in names(patterns)) {
  pat <- patterns[[pattern_name]]
  message("Running UCI path for pattern: ", pat$name)
  uci <- lelkes_run_sensitivity_path(
    dat = dat,
    pattern = pat,
    method = "uci",
    delta_grid_y = delta_grid_y,
    uci_steps = uci_steps
  )
  uci$delta_outcome_units <- uci$delta
  uci$delta_share_of_residual_y_sd <- uci$delta / metrics$residual_y_sd
  uci$residual_y_sd <- metrics$residual_y_sd
  uci_rows[[pattern_name]] <- uci

  message("Running LTZ path for pattern: ", pat$name)
  ltz <- lelkes_run_sensitivity_path(
    dat = dat,
    pattern = pat,
    method = "ltz",
    delta_grid_y = delta_grid_y,
    uci_steps = uci_steps
  )
  ltz$delta_outcome_units <- ltz$delta
  ltz$delta_share_of_residual_y_sd <- ltz$delta / metrics$residual_y_sd
  ltz$residual_y_sd <- metrics$residual_y_sd
  ltz_rows[[pattern_name]] <- ltz
}

uci_dt <- do.call(rbind, uci_rows)
row.names(uci_dt) <- NULL
ltz_dt <- do.call(rbind, ltz_rows)
row.names(ltz_dt) <- NULL
paths_dt <- rbind(uci_dt, ltz_dt)
row.names(paths_dt) <- NULL

sens_summary <- lelkes_sensitivity_summary(paths_dt, baseline_dt)
app_summary <- lelkes_application_summary(baseline_dt, sens_summary, bpe_status, bpe_diag = bpe_dt)

lelkes_write_csv(baseline_dt, file.path(out_dirs$tables, "lelkes_baseline_iv.csv"))
lelkes_write_csv(uci_dt, file.path(out_dirs$tables, "lelkes_uci_paths.csv"))
lelkes_write_csv(ltz_dt, file.path(out_dirs$tables, "lelkes_ltz_paths.csv"))
lelkes_write_csv(sens_summary, file.path(out_dirs$tables, "lelkes_sensitivity_summary.csv"))
lelkes_write_csv(app_summary, file.path(out_dirs$tables, "lelkes_application_summary.csv"))
lelkes_write_csv(bpe_status, file.path(out_dirs$tables, "lelkes_bpe_status.csv"))
lelkes_write_csv(bpe_dt, file.path(out_dirs$tables, "lelkes_bpe_diagnostics.csv"))

lelkes_plot_paths(uci_dt, file.path(out_dirs$figures, "lelkes_uci_paths.pdf"), "pdf", method_label = "UCI")
lelkes_plot_paths(uci_dt, file.path(out_dirs$figures, "lelkes_uci_paths.png"), "png", method_label = "UCI")
lelkes_plot_paths(ltz_dt, file.path(out_dirs$figures, "lelkes_ltz_paths.pdf"), "pdf", method_label = "LTZ")
lelkes_plot_paths(ltz_dt, file.path(out_dirs$figures, "lelkes_ltz_paths.png"), "png", method_label = "LTZ")
lelkes_plot_combined_paths(uci_dt, ltz_dt, file.path(out_dirs$figures, "lelkes_combined_sensitivity_paths.pdf"), "pdf")
lelkes_plot_combined_paths(uci_dt, ltz_dt, file.path(out_dirs$figures, "lelkes_combined_sensitivity_paths.png"), "png")
lelkes_plot_bpe(bpe_dt, file.path(out_dirs$figures, "lelkes_bpe_diagnostics.pdf"), "pdf")
lelkes_plot_bpe(bpe_dt, file.path(out_dirs$figures, "lelkes_bpe_diagnostics.png"), "png")

rplots_path <- file.path(script_dir, "Rplots.pdf")
if (file.exists(rplots_path)) {
  unlink(rplots_path)
}

run_end <- Sys.time()
log_lines <- c(
  "# Lelkes, Sood, and Iyengar Table 1 Updated SPLIV Replication Log",
  "",
  paste0("Start time: ", format(run_start, "%Y-%m-%d %H:%M:%S %Z")),
  paste0("End time: ", format(run_end, "%Y-%m-%d %H:%M:%S %Z")),
  paste0("FAST mode: ", fast_mode),
  "",
  "## Package",
  paste0("- Path: ", pkg_info$path),
  paste0("- Version: ", pkg_info$version),
  paste0("- Loader: ", pkg_info$loader),
  "",
  "## Data",
  paste0("- Individual merged file: ", data_files$merged),
  paste0("- County merged file: ", data_files$county),
  paste0("- Raw individual rows: ", nrow(raw$merged)),
  paste0("- Analysis rows after Table 1 construction and complete cases: ", nrow(dat)),
  paste0("- Required variables: ", paste(lelkes_required_variables(), collapse = ", ")),
  "",
  "## Specification",
  "- Outcome: affective_polarization = zero1(infeels - outfeels)",
  "- Endogenous treatment: log_providers = log(providers)",
  "- Instrument: log_Total = log(Total), the ROW index",
  "- Controls: as.factor(year), region, percent_black, percent_white, percent_male, lowed, unemploymentrate, density, log_HHINC",
  "- Cluster: state",
  "",
  "## Delta Grid",
  paste0("- Delta shares of residual outcome SD: ", paste(delta_share_grid, collapse = ", ")),
  paste0("- UCI theta grid steps per delta: ", uci_steps),
  paste0("- Residual outcome SD: ", signif(metrics$residual_y_sd, 6)),
  "- Delta outcome-unit grid is delta_share * residual_y_sd.",
  "- UCI remains the conservative main sensitivity benchmark.",
  "- LTZ is added for comparability with the original plausibly exogenous analysis.",
  "- UCI and LTZ use the same outcome-unit scaling and direct-effect patterns.",
  "",
  "## Patterns",
  paste0("- ", patterns$uniform$name, ": ", patterns$uniform$rationale),
  paste0("- ", patterns$density$name, ": ", patterns$density$rationale),
  "",
  "## BPE Status",
  paste0("- Status: ", bpe_status$bpe_status),
  paste0("- Reason: ", bpe_status$reason),
  paste0("- Design: ", bpe_design_info$design$name),
  "- Rule: density <= 10th percentile in the Table 1 analysis sample.",
  paste0("- Density threshold: ", signif(bpe_design_info$threshold, 6)),
  "- This threshold is pre-specified in code and is not tuned by first-stage strength, equivalence results, reduced-form estimates, or outcomes.",
  "- Minimum subset size: 200",
  "- Minimum state clusters in S: 30",
  "- Main equivalence margin: 0.05 * residualized treatment SD.",
  "- Appendix diagnostic margin: 0.10 * residualized treatment SD; not used to rescue main eligibility.",
  paste0("- Main BPE n_S: ", bpe_dt$n_S[bpe_dt$margin_share == 0.05]),
  paste0("- Main BPE G_S: ", bpe_dt$G_S[bpe_dt$margin_share == 0.05]),
  paste0("- Main BPE equivalence passed: ", bpe_dt$equivalence_passed[bpe_dt$margin_share == 0.05]),
  paste0("- Main BPE eligibility passed: ", bpe_dt$bpe_eligibility_passed[bpe_dt$margin_share == 0.05]),
  paste0("- BPE estimation source: ", bpe_dt$bpe_estimation_source[bpe_dt$margin_share == 0.05]),
  "",
  "## Baseline",
  paste0("- Estimate: ", signif(baseline_dt$estimate, 6)),
  paste0("- 95% CI: [", signif(baseline_dt$conf_low, 6), ", ", signif(baseline_dt$conf_high, 6), "]"),
  paste0("- First-stage coefficient on log_Total: ", signif(baseline_dt$first_stage_coefficient_full, 6)),
  paste0("- First-stage F diagnostic: ", signif(baseline_dt$first_stage_f_full_diagnostic_only, 6)),
  "",
  "## Outputs",
  paste0("- ", file.path(out_dirs$tables, "lelkes_baseline_iv.csv")),
  paste0("- ", file.path(out_dirs$tables, "lelkes_uci_paths.csv")),
  paste0("- ", file.path(out_dirs$tables, "lelkes_ltz_paths.csv")),
  paste0("- ", file.path(out_dirs$tables, "lelkes_sensitivity_summary.csv")),
  paste0("- ", file.path(out_dirs$tables, "lelkes_application_summary.csv")),
  paste0("- ", file.path(out_dirs$tables, "lelkes_bpe_status.csv")),
  paste0("- ", file.path(out_dirs$tables, "lelkes_bpe_diagnostics.csv")),
  paste0("- ", file.path(out_dirs$figures, "lelkes_uci_paths.pdf")),
  paste0("- ", file.path(out_dirs$figures, "lelkes_uci_paths.png")),
  paste0("- ", file.path(out_dirs$figures, "lelkes_ltz_paths.pdf")),
  paste0("- ", file.path(out_dirs$figures, "lelkes_ltz_paths.png")),
  paste0("- ", file.path(out_dirs$figures, "lelkes_combined_sensitivity_paths.pdf")),
  paste0("- ", file.path(out_dirs$figures, "lelkes_combined_sensitivity_paths.png")),
  paste0("- ", file.path(out_dirs$figures, "lelkes_bpe_diagnostics.pdf")),
  paste0("- ", file.path(out_dirs$figures, "lelkes_bpe_diagnostics.png"))
)
lelkes_write_log(file.path(out_dirs$logs, "lelkes_run_log.md"), log_lines)

cat("\nLelkes Table 1 updated replication complete.\n")
cat("Data files:\n")
cat("  ", data_files$merged, "\n", sep = "")
cat("  ", data_files$county, "\n", sep = "")
cat("FAST mode:", fast_mode, "\n")
cat("Delta share grid:", paste(delta_share_grid, collapse = ", "), "\n")
cat("\nBaseline estimate:\n")
print(baseline_dt[, c("estimate", "std_error", "conf_low", "conf_high", "residual_y_sd", "first_stage_f_full_diagnostic_only")])
cat("\nSensitivity paths:\n")
print(table(paths_dt$method, paths_dt$pattern_name))
cat("\nBPE status:\n")
print(bpe_status)
cat("\nBPE diagnostics, main 5% margin:\n")
print(bpe_dt[
  bpe_dt$margin_share == 0.05,
  c(
    "density_threshold", "n_S", "G_S", "equivalence_margin",
    "first_stage_coefficient", "first_stage_ci_low", "first_stage_ci_high",
    "equivalence_passed", "bpe_eligibility_passed", "bpe_estimation_source",
    "bpe_beta_estimate", "bpe_beta_conf_low", "bpe_beta_conf_high",
    "reason_if_not_applicable"
  )
])
cat("\nOutputs written under:", out_dirs$root, "\n")

invisible(list(
  baseline = baseline_dt,
  uci_paths = uci_dt,
  ltz_paths = ltz_dt,
  sensitivity_summary = sens_summary,
  application_summary = app_summary,
  bpe_status = bpe_status,
  bpe_diagnostics = bpe_dt,
  package = pkg_info,
  data_files = data_files,
  fast_mode = fast_mode
))
