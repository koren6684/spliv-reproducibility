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

source(file.path(script_dir, "R", "helpers_koren_paths.R"))
source(file.path(script_dir, "R", "helpers_koren_data.R"))
source(file.path(script_dir, "R", "helpers_koren_models.R"))
source(file.path(script_dir, "R", "helpers_koren_outputs.R"))

run_start <- Sys.time()
out_dirs <- koren_prepare_dirs(script_dir)
pkg_info <- koren_load_spliv(script_dir)
source(file.path(script_dir, "koren_design.R"))

fast_mode <- koren_bool_env("SPLIV_REPL_FAST_MODE", default = TRUE)
delta_share_grid <- if (isTRUE(fast_mode)) {
  seq(0, 0.20, by = 0.05)
} else {
  seq(0, 0.20, by = 0.01)
}
uci_steps_default <- if (isTRUE(fast_mode)) 5L else 21L
uci_steps <- as.integer(Sys.getenv("KOREN_UCI_STEPS", as.character(uci_steps_default)))
if (!is.finite(uci_steps) || uci_steps < 2) {
  stop("`KOREN_UCI_STEPS` must be an integer >= 2.", call. = FALSE)
}

message("FAST mode: ", fast_mode)
message("Delta share grid: ", paste(delta_share_grid, collapse = ", "))
message("UCI theta grid steps: ", uci_steps)

data_file <- koren_find_panel_data(script_dir)
message("Reading Koren panel data: ", data_file)
dat_raw <- koren_read_panel_data(data_file)
koren_validate_variables(dat_raw)

patterns <- make_koren_patterns()
specs <- koren_table3_specs()

baseline_rows <- list()
path_rows <- list()
ltz_rows <- list()
bpe_rows <- list()
run_notes <- character()

for (i in seq_len(nrow(specs))) {
  spec <- specs[i, , drop = FALSE]
  crop <- spec$crop
  treatment <- spec$treatment

  message("Running Koren Table 3 updated workflow for ", crop, " (", treatment, ")")
  dat <- koren_analysis_data(dat_raw, treatment)
  bpe_design_info <- make_koren_bpe_design_main(dat)

  fit0 <- koren_baseline_fit(dat, treatment)
  beta0 <- koren_extract_term(fit0, treatment)
  metrics <- koren_residual_metrics(fit0, treatment)
  delta_grid_y <- delta_share_grid * metrics$residual_y_sd

  baseline_rows[[crop]] <- data.frame(
    crop = crop,
    treatment = treatment,
    outcome = "acled_inc_sum",
    instrument = "spi6",
    fixed_effects = "gid + year",
    cluster = "gid",
    nobs = metrics$nobs,
    n_clusters = metrics$n_clusters,
    estimate = koren_scalar(beta0$estimate),
    std_error = koren_scalar(beta0$std.error),
    conf_low = koren_scalar(beta0$conf.low),
    conf_high = koren_scalar(beta0$conf.high),
    residual_y_sd = metrics$residual_y_sd,
    residual_x_sd = metrics$residual_x_sd,
    residual_z_sd = metrics$residual_z_sd,
    first_stage_coefficient_full = metrics$first_stage_full[["coefficient"]],
    first_stage_se_full = metrics$first_stage_full[["se"]],
    first_stage_f_full_diagnostic_only = metrics$first_stage_full[["f_statistic"]],
    stringsAsFactors = FALSE
  )

  for (pattern_name in names(patterns)) {
    pat <- patterns[[pattern_name]]
    path <- koren_run_sensitivity_path(
      dat = dat,
      treatment = treatment,
      pattern = pat,
      delta_grid_y = delta_grid_y,
      method = "uci",
      uci_steps = uci_steps
    )
    path$crop <- crop
    path$treatment <- treatment
    path$delta_outcome_units <- path$delta
    path$delta_share_of_residual_y_sd <- path$delta / metrics$residual_y_sd
    path$residual_y_sd <- metrics$residual_y_sd
    path_rows[[paste(crop, pattern_name, sep = "_")]] <- path

    ltz_path <- koren_run_sensitivity_path(
      dat = dat,
      treatment = treatment,
      pattern = pat,
      delta_grid_y = delta_grid_y,
      method = "ltz",
      uci_steps = uci_steps
    )
    ltz_path$crop <- crop
    ltz_path$treatment <- treatment
    ltz_path$delta_outcome_units <- ltz_path$delta
    ltz_path$delta_share_of_residual_y_sd <- ltz_path$delta / metrics$residual_y_sd
    ltz_path$residual_y_sd <- metrics$residual_y_sd
    ltz_rows[[paste(crop, pattern_name, sep = "_")]] <- ltz_path
  }

  bpe <- koren_run_bpe(
    dat = dat,
    crop = crop,
    treatment = treatment,
    design_info = bpe_design_info,
    residual_x_sd = metrics$residual_x_sd,
    min_n_S = 200,
    min_clusters_S = 30
  )
  bpe_diag <- bpe$diagnostics
  bpe_diag$residualized_treatment_sd_full <- metrics$residual_x_sd
  bpe_rows[[crop]] <- bpe_diag

  run_notes <- c(
    run_notes,
    paste0(
      "- ", crop, ": n=", metrics$nobs,
      ", residual_y_sd=", signif(metrics$residual_y_sd, 5),
      ", residual_x_sd=", signif(metrics$residual_x_sd, 5),
      ", primary BPE rule=", bpe_design_info$subset_rule,
      ", BPE main eligible=", bpe_diag$bpe_eligibility_passed[bpe_diag$margin_share == 0.05]
    )
  )
}

baseline_dt <- do.call(rbind, baseline_rows)
paths_dt <- do.call(rbind, path_rows)
row.names(paths_dt) <- NULL
ltz_dt <- do.call(rbind, ltz_rows)
row.names(ltz_dt) <- NULL
bpe_dt <- do.call(rbind, bpe_rows)
row.names(bpe_dt) <- NULL

sens_summary <- koren_sensitivity_summary(paths_dt, baseline_dt)
ltz_summary <- koren_sensitivity_summary(ltz_dt, baseline_dt)
app_summary <- koren_application_summary(baseline_dt, sens_summary, bpe_dt, ltz_summary = ltz_summary)

koren_write_csv(baseline_dt, file.path(out_dirs$tables, "koren_baseline_iv.csv"))
koren_write_csv(paths_dt, file.path(out_dirs$tables, "koren_sensitivity_paths.csv"))
koren_write_csv(sens_summary, file.path(out_dirs$tables, "koren_sensitivity_summary.csv"))
koren_write_csv(ltz_dt, file.path(out_dirs$tables, "koren_ltz_paths.csv"))
koren_write_csv(ltz_summary, file.path(out_dirs$tables, "koren_ltz_summary.csv"))
koren_write_csv(bpe_dt, file.path(out_dirs$tables, "koren_bpe_diagnostics.csv"))
koren_write_csv(app_summary, file.path(out_dirs$tables, "koren_application_summary.csv"))

koren_plot_sensitivity(paths_dt, file.path(out_dirs$figures, "koren_sensitivity_paths.pdf"), "pdf")
koren_plot_sensitivity(paths_dt, file.path(out_dirs$figures, "koren_sensitivity_paths.png"), "png")
koren_plot_sensitivity(
  ltz_dt,
  file.path(out_dirs$figures, "koren_ltz_paths.pdf"),
  "pdf",
  method_label = "LTZ",
  overall_title = "LTZ sensitivity paths"
)
koren_plot_sensitivity(
  ltz_dt,
  file.path(out_dirs$figures, "koren_ltz_paths.png"),
  "png",
  method_label = "LTZ",
  overall_title = "LTZ sensitivity paths"
)
koren_plot_bpe(bpe_dt, file.path(out_dirs$figures, "koren_bpe_diagnostics.pdf"), "pdf")
koren_plot_bpe(bpe_dt, file.path(out_dirs$figures, "koren_bpe_diagnostics.png"), "png")

rplots_path <- file.path(script_dir, "Rplots.pdf")
if (file.exists(rplots_path)) {
  unlink(rplots_path)
}

run_end <- Sys.time()
main_bpe <- bpe_dt[bpe_dt$margin_share == 0.05, , drop = FALSE]
log_lines <- c(
  "# Koren 2018 Table 3 Updated SPLIV Replication Log",
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
  paste0("- File: ", data_file),
  paste0("- Raw rows: ", nrow(dat_raw)),
  paste0("- Raw columns: ", ncol(dat_raw)),
  paste0("- Required variables: ", paste(koren_required_variables(), collapse = ", ")),
  "",
  "## Specification",
  "- Outcome: acled_inc_sum",
  "- Endogenous treatments: maize_yield, wheat_yield",
  "- Instrument: spi6",
  "- Fixed effects: gid + year",
  "- Cluster: gid",
  "- Controls: none in this Table 3 first pass",
  "",
  "## Delta Grid",
  paste0("- Delta shares of residual outcome SD: ", paste(delta_share_grid, collapse = ", ")),
  paste0("- UCI theta grid steps per delta: ", uci_steps),
  "- Delta outcome-unit grids are crop-specific: delta_share * residual_y_sd.",
  "- LTZ is run on the same delta shares and outcome-unit grids for comparability with the original plausibly exogenous analysis.",
  "- UCI remains the conservative main sensitivity benchmark.",
  "- LTZ uses the same outcome-unit scaling and same direct-effect patterns as UCI.",
  "",
  "## Patterns",
  paste0("- ", patterns$uniform$name, ": ", patterns$uniform$rationale),
  paste0("- ", patterns$sparsebare$name, ": ", patterns$sparsebare$rationale),
  "- `sparsebare` remains the direct-effect pattern for sensitivity and is present in the original Koren panel file.",
  "",
  "## Confirmatory BPE Design",
  "- Main rule: sparsebare >= 90th percentile in the crop-specific analysis sample.",
  "- Rationale: the crop-yield treatment channel should be weakest where sparse/bare land-cover share is highest and ordinary crop production is least central.",
  "- `sparsebare` exists in the original Koren panel file and is used as both the direct-effect pattern and the main confirmatory BPE inactive subset.",
  "- `sparsebare` is a SPLIV extension design variable, not part of the original Koren estimating equation.",
  "- No BPE threshold is chosen by first-stage strength, equivalence results, reduced-form estimates, or outcomes.",
  "- No BPE subgroup search was used.",
  "- Minimum subset size: 200",
  "- Minimum clusters in S: 30",
  "- Main equivalence margin: 0.05 * residualized treatment SD.",
  "- Appendix diagnostic margin: 0.10 * residualized treatment SD; not used to rescue main eligibility.",
  "",
  "## Crop Notes",
  run_notes,
  "",
  "## BPE Eligibility",
  paste0(
    "- ", main_bpe$crop,
    ": eligible=", main_bpe$bpe_eligibility_passed,
    ", equivalence_passed=", main_bpe$equivalence_passed,
    ", estimation_source=", main_bpe$bpe_estimation_source,
    ", reason=", ifelse(nchar(main_bpe$reason_if_not_applicable), main_bpe$reason_if_not_applicable, "eligible")
  ),
  "",
  "## BPE Estimation Notes",
  paste0(
    "- ", main_bpe$crop,
    ": ", ifelse(nchar(main_bpe$bpe_fit_error), main_bpe$bpe_fit_error, "No BPE estimation error.")
  ),
  "",
  "## Outputs",
  paste0("- ", file.path(out_dirs$tables, "koren_baseline_iv.csv")),
  paste0("- ", file.path(out_dirs$tables, "koren_sensitivity_paths.csv")),
  paste0("- ", file.path(out_dirs$tables, "koren_sensitivity_summary.csv")),
  paste0("- ", file.path(out_dirs$tables, "koren_ltz_paths.csv")),
  paste0("- ", file.path(out_dirs$tables, "koren_ltz_summary.csv")),
  paste0("- ", file.path(out_dirs$tables, "koren_bpe_diagnostics.csv")),
  paste0("- ", file.path(out_dirs$tables, "koren_application_summary.csv")),
  paste0("- ", file.path(out_dirs$figures, "koren_sensitivity_paths.pdf")),
  paste0("- ", file.path(out_dirs$figures, "koren_sensitivity_paths.png")),
  paste0("- ", file.path(out_dirs$figures, "koren_ltz_paths.pdf")),
  paste0("- ", file.path(out_dirs$figures, "koren_ltz_paths.png")),
  paste0("- ", file.path(out_dirs$figures, "koren_bpe_diagnostics.pdf")),
  paste0("- ", file.path(out_dirs$figures, "koren_bpe_diagnostics.png"))
)
koren_write_log(file.path(out_dirs$logs, "koren_run_log.md"), log_lines)

cat("\nKoren 2018 Table 3 updated replication complete.\n")
cat("Data file:", data_file, "\n")
cat("FAST mode:", fast_mode, "\n")
cat("Delta share grid:", paste(delta_share_grid, collapse = ", "), "\n")
cat("\nBaseline estimates:\n")
print(baseline_dt[, c("crop", "treatment", "estimate", "std_error", "conf_low", "conf_high", "residual_y_sd")])
cat("\nBPE main eligibility:\n")
print(main_bpe[, c("crop", "design_name", "subset_definition", "n_S", "G_S", "equivalence_margin", "equivalence_passed", "bpe_eligibility_passed", "reason_if_not_applicable")])
cat("\nBPE estimation source/errors:\n")
print(main_bpe[, c("crop", "bpe_estimation_source", "bpe_beta_estimate", "bpe_beta_conf_low", "bpe_beta_conf_high", "bpe_fit_error")])
cat("\nOutputs written under:", out_dirs$root, "\n")

invisible(list(
  baseline = baseline_dt,
  sensitivity_paths = paths_dt,
  sensitivity_summary = sens_summary,
  ltz_paths = ltz_dt,
  ltz_summary = ltz_summary,
  bpe_diagnostics = bpe_dt,
  application_summary = app_summary,
  package = pkg_info,
  data_file = data_file,
  fast_mode = fast_mode
))
