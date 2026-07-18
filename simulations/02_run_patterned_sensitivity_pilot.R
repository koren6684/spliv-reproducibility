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
config <- spliv_sim_resolve_config(paths, profile = Sys.getenv("SPLIV_SIM_PROFILE", "pilot"))
spliv_sim_print_config(config)

log_file <- spliv_sim_new_log_file(paths, "02_run_patterned_sensitivity_pilot", config$profile)
spliv_sim_log_header(log_file, "02_run_patterned_sensitivity_pilot", config, package_info, make_patterned_scenarios(config))

family_run <- run_patterned_family(config, paths, package_info, log_file = log_file)
summary_out <- summarize_patterned_family(family_run$results, paths, config)

spliv_sim_log(log_file, "Patterned summary written to: ", summary_out$file)
spliv_sim_log(log_file, "Patterned figure files created in: ", paths$figures_dir)
print(family_run$counts)
