#!/usr/bin/env Rscript

cat("This is the full subgroup-search stress-test run and may take many hours.\n")

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
config <- spliv_sim_resolve_config(paths, profile = Sys.getenv("SPLIV_SIM_PROFILE", "full"))
if (!identical(config$profile, "full")) {
  stop("`92_run_full_subgroup_search_stress.R` requires `SPLIV_SIM_PROFILE=full`.")
}
spliv_sim_print_config(config)

log_file <- spliv_sim_new_log_file(paths, "92_run_full_subgroup_search_stress", config$profile)
spliv_sim_log_header(
  log_file,
  "92_run_full_subgroup_search_stress",
  config,
  package_info,
  make_search_stress_scenarios(config)
)

family_run <- run_search_stress_family(config, paths, package_info, log_file = log_file)
summary_out <- summarize_search_stress_family(family_run$results, paths, config)

spliv_sim_log(log_file, "Subgroup-search stress summary written to: ", summary_out$file)
spliv_sim_log(log_file, "Subgroup-search stress figure written to: ", summary_out$figure)
print(family_run$counts)
