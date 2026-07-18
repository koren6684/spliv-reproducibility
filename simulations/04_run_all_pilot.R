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

log_file <- spliv_sim_new_log_file(paths, "04_run_all_pilot", config$profile)
spliv_sim_log_header(log_file, "04_run_all_pilot", config, package_info)

patterned_run <- run_patterned_family(config, paths, package_info, log_file = log_file)
patterned_summary <- summarize_patterned_family(patterned_run$results, paths, config)
spliv_sim_log(log_file, "Patterned summary file: ", patterned_summary$file)

bpe_run <- run_bpe_family(config, paths, package_info, log_file = log_file)
bpe_summary <- summarize_bpe_family(bpe_run$results, paths, config)
spliv_sim_log(log_file, "BPE summary file: ", bpe_summary$file)

search_run <- run_search_family(config, paths, package_info, log_file = log_file)
search_summary <- summarize_search_family(search_run$results, paths, config)
spliv_sim_log(log_file, "Subgroup-search summary file: ", search_summary$file)

spliv_sim_log(log_file, "Pilot tables written in: ", paths$tables_dir)
spliv_sim_log(log_file, "Pilot figures written in: ", paths$figures_dir)

print(list(
  patterned = patterned_run$counts,
  bpe = bpe_run$counts,
  subgroup_search = search_run$counts
))
