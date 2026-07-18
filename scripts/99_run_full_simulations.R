#!/usr/bin/env Rscript

cat("WARNING: the full simulation can take many hours or days and writes large checkpoint trees.\n")
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)) else getwd()
repo_root <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
sim_root <- file.path(repo_root, "simulations")
pkg_path <- Sys.getenv("SPLIV_PACKAGE_PATH", file.path(repo_root, "..", "spliv"))
Sys.setenv(SPLIV_PACKAGE_PATH = normalizePath(pkg_path, mustWork = TRUE), SPLIV_SIM_PROFILE = "full")
setwd(sim_root)
status <- system2(Sys.which("Rscript"), c("--vanilla", shQuote(file.path(sim_root, "99_run_all_full.R"))))
if (!identical(as.integer(status), 0L)) quit(status = status)
