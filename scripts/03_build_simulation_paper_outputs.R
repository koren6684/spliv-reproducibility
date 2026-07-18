#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)) else getwd()
repo_root <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
sim_root <- file.path(repo_root, "simulations")
setwd(sim_root)
status <- system2(Sys.which("Rscript"), c("--vanilla", shQuote(file.path(sim_root, "05_make_paper_sim_outputs.R"))))
if (!identical(as.integer(status), 0L)) quit(status = status)
