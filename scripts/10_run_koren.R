#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)) else getwd()
repo_root <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
rep_root <- file.path(repo_root, "replications", "koren2018")
source(file.path(script_dir, "helpers_spliv_package.R"), local = FALSE)
spliv_load_package(report = FALSE)
setwd(repo_root)
status <- system2(Sys.which("Rscript"), c(shQuote(file.path(rep_root, "run_replication.R"))))
if (!identical(as.integer(status), 0L)) quit(status = status)
