#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)) else getwd()
repo_root <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
source(file.path(script_dir, "helpers_spliv_package.R"), local = FALSE)
package_info <- spliv_load_package(report = TRUE)

cat("R:", R.version.string, "\n")
cat("Repository:", repo_root, "\n")
cat("Package path:", package_info$path, "\n")

required <- c("spliv", "fixest", "testthat", "ggplot2", "knitr", "rmarkdown")
ok <- vapply(required, requireNamespace, logical(1), quietly = TRUE)
if (any(!ok)) {
  stop("Missing required R package(s): ", paste(names(ok)[!ok], collapse = ", "), call. = FALSE)
}
if (as.integer(parallel::detectCores(logical = TRUE)) > 2L) {
  cat("Note: validation commands should use at most two cores; simulation defaults are explicitly bounded in pilot mode.\n")
}
cat("Environment check passed.\n")
