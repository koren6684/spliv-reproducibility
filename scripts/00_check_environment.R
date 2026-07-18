#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)) else getwd()
repo_root <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
pkg_path <- Sys.getenv("SPLIV_PACKAGE_PATH", file.path(repo_root, "..", "spliv"))

cat("R:", R.version.string, "\n")
cat("Repository:", repo_root, "\n")
cat("Package path:", normalizePath(pkg_path, mustWork = FALSE), "\n")
if (!file.exists(file.path(pkg_path, "DESCRIPTION"))) {
  stop("The staged spliv package is missing. Set SPLIV_PACKAGE_PATH to its project-relative or absolute path.", call. = FALSE)
}

required <- c("devtools", "fixest", "testthat", "ggplot2", "knitr", "rmarkdown")
ok <- vapply(required, requireNamespace, logical(1), quietly = TRUE)
if (any(!ok)) {
  stop("Missing required R package(s): ", paste(names(ok)[!ok], collapse = ", "), call. = FALSE)
}
if (as.integer(parallel::detectCores(logical = TRUE)) > 2L) {
  cat("Note: validation commands should use at most two cores; simulation defaults are explicitly bounded in pilot mode.\n")
}
cat("Environment check passed.\n")
