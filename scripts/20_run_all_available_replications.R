#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)) else getwd()
repo_root <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
pkg_path <- Sys.getenv("SPLIV_PACKAGE_PATH", file.path(repo_root, "..", "spliv"))
Sys.setenv(SPLIV_PACKAGE_PATH = normalizePath(pkg_path, mustWork = TRUE))

run_if_configured <- function(label, env_name, script) {
  value <- Sys.getenv(env_name, "")
  if (!nzchar(value)) {
    cat(label, " skipped: set ", env_name, " to run this replication.\n", sep = "")
    return(invisible(TRUE))
  }
  cat(label, " configured at ", value, "; running.\n", sep = "")
  status <- system2(Sys.which("Rscript"), c("--vanilla", shQuote(script)))
  if (!identical(as.integer(status), 0L)) quit(status = status)
  invisible(TRUE)
}

run_if_configured("Koren", "SPLIV_KOREN_DATA", file.path(repo_root, "scripts", "10_run_koren.R"))
run_if_configured("Lelkes", "SPLIV_LELKES_DATA", file.path(repo_root, "scripts", "11_run_lelkes.R"))
cat("Available-replications dispatch complete; no data were substituted.\n")
