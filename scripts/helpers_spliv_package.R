# Shared SPLIV package loader for all public reproducibility scripts.

spliv_load_package <- function(report = TRUE) {
  source_path <- Sys.getenv("SPLIV_PACKAGE_PATH", unset = "")

  if (nzchar(source_path)) {
    source_path <- normalizePath(source_path, mustWork = TRUE)
    if (!requireNamespace("devtools", quietly = TRUE)) {
      stop(
        "SPLIV_PACKAGE_PATH is set, but `devtools` is not installed. ",
        "Install devtools or unset SPLIV_PACKAGE_PATH to use the installed package.",
        call. = FALSE
      )
    }
    devtools::load_all(source_path, quiet = TRUE)
    loader <- "devtools::load_all"
    package_path <- source_path
    source_project_path <- source_path
  } else {
    if (!requireNamespace("spliv", quietly = TRUE)) {
      stop(
        "The spliv package is not installed. Run renv::restore() from the reproducibility repository root.",
        call. = FALSE
      )
    }
    suppressPackageStartupMessages(library(spliv))
    package_path <- normalizePath(find.package("spliv"), mustWork = TRUE)
    loader <- "installed package"
    source_project_path <- NA_character_
  }

  version <- packageVersion("spliv")
  if (version < package_version("0.1.0")) {
    stop(
      "The installed spliv package is version ", as.character(version),
      "; version 0.1.0 or newer is required.",
      call. = FALSE
    )
  }

  info <- list(
    path = package_path,
    version = as.character(version),
    loader = loader,
    source_project_path = source_project_path
  )

  if (isTRUE(report)) {
    cat("Loaded SPLIV package:\n")
    cat("  path:", info$path, "\n")
    cat("  version:", info$version, "\n")
    cat("  loader:", info$loader, "\n")
    if (!is.na(info$source_project_path)) {
      cat("  source project:", info$source_project_path, "\n")
    }
  }

  invisible(info)
}
