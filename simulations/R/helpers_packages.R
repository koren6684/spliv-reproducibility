spliv_sim_load_package <- function(paths) {
  pkg_path <- Sys.getenv("SPLIV_PACKAGE_PATH", "")
  if (!nzchar(pkg_path)) {
    pkg_path <- spliv_sim_relative_package_path(paths)
  }
  pkg_path <- .spliv_sim_normalize_path(pkg_path)

  if (!file.exists(file.path(pkg_path, "DESCRIPTION"))) {
    stop(
      "Could not find the local SPLIV package at `", pkg_path,
      "`. Set `SPLIV_PACKAGE_PATH` to override the default `../spliv` path."
    )
  }

  loader <- NULL
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(pkg_path, quiet = TRUE, export_all = FALSE)
    loader <- "devtools::load_all"
  } else if (requireNamespace("remotes", quietly = TRUE)) {
    remotes::install_local(
      pkg_path,
      upgrade = "never",
      dependencies = FALSE,
      quiet = TRUE
    )
    loader <- "remotes::install_local"
  } else {
    install.packages(pkg_path, repos = NULL, type = "source")
    loader <- "install.packages(type = 'source')"
  }

  suppressPackageStartupMessages(library(spliv))

  version <- tryCatch(
    as.character(utils::packageVersion("spliv")),
    error = function(e) NA_character_
  )

  info <- list(
    path = pkg_path,
    version = version,
    loader = loader
  )

  cat("Loaded SPLIV package:\n")
  cat("  path:", info$path, "\n")
  cat("  version:", info$version, "\n")
  cat("  loader:", info$loader, "\n")

  invisible(info)
}
