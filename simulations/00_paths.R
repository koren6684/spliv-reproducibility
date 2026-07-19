"%||%" <- function(x, y) {
  if (is.null(x)) {
    y
  } else {
    x
  }
}

.spliv_sim_normalize_path <- function(path) {
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

spliv_sim_script_dir <- function(default = getwd()) {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(.spliv_sim_normalize_path(sub("^--file=", "", file_arg[1]))))
  }

  frames <- sys.frames()
  for (idx in rev(seq_along(frames))) {
    ofile <- frames[[idx]]$ofile
    if (!is.null(ofile)) {
      return(dirname(.spliv_sim_normalize_path(ofile)))
    }
  }

  .spliv_sim_normalize_path(default)
}

spliv_sim_resolve_root <- function(root = NULL) {
  env_root <- Sys.getenv("SPLIV_SIMS_ROOT", "")
  candidates <- unique(Filter(
    nzchar,
    c(root, env_root, getwd(), spliv_sim_script_dir())
  ))

  for (candidate in candidates) {
    candidate <- .spliv_sim_normalize_path(candidate)
    if (file.exists(file.path(candidate, "00_config.R")) &&
        dir.exists(file.path(candidate, "R"))) {
      return(candidate)
    }
  }

  stop(
    "Could not resolve the sims root. Set `SPLIV_SIMS_ROOT` or run from the sims directory."
  )
}

spliv_sim_paths <- function(root = NULL) {
  root <- spliv_sim_resolve_root(root)
  paths <- list(
    root = root,
    r_dir = file.path(root, "R"),
    output_root = file.path(root, "output"),
    logs_dir = file.path(root, "logs"),
    figures_dir = file.path(root, "figures"),
    tables_dir = file.path(root, "tables")
  )

  dir.create(paths$output_root, recursive = TRUE, showWarnings = FALSE)
  dir.create(paths$logs_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(paths$figures_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(paths$tables_dir, recursive = TRUE, showWarnings = FALSE)
  paths
}

spliv_sim_source_helpers <- function(paths) {
  shared_loader <- file.path(dirname(paths$root), "scripts", "helpers_spliv_package.R")
  if (!file.exists(shared_loader)) {
    stop("Could not find the shared SPLIV package loader at `", shared_loader, ".`", call. = FALSE)
  }
  source(shared_loader, local = FALSE)

  helper_files <- c(
    "helpers_packages.R",
    "helpers_parallel.R",
    "helpers_checkpointing.R",
    "helpers_dgp.R",
    "helpers_patterns.R",
    "helpers_estimators.R",
    "helpers_plots.R",
    "helpers_summaries.R"
  )

  for (helper_file in helper_files) {
    source(file.path(paths$r_dir, helper_file), local = FALSE)
  }

  invisible(paths)
}
