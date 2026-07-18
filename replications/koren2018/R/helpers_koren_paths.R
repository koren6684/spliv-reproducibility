`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

koren_bool_env <- function(name, default = TRUE) {
  val <- Sys.getenv(name, if (isTRUE(default)) "TRUE" else "FALSE")
  tolower(trimws(val)) %in% c("1", "true", "t", "yes", "y")
}

koren_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)))
  }
  normalizePath(getwd(), mustWork = TRUE)
}

koren_prepare_dirs <- function(script_dir) {
  out_root <- file.path(script_dir, "output")
  dirs <- list(
    root = out_root,
    tables = file.path(out_root, "tables"),
    figures = file.path(out_root, "figures"),
    logs = file.path(out_root, "logs")
  )
  for (d in dirs) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
  dirs
}

koren_update_root <- function(script_dir) {
  normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
}

koren_project_root <- function(script_dir) {
  dirname(koren_update_root(script_dir))
}

koren_resolve_package_path <- function(script_dir) {
  env_path <- Sys.getenv("SPLIV_PACKAGE_PATH", "")
  if (nzchar(env_path)) {
    return(normalizePath(env_path, mustWork = TRUE))
  }
  normalizePath(file.path(script_dir, "..", "..", "..", "spliv"), mustWork = TRUE)
}

koren_load_spliv <- function(script_dir) {
  pkg_path <- koren_resolve_package_path(script_dir)
  desc_path <- file.path(pkg_path, "DESCRIPTION")
  version <- if (file.exists(desc_path)) {
    as.character(read.dcf(desc_path)[1, "Version"])
  } else {
    NA_character_
  }

  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(pkg_path, quiet = TRUE, export_all = FALSE)
    loader <- "devtools::load_all"
  } else if (requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(pkg_path, quiet = TRUE, export_all = FALSE)
    loader <- "pkgload::load_all"
  } else {
    stop(
      "Neither `devtools` nor `pkgload` is installed. Install one of them or set up `spliv` on `.libPaths()`.",
      call. = FALSE
    )
  }

  message("Loaded spliv from: ", pkg_path)
  message("spliv version: ", version)
  message("spliv loader: ", loader)
  list(path = pkg_path, version = version, loader = loader)
}

koren_find_panel_data <- function(script_dir) {
  env_path <- Sys.getenv("SPLIV_KOREN_DATA", "")
  if (!nzchar(env_path)) {
    env_path <- Sys.getenv("KOREN_PANEL_DATA_PATH", "")
  }
  if (nzchar(env_path)) {
    if (!file.exists(env_path)) {
      stop("`SPLIV_KOREN_DATA`/`KOREN_PANEL_DATA_PATH` was set, but the Koren file does not exist: ", env_path, call. = FALSE)
    }
    return(normalizePath(env_path, mustWork = TRUE))
  }
  stop(
    "Koren replication data are not configured. Obtain the restricted study data, " ,
    "set `SPLIV_KOREN_DATA` to the full path of `crop.dat.af.rnr3.dta`, and run `Rscript check_data.R`.",
    call. = FALSE
  )
}
