lelkes_bool_env <- function(name, default = TRUE) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) {
    return(default)
  }
  tolower(value) %in% c("1", "true", "t", "yes", "y")
}

lelkes_prepare_dirs <- function(script_dir) {
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

lelkes_package_path <- function(script_dir) {
  env_path <- Sys.getenv("SPLIV_PACKAGE_PATH", unset = "")
  if (nzchar(env_path)) {
    return(normalizePath(env_path, mustWork = TRUE))
  }
  normalizePath(file.path(script_dir, "..", "..", "..", "spliv"), mustWork = TRUE)
}

lelkes_load_spliv <- function(script_dir) {
  pkg_path <- lelkes_package_path(script_dir)
  loader <- NA_character_
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(pkg_path, quiet = TRUE)
    loader <- "devtools::load_all"
  } else if (requireNamespace("remotes", quietly = TRUE)) {
    remotes::install_local(pkg_path, quiet = TRUE, upgrade = "never")
    library(spliv)
    loader <- "remotes::install_local"
  } else {
    install.packages(pkg_path, repos = NULL, type = "source")
    library(spliv)
    loader <- "install.packages(source)"
  }
  version <- as.character(utils::packageVersion("spliv"))
  message("Loaded spliv from: ", pkg_path)
  message("spliv version: ", version)
  message("spliv loader: ", loader)
  list(path = pkg_path, version = version, loader = loader)
}
