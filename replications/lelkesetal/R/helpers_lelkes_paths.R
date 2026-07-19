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

lelkes_load_spliv <- function(script_dir) {
  if (!exists("spliv_load_package", mode = "function", inherits = FALSE)) {
    shared_loader <- file.path(dirname(dirname(script_dir)), "scripts", "helpers_spliv_package.R")
    source(shared_loader, local = FALSE)
  }
  spliv_load_package(report = TRUE)
}
