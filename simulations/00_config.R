spliv_sim_detect_cores <- function() {
  cores <- suppressWarnings(parallel::detectCores(logical = TRUE))
  if (!is.finite(cores) || is.na(cores)) {
    return(1L)
  }
  as.integer(max(1L, cores - 1L))
}

.spliv_sim_env_int <- function(name, default = NULL) {
  value <- Sys.getenv(name, "")
  if (!nzchar(value)) {
    return(default)
  }
  parsed <- suppressWarnings(as.integer(value))
  if (is.na(parsed)) {
    stop("Environment variable `", name, "` must be an integer.")
  }
  parsed
}

.spliv_sim_env_num <- function(name, default = NULL) {
  value <- Sys.getenv(name, "")
  if (!nzchar(value)) {
    return(default)
  }
  parsed <- suppressWarnings(as.numeric(value))
  if (is.na(parsed)) {
    stop("Environment variable `", name, "` must be numeric.")
  }
  parsed
}

.spliv_sim_env_flag <- function(name, default = FALSE) {
  value <- Sys.getenv(name, "")
  if (!nzchar(value)) {
    return(default)
  }
  value <- tolower(trimws(value))
  if (value %in% c("1", "true", "t", "yes", "y")) {
    return(TRUE)
  }
  if (value %in% c("0", "false", "f", "no", "n")) {
    return(FALSE)
  }
  stop("Environment variable `", name, "` must be TRUE/FALSE-like.")
}

spliv_sim_profiles <- function() {
  list(
    pilot = list(
      profile = "pilot",
      R = 10L,
      nx = 8L,
      ny = 8L,
      T_periods = 8L,
      delta_grid = seq(0, 0.20, by = 0.05),
      n_cores = 1L,
      base_seed = 20260424L
    ),
    moderate = list(
      profile = "moderate",
      R = 100L,
      nx = 12L,
      ny = 12L,
      T_periods = 20L,
      delta_grid = seq(0, 0.20, by = 0.025),
      n_cores = spliv_sim_detect_cores(),
      base_seed = 20260424L
    ),
    full = list(
      profile = "full",
      R = 1000L,
      nx = 20L,
      ny = 20L,
      T_periods = 50L,
      delta_grid = seq(0, 0.20, by = 0.01),
      n_cores = spliv_sim_detect_cores(),
      base_seed = 20260424L
    )
  )
}

spliv_sim_resolve_config <- function(paths,
                                     profile = Sys.getenv("SPLIV_SIM_PROFILE", "pilot")) {
  profiles <- spliv_sim_profiles()
  if (!profile %in% names(profiles)) {
    stop(
      "`SPLIV_SIM_PROFILE` must be one of: ",
      paste(names(profiles), collapse = ", ")
    )
  }

  config <- profiles[[profile]]
  config$R <- .spliv_sim_env_int("SPLIV_SIM_R", config$R)
  config$nx <- .spliv_sim_env_int("SPLIV_SIM_NX", config$nx)
  config$ny <- .spliv_sim_env_int("SPLIV_SIM_NY", config$ny)
  config$T_periods <- .spliv_sim_env_int("SPLIV_SIM_T", config$T_periods)
  config$n_cores <- .spliv_sim_env_int("SPLIV_SIM_CORES", config$n_cores)
  config$base_seed <- .spliv_sim_env_int("SPLIV_SIM_SEED", config$base_seed)
  config$overwrite <- .spliv_sim_env_flag("SPLIV_SIM_OVERWRITE", FALSE)

  delta_max <- .spliv_sim_env_num(
    "SPLIV_SIM_DELTA_MAX",
    max(config$delta_grid)
  )
  delta_step <- .spliv_sim_env_num(
    "SPLIV_SIM_DELTA_STEP",
    if (length(config$delta_grid) > 1L) {
      config$delta_grid[2] - config$delta_grid[1]
    } else {
      0.05
    }
  )
  if (!is.finite(delta_max) || delta_max < 0) {
    stop("Resolved `delta_max` must be non-negative.")
  }
  if (!is.finite(delta_step) || delta_step <= 0) {
    stop("Resolved `delta_step` must be strictly positive.")
  }
  config$delta_grid <- round(seq(0, delta_max, by = delta_step), 10)

  output_root <- Sys.getenv("SPLIV_SIM_OUTPUT_DIR", "")
  if (!nzchar(output_root)) {
    output_root <- paths$output_root
  } else if (!grepl("^(/|[A-Za-z]:)", output_root)) {
    output_root <- file.path(paths$root, output_root)
  }

  config$output_root <- .spliv_sim_normalize_path(output_root)
  config$output_dir <- file.path(config$output_root, config$profile)
  config$raw_dir <- file.path(config$output_dir, "raw")
  config$run_timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")

  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(config$raw_dir, recursive = TRUE, showWarnings = FALSE)

  config
}

spliv_sim_print_config <- function(config) {
  cat("Resolved simulation config:\n")
  cat("  profile:", config$profile, "\n")
  cat("  R:", config$R, "\n")
  cat("  nx:", config$nx, "\n")
  cat("  ny:", config$ny, "\n")
  cat("  T_periods:", config$T_periods, "\n")
  cat("  delta_grid:", paste(config$delta_grid, collapse = ", "), "\n")
  cat("  n_cores:", config$n_cores, "\n")
  cat("  overwrite:", config$overwrite, "\n")
  cat("  base_seed:", config$base_seed, "\n")
  cat("  output_dir:", config$output_dir, "\n")
  invisible(config)
}
