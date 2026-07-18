spliv_sim_now <- function() {
  format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
}

spliv_sim_sanitize_name <- function(x) {
  x <- tolower(as.character(x))
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}

spliv_sim_task_seed <- function(base_seed, scenario_id, replicate_id) {
  seed <- as.integer((base_seed + scenario_id * 100000L + replicate_id) %% .Machine$integer.max)
  if (!is.finite(seed) || seed <= 0L) {
    seed <- as.integer(abs(base_seed + scenario_id + replicate_id) + 1L)
  }
  seed
}

spliv_sim_result_file <- function(config, family, scenario_name, replicate_id) {
  file.path(
    config$raw_dir,
    sprintf(
      "%s__%s__rep%04d__%s.rds",
      spliv_sim_sanitize_name(family),
      spliv_sim_sanitize_name(scenario_name),
      as.integer(replicate_id),
      config$profile
    )
  )
}

spliv_sim_save_result <- function(result, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(result, file = file)
  invisible(file)
}

spliv_sim_collect_files <- function(config, family) {
  pattern <- paste0("^", spliv_sim_sanitize_name(family), "__.*__", config$profile, "\\.rds$")
  sort(list.files(config$raw_dir, pattern = pattern, full.names = TRUE))
}

spliv_sim_collect_results <- function(config, family) {
  files <- spliv_sim_collect_files(config, family)
  lapply(files, readRDS)
}

spliv_sim_rbind_fill <- function(dfs) {
  dfs <- Filter(function(x) is.data.frame(x) && nrow(x) > 0L, dfs)
  if (length(dfs) == 0L) {
    return(data.frame())
  }

  all_names <- unique(unlist(lapply(dfs, names), use.names = FALSE))
  filled <- lapply(dfs, function(df) {
    missing_names <- setdiff(all_names, names(df))
    for (name in missing_names) {
      df[[name]] <- NA
    }
    df[all_names]
  })
  out <- do.call(rbind, filled)
  rownames(out) <- NULL
  out
}

spliv_sim_nearest_delta <- function(delta_grid, value) {
  delta_grid[[which.min(abs(delta_grid - value))]]
}

spliv_sim_make_tasks <- function(scenarios, R) {
  tasks <- vector("list", length = nrow(scenarios) * R)
  idx <- 1L
  for (scenario_idx in seq_len(nrow(scenarios))) {
    for (rep_id in seq_len(R)) {
      tasks[[idx]] <- list(
        scenario_idx = scenario_idx,
        replicate_id = rep_id
      )
      idx <- idx + 1L
    }
  }
  tasks
}

spliv_sim_new_log_file <- function(paths, run_name, profile) {
  file.path(
    paths$logs_dir,
    sprintf(
      "%s_%s_%s.log",
      spliv_sim_sanitize_name(run_name),
      spliv_sim_sanitize_name(profile),
      format(Sys.time(), "%Y%m%d_%H%M%S")
    )
  )
}

spliv_sim_log <- function(log_file, ...) {
  message_text <- paste(..., collapse = "")
  line <- paste0("[", spliv_sim_now(), "] ", message_text)
  cat(line, "\n")
  cat(line, "\n", file = log_file, append = TRUE)
  invisible(line)
}

spliv_sim_log_header <- function(log_file, run_name, config, package_info, scenarios = NULL) {
  spliv_sim_log(log_file, "Run: ", run_name)
  spliv_sim_log(log_file, "Profile: ", config$profile)
  spliv_sim_log(log_file, "Package path: ", package_info$path)
  spliv_sim_log(log_file, "Package version: ", package_info$version)
  spliv_sim_log(log_file, "Package loader: ", package_info$loader)
  spliv_sim_log(log_file, "Output dir: ", config$output_dir)
  spliv_sim_log(log_file, "Overwrite existing results: ", config$overwrite)
  spliv_sim_log(
    log_file,
    "Resolved config: ",
    paste(
      c(
        paste0("R=", config$R),
        paste0("nx=", config$nx),
        paste0("ny=", config$ny),
        paste0("T_periods=", config$T_periods),
        paste0("delta_grid=[", paste(config$delta_grid, collapse = ", "), "]"),
        paste0("n_cores=", config$n_cores),
        paste0("base_seed=", config$base_seed)
      ),
      collapse = "; "
    )
  )
  if (!is.null(scenarios) && nrow(scenarios) > 0L) {
    spliv_sim_log(
      log_file,
      "Scenarios: ",
      paste(scenarios$scenario_name, collapse = ", ")
    )
  }
}

spliv_sim_log_footer <- function(log_file, counters) {
  spliv_sim_log(
    log_file,
    "Completed: ", counters$completed,
    "; skipped: ", counters$skipped,
    "; failed: ", counters$failed
  )
}

spliv_sim_count_status <- function(status_rows) {
  status_values <- vapply(status_rows, function(x) x$status %||% "failed", character(1))
  list(
    completed = sum(status_values == "completed"),
    skipped = sum(status_values == "skipped"),
    failed = sum(status_values == "failed")
  )
}
