lelkes_old_source_roots <- function(script_dir) {
  data_root <- Sys.getenv("SPLIV_LELKES_DATA", "")
  if (!nzchar(data_root)) character(0) else normalizePath(data_root, mustWork = FALSE)
}

lelkes_find_file_case_insensitive <- function(root, filename) {
  if (!dir.exists(root)) {
    return(NA_character_)
  }
  files <- list.files(root, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  hit <- files[tolower(basename(files)) == tolower(filename)]
  if (length(hit)) {
    normalizePath(hit[[1]], mustWork = TRUE)
  } else {
    NA_character_
  }
}

lelkes_find_data_files <- function(script_dir) {
  roots <- lelkes_old_source_roots(script_dir)
  merged_file <- county_file <- NA_character_
  for (root in roots) {
    if (is.na(merged_file)) {
      merged_file <- lelkes_find_file_case_insensitive(root, "mergeddataset.RData")
    }
    if (is.na(county_file)) {
      county_file <- lelkes_find_file_case_insensitive(root, "mergedcounty.RData")
    }
  }
  if (is.na(merged_file) || is.na(county_file)) {
    stop(
      "Lelkes replication data are not configured. Obtain the restricted study data, " ,
      "set `SPLIV_LELKES_DATA` to the directory containing `mergeddataset.RData` and `mergedcounty.RData`, " ,
      "and run `Rscript check_data.R`.",
      call. = FALSE
    )
  }
  list(merged = merged_file, county = county_file)
}

zero1 <- function(x, minx = NA, maxx = NA) {
  if (is.na(minx)) {
    return((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
  }
  (x - minx) / (maxx - minx)
}

lelkes_required_variables <- function() {
  c(
    "infeels", "outfeels", "providers", "Total", "year", "region",
    "percent_black", "percent_white", "percent_male", "lowed",
    "unemploymentrate", "density", "HHINC", "state"
  )
}

lelkes_load_raw_data <- function(data_files) {
  env <- new.env(parent = emptyenv())
  load(data_files$merged, envir = env)
  load(data_files$county, envir = env)
  if (!exists("merged", envir = env, inherits = FALSE)) {
    stop("`mergeddataset.RData` did not contain object `merged`.", call. = FALSE)
  }
  if (!exists("countyproviders", envir = env, inherits = FALSE)) {
    stop("`mergedcounty.RData` did not contain object `countyproviders`.", call. = FALSE)
  }
  list(
    merged = get("merged", envir = env),
    countyproviders = get("countyproviders", envir = env)
  )
}

lelkes_prepare_analysis_data <- function(raw) {
  merged <- raw$merged
  missing <- setdiff(lelkes_required_variables(), names(merged))
  if (length(missing)) {
    stop(
      "Lelkes merged data are missing required variables: ",
      paste(missing, collapse = ", "),
      ". Available variables: ",
      paste(names(merged), collapse = ", "),
      call. = FALSE
    )
  }

  dat <- merged[, c(1:16, 21:22)]
  dat$affective_polarization <- zero1(dat$infeels - dat$outfeels)
  dat$log_providers <- log(dat$providers)
  dat$log_Total <- log(dat$Total)
  dat$log_HHINC <- log(dat$HHINC)

  analysis_vars <- c(
    "affective_polarization", "log_providers", "log_Total",
    "year", "region", "percent_black", "percent_white", "percent_male",
    "lowed", "unemploymentrate", "density", "log_HHINC", "state",
    "infeels", "outfeels", "providers", "Total", "HHINC"
  )
  dat <- stats::na.omit(dat[, analysis_vars])
  dat$state <- droplevels(as.factor(dat$state))
  dat$region <- droplevels(as.factor(dat$region))
  dat
}
