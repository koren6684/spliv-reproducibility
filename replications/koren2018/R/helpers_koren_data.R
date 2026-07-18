koren_read_panel_data <- function(path) {
  if (requireNamespace("foreign", quietly = TRUE)) {
    return(foreign::read.dta(path))
  }
  if (requireNamespace("haven", quietly = TRUE)) {
    return(as.data.frame(haven::read_dta(path)))
  }
  stop("Reading Stata data requires either the `foreign` or `haven` package.", call. = FALSE)
}

koren_required_variables <- function() {
  c(
    "acled_inc_sum", "maize_yield", "wheat_yield", "spi6", "gid", "year",
    "sparsebare"
  )
}

koren_validate_variables <- function(dat) {
  missing <- setdiff(koren_required_variables(), names(dat))
  if (length(missing) > 0) {
    if ("sparsebare" %in% missing) {
      land_like <- grep(
        "sparse|bare|crop|land|cover|forest|water|irri|mnt|grass|pasture",
        names(dat),
        ignore.case = TRUE,
        value = TRUE
      )
      stop(
        "Required Koren land-cover variable missing: sparsebare",
        ". Possible land-cover variables available are:\n",
        paste0("- ", land_like, collapse = "\n"),
        call. = FALSE
      )
    }
    stop(
      "Missing required Koren Table 3 variable(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

koren_table3_specs <- function() {
  data.frame(
    crop = c("maize", "wheat"),
    treatment = c("maize_yield", "wheat_yield"),
    outcome = "acled_inc_sum",
    instrument = "spi6",
    fe = "gid + year",
    cluster = "gid",
    stringsAsFactors = FALSE
  )
}

koren_analysis_data <- function(dat, treatment) {
  vars <- c("acled_inc_sum", treatment, "spi6", "gid", "year", "sparsebare")
  keep <- stats::complete.cases(dat[, vars])
  out <- dat[keep, , drop = FALSE]
  if (!nrow(out)) {
    stop("No complete observations remain for treatment `", treatment, "`.", call. = FALSE)
  }
  out
}

koren_iv_formula <- function(treatment) {
  stats::as.formula(paste0("acled_inc_sum ~ ", treatment, " | spi6"))
}
