make_lelkes_patterns <- function() {
  list(
    uniform = spliv_pattern(
      name = "Uniform direct effect",
      pattern = ~ 1,
      rationale = "The instrument is allowed to have the same direct effect in every observation.",
      pattern_type = "uniform",
      normalize = "max_abs"
    ),
    density = spliv_pattern(
      name = "Density direct-effect pattern",
      pattern = ~ density,
      rationale = paste(
        "The possible direct effect of the ROW instrument on affective polarization",
        "is allowed to be larger in denser counties, where broadband markets,",
        "media environments, and political communication exposure are more developed."
      ),
      variables_used = "density",
      pattern_type = "theory_defined",
      normalize = "max_abs"
    )
  )
}

make_lelkes_bpe_design <- function(dat) {
  if (!"density" %in% names(dat)) {
    stop("BPE design requires variable `density`, but it was not found.", call. = FALSE)
  }

  density_threshold <- stats::quantile(dat$density, 0.10, na.rm = TRUE)

  design <- bpe_design(
    name = "Bottom-decile population density subset",
    subset = function(data) data$density <= density_threshold,
    rationale = paste(
      "The ROW instrument should have weaker or no first-stage effect on provider entry",
      "in very low-density counties, where market size and infrastructure economics",
      "plausibly dominate right-of-way regulatory costs. This subset is defined before",
      "BPE diagnostics and is not selected by first-stage strength."
    ),
    variables_used = "density",
    subset_type = "theory_defined",
    pre_specified = TRUE,
    transportability_rationale = paste(
      "The direct ROW-polarization relationship in low-density counties is assumed",
      "informative about possible direct ROW-polarization channels in the target sample,",
      "but this transportability assumption is substantively stronger than in the Koren",
      "application and should be interpreted cautiously."
    ),
    notes = paste0(
      "Pre-specified rule: density <= 10th percentile. Threshold = ",
      signif(density_threshold, 4),
      ". This threshold is not tuned using first-stage diagnostics."
    )
  )

  list(design = design, threshold = density_threshold)
}
