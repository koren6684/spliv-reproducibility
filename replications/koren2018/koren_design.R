make_koren_patterns <- function() {
  pat_uniform <- spliv_pattern(
    name = "Uniform direct effect",
    pattern = ~ 1,
    rationale = "The instrument is allowed to have the same direct effect in every observation.",
    pattern_type = "uniform",
    normalize = "max_abs"
  )

  pat_sparsebare <- spliv_pattern(
    name = "Sparse/bare direct-effect pattern",
    pattern = ~ sparsebare,
    rationale = paste(
      "The possible direct effect of drought on conflict is allowed to be larger",
      "where sparse or bare land cover is greater, capturing non-crop drought channels",
      "such as mobility, pastoral livelihoods, market access, or local resource stress."
    ),
    variables_used = "sparsebare",
    pattern_type = "theory_defined",
    normalize = "max_abs"
  )

  list(
    uniform = pat_uniform,
    sparsebare = pat_sparsebare
  )
}

make_koren_bpe_design_main <- function(dat) {
  if (!"sparsebare" %in% names(dat)) {
    stop("Main BPE design requires `sparsebare`, but it was not found.", call. = FALSE)
  }
  sparse_threshold <- stats::quantile(dat$sparsebare, 0.90, na.rm = TRUE)
  if (!is.finite(sparse_threshold)) {
    stop("Could not compute the 90th percentile of `sparsebare`.", call. = FALSE)
  }
  design <- bpe_design(
    name = "Top-decile sparse/bare land subset",
    subset = function(data) data$sparsebare >= sparse_threshold,
    rationale = paste(
      "The crop-yield treatment channel should be weak or absent in cells",
      "with very high sparse/bare land cover. This subset is defined before",
      "BPE diagnostics and is not selected by first-stage strength."
    ),
    variables_used = "sparsebare",
    subset_type = "theory_defined",
    pre_specified = TRUE,
    transportability_rationale = paste(
      "The direct drought-conflict relationship in high sparse/bare cells is assumed",
      "informative about possible direct drought-conflict channels in the target sample,",
      "with uncertainty carried forward."
    ),
    notes = paste0(
      "Pre-specified rule: sparsebare >= 90th percentile in the analysis sample. Threshold = ",
      signif(sparse_threshold, 4),
      ". This threshold is not tuned using first-stage diagnostics."
    )
  )

  list(
    design = design,
    design_role = "primary",
    design_id = "primary_top_decile_sparsebare",
    design_name = "Top-decile sparse/bare land subset",
    subset_rule = paste0("sparsebare >= 90th percentile (", signif(sparse_threshold, 6), ")"),
    design_variable = "sparsebare",
    threshold = as.numeric(sparse_threshold),
    threshold_label = "90th percentile"
  )
}

make_koren_bpe_design_primary <- function(dat) {
  make_koren_bpe_design_main(dat)
}

make_koren_bpe_design <- function(dat) {
  make_koren_bpe_design_main(dat)
}
