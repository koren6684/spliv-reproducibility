spliv_sim_pattern_objects <- function(data) {
  list(
    uniform = spliv_pattern(
      name = "Uniform pattern",
      pattern = ~ 1,
      rationale = "Uniform direct-effect benchmark.",
      pattern_type = "uniform"
    ),
    zone = spliv_pattern(
      name = "Zone pattern",
      pattern = ~ zone_pattern,
      rationale = "Direct effects are concentrated in a theory-defined high-exposure zone.",
      variables_used = "zone_pattern",
      pattern_type = "zone"
    ),
    gradient = spliv_pattern(
      name = "Exposure gradient",
      pattern = ~ exposure_gradient,
      rationale = "Direct effects increase with a generic alternative-channel exposure gradient.",
      variables_used = "exposure_gradient",
      pattern_type = "continuous"
    )
  )
}

spliv_sim_incorrect_pattern_key <- function(truth_pattern) {
  truth_pattern <- match.arg(truth_pattern, c("uniform", "zone", "gradient"))
  switch(
    truth_pattern,
    uniform = "zone",
    zone = "gradient",
    gradient = "zone"
  )
}

spliv_sim_bpe_design <- function() {
  bpe_design(
    name = "Theory-defined inactive subset",
    subset = ~ inactive_condition == 1,
    rationale = paste(
      "The instrument is expected not to shift the endogenous treatment in",
      "this subset because the simulated treatment channel is absent by design."
    ),
    variables_used = "inactive_condition",
    subset_type = "theory_defined",
    pre_specified = TRUE,
    transportability_rationale = paste(
      "The direct instrument-outcome path in the inactive subset is simulated",
      "to be informative about the target sample, subject to any transport gap",
      "built into the scenario."
    )
  )
}
