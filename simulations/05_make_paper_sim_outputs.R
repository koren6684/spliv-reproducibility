#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) > 0L) {
  raw_path <- sub("^--file=", "", file_arg[1])
  raw_path <- gsub("~\\+~", " ", raw_path)
  script_path <- tryCatch(normalizePath(raw_path, mustWork = TRUE), error = function(e) raw_path)
  setwd(dirname(script_path))
}

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("The paper-output script requires ggplot2 for figures.")
}

paper_tables_dir <- file.path("tables", "paper")
paper_figures_dir <- file.path("figures", "paper")
dir.create(paper_tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(paper_figures_dir, recursive = TRUE, showWarnings = FALSE)

source_files <- c(
  patterned = file.path("tables", "full_patterned_summary.csv"),
  bpe = file.path("tables", "full_bpe_summary.csv"),
  stress = file.path("tables", "full_subgroup_search_stress.csv")
)
missing_sources <- source_files[!file.exists(source_files)]
if (length(missing_sources) > 0L) {
  stop("Missing required source CSV(s): ", paste(missing_sources, collapse = ", "))
}

ggplot2 <- asNamespace("ggplot2")

theme_paper <- function(base_size = 11) {
  ggplot2$theme_minimal(base_size = base_size) +
    ggplot2$theme(
      panel.grid.minor = ggplot2$element_blank(),
      legend.position = "bottom",
      strip.text = ggplot2$element_text(face = "bold"),
      plot.title = ggplot2$element_text(face = "bold"),
      plot.subtitle = ggplot2$element_text(color = "gray30")
    )
}

parse_p_number <- function(x, pattern) {
  token <- sub(pattern, "\\1", x)
  suppressWarnings(as.numeric(gsub("p", ".", token, fixed = TRUE)))
}

round_df <- function(df, digits = 3) {
  out <- df
  num_cols <- vapply(out, is.numeric, logical(1))
  out[num_cols] <- lapply(out[num_cols], function(x) round(x, digits))
  out
}

latex_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("_", "\\\\_", x)
  x <- gsub("%", "\\\\%", x)
  x <- gsub("&", "\\\\&", x)
  x
}

write_latex_table <- function(df, file, caption, label) {
  display <- round_df(df, digits = 3)
  display[] <- lapply(display, latex_escape)
  align <- paste0("l", paste(rep("r", max(0L, ncol(display) - 1L)), collapse = ""))
  lines <- c(
    "\\begin{table}[!htbp]",
    "\\centering",
    paste0("\\caption{", latex_escape(caption), "}"),
    paste0("\\label{", latex_escape(label), "}"),
    paste0("\\begin{tabular}{", align, "}"),
    "\\hline",
    paste(latex_escape(names(display)), collapse = " & "),
    "\\\\",
    "\\hline"
  )
  body <- apply(display, 1L, function(row) paste(row, collapse = " & "))
  lines <- c(lines, paste0(body, " \\\\"), "\\hline", "\\end{tabular}", "\\end{table}")
  writeLines(lines, file)
  invisible(file)
}

write_table_outputs <- function(df, stem, caption, label) {
  csv_file <- file.path(paper_tables_dir, paste0(stem, ".csv"))
  tex_file <- file.path(paper_tables_dir, paste0(stem, ".tex"))
  utils::write.csv(round_df(df, digits = 4), csv_file, row.names = FALSE)
  write_latex_table(df, tex_file, caption = caption, label = label)
  c(csv = csv_file, tex = tex_file)
}

save_plot <- function(plot, stem, width = 8, height = 5) {
  png_file <- file.path(paper_figures_dir, paste0(stem, ".png"))
  pdf_file <- file.path(paper_figures_dir, paste0(stem, ".pdf"))
  ggplot2$ggsave(png_file, plot, width = width, height = height, dpi = 300)
  ggplot2$ggsave(pdf_file, plot, width = width, height = height, device = grDevices::pdf)
  c(png = png_file, pdf = pdf_file)
}

patterned <- utils::read.csv(source_files[["patterned"]], stringsAsFactors = FALSE)
bpe <- utils::read.csv(source_files[["bpe"]], stringsAsFactors = FALSE)
stress <- utils::read.csv(source_files[["stress"]], stringsAsFactors = FALSE)

patterned$pi <- parse_p_number(patterned$scenario_name, ".*_pi_([0-9]p[0-9]+)_theta_.*")
patterned$theta_true <- parse_p_number(patterned$scenario_name, ".*_theta_([0-9]p[0-9]+)$")
patterned_at_truth <- patterned[abs(patterned$delta - patterned$theta_true) < 1e-10, , drop = FALSE]

patterned_wide <- reshape(
  patterned_at_truth[, c(
    "truth_pattern", "pi", "theta_true", "pattern_label",
    "coverage_rate", "mean_interval_width", "zero_inclusion_rate"
  )],
  idvar = c("truth_pattern", "pi", "theta_true"),
  timevar = "pattern_label",
  direction = "wide"
)
names(patterned_wide) <- sub("coverage_rate\\.", "coverage_", names(patterned_wide))
names(patterned_wide) <- sub("mean_interval_width\\.", "width_", names(patterned_wide))
names(patterned_wide) <- sub("zero_inclusion_rate\\.", "zero_", names(patterned_wide))

patterned_table <- patterned_wide[order(patterned_wide$truth_pattern, patterned_wide$pi, patterned_wide$theta_true), ]
patterned_table$width_ratio_correct_vs_default <- patterned_table$width_correct / patterned_table$width_default
patterned_table$width_ratio_correct_vs_uniform <- patterned_table$width_correct / patterned_table$width_uniform
patterned_table <- patterned_table[, c(
  "truth_pattern", "pi", "theta_true",
  "coverage_default", "coverage_uniform", "coverage_correct",
  "width_default", "width_uniform", "width_correct",
  "width_ratio_correct_vs_default", "width_ratio_correct_vs_uniform"
)]

patterned_table_files <- write_table_outputs(
  patterned_table,
  stem = "table_1_patterned_sensitivity_at_true_delta",
  caption = "Patterned sensitivity coverage and interval widths at the true direct-effect magnitude.",
  label = "tab:patterned-sensitivity"
)

ratio_df <- rbind(
  data.frame(
    truth_pattern = patterned_table$truth_pattern,
    pi = patterned_table$pi,
    theta_true = patterned_table$theta_true,
    reference = "Default",
    width_ratio = patterned_table$width_ratio_correct_vs_default,
    coverage_correct = patterned_table$coverage_correct,
    stringsAsFactors = FALSE
  ),
  data.frame(
    truth_pattern = patterned_table$truth_pattern,
    pi = patterned_table$pi,
    theta_true = patterned_table$theta_true,
    reference = "Uniform",
    width_ratio = patterned_table$width_ratio_correct_vs_uniform,
    coverage_correct = patterned_table$coverage_correct,
    stringsAsFactors = FALSE
  )
)
ratio_df$coverage_band <- ifelse(ratio_df$coverage_correct >= 0.95, "Coverage >= 0.95", "Coverage < 0.95")

patterned_plot <- ggplot2$ggplot(
  ratio_df,
  ggplot2$aes(
    x = factor(theta_true),
    y = width_ratio,
    shape = reference,
    color = coverage_band
  )
) +
  ggplot2$geom_hline(yintercept = 1, linetype = "dashed", color = "gray35") +
  ggplot2$geom_point(
    size = 3.1,
    position = ggplot2$position_dodge(width = 0.45)
  ) +
  ggplot2$geom_text(
    ggplot2$aes(label = sprintf("%.2f", coverage_correct)),
    position = ggplot2$position_dodge(width = 0.45),
    vjust = -0.85,
    size = 2.7,
    show.legend = FALSE
  ) +
  ggplot2$facet_grid(truth_pattern ~ pi, labeller = ggplot2$label_both) +
  ggplot2$scale_color_manual(values = c("Coverage >= 0.95" = "#2b8cbe", "Coverage < 0.95" = "#e34a33")) +
  ggplot2$scale_shape_manual(values = c("Default" = 16, "Uniform" = 17)) +
  ggplot2$coord_cartesian(ylim = c(0, max(1.15, max(ratio_df$width_ratio, na.rm = TRUE) * 1.15))) +
  ggplot2$labs(
    title = "Correct-pattern sensitivity tightens intervals relative to scalar benchmarks",
    subtitle = "Points show interval-width ratios at delta = theta_true; labels show coverage of the correct-pattern interval.",
    x = "True direct-effect magnitude",
    y = "Interval width ratio",
    shape = "Reference",
    color = "Coverage"
  ) +
  theme_paper()
patterned_plot_files <- save_plot(
  patterned_plot,
  stem = "figure_1_patterned_width_ratios",
  width = 9,
  height = 6.8
)

bpe$theta_true <- bpe$theta_inactive
bpe$gamma_bias <- bpe$mean_gamma_hat - bpe$theta_true
bpe$gamma_abs_error <- abs(bpe$gamma_bias)
bpe_table <- bpe[order(bpe$transport_gap, bpe$theta_true, bpe$pi_active, bpe$inactive_share), c(
  "transport_gap", "theta_true", "pi_active", "inactive_share",
  "eligibility_rate", "equivalence_rate", "mean_first_stage",
  "mean_first_stage_f", "mean_gamma_hat", "gamma_bias", "gamma_abs_error",
  "coverage_rate", "mean_interval_width"
)]
names(bpe_table)[names(bpe_table) == "pi_active"] <- "pi"

bpe_table_files <- write_table_outputs(
  bpe_table,
  stem = "table_2_bpe_design_performance",
  caption = "Confirmatory BPE performance by inactive-subset design scenario.",
  label = "tab:bpe-design-performance"
)

bpe_plot_df <- bpe
bpe_plot_df$theta_label <- factor(sprintf("%.2f", bpe_plot_df$theta_true))
bpe_plot_df$pi_label <- factor(sprintf("pi = %.1f", bpe_plot_df$pi_active))
bpe_plot_df$inactive_label <- factor(sprintf("inactive share = %.2f", bpe_plot_df$inactive_share))
bpe_y_min <- max(
  0,
  floor((min(bpe_plot_df$coverage_rate, na.rm = TRUE) - 0.02) * 20) / 20
)
bpe_plot <- ggplot2$ggplot(
  bpe_plot_df,
  ggplot2$aes(
    x = theta_label,
    y = coverage_rate,
    color = pi_label,
    shape = inactive_label,
    group = interaction(pi_label, inactive_label)
  )
) +
  ggplot2$geom_hline(yintercept = 0.95, linetype = "dashed", color = "gray35") +
  ggplot2$geom_line(linewidth = 0.55, alpha = 0.8) +
  ggplot2$geom_point(size = 2.7) +
  ggplot2$facet_wrap(~ transport_gap, labeller = ggplot2$label_both) +
  ggplot2$coord_cartesian(ylim = c(bpe_y_min, 1.01)) +
  ggplot2$labs(
    title = "BPE coverage across transport-gap and direct-effect scenarios",
    subtitle = "Dashed line marks 95 percent coverage.",
    x = "Theta in the inactive subset",
    y = "Coverage",
    color = "First stage",
    shape = "Subset size"
  ) +
  theme_paper()
bpe_plot_files <- save_plot(
  bpe_plot,
  stem = "figure_2_bpe_coverage",
  width = 8,
  height = 4.8
)

stress_table <- stress[order(stress$noise_multiplier, stress$pi, stress$subgroup_share, stress$K), c(
  "pi", "K", "subgroup_share", "noise_multiplier", "R",
  "false_f_rate", "false_ci_rate", "false_equivalence_rate",
  "mean_min_F", "median_min_F", "mean_selected_group_share",
  "mean_selected_group_pi_hat", "mean_selected_group_pi_se",
  "mean_selected_group_abs_t"
)]

stress_table_files <- write_table_outputs(
  stress_table,
  stem = "table_3_subgroup_search_stress",
  caption = "False inactive-subgroup selection under data-driven subgroup search.",
  label = "tab:subgroup-search-stress"
)

stress_plot_df <- stress
stress_plot_df$pi_label <- factor(sprintf("pi = %.2f", stress_plot_df$pi))
stress_plot_df$share_label <- factor(sprintf("share = %.2f", stress_plot_df$subgroup_share))
stress_plot_df$noise_label <- factor(sprintf("noise x %.0f", stress_plot_df$noise_multiplier))
stress_f_plot <- ggplot2$ggplot(
  stress_plot_df,
  ggplot2$aes(
    x = factor(K),
    y = false_f_rate,
    color = share_label,
    group = share_label
  )
) +
  ggplot2$geom_line(linewidth = 0.75) +
  ggplot2$geom_point(size = 2.4) +
  ggplot2$facet_grid(noise_label ~ pi_label) +
  ggplot2$coord_cartesian(ylim = c(0, 1)) +
  ggplot2$labs(
    title = "Data-driven subgroup search can find false weak-first-stage groups",
    subtitle = "No true inactive subset exists; each point is the probability that any searched group has F <= 5.",
    x = "Number of searched candidate groups",
    y = "False F <= 5 selection rate",
    color = "Candidate subgroup share"
  ) +
  theme_paper()
stress_f_plot_files <- save_plot(
  stress_f_plot,
  stem = "figure_3_subgroup_search_false_f",
  width = 8.5,
  height = 5.8
)

stress_long <- rbind(
  data.frame(
    pi = stress$pi,
    K = stress$K,
    subgroup_share = stress$subgroup_share,
    noise_multiplier = stress$noise_multiplier,
    rule = "CI includes zero",
    rate = stress$false_ci_rate,
    stringsAsFactors = FALSE
  ),
  data.frame(
    pi = stress$pi,
    K = stress$K,
    subgroup_share = stress$subgroup_share,
    noise_multiplier = stress$noise_multiplier,
    rule = "Loose equivalence pass",
    rate = stress$false_equivalence_rate,
    stringsAsFactors = FALSE
  )
)
stress_long$pi_label <- factor(sprintf("pi = %.2f", stress_long$pi))
stress_long$share_label <- factor(sprintf("share = %.2f", stress_long$subgroup_share))
stress_long$noise_label <- factor(sprintf("noise x %.0f", stress_long$noise_multiplier))
stress_appendix_plot <- ggplot2$ggplot(
  stress_long,
  ggplot2$aes(
    x = factor(K),
    y = rate,
    color = share_label,
    group = share_label
  )
) +
  ggplot2$geom_line(linewidth = 0.7) +
  ggplot2$geom_point(size = 2.1) +
  ggplot2$facet_grid(rule + noise_label ~ pi_label) +
  ggplot2$coord_cartesian(ylim = c(0, 1)) +
  ggplot2$labs(
    title = "Additional subgroup-search false-selection diagnostics",
    x = "Number of searched candidate groups",
    y = "False-selection rate",
    color = "Candidate subgroup share"
  ) +
  theme_paper(base_size = 10)
stress_appendix_plot_files <- save_plot(
  stress_appendix_plot,
  stem = "figure_A1_subgroup_search_ci_equivalence",
  width = 8.5,
  height = 8.5
)

readme_file <- file.path(paper_tables_dir, "README_paper_sim_outputs.md")
source_info <- data.frame(
  source = names(source_files),
  path = unname(source_files),
  rows = c(nrow(patterned), nrow(bpe), nrow(stress)),
  modified = vapply(source_files, function(path) {
    format(file.info(path)$mtime, "%Y-%m-%d %H:%M:%S %Z")
  }, character(1)),
  stringsAsFactors = FALSE
)

produced_files <- c(
  patterned_table_files,
  patterned_plot_files,
  bpe_table_files,
  bpe_plot_files,
  stress_table_files,
  stress_f_plot_files,
  stress_appendix_plot_files,
  readme = readme_file
)

readme_lines <- c(
  "# Paper Simulation Outputs",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "This is a post-processing pass only. It reads existing CSV summaries and does not run or modify any simulations.",
  "",
  "## Source CSVs",
  "",
  paste0("- `", source_info$path, "`: ", source_info$rows, " rows; modified ", source_info$modified),
  "",
  "## Produced Tables",
  "",
  paste0("- `", patterned_table_files[["csv"]], "` and `", patterned_table_files[["tex"]], "`"),
  paste0("- `", bpe_table_files[["csv"]], "` and `", bpe_table_files[["tex"]], "`"),
  paste0("- `", stress_table_files[["csv"]], "` and `", stress_table_files[["tex"]], "`"),
  "",
  "## Produced Figures",
  "",
  paste0("- `", patterned_plot_files[["png"]], "` and `", patterned_plot_files[["pdf"]], "`"),
  paste0("- `", bpe_plot_files[["png"]], "` and `", bpe_plot_files[["pdf"]], "`"),
  paste0("- `", stress_f_plot_files[["png"]], "` and `", stress_f_plot_files[["pdf"]], "`"),
  paste0("- `", stress_appendix_plot_files[["png"]], "` and `", stress_appendix_plot_files[["pdf"]], "`"),
  "",
  "## Notes",
  "",
  "- Patterned-sensitivity table and figure use rows where `delta = theta_true`.",
  "- Patterned figure compares the correct-pattern interval width against the default and uniform sensitivity benchmarks; labels report correct-pattern coverage.",
  "- BPE table reports gamma recovery as `mean_gamma_hat - theta_true`.",
  "- Subgroup-search stress outputs use `full_subgroup_search_stress.csv`; the older blank subgroup-search summary is not used."
)
writeLines(readme_lines, readme_file)

cat("Paper simulation outputs generated.\n")
cat("Source CSVs:\n")
print(source_info, row.names = FALSE)
cat("\nProduced files:\n")
print(unname(produced_files))
