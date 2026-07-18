# Paper Simulation Outputs

Generated: 2026-04-27 08:59:37 EDT

This is a post-processing pass only. It reads existing CSV summaries and does not run or modify any simulations.

The final paper simulation design uses **R = 1,000 Monte Carlo replications**.
The reference tables and figures below are the small, post-processed outputs
for that design; `table_3_subgroup_search_stress.csv` records the replication
count explicitly in its `R` column.

## Source CSVs

- `tables/full_patterned_summary.csv`: 1512 rows; modified 2026-04-25 03:43:21 EDT
- `tables/full_bpe_summary.csv`: 24 rows; modified 2026-04-26 10:17:27 EDT
- `tables/full_subgroup_search_stress.csv`: 72 rows; modified 2026-04-26 13:40:14 EDT

## Produced Tables

- `tables/paper/table_1_patterned_sensitivity_at_true_delta.csv` and `tables/paper/table_1_patterned_sensitivity_at_true_delta.tex`
- `tables/paper/table_2_bpe_design_performance.csv` and `tables/paper/table_2_bpe_design_performance.tex`
- `tables/paper/table_3_subgroup_search_stress.csv` and `tables/paper/table_3_subgroup_search_stress.tex`

## Produced Figures

- `figures/paper/figure_1_patterned_width_ratios.png` and `figures/paper/figure_1_patterned_width_ratios.pdf`
- `figures/paper/figure_2_bpe_coverage.png` and `figures/paper/figure_2_bpe_coverage.pdf`
- `figures/paper/figure_3_subgroup_search_false_f.png` and `figures/paper/figure_3_subgroup_search_false_f.pdf`
- `figures/paper/figure_A1_subgroup_search_ci_equivalence.png` and `figures/paper/figure_A1_subgroup_search_ci_equivalence.pdf`

## Notes

- Patterned-sensitivity table and figure use rows where `delta = theta_true`.
- Patterned figure compares the correct-pattern interval width against the default and uniform sensitivity benchmarks; labels report correct-pattern coverage.
- BPE table reports gamma recovery as `mean_gamma_hat - theta_true`.
- Subgroup-search stress outputs use `full_subgroup_search_stress.csv`; the older blank subgroup-search summary is not used.
