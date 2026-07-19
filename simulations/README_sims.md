# SPLIV simulations

This directory contains portable simulation code for patterned sensitivity,
confirmatory BPE validation, and the negative subgroup-search stress test. All
randomness is seeded; pilot mode uses `R = 10`, an 8-by-8 grid, eight periods,
one core, and base seed `20260424`. The full profile uses `R = 1000`, a 20-by-
20 grid, 50 periods, and a delta step of `0.01`; it can take many hours or
days and writes checkpoint `.rds` files under ignored output directories.

`00_paths.R` and `00_config.R` resolve project-relative paths and environment
overrides. The `R/` directory contains DGP, estimator, summary, plotting,
parallel, and checkpoint helpers. `01_smoke_test.R` is the smallest end-to-end
check; `04_run_all_pilot.R` runs all pilot families; `05_make_paper_sim_outputs.R`
builds publication tables/figures after full summaries exist; and
`99_run_all_full.R` is the production runner.

Run from this directory, or use the wrappers in the repository-level `scripts/`
directory:

```bash
Rscript 01_smoke_test.R
SPLIV_SIM_PROFILE=pilot SPLIV_SIM_CORES=1 Rscript 04_run_all_pilot.R
# Long-running production command:
SPLIV_SIM_PROFILE=full Rscript 99_run_all_full.R
```

The simulation scripts load the installed `spliv` package by default. Set
`SPLIV_PACKAGE_PATH` only as a development override for a source checkout; it
must point to an existing package project and requires `devtools`. Set
`SPLIV_SIM_OUTPUT_DIR` to a project-relative output directory to keep generated
artifacts out of version control.
