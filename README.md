# SPLIV reproducibility repository

This repository is the companion to the [`spliv`](https://github.com/koren6684/spliv) R package. It
contains synthetic simulations, the small reference tables/figures needed to
check the paper workflows, and data-gated runners for the Koren (2018) and
Lelkes, Sood, and Iyengar (2017) applications. The package repository contains
the reusable estimators; this repository contains computational workflows and
replication specifications.

## Quick environment and smoke test

From the repository root, the standard public workflow is:

```bash
Rscript -e 'renv::restore(prompt = FALSE)'
Rscript scripts/00_check_environment.R
Rscript scripts/01_run_smoke.R
```

The smoke test uses a tiny synthetic panel, runs baseline IV, a patterned UCI
path, and confirmatory BPE diagnostics, and writes only ignored generated
outputs under `simulations/`.

`renv::restore()` installs `spliv` into the reproducibility environment; a
sibling `spliv` source checkout is not required. The equivalent direct install
after the package tag exists is:

```r
remotes::install_github("koren6684/spliv@v0.1.0")
```

`SPLIV_PACKAGE_PATH` is an optional development override. If set, the scripts
require `devtools` and load that source checkout; when unset, they load the
installed package from the active renv library.

```bash
SPLIV_PACKAGE_PATH=/absolute/path/to/spliv Rscript scripts/01_run_smoke.R
```

The project is pinned to `spliv` version 0.1.0 and the package versions
recorded in `renv.lock`.

## Pilot and full simulations

```bash
Rscript scripts/02_run_simulation_pilot.R
Rscript scripts/03_build_simulation_paper_outputs.R
```

The pilot is intentionally quick (10 replicates, one core, 8-by-8 grid, eight
periods, seed `20260424`). The full production workflow is explicitly long:

```bash
Rscript scripts/99_run_full_simulations.R
```

It uses 1,000 replicates, a 20-by-20 grid, 50 periods, and a delta step of
0.01; runtime and storage depend on hardware and may be measured in hours or
days. Set `SPLIV_SIM_CORES` explicitly when running production simulations.
Checkpoint files, logs, and generated tables/figures are ignored and are not
part of the public repository.

## Koren (2018) data and replication

The Koren runner targets the original crop-yield/conflict panel. Obtain the
authorized file from the original study archive and set:

```bash
export SPLIV_KOREN_DATA=/absolute/path/to/crop.dat.af.rnr3.dta
Rscript replications/koren2018/check_data.R
Rscript scripts/10_run_koren.R
```

The staged code requires `acled_inc_sum`, `maize_yield`, `wheat_yield`, `spi6`,
`gid`, `year`, and `sparsebare`; see
[`replications/koren2018/DATA.md`](replications/koren2018/DATA.md) for the
expected 72,169 complete rows/6,680 clusters and checksum notes. Missing data
produce an acquisition instruction, not a generic file-not-found error.

## Lelkes et al (2017) data and replication

Obtain the authorized archive for the Lelkes, Sood, and Iyengar application and
set the directory containing both required RData files:

```bash
export SPLIV_LELKES_DATA=/absolute/path/to/lelkes-data-directory
Rscript replications/lelkesetal/check_data.R
Rscript scripts/11_run_lelkes.R
```

The required files are `mergeddataset.RData` (object `merged`) and
`mergedcounty.RData` (object `countyproviders`). See
[`replications/lelkesetal/DATA.md`](replications/lelkesetal/DATA.md) for the
required variables, expected 114,803 complete Table 1 rows/48 state clusters,
and checksum notes.

Run all configured replications with:

```bash
Rscript scripts/20_run_all_available_replications.R
```

This dispatcher skips unconfigured data with a clear message and never
substitutes another dataset or specification.

## Recreating tables and figures

Simulation paper outputs are regenerated with
`scripts/03_build_simulation_paper_outputs.R` after the full summary CSVs exist.
The repository includes only the small final reference tables and publication
figures under [`results/reference`](results/reference). Generated outputs must
be compared with those references; they are never overwritten silently.
Replication expected outputs are kept per application under
`replications/*/expected_outputs/`.

## Data and citation

Third-party empirical data are not redistributed because the source archive and
provider terms have not been confirmed for public redistribution. The exact
study citation, DOI, and acquisition URLs in each `DATA.md` require author
confirmation before publication. Cite the `spliv` package and this repository;
see [`CITATION.cff`](CITATION.cff).
