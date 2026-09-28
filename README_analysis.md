# MPSE analysis framework

This framework separates raw-data import, validation, cohort derivation, analysis,
and reporting. Raw files in `mpse_data/` and `mosaic_data/` are read but never
modified.

## Current structure

- `R/00_config.R`: file paths, study dates, implementation rules, and category definitions
- `R/01_utils.R`: shared validation and date-parsing helpers
- `R/02_import.R`: typed imports and MRN canonicalization
- `R/03_validate.R`: source counts, identifier overlap, duplicate-key checks, and missingness
- `R/04_admissions.R`: selects one admission record per screened patient
- `R/05_scores.R`: reconstructs raw MPSE score at clinical review
- `R/06_testing.R`: collapses test records to patient-level testing and result summaries
- `R/07_cohorts.R`: creates the screened and reviewed analytic cohorts
- `R/08_qc_assertions.R`: required checks and non-blocking QC warnings
- `R/09_descriptive_tables.R`: cohort, testing, and missingness summaries
- `R/10_cohort_flow.R`: cohort-flow table and PNG figure
- `analysis/01_descriptive.R`: descriptive-analysis entry point
- `run_analysis.R`: entry point
- `output/qc/`: machine-readable validation summaries
- `output/data/`, `output/tables/`, `output/figures/`: reserved for later stages

## Run with de-identified development data

From the `GK_MPSE_analysis` R project directory:

```r
source("run_analysis.R")
```

From a terminal:

```sh
Rscript run_analysis.R --mode=synthetic
```

The command writes QC checks to `output/qc/`, analytic extracts to
`output/data/`, descriptive CSV/HTML tables to `output/tables/`, and the
cohort-flow figure to `output/figures/`.

The complete de-identified data are the default while the framework is under
development. The configuration value remains `synthetic` because the supplied
filenames retain that prefix.

## Run with real data

In RStudio:

```r
options(gk_mpse.data_mode = "real")
source("run_analysis.R")
```

Or from a terminal:

```sh
Rscript run_analysis.R --mode=real
```

Historical score percentiles are not required. The contemporaneous percentile
recorded at clinical review (`mpse_cp_percent`) is the authoritative percentile
for review-day grouping.

## Planned next stages

1. Review the automated checks and descriptive outputs.
2. Resolve or document remaining warning-level exceptions.
3. Add exploratory percentile-group testing and turnaround-time figures.

The de-identified dates preserve within-patient intervals but can distort calendar
phase, backlog, and monthly-volume assignments. Those fields are retained for code
testing and explicitly marked provisional until the original dates are used.

Admission matching ranks candidates by calendar-day distance from the target
admission date, then by the configured unit priority, and finally by whether the
recorded arrival-discharge interval contains the target date. Missing recommendation
results are categorized as `Other` when the corresponding recommendation is
`Other`. Project 8146 is also categorized as `Other` because the patient died
before a testing decision could occur. The original source value is retained in
`mpse_rec_result_raw`.
