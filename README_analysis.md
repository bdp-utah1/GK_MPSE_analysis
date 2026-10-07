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
- `R/11_score_regressions.R`: exploratory score–outcome models, model-population audits, and figures
- `analysis/01_descriptive.R`: descriptive-analysis entry point
- `analysis/02_score_regressions.R`: score-regression entry point
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
`output/data/`, descriptive and regression CSV/HTML tables to `output/tables/`,
and cohort-flow and regression figures to `output/figures/`.

`output/data/admission_match_manifest.csv` contains one row per screened patient
with the admission-matching target, selected admission, ranking fields, candidate
counts, and review flag. It intentionally retains patient linkage so the selected
admission can be joined back into patient-level analyses.

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

## Exploratory score regressions

The regression stage treats the recorded percentile at clinical review as the
authoritative percentile. It reports the clinical implementation comparison
(`>=50th` versus `<50th` percentile), a continuous percentile analysis per 10
percentile points, and a parallel reconstructed raw-score analysis per 10 score
points. All models are unadjusted and exploratory.

Binary testing and result outcomes use logistic regression. Testing outcomes are
post-review outcomes: same-day orders count as post-review, the any-testing risk
set excludes patients with any pre-review test, and the WGS risk set excludes only
patients with pre-review WGS. Orders with missing dates remain temporally unknown.
Observed nonnegative consultation, referral, and test-order intervals use linear
regression. Patients without a documented event/date are not assigned a censoring
time, so these are not survival models. Pre-existing events are excluded from
timing regressions and tabulated separately.

A nonmissing laboratory classification is an official result even when the result
date is missing; the missing date only prevents result-turnaround analysis.
Patients with unmatched result values are excluded from the relevant diagnostic-
result models until those values are adjudicated. Models with a zero cell or
insufficient outcome variation are reported as not estimable by ordinary logistic
regression rather than displaying unstable odds ratios.

June 2, 2026 is the inclusive admission and clinical-review cutoff. Follow-up
outcomes through September 1, 2026 are included. Results are reported overall and
by implementation phase. Admission-based timing outcomes also include a
sensitivity population excluding the startup backlog.

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
