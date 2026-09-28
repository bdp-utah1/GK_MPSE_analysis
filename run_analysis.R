suppressPackageStartupMessages({
  library(tidyverse)
  library(fs)
  library(gt)
})

source(file.path("R", "00_config.R"))
source(file.path("R", "01_utils.R"))
source(file.path("R", "02_import.R"))
source(file.path("R", "03_validate.R"))
source(file.path("R", "04_admissions.R"))
source(file.path("R", "05_scores.R"))
source(file.path("R", "06_testing.R"))
source(file.path("R", "07_cohorts.R"))

arguments <- commandArgs(trailingOnly = TRUE)
mode_argument <- grep("^--mode=", arguments, value = TRUE)
data_mode <- if (length(mode_argument) > 0L) {
  sub("^--mode=", "", mode_argument[[1]])
} else {
  NULL
}

config <- analysis_config(data_mode = data_mode)
assert_project_root(config$project_root)
make_output_dirs(config)

message("Importing ", config$data_mode, " data...")
sources <- import_analysis_sources(config)
qc <- run_import_qc(sources, config, write_files = TRUE)

admission_targets <- make_admission_targets(sources)
admission_matches <- match_index_admissions(
  admission_targets,
  sources$mpse_admissions,
  config
)
score_matches <- match_scores_at_review(
  sources$reviewed_patients,
  sources$mpse_scores,
  config
)
test_summary <- summarize_test_results(sources$test_results, config)
cohorts <- build_analysis_cohorts(
  sources,
  admission_matches,
  score_matches,
  test_summary,
  config
)

derived_qc <- list(
  admission_matches = make_admission_qc(admission_matches),
  score_discrepancies = make_score_qc(score_matches),
  unmapped_result_values = test_summary$result_value_qc,
  missing_test_names = test_summary$test_name_qc,
  test_record_count_mismatches = reconcile_test_record_counts(
    sources$reviewed_patients,
    test_summary$patient_summary
  )
)
purrr::iwalk(derived_qc, function(data, name) {
  readr::write_csv(
    data,
    file.path(config$output_dirs$qc, paste0(name, ".csv")),
    na = ""
  )
})

write_cohort_outputs(cohorts, config)

analysis_state <- list(
  config = config,
  sources = sources,
  qc = c(qc, derived_qc),
  admission_matches = admission_matches,
  score_matches = score_matches,
  test_summary = test_summary,
  cohorts = cohorts
)

message("Import, derivation, and validation complete.")
message("QC files: output/qc/")
message("De-identified analytic extracts: output/data/")

invisible(analysis_state)
