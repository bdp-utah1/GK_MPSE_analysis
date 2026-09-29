analysis_config <- function(data_mode = NULL, project_root = getwd()) {
  if (is.null(data_mode)) {
    data_mode <- getOption(
      "gk_mpse.data_mode",
      Sys.getenv("GK_MPSE_DATA_MODE", unset = "synthetic")
    )
  }

  data_mode <- tolower(data_mode)
  if (!data_mode %in% c("synthetic", "real")) {
    stop("data_mode must be either 'synthetic' or 'real'.", call. = FALSE)
  }

  project_root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)

  file_sets <- list(
    synthetic = list(
      mpse_patients = "mpse_data/synthetic_mpse_patients.tsv",
      mpse_admissions = "mpse_data/synthetic_mpse_admissions.tsv",
      mpse_scores = "mpse_data/synthetic_mpse_scores.tsv",
      reviewed_patients = "mosaic_data/synthetic_patient_data.csv",
      test_results = "mosaic_data/synthetic_test_results_data.csv"
    ),
    real = list(
      mpse_patients = "mpse_data/patient.tsv",
      mpse_admissions = "mpse_data/patient_stay.tsv",
      mpse_scores = "mpse_data/patient_mpse_score.tsv",
      reviewed_patients = "mosaic_data/patient_level_data.csv",
      test_results = "mosaic_data/test_results_data.csv"
    )
  )

  relative_files <- file_sets[[data_mode]]
  input_files <- lapply(
    relative_files,
    function(path) file.path(project_root, path)
  )

  list(
    project_root = project_root,
    data_mode = data_mode,
    input_files = input_files,
    output_dirs = list(
      data = file.path(project_root, "output", "data"),
      qc = file.path(project_root, "output", "qc"),
      tables = file.path(project_root, "output", "tables"),
      figures = file.path(project_root, "output", "figures")
    ),
    dates = list(
      study_start = as.Date("2025-08-11"),
      study_end = as.Date("2026-06-02"),
      implementation_change = as.Date("2025-10-30"),
      data_freeze = as.Date("2026-09-01")
    ),
    rules = list(
      score_carry_forward_days = 3L,
      score_discrepancy_warning_threshold = 1.0,
      hpo_source = "ClinPhen_bulk",
      percentile_cutoff = 50,
      regression_confidence_level = 0.95,
      percentile_regression_scale = 10,
      raw_score_regression_scale = 10,
      expected_screened_patients = 314L,
      expected_reviewed_patients = 254L,
      other_recommendation_result_project_ids = c("8146"),
      wgs_tests = c(
        "WGS", "WGS reanalysis", "rWGS", "rWGS Reanalysis", "urWGS"
      ),
      result_hierarchy = c(
        "Negative" = 1L,
        "Uncertain" = 2L,
        "Positive" = 3L
      ),
      unit_priority = c(
        "CPH 4E NEONATAL ICU" = 1L,
        "PC_4E NICU" = 2L,
        "CPL NEONATAL ICU" = 3L,
        "PL_NICU" = 4L
      )
    )
  )
}
