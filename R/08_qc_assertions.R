make_pipeline_check <- function(
    check_id,
    domain,
    severity,
    passed,
    observed,
    expected,
    issue_message) {
  status <- dplyr::case_when(
    isTRUE(passed) ~ "PASS",
    severity == "ERROR" ~ "FAIL",
    TRUE ~ "WARN"
  )

  tibble::tibble(
    check_id = check_id,
    domain = domain,
    severity = severity,
    status = status,
    observed = as.character(observed),
    expected = as.character(expected),
    message = dplyr::if_else(
      status == "PASS",
      "Check passed.",
      issue_message
    )
  )
}

build_pipeline_checks <- function(
    sources,
    cohorts,
    admission_matches,
    score_matches,
    test_summary,
    config) {
  screened <- cohorts$screened
  reviewed <- cohorts$reviewed
  test_count_mismatches <- reconcile_test_record_counts(
    sources$reviewed_patients,
    test_summary$patient_summary
  )

  invalid_percentiles <- sum(
    is.na(reviewed$recorded_review_percentile) |
      reviewed$recorded_review_percentile < 0 |
      reviewed$recorded_review_percentile > 100
  )
  missing_selected_admissions <- sum(is.na(screened$selected_unit))
  negative_ages <- sum(
    !is.na(screened$age_at_admission_days) &
      screened$age_at_admission_days < 0
  )
  missing_rec_results <- sum(is.na(reviewed$mpse_rec_result))
  missing_review_scores <- sum(is.na(score_matches$reconstructed_review_score))
  large_score_discrepancies <- sum(
    abs(score_matches$score_difference) >=
      config$rules$score_discrepancy_warning_threshold,
    na.rm = TRUE
  )
  reviewed_nearest_admissions <- admission_matches$selected |>
    dplyr::filter(
      .data$reviewed,
      .data$admission_match_method == "Nearest arrival date"
    ) |>
    nrow()
  unmapped_result_records <- test_summary$records |>
    dplyr::filter(
      !is.na(.data$diagnosis_by_lab),
      !.data$result_is_classifiable
    ) |>
    nrow()
  tested_without_primary_order_date <- reviewed |>
    dplyr::filter(.data$has_any_test, is.na(.data$primary_order_date)) |>
    nrow()
  official_results_without_date <- reviewed |>
    dplyr::filter(.data$has_official_result_missing_date) |>
    nrow()
  enrollment_date_inconsistencies <- reviewed |>
    dplyr::filter(
      dplyr::coalesce(.data$gk_enrollment_status == "Enrolled", FALSE) !=
        !is.na(.data$gk_enrollment_date_date)
    ) |>
    nrow()
  reviewed_after_cutoff <- if (config$data_mode == "real") {
    sources$reviewed_patients |>
      dplyr::filter(
        (!is.na(.data$mpse_admit_date_date) &
          .data$mpse_admit_date_date > config$dates$admission_cutoff) |
          (!is.na(.data$mpse_review_date_date) &
            .data$mpse_review_date_date > config$dates$last_review_date)
      ) |>
      nrow()
  } else {
    0L
  }
  first_week_review_count <- reviewed |>
    dplyr::filter(
      .data$mpse_review_date_date >= config$dates$study_start,
      .data$mpse_review_date_date <= config$dates$startup_backlog_end
    ) |>
    nrow()
  startup_backlog_count <- sum(reviewed$startup_backlog_patient)
  after_freeze_record_count <- nrow(test_summary$after_freeze_qc)
  reviewed_not_in_screened <- sum(
    !sources$reviewed_patients$mrn %in% sources$mpse_patients$patient_id
  )

  dplyr::bind_rows(
    make_pipeline_check(
      "screened_cohort_count",
      "Cohort",
      "ERROR",
      nrow(screened) == config$rules$expected_screened_patients,
      nrow(screened),
      config$rules$expected_screened_patients,
      "The screened cohort count differs from the expected study cohort."
    ),
    make_pipeline_check(
      "reviewed_cohort_count",
      "Cohort",
      "ERROR",
      nrow(reviewed) == config$rules$expected_reviewed_patients,
      nrow(reviewed),
      config$rules$expected_reviewed_patients,
      "The reviewed cohort count differs from the expected study cohort."
    ),
    make_pipeline_check(
      "unique_screened_patient_ids",
      "Identifiers",
      "ERROR",
      dplyr::n_distinct(sources$mpse_patients$patient_id) ==
        nrow(sources$mpse_patients),
      dplyr::n_distinct(sources$mpse_patients$patient_id),
      nrow(sources$mpse_patients),
      "MPSE patient identifiers are not unique."
    ),
    make_pipeline_check(
      "unique_reviewed_project_ids",
      "Identifiers",
      "ERROR",
      dplyr::n_distinct(sources$reviewed_patients$project_id) ==
        nrow(sources$reviewed_patients),
      dplyr::n_distinct(sources$reviewed_patients$project_id),
      nrow(sources$reviewed_patients),
      "Reviewed project identifiers are not unique."
    ),
    make_pipeline_check(
      "reviewed_patients_in_screened_cohort",
      "Identifiers",
      "ERROR",
      reviewed_not_in_screened == 0L,
      reviewed_not_in_screened,
      0L,
      "One or more reviewed patients are absent from the MPSE patient file."
    ),
    make_pipeline_check(
      "selected_admission_available",
      "Admissions",
      "ERROR",
      missing_selected_admissions == 0L,
      missing_selected_admissions,
      0L,
      "One or more screened patients do not have a selected admission."
    ),
    make_pipeline_check(
      "nonnegative_age_at_admission",
      "Demographics",
      "ERROR",
      negative_ages == 0L,
      negative_ages,
      0L,
      "One or more calculated ages at admission are negative."
    ),
    make_pipeline_check(
      "valid_recorded_review_percentile",
      "MPSE scores",
      "ERROR",
      invalid_percentiles == 0L,
      invalid_percentiles,
      0L,
      "Recorded review percentiles are missing or outside 0–100."
    ),
    make_pipeline_check(
      "recommendation_result_complete",
      "Clinical review",
      "ERROR",
      missing_rec_results == 0L,
      missing_rec_results,
      0L,
      "One or more reviewed patients lack a derived recommendation result."
    ),
    make_pipeline_check(
      "enrollment_status_date_consistent",
      "Gene Kids enrollment",
      "ERROR",
      enrollment_date_inconsistencies == 0L,
      enrollment_date_inconsistencies,
      0L,
      paste0(
        "One or more patients have inconsistent GK enrollment status and ",
        "enrollment-date values."
      )
    ),
    make_pipeline_check(
      "reviewed_records_within_study_cutoffs",
      "Cohort",
      "ERROR",
      reviewed_after_cutoff == 0L,
      reviewed_after_cutoff,
      0L,
      "One or more reviewed records occur after the admission/review cutoff."
    ),
    make_pipeline_check(
      "test_record_counts_reconcile",
      "Testing",
      "ERROR",
      nrow(test_count_mismatches) == 0L,
      nrow(test_count_mismatches),
      0L,
      "Recorded and observed test-result row counts do not reconcile."
    ),
    make_pipeline_check(
      "review_score_within_carry_window",
      "MPSE scores",
      "WARNING",
      missing_review_scores == 0L,
      missing_review_scores,
      0L,
      "One or more review-day raw scores could not be reconstructed."
    ),
    make_pipeline_check(
      "review_score_difference_below_warning_threshold",
      "MPSE scores",
      "WARNING",
      large_score_discrepancies == 0L,
      large_score_discrepancies,
      paste0(
        "0 with absolute difference >= ",
        config$rules$score_discrepancy_warning_threshold
      ),
      paste0(
        "One or more recorded versus reconstructed scores differ by at least ",
        config$rules$score_discrepancy_warning_threshold,
        "."
      )
    ),
    make_pipeline_check(
      "reviewed_admission_requires_nearest_date",
      "Admissions",
      "WARNING",
      reviewed_nearest_admissions == 0L,
      reviewed_nearest_admissions,
      0L,
      "One or more reviewed patients required nearest-date admission matching."
    ),
    make_pipeline_check(
      "all_result_values_classifiable",
      "Testing",
      "WARNING",
      unmapped_result_records == 0L,
      unmapped_result_records,
      0L,
      "One or more laboratory-result values are outside the primary hierarchy."
    ),
    make_pipeline_check(
      "tested_patients_have_primary_order_date",
      "Testing",
      "WARNING",
      tested_without_primary_order_date == 0L,
      tested_without_primary_order_date,
      0L,
      "One or more tested patients lack a usable primary order date."
    ),
    make_pipeline_check(
      "official_results_have_result_date",
      "Testing",
      "WARNING",
      official_results_without_date == 0L,
      official_results_without_date,
      0L,
      paste0(
        "One or more official laboratory classifications lack a result date; ",
        "they remain eligible for diagnostic analysis but not result-timing ",
        "analysis."
      )
    ),
    make_pipeline_check(
      "test_records_within_data_freeze",
      "Testing",
      "WARNING",
      after_freeze_record_count == 0L,
      after_freeze_record_count,
      0L,
      "One or more testing records occur after the configured data freeze."
    ),
    make_pipeline_check(
      "first_week_review_count_matches_methods",
      "Implementation",
      "WARNING",
      config$data_mode != "real" || first_week_review_count == 25L,
      first_week_review_count,
      25L,
      "The first-week review count does not match the methods description."
    ),
    make_pipeline_check(
      "startup_backlog_count_matches_methods",
      "Implementation",
      "WARNING",
      config$data_mode != "real" || startup_backlog_count == 23L,
      startup_backlog_count,
      23L,
      "The startup-backlog count does not match the methods description."
    ),
    make_pipeline_check(
      "calendar_fields_use_original_dates",
      "Calendar",
      "WARNING",
      config$data_mode == "real",
      config$data_mode,
      "real",
      paste0(
        "Implementation phase, backlog status, and calendar summaries are ",
        "provisional with shifted dates."
      )
    )
  )
}

run_pipeline_checks <- function(
    sources,
    cohorts,
    admission_matches,
    score_matches,
    test_summary,
    config,
    stop_on_failure = TRUE) {
  checks <- build_pipeline_checks(
    sources,
    cohorts,
    admission_matches,
    score_matches,
    test_summary,
    config
  )

  readr::write_csv(
    checks,
    file.path(config$output_dirs$qc, "pipeline_checks.csv"),
    na = ""
  )

  failure_count <- sum(checks$status == "FAIL")
  warning_count <- sum(checks$status == "WARN")

  if (warning_count > 0L) {
    message(
      warning_count,
      " QC warning(s); see output/qc/pipeline_checks.csv."
    )
  }

  if (stop_on_failure && failure_count > 0L) {
    stop(
      failure_count,
      " required QC check(s) failed; see output/qc/pipeline_checks.csv.",
      call. = FALSE
    )
  }

  checks
}
