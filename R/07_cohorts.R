timing_category <- function(interval_days) {
  dplyr::case_when(
    is.na(interval_days) ~ "Event/date unavailable",
    interval_days < 0L ~ "Pre-existing event",
    interval_days == 0L ~ "Same-day event",
    TRUE ~ "Event after reference date"
  )
}

build_analysis_cohorts <- function(
    sources,
    admission_matches,
    score_matches,
    test_summary,
    config) {
  admission_fields <- admission_matches$selected |>
    dplyr::select(
      patient_id,
      selected_unit,
      selected_arrive_date,
      selected_dischg_date,
      admission_match_method,
      arrival_distance_days,
      admission_candidate_count,
      admission_match_review_flag
    )

  reviewed_lookup <- sources$reviewed_patients |>
    dplyr::select(
      project_id,
      mrn,
      reviewed_dob_date = dob_date,
      proband_gender,
      mpse_admit_date_date
    )

  screened <- sources$mpse_patients |>
    dplyr::transmute(
      patient_id = .data$patient_id,
      mpse_sex = .data$sex,
      mpse_dob_date = .data$dob_date,
      first_seen_date = .data$first_seen_date,
      last_seen_date = .data$last_seen_date
    ) |>
    dplyr::left_join(reviewed_lookup, by = c("patient_id" = "mrn")) |>
    dplyr::left_join(admission_fields, by = "patient_id") |>
    dplyr::mutate(
      reviewed = !is.na(.data$project_id),
      analysis_dob_date = dplyr::coalesce(
        .data$reviewed_dob_date, .data$mpse_dob_date
      ),
      analysis_admit_date = dplyr::coalesce(
        .data$mpse_admit_date_date, .data$selected_arrive_date
      ),
      age_at_admission_days = as.integer(
        .data$analysis_admit_date - .data$analysis_dob_date
      ),
      sex = dplyr::coalesce(.data$proband_gender, .data$mpse_sex),
      date_shifted_calendar_fields_provisional = TRUE
    )

  mpse_demographics <- sources$mpse_patients |>
    dplyr::transmute(
      mrn = .data$patient_id,
      mpse_sex = .data$sex,
      mpse_dob_date = .data$dob_date,
      first_seen_date = .data$first_seen_date,
      last_seen_date = .data$last_seen_date
    )

  reviewed <- sources$reviewed_patients |>
    dplyr::left_join(mpse_demographics, by = "mrn") |>
    dplyr::left_join(admission_fields, by = c("mrn" = "patient_id")) |>
    dplyr::left_join(
      score_matches |>
        dplyr::select(-mrn, -review_date),
      by = "project_id"
    ) |>
    dplyr::left_join(test_summary$patient_summary, by = "project_id") |>
    dplyr::mutate(
      has_any_test = dplyr::coalesce(.data$has_any_test, FALSE),
      has_wgs = dplyr::coalesce(.data$has_wgs, FALSE),
      has_non_wgs = dplyr::coalesce(.data$has_non_wgs, FALSE),
      has_test_with_missing_name = dplyr::coalesce(
        .data$has_test_with_missing_name, FALSE
      ),
      test_record_count = dplyr::coalesce(.data$test_record_count, 0L),
      primary_test_type = dplyr::coalesce(
        .data$primary_test_type, "No testing"
      ),
      has_completed_result_record = dplyr::coalesce(
        .data$has_completed_result_record, FALSE
      ),
      has_classifiable_completed_result = dplyr::coalesce(
        .data$has_classifiable_completed_result, FALSE
      ),
      has_completed_wgs_result_record = dplyr::coalesce(
        .data$has_completed_wgs_result_record, FALSE
      ),
      has_classifiable_completed_wgs_result = dplyr::coalesce(
        .data$has_classifiable_completed_wgs_result, FALSE
      ),
      has_unmapped_result = dplyr::coalesce(
        .data$has_unmapped_result, FALSE
      ),
      has_unmapped_wgs_result = dplyr::coalesce(
        .data$has_unmapped_wgs_result, FALSE
      ),
      gk_enrolled = dplyr::coalesce(
        .data$gk_enrollment_status == "Enrolled", FALSE
      ),
      age_at_admission_days = as.integer(
        .data$mpse_admit_date_date - .data$dob_date
      ),
      implementation_phase = dplyr::case_when(
        .data$mpse_review_date_date < config$dates$implementation_change ~
          "Initial threshold-based phase",
        TRUE ~ "Expanded review phase"
      ),
      backlog_patient = .data$mpse_admit_date_date < config$dates$study_start,
      date_shifted_calendar_fields_provisional = TRUE,
      days_admission_to_review = as.integer(
        .data$mpse_review_date_date - .data$mpse_admit_date_date
      ),
      days_admission_to_referral = as.integer(
        .data$gk_referral_received_date_date - .data$mpse_admit_date_date
      ),
      days_review_to_referral = as.integer(
        .data$gk_referral_received_date_date - .data$mpse_review_date_date
      ),
      days_admission_to_test_order = as.integer(
        .data$primary_order_date - .data$mpse_admit_date_date
      ),
      days_review_to_test_order = as.integer(
        .data$primary_order_date - .data$mpse_review_date_date
      ),
      referral_timing_from_review = timing_category(
        .data$days_review_to_referral
      ),
      test_timing_from_review = timing_category(
        .data$days_review_to_test_order
      ),
      demographic_dob_discrepancy = !is.na(.data$mpse_dob_date) &
        !is.na(.data$dob_date) & .data$mpse_dob_date != .data$dob_date,
      demographic_sex_discrepancy = !is.na(.data$mpse_sex) &
        !is.na(.data$proband_gender) & .data$mpse_sex != .data$proband_gender
    )

  list(screened = screened, reviewed = reviewed)
}

write_cohort_outputs <- function(cohorts, config) {
  reviewed_output <- cohorts$reviewed |>
    dplyr::select(
      project_id,
      proband_gender,
      age_at_admission_days,
      selected_unit,
      admission_match_method,
      arrival_distance_days,
      implementation_phase,
      backlog_patient,
      recorded_review_score,
      reconstructed_review_score,
      score_difference,
      score_discrepancy,
      score_age_days,
      recorded_review_percentile,
      review_percentile_group,
      mpse_rec,
      mpse_rec_result,
      gk_enrollment_status,
      gk_enrollment_type,
      gk_enrolled,
      has_any_test,
      has_wgs,
      has_non_wgs,
      primary_test_type,
      has_completed_result_record,
      has_classifiable_completed_result,
      has_completed_wgs_result_record,
      has_classifiable_completed_wgs_result,
      has_unmapped_result,
      has_unmapped_wgs_result,
      patient_result,
      wgs_patient_result,
      days_admission_to_review,
      days_admission_to_referral,
      days_review_to_referral,
      days_admission_to_test_order,
      days_review_to_test_order,
      referral_timing_from_review,
      test_timing_from_review,
      date_shifted_calendar_fields_provisional
    )

  readr::write_csv(
    reviewed_output,
    file.path(config$output_dirs$data, "reviewed_cohort_analysis.csv"),
    na = ""
  )

  screened_output <- cohorts$screened |>
    dplyr::transmute(
      screened_record = dplyr::row_number(),
      project_id = .data$project_id,
      reviewed = .data$reviewed,
      sex = .data$sex,
      age_at_admission_days = .data$age_at_admission_days,
      selected_unit = .data$selected_unit,
      admission_match_method = .data$admission_match_method,
      arrival_distance_days = .data$arrival_distance_days,
      date_shifted_calendar_fields_provisional =
        .data$date_shifted_calendar_fields_provisional
    )

  readr::write_csv(
    screened_output,
    file.path(config$output_dirs$data, "screened_cohort_analysis.csv"),
    na = ""
  )
}
