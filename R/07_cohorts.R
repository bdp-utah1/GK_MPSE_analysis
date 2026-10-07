timing_category <- function(interval_days) {
  dplyr::case_when(
    is.na(interval_days) ~ "Event/date unavailable",
    interval_days < 0L ~ "Pre-existing event",
    interval_days == 0L ~ "Same-day event",
    TRUE ~ "Event after reference date"
  )
}

documented_event_timing <- function(event_date, review_date, config) {
  dplyr::case_when(
    is.na(event_date) ~ "Unknown/not documented",
    config$data_mode == "real" & event_date > config$dates$data_freeze ~
      "After data freeze",
    is.na(review_date) ~ "Unknown review date",
    event_date < review_date ~ "Pre-review",
    TRUE ~ "Post-review"
  )
}

event_date_for_analysis <- function(event_date, config) {
  if (config$data_mode == "real") {
    event_date[event_date > config$dates$data_freeze] <- as.Date(NA)
  }
  event_date
}

build_analysis_cohorts <- function(
    sources,
    admission_matches,
    score_matches,
    test_summary,
    config) {
  reviewed_source <- sources$reviewed_patients
  if (config$data_mode == "real") {
    reviewed_source <- reviewed_source |>
      dplyr::filter(
        is.na(.data$mpse_admit_date_date) |
          .data$mpse_admit_date_date <= config$dates$admission_cutoff,
        is.na(.data$mpse_review_date_date) |
          .data$mpse_review_date_date <= config$dates$last_review_date
      )
  }

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

  reviewed_lookup <- reviewed_source |>
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
      date_shifted_calendar_fields_provisional = config$data_mode != "real"
    )

  if (config$data_mode == "real") {
    screened <- screened |>
      dplyr::filter(
        is.na(.data$analysis_admit_date) |
          .data$analysis_admit_date <= config$dates$admission_cutoff
      )
  }

  mpse_demographics <- sources$mpse_patients |>
    dplyr::transmute(
      mrn = .data$patient_id,
      mpse_sex = .data$sex,
      mpse_dob_date = .data$dob_date,
      first_seen_date = .data$first_seen_date,
      last_seen_date = .data$last_seen_date
    )

  reviewed <- reviewed_source |>
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
      any_pre_review_test = dplyr::coalesce(
        .data$any_pre_review_test, FALSE
      ),
      any_post_review_test = dplyr::coalesce(
        .data$any_post_review_test, FALSE
      ),
      any_unknown_timing_test = dplyr::coalesce(
        .data$any_unknown_timing_test, FALSE
      ),
      wgs_pre_review = dplyr::coalesce(.data$wgs_pre_review, FALSE),
      wgs_post_review = dplyr::coalesce(.data$wgs_post_review, FALSE),
      wgs_unknown_timing = dplyr::coalesce(
        .data$wgs_unknown_timing, FALSE
      ),
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
      has_official_result_missing_date = dplyr::coalesce(
        .data$has_official_result_missing_date, FALSE
      ),
      has_official_wgs_result_missing_date = dplyr::coalesce(
        .data$has_official_wgs_result_missing_date, FALSE
      ),
      gk_enrolled = .data$gk_enrollment_status == "Enrolled" &
        !is.na(.data$gk_enrollment_date_date) &
        !(
          config$data_mode == "real" &
            .data$gk_enrollment_date_date > config$dates$data_freeze
        ),
      gk_enrolled = dplyr::coalesce(.data$gk_enrolled, FALSE),
      gk_enrolled_before_review = .data$gk_enrolled &
        .data$gk_enrollment_date_date < .data$mpse_review_date_date,
      eligible_for_post_review_gk_enrollment =
        !.data$gk_enrolled_before_review,
      post_review_gk_enrollment =
        .data$eligible_for_post_review_gk_enrollment &
          .data$gk_enrolled &
          .data$gk_enrollment_date_date >= .data$mpse_review_date_date,
      eligible_for_post_review_any_test = !.data$any_pre_review_test,
      post_review_any_test = dplyr::case_when(
        !.data$eligible_for_post_review_any_test ~ NA,
        .data$any_post_review_test ~ TRUE,
        .data$any_unknown_timing_test ~ NA,
        TRUE ~ FALSE
      ),
      eligible_for_post_review_wgs = !.data$wgs_pre_review,
      post_review_wgs = dplyr::case_when(
        !.data$eligible_for_post_review_wgs ~ NA,
        .data$wgs_post_review ~ TRUE,
        .data$wgs_unknown_timing ~ NA,
        TRUE ~ FALSE
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
      startup_backlog_patient =
        .data$mpse_admit_date_date < config$dates$study_start &
          .data$mpse_review_date_date >= config$dates$study_start &
          .data$mpse_review_date_date <= config$dates$startup_backlog_end,
      prestudy_admission_reviewed_later =
        .data$mpse_admit_date_date < config$dates$study_start &
          .data$mpse_review_date_date > config$dates$startup_backlog_end,
      date_shifted_calendar_fields_provisional = config$data_mode != "real",
      days_admission_to_review = as.integer(
        .data$mpse_review_date_date - .data$mpse_admit_date_date
      ),
      days_admission_to_referral = as.integer(
        event_date_for_analysis(
          .data$gk_referral_received_date_date,
          config
        ) - .data$mpse_admit_date_date
      ),
      days_review_to_referral = as.integer(
        event_date_for_analysis(
          .data$gk_referral_received_date_date,
          config
        ) - .data$mpse_review_date_date
      ),
      days_admission_to_genetics_consult = as.integer(
        event_date_for_analysis(
          .data$mpse_genetic_consult_date_date,
          config
        ) - .data$mpse_admit_date_date
      ),
      days_review_to_genetics_consult = as.integer(
        event_date_for_analysis(
          .data$mpse_genetic_consult_date_date,
          config
        ) - .data$mpse_review_date_date
      ),
      days_admission_to_test_order = as.integer(
        .data$primary_order_date - .data$mpse_admit_date_date
      ),
      days_review_to_test_order = as.integer(
        .data$primary_order_date - .data$mpse_review_date_date
      ),
      days_admission_to_post_review_any_order = as.integer(
        .data$first_post_review_any_order_date - .data$mpse_admit_date_date
      ),
      days_review_to_post_review_any_order = as.integer(
        .data$first_post_review_any_order_date - .data$mpse_review_date_date
      ),
      days_admission_to_post_review_wgs_order = as.integer(
        .data$first_post_review_wgs_order_date - .data$mpse_admit_date_date
      ),
      days_review_to_post_review_wgs_order = as.integer(
        .data$first_post_review_wgs_order_date - .data$mpse_review_date_date
      ),
      genetics_consult_timing = documented_event_timing(
        .data$mpse_genetic_consult_date_date,
        .data$mpse_review_date_date,
        config
      ),
      gk_referral_timing = documented_event_timing(
        .data$gk_referral_received_date_date,
        .data$mpse_review_date_date,
        config
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
      startup_backlog_patient,
      prestudy_admission_reviewed_later,
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
      gk_enrolled_before_review,
      eligible_for_post_review_gk_enrollment,
      post_review_gk_enrollment,
      has_any_test,
      has_wgs,
      has_non_wgs,
      any_pre_review_test,
      any_post_review_test,
      any_unknown_timing_test,
      wgs_pre_review,
      wgs_post_review,
      wgs_unknown_timing,
      eligible_for_post_review_any_test,
      post_review_any_test,
      eligible_for_post_review_wgs,
      post_review_wgs,
      primary_test_type,
      has_completed_result_record,
      has_classifiable_completed_result,
      has_completed_wgs_result_record,
      has_classifiable_completed_wgs_result,
      has_unmapped_result,
      has_unmapped_wgs_result,
      has_official_result_missing_date,
      has_official_wgs_result_missing_date,
      patient_result,
      wgs_patient_result,
      days_admission_to_review,
      days_admission_to_referral,
      days_review_to_referral,
      days_admission_to_genetics_consult,
      days_review_to_genetics_consult,
      days_admission_to_test_order,
      days_review_to_test_order,
      days_admission_to_post_review_any_order,
      days_review_to_post_review_any_order,
      days_admission_to_post_review_wgs_order,
      days_review_to_post_review_wgs_order,
      genetics_consult_timing,
      gk_referral_timing,
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
