min_date_or_na <- function(x) {
  if (length(x) == 0L || all(is.na(x))) {
    return(as.Date(NA))
  }
  min(x, na.rm = TRUE)
}

highest_classifiable_result <- function(result, eligible, hierarchy) {
  values <- result[eligible & result %in% names(hierarchy)]
  if (length(values) == 0L) {
    return(NA_character_)
  }
  ranks <- unname(hierarchy[values])
  values[[which.max(ranks)]]
}

date_within_outcome_window <- function(date, config) {
  if (config$data_mode == "real") {
    is.na(date) | date <= config$dates$data_freeze
  } else {
    rep(TRUE, length(date))
  }
}

prepare_test_records <- function(test_results, reviewed_patients, config) {
  normalized_wgs <- tolower(config$rules$wgs_tests)
  classifiable_results <- names(config$rules$result_hierarchy)

  test_results |>
    dplyr::left_join(
      reviewed_patients |>
        dplyr::select(
          project_id,
          review_date = mpse_review_date_date
        ),
      by = "project_id"
    ) |>
    dplyr::mutate(
      genetic_test = blank_to_na(.data$genetic_test),
      diagnosis_by_lab = blank_to_na(.data$diagnosis_by_lab),
      normalized_test = tolower(stringr::str_squish(.data$genetic_test)),
      is_wgs = !is.na(.data$normalized_test) &
        .data$normalized_test %in% normalized_wgs,
      result_is_classifiable = .data$diagnosis_by_lab %in% classifiable_results,
      test_in_analysis_window = date_within_outcome_window(
        .data$order_submitted_date,
        config
      ),
      result_in_analysis_window = .data$test_in_analysis_window &
        date_within_outcome_window(.data$result_date_date, config),
      official_result_available = !is.na(.data$diagnosis_by_lab) &
        .data$result_in_analysis_window,
      classifiable_result_available = .data$result_is_classifiable &
        .data$result_in_analysis_window,
      order_timing_from_review = dplyr::case_when(
        !.data$test_in_analysis_window ~ "After data freeze",
        is.na(.data$order_submitted_date) ~ "Unknown order date",
        is.na(.data$review_date) ~ "Unknown review date",
        .data$order_submitted_date < .data$review_date ~ "Pre-review",
        TRUE ~ "Post-review"
      ),
      result_completion_status = dplyr::case_when(
        .data$official_result_available & !is.na(.data$result_date_date) ~
          "Result date and classification recorded",
        .data$official_result_available ~ "Classification only",
        .data$result_in_analysis_window & !is.na(.data$result_date_date) ~
          "Result date only",
        !.data$result_in_analysis_window ~ "After data freeze",
        TRUE ~ "Pending or not recorded"
      ),
      completed_and_classifiable = .data$classifiable_result_available
    )
}

summarize_test_results <- function(test_results, reviewed_patients, config) {
  prepared <- prepare_test_records(test_results, reviewed_patients, config)
  hierarchy <- config$rules$result_hierarchy

  patient_summary <- prepared |>
    dplyr::group_by(.data$project_id) |>
    dplyr::summarise(
      test_record_count = dplyr::n(),
      analysis_test_record_count = sum(.data$test_in_analysis_window),
      has_any_test = any(.data$test_in_analysis_window),
      has_wgs = any(.data$is_wgs & .data$test_in_analysis_window),
      has_non_wgs = any(
        !.data$is_wgs &
          !is.na(.data$genetic_test) &
          .data$test_in_analysis_window
      ),
      has_test_with_missing_name = any(is.na(.data$genetic_test)),
      primary_test_type = dplyr::if_else(
        any(.data$is_wgs & .data$test_in_analysis_window),
        "WGS",
        "Non-WGS"
      ),
      first_any_order_date = min_date_or_na(
        .data$order_submitted_date[.data$test_in_analysis_window]
      ),
      first_wgs_order_date = min_date_or_na(
        .data$order_submitted_date[
          .data$is_wgs & .data$test_in_analysis_window
        ]
      ),
      first_post_review_any_order_date = min_date_or_na(
        .data$order_submitted_date[
          .data$order_timing_from_review == "Post-review"
        ]
      ),
      first_post_review_wgs_order_date = min_date_or_na(
        .data$order_submitted_date[
          .data$is_wgs & .data$order_timing_from_review == "Post-review"
        ]
      ),
      primary_order_date = if (
        any(.data$is_wgs & .data$test_in_analysis_window)
      ) {
        min_date_or_na(
          .data$order_submitted_date[
            .data$is_wgs & .data$test_in_analysis_window
          ]
        )
      } else {
        min_date_or_na(
          .data$order_submitted_date[.data$test_in_analysis_window]
        )
      },
      any_pre_review_test = any(
        .data$order_timing_from_review == "Pre-review"
      ),
      any_post_review_test = any(
        .data$order_timing_from_review == "Post-review"
      ),
      any_unknown_timing_test = any(
        .data$order_timing_from_review %in%
          c("Unknown order date", "Unknown review date")
      ),
      wgs_pre_review = any(
        .data$is_wgs & .data$order_timing_from_review == "Pre-review"
      ),
      wgs_post_review = any(
        .data$is_wgs & .data$order_timing_from_review == "Post-review"
      ),
      wgs_unknown_timing = any(
        .data$is_wgs &
          .data$order_timing_from_review %in%
            c("Unknown order date", "Unknown review date")
      ),
      has_completed_result_record = any(
        .data$official_result_available
      ),
      has_classifiable_completed_result = any(
        .data$classifiable_result_available
      ),
      has_completed_wgs_result_record = any(
        .data$is_wgs & .data$official_result_available
      ),
      has_classifiable_completed_wgs_result = any(
        .data$is_wgs & .data$classifiable_result_available
      ),
      has_unmapped_result = any(
        .data$official_result_available & !.data$result_is_classifiable
      ),
      has_unmapped_wgs_result = any(
        .data$official_result_available &
          !.data$result_is_classifiable &
          .data$is_wgs
      ),
      has_official_result_missing_date = any(
        .data$official_result_available & is.na(.data$result_date_date)
      ),
      has_official_wgs_result_missing_date = any(
        .data$is_wgs &
          .data$official_result_available &
          is.na(.data$result_date_date)
      ),
      patient_result = highest_classifiable_result(
        .data$diagnosis_by_lab,
        .data$classifiable_result_available,
        hierarchy
      ),
      wgs_patient_result = highest_classifiable_result(
        .data$diagnosis_by_lab,
        .data$classifiable_result_available & .data$is_wgs,
        hierarchy
      ),
      .groups = "drop"
    )

  result_value_qc <- prepared |>
    dplyr::filter(
      .data$official_result_available,
      !.data$result_is_classifiable
    ) |>
    dplyr::count(.data$diagnosis_by_lab, name = "test_record_count") |>
    dplyr::arrange(dplyr::desc(.data$test_record_count))

  test_name_qc <- prepared |>
    dplyr::filter(.data$test_in_analysis_window, is.na(.data$genetic_test)) |>
    dplyr::count(.data$genetic_test, name = "test_record_count")

  after_freeze_qc <- prepared |>
    dplyr::filter(
      !.data$test_in_analysis_window |
        (!is.na(.data$result_date_date) &
          config$data_mode == "real" &
          .data$result_date_date > config$dates$data_freeze)
    ) |>
    dplyr::select(
      project_id,
      genetic_test,
      order_submitted_date,
      result_date_date,
      diagnosis_by_lab,
      test_in_analysis_window,
      result_in_analysis_window
    )

  list(
    records = prepared,
    patient_summary = patient_summary,
    result_value_qc = result_value_qc,
    test_name_qc = test_name_qc,
    after_freeze_qc = after_freeze_qc
  )
}

reconcile_test_record_counts <- function(reviewed_patients, patient_summary) {
  reviewed_patients |>
    dplyr::transmute(
      project_id = .data$project_id,
      recorded_test_count = suppressWarnings(
        as.integer(.data$test_results_record_count)
      )
    ) |>
    dplyr::left_join(
      patient_summary |>
        dplyr::select(project_id, observed_test_count = test_record_count),
      by = "project_id"
    ) |>
    dplyr::mutate(
      recorded_test_count = dplyr::coalesce(.data$recorded_test_count, 0L),
      observed_test_count = dplyr::coalesce(.data$observed_test_count, 0L),
      count_matches = .data$recorded_test_count == .data$observed_test_count
    ) |>
    dplyr::filter(!.data$count_matches)
}
