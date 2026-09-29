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

prepare_test_records <- function(test_results, config) {
  normalized_wgs <- tolower(config$rules$wgs_tests)
  classifiable_results <- names(config$rules$result_hierarchy)

  test_results |>
    dplyr::mutate(
      genetic_test = blank_to_na(.data$genetic_test),
      diagnosis_by_lab = blank_to_na(.data$diagnosis_by_lab),
      normalized_test = tolower(stringr::str_squish(.data$genetic_test)),
      is_wgs = !is.na(.data$normalized_test) &
        .data$normalized_test %in% normalized_wgs,
      result_is_classifiable = .data$diagnosis_by_lab %in% classifiable_results,
      result_completion_status = dplyr::case_when(
        !is.na(.data$result_date_date) & !is.na(.data$diagnosis_by_lab) ~
          "Result date and classification recorded",
        !is.na(.data$result_date_date) ~ "Result date only",
        !is.na(.data$diagnosis_by_lab) ~ "Classification only",
        TRUE ~ "Pending or not recorded"
      ),
      completed_and_classifiable = !is.na(.data$result_date_date) &
        .data$result_is_classifiable
    )
}

summarize_test_results <- function(test_results, config) {
  prepared <- prepare_test_records(test_results, config)
  hierarchy <- config$rules$result_hierarchy

  patient_summary <- prepared |>
    dplyr::group_by(.data$project_id) |>
    dplyr::summarise(
      test_record_count = dplyr::n(),
      has_any_test = TRUE,
      has_wgs = any(.data$is_wgs, na.rm = TRUE),
      has_non_wgs = any(!.data$is_wgs & !is.na(.data$genetic_test)),
      has_test_with_missing_name = any(is.na(.data$genetic_test)),
      primary_test_type = dplyr::if_else(
        any(.data$is_wgs, na.rm = TRUE), "WGS", "Non-WGS"
      ),
      primary_order_date = if (any(.data$is_wgs, na.rm = TRUE)) {
        min_date_or_na(.data$order_submitted_date[.data$is_wgs])
      } else {
        min_date_or_na(.data$order_submitted_date)
      },
      has_completed_result_record = any(
        !is.na(.data$result_date_date) & !is.na(.data$diagnosis_by_lab)
      ),
      has_classifiable_completed_result = any(
        .data$completed_and_classifiable
      ),
      has_completed_wgs_result_record = any(
        .data$is_wgs &
          !is.na(.data$result_date_date) &
          !is.na(.data$diagnosis_by_lab)
      ),
      has_classifiable_completed_wgs_result = any(
        .data$is_wgs & .data$completed_and_classifiable
      ),
      has_unmapped_result = any(
        !is.na(.data$diagnosis_by_lab) & !.data$result_is_classifiable
      ),
      has_unmapped_wgs_result = any(
        !is.na(.data$diagnosis_by_lab) &
          !.data$result_is_classifiable &
          .data$is_wgs
      ),
      patient_result = highest_classifiable_result(
        .data$diagnosis_by_lab,
        .data$completed_and_classifiable,
        hierarchy
      ),
      wgs_patient_result = highest_classifiable_result(
        .data$diagnosis_by_lab,
        .data$completed_and_classifiable & .data$is_wgs,
        hierarchy
      ),
      .groups = "drop"
    )

  result_value_qc <- prepared |>
    dplyr::filter(
      !is.na(.data$diagnosis_by_lab),
      !.data$result_is_classifiable
    ) |>
    dplyr::count(.data$diagnosis_by_lab, name = "test_record_count") |>
    dplyr::arrange(dplyr::desc(.data$test_record_count))

  test_name_qc <- prepared |>
    dplyr::filter(is.na(.data$genetic_test)) |>
    dplyr::count(.data$genetic_test, name = "test_record_count")

  list(
    records = prepared,
    patient_summary = patient_summary,
    result_value_qc = result_value_qc,
    test_name_qc = test_name_qc
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
