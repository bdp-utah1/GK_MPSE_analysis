unit_priority_value <- function(unit, config) {
  priority <- unname(config$rules$unit_priority[unit])
  fallback <- length(config$rules$unit_priority) + 1L
  dplyr::coalesce(as.integer(priority), fallback)
}

make_admission_targets <- function(sources) {
  reviewed_targets <- sources$reviewed_patients |>
    dplyr::transmute(
      patient_id = .data$mrn,
      project_id = .data$project_id,
      reviewed = TRUE,
      target_date = .data$mpse_admit_date_date,
      target_date_source = "recorded MPSE admission date"
    )

  sources$mpse_patients |>
    dplyr::transmute(
      patient_id = .data$patient_id,
      first_seen_date = .data$first_seen_date
    ) |>
    dplyr::left_join(reviewed_targets, by = "patient_id") |>
    dplyr::mutate(
      reviewed = dplyr::coalesce(.data$reviewed, FALSE),
      target_date = dplyr::coalesce(.data$target_date, .data$first_seen_date),
      target_date_source = dplyr::coalesce(
        .data$target_date_source,
        "first MPSE-seen date"
      )
    ) |>
    dplyr::select(-first_seen_date)
}

match_index_admissions <- function(targets, admissions, config) {
  candidates <- targets |>
    dplyr::left_join(admissions, by = "patient_id", relationship = "many-to-many") |>
    dplyr::mutate(
      unit_priority = unit_priority_value(.data$unit, config),
      same_start_date = !is.na(.data$arrive_date) &
        .data$arrive_date == .data$target_date,
      recorded_interval_overlap = !is.na(.data$arrive_date) &
        !is.na(.data$dischg_date) &
        .data$arrive_date <= .data$target_date &
        .data$dischg_date >= .data$target_date,
      arrival_distance_days = abs(as.integer(.data$arrive_date - .data$target_date))
    )

  selected <- candidates |>
    dplyr::group_by(.data$patient_id) |>
    dplyr::mutate(
      admission_candidate_count = sum(!is.na(.data$arrive_date)),
      recorded_interval_overlap_count = sum(
        .data$recorded_interval_overlap,
        na.rm = TRUE
      ),
      nearest_arrival_distance_days = if (all(is.na(.data$arrival_distance_days))) {
        NA_integer_
      } else {
        min(.data$arrival_distance_days, na.rm = TRUE)
      },
      nearest_arrival_candidate_count = sum(
        .data$arrival_distance_days == .data$nearest_arrival_distance_days,
        na.rm = TRUE
      )
    ) |>
    dplyr::arrange(
      .data$arrival_distance_days,
      .data$unit_priority,
      dplyr::desc(.data$recorded_interval_overlap),
      .data$arrive_date,
      .by_group = TRUE
    ) |>
    dplyr::slice_head(n = 1L) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      admission_match_method = dplyr::case_when(
        is.na(.data$target_date) ~ "Missing target date",
        is.na(.data$arrive_date) ~ "No admission record",
        .data$same_start_date ~ "Same start date",
        .data$recorded_interval_overlap ~ "Recorded interval overlap",
        TRUE ~ "Nearest arrival date"
      ),
      admission_match_review_flag = .data$admission_match_method %in% c(
        "Missing target date", "No admission record", "Nearest arrival date"
      )
    ) |>
    dplyr::rename(
      selected_unit = unit,
      selected_arrive_dtm = arrive_dtm,
      selected_dischg_dtm = dischg_dtm,
      selected_arrive_date = arrive_date,
      selected_dischg_date = dischg_date
    )

  list(selected = selected, candidates = candidates)
}

make_admission_qc <- function(admission_matches) {
  admission_matches$selected |>
    dplyr::filter(
      .data$admission_match_review_flag |
        .data$admission_candidate_count > 1L
    ) |>
    dplyr::transmute(
      project_id = .data$project_id,
      reviewed = .data$reviewed,
      target_date_source = .data$target_date_source,
      selected_unit = .data$selected_unit,
      admission_match_method = .data$admission_match_method,
      arrival_distance_days = .data$arrival_distance_days,
      admission_candidate_count = .data$admission_candidate_count,
      nearest_arrival_candidate_count = .data$nearest_arrival_candidate_count,
      recorded_interval_overlap_count =
        .data$recorded_interval_overlap_count,
      review_flag = .data$admission_match_review_flag
    )
}
