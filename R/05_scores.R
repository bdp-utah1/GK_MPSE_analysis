match_scores_at_review <- function(reviewed_patients, scores, config) {
  review_fields <- reviewed_patients |>
    dplyr::transmute(
      project_id = .data$project_id,
      mrn = .data$mrn,
      review_date = .data$mpse_review_date_date,
      recorded_review_score = .data$mpse_cp_score,
      recorded_review_percentile = .data$mpse_cp_percent
    )

  eligible_scores <- review_fields |>
    dplyr::inner_join(
      scores |>
        dplyr::select(
          patient_id,
          score_date,
          days_seen,
          mpse_score
        ),
      by = c("mrn" = "patient_id"),
      relationship = "many-to-many"
    ) |>
    dplyr::mutate(
      score_age_days = as.integer(.data$review_date - .data$score_date)
    ) |>
    dplyr::filter(
      !is.na(.data$score_age_days),
      .data$score_age_days >= 0L,
      .data$score_age_days <= config$rules$score_carry_forward_days
    )

  selected_scores <- eligible_scores |>
    dplyr::arrange(
      .data$project_id,
      .data$score_age_days,
      dplyr::desc(.data$score_date)
    ) |>
    dplyr::group_by(.data$project_id) |>
    dplyr::slice_head(n = 1L) |>
    dplyr::ungroup() |>
    dplyr::transmute(
      project_id = .data$project_id,
      reconstructed_score_date = .data$score_date,
      reconstructed_days_seen = .data$days_seen,
      reconstructed_review_score = .data$mpse_score,
      score_age_days = .data$score_age_days
    )

  review_fields |>
    dplyr::left_join(selected_scores, by = "project_id") |>
    dplyr::mutate(
      score_match_status = dplyr::case_when(
        is.na(.data$reconstructed_review_score) ~ "No score within carry window",
        .data$score_age_days == 0L ~ "Score computed on review date",
        TRUE ~ paste0("Carried forward ", .data$score_age_days, " day(s)")
      ),
      score_difference = .data$reconstructed_review_score -
        .data$recorded_review_score,
      score_discrepancy = dplyr::case_when(
        is.na(.data$reconstructed_review_score) |
          is.na(.data$recorded_review_score) ~ NA,
        TRUE ~ !dplyr::near(
          .data$reconstructed_review_score,
          .data$recorded_review_score
        )
      ),
      review_percentile_group = dplyr::case_when(
        is.na(.data$recorded_review_percentile) ~ "Missing",
        .data$recorded_review_percentile < config$rules$percentile_cutoff ~
          "<50th percentile",
        TRUE ~ ">=50th percentile"
      )
    )
}

make_score_qc <- function(score_matches) {
  score_matches |>
    dplyr::filter(
      is.na(.data$reconstructed_review_score) |
        dplyr::coalesce(.data$score_discrepancy, FALSE)
    ) |>
    dplyr::select(
      project_id,
      recorded_review_score,
      reconstructed_review_score,
      score_difference,
      score_age_days,
      score_match_status,
      score_discrepancy
    )
}
