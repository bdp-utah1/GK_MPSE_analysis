format_n_percent <- function(n, denominator) {
  if (is.na(denominator) || denominator == 0L) {
    return(as.character(n))
  }
  sprintf("%d/%d (%.1f%%)", n, denominator, 100 * n / denominator)
}

format_median_iqr <- function(median, q1, q3) {
  if (is.na(median)) {
    return("Not available")
  }
  sprintf("%.1f [%.1f, %.1f]", median, q1, q3)
}

make_count_row <- function(
    section,
    metric,
    numerator,
    denominator = NA_integer_,
    denominator_label = NA_character_) {
  tibble::tibble(
    section = section,
    metric = metric,
    numerator = as.integer(numerator),
    denominator = as.integer(denominator),
    denominator_label = denominator_label,
    percent = dplyr::if_else(
      !is.na(denominator) & denominator > 0L,
      100 * numerator / denominator,
      NA_real_
    ),
    display = format_n_percent(numerator, denominator)
  )
}

make_cohort_accounting_table <- function(cohorts) {
  screened <- cohorts$screened
  reviewed <- cohorts$reviewed

  n_screened <- nrow(screened)
  n_reviewed <- nrow(reviewed)
  n_tested <- sum(reviewed$has_any_test)
  n_wgs <- sum(reviewed$has_wgs)
  n_gk <- sum(reviewed$gk_enrolled)
  n_completed <- sum(reviewed$has_completed_result_record)
  adjudicated_result_population <- reviewed |>
    dplyr::filter(
      .data$has_classifiable_completed_result,
      !.data$has_unmapped_result
    )
  n_classifiable <- nrow(adjudicated_result_population)
  n_positive <- sum(
    adjudicated_result_population$patient_result == "Positive",
    na.rm = TRUE
  )
  n_uncertain <- sum(
    adjudicated_result_population$patient_result == "Uncertain",
    na.rm = TRUE
  )
  n_negative <- sum(
    adjudicated_result_population$patient_result == "Negative",
    na.rm = TRUE
  )

  dplyr::bind_rows(
    make_count_row("Cohort", "Screened by MPSE", n_screened),
    make_count_row(
      "Cohort", "Clinician reviewed", n_reviewed, n_screened, "Screened"
    ),
    make_count_row(
      "Testing", "Any genetic testing", n_tested, n_reviewed, "Reviewed"
    ),
    make_count_row(
      "Testing", "WGS", n_wgs, n_reviewed, "Reviewed"
    ),
    make_count_row(
      "Testing", "WGS among tested patients", n_wgs, n_tested, "Tested"
    ),
    make_count_row(
      "Gene Kids", "GK enrolled", n_gk, n_reviewed, "Reviewed"
    ),
    make_count_row(
      "Results",
      "Official result classification",
      n_completed,
      n_tested,
      "Tested"
    ),
    make_count_row(
      "Results",
      "Classifiable result without unmatched values",
      n_classifiable,
      n_tested,
      "Tested"
    ),
    make_count_row(
      "Diagnostic result",
      "Positive",
      n_positive,
      n_classifiable,
      "Classifiable results without unmatched values"
    ),
    make_count_row(
      "Diagnostic result",
      "Uncertain",
      n_uncertain,
      n_classifiable,
      "Classifiable results without unmatched values"
    ),
    make_count_row(
      "Diagnostic result",
      "Negative",
      n_negative,
      n_classifiable,
      "Classifiable results without unmatched values"
    )
  )
}

expand_summary_groups <- function(data, group_variable, overall_label = "Overall") {
  grouped <- data |>
    dplyr::mutate(summary_group = as.character(.data[[group_variable]]))

  overall <- data |>
    dplyr::mutate(summary_group = overall_label)

  dplyr::bind_rows(overall, grouped)
}

summarize_continuous_by_group <- function(
    data,
    variable,
    characteristic,
    units = NA_character_) {
  data |>
    dplyr::group_by(.data$summary_group) |>
    dplyr::summarise(
      characteristic = characteristic,
      level = NA_character_,
      units = units,
      n = sum(!is.na(.data[[variable]])),
      denominator = dplyr::n(),
      percent = NA_real_,
      median = stats::median(.data[[variable]], na.rm = TRUE),
      q1 = stats::quantile(.data[[variable]], 0.25, na.rm = TRUE, names = FALSE),
      q3 = stats::quantile(.data[[variable]], 0.75, na.rm = TRUE, names = FALSE),
      minimum = min(.data[[variable]], na.rm = TRUE),
      maximum = max(.data[[variable]], na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(c("median", "q1", "q3", "minimum", "maximum")),
        ~ dplyr::if_else(is.infinite(.x), NA_real_, as.numeric(.x))
      ),
      display = purrr::pmap_chr(
        list(.data$median, .data$q1, .data$q3),
        format_median_iqr
      )
    )
}

summarize_categorical_by_group <- function(
    data,
    variable,
    characteristic,
    missing_label = "Not recorded") {
  data |>
    dplyr::mutate(
      summary_level = dplyr::coalesce(
        as.character(.data[[variable]]),
        missing_label
      )
    ) |>
    dplyr::count(
      .data$summary_group,
      .data$summary_level,
      name = "n"
    ) |>
    dplyr::group_by(.data$summary_group) |>
    dplyr::mutate(
      denominator = sum(.data$n),
      percent = 100 * .data$n / .data$denominator
    ) |>
    dplyr::ungroup() |>
    dplyr::transmute(
      summary_group = .data$summary_group,
      characteristic = characteristic,
      level = .data$summary_level,
      units = NA_character_,
      n = .data$n,
      denominator = .data$denominator,
      percent = .data$percent,
      median = NA_real_,
      q1 = NA_real_,
      q3 = NA_real_,
      minimum = NA_real_,
      maximum = NA_real_,
      display = sprintf(
        "%d/%d (%.1f%%)",
        .data$n,
        .data$denominator,
        .data$percent
      )
    )
}

make_screened_characteristics_table <- function(cohorts) {
  data <- cohorts$screened |>
    dplyr::mutate(
      review_group = dplyr::if_else(
        .data$reviewed,
        "Reviewed",
        "Not reviewed"
      )
    ) |>
    expand_summary_groups("review_group")

  dplyr::bind_rows(
    summarize_continuous_by_group(
      data,
      "age_at_admission_days",
      "Age at admission",
      "days"
    ),
    summarize_categorical_by_group(data, "sex", "Gender"),
    summarize_categorical_by_group(
      data,
      "selected_unit",
      "Selected admission unit"
    ),
    summarize_categorical_by_group(
      data,
      "admission_match_method",
      "Admission matching method"
    )
  ) |>
    dplyr::arrange(
      factor(
        .data$summary_group,
        levels = c("Overall", "Reviewed", "Not reviewed")
      ),
      .data$characteristic,
      .data$level
    )
}

make_reviewed_characteristics_table <- function(cohorts) {
  data <- cohorts$reviewed |>
    dplyr::mutate(
      genetics_consult_recorded = dplyr::if_else(
        !is.na(.data$mpse_genetic_consult_date_date),
        "Genetics consultation recorded",
        "No genetics consultation date"
      ),
      gk_status_summary = dplyr::coalesce(
        .data$gk_enrollment_status,
        "No documented GK enrollment status"
      ),
      testing_summary = dplyr::case_when(
        .data$has_wgs ~ "WGS",
        .data$has_any_test ~ "Non-WGS testing only",
        TRUE ~ "Untested"
      ),
      result_summary = dplyr::case_when(
        .data$has_unmapped_result ~ "Unmatched result value",
        !is.na(.data$patient_result) ~ .data$patient_result,
        TRUE ~ "Untested"
      )
    ) |>
    expand_summary_groups("review_percentile_group")

  dplyr::bind_rows(
    summarize_continuous_by_group(
      data,
      "age_at_admission_days",
      "Age at admission",
      "days"
    ),
    summarize_continuous_by_group(
      data,
      "reconstructed_review_score",
      "Reconstructed raw MPSE score at review",
      "score"
    ),
    summarize_continuous_by_group(
      data,
      "recorded_review_percentile",
      "Recorded MPSE percentile at review",
      "percentile"
    ),
    summarize_continuous_by_group(
      data,
      "days_admission_to_review",
      "Admission to clinical review",
      "days"
    ),
    summarize_categorical_by_group(data, "proband_gender", "Gender"),
    summarize_categorical_by_group(
      data,
      "mpse_rec",
      "MPSE recommendation"
    ),
    summarize_categorical_by_group(
      data,
      "mpse_rec_result",
      "MPSE recommendation result"
    ),
    summarize_categorical_by_group(
      data,
      "genetics_consult_recorded",
      "Genetics consultation"
    ),
    summarize_categorical_by_group(
      data,
      "gk_status_summary",
      "GK enrollment status"
    ),
    summarize_categorical_by_group(
      data,
      "testing_summary",
      "Testing"
    ),
    summarize_categorical_by_group(
      data,
      "result_summary",
      "Patient-level result"
    ),
    summarize_categorical_by_group(
      data,
      "implementation_phase",
      "Implementation phase"
    ),
    summarize_categorical_by_group(
      data,
      "genetics_consult_timing",
      "Documented Genetics consultation timing"
    ),
    summarize_categorical_by_group(
      data,
      "gk_referral_timing",
      "Documented GK referral timing"
    )
  ) |>
    dplyr::arrange(
      factor(
        .data$summary_group,
        levels = c(
          "Overall",
          "<50th percentile",
          ">=50th percentile",
          "Missing"
        )
      ),
      .data$characteristic,
      .data$level
    )
}

make_testing_diagnostic_table <- function(cohorts) {
  reviewed <- cohorts$reviewed
  n_reviewed <- nrow(reviewed)
  n_tested <- sum(reviewed$has_any_test)
  diagnostic_population <- reviewed |>
    dplyr::filter(
      .data$has_classifiable_completed_result,
      !.data$has_unmapped_result
    )
  n_classifiable <- nrow(diagnostic_population)
  any_test_at_risk <- reviewed |>
    dplyr::filter(.data$eligible_for_post_review_any_test)
  any_test_known <- any_test_at_risk |>
    dplyr::filter(!is.na(.data$post_review_any_test))
  wgs_at_risk <- reviewed |>
    dplyr::filter(.data$eligible_for_post_review_wgs)
  wgs_known <- wgs_at_risk |>
    dplyr::filter(!is.na(.data$post_review_wgs))
  gk_at_risk <- reviewed |>
    dplyr::filter(.data$eligible_for_post_review_gk_enrollment)

  dplyr::bind_rows(
    make_count_row(
      "Testing uptake", "Any genetic testing", n_tested, n_reviewed, "Reviewed"
    ),
    make_count_row(
      "Testing uptake", "WGS", sum(reviewed$has_wgs), n_reviewed, "Reviewed"
    ),
    make_count_row(
      "Testing uptake",
      "Any non-WGS testing",
      sum(reviewed$has_non_wgs),
      n_reviewed,
      "Reviewed"
    ),
    make_count_row(
      "Post-review uptake",
      "Any testing after review",
      sum(any_test_known$post_review_any_test),
      nrow(any_test_known),
      "No pre-review testing and known order timing"
    ),
    make_count_row(
      "Post-review uptake",
      "Any-test timing unknown",
      sum(is.na(any_test_at_risk$post_review_any_test)),
      nrow(any_test_at_risk),
      "No pre-review testing"
    ),
    make_count_row(
      "Post-review uptake",
      "WGS after review",
      sum(wgs_known$post_review_wgs),
      nrow(wgs_known),
      "No pre-review WGS and known WGS order timing"
    ),
    make_count_row(
      "Post-review uptake",
      "WGS timing unknown",
      sum(is.na(wgs_at_risk$post_review_wgs)),
      nrow(wgs_at_risk),
      "No pre-review WGS"
    ),
    make_count_row(
      "Post-review uptake",
      "GK enrollment after review",
      sum(gk_at_risk$post_review_gk_enrollment),
      nrow(gk_at_risk),
      "Not enrolled before review"
    ),
    make_count_row(
      "Testing uptake",
      "GK enrolled",
      sum(reviewed$gk_enrolled),
      n_reviewed,
      "Reviewed"
    ),
    make_count_row(
      "Result availability",
      "Official result classification",
      sum(reviewed$has_completed_result_record),
      n_tested,
      "Tested"
    ),
    make_count_row(
      "Result availability",
      "Classifiable result without unmatched values",
      n_classifiable,
      n_tested,
      "Tested"
    ),
    make_count_row(
      "Result availability",
      "Unmatched result value",
      sum(reviewed$has_unmapped_result),
      n_tested,
      "Tested"
    ),
    make_count_row(
      "Diagnostic result",
      "Positive",
      sum(diagnostic_population$patient_result == "Positive", na.rm = TRUE),
      n_classifiable,
      "Classifiable results without unmatched values"
    ),
    make_count_row(
      "Diagnostic result",
      "Uncertain",
      sum(diagnostic_population$patient_result == "Uncertain", na.rm = TRUE),
      n_classifiable,
      "Classifiable results without unmatched values"
    ),
    make_count_row(
      "Diagnostic result",
      "Negative",
      sum(diagnostic_population$patient_result == "Negative", na.rm = TRUE),
      n_classifiable,
      "Classifiable results without unmatched values"
    )
  )
}

make_missingness_table <- function(cohorts) {
  reviewed <- cohorts$reviewed
  tested <- reviewed |>
    dplyr::filter(.data$has_any_test)
  official_results <- reviewed |>
    dplyr::filter(.data$has_completed_result_record)
  status_enrolled <- reviewed |>
    dplyr::filter(.data$gk_enrollment_status == "Enrolled")

  missingness_row <- function(data, variable, label, population) {
    n_missing <- sum(is.na(data[[variable]]))
    denominator <- nrow(data)
    tibble::tibble(
      population = population,
      variable = variable,
      label = label,
      n_missing = n_missing,
      denominator = denominator,
      percent_missing = dplyr::if_else(
        denominator > 0L,
        100 * n_missing / denominator,
        NA_real_
      ),
      display = format_n_percent(n_missing, denominator)
    )
  }

  missingness_count_row <- function(
      variable,
      label,
      population,
      n_missing,
      denominator) {
    tibble::tibble(
      population = population,
      variable = variable,
      label = label,
      n_missing = as.integer(n_missing),
      denominator = as.integer(denominator),
      percent_missing = dplyr::if_else(
        denominator > 0L,
        100 * n_missing / denominator,
        NA_real_
      ),
      display = format_n_percent(n_missing, denominator)
    )
  }

  dplyr::bind_rows(
    missingness_row(
      reviewed,
      "reconstructed_review_score",
      "Reconstructed review-day raw score",
      "Reviewed"
    ),
    missingness_row(
      reviewed,
      "recorded_review_percentile",
      "Recorded review percentile",
      "Reviewed"
    ),
    missingness_row(
      status_enrolled,
      "gk_enrollment_date_date",
      "Enrollment date among records marked Enrolled",
      "Status marked Enrolled"
    ),
    missingness_row(
      tested,
      "primary_order_date",
      "Primary test order date",
      "Tested"
    ),
    missingness_count_row(
      "result_date_date",
      "Patients with at least one official result missing its result date",
      "Official result available",
      sum(official_results$has_official_result_missing_date),
      nrow(official_results)
    )
  )
}

make_event_documentation_table <- function(cohorts) {
  reviewed <- cohorts$reviewed |>
    dplyr::mutate(
      enrollment_timing = dplyr::case_when(
        .data$gk_enrolled_before_review ~ "Pre-review",
        .data$post_review_gk_enrollment ~ "Post-review",
        TRUE ~ "Not enrolled"
      )
    )

  event_rows <- list(
    genetics_consult_timing = "Documented Genetics consultation",
    gk_referral_timing = "Documented GK referral",
    enrollment_timing = "GK enrollment"
  )

  purrr::imap_dfr(event_rows, function(label, variable) {
    reviewed |>
      dplyr::count(.data[[variable]], name = "n") |>
      dplyr::mutate(
        event = label,
        status = .data[[variable]],
        denominator = nrow(reviewed),
        percent = 100 * .data$n / .data$denominator,
        .before = 1
      ) |>
      dplyr::select(-dplyr::all_of(variable))
  })
}

make_phase_backlog_table <- function(cohorts) {
  cohorts$reviewed |>
    dplyr::mutate(
      backlog_category = dplyr::case_when(
        .data$startup_backlog_patient ~ "Startup backlog (first review week)",
        .data$prestudy_admission_reviewed_later ~
          "Pre-study admission reviewed later",
        TRUE ~ "Admitted on or after study start"
      )
    ) |>
    dplyr::count(
      .data$implementation_phase,
      .data$backlog_category,
      .data$review_percentile_group,
      name = "n"
    )
}

build_descriptive_tables <- function(cohorts) {
  list(
    cohort_accounting = make_cohort_accounting_table(cohorts),
    screened_characteristics = make_screened_characteristics_table(cohorts),
    reviewed_characteristics = make_reviewed_characteristics_table(cohorts),
    testing_diagnostic_summary = make_testing_diagnostic_table(cohorts),
    key_field_missingness = make_missingness_table(cohorts),
    event_documentation = make_event_documentation_table(cohorts),
    phase_backlog_summary = make_phase_backlog_table(cohorts)
  )
}

write_descriptive_tables <- function(tables, config) {
  purrr::iwalk(tables, function(data, name) {
    csv_path <- file.path(config$output_dirs$tables, paste0(name, ".csv"))
    html_path <- file.path(config$output_dirs$tables, paste0(name, ".html"))

    readr::write_csv(data, csv_path, na = "")

    gt_table <- gt::gt(data) |>
      gt::opt_row_striping() |>
      gt::tab_options(table.font.size = gt::px(12))
    gt::gtsave(gt_table, html_path)
  })
}
