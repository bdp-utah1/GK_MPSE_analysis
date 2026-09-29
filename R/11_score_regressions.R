prepare_score_regression_data <- function(reviewed, config) {
  reviewed |>
    dplyr::mutate(
      score_above_50 = as.integer(
        .data$recorded_review_percentile >= config$rules$percentile_cutoff
      ),
      percentile_per_10 = .data$recorded_review_percentile /
        config$rules$percentile_regression_scale,
      raw_score_per_10 = .data$reconstructed_review_score /
        config$rules$raw_score_regression_scale
    )
}

make_model_spec <- function(
    model_id,
    outcome_label,
    population_label,
    model_type,
    data,
    source_population_n,
    unavailable_outcome_n = 0L,
    preexisting_excluded_n = 0L,
    unmatched_result_excluded_n = 0L) {
  list(
    model_id = model_id,
    outcome_label = outcome_label,
    population_label = population_label,
    model_type = model_type,
    data = data,
    source_population_n = as.integer(source_population_n),
    unavailable_outcome_n = as.integer(unavailable_outcome_n),
    preexisting_excluded_n = as.integer(preexisting_excluded_n),
    unmatched_result_excluded_n = as.integer(
      unmatched_result_excluded_n
    )
  )
}

make_binary_spec <- function(
    data,
    model_id,
    outcome_label,
    population_label,
    outcome,
    source_population_n = nrow(data),
    unavailable_outcome_n = 0L,
    unmatched_result_excluded_n = 0L) {
  model_data <- data |>
    dplyr::mutate(model_outcome = as.integer(.data[[outcome]])) |>
    dplyr::filter(!is.na(.data$model_outcome))

  make_model_spec(
    model_id = model_id,
    outcome_label = outcome_label,
    population_label = population_label,
    model_type = "Binary",
    data = model_data,
    source_population_n = source_population_n,
    unavailable_outcome_n = unavailable_outcome_n,
    unmatched_result_excluded_n = unmatched_result_excluded_n
  )
}

make_timing_spec <- function(
    data,
    model_id,
    outcome_label,
    population_label,
    interval_variable) {
  observed <- data |>
    dplyr::filter(!is.na(.data[[interval_variable]]))
  model_data <- observed |>
    dplyr::filter(.data[[interval_variable]] >= 0L) |>
    dplyr::mutate(model_outcome = as.numeric(.data[[interval_variable]]))

  make_model_spec(
    model_id = model_id,
    outcome_label = outcome_label,
    population_label = population_label,
    model_type = "Continuous",
    data = model_data,
    source_population_n = nrow(data),
    unavailable_outcome_n = sum(is.na(data[[interval_variable]])),
    preexisting_excluded_n = sum(
      observed[[interval_variable]] < 0L,
      na.rm = TRUE
    )
  )
}

build_score_model_specs <- function(reviewed, config) {
  data <- prepare_score_regression_data(reviewed, config)

  all_completed <- data |>
    dplyr::filter(.data$has_completed_result_record)
  all_diagnostic <- all_completed |>
    dplyr::filter(
      !.data$has_unmapped_result,
      !is.na(.data$patient_result)
    ) |>
    dplyr::mutate(
      positive_result = .data$patient_result == "Positive",
      uncertain_result = .data$patient_result == "Uncertain"
    )

  wgs_completed <- data |>
    dplyr::filter(.data$has_completed_wgs_result_record)
  wgs_diagnostic <- wgs_completed |>
    dplyr::filter(
      !.data$has_unmapped_wgs_result,
      !is.na(.data$wgs_patient_result)
    ) |>
    dplyr::mutate(
      positive_wgs_result = .data$wgs_patient_result == "Positive",
      uncertain_wgs_result = .data$wgs_patient_result == "Uncertain"
    )

  gk_completed <- data |>
    dplyr::filter(
      .data$gk_enrolled,
      .data$has_completed_wgs_result_record
    )
  gk_diagnostic <- gk_completed |>
    dplyr::filter(
      !.data$has_unmapped_wgs_result,
      !is.na(.data$wgs_patient_result)
    ) |>
    dplyr::mutate(
      positive_gk_result = .data$wgs_patient_result == "Positive",
      uncertain_gk_result = .data$wgs_patient_result == "Uncertain"
    )

  all_unavailable <- sum(
    !all_completed$has_unmapped_result & is.na(all_completed$patient_result)
  )
  wgs_unavailable <- sum(
    !wgs_completed$has_unmapped_wgs_result &
      is.na(wgs_completed$wgs_patient_result)
  )
  gk_unavailable <- sum(
    !gk_completed$has_unmapped_wgs_result &
      is.na(gk_completed$wgs_patient_result)
  )

  test_timing_population <- data |>
    dplyr::filter(.data$has_any_test)

  list(
    make_binary_spec(
      data, "any_testing", "Any genetic testing", "All reviewed patients",
      "has_any_test"
    ),
    make_binary_spec(
      data, "wgs_testing", "WGS", "All reviewed patients", "has_wgs"
    ),
    make_binary_spec(
      data, "gk_enrollment", "GK enrollment", "All reviewed patients",
      "gk_enrolled"
    ),
    make_binary_spec(
      all_diagnostic,
      "positive_any_result",
      "Positive diagnostic result (any test)",
      "Completed-result patients without unmatched values",
      "positive_result",
      source_population_n = nrow(all_completed),
      unavailable_outcome_n = all_unavailable,
      unmatched_result_excluded_n = sum(all_completed$has_unmapped_result)
    ),
    make_binary_spec(
      all_diagnostic,
      "uncertain_any_result",
      "Uncertain result (any test)",
      "Completed-result patients without unmatched values",
      "uncertain_result",
      source_population_n = nrow(all_completed),
      unavailable_outcome_n = all_unavailable,
      unmatched_result_excluded_n = sum(all_completed$has_unmapped_result)
    ),
    make_binary_spec(
      wgs_diagnostic,
      "positive_wgs_result",
      "Positive WGS result",
      "Completed-WGS-result patients without unmatched WGS values",
      "positive_wgs_result",
      source_population_n = nrow(wgs_completed),
      unavailable_outcome_n = wgs_unavailable,
      unmatched_result_excluded_n = sum(
        wgs_completed$has_unmapped_wgs_result
      )
    ),
    make_binary_spec(
      wgs_diagnostic,
      "uncertain_wgs_result",
      "Uncertain WGS result",
      "Completed-WGS-result patients without unmatched WGS values",
      "uncertain_wgs_result",
      source_population_n = nrow(wgs_completed),
      unavailable_outcome_n = wgs_unavailable,
      unmatched_result_excluded_n = sum(
        wgs_completed$has_unmapped_wgs_result
      )
    ),
    make_binary_spec(
      gk_diagnostic,
      "positive_gk_result",
      "Positive GK WGS result",
      "GK-enrolled patients with completed WGS results",
      "positive_gk_result",
      source_population_n = nrow(gk_completed),
      unavailable_outcome_n = gk_unavailable,
      unmatched_result_excluded_n = sum(
        gk_completed$has_unmapped_wgs_result
      )
    ),
    make_binary_spec(
      gk_diagnostic,
      "uncertain_gk_result",
      "Uncertain GK WGS result",
      "GK-enrolled patients with completed WGS results",
      "uncertain_gk_result",
      source_population_n = nrow(gk_completed),
      unavailable_outcome_n = gk_unavailable,
      unmatched_result_excluded_n = sum(
        gk_completed$has_unmapped_wgs_result
      )
    ),
    make_timing_spec(
      data,
      "admission_to_referral",
      "Admission to GK referral",
      "All reviewed patients",
      "days_admission_to_referral"
    ),
    make_timing_spec(
      data,
      "review_to_referral",
      "Clinical review to GK referral",
      "All reviewed patients",
      "days_review_to_referral"
    ),
    make_timing_spec(
      test_timing_population,
      "admission_to_test_order",
      "Admission to primary test order",
      "Patients with any testing record",
      "days_admission_to_test_order"
    ),
    make_timing_spec(
      test_timing_population,
      "review_to_test_order",
      "Clinical review to primary test order",
      "Patients with any testing record",
      "days_review_to_test_order"
    )
  )
}

score_predictor_specs <- function() {
  tibble::tribble(
    ~predictor, ~predictor_label,
    "score_above_50", ">=50th versus <50th recorded percentile",
    "percentile_per_10", "Per 10-percentile-point increase",
    "raw_score_per_10", "Per 10-point reconstructed raw-score increase"
  )
}

fit_score_model <- function(spec, predictor_spec, config) {
  predictor <- predictor_spec$predictor[[1]]
  predictor_label <- predictor_spec$predictor_label[[1]]
  data <- spec$data |>
    dplyr::filter(
      !is.na(.data$model_outcome),
      !is.na(.data[[predictor]])
    )
  analysis_n <- nrow(data)
  event_n <- if (spec$model_type == "Binary") {
    sum(data$model_outcome == 1L)
  } else {
    NA_integer_
  }
  non_event_n <- if (spec$model_type == "Binary") {
    sum(data$model_outcome == 0L)
  } else {
    NA_integer_
  }

  warning_messages <- character()
  estimable <- analysis_n > 1L &&
    dplyr::n_distinct(data[[predictor]]) > 1L &&
    dplyr::n_distinct(data$model_outcome) > 1L

  if (spec$model_type == "Binary" && predictor == "score_above_50") {
    cells <- table(
      factor(data[[predictor]], levels = c(0, 1)),
      factor(data$model_outcome, levels = c(0, 1))
    )
    if (any(cells == 0L)) {
      estimable <- FALSE
      warning_messages <- c(
        warning_messages,
        paste0(
          "Ordinary logistic regression not estimated because a 2x2 cell ",
          "is zero (complete or quasi-complete separation)."
        )
      )
    } else if (any(cells < 5L)) {
      warning_messages <- c(
        warning_messages,
        "Sparse outcome cells; interpret the odds ratio cautiously."
      )
    }
  }

  if (!estimable && length(warning_messages) == 0L) {
    warning_messages <- c(
      warning_messages,
      "Outcome or predictor variation is insufficient for estimation."
    )
  }

  estimate <- std_error <- conf_low <- conf_high <- p_value <- NA_real_
  model_converged <- NA
  estimate_type <- if (spec$model_type == "Binary") {
    "Odds ratio"
  } else {
    "Mean difference in days"
  }

  if (estimable) {
    formula <- stats::reformulate(predictor, response = "model_outcome")
    fit_warnings <- character()
    fit <- withCallingHandlers(
      if (spec$model_type == "Binary") {
        stats::glm(formula, data = data, family = stats::binomial())
      } else {
        stats::lm(formula, data = data)
      },
      warning = function(warning) {
        fit_warnings <<- c(fit_warnings, conditionMessage(warning))
        invokeRestart("muffleWarning")
      }
    )
    coefficients <- summary(fit)$coefficients
    if (nrow(coefficients) >= 2L) {
      beta <- coefficients[2, 1]
      std_error <- coefficients[2, 2]
      p_value <- coefficients[2, ncol(coefficients)]
      critical_value <- stats::qnorm(
        1 - (1 - config$rules$regression_confidence_level) / 2
      )
      beta_low <- beta - critical_value * std_error
      beta_high <- beta + critical_value * std_error
      if (spec$model_type == "Binary") {
        estimate <- exp(beta)
        conf_low <- exp(beta_low)
        conf_high <- exp(beta_high)
        model_converged <- isTRUE(fit$converged)
        if (!model_converged || !all(is.finite(c(estimate, conf_low, conf_high)))) {
          warning_messages <- c(
            warning_messages,
            "Logistic model did not produce a stable finite estimate."
          )
        }
      } else {
        estimate <- beta
        conf_low <- beta_low
        conf_high <- beta_high
        model_converged <- TRUE
      }
    }
    warning_messages <- c(warning_messages, fit_warnings)
  }

  tibble::tibble(
    model_id = spec$model_id,
    outcome = spec$outcome_label,
    population = spec$population_label,
    model_type = spec$model_type,
    predictor = predictor,
    predictor_label = predictor_label,
    analysis_n = analysis_n,
    event_n = event_n,
    non_event_n = non_event_n,
    estimate_type = estimate_type,
    estimate = estimate,
    std_error = std_error,
    conf_low = conf_low,
    conf_high = conf_high,
    p_value = p_value,
    model_converged = model_converged,
    status = dplyr::if_else(estimable, "Estimated", "Not estimable"),
    model_note = dplyr::if_else(
      length(warning_messages) == 0L,
      "",
      paste(unique(warning_messages), collapse = " ")
    )
  )
}

make_score_regression_results <- function(specs, config) {
  predictors <- score_predictor_specs()
  purrr::map_dfr(specs, function(spec) {
    purrr::map_dfr(seq_len(nrow(predictors)), function(index) {
      fit_score_model(spec, predictors[index, ], config)
    })
  }) |>
    dplyr::mutate(
      estimate_display = dplyr::if_else(
        is.na(.data$estimate),
        "Not estimable",
        sprintf(
          "%.2f (%.2f, %.2f)",
          .data$estimate,
          .data$conf_low,
          .data$conf_high
        )
      )
    )
}

make_score_model_population_audit <- function(specs) {
  purrr::map_dfr(specs, function(spec) {
    tibble::tibble(
      model_id = spec$model_id,
      outcome = spec$outcome_label,
      source_population = spec$population_label,
      source_population_n = spec$source_population_n,
      outcome_or_date_unavailable_n = spec$unavailable_outcome_n,
      preexisting_event_excluded_n = spec$preexisting_excluded_n,
      unmatched_result_excluded_n = spec$unmatched_result_excluded_n,
      analysis_n = nrow(spec$data)
    )
  })
}

summarize_score_group <- function(data, spec, group_label) {
  if (spec$model_type == "Binary") {
    events <- sum(data$model_outcome == 1L)
    denominator <- nrow(data)
    tibble::tibble(
      model_id = spec$model_id,
      outcome = spec$outcome_label,
      model_type = spec$model_type,
      percentile_group = group_label,
      n = denominator,
      events = events,
      percent = dplyr::if_else(
        denominator > 0L,
        100 * events / denominator,
        NA_real_
      ),
      mean = NA_real_,
      standard_deviation = NA_real_,
      median = NA_real_,
      q1 = NA_real_,
      q3 = NA_real_
    )
  } else {
    values <- data$model_outcome
    tibble::tibble(
      model_id = spec$model_id,
      outcome = spec$outcome_label,
      model_type = spec$model_type,
      percentile_group = group_label,
      n = length(values),
      events = NA_integer_,
      percent = NA_real_,
      mean = mean(values),
      standard_deviation = stats::sd(values),
      median = stats::median(values),
      q1 = stats::quantile(values, 0.25, names = FALSE),
      q3 = stats::quantile(values, 0.75, names = FALSE)
    )
  }
}

make_score_group_descriptives <- function(specs) {
  purrr::map_dfr(specs, function(spec) {
    overall <- summarize_score_group(spec$data, spec, "Overall")
    grouped <- spec$data |>
      dplyr::group_split(.data$review_percentile_group, .keep = TRUE)
    group_rows <- purrr::map_dfr(grouped, function(group_data) {
      summarize_score_group(
        group_data,
        spec,
        as.character(group_data$review_percentile_group[[1]])
      )
    })
    dplyr::bind_rows(overall, group_rows)
  }) |>
    dplyr::mutate(
      display = dplyr::case_when(
        .data$model_type == "Binary" ~ sprintf(
          "%d/%d (%.1f%%)",
          .data$events,
          .data$n,
          .data$percent
        ),
        TRUE ~ sprintf(
          "%.1f [%.1f, %.1f]",
          .data$median,
          .data$q1,
          .data$q3
        )
      )
    )
}

make_timing_status_table <- function(reviewed) {
  timing_variables <- tibble::tribble(
    ~interval_variable, ~outcome,
    "days_admission_to_referral", "Admission to GK referral",
    "days_review_to_referral", "Clinical review to GK referral",
    "days_admission_to_test_order", "Admission to primary test order",
    "days_review_to_test_order", "Clinical review to primary test order"
  )

  purrr::map_dfr(seq_len(nrow(timing_variables)), function(index) {
    variable <- timing_variables$interval_variable[[index]]
    outcome <- timing_variables$outcome[[index]]
    reviewed |>
      dplyr::mutate(
        event_status = timing_category(.data[[variable]]),
        score_group = as.character(.data$review_percentile_group)
      ) |>
      dplyr::bind_rows(
        reviewed |>
          dplyr::mutate(
            event_status = timing_category(.data[[variable]]),
            score_group = "Overall"
          )
      ) |>
      dplyr::count(.data$score_group, .data$event_status, name = "n") |>
      dplyr::group_by(.data$score_group) |>
      dplyr::mutate(
        denominator = sum(.data$n),
        percent = 100 * .data$n / .data$denominator
      ) |>
      dplyr::ungroup() |>
      dplyr::mutate(
        interval_variable = variable,
        outcome = outcome,
        .before = 1
      )
  })
}

make_score_regression_plots <- function(results) {
  binary_data <- results |>
    dplyr::filter(
      .data$model_type == "Binary",
      .data$status == "Estimated",
      is.finite(.data$estimate),
      is.finite(.data$conf_low),
      is.finite(.data$conf_high)
    )
  timing_data <- results |>
    dplyr::filter(
      .data$model_type == "Continuous",
      .data$status == "Estimated",
      is.finite(.data$estimate),
      is.finite(.data$conf_low),
      is.finite(.data$conf_high)
    )

  binary_plot <- ggplot2::ggplot(
    binary_data,
    ggplot2::aes(
      x = .data$estimate,
      y = .data$outcome,
      xmin = .data$conf_low,
      xmax = .data$conf_high
    )
  ) +
    ggplot2::geom_vline(xintercept = 1, color = "#9AA5B1", linetype = 2) +
    ggplot2::geom_errorbarh(height = 0.18, color = "#4B6B88") +
    ggplot2::geom_point(size = 2.4, color = "#1F4E79") +
    ggplot2::scale_x_log10() +
    ggplot2::facet_wrap(
      ggplot2::vars(.data$predictor_label),
      scales = "free_y",
      ncol = 1
    ) +
    ggplot2::labs(
      title = "Exploratory score associations with binary outcomes",
      subtitle = "Unadjusted logistic regression; points are odds ratios with 95% confidence intervals",
      x = "Odds ratio (log scale)",
      y = NULL,
      caption = "Models with separation or insufficient variation are omitted from the figure and retained in the results table."
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold")
    )

  timing_plot <- ggplot2::ggplot(
    timing_data,
    ggplot2::aes(
      x = .data$estimate,
      y = .data$outcome,
      xmin = .data$conf_low,
      xmax = .data$conf_high
    )
  ) +
    ggplot2::geom_vline(xintercept = 0, color = "#9AA5B1", linetype = 2) +
    ggplot2::geom_errorbarh(height = 0.18, color = "#6A5A8C") +
    ggplot2::geom_point(size = 2.4, color = "#4D3975") +
    ggplot2::facet_wrap(
      ggplot2::vars(.data$predictor_label),
      scales = "free_y",
      ncol = 1
    ) +
    ggplot2::labs(
      title = "Exploratory score associations with observed timing outcomes",
      subtitle = "Unadjusted linear regression; points are mean differences in days with 95% confidence intervals",
      x = "Mean difference (days)",
      y = NULL,
      caption = "Only observed intervals on or after the reference date are modeled; this is not a time-to-event survival analysis."
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold")
    )

  list(binary = binary_plot, timing = timing_plot)
}

make_score_regression_checks <- function(results, population_audit) {
  expected_analysis_n <- population_audit$analysis_n[
    match(results$model_id, population_audit$model_id)
  ]
  population_reconciliation_failures <- population_audit |>
    dplyr::filter(
      .data$source_population_n !=
        .data$outcome_or_date_unavailable_n +
          .data$preexisting_event_excluded_n +
          .data$unmatched_result_excluded_n +
          .data$analysis_n
    ) |>
    nrow()
  oversized_model_populations <- sum(
    results$analysis_n > expected_analysis_n,
    na.rm = TRUE
  )
  nonfinite_estimates <- results |>
    dplyr::filter(
      .data$status == "Estimated",
      !is.finite(.data$estimate) |
        !is.finite(.data$conf_low) |
        !is.finite(.data$conf_high)
    ) |>
    nrow()
  undocumented_unestimated <- results |>
    dplyr::filter(
      .data$status == "Not estimable",
      is.na(.data$model_note) | .data$model_note == ""
    ) |>
    nrow()
  not_estimable_count <- sum(results$status == "Not estimable")

  dplyr::bind_rows(
    make_pipeline_check(
      "score_model_populations_reconcile",
      "Score regressions",
      "ERROR",
      population_reconciliation_failures == 0L,
      population_reconciliation_failures,
      0L,
      "One or more score-model populations do not reconcile."
    ),
    make_pipeline_check(
      "score_model_complete_case_counts_valid",
      "Score regressions",
      "ERROR",
      oversized_model_populations == 0L,
      oversized_model_populations,
      0L,
      "A fitted model contains more observations than its eligible population."
    ),
    make_pipeline_check(
      "score_model_estimates_finite",
      "Score regressions",
      "ERROR",
      nonfinite_estimates == 0L,
      nonfinite_estimates,
      0L,
      "One or more estimated models produced nonfinite estimates."
    ),
    make_pipeline_check(
      "unestimated_score_models_documented",
      "Score regressions",
      "ERROR",
      undocumented_unestimated == 0L,
      undocumented_unestimated,
      0L,
      "One or more unestimated models lack an explanatory note."
    ),
    make_pipeline_check(
      "score_models_not_estimable",
      "Score regressions",
      "WARNING",
      not_estimable_count == 0L,
      not_estimable_count,
      0L,
      "Some exploratory score models were not estimable; see the results table."
    )
  )
}

build_score_regression_analysis <- function(cohorts, config) {
  specs <- build_score_model_specs(cohorts$reviewed, config)
  results <- make_score_regression_results(specs, config)
  population_audit <- make_score_model_population_audit(specs)
  list(
    results = results,
    group_descriptives = make_score_group_descriptives(specs),
    population_audit = population_audit,
    timing_status = make_timing_status_table(cohorts$reviewed),
    plots = make_score_regression_plots(results),
    checks = make_score_regression_checks(results, population_audit)
  )
}

write_score_regression_analysis <- function(analysis, config) {
  tables <- list(
    score_regression_results = analysis$results,
    score_group_descriptives = analysis$group_descriptives,
    score_model_population_audit = analysis$population_audit,
    score_timing_status = analysis$timing_status
  )

  purrr::iwalk(tables, function(data, name) {
    csv_path <- file.path(config$output_dirs$tables, paste0(name, ".csv"))
    html_path <- file.path(config$output_dirs$tables, paste0(name, ".html"))
    readr::write_csv(data, csv_path, na = "")
    gt_table <- gt::gt(data) |>
      gt::opt_row_striping() |>
      gt::tab_options(table.font.size = gt::px(12))
    gt::gtsave(gt_table, html_path)
  })

  readr::write_csv(
    analysis$checks,
    file.path(config$output_dirs$qc, "score_model_checks.csv"),
    na = ""
  )

  failure_count <- sum(analysis$checks$status == "FAIL")
  if (failure_count > 0L) {
    stop(
      failure_count,
      " required score-model check(s) failed; see ",
      "output/qc/score_model_checks.csv.",
      call. = FALSE
    )
  }

  ggplot2::ggsave(
    file.path(config$output_dirs$figures, "score_binary_models.png"),
    analysis$plots$binary,
    width = 11,
    height = 12,
    units = "in",
    dpi = 300,
    bg = "white"
  )
  ggplot2::ggsave(
    file.path(config$output_dirs$figures, "score_timing_models.png"),
    analysis$plots$timing,
    width = 11,
    height = 9,
    units = "in",
    dpi = 300,
    bg = "white"
  )
}
