make_import_summary <- function(sources) {
  dplyr::bind_rows(
    tibble::tibble(
      dataset = "reviewed_patients",
      rows = nrow(sources$reviewed_patients),
      unique_patients = dplyr::n_distinct(sources$reviewed_patients$mrn),
      unique_projects = dplyr::n_distinct(sources$reviewed_patients$project_id)
    ),
    tibble::tibble(
      dataset = "test_results",
      rows = nrow(sources$test_results),
      unique_patients = NA_integer_,
      unique_projects = dplyr::n_distinct(sources$test_results$project_id)
    ),
    tibble::tibble(
      dataset = "mpse_patients",
      rows = nrow(sources$mpse_patients),
      unique_patients = dplyr::n_distinct(sources$mpse_patients$patient_id),
      unique_projects = NA_integer_
    ),
    tibble::tibble(
      dataset = "mpse_admissions",
      rows = nrow(sources$mpse_admissions),
      unique_patients = dplyr::n_distinct(sources$mpse_admissions$patient_id),
      unique_projects = NA_integer_
    ),
    tibble::tibble(
      dataset = "mpse_scores",
      rows = nrow(sources$mpse_scores),
      unique_patients = dplyr::n_distinct(sources$mpse_scores$patient_id),
      unique_projects = NA_integer_
    )
  )
}

make_identifier_overlap <- function(sources) {
  reviewed_mrns <- unique(stats::na.omit(sources$reviewed_patients$mrn))
  reviewed_projects <- unique(stats::na.omit(sources$reviewed_patients$project_id))

  tibble::tibble(
    comparison = c(
      "reviewed MRNs found in MPSE patients",
      "reviewed MRNs found in MPSE admissions",
      "reviewed MRNs found in MPSE scores",
      "test-result project IDs found in reviewed patients"
    ),
    reference_n = c(
      length(reviewed_mrns),
      length(reviewed_mrns),
      length(reviewed_mrns),
      dplyr::n_distinct(sources$test_results$project_id, na.rm = TRUE)
    ),
    matched_n = c(
      sum(reviewed_mrns %in% sources$mpse_patients$patient_id),
      sum(reviewed_mrns %in% sources$mpse_admissions$patient_id),
      sum(reviewed_mrns %in% sources$mpse_scores$patient_id),
      sum(unique(stats::na.omit(sources$test_results$project_id)) %in% reviewed_projects)
    )
  ) |>
    dplyr::mutate(unmatched_n = .data$reference_n - .data$matched_n)
}

make_duplicate_key_summary <- function(sources) {
  tibble::tibble(
    check = c(
      "reviewed project_id",
      "reviewed MRN",
      "MPSE patient_id",
      "MPSE patient_id + score_date"
    ),
    duplicate_key_groups = c(
      count_duplicate_keys(sources$reviewed_patients, "project_id"),
      count_duplicate_keys(sources$reviewed_patients, "mrn"),
      count_duplicate_keys(sources$mpse_patients, "patient_id"),
      count_duplicate_keys(sources$mpse_scores, c("patient_id", "score_date"))
    )
  )
}

make_missingness_summary <- function(sources) {
  key_fields <- list(
    reviewed_patients = c(
      "project_id", "mrn", "dob_date", "proband_gender",
      "mpse_admit_date_date", "mpse_review_date_date", "mpse_cp_score",
      "mpse_cp_percent", "mpse_rec", "mpse_rec_result",
      "gk_referral_received_date_date", "gk_enrollment_status",
      "gk_enrollment_date_date"
    ),
    test_results = c(
      "project_id", "test_status", "genetic_test", "order_submitted_date",
      "result_date_date", "diagnosis_by_lab"
    ),
    mpse_patients = c(
      "patient_id", "sex", "dob_date", "first_seen_date", "last_seen_date"
    ),
    mpse_admissions = c(
      "patient_id", "unit", "arrive_date", "dischg_date"
    ),
    mpse_scores = c(
      "patient_id", "score_date", "days_seen", "mpse_score"
    )
  )

  purrr::imap_dfr(key_fields, function(fields, dataset_name) {
    data <- sources[[dataset_name]]
    fields <- intersect(fields, names(data))

    purrr::map_dfr(fields, function(field) {
      tibble::tibble(
        dataset = dataset_name,
        variable = field,
        n_missing = sum(is.na(data[[field]])),
        n_total = nrow(data),
        percent_missing = if (nrow(data) == 0L) {
          NA_real_
        } else {
          100 * sum(is.na(data[[field]])) / nrow(data)
        }
      )
    })
  })
}

run_import_qc <- function(sources, config, write_files = TRUE) {
  qc <- list(
    import_summary = make_import_summary(sources),
    identifier_overlap = make_identifier_overlap(sources),
    duplicate_key_summary = make_duplicate_key_summary(sources),
    missingness = make_missingness_summary(sources)
  )

  if (write_files) {
    make_output_dirs(config)
    purrr::iwalk(qc, function(data, name) {
      readr::write_csv(
        data,
        file.path(config$output_dirs$qc, paste0(name, ".csv")),
        na = ""
      )
    })
  }

  unmatched <- qc$identifier_overlap |>
    dplyr::filter(.data$unmatched_n > 0L)
  if (nrow(unmatched) > 0L) {
    warning(
      "Some identifiers do not match across sources; see output/qc/identifier_overlap.csv.",
      call. = FALSE
    )
  }

  qc
}
