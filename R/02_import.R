read_all_character <- function(path, delimiter) {
  if (identical(delimiter, "\t")) {
    utils::read.delim(
      path,
      header = TRUE,
      sep = "\t",
      quote = "",
      fill = TRUE,
      comment.char = "",
      colClasses = "character",
      na.strings = c("", "NA"),
      check.names = FALSE,
      stringsAsFactors = FALSE
    ) |>
      tibble::as_tibble()
  } else {
    readr::read_csv(
      path,
      col_types = readr::cols(.default = readr::col_character()),
      na = c("", "NA"),
      show_col_types = FALSE,
      progress = FALSE
    )
  }
}

read_scores_file <- function(path) {
  header <- readLines(path, n = 1L, warn = FALSE)

  if (grepl("mpse_score  mpse_percentile", header, fixed = TRUE)) {
    warning(
      paste0(
        "The synthetic score header contains spaces instead of a tab between ",
        "mpse_score and mpse_percentile. The import repaired that known ",
        "synthetic-file issue without changing the source file."
      ),
      call. = FALSE
    )

    scores <- readr::read_tsv(
      path,
      skip = 1L,
      col_names = c(
        "patient_id",
        "date",
        "days_seen",
        "mpse_score",
        "mpse_percentile"
      ),
      col_types = readr::cols(.default = readr::col_character()),
      na = c("", "NA"),
      show_col_types = FALSE,
      progress = FALSE
    )
  } else {
    scores <- read_all_character(path, "\t")
  }

  scores
}

make_mrn_crosswalk <- function(reviewed_patients) {
  identity_map <- reviewed_patients |>
    dplyr::transmute(
      source_mrn = blank_to_na(.data$mrn),
      canonical_mrn = blank_to_na(.data$mrn),
      mapping_source = "current_mrn"
    )

  if (!"other_mrn" %in% names(reviewed_patients)) {
    reviewed_patients$other_mrn <- NA_character_
  }

  historical_map <- reviewed_patients |>
    dplyr::transmute(
      source_mrn = blank_to_na(.data$other_mrn),
      canonical_mrn = blank_to_na(.data$mrn),
      mapping_source = "other_mrn"
    ) |>
    dplyr::filter(!is.na(.data$source_mrn))

  crosswalk <- dplyr::bind_rows(identity_map, historical_map) |>
    dplyr::filter(!is.na(.data$source_mrn), !is.na(.data$canonical_mrn)) |>
    dplyr::distinct()

  ambiguous <- crosswalk |>
    dplyr::distinct(.data$source_mrn, .data$canonical_mrn) |>
    dplyr::count(.data$source_mrn, name = "n_canonical") |>
    dplyr::filter(.data$n_canonical > 1L)

  if (nrow(ambiguous) > 0L) {
    stop(
      "At least one source MRN maps to more than one canonical MRN.",
      call. = FALSE
    )
  }

  crosswalk |>
    dplyr::arrange(.data$source_mrn, .data$mapping_source) |>
    dplyr::distinct(.data$source_mrn, .keep_all = TRUE)
}

apply_mrn_crosswalk <- function(data, id_column, crosswalk) {
  if (!id_column %in% names(data)) {
    stop(paste("Missing identifier column:", id_column), call. = FALSE)
  }

  data$.source_mrn_for_join <- blank_to_na(data[[id_column]])

  data |>
    dplyr::left_join(
      crosswalk |>
        dplyr::select(dplyr::all_of(c("source_mrn", "canonical_mrn"))),
      by = c(".source_mrn_for_join" = "source_mrn")
    ) |>
    dplyr::mutate(
      "{id_column}" := dplyr::coalesce(
        .data$canonical_mrn,
        .data$.source_mrn_for_join
      )
    ) |>
    dplyr::select(
      -dplyr::all_of(c(".source_mrn_for_join", "canonical_mrn"))
    )
}

import_analysis_sources <- function(config) {
  assert_files_exist(config$input_files)

  reviewed <- read_all_character(config$input_files$reviewed_patients, ",")
  tests <- read_all_character(config$input_files$test_results, ",")
  patients <- read_all_character(config$input_files$mpse_patients, "\t")
  admissions <- read_all_character(config$input_files$mpse_admissions, "\t")
  scores <- read_scores_file(config$input_files$mpse_scores)

  assert_required_columns(
    reviewed,
    c(
      "project_id", "mrn", "dob", "proband_gender", "mpse_admit_date",
      "mpse_review_date", "mpse_cp_score", "mpse_cp_percent", "mpse_rec",
      "mpse_rec_result", "gk_referral_received_date", "gk_enrollment_status",
      "gk_enrollment_date", "gk_enrollment_type"
    ),
    "reviewed patient data"
  )
  assert_required_columns(
    tests,
    c(
      "project_id", "genetic_test", "order_submitted", "result_date",
      "diagnosis_by_lab"
    ),
    "test results data"
  )
  assert_required_columns(
    patients,
    c("patient_id", "sex", "dob", "first_seen", "last_seen"),
    "MPSE patient data"
  )
  assert_required_columns(
    admissions,
    c("patient_id", "unit", "arrive_dtm", "dischg_dtm"),
    "MPSE admission data"
  )
  assert_required_columns(
    scores,
    c("patient_id", "date", "days_seen", "mpse_score"),
    "MPSE score data"
  )

  if ("hpo_source" %in% names(scores)) {
    scores <- scores |>
      dplyr::filter(.data$hpo_source == config$rules$hpo_source)
  }

  reviewed <- reviewed |>
    dplyr::mutate(
      project_id = blank_to_na(.data$project_id),
      mrn = blank_to_na(.data$mrn),
      mpse_rec_result_raw = blank_to_na(.data$mpse_rec_result),
      mpse_rec_result = dplyr::case_when(
        is.na(.data$mpse_rec_result_raw) &
          (
            .data$mpse_rec == "Other" |
              .data$project_id %in%
                config$rules$other_recommendation_result_project_ids
          ) ~ "Other",
        TRUE ~ .data$mpse_rec_result_raw
      ),
      mpse_cp_score = parse_number_quietly(.data$mpse_cp_score),
      mpse_cp_percent = parse_number_quietly(.data$mpse_cp_percent)
    )

  if (!"other_mrn" %in% names(reviewed)) {
    reviewed$other_mrn <- NA_character_
  }

  reviewed_date_columns <- c(
    "dob", "mpse_admit_date", "mpse_genetic_consult_date", "mpse_review_date",
    "gk_approval_date", "gk_referral_received_date", "gk_enrollment_date"
  )
  for (column in intersect(reviewed_date_columns, names(reviewed))) {
    reviewed <- add_date_copy(reviewed, column, paste0(column, "_date"))
  }

  crosswalk <- make_mrn_crosswalk(reviewed)

  patients <- patients |>
    dplyr::mutate(patient_id = blank_to_na(.data$patient_id)) |>
    apply_mrn_crosswalk("patient_id", crosswalk) |>
    add_date_copy("dob", "dob_date") |>
    add_date_copy("first_seen", "first_seen_date") |>
    add_date_copy("last_seen", "last_seen_date")

  admissions <- admissions |>
    dplyr::mutate(patient_id = blank_to_na(.data$patient_id)) |>
    apply_mrn_crosswalk("patient_id", crosswalk) |>
    add_date_copy("arrive_dtm", "arrive_date") |>
    add_date_copy("dischg_dtm", "dischg_date")

  scores <- scores |>
    dplyr::mutate(
      patient_id = blank_to_na(.data$patient_id),
      score_date = parse_calendar_date(.data$date),
      days_seen = suppressWarnings(as.integer(.data$days_seen)),
      mpse_score = parse_number_quietly(.data$mpse_score)
    ) |>
    apply_mrn_crosswalk("patient_id", crosswalk)

  tests <- tests |>
    dplyr::mutate(project_id = blank_to_na(.data$project_id)) |>
    add_date_copy("data_group_record_date", "data_group_record_date_date") |>
    add_date_copy("order_submitted", "order_submitted_date") |>
    add_date_copy("result_date", "result_date_date")

  list(
    reviewed_patients = reviewed,
    test_results = tests,
    mpse_patients = patients,
    mpse_admissions = admissions,
    mpse_scores = scores,
    mrn_crosswalk = crosswalk
  )
}
