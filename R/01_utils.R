assert_project_root <- function(project_root) {
  rproj_files <- list.files(project_root, pattern = "[.]Rproj$", full.names = TRUE)
  if (length(rproj_files) == 0L) {
    stop(
      "Run the analysis from the GK_MPSE_analysis R project directory.",
      call. = FALSE
    )
  }
  invisible(project_root)
}

assert_files_exist <- function(paths) {
  exists <- vapply(paths, file.exists, logical(1))
  if (any(!exists)) {
    missing <- unlist(paths[!exists], use.names = FALSE)
    stop(
      paste("Missing required input file(s):", paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }
  invisible(paths)
}

assert_required_columns <- function(data, required, dataset_name) {
  missing <- setdiff(required, names(data))
  if (length(missing) > 0L) {
    stop(
      paste0(
        dataset_name,
        " is missing required column(s): ",
        paste(missing, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(data)
}

blank_to_na <- function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  x
}

parse_calendar_date <- function(x) {
  x <- blank_to_na(x)
  as.Date(substr(x, 1L, 10L), format = "%Y-%m-%d")
}

parse_number_quietly <- function(x) {
  suppressWarnings(readr::parse_double(blank_to_na(x), na = c("", "NA")))
}

make_output_dirs <- function(config) {
  invisible(lapply(config$output_dirs, fs::dir_create, recurse = TRUE))
}

add_date_copy <- function(data, source_column, destination_column) {
  if (source_column %in% names(data)) {
    data[[destination_column]] <- parse_calendar_date(data[[source_column]])
  }
  data
}

count_duplicate_keys <- function(data, keys) {
  if (!all(keys %in% names(data))) {
    return(NA_integer_)
  }

  data |>
    dplyr::count(dplyr::across(dplyr::all_of(keys)), name = "n") |>
    dplyr::filter(.data$n > 1L) |>
    nrow()
}
