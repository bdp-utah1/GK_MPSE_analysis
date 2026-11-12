library(tidyverse)


# Analysis settings --------------------------------------------------------

analysis_output_dir <- Sys.getenv(
  "PHENOTYPE_ANALYSIS_OUTPUT_DIR",
  unset = "output/tables/phenotype_analysis"
)
min_term_patients_for_model <- 5L
hpo_ontology_release <- "2025-05-06"
hpo_ontology_path <- "docs/hp_v2025-05-06.json"
training_hpo_release <- "2018-12-21"
training_hpo_path <- "mpse_data/HPO_versions/v2018-12-21.obo"
training_raw_hpo_codes_path <- paste0(
  "mpse_data/HPO_versions/training_dataset_raw_codes.txt"
)
training_clean_hpo_codes_path <- paste0(
  "mpse_data/HPO_versions/training_dataset_clean_codes.txt"
)
hpo_modifier_roots <- c(
  "HP:0000005", # Mode of inheritance
  "HP:0012823", # Clinical modifier
  "HP:0040279", # Frequency
  "HP:0031797", # Clinical course
  "HP:0032223", # Blood group
  "HP:0032443"  # Past medical history
)

dir.create(analysis_output_dir, recursive = TRUE, showWarnings = FALSE)


# Import and assemble the analysis data -----------------------------------

modl_feat <- read_tsv(
  "mpse_data/mpse_model_feature_coefficients.tsv",
  show_col_types = FALSE
)


# Build the cohort-release HPO graph --------------------------------------

uri_to_hpo <- function(value) {
  str_replace(
    value,
    "^http://purl\\.obolibrary\\.org/obo/HP_",
    "HP:"
  )
}

extract_obsolete_replacement <- function(node) {
  properties <- node$meta$basicPropertyValues
  if (is.null(properties)) {
    return(NA_character_)
  }

  replacement_properties <- keep(
    properties,
    ~ identical(
      .x$pred,
      "http://purl.obolibrary.org/obo/IAO_0100001"
    )
  )

  if (length(replacement_properties) == 0L) {
    return(NA_character_)
  }

  uri_to_hpo(replacement_properties[[1]]$val)
}

extract_node_property_values <- function(node, predicate) {
  properties <- node$meta$basicPropertyValues
  if (is.null(properties)) {
    return(character())
  }

  properties %>%
    keep(~ identical(.x$pred, predicate)) %>%
    map_chr("val")
}

extract_node_creation_date <- function(node) {
  properties <- node$meta$basicPropertyValues
  if (is.null(properties)) {
    return(NA_character_)
  }

  date_values <- properties %>%
    keep(~ str_detect(.x$pred, "(creation_date|/date)$")) %>%
    map_chr("val")

  first(date_values, default = NA_character_)
}

extract_hpo_node_table <- function(ontology_graph) {
  ontology_graph$nodes %>%
    keep(~ str_detect(.x$id, "/HP_[0-9]{7}$")) %>%
    map_dfr(~ tibble(
      hpo = uri_to_hpo(.x$id),
      hpo_description = if (is.null(.x$lbl)) NA_character_ else .x$lbl,
      is_obsolete = isTRUE(.x$meta$deprecated),
      replacement_hpo = extract_obsolete_replacement(.x),
      alternative_hpo = paste(
        extract_node_property_values(
          .x,
          paste0(
            "http://www.geneontology.org/formats/",
            "oboInOwl#hasAlternativeId"
          )
        ),
        collapse = "; "
      ),
      term_creation_date = extract_node_creation_date(.x)
    ))
}

ontology_graph <- jsonlite::fromJSON(
  hpo_ontology_path,
  simplifyVector = FALSE
)$graphs[[1]]

ontology_version <- str_extract(
  ontology_graph$meta$version,
  "[0-9]{4}-[0-9]{2}-[0-9]{2}"
)

if (!identical(ontology_version, hpo_ontology_release)) {
  stop(
    "The HPO JSON release does not match the configured ontology release."
  )
}

ontology_node_table <- extract_hpo_node_table(ontology_graph)

ontology_edge_table <- ontology_graph$edges %>%
  map_dfr(~ tibble(
    child = uri_to_hpo(.x$sub),
    parent = uri_to_hpo(.x$obj),
    predicate = .x$pred
  )) %>%
  filter(
    predicate == "is_a",
    child %in% ontology_node_table$hpo,
    parent %in% ontology_node_table$hpo
  )

ontology_parents <- set_names(
  vector("list", nrow(ontology_node_table)),
  ontology_node_table$hpo
)
edge_parents <- split(
  ontology_edge_table$parent,
  ontology_edge_table$child
)
ontology_parents[names(edge_parents)] <- edge_parents

hpo_ontology <- ontologyIndex::ontology_index(
  parents = ontology_parents,
  id = ontology_node_table$hpo,
  name = set_names(
    coalesce(
      ontology_node_table$hpo_description,
      ontology_node_table$hpo
    ),
    ontology_node_table$hpo
  ),
  obsolete = set_names(
    ontology_node_table$is_obsolete,
    ontology_node_table$hpo
  ),
  version = hpo_ontology_release
)

if (any(!hpo_modifier_roots %in% hpo_ontology$id)) {
  stop("At least one PyHPO modifier root is absent from the HPO release.")
}

invalid_replacements <- ontology_node_table %>%
  filter(
    is_obsolete,
    !is.na(replacement_hpo),
    !replacement_hpo %in% hpo
  )

if (nrow(invalid_replacements) > 0L) {
  stop("At least one obsolete HPO replacement is absent from the ontology.")
}

# PyHPO treats each modifier root and all of its descendants as modifiers.
hpo_modifier_terms <- ontologyIndex::get_descendants(
  hpo_ontology,
  roots = hpo_modifier_roots,
  exclude_roots = FALSE
)


# Build the confirmed training-release HPO graph -------------------------

training_hpo_header <- read_lines(training_hpo_path, n_max = 50L)
training_ontology_version <- training_hpo_header %>%
  str_subset("^data-version: releases/[0-9]{4}-[0-9]{2}-[0-9]{2}$") %>%
  str_extract("[0-9]{4}-[0-9]{2}-[0-9]{2}") %>%
  first(default = NA_character_)

if (!identical(training_ontology_version, training_hpo_release)) {
  stop("The training HPO OBO file does not match the configured release.")
}

training_hpo_ontology <- ontologyIndex::get_ontology(
  training_hpo_path,
  propagate_relationships = "is_a",
  extract_tags = "everything",
  merge_equivalent_terms = FALSE
)

training_hpo_term_table <- tibble(
  training_hpo = unname(training_hpo_ontology$id),
  training_hpo_description = unname(
    training_hpo_ontology$name[training_hpo_ontology$id]
  ),
  training_hpo_is_obsolete = unname(
    training_hpo_ontology$obsolete[training_hpo_ontology$id]
  )
)

training_hpo_alternative_id_lookup <- imap_dfr(
  training_hpo_ontology$alt_id,
  ~ {
    if (length(.x) == 0L) {
      return(tibble(
        training_hpo_alternative_id = character(),
        training_hpo = character()
      ))
    }

    tibble(
      training_hpo_alternative_id = .x,
      training_hpo = .y
    )
  }
) %>%
  filter(training_hpo_alternative_id != "") %>%
  distinct()

read_hpo_code_set <- function(path) {
  codes <- read_lines(path) %>%
    str_trim() %>%
    discard(~ .x == "")

  if (anyDuplicated(codes)) {
    stop("Duplicate HPO codes were found in: ", path)
  }
  if (any(!str_detect(codes, "^HP:[0-9]{7}$"))) {
    stop("Invalid HPO codes were found in: ", path)
  }

  codes
}

training_raw_hpo_codes <- read_hpo_code_set(training_raw_hpo_codes_path)
training_clean_hpo_codes <- read_hpo_code_set(training_clean_hpo_codes_path)

if (length(setdiff(training_clean_hpo_codes, training_raw_hpo_codes)) > 0L) {
  stop("Clean training HPO codes are not a subset of the raw codes.")
}

if (!setequal(training_clean_hpo_codes, modl_feat$codes)) {
  stop("Clean training HPO codes do not match the model feature vocabulary.")
}

unknown_training_hpo_codes <- setdiff(
  training_raw_hpo_codes,
  c(
    training_hpo_term_table$training_hpo,
    training_hpo_alternative_id_lookup$training_hpo_alternative_id
  )
)

if (length(unknown_training_hpo_codes) > 0L) {
  stop("Training code sets contain IDs absent from the 2018 HPO release.")
}

training_raw_only_hpo_codes <- setdiff(
  training_raw_hpo_codes,
  training_clean_hpo_codes
)

training_clean_ancestor_codes <- training_hpo_ontology$ancestors[
  training_clean_hpo_codes
] %>%
  unlist(use.names = FALSE) %>%
  unique()

training_parent_removed_hpo_codes <- intersect(
  training_raw_only_hpo_codes,
  training_clean_ancestor_codes
)

training_other_removed_hpo_codes <- setdiff(
  training_raw_only_hpo_codes,
  training_parent_removed_hpo_codes
)

training_retained_descendant_lookup <- map_dfr(
  training_parent_removed_hpo_codes,
  function(parent_hpo) {
    retained_descendants <- training_clean_hpo_codes[
      map_lgl(
        training_hpo_ontology$ancestors[training_clean_hpo_codes],
        ~ parent_hpo %in% .x
      )
    ]

    tibble(
      training_hpo = parent_hpo,
      retained_descendant_count = length(retained_descendants),
      retained_descendant_hpo = paste(
        sort(retained_descendants),
        collapse = "; "
      ),
      retained_descendant_description = paste(
        sort(unique(training_hpo_ontology$name[retained_descendants])),
        collapse = "; "
      )
    )
  }
)


defunct_mrn_4879 <- "X014455166"
current_mrn_4879 <- "X014453675"

cern_unmerged_mrn_4048 <- "650480153"
epic_unmerged_mrn_4048 <- "X014183906"

mosc_patnt <- read_csv(
  "~/Documents/phd/yandell/gene_kids/data/GK_MPSE_analysis_copy/mosaic_data/patient_level_data.csv",
  show_col_types = FALSE
) %>%
  select(
    project_id,
    mrn,
    dob,
    mpse_admit_date,
    mpse_review_date,
    mpse_cp_percent
  )

mosc_outcm <- read_csv(
  "~/Documents/phd/yandell/gene_kids/data/GK_MPSE_analysis_copy/output/data/reviewed_cohort_analysis.csv",
  show_col_types = FALSE
) %>%
  select(
    project_id,
    gk_enrolled,
    has_any_test,
    has_classifiable_completed_result,
    patient_result
  )

mosc_data <- inner_join(mosc_patnt, mosc_outcm, by = "project_id")


mpse_patnt <- read_tsv(
  "~/Documents/phd/yandell/gene_kids/data/GK_MPSE_analysis_copy/mpse_data/filtered_patient.tsv",
  show_col_types = FALSE
)

mpse_pheno <- read_tsv(
  "~/Documents/phd/yandell/gene_kids/data/GK_MPSE_analysis_copy/mpse_data/filtered_patient_hpo.tsv",
  show_col_types = FALSE
) %>%
  select(-note_dtm)

mpse_data <- inner_join(mpse_patnt, mpse_pheno, by = "patient_id") %>%
  mutate(
    patient_id = str_replace(
      patient_id,
      defunct_mrn_4879,
      current_mrn_4879
    ),
    patient_id = str_replace(
      patient_id,
      cern_unmerged_mrn_4048,
      epic_unmerged_mrn_4048
    )
  ) %>%
  rename(mrn = patient_id)

df <- full_join(mosc_data, mpse_data, by = join_by(mrn, dob)) %>%
  mutate(
    subgrp = factor(
      case_when(
        is.na(project_id) ~ "Unreviewed",
        !has_any_test ~ "Reviewed & Untested",
        has_classifiable_completed_result ~ patient_result,
        TRUE ~ NA_character_
      ),
      levels = c(
        "Positive",
        "Uncertain",
        "Negative",
        "Reviewed & Untested",
        "Unreviewed"
      )
    )
  ) %>%
  select(
    mrn,
    dob,
    mpse_admit_date,
    mpse_review_date,
    mpse_cp_percent,
    subgrp,
    first_seen,
    load_date,
    note_id,
    note_type,
    hpo
  )


# Define the reviewed cohort at the patient level -------------------------

reviewed_patients <- df %>%
  filter(subgrp != "Unreviewed") %>%
  distinct(
    mrn,
    dob,
    mpse_admit_date,
    mpse_review_date,
    mpse_cp_percent,
    subgrp
  ) %>%
  mutate(
    patient_key = paste(mrn, dob, sep = "|"),
    days_admission_to_review = as.integer(
      mpse_review_date - mpse_admit_date
    ),
    review_weeks = days_admission_to_review / 7
  )

if (anyDuplicated(reviewed_patients$patient_key)) {
  stop("Reviewed cohort contains duplicate MRN/date-of-birth patient keys.")
}

required_reviewed_fields <- c(
  "subgrp",
  "mpse_admit_date",
  "mpse_review_date",
  "mpse_cp_percent",
  "days_admission_to_review"
)

if (anyNA(reviewed_patients[required_reviewed_fields])) {
  stop("Reviewed patients have missing subgroup, date, or percentile fields.")
}

if (any(reviewed_patients$days_admission_to_review < 0)) {
  stop("At least one clinical review date precedes the admission date.")
}

if (any(!reviewed_patients$subgrp %in% c(
  "Positive",
  "Uncertain",
  "Negative",
  "Reviewed & Untested"
))) {
  stop("The reviewed cohort contains an unexpected diagnostic category.")
}


# Restrict phenotype data to terms available by clinical review -----------

# MPSE uses cumulative phenotype data. There is intentionally no lower date
# bound: every term loaded on or before clinical review is considered present.
reviewed_hpo <- df %>%
  filter(
    subgrp != "Unreviewed",
    !is.na(hpo),
    !is.na(load_date),
    load_date <= mpse_review_date
  ) %>%
  mutate(patient_key = paste(mrn, dob, sep = "|")) %>%
  select(patient_key, load_date, note_id, note_type, hpo)

if (any(!str_detect(reviewed_hpo$hpo, "^HP:[0-9]{7}$"))) {
  stop("At least one phenotype value is not a valid seven-digit HPO code.")
}

# One row per patient-term implements cumulative binary term presence before
# applying the BasicHPOSet-compatible preprocessing used by MPSE.
raw_term_presence <- reviewed_hpo %>%
  distinct(patient_key, hpo)

unknown_hpo_terms <- setdiff(
  unique(raw_term_presence$hpo),
  hpo_ontology$id
)

if (length(unknown_hpo_terms) > 0L) {
  stop(
    paste0(
      "Pre-review phenotype data contain HPO terms absent from release ",
      hpo_ontology_release,
      ": ",
      paste(unknown_hpo_terms, collapse = ", ")
    )
  )
}

# BasicHPOSet behavior is reproduced in three explicit stages:
#   1. replace obsolete terms, dropping those without a replacement;
#   2. remove modifier terms and descendants of modifier roots;
#   3. remove a parent whenever a descendant is present for the patient.
term_resolution <- raw_term_presence %>%
  left_join(
    ontology_node_table %>%
      select(hpo, is_obsolete, replacement_hpo),
    by = "hpo"
  ) %>%
  mutate(
    resolved_hpo = case_when(
      !is_obsolete ~ hpo,
      !is.na(replacement_hpo) ~ replacement_hpo,
      TRUE ~ NA_character_
    )
  )

resolved_term_presence <- term_resolution %>%
  filter(!is.na(resolved_hpo)) %>%
  distinct(patient_key, hpo = resolved_hpo)

non_modifier_term_presence <- resolved_term_presence %>%
  filter(!hpo %in% hpo_modifier_terms)

term_presence <- non_modifier_term_presence %>%
  group_by(patient_key) %>%
  summarise(
    hpo = list(ontologyIndex::minimal_set(hpo_ontology, hpo)),
    .groups = "drop"
  ) %>%
  unnest_longer(hpo)

if (any(term_presence$hpo %in% hpo_modifier_terms)) {
  stop("Modifier HPO terms remain after preprocessing.")
}

if (any(
  ontology_node_table$is_obsolete[
    match(term_presence$hpo, ontology_node_table$hpo)
  ]
)) {
  stop("Obsolete HPO terms remain after preprocessing.")
}

raw_patient_counts <- reviewed_hpo %>%
  group_by(patient_key) %>%
  summarise(
    term_yielding_note_count = n_distinct(note_id),
    raw_unique_hpo_count = n_distinct(hpo),
    .groups = "drop"
  )

obsolete_patient_counts <- term_resolution %>%
  group_by(patient_key) %>%
  summarise(
    obsolete_hpo_input_count = sum(is_obsolete),
    obsolete_hpo_replaced_count = sum(
      is_obsolete & !is.na(replacement_hpo)
    ),
    obsolete_hpo_removed_without_replacement_count = sum(
      is_obsolete & is.na(replacement_hpo)
    ),
    .groups = "drop"
  )

patient_phenotype_counts <- reviewed_patients %>%
  left_join(raw_patient_counts, by = "patient_key") %>%
  left_join(obsolete_patient_counts, by = "patient_key") %>%
  left_join(
    resolved_term_presence %>%
      count(patient_key, name = "post_obsolete_unique_hpo_count"),
    by = "patient_key"
  ) %>%
  left_join(
    non_modifier_term_presence %>%
      count(patient_key, name = "post_modifier_unique_hpo_count"),
    by = "patient_key"
  ) %>%
  left_join(
    term_presence %>%
      count(patient_key, name = "model_retained_hpo_count"),
    by = "patient_key"
  ) %>%
  mutate(
    across(
      c(
        term_yielding_note_count,
        raw_unique_hpo_count,
        obsolete_hpo_input_count,
        obsolete_hpo_replaced_count,
        obsolete_hpo_removed_without_replacement_count,
        post_obsolete_unique_hpo_count,
        post_modifier_unique_hpo_count,
        model_retained_hpo_count
      ),
      ~ as.integer(replace_na(.x, 0L))
    ),
    modifier_hpo_count_removed =
      post_obsolete_unique_hpo_count - post_modifier_unique_hpo_count,
    parent_hpo_count_removed =
      post_modifier_unique_hpo_count - model_retained_hpo_count,
    hpo_terms_removed_by_preprocessing =
      raw_unique_hpo_count - model_retained_hpo_count
  )

model_term_lookup <- modl_feat %>%
  transmute(
    hpo = codes,
    hpo_description = code_descrip,
    model_coefficient = coef,
    training_control_patients_with_term = as.integer(ctrl_cnt),
    training_case_patients_with_term = as.integer(case_cnt),
    training_patients_with_term = as.integer(ctrl_cnt + case_cnt)
  )

term_annotation_lookup <- ontology_node_table %>%
  select(hpo, hpo_description) %>%
  left_join(
    model_term_lookup %>% select(-hpo_description),
    by = "hpo"
  )


# Reusable statistical helpers -------------------------------------------

summarize_numeric_metrics <- function(data, group_variable, metrics) {
  overall <- data %>%
    mutate(summary_group = "Overall")

  grouped <- data %>%
    mutate(summary_group = as.character(.data[[group_variable]]))

  bind_rows(overall, grouped) %>%
    pivot_longer(
      cols = all_of(metrics),
      names_to = "metric",
      values_to = "value"
    ) %>%
    group_by(summary_group, metric) %>%
    summarise(
      n = sum(!is.na(value)),
      missing = sum(is.na(value)),
      mean = mean(value, na.rm = TRUE),
      sd = sd(value, na.rm = TRUE),
      median = median(value, na.rm = TRUE),
      q1 = quantile(value, 0.25, na.rm = TRUE, names = FALSE),
      q3 = quantile(value, 0.75, na.rm = TRUE, names = FALSE),
      minimum = min(value, na.rm = TRUE),
      maximum = max(value, na.rm = TRUE),
      .groups = "drop"
    )
}

wilcoxon_group_comparison <- function(data, metric) {
  positive <- data %>%
    filter(subgrp == "Positive") %>%
    pull(all_of(metric))
  negative <- data %>%
    filter(subgrp == "Negative") %>%
    pull(all_of(metric))

  test <- wilcox.test(positive, negative, exact = FALSE)

  tibble(
    metric = metric,
    positive_n = sum(!is.na(positive)),
    negative_n = sum(!is.na(negative)),
    positive_median = median(positive, na.rm = TRUE),
    positive_q1 = quantile(positive, 0.25, na.rm = TRUE, names = FALSE),
    positive_q3 = quantile(positive, 0.75, na.rm = TRUE, names = FALSE),
    negative_median = median(negative, na.rm = TRUE),
    negative_q1 = quantile(negative, 0.25, na.rm = TRUE, names = FALSE),
    negative_q3 = quantile(negative, 0.75, na.rm = TRUE, names = FALSE),
    median_difference_positive_minus_negative =
      median(positive, na.rm = TRUE) - median(negative, na.rm = TRUE),
    wilcoxon_w = unname(test$statistic),
    p_value = test$p.value
  )
}

fit_adjusted_count_model <- function(data, outcome) {
  formula <- reformulate(
    c("result_positive", "review_weeks"),
    response = outcome
  )

  negative_binomial_fit <- tryCatch(
    suppressWarnings(MASS::glm.nb(formula, data = data)),
    error = function(e) NULL
  )

  if (is.null(negative_binomial_fit)) {
    fit <- glm(formula, data = data, family = quasipoisson())
    model_family <- "Quasi-Poisson fallback"
  } else {
    fit <- negative_binomial_fit
    model_family <- "Negative binomial"
  }

  coefficient_table <- coef(summary(fit))
  estimate <- coefficient_table["result_positive", "Estimate"]
  standard_error <- coefficient_table["result_positive", "Std. Error"]
  p_value <- coefficient_table[
    "result_positive",
    ncol(coefficient_table)
  ]

  tibble(
    metric = outcome,
    model_family = model_family,
    comparison = "Positive vs Negative",
    adjusted_incidence_rate_ratio = exp(estimate),
    confidence_lower_95 = exp(estimate - 1.96 * standard_error),
    confidence_upper_95 = exp(estimate + 1.96 * standard_error),
    p_value = p_value,
    adjustment = "Admission-to-review time (weeks)"
  )
}

safe_fisher_test <- function(
    positive_with,
    positive_without,
    negative_with,
    negative_without) {
  contingency_table <- matrix(
    c(
      positive_with,
      positive_without,
      negative_with,
      negative_without
    ),
    nrow = 2,
    byrow = TRUE
  )

  test <- fisher.test(contingency_table, conf.int = TRUE)

  tibble(
    odds_ratio = unname(test$estimate),
    confidence_lower_95 = test$conf.int[[1]],
    confidence_upper_95 = test$conf.int[[2]],
    p_value = test$p.value
  )
}

fit_term_result_model <- function(hpo_code, patient_data, term_patient_keys) {
  term_present <- as.integer(
    patient_data$patient_key %in% term_patient_keys[[hpo_code]]
  )
  result_positive <- patient_data$result_positive
  contingency_table <- table(
    factor(result_positive, levels = c(0, 1)),
    factor(term_present, levels = c(0, 1))
  )

  model_is_supported <-
    sum(term_present) >= min_term_patients_for_model &&
    sum(term_present == 0) >= min_term_patients_for_model &&
    all(contingency_table > 0)

  if (!model_is_supported) {
    return(tibble(
      adjusted_odds_ratio = NA_real_,
      adjusted_confidence_lower_95 = NA_real_,
      adjusted_confidence_upper_95 = NA_real_,
      adjusted_p_value = NA_real_,
      adjusted_model_status = "Sparse/separated; Fisher exact test only"
    ))
  }

  model_data <- patient_data %>%
    mutate(term_present = term_present)

  fit <- tryCatch(
    suppressWarnings(glm(
      term_present ~ result_positive + review_weeks,
      data = model_data,
      family = binomial()
    )),
    error = function(e) NULL
  )

  if (is.null(fit) || !fit$converged) {
    return(tibble(
      adjusted_odds_ratio = NA_real_,
      adjusted_confidence_lower_95 = NA_real_,
      adjusted_confidence_upper_95 = NA_real_,
      adjusted_p_value = NA_real_,
      adjusted_model_status = "Adjusted logistic model did not converge"
    ))
  }

  coefficient_table <- coef(summary(fit))
  estimate <- coefficient_table["result_positive", "Estimate"]
  standard_error <- coefficient_table["result_positive", "Std. Error"]
  p_value <- coefficient_table["result_positive", "Pr(>|z|)"]

  if (any(!is.finite(c(estimate, standard_error, p_value)))) {
    return(tibble(
      adjusted_odds_ratio = NA_real_,
      adjusted_confidence_lower_95 = NA_real_,
      adjusted_confidence_upper_95 = NA_real_,
      adjusted_p_value = NA_real_,
      adjusted_model_status = "Adjusted logistic estimate was unstable"
    ))
  }

  tibble(
    adjusted_odds_ratio = exp(estimate),
    adjusted_confidence_lower_95 = exp(estimate - 1.96 * standard_error),
    adjusted_confidence_upper_95 = exp(estimate + 1.96 * standard_error),
    adjusted_p_value = p_value,
    adjusted_model_status = "Fitted"
  )
}

fit_term_percentile_model <- function(
    hpo_code,
    patient_data,
    term_patient_keys) {
  term_present <- as.integer(
    patient_data$patient_key %in% term_patient_keys[[hpo_code]]
  )
  with_term <- patient_data$mpse_cp_percent[term_present == 1]
  without_term <- patient_data$mpse_cp_percent[term_present == 0]
  with_term_days <- patient_data$days_admission_to_review[term_present == 1]
  without_term_days <- patient_data$days_admission_to_review[term_present == 0]

  model_is_supported <-
    length(with_term) >= min_term_patients_for_model &&
    length(without_term) >= min_term_patients_for_model

  descriptive_result <- tibble(
    n_with_term = length(with_term),
    n_without_term = length(without_term),
    patient_percent_with_term = 100 * length(with_term) / nrow(patient_data),
    percentile_median_with_term = median(with_term),
    percentile_q1_with_term = quantile(
      with_term,
      0.25,
      names = FALSE
    ),
    percentile_q3_with_term = quantile(
      with_term,
      0.75,
      names = FALSE
    ),
    percentile_median_without_term = median(without_term),
    percentile_q1_without_term = quantile(
      without_term,
      0.25,
      names = FALSE
    ),
    percentile_q3_without_term = quantile(
      without_term,
      0.75,
      names = FALSE
    ),
    percentile_median_difference =
      median(with_term) - median(without_term),
    review_days_median_with_term = median(with_term_days),
    review_days_q1_with_term = quantile(
      with_term_days,
      0.25,
      names = FALSE
    ),
    review_days_q3_with_term = quantile(
      with_term_days,
      0.75,
      names = FALSE
    ),
    review_days_median_without_term = median(without_term_days),
    review_days_q1_without_term = quantile(
      without_term_days,
      0.25,
      names = FALSE
    ),
    review_days_q3_without_term = quantile(
      without_term_days,
      0.75,
      names = FALSE
    )
  )

  if (!model_is_supported) {
    return(descriptive_result %>%
      mutate(
        wilcoxon_w = NA_real_,
        p_value = NA_real_,
        adjusted_percentile_difference = NA_real_,
        adjusted_confidence_lower_95 = NA_real_,
        adjusted_confidence_upper_95 = NA_real_,
        adjusted_p_value = NA_real_,
        adjusted_model_status = paste0(
          "Fewer than ",
          min_term_patients_for_model,
          " patients in one comparison group"
        )
      ))
  }

  rank_test <- wilcox.test(with_term, without_term, exact = FALSE)
  model_data <- patient_data %>%
    mutate(term_present = term_present)
  fit <- lm(
    mpse_cp_percent ~ term_present + review_weeks,
    data = model_data
  )
  coefficient_table <- coef(summary(fit))
  estimate <- coefficient_table["term_present", "Estimate"]
  standard_error <- coefficient_table["term_present", "Std. Error"]
  critical_value <- qt(0.975, df = fit$df.residual)

  descriptive_result %>%
    mutate(
      wilcoxon_w = unname(rank_test$statistic),
      p_value = rank_test$p.value,
      adjusted_percentile_difference = estimate,
      adjusted_confidence_lower_95 =
        estimate - critical_value * standard_error,
      adjusted_confidence_upper_95 =
        estimate + critical_value * standard_error,
      adjusted_p_value = coefficient_table[
        "term_present",
        "Pr(>|t|)"
      ],
      adjusted_model_status = "Fitted"
    )
}


# 1. Cohort terms that were absent from model training --------------------

cohort_term_frequency <- term_presence %>%
  count(hpo, name = "cohort_patients_with_term") %>%
  mutate(
    reviewed_cohort_n = nrow(reviewed_patients),
    cohort_patient_percent =
      100 * cohort_patients_with_term / reviewed_cohort_n
  )

raw_terms_absent_from_training <- raw_term_presence %>%
  distinct(hpo) %>%
  anti_join(model_term_lookup, by = "hpo")

terms_absent_from_training_base <- cohort_term_frequency %>%
  anti_join(model_term_lookup, by = "hpo") %>%
  left_join(
    ontology_node_table %>% select(hpo, hpo_description),
    by = "hpo"
  ) %>%
  relocate(hpo_description, .after = hpo) %>%
  arrange(desc(cohort_patients_with_term), hpo)

# Map current HPO identifiers to their 2018 primary identifiers. This checks
# direct primary IDs, 2018 alternative IDs, and historical IDs retained as
# alternatives on the current term.
current_term_alternative_ids <- ontology_node_table %>%
  filter(alternative_hpo != "") %>%
  separate_longer_delim(alternative_hpo, delim = "; ") %>%
  transmute(
    hpo,
    current_release_alternative_hpo = alternative_hpo
  )

direct_training_primary_mappings <- terms_absent_from_training_base %>%
  inner_join(
    training_hpo_term_table,
    by = c("hpo" = "training_hpo")
  ) %>%
  transmute(
    hpo,
    training_hpo = hpo,
    training_hpo_description,
    training_id_mapping_method = "Same primary HPO identifier"
  )

training_alternative_id_mappings <- terms_absent_from_training_base %>%
  inner_join(
    training_hpo_alternative_id_lookup,
    by = c("hpo" = "training_hpo_alternative_id")
  ) %>%
  left_join(training_hpo_term_table, by = "training_hpo") %>%
  transmute(
    hpo,
    training_hpo,
    training_hpo_description,
    training_id_mapping_method = "Current ID was a 2018 alternative ID"
  )

current_alternative_to_training_primary_mappings <-
  current_term_alternative_ids %>%
  inner_join(
    training_hpo_term_table,
    by = c("current_release_alternative_hpo" = "training_hpo")
  ) %>%
  transmute(
    hpo,
    training_hpo = current_release_alternative_hpo,
    training_hpo_description,
    training_id_mapping_method = paste0(
      "Current term contains the 2018 primary ID as an alternative ID"
    )
  )

current_alternative_to_training_alternative_mappings <-
  current_term_alternative_ids %>%
  inner_join(
    training_hpo_alternative_id_lookup,
    by = c(
      "current_release_alternative_hpo" =
        "training_hpo_alternative_id"
    )
  ) %>%
  left_join(training_hpo_term_table, by = "training_hpo") %>%
  transmute(
    hpo,
    training_hpo,
    training_hpo_description,
    training_id_mapping_method = paste0(
      "Current term contains a 2018 alternative ID"
    )
  )

training_concept_mappings <- bind_rows(
  direct_training_primary_mappings,
  training_alternative_id_mappings,
  current_alternative_to_training_primary_mappings,
  current_alternative_to_training_alternative_mappings
) %>%
  group_by(hpo, training_hpo, training_hpo_description) %>%
  summarise(
    training_id_mapping_method = paste(
      sort(unique(training_id_mapping_method)),
      collapse = "; "
    ),
    .groups = "drop"
  ) %>%
  mutate(
    present_in_training_raw_codes =
      training_hpo %in% training_raw_hpo_codes,
    present_in_training_clean_codes =
      training_hpo %in% training_clean_hpo_codes,
    removed_as_parent_during_training =
      training_hpo %in% training_parent_removed_hpo_codes,
    removed_by_other_training_preprocessing =
      training_hpo %in% training_other_removed_hpo_codes
  ) %>%
  left_join(training_retained_descendant_lookup, by = "training_hpo")

collapse_character_values <- function(values) {
  values <- values[!is.na(values) & values != ""] %>%
    unique() %>%
    sort()
  if (length(values) == 0L) {
    return(NA_character_)
  }
  paste(values, collapse = "; ")
}

training_concept_status <- training_concept_mappings %>%
  group_by(hpo) %>%
  summarise(
    training_hpo_ids = collapse_character_values(training_hpo),
    training_hpo_descriptions = collapse_character_values(
      training_hpo_description
    ),
    training_id_mapping_methods = collapse_character_values(
      training_id_mapping_method
    ),
    training_term_existed = TRUE,
    matching_training_raw_hpo = collapse_character_values(
      training_hpo[training_hpo %in% training_raw_hpo_codes]
    ),
    matching_training_clean_hpo = collapse_character_values(
      training_hpo[training_hpo %in% training_clean_hpo_codes]
    ),
    present_in_training_raw_codes = any(present_in_training_raw_codes),
    present_in_training_clean_codes = any(present_in_training_clean_codes),
    removed_as_parent_during_training = any(
      removed_as_parent_during_training
    ),
    removed_by_other_training_preprocessing = any(
      removed_by_other_training_preprocessing
    ),
    retained_descendant_count = sum(
      replace_na(retained_descendant_count, 0L)
    ),
    retained_descendant_hpo = collapse_character_values(
      retained_descendant_hpo
    ),
    retained_descendant_description = collapse_character_values(
      retained_descendant_description
    ),
    .groups = "drop"
  )

terms_absent_from_training <- terms_absent_from_training_base %>%
  left_join(training_concept_status, by = "hpo") %>%
  left_join(
    ontology_node_table %>%
      select(
        hpo,
        current_release_alternative_hpo = alternative_hpo,
        hpo_term_creation_date = term_creation_date
      ),
    by = "hpo"
  ) %>%
  mutate(
    training_hpo_release = training_hpo_release,
    across(
      c(
        training_term_existed,
        present_in_training_raw_codes,
        present_in_training_clean_codes,
        removed_as_parent_during_training,
        removed_by_other_training_preprocessing
      ),
      ~ replace_na(.x, FALSE)
    ),
    retained_descendant_count = as.integer(
      replace_na(retained_descendant_count, 0)
    ),
    classification_code = case_when(
      present_in_training_clean_codes ~
        "represented_under_historical_identifier",
      present_in_training_raw_codes & removed_as_parent_during_training ~
        "removed_by_training_parent_trimming",
      present_in_training_raw_codes &
        removed_by_other_training_preprocessing ~
        "removed_by_other_training_preprocessing",
      training_term_existed ~ "not_observed_in_training_raw_data",
      TRUE ~ "new_hpo_term_after_training_release"
    ),
    classification = case_when(
      classification_code == "represented_under_historical_identifier" ~
        paste0(
          "Represented in the model under a historical HPO identifier"
        ),
      classification_code == "removed_by_training_parent_trimming" ~
        paste0(
          "Observed in raw training data but removed because a more ",
          "specific descendant was retained"
        ),
      classification_code ==
        "removed_by_other_training_preprocessing" ~
        "Observed in raw training data but removed by another preprocessing step",
      classification_code == "not_observed_in_training_raw_data" ~
        paste0(
          "Present in the 2018 HPO release but not observed in raw ",
          "training data"
        ),
      TRUE ~ "HPO concept added after the 2018 training release"
    )
  ) %>%
  relocate(
    classification_code,
    classification,
    training_hpo_release,
    .after = hpo_description
  ) %>%
  arrange(
    factor(
      classification_code,
      levels = c(
        "new_hpo_term_after_training_release",
        "represented_under_historical_identifier",
        "removed_by_training_parent_trimming",
        "removed_by_other_training_preprocessing",
        "not_observed_in_training_raw_data"
      )
    ),
    desc(cohort_patients_with_term),
    hpo
  )

term_classification_summary <- terms_absent_from_training %>%
  group_by(
    training_hpo_release,
    classification_code,
    classification
  ) %>%
  summarise(
    term_count = n(),
    term_percent = 100 * n() / nrow(terms_absent_from_training),
    cohort_patient_term_pairs = sum(cohort_patients_with_term),
    .groups = "drop"
  ) %>%
  arrange(desc(term_count), classification_code)

training_reference_validation <- tibble(
  metric = c(
    "Confirmed training HPO release",
    "Training HPO OBO file",
    "Primary terms in training HPO release",
    "Unique raw training HPO codes",
    "Unique cleaned training HPO codes",
    "Model feature codes",
    "Clean training codes match model features",
    "Training codes removed during preprocessing",
    "Removed codes that are parents of retained codes",
    "Removed codes not explained by parent trimming"
  ),
  value = c(
    training_hpo_release,
    training_hpo_path,
    length(training_hpo_ontology$id),
    length(training_raw_hpo_codes),
    length(training_clean_hpo_codes),
    nrow(model_term_lookup),
    setequal(training_clean_hpo_codes, model_term_lookup$hpo),
    length(training_raw_only_hpo_codes),
    length(training_parent_removed_hpo_codes),
    length(training_other_removed_hpo_codes)
  )
)

training_preprocessing_removed_terms <- tibble(
  training_hpo = sort(training_raw_only_hpo_codes)
) %>%
  left_join(training_hpo_term_table, by = "training_hpo") %>%
  mutate(
    removal_mechanism = if_else(
      training_hpo %in% training_parent_removed_hpo_codes,
      "Parent trimming",
      "Other preprocessing"
    )
  ) %>%
  left_join(training_retained_descendant_lookup, by = "training_hpo") %>%
  left_join(
    training_concept_mappings %>%
      filter(hpo %in% terms_absent_from_training$hpo) %>%
      group_by(training_hpo) %>%
      summarise(
        matching_current_cohort_hpo = collapse_character_values(hpo),
        .groups = "drop"
      ),
    by = "training_hpo"
  ) %>%
  arrange(removal_mechanism, training_hpo)

hpo_preprocessing_qc <- tibble(
  metric = c(
    "Reviewed patients",
    "Patients with at least one raw pre-review HPO term",
    "Raw patient-HPO pairs",
    "Obsolete input patient-HPO pairs",
    "Obsolete input pairs replaced",
    "Obsolete input pairs removed without replacement",
    "Patient-HPO pairs after obsolete-term resolution",
    "Modifier patient-HPO pairs removed",
    "Patient-HPO pairs after modifier removal",
    "Parent patient-HPO pairs removed",
    "Final model-retained patient-HPO pairs",
    "Raw distinct cohort HPO terms",
    "Final model-retained distinct cohort HPO terms",
    "Raw distinct cohort terms absent from training",
    "Final model-retained distinct cohort terms absent from training"
  ),
  value = c(
    nrow(reviewed_patients),
    n_distinct(raw_term_presence$patient_key),
    nrow(raw_term_presence),
    sum(term_resolution$is_obsolete),
    sum(
      term_resolution$is_obsolete &
        !is.na(term_resolution$replacement_hpo)
    ),
    sum(
      term_resolution$is_obsolete &
        is.na(term_resolution$replacement_hpo)
    ),
    nrow(resolved_term_presence),
    nrow(resolved_term_presence) - nrow(non_modifier_term_presence),
    nrow(non_modifier_term_presence),
    nrow(non_modifier_term_presence) - nrow(term_presence),
    nrow(term_presence),
    n_distinct(raw_term_presence$hpo),
    n_distinct(term_presence$hpo),
    nrow(raw_terms_absent_from_training),
    nrow(terms_absent_from_training)
  )
)


# 2. Term-yielding notes and unique HPO counts per patient ----------------

cohort_category_counts <- bind_rows(
  tibble(
    summary_group = "Overall",
    n = nrow(patient_phenotype_counts)
  ),
  patient_phenotype_counts %>%
    count(subgrp, name = "n") %>%
    transmute(summary_group = as.character(subgrp), n)
) %>%
  mutate(
    reviewed_cohort_n = nrow(patient_phenotype_counts),
    percent_of_reviewed_cohort = 100 * n / reviewed_cohort_n
  )

patient_count_summaries <- summarize_numeric_metrics(
  patient_phenotype_counts,
  group_variable = "subgrp",
  metrics = c(
    "term_yielding_note_count",
    "raw_unique_hpo_count",
    "model_retained_hpo_count",
    "hpo_terms_removed_by_preprocessing",
    "days_admission_to_review"
  )
)

positive_negative_patient_comparisons <- map_dfr(
  c(
    "term_yielding_note_count",
    "raw_unique_hpo_count",
    "model_retained_hpo_count",
    "hpo_terms_removed_by_preprocessing",
    "days_admission_to_review"
  ),
  ~ wilcoxon_group_comparison(patient_phenotype_counts, .x)
) %>%
  mutate(p_value_bh = p.adjust(p_value, method = "BH"))

diagnostic_patients <- patient_phenotype_counts %>%
  filter(subgrp %in% c("Positive", "Negative")) %>%
  mutate(result_positive = as.integer(subgrp == "Positive"))

positive_negative_adjusted_count_models <- map_dfr(
  c(
    "term_yielding_note_count",
    "raw_unique_hpo_count",
    "model_retained_hpo_count"
  ),
  ~ fit_adjusted_count_model(diagnostic_patients, .x)
) %>%
  mutate(p_value_bh = p.adjust(p_value, method = "BH"))


# 3. Term frequency comparison: Positive vs Negative ---------------------

diagnostic_term_presence <- term_presence %>%
  inner_join(
    diagnostic_patients %>% select(patient_key, subgrp),
    by = "patient_key"
  )

n_positive <- sum(diagnostic_patients$subgrp == "Positive")
n_negative <- sum(diagnostic_patients$subgrp == "Negative")

term_result_counts <- diagnostic_term_presence %>%
  group_by(hpo) %>%
  summarise(
    positive_with_term = sum(subgrp == "Positive"),
    negative_with_term = sum(subgrp == "Negative"),
    .groups = "drop"
  ) %>%
  mutate(
    positive_n = n_positive,
    negative_n = n_negative,
    positive_without_term = positive_n - positive_with_term,
    negative_without_term = negative_n - negative_with_term,
    positive_percent = 100 * positive_with_term / positive_n,
    negative_percent = 100 * negative_with_term / negative_n,
    prevalence_difference_percent = positive_percent - negative_percent
  )

term_patient_keys <- split(term_presence$patient_key, term_presence$hpo)

term_frequency_positive_vs_negative <- term_result_counts %>%
  mutate(
    fisher_result = pmap(
      list(
        positive_with_term,
        positive_without_term,
        negative_with_term,
        negative_without_term
      ),
      safe_fisher_test
    ),
    adjusted_result = map(
      hpo,
      ~ fit_term_result_model(
        .x,
        diagnostic_patients,
        term_patient_keys
      )
    )
  ) %>%
  unnest(fisher_result) %>%
  unnest(adjusted_result) %>%
  left_join(term_annotation_lookup, by = "hpo") %>%
  mutate(
    p_value_bh = p.adjust(p_value, method = "BH"),
    adjusted_p_value_bh = p.adjust(
      adjusted_p_value,
      method = "BH"
    )
  ) %>%
  arrange(p_value_bh, p_value, hpo)


# 4. HPO term associations with review-date MPSE percentile --------------

term_percentile_associations <- tibble(
  hpo = sort(unique(term_presence$hpo))
) %>%
  mutate(
    association_result = map(
      hpo,
      ~ fit_term_percentile_model(
        .x,
        patient_phenotype_counts,
        term_patient_keys
      )
    )
  ) %>%
  unnest(association_result) %>%
  left_join(term_annotation_lookup, by = "hpo") %>%
  mutate(
    p_value_bh = p.adjust(p_value, method = "BH"),
    adjusted_p_value_bh = p.adjust(
      adjusted_p_value,
      method = "BH"
    )
  ) %>%
  arrange(p_value_bh, p_value, hpo)


# Write analysis-ready tables --------------------------------------------

analysis_parameters <- tribble(
  ~parameter, ~value,
  "Phenotype cutoff", "load_date <= mpse_review_date",
  "Cutoff inclusivity", "Terms loaded on the review date are included",
  "Lower phenotype date bound", "None; MPSE phenotype data are cumulative",
  "HPO ontology release", hpo_ontology_release,
  "HPO ontology file", hpo_ontology_path,
  "Confirmed training HPO release", training_hpo_release,
  "Training HPO ontology file", training_hpo_path,
  "Raw training HPO code file", training_raw_hpo_codes_path,
  "Clean training HPO code file", training_clean_hpo_codes_path,
  "Training vocabulary validation", paste0(
    "The clean training code set must exactly equal the model feature ",
    "code set"
  ),
  "Training parent-trimming evidence", paste0(
    "A raw-only training code is classified as parent-trimmed when it is ",
    "an ancestor of at least one retained clean training code"
  ),
  "Training obsolete-term handling", paste0(
    "No obsolete-term filtering is assumed; classification follows the ",
    "observed raw and clean training code sets"
  ),
  "Cohort-only term classification policy", paste0(
    "Compare current IDs and alternative IDs with the confirmed 2018 HPO ",
    "release, then use raw and clean training-code membership"
  ),
  "HPO preprocessing", paste0(
    "BasicHPOSet-compatible: replace obsolete terms, remove modifier ",
    "branches, and retain only the most specific patient terms"
  ),
  "Modifier roots", paste(hpo_modifier_roots, collapse = "; "),
  "Obsolete terms without replacement", "Removed",
  "Term definition", paste0(
    "Present at least once per patient by clinical review and retained ",
    "after BasicHPOSet-compatible preprocessing"
  ),
  "Note definition", "Distinct note_id yielding at least one HPO term",
  "Analysis cohort", "All clinically reviewed patients",
  "Direct diagnostic comparison", "Positive vs Negative",
  "HPO hierarchy relation", "is_a",
  "Admission-to-review definition", "mpse_review_date - mpse_admit_date",
  "Minimum patients for regression", as.character(
    min_term_patients_for_model
  ),
  "Multiple-testing correction", "Benjamini-Hochberg false-discovery rate"
)

write_csv(
  analysis_parameters,
  file.path(analysis_output_dir, "00_analysis_parameters.csv"),
  na = ""
)
write_csv(
  cohort_category_counts,
  file.path(analysis_output_dir, "00_reviewed_cohort_counts.csv"),
  na = ""
)
write_csv(
  hpo_preprocessing_qc,
  file.path(analysis_output_dir, "00_hpo_preprocessing_qc.csv"),
  na = ""
)
write_csv(
  training_reference_validation,
  file.path(
    analysis_output_dir,
    "01a_training_reference_validation.csv"
  ),
  na = ""
)
write_csv(
  terms_absent_from_training,
  file.path(analysis_output_dir, "01_terms_absent_from_training.csv"),
  na = ""
)
write_csv(
  term_classification_summary,
  file.path(
    analysis_output_dir,
    "01b_term_classification_summary.csv"
  ),
  na = ""
)
write_csv(
  training_preprocessing_removed_terms,
  file.path(
    analysis_output_dir,
    "01c_training_preprocessing_removed_terms.csv"
  ),
  na = ""
)
write_csv(
  patient_count_summaries,
  file.path(analysis_output_dir, "02a_patient_counts_by_category.csv"),
  na = ""
)
write_csv(
  positive_negative_patient_comparisons,
  file.path(
    analysis_output_dir,
    "02b_positive_vs_negative_patient_comparisons.csv"
  ),
  na = ""
)
write_csv(
  positive_negative_adjusted_count_models,
  file.path(
    analysis_output_dir,
    "02c_positive_vs_negative_adjusted_count_models.csv"
  ),
  na = ""
)
write_csv(
  term_frequency_positive_vs_negative,
  file.path(
    analysis_output_dir,
    "03_term_frequencies_positive_vs_negative.csv"
  ),
  na = ""
)
write_csv(
  term_percentile_associations,
  file.path(
    analysis_output_dir,
    "04_term_associations_with_mpse_percentile.csv"
  ),
  na = ""
)


message("Phenotype analysis complete.")
message("Reviewed patients: ", nrow(reviewed_patients))
message(
  "Raw pre-review HPO terms: ",
  n_distinct(raw_term_presence$hpo)
)
message(
  "Model-retained pre-review HPO terms: ",
  n_distinct(term_presence$hpo)
)
message(
  "Model-retained terms absent from model training: ",
  nrow(terms_absent_from_training)
)
message(
  "Confirmed training HPO release: ",
  training_hpo_release
)
walk2(
  term_classification_summary$classification_code,
  term_classification_summary$term_count,
  ~ message("  ", .x, ": ", .y)
)
message("Tables written to: ", normalizePath(analysis_output_dir))
