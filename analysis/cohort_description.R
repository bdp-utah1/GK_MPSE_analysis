library(tidyverse)
library(gt)
library(gtsummary)


# File paths ------------------------------------------------------------------

project_dir <- "/home/bennet/Documents/phd/yandell/gene_kids/data/GK_MPSE_analysis"

cohort_file <- file.path(
  project_dir,
  "output/data/master_cohort_analysis.csv"
)

# Set GK_MPSE_FIGURE_DIR to write figures somewhere other than the project
# output directory (useful when testing the script).
figure_dir <- Sys.getenv(
  "GK_MPSE_FIGURE_DIR",
  unset = file.path(project_dir, "output/figures")
)

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)


# Editable demographic mappings ----------------------------------------------

# The names on the left are values in the source data. The values on the right
# are the categories displayed in Table 1. Edit the right-hand side to change
# how a source category is grouped.
race_map <- c(
  "White" = "White",
  "White or Caucasian" = "White",
  "Black or African American" = "Black or African American",
  "Asian" = "Asian",
  "Filipino American" = "Asian",
  "American Indian or Alaska Native" = "American Indian or Alaska Native",
  "Alaskan Native American" = "American Indian or Alaska Native",
  "Native American" = "American Indian or Alaska Native",
  "Other Pacific Islander" = "Native Hawaiian or Other Pacific Islander",
  "Multiple" = "Other or Multiple",
  "Other" = "Other or Multiple",
  "Unknown" = "Unknown or declined",
  "Patient Declined" = "Unknown or declined"
)

ethnicity_map <- c(
  "Non-Hispanic" = "Non-Hispanic",
  "Not Hispanic, Latino, or Spanish Origin" = "Non-Hispanic",
  "Hispanic" = "Hispanic",
  "Hispanic, Latino, or Spanish Origin" = "Hispanic",
  "Mexican, Mexican American, or Chicano/a" = "Hispanic",
  "Other Hispanic, Latino/a, or Spanish origin" = "Hispanic",
  "Puerto Rican" = "Hispanic",
  "Unknown" = "Unknown or declined",
  "Patient Declined" = "Unknown or declined"
)


# Prepare analysis data -------------------------------------------------------

cohort_manifest <- read_csv(cohort_file, show_col_types = FALSE)

unmapped_race <- setdiff(
  unique(na.omit(cohort_manifest$race)),
  names(race_map)
)
unmapped_ethnicity <- setdiff(
  unique(na.omit(cohort_manifest$ethnicity)),
  names(ethnicity_map)
)

if (length(unmapped_race) > 0) {
  stop(
    "Add the following source values to race_map: ",
    paste(unmapped_race, collapse = ", ")
  )
}

if (length(unmapped_ethnicity) > 0) {
  stop(
    "Add the following source values to ethnicity_map: ",
    paste(unmapped_ethnicity, collapse = ", ")
  )
}

cohort <- cohort_manifest %>%
  mutate(
    is_reviewed = !is.na(project_id),

    # Testing status is only defined for reviewed patients. Within reviewed
    # patients, this is equivalent to has_any_test in the current dataset.
    is_tested = case_when(
      !is_reviewed ~ NA,
      TRUE ~ !is.na(patient_result)
    ),

    age_at_admission_days = as.integer(selected_arrive_date - dob),
    days_admission_to_review = as.integer(
      reviewed_raw_mpse_review_date - selected_arrive_date
    ),

    # Convert strings such as "39 3/7" to decimal weeks.
    gestational_age_weeks =
      as.numeric(str_extract(gest_age, "^\\d+")) +
      as.numeric(str_match(gest_age, "(\\d+)/7")[, 2]) / 7,

    gender = factor(gender, levels = c("Female", "Male")),
    race = recode(
      race,
      !!!race_map,
      .default = NA_character_,
      .missing = "Unknown or declined"
    ) %>%
      factor(levels = c(
        "White",
        "Black or African American",
        "Asian",
        "American Indian or Alaska Native",
        "Native Hawaiian or Other Pacific Islander",
        "Other or Multiple",
        "Unknown or declined"
      )),
    ethnicity = recode(
      ethnicity,
      !!!ethnicity_map,
      .default = NA_character_,
      .missing = "Unknown or declined"
    ) %>%
      factor(levels = c(
        "Non-Hispanic",
        "Hispanic",
        "Unknown or declined"
      )),
    coi_natrnk = factor(
      coi_natrnk,
      levels = c("Very Low", "Low", "Moderate", "High", "Very High", "Unknown")
    ),
    review_status = factor(
      if_else(is_reviewed, "Reviewed", "Unreviewed"),
      levels = c("Reviewed", "Unreviewed")
    ),
    testing_status = factor(
      case_when(
        is_tested %in% TRUE ~ "Tested",
        is_tested %in% FALSE ~ "Untested",
        TRUE ~ NA_character_
      ),
      levels = c("Tested", "Untested")
    ),
    diagnosis = factor(
      patient_result,
      levels = c("Positive", "Uncertain", "Negative")
    )
  ) %>%
  select(
    gender,
    race,
    ethnicity,
    gestational_age_weeks,
    age_at_admission_days,
    coi_natrnk,
    days_admission_to_review,
    is_reviewed,
    is_tested,
    has_any_test,
    review_status,
    testing_status,
    diagnosis
  )

if (!identical(cohort$is_tested, cohort$has_any_test)) {
  stop("is_tested and has_any_test are not equivalent in the current dataset.")
}

cohort <- cohort %>%
  select(-has_any_test)


# Shared table settings -------------------------------------------------------

table_variables <- c(
  "gender",
  "race",
  "ethnicity",
  "gestational_age_weeks",
  "age_at_admission_days",
  "coi_natrnk",
  "days_admission_to_review"
)

table_labels <- list(
  gender ~ "Gender",
  race ~ "Race",
  ethnicity ~ "Ethnicity",
  gestational_age_weeks ~ "Gestational age at birth, weeks",
  age_at_admission_days ~ "Age at admission, days",
  coi_natrnk ~ "COI national rank",
  days_admission_to_review ~ "Time from admission to MPSE review, days"
)

summary_statistics <- list(
  all_continuous() ~ "{median} ({p25}, {p75})",
  all_categorical() ~ "{n} ({p}%)"
)

summary_digits <- list(
  gestational_age_weeks ~ 1,
  age_at_admission_days ~ 0,
  days_admission_to_review ~ 0,
  all_categorical() ~ c(0, 1)
)

# Fisher's exact tests use Monte Carlo simulation because several demographic
# categories are sparse. Setting the seed makes the reported p-values
# reproducible when the input data are unchanged.
set.seed(20261124)

pvalue_format <- function(x) {
  ifelse(is.na(x), "\u2014", style_pvalue(x, digits = 3))
}

categorical_test_args <- list(simulate.p.value = TRUE, B = 100000)

table_footnote <- paste(
  "Continuous variables were compared using the Wilcoxon rank-sum test",
  "(two groups) or Kruskal-Wallis test (three groups).",
  "Categorical variables were compared using Fisher's exact test with",
  "Monte Carlo simulated p-values (100,000 replicates).",
  "Unknown and declined responses are shown explicitly.",
  "In Panel A, the full-cohort time-to-review estimate includes reviewed",
  "patients only; this measure is undefined for unreviewed patients."
)


style_table <- function(table, subtitle) {
  table %>%
    as_gt() %>%
    tab_header(
      title = md("**Table 1. Patient cohort characteristics**"),
      subtitle = md(subtitle)
    ) %>%
    tab_source_note(source_note = table_footnote) %>%
    tab_options(
      table.width = pct(100),
      table.font.size = px(12),
      data_row.padding = px(4),
      column_labels.padding = px(5),
      heading.align = "left"
    )
}


# Panel A: full cohort and review status -------------------------------------

# The overall column is descriptive only. P-values compare Reviewed with
# Unreviewed. Time to review is undefined for unreviewed patients and is not
# tested in this panel.
table1a <- cohort %>%
  select(all_of(table_variables), review_status) %>%
  tbl_summary(
    by = review_status,
    statistic = summary_statistics,
    digits = summary_digits,
    missing = "no",
    label = table_labels
  ) %>%
  add_overall(
    last = FALSE,
    col_label = "**Full cohort**  \nN = {N}"
  ) %>%
  add_p(
    include = -days_admission_to_review,
    test = list(
      all_continuous() ~ "wilcox.test",
      all_categorical() ~ "fisher.test"
    ),
    test.args = all_categorical() ~ categorical_test_args,
    pvalue_fun = pvalue_format
  ) %>%
  modify_table_body(
    ~ .x %>%
      mutate(
        stat_2 = if_else(
          variable == "days_admission_to_review" & row_type == "label",
          "\u2014",
          stat_2
        ),
        p.value = if_else(
          variable == "days_admission_to_review" & row_type == "label",
          -1,
          p.value
        )
      )
  ) %>%
  modify_header(
    label = "**Characteristic**",
    p.value = "**p-value**"
  ) %>%
  remove_footnote_header(columns = p.value) %>%
  modify_fmt_fun(
    p.value ~ function(x) rep("\u2014", length(x)),
    rows = variable == "days_admission_to_review" & row_type == "label"
  ) %>%
  bold_labels()

table1a_gt <- style_table(
  table1a,
  "**A. Full cohort by MPSE review status**"
)


# Panel B: reviewed cohort and testing status --------------------------------

# The overall column is descriptive only. P-values compare Tested with Untested
# among reviewed patients.
table1b <- cohort %>%
  filter(is_reviewed) %>%
  select(all_of(table_variables), testing_status) %>%
  tbl_summary(
    by = testing_status,
    statistic = summary_statistics,
    digits = summary_digits,
    missing = "no",
    label = table_labels
  ) %>%
  add_overall(
    last = FALSE,
    col_label = "**Reviewed patients**  \nN = {N}"
  ) %>%
  add_p(
    test = list(
      all_continuous() ~ "wilcox.test",
      all_categorical() ~ "fisher.test"
    ),
    test.args = all_categorical() ~ categorical_test_args,
    pvalue_fun = pvalue_format
  ) %>%
  modify_header(
    label = "**Characteristic**",
    p.value = "**p-value**"
  ) %>%
  remove_footnote_header(columns = p.value) %>%
  bold_labels()

table1b_gt <- style_table(
  table1b,
  "**B. Reviewed cohort by testing status**"
)


# Panel C: tested cohort and diagnostic result -------------------------------

# P-values are overall three-group comparisons across Positive, Uncertain, and
# Negative diagnostic results.
table1c <- cohort %>%
  filter(is_reviewed, is_tested %in% TRUE) %>%
  select(all_of(table_variables), diagnosis) %>%
  tbl_summary(
    by = diagnosis,
    statistic = summary_statistics,
    digits = summary_digits,
    missing = "no",
    label = table_labels
  ) %>%
  add_p(
    test = list(
      all_continuous() ~ "kruskal.test",
      all_categorical() ~ "fisher.test"
    ),
    test.args = all_categorical() ~ categorical_test_args,
    pvalue_fun = pvalue_format
  ) %>%
  modify_header(
    label = "**Characteristic**",
    p.value = "**p-value**"
  ) %>%
  remove_footnote_header(columns = p.value) %>%
  bold_labels()

table1c_gt <- style_table(
  table1c,
  "**C. Tested cohort by diagnostic result**"
)


# Save output -----------------------------------------------------------------

output_tables <- list(
  cohort_characteristics_part_a_review_status = table1a_gt,
  cohort_characteristics_part_b_testing_status = table1b_gt,
  cohort_characteristics_part_c_diagnosis = table1c_gt
)

if (Sys.getenv("GK_MPSE_SKIP_SAVE", unset = "false") != "true") {
  iwalk(
    output_tables,
    ~ gtsave(
      .x,
      filename = file.path(figure_dir, paste0(.y, ".png")),
      zoom = 3,
      expand = 5
    )
  )
}
