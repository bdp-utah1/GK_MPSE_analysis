library(tidyverse)
library(gt)
library(gtsummary)

screened_cohort_manifest <- read_csv("/home/bennet/Documents/phd/yandell/gene_kids/data/GK_MPSE_analysis_copy/output/data/screened_cohort_analysis.csv") %>% 
  filter(!reviewed) %>% 
  select(sex, age_at_admission_days)
reviewed_cohort_manifest <- read_csv("/home/bennet/Documents/phd/yandell/gene_kids/data/GK_MPSE_analysis_copy/output/data/reviewed_cohort_analysis.csv") %>% 
  rename(sex=proband_gender)

cohort <- bind_rows(screened_cohort_manifest, reviewed_cohort_manifest) %>% 
  mutate(is_reviewed = !is.na(project_id),
         is_genekids = replace_na(gk_enrolled, FALSE),
         is_clinical = (!is_genekids & !is.na(patient_result))) %>% 
  select(sex, age_at_admission_days, days_admission_to_review, is_reviewed, is_clinical, is_genekids)





make_table1 <- function(data, header_name) {
  
  n <- nrow(data)
  
  data %>%
    select(
      age_at_admission_days,
      sex,
      days_admission_to_review
      #time_review_enroll
    ) %>%
    tbl_summary(
      statistic = list(
        all_continuous() ~ "{median} ({p25}, {p75})",
        all_categorical() ~ "{n} ({p}%)"
      ),
      digits = list(
        all_continuous() ~ 0,
        all_categorical() ~ c(0, 1)
      ),
      missing = "no",
      label = list(
        age_at_admission_days ~ "Age at admission, days",
        sex ~ "Sex",
        days_admission_to_review ~ "Time from admission to MPSE review, days"
        #time_review_enroll ~ "Time from MPSE review to GK enrollment, days"
      )
    ) %>%
    
    # Remove the redundant "Sex" parent row.
    # Keep the Female and Male rows.
    modify_table_body(
      ~ .x %>%
        filter(!(variable == "sex" & row_type == "label"))
    ) %>%
    
    # Put cohort name and N in the column header
    modify_header(
      stat_0 ~ paste0("**", header_name, "**<br>N = ", n)
    )
}

tbl_all <- make_table1(
  cohort,
  "All patients"
)

tbl_reviewed <- make_table1(
  cohort %>% filter(is_reviewed),
  "Reviewed patients"
)

tbl_clinical <- make_table1(
  cohort %>% filter(is_clinical),
  "Clinical testing patients"
)

tbl_genekids <- make_table1(
  cohort %>% filter(is_genekids),
  "GK testing patients"
)

tbl_all <- tbl_all %>%
  modify_table_body(
    ~ .x %>%
      mutate(
        stat_0 = case_when(
          variable == "days_admission_to_review" & row_type == "label"
          ~ "—",
          TRUE
          ~ stat_0
        )
      )
  )

# tbl_all <- tbl_all %>%
#   modify_table_body(
#     ~ .x %>%
#       mutate(
#         stat_0 = case_when(
#           variable == "time_review_enroll" & row_type == "label"
#           ~ "—",
#           TRUE ~ stat_0
#         )
#       )
#   )

tbl_reviewed <- tbl_reviewed %>%
  modify_table_body(
    ~ .x %>%
      mutate(
        stat_0 = case_when(
          variable == "time_review_enroll" & row_type == "label"
          ~ "—",
          TRUE ~ stat_0
        )
      )
  )

# tbl_clinical <- tbl_clinical %>%
#   modify_table_body(
#     ~ .x %>%
#       mutate(
#         stat_0 = case_when(
#           variable == "time_review_enroll" & row_type == "label"
#           ~ "—",
#           TRUE ~ stat_0
#         )
#       )
#   )

table1 <- tbl_merge(
  tbls = list(
    tbl_all,
    tbl_reviewed,
    tbl_clinical,
    tbl_genekids
  ),
  tab_spanner = FALSE
)

table1_gt <- table1 %>%
  as_gt() %>%
  tab_header(
    title = md("**Table 1. Patient cohort characteristics**")
  ) %>%
  tab_options(
    table.width = px(450),
    table.font.size = px(12),
    data_row.padding = px(4),
    column_labels.padding = px(4)
  )


gtsave(
  table1_gt,
  filename = "/home/bennet/Documents/phd/yandell/gene_kids/data/GK_MPSE_analysis_copy/output/figures/cohort_characteristics.png",
  zoom = 3,
  expand = 5
)
