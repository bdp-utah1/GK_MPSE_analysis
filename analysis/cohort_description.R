library(tidyverse)
library(gt)
library(gtsummary)


mosc_data <- read_csv("/home/bennet/Documents/phd/yandell/gene_kids/data/GK_MPSE_analysis_copy/mosaic_data/patient_level_data.csv") %>% 
  select(mrn, dob, proband_gender, mpse_admit_date, mpse_review_date, gk_enrollment_date, test_results_record_count) %>% 
  rename(sex = proband_gender, admit_date = mpse_admit_date)


mpse_patnt <- read_tsv("/home/bennet/Documents/phd/yandell/gene_kids/data/GK_MPSE_analysis_copy/mpse_data/filtered_patient.tsv")
mpse_visit <- read_tsv("/home/bennet/Documents/phd/yandell/gene_kids/data/GK_MPSE_analysis_copy/mpse_data/filtered_patient_stay.tsv") %>% 
  arrange(patient_id, arrive_dtm) %>% 
  distinct(patient_id, .keep_all = TRUE)

mpse_data <- inner_join(mpse_patnt, mpse_visit, by=join_by(patient_id)) %>% 
  select(patient_id, dob, sex, arrive_dtm) %>% 
  mutate(arrive_dtm = as.Date(arrive_dtm)) %>% 
  rename(mrn = patient_id, admit_date = arrive_dtm) %>% 
  filter(!mrn %in% mosc_data$mrn)

cohort <- bind_rows(mosc_data, mpse_data) %>% 
  mutate(age_at_admit = as.integer(admit_date - dob),
         time_admit_review = as.integer(mpse_review_date - admit_date),
         time_review_enroll = as.integer(gk_enrollment_date - mpse_review_date)) %>% 
  mutate(is_reviewed = !is.na(mpse_review_date),
         is_genekids = !is.na(gk_enrollment_date),
         is_clinical = (!is_genekids & test_results_record_count != 0)) %>% 
  select(sex, age_at_admit, time_admit_review, time_review_enroll, is_reviewed, is_clinical, is_genekids)





make_table1 <- function(data, header_name) {
  
  n <- nrow(data)
  
  data %>%
    select(
      age_at_admit,
      sex,
      time_admit_review,
      time_review_enroll
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
        age_at_admit ~ "Age at admission, days",
        sex ~ "Sex",
        time_admit_review ~ "Time from admission to MPSE review, days",
        time_review_enroll ~ "Time from MPSE review to GK enrollment, days"
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
          variable == "time_admit_review" & row_type == "label"
          ~ "—",
          TRUE
          ~ stat_0
        )
      )
  )

tbl_all <- tbl_all %>%
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

tbl_clinical <- tbl_clinical %>%
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
  )
