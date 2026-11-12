library(tidyverse)


scoring_process_change <- as.Date("2025-05-16")
frst_review_date <- as.Date("2025-08-11")
last_review_date <- as.Date("2026-06-02")

defunct_mrn_4879 <- "X014455166"
current_mrn_4879 <- "X014453675"

cern_unmerged_mrn_4048 <- "650480153"
epic_unmerged_mrn_4048 <- "X014183906"

mrn_dob_4104 <- "X014208664"
bad_dob_4104 <- as.Date("2025-07-01")
fix_dob_4104 <- as.Date("2025-07-02")

cohort_ids <- read_csv("/home/bennet/Documents/phd/yandell/gene_kids/data/GK_MPSE_analysis_copy/mosaic_data/patient_level_data.csv") %>% 
  select(project_id, mrn, other_mrn, dob, proband_gender, mpse_admit_date)

study_period_total_screened <- read_tsv("/home/bennet/Documents/phd/yandell/gene_kids/data/GK_MPSE_analysis_copy/mpse_data/patient_stay.tsv") %>% 
  mutate(patient_id = str_replace(patient_id, defunct_mrn_4879, current_mrn_4879)) %>% 
  filter((unit %in% c("CPH 4E NEONATAL ICU", "PC_4E NICU")) | (patient_id %in% cohort_ids$mrn)) %>% 
  filter((dischg_dtm > frst_review_date) | (patient_id %in% cohort_ids$mrn)) %>% 
  filter((arrive_dtm < last_review_date) | (patient_id %in% cohort_ids$mrn)) %>% 
  filter((arrive_dtm >= scoring_process_change) | (patient_id %in% cohort_ids$mrn)) %>% 
  filter(patient_id != epic_unmerged_mrn_4048) %>% 
  distinct(patient_id, .keep_all = TRUE)



mpse_patnt <- read_tsv("/home/bennet/Documents/phd/yandell/gene_kids/data/GK_MPSE_analysis_copy/mpse_data/patient.tsv") %>% 
  mutate(patient_id = str_replace(patient_id, defunct_mrn_4879, current_mrn_4879)) %>% 
  mutate(dob = if_else((patient_id==mrn_dob_4104 & dob==bad_dob_4104), fix_dob_4104, dob)) %>% 
  select(patient_id, sex, dob, first_seen, last_seen) %>% 
  filter(patient_id %in% study_period_total_screened$patient_id) %>% 
  mutate(patient_id = str_replace(patient_id, cern_unmerged_mrn_4048, epic_unmerged_mrn_4048))
write_tsv(mpse_patnt, "/home/bennet/Documents/phd/yandell/gene_kids/data/GK_MPSE_analysis_copy/mpse_data/filtered_patient.tsv", na="")

mpse_visit <- read_tsv("/home/bennet/Documents/phd/yandell/gene_kids/data/GK_MPSE_analysis_copy/mpse_data/patient_stay.tsv") %>% 
  mutate(patient_id = str_replace(patient_id, defunct_mrn_4879, current_mrn_4879)) %>% 
  select(patient_id, unit, arrive_dtm, dischg_dtm) %>% 
  filter(patient_id %in% study_period_total_screened$patient_id) %>% 
  mutate(patient_id = str_replace(patient_id, cern_unmerged_mrn_4048, epic_unmerged_mrn_4048))
write_tsv(mpse_visit, "/home/bennet/Documents/phd/yandell/gene_kids/data/GK_MPSE_analysis_copy/mpse_data/filtered_patient_stay.tsv", na="")

mpse_score <- read_tsv("/home/bennet/Documents/phd/yandell/gene_kids/data/GK_MPSE_analysis_copy/mpse_data/patient_mpse_score.tsv") %>% 
  mutate(patient_id = str_replace(patient_id, defunct_mrn_4879, current_mrn_4879)) %>% 
  filter(hpo_source == "ClinPhen_bulk") %>% 
  select(patient_id, date, days_seen, mpse_score) %>% 
  filter(patient_id %in% study_period_total_screened$patient_id) %>% 
  mutate(patient_id = str_replace(patient_id, cern_unmerged_mrn_4048, epic_unmerged_mrn_4048))
write_tsv(mpse_score, "/home/bennet/Documents/phd/yandell/gene_kids/data/GK_MPSE_analysis_copy/mpse_data/filtered_patient_mpse_score.tsv", na="")

mpse_pheno <- read_tsv("/home/bennet/Documents/phd/yandell/gene_kids/data/GK_MPSE_analysis_copy/mpse_data/patient_hpo.tsv") %>% 
  mutate(patient_id = str_replace(patient_id, defunct_mrn_4879, current_mrn_4879)) %>% 
  filter(hpo_source == "ClinPhen_bulk") %>% 
  select(patient_id, load_date, note_id, note_type, note_dtm, hpo) %>% 
  filter(patient_id %in% study_period_total_screened$patient_id) %>% 
  mutate(patient_id = str_replace(patient_id, cern_unmerged_mrn_4048, epic_unmerged_mrn_4048))
write_tsv(mpse_pheno, "/home/bennet/Documents/phd/yandell/gene_kids/data/GK_MPSE_analysis_copy/mpse_data/filtered_patient_hpo.tsv", na="")




