library(tidyverse)


defunct_mrn_4879 <- "X014455166"
current_mrn_4879 <- "X014453675"

cern_unmerged_mrn_4048 <- "650480153"
epic_unmerged_mrn_4048 <- "X014183906"



mosc_patnt <- read_csv("/home/bennet/Documents/phd/yandell/gene_kids/data/GK_MPSE_analysis_copy/mosaic_data/patient_level_data.csv") %>% 
  select(project_id, mrn, dob, mpse_admit_date, mpse_review_date)

mosc_outcm <- read_csv("/home/bennet/Documents/phd/yandell/gene_kids/data/GK_MPSE_analysis_copy/output/data/reviewed_cohort_analysis.csv") %>% 
  select(project_id, gk_enrolled, has_any_test, has_classifiable_completed_result, patient_result)

mosc_data <- inner_join(mosc_patnt, mosc_outcm, by="project_id")

mpse_patnt <- read_tsv("/home/bennet/Documents/phd/yandell/gene_kids/data/GK_MPSE_analysis_copy/mpse_data/filtered_patient.tsv")
mpse_score <- read_tsv("/home/bennet/Documents/phd/yandell/gene_kids/data/GK_MPSE_analysis_copy/mpse_data/filtered_patient_mpse_score.tsv")

mpse_data <- inner_join(mpse_patnt, mpse_score, by="patient_id") %>% 
  mutate(patient_id = str_replace(patient_id, defunct_mrn_4879, current_mrn_4879),
         patient_id = str_replace(patient_id, cern_unmerged_mrn_4048, epic_unmerged_mrn_4048)) %>% 
  rename(mrn=patient_id) %>% 
  select(-days_seen)

df <- full_join(mosc_data, mpse_data, by=join_by(mrn, dob)) %>%
  arrange(mrn, date) %>% 
  group_by(mrn) %>% 
  complete(date = seq(min(date), max(date), by = "day")) %>%
  fill(everything(), .direction = "down") %>% 
  mutate(manifest_rank = rank(date)) %>% 
  ungroup()




score_trajectories <- df %>% 
  filter(manifest_rank <= 30) %>% 
  mutate(subgrp = factor(case_when(
    is.na(project_id) ~ "Unreviewed",
    project_id & !has_any_test ~ "Reviewed & Untested",
    has_classifiable_completed_result ~ patient_result
  ), levels=c("Positive", "Uncertain", "Negative", "Reviewed & Untested", "Unreviewed"))) %>% 
  ggplot(aes(x=manifest_rank, y=mpse_score, group=mrn)) + 
  geom_line(color="#2C3E50", alpha=0.5) +
  geom_smooth(inherit.aes = FALSE, aes(x=manifest_rank, y=mpse_score, group=subgrp), method="lm", formula=y~poly(x, 8), se=FALSE, color="#E74C3C", linewidth=1.2) +
  facet_wrap(~subgrp, nrow=1) + 
  scale_x_continuous(name="Days After Admission",
                     expand=c(0,0),
                     limits=c(0,30),
                     minor_breaks=FALSE) +
  scale_y_continuous(name="MPSE Score",
                     breaks=seq(-50,200,50),
                     labels=seq(-50,200,50),
                     minor_breaks=FALSE) + 
  theme_bw() +
  theme(panel.spacing = unit(0.7, "cm"),
        panel.grid.major = element_line(color = "grey70", linewidth = 0.4))

ggsave(
  filename = "/home/bennet/Documents/phd/yandell/gene_kids/data/GK_MPSE_analysis/output/figures/grouped_score_trajectories.svg",
  plot = score_trajectories,
  width = 8,
  height = 4,
  units = "in",
  bg = "white"
)


df %>% 
  #filter(is.na(mpse_review_date) | (date <= mpse_review_date)) %>% 
  #filter(!is.na(mpse_review_date) | (manifest_rank <= 30)) %>% 
  filter(manifest_rank <= 30) %>% 
  mutate(subgrp = factor(case_when(
    is.na(project_id) ~ "Unreviewed",
    project_id & !has_any_test ~ "Reviewed & Untested",
    has_any_test & !has_classifiable_completed_result ~ "Tested & Incomplete",
    has_classifiable_completed_result ~ patient_result
  ), levels=c("Positive", "Uncertain", "Negative", "Reviewed & Untested", "Unreviewed"))) %>% 
  ggplot(aes(x=manifest_rank, y=mpse_score, group=subgrp, color=subgrp)) + 
  geom_smooth()
