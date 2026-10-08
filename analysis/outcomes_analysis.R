library(tidyverse)
library(ggsignif)


mosc_ptnt_data <- read_csv("/home/bennet/Documents/phd/yandell/gene_kids/data/GK_MPSE_analysis_copy/mosaic_data/patient_level_data.csv") %>% 
  select(mrn, dob, mpse_admit_date, mpse_review_date, mpse_cp_score, mpse_cp_percent, gk_enrollment_date, test_results_record_count) %>% 
  rename(sex = proband_gender, admit_date = mpse_admit_date)

mosc_test_data <- read_csv("/home/bennet/Documents/phd/yandell/gene_kids/data/GK_MPSE_analysis_copy/mosaic_data/test_results_data.csv")
admission_matching <- read_csv("/home/bennet/Documents/phd/yandell/gene_kids/data/GK_MPSE_analysis_copy/output/data/admission_match_manifest.csv")
  

mpse_patnt <- read_tsv("/home/bennet/Documents/phd/yandell/gene_kids/data/GK_MPSE_analysis_copy/mpse_data/filtered_patient.tsv")
mpse_visit <- read_tsv("/home/bennet/Documents/phd/yandell/gene_kids/data/GK_MPSE_analysis_copy/mpse_data/filtered_patient_stay.tsv") %>% 
  arrange(patient_id, arrive_dtm) %>% 
  distinct(patient_id, .keep_all = TRUE)




reviewed_cohort_manifest <- read_csv("/home/bennet/Documents/phd/yandell/gene_kids/data/GK_MPSE_analysis_copy/output/data/reviewed_cohort_analysis.csv")

df <- reviewed_cohort_manifest %>% 
  select(project_id, recorded_review_score, reconstructed_review_score, score_difference, 
         recorded_review_percentile, review_percentile_group, 
         gk_enrolled, has_any_test, primary_test_type, has_classifiable_completed_result, patient_result)


# MPSE score association with Pos vs Neg Dx
pos_vs_neg_logit_data <- df %>% 
  select(recorded_review_score, reconstructed_review_score, recorded_review_percentile, patient_result) %>% 
  drop_na() %>% 
  filter(patient_result == "Positive" | patient_result == "Negative") %>% 
  mutate(outcome = if_else(patient_result=="Positive", 1, 0))

summary(glm(outcome ~ recorded_review_percentile, data = pos_vs_neg_logit_data, family = binomial))
t.test(recorded_review_percentile ~ outcome, data=pos_vs_neg_logit_data)



# MPSE score association with Pos vs Unc Dx
pos_vs_unc_logit_data <- df %>% 
  select(recorded_review_score, reconstructed_review_score, recorded_review_percentile, patient_result) %>% 
  drop_na() %>% 
  filter(patient_result == "Positive" | patient_result == "Uncertain") %>% 
  mutate(outcome = if_else(patient_result=="Positive", 1, 0))

summary(glm(outcome ~ recorded_review_percentile, data = pos_vs_unc_logit_data, family = binomial))
t.test(recorded_review_percentile ~ outcome, data=pos_vs_unc_logit_data)



# MPSE score association with Pos vs Inc Dx
pos_vs_inc_logit_data <- df %>% 
  mutate(patient_result = if_else((has_any_test & is.na(patient_result)), "Incomplete", patient_result)) %>% 
  select(recorded_review_score, reconstructed_review_score, recorded_review_percentile, patient_result) %>% 
  drop_na() %>% 
  filter(patient_result == "Positive" | patient_result == "Incomplete") %>% 
  mutate(outcome = if_else(patient_result=="Positive", 1, 0))

summary(glm(outcome ~ recorded_review_percentile, data = pos_vs_inc_logit_data, family = binomial))
t.test(recorded_review_percentile ~ outcome, data=pos_vs_inc_logit_data)



# MPSE score association with Pos vs Neg Dx (WGS only)
pos_vs_neg_wgs_logit_data <- df %>% 
  select(recorded_review_score, reconstructed_review_score, recorded_review_percentile, primary_test_type, patient_result) %>% 
  drop_na() %>% 
  filter(primary_test_type == "WGS") %>% 
  filter(patient_result == "Positive" | patient_result == "Negative") %>% 
  mutate(outcome = if_else(patient_result=="Positive", 1, 0))

summary(glm(outcome ~ recorded_review_percentile, data = pos_vs_neg_wgs_logit_data, family = binomial))
t.test(recorded_review_percentile ~ outcome, data=pos_vs_neg_wgs_logit_data)



# MPSE score association with Pos vs Unc Dx (WGS only)
pos_vs_unc_wgs_logit_data <- df %>% 
  select(recorded_review_score, reconstructed_review_score, recorded_review_percentile, primary_test_type, patient_result) %>% 
  drop_na() %>% 
  filter(primary_test_type == "WGS") %>% 
  filter(patient_result == "Positive" | patient_result == "Uncertain") %>% 
  mutate(outcome = if_else(patient_result=="Positive", 1, 0))

summary(glm(outcome ~ recorded_review_percentile, data = pos_vs_unc_wgs_logit_data, family = binomial))
t.test(recorded_review_percentile ~ outcome, data=pos_vs_unc_wgs_logit_data)



# MPSE score association with Pos vs Inc Dx (WGS only)
pos_vs_inc_wgs_logit_data <- df %>% 
  mutate(patient_result = if_else((has_any_test & is.na(patient_result)), "Incomplete", patient_result)) %>% 
  select(recorded_review_score, reconstructed_review_score, recorded_review_percentile, primary_test_type, patient_result) %>% 
  drop_na() %>% 
  filter(primary_test_type == "WGS") %>% 
  filter(patient_result == "Positive" | patient_result == "Incomplete") %>% 
  mutate(outcome = if_else(patient_result=="Positive", 1, 0))

summary(glm(outcome ~ recorded_review_percentile, data = pos_vs_inc_wgs_logit_data, family = binomial))
t.test(recorded_review_percentile ~ outcome, data=pos_vs_inc_wgs_logit_data)



all_tested_percentile <- reviewed_cohort_manifest %>% 
  mutate(patient_result = if_else((has_any_test & is.na(patient_result)), "Incomplete", patient_result)) %>% 
  drop_na(patient_result) %>% 
  select(recorded_review_percentile, patient_result) %>% 
  mutate(subset="All tested")

wgs_only_percentile <- reviewed_cohort_manifest %>% 
  filter(primary_test_type == "WGS") %>% 
  mutate(patient_result = if_else((has_any_test & is.na(patient_result)), "Incomplete", patient_result)) %>% 
  drop_na(patient_result) %>% 
  select(recorded_review_percentile, patient_result) %>% 
  mutate(subset="WGS only")

untested_percentile <- reviewed_cohort_manifest %>% 
  filter(!has_any_test) %>% 
  select(recorded_review_percentile, patient_result) %>% 
  mutate(patient_result = "Untested",
         subset="Untested")

grouped_percentiles_data <- bind_rows(all_tested_percentile, wgs_only_percentile, untested_percentile) %>%
  mutate(
    subset = factor(
      subset,
      levels = c("All tested", "WGS only", "Untested")
    ),
    patient_result = factor(
      patient_result,
      levels = c("Positive", "Uncertain", "Negative", "Incomplete", "Untested")
    )
  ) %>%
  group_by(subset, patient_result) %>%
  mutate(
    x_label = paste0(
      patient_result,
      "\n(n = ", n(), ")"
    )
  ) %>%
  ungroup() %>%
  mutate(
    x_label = factor(
      x_label,
      levels = unique(x_label[
        order(patient_result)
      ])
    )
  )

# Sample sizes for x-axis labels
group_counts <- grouped_percentiles_data %>%
  count(subset, patient_result, .drop = FALSE) %>%
  filter(n > 0)

# Custom color palette
result_colors <- c(
  "Positive"   = "#D55E00",
  "Uncertain"  = "#009E73",
  "Negative"   = "#0072B2",
  "Incomplete" = "#8E44AD",
  "Untested"   = "#999999"
)

significance_data <- tibble(
  subset = factor(
    c("All tested", "WGS only", "All tested", "WGS only", "All tested", "WGS only"),
    levels = c("All tested", "WGS only", "Untested")
  ),
  xmin = c(1, 1, 1, 1, 1, 1),
  xmax = c(2, 2, 3, 3, 4, 4),
  y_position = c(102, 102, 106, 106, 110, 110),
  annotations = c("p=0.15", "p=0.29", "p=0.0047", "p=0.0036", "p=0.022", "p=0.076")
)

grouped_percentiles_plot <- grouped_percentiles_data %>% 
  ggplot(
    aes(
      x = x_label,
      y = recorded_review_percentile,
      fill = patient_result,
      color = patient_result
    )
  ) +
  geom_boxplot(
    width = 0.55,
    alpha = 0.08,
    outlier.shape = NA,
    linewidth = 0.7
  ) +
  geom_jitter(
    width = 0.15,
    height = 0,
    size = 1.9,
    alpha = 0.3
  ) +
  geom_signif(
    data = significance_data,
    aes(
      xmin = xmin,
      xmax = xmax,
      annotations = annotations,
      y_position = y_position
    ),
    manual = TRUE,
    inherit.aes = FALSE,
    textsize = 5,
    tip_length = 0,
    linewidth = 0.6
  ) +
  facet_grid(
    ~ subset,
    scales = "free_x",
    space = "free_x"
  ) +
  scale_fill_manual(values = result_colors) +
  scale_color_manual(values = result_colors) +
  scale_y_continuous(
    breaks = seq(0, 100, 10),
    expand = expansion(mult = c(0, 0.02))
  ) +
  coord_cartesian(ylim = c(0, 100), clip = "off") +
  labs(
    x = "Diagnostic category",
    y = "MPSE score percentile"
  ) +
  theme_classic(base_size = 16) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 17,
      margin = margin(b = 80)
    ),
    axis.text.x = element_text(size = 16),
    axis.title.x = element_text(size = 17),
    axis.text.y = element_text(size = 16),
    axis.title.y = element_text(size = 17),
    panel.spacing.x = unit(1.4, "lines"),
    legend.position = "none"
  )


ggsave(
  filename = "/home/bennet/Documents/phd/yandell/gene_kids/data/GK_MPSE_analysis/output/figures/grouped_percentile_boxplots.svg",
  plot = grouped_percentiles_plot,
  width = 12,
  height = 11,
  units = "in",
  bg = "white"
)
