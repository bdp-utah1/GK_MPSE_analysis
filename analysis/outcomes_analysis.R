library(tidyverse)
library(ggsignif)



cohort_manifest <- read_csv("/home/bennet/Documents/phd/yandell/gene_kids/data/GK_MPSE_analysis_copy/output/data/master_cohort_analysis.csv")


# MPSE score association with Pos vs Neg Dx
pos_vs_neg_logit_data <- cohort_manifest %>% 
  select(recorded_review_percentile, patient_result) %>% 
  drop_na() %>% 
  filter(patient_result == "Positive" | patient_result == "Negative") %>% 
  mutate(outcome = if_else(patient_result=="Positive", 1, 0))

summary(glm(outcome ~ recorded_review_percentile, data = pos_vs_neg_logit_data, family = binomial))
t.test(recorded_review_percentile ~ outcome, data=pos_vs_neg_logit_data)


# MPSE score association with Pos vs Unc Dx
pos_vs_unc_logit_data <- cohort_manifest %>% 
  select(recorded_review_score, reconstructed_review_score, recorded_review_percentile, patient_result) %>% 
  drop_na() %>% 
  filter(patient_result == "Positive" | patient_result == "Uncertain") %>% 
  mutate(outcome = if_else(patient_result=="Positive", 1, 0))

summary(glm(outcome ~ recorded_review_percentile, data = pos_vs_unc_logit_data, family = binomial))
t.test(recorded_review_percentile ~ outcome, data=pos_vs_unc_logit_data)






tested_percentile <- cohort_manifest %>% 
  drop_na(patient_result) %>% 
  select(recorded_review_percentile, patient_result) %>% 
  mutate(subset="Tested")

# Map the original Gene Kids enrollment statuses to broader reasons for why
# testing was not completed. To change a grouping, edit the category on the
# right; any status not listed here is assigned to "Reason Unknown".
untested_status_map <- c(
  "Declined scheduling" = "Declined Testing",
  "Declined post-consent" = "Declined Testing",
  "Declined during consent" = "Declined Testing",
  "Family will contact us to schedule if interested" = "Declined Testing",
  "Family will contact us to consent if interested (consent visit completed)" = "Declined Testing",
  "UTC for scheduling" = "Unable to Contact",
  "UTC for consent forms" = "Unable to Contact",
  "Reschedule (no show/cancel)" = "Unable to Contact",
  "Deceased (no contact at this time)" = "Deceased"
)

# Change the numbers here to use different plotting symbols.
untested_reason_shapes <- c(
  "Declined Testing" = 17,
  "Unable to Contact" = 15,
  "Deceased" = 4,
  "Reason Unknown" = 1
)

untested_percentile <- cohort_manifest %>% 
  filter(is.na(patient_result)) %>% 
  drop_na(recorded_review_percentile) %>% 
  select(recorded_review_percentile, patient_result, gk_enrollment_status) %>% 
  mutate(patient_result = "Untested",
         subset = "Untested",
         untested_reason = recode(
           gk_enrollment_status,
           !!!untested_status_map,
           .default = "Reason Unknown",
           .missing = "Reason Unknown"
         ),
         untested_reason = factor(
           untested_reason,
           levels = names(untested_reason_shapes)
         ))

grouped_percentiles_data <- bind_rows(tested_percentile, untested_percentile) %>%
  mutate(
    subset = factor(
      subset,
      levels = c("Tested", "Untested")
    ),
    patient_result = factor(
      patient_result,
      levels = c("Positive", "Uncertain", "Negative", "Untested")
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
  # "Incomplete" = "#8E44AD",
  "Untested"   = "#999999"
)

significance_data <- tibble(
  subset = factor(
    c("Tested", "Tested"),
    levels = c("Tested", "Untested")
  ),
  xmin = c(1, 1),
  xmax = c(2, 3),
  y_position = c(102, 106),
  annotations = c("p=0.10", "p=0.0012")
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
    linewidth = 0.9
  ) +
  geom_jitter(
    data = grouped_percentiles_data %>% filter(subset == "Tested"),
    width = 0.15,
    height = 0,
    size = 1.9,
    alpha = 0.5
  ) +
  geom_jitter(
    data = grouped_percentiles_data %>% filter(subset == "Untested"),
    aes(shape = untested_reason),
    width = 0.15,
    height = 0,
    size = 1.9,
    alpha = 0.5
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
  scale_shape_manual(values = untested_reason_shapes, drop = FALSE) +
  guides(
    fill = "none",
    color = "none",
    shape = guide_legend(
      title = "Reason not tested",
      override.aes = list(color = "#555555", alpha = 0.8, size = 3)
    )
  ) +
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
    strip.text = element_blank(),
    axis.text.x = element_text(size = 16),
    axis.title.x = element_text(size = 17),
    axis.text.y = element_text(size = 16),
    axis.title.y = element_text(size = 17),
    plot.margin = margin(t=45, r=10, b=10, l=10, unit="pt"),
    panel.spacing.x = unit(1.4, "lines"),
    legend.position = "right"
  )


ggsave(
  filename = "/home/bennet/Documents/phd/yandell/gene_kids/data/GK_MPSE_analysis/output/figures/grouped_percentile_boxplots_small.svg",
  plot = grouped_percentiles_plot,
  width = 8,
  height = 10,
  units = "in",
  bg = "white"
)
