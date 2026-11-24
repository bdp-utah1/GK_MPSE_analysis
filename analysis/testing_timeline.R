library(tidyverse)

master <- read_csv("output/data/master_cohort_analysis.csv")


df <- master %>% 
  filter(reviewed == TRUE & !is.na(patient_result)) %>% 
  select(project_id, 
         dob, 
         reviewed_raw_mpse_admit_date, 
         reviewed_raw_mpse_review_date,
         reviewed_raw_mpse_genetic_consult_date, 
         reviewed_raw_gk_enrollment_date,
         #primary_test_raw_genetic_test,
         primary_test_raw_order_submitted,
         #primary_test_raw_sample_collection_date,
         #primary_test_raw_sample_received_date,
         primary_test_raw_result_date,
         patient_result) %>% 
  rename(birth_date = dob,
         admit_date = reviewed_raw_mpse_admit_date,
         review_date = reviewed_raw_mpse_review_date,
         genetic_consult_date = reviewed_raw_mpse_genetic_consult_date, 
         GK_enrollment_date = reviewed_raw_gk_enrollment_date,
         #test_type = primary_test_raw_genetic_test,
         test_order_date = primary_test_raw_order_submitted,
         #test_sample_collection_date = primary_test_raw_sample_collection_date,
         #test_sample_received_date = primary_test_raw_sample_received_date,
         test_result_date = primary_test_raw_result_date,
         test_diagnostic_category = patient_result) %>% 
  mutate(
    test_diagnostic_category = factor(
      test_diagnostic_category,
      levels = c("Positive", "Uncertain", "Negative")
    )
  )


# Patient-level diagnostic timeline ---------------------------------------

# Birth date is retained in df for other analyses but is intentionally not
# included as an event in the diagnostic timeline. Admission is represented by
# the vertical day-zero reference line rather than a repeated point.
event_key <- tribble(
  ~milestone,                  ~label,
  "genetic_consult_date",      "Genetics consult",
  "review_date",               "MPSE review",
  "GK_enrollment_date",        "GK enrollment",
  "test_order_date",           "Test order",
  "test_result_date",          "Result date"
)

# Sort patients within each diagnostic category by time from admission to
# result. Patients without a recorded result appear last in their facet.
timeline_base <- df %>%
  mutate(result_days = as.integer(test_result_date - admit_date)) %>%
  group_by(test_diagnostic_category) %>%
  arrange(
    is.na(result_days),
    result_days,
    project_id,
    .by_group = TRUE
  ) %>%
  mutate(patient_rank = row_number()) %>%
  ungroup()

timeline_long <- timeline_base %>%
  mutate(admission_anchor = admit_date) %>%
  select(
    project_id,
    patient_rank,
    test_diagnostic_category,
    admission_anchor,
    all_of(event_key$milestone)
  ) %>%
  pivot_longer(
    cols = all_of(event_key$milestone),
    names_to = "milestone",
    values_to = "event_date"
  ) %>%
  filter(!is.na(event_date)) %>%
  left_join(event_key, by = "milestone") %>%
  mutate(
    days = as.integer(event_date - admission_anchor),
    label = factor(label, levels = event_key$label),
    timeline_panel = factor(
      as.character(test_diagnostic_category),
      levels = c("Positive", "Uncertain", "Negative", "Category means")
    )
  )

# Use the same x-axis range and breaks for the timeline and boxplot figures.
shared_x_limits <- range(timeline_long$days, na.rm = TRUE)
shared_x_breaks <- scales::breaks_pretty(n = 8)(shared_x_limits)

# The thin line is a visual guide spanning all observed milestones for a
# patient; it does not imply continuous treatment or follow-up.
timeline_spans <- timeline_long %>%
  group_by(
    timeline_panel,
    test_diagnostic_category,
    project_id,
    patient_rank
  ) %>%
  summarise(
    first_day = min(c(0, days)),
    last_day = max(c(0, days)),
    .groups = "drop"
  )

category_labels <- timeline_base %>%
  count(test_diagnostic_category, name = "n") %>%
  transmute(
    timeline_panel = as.character(test_diagnostic_category),
    facet_label = paste0(test_diagnostic_category, " (n = ", n, ")")
  ) %>%
  bind_rows(
    tibble(
      timeline_panel = "Category means",
      facet_label = "Category means"
    )
  ) %>%
  deframe()

# Add a compact fourth facet containing the mean elapsed time for every
# milestone within each diagnostic category.
mean_y_positions <- tibble(
  test_diagnostic_category = factor(
    c("Positive", "Uncertain", "Negative"),
    levels = c("Positive", "Uncertain", "Negative")
  ),
  mean_y = c(1.25, 2.5, 3.75)
)

timeline_means <- timeline_long %>%
  group_by(test_diagnostic_category, label) %>%
  summarise(mean_days = mean(days), .groups = "drop") %>%
  left_join(mean_y_positions, by = "test_diagnostic_category") %>%
  mutate(
    timeline_panel = factor(
      "Category means",
      levels = c("Positive", "Uncertain", "Negative", "Category means")
    )
  )

timeline_mean_spans <- timeline_means %>%
  group_by(timeline_panel, test_diagnostic_category, mean_y) %>%
  summarise(
    first_day = min(c(0, mean_days)),
    last_day = max(c(0, mean_days)),
    .groups = "drop"
  )

timeline_mean_labels <- timeline_means %>%
  distinct(timeline_panel, test_diagnostic_category, mean_y) %>%
  mutate(label_x = -6)

timeline_mean_panel_extent <- tibble(
  timeline_panel = factor(
    c("Category means", "Category means"),
    levels = c("Positive", "Uncertain", "Negative", "Category means")
  ),
  extent_y = c(0, 5)
)

milestone_colours <- c(
  "Genetics consult" = "#CC79A7",
  "MPSE review" = "#0072B2",
  "GK enrollment" = "#009E73",
  "Test order" = "#E69F00",
  "Result date" = "#737373"
)

timeline_plot <- ggplot() +
  geom_vline(
    xintercept = 0,
    colour = "grey45",
    linewidth = 0.4
  ) +
  geom_segment(
    data = timeline_spans,
    aes(
      x = first_day,
      xend = last_day,
      y = patient_rank,
      yend = patient_rank
    ),
    colour = "grey82",
    linewidth = 0.45
  ) +
  geom_point(
    data = timeline_long,
    aes(
      x = days,
      y = patient_rank,
      colour = label,
      shape = label
    ),
    size = 1.7,
    stroke = 0.5,
    position = position_jitter(
      width = 0,
      height = 0.12,
      seed = 123
    )
  ) +
  geom_segment(
    data = timeline_mean_spans,
    aes(
      x = first_day,
      xend = last_day,
      y = mean_y,
      yend = mean_y
    ),
    colour = "grey55",
    linewidth = 0.65
  ) +
  geom_point(
    data = timeline_means,
    aes(
      x = mean_days,
      y = mean_y,
      colour = label,
      shape = label
    ),
    size = 2.1,
    stroke = 0.6
  ) +
  geom_text(
    data = timeline_mean_labels,
    aes(
      x = label_x,
      y = mean_y,
      label = test_diagnostic_category
    ),
    hjust = 1,
    colour = "grey20",
    size = 3
  ) +
  geom_blank(
    data = timeline_mean_panel_extent,
    aes(x = 0, y = extent_y)
  ) +
  facet_grid(
    rows = vars(timeline_panel),
    scales = "free_y",
    space = "free_y",
    labeller = labeller(timeline_panel = category_labels)
  ) +
  scale_y_reverse(breaks = NULL) +
  scale_x_continuous(
    breaks = shared_x_breaks,
    limits = shared_x_limits,
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_colour_manual(
    values = milestone_colours,
    name = NULL
  ) +
  scale_shape_manual(
    values = c(17, 15, 4, 8, 7),
    name = NULL
  ) +
  guides(
    colour = guide_legend(nrow = 1, byrow = TRUE),
    shape = guide_legend(nrow = 1, byrow = TRUE)
  ) +
  labs(
    x = "Days relative to admission",
    y = "Patients, sorted by admission-to-result time"
  ) +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    strip.text.y = element_text(angle = 0)
  )

timeline_plot

figure_dir <- "output/figures"
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(
  filename = file.path(figure_dir, "testing_timeline.pdf"),
  plot = timeline_plot,
  width = 7.5,
  height = 9.5
)

ggsave(
  filename = file.path(figure_dir, "testing_timeline.png"),
  plot = timeline_plot,
  width = 7.5,
  height = 9.5,
  dpi = 600,
  bg = "white"
)


# Elapsed-time data by diagnostic category -------------------------------

# Signed elapsed times are retained. Negative values indicate that testing
# occurred before the MPSE admission used as the time origin.
elapsed_time_data <- df %>%
  transmute(
    project_id,
    test_diagnostic_category,
    mpse_review_days = as.integer(review_date - admit_date),
    `MPSE review` = as.integer(review_date - admit_date),
    `Test order` = as.integer(test_order_date - admit_date),
    `Result date` = as.integer(test_result_date - admit_date)
  ) %>%
  pivot_longer(
    cols = c(`MPSE review`, `Test order`, `Result date`),
    names_to = "event",
    values_to = "days"
  ) %>%
  filter(!is.na(days)) %>%
  mutate(
    event = factor(
      event,
      levels = c("MPSE review", "Test order", "Result date")
    )
  )

category_colours <- c(
  "Negative" = "#0072B2",
  "Positive" = "#D55E00",
  "Uncertain" = "#009E73"
)

# Boxplot and mean comparisons -------------------------------------------

# The boxplots show the full patient-level distributions. Means are marked
# with white diamonds because the accompanying tests compare means, while the
# box and centre line show the interquartile range and median, respectively.
boxplot_data <- elapsed_time_data %>%
  mutate(
    test_diagnostic_category = factor(
      test_diagnostic_category,
      levels = c("Positive", "Uncertain", "Negative")
    ))

elapsed_time_summary <- boxplot_data %>%
  group_by(event, test_diagnostic_category) %>%
  summarise(
    n = n(),
    mean_days = mean(days),
    sd_days = sd(days),
    median_days = median(days),
    q1_days = quantile(days, 0.25),
    q3_days = quantile(days, 0.75),
    min_days = min(days),
    max_days = max(days),
    .groups = "drop"
  )

# Pairwise Welch tests compare means without assuming equal variances or equal
# group sizes. P-values are reported without multiplicity adjustment.
category_pairs <- combn(
  levels(boxplot_data$test_diagnostic_category),
  2,
  simplify = FALSE
)

mean_ci_95 <- function(values) {
  n <- length(values)
  mean_value <- mean(values)
  margin <- qt(0.975, df = n - 1) * sd(values) / sqrt(n)

  c(
    mean = mean_value,
    conf_low = mean_value - margin,
    conf_high = mean_value + margin
  )
}

run_pairwise_tests <- function(data) {
  data %>%
    group_by(event) %>%
    group_modify(~ {
      map_dfr(category_pairs, function(pair) {
        group_1_values <- .x %>%
          filter(test_diagnostic_category == pair[1]) %>%
          pull(days)
        group_2_values <- .x %>%
          filter(test_diagnostic_category == pair[2]) %>%
          pull(days)

        group_1_summary <- mean_ci_95(group_1_values)
        group_2_summary <- mean_ci_95(group_2_values)

        welch_fit <- t.test(
          group_1_values,
          group_2_values,
          var.equal = FALSE
        )

        tibble(
          group_1 = pair[1],
          group_2 = pair[2],
          n_group_1 = length(group_1_values),
          n_group_2 = length(group_2_values),
          mean_group_1 = group_1_summary["mean"],
          mean_group_1_conf_low = group_1_summary["conf_low"],
          mean_group_1_conf_high = group_1_summary["conf_high"],
          mean_group_2 = group_2_summary["mean"],
          mean_group_2_conf_low = group_2_summary["conf_low"],
          mean_group_2_conf_high = group_2_summary["conf_high"],
          mean_difference = group_1_summary["mean"] - group_2_summary["mean"],
          welch_conf_low = welch_fit$conf.int[1],
          welch_conf_high = welch_fit$conf.int[2],
          welch_statistic = unname(welch_fit$statistic),
          welch_degrees_freedom = unname(welch_fit$parameter),
          welch_p_value = welch_fit$p.value
        )
      })
    }) %>%
    ungroup()
}

pairwise_results <- run_pairwise_tests(boxplot_data)

# Sensitivity analysis: for test order and result date, exclude observations
# that occurred before the patient's MPSE review. MPSE review observations are
# retained because the event and reference date are identical.
sensitivity_data <- boxplot_data %>%
  filter(days >= mpse_review_days)

pairwise_sensitivity_results <- run_pairwise_tests(
  sensitivity_data
)

# Sensitivity analysis: retain only observations that occurred on or after
# admission. Same-day observations are retained.
admission_sensitivity_data <- boxplot_data %>%
  filter(days >= 0)

pairwise_admission_sensitivity_results <- run_pairwise_tests(
  admission_sensitivity_data
)

# Publication-oriented table of pairwise Welch comparisons.
format_estimate_ci <- function(estimate, conf_low, conf_high) {
  sprintf("%.1f (%.1f to %.1f)", estimate, conf_low, conf_high)
}

format_p_value <- function(p_value) {
  if_else(
    p_value < 0.001,
    "<0.001",
    sprintf("%.3f", p_value)
  )
}

pairwise_comparison_table <- bind_rows(
  pairwise_results %>%
    mutate(analysis = "Primary analysis"),
  pairwise_sensitivity_results %>%
    mutate(
      analysis = "Sensitivity: on/after MPSE review"
    ),
  pairwise_admission_sensitivity_results %>%
    mutate(
      analysis = "Sensitivity: on/after admission"
    )
) %>%
  transmute(
    Analysis = analysis,
    Event = as.character(event),
    Comparison = paste(group_1, "vs", group_2),
    `Group 1 mean days (95% CI)` = paste0(
      group_1,
      ", n = ",
      n_group_1,
      ": ",
      format_estimate_ci(
        mean_group_1,
        mean_group_1_conf_low,
        mean_group_1_conf_high
      )
    ),
    `Group 2 mean days (95% CI)` = paste0(
      group_2,
      ", n = ",
      n_group_2,
      ": ",
      format_estimate_ci(
        mean_group_2,
        mean_group_2_conf_low,
        mean_group_2_conf_high
      )
    ),
    `Mean difference days (95% CI)` = format_estimate_ci(
      mean_difference,
      welch_conf_low,
      welch_conf_high
    ),
    `Welch p` = format_p_value(welch_p_value)
  )

boxplot_plot <- ggplot(
  boxplot_data,
  aes(
    x = days,
    y = event,
    group = interaction(event, fct_rev(test_diagnostic_category)),
    fill = test_diagnostic_category,
    colour = test_diagnostic_category
  )
) +
  geom_vline(
    xintercept = 0,
    colour = "grey55",
    linewidth = 0.4
  ) +
  geom_boxplot(
    width = 0.48,
    alpha = 0.35,
    linewidth = 0.55,
    outlier.shape = NA,
    orientation = "y",
    position = position_dodge(width = 0.62)
  ) +
  geom_jitter(
    alpha = 0.5,
    size = 1.25,
    stroke = 0,
    position = position_jitterdodge(
      jitter.width = 0,
      jitter.height = 0.05,
      dodge.width = 0.62,
      seed = 123
    )
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 23,
    size = 2.8,
    stroke = 0.65,
    colour = "grey20",
    fill = "white",
    orientation = "y",
    position = position_dodge(width = 0.62)
  ) +
  scale_fill_manual(
    values = category_colours,
    name = NULL
  ) +
  scale_colour_manual(
    values = category_colours,
    name = NULL
  ) +
  scale_x_continuous(
    breaks = shared_x_breaks,
    limits = shared_x_limits,
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_y_discrete(
    limits = rev(levels(boxplot_data$event))
  ) +
  guides(
    fill = guide_legend(nrow = 1, byrow = TRUE),
    colour = guide_legend(nrow = 1, byrow = TRUE)
  ) +
  labs(
    x = "Days relative to admission",
    y = "Event"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom"
  )

boxplot_plot

ggsave(
  filename = file.path(figure_dir, "testing_timeline_boxplot.pdf"),
  plot = boxplot_plot,
  width = 9,
  height = 4
)

ggsave(
  filename = file.path(figure_dir, "testing_timeline_boxplot.png"),
  plot = boxplot_plot,
  width = 9,
  height = 4,
  dpi = 600,
  bg = "white"
)

statistics_dir <- "output/statistics"
dir.create(statistics_dir, recursive = TRUE, showWarnings = FALSE)

write_csv(
  elapsed_time_summary,
  file.path(statistics_dir, "testing_timeline_descriptive_statistics.csv")
)

write_csv(
  pairwise_results,
  file.path(statistics_dir, "testing_timeline_pairwise_tests.csv")
)

write_csv(
  pairwise_sensitivity_results,
  file.path(
    statistics_dir,
    "testing_timeline_pairwise_tests_after_mpse_review.csv"
  )
)

write_csv(
  pairwise_admission_sensitivity_results,
  file.path(
    statistics_dir,
    "testing_timeline_pairwise_tests_after_admission.csv"
  )
)

write_csv(
  pairwise_comparison_table,
  file.path(
    statistics_dir,
    "testing_timeline_pairwise_comparison_table.csv"
  )
)
