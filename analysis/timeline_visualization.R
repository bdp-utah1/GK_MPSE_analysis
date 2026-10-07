library(tidyverse)

mosc_data <- read_csv("/home/bennet/Documents/phd/yandell/gene_kids/data/GK_MPSE_analysis_copy/mosaic_data/patient_level_data.csv") %>% 
  select(mrn, mpse_admit_date, mpse_review_date, gk_enrollment_date)

mpse_data <- read_tsv("/home/bennet/Documents/phd/yandell/gene_kids/data/GK_MPSE_analysis_copy/mpse_data/filtered_patient_stay.tsv") %>% 
  filter(!patient_id %in% mosc_data$mrn) %>% 
  arrange(patient_id, arrive_dtm) %>% 
  distinct(patient_id, .keep_all = TRUE) %>% 
  select(patient_id, arrive_dtm) %>% 
  rename(mrn="patient_id", mpse_admit_date="arrive_dtm") %>% 
  mutate(mpse_admit_date = as.Date(mpse_admit_date))

df <- bind_rows(mosc_data, mpse_data)


study_start <- as.Date("2025-08-01")
study_end <- as.Date("2026-07-01")

counts_2wk <- df %>%
  pivot_longer(
    cols = c(mpse_admit_date, mpse_review_date, gk_enrollment_date),
    names_to = "event_type",
    values_to = "event_date"
  ) %>%
  mutate(event_date = as.Date(event_date)) %>%
  filter(!is.na(event_date)) %>%  # removes patients without enrollment dates
  filter(event_date >= study_start,
         event_date <= study_end) %>%
  mutate(
    interval_start = study_start +
      14 * floor(as.numeric(event_date - study_start) / 14),
    event_type = recode(
      event_type,
      mpse_admit_date = "NICU Admits",
      mpse_review_date = "MPSE Reviews",
      gk_enrollment_date = "GK Enrollments"
    )
  ) %>%
  count(interval_start, event_type, name = "patient_count") %>%
  complete(
    interval_start = seq(
      study_start,
      max(interval_start),
      by = "14 days"
    ),
    event_type = c("NICU Admits", "MPSE Reviews", "GK Enrollments"),
    fill = list(patient_count = 0)
  )

counts_monthly <- df %>%
  pivot_longer(
    cols = c(mpse_admit_date, mpse_review_date, gk_enrollment_date),
    names_to = "event_type",
    values_to = "event_date"
  ) %>%
  mutate(event_date = as.Date(event_date)) %>%
  filter(!is.na(event_date), event_date >= study_start) %>%
  mutate(
    month_number = floor(
      time_length(interval(study_start, event_date), unit = "months")
    ),
    interval_start = study_start %m+% months(month_number),
    event_type = recode(
      event_type,
      mpse_admit_date = "NICU Admits",
      mpse_review_date = "MPSE Reviews",
      gk_enrollment_date = "GK Enrollments"
    )
  ) %>%
  count(interval_start, event_type, name = "patient_count") %>%
  complete(
    interval_start = seq(
      study_start,
      max(interval_start),
      by = "1 month"
    ),
    event_type = c("NICU Admits", "MPSE Reviews", "GK Enrollments"),
    fill = list(patient_count = 0)
  ) %>% 
  mutate(event_type = factor(event_type, levels=c("NICU Admits","MPSE Reviews","GK Enrollments")))


ggplot(
  counts_2wk,
  aes(
    x = interval_start,
    y = patient_count,
    color = event_type,
    group = event_type
  )
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b\n%Y") +
  scale_color_manual(
    values = c(
      "NICU Admits" = "#0072B2",
      "MPSE Reviews" = "#D55E00",
      "GK Enrollments" = "#009E73"
    )
  ) +
  labs(
    x = "Two-week interval beginning",
    y = "Patient count",
    color = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank()
  )




plot_start <- as.Date("2025-08-01")
plot_end   <- as.Date("2026-06-02")
month_barplot <- ggplot(
  counts_monthly,
  aes(
    x = interval_start,
    y = patient_count,
    fill = event_type
  )
) +
  geom_col(
    position = position_dodge(width = 24),
    width = 22
  ) +
  scale_x_date(
    breaks = seq(plot_start, plot_end, by = "1 month"),
    date_labels = "%b\n%Y",
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  coord_cartesian(xlim = c(plot_start, plot_end)) +
  scale_fill_manual(
    values = c(
      "NICU Admits" = "#0072B2",
      "MPSE Reviews" = "#D55E00",
      "GK Enrollments" = "#009E73"
    )
  ) +
  labs(
    x = "One-month interval",
    y = "Patient count",
    fill = NULL
  ) +
  theme_minimal(base_size = 22) +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 20),
    legend.key.height = grid::unit(0.28, "in"),
    legend.key.width = grid::unit(0.28, "in"),
    
    axis.title.x = element_text(
      size = 26,
      face = "bold",
      margin = margin(t = 12)
    ),
    axis.title.y = element_text(
      size = 26,
      face = "bold",
      margin = margin(r = 12)
    ),
    axis.text.x = element_text(
      size = 18,
      lineheight = 0.9,
      margin = margin(t = 8)
    ),
    axis.text.y = element_text(size = 18),
    
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.margin = margin(12, 18, 12, 12)
  )

ggsave(
  filename = "/home/bennet/Documents/phd/yandell/gene_kids/data/GK_MPSE_analysis/implementation_data/MPSE_study_timeline_vis_large.svg",
  plot = month_barplot,
  width = 18,
  height = 11,
  units = "in",
  bg = "white"
)
