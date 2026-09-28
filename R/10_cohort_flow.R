build_cohort_flow <- function(cohorts) {
  screened <- cohorts$screened
  reviewed <- cohorts$reviewed

  counts <- list(
    screened = nrow(screened),
    reviewed = nrow(reviewed),
    not_reviewed = nrow(screened) - nrow(reviewed),
    tested = sum(reviewed$has_any_test),
    not_tested = sum(!reviewed$has_any_test),
    classifiable = sum(reviewed$has_classifiable_completed_result),
    not_classifiable = sum(
      reviewed$has_any_test & !reviewed$has_classifiable_completed_result
    ),
    positive = sum(reviewed$patient_result == "Positive", na.rm = TRUE),
    uncertain = sum(reviewed$patient_result == "Uncertain", na.rm = TRUE),
    negative = sum(reviewed$patient_result == "Negative", na.rm = TRUE),
    wgs = sum(reviewed$has_wgs),
    gk_enrolled = sum(reviewed$gk_enrolled)
  )

  flow_table <- dplyr::bind_rows(
    make_count_row("Screening", "Screened by MPSE", counts$screened),
    make_count_row(
      "Clinical review",
      "Reviewed",
      counts$reviewed,
      counts$screened,
      "Screened"
    ),
    make_count_row(
      "Clinical review",
      "Not reviewed",
      counts$not_reviewed,
      counts$screened,
      "Screened"
    ),
    make_count_row(
      "Testing",
      "Any genetic testing",
      counts$tested,
      counts$reviewed,
      "Reviewed"
    ),
    make_count_row(
      "Testing",
      "No testing record",
      counts$not_tested,
      counts$reviewed,
      "Reviewed"
    ),
    make_count_row(
      "Results",
      "Completed and classifiable result",
      counts$classifiable,
      counts$tested,
      "Tested"
    ),
    make_count_row(
      "Results",
      "No completed classifiable result",
      counts$not_classifiable,
      counts$tested,
      "Tested"
    ),
    make_count_row(
      "Diagnostic result",
      "Positive",
      counts$positive,
      counts$classifiable,
      "Completed and classifiable results"
    ),
    make_count_row(
      "Diagnostic result",
      "Uncertain",
      counts$uncertain,
      counts$classifiable,
      "Completed and classifiable results"
    ),
    make_count_row(
      "Diagnostic result",
      "Negative",
      counts$negative,
      counts$classifiable,
      "Completed and classifiable results"
    ),
    make_count_row(
      "Overlapping subset",
      "WGS",
      counts$wgs,
      counts$reviewed,
      "Reviewed"
    ),
    make_count_row(
      "Overlapping subset",
      "GK enrolled",
      counts$gk_enrolled,
      counts$reviewed,
      "Reviewed"
    )
  )

  nodes <- tibble::tribble(
    ~node, ~label, ~x, ~y, ~count,
    "screened", "Screened by MPSE", 1, 0, counts$screened,
    "reviewed", "Clinician reviewed", 2, 1.2, counts$reviewed,
    "not_reviewed", "Not reviewed", 2, -1.2, counts$not_reviewed,
    "tested", "Any genetic testing", 3, 2.0, counts$tested,
    "not_tested", "No testing record", 3, 0.3, counts$not_tested,
    "classifiable", "Completed and\nclassifiable result", 4, 2.8,
      counts$classifiable,
    "not_classifiable", "No completed\nclassifiable result", 4, 1.2,
      counts$not_classifiable,
    "positive", "Positive", 5, 3.6, counts$positive,
    "uncertain", "Uncertain", 5, 2.8, counts$uncertain,
    "negative", "Negative", 5, 2.0, counts$negative
  ) |>
    dplyr::mutate(node_label = paste0(.data$label, "\n", "n = ", .data$count))

  edges <- tibble::tribble(
    ~from, ~to,
    "screened", "reviewed",
    "screened", "not_reviewed",
    "reviewed", "tested",
    "reviewed", "not_tested",
    "tested", "classifiable",
    "tested", "not_classifiable",
    "classifiable", "positive",
    "classifiable", "uncertain",
    "classifiable", "negative"
  ) |>
    dplyr::left_join(
      nodes |>
        dplyr::select(from = node, x_from = x, y_from = y),
      by = "from"
    ) |>
    dplyr::left_join(
      nodes |>
        dplyr::select(to = node, x_to = x, y_to = y),
      by = "to"
    )

  plot <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = edges,
      ggplot2::aes(
        x = .data$x_from + 0.18,
        y = .data$y_from,
        xend = .data$x_to - 0.18,
        yend = .data$y_to
      ),
      linewidth = 0.6,
      color = "#7A8694",
      arrow = grid::arrow(length = grid::unit(0.12, "inches"))
    ) +
    ggplot2::geom_label(
      data = nodes,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$node_label),
      size = 3.4,
      lineheight = 1.05,
      label.size = 0.35,
      label.padding = grid::unit(0.18, "lines"),
      color = "#18212B",
      fill = "#F7F9FB"
    ) +
    ggplot2::annotate(
      "label",
      x = 3.5,
      y = -1.25,
      label = paste0(
        "Overlapping reviewed-cohort subsets\n",
        "WGS: n = ", counts$wgs,
        "    GK enrolled: n = ", counts$gk_enrolled
      ),
      size = 3.3,
      label.size = 0.3,
      fill = "#FFF7E6",
      color = "#5A451A"
    ) +
    ggplot2::annotate(
      "text",
      x = 3.5,
      y = -1.85,
      label = paste0(
        "WGS and GK enrollment are not treated as mutually exclusive ",
        "sequential stages."
      ),
      size = 3.1,
      color = "#5A6470"
    ) +
    ggplot2::coord_cartesian(
      xlim = c(0.6, 5.4),
      ylim = c(-2.1, 4.0),
      clip = "off"
    ) +
    ggplot2::labs(
      title = "MPSE implementation cohort flow",
      subtitle = paste0(
        "Diagnostic-result categories use completed, classifiable results only"
      )
    ) +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 15),
      plot.subtitle = ggplot2::element_text(color = "#4D5966"),
      plot.margin = ggplot2::margin(15, 25, 15, 25)
    )

  list(counts = counts, table = flow_table, plot = plot)
}

write_cohort_flow <- function(flow, config) {
  readr::write_csv(
    flow$table,
    file.path(config$output_dirs$tables, "cohort_flow.csv"),
    na = ""
  )

  gt_table <- gt::gt(flow$table) |>
    gt::opt_row_striping() |>
    gt::tab_header(title = "MPSE implementation cohort flow")
  gt::gtsave(
    gt_table,
    file.path(config$output_dirs$tables, "cohort_flow.html")
  )

  ggplot2::ggsave(
    filename = file.path(config$output_dirs$figures, "cohort_flow.png"),
    plot = flow$plot,
    width = 12,
    height = 7,
    units = "in",
    dpi = 300,
    bg = "white"
  )
}
