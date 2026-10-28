#!/usr/bin/env Rscript

# N-Focused Energy Score Term Decomposition Figure (standalone)
#
# Companion figure to n_focused_summary_energy_decomposition_plot.R, showing
# just the term decomposition panel (previously the RHS of that figure) on
# its own, without the mean energy score log-scale plot.
#
# For the same symmetric alpha=(10,10,10) setting and the same two
# wrong-alpha scenarios (scenario3 extreme_dispersed, scenario5
# extreme_narrow): energy score term decomposition vs n_counts
# (ES = term1 + (-term2)), following the style of
# results/term_decomposition_analysis/plot1_*, restricted to just these two
# scenarios.
#
# Output: one PNG (both scenarios shown together, faceted by column, term
# type faceted by row), no focal-N marker since this figure spans all N.

library(ggplot2)
library(dplyr)
library(tidyr)

# ============================================================================
# Config
# ============================================================================

data_dir   <- "results/data/alpha_10_10_10"
output_dir <- "results/n_focused_summary"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

scenario_files <- c(
  dispersed = "scenario3_extreme_dispersed_alpha_10_10_10.rds",
  narrow    = "scenario5_extreme_narrow_alpha_10_10_10.rds"
)

focal_n_values <- c(2, 5)

# ============================================================================
# Load data
# ============================================================================

cat("Loading scenario data...\n")

scenario_data <- lapply(names(scenario_files), function(key) {
  df <- readRDS(file.path(data_dir, scenario_files[key]))
  df$scenario_key <- key
  df
})
combined_df <- do.call(rbind, scenario_data)

combined_df <- combined_df %>%
  mutate(
    scenario_label = case_when(
      scenario_key == "dispersed" ~ "Dispersed (α=0.5,0.5,0.5)",
      scenario_key == "narrow"    ~ "Narrow (α=200,200,200)"
    ),
    method_label = case_when(
      method == "mn" ~ "Multinomial",
      method == "non-mn" ~ "Non-Multinomial"
    )
  )

# ============================================================================
# Term decomposition vs n_counts (ES = term1 + (-term2))
# ============================================================================

term_summary <- combined_df %>%
  group_by(scenario_label, method, n_counts) %>%
  summarise(
    mean_es_correct = mean(es_correct, na.rm = TRUE),
    mean_es_incorrect = mean(es_incorrect, na.rm = TRUE),
    mean_term1_correct = mean(term1_correct, na.rm = TRUE),
    mean_term1_incorrect = mean(term1_incorrect, na.rm = TRUE),
    mean_term2_correct = mean(term2_correct, na.rm = TRUE),
    mean_term2_incorrect = mean(term2_incorrect, na.rm = TRUE),
    .groups = "drop"
  )

term_long <- term_summary %>%
  pivot_longer(
    cols = c(mean_es_correct, mean_es_incorrect,
             mean_term1_correct, mean_term1_incorrect,
             mean_term2_correct, mean_term2_incorrect),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    # Flip term2 by multiplying by -1 so that ES = term1 + (-term2) visually
    value = ifelse(grepl("term2", metric), -value, value),
    term_type = factor(case_when(
      grepl("es", metric) ~ "ES",
      grepl("term1", metric) ~ "Term 1",
      grepl("term2", metric) ~ "-Term 2"
    ), levels = c("ES", "Term 1", "-Term 2")),
    forecast_type = ifelse(grepl("_correct$", metric), "Correct Forecast", "Incorrect Forecast")
  )

term_plot <- ggplot(term_long,
       aes(x = n_counts, y = value, color = forecast_type, linetype = method,
           group = interaction(forecast_type, method))) +
  # geom_vline(xintercept = focal_n_values, linetype = "dashed", color = "#984EA3", linewidth = 1.1) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  facet_grid(term_type ~ scenario_label, scales = "free_y") +
  scale_x_log10(breaks = c(1, 5, 25, 100, 500)) +
  scale_color_manual(
    values = c("Correct Forecast" = "#2166AC", "Incorrect Forecast" = "#E41A1C"),
    name = "Forecast Type"
  ) +
  scale_linetype_manual(
    values = c("mn" = "solid", "non-mn" = "dashed"),
    name = "Method",
    labels = c("mn" = "Multinomial", "non-mn" = "Non-Multinomial")
  ) +
  labs(
    title = "Symmetric α=(10,10,10): Energy Score Term Decomposition vs n_counts",
    subtitle = "ES = Term 1 + (-Term 2)",
    x = "N (log scale)",
    y = "Mean Value"
  ) +
  theme_bw(base_size = 22) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 19),
    plot.title = element_text(face = "bold", size = 23),
    plot.subtitle = element_text(size = 19),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 18),
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(face = "bold", size = 19)
  )

# ============================================================================
# Save
# ============================================================================

cat("Assembling term decomposition figure...\n")

out_path <- file.path(output_dir, "n_focused_term_decomposition_symmetric.png")
ggsave(out_path, term_plot, width = 14, height = 13, dpi = 400)

cat("  Saved:", out_path, "\n")
cat("\nDone!\n")
