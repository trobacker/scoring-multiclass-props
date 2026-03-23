#!/usr/bin/env Rscript
# Comprehensive plotting for ASYMMETRIC alpha scenarios (2, 3, 15)
# Creates energy score plots and discrimination plots across all n_counts

library(ggplot2)
library(tidyr)
library(dplyr)

# Define paths
data_dir <- "./results/data/alpha_2_3_15"
output_dir <- "./results/plotting_code"

cat("=========================================\n")
cat("ASYMMETRIC ALPHA PLOTTING\n")
cat("Alpha = (2, 3, 15)\n")
cat("=========================================\n\n")

# Scenario file names and labels
scenarios <- data.frame(
  file = c(
    "scenario1_misspecified_alpha_2_3_15.rds",
    "scenario2_narrow_calibrated_alpha_2_3_15.rds",
    "scenario3_dispersed_calibrated_alpha_2_3_15.rds",
    "scenario4_less_misspecified_alpha_2_3_15.rds"
  ),
  scenario_label = c(
    "Scenario 1: Misspecified",
    "Scenario 2: Narrow/Calibrated",
    "Scenario 3: Dispersed/Calibrated",
    "Scenario 4: Less Misspecified"
  ),
  stringsAsFactors = FALSE
)

cat("Loading and combining scenario data...\n")

# Load all scenarios
all_data <- list()
for (i in 1:nrow(scenarios)) {
  file_path <- file.path(data_dir, scenarios$file[i])
  cat("  Loading:", scenarios$file[i], "\n")
  df <- readRDS(file_path)
  df$scenario_label <- scenarios$scenario_label[i]
  all_data[[i]] <- df
}

# Combine all scenarios
combined_df <- do.call(rbind, all_data)

cat("Total rows loaded:", nrow(combined_df), "\n\n")

# Add method labels
combined_df <- combined_df %>%
  mutate(
    method_label = case_when(
      method == "mn" ~ "Multinomial",
      method == "non-mn" ~ "Non-Multinomial"
    )
  )

# Calculate summary statistics
cat("Computing summary statistics...\n")
summary_stats <- combined_df %>%
  group_by(scenario_label, method_label, n_counts) %>%
  summarise(
    mean_es_correct = mean(es_correct),
    mean_es_incorrect = mean(es_incorrect),
    mean_es_diff = mean(es_diff),
    prop_correct_wins = mean(es_diff < 0),
    .groups = "drop"
  )

# Save summary statistics
write.csv(summary_stats,
          file.path(output_dir, "asymmetric_alpha_summary_stats.csv"),
          row.names = FALSE)
cat("  Summary statistics saved\n\n")

# ============================================================================
# PLOT 1: Energy Scores vs n_counts (Mean values)
# ============================================================================
cat("Creating Plot 1: Mean Energy Scores vs n_counts...\n")

# Reshape for plotting
es_long <- summary_stats %>%
  pivot_longer(
    cols = c(mean_es_correct, mean_es_incorrect),
    names_to = "forecast_type",
    values_to = "mean_energy_score"
  ) %>%
  mutate(
    forecast_type = case_when(
      forecast_type == "mean_es_correct" ~ "Correct Forecast",
      forecast_type == "mean_es_incorrect" ~ "Incorrect Forecast"
    )
  )

p1 <- ggplot(es_long,
             aes(x = n_counts, y = mean_energy_score,
                 color = scenario_label,
                 linetype = forecast_type,
                 shape = method_label)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  facet_wrap(~ scenario_label, ncol = 2) +
  scale_x_log10(
    breaks = c(1, 2, 3, 4, 5, 10, 25, 50, 100, 500),
    labels = c(1, 2, 3, 4, 5, 10, 25, 50, 100, 500)
  ) +
  scale_color_manual(
    values = c(
      "Scenario 1: Misspecified" = "#E41A1C",
      "Scenario 2: Narrow/Calibrated" = "#377EB8",
      "Scenario 3: Dispersed/Calibrated" = "#4DAF4A",
      "Scenario 4: Less Misspecified" = "#984EA3"
    ),
    name = "Scenario"
  ) +
  scale_linetype_manual(
    values = c("Correct Forecast" = "solid",
               "Incorrect Forecast" = "dashed"),
    name = "Forecast Type"
  ) +
  scale_shape_manual(
    values = c("Multinomial" = 16, "Non-Multinomial" = 17),
    name = "Method"
  ) +
  labs(
    title = "Mean Energy Scores by Scenario: Asymmetric Alpha (2,3,15)",
    subtitle = "Lower scores indicate better forecasts",
    x = "n_counts (log scale)",
    y = "Mean Energy Score"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", size = 13),
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(face = "bold", size = 10)
  ) +
  guides(
    color = guide_legend(order = 1, ncol = 2),
    linetype = guide_legend(order = 2),
    shape = guide_legend(order = 3)
  )

ggsave(file.path(output_dir, "asymmetric_alpha_energy_scores_vs_n.png"),
       p1, width = 11, height = 8, dpi = 300)
cat("  Plot saved: asymmetric_alpha_energy_scores_vs_n.png\n\n")

# ============================================================================
# PLOT 2: Discrimination (Win %) vs n_counts
# ============================================================================
cat("Creating Plot 2: Discrimination (Win %) vs n_counts...\n")

p2 <- ggplot(summary_stats,
             aes(x = n_counts, y = prop_correct_wins,
                 color = scenario_label,
                 linetype = method_label,
                 shape = method_label)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "black", alpha = 0.5) +
  facet_wrap(~ scenario_label, ncol = 2) +
  scale_x_log10(
    breaks = c(1, 2, 3, 4, 5, 10, 25, 50, 100, 500),
    labels = c(1, 2, 3, 4, 5, 10, 25, 50, 100, 500)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    labels = scales::percent
  ) +
  scale_color_manual(
    values = c(
      "Scenario 1: Misspecified" = "#E41A1C",
      "Scenario 2: Narrow/Calibrated" = "#377EB8",
      "Scenario 3: Dispersed/Calibrated" = "#4DAF4A",
      "Scenario 4: Less Misspecified" = "#984EA3"
    ),
    name = "Scenario"
  ) +
  scale_linetype_manual(
    values = c("Multinomial" = "solid", "Non-Multinomial" = "dashed"),
    name = "Method"
  ) +
  scale_shape_manual(
    values = c("Multinomial" = 16, "Non-Multinomial" = 17),
    name = "Method"
  ) +
  labs(
    title = "Discrimination Power by Scenario: Asymmetric Alpha (2,3,15)",
    subtitle = "Proportion of simulations where correct forecast scores better (lower ES)",
    x = "n_counts (log scale)",
    y = "Proportion Correct Wins"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", size = 13),
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(face = "bold", size = 10)
  ) +
  guides(
    color = guide_legend(order = 1, ncol = 2),
    linetype = guide_legend(order = 2),
    shape = guide_legend(order = 2)
  )

ggsave(file.path(output_dir, "asymmetric_alpha_discrimination_vs_n.png"),
       p2, width = 11, height = 8, dpi = 300)
cat("  Plot saved: asymmetric_alpha_discrimination_vs_n.png\n\n")

# ============================================================================
# PLOT 3: Energy Score Differences vs n_counts
# ============================================================================
cat("Creating Plot 3: Energy Score Differences vs n_counts...\n")

p3 <- ggplot(summary_stats,
             aes(x = n_counts, y = mean_es_diff,
                 color = scenario_label,
                 linetype = method_label,
                 shape = method_label)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", alpha = 0.5) +
  facet_wrap(~ scenario_label, ncol = 2, scales = "free_y") +
  scale_x_log10(
    breaks = c(1, 2, 3, 4, 5, 10, 25, 50, 100, 500),
    labels = c(1, 2, 3, 4, 5, 10, 25, 50, 100, 500)
  ) +
  scale_color_manual(
    values = c(
      "Scenario 1: Misspecified" = "#E41A1C",
      "Scenario 2: Narrow/Calibrated" = "#377EB8",
      "Scenario 3: Dispersed/Calibrated" = "#4DAF4A",
      "Scenario 4: Less Misspecified" = "#984EA3"
    ),
    name = "Scenario"
  ) +
  scale_linetype_manual(
    values = c("Multinomial" = "solid", "Non-Multinomial" = "dashed"),
    name = "Method"
  ) +
  scale_shape_manual(
    values = c("Multinomial" = 16, "Non-Multinomial" = 17),
    name = "Method"
  ) +
  labs(
    title = "Mean Energy Score Differences by Scenario: Asymmetric Alpha (2,3,15)",
    subtitle = "ES Difference = Correct - Incorrect (negative values = correct forecast is better)",
    x = "n_counts (log scale)",
    y = "Mean ES Difference"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", size = 13),
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(face = "bold", size = 10)
  ) +
  guides(
    color = guide_legend(order = 1, ncol = 2),
    linetype = guide_legend(order = 2),
    shape = guide_legend(order = 2)
  )

ggsave(file.path(output_dir, "asymmetric_alpha_es_differences_vs_n.png"),
       p3, width = 11, height = 8, dpi = 300)
cat("  Plot saved: asymmetric_alpha_es_differences_vs_n.png\n\n")

# ============================================================================
# Print summary
# ============================================================================
cat("=========================================\n")
cat("SUMMARY\n")
cat("=========================================\n")
cat("Scenarios analyzed:", nrow(scenarios), "\n")
cat("Total simulations:", nrow(combined_df), "\n")
cat("\nFiles created:\n")
cat("  - asymmetric_alpha_energy_scores_vs_n.png\n")
cat("  - asymmetric_alpha_discrimination_vs_n.png\n")
cat("  - asymmetric_alpha_es_differences_vs_n.png\n")
cat("  - asymmetric_alpha_summary_stats.csv\n\n")

cat("Done!\n")
