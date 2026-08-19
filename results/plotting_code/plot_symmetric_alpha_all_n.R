#!/usr/bin/env Rscript
# Comprehensive plotting for SYMMETRIC alpha scenarios (10, 10, 10)
# Creates energy score plots and discrimination plots across all n_counts
# Includes all 7 scenarios (including ultra-dispersed and ultra-narrow)

library(ggplot2)
library(tidyr)
library(dplyr)

# Define paths
data_dir <- "./results/data/alpha_10_10_10"
output_dir <- "./results/plotting_code"

cat("=========================================\n")
cat("SYMMETRIC ALPHA PLOTTING\n")
cat("Alpha = (10, 10, 10)\n")
cat("=========================================\n\n")

# Scenario file names and labels
scenarios <- data.frame(
  file = c(
    "scenario1_mislocation_symmetric_alpha_10_10_10.rds",
    "scenario2_moderate_dispersed_alpha_10_10_10.rds",
    "scenario3_extreme_dispersed_alpha_10_10_10.rds",
    "scenario4_moderate_narrow_alpha_10_10_10.rds",
    "scenario5_extreme_narrow_alpha_10_10_10.rds",
    "scenario6_ultra_dispersed_alpha_10_10_10.rds",
    "scenario7_ultra_narrow_alpha_10_10_10.rds"
  ),
  scenario_label = c(
    "S1: Mislocation (α=5,10,15)",
    "S2: Moderate Dispersed (α=2,2,2)",
    "S3: Extreme Dispersed (α=0.5,0.5,0.5)",
    "S4: Moderate Narrow (α=50,50,50)",
    "S5: Extreme Narrow (α=200,200,200)",
    "S6: Ultra Dispersed (α=0.1,0.1,0.1)",
    "S7: Ultra Narrow (α=1000,1000,1000)"
  ),
  stringsAsFactors = FALSE
)

cat("Loading and combining scenario data...\n")

# Load all scenarios
all_data <- list()
for (i in 1:nrow(scenarios)) {
  file_path <- file.path(data_dir, scenarios$file[i])
  if (file.exists(file_path)) {
    cat("  Loading:", scenarios$file[i], "\n")
    df <- readRDS(file_path)
    df$scenario_label <- scenarios$scenario_label[i]
    all_data[[i]] <- df
  } else {
    cat("  WARNING: File not found:", scenarios$file[i], "\n")
  }
}

# Combine all scenarios
combined_df <- do.call(rbind, all_data)

cat("Total rows loaded:", nrow(combined_df), "\n")
cat("Scenarios loaded:", length(unique(combined_df$scenario_label)), "\n\n")

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
          file.path(output_dir, "symmetric_alpha_summary_stats.csv"),
          row.names = FALSE)
cat("  Summary statistics saved\n\n")

# Define color palette for 7 scenarios
scenario_colors <- c(
  "S1: Mislocation (α=5,10,15)" = "#984EA3",                    # purple (location error)
  "S2: Moderate Dispersed (α=2,2,2)" = "#A6D854",        # light green
  "S3: Extreme Dispersed (α=0.5,0.5,0.5)" = "#4DAF4A",        # green
  "S6: Ultra Dispersed (α=0.1,0.1,0.1)" = "#1B7837",         # dark green
  "S4: Moderate Narrow (α=50,50,50)" = "#FFD92F",           # light orange
  "S5: Extreme Narrow (α=200,200,200)" = "#E78AC3",           # pink
  "S7: Ultra Narrow (α=1000,1000,1000)" = "#E41A1C"             # red
)

# ============================================================================
# PLOT 1: Energy Scores vs n_counts (Mean values) - Faceted by method
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
                 color = forecast_type,
                 linetype = forecast_type,
                 shape = method_label)) +
  geom_line(linewidth = 0.9, alpha = 0.7) +
  geom_point(size = 2.5, alpha = 0.7) +
  facet_wrap(~ scenario_label, ncol = 3, scales = "free_y") +
  scale_x_log10(
    breaks = c(1, 5, 25, 100, 500),
    labels = c(1, 5, 25, 100, 500)
  ) +
  scale_y_log10() +
  scale_color_manual(
    values = c("Correct Forecast" = "#2166AC", "Incorrect Forecast" = "#E41A1C"),
    name = "Forecast Type"
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
    title = "Mean Energy Scores by Scenario: Symmetric Alpha (10,10,10)",
    subtitle = "Lower scores indicate better forecasts. Colors distinguish forecast types within each scenario.",
    x = "n_counts (log scale)",
    y = "Mean Energy Score (log scale)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    plot.title = element_text(face = "bold", size = 13),
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(face = "bold", size = 9)
  ) +
  guides(
    color = guide_legend(order = 1, override.aes = list(alpha = 1)),
    linetype = guide_legend(order = 1),
    shape = guide_legend(order = 2, override.aes = list(alpha = 1))
  )

ggsave(file.path(output_dir, "symmetric_alpha_energy_scores_vs_n.png"),
       p1, width = 13, height = 9, dpi = 300)
cat("  Plot saved: symmetric_alpha_energy_scores_vs_n.png\n\n")

# ============================================================================
# PLOT 2: Discrimination (Win %) vs n_counts
# ============================================================================
cat("Creating Plot 2: Discrimination (Win %) vs n_counts...\n")

p2 <- ggplot(summary_stats,
             aes(x = n_counts, y = prop_correct_wins,
                 color = method_label,
                 linetype = method_label,
                 shape = method_label)) +
  geom_line(linewidth = 1.2, alpha = 0.7) +
  geom_point(size = 3, alpha = 0.7) +
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "black", alpha = 0.5) +
  facet_wrap(~ scenario_label, ncol = 3) +
  scale_x_log10(
    breaks = c(1, 5, 25, 100, 500),
    labels = c(1, 5, 25, 100, 500)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    labels = scales::percent
  ) +
  scale_color_manual(
    values = c("Multinomial" = "#2166AC", "Non-Multinomial" = "#FF7F00"),
    name = "Method"
  ) +
  scale_linetype_manual(
    values = c("Multinomial" = "solid", "Non-Multinomial" = "solid"),
    name = "Method"
  ) +
  scale_shape_manual(
    values = c("Multinomial" = 16, "Non-Multinomial" = 17),
    name = "Method"
  ) +
  labs(
    title = "Discrimination Power by Scenario: Symmetric Alpha (10,10,10)",
    subtitle = "Proportion of simulations where correct forecast scores better (lower ES). Colors distinguish methods.",
    x = "n_counts (log scale)",
    y = "Proportion Correct Wins"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    plot.title = element_text(face = "bold", size = 13),
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(face = "bold", size = 9)
  ) +
  guides(
    color = guide_legend(order = 1, override.aes = list(alpha = 1)),
    linetype = guide_legend(order = 1),
    shape = guide_legend(order = 1)
  )

ggsave(file.path(output_dir, "symmetric_alpha_discrimination_vs_n.png"),
       p2, width = 13, height = 9, dpi = 300)
cat("  Plot saved: symmetric_alpha_discrimination_vs_n.png\n\n")

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
  geom_point(size = 2.5) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", alpha = 0.5) +
  facet_wrap(~ scenario_label, ncol = 3, scales = "free_y") +
  scale_x_log10(
    breaks = c(1, 5, 25, 100, 500),
    labels = c(1, 5, 25, 100, 500)
  ) +
  scale_color_manual(
    values = scenario_colors,
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
    title = "Mean Energy Score Differences by Scenario: Symmetric Alpha (10,10,10)",
    subtitle = "ES Difference = Correct - Incorrect (negative values = correct forecast is better)",
    x = "n_counts (log scale)",
    y = "Mean ES Difference"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    plot.title = element_text(face = "bold", size = 13),
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(face = "bold", size = 9)
  ) +
  guides(
    color = guide_legend(order = 1, ncol = 3),
    linetype = guide_legend(order = 2),
    shape = guide_legend(order = 2)
  )

ggsave(file.path(output_dir, "symmetric_alpha_es_differences_vs_n.png"),
       p3, width = 13, height = 9, dpi = 300)
cat("  Plot saved: symmetric_alpha_es_differences_vs_n.png\n\n")

# ============================================================================

# ============================================================================
# Print summary
# ============================================================================
cat("=========================================\n")
cat("SUMMARY\n")
cat("=========================================\n")
cat("Scenarios analyzed:", length(unique(combined_df$scenario_label)), "\n")
cat("Total simulations:", nrow(combined_df), "\n")
cat("\nFiles created:\n")
cat("  Energy Score Plots:\n")
cat("    - symmetric_alpha_energy_scores_vs_n.png\n")
cat("    - symmetric_alpha_discrimination_vs_n.png\n")
cat("    - symmetric_alpha_es_differences_vs_n.png\n")
cat("  Summary Statistics:\n")
cat("    - symmetric_alpha_summary_stats.csv\n\n")

cat("Done!\n")
