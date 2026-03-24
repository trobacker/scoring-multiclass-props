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
    "S1: Mislocation",
    "S2: Moderate Dispersed (5x)",
    "S3: Extreme Dispersed (20x)",
    "S4: Moderate Narrow (5x)",
    "S5: Extreme Narrow (20x)",
    "S6: Ultra Dispersed (100x)",
    "S7: Ultra Narrow (100x)"
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
  "S1: Mislocation" = "#984EA3",                    # purple (location error)
  "S2: Moderate Dispersed (5x)" = "#A6D854",        # light green
  "S3: Extreme Dispersed (20x)" = "#4DAF4A",        # green
  "S6: Ultra Dispersed (100x)" = "#1B7837",         # dark green
  "S4: Moderate Narrow (5x)" = "#FFD92F",           # light orange
  "S5: Extreme Narrow (20x)" = "#E78AC3",           # orange
  "S7: Ultra Narrow (100x)" = "#E41A1C"             # red
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
                 color = scenario_label,
                 linetype = forecast_type,
                 shape = method_label)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_wrap(~ scenario_label, ncol = 3) +
  scale_x_log10(
    breaks = c(1, 5, 25, 100, 500),
    labels = c(1, 5, 25, 100, 500)
  ) +
  scale_color_manual(
    values = scenario_colors,
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
    title = "Mean Energy Scores by Scenario: Symmetric Alpha (10,10,10)",
    subtitle = "Lower scores indicate better forecasts",
    x = "n_counts (log scale)",
    y = "Mean Energy Score"
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
    shape = guide_legend(order = 3)
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
                 color = scenario_label,
                 linetype = method_label,
                 shape = method_label)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
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
    title = "Discrimination Power by Scenario: Symmetric Alpha (10,10,10)",
    subtitle = "Proportion of simulations where correct forecast scores better (lower ES)",
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
    color = guide_legend(order = 1, ncol = 3),
    linetype = guide_legend(order = 2),
    shape = guide_legend(order = 2)
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
# PLOT 4: Dose-Response (Dispersion scenarios only)
# ============================================================================
cat("Creating Plot 4: Dose-Response for Dispersion Errors...\n")

# Filter to dispersion scenarios only (exclude mislocation)
dispersion_scenarios <- summary_stats %>%
  filter(scenario_label != "S1: Mislocation") %>%
  mutate(
    error_type = case_when(
      grepl("Dispersed", scenario_label) ~ "Underconfident",
      grepl("Narrow", scenario_label) ~ "Overconfident"
    ),
    concentration_ratio = case_when(
      grepl("5x", scenario_label) ~ 5,
      grepl("20x", scenario_label) ~ 20,
      grepl("100x", scenario_label) ~ 100
    )
  )

p4 <- ggplot(dispersion_scenarios %>% filter(n_counts %in% c(1, 5, 25, 100)),
             aes(x = concentration_ratio, y = prop_correct_wins,
                 color = factor(n_counts),
                 linetype = method_label,
                 shape = method_label)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "black", alpha = 0.5) +
  scale_x_log10(
    breaks = c(5, 20, 100),
    labels = c("5x", "20x", "100x")
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    labels = scales::percent
  ) +
  scale_color_brewer(
    palette = "Set1",
    name = "n_counts"
  ) +
  scale_linetype_manual(
    values = c("Multinomial" = "solid", "Non-Multinomial" = "dashed"),
    name = "Method"
  ) +
  scale_shape_manual(
    values = c("Multinomial" = 16, "Non-Multinomial" = 17),
    name = "Method"
  ) +
  facet_wrap(~ error_type, ncol = 2) +
  labs(
    title = "Dose-Response: Discrimination vs. Concentration Error Magnitude",
    subtitle = "How does discrimination change with increasing misspecification?",
    x = "Concentration Ratio (log scale)",
    y = "Proportion Correct Wins"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    plot.title = element_text(face = "bold", size = 13),
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(face = "bold")
  ) +
  guides(
    color = guide_legend(order = 1),
    linetype = guide_legend(order = 2),
    shape = guide_legend(order = 2)
  )

ggsave(file.path(output_dir, "symmetric_alpha_dose_response.png"),
       p4, width = 12, height = 6, dpi = 300)
cat("  Plot saved: symmetric_alpha_dose_response.png\n\n")

# ============================================================================
# TERM DECOMPOSITION ANALYSIS
# ============================================================================
cat("Creating term decomposition plots...\n")

# Calculate term summary statistics
term_summary_stats <- combined_df %>%
  group_by(scenario_label, method_label, n_counts) %>%
  summarise(
    mean_term1_correct = mean(term1_correct),
    mean_term1_incorrect = mean(term1_incorrect),
    mean_term2_correct = mean(term2_correct),
    mean_term2_incorrect = mean(term2_incorrect),
    mean_term1_diff = mean(term1_diff),
    mean_term2_diff = mean(term2_diff),
    .groups = "drop"
  )

# PLOT 5: Term Decomposition - Side-by-side bars for selected scenarios
cat("  Plot 5: Term decomposition bar chart\n")

# Filter to a few key n_counts values for clarity
decomp_data <- term_summary_stats %>%
  filter(n_counts %in% c(5, 25, 100)) %>%
  select(scenario_label, method_label, n_counts, mean_term1_correct, mean_term2_correct) %>%
  pivot_longer(cols = c(mean_term1_correct, mean_term2_correct),
               names_to = "term",
               values_to = "value") %>%
  mutate(term = factor(term,
                      levels = c("mean_term1_correct", "mean_term2_correct"),
                      labels = c("Term 1: Distance to observation",
                                "Term 2: 0.5 × Pairwise distance")))

p5 <- ggplot(decomp_data, aes(x = factor(n_counts), y = value, fill = term)) +
  geom_col(position = "dodge", alpha = 0.8, width = 0.7) +
  facet_wrap(~ scenario_label, ncol = 3, scales = "free_y") +
  labs(
    title = "Energy Score Decomposition: Correct Forecasts",
    subtitle = "ES = Term 1 - Term 2. Selected n_counts values shown.",
    x = "n_counts",
    y = "Term Value",
    fill = "Component"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(face = "bold", size = 9)
  ) +
  scale_fill_manual(values = c("#E41A1C", "#377EB8"))

ggsave(file.path(output_dir, "symmetric_alpha_term_decomposition_bars.png"),
       p5, width = 13, height = 9, dpi = 300)

# PLOT 6: Term1 vs Term2 relationship scatter
cat("  Plot 6: Term1 vs Term2 scatter\n")

# Sample data for clarity (use n_counts = 100)
scatter_data <- combined_df %>%
  filter(n_counts == 100) %>%
  select(scenario_label, method_label, term1_correct, term2_correct)

p6 <- ggplot(scatter_data, aes(x = term2_correct, y = term1_correct, color = scenario_label)) +
  geom_point(alpha = 0.3, size = 1) +
  stat_ellipse(level = 0.95, linewidth = 1) +
  facet_wrap(~ method_label, ncol = 2) +
  labs(
    title = "Term 1 vs Term 2 Relationship (n_counts = 100)",
    subtitle = "95% confidence ellipses by scenario. ES = Term 1 - Term 2",
    x = "Term 2: 0.5 × Mean Pairwise Distance",
    y = "Term 1: Mean Distance to Observation",
    color = "Scenario"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    legend.position = "bottom",
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(face = "bold")
  ) +
  scale_color_manual(values = scenario_colors) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 3), ncol = 3))

ggsave(file.path(output_dir, "symmetric_alpha_term_scatter.png"),
       p6, width = 12, height = 6, dpi = 300)

# PLOT 7: Term differences vs n_counts
cat("  Plot 7: Term differences vs n_counts\n")

# Reshape for plotting
term_diff_long <- term_summary_stats %>%
  pivot_longer(cols = c(mean_term1_diff, mean_term2_diff),
               names_to = "term",
               values_to = "mean_diff") %>%
  mutate(term = factor(term,
                      levels = c("mean_term1_diff", "mean_term2_diff"),
                      labels = c("Term 1 difference", "Term 2 difference")))

p7 <- ggplot(term_diff_long, aes(x = n_counts, y = mean_diff,
                                 color = scenario_label,
                                 linetype = method_label,
                                 shape = method_label)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", alpha = 0.5) +
  facet_wrap(~ term, ncol = 2, scales = "free_y") +
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
    title = "Term Differences by Scenario (Correct - Incorrect)",
    subtitle = "How term1 and term2 contribute to overall ES difference",
    x = "n_counts (log scale)",
    y = "Mean Term Difference"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    legend.position = "bottom",
    legend.box = "vertical",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(face = "bold")
  ) +
  guides(
    color = guide_legend(order = 1, ncol = 3),
    linetype = guide_legend(order = 2),
    shape = guide_legend(order = 2)
  )

ggsave(file.path(output_dir, "symmetric_alpha_term_differences_vs_n.png"),
       p7, width = 13, height = 6, dpi = 300)

cat("  Term decomposition plots saved\n\n")

# Save term summary statistics
write.csv(term_summary_stats,
          file.path(output_dir, "symmetric_alpha_term_summary_stats.csv"),
          row.names = FALSE)
cat("  Term summary statistics saved\n\n")

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
cat("    - symmetric_alpha_dose_response.png\n")
cat("  Term Decomposition Plots:\n")
cat("    - symmetric_alpha_term_decomposition_bars.png\n")
cat("    - symmetric_alpha_term_scatter.png\n")
cat("    - symmetric_alpha_term_differences_vs_n.png\n")
cat("  Summary Statistics:\n")
cat("    - symmetric_alpha_summary_stats.csv\n")
cat("    - symmetric_alpha_term_summary_stats.csv\n\n")

cat("Done!\n")
