#!/usr/bin/env Rscript

# Term Decomposition Analysis Visualizations
# This script creates plots to highlight the roles of term1 and term2 in energy score discrimination

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(viridisLite)

# Set output directory
output_dir <- "results/term_decomposition_analysis"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# Load all simulation data
# ============================================================================

cat("Loading simulation data...\n")

# Symmetric scenarios (α = 10, 10, 10)
symmetric_scenarios <- list(
  scenario1_mislocation = "results/data/alpha_10_10_10/scenario1_mislocation_symmetric_alpha_10_10_10.rds",
  scenario2_moderate_dispersed = "results/data/alpha_10_10_10/scenario2_moderate_dispersed_alpha_10_10_10.rds",
  scenario3_extreme_dispersed = "results/data/alpha_10_10_10/scenario3_extreme_dispersed_alpha_10_10_10.rds",
  scenario4_moderate_narrow = "results/data/alpha_10_10_10/scenario4_moderate_narrow_alpha_10_10_10.rds",
  scenario5_extreme_narrow = "results/data/alpha_10_10_10/scenario5_extreme_narrow_alpha_10_10_10.rds",
  scenario6_ultra_dispersed = "results/data/alpha_10_10_10/scenario6_ultra_dispersed_alpha_10_10_10.rds",
  scenario7_ultra_narrow = "results/data/alpha_10_10_10/scenario7_ultra_narrow_alpha_10_10_10.rds"
)

# Asymmetric scenarios (α = 2, 3, 15)
asymmetric_scenarios <- list(
  scenario1_misspecified = "results/data/alpha_2_3_15/scenario1_misspecified_alpha_2_3_15.rds",
  scenario2_narrow_calibrated = "results/data/alpha_2_3_15/scenario2_narrow_calibrated_alpha_2_3_15.rds",
  scenario3_dispersed_calibrated = "results/data/alpha_2_3_15/scenario3_dispersed_calibrated_alpha_2_3_15.rds",
  scenario4_less_misspecified = "results/data/alpha_2_3_15/scenario4_less_misspecified_alpha_2_3_15.rds"
)

# Load all data
all_data <- bind_rows(
  lapply(names(symmetric_scenarios), function(name) {
    readRDS(symmetric_scenarios[[name]]) %>%
      mutate(
        scenario_short = name,
        alpha_type = "symmetric"
      )
  }),
  lapply(names(asymmetric_scenarios), function(name) {
    readRDS(asymmetric_scenarios[[name]]) %>%
      mutate(
        scenario_short = name,
        alpha_type = "asymmetric"
      )
  })
)

cat("Data loaded. Total rows:", nrow(all_data), "\n")

# ============================================================================
# PLOT 1: Term Decomposition Line Plots
# ============================================================================

cat("\n=== Creating Plot 1: Term Decomposition Line Plots ===\n")

# Calculate summary statistics by scenario, N, and method
term_summary <- all_data %>%
  group_by(scenario_short, alpha_type, n_counts, method, alpha_true, alpha_wrong) %>%
  summarise(
    mean_es_correct = mean(es_correct, na.rm = TRUE),
    mean_es_incorrect = mean(es_incorrect, na.rm = TRUE),
    mean_term1_correct = mean(term1_correct, na.rm = TRUE),
    mean_term1_incorrect = mean(term1_incorrect, na.rm = TRUE),
    mean_term2_correct = mean(term2_correct, na.rm = TRUE),
    mean_term2_incorrect = mean(term2_incorrect, na.rm = TRUE),
    .groups = 'drop'
  )

# Reshape for plotting
term_long <- term_summary %>%
  pivot_longer(
    cols = c(mean_es_correct, mean_es_incorrect,
             mean_term1_correct, mean_term1_incorrect,
             mean_term2_correct, mean_term2_incorrect),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    term_type = case_when(
      grepl("term1", metric) ~ "Term 1: E[||y - Y||] (distance from observation)",
      grepl("term2", metric) ~ "Term 2: 0.5*E[||Y - Y'||] (distance between forecast samples)",
      grepl("es", metric) ~ "Energy Score (ES = term1 - term2)"
    ),
    forecast_type = ifelse(grepl("_correct$", metric), "Correct Forecast", "Incorrect Forecast")
  )

# Plot for symmetric scenarios
symmetric_plot <- term_long %>%
  filter(alpha_type == "symmetric") %>%
  mutate(scenario_label = gsub("_calibrated", "", scenario_short)) %>%
  ggplot(aes(x = n_counts, y = value, color = forecast_type, linetype = method,
             group = interaction(forecast_type, method))) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  facet_grid(term_type ~ scenario_label, scales = "free_y") +
  scale_x_log10(breaks = c(1, 2, 5, 10, 25, 50, 100, 500)) +
  scale_color_manual(values = c("Correct Forecast" = "#1f77b4", "Incorrect Forecast" = "#d62728")) +
  scale_linetype_manual(values = c("mn" = "solid", "non-mn" = "dashed")) +
  labs(
    title = "Energy Score Term Decomposition: Symmetric Scenarios (α = 10, 10, 10)",
    subtitle = "ES = term1 - term2. Note: different y-axis scales per row. Lower ES is better.",
    x = "N (sample size, log scale)",
    y = "Mean Value",
    color = "Forecast Type",
    linetype = "Method"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  file.path(output_dir, "plot1_term_decomposition_symmetric.png"),
  symmetric_plot,
  width = 16,
  height = 10,
  dpi = 300
)

cat("  Saved: plot1_term_decomposition_symmetric.png\n")

# Plot for asymmetric scenarios
asymmetric_plot <- term_long %>%
  filter(alpha_type == "asymmetric") %>%
  mutate(scenario_label = gsub("_calibrated", "", scenario_short)) %>%
  ggplot(aes(x = n_counts, y = value, color = forecast_type, linetype = method,
             group = interaction(forecast_type, method))) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  facet_grid(term_type ~ scenario_label, scales = "free_y") +
  scale_x_log10(breaks = c(1, 2, 5, 10, 25, 50, 100, 500)) +
  scale_color_manual(values = c("Correct Forecast" = "#1f77b4", "Incorrect Forecast" = "#d62728")) +
  scale_linetype_manual(values = c("mn" = "solid", "non-mn" = "dashed")) +
  labs(
    title = "Energy Score Term Decomposition: Asymmetric Scenarios (α = 2, 3, 15)",
    subtitle = "ES = term1 - term2. Note: different y-axis scales per row. Lower ES is better.",
    x = "N (sample size, log scale)",
    y = "Mean Value",
    color = "Forecast Type",
    linetype = "Method"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  file.path(output_dir, "plot1_term_decomposition_asymmetric.png"),
  asymmetric_plot,
  width = 14,
  height = 10,
  dpi = 300
)

cat("  Saved: plot1_term_decomposition_asymmetric.png\n")

# ============================================================================
# PLOT 3: Term Scatter Plots (Individual simulation level)
# ============================================================================

cat("\n=== Creating Plot 3: Term Scatter Plots ===\n")

# Focus on key scenarios and N values
scatter_scenarios <- c("scenario5_extreme_narrow", "scenario3_extreme_dispersed", "scenario1_mislocation")
scatter_n_values <- c(1, 5, 10, 100, 500)

scatter_data <- all_data %>%
  filter(
    scenario_short %in% scatter_scenarios,
    n_counts %in% scatter_n_values
  ) %>%
  mutate(
    correct_wins = es_correct < es_incorrect,
    n_label = paste0("N = ", n_counts)
  )

# Create scatter plot
scatter_plot <- scatter_data %>%
  ggplot(aes(x = term1_diff, y = term2_diff, color = correct_wins)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_point(alpha = 0.3, size = 0.8) +
  facet_grid(
    factor(n_label, levels = paste0("N = ", scatter_n_values)) ~
    factor(scenario_short, levels = scatter_scenarios) + method,
    scales = "free"
  ) +
  scale_color_manual(
    values = c("TRUE" = "#2ca02c", "FALSE" = "#d62728"),
    labels = c("TRUE" = "Correct Wins", "FALSE" = "Incorrect Wins")
  ) +
  labs(
    title = "Term1 vs Term2 Differences: Individual Simulation Runs",
    subtitle = "term_diff = correct - incorrect. Positive term2_diff favors correct forecast",
    x = "term1_diff (correct - incorrect)",
    y = "term2_diff (correct - incorrect)",
    color = "Outcome"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 7),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y = element_text(size = 7)
  )

ggsave(
  file.path(output_dir, "plot3_term_scatter.png"),
  scatter_plot,
  width = 14,
  height = 12,
  dpi = 300
)

cat("  Saved: plot3_term_scatter.png\n")

# ============================================================================
# PLOT 6: Variance vs Term2 Validation Plot
# ============================================================================

cat("\n=== Creating Plot 6: Variance vs Term2 Validation ===\n")

# For this plot, we need to understand that term2 should be proportional to
# the spread of the forecast distribution. Since term2 is constant for each
# forecast type × N × method combination, we'll plot those aggregate values

variance_term2_data <- all_data %>%
  group_by(scenario_short, alpha_type, n_counts, method, alpha_true, alpha_wrong) %>%
  summarise(
    term2_correct = first(term2_correct),
    term2_incorrect = first(term2_incorrect),
    .groups = 'drop'
  ) %>%
  pivot_longer(
    cols = c(term2_correct, term2_incorrect),
    names_to = "forecast_type",
    values_to = "term2"
  ) %>%
  mutate(
    forecast_type = ifelse(forecast_type == "term2_correct", "Correct", "Incorrect"),
    # Parse alpha values to estimate "concentration" (sum of alphas)
    alpha_str = ifelse(forecast_type == "Correct", alpha_true, alpha_wrong),
    scenario_type = case_when(
      grepl("dispersed", scenario_short) ~ "Dispersed",
      grepl("narrow", scenario_short) ~ "Narrow",
      TRUE ~ "Mislocation"
    )
  )

# Plot term2 vs N, colored by scenario type and forecast type
variance_plot <- variance_term2_data %>%
  ggplot(aes(x = n_counts, y = term2, color = interaction(scenario_type, forecast_type),
             shape = method)) +
  geom_line(aes(group = interaction(scenario_short, forecast_type, method)), alpha = 0.6) +
  geom_point(size = 2.5, alpha = 0.8) +
  facet_wrap(~ alpha_type, scales = "free_y", ncol = 1) +
  scale_x_log10(breaks = c(1, 2, 5, 10, 25, 50, 100, 500)) +
  scale_shape_manual(values = c("mn" = 16, "non-mn" = 17)) +
  labs(
    title = "Term2 Values Across Scenarios and Sample Sizes",
    subtitle = "Term2 = -0.5*E[||Y - Y'||] reflects forecast spread/variance\nMore negative = higher variance. Shows convergence of MN → non-MN as N increases",
    x = "N (sample size, log scale)",
    y = "Term2 Value",
    color = "Scenario × Forecast",
    shape = "Method"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 8)
  )

ggsave(
  file.path(output_dir, "plot6_variance_term2_validation.png"),
  variance_plot,
  width = 12,
  height = 10,
  dpi = 300
)

cat("  Saved: plot6_variance_term2_validation.png\n")

# ============================================================================
# PLOT 7: Term Convergence Plot (Multinomial vs Non-Multinomial)
# ============================================================================

cat("\n=== Creating Plot 7: Term Convergence Plot ===\n")

# Calculate the difference between MN and non-MN term2 values
convergence_data <- all_data %>%
  group_by(scenario_short, alpha_type, n_counts, method) %>%
  summarise(
    term2_correct = first(term2_correct),
    term2_incorrect = first(term2_incorrect),
    .groups = 'drop'
  ) %>%
  pivot_wider(
    names_from = method,
    values_from = c(term2_correct, term2_incorrect)
  ) %>%
  mutate(
    diff_correct = abs(term2_correct_mn - `term2_correct_non-mn`),
    diff_incorrect = abs(term2_incorrect_mn - `term2_incorrect_non-mn`)
  ) %>%
  select(scenario_short, alpha_type, n_counts, diff_correct, diff_incorrect) %>%
  pivot_longer(
    cols = c(diff_correct, diff_incorrect),
    names_to = "forecast_type",
    values_to = "mn_nonmn_diff"
  ) %>%
  mutate(
    forecast_type = ifelse(forecast_type == "diff_correct", "Correct Forecast", "Incorrect Forecast")
  )

convergence_plot <- convergence_data %>%
  ggplot(aes(x = n_counts, y = mn_nonmn_diff, color = scenario_short, linetype = forecast_type)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_wrap(~ alpha_type, ncol = 1, scales = "free_y") +
  scale_x_log10(breaks = c(1, 2, 5, 10, 25, 50, 100, 500)) +
  scale_y_log10() +
  scale_linetype_manual(values = c("Correct Forecast" = "solid", "Incorrect Forecast" = "dashed")) +
  labs(
    title = "Convergence of Multinomial to Non-Multinomial Term2",
    subtitle = "Showing |term2_MN - term2_nonMN| decreasing as N increases",
    x = "N (sample size, log scale)",
    y = "|term2_MN - term2_nonMN| (log scale)",
    color = "Scenario",
    linetype = "Forecast Type"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 8)
  )

ggsave(
  file.path(output_dir, "plot7_term_convergence.png"),
  convergence_plot,
  width = 12,
  height = 10,
  dpi = 300
)

cat("  Saved: plot7_term_convergence.png\n")

# ============================================================================
# PLOT 8: Discrimination Decomposition Table
# ============================================================================

cat("\n=== Creating Plot 8: Discrimination Decomposition Table ===\n")

# Create comprehensive table
decomp_table <- all_data %>%
  group_by(scenario_short, alpha_type, n_counts, method, alpha_true, alpha_wrong) %>%
  summarise(
    discrimination_pct = mean(es_correct < es_incorrect) * 100,
    mean_es_diff = mean(es_diff),
    mean_term1_diff = mean(term1_diff),
    mean_term2_diff = mean(term2_diff),
    .groups = 'drop'
  ) %>%
  mutate(
    # Calculate contribution ratio: how much of es_diff comes from term2
    term2_contrib_ratio = abs(mean_term2_diff) / (abs(mean_term1_diff) + abs(mean_term2_diff)),
    # Format for readability
    across(c(mean_es_diff, mean_term1_diff, mean_term2_diff), ~round(., 4)),
    across(c(discrimination_pct, term2_contrib_ratio), ~round(., 2))
  ) %>%
  arrange(alpha_type, scenario_short, method, n_counts)

# Save as CSV
write.csv(
  decomp_table,
  file.path(output_dir, "plot8_decomposition_table.csv"),
  row.names = FALSE
)

cat("  Saved: plot8_decomposition_table.csv\n")

# Also create a visual table for key scenarios and N values
key_table <- decomp_table %>%
  filter(
    scenario_short %in% c("scenario5_extreme_narrow", "scenario3_extreme_dispersed"),
    n_counts %in% c(1, 5, 10, 100, 500)
  ) %>%
  select(scenario_short, n_counts, method, discrimination_pct, mean_es_diff,
         mean_term1_diff, mean_term2_diff, term2_contrib_ratio)

# Create a heatmap-style visualization of the table
table_long <- key_table %>%
  pivot_longer(
    cols = c(discrimination_pct, term2_contrib_ratio),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric_label = ifelse(metric == "discrimination_pct", "Discrimination %", "Term2 Contribution")
  )

table_plot <- table_long %>%
  ggplot(aes(x = factor(n_counts), y = interaction(scenario_short, method), fill = value)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(value, 1)), size = 3) +
  facet_wrap(~ metric_label, ncol = 2, scales = "free_x") +
  scale_fill_viridis_c(option = "plasma") +
  labs(
    title = "Discrimination and Term2 Contribution: Key Scenarios",
    subtitle = "Term2 contribution = |term2_diff| / (|term1_diff| + |term2_diff|)",
    x = "N (sample size)",
    y = "Scenario × Method",
    fill = "Value"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0),
    axis.text.y = element_text(size = 8),
    legend.position = "right"
  )

ggsave(
  file.path(output_dir, "plot8_decomposition_table_heatmap.png"),
  table_plot,
  width = 10,
  height = 6,
  dpi = 300
)

cat("  Saved: plot8_decomposition_table_heatmap.png\n")

# ============================================================================
# PLOT 9: Quantization Effect Visualization
# ============================================================================

cat("\n=== Creating Plot 9: Quantization Effect Visualization ===\n")

# This requires us to generate actual Dirichlet samples and their multinomial
# counterparts. We'll do this for the narrow scenario as an example.

source("src/helper-functions.R")

set.seed(42)

# Parameters
n_samples <- 100
n_values_to_show <- c(1, 5, 10, 100)
alpha_correct <- c(10, 10, 10)
alpha_narrow <- c(200, 200, 200)

# Function to create quantization plot for one N value
create_quantization_plot <- function(N, alpha_correct, alpha_narrow, n_samples) {

  # Generate Dirichlet samples
  correct_samples <- brms::rdirichlet(n_samples, alpha_correct)
  narrow_samples <- brms::rdirichlet(n_samples, alpha_narrow)

  # For multinomial resampling, we need to pick one observation
  # Let's use the mean as the "true" observation
  true_theta <- alpha_correct / sum(alpha_correct)

  # Generate multinomial samples from each Dirichlet sample
  correct_mn_samples <- t(sapply(1:n_samples, function(i) {
    counts <- mc2d::rmultinomial(1, N, correct_samples[i, ])
    counts / N  # Convert to proportions
  }))

  narrow_mn_samples <- t(sapply(1:n_samples, function(i) {
    counts <- mc2d::rmultinomial(1, N, narrow_samples[i, ])
    counts / N  # Convert to proportions
  }))

  # Combine into data frames
  df_correct <- bind_rows(
    data.frame(correct_samples, type = "Original Dirichlet", forecast = "Correct"),
    data.frame(correct_mn_samples, type = "Multinomial Resampled", forecast = "Correct")
  )
  colnames(df_correct)[1:3] <- c("p1", "p2", "p3")

  df_narrow <- bind_rows(
    data.frame(narrow_samples, type = "Original Dirichlet", forecast = "Narrow (α=200,200,200)"),
    data.frame(narrow_mn_samples, type = "Multinomial Resampled", forecast = "Narrow (α=200,200,200)")
  )
  colnames(df_narrow)[1:3] <- c("p1", "p2", "p3")

  df_all <- bind_rows(df_correct, df_narrow)

  # Create ternary plot using basic plotting
  # We'll convert to 2D coordinates for ggplot
  df_all$x <- df_all$p2 + df_all$p3 / 2
  df_all$y <- df_all$p3 * sqrt(3) / 2

  p <- ggplot(df_all, aes(x = x, y = y, color = type, shape = type)) +
    geom_point(alpha = 0.5, size = 2) +
    facet_wrap(~ forecast, ncol = 2) +
    scale_color_manual(values = c("Original Dirichlet" = "#1f77b4", "Multinomial Resampled" = "#ff7f0e")) +
    scale_shape_manual(values = c("Original Dirichlet" = 16, "Multinomial Resampled" = 17)) +
    coord_fixed() +
    labs(
      title = sprintf("Quantization Effect at N = %d", N),
      subtitle = "Blue circles: Original Dirichlet samples | Orange triangles: After multinomial resampling",
      color = "Sample Type",
      shape = "Sample Type"
    ) +
    theme_void() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, size = 9)
    )

  return(p)
}

# Create plots for each N value
quantization_plots <- lapply(n_values_to_show, function(N) {
  create_quantization_plot(N, alpha_correct, alpha_narrow, n_samples)
})

# Combine into a single figure
combined_quantization <- wrap_plots(quantization_plots, ncol = 2)

ggsave(
  file.path(output_dir, "plot9_quantization_effect.png"),
  combined_quantization,
  width = 14,
  height = 12,
  dpi = 300
)

cat("  Saved: plot9_quantization_effect.png\n")

# ============================================================================
# Summary
# ============================================================================

cat("\n=========================================\n")
cat("SUMMARY\n")
cat("=========================================\n\n")

cat("All plots saved to:", output_dir, "\n\n")

cat("Files created:\n")
cat("  1. plot1_term_decomposition_symmetric.png - Term decomposition for symmetric scenarios\n")
cat("  2. plot1_term_decomposition_asymmetric.png - Term decomposition for asymmetric scenarios\n")
cat("  3. plot3_term_scatter.png - Scatter plots of term1_diff vs term2_diff\n")
cat("  4. plot6_variance_term2_validation.png - Term2 values across scenarios\n")
cat("  5. plot7_term_convergence.png - Convergence of MN to non-MN term2\n")
cat("  6. plot8_decomposition_table.csv - Full decomposition table (CSV)\n")
cat("  7. plot8_decomposition_table_heatmap.png - Visual heatmap of key metrics\n")
cat("  8. plot9_quantization_effect.png - Visualization of multinomial quantization\n\n")

cat("Done!\n")
