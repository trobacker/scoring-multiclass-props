#!/usr/bin/env Rscript
# Combine scenarios for small n_counts (1, 2, 3, 4)
# Creates a single plot comparing all 4 scenarios

library(ggplot2)
library(tidyr)
library(dplyr)

# Define paths
data_dir <- "./results/data/alpha_2_3_15"
output_dir <- "./results/plotting_code"

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

cat("Total rows loaded:", nrow(combined_df), "\n")

# Filter to small n_counts only
combined_df <- combined_df %>%
  filter(n_counts %in% c(1, 2, 3, 4)) %>%
  mutate(
    method_label = case_when(
      method == "mn" ~ "Multinomial Sampling",
      method == "non-mn" ~ "Non-Multinomial"
    )
  )

cat("Rows after filtering to n_counts 1-4:", nrow(combined_df), "\n")

# Pivot longer for plotting (es_correct vs es_incorrect)
df_long <- combined_df %>%
  pivot_longer(
    cols = c(es_correct, es_incorrect),
    names_to = "forecast_type",
    values_to = "energy_score"
  ) %>%
  mutate(
    # Clean up labels
    forecast_type = case_when(
      forecast_type == "es_correct" ~ "Correct Forecast",
      forecast_type == "es_incorrect" ~ "Incorrect Forecast"
    )
  )

cat("Creating plot...\n")

# Create the combined plot
p <- ggplot(df_long, aes(x = factor(n_counts), y = energy_score,
                         color = scenario_label, shape = method_label)) +
  geom_point(alpha = 0.4, size = 1.5,
             position = position_jitterdodge(jitter.width = 0.2,
                                            dodge.width = 0.7)) +
  facet_wrap(~ forecast_type, ncol = 2) +
  scale_shape_manual(
    values = c("Multinomial Sampling" = 16,      # filled circle
               "Non-Multinomial" = 4),            # x shape
    name = "Method"
  ) +
  scale_color_manual(
    values = c(
      "Scenario 1: Misspecified" = "#E41A1C",           # red
      "Scenario 2: Narrow/Calibrated" = "#377EB8",      # blue
      "Scenario 3: Dispersed/Calibrated" = "#4DAF4A",   # green
      "Scenario 4: Less Misspecified" = "#984EA3"       # purple
    ),
    name = "Scenario"
  ) +
  labs(
    title = "Energy Scores Across Scenarios: Small Sample Sizes",
    subtitle = "Comparing Multinomial vs Non-Multinomial sampling for n_counts = 1-4",
    x = "n_counts (number of sequences)",
    y = "Energy Score"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.margin = margin(t = 5),
    plot.title = element_text(face = "bold", size = 14),
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(face = "bold", size = 11)
  ) +
  guides(
    color = guide_legend(ncol = 2, override.aes = list(size = 3, alpha = 1)),
    shape = guide_legend(ncol = 2, override.aes = list(size = 3, alpha = 1))
  )

# Save the plot
output_file <- file.path(output_dir, "combined_scenarios_small_n_comparison.png")
ggsave(output_file, plot = p, width = 12, height = 7, dpi = 300)

cat("\nPlot saved to:", output_file, "\n")

# Also create a version showing differences (es_diff)
cat("\nCreating difference plot...\n")

# Calculate means for each combination
diff_means <- combined_df %>%
  group_by(scenario_label, method_label, n_counts) %>%
  summarise(
    mean_diff = mean(es_diff),
    .groups = "drop"
  )

p_diff <- ggplot(combined_df, aes(x = factor(n_counts), y = es_diff,
                                   color = scenario_label, shape = method_label)) +
  geom_point(alpha = 0.4, size = 1.5,
             position = position_jitterdodge(jitter.width = 0.2,
                                            dodge.width = 0.7)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
  # Add mean points - make them stand out
  geom_point(data = diff_means,
             aes(x = factor(n_counts), y = mean_diff,
                 color = scenario_label, shape = method_label),
             size = 5, alpha = 1, stroke = 1.5,
             position = position_dodge(width = 0.7)) +
  scale_shape_manual(
    values = c("Multinomial Sampling" = 16,
               "Non-Multinomial" = 4),
    name = "Method"
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
  labs(
    title = "Energy Score Differences Across Scenarios: Small Sample Sizes",
    subtitle = "ES Difference = Correct - Incorrect (negative = correct better). Large symbols = mean difference.",
    x = "n_counts (number of sequences)",
    y = "ES Difference (Correct - Incorrect)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.margin = margin(t = 5),
    plot.title = element_text(face = "bold", size = 14)
  ) +
  guides(
    color = guide_legend(ncol = 2, override.aes = list(size = 3, alpha = 1)),
    shape = guide_legend(ncol = 2, override.aes = list(size = 3, alpha = 1))
  )

output_file_diff <- file.path(output_dir, "combined_scenarios_small_n_differences.png")
ggsave(output_file_diff, plot = p_diff, width = 10, height = 7, dpi = 300)

cat("Difference plot saved to:", output_file_diff, "\n")

# Summary statistics
cat("\n=== Summary Statistics ===\n")
summary_stats <- df_long %>%
  group_by(scenario_label, method_label, n_counts, forecast_type) %>%
  summarise(
    mean_es = mean(energy_score),
    median_es = median(energy_score),
    sd_es = sd(energy_score),
    .groups = "drop"
  )

print(summary_stats)

cat("\n=== Difference Statistics ===\n")
diff_stats <- combined_df %>%
  group_by(scenario_label, method_label, n_counts) %>%
  summarise(
    mean_diff = mean(es_diff),
    median_diff = median(es_diff),
    prop_correct_better = mean(es_diff < 0),
    .groups = "drop"
  ) %>%
  arrange(scenario_label, method_label, n_counts)

print(diff_stats)

cat("\nDone!\n")
