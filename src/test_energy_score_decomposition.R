#!/usr/bin/env Rscript
# Comprehensive test of manual_energy_score() function
# Tests ES decomposition across multiple forecast scenarios
# Compares custom implementation with scoringRules::es_sample()

library(brms)
library(dplyr)
library(ggplot2)
library(tidyr)
library(scoringRules)

# Source helper functions
source("src/helper-functions.R")

cat("=========================================\n")
cat("ENERGY SCORE DECOMPOSITION ANALYSIS\n")
cat("=========================================\n\n")

# Set seed for reproducibility
set.seed(42)

# Simulation parameters
n_sims <- 500  # Number of simulations per scenario
n_samples <- 100  # Number of forecast samples (thetas)
alpha_true <- c(5, 2, 3)  # True distribution

cat("Simulation settings:\n")
cat("  Simulations per scenario:", n_sims, "\n")
cat("  Forecast samples:", n_samples, "\n")
cat("  True alpha:", paste(alpha_true, collapse = ", "), "\n\n")

# Define forecast scenarios to test
scenarios <- list(
  list(
    name = "Ground Truth",
    alpha = alpha_true,
    description = "Perfect forecast (same as true distribution)"
  ),
  list(
    name = "Slight Mislocation",
    alpha = c(4.5, 2.5, 3),
    description = "Small location error"
  ),
  list(
    name = "Moderate Mislocation",
    alpha = c(3, 4, 3),
    description = "Moderate location error"
  ),
  list(
    name = "Severe Mislocation",
    alpha = c(2, 8, 1),
    description = "Large location error"
  ),
  list(
    name = "Overconfident (5x)",
    alpha = alpha_true * 5,  # Narrower distribution
    description = "Too narrow/concentrated"
  ),
  list(
    name = "Underconfident (5x)",
    alpha = alpha_true / 5,  # Wider distribution
    description = "Too dispersed"
  )
)

cat("Testing", length(scenarios), "forecast scenarios:\n")
for (i in seq_along(scenarios)) {
  cat(sprintf("  %d. %s: %s\n", i, scenarios[[i]]$name, scenarios[[i]]$description))
}
cat("\n")

# Run simulations for all scenarios
cat("Running simulations...\n")
results_list <- list()

for (s in seq_along(scenarios)) {
  scenario <- scenarios[[s]]
  cat(sprintf("  Scenario %d/%d: %s\n", s, length(scenarios), scenario$name))

  # Pre-generate forecast samples for this scenario
  # rdirichlet(n, alpha) returns n×d matrix, we need d×n for manual_energy_score
  forecast_matrix <- t(rdirichlet(n_samples, alpha = scenario$alpha))

  # Storage for this scenario
  scenario_results <- data.frame(
    scenario = character(n_sims),
    sim_id = integer(n_sims),
    es_custom = numeric(n_sims),
    es_scoringRules = numeric(n_sims),
    term1 = numeric(n_sims),
    term2 = numeric(n_sims),
    es_diff = numeric(n_sims),
    stringsAsFactors = FALSE
  )

  # Run simulations
  for (i in 1:n_sims) {
    # Generate observation from true distribution
    y_obs <- as.vector(rdirichlet(1, alpha = alpha_true))

    # Calculate ES using both methods
    custom_result <- manual_energy_score(y = y_obs, dat = forecast_matrix)
    scoringRules_es <- es_sample(y = y_obs, dat = forecast_matrix)

    # Store results
    scenario_results$scenario[i] <- scenario$name
    scenario_results$sim_id[i] <- i
    scenario_results$es_custom[i] <- custom_result$ES
    scenario_results$es_scoringRules[i] <- scoringRules_es
    scenario_results$term1[i] <- custom_result$term1
    scenario_results$term2[i] <- custom_result$term2
    scenario_results$es_diff[i] <- custom_result$ES - scoringRules_es
  }

  results_list[[s]] <- scenario_results
}

# Combine all results
all_results <- bind_rows(results_list)

cat("  Complete!\n\n")

# Add scenario factor with ordered levels
all_results$scenario <- factor(all_results$scenario,
                               levels = sapply(scenarios, function(x) x$name))

# ============================================================================
# ANALYSIS 1: Verify custom ES matches scoringRules::es_sample
# ============================================================================

cat("=========================================\n")
cat("VERIFICATION: Custom ES vs scoringRules\n")
cat("=========================================\n\n")

verification_summary <- all_results %>%
  group_by(scenario) %>%
  summarise(
    mean_custom = mean(es_custom),
    mean_scoringRules = mean(es_scoringRules),
    mean_diff = mean(es_diff),
    max_abs_diff = max(abs(es_diff)),
    .groups = "drop"
  )

print(verification_summary)
cat("\n")

# Check if differences are within floating-point tolerance
tolerance <- 1e-10
max_diff <- max(abs(all_results$es_diff))

if (max_diff < tolerance) {
  cat(sprintf("✓ SUCCESS: Maximum difference (%.2e) is within tolerance (%.2e)\n",
              max_diff, tolerance))
  cat("  Custom ES function matches scoringRules::es_sample()\n\n")
} else {
  cat(sprintf("⚠ WARNING: Maximum difference (%.2e) exceeds tolerance (%.2e)\n",
              max_diff, tolerance))
  cat("  This may indicate a problem with the custom implementation\n\n")
}

# ============================================================================
# ANALYSIS 2: Energy Score Decomposition
# ============================================================================

cat("=========================================\n")
cat("ENERGY SCORE DECOMPOSITION ANALYSIS\n")
cat("=========================================\n\n")

decomposition_summary <- all_results %>%
  group_by(scenario) %>%
  summarise(
    mean_es = mean(es_custom),
    mean_term1 = mean(term1),
    mean_term2 = mean(term2),
    sd_es = sd(es_custom),
    term1_pct = mean(term1) / (mean(term1) + mean(term2)) * 100,
    term2_pct = mean(term2) / (mean(term1) + mean(term2)) * 100,
    .groups = "drop"
  )

cat("Mean values by scenario:\n")
print(decomposition_summary)
cat("\n")

# ============================================================================
# PLOTTING
# ============================================================================

cat("Creating visualizations...\n")

# Create output directory
output_dir <- "results/test_results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# PLOT 1: Verification - Custom ES vs scoringRules ES
cat("  Plot 1: ES comparison (custom vs scoringRules)\n")

p1 <- ggplot(all_results, aes(x = es_scoringRules, y = es_custom)) +
  geom_point(alpha = 0.3, color = "steelblue") +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
  facet_wrap(~ scenario, scales = "free", ncol = 3) +
  labs(
    title = "Verification: Custom ES vs scoringRules::es_sample()",
    subtitle = "Points should fall on red diagonal line (y = x)",
    x = "scoringRules::es_sample()",
    y = "manual_energy_score()$ES"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold"),
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(face = "bold")
  )

ggsave(file.path(output_dir, "es_verification_plot.png"),
       p1, width = 12, height = 8, dpi = 300)

# PLOT 2: Energy Score by Scenario
cat("  Plot 2: ES distribution by scenario\n")

p2 <- ggplot(all_results, aes(x = scenario, y = es_custom, fill = scenario)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.3) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
  labs(
    title = "Energy Score Distribution by Forecast Scenario",
    subtitle = "Lower scores indicate better forecasts. White diamonds show mean values.",
    x = "Forecast Scenario",
    y = "Energy Score"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  scale_fill_brewer(palette = "Set2")

ggsave(file.path(output_dir, "es_by_scenario_boxplot.png"),
       p2, width = 10, height = 6, dpi = 300)

# PLOT 3: Term Decomposition - Stacked Bar
cat("  Plot 3: Term decomposition (stacked bars)\n")

decomp_long <- decomposition_summary %>%
  select(scenario, mean_term1, mean_term2) %>%
  pivot_longer(cols = c(mean_term1, mean_term2),
               names_to = "term",
               values_to = "value") %>%
  mutate(term = factor(term,
                      levels = c("mean_term1", "mean_term2"),
                      labels = c("Term 1: Mean distance to observation",
                                "Term 2: 0.5 × Mean pairwise distance")))

p3 <- ggplot(decomp_long, aes(x = scenario, y = value, fill = term)) +
  geom_col(position = "dodge", alpha = 0.8, width = 0.7) +
  labs(
    title = "Energy Score Decomposition by Scenario",
    subtitle = "ES = Term 1 - Term 2. Bars show raw values, not percentages.",
    x = "Forecast Scenario",
    y = "Term Value",
    fill = "Component"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  ) +
  scale_fill_manual(values = c("#E41A1C", "#377EB8"))

ggsave(file.path(output_dir, "es_decomposition_barplot.png"),
       p3, width = 10, height = 6, dpi = 300)

# PLOT 4: Term Decomposition - Individual Terms
cat("  Plot 4: Individual term distributions\n")

terms_long <- all_results %>%
  select(scenario, sim_id, term1, term2) %>%
  pivot_longer(cols = c(term1, term2),
               names_to = "term",
               values_to = "value") %>%
  mutate(term = factor(term,
                      levels = c("term1", "term2"),
                      labels = c("Term 1", "Term 2")))

p4 <- ggplot(terms_long, aes(x = scenario, y = value, fill = term)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.3) +
  facet_wrap(~ term, ncol = 1, scales = "free_y") +
  labs(
    title = "Energy Score Term Distributions by Scenario",
    subtitle = "ES = Term 1 - Term 2 (sharpness penalty minus diversity credit)",
    x = "Forecast Scenario",
    y = "Term Value",
    fill = "Term"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(face = "bold")
  ) +
  scale_fill_manual(values = c("#E41A1C", "#377EB8"))

ggsave(file.path(output_dir, "es_terms_boxplot.png"),
       p4, width = 10, height = 8, dpi = 300)

# PLOT 5: Term Relationship Scatter
cat("  Plot 5: Term1 vs Term2 relationship\n")

p5 <- ggplot(all_results, aes(x = term2, y = term1, color = scenario)) +
  geom_point(alpha = 0.4, size = 1.5) +
  stat_ellipse(aes(group = scenario), level = 0.95, linewidth = 1) +
  labs(
    title = "Relationship Between Energy Score Terms",
    subtitle = "95% confidence ellipses by scenario. ES = Term 1 - Term 2",
    x = "Term 2: 0.5 × Mean Pairwise Distance (diversity credit)",
    y = "Term 1: Mean Distance to Observation (sharpness penalty)",
    color = "Scenario"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  ) +
  scale_color_brewer(palette = "Set2") +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 3)))

ggsave(file.path(output_dir, "es_terms_scatter.png"),
       p5, width = 10, height = 7, dpi = 300)

# PLOT 6: ES vs Terms
cat("  Plot 6: ES vs individual terms\n")

es_vs_terms <- all_results %>%
  select(scenario, es_custom, term1, term2) %>%
  pivot_longer(cols = c(term1, term2),
               names_to = "term",
               values_to = "term_value")

p6 <- ggplot(es_vs_terms, aes(x = term_value, y = es_custom, color = scenario)) +
  geom_point(alpha = 0.3, size = 1) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1) +
  facet_wrap(~ term, scales = "free_x", ncol = 2,
             labeller = labeller(term = c(
               term1 = "Term 1: Distance to Observation",
               term2 = "Term 2: 0.5 × Pairwise Distance"
             ))) +
  labs(
    title = "Energy Score vs Individual Terms",
    subtitle = "ES = Term 1 - Term 2. Term 1 positively contributes, Term 2 negatively contributes.",
    x = "Term Value",
    y = "Energy Score",
    color = "Scenario"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom",
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(face = "bold")
  ) +
  scale_color_brewer(palette = "Set2") +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 3)))

ggsave(file.path(output_dir, "es_vs_terms_scatter.png"),
       p6, width = 12, height = 6, dpi = 300)

cat("  All plots saved!\n\n")

# ============================================================================
# SAVE RESULTS
# ============================================================================

cat("Saving results to results/test_results/...\n")

# Save full results
saveRDS(all_results, file.path(output_dir, "es_decomposition_results.rds"))
cat("  Full results:", file.path(output_dir, "es_decomposition_results.rds"), "\n")

# Save summary statistics
write.csv(decomposition_summary,
          file.path(output_dir, "es_decomposition_summary.csv"),
          row.names = FALSE)
cat("  Summary statistics:", file.path(output_dir, "es_decomposition_summary.csv"), "\n")

# Save verification results
write.csv(verification_summary,
          file.path(output_dir, "es_verification_summary.csv"),
          row.names = FALSE)
cat("  Verification summary:", file.path(output_dir, "es_verification_summary.csv"), "\n\n")

# ============================================================================
# FINAL SUMMARY
# ============================================================================

cat("=========================================\n")
cat("SUMMARY\n")
cat("=========================================\n\n")

cat("Simulations completed:\n")
cat(sprintf("  Total simulations: %d (6 scenarios × %d iterations)\n",
            nrow(all_results), n_sims))
cat(sprintf("  Forecast samples per iteration: %d\n", n_samples))
cat("\n")

cat("Key findings:\n")
cat(sprintf("  1. ES Verification: Max difference = %.2e (tolerance: %.2e)\n",
            max_diff, tolerance))

# Find best and worst scenarios
best_scenario <- decomposition_summary$scenario[which.min(decomposition_summary$mean_es)]
worst_scenario <- decomposition_summary$scenario[which.max(decomposition_summary$mean_es)]

cat(sprintf("  2. Best performing forecast: %s (mean ES = %.4f)\n",
            best_scenario, min(decomposition_summary$mean_es)))
cat(sprintf("  3. Worst performing forecast: %s (mean ES = %.4f)\n",
            worst_scenario, max(decomposition_summary$mean_es)))

cat("\nFiles created in results/test_results/:\n")
cat("  Plots (6):\n")
cat("    - es_verification_plot.png\n")
cat("    - es_by_scenario_boxplot.png\n")
cat("    - es_decomposition_barplot.png\n")
cat("    - es_terms_boxplot.png\n")
cat("    - es_terms_scatter.png\n")
cat("    - es_vs_terms_scatter.png\n")
cat("  Data (3):\n")
cat("    - es_decomposition_results.rds\n")
cat("    - es_decomposition_summary.csv\n")
cat("    - es_verification_summary.csv\n\n")

cat("Done!\n")
