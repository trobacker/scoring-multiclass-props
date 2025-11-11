#!/usr/bin/env Rscript
# Parallel execution of higher-dimensional simulation scenarios (K > 3 classes)
# This script extends the 3-class simulations to arbitrary dimensions
# Usage: Rscript run_simulations_highdim.R
# Note: Run from ./src directory (paths assume execution from src/)

library(ggplot2)
library(tidyr)
library(dplyr)
library(scoringRules)
library(brms)
library(parallel)

# Source helper functions
source("./src/helper-functions.R")

cat("=========================================\n")
cat("HIGH-DIMENSIONAL SIMULATION RUNNER\n")
cat("=========================================\n\n")

# Detect cores
n_cores_total <- detectCores()
n_cores_inner <- max(1, n_cores_total - 2)  # Cores for nsim parallelization
n_cores_outer <- min(4, n_cores_total)      # Cores for scenario parallelization

cat("System has", n_cores_total, "cores\n")
cat("Using", n_cores_outer, "cores for scenario-level parallelization\n")
cat("Using", n_cores_inner, "cores per scenario for iteration parallelization\n\n")

# es_sampling function (dimension-agnostic, parallelized)
es_sampling <- function(alpha,
                        thetas = NULL,
                        thetas_wrong = NULL,
                        n_counts = 5,
                        nsim = 1000,
                        seed = 42,
                        parallel = TRUE,
                        n_cores = parallel::detectCores() - 2){

  K <- length(alpha)  # Number of classes (dimension)
  p_matrix <- matrix(unlist(thetas), nrow = K, byrow = FALSE)
  p_matrix_wrong <- matrix(unlist(thetas_wrong), nrow = K, byrow = FALSE)

  # True mean of Dirichlet(alpha) is alpha/sum(alpha)
  true_mean <- alpha / sum(alpha)

  if(parallel && n_cores > 1) {
    results <- parallel::mclapply(1:nsim, function(i) {
      set.seed(seed + i)
      prob_vec <- as.numeric(brms::rdirichlet(1, alpha))
      C <- rmultinom(n = 1, size = n_counts, prob = prob_vec)
      p <- C / n_counts
      es_correct <- scoringRules::es_sample(y = as.numeric(p), dat = p_matrix)
      es_incorrect <- scoringRules::es_sample(y = as.numeric(p), dat = p_matrix_wrong)
      # Score forecasts against true mean as observation
      es_correct_vs_true_mean <- scoringRules::es_sample(y = true_mean, dat = p_matrix)
      es_incorrect_vs_true_mean <- scoringRules::es_sample(y = true_mean, dat = p_matrix_wrong)
      return(c(es_correct = es_correct, es_incorrect = es_incorrect,
               es_correct_vs_true_mean = es_correct_vs_true_mean,
               es_incorrect_vs_true_mean = es_incorrect_vs_true_mean))
    }, mc.cores = n_cores, mc.set.seed = FALSE)

    results_matrix <- do.call(rbind, results)
    es_correct <- results_matrix[, 1]
    es_incorrect <- results_matrix[, 2]
    es_correct_vs_true_mean <- results_matrix[, 3]
    es_incorrect_vs_true_mean <- results_matrix[, 4]
  } else {
    set.seed(seed)
    es_correct <- rep(NA, nsim)
    es_incorrect <- rep(NA, nsim)
    es_correct_vs_true_mean <- rep(NA, nsim)
    es_incorrect_vs_true_mean <- rep(NA, nsim)
    for(i in 1:nsim){
      prob_vec <- as.numeric(brms::rdirichlet(1, alpha))
      C <- rmultinom(n = 1, size = n_counts, prob = prob_vec)
      p <- C / n_counts
      es_correct[i] <- scoringRules::es_sample(y = as.numeric(p), dat = p_matrix)
      es_incorrect[i] <- scoringRules::es_sample(y = as.numeric(p), dat = p_matrix_wrong)
      # Score forecasts against true mean as observation
      es_correct_vs_true_mean[i] <- scoringRules::es_sample(y = true_mean, dat = p_matrix)
      es_incorrect_vs_true_mean[i] <- scoringRules::es_sample(y = true_mean, dat = p_matrix_wrong)
    }
  }
  return(list(es_correct = es_correct, es_incorrect = es_incorrect,
              es_correct_vs_true_mean = es_correct_vs_true_mean,
              es_incorrect_vs_true_mean = es_incorrect_vs_true_mean))
}

# es_mn_sampling function (dimension-agnostic, parallelized)
es_mn_sampling <- function(alpha,
                           thetas = NULL,
                           thetas_wrong = NULL,
                           n_counts = 5,
                           nsim = 1000,
                           N_multinomial = 100,
                           seed = 42,
                           parallel = TRUE,
                           n_cores = parallel::detectCores() - 2){

  # True mean of Dirichlet(alpha) is alpha/sum(alpha)
  K <- length(alpha)
  true_mean <- alpha / sum(alpha)

  if(parallel && n_cores > 1) {
    results <- parallel::mclapply(1:nsim, function(i) {
      set.seed(seed + i)
      prob_vec <- as.numeric(brms::rdirichlet(1, alpha))
      C <- rmultinom(n = 1, size = n_counts, prob = prob_vec)
      p <- C / n_counts

      samp_multinomial_counts <- do.call(cbind,
        lapply(thetas, function(theta) rmultinom(n = N_multinomial, size = n_counts, prob = theta))
      )
      samp_multinomial_counts_wrong <- do.call(cbind,
        lapply(thetas_wrong, function(theta) rmultinom(n = N_multinomial, size = n_counts, prob = theta))
      )

      p_matrix <- samp_multinomial_counts / n_counts
      p_matrix_wrong <- samp_multinomial_counts_wrong / n_counts

      es_correct <- scoringRules::es_sample(y = as.numeric(p), dat = p_matrix)
      es_incorrect <- scoringRules::es_sample(y = as.numeric(p), dat = p_matrix_wrong)
      # Score forecasts against true mean as observation
      es_correct_vs_true_mean <- scoringRules::es_sample(y = true_mean, dat = p_matrix)
      es_incorrect_vs_true_mean <- scoringRules::es_sample(y = true_mean, dat = p_matrix_wrong)
      return(c(es_correct = es_correct, es_incorrect = es_incorrect,
               es_correct_vs_true_mean = es_correct_vs_true_mean,
               es_incorrect_vs_true_mean = es_incorrect_vs_true_mean))
    }, mc.cores = n_cores, mc.set.seed = FALSE)

    results_matrix <- do.call(rbind, results)
    es_correct <- results_matrix[, 1]
    es_incorrect <- results_matrix[, 2]
    es_correct_vs_true_mean <- results_matrix[, 3]
    es_incorrect_vs_true_mean <- results_matrix[, 4]
  } else {
    set.seed(seed)
    es_correct <- rep(NA, nsim)
    es_incorrect <- rep(NA, nsim)
    es_correct_vs_true_mean <- rep(NA, nsim)
    es_incorrect_vs_true_mean <- rep(NA, nsim)
    for(i in 1:nsim){
      prob_vec <- as.numeric(brms::rdirichlet(1, alpha))
      C <- rmultinom(n = 1, size = n_counts, prob = prob_vec)
      p <- C / n_counts
      samp_multinomial_counts <- do.call(cbind,
        lapply(thetas, function(theta) rmultinom(n = N_multinomial, size = n_counts, prob = theta))
      )
      samp_multinomial_counts_wrong <- do.call(cbind,
        lapply(thetas_wrong, function(theta) rmultinom(n = N_multinomial, size = n_counts, prob = theta))
      )
      p_matrix <- samp_multinomial_counts / n_counts
      p_matrix_wrong <- samp_multinomial_counts_wrong / n_counts

      es_correct[i] <- scoringRules::es_sample(y = as.numeric(p), dat = p_matrix)
      es_incorrect[i] <- scoringRules::es_sample(y = as.numeric(p), dat = p_matrix_wrong)
      # Score forecasts against true mean as observation
      es_correct_vs_true_mean[i] <- scoringRules::es_sample(y = true_mean, dat = p_matrix)
      es_incorrect_vs_true_mean[i] <- scoringRules::es_sample(y = true_mean, dat = p_matrix_wrong)
    }
  }
  return(list(es_correct = es_correct, es_incorrect = es_incorrect,
              es_correct_vs_true_mean = es_correct_vs_true_mean,
              es_incorrect_vs_true_mean = es_incorrect_vs_true_mean))
}

# Helper functions
result_to_df <- function(result, n_counts, method = "mn") {
  data.frame(
    es_correct = result$es_correct,
    es_incorrect = result$es_incorrect,
    es_correct_vs_true_mean = result$es_correct_vs_true_mean,
    es_incorrect_vs_true_mean = result$es_incorrect_vs_true_mean,
    es_diff = result$es_correct - result$es_incorrect,
    es_diff_true_mean_obs = result$es_correct_vs_true_mean - result$es_incorrect_vs_true_mean,
    n_counts = n_counts,
    method = method,
    sim_label = paste0(method, "_n", n_counts)
  )
}

run_simulation_experiment <- function(alpha, thetas, thetas_wrong,
                                      n_counts_vec = c(1, 2, 3, 4, 5, 10, 25, 50, 100, 500),
                                      nsim = 1000, base_seed = 42,
                                      parallel = TRUE, n_cores = parallel::detectCores() - 2) {
  results_list <- list()
  for (i in seq_along(n_counts_vec)) {
    n <- n_counts_vec[i]
    cat("  n_counts =", n, "...")
    result_nonmn <- es_sampling(thetas = thetas, thetas_wrong = thetas_wrong, alpha = alpha,
                                seed = base_seed + i, nsim = nsim, n_counts = n,
                                parallel = parallel, n_cores = n_cores)
    result_mn <- es_mn_sampling(thetas = thetas, thetas_wrong = thetas_wrong, alpha = alpha,
                               seed = base_seed + 100 + i, nsim = nsim, n_counts = n,
                               parallel = parallel, n_cores = n_cores)
    results_list[[paste0("nonmn_", i)]] <- result_to_df(result_nonmn, n, "non-mn")
    results_list[[paste0("mn_", i)]] <- result_to_df(result_mn, n, "mn")
    cat(" done\n")
  }
  df_combined <- do.call(rbind, results_list)
  rownames(df_combined) <- NULL
  return(df_combined)
}

save_results <- function(df, scenario_name, alpha, results_dir = "./results/data") {
  alpha_str <- paste(alpha, collapse = "_")
  alpha_dir <- file.path(results_dir, paste0("alpha_", alpha_str))
  if (!dir.exists(alpha_dir)) {
    dir.create(alpha_dir, recursive = TRUE)
  }
  filename <- file.path(alpha_dir, paste0(scenario_name, "_alpha_", alpha_str, ".rds"))
  saveRDS(df, filename)
  cat("Results saved to:", filename, "\n")
  return(filename)
}

# Save plot to disk (organized by alpha subdirectory)
save_plot <- function(plot_obj, plot_name, scenario_name, alpha,
                     results_dir = "./results/plots",
                     width = 10, height = 6) {
  alpha_str <- paste(alpha, collapse = "_")
  alpha_dir <- file.path(results_dir, paste0("alpha_", alpha_str))

  if (!dir.exists(alpha_dir)) {
    dir.create(alpha_dir, recursive = TRUE)
  }

  filename <- file.path(alpha_dir, paste0(scenario_name, "_", plot_name, "_alpha_", alpha_str, ".png"))
  ggsave(filename, plot = plot_obj, width = width, height = height, dpi = 300)
  cat("Plot saved to:", filename, "\n")
  return(filename)
}

# Generate all plots for a simulation result
generate_all_plots <- function(df_results) {
  plots <- list()

  # Plot 1: Faceted boxplots comparing correct vs incorrect (random observations)
  # COMMENTED OUT - Dotplots are preferred for visualization
  # df_long <- df_results %>%
  #   pivot_longer(
  #     cols = c(es_correct, es_incorrect),
  #     names_to = "forecast_type",
  #     values_to = "energy_score"
  #   )
  #
  # plots$comparison <- ggplot(df_long, aes(x = factor(n_counts), y = energy_score, fill = forecast_type)) +
  #   geom_boxplot() +
  #   facet_wrap(~ method, labeller = labeller(method = c("mn" = "Multinomial Sampling",
  #                                                         "non-mn" = "No Multinomial Sampling"))) +
  #   labs(
  #     title = "Energy Scores Comparison by Method and Sample Size",
  #     subtitle = "Observations drawn from Dir(alpha)",
  #     x = "n_counts (number of sequences)",
  #     y = "Energy Score",
  #     fill = "Forecast Type"
  #   ) +
  #   scale_fill_manual(
  #     values = c(es_correct = rgb(0, 0, 1, 0.4),
  #                es_incorrect = rgb(1, 0, 0, 0.5)),
  #     labels = c("Correct (samples from true alpha)",
  #                "Incorrect (samples from wrong alpha)")
  #   ) +
  #   theme_bw() +
  #   theme(legend.position = "bottom")

  # Plot 1b: Scores when observation is true mean
  # COMMENTED OUT - Dotplots are preferred for visualization
  df_long_true_mean <- df_results %>%
    pivot_longer(
      cols = c(es_correct_vs_true_mean, es_incorrect_vs_true_mean),
      names_to = "forecast_type",
      values_to = "energy_score"
    )

  # plots$comparison_true_mean_obs <- ggplot(df_long_true_mean, aes(x = factor(n_counts), y = energy_score, fill = forecast_type)) +
  #   geom_boxplot() +
  #   facet_wrap(~ method, labeller = labeller(method = c("mn" = "Multinomial Sampling",
  #                                                         "non-mn" = "No Multinomial Sampling"))) +
  #   labs(
  #     title = "Energy Scores When Observation = True Mean",
  #     subtitle = "Observation is alpha/sum(alpha) (deterministic)",
  #     x = "n_counts (number of sequences)",
  #     y = "Energy Score",
  #     fill = "Forecast Type"
  #   ) +
  #   scale_fill_manual(
  #     values = c(es_correct_vs_true_mean = rgb(0, 0, 1, 0.4),
  #                es_incorrect_vs_true_mean = rgb(1, 0, 0, 0.5)),
  #     labels = c("Correct forecast vs true mean",
  #                "Incorrect forecast vs true mean")
  #   ) +
  #   theme_bw() +
  #   theme(legend.position = "bottom")

  # Plot 2: Side-by-side comparison by method
  # COMMENTED OUT - Dotplots are preferred for visualization
  # plots$differences_method <- ggplot(df_results, aes(x = factor(n_counts), y = es_diff, fill = method)) +
  #   geom_boxplot() +
  #   geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.5) +
  #   labs(
  #     title = "Energy Score Differences by Method and Sample Size",
  #     subtitle = "Comparing multinomial vs non-multinomial sampling",
  #     x = "n_counts (number of sequences)",
  #     y = "ES Difference (Correct - Incorrect)",
  #     fill = "Method"
  #   ) +
  #   scale_fill_manual(
  #     values = c("mn" = "steelblue", "non-mn" = "coral"),
  #     labels = c("Multinomial Sampling", "No Multinomial Sampling")
  #   ) +
  #   theme_bw() +
  #   theme(legend.position = "bottom")

  # Prepare data for dotplots
  df_long <- df_results %>%
    pivot_longer(
      cols = c(es_correct, es_incorrect),
      names_to = "forecast_type",
      values_to = "energy_score"
    )

  # Plot 3: Dotplot version of comparison (shows discreteness better)
  plots$comparison_dotplot <- ggplot(df_long, aes(x = factor(n_counts), y = energy_score, color = forecast_type)) +
    geom_point(alpha = 0.3, size = 0.8, position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.6)) +
    facet_wrap(~ method, labeller = labeller(method = c("mn" = "Multinomial Sampling",
                                                          "non-mn" = "No Multinomial Sampling"))) +
    labs(
      title = "Energy Scores Comparison by Method and Sample Size (Dotplot)",
      subtitle = "Observations drawn from Dir(alpha) - reveals discreteness at small n_counts",
      x = "n_counts (number of sequences)",
      y = "Energy Score",
      color = "Forecast Type"
    ) +
    scale_color_manual(
      values = c(es_correct = "blue", es_incorrect = "red"),
      labels = c("Correct (samples from true alpha)",
                 "Incorrect (samples from wrong alpha)")
    ) +
    theme_bw() +
    theme(legend.position = "bottom")

  # Plot 3b: Dotplot when observation is true mean
  plots$comparison_true_mean_obs_dotplot <- ggplot(df_long_true_mean, aes(x = factor(n_counts), y = energy_score, color = forecast_type)) +
    geom_point(alpha = 0.3, size = 0.8, position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.6)) +
    facet_wrap(~ method, labeller = labeller(method = c("mn" = "Multinomial Sampling",
                                                          "non-mn" = "No Multinomial Sampling"))) +
    labs(
      title = "Energy Scores When Observation = True Mean (Dotplot)",
      subtitle = "Observation is alpha/sum(alpha) (deterministic)",
      x = "n_counts (number of sequences)",
      y = "Energy Score",
      color = "Forecast Type"
    ) +
    scale_color_manual(
      values = c(es_correct_vs_true_mean = "blue", es_incorrect_vs_true_mean = "red"),
      labels = c("Correct forecast vs true mean",
                 "Incorrect forecast vs true mean")
    ) +
    theme_bw() +
    theme(legend.position = "bottom")

  # Plot 4: Dotplot version of differences by method (shows discreteness better)
  plots$differences_method_dotplot <- ggplot(df_results, aes(x = factor(n_counts), y = es_diff, color = method)) +
    geom_point(alpha = 0.3, size = 0.8, position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.6)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
    labs(
      title = "Energy Score Differences by Method and Sample Size (Dotplot)",
      subtitle = "Each dot represents one simulation run - reveals discreteness at small n_counts",
      x = "n_counts (number of sequences)",
      y = "ES Difference (Correct - Incorrect)",
      color = "Method"
    ) +
    scale_color_manual(
      values = c("mn" = "steelblue", "non-mn" = "coral"),
      labels = c("Multinomial Sampling", "No Multinomial Sampling")
    ) +
    theme_bw() +
    theme(legend.position = "bottom")

  # Plot 5: Dotplot comparing forecast performance when obs = true mean
  plots$diff_true_mean_obs_dotplot <- ggplot(df_results, aes(x = factor(n_counts), y = es_diff_true_mean_obs, color = method)) +
    geom_point(alpha = 0.3, size = 0.8, position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.6)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
    labs(
      title = "Forecast Comparison When Observation = True Mean (Dotplot)",
      subtitle = "Negative values indicate wrong forecast performs better (surprising!)",
      x = "n_counts (number of sequences)",
      y = "ES Difference (Correct - Incorrect) at True Mean Observation",
      color = "Method"
    ) +
    scale_color_manual(
      values = c("mn" = "steelblue", "non-mn" = "coral"),
      labels = c("Multinomial Sampling", "No Multinomial Sampling")
    ) +
    theme_bw() +
    theme(legend.position = "bottom")

  return(plots)
}

# Define higher-dimensional scenarios (K = 5 classes)
scenarios <- list(
  list(
    name = "scenario1_5class_ascending",
    alpha = c(2, 5, 10, 20, 40),
    alpha_wrong = c(40, 20, 10, 5, 2),  # Reversed proportions
    base_seed = 42
  ),
  list(
    name = "scenario2_5class_narrow_calibrated",
    alpha = c(2, 5, 10, 20, 40),
    alpha_wrong = c(2, 5, 10, 20, 40) * 10,  # Same proportions, less dispersed
    base_seed = 52
  ),
  list(
    name = "scenario3_5class_dispersed_calibrated",
    alpha = c(2, 5, 10, 20, 40),
    alpha_wrong = c(2, 5, 10, 20, 40) / 5,  # Same proportions, more dispersed
    base_seed = 62
  ),
  list(
    name = "scenario4_5class_uniform_vs_peaked",
    alpha = c(10, 10, 10, 10, 10),  # Uniform
    alpha_wrong = c(2, 2, 50, 2, 2),  # Peaked at middle class
    base_seed = 72
  ),
  list(
    name = "scenario5_5class_moderate_misspec",
    alpha = c(2, 5, 10, 20, 40),
    alpha_wrong = c(3, 7, 12, 18, 35),  # Moderately different
    base_seed = 82
  )
)

# Function to run a single scenario
run_scenario <- function(scenario_config, n_samp = 100,
                        n_counts_vec = c(1, 2, 3, 4, 5, 10, 25, 50, 100, 500),
                        nsim = 1000, n_cores_inner = 12) {

  cat("\n========================================\n")
  cat("Running:", scenario_config$name, "\n")
  cat("Alpha (true):", paste(scenario_config$alpha, collapse = ", "), "\n")
  cat("Alpha (wrong):", paste(scenario_config$alpha_wrong, collapse = ", "), "\n")
  cat("K (number of classes):", length(scenario_config$alpha), "\n")
  cat("========================================\n")

  # Pre-allocate distributions
  thetas <- lapply(FUN = draw_one_dirichlet, X = rep(1, n_samp), alpha = scenario_config$alpha)
  thetas_wrong <- lapply(FUN = draw_one_dirichlet, X = rep(1, n_samp), alpha = scenario_config$alpha_wrong)

  # Run simulations
  df_results <- run_simulation_experiment(
    alpha = scenario_config$alpha,
    thetas = thetas,
    thetas_wrong = thetas_wrong,
    n_counts_vec = n_counts_vec,
    nsim = nsim,
    base_seed = scenario_config$base_seed,
    parallel = TRUE,
    n_cores = n_cores_inner
  )

  # Add metadata
  df_results$scenario <- scenario_config$name
  df_results$alpha_true <- paste(scenario_config$alpha, collapse = "_")
  df_results$alpha_wrong <- paste(scenario_config$alpha_wrong, collapse = "_")
  df_results$K <- length(scenario_config$alpha)

  # Generate plots
  plots <- generate_all_plots(df_results)

  # Save results
  save_results(df_results, scenario_config$name, scenario_config$alpha)

  # Save plots
  # COMMENTED OUT - Boxplots not needed, dotplots are preferred
  # save_plot(plots$comparison, "comparison", scenario_config$name, scenario_config$alpha)
  # save_plot(plots$comparison_true_mean_obs, "comparison_true_mean_obs", scenario_config$name, scenario_config$alpha)
  # save_plot(plots$differences_method, "differences_method", scenario_config$name, scenario_config$alpha)
  save_plot(plots$comparison_dotplot, "comparison_dotplot", scenario_config$name, scenario_config$alpha)
  save_plot(plots$comparison_true_mean_obs_dotplot, "comparison_true_mean_obs_dotplot", scenario_config$name, scenario_config$alpha)
  save_plot(plots$differences_method_dotplot, "differences_method_dotplot", scenario_config$name, scenario_config$alpha)
  save_plot(plots$diff_true_mean_obs_dotplot, "diff_true_mean_obs_dotplot", scenario_config$name, scenario_config$alpha)

  # Note: No ternary plots for K > 3

  cat("Completed:", scenario_config$name, "\n")
  return(df_results)
}

# Record start time
start_time <- Sys.time()

# Detect if running from RStudio (fork safety issues on macOS)
in_rstudio <- Sys.getenv("RSTUDIO") == "1"

if (in_rstudio) {
  cat("\nDetected RStudio environment - running scenarios SEQUENTIALLY\n")
  cat("(Using parallelization within each scenario for ~10x speedup)\n")
  cat("For maximum speed, run from terminal: Rscript src/run_simulations_highdim.R\n\n")

  # Run scenarios sequentially (but each uses parallel processing internally)
  all_results <- lapply(scenarios, function(scenario_config) {
    run_scenario(scenario_config,
                 n_samp = 100,
                 nsim = 1000,
                 n_cores_inner = n_cores_inner)
  })
} else {
  # Run all scenarios in parallel (scenario-level parallelization)
  cat("\nRunning", length(scenarios), "scenarios in parallel...\n")

  all_results <- parallel::mclapply(scenarios, function(scenario_config) {
    run_scenario(scenario_config,
                 n_samp = 100,
                 nsim = 1000,
                 n_cores_inner = n_cores_inner)
  }, mc.cores = n_cores_outer)
}

# Record end time
end_time <- Sys.time()
elapsed_time <- difftime(end_time, start_time, units = "mins")

cat("\n=========================================\n")
cat("ALL SIMULATIONS COMPLETE!\n")
cat("Total time:", round(elapsed_time, 2), "minutes\n")
cat("=========================================\n")

# Combine all results
combined_results <- do.call(rbind, all_results)

# Save combined results
combined_file <- "./results/all_scenarios_highdim_combined.rds"
saveRDS(combined_results, combined_file)
cat("Combined results saved to:", combined_file, "\n")

cat("\nDone! Check the results/ directory for output files.\n")
