#!/usr/bin/env Rscript
# ADDITIONAL extreme scenarios for SYMMETRIC alpha
# True alpha = c(10, 10, 10) - symmetric/uniform distribution
# Two new scenarios: ultra-dispersed and ultra-narrow
# Usage: Rscript run_simulations_alpha_symmetric_additional.R

library(ggplot2)
library(tidyr)
library(dplyr)
library(scoringRules)
library(brms)
library(Ternary)
library(parallel)

# Source helper functions - robust path resolution for different execution contexts
if (file.exists("./src/helper-functions.R")) {
  # Running from project root
  source("./src/helper-functions.R")
} else if (file.exists("helper-functions.R")) {
  # Running from src/ directory
  source("helper-functions.R")
} else {
  stop("Cannot find helper-functions.R. Please run from project root or src/ directory.")
}

cat("=========================================\n")
cat("PARALLEL SIMULATION RUNNER\n")
cat("=========================================\n\n")

# Detect cores
n_cores_total <- detectCores()
n_cores_inner <- max(1, n_cores_total - 2)  # Cores for nsim parallelization
n_cores_outer <- min(4, n_cores_total)      # Cores for scenario parallelization

cat("System has", n_cores_total, "cores\n")
cat("Using", n_cores_outer, "cores for scenario-level parallelization\n")
cat("Using", n_cores_inner, "cores per scenario for iteration parallelization\n\n")

# Load all functions from simulation_1.qmd
# These are copied here for standalone execution

# es_sampling function (parallelized)
# Now generates fresh forecast samples for each simulation iteration
es_sampling <- function(alpha = c(2,4,4),
                           alpha_wrong = NULL,
                           n_counts = 5,
                           n_samp = 100,
                           nsim = 1000,
                           seed = 42,
                           parallel = TRUE,
                           n_cores = parallel::detectCores() - 2){

  K <- length(alpha)

  # True mean of Dirichlet(alpha) is alpha/sum(alpha)
  true_mean <- alpha / sum(alpha)

  if(parallel && n_cores > 1) {
    results <- parallel::mclapply(1:nsim, function(i) {
      set.seed(seed + i)

      # Generate fresh forecast samples for this iteration
      thetas <- lapply(FUN = draw_one_dirichlet, X = rep(1, n_samp), alpha = alpha)
      thetas_wrong <- lapply(FUN = draw_one_dirichlet, X = rep(1, n_samp), alpha = alpha_wrong)
      p_matrix <- matrix(unlist(thetas), nrow = K, byrow = FALSE)
      p_matrix_wrong <- matrix(unlist(thetas_wrong), nrow = K, byrow = FALSE)

      # Generate observation
      prob_vec <- as.numeric(brms::rdirichlet(1, alpha))
      C <- rmultinom(n = 1, size = n_counts, prob = prob_vec)
      p <- C / n_counts

      # Score forecasts using manual energy score (returns ES, term1, term2)
      es_result_correct <- manual_energy_score(y = as.numeric(p), dat = p_matrix)
      es_result_incorrect <- manual_energy_score(y = as.numeric(p), dat = p_matrix_wrong)
      # Score forecasts against true mean as observation
      es_result_correct_vs_true_mean <- manual_energy_score(y = true_mean, dat = p_matrix)
      es_result_incorrect_vs_true_mean <- manual_energy_score(y = true_mean, dat = p_matrix_wrong)

      return(c(
        es_correct = es_result_correct$ES,
        es_incorrect = es_result_incorrect$ES,
        es_correct_vs_true_mean = es_result_correct_vs_true_mean$ES,
        es_incorrect_vs_true_mean = es_result_incorrect_vs_true_mean$ES,
        term1_correct = es_result_correct$term1,
        term1_incorrect = es_result_incorrect$term1,
        term1_correct_vs_true_mean = es_result_correct_vs_true_mean$term1,
        term1_incorrect_vs_true_mean = es_result_incorrect_vs_true_mean$term1,
        term2_correct = es_result_correct$term2,
        term2_incorrect = es_result_incorrect$term2,
        term2_correct_vs_true_mean = es_result_correct_vs_true_mean$term2,
        term2_incorrect_vs_true_mean = es_result_incorrect_vs_true_mean$term2
      ))
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
    term1_correct <- rep(NA, nsim)
    term1_incorrect <- rep(NA, nsim)
    term1_correct_vs_true_mean <- rep(NA, nsim)
    term1_incorrect_vs_true_mean <- rep(NA, nsim)
    term2_correct <- rep(NA, nsim)
    term2_incorrect <- rep(NA, nsim)
    term2_correct_vs_true_mean <- rep(NA, nsim)
    term2_incorrect_vs_true_mean <- rep(NA, nsim)

    for(i in 1:nsim){
      # Generate fresh forecast samples for this iteration
      thetas <- lapply(FUN = draw_one_dirichlet, X = rep(1, n_samp), alpha = alpha)
      thetas_wrong <- lapply(FUN = draw_one_dirichlet, X = rep(1, n_samp), alpha = alpha_wrong)
      p_matrix <- matrix(unlist(thetas), nrow = K, byrow = FALSE)
      p_matrix_wrong <- matrix(unlist(thetas_wrong), nrow = K, byrow = FALSE)

      # Generate observation
      prob_vec <- as.numeric(brms::rdirichlet(1, alpha))
      C <- rmultinom(n = 1, size = n_counts, prob = prob_vec)
      p <- C / n_counts

      # Score forecasts using manual energy score (returns ES, term1, term2)
      es_result_correct <- manual_energy_score(y = as.numeric(p), dat = p_matrix)
      es_result_incorrect <- manual_energy_score(y = as.numeric(p), dat = p_matrix_wrong)
      # Score forecasts against true mean as observation
      es_result_correct_vs_true_mean <- manual_energy_score(y = true_mean, dat = p_matrix)
      es_result_incorrect_vs_true_mean <- manual_energy_score(y = true_mean, dat = p_matrix_wrong)

      # Extract ES values
      es_correct[i] <- es_result_correct$ES
      es_incorrect[i] <- es_result_incorrect$ES
      es_correct_vs_true_mean[i] <- es_result_correct_vs_true_mean$ES
      es_incorrect_vs_true_mean[i] <- es_result_incorrect_vs_true_mean$ES

      # Extract term1 values
      term1_correct[i] <- es_result_correct$term1
      term1_incorrect[i] <- es_result_incorrect$term1
      term1_correct_vs_true_mean[i] <- es_result_correct_vs_true_mean$term1
      term1_incorrect_vs_true_mean[i] <- es_result_incorrect_vs_true_mean$term1

      # Extract term2 values
      term2_correct[i] <- es_result_correct$term2
      term2_incorrect[i] <- es_result_incorrect$term2
      term2_correct_vs_true_mean[i] <- es_result_correct_vs_true_mean$term2
      term2_incorrect_vs_true_mean[i] <- es_result_incorrect_vs_true_mean$term2
    }
  }
  return(list(
    es_correct = es_correct, es_incorrect = es_incorrect,
    es_correct_vs_true_mean = es_correct_vs_true_mean,
    es_incorrect_vs_true_mean = es_incorrect_vs_true_mean,
    term1_correct = term1_correct, term1_incorrect = term1_incorrect,
    term1_correct_vs_true_mean = term1_correct_vs_true_mean,
    term1_incorrect_vs_true_mean = term1_incorrect_vs_true_mean,
    term2_correct = term2_correct, term2_incorrect = term2_incorrect,
    term2_correct_vs_true_mean = term2_correct_vs_true_mean,
    term2_incorrect_vs_true_mean = term2_incorrect_vs_true_mean
  ))
}

# es_mn_sampling function (parallelized)
# Now generates fresh forecast samples for each simulation iteration
# IMPORTANT: Uses same seed offset as es_sampling to ensure same base thetas
es_mn_sampling <- function(alpha = c(2,4,4),
                           alpha_wrong = NULL,
                           n_counts = 5,
                           n_samp = 100,
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

      # Generate fresh forecast samples for this iteration (same as es_sampling with same seed)
      thetas <- lapply(FUN = draw_one_dirichlet, X = rep(1, n_samp), alpha = alpha)
      thetas_wrong <- lapply(FUN = draw_one_dirichlet, X = rep(1, n_samp), alpha = alpha_wrong)

      # Generate observation
      prob_vec <- as.numeric(brms::rdirichlet(1, alpha))
      C <- rmultinom(n = 1, size = n_counts, prob = prob_vec)
      p <- C / n_counts

      # Apply multinomial sampling to forecasts
      samp_multinomial_counts <- do.call(cbind,
        lapply(thetas, function(theta) rmultinom(n = N_multinomial, size = n_counts, prob = theta))
      )
      samp_multinomial_counts_wrong <- do.call(cbind,
        lapply(thetas_wrong, function(theta) rmultinom(n = N_multinomial, size = n_counts, prob = theta))
      )

      p_matrix <- samp_multinomial_counts / n_counts
      p_matrix_wrong <- samp_multinomial_counts_wrong / n_counts

      # Score forecasts using manual energy score (returns ES, term1, term2)
      es_result_correct <- manual_energy_score(y = as.numeric(p), dat = p_matrix)
      es_result_incorrect <- manual_energy_score(y = as.numeric(p), dat = p_matrix_wrong)
      # Score forecasts against true mean as observation
      es_result_correct_vs_true_mean <- manual_energy_score(y = true_mean, dat = p_matrix)
      es_result_incorrect_vs_true_mean <- manual_energy_score(y = true_mean, dat = p_matrix_wrong)

      return(c(
        es_correct = es_result_correct$ES,
        es_incorrect = es_result_incorrect$ES,
        es_correct_vs_true_mean = es_result_correct_vs_true_mean$ES,
        es_incorrect_vs_true_mean = es_result_incorrect_vs_true_mean$ES,
        term1_correct = es_result_correct$term1,
        term1_incorrect = es_result_incorrect$term1,
        term1_correct_vs_true_mean = es_result_correct_vs_true_mean$term1,
        term1_incorrect_vs_true_mean = es_result_incorrect_vs_true_mean$term1,
        term2_correct = es_result_correct$term2,
        term2_incorrect = es_result_incorrect$term2,
        term2_correct_vs_true_mean = es_result_correct_vs_true_mean$term2,
        term2_incorrect_vs_true_mean = es_result_incorrect_vs_true_mean$term2
      ))
    }, mc.cores = n_cores, mc.set.seed = FALSE)

    results_matrix <- do.call(rbind, results)
    es_correct <- results_matrix[, 1]
    es_incorrect <- results_matrix[, 2]
    es_correct_vs_true_mean <- results_matrix[, 3]
    es_incorrect_vs_true_mean <- results_matrix[, 4]
    term1_correct <- results_matrix[, 5]
    term1_incorrect <- results_matrix[, 6]
    term1_correct_vs_true_mean <- results_matrix[, 7]
    term1_incorrect_vs_true_mean <- results_matrix[, 8]
    term2_correct <- results_matrix[, 9]
    term2_incorrect <- results_matrix[, 10]
    term2_correct_vs_true_mean <- results_matrix[, 11]
    term2_incorrect_vs_true_mean <- results_matrix[, 12]
  } else {
    set.seed(seed)
    es_correct <- rep(NA, nsim)
    es_incorrect <- rep(NA, nsim)
    es_correct_vs_true_mean <- rep(NA, nsim)
    es_incorrect_vs_true_mean <- rep(NA, nsim)
    term1_correct <- rep(NA, nsim)
    term1_incorrect <- rep(NA, nsim)
    term1_correct_vs_true_mean <- rep(NA, nsim)
    term1_incorrect_vs_true_mean <- rep(NA, nsim)
    term2_correct <- rep(NA, nsim)
    term2_incorrect <- rep(NA, nsim)
    term2_correct_vs_true_mean <- rep(NA, nsim)
    term2_incorrect_vs_true_mean <- rep(NA, nsim)

    for(i in 1:nsim){
      # Generate fresh forecast samples for this iteration
      thetas <- lapply(FUN = draw_one_dirichlet, X = rep(1, n_samp), alpha = alpha)
      thetas_wrong <- lapply(FUN = draw_one_dirichlet, X = rep(1, n_samp), alpha = alpha_wrong)

      # Generate observation
      prob_vec <- as.numeric(brms::rdirichlet(1, alpha))
      C <- rmultinom(n = 1, size = n_counts, prob = prob_vec)
      p <- C / n_counts

      # Apply multinomial sampling to forecasts
      samp_multinomial_counts <- do.call(cbind,
        lapply(thetas, function(theta) rmultinom(n = N_multinomial, size = n_counts, prob = theta))
      )
      samp_multinomial_counts_wrong <- do.call(cbind,
        lapply(thetas_wrong, function(theta) rmultinom(n = N_multinomial, size = n_counts, prob = theta))
      )
      p_matrix <- samp_multinomial_counts / n_counts
      p_matrix_wrong <- samp_multinomial_counts_wrong / n_counts

      # Score forecasts using manual energy score (returns ES, term1, term2)
      es_result_correct <- manual_energy_score(y = as.numeric(p), dat = p_matrix)
      es_result_incorrect <- manual_energy_score(y = as.numeric(p), dat = p_matrix_wrong)
      # Score forecasts against true mean as observation
      es_result_correct_vs_true_mean <- manual_energy_score(y = true_mean, dat = p_matrix)
      es_result_incorrect_vs_true_mean <- manual_energy_score(y = true_mean, dat = p_matrix_wrong)

      # Extract ES values
      es_correct[i] <- es_result_correct$ES
      es_incorrect[i] <- es_result_incorrect$ES
      es_correct_vs_true_mean[i] <- es_result_correct_vs_true_mean$ES
      es_incorrect_vs_true_mean[i] <- es_result_incorrect_vs_true_mean$ES

      # Extract term1 values
      term1_correct[i] <- es_result_correct$term1
      term1_incorrect[i] <- es_result_incorrect$term1
      term1_correct_vs_true_mean[i] <- es_result_correct_vs_true_mean$term1
      term1_incorrect_vs_true_mean[i] <- es_result_incorrect_vs_true_mean$term1

      # Extract term2 values
      term2_correct[i] <- es_result_correct$term2
      term2_incorrect[i] <- es_result_incorrect$term2
      term2_correct_vs_true_mean[i] <- es_result_correct_vs_true_mean$term2
      term2_incorrect_vs_true_mean[i] <- es_result_incorrect_vs_true_mean$term2
    }
  }
  return(list(
    es_correct = es_correct, es_incorrect = es_incorrect,
    es_correct_vs_true_mean = es_correct_vs_true_mean,
    es_incorrect_vs_true_mean = es_incorrect_vs_true_mean,
    term1_correct = term1_correct, term1_incorrect = term1_incorrect,
    term1_correct_vs_true_mean = term1_correct_vs_true_mean,
    term1_incorrect_vs_true_mean = term1_incorrect_vs_true_mean,
    term2_correct = term2_correct, term2_incorrect = term2_incorrect,
    term2_correct_vs_true_mean = term2_correct_vs_true_mean,
    term2_incorrect_vs_true_mean = term2_incorrect_vs_true_mean
  ))
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
    # Add term1 fields
    term1_correct = result$term1_correct,
    term1_incorrect = result$term1_incorrect,
    term1_correct_vs_true_mean = result$term1_correct_vs_true_mean,
    term1_incorrect_vs_true_mean = result$term1_incorrect_vs_true_mean,
    term1_diff = result$term1_correct - result$term1_incorrect,
    term1_diff_true_mean_obs = result$term1_correct_vs_true_mean - result$term1_incorrect_vs_true_mean,
    # Add term2 fields
    term2_correct = result$term2_correct,
    term2_incorrect = result$term2_incorrect,
    term2_correct_vs_true_mean = result$term2_correct_vs_true_mean,
    term2_incorrect_vs_true_mean = result$term2_incorrect_vs_true_mean,
    term2_diff = result$term2_correct - result$term2_incorrect,
    term2_diff_true_mean_obs = result$term2_correct_vs_true_mean - result$term2_incorrect_vs_true_mean,
    n_counts = n_counts,
    method = method,
    sim_label = paste0(method, "_n", n_counts)
  )
}

run_simulation_experiment <- function(alpha, alpha_wrong, n_samp = 100,
                                      n_counts_vec = c(1, 2, 3, 4, 5, 10, 25, 50, 100, 500),
                                      nsim = 1000, base_seed = 42,
                                      parallel = TRUE, n_cores = parallel::detectCores() - 2) {
  results_list <- list()
  for (i in seq_along(n_counts_vec)) {
    n <- n_counts_vec[i]
    cat("  n_counts =", n, "...")
    result_nonmn <- es_sampling(alpha = alpha, alpha_wrong = alpha_wrong, n_samp = n_samp,
                                seed = base_seed + i, nsim = nsim, n_counts = n,
                                parallel = parallel, n_cores = n_cores)
    result_mn <- es_mn_sampling(alpha = alpha, alpha_wrong = alpha_wrong, n_samp = n_samp,
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

# Visualize alpha distributions on ternary plot
plot_alpha_ternary <- function(alpha, alpha_wrong, n_samp = 100) {
  # Generate samples
  thetas <- lapply(FUN = draw_one_dirichlet, X = rep(1, n_samp), alpha = alpha)
  thetas_wrong <- lapply(FUN = draw_one_dirichlet, X = rep(1, n_samp), alpha = alpha_wrong)

  # Calculate mean and mode
  theta_true <- alpha / sum(alpha)
  theta_mode <- (alpha - 1) / (sum(alpha) - length(alpha))

  par(mar = c(0.3, 0.3, 1.3, 0.3))

  # Plot true alpha distribution
  TernaryPlot(alab = "Variant A prevalence \u2192",
              blab = "\u2190 Variant B prevalence",
              clab = "Variant C prevalence \u2192",
              main = paste0("True: Dir(", paste(alpha, collapse = ", "), ") vs Wrong: Dir(", paste(alpha_wrong, collapse = ", "), ")"),
              region = Ternary:::ternRegionDefault / 100,
              point = "right", lab.cex = 0.8, grid.minor.lines = 0,
              grid.lty = "solid", col = rgb(0.9, 0.9, 0.9), grid.col = "white",
              axis.col = rgb(0.6, 0.6, 0.6), ticks.col = rgb(0.6, 0.6, 0.6),
              axis.rotate = FALSE,
              padding = 0.08)

  # Color by true alpha density
  cols <- TernaryPointValues(Func = dd_func, alpha = alpha)
  ColourTernary(cols, spectrum = rev(hcl.colors(10, palette = "viridis", alpha = 0.6)))

  # Add points
  AddToTernary(graphics::points, thetas, pch = 20, cex = 0.6, col = rgb(0, 0, 1, 0.3))
  AddToTernary(graphics::points, thetas_wrong, pch = 20, cex = 0.6, col = rgb(1, 0, 0, 0.3))
  AddToTernary(graphics::points, theta_true, pch = 21, cex = 1.0, col = "blue", bg = "lightblue", lwd = 2)

  # Add legend
  legend("topright",
         legend = c("True alpha samples", "Wrong alpha samples", "True mean"),
         pch = c(20, 20, 21),
         col = c(rgb(0, 0, 1, 0.5), rgb(1, 0, 0, 0.5), "blue"),
         pt.bg = c(NA, NA, "lightblue"),
         pt.cex = c(1, 1, 1.5),
         bty = "n")
}

# Save base R graphics plot (for ternary plots)
save_ternary_plot <- function(plot_func, scenario_name, alpha, alpha_wrong,
                              n_samp = 100,
                              results_dir = "./results/plots",
                              width = 800, height = 800) {
  # Create alpha-specific subdirectory
  alpha_str <- paste(alpha, collapse = "_")
  alpha_dir <- file.path(results_dir, paste0("alpha_", alpha_str))

  if (!dir.exists(alpha_dir)) {
    dir.create(alpha_dir, recursive = TRUE)
  }

  # Create filename with alpha (no timestamp - overwrites existing)
  filename <- file.path(alpha_dir, paste0(scenario_name, "_ternary_alpha_", alpha_str, ".png"))

  # Save using base R graphics
  png(filename, width = width, height = height, res = 150)
  plot_func(alpha, alpha_wrong, n_samp)
  dev.off()

  cat("Ternary plot saved to:", filename, "\n")
  return(filename)
}

# Define ADDITIONAL extreme scenarios with SYMMETRIC true alpha
# True alpha = c(10, 10, 10): symmetric, centered at (1/3, 1/3, 1/3), sum=30
# These extend beyond the previous 5 scenarios with even more extreme dispersion
scenarios <- list(
  # Scenario 6: Ultra-dispersion (extremely underconfident, 100x less concentrated)
  list(
    name = "scenario6_ultra_dispersed",
    alpha = c(10, 10, 10),
    alpha_wrong = c(0.1, 0.1, 0.1),  # sum=0.3 (100x less concentrated), same mean
    base_seed = 132
  ),
  # Scenario 7: Ultra-narrow (extremely overconfident, 100x more concentrated)
  list(
    name = "scenario7_ultra_narrow",
    alpha = c(10, 10, 10),
    alpha_wrong = c(1000, 1000, 1000),  # sum=3000 (100x more concentrated), same mean
    base_seed = 142
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
  cat("========================================\n")

  # Run simulations (forecasts now generated fresh for each iteration)
  df_results <- run_simulation_experiment(
    alpha = scenario_config$alpha,
    alpha_wrong = scenario_config$alpha_wrong,
    n_samp = n_samp,
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

  # Save ternary plot
  save_ternary_plot(plot_alpha_ternary, scenario_config$name,
                   scenario_config$alpha, scenario_config$alpha_wrong,
                   n_samp = n_samp)

  cat("Completed:", scenario_config$name, "\n")
  return(df_results)
}

# Record start time
start_time <- Sys.time()

# Detect if running from RStudio or Positron (fork safety issues on macOS)
in_rstudio <- Sys.getenv("RSTUDIO") == "1"
in_positron <- Sys.getenv("POSITRON") == "1"
in_gui <- in_rstudio || in_positron

if (in_gui) {
  ide_name <- if (in_positron) "Positron" else "RStudio"
  cat("\n=========================================\n")
  cat("DETECTED", ide_name, "ENVIRONMENT\n")
  cat("=========================================\n")
  cat("Running scenarios SEQUENTIALLY to avoid forking issues on macOS\n")
  cat("(Using parallelization within each scenario for speedup)\n")
  cat("\nFor maximum speed, run from terminal:\n")
  cat("  Rscript src/run_simulations_alpha_symmetric.R\n\n")

  # Run scenarios sequentially (but each uses parallel processing internally)
  # Use fewer cores to be conservative in GUI environments
  n_cores_safe <- max(1, min(n_cores_inner, 4))
  cat("Using", n_cores_safe, "cores per scenario (conservative setting for GUI)\n\n")

  all_results <- lapply(scenarios, function(scenario_config) {
    run_scenario(scenario_config,
                 n_samp = 100,
                 nsim = 1000,
                 n_cores_inner = n_cores_safe)
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
combined_file <- "./results/data/alpha_10_10_10/additional_scenarios_combined.rds"
saveRDS(combined_results, combined_file)
cat("Combined results saved to:", combined_file, "\n")

cat("\nDone! Check the results/ directory for output files.\n")
cat("Individual scenario files saved in: ./results/data/alpha_10_10_10/\n")
