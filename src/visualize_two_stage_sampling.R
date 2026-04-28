#!/usr/bin/env Rscript
# Two-Stage Sampling Visualization: Theoretical DM vs Empirical Outcomes
# Shows theoretical Dirichlet-Multinomial distributions alongside
# empirical outcomes from actual two-stage sampling (Dirichlet -> Multinomial)
#
# Purpose: Explains why discrimination is poor at small N
# - At small N: Empirical outcomes from correct and wrong forecasts overlap heavily
# - At large N: Empirical outcomes become more distinguishable
#
# Output: 2x2 panel plots for each scenario and N value
# - TOP LEFT: Theoretical DM (correct alpha)
# - TOP RIGHT: Empirical two-stage sampling (correct alpha)
# - BOTTOM LEFT: Theoretical DM (wrong alpha)
# - BOTTOM RIGHT: Empirical two-stage sampling (wrong alpha)

library(extraDistr)  # For ddirmnom, rdirmnom
library(brms)        # For rdirichlet
library(Ternary)     # For ternary plots
library(viridisLite) # For color palettes

# Source helper functions
if (file.exists("./src/helper-functions.R")) {
  source("./src/helper-functions.R")
} else if (file.exists("helper-functions.R")) {
  source("helper-functions.R")
} else {
  stop("Cannot find helper-functions.R")
}

# Configuration
N_theta_samples <- 100   # Number of theta samples from Dirichlet
N_multinomial <- 100     # Number of multinomial samples per theta
N_values <- c(1, 2, 3, 4, 5, 10, 25, 50, 100, 500)
output_dir <- "./results/plots/data_dists"

# Ensure output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created directory:", output_dir, "\n")
}

# Function to plot a single ternary with adaptive point sizing
plot_one_ternary <- function(dat_outcomes, point_colors, title, grid_coords = NULL) {
  # Convert outcomes to proportions
  dat_props <- t(apply(dat_outcomes, 2, function(col) col / sum(col)))

  # Adjust point size based on N
  point_size <- if (n >= 100) {
    0.8  # Much smaller for large N
  } else if (n >= 50) {
    1.2  # Medium-small for medium-large N
  } else {
    2.0  # Normal size for small N
  }

  line_width <- if (n >= 100) {
    0.8
  } else if (n >= 50) {
    1.0
  } else {
    1.5
  }

  # Create ternary plot
  TernaryPlot(
    alab = "A \u2192",
    blab = "\u2190 B",
    clab = "C \u2192",
    main = title,
    region = Ternary:::ternRegionDefault / 100,
    point = "right",
    lab.cex = 0.8,
    grid.minor.lines = 0,
    grid.lty = "solid",
    col = rgb(0.9, 0.9, 0.9),
    grid.col = "white",
    axis.col = rgb(0.6, 0.6, 0.6),
    ticks.col = rgb(0.6, 0.6, 0.6),
    axis.rotate = FALSE,
    padding = 0.08
  )

  # Add Dirichlet density blanket if provided
  if (!is.null(grid_coords)) {
    AddToTernary(graphics::points,
                 cbind(grid_coords$a, grid_coords$b, grid_coords$c),
                 pch = 15,
                 col = grid_coords$color,
                 cex = 0.5)
  }

  # Add empirical points
  AddToTernary(
    graphics::points,
    dat_props,
    pch = 21,
    cex = point_size,
    col = point_colors,
    bg = point_colors,
    lwd = line_width
  )
}

# Function to create 2x2 comparison plot
plot_two_stage_correct_vs_wrong <- function(n, alpha_correct, alpha_wrong,
                                            correct_label, wrong_label, filename) {

  cat("Creating 2x2 comparison for N=", n, "\n", sep="")

  # Stage 1: Sample thetas from Dirichlet
  set.seed(42)
  thetas_correct <- lapply(1:N_theta_samples, function(i) {
    as.vector(rdirichlet(1, alpha = alpha_correct))
  })
  thetas_wrong <- lapply(1:N_theta_samples, function(i) {
    as.vector(rdirichlet(1, alpha = alpha_wrong))
  })

  # Stage 2: Sample counts from Multinomial for each theta
  # Correct forecasts
  empirical_counts_correct <- do.call(cbind,
    lapply(thetas_correct, function(theta) {
      rmultinom(n = N_multinomial, size = n, prob = theta)
    })
  )

  # Wrong forecasts
  empirical_counts_wrong <- do.call(cbind,
    lapply(thetas_wrong, function(theta) {
      rmultinom(n = N_multinomial, size = n, prob = theta)
    })
  )

  # Get all possible outcomes for theoretical plots
  # enumerate returns outcomes as rows (n_outcomes × 3), transpose for plotting (3 × n_outcomes)
  all_outcomes_rows <- enumerate_multinomial_sample_space(n)
  all_outcomes <- t(all_outcomes_rows)

  # Calculate theoretical DM probabilities (iterate over rows of original matrix)
  prob_correct <- sapply(1:nrow(all_outcomes_rows), function(i) {
    ddirmnom(x = all_outcomes_rows[i, ], size = n, alpha = alpha_correct)
  })
  prob_wrong <- sapply(1:nrow(all_outcomes_rows), function(i) {
    ddirmnom(x = all_outcomes_rows[i, ], size = n, alpha = alpha_wrong)
  })

  # Create Dirichlet density grid for correct alpha
  grid_res <- 150
  grid_coords_correct <- expand.grid(
    a = seq(0, 1, length.out = grid_res),
    b = seq(0, 1, length.out = grid_res)
  )
  grid_coords_correct$c <- 1 - grid_coords_correct$a - grid_coords_correct$b
  grid_coords_correct <- grid_coords_correct[
    grid_coords_correct$c >= 0.001 & grid_coords_correct$c <= 0.999 &
    grid_coords_correct$a >= 0.001 & grid_coords_correct$a <= 0.999 &
    grid_coords_correct$b >= 0.001 & grid_coords_correct$b <= 0.999,
  ]

  # Renormalize
  row_sums <- grid_coords_correct$a + grid_coords_correct$b + grid_coords_correct$c
  grid_coords_correct$a <- grid_coords_correct$a / row_sums
  grid_coords_correct$b <- grid_coords_correct$b / row_sums
  grid_coords_correct$c <- grid_coords_correct$c / row_sums

  # Calculate Dirichlet density using proper function from helper-functions.R
  grid_coords_correct$density <- dd_func(grid_coords_correct$a, grid_coords_correct$b,
                                         grid_coords_correct$c, alpha_correct)

  # Map to colors (mako palette)
  n_dens_colors <- 100
  dens_colors <- viridisLite::mako(n_dens_colors, alpha = 0.8, begin = 0, end = 1, direction = -1)

  # Handle edge case where density values are too similar
  if (max(grid_coords_correct$density) - min(grid_coords_correct$density) < 1e-10) {
    grid_coords_correct$color <- rep(dens_colors[50], nrow(grid_coords_correct))
  } else {
    dens_breaks <- seq(min(grid_coords_correct$density),
                       max(grid_coords_correct$density),
                       length.out = n_dens_colors + 1)
    grid_coords_correct$color <- dens_colors[cut(grid_coords_correct$density,
                                                 breaks = dens_breaks,
                                                 labels = FALSE,
                                                 include.lowest = TRUE)]
  }

  # Create Dirichlet density grid for wrong alpha
  grid_coords_wrong <- expand.grid(
    a = seq(0, 1, length.out = grid_res),
    b = seq(0, 1, length.out = grid_res)
  )
  grid_coords_wrong$c <- 1 - grid_coords_wrong$a - grid_coords_wrong$b
  grid_coords_wrong <- grid_coords_wrong[
    grid_coords_wrong$c >= 0.001 & grid_coords_wrong$c <= 0.999 &
    grid_coords_wrong$a >= 0.001 & grid_coords_wrong$a <= 0.999 &
    grid_coords_wrong$b >= 0.001 & grid_coords_wrong$b <= 0.999,
  ]

  # Renormalize
  row_sums <- grid_coords_wrong$a + grid_coords_wrong$b + grid_coords_wrong$c
  grid_coords_wrong$a <- grid_coords_wrong$a / row_sums
  grid_coords_wrong$b <- grid_coords_wrong$b / row_sums
  grid_coords_wrong$c <- grid_coords_wrong$c / row_sums

  # Calculate Dirichlet density using proper function from helper-functions.R
  grid_coords_wrong$density <- dd_func(grid_coords_wrong$a, grid_coords_wrong$b,
                                       grid_coords_wrong$c, alpha_wrong)

  # Map to colors (handle edge case)
  if (max(grid_coords_wrong$density) - min(grid_coords_wrong$density) < 1e-10) {
    grid_coords_wrong$color <- rep(dens_colors[50], nrow(grid_coords_wrong))
  } else {
    dens_breaks <- seq(min(grid_coords_wrong$density),
                       max(grid_coords_wrong$density),
                       length.out = n_dens_colors + 1)
    grid_coords_wrong$color <- dens_colors[cut(grid_coords_wrong$density,
                                               breaks = dens_breaks,
                                               labels = FALSE,
                                               include.lowest = TRUE)]
  }

  # Color points by probability
  n_colors <- 100
  color_palette <- viridisLite::viridis(n_colors, begin = 0.1, end = 0.95)

  # Handle case where probabilities are too similar
  if (max(prob_correct) - min(prob_correct) < 1e-10) {
    point_colors_correct <- rep(color_palette[50], length(prob_correct))
  } else {
    prob_breaks_correct <- seq(min(prob_correct), max(prob_correct), length.out = n_colors + 1)
    point_colors_correct <- color_palette[cut(prob_correct, breaks = prob_breaks_correct,
                                             labels = FALSE, include.lowest = TRUE)]
  }

  if (max(prob_wrong) - min(prob_wrong) < 1e-10) {
    point_colors_wrong <- rep(color_palette[50], length(prob_wrong))
  } else {
    prob_breaks_wrong <- seq(min(prob_wrong), max(prob_wrong), length.out = n_colors + 1)
    point_colors_wrong <- color_palette[cut(prob_wrong, breaks = prob_breaks_wrong,
                                           labels = FALSE, include.lowest = TRUE)]
  }

  # For empirical outcomes, use density-based coloring
  # Count unique outcomes and their frequencies
  unique_outcomes_correct <- unique(empirical_counts_correct, MARGIN = 2)
  freqs_correct <- apply(unique_outcomes_correct, 2, function(outcome) {
    sum(apply(empirical_counts_correct, 2, function(x) all(x == outcome)))
  })

  unique_outcomes_wrong <- unique(empirical_counts_wrong, MARGIN = 2)
  freqs_wrong <- apply(unique_outcomes_wrong, 2, function(outcome) {
    sum(apply(empirical_counts_wrong, 2, function(x) all(x == outcome)))
  })

  # Calculate overlap for diagnostics
  overlap_count <- sum(apply(unique_outcomes_correct, 2, function(correct_outcome) {
    any(apply(unique_outcomes_wrong, 2, function(wrong_outcome) {
      all(correct_outcome == wrong_outcome)
    }))
  }))
  overlap_pct <- round(100 * overlap_count / ncol(unique_outcomes_correct), 1)

  # Calculate coverage (percentage of theoretical sample space represented)
  n_theoretical_outcomes <- ncol(all_outcomes)
  coverage_correct_pct <- round(100 * ncol(unique_outcomes_correct) / n_theoretical_outcomes, 1)
  coverage_wrong_pct <- round(100 * ncol(unique_outcomes_wrong) / n_theoretical_outcomes, 1)

  cat("  Empirical (correct):", ncol(unique_outcomes_correct), "outcomes (",
      coverage_correct_pct, "% coverage)\n", sep="")
  cat("  Empirical (wrong):", ncol(unique_outcomes_wrong), "outcomes (",
      coverage_wrong_pct, "% coverage)\n", sep="")
  cat("  Overlap:", overlap_count, "outcomes appear in both (",
      overlap_pct, "% of correct outcomes)\n", sep="")

  # Color empirical points by frequency
  # Handle case where frequencies are too similar
  if (max(freqs_correct) - min(freqs_correct) < 1e-10) {
    empirical_colors_correct <- rep(color_palette[50], length(freqs_correct))
  } else {
    freq_breaks_correct <- seq(min(freqs_correct), max(freqs_correct), length.out = n_colors + 1)
    empirical_colors_correct <- color_palette[cut(freqs_correct, breaks = freq_breaks_correct,
                                                 labels = FALSE, include.lowest = TRUE)]
  }

  if (max(freqs_wrong) - min(freqs_wrong) < 1e-10) {
    empirical_colors_wrong <- rep(color_palette[50], length(freqs_wrong))
  } else {
    freq_breaks_wrong <- seq(min(freqs_wrong), max(freqs_wrong), length.out = n_colors + 1)
    empirical_colors_wrong <- color_palette[cut(freqs_wrong, breaks = freq_breaks_wrong,
                                               labels = FALSE, include.lowest = TRUE)]
  }

  # Create 2x2 panel plot
  png(file.path(output_dir, filename),
      width = 1600, height = 1600, res = 150)

  par(mfrow = c(2, 2), mar = c(1, 1, 3, 1), oma = c(0, 0, 3, 0))

  # TOP LEFT: Theoretical DM (correct)
  plot_one_ternary(
    all_outcomes,
    point_colors_correct,
    paste0("Theoretical DM\n", correct_label),
    grid_coords_correct
  )

  # TOP RIGHT: Empirical two-stage (correct)
  plot_one_ternary(
    unique_outcomes_correct,
    empirical_colors_correct,
    paste0("Empirical: Two-Stage\n(", N_theta_samples, " \u03B8 × ", N_multinomial, " MN)"),
    NULL
  )

  # BOTTOM LEFT: Theoretical DM (wrong)
  plot_one_ternary(
    all_outcomes,
    point_colors_wrong,
    paste0("Theoretical DM\nMisspecified: ", wrong_label),
    grid_coords_wrong
  )

  # BOTTOM RIGHT: Empirical two-stage (wrong)
  plot_one_ternary(
    unique_outcomes_wrong,
    empirical_colors_wrong,
    paste0("Empirical: Two-Stage\n(", N_theta_samples, " \u03B8 × ", N_multinomial, " MN)"),
    NULL
  )

  # Add overall title
  mtext(paste0("N=", n, ": ", correct_label, " vs ", wrong_label,
               " (Overlap: ", overlap_pct, "%)"),
        side = 3, line = 1, outer = TRUE, cex = 1.3, font = 2)

  dev.off()

  cat("  Saved:", filename, "\n")
}

# Define all simulation scenarios
cat("\n=========================================\n")
cat("TWO-STAGE SAMPLING VISUALIZATION\n")
cat("=========================================\n\n")

# SYMMETRIC SCENARIOS (true alpha = 10, 10, 10)
symmetric_scenarios <- list(
  list(
    name = "scenario1_mislocation_symmetric",
    alpha_correct = c(10, 10, 10),
    alpha_wrong = c(5, 10, 15),
    correct_label = "Correct: \u03B1=(10, 10, 10)",
    wrong_label = "\u03B1=(5, 10, 15)"
  ),
  list(
    name = "scenario2_moderate_dispersed",
    alpha_correct = c(10, 10, 10),
    alpha_wrong = c(2, 2, 2),
    correct_label = "Correct: \u03B1=(10, 10, 10)",
    wrong_label = "\u03B1=(2, 2, 2)"
  ),
  list(
    name = "scenario3_extreme_dispersed",
    alpha_correct = c(10, 10, 10),
    alpha_wrong = c(0.5, 0.5, 0.5),
    correct_label = "Correct: \u03B1=(10, 10, 10)",
    wrong_label = "\u03B1=(0.5, 0.5, 0.5)"
  ),
  list(
    name = "scenario4_moderate_narrow",
    alpha_correct = c(10, 10, 10),
    alpha_wrong = c(50, 50, 50),
    correct_label = "Correct: \u03B1=(10, 10, 10)",
    wrong_label = "\u03B1=(50, 50, 50)"
  ),
  list(
    name = "scenario5_extreme_narrow",
    alpha_correct = c(10, 10, 10),
    alpha_wrong = c(200, 200, 200),
    correct_label = "Correct: \u03B1=(10, 10, 10)",
    wrong_label = "\u03B1=(200, 200, 200)"
  ),
  list(
    name = "scenario6_ultra_dispersed",
    alpha_correct = c(10, 10, 10),
    alpha_wrong = c(0.1, 0.1, 0.1),
    correct_label = "Correct: \u03B1=(10, 10, 10)",
    wrong_label = "\u03B1=(0.1, 0.1, 0.1)"
  ),
  list(
    name = "scenario7_ultra_narrow",
    alpha_correct = c(10, 10, 10),
    alpha_wrong = c(1000, 1000, 1000),
    correct_label = "Correct: \u03B1=(10, 10, 10)",
    wrong_label = "\u03B1=(1000, 1000, 1000)"
  )
)

# ASYMMETRIC SCENARIOS (true alpha = 2, 3, 15)
asymmetric_scenarios <- list(
  list(
    name = "scenario1_misspecified",
    alpha_correct = c(2, 3, 15),
    alpha_wrong = c(15, 3, 2),
    correct_label = "Correct: \u03B1=(2, 3, 15)",
    wrong_label = "\u03B1=(15, 3, 2)"
  ),
  list(
    name = "scenario2_narrow_calibrated",
    alpha_correct = c(2, 3, 15),
    alpha_wrong = c(200, 300, 1500),
    correct_label = "Correct: \u03B1=(2, 3, 15)",
    wrong_label = "\u03B1=(200, 300, 1500)"
  ),
  list(
    name = "scenario3_dispersed_calibrated",
    alpha_correct = c(2, 3, 15),
    alpha_wrong = c(0.4, 0.6, 3),
    correct_label = "Correct: \u03B1=(2, 3, 15)",
    wrong_label = "\u03B1=(0.4, 0.6, 3)"
  ),
  list(
    name = "scenario4_less_misspecified",
    alpha_correct = c(2, 3, 15),
    alpha_wrong = c(4, 5, 11),
    correct_label = "Correct: \u03B1=(2, 3, 15)",
    wrong_label = "\u03B1=(4, 5, 11)"
  )
)

# Generate plots for all scenarios
cat("Generating plots for", length(symmetric_scenarios) + length(asymmetric_scenarios), "scenarios\n")
cat("N values:", paste(N_values, collapse = ", "), "\n\n")

# Process symmetric scenarios
cat("SYMMETRIC SCENARIOS:\n")
for (scenario in symmetric_scenarios) {
  cat("\n", scenario$name, ":\n", sep="")
  for (n in N_values) {
    filename <- paste0("two_stage_2x2_symmetric_", scenario$name, "_N", n, ".png")
    plot_two_stage_correct_vs_wrong(
      n = n,
      alpha_correct = scenario$alpha_correct,
      alpha_wrong = scenario$alpha_wrong,
      correct_label = scenario$correct_label,
      wrong_label = scenario$wrong_label,
      filename = filename
    )
  }
}

# Process asymmetric scenarios
cat("\nASYMMETRIC SCENARIOS:\n")
for (scenario in asymmetric_scenarios) {
  cat("\n", scenario$name, ":\n", sep="")
  for (n in N_values) {
    filename <- paste0("two_stage_2x2_asymmetric_", scenario$name, "_N", n, ".png")
    plot_two_stage_correct_vs_wrong(
      n = n,
      alpha_correct = scenario$alpha_correct,
      alpha_wrong = scenario$alpha_wrong,
      correct_label = scenario$correct_label,
      wrong_label = scenario$wrong_label,
      filename = filename
    )
  }
}

cat("\n=========================================\n")
cat("COMPLETE!\n")
cat("Generated", (length(symmetric_scenarios) + length(asymmetric_scenarios)) * length(N_values), "plots\n")
cat("Output directory:", output_dir, "\n")
cat("=========================================\n")
