#!/usr/bin/env Rscript
# Data Distribution Visualizations for All Simulation Scenarios
# Creates comprehensive visualizations of Dirichlet and Dirichlet-Multinomial distributions
#
# Plot types:
# 1. Single ternary with Dirichlet density overlay
# 2. Side-by-side comparison ternaries (correct vs wrong)
# 3. Three-panel dispersion comparisons (dispersed/correct/narrow)
# 4. Multinomial outcome distributions for various N

library(extraDistr)  # For ddirmnom, rdirmnom
library(brms)        # For rdirichlet, ddirichlet
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
N_values <- c(1, 2, 3, 4, 5, 10, 25, 50, 100, 500)
n_samp <- 100  # Number of samples for visualization
output_dir <- "./results/plots/data_dists"

# Ensure output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

cat("\n=========================================\n")
cat("DATA DISTRIBUTION VISUALIZATIONS\n")
cat("=========================================\n\n")

# ============================================================================
# Helper Functions
# ============================================================================

# Create Dirichlet density grid overlay
create_dirichlet_grid <- function(alpha, grid_res = 150) {
  grid_coords <- expand.grid(
    a = seq(0, 1, length.out = grid_res),
    b = seq(0, 1, length.out = grid_res)
  )
  grid_coords$c <- 1 - grid_coords$a - grid_coords$b

  # Filter to simplex interior
  grid_coords <- grid_coords[
    grid_coords$c >= 0.001 & grid_coords$c <= 0.999 &
    grid_coords$a >= 0.001 & grid_coords$a <= 0.999 &
    grid_coords$b >= 0.001 & grid_coords$b <= 0.999,
  ]

  # Renormalize
  row_sums <- grid_coords$a + grid_coords$b + grid_coords$c
  grid_coords$a <- grid_coords$a / row_sums
  grid_coords$b <- grid_coords$b / row_sums
  grid_coords$c <- grid_coords$c / row_sums

  # Calculate Dirichlet density using proper function from helper-functions.R
  grid_coords$density <- dd_func(grid_coords$a, grid_coords$b, grid_coords$c, alpha)

  # Map to colors (mako palette)
  n_dens_colors <- 100
  dens_colors <- viridisLite::mako(n_dens_colors, alpha = 0.8, begin = 0, end = 1, direction = -1)

  # Handle edge case where density values are too similar (e.g., very large alpha)
  if (max(grid_coords$density) - min(grid_coords$density) < 1e-10) {
    grid_coords$color <- rep(dens_colors[50], nrow(grid_coords))
  } else {
    dens_breaks <- seq(min(grid_coords$density),
                       max(grid_coords$density),
                       length.out = n_dens_colors + 1)
    grid_coords$color <- dens_colors[cut(grid_coords$density,
                                         breaks = dens_breaks,
                                         labels = FALSE,
                                         include.lowest = TRUE)]
  }

  return(grid_coords)
}

# Plot single ternary with Dirichlet overlay and multinomial outcomes
plot_single_ternary_with_dm <- function(alpha, n, alpha_label, filename) {
  cat("  Creating ternary for", alpha_label, "N=", n, "\n")

  # Get multinomial sample space
  all_outcomes_rows <- enumerate_multinomial_sample_space(n)

  # Calculate DM probabilities
  probs <- sapply(1:nrow(all_outcomes_rows), function(i) {
    ddirmnom(x = all_outcomes_rows[i, ], size = n, alpha = alpha)
  })

  # Create Dirichlet density grid
  grid_coords <- create_dirichlet_grid(alpha)

  # Color points by probability
  n_colors <- 100
  color_palette <- viridisLite::viridis(n_colors, begin = 0.1, end = 0.95)

  if (max(probs) - min(probs) < 1e-10) {
    point_colors <- rep(color_palette[50], length(probs))
  } else {
    prob_breaks <- seq(min(probs), max(probs), length.out = n_colors + 1)
    point_colors <- color_palette[cut(probs, breaks = prob_breaks,
                                     labels = FALSE, include.lowest = TRUE)]
  }

  # Create plot
  png(file.path(output_dir, filename), width = 1100, height = 800, res = 150)

  layout(matrix(c(1, 2), nrow = 1), widths = c(3.0, 1.5))
  par(mar = c(2, 2, 3, 4))

  TernaryPlot(
    alab = "Category A \u2192",
    blab = "\u2190 Category B",
    clab = "Category C \u2192",
    main = paste0(alpha_label, " (N=", n, ")"),
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

  # Add Dirichlet density blanket
  AddToTernary(graphics::points,
               cbind(grid_coords$a, grid_coords$b, grid_coords$c),
               pch = 15,
               col = grid_coords$color,
               cex = 0.5)

  # Add multinomial outcomes
  outcomes_props <- t(apply(all_outcomes_rows, 1, function(x) x / sum(x)))
  AddToTernary(
    graphics::points,
    outcomes_props,
    pch = 21,
    cex = 1.5,
    col = point_colors,
    bg = point_colors,
    lwd = 1.2
  )

  # Add legend
  par(mar = c(2, 0, 3, 2))
  plot.new()
  legend("left",
         legend = c("Dirichlet Prior", "Multinomial Outcomes"),
         pch = c(15, 21),
         pt.cex = c(1, 1.5),
         col = c(viridisLite::mako(1, alpha = 0.8), color_palette[50]),
         pt.bg = c(NA, color_palette[50]),
         bty = "n",
         cex = 1.0,
         title = "Distribution")

  dev.off()
}

# Plot comparison ternaries (correct vs wrong)
plot_comparison_ternaries <- function(alpha_correct, alpha_wrong, n,
                                     correct_label, wrong_label, filename) {
  cat("  Creating comparison:", correct_label, "vs", wrong_label, "N=", n, "\n")

  # Get multinomial sample space
  all_outcomes_rows <- enumerate_multinomial_sample_space(n)

  # Calculate DM probabilities
  probs_correct <- sapply(1:nrow(all_outcomes_rows), function(i) {
    ddirmnom(x = all_outcomes_rows[i, ], size = n, alpha = alpha_correct)
  })
  probs_wrong <- sapply(1:nrow(all_outcomes_rows), function(i) {
    ddirmnom(x = all_outcomes_rows[i, ], size = n, alpha = alpha_wrong)
  })

  # Create Dirichlet density grids
  grid_correct <- create_dirichlet_grid(alpha_correct)
  grid_wrong <- create_dirichlet_grid(alpha_wrong)

  # Color points
  n_colors <- 100
  color_palette <- viridisLite::viridis(n_colors, begin = 0.1, end = 0.95)

  if (max(probs_correct) - min(probs_correct) < 1e-10) {
    colors_correct <- rep(color_palette[50], length(probs_correct))
  } else {
    prob_breaks <- seq(min(probs_correct), max(probs_correct), length.out = n_colors + 1)
    colors_correct <- color_palette[cut(probs_correct, breaks = prob_breaks,
                                       labels = FALSE, include.lowest = TRUE)]
  }

  if (max(probs_wrong) - min(probs_wrong) < 1e-10) {
    colors_wrong <- rep(color_palette[50], length(probs_wrong))
  } else {
    prob_breaks <- seq(min(probs_wrong), max(probs_wrong), length.out = n_colors + 1)
    colors_wrong <- color_palette[cut(probs_wrong, breaks = prob_breaks,
                                     labels = FALSE, include.lowest = TRUE)]
  }

  # Create 2-panel plot
  png(file.path(output_dir, filename), width = 1600, height = 800, res = 150)

  par(mfrow = c(1, 2), mar = c(2, 2, 3, 2))

  outcomes_props <- t(apply(all_outcomes_rows, 1, function(x) x / sum(x)))

  # Left panel: Correct
  TernaryPlot(
    alab = "A \u2192",
    blab = "\u2190 B",
    clab = "C \u2192",
    main = paste0("Correct: ", correct_label, " (N=", n, ")"),
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

  AddToTernary(graphics::points,
               cbind(grid_correct$a, grid_correct$b, grid_correct$c),
               pch = 15, col = grid_correct$color, cex = 0.5)

  AddToTernary(graphics::points, outcomes_props,
               pch = 21, cex = 2.2, col = colors_correct,
               bg = colors_correct, lwd = 1.2)

  # Right panel: Wrong
  TernaryPlot(
    alab = "A \u2192",
    blab = "\u2190 B",
    clab = "C \u2192",
    main = paste0("Misspecified: ", wrong_label, " (N=", n, ")"),
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

  AddToTernary(graphics::points,
               cbind(grid_wrong$a, grid_wrong$b, grid_wrong$c),
               pch = 15, col = grid_wrong$color, cex = 0.5)

  AddToTernary(graphics::points, outcomes_props,
               pch = 21, cex = 2.2, col = colors_wrong,
               bg = colors_wrong, lwd = 1.2)

  dev.off()
}

# Plot 3-panel dispersion comparison
plot_dispersion_comparison <- function(alpha_dispersed, alpha_correct, alpha_narrow, n,
                                      scenario_name, filename) {
  cat("  Creating 3-panel dispersion for", scenario_name, "N=", n, "\n")

  # Get multinomial sample space
  all_outcomes_rows <- enumerate_multinomial_sample_space(n)
  outcomes_props <- t(apply(all_outcomes_rows, 1, function(x) x / sum(x)))

  # Calculate probabilities for all three
  probs_dispersed <- sapply(1:nrow(all_outcomes_rows), function(i) {
    ddirmnom(x = all_outcomes_rows[i, ], size = n, alpha = alpha_dispersed)
  })
  probs_correct <- sapply(1:nrow(all_outcomes_rows), function(i) {
    ddirmnom(x = all_outcomes_rows[i, ], size = n, alpha = alpha_correct)
  })
  probs_narrow <- sapply(1:nrow(all_outcomes_rows), function(i) {
    ddirmnom(x = all_outcomes_rows[i, ], size = n, alpha = alpha_narrow)
  })

  # Create grids
  grid_dispersed <- create_dirichlet_grid(alpha_dispersed)
  grid_correct <- create_dirichlet_grid(alpha_correct)
  grid_narrow <- create_dirichlet_grid(alpha_narrow)

  # Color points
  n_colors <- 100
  color_palette <- viridisLite::viridis(n_colors, begin = 0.1, end = 0.95)

  # Helper to color points
  color_probs <- function(probs) {
    if (max(probs) - min(probs) < 1e-10) {
      rep(color_palette[50], length(probs))
    } else {
      prob_breaks <- seq(min(probs), max(probs), length.out = n_colors + 1)
      color_palette[cut(probs, breaks = prob_breaks, labels = FALSE, include.lowest = TRUE)]
    }
  }

  colors_dispersed <- color_probs(probs_dispersed)
  colors_correct <- color_probs(probs_correct)
  colors_narrow <- color_probs(probs_narrow)

  # Create 3-panel plot
  png(file.path(output_dir, filename), width = 2200, height = 800, res = 150)

  par(mfrow = c(1, 3), mar = c(2, 2, 3, 2))

  # Left: Dispersed
  TernaryPlot(
    alab = "A \u2192", blab = "\u2190 B", clab = "C \u2192",
    main = paste0("Dispersed: \u03B1=(", paste(alpha_dispersed, collapse = ", "), ")"),
    region = Ternary:::ternRegionDefault / 100,
    point = "right", lab.cex = 0.8, grid.minor.lines = 0,
    grid.lty = "solid", col = rgb(0.9, 0.9, 0.9), grid.col = "white",
    axis.col = rgb(0.6, 0.6, 0.6), ticks.col = rgb(0.6, 0.6, 0.6),
    axis.rotate = FALSE, padding = 0.08
  )
  AddToTernary(graphics::points,
               cbind(grid_dispersed$a, grid_dispersed$b, grid_dispersed$c),
               pch = 15, col = grid_dispersed$color, cex = 0.5)
  AddToTernary(graphics::points, outcomes_props,
               pch = 21, cex = 2.2, col = colors_dispersed,
               bg = colors_dispersed, lwd = 1.2)

  # Middle: Correct
  TernaryPlot(
    alab = "A \u2192", blab = "\u2190 B", clab = "C \u2192",
    main = paste0("Correct: \u03B1=(", paste(alpha_correct, collapse = ", "), ")"),
    region = Ternary:::ternRegionDefault / 100,
    point = "right", lab.cex = 0.8, grid.minor.lines = 0,
    grid.lty = "solid", col = rgb(0.9, 0.9, 0.9), grid.col = "white",
    axis.col = rgb(0.6, 0.6, 0.6), ticks.col = rgb(0.6, 0.6, 0.6),
    axis.rotate = FALSE, padding = 0.08
  )
  AddToTernary(graphics::points,
               cbind(grid_correct$a, grid_correct$b, grid_correct$c),
               pch = 15, col = grid_correct$color, cex = 0.5)
  AddToTernary(graphics::points, outcomes_props,
               pch = 21, cex = 2.2, col = colors_correct,
               bg = colors_correct, lwd = 1.2)

  # Right: Narrow
  TernaryPlot(
    alab = "A \u2192", blab = "\u2190 B", clab = "C \u2192",
    main = paste0("Narrow: \u03B1=(", paste(alpha_narrow, collapse = ", "), ")"),
    region = Ternary:::ternRegionDefault / 100,
    point = "right", lab.cex = 0.8, grid.minor.lines = 0,
    grid.lty = "solid", col = rgb(0.9, 0.9, 0.9), grid.col = "white",
    axis.col = rgb(0.6, 0.6, 0.6), ticks.col = rgb(0.6, 0.6, 0.6),
    axis.rotate = FALSE, padding = 0.08
  )
  AddToTernary(graphics::points,
               cbind(grid_narrow$a, grid_narrow$b, grid_narrow$c),
               pch = 15, col = grid_narrow$color, cex = 0.5)
  AddToTernary(graphics::points, outcomes_props,
               pch = 21, cex = 2.2, col = colors_narrow,
               bg = colors_narrow, lwd = 1.2)

  # Overall title
  mtext(paste0("N=", n, ": Dispersion Comparison"),
        side = 3, line = -2, outer = TRUE, cex = 1.3, font = 2)

  dev.off()
}

# ============================================================================
# Define Scenarios
# ============================================================================

# SYMMETRIC SCENARIOS
symmetric_base <- c(10, 10, 10)
symmetric_scenarios <- list(
  list(
    name = "scenario1_mislocation_symmetric",
    alpha_correct = symmetric_base,
    alpha_wrong = c(5, 10, 15),
    alpha_dispersed = c(2, 2, 2),
    alpha_narrow = c(50, 50, 50)
  ),
  list(
    name = "scenario2_moderate_dispersed",
    alpha_correct = symmetric_base,
    alpha_wrong = c(2, 2, 2),
    alpha_dispersed = c(0.5, 0.5, 0.5),
    alpha_narrow = c(50, 50, 50)
  ),
  list(
    name = "scenario3_extreme_dispersed",
    alpha_correct = symmetric_base,
    alpha_wrong = c(0.5, 0.5, 0.5),
    alpha_dispersed = c(0.1, 0.1, 0.1),
    alpha_narrow = c(50, 50, 50)
  ),
  list(
    name = "scenario4_moderate_narrow",
    alpha_correct = symmetric_base,
    alpha_wrong = c(50, 50, 50),
    alpha_dispersed = c(2, 2, 2),
    alpha_narrow = c(200, 200, 200)
  ),
  list(
    name = "scenario5_extreme_narrow",
    alpha_correct = symmetric_base,
    alpha_wrong = c(200, 200, 200),
    alpha_dispersed = c(2, 2, 2),
    alpha_narrow = c(1000, 1000, 1000)
  ),
  list(
    name = "scenario6_ultra_dispersed",
    alpha_correct = symmetric_base,
    alpha_wrong = c(0.1, 0.1, 0.1),
    alpha_dispersed = c(0.01, 0.01, 0.01),
    alpha_narrow = c(50, 50, 50)
  ),
  list(
    name = "scenario7_ultra_narrow",
    alpha_correct = symmetric_base,
    alpha_wrong = c(1000, 1000, 1000),
    alpha_dispersed = c(2, 2, 2),
    alpha_narrow = c(5000, 5000, 5000)
  )
)

# ASYMMETRIC SCENARIOS
asymmetric_base <- c(2, 3, 15)
asymmetric_scenarios <- list(
  list(
    name = "scenario1_misspecified",
    alpha_correct = asymmetric_base,
    alpha_wrong = c(15, 3, 2),
    alpha_dispersed = c(0.4, 0.6, 3),
    alpha_narrow = c(20, 30, 150)
  ),
  list(
    name = "scenario2_narrow_calibrated",
    alpha_correct = asymmetric_base,
    alpha_wrong = c(200, 300, 1500),
    alpha_dispersed = c(0.4, 0.6, 3),
    alpha_narrow = c(1000, 1500, 7500)
  ),
  list(
    name = "scenario3_dispersed_calibrated",
    alpha_correct = asymmetric_base,
    alpha_wrong = c(0.4, 0.6, 3),
    alpha_dispersed = c(0.2, 0.3, 1.5),
    alpha_narrow = c(20, 30, 150)
  ),
  list(
    name = "scenario4_less_misspecified",
    alpha_correct = asymmetric_base,
    alpha_wrong = c(4, 5, 11),
    alpha_dispersed = c(0.4, 0.6, 3),
    alpha_narrow = c(20, 30, 150)
  )
)

# ============================================================================
# Generate All Plots
# ============================================================================

cat("Generating visualizations for",
    length(symmetric_scenarios) + length(asymmetric_scenarios), "scenarios\n")
cat("N values:", paste(N_values, collapse = ", "), "\n\n")

# Process symmetric scenarios
cat("SYMMETRIC SCENARIOS:\n")
for (scenario in symmetric_scenarios) {
  cat("\n", scenario$name, ":\n", sep="")

  for (n in N_values) {
    # Single ternary plots
    plot_single_ternary_with_dm(
      scenario$alpha_correct, n,
      paste0("\u03B1=(", paste(scenario$alpha_correct, collapse = ", "), ")"),
      paste0("single_ternary_symmetric_", scenario$name, "_N", n, ".png")
    )

    # Comparison ternaries
    plot_comparison_ternaries(
      scenario$alpha_correct, scenario$alpha_wrong, n,
      paste0("\u03B1=(", paste(scenario$alpha_correct, collapse = ", "), ")"),
      paste0("\u03B1=(", paste(scenario$alpha_wrong, collapse = ", "), ")"),
      paste0("comparison_ternary_symmetric_", scenario$name, "_N", n, ".png")
    )

    # 3-panel dispersion comparison
    plot_dispersion_comparison(
      scenario$alpha_dispersed, scenario$alpha_correct, scenario$alpha_narrow, n,
      scenario$name,
      paste0("dispersion_comparison_symmetric_", scenario$name, "_N", n, ".png")
    )
  }
}

# Process asymmetric scenarios
cat("\nASYMMETRIC SCENARIOS:\n")
for (scenario in asymmetric_scenarios) {
  cat("\n", scenario$name, ":\n", sep="")

  for (n in N_values) {
    # Single ternary plots
    plot_single_ternary_with_dm(
      scenario$alpha_correct, n,
      paste0("\u03B1=(", paste(scenario$alpha_correct, collapse = ", "), ")"),
      paste0("single_ternary_asymmetric_", scenario$name, "_N", n, ".png")
    )

    # Comparison ternaries
    plot_comparison_ternaries(
      scenario$alpha_correct, scenario$alpha_wrong, n,
      paste0("\u03B1=(", paste(scenario$alpha_correct, collapse = ", "), ")"),
      paste0("\u03B1=(", paste(scenario$alpha_wrong, collapse = ", "), ")"),
      paste0("comparison_ternary_asymmetric_", scenario$name, "_N", n, ".png")
    )

    # 3-panel dispersion comparison
    plot_dispersion_comparison(
      scenario$alpha_dispersed, scenario$alpha_correct, scenario$alpha_narrow, n,
      scenario$name,
      paste0("dispersion_comparison_asymmetric_", scenario$name, "_N", n, ".png")
    )
  }
}

cat("\n=========================================\n")
cat("COMPLETE!\n")
total_plots <- (length(symmetric_scenarios) + length(asymmetric_scenarios)) * length(N_values) * 3
cat("Generated", total_plots, "plots (3 types per scenario/N combination)\n")
cat("Output directory:", output_dir, "\n")
cat("=========================================\n")
