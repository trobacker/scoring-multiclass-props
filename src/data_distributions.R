#!/usr/bin/env Rscript
# Visualizations of Dirichlet-Multinomial Count Distributions
# Shows the induced distribution of realized counts N under Dirichlet-Multinomial model

library(ggplot2)
library(extraDistr)
library(brms)
library(Ternary)
library(tidyr)
library(dplyr)

# Source helper functions
if (file.exists("./src/helper-functions.R")) {
  source("./src/helper-functions.R")
} else if (file.exists("helper-functions.R")) {
  source("helper-functions.R")
} else {
  stop("Cannot find helper-functions.R")
}

# ============================================================================
# CONFIGURATION: Change these parameters to explore different settings
# ============================================================================

# where to save plots
output_dir <- "./results/plots/data_dists"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Alpha values from simulation scenarios
alpha_symmetric <- c(10, 10, 10)     # Symmetric from simulations
alpha_asymmetric <- c(2, 3, 15)      # Asymmetric from simulations

# N values to visualize (matching simulations)
N_values <- c(1, 2, 3, 4, 5, 10, 25, 50, 100, 500)

# Number of samples to draw for heatmap visualization
n_samples <- 100

# ============================================================================
# PART 1: MULTIPLE TERNARY PLOTS with Color-Coded Densities
# ============================================================================

cat("Creating ternary plots for K=3 Dirichlet-Multinomial with color legends...\n\n")

# Function to create ternary plot with color-coded densities and legend
plot_ternary_with_density_legend <- function(n, alpha, plot_title, filename) {
  require(extraDistr)
  require(Ternary)

  # Enumerate all possible outcomes
  dat <- enumerate_multinomial_sample_space(n)

  # Calculate densities for all outcomes
  densities <- ddirmnom(dat, size = n, alpha = alpha)

  # Create color palette based on densities
  # Use equal-width binning for better visualization
  n_colors <- 10
  # Use viridis palette (color-blind friendly)
  color_palette <- viridisLite::viridis(n_colors, begin = 0.1, end = 0.95)

  # Get range of non-zero densities
  dens_nonzero <- densities[densities > 0]
  dens_min <- min(dens_nonzero)
  dens_max <- max(dens_nonzero)

  # Create bins with a small epsilon to ensure uniqueness
  breaks <- seq(dens_min - 1e-10, dens_max + 1e-10, length.out = n_colors + 1)

  # Assign colors based on density bins
  point_colors <- rep(rgb(0.9, 0.9, 0.9, 0.3), length(densities))  # Default for zero
  for (i in seq_along(densities)) {
    if (densities[i] > 0) {
      # Find which bin this density falls into
      bin_idx <- findInterval(densities[i], breaks, all.inside = TRUE)
      bin_idx <- min(bin_idx, n_colors)  # Ensure within range
      point_colors[i] <- color_palette[bin_idx]
    }
  }

  # Create plot
  png(file.path(output_dir, filename),
      width = 1100, height = 800, res = 150)

  # Set up layout for plot + legend
  layout(matrix(c(1, 2), nrow = 1), widths = c(3.0, 1.5))

  # Main ternary plot
  par(mar = c(1, 1, 3, 1))
  TernaryPlot(alab = "Variant A count →",
              blab = "← Variant B count ",
              clab = "Variant C count →",
              region = Ternary:::ternRegionDefault/(100/n),
              grid.lines = 5,
              point = "right", lab.cex = 0.8, grid.minor.lines = 0,
              grid.lty = "solid", col = rgb(0.9, 0.9, 0.9),
              grid.col = "white",
              axis.col = rgb(0.6, 0.6, 0.6),
              ticks.col = rgb(0.6, 0.6, 0.6),
              axis.rotate = FALSE,
              padding = 0.08)

  # Add points with color based on density
  point_size <- 2
  AddToTernary(graphics::points,
               dat,
               pch = 19,
               cex = point_size,
               col = point_colors)

  # Add title
  mtext(plot_title, side = 3, line = 1, cex = 1.2, font = 2)

  # Add legend
  par(mar = c(2, 2, 3, 4))
  plot.new()

  # Get density range for legend
  dens_range <- range(densities[densities > 0])
  legend_vals <- seq(dens_range[1], dens_range[2], length.out = n_colors)

  # Create color bar
  legend_y <- seq(0.15, 0.85, length.out = n_colors)

  for (i in 1:n_colors) {
    rect(0.20, legend_y[i], 0.40,
         legend_y[min(i+1, n_colors)],
         col = color_palette[i], border = NA)
  }

  # Add legend labels - positioned to avoid cutoff
  text(0.30, 0.92, "Probability", cex = 0.9, font = 2)
  text(0.50, 0.85, sprintf("%.4f", dens_range[2]), pos = 4, cex = 0.75)
  text(0.50, 0.5, sprintf("%.4f", mean(dens_range)), pos = 4, cex = 0.75)
  text(0.50, 0.15, sprintf("%.4f", dens_range[1]), pos = 4, cex = 0.75)

  dev.off()

  cat("  Saved:", filename, "\n")
}

# Function to create ternary plot with Dirichlet density overlay
plot_ternary_with_dirichlet_overlay <- function(n, alpha, plot_title, filename) {
  require(extraDistr)
  require(Ternary)
  require(brms)

  # Enumerate all possible multinomial outcomes
  dat <- enumerate_multinomial_sample_space(n)

  # Calculate Dirichlet-Multinomial densities for outcomes
  dm_densities <- ddirmnom(dat, size = n, alpha = alpha)

  # Create color palette for points
  n_colors <- 10
  # Use viridis palette (color-blind friendly, matching regular ternary plots)
  color_palette <- viridisLite::viridis(n_colors, begin = 0.1, end = 0.95)

  # Get range of non-zero densities
  dens_nonzero <- dm_densities[dm_densities > 0]
  dens_min <- min(dens_nonzero)
  dens_max <- max(dens_nonzero)

  # Create bins for point colors
  breaks <- seq(dens_min - 1e-10, dens_max + 1e-10, length.out = n_colors + 1)

  # Assign colors to points based on DM density
  point_colors <- rep(rgb(0.9, 0.9, 0.9, 0.3), length(dm_densities))
  for (i in seq_along(dm_densities)) {
    if (dm_densities[i] > 0) {
      bin_idx <- findInterval(dm_densities[i], breaks, all.inside = TRUE)
      bin_idx <- min(bin_idx, n_colors)
      point_colors[i] <- color_palette[bin_idx]
    }
  }

  # Create plot
  png(file.path(output_dir, filename),
      width = 1100, height = 800, res = 150)

  # Set up layout for plot + legend
  layout(matrix(c(1, 2), nrow = 1), widths = c(3.0, 1.5))

  # Main ternary plot
  par(mar = c(1, 1, 3, 1))
  TernaryPlot(alab = "Variant A proportion →",
              blab = "← Variant B proportion ",
              clab = "Variant C proportion →",
              grid.lines = 5,
              point = "right", lab.cex = 0.8, grid.minor.lines = 0,
              grid.lty = "solid", col = rgb(0.9, 0.9, 0.9),
              grid.col = "white",
              axis.col = rgb(0.6, 0.6, 0.6),
              ticks.col = rgb(0.6, 0.6, 0.6),
              axis.rotate = FALSE,
              padding = 0.08)

  # Add Dirichlet density as colored blanket using dense grid
  # Create a fine grid over the simplex
  grid_res <- 150
  grid_coords <- expand.grid(
    a = seq(0, 1, length.out = grid_res),
    b = seq(0, 1, length.out = grid_res)
  )
  # Filter to valid simplex coordinates
  grid_coords$c <- 1 - grid_coords$a - grid_coords$b
  grid_coords <- grid_coords[grid_coords$c >= 0 & grid_coords$c <= 1, ]

  # Calculate Dirichlet densities
  grid_coords$density <- brms::ddirichlet(
    cbind(grid_coords$a, grid_coords$b, grid_coords$c),
    alpha
  )

  # Create color mapping using mako palette
  n_dens_colors <- 100
  dens_colors <- viridisLite::mako(n_dens_colors, alpha = 0.8, begin = 0, end = 1, direction = -1)
  dens_min <- min(grid_coords$density)
  dens_max <- max(grid_coords$density)
  dens_breaks <- seq(dens_min, dens_max, length.out = n_dens_colors + 1)

  grid_coords$color_idx <- cut(grid_coords$density,
                               breaks = dens_breaks,
                               labels = FALSE,
                               include.lowest = TRUE)
  grid_coords$color_idx[is.na(grid_coords$color_idx)] <- 1
  grid_coords$color <- dens_colors[grid_coords$color_idx]

  # Add background density as fine points to create blanket effect
  AddToTernary(graphics::points,
               cbind(grid_coords$a, grid_coords$b, grid_coords$c),
               pch = 15,  # Filled square
               col = grid_coords$color,
               cex = 0.5)

  # Convert counts to proportions for plotting
  dat_props <- dat / n

  # Add points for discrete multinomial outcomes
  AddToTernary(graphics::points,
               dat_props,
               pch = 21,  # Circle with border
               cex = 1.8,
               col = "black",  # Black border for visibility
               bg = point_colors,  # Fill color
               lwd = 1)

  # Add title
  mtext(plot_title, side = 3, line = 1, cex = 1.2, font = 2)

  # Add legend
  par(mar = c(2, 2, 3, 4))
  plot.new()

  # Legend for background
  text(0.5, 0.92, "Background:", cex = 0.85, font = 2)
  text(0.5, 0.85, "Dirichlet density", cex = 0.75)

  # Legend for points
  text(0.5, 0.70, "Points:", cex = 0.85, font = 2)
  text(0.5, 0.63, "DM count outcomes", cex = 0.75)

  # Color scale for DM points
  dm_dens_range <- range(dm_densities[dm_densities > 0])
  legend_y <- seq(0.40, 0.55, length.out = 6)

  for (i in 1:(length(legend_y)-1)) {
    idx <- round(seq(1, n_colors, length.out = length(legend_y)-1))[i]
    rect(0.25, legend_y[i], 0.40,
         legend_y[i+1],
         col = color_palette[idx], border = NA)
  }

  text(0.22, 0.56, sprintf("%.3f", dm_dens_range[2]), pos = 2, cex = 0.65)
  text(0.22, 0.39, sprintf("%.3f", dm_dens_range[1]), pos = 2, cex = 0.65)

  dev.off()

  cat("  Saved:", filename, "\n")
}

# Generate plots for main symmetric and asymmetric alphas
cat("BASIC TERNARY PLOTS:\n")
cat("Symmetric alpha (10, 10, 10) from simulations:\n")
for (n_val in N_values) {
  plot_title <- sprintf("Dirichlet-Multinomial: N=%d, α=(10, 10, 10)", n_val)
  filename <- sprintf("ternary_symmetric_correct_N%d.png", n_val)
  plot_ternary_with_density_legend(n_val, alpha_symmetric, plot_title, filename)
}
cat("\n")

cat("Asymmetric alpha (2, 3, 15) from simulations:\n")
for (n_val in N_values) {
  plot_title <- sprintf("Dirichlet-Multinomial: N=%d, α=(2, 3, 15)", n_val)
  filename <- sprintf("ternary_asymmetric_correct_N%d.png", n_val)
  plot_ternary_with_density_legend(n_val, alpha_asymmetric, plot_title, filename)
}
cat("\n")

# ============================================================================
# PART 1b: Ternary plots with Dirichlet density overlay
# ============================================================================

cat("Creating ternary plots WITH Dirichlet overlay...\n\n")

# Generate overlay plots for symmetric alpha (correct)
cat("OVERLAY TERNARY PLOTS:\n")
cat("Symmetric alpha (10, 10, 10) with Dirichlet overlay:\n")
for (n_val in N_values) {
  plot_title <- sprintf("DM + Dirichlet Overlay: N=%d, α=(10, 10, 10)", n_val)
  filename <- sprintf("ternary_overlay_symmetric_correct_N%d.png", n_val)
  plot_ternary_with_dirichlet_overlay(n_val, alpha_symmetric, plot_title, filename)
}
cat("\n")

# Generate overlay plots for asymmetric alpha (correct)
cat("Asymmetric alpha (2, 3, 15) with Dirichlet overlay:\n")
for (n_val in N_values) {
  plot_title <- sprintf("DM + Dirichlet Overlay: N=%d, α=(2, 3, 15)", n_val)
  filename <- sprintf("ternary_overlay_asymmetric_correct_N%d.png", n_val)
  plot_ternary_with_dirichlet_overlay(n_val, alpha_asymmetric, plot_title, filename)
}
cat("\n")

# ============================================================================
# PART 1c: Side-by-side comparison: Correct vs Misspecified alpha
# ============================================================================

# Function to create side-by-side ternary comparison
plot_ternary_comparison <- function(n, alpha_correct, alpha_wrong,
                                   correct_label, wrong_label, filename) {
  require(extraDistr)
  require(Ternary)
  require(brms)

  # Helper function to plot one ternary with overlay
  plot_one_ternary <- function(alpha, title_text) {
    # Enumerate outcomes
    dat <- enumerate_multinomial_sample_space(n)
    dm_densities <- ddirmnom(dat, size = n, alpha = alpha)

    # Create color palette for points
    n_colors <- 10
    color_palette <- viridisLite::viridis(n_colors, begin = 0.1, end = 0.95)

    # Get range and assign colors
    dens_nonzero <- dm_densities[dm_densities > 0]
    dens_min <- min(dens_nonzero)
    dens_max <- max(dens_nonzero)
    breaks <- seq(dens_min - 1e-10, dens_max + 1e-10, length.out = n_colors + 1)

    point_colors <- rep(rgb(0.9, 0.9, 0.9, 0.3), length(dm_densities))
    for (i in seq_along(dm_densities)) {
      if (dm_densities[i] > 0) {
        bin_idx <- findInterval(dm_densities[i], breaks, all.inside = TRUE)
        bin_idx <- min(bin_idx, n_colors)
        point_colors[i] <- color_palette[bin_idx]
      }
    }

    # Plot ternary
    TernaryPlot(alab = "Variant A →",
                blab = "← Variant B",
                clab = "Variant C →",
                grid.lines = 5,
                point = "right", lab.cex = 0.7, grid.minor.lines = 0,
                grid.lty = "solid", col = rgb(0.9, 0.9, 0.9),
                grid.col = "white",
                axis.col = rgb(0.6, 0.6, 0.6),
                ticks.col = rgb(0.6, 0.6, 0.6),
                axis.rotate = FALSE,
                padding = 0.08)

    # Add Dirichlet density blanket
    grid_res <- 100
    grid_coords <- expand.grid(
      a = seq(0, 1, length.out = grid_res),
      b = seq(0, 1, length.out = grid_res)
    )
    grid_coords$c <- 1 - grid_coords$a - grid_coords$b
    # Filter to valid simplex coordinates with stricter bounds to avoid edge issues
    grid_coords <- grid_coords[grid_coords$c >= 0.001 & grid_coords$c <= 0.999 &
                               grid_coords$a >= 0.001 & grid_coords$a <= 0.999 &
                               grid_coords$b >= 0.001 & grid_coords$b <= 0.999, ]

    # Renormalize to ensure sum = 1 exactly
    row_sums <- grid_coords$a + grid_coords$b + grid_coords$c
    grid_coords$a <- grid_coords$a / row_sums
    grid_coords$b <- grid_coords$b / row_sums
    grid_coords$c <- grid_coords$c / row_sums

    grid_coords$density <- brms::ddirichlet(
      cbind(grid_coords$a, grid_coords$b, grid_coords$c),
      alpha
    )

    # Handle extreme/infinite density values
    grid_coords$density[is.infinite(grid_coords$density) | is.na(grid_coords$density)] <- NA

    n_dens_colors <- 100
    dens_colors <- viridisLite::mako(n_dens_colors, alpha = 0.7, begin = 0, end = 1, direction = -1)
    # Use robust min/max that ignores NA values
    dens_min_bg <- min(grid_coords$density, na.rm = TRUE)
    dens_max_bg <- max(grid_coords$density, na.rm = TRUE)

    # Ensure finite values
    if (!is.finite(dens_min_bg) || !is.finite(dens_max_bg)) {
      dens_min_bg <- 0
      dens_max_bg <- 1
    }

    dens_breaks <- seq(dens_min_bg, dens_max_bg, length.out = n_dens_colors + 1)

    grid_coords$color_idx <- cut(grid_coords$density,
                                 breaks = dens_breaks,
                                 labels = FALSE,
                                 include.lowest = TRUE)
    grid_coords$color_idx[is.na(grid_coords$color_idx)] <- 1
    grid_coords$color <- dens_colors[grid_coords$color_idx]

    AddToTernary(graphics::points,
                 cbind(grid_coords$a, grid_coords$b, grid_coords$c),
                 pch = 15,
                 col = grid_coords$color,
                 cex = 0.4)

    # Add DM outcome points
    dat_props <- dat / n
    AddToTernary(graphics::points,
                 dat_props,
                 pch = 21,
                 cex = 2.2,
                 col = "white",
                 bg = point_colors,
                 lwd = 1.5)

    # Add title
    mtext(title_text, side = 3, line = 0.5, cex = 0.9, font = 2)
  }

  # Create plot
  png(file.path(output_dir, filename),
      width = 1400, height = 700, res = 150)

  par(mfrow = c(1, 2), mar = c(1, 1, 3, 1))

  # Left: Correct alpha
  plot_one_ternary(alpha_correct,
                   sprintf("Correct: α=(%s)", paste(alpha_correct, collapse=", ")))

  # Right: Misspecified alpha
  plot_one_ternary(alpha_wrong,
                   sprintf("Misspecified: α=(%s)", paste(alpha_wrong, collapse=", ")))

  # Add overall title
  mtext(sprintf("N=%d: %s vs %s", n, correct_label, wrong_label),
        side = 3, line = -1.5, outer = TRUE, cex = 1.2, font = 2)

  dev.off()

  cat("  Saved:", filename, "\n")
}

cat("Creating side-by-side comparison plots...\n\n")

# Define all simulation scenarios
symmetric_scenarios <- list(
  list(name = "scenario1_mislocation_symmetric", alpha_correct = c(10, 10, 10), alpha_wrong = c(5, 10, 15)),
  list(name = "scenario2_moderate_dispersed", alpha_correct = c(10, 10, 10), alpha_wrong = c(2, 2, 2)),
  list(name = "scenario3_extreme_dispersed", alpha_correct = c(10, 10, 10), alpha_wrong = c(0.5, 0.5, 0.5)),
  list(name = "scenario4_moderate_narrow", alpha_correct = c(10, 10, 10), alpha_wrong = c(50, 50, 50)),
  list(name = "scenario5_extreme_narrow", alpha_correct = c(10, 10, 10), alpha_wrong = c(200, 200, 200)),
  list(name = "scenario6_ultra_dispersed", alpha_correct = c(10, 10, 10), alpha_wrong = c(0.1, 0.1, 0.1)),
  list(name = "scenario7_ultra_narrow", alpha_correct = c(10, 10, 10), alpha_wrong = c(1000, 1000, 1000))
)

asymmetric_scenarios <- list(
  list(name = "scenario1_misspecified", alpha_correct = c(2, 3, 15), alpha_wrong = c(15, 3, 2)),
  list(name = "scenario2_narrow_calibrated", alpha_correct = c(2, 3, 15), alpha_wrong = c(200, 300, 1500)),
  list(name = "scenario3_dispersed_calibrated", alpha_correct = c(2, 3, 15), alpha_wrong = c(0.4, 0.6, 3)),
  list(name = "scenario4_less_misspecified", alpha_correct = c(2, 3, 15), alpha_wrong = c(4, 5, 11))
)

# Generate comparison plots for all symmetric scenarios
cat("SYMMETRIC SCENARIOS:\n")
for (scenario in symmetric_scenarios) {
  cat(sprintf("  %s:\n", scenario$name))
  for (n_val in N_values) {
    filename <- sprintf("comparison_symmetric_%s_N%d.png", scenario$name, n_val)
    plot_ternary_comparison(n_val, scenario$alpha_correct, scenario$alpha_wrong,
                           paste0("α=(", paste(scenario$alpha_correct, collapse=", "), ")"),
                           paste0("α=(", paste(scenario$alpha_wrong, collapse=", "), ")"),
                           filename)
  }
}
cat("\n")

# Generate comparison plots for all asymmetric scenarios
cat("ASYMMETRIC SCENARIOS:\n")
for (scenario in asymmetric_scenarios) {
  cat(sprintf("  %s:\n", scenario$name))
  for (n_val in N_values) {
    filename <- sprintf("comparison_asymmetric_%s_N%d.png", scenario$name, n_val)
    plot_ternary_comparison(n_val, scenario$alpha_correct, scenario$alpha_wrong,
                           paste0("α=(", paste(scenario$alpha_correct, collapse=", "), ")"),
                           paste0("α=(", paste(scenario$alpha_wrong, collapse=", "), ")"),
                           filename)
  }
}
cat("\n")

# ============================================================================
# PART 1d: Three-panel comparison: Dispersed vs Correct vs Narrow
# ============================================================================

# Function to create three-panel ternary comparison
plot_ternary_dispersion_comparison <- function(n, alpha_dispersed, alpha_correct, alpha_narrow,
                                              filename) {
  require(extraDistr)
  require(Ternary)
  require(brms)

  # Helper function to plot one ternary with overlay
  plot_one_ternary_small <- function(alpha, title_text) {
    # Enumerate outcomes
    dat <- enumerate_multinomial_sample_space(n)
    dm_densities <- ddirmnom(dat, size = n, alpha = alpha)

    # Create color palette for points
    n_colors <- 10
    color_palette <- viridisLite::viridis(n_colors, begin = 0.1, end = 0.95)

    # Get range and assign colors
    dens_nonzero <- dm_densities[dm_densities > 0]
    dens_min <- min(dens_nonzero)
    dens_max <- max(dens_nonzero)
    breaks <- seq(dens_min - 1e-10, dens_max + 1e-10, length.out = n_colors + 1)

    point_colors <- rep(rgb(0.9, 0.9, 0.9, 0.3), length(dm_densities))
    for (i in seq_along(dm_densities)) {
      if (dm_densities[i] > 0) {
        bin_idx <- findInterval(dm_densities[i], breaks, all.inside = TRUE)
        bin_idx <- min(bin_idx, n_colors)
        point_colors[i] <- color_palette[bin_idx]
      }
    }

    # Plot ternary
    TernaryPlot(alab = "A →",
                blab = "← B",
                clab = "C →",
                grid.lines = 5,
                point = "right", lab.cex = 0.6, grid.minor.lines = 0,
                grid.lty = "solid", col = rgb(0.9, 0.9, 0.9),
                grid.col = "white",
                axis.col = rgb(0.6, 0.6, 0.6),
                ticks.col = rgb(0.6, 0.6, 0.6),
                axis.rotate = FALSE,
                padding = 0.08)

    # Add Dirichlet density blanket (matching style from correct vs wrong plots)
    grid_res <- 150
    grid_coords <- expand.grid(
      a = seq(0, 1, length.out = grid_res),
      b = seq(0, 1, length.out = grid_res)
    )
    grid_coords$c <- 1 - grid_coords$a - grid_coords$b
    # Filter to valid simplex coordinates with stricter bounds to avoid edge issues
    grid_coords <- grid_coords[grid_coords$c >= 0.001 & grid_coords$c <= 0.999 &
                               grid_coords$a >= 0.001 & grid_coords$a <= 0.999 &
                               grid_coords$b >= 0.001 & grid_coords$b <= 0.999, ]

    # Renormalize to ensure sum = 1 exactly
    row_sums <- grid_coords$a + grid_coords$b + grid_coords$c
    grid_coords$a <- grid_coords$a / row_sums
    grid_coords$b <- grid_coords$b / row_sums
    grid_coords$c <- grid_coords$c / row_sums

    grid_coords$density <- brms::ddirichlet(
      cbind(grid_coords$a, grid_coords$b, grid_coords$c),
      alpha
    )

    # Handle extreme/infinite density values
    grid_coords$density[is.infinite(grid_coords$density) | is.na(grid_coords$density)] <- NA

    n_dens_colors <- 100
    dens_colors <- viridisLite::mako(n_dens_colors, alpha = 0.8, begin = 0, end = 1, direction = -1)
    # Use robust min/max that ignores NA values
    dens_min_bg <- min(grid_coords$density, na.rm = TRUE)
    dens_max_bg <- max(grid_coords$density, na.rm = TRUE)

    # Ensure finite values
    if (!is.finite(dens_min_bg) || !is.finite(dens_max_bg)) {
      dens_min_bg <- 0
      dens_max_bg <- 1
    }

    dens_breaks <- seq(dens_min_bg, dens_max_bg, length.out = n_dens_colors + 1)

    grid_coords$color_idx <- cut(grid_coords$density,
                                 breaks = dens_breaks,
                                 labels = FALSE,
                                 include.lowest = TRUE)
    grid_coords$color_idx[is.na(grid_coords$color_idx)] <- 1
    grid_coords$color <- dens_colors[grid_coords$color_idx]

    AddToTernary(graphics::points,
                 cbind(grid_coords$a, grid_coords$b, grid_coords$c),
                 pch = 15,
                 col = grid_coords$color,
                 cex = 0.5)

    # Add DM outcome points
    dat_props <- dat / n
    AddToTernary(graphics::points,
                 dat_props,
                 pch = 21,
                 cex = 2.0,
                 col = "white",
                 bg = point_colors,
                 lwd = 1.5)

    # Add title
    mtext(title_text, side = 3, line = 0.3, cex = 0.75, font = 2)
  }

  # Create plot
  png(file.path(output_dir, filename),
      width = 1800, height = 600, res = 150)

  par(mfrow = c(1, 3), mar = c(1, 1, 2.5, 1))

  # Left: Dispersed (smaller alpha = more spread)
  plot_one_ternary_small(alpha_dispersed,
                         sprintf("Dispersed\nα=(%s)", paste(alpha_dispersed, collapse=", ")))

  # Center: Correct
  plot_one_ternary_small(alpha_correct,
                         sprintf("Correct\nα=(%s)", paste(alpha_correct, collapse=", ")))

  # Right: Narrow (larger alpha = more concentrated)
  plot_one_ternary_small(alpha_narrow,
                         sprintf("Narrow\nα=(%s)", paste(alpha_narrow, collapse=", ")))

  # Add overall title
  mtext(sprintf("N=%d: Dispersion Comparison", n),
        side = 3, line = -1.5, outer = TRUE, cex = 1.2, font = 2)

  dev.off()

  cat("  Saved:", filename, "\n")
}

cat("Creating dispersion comparison plots...\n\n")

# Define dispersion comparison sets for each scenario type
symmetric_dispersion_sets <- list(
  list(dispersed = c(2, 2, 2), correct = c(10, 10, 10), narrow = c(50, 50, 50)),
  list(dispersed = c(0.5, 0.5, 0.5), correct = c(10, 10, 10), narrow = c(200, 200, 200)),
  list(dispersed = c(0.1, 0.1, 0.1), correct = c(10, 10, 10), narrow = c(1000, 1000, 1000))
)

asymmetric_dispersion_sets <- list(
  list(dispersed = c(0.4, 0.6, 3), correct = c(2, 3, 15), narrow = c(20, 30, 150)),
  list(dispersed = c(1, 1.5, 7.5), correct = c(2, 3, 15), narrow = c(10, 15, 75))
)

# Generate dispersion comparison plots for symmetric cases
cat("SYMMETRIC DISPERSION COMPARISONS:\n")
for (i in seq_along(symmetric_dispersion_sets)) {
  dset <- symmetric_dispersion_sets[[i]]
  cat(sprintf("  Set %d: Dispersed=(%s), Correct=(%s), Narrow=(%s)\n",
              i, paste(dset$dispersed, collapse=", "),
              paste(dset$correct, collapse=", "),
              paste(dset$narrow, collapse=", ")))
  for (n_val in N_values) {
    filename <- sprintf("dispersion_symmetric_set%d_N%d.png", i, n_val)
    plot_ternary_dispersion_comparison(n_val,
                                      dset$dispersed,
                                      dset$correct,
                                      dset$narrow,
                                      filename)
  }
}
cat("\n")

# Generate dispersion comparison plots for asymmetric cases
cat("ASYMMETRIC DISPERSION COMPARISONS:\n")
for (i in seq_along(asymmetric_dispersion_sets)) {
  dset <- asymmetric_dispersion_sets[[i]]
  cat(sprintf("  Set %d: Dispersed=(%s), Correct=(%s), Narrow=(%s)\n",
              i, paste(dset$dispersed, collapse=", "),
              paste(dset$correct, collapse=", "),
              paste(dset$narrow, collapse=", ")))
  for (n_val in N_values) {
    filename <- sprintf("dispersion_asymmetric_set%d_N%d.png", i, n_val)
    plot_ternary_dispersion_comparison(n_val,
                                      dset$dispersed,
                                      dset$correct,
                                      dset$narrow,
                                      filename)
  }
}
cat("\n")

# ============================================================================
# PART 2: Additional plots for K=3 and K=5 dimensions
# ============================================================================

# Function to create all plot types for a given alpha and N
create_additional_plots <- function(n_val, alpha, alpha_label, n_samples = 100) {
  K <- length(alpha)

  cat(sprintf("Creating additional plots for N=%d, α=%s...\n", n_val, alpha_label))

  # Draw samples from Dirichlet-Multinomial
  set.seed(42 + n_val)  # Different seed for each N
  samples <- rdirmnom(n = n_samples, size = n_val, alpha = alpha)

  # Convert to data frame for plotting
  samples_df <- as.data.frame(samples)
  colnames(samples_df) <- paste0("Category_", 1:K)
  samples_df$sample_id <- 1:n_samples

  # Reshape to long format
  samples_long <- samples_df %>%
    pivot_longer(cols = starts_with("Category"),
                 names_to = "Category",
                 values_to = "Count")

  # Expected proportions from Dirichlet
  expected_props <- alpha / sum(alpha)
  expected_counts <- expected_props * n_val

  # ============================================================================
  # Plot 1: Density Heatmap showing probability of all outcomes
  # ============================================================================

  if (K == 3) {
    # For K=3, create a 2D grid (Cat1 vs Cat2, Cat3 is implied)
    outcome_space <- enumerate_multinomial_sample_space(n_val)
    densities <- ddirmnom(outcome_space, size = n_val, alpha = alpha)

    density_df <- data.frame(
      Cat1 = outcome_space[, 1],
      Cat2 = outcome_space[, 2],
      Cat3 = outcome_space[, 3],
      Density = densities
    )

    p1 <- ggplot(density_df, aes(x = Cat1, y = Cat2, fill = Density)) +
      geom_tile(color = "white", linewidth = 0.3) +
      scale_fill_viridis_c(option = "plasma", begin = 0.1, end = 0.95) +
      labs(
        title = "Outcome Density Heatmap",
        subtitle = paste0("N=", n_val, ", α=(", paste(alpha, collapse=", "),
                         "). Each cell = P(outcome) under Dirichlet-Multinomial"),
        x = "Category 1 Count",
        y = "Category 2 Count",
        fill = "Probability"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", size = 14),
        panel.grid = element_blank()
      ) +
      coord_equal()
  } else {
    # For K>3, show outcomes as rows sorted by density
    outcome_space <- matrix(0, nrow = 0, ncol = K)
    # Generate all possible outcomes (simplified for K>3)
    samples_for_heatmap <- samples[1:min(50, nrow(samples)), ]
    densities_for_heatmap <- apply(samples_for_heatmap, 1, function(x) {
      ddirmnom(matrix(x, nrow = 1), size = n_val, alpha = alpha)
    })

    # Create long format for heatmap
    heatmap_df <- data.frame(samples_for_heatmap)
    colnames(heatmap_df) <- paste0("Category_", 1:K)
    heatmap_df$Density <- densities_for_heatmap
    heatmap_df$Outcome_ID <- 1:nrow(heatmap_df)
    heatmap_df <- heatmap_df %>% arrange(desc(Density))
    heatmap_df$Outcome_ID <- 1:nrow(heatmap_df)

    heatmap_long <- heatmap_df %>%
      pivot_longer(cols = starts_with("Category"),
                   names_to = "Category",
                   values_to = "Count")

    p1 <- ggplot(heatmap_long,
                 aes(x = Category, y = factor(Outcome_ID), fill = Count)) +
      geom_tile(color = "white", linewidth = 0.2) +
      scale_fill_viridis_c(option = "plasma", begin = 0.1, end = 0.95) +
      labs(
        title = "Outcome Density (Top 50 most likely)",
        subtitle = paste0("N=", n_val, ", α=(", paste(alpha, collapse=", "), ")"),
        x = "Category",
        y = "Outcome (sorted by density)",
        fill = "Count"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", size = 14),
        panel.grid = element_blank(),
        axis.text.y = element_blank()
      )
  }

  filename <- sprintf("heatmap_density_%s_N%d.png", alpha_label, n_val)
  ggsave(file.path(output_dir, filename),
         p1, width = 10, height = 8, dpi = 300)

  # ============================================================================
  # Plot 2: Distribution of counts per category (violin plots)
  # ============================================================================

  samples_long <- samples_long %>%
    mutate(Category = factor(Category,
                             levels = paste0("Category_", 1:K)))

  # Add expected values for reference
  expected_df <- data.frame(
    Category = factor(paste0("Category_", 1:K),
                     levels = paste0("Category_", 1:K)),
    Expected = expected_counts
  )

  p2 <- ggplot(samples_long, aes(x = Category, y = Count)) +
    geom_violin(fill = "steelblue", alpha = 0.6) +
    geom_boxplot(width = 0.2, alpha = 0.8, outlier.alpha = 0.3) +
    geom_point(data = expected_df,
               aes(y = Expected),
               color = "red", size = 4, shape = 18) +
    geom_hline(data = expected_df,
               aes(yintercept = Expected),
               color = "red", linetype = "dashed", alpha = 0.5) +
    labs(
      title = "Distribution of Counts per Category",
      subtitle = paste0("Dirichlet-Multinomial: N=", n_val, ", α=(",
                       paste(alpha, collapse=", "), "). Red diamonds = expected counts"),
      x = "Category",
      y = "Count"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  filename <- sprintf("distributions_%s_N%d.png", alpha_label, n_val)
  ggsave(file.path(output_dir, filename),
         p2, width = 10, height = 6, dpi = 300)

  # ============================================================================
  # Plot 3: Correlation heatmap between categories
  # ============================================================================

  # Calculate correlation matrix
  cor_matrix <- cor(samples[, 1:K])
  colnames(cor_matrix) <- paste0("Cat_", 1:K)
  rownames(cor_matrix) <- paste0("Cat_", 1:K)

  # Convert to long format for ggplot
  cor_long <- cor_matrix %>%
    as.data.frame() %>%
    mutate(Category1 = rownames(.)) %>%
    pivot_longer(cols = -Category1,
                 names_to = "Category2",
                 values_to = "Correlation")

  p3 <- ggplot(cor_long,
               aes(x = Category1, y = Category2, fill = Correlation)) +
    geom_tile(color = "white") +
    geom_text(aes(label = round(Correlation, 2)),
              color = "white", size = 4, fontface = "bold") +
    scale_fill_gradient2(
      low = "blue", mid = "white", high = "red",
      midpoint = 0, limits = c(-1, 1)
    ) +
    labs(
      title = "Correlation Between Category Counts",
      subtitle = paste0("Negative correlations expected (counts constrained to sum to N=", n_val, ")"),
      x = NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    coord_fixed()

  filename <- sprintf("correlations_%s_N%d.png", alpha_label, n_val)
  ggsave(file.path(output_dir, filename),
         p3, width = 8, height = 7, dpi = 300)

  # ============================================================================
  # Plot 4: Empirical vs Expected Proportions
  # ============================================================================

  empirical_props <- samples_long %>%
    group_by(Category) %>%
    summarise(
      Mean_Count = mean(Count),
      .groups = "drop"
    ) %>%
    mutate(
      Empirical_Proportion = Mean_Count / n_val,
      Expected_Proportion = expected_props,
      Category_Num = 1:K
    )

  p4 <- ggplot(empirical_props,
               aes(x = Category_Num)) +
    geom_bar(aes(y = Empirical_Proportion, fill = "Empirical"),
             stat = "identity", alpha = 0.7, width = 0.7) +
    geom_point(aes(y = Expected_Proportion, color = "Expected"),
               size = 4, shape = 18) +
    geom_line(aes(y = Expected_Proportion, color = "Expected", group = 1),
              linewidth = 1, linetype = "dashed") +
    scale_fill_manual(values = c("Empirical" = "steelblue")) +
    scale_color_manual(values = c("Expected" = "red")) +
    scale_x_continuous(breaks = 1:K, labels = paste0("Cat_", 1:K)) +
    labs(
      title = "Empirical vs Expected Proportions",
      subtitle = paste0("Based on ", n_samples, " samples from Dirichlet-Multinomial"),
      x = "Category",
      y = "Proportion",
      fill = NULL,
      color = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      legend.position = "bottom"
    )

  filename <- sprintf("proportions_%s_N%d.png", alpha_label, n_val)
  ggsave(file.path(output_dir, filename),
         p4, width = 10, height = 6, dpi = 300)

  cat(sprintf("  Saved 4 plots for N=%d, α=%s\n", n_val, alpha_label))
}

# Generate additional plots for K=3 symmetric and asymmetric alpha (correct only)
cat("\nADDITIONAL PLOTS (heatmaps, distributions, correlations, proportions):\n")
cat("Symmetric alpha (10, 10, 10):\n")
for (n_val in N_values) {
  create_additional_plots(n_val, alpha_symmetric, "symmetric_correct", n_samples)
}
cat("\n")

cat("Asymmetric alpha (2, 3, 15):\n")
for (n_val in N_values) {
  create_additional_plots(n_val, alpha_asymmetric, "asymmetric_correct", n_samples)
}
cat("\n")

# ============================================================================
# PART 3: K=5 example
# ============================================================================

# Use one of the N values for K=5 visualization
N <- 10
alpha_5d <- c(2, 5, 10, 20, 40)  # Increasing probabilities

create_additional_plots(N, alpha_5d, "K5", n_samples)
cat("\n")

# ============================================================================
# SUMMARY
# ============================================================================

cat("=========================================\n")
cat("SUMMARY\n")
cat("=========================================\n\n")
cat("Configuration:\n")
cat("  Alpha values from simulations:\n")
cat("    - Symmetric correct (K=3): (10, 10, 10)\n")
cat("    - Asymmetric correct (K=3): (2, 3, 15)\n")
cat("    - 7 symmetric scenarios with different misspecifications\n")
cat("    - 4 asymmetric scenarios with different misspecifications\n")
cat("    - 3 symmetric dispersion sets\n")
cat("    - 2 asymmetric dispersion sets\n")
cat("    - K=5 example: (2, 5, 10, 20, 40)\n")
cat("  N values:", paste(N_values, collapse=", "), "\n")
cat("  Samples drawn:", n_samples, "\n\n")

cat("Files created in", output_dir, ":\n")
cat("  Basic ternary plots (2 alpha types × ", length(N_values), " N values):\n")
cat("    - ternary_symmetric_correct_N*.png\n")
cat("    - ternary_asymmetric_correct_N*.png\n")
cat("  Overlay ternary plots (2 alpha types × ", length(N_values), " N values):\n")
cat("    - ternary_overlay_symmetric_correct_N*.png\n")
cat("    - ternary_overlay_asymmetric_correct_N*.png\n")
cat("  Comparison plots (11 scenarios × ", length(N_values), " N values):\n")
cat("    - comparison_symmetric_scenario*_N*.png (7 scenarios)\n")
cat("    - comparison_asymmetric_scenario*_N*.png (4 scenarios)\n")
cat("  Dispersion comparison plots (5 sets × ", length(N_values), " N values):\n")
cat("    - dispersion_symmetric_set*_N*.png (3 sets)\n")
cat("    - dispersion_asymmetric_set*_N*.png (2 sets)\n")
cat("  Additional plots (2 alpha types × ", length(N_values), " N values × 4 plot types):\n")
cat("    - heatmap_density_*_correct_N*.png\n")
cat("    - distributions_*_correct_N*.png\n")
cat("    - correlations_*_correct_N*.png\n")
cat("    - proportions_*_correct_N*.png\n")
cat("  K=5 example (N=10, 4 plot types):\n")
cat("    - heatmap_density_K5_N10.png\n")
cat("    - distributions_K5_N10.png\n")
cat("    - correlations_K5_N10.png\n")
cat("    - proportions_K5_N10.png\n\n")

n_basic <- length(N_values) * 2
n_overlay <- length(N_values) * 2
n_comparison <- length(N_values) * (7 + 4)  # 11 scenarios
n_dispersion <- length(N_values) * (3 + 2)  # 5 dispersion sets
n_additional <- length(N_values) * 2 * 4    # 2 alpha types × 4 plot types
n_k5 <- 4

total_plots <- n_basic + n_overlay + n_comparison + n_dispersion + n_additional + n_k5
cat("Total plots generated:", total_plots, "\n")
cat("  Basic ternaries:", n_basic, "\n")
cat("  Overlay ternaries:", n_overlay, "\n")
cat("  Comparison plots:", n_comparison, "\n")
cat("  Dispersion plots:", n_dispersion, "\n")
cat("  Additional plots:", n_additional, "\n")
cat("  K=5 plots:", n_k5, "\n\n")

cat("Key insights:\n")
cat("  - Ternary plots show probability mass with color-coded densities (viridis palette)\n")
cat("  - Overlay plots show Dirichlet density (continuous, background) + DM outcomes (discrete, points)\n")
cat("  - Comparison plots show correct vs misspecified alpha side-by-side\n")
cat("  - Dispersion plots show effect of alpha magnitude: dispersed (small α) vs narrow (large α)\n")
cat("  - Symmetric alpha (10,10,10) shows uniform distribution\n")
cat("  - Asymmetric alpha (2,3,15) shows concentration toward variant C\n")
cat("  - Small N values show discrete nature and impact of misspecification\n")
cat("  - Dispersed priors (small α) = more uncertainty, spread across simplex\n")
cat("  - Narrow priors (large α) = more certainty, concentrated in center/mode\n")
cat("  - Density heatmaps show probability of each outcome (Cat1 vs Cat2 grid for K=3)\n")
cat("  - Violin plots reveal the variability in each category\n")
cat("  - Correlations show the constraint structure (sum to N)\n")
cat("  - Proportion comparisons validate sampling accuracy\n")
cat("  - Color-blind friendly viridis palette used throughout\n\n")

cat("To modify:\n")
cat("  - Change N_values to see different sample sizes\n")
cat("  - Adjust alpha_symmetric or alpha_asymmetric for different priors\n")
cat("  - Increase n_samples for more stable empirical estimates\n\n")

cat("Done!\n")
